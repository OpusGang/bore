// This file is part of bore.
// Copyright (C) 2024 OpusGang
//
// bore is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <stdint.h>

#include "common.h"
#include "VapourSynth4.h"
#include "VSHelper4.h"

#ifndef BORE_VERSION_MAJOR
#define BORE_VERSION_MAJOR 1
#define BORE_VERSION_MINOR 0
#define BORE_VERSION_PATCH 1
#endif

// VapourSynth plugin versions only have major and minor parts: encode minor.patch as
// minor * 100 + patch (1.0.1 -> 1.1, 1.2.0 -> 1.200), which keeps versions ordered.
#define BORE_VS_VERSION VS_MAKE_VERSION(BORE_VERSION_MAJOR, BORE_VERSION_MINOR * 100 + BORE_VERSION_PATCH)

typedef struct {
    LinearRegressionData data;
    VSNode *node;
    VSNode *weight_mask;
} VS_LinearRegressionData;

static void set_frame_props(VSFrame* frame, const double* props, const VSAPI* vsapi)
{
    VSMap* __restrict dst_props = vsapi->getFramePropertiesRW(frame);
    vsapi->mapSetFloat(dst_props, "BoreAdjustment", props[0], maAppend);
    vsapi->mapSetFloat(dst_props, "BoreCovariance", props[1], maAppend);
    vsapi->mapSetFloat(dst_props, "BoreSumSquares", props[2], maAppend);
}

// Process the border lines of one plane, innermost to outermost. If props_frame is set,
// the fit results of every line are appended to its frame properties (debug mode).
static void processBorders(const LinearRegressionData* data, int w, int h, ptrdiff_t stride, float* __restrict dstp,
    BoreLineParams* params, VSFrame* props_frame, const VSAPI* vsapi)
{
    for (int row = data->top - 1; row > -1; --row) {
        params->mask_dist = data->top - row;
        data->processRow(row, w, h, stride, dstp, params);
        if (props_frame)
            set_frame_props(props_frame, params->c1_cov11_sumsq, vsapi);
    }
    for (int row = h - data->bottom; row < h; ++row) {
        params->mask_dist = data->bottom + row - h + 1;
        data->processRow(row, w, h, stride, dstp, params);
        if (props_frame)
            set_frame_props(props_frame, params->c1_cov11_sumsq, vsapi);
    }
    for (int column = data->left - 1; column > -1; --column) {
        params->mask_dist = data->left - column;
        data->processColumn(column, w, h, stride, dstp, params);
        if (props_frame)
            set_frame_props(props_frame, params->c1_cov11_sumsq, vsapi);
    }
    for (int column = w - data->right; column < w; ++column) {
        params->mask_dist = data->right + column - w + 1;
        data->processColumn(column, w, h, stride, dstp, params);
        if (props_frame)
            set_frame_props(props_frame, params->c1_cov11_sumsq, vsapi);
    }
}

static const VSFrame *VS_CC singlePlaneGetFrame(int n, int activationReason, void *instanceData, void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) {
    (void)frameData;
    VS_LinearRegressionData *d = (VS_LinearRegressionData *)instanceData;

    if (activationReason == arInitial) {
        vsapi->requestFrameFilter(n, d->node, frameCtx);
        if (d->weight_mask)
            vsapi->requestFrameFilter(n, d->weight_mask, frameCtx);
    } else if (activationReason == arAllFramesReady) {
        const LinearRegressionData* data = &d->data;
        const VSFrame *src = vsapi->getFrameFilter(n, d->node, frameCtx);
        const VSFrame *weight_mask = NULL;

        VSFrame *dst = vsapi->copyFrame(src, core);

        const int plane = data->plane;
        ptrdiff_t stride = vsapi->getStride(dst, plane) / 4;
        int w = vsapi->getFrameWidth(src, plane);
        int h = vsapi->getFrameHeight(src, plane);
        float* __restrict dstp = (float *) vsapi->getWritePtr(dst, plane);
        double c1_cov11_sumsq[3] = { 0.0 };

        BoreScratch* scratch = boreScratchAlloc(data, w, h);
        if (!scratch) {
            vsapi->setFilterError("bore: out of memory", frameCtx);
            vsapi->freeFrame(dst);
            vsapi->freeFrame(src);
            return NULL;
        }

        BoreLineParams params = {
            .ref_line_size = data->ref_line_size,
            .sigmaS = data->sigmaS,
            .sigmaR = data->sigmaR,
            .sigmaD = data->sigmaD,
            .c1_cov11_sumsq = c1_cov11_sumsq,
            .scratch = scratch
        };

        if (d->weight_mask) {
            weight_mask = vsapi->getFrameFilter(n, d->weight_mask, frameCtx);
            params.wmaskstride = vsapi->getStride(weight_mask, 0) / 4;
            params.wmaskp = (const float*) vsapi->getReadPtr(weight_mask, 0);
        }

        processBorders(data, w, h, stride, dstp, &params, data->mode == LINREG_MODE_SINGLE_DEBUG ? dst : NULL, vsapi);

        boreScratchFree(scratch);
        vsapi->freeFrame(weight_mask);
        vsapi->freeFrame(src);

        return dst;
    }

    return NULL;
}

static const VSFrame *VS_CC multiPlaneGetFrame(int n, int activationReason, void *instanceData, void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) {
    (void)frameData;
    VS_LinearRegressionData *d = (VS_LinearRegressionData *)instanceData;

    if (activationReason == arInitial) {
        vsapi->requestFrameFilter(n, d->node, frameCtx);
    } else if (activationReason == arAllFramesReady) {
        const LinearRegressionData* data = &d->data;
        const VSFrame *src = vsapi->getFrameFilter(n, d->node, frameCtx);

        VSFrame *dst = vsapi->copyFrame(src, core);

        const int plane = data->plane;
        ptrdiff_t stride = vsapi->getStride(dst, plane) / 4;
        int w = vsapi->getFrameWidth(src, plane);
        int h = vsapi->getFrameHeight(src, plane);
        float* dstp = (float*)vsapi->getWritePtr(dst, plane);
        float* dstp1 = (float*)vsapi->getWritePtr(dst, 0);
        float* dstp2 = (float*)vsapi->getWritePtr(dst, 1);
        float* dstp3 = (float*)vsapi->getWritePtr(dst, 2);

        BoreScratch* scratch = boreScratchAlloc(data, w, h);
        if (!scratch) {
            vsapi->setFilterError("bore: out of memory", frameCtx);
            vsapi->freeFrame(dst);
            vsapi->freeFrame(src);
            return NULL;
        }

        for (int row = data->top - 1; row > -1; --row)
            processRowMLR(row, w, h, stride, dstp, dstp1, dstp2, dstp3, scratch);
        for (int row = h - data->bottom; row < h; ++row)
            processRowMLR(row, w, h, stride, dstp, dstp1, dstp2, dstp3, scratch);
        for (int column = data->left - 1; column > -1; --column)
            processColumnMLR(column, w, h, stride, dstp, dstp1, dstp2, dstp3, scratch);
        for (int column = w - data->right; column < w; ++column)
            processColumnMLR(column, w, h, stride, dstp, dstp1, dstp2, dstp3, scratch);

        boreScratchFree(scratch);
        vsapi->freeFrame(src);

        return dst;
    }

    return NULL;
}

static void VS_CC linearRegressionFree(void *instanceData, VSCore *core, const VSAPI *vsapi) {
    (void)core;
    VS_LinearRegressionData *d = (VS_LinearRegressionData *)instanceData;
    vsapi->freeNode(d->node);
    vsapi->freeNode(d->weight_mask);
    free(d);
}

static void VS_CC linearRegressionCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi) {
    LinRegMode mode = (LinRegMode)(intptr_t)userData;
    int err, w, h;
    const char *error = NULL;
    VSNode *weight_mask = NULL;
    VS_LinearRegressionData *d = NULL;

    VSNode* node = vsapi->mapGetNode(in, "clip", 0, 0);
    const VSVideoInfo *vi = vsapi->getVideoInfo(node);

    int plane = vsapi->mapGetIntSaturated(in, "plane", 0, &err);
    if (err)
        plane = 0;

    if (!vsh_isConstantVideoFormat(vi) || vi->format.sampleType != stFloat || vi->format.bitsPerSample != 32) {
        error = "bore: only constant format single float clip input is supported";
        goto fail;
    }

    if (vi->width == 0 || vi->height == 0) {
        error = "bore: only constant resolution clip input is supported";
        goto fail;
    }

    if (vi->format.numPlanes == 1 && plane > 0) {
        error = "bore: plane cannot be bigger than 0";
        goto fail;
    }

    if (plane < 0 || plane >= vi->format.numPlanes) {
        error = "bore: plane must be in [0, number of planes - 1]";
        goto fail;
    }

    if (mode == LINREG_MODE_MULTI) {
        if (vi->format.numPlanes == 1) {
            error = "bore: clip must have 3 planes";
            goto fail;
        }
        if (vi->format.colorFamily != cfRGB && (vi->format.colorFamily != cfYUV || vi->format.subSamplingW != 0 || vi->format.subSamplingH != 0)) {
            error = "bore: only RGB or YUV444 format is supported";
            goto fail;
        }
    }

    w = vi->width;
    h = vi->height;

    if (plane > 0) {
        w = w >> vi->format.subSamplingW;
        h = h >> vi->format.subSamplingH;
    }

    weight_mask = vsapi->mapGetNode(in, "weight_mask", 0, &err);
    if (err)
        weight_mask = NULL;
    else {
        const VSVideoInfo *ivi = vsapi->getVideoInfo(weight_mask);
        if (!vsh_isConstantVideoFormat(ivi) || ivi->format.sampleType != stFloat || ivi->format.bitsPerSample != 32) {
            error = "bore: only constant format single float weight_mask input is supported";
            goto fail;
        }

        if (ivi->format.numPlanes > 1) {
            error = "bore: weight_mask must be in gray format";
            goto fail;
        }

        if (ivi->width != w || ivi->height != h) {
            error = "bore: clip and weight_mask must have matching dimensions";
            goto fail;
        }

        if (ivi->width == 0 || ivi->height == 0) {
            error = "bore: only constant resolution weight_mask input is supported";
            goto fail;
        }
    }

    int top = vsapi->mapGetIntSaturated(in, "top", 0, &err);
    if (err)
        top = 0;
    else if (top < 0 || top > h / 2) {
        error = "bore: top must be in [0, height / 2]";
        goto fail;
    }

    int bottom = vsapi->mapGetIntSaturated(in, "bottom", 0, &err);
    if (err)
        bottom = 0;
    else if (bottom < 0 || bottom > h / 2) {
        error = "bore: bottom must be in [0, height / 2]";
        goto fail;
    }

    int left = vsapi->mapGetIntSaturated(in, "left", 0, &err);
    if (err)
        left = 0;
    else if (left < 0 || left > w / 2) {
        error = "bore: left must be in [0, width / 2]";
        goto fail;
    }

    int right = vsapi->mapGetIntSaturated(in, "right", 0, &err);
    if (err)
        right = 0;
    else if (right < 0 || right > w / 2) {
        error = "bore: right must be in [0, width / 2]";
        goto fail;
    }

    int ref_line_size = 0;
    if (mode == LINREG_MODE_SINGLE_LIMITED || mode == LINREG_MODE_SINGLE_WEIGHTED) {
        ref_line_size = vsapi->mapGetIntSaturated(in, "ref_line_size", 0, &err);
        if (err)
            ref_line_size = 100;
        else if (ref_line_size < 0) {
            error = "bore: ref_line_size cannot be negative";
            goto fail;
        }
    }

    double sigmaS = 0.0, sigmaR = 0.0, sigmaD = 0.0;
    if (mode == LINREG_MODE_SINGLE_WEIGHTED) {
        sigmaS = vsapi->mapGetFloat(in, "sigmaS", 0, &err);
        if (err)
            sigmaS = 50.0;
        else if (!(sigmaS > 0.0)) {
            error = "bore: sigmaS must be greater than 0";
            goto fail;
        }

        sigmaR = vsapi->mapGetFloat(in, "sigmaR", 0, &err);
        if (err)
            sigmaR = 0.5;
        else if (!(sigmaR > 0.0)) {
            error = "bore: sigmaR must be greater than 0";
            goto fail;
        }

        sigmaD = vsapi->mapGetFloat(in, "sigmaD", 0, &err);
        if (err)
            sigmaD = 1.5;
        else if (!(sigmaD > 0.0)) {
            error = "bore: sigmaD must be greater than 0";
            goto fail;
        }
    }

    d = (VS_LinearRegressionData*)calloc(1, sizeof(VS_LinearRegressionData));
    if (!d) {
        error = "bore: out of memory";
        goto fail;
    }
    d->node = node;
    d->weight_mask = weight_mask;
    d->data.mode = mode;
    d->data.plane = plane;
    d->data.top = top;
    d->data.bottom = bottom;
    d->data.left = left;
    d->data.right = right;
    d->data.ref_line_size = ref_line_size;
    d->data.sigmaS = sigmaS;
    d->data.sigmaR = sigmaR;
    d->data.sigmaD = sigmaD;

    switch (mode) {
        case LINREG_MODE_SINGLE:
            if (weight_mask) {
                d->data.processRow = &processRowSLRMasked;
                d->data.processColumn = &processColumnSLRMasked;
            } else {
                d->data.processRow = &processRowSLR;
                d->data.processColumn = &processColumnSLR;
            }
            break;
        case LINREG_MODE_SINGLE_LIMITED:
            if (weight_mask) {
                d->data.processRow = &processRowSLRRefMasked;
                d->data.processColumn = &processColumnSLRRefMasked;
            } else {
                d->data.processRow = &processRowSLRRef;
                d->data.processColumn = &processColumnSLRRef;
            }
            break;
        case LINREG_MODE_SINGLE_WEIGHTED:
            if (weight_mask) {
                d->data.processRow = &processRowWSLRMasked;
                d->data.processColumn = &processColumnWSLRMasked;
            } else {
                d->data.processRow = &processRowWSLR;
                d->data.processColumn = &processColumnWSLR;
            }
            break;
        case LINREG_MODE_SINGLE_DEBUG:
            if (weight_mask) {
                d->data.processRow = &debugRowSLRMasked;
                d->data.processColumn = &debugColumnSLRMasked;
            } else {
                d->data.processRow = &debugRowSLR;
                d->data.processColumn = &debugColumnSLR;
            }
            break;
        case LINREG_MODE_MULTI:
            break;
    }

    // MultiPlane accepts weight_mask but does not use it, so it does not request its frames.
    VSFilterDependency deps[] = {{node, rpStrictSpatial}, {weight_mask, rpStrictSpatial}};
    const int num_deps = (weight_mask && mode != LINREG_MODE_MULTI) ? 2 : 1;

    switch (mode) {
        case LINREG_MODE_SINGLE:
            vsapi->createVideoFilter(out, "SinglePlane", vi, singlePlaneGetFrame, linearRegressionFree, fmParallel, deps, num_deps, d, core);
            break;
        case LINREG_MODE_MULTI:
            vsapi->createVideoFilter(out, "MultiPlane", vi, multiPlaneGetFrame, linearRegressionFree, fmParallel, deps, num_deps, d, core);
            break;
        case LINREG_MODE_SINGLE_LIMITED:
            vsapi->createVideoFilter(out, "SinglePlaneLimited", vi, singlePlaneGetFrame, linearRegressionFree, fmParallel, deps, num_deps, d, core);
            break;
        case LINREG_MODE_SINGLE_WEIGHTED:
            vsapi->createVideoFilter(out, "SinglePlaneWeighted", vi, singlePlaneGetFrame, linearRegressionFree, fmParallel, deps, num_deps, d, core);
            break;
        case LINREG_MODE_SINGLE_DEBUG:
            vsapi->createVideoFilter(out, "SinglePlaneDebug", vi, singlePlaneGetFrame, linearRegressionFree, fmParallel, deps, num_deps, d, core);
            break;
    }
    return;

fail:
    vsapi->mapSetError(out, error);
    vsapi->freeNode(node);
    vsapi->freeNode(weight_mask);
}

VS_EXTERNAL_API(void) VapourSynthPluginInit2(VSPlugin *plugin, const VSPLUGINAPI *vspapi) {
    boreInit();
    vspapi->configPlugin("ng.opusga.bore", "bore", "bore plugin", BORE_VS_VERSION, VAPOURSYNTH_API_VERSION, 0, plugin);
    vspapi->registerFunction("SinglePlane", "clip:vnode;left:int:opt;right:int:opt;top:int:opt;bottom:int:opt;weight_mask:vnode:opt;plane:int:opt;", "clip:vnode;", linearRegressionCreate, (void *)(LINREG_MODE_SINGLE), plugin);
    vspapi->registerFunction("MultiPlane", "clip:vnode;left:int:opt;right:int:opt;top:int:opt;bottom:int:opt;weight_mask:vnode:opt;plane:int:opt;", "clip:vnode;", linearRegressionCreate, (void *)(LINREG_MODE_MULTI), plugin);
    vspapi->registerFunction("SinglePlaneLimited", "clip:vnode;left:int:opt;right:int:opt;top:int:opt;bottom:int:opt;weight_mask:vnode:opt;ref_line_size:int:opt;plane:int:opt;", "clip:vnode;", linearRegressionCreate, (void *)(LINREG_MODE_SINGLE_LIMITED), plugin);
    vspapi->registerFunction("SinglePlaneWeighted", "clip:vnode;left:int:opt;right:int:opt;top:int:opt;bottom:int:opt;weight_mask:vnode:opt;sigmaS:float:opt;sigmaR:float:opt;sigmaD:float:opt;ref_line_size:int:opt;plane:int:opt;", "clip:vnode;", linearRegressionCreate, (void *)(LINREG_MODE_SINGLE_WEIGHTED), plugin);
    vspapi->registerFunction("SinglePlaneDebug", "clip:vnode;left:int:opt;right:int:opt;top:int:opt;bottom:int:opt;weight_mask:vnode:opt;plane:int:opt;", "clip:vnode;", linearRegressionCreate, (void *)(LINREG_MODE_SINGLE_DEBUG), plugin);
}
