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

#include "common.h"
#include "VapourSynth4.h"
#include "VSHelper4.h"

typedef struct {
    union
    {
        LinearRegressionData data;
    } shared;
    VSNode *node;
    VSNode *ignore_mask;
} VS_LinearRegressionData;

static const VSFrame *VS_CC singlePlaneGetFrame(int n, int activationReason, void *instanceData, void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) {
    VS_LinearRegressionData*d = (VS_LinearRegressionData*)instanceData;

    if (activationReason == arInitial) {
        vsapi->requestFrameFilter(n, d->node, frameCtx);
        if (d->ignore_mask)
            vsapi->requestFrameFilter(n, d->ignore_mask, frameCtx);
    } else if (activationReason == arAllFramesReady) {
        const VSFrame *src = vsapi->getFrameFilter(n, d->node, frameCtx);
        const VSFrame *ignore_mask = NULL;
        const unsigned char *imaskp = NULL;

        VSFrame *dst = vsapi->copyFrame(src, core);

        const int plane = d->shared.data.plane;
        const int top = d->shared.data.top;
        const int bottom = d->shared.data.bottom;
        const int left = d->shared.data.left;
        const int right = d->shared.data.right;
        const int ref_line_size = d->shared.data.ref_line_size;
        const double sigmaS = d->shared.data.sigmaS;
        const double sigmaR = d->shared.data.sigmaR;
        const double sigmaD = d->shared.data.sigmaD;
        ptrdiff_t stride = vsapi->getStride(dst, plane) / 4;
        int w = vsapi->getFrameWidth(src, plane);
        int h = vsapi->getFrameHeight(src, plane);
        float* __restrict dstp = (float *) vsapi->getWritePtr(dst, plane);

        if (d->ignore_mask) {
            ignore_mask = vsapi->getFrameFilter(n, d->ignore_mask, frameCtx);
            ptrdiff_t imaskstride = vsapi->getStride(ignore_mask, 0);
            imaskp = vsapi->getReadPtr(ignore_mask, 0);
            if (top != 0) {
                for (int row = top - 1; row > -1; --row)
                    d->shared.data.processRow(row, w, h, stride, dstp, ref_line_size, sigmaS, sigmaR, sigmaD, imaskp, imaskstride, top - row);
            }
            if (bottom != 0) {
                for (int row = h - bottom; row < h; ++row)
                    d->shared.data.processRow(row, w, h, stride, dstp, ref_line_size, sigmaS, sigmaR, sigmaD, imaskp, imaskstride, bottom + row - h + 1);
            }
            if (left != 0) {
                for (int column = left - 1; column > -1; --column)
                    d->shared.data.processColumn(column, w, h, stride, dstp, ref_line_size, sigmaS, sigmaR, sigmaD, imaskp, imaskstride,left - column);
            }
            if (right != 0) {
                for (int column = w - right; column < w; ++column)
                    d->shared.data.processColumn(column, w, h, stride, dstp, ref_line_size, sigmaS, sigmaR, sigmaD, imaskp, imaskstride, right + column - w + 1);
            }
        } else {
            if (top != 0) {
                for (int row = top - 1; row > -1; --row)
                    d->shared.data.processRow(row, w, h, stride, dstp, ref_line_size, sigmaS, sigmaR, sigmaD);
            }
            if (bottom != 0) {
                for (int row = h - bottom; row < h; ++row)
                    d->shared.data.processRow(row, w, h, stride, dstp, ref_line_size, sigmaS, sigmaR, sigmaD);
            }
            if (left != 0) {
                for (int column = left - 1; column > -1; --column)
                    d->shared.data.processColumn(column, w, h, stride, dstp, ref_line_size, sigmaS, sigmaR, sigmaD);
            }
            if (right != 0) {
                for (int column = w - right; column < w; ++column)
                    d->shared.data.processColumn(column, w, h, stride, dstp, ref_line_size, sigmaS, sigmaR, sigmaD);
            }
        }

        vsapi->freeFrame(ignore_mask);
        vsapi->freeFrame(src);

        return dst;
    }

    return NULL;
}

static const VSFrame *VS_CC multiPlaneGetFrame(int n, int activationReason, void *instanceData, void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) {
    VS_LinearRegressionData *d = (VS_LinearRegressionData *)instanceData;

    if (activationReason == arInitial) {
        vsapi->requestFrameFilter(n, d->node, frameCtx);
    } else if (activationReason == arAllFramesReady) {
        const VSFrame *src = vsapi->getFrameFilter(n, d->node, frameCtx);

        VSFrame *dst = vsapi->copyFrame(src, core);

        const int plane = d->shared.data.plane;
        const int top = d->shared.data.top;
        const int bottom = d->shared.data.bottom;
        const int left = d->shared.data.left;
        const int right = d->shared.data.right;
        ptrdiff_t stride = vsapi->getStride(dst, plane) / 4;
        int w = vsapi->getFrameWidth(src, plane);
        int h = vsapi->getFrameHeight(src, plane);
        float* dstp = (float*)vsapi->getWritePtr(dst, plane);
        float* dstp1 = (float*)vsapi->getWritePtr(dst, 0);
        float* dstp2 = (float*)vsapi->getWritePtr(dst, 1);
        float* dstp3 = (float*)vsapi->getWritePtr(dst, 2);

        if (top != 0) {
            for (int row = top - 1; row > -1; --row)
                processRowMLR(row, w, h, stride, dstp, dstp1, dstp2, dstp3);
        }
        if (bottom != 0) {
            for (int row = h - bottom; row < h; ++row)
                processRowMLR(row, w, h, stride, dstp, dstp1, dstp2, dstp3);
        }
        if (left != 0) {
            for (int column = left - 1; column > -1; --column)
                processColumnMLR(column, w, h, stride, dstp, dstp1, dstp2, dstp3);
        }
        if (right != 0) {
            for (int column = w - right; column < w; ++column)
                processColumnMLR(column, w, h, stride, dstp, dstp1, dstp2, dstp3);
        }

        vsapi->freeFrame(src);

        return dst;
    }

    return NULL;
}

static void set_frame_props(VSFrame* frame, double* props, const VSAPI* vsapi)
{
    VSMap* __restrict dst_props = vsapi->getFramePropertiesRW(frame);
    vsapi->mapSetFloat(dst_props, "BoreAdjustment", props[0], maAppend);
    vsapi->mapSetFloat(dst_props, "BoreCovariance", props[1], maAppend);
    vsapi->mapSetFloat(dst_props, "BoreSumSquares", props[2], maAppend);
}

static const VSFrame *VS_CC singlePlaneDebugGetFrame(int n, int activationReason, void *instanceData, void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) {
    VS_LinearRegressionData *d = (VS_LinearRegressionData *)instanceData;

    if (activationReason == arInitial) {
        vsapi->requestFrameFilter(n, d->node, frameCtx);
        if (d->ignore_mask)
            vsapi->requestFrameFilter(n, d->ignore_mask, frameCtx);
    } else if (activationReason == arAllFramesReady) {
        const VSFrame *src = vsapi->getFrameFilter(n, d->node, frameCtx);
        const VSFrame *ignore_mask = NULL;
        const unsigned char *imaskp = NULL;

        VSFrame *dst = vsapi->copyFrame(src, core);

        const int plane = d->shared.data.plane;
        const int top = d->shared.data.top;
        const int bottom = d->shared.data.bottom;
        const int left = d->shared.data.left;
        const int right = d->shared.data.right;
        ptrdiff_t stride = vsapi->getStride(dst, plane) / 4;
        int w = vsapi->getFrameWidth(src, plane);
        int h = vsapi->getFrameHeight(src, plane);
        float* __restrict dstp = (float *) vsapi->getWritePtr(dst, plane);
        double c1_cov11_sumsq[3] = { 0.0 };
        double* props = c1_cov11_sumsq;

        if (d->ignore_mask) {
            ignore_mask = vsapi->getFrameFilter(n, d->ignore_mask, frameCtx);
            ptrdiff_t imaskstride = vsapi->getStride(ignore_mask, 0);
            imaskp = vsapi->getReadPtr(ignore_mask, 0);
            if (top != 0) {
                for (int row = top - 1; row > -1; --row)
                {
                    debugRowSLRMasked(row, w, h, stride, dstp, &props, imaskp, imaskstride, top - row);
                    set_frame_props(dst, props, vsapi);
                }
            }
            if (bottom != 0) {
                for (int row = h - bottom; row < h; ++row)
                {
                    debugRowSLRMasked(row, w, h, stride, dstp, &props, imaskp, imaskstride, bottom + row - h + 1);
                    set_frame_props(dst, props, vsapi);
                }
            }
            if (left != 0) {
                for (int column = left - 1; column > -1; --column)
                {
                    debugColumnSLRMasked(column, w, h, stride, dstp, &props, imaskp, imaskstride, left - column);
                    set_frame_props(dst, props, vsapi);
                }
            }
            if (right != 0) {
                for (int column = w - right; column < w; ++column)
                {
                    debugColumnSLRMasked(column, w, h, stride, dstp, &props, imaskp, imaskstride, right + column - w + 1);
                    set_frame_props(dst, props, vsapi);
                }
            }
            vsapi->freeFrame(ignore_mask);
        } else {
            if (top != 0) {
                for (int row = top - 1; row > -1; --row)
                {
                    debugRowSLR(row, w, h, stride, dstp, &props);
                    set_frame_props(dst, props, vsapi);
                }
            }
            if (bottom != 0) {
                for (int row = h - bottom; row < h; ++row)
                {
                    debugRowSLR(row, w, h, stride, dstp, &props);
                    set_frame_props(dst, props, vsapi);
                }
            }
            if (left != 0) {
                for (int column = left - 1; column > -1; --column)
                {
                    debugColumnSLR(column, w, h, stride, dstp, &props);
                    set_frame_props(dst, props, vsapi);
                }
            }
            if (right != 0) {
                for (int column = w - right; column < w; ++column)
                {
                    debugColumnSLR(column, w, h, stride, dstp, &props);
                    set_frame_props(dst, props, vsapi);
                }
            }
        }
        
        vsapi->freeFrame(src);

        return dst;
    }

    return NULL;
}

static void VS_CC linearRegressionFree(void *instanceData, VSCore *core, const VSAPI *vsapi) {
    VS_LinearRegressionData *d = (VS_LinearRegressionData *)instanceData;
    vsapi->freeNode(d->node);
    vsapi->freeNode(d->ignore_mask);
    free(d);
}

static void VS_CC linearRegressionCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi) {
    LinRegMode mode = (LinRegMode) userData;
    int err, w, h;

    VSNode* node = vsapi->mapGetNode(in, "clip", 0, 0);
    const VSVideoInfo *vi = vsapi->getVideoInfo(node);

    int plane = (int)vsapi->mapGetInt(in, "plane", 0, &err);
    if (err)
        plane = 0;

    if (!vsh_isConstantVideoFormat(vi) || vi->format.sampleType != stFloat) {
        vsapi->mapSetError(out, "bore: only constant format single float clip input is supported");
        vsapi->freeNode(node);
        return;
    }
    
    if (vi->width == 0 || vi->height == 0) {
        vsapi->mapSetError(out, "bore: only constant resolution clip input is supported");
        vsapi->freeNode(node);
        return;
    }

    if (vi->format.numPlanes == 1 && plane > 0) {
        vsapi->mapSetError(out, "bore: plane cannot be bigger than 0");
        vsapi->freeNode(node);
        return;
    }

    w = vi->width;
    h = vi->height;
    
    if (plane > 0) {
        w = w >> vi->format.subSamplingW;
        h = h >> vi->format.subSamplingH;
    }
    
    VSNode* ignore_mask = vsapi->mapGetNode(in, "ignore_mask", 0, &err);
    if (err)
        ignore_mask = NULL;
    else {
        const VSVideoInfo *ivi = vsapi->getVideoInfo(ignore_mask);
        if (!vsh_isConstantVideoFormat(ivi) || (ivi->format.sampleType != stInteger && ivi->format.bitsPerSample != 8)) {
            vsapi->mapSetError(out, "bore: only constant format 8-bit ignore_mask input is supported");
            vsapi->freeNode(node);
            vsapi->freeNode(ignore_mask);
            return;
        }

        if (ivi->format.numPlanes > 1) {
            vsapi->mapSetError(out, "bore: ignore_mask must be in gray format");
            vsapi->freeNode(node);
            vsapi->freeNode(ignore_mask);
        return;
        }

        if (ivi->width != w || ivi->height != h) {
            vsapi->mapSetError(out, "bore: clip and ignore_mask must have matching dimensions");
            vsapi->freeNode(node);
            vsapi->freeNode(ignore_mask);
            return;
        }

        if (ivi->width == 0 || ivi->height == 0) {
            vsapi->mapSetError(out, "bore: only constant resolution ignore_mask input is supported");
            vsapi->freeNode(node);
            vsapi->freeNode(ignore_mask);
            return;
        }
    }

    int top = (int)vsapi->mapGetInt(in, "top", 0, &err);
    if (err)
        top = 0;
    else if (top > h / 2) {
        vsapi->mapSetError(out, "bore: top must be in [0, height / 2]");
        vsapi->freeNode(node);
        if (ignore_mask)
            vsapi->freeNode(ignore_mask);
        return;
    }

    int bottom = (int)vsapi->mapGetInt(in, "bottom", 0, &err);
    if (err)
        bottom = 0;
    else if (bottom > h / 2) {
        vsapi->mapSetError(out, "bore: bottom must be in [0, height / 2]");
        vsapi->freeNode(node);
        if (ignore_mask)
            vsapi->freeNode(ignore_mask);
        return;
    }

    int left = (int)vsapi->mapGetInt(in, "left", 0, &err);
    if (err)
        left = 0;
    else if (left > w / 2) {
        vsapi->mapSetError(out, "bore: left must be in [0, width / 2]");
        vsapi->freeNode(node);
        if (ignore_mask)
            vsapi->freeNode(ignore_mask);
        return;
    }

    int right = (int)vsapi->mapGetInt(in, "right", 0, &err);
    if (err)
        right = 0;
    else if (right > w / 2) {
        vsapi->mapSetError(out, "bore: right must be in [0, width / 2]");
        vsapi->freeNode(node);
        if (ignore_mask)
            vsapi->freeNode(ignore_mask);
        return;
    }

    int ref_line_size = 0;
    if (mode == LINREG_MODE_SINGLE_LIMITED || mode == LINREG_MODE_SINGLE_WEIGHTED) {
        ref_line_size = (int)vsapi->mapGetInt(in, "ref_line_size", 0, &err);
        if (err)
            ref_line_size = 100;
    }

    VS_LinearRegressionData* d = (VS_LinearRegressionData*)malloc(sizeof(VS_LinearRegressionData));
    d->node = node;
    d->ignore_mask = ignore_mask;
    d->shared.data.plane = plane;
    d->shared.data.top = top;
    d->shared.data.bottom = bottom;
    d->shared.data.left = left;
    d->shared.data.right = right;
    d->shared.data.ref_line_size = ref_line_size;

    switch (mode) {
        case LINREG_MODE_SINGLE:
            if (ignore_mask) {
                d->shared.data.processRow = &processRowSLRMasked;
                d->shared.data.processColumn = &processColumnSLRMasked;
            } else {
                d->shared.data.processRow = &processRowSLR;
                d->shared.data.processColumn = &processColumnSLR;
            }
            break;
        case LINREG_MODE_MULTI:
            if (vi->format.numPlanes == 1) {
                vsapi->mapSetError(out, "bore: clip must have 3 planes");
                vsapi->freeNode(node);
                if (ignore_mask)
                    vsapi->freeNode(ignore_mask);
                free(d);
                return;
            }
            break;
        case LINREG_MODE_SINGLE_LIMITED:
            if (ignore_mask) {
                d->shared.data.processRow = &processRowSLRRefMasked;
                d->shared.data.processColumn = &processColumnSLRRefMasked;
            } else {
                d->shared.data.processRow = &processRowSLRRef;
                d->shared.data.processColumn = &processColumnSLRRef;
            }
            break;
        case LINREG_MODE_SINGLE_WEIGHTED:
            d->shared.data.sigmaS = vsapi->mapGetFloat(in, "sigmaS", 0, &err);
            if (err)
                d->shared.data.sigmaS = 50.0;

            d->shared.data.sigmaR = vsapi->mapGetFloat(in, "sigmaR", 0, &err);
            if (err)
                d->shared.data.sigmaR = 0.5;

            d->shared.data.sigmaD = vsapi->mapGetFloat(in, "sigmaD", 0, &err);
            if (err)
                d->shared.data.sigmaD = 1.5;

            if (ignore_mask) {
                d->shared.data.processRow = &processRowWSLRMasked;
                d->shared.data.processColumn = &processColumnWSLRMasked;
            } else {
                d->shared.data.processRow = &processRowWSLR;
                d->shared.data.processColumn = &processColumnWSLR;
            }
            break;
        case LINREG_MODE_SINGLE_DEBUG:
            break;
    }

    VSFilterDependency deps[] = {{node, rpStrictSpatial}};

    switch (mode) {
        case LINREG_MODE_SINGLE:
            vsapi->createVideoFilter(out, "SinglePlane", vi, singlePlaneGetFrame, linearRegressionFree, fmParallel, deps, 1, d, core);
            break;
        case LINREG_MODE_MULTI:
            vsapi->createVideoFilter(out, "MultiPlane", vi, multiPlaneGetFrame, linearRegressionFree, fmParallel, deps, 1, d, core);
            break;
        case LINREG_MODE_SINGLE_LIMITED:
            vsapi->createVideoFilter(out, "SinglePlaneLimited", vi, singlePlaneGetFrame, linearRegressionFree, fmParallel, deps, 1, d, core);
            break;
        case LINREG_MODE_SINGLE_WEIGHTED:
            vsapi->createVideoFilter(out, "SinglePlaneWeighted", vi, singlePlaneGetFrame, linearRegressionFree, fmParallel, deps, 1, d, core);
            break;
        case LINREG_MODE_SINGLE_DEBUG:
            vsapi->createVideoFilter(out, "SinglePlaneDebug", vi, singlePlaneDebugGetFrame, linearRegressionFree, fmParallel, deps, 1, d, core);
            break;
    }
}

VS_EXTERNAL_API(void) VapourSynthPluginInit2(VSPlugin *plugin, const VSPLUGINAPI *vspapi) {
    vspapi->configPlugin("ng.opusga.bore", "bore", "bore plugin", VS_MAKE_VERSION(1, 0), VAPOURSYNTH_API_VERSION, 0, plugin);
    vspapi->registerFunction("SinglePlane", "clip:vnode;left:int:opt;right:int:opt;top:int:opt;bottom:int:opt;ignore_mask:vnode:opt;plane:int:opt;", "clip:vnode;", linearRegressionCreate, (void *)(LINREG_MODE_SINGLE), plugin);
    vspapi->registerFunction("MultiPlane", "clip:vnode;left:int:opt;right:int:opt;top:int:opt;bottom:int:opt;ignore_mask:vnode:opt;plane:int:opt;", "clip:vnode;", linearRegressionCreate, (void *)(LINREG_MODE_MULTI), plugin);
    vspapi->registerFunction("SinglePlaneLimited", "clip:vnode;left:int:opt;right:int:opt;top:int:opt;bottom:int:opt;ignore_mask:vnode:opt;ref_line_size:int:opt;plane:int:opt;", "clip:vnode;", linearRegressionCreate, (void *)(LINREG_MODE_SINGLE_LIMITED), plugin);
    vspapi->registerFunction("SinglePlaneWeighted", "clip:vnode;left:int:opt;right:int:opt;top:int:opt;bottom:int:opt;ignore_mask:vnode:opt;sigmaS:float:opt;sigmaR:float:opt;sigmaD:float:opt;ref_line_size:int:opt;plane:int:opt;", "clip:vnode;", linearRegressionCreate, (void *)(LINREG_MODE_SINGLE_WEIGHTED), plugin);
    vspapi->registerFunction("SinglePlaneDebug", "clip:vnode;left:int:opt;right:int:opt;top:int:opt;bottom:int:opt;ignore_mask:vnode:opt;plane:int:opt;", "clip:vnode;", linearRegressionCreate, (void *)(LINREG_MODE_SINGLE_DEBUG), plugin);
}
