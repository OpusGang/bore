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

#include <stdlib.h>
#include <vapoursynth/VapourSynth4.h>
#include <vapoursynth/VSHelper4.h>
#include <gsl/gsl_multifit.h>

typedef struct {
    VSNode *node;
    float lower;
    float upper;
    float thrlo;
    float thrhi;
    int plane;
    int top;
    int bottom;
    int left;
    int right;
    int step;
} FixBrightnessData;

static void processRow(int row, int w, int h, ptrdiff_t stride, float *dstp, FixBrightnessData *d) {
    int x;
    float cur;
    float ref;

    int sign = 1;
    if (row > h / 2)
        sign = -1;

    float tmp;
    float sum = 0.f;
    int div = 0;

    // go through row and get current and reference
    dstp += stride * row;
    for (x = 0; x < w; x += d->step) {
        cur = (dstp[x]);
        ref = (dstp[sign * stride + x]);
        tmp = ref / cur;
        // add to sum if current and reference are within (lower, upper), quotient is within thresholds and there's no sign change
        if (cur < d->upper && cur > d->lower && ref < d->upper && ref > d->lower && (tmp > d->thrlo) && (tmp < d->thrhi) && ((cur > 0 && ref > 0) || (cur < 0 && ref < 0))) {
            sum += tmp;
            div += 1;
        }
    }

    // division by zero => neutral
    if (div == 0)
        sum = 1.f;
    // otherwise just mean of the quotient of the two rows
    else
        sum /= div;

    // adjust each pixel
    for (x = 0; x < w; x++) {
        if (dstp[x] < d->upper && dstp[x] > d->lower)
            dstp[x] *= sum;
    }
}

static void processColumn(int column, int w, int h, ptrdiff_t stride, float *dstp, FixBrightnessData *d) {
    int x;
    float cur;
    float ref;

    int sign = 1;
    if (column > w / 2)
        sign = -1;

    float tmp;
    float sum = 0.f;
    int div = 0;

    // go through column and get current and reference
    for (x = 0; x < h; x += d->step) {
        cur = (dstp[x * stride + column]);
        ref = (dstp[x * stride + column + sign]);
        tmp = ref / cur;
        // add to sum if current and reference are within (lower, upper), quotient is within thresholds and there's no sign change
        if (cur < d->upper && cur > d->lower && ref < d->upper && ref > d->lower && (tmp > d->thrlo) && (tmp < d->thrhi) && ((cur > 0 && ref > 0) || (cur < 0 && ref < 0))) {
            sum += tmp;
            div += 1;
        }
    }

    // division by zero => neutral
    if (div == 0)
        sum = 1.f;
    // otherwise just mean of the quotient of the two columns
    else
        sum /= div;

    // adjust each pixel
    for (x = 0; x < h; x++) {
        if (dstp[x * stride + column] < d->upper && dstp[x * stride + column] > d->lower)
            dstp[x * stride + column] *= sum;
    }
}

static const VSFrame *VS_CC fixBrightnessGetFrame(int n, int activationReason, void *instanceData, void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) {
    FixBrightnessData *d = (FixBrightnessData *)instanceData;

    if (activationReason == arInitial) {
        vsapi->requestFrameFilter(n, d->node, frameCtx);
    } else if (activationReason == arAllFramesReady) {
        const VSFrame *src = vsapi->getFrameFilter(n, d->node, frameCtx);

        VSFrame *dst = vsapi->copyFrame(src, core);

        ptrdiff_t stride = vsapi->getStride(dst, d->plane) / 4;
        int w = vsapi->getFrameWidth(src, d->plane);
        int h = vsapi->getFrameHeight(src, d->plane);
        float *dstp = (float *) vsapi->getWritePtr(dst, d->plane);

        if (d->top != 0) {
            for (int row = d->top - 1; row > -1; --row)
                processRow(row, w, h, stride, dstp, d);
        }
        if (d->bottom != 0) {
            for (int row = h - d->bottom; row < h; ++row)
                processRow(row, w, h, stride, dstp, d);
        }
        if (d->left != 0) {
            for (int column = d->left - 1; column > -1; --column)
                processColumn(column, w, h, stride, dstp, d);
        }
        if (d->right != 0) {
            for (int column = w - d->right; column < w; ++column)
                processColumn(column, w, h, stride, dstp, d);
        }


        vsapi->freeFrame(src);

        return dst;
    }

    return NULL;
}

static void VS_CC fixBrightnessFree(void *instanceData, VSCore *core, const VSAPI *vsapi) {
    FixBrightnessData *d = (FixBrightnessData *)instanceData;
    vsapi->freeNode(d->node);
    free(d);
}

static void VS_CC fixBrightnessCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi) {
    FixBrightnessData d;
    FixBrightnessData *data;
    int err;

    d.node = vsapi->mapGetNode(in, "clip", 0, 0);
    const VSVideoInfo *vi = vsapi->getVideoInfo(d.node);

    if (!vsh_isConstantVideoFormat(vi) || vi->format.sampleType != stFloat) {
        vsapi->mapSetError(out, "FixBrightness: only constant format single float input supported");
        vsapi->freeNode(d.node);
        return;
    }

    d.lower = (float)vsapi->mapGetFloat(in, "lower", 0, &err);
    if (err)
        d.lower = 0.f;

    d.upper = (float)vsapi->mapGetFloat(in, "upper", 0, &err);
    if (err)
        d.upper = 1.f;

    d.thrlo = (float)vsapi->mapGetFloat(in, "thrlo", 0, &err);
    if (err)
        d.thrlo = 0.1;

    d.thrhi = (float)vsapi->mapGetFloat(in, "thrhi", 0, &err);
    if (err)
        d.thrhi = 8.f;

    d.top = vsapi->mapGetInt(in, "top", 0, &err);
    if (err)
        d.top = 0;

    else if (d.top > vi->height / 2) {
        vsapi->mapSetError(out, "FixBrightness: top must be in [0, height / 2]");
        vsapi->freeNode(d.node);
        return;
    }

    d.bottom = vsapi->mapGetInt(in, "bottom", 0, &err);
    if (err)
        d.bottom = 0;

    else if (d.bottom > vi->height / 2) {
        vsapi->mapSetError(out, "FixBrightness: bottom must be in [0, height / 2]");
        vsapi->freeNode(d.node);
        return;
    }

    d.left = vsapi->mapGetInt(in, "left", 0, &err);
    if (err)
        d.left = 0;

    else if (d.left > vi->width / 2) {
        vsapi->mapSetError(out, "FixBrightness: left must be in [0, width / 2]");
        vsapi->freeNode(d.node);
        return;
    }

    d.right = vsapi->mapGetInt(in, "right", 0, &err);
    if (err)
        d.right = 0;

    else if (d.right > vi->width / 2) {
        vsapi->mapSetError(out, "FixBrightness: right must be in [0, width / 2]");
        vsapi->freeNode(d.node);
        return;
    }

    d.step = vsapi->mapGetInt(in, "step", 0, &err);
    if (err)
        d.step = 1;

    d.plane = vsapi->mapGetInt(in, "plane", 0, &err);
    if (err)
        d.plane = 0;

    data = (FixBrightnessData *)malloc(sizeof(d));
    *data = d;

    VSFilterDependency deps[] = {{d.node, rpStrictSpatial}};
    vsapi->createVideoFilter(out, "FixBrightness", vi, fixBrightnessGetFrame, fixBrightnessFree, fmParallel, deps, 1, data, core);
}

typedef struct {
    VSNode *node;
    float lower;
    float upper;
    int plane;
    int top;
    int bottom;
    int left;
    int right;
} MultipleLinearRegressionData;

static void processRowMLR(int row, int w, int h, ptrdiff_t stride, float *dstp, float *dstp1, float *dstp2, float *dstp3, MultipleLinearRegressionData *d) {
    int i;
    int sign = 1;
    if (row > h / 2)
        sign = -1;

    dstp += stride * row;
    dstp1 += stride * row;
    dstp2 += stride * row;
    dstp3 += stride * row;
    gsl_vector *y = gsl_vector_alloc(w);
    gsl_matrix *x = gsl_matrix_alloc(w, 3);
    for (i = 0; i < w; i++) {
        gsl_vector_set(y, i, dstp[i + stride * sign]);
        gsl_matrix_set(x, i, 0, dstp1[i]);
        gsl_matrix_set(x, i, 1, dstp2[i]);
        gsl_matrix_set(x, i, 2, dstp3[i]);
    }

    gsl_multifit_linear_workspace *ws = gsl_multifit_linear_alloc(w, 3);
    double chisq;
    gsl_matrix *cov = gsl_matrix_alloc(3, 3);
    gsl_vector *b = gsl_vector_alloc(3);
    gsl_multifit_linear(x, y, b, cov, &chisq, ws);

    // adjust each pixel
    for (i = 0; i < w; i++) {
        if (dstp[i] < d->upper && dstp[i] > d->lower) {
            dstp[i] = gsl_vector_get(b, 0) * dstp1[i] + gsl_vector_get(b, 1) * dstp2[i] + gsl_vector_get(b, 2) * dstp3[i];
        }
    }
}

static void processColumnMLR(int column, int w, int h, ptrdiff_t stride, float *dstp, float *dstp1, float *dstp2, float *dstp3, MultipleLinearRegressionData *d) {
    int i;
    int j;
    int sign = 1;
    if (column > w / 2)
        sign = -1;


    gsl_vector *y = gsl_vector_alloc(h);
    gsl_matrix *x = gsl_matrix_alloc(h, 3);
    for (i = 0; i < h; i++) {
        gsl_vector_set(y, i, dstp[i * stride + column + sign]);
        gsl_matrix_set(x, i, 0, dstp1[i * stride + column]);
        gsl_matrix_set(x, i, 1, dstp2[i * stride + column]);
        gsl_matrix_set(x, i, 2, dstp3[i * stride + column]);
    }
    gsl_multifit_linear_workspace *ws = gsl_multifit_linear_alloc(h, 3);
    double chisq;
    gsl_matrix *cov = gsl_matrix_alloc(3, 3);
    gsl_vector *b = gsl_vector_alloc(3);
    gsl_multifit_linear(x, y, b, cov, &chisq, ws);

    // adjust each pixel
    for (i = 0; i < h; i++) {
        j = i * stride + column;
        if (dstp[j] < d->upper && dstp[j] > d->lower) {
            dstp[j] = gsl_vector_get(b, 0) * dstp1[j] + gsl_vector_get(b, 1) * dstp2[j] + gsl_vector_get(b, 2) * dstp3[j];
        }
    }
}

static const VSFrame *VS_CC multipleLinearRegressionGetFrame(int n, int activationReason, void *instanceData, void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) {
    MultipleLinearRegressionData *d = (MultipleLinearRegressionData *)instanceData;

    if (activationReason == arInitial) {
        vsapi->requestFrameFilter(n, d->node, frameCtx);
    } else if (activationReason == arAllFramesReady) {
        const VSFrame *src = vsapi->getFrameFilter(n, d->node, frameCtx);

        VSFrame *dst = vsapi->copyFrame(src, core);

        ptrdiff_t stride = vsapi->getStride(dst, d->plane) / 4;
        int w = vsapi->getFrameWidth(src, d->plane);
        int h = vsapi->getFrameHeight(src, d->plane);
        float *dstp = (float *) vsapi->getWritePtr(dst, d->plane);
        float *dstp1 = (float *) vsapi->getWritePtr(dst, 0);
        float *dstp2 = (float *) vsapi->getWritePtr(dst, 1);
        float *dstp3 = (float *) vsapi->getWritePtr(dst, 2);

        if (d->top != 0) {
            for (int row = d->top - 1; row > -1; --row)
                processRowMLR(row, w, h, stride, dstp, dstp1, dstp2, dstp3, d);
        }
        if (d->bottom != 0) {
            for (int row = h - d->bottom; row < h; ++row)
                processRowMLR(row, w, h, stride, dstp, dstp1, dstp2, dstp3, d);
        }
        if (d->left != 0) {
            for (int column = d->left - 1; column > -1; --column)
                processColumnMLR(column, w, h, stride, dstp, dstp1, dstp2, dstp3, d);
        }
        if (d->right != 0) {
            for (int column = w - d->right; column < w; ++column)
                processColumnMLR(column, w, h, stride, dstp, dstp1, dstp2, dstp3, d);
        }


        vsapi->freeFrame(src);

        return dst;
    }

    return NULL;
}

static void VS_CC multipleLinearRegressionFree(void *instanceData, VSCore *core, const VSAPI *vsapi) {
    MultipleLinearRegressionData *d = (MultipleLinearRegressionData *)instanceData;
    vsapi->freeNode(d->node);
    free(d);
}

static void VS_CC multipleLinearRegressionCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi) {
    MultipleLinearRegressionData d;
    MultipleLinearRegressionData *data;
    int err;

    d.node = vsapi->mapGetNode(in, "clip", 0, 0);
    const VSVideoInfo *vi = vsapi->getVideoInfo(d.node);

    if (!vsh_isConstantVideoFormat(vi) || vi->format.sampleType != stFloat || vi->format.subSamplingH > 0 || vi->format.subSamplingW > 0 || vi->format.numPlanes != 3) {
        vsapi->mapSetError(out, "MultipleLinearRegression: only constant format single float input with three planes and no subsampling are supported");
        vsapi->freeNode(d.node);
        return;
    }

    d.lower = (float)vsapi->mapGetFloat(in, "lower", 0, &err);
    if (err)
        d.lower = 0.f;

    d.upper = (float)vsapi->mapGetFloat(in, "upper", 0, &err);
    if (err)
        d.upper = 1.f;

    d.top = vsapi->mapGetInt(in, "top", 0, &err);
    if (err)
        d.top = 0;

    else if (d.top > vi->height / 2) {
        vsapi->mapSetError(out, "MultipleLinearRegression: top must be in [0, height / 2]");
        vsapi->freeNode(d.node);
        return;
    }

    d.bottom = vsapi->mapGetInt(in, "bottom", 0, &err);
    if (err)
        d.bottom = 0;

    else if (d.bottom > vi->height / 2) {
        vsapi->mapSetError(out, "MultipleLinearRegression: bottom must be in [0, height / 2]");
        vsapi->freeNode(d.node);
        return;
    }


    d.left = vsapi->mapGetInt(in, "left", 0, &err);
    if (err)
        d.left = 0;

    else if (d.left > vi->width / 2) {
        vsapi->mapSetError(out, "MultipleLinearRegression: left must be in [0, width / 2]");
        vsapi->freeNode(d.node);
        return;
    }


    d.right = vsapi->mapGetInt(in, "right", 0, &err);
    if (err)
        d.right = 0;

    else if (d.right > vi->width / 2) {
        vsapi->mapSetError(out, "MultipleLinearRegression: right must be in [0, width / 2]");
        vsapi->freeNode(d.node);
        return;
    }

    d.plane = vsapi->mapGetInt(in, "plane", 0, &err);
    if (err)
        d.plane = 0;

    data = (MultipleLinearRegressionData *)malloc(sizeof(d));
    *data = d;

    VSFilterDependency deps[] = {{d.node, rpStrictSpatial}};
    vsapi->createVideoFilter(out, "MultipleLinearRegression", vi, multipleLinearRegressionGetFrame, multipleLinearRegressionFree, fmParallel, deps, 1, data, core);
}

VS_EXTERNAL_API(void) VapourSynthPluginInit2(VSPlugin *plugin, const VSPLUGINAPI *vspapi) {
    vspapi->configPlugin("ng.opusga.bore", "bore", "bore plugin", VS_MAKE_VERSION(1, 0), VAPOURSYNTH_API_VERSION, 0, plugin);
    vspapi->registerFunction("FixBrightness", "clip:vnode;top:int:opt;bottom:int:opt;left:int:opt;right:int:opt;lower:float:opt;upper:float:opt;thrlo:float:opt;thrhi:float:opt;step:int:opt;plane:int:opt;", "clip:vnode;", fixBrightnessCreate, NULL, plugin);
    vspapi->registerFunction("Balance", "clip:vnode;top:int:opt;bottom:int:opt;left:int:opt;right:int:opt;lower:float:opt;upper:float:opt;plane:int:opt;", "clip:vnode;", multipleLinearRegressionCreate, NULL, plugin);
}

