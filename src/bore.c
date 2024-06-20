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

#include <math.h>
#include <stdlib.h>
#include <vapoursynth/VapourSynth4.h>
#include <vapoursynth/VSHelper4.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>

typedef struct {
    VSNode *node;
    VSNode *ignore_mask;
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
        // add to sum if quotient is within thresholds and there's no sign change
        if ((tmp > d->thrlo) && (tmp < d->thrhi) && ((cur > 0 && ref > 0) || (cur < 0 && ref < 0))) {
            sum += tmp;
            div += 1;
        }
    }

    // division by zero => neutral
    if (div == 0)
        return;
    // otherwise just mean of the quotient of the two rows
    else
        sum /= div;

    // adjust each pixel
    for (x = 0; x < w; x++) {
        dstp[x] *= sum;
    }
}

static void processRowMasked(int row, int w, int h, ptrdiff_t stride, float *dstp, const unsigned char * restrict imaskp, FixBrightnessData *d) {
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
    imaskp += stride * row;
    for (x = 0; x < w; x += d->step) {
        if (imaskp[x] < 128 && imaskp[sign * stride + x] < 128) {
            cur = (dstp[x]);
            ref = (dstp[sign * stride + x]);
            tmp = ref / cur;
            // add to sum if quotient is within thresholds and there's no sign change
            if ((tmp > d->thrlo) && (tmp < d->thrhi) && ((cur > 0 && ref > 0) || (cur < 0 && ref < 0))) {
                sum += tmp;
                div += 1;
            }
        }
    }

    // division by zero => neutral
    if (div == 0)
        return;
    // otherwise just mean of the quotient of the two rows
    else
        sum /= div;

    // adjust each pixel
    for (x = 0; x < w; x++) {
        if (imaskp[x] < 128 && imaskp[sign * stride + x] < 128) {
            dstp[x] *= sum;
        }
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
        // add to sum if quotient is within thresholds and there's no sign change
        if ((tmp > d->thrlo) && (tmp < d->thrhi) && ((cur > 0 && ref > 0) || (cur < 0 && ref < 0))) {
            sum += tmp;
            div += 1;
        }
    }

    // division by zero => neutral
    if (div == 0)
        return;
    // otherwise just mean of the quotient of the two columns
    else
        sum /= div;

    // adjust each pixel
    for (x = 0; x < h; x++) {
        dstp[x * stride + column] *= sum;
    }
}

static void processColumnMasked(int column, int w, int h, ptrdiff_t stride, float *dstp, const unsigned char * restrict imaskp, FixBrightnessData *d) {
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
        if (imaskp[x * stride + column] < 128 && imaskp[x * stride + column + sign] < 128) {
            cur = (dstp[x * stride + column]);
            ref = (dstp[x * stride + column + sign]);
            tmp = ref / cur;
            // add to sum if quotient is within thresholds and there's no sign change
            if ((tmp > d->thrlo) && (tmp < d->thrhi) && ((cur > 0 && ref > 0) || (cur < 0 && ref < 0))) {
                sum += tmp;
                div += 1;
            }
        }
    }

    // division by zero => neutral
    if (div == 0)
        return;
    // otherwise just mean of the quotient of the two columns
    else
        sum /= div;

    // adjust each pixel
    for (x = 0; x < h; x++) {
        if (imaskp[x * stride + column] < 128 && imaskp[x * stride + column + sign] < 128) {
            dstp[x * stride + column] *= sum;
        }
    }
}

static const VSFrame *VS_CC fixBrightnessGetFrame(int n, int activationReason, void *instanceData, void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) {
    FixBrightnessData *d = (FixBrightnessData *)instanceData;

    if (activationReason == arInitial) {
        vsapi->requestFrameFilter(n, d->node, frameCtx);
        if (d->ignore_mask)
            vsapi->requestFrameFilter(n, d->ignore_mask, frameCtx);
    } else if (activationReason == arAllFramesReady) {
        const VSFrame *src = vsapi->getFrameFilter(n, d->node, frameCtx);
        const VSFrame *ignore_mask = NULL;
        const unsigned char *imaskp = NULL;

        VSFrame *dst = vsapi->copyFrame(src, core);

        ptrdiff_t stride = vsapi->getStride(dst, d->plane) / 4;
        int w = vsapi->getFrameWidth(src, d->plane);
        int h = vsapi->getFrameHeight(src, d->plane);
        float *dstp = (float *) vsapi->getWritePtr(dst, d->plane);

        if (d->ignore_mask) {
            ignore_mask = vsapi->getFrameFilter(n, d->ignore_mask, frameCtx);
            imaskp = vsapi->getReadPtr(ignore_mask, d->plane);

            if (d->top != 0) {
                for (int row = d->top - 1; row > -1; --row)
                    processRowMasked(row, w, h, stride, dstp, imaskp, d);
            }
            if (d->bottom != 0) {
                for (int row = h - d->bottom; row < h; ++row)
                    processRowMasked(row, w, h, stride, dstp, imaskp, d);
            }
            if (d->left != 0) {
                for (int column = d->left - 1; column > -1; --column)
                    processColumnMasked(column, w, h, stride, dstp, imaskp, d);
            }
            if (d->right != 0) {
                for (int column = w - d->right; column < w; ++column)
                    processColumnMasked(column, w, h, stride, dstp, imaskp, d);
            }
        } else {
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
        }


        vsapi->freeFrame(src);
        vsapi->freeFrame(ignore_mask);

        return dst;
    }

    return NULL;
}

static void VS_CC fixBrightnessFree(void *instanceData, VSCore *core, const VSAPI *vsapi) {
    FixBrightnessData *d = (FixBrightnessData *)instanceData;
    vsapi->freeNode(d->node);
    vsapi->freeNode(d->ignore_mask);
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

    d.ignore_mask = vsapi->mapGetNode(in, "ignore_mask", 0, &err);
    if (err)
        d.ignore_mask = NULL;
    else {
        vsapi->freeNode(d.ignore_mask);
    }
    
    d.thrlo = (float)vsapi->mapGetFloat(in, "thrlo", 0, &err);
    if (err)
        d.thrlo = 0.1;

    d.thrhi = (float)vsapi->mapGetFloat(in, "thrhi", 0, &err);
    if (err)
        d.thrhi = 8.f;

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
        vsapi->freeNode(d.ignore_mask);
        return;
    }

    d.bottom = vsapi->mapGetInt(in, "bottom", 0, &err);
    if (err)
        d.bottom = 0;

    else if (d.bottom > vi->height / 2) {
        vsapi->mapSetError(out, "FixBrightness: bottom must be in [0, height / 2]");
        vsapi->freeNode(d.node);
        vsapi->freeNode(d.ignore_mask);
        return;
    }

    d.left = vsapi->mapGetInt(in, "left", 0, &err);
    if (err)
        d.left = 0;

    else if (d.left > vi->width / 2) {
        vsapi->mapSetError(out, "FixBrightness: left must be in [0, width / 2]");
        vsapi->freeNode(d.node);
        vsapi->freeNode(d.ignore_mask);
        return;
    }

    d.right = vsapi->mapGetInt(in, "right", 0, &err);
    if (err)
        d.right = 0;

    else if (d.right > vi->width / 2) {
        vsapi->mapSetError(out, "FixBrightness: right must be in [0, width / 2]");
        vsapi->freeNode(d.node);
        vsapi->freeNode(d.ignore_mask);
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
    int plane;
    int top;
    int bottom;
    int left;
    int right;
    int mode;
    int ref_line_size;
    double sigmaS;
    double sigmaR;
} LinearRegressionData;

static void processRowSLR(int row, int w, int h, ptrdiff_t stride, float *dstp) {
    int sign = 1;
    if (row > h / 2)
        sign = -1;

    double *cur, *ref;
    int i;

    cur = malloc(sizeof(double) * w);
    ref = malloc(sizeof(double) * w);

    dstp += row * stride;
    for (i = 0; i < w; i++) {
        cur[i] = dstp[i];
        ref[i] = dstp[sign * stride + i];
    }

    double cov11, sumsq;

    double c1;

    const double *const_cur = cur;
    const double *const_ref = ref;

    int status = gsl_fit_mul(const_cur, 1, const_ref, 1, w, &c1, &cov11, &sumsq);

    if (!status && isfinite(c1)) {
        // adjust each pixel
        for (i = 0; i < w; i++) {
            dstp[i] *= c1;
        }
    }

    free(cur);
    free(ref);
}

static void processColumnSLR(int column, int w, int h, ptrdiff_t stride, float *dstp) {
    int sign = 1;
    if (column > w / 2)
        sign = -1;

    double *cur, *ref;
    int i;

    cur = malloc(sizeof(double) * h);
    ref = malloc(sizeof(double) * h);

    dstp += column;
    for (i = 0; i < h; i++) {
        cur[i] = dstp[i * stride];
        ref[i] = dstp[sign + stride * i];
    }

    double cov11, sumsq;

    double c1;

    const double *const_cur = cur;
    const double *const_ref = ref;

    int status = gsl_fit_mul(const_cur, 1, const_ref, 1, h, &c1, &cov11, &sumsq);

    if (!status && isfinite(c1)) {
        // adjust each pixel
        for (i = 0; i < h; i++) {
            dstp[i * stride] *= c1;
        }
    }

    free(cur);
    free(ref);
}

static void debugRowSLR(int row, int w, int h, ptrdiff_t stride, float *dstp, VSFrame *dst, const VSAPI *vsapi) {
    int sign = 1;
    if (row > h / 2)
        sign = -1;

    double *cur, *ref;
    int i;

    cur = malloc(sizeof(double) * w);
    ref = malloc(sizeof(double) * w);

    dstp += row * stride;
    for (i = 0; i < w; i++) {
        cur[i] = dstp[i];
        ref[i] = dstp[sign * stride + i];
    }

    double cov11, sumsq;

    double c1;

    const double *const_cur = cur;
    const double *const_ref = ref;

    int status = gsl_fit_mul(const_cur, 1, const_ref, 1, w, &c1, &cov11, &sumsq);

    /* if (!status || isfinite(c1)) { */
    /*     c1 = 0; */
    /* } */
    vsapi->mapSetFloat(vsapi->getFramePropertiesRW(dst), "BalanceAdjustment", c1, maAppend);
    /* vsapi->mapSetFloat(vsapi->getFramePropertiesRW(dst), "BalanceCovariance", cov11, maAppend); */
    /* vsapi->mapSetFloat(vsapi->getFramePropertiesRW(dst), "BalanceSumSquares", sumsq, maAppend); */

    free(cur);
    free(ref);
}

static void debugColumnSLR(int column, int w, int h, ptrdiff_t stride, float *dstp, VSFrame *dst, const VSAPI *vsapi) {
    int sign = 1;
    if (column > w / 2)
        sign = -1;

    double *cur, *ref;
    int i;

    cur = malloc(sizeof(double) * h);
    ref = malloc(sizeof(double) * h);

    dstp += column;
    for (i = 0; i < h; i++) {
        cur[i] = dstp[i * stride];
        ref[i] = dstp[sign + stride * i];
    }

    double cov11, sumsq;

    double c1;

    const double *const_cur = cur;
    const double *const_ref = ref;

    int status = gsl_fit_mul(const_cur, 1, const_ref, 1, h, &c1, &cov11, &sumsq);

    if (!status || isfinite(c1)) {
        c1 = 0;
    }
    vsapi->mapSetFloat(vsapi->getFramePropertiesRW(dst), "BalanceAdjustment", c1, maAppend);
    vsapi->mapSetFloat(vsapi->getFramePropertiesRW(dst), "BalanceCovariance", cov11, maAppend);
    vsapi->mapSetFloat(vsapi->getFramePropertiesRW(dst), "BalanceSumSquares", sumsq, maAppend);

    free(cur);
    free(ref);
}

static void processRowMLR(int row, int w, int h, ptrdiff_t stride, float *dstp, float *dstp1, float *dstp2, float *dstp3) {
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
    int status = gsl_multifit_linear(x, y, b, cov, &chisq, ws);

    if (!status) {
        // adjust each pixel
        for (i = 0; i < w; i++) {
            dstp[i] = gsl_vector_get(b, 0) * dstp1[i] + gsl_vector_get(b, 1) * dstp2[i] + gsl_vector_get(b, 2) * dstp3[i];
        }
    }

    gsl_vector_free(b);
    gsl_matrix_free(cov);
    gsl_multifit_linear_free(ws);
    gsl_matrix_free(x);
    gsl_vector_free(y);
}

static void processColumnMLR(int column, int w, int h, ptrdiff_t stride, float *dstp, float *dstp1, float *dstp2, float *dstp3) {
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
    int status = gsl_multifit_linear(x, y, b, cov, &chisq, ws);

    if (!status) {
        // adjust each pixel
        for (i = 0; i < h; i++) {
            j = i * stride + column;
            dstp[j] = gsl_vector_get(b, 0) * dstp1[j] + gsl_vector_get(b, 1) * dstp2[j] + gsl_vector_get(b, 2) * dstp3[j];
        }
    }

    gsl_vector_free(b);
    gsl_matrix_free(cov);
    gsl_multifit_linear_free(ws);
    gsl_matrix_free(x);
    gsl_vector_free(y);
}

static void processRowSLRRef(int row, int w, int h, ptrdiff_t stride, float *dstp, int ref_line_size) {
    int sign = 1;
    if (row > h / 2)
        sign = -1;

    int i;
    int start, stop;
    int status;

    double c1;
    double cov11, sumsq;

    double *cur, *ref;
    cur = malloc(sizeof(double) * w);
    ref = malloc(sizeof(double) * w);

    dstp += row * stride;
    for (i = 0; i < w; i++) {
        cur[i] = dstp[i];
        ref[i] = dstp[sign * stride + i];
    }

    for (i = 0; i < w; i++) {
        start = i - ref_line_size;
        stop = i + ref_line_size;
        if (start < 0)
            start = 0;
        if (stop >= w)
            stop = w - 1;

        const double *const_cur = cur + start;
        const double *const_ref = ref + start;

        status = gsl_fit_mul(const_cur, 1, const_ref, 1, stop - start, &c1, &cov11, &sumsq);

        if (!status && isfinite(c1)) 
            dstp[i] *= c1;
    }

    free(cur);
    free(ref);
}

static void processColumnSLRRef(int column, int w, int h, ptrdiff_t stride, float *dstp, int ref_line_size) {
    int sign = 1;
    if (column > w / 2)
        sign = -1;

    int i;
    int start, stop;
    int status;

    double c1;
    double cov11, sumsq;

    double *cur, *ref;
    cur = malloc(sizeof(double) * ref_line_size);
    ref = malloc(sizeof(double) * ref_line_size);

    dstp += column;
    for (i = 0; i < h; i++) {
        cur[i] = dstp[i * stride];
        ref[i] = dstp[sign + stride * i];
    }

    for (i = 0; i < h; i++) {
        start = i - ref_line_size;
        stop = i + ref_line_size;
        if (start < 0)
            start = 0;
        if (stop >= h)
            stop = h - 1;

        const double *const_cur = cur + start;
        const double *const_ref = ref + start;

        status = gsl_fit_mul(const_cur, 1, const_ref, 1, stop - start, &c1, &cov11, &sumsq);

        if (!status && isfinite(c1)) 
            dstp[i * stride] *= c1;
    }

    free(cur);
    free(ref);
}

static void processRowWSLR(int row, int w, int h, ptrdiff_t stride, float *dstp, int ref_line_size, double sigmaS, double sigmaR) {
    int sign = 1;
    if (row > h / 2)
        sign = -1;

    int i;
    int j;
    int k;
    int start, stop;
    int status;

    double c1;
    double cov11, sumsq;

    double *cur, *ref;
    cur = malloc(sizeof(double) * w);
    ref = malloc(sizeof(double) * w);

    dstp += row * stride;
    for (i = 0; i < w; i++) {
        cur[i] = dstp[i];
        ref[i] = dstp[sign * stride + i];
    }

    for (i = 0; i < w; i++) {
        start = i - ref_line_size;
        stop = i + ref_line_size;
        if (start < 0)
            start = 0;
        if (stop >= w)
            stop = w - 1;

        const double *const_cur = cur + start;
        const double *const_ref = ref + start;

        double *weights;
        weights = malloc(sizeof(double) * (stop - start));

        for (j = start; j < stop; j++) {
            k = j - start;
            weights[k] = exp(-((j - i) * (j - i)) / (sigmaS * sigmaS) - ((dstp[j] - dstp[i]) * (dstp[j] - dstp[i])) / (sigmaR * sigmaR));
        }

        status = gsl_fit_wmul(const_cur, 1, weights, 1, const_ref, 1, stop - start, &c1, &cov11, &sumsq);

        if (!status && isfinite(c1)) 
            dstp[i] *= c1;

        free(weights);
    }

    free(cur);
    free(ref);
}

static void processColumnWSLR(int column, int w, int h, ptrdiff_t stride, float *dstp, int ref_line_size, double sigmaS, double sigmaR) {
    int sign = 1;
    if (column > w / 2)
        sign = -1;

    int i;
    int j;
    int k;
    int start, stop;
    int status;

    double c1;
    double cov11, sumsq;

    double *cur, *ref;
    cur = malloc(sizeof(double) * w);
    ref = malloc(sizeof(double) * w);

    dstp += column;
    for (i = 0; i < h; i++) {
        cur[i] = dstp[i];
        ref[i] = dstp[sign + i];
    }
    for (i = 0; i < h; i++) {
        start = i - ref_line_size;
        stop = i + ref_line_size;
        if (start < 0)
            start = 0;
        if (stop >= h)
            stop = h - 1;

        const double *const_cur = cur + start;
        const double *const_ref = ref + start;
        
        double *weights;
        weights = malloc(sizeof(double) * (stop - start));
        for (j = start; j < stop; j++) {
            k = j - start;
            weights[k] = exp(-((j - i) * (j - i)) / (sigmaS * sigmaS) - ((dstp[j] - dstp[i]) * (dstp[j] - dstp[i])) / (sigmaR * sigmaR));
        }

        status = gsl_fit_wmul(const_cur, 1, weights, 1, const_ref, 1, stop - start, &c1, &cov11, &sumsq);

        if (!status && isfinite(c1)) 
            dstp[i * stride] *= c1;

        free(weights);
    }

    free(cur);
    free(ref);
}

static const VSFrame *VS_CC singlePlaneGetFrame(int n, int activationReason, void *instanceData, void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) {
    LinearRegressionData *d = (LinearRegressionData *)instanceData;

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
                processRowSLR(row, w, h, stride, dstp);
        }
        if (d->bottom != 0) {
            for (int row = h - d->bottom; row < h; ++row)
                processRowSLR(row, w, h, stride, dstp);
        }
        if (d->left != 0) {
            for (int column = d->left - 1; column > -1; --column)
                processColumnSLR(column, w, h, stride, dstp);
        }
        if (d->right != 0) {
            for (int column = w - d->right; column < w; ++column)
                processColumnSLR(column, w, h, stride, dstp);
        }

        vsapi->freeFrame(src);

        return dst;
    }

    return NULL;
}

static const VSFrame *VS_CC multiPlaneGetFrame(int n, int activationReason, void *instanceData, void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) {
    LinearRegressionData *d = (LinearRegressionData *)instanceData;

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
                processRowMLR(row, w, h, stride, dstp, dstp1, dstp2, dstp3);
        }
        if (d->bottom != 0) {
            for (int row = h - d->bottom; row < h; ++row)
                processRowMLR(row, w, h, stride, dstp, dstp1, dstp2, dstp3);
        }
        if (d->left != 0) {
            for (int column = d->left - 1; column > -1; --column)
                processColumnMLR(column, w, h, stride, dstp, dstp1, dstp2, dstp3);
        }
        if (d->right != 0) {
            for (int column = w - d->right; column < w; ++column)
                processColumnMLR(column, w, h, stride, dstp, dstp1, dstp2, dstp3);
        }

        vsapi->freeFrame(src);

        return dst;
    }

    return NULL;
}

static const VSFrame *VS_CC singlePlaneLimitedGetFrame(int n, int activationReason, void *instanceData, void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) {
    LinearRegressionData *d = (LinearRegressionData *)instanceData;

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
                processRowSLRRef(row, w, h, stride, dstp, d->ref_line_size);
        }
        if (d->bottom != 0) {
            for (int row = h - d->bottom; row < h; ++row)
                processRowSLRRef(row, w, h, stride, dstp, d->ref_line_size);
        }
        if (d->left != 0) {
            for (int column = d->left - 1; column > -1; --column)
                processColumnSLRRef(column, w, h, stride, dstp, d->ref_line_size);
        }
        if (d->right != 0) {
            for (int column = w - d->right; column < w; ++column)
                processColumnSLRRef(column, w, h, stride, dstp, d->ref_line_size);
        }

        vsapi->freeFrame(src);

        return dst;
    }

    return NULL;
}

static const VSFrame *VS_CC singlePlaneWeightedGetFrame(int n, int activationReason, void *instanceData, void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) {
    LinearRegressionData *d = (LinearRegressionData *)instanceData;

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
                processRowWSLR(row, w, h, stride, dstp, d->ref_line_size, d->sigmaS, d->sigmaR);
        }
        if (d->bottom != 0) {
            for (int row = h - d->bottom; row < h; ++row)
                processRowWSLR(row, w, h, stride, dstp, d->ref_line_size, d->sigmaS, d->sigmaR);
        }
        if (d->left != 0) {
            for (int column = d->left - 1; column > -1; --column)
                processColumnWSLR(column, w, h, stride, dstp, d->ref_line_size, d->sigmaS, d->sigmaR);
        }
        if (d->right != 0) {
            for (int column = w - d->right; column < w; ++column)
                processColumnWSLR(column, w, h, stride, dstp, d->ref_line_size, d->sigmaS, d->sigmaR);
        }

        vsapi->freeFrame(src);

        return dst;
    }

    return NULL;
}

static const VSFrame *VS_CC singlePlaneDebugGetFrame(int n, int activationReason, void *instanceData, void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) {
    LinearRegressionData *d = (LinearRegressionData *)instanceData;

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
                debugRowSLR(row, w, h, stride, dstp, dst, vsapi);
        }
        if (d->bottom != 0) {
            for (int row = h - d->bottom; row < h; ++row)
                debugRowSLR(row, w, h, stride, dstp, dst, vsapi);
        }
        if (d->left != 0) {
            for (int column = d->left - 1; column > -1; --column)
                debugColumnSLR(column, w, h, stride, dstp, dst, vsapi);
        }
        if (d->right != 0) {
            for (int column = w - d->right; column < w; ++column)
                debugColumnSLR(column, w, h, stride, dstp, dst, vsapi);
        }
        
        vsapi->freeFrame(src);

        return dst;
    }

    return NULL;
}

static void VS_CC linearRegressionFree(void *instanceData, VSCore *core, const VSAPI *vsapi) {
    LinearRegressionData *d = (LinearRegressionData *)instanceData;
    vsapi->freeNode(d->node);
    free(d);
}

static void VS_CC singlePlaneCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi) {
    LinearRegressionData d;
    LinearRegressionData *data;
    int err;

    d.node = vsapi->mapGetNode(in, "clip", 0, 0);
    const VSVideoInfo *vi = vsapi->getVideoInfo(d.node);

    d.top = vsapi->mapGetInt(in, "top", 0, &err);
    if (err)
        d.top = 0;

    else if (d.top > vi->height / 2) {
        vsapi->mapSetError(out, "Balance: top must be in [0, height / 2]");
        vsapi->freeNode(d.node);
        return;
    }

    d.bottom = vsapi->mapGetInt(in, "bottom", 0, &err);
    if (err)
        d.bottom = 0;

    else if (d.bottom > vi->height / 2) {
        vsapi->mapSetError(out, "Balance: bottom must be in [0, height / 2]");
        vsapi->freeNode(d.node);
        return;
    }

    d.left = vsapi->mapGetInt(in, "left", 0, &err);
    if (err)
        d.left = 0;

    else if (d.left > vi->width / 2) {
        vsapi->mapSetError(out, "Balance: left must be in [0, width / 2]");
        vsapi->freeNode(d.node);
        return;
    }

    d.right = vsapi->mapGetInt(in, "right", 0, &err);
    if (err)
        d.right = 0;

    else if (d.right > vi->width / 2) {
        vsapi->mapSetError(out, "Balance: right must be in [0, width / 2]");
        vsapi->freeNode(d.node);
        return;
    }

    d.plane = vsapi->mapGetInt(in, "plane", 0, &err);
    if (err)
        d.plane = 0;

    if (!vsh_isConstantVideoFormat(vi) || vi->format.sampleType != stFloat) {
        vsapi->mapSetError(out, "Balance: only constant format single float input");
        vsapi->freeNode(d.node);
        return;
    }

    data = (LinearRegressionData *)malloc(sizeof(d));
    *data = d;

    VSFilterDependency deps[] = {{d.node, rpStrictSpatial}};
    vsapi->createVideoFilter(out, "SinglePlane", vi, singlePlaneGetFrame, linearRegressionFree, fmParallel, deps, 1, data, core);
}

static void VS_CC multiPlaneCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi) {
    LinearRegressionData d;
    LinearRegressionData *data;
    int err;

    d.node = vsapi->mapGetNode(in, "clip", 0, 0);
    const VSVideoInfo *vi = vsapi->getVideoInfo(d.node);

    d.top = vsapi->mapGetInt(in, "top", 0, &err);
    if (err)
        d.top = 0;

    else if (d.top > vi->height / 2) {
        vsapi->mapSetError(out, "Balance: top must be in [0, height / 2]");
        vsapi->freeNode(d.node);
        return;
    }

    d.bottom = vsapi->mapGetInt(in, "bottom", 0, &err);
    if (err)
        d.bottom = 0;

    else if (d.bottom > vi->height / 2) {
        vsapi->mapSetError(out, "Balance: bottom must be in [0, height / 2]");
        vsapi->freeNode(d.node);
        return;
    }

    d.left = vsapi->mapGetInt(in, "left", 0, &err);
    if (err)
        d.left = 0;

    else if (d.left > vi->width / 2) {
        vsapi->mapSetError(out, "Balance: left must be in [0, width / 2]");
        vsapi->freeNode(d.node);
        return;
    }

    d.right = vsapi->mapGetInt(in, "right", 0, &err);
    if (err)
        d.right = 0;

    else if (d.right > vi->width / 2) {
        vsapi->mapSetError(out, "Balance: right must be in [0, width / 2]");
        vsapi->freeNode(d.node);
        return;
    }

    d.plane = vsapi->mapGetInt(in, "plane", 0, &err);
    if (err)
        d.plane = 0;

    if (!vsh_isConstantVideoFormat(vi) || vi->format.sampleType != stFloat) {
        vsapi->mapSetError(out, "Balance: only constant format single float input");
        vsapi->freeNode(d.node);
        return;
    }
    if (vi->format.subSamplingH > 0 || vi->format.subSamplingW > 0 || vi->format.numPlanes != 3) {
        vsapi->mapSetError(out, "Balance: mode 1 (multiple) requires three planes with no subsampling");
        vsapi->freeNode(d.node);
        return;
    }

    data = (LinearRegressionData *)malloc(sizeof(d));
    *data = d;

    VSFilterDependency deps[] = {{d.node, rpStrictSpatial}};
    vsapi->createVideoFilter(out, "MultiPlane", vi, multiPlaneGetFrame, linearRegressionFree, fmParallel, deps, 1, data, core);
}

static void VS_CC singlePlaneLimitedCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi) {
    LinearRegressionData d;
    LinearRegressionData *data;
    int err;

    d.node = vsapi->mapGetNode(in, "clip", 0, 0);
    const VSVideoInfo *vi = vsapi->getVideoInfo(d.node);

    d.top = vsapi->mapGetInt(in, "top", 0, &err);
    if (err)
        d.top = 0;

    else if (d.top > vi->height / 2) {
        vsapi->mapSetError(out, "Balance: top must be in [0, height / 2]");
        vsapi->freeNode(d.node);
        return;
    }

    d.bottom = vsapi->mapGetInt(in, "bottom", 0, &err);
    if (err)
        d.bottom = 0;

    else if (d.bottom > vi->height / 2) {
        vsapi->mapSetError(out, "Balance: bottom must be in [0, height / 2]");
        vsapi->freeNode(d.node);
        return;
    }

    d.left = vsapi->mapGetInt(in, "left", 0, &err);
    if (err)
        d.left = 0;

    else if (d.left > vi->width / 2) {
        vsapi->mapSetError(out, "Balance: left must be in [0, width / 2]");
        vsapi->freeNode(d.node);
        return;
    }

    d.right = vsapi->mapGetInt(in, "right", 0, &err);
    if (err)
        d.right = 0;

    else if (d.right > vi->width / 2) {
        vsapi->mapSetError(out, "Balance: right must be in [0, width / 2]");
        vsapi->freeNode(d.node);
        return;
    }

    d.plane = vsapi->mapGetInt(in, "plane", 0, &err);
    if (err)
        d.plane = 0;

    if (!vsh_isConstantVideoFormat(vi) || vi->format.sampleType != stFloat) {
        vsapi->mapSetError(out, "Balance: only constant format single float input");
        vsapi->freeNode(d.node);
        return;
    }

    d.ref_line_size = vsapi->mapGetInt(in, "ref_line_size", 0, &err) / 2;
    if (err)
        d.ref_line_size = 5;

    data = (LinearRegressionData *)malloc(sizeof(d));
    *data = d;

    VSFilterDependency deps[] = {{d.node, rpStrictSpatial}};
    vsapi->createVideoFilter(out, "SinglePlaneLimited", vi, singlePlaneLimitedGetFrame, linearRegressionFree, fmParallel, deps, 1, data, core);
}

static void VS_CC singlePlaneWeightedCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi) {
    LinearRegressionData d;
    LinearRegressionData *data;
    int err;

    d.node = vsapi->mapGetNode(in, "clip", 0, 0);
    const VSVideoInfo *vi = vsapi->getVideoInfo(d.node);

    d.top = vsapi->mapGetInt(in, "top", 0, &err);
    if (err)
        d.top = 0;

    else if (d.top > vi->height / 2) {
        vsapi->mapSetError(out, "Balance: top must be in [0, height / 2]");
        vsapi->freeNode(d.node);
        return;
    }

    d.bottom = vsapi->mapGetInt(in, "bottom", 0, &err);
    if (err)
        d.bottom = 0;

    else if (d.bottom > vi->height / 2) {
        vsapi->mapSetError(out, "Balance: bottom must be in [0, height / 2]");
        vsapi->freeNode(d.node);
        return;
    }

    d.left = vsapi->mapGetInt(in, "left", 0, &err);
    if (err)
        d.left = 0;

    else if (d.left > vi->width / 2) {
        vsapi->mapSetError(out, "Balance: left must be in [0, width / 2]");
        vsapi->freeNode(d.node);
        return;
    }

    d.right = vsapi->mapGetInt(in, "right", 0, &err);
    if (err)
        d.right = 0;

    else if (d.right > vi->width / 2) {
        vsapi->mapSetError(out, "Balance: right must be in [0, width / 2]");
        vsapi->freeNode(d.node);
        return;
    }

    d.plane = vsapi->mapGetInt(in, "plane", 0, &err);
    if (err)
        d.plane = 0;

    if (!vsh_isConstantVideoFormat(vi) || vi->format.sampleType != stFloat) {
        vsapi->mapSetError(out, "Balance: only constant format single float input");
        vsapi->freeNode(d.node);
        return;
    }

    d.sigmaS = vsapi->mapGetFloat(in, "sigmaS", 0, &err);
    if (err)
        d.sigmaS = 3.0;

    d.sigmaR = vsapi->mapGetFloat(in, "sigmaR", 0, &err);
    if (err)
        d.sigmaR = 0.5;

    d.ref_line_size = vsapi->mapGetInt(in, "ref_line_size", 0, &err) / 2;
    if (err)
        d.ref_line_size = 5;

    data = (LinearRegressionData *)malloc(sizeof(d));
    *data = d;

    VSFilterDependency deps[] = {{d.node, rpStrictSpatial}};
    vsapi->createVideoFilter(out, "SinglePlaneWeighted", vi, singlePlaneWeightedGetFrame, linearRegressionFree, fmParallel, deps, 1, data, core);
}

static void VS_CC singlePlaneDebugCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi) {
    LinearRegressionData d;
    LinearRegressionData *data;
    int err;

    d.node = vsapi->mapGetNode(in, "clip", 0, 0);
    const VSVideoInfo *vi = vsapi->getVideoInfo(d.node);

    d.top = vsapi->mapGetInt(in, "top", 0, &err);
    if (err)
        d.top = 0;

    else if (d.top > vi->height / 2) {
        vsapi->mapSetError(out, "Balance: top must be in [0, height / 2]");
        vsapi->freeNode(d.node);
        return;
    }

    d.bottom = vsapi->mapGetInt(in, "bottom", 0, &err);
    if (err)
        d.bottom = 0;

    else if (d.bottom > vi->height / 2) {
        vsapi->mapSetError(out, "Balance: bottom must be in [0, height / 2]");
        vsapi->freeNode(d.node);
        return;
    }

    d.left = vsapi->mapGetInt(in, "left", 0, &err);
    if (err)
        d.left = 0;

    else if (d.left > vi->width / 2) {
        vsapi->mapSetError(out, "Balance: left must be in [0, width / 2]");
        vsapi->freeNode(d.node);
        return;
    }

    d.right = vsapi->mapGetInt(in, "right", 0, &err);
    if (err)
        d.right = 0;

    else if (d.right > vi->width / 2) {
        vsapi->mapSetError(out, "Balance: right must be in [0, width / 2]");
        vsapi->freeNode(d.node);
        return;
    }

    d.plane = vsapi->mapGetInt(in, "plane", 0, &err);
    if (err)
        d.plane = 0;

    if (!vsh_isConstantVideoFormat(vi) || vi->format.sampleType != stFloat) {
        vsapi->mapSetError(out, "Balance: only constant format single float input");
        vsapi->freeNode(d.node);
        return;
    }

    data = (LinearRegressionData *)malloc(sizeof(d));
    *data = d;

    VSFilterDependency deps[] = {{d.node, rpStrictSpatial}};
    vsapi->createVideoFilter(out, "SinglePlaneDebug", vi, singlePlaneDebugGetFrame, linearRegressionFree, fmParallel, deps, 1, data, core);
}

VS_EXTERNAL_API(void) VapourSynthPluginInit2(VSPlugin *plugin, const VSPLUGINAPI *vspapi) {
    vspapi->configPlugin("ng.opusga.bore", "bore", "bore plugin", VS_MAKE_VERSION(1, 0), VAPOURSYNTH_API_VERSION, 0, plugin);
    vspapi->registerFunction("FixBrightness", "clip:vnode;top:int:opt;bottom:int:opt;left:int:opt;right:int:opt;ignore_mask:vnode:opt;thrlo:float:opt;thrhi:float:opt;step:int:opt;plane:int:opt;", "clip:vnode;", fixBrightnessCreate, NULL, plugin);
    vspapi->registerFunction("SinglePlane", "clip:vnode;top:int:opt;bottom:int:opt;left:int:opt;right:int:opt;plane:int:opt;", "clip:vnode;", singlePlaneCreate, NULL, plugin);
    vspapi->registerFunction("MultiPlane", "clip:vnode;top:int:opt;bottom:int:opt;left:int:opt;right:int:opt;plane:int:opt;", "clip:vnode;", multiPlaneCreate, NULL, plugin);
    vspapi->registerFunction("SinglePlaneLimited", "clip:vnode;top:int:opt;bottom:int:opt;left:int:opt;right:int:opt;ref_line_size:int:opt;plane:int:opt;", "clip:vnode;", singlePlaneLimitedCreate, NULL, plugin);
    vspapi->registerFunction("SinglePlaneWeighted", "clip:vnode;top:int:opt;bottom:int:opt;left:int:opt;right:int:opt;sigmaS:float:opt;sigmaR:float:opt;ref_line_size:int:opt;plane:int:opt;", "clip:vnode;", singlePlaneWeightedCreate, NULL, plugin);
    vspapi->registerFunction("SinglePlaneDebug", "clip:vnode;top:int:opt;bottom:int:opt;left:int:opt;right:int:opt;plane:int:opt;", "clip:vnode;", singlePlaneDebugCreate, NULL, plugin);
}

