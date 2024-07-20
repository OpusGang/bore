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

typedef enum LinRegMode {
    LINREG_MODE_SINGLE = 1,
    LINREG_MODE_MULTI = 2,
    LINREG_MODE_SINGLE_LIMITED = 3,
    LINREG_MODE_SINGLE_WEIGHTED = 4,
    LINREG_MODE_SINGLE_DEBUG = 5
} LinRegMode;

typedef struct {
    VSNode *node;
    VSNode *ignore_mask;
    int plane;
    int top;
    int bottom;
    int left;
    int right;
    int ref_line_size;
    double sigmaS;
    double sigmaR;
    double sigmaD;
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
static void processRowSLRMasked(int row, int w, int h, ptrdiff_t stride, float *dstp, const unsigned char * restrict imaskp, ptrdiff_t imaskstride, int mask_dist) {
    int sign = 1;
    if (row > h / 2)
        sign = -1;

    double *cur, *ref;
    double *weights;
    int i;

    cur = malloc(sizeof(double) * w);
    ref = malloc(sizeof(double) * w);
    weights = malloc(sizeof(double) * w);

    dstp += row * stride;
    imaskp += row * imaskstride;
    for (i = 0; i < w; i++) {
        cur[i] = dstp[i];
        ref[i] = dstp[sign * stride + i];
        if (imaskp[i] < 128 && imaskp[sign * mask_dist * imaskstride + i] < 128)
            weights[i] = 1.0;
        else
            weights[i] = 0.0;
    }
    
    double cov11, sumsq;

    double c1;

    const double *const_cur = cur;
    const double *const_ref = ref;
    const double *const_weights = weights;

    int status = gsl_fit_wmul(const_cur, 1, const_weights, 1, const_ref, 1, w, &c1, &cov11, &sumsq);

    if (!status && isfinite(c1)) {
        // adjust each pixel
        for (i = 0; i < w; i++) {
            dstp[i] *= c1;
        }
    }
    
    free(cur);
    free(ref);
    free(weights);
}

static void processColumnSLRMasked(int column, int w, int h, ptrdiff_t stride, float *dstp, const unsigned char * restrict imaskp, ptrdiff_t imaskstride, int mask_dist) {
    int sign = 1;
    if (column > w / 2)
        sign = -1;

    double *cur, *ref;
    double *weights;
    int i;

    cur = malloc(sizeof(double) * h);
    ref = malloc(sizeof(double) * h);
    weights = malloc(sizeof(double) * h);

    imaskp += column;
    dstp += column;
    for (i = 0; i < h; i++) {
        cur[i] = dstp[i * stride];
        ref[i] = dstp[sign + stride * i];
        if (imaskp[i * imaskstride] < 128 && imaskp[sign * mask_dist + imaskstride * i] < 128)
            weights[i] = 1.0;
        else
            weights[i] = 0.0;
    }

    double cov11, sumsq;

    double c1;

    const double *const_cur = cur;
    const double *const_ref = ref;
    const double *const_weights = weights;

    int status = gsl_fit_wmul(const_cur, 1, const_weights, 1, const_ref, 1, h, &c1, &cov11, &sumsq);

    if (!status && isfinite(c1)) {
        // adjust each pixel
        for (i = 0; i < h; i++) {
            dstp[i * stride] *= c1;
        }
    }
    
    free(cur);
    free(ref);
    free(weights);
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

    if (!status || isfinite(c1)) {
        c1 = 0;
    }
    vsapi->mapSetFloat(vsapi->getFramePropertiesRW(dst), "BoreAdjustment", c1, maAppend);
    vsapi->mapSetFloat(vsapi->getFramePropertiesRW(dst), "BoreCovariance", cov11, maAppend);
    vsapi->mapSetFloat(vsapi->getFramePropertiesRW(dst), "BoreSumSquares", sumsq, maAppend);

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
    vsapi->mapSetFloat(vsapi->getFramePropertiesRW(dst), "BoreAdjustment", c1, maAppend);
    vsapi->mapSetFloat(vsapi->getFramePropertiesRW(dst), "BoreCovariance", cov11, maAppend);
    vsapi->mapSetFloat(vsapi->getFramePropertiesRW(dst), "BoreSumSquares", sumsq, maAppend);

    free(cur);
    free(ref);
}

static void debugRowSLRMasked(int row, int w, int h, ptrdiff_t stride, float *dstp, VSFrame *dst, const VSAPI *vsapi, const unsigned char * restrict imaskp, ptrdiff_t imaskstride, int mask_dist) {
    int sign = 1;
    if (row > h / 2)
        sign = -1;

    double *cur, *ref;
    double *weights;
    int i;

    cur = malloc(sizeof(double) * w);
    ref = malloc(sizeof(double) * w);
    weights = malloc(sizeof(double) * w);

    dstp += row * stride;
    imaskp += row * imaskstride;
    for (i = 0; i < w; i++) {
        cur[i] = dstp[i];
        ref[i] = dstp[sign * stride + i];
        if (imaskp[i] < 128 && imaskp[sign * mask_dist * imaskstride + i] < 128)
            weights[i] = 1.0;
        else
            weights[i] = 0.0;
    }

    double cov11, sumsq;

    double c1;

    const double *const_cur = cur;
    const double *const_ref = ref;
    const double *const_weights = weights;

    int status = gsl_fit_wmul(const_cur, 1, const_weights, 1, const_ref, 1, w, &c1, &cov11, &sumsq);

    if (!status || isfinite(c1)) {
        c1 = 0;
    }
    vsapi->mapSetFloat(vsapi->getFramePropertiesRW(dst), "BoreAdjustment", c1, maAppend);
    vsapi->mapSetFloat(vsapi->getFramePropertiesRW(dst), "BoreCovariance", cov11, maAppend);
    vsapi->mapSetFloat(vsapi->getFramePropertiesRW(dst), "BoreSumSquares", sumsq, maAppend);

    free(cur);
    free(ref);
    free(weights);
}

static void debugColumnSLRMasked(int column, int w, int h, ptrdiff_t stride, float *dstp, VSFrame *dst, const VSAPI *vsapi, const unsigned char * restrict imaskp, ptrdiff_t imaskstride, int mask_dist) {
    int sign = 1;
    if (column > w / 2)
        sign = -1;

    double *cur, *ref;
    double *weights;
    int i;

    cur = malloc(sizeof(double) * h);
    ref = malloc(sizeof(double) * h);
    weights = malloc(sizeof(double) * h);

    imaskp += column;
    dstp += column;
    for (i = 0; i < h; i++) {
        cur[i] = dstp[i * stride];
        ref[i] = dstp[sign + stride * i];
        if (imaskp[i * imaskstride] < 128 && imaskp[sign * mask_dist + imaskstride * i] < 128)
            weights[i] = 1.0;
        else
            weights[i] = 0.0;
    }

    double cov11, sumsq;

    double c1;

    const double *const_cur = cur;
    const double *const_ref = ref;
    const double *const_weights = weights;

    int status = gsl_fit_wmul(const_cur, 1, const_weights, 1, const_ref, 1, h, &c1, &cov11, &sumsq);

    if (!status || isfinite(c1)) {
        c1 = 0;
    }
    vsapi->mapSetFloat(vsapi->getFramePropertiesRW(dst), "BoreAdjustment", c1, maAppend);
    vsapi->mapSetFloat(vsapi->getFramePropertiesRW(dst), "BoreCovariance", cov11, maAppend);
    vsapi->mapSetFloat(vsapi->getFramePropertiesRW(dst), "BoreSumSquares", sumsq, maAppend);

    free(cur);
    free(ref);
    free(weights);
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

    const double *const_cur = cur;
    const double *const_ref = ref;

    for (i = 0; i < w; i++) {
        start = i - ref_line_size;
        stop = i + ref_line_size + 1;
        if (start < 0)
            start = 0;
        if (stop >= w)
            stop = w - 1;

        status = gsl_fit_mul(const_cur + start, 1, const_ref + start, 1, stop - start, &c1, &cov11, &sumsq);

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
    cur = malloc(sizeof(double) * h);
    ref = malloc(sizeof(double) * h);

    dstp += column;
    for (i = 0; i < h; i++) {
        cur[i] = dstp[i * stride];
        ref[i] = dstp[sign + stride * i];
    }

    const double *const_cur = cur;
    const double *const_ref = ref;

    for (i = 0; i < h; i++) {
        start = i - ref_line_size;
        stop = i + ref_line_size + 1;
        if (start < 0)
            start = 0;
        if (stop >= h)
            stop = h - 1;

        status = gsl_fit_mul(const_cur + start, 1, const_ref + start, 1, stop - start, &c1, &cov11, &sumsq);

        if (!status && isfinite(c1)) 
            dstp[i * stride] *= c1;
    }

    free(cur);
    free(ref);
}

static void processRowSLRRefMasked(int row, int w, int h, ptrdiff_t stride, float *dstp, int ref_line_size, const unsigned char * restrict imaskp, ptrdiff_t imaskstride, int mask_dist) {
    int sign = 1;
    if (row > h / 2)
        sign = -1;

    int i;
    int start, stop;
    int status;

    double c1;
    double cov11, sumsq;

    double *cur, *ref;
    double *weights;

    cur = malloc(sizeof(double) * w);
    ref = malloc(sizeof(double) * w);
    weights = malloc(sizeof(double) * w);

    dstp += row * stride;
    imaskp += row * imaskstride;
    for (i = 0; i < w; i++) {
        cur[i] = dstp[i];
        ref[i] = dstp[sign * stride + i];
        if (imaskp[i] < 128 && imaskp[sign * mask_dist * imaskstride + i] < 128)
            weights[i] = 1.0;
        else
            weights[i] = 0.0;
    }

    const double *const_cur = cur;
    const double *const_ref = ref;
    const double *const_weights = weights;

    for (i = 0; i < w; i++) {
        start = i - ref_line_size;
        stop = i + ref_line_size + 1;
        if (start < 0)
            start = 0;
        if (stop >= w)
            stop = w - 1;

        status = gsl_fit_wmul(const_cur + start, 1, const_weights + start, 1, const_ref + start, 1, stop - start, &c1, &cov11, &sumsq);

        if (!status && isfinite(c1)) 
            dstp[i] *= c1;
    }

    free(cur);
    free(ref);
    free(weights);
}

static void processColumnSLRRefMasked(int column, int w, int h, ptrdiff_t stride, float *dstp, int ref_line_size, const unsigned char * restrict imaskp, ptrdiff_t imaskstride, int mask_dist) {
    int sign = 1;
    if (column > w / 2)
        sign = -1;

    int i;
    int start, stop;
    int status;

    double c1;
    double cov11, sumsq;

    double *cur, *ref;
    double *weights;

    cur = malloc(sizeof(double) * h);
    ref = malloc(sizeof(double) * h);
    weights = malloc(sizeof(double) * h);

    imaskp += column;
    dstp += column;
    for (i = 0; i < h; i++) {
        cur[i] = dstp[i * stride];
        ref[i] = dstp[sign + stride * i];
        if (imaskp[i * imaskstride] < 128 && imaskp[sign * mask_dist + imaskstride * i] < 128)
            weights[i] = 1.0;
        else
            weights[i] = 1.0;
    }

    const double *const_cur = cur;
    const double *const_ref = ref;
    const double *const_weights = weights;

    for (i = 0; i < h; i++) {
        start = i - ref_line_size;
        stop = i + ref_line_size + 1;
        if (start < 0)
            start = 0;
        if (stop >= h)
            stop = h - 1;

        status = gsl_fit_wmul(const_cur + start, 1, const_weights + start, 1, const_ref + start, 1, stop - start, &c1, &cov11, &sumsq);

        if (!status && isfinite(c1)) 
            dstp[i * stride] *= c1;
    }

    free(cur);
    free(ref);
    free(weights);
}

static void processRowWSLR(int row, int w, int h, ptrdiff_t stride, float *dstp, int ref_line_size, double sigmaS, double sigmaR, double sigmaD) {
    int sign = 1;
    if (row > h / 2)
        sign = -1;

    int i;
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

    double w_s;
    double w_c;
    double w_d;
    const double *const_cur = cur;
    const double *const_ref = ref;

    for (i = 0; i < w; i++) {
        start = i - ref_line_size;
        stop = i + ref_line_size + 1;
        if (start < 0)
            start = 0;
        if (stop >= w)
            stop = w - 1;

        double *weights;
        weights = malloc(sizeof(double) * (stop - start));

        for (k = start; k < stop; k++) {
            w_s = exp(-((k - i) * (k - i)) / (sigmaS * sigmaS));
            w_c = exp(-((const_cur[k] - const_cur[i]) * (const_cur[k] - const_cur[i]) / (sigmaR * sigmaR)));
            w_d = exp(-(const_ref[k] / const_cur[k] - const_ref[i] / const_cur[i]) * (const_ref[k] / const_cur[k] - const_ref[i] / const_cur[i]) / (sigmaD * sigmaD));
            weights[k - start] = w_s * w_c * w_d;
            /* weights[k] = exp(-((j - i) * (j - i)) / (sigmaS * sigmaS) - ((const_cur[k] - const_cur[ref_line_size]) * (const_cur[k] - const_cur[ref_line_size])) / (sigmaR * sigmaR) - (const_ref[k] / const_cur[k] - const_cur[ref_line_size] / dstp[sign * stride + i]) * (const_ref[k] / const_cur[k] - const_cur[ref_line_size] / dstp[sign * stride + i]) / (sigmaD * sigmaD)); */
            /* weights[k] = exp(-((j - i) * (j - i)) / (sigmaS * sigmaS) - ((const_cur[k] - const_cur[ref_line_size]) * (const_cur[k] - const_cur[ref_line_size])) / (sigmaR * sigmaR)); */
            /* if (i == 1370) { */
            /*     dstp[j] = w_d;//(dstp[j] - dstp[i]) * (dstp[j] - dstp[i]); */
            /*     dstp[i] = 1.0; */
            /* } */
        }

        const double *const_weights = weights;

        status = gsl_fit_wmul(const_cur + start, 1, const_weights, 1, const_ref + start, 1, stop - start, &c1, &cov11, &sumsq);

        if (!status && isfinite(c1)) 
            dstp[i] *= c1;

        free(weights);
    }

    free(cur);
    free(ref);
}

static void processColumnWSLR(int column, int w, int h, ptrdiff_t stride, float *dstp, int ref_line_size, double sigmaS, double sigmaR, double sigmaD) {
    int sign = 1;
    if (column > w / 2)
        sign = -1;

    int i;
    int k;
    int start, stop;
    int status;

    double c1;
    double cov11, sumsq;

    double *cur, *ref;
    cur = malloc(sizeof(double) * h);
    ref = malloc(sizeof(double) * h);

    dstp += column;
    for (i = 0; i < h; i++) {
        cur[i] = dstp[i * stride];
        ref[i] = dstp[sign + stride * i];
    }
    
    const double *const_cur = cur;
    const double *const_ref = ref;

    for (i = 0; i < h; i++) {
        start = i - ref_line_size;
        stop = i + ref_line_size + 1;
        if (start < 0)
            start = 0;
        if (stop >= h)
            stop = h - 1;
        
        double *weights;
        weights = malloc(sizeof(double) * (stop - start));

        for (k = start; k < stop; k++) {
            weights[k - start] = exp(-((k - i) * (k - i)) / (sigmaS * sigmaS) - ((const_cur[k] - dstp[i * stride]) * (const_cur[k] - dstp[i * stride])) / (sigmaR * sigmaR) - (const_ref[k] / const_cur[k] - dstp[i] / dstp[sign * stride + i]) * (const_ref[k] / const_cur[k] - dstp[i] / dstp[sign * stride + i]) / (sigmaD * sigmaD));
        }

        const double *const_weights = weights;

        status = gsl_fit_wmul(const_cur + start, 1, const_weights, 1, const_ref + start, 1, stop - start, &c1, &cov11, &sumsq);

        if (!status && isfinite(c1)) 
            dstp[i * stride] *= c1;

        free(weights);
    }

    free(cur);
    free(ref);
}


static void processRowWSLRMasked(int row, int w, int h, ptrdiff_t stride, float *dstp, int ref_line_size, double sigmaS, double sigmaR, double sigmaD, const unsigned char * restrict imaskp, ptrdiff_t imaskstride, int mask_dist) {
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
    imaskp += row * imaskstride;
    for (i = 0; i < w; i++) {
        cur[i] = dstp[i];
        ref[i] = dstp[sign * stride + i];
    }

    double w_s;
    double w_c;
    double w_d;

    const double *const_cur = cur;
    const double *const_ref = ref;
    double *weights;
    weights = malloc(sizeof(double) * (2 * ref_line_size));

    sign *= mask_dist;

    for (i = 0; i < w; i++) {
        start = i - ref_line_size;
        stop = i + ref_line_size + 1;
        if (start < 0)
            start = 0;
        if (stop >= w)
            stop = w - 1;

        for (k = 0; k < stop - start; k++) {
            j = k + start;
            if (imaskp[j] < 128 && imaskp[sign * imaskstride + j] < 128) {
                w_s = exp(-((j - i) * (j - i)) / (sigmaS * sigmaS));
                w_c = exp(-((const_ref[j] - const_ref[i]) * (const_ref[j] - const_ref[i]) / (sigmaR * sigmaR)));
                w_d = exp(-(const_ref[j] / const_cur[j] - const_ref[i] / const_cur[i]) * (const_ref[j] / const_cur[j] - const_ref[i] / const_cur[i]) / (sigmaD * sigmaD));
                weights[k] = w_s * w_c * w_d;
            } else
                weights[k] = 0.0;
            /* if (i == 1370) { */
            /*     dstp[j] = w_d;//(dstp[j] - dstp[i]) * (dstp[j] - dstp[i]); */
            /*     dstp[i] = 1.0; */
            /* } */
        }

        status = gsl_fit_wmul(const_cur + start, 1, weights, 1, const_ref + start, 1, stop - start, &c1, &cov11, &sumsq);

        if (!status && isfinite(c1)) 
            dstp[i] *= c1;

    }

    free(weights);
    free(cur);
    free(ref);
}

static void processColumnWSLRMasked(int column, int w, int h, ptrdiff_t stride, float *dstp, int ref_line_size, double sigmaS, double sigmaR, double sigmaD, const unsigned char * restrict imaskp, ptrdiff_t imaskstride, int mask_dist) {
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
    cur = malloc(sizeof(double) * h);
    ref = malloc(sizeof(double) * h);
        
    double *weights;
    weights = malloc(sizeof(double) * (2 * ref_line_size));

    const double *const_cur = cur;
    const double *const_ref = ref;

    imaskp += column;
    dstp += column;
    for (i = 0; i < h; i++) {
        cur[i] = dstp[i * stride];
        ref[i] = dstp[sign + stride * i];
    }

    double w_s, w_c, w_d;

    float cur_cur, cur_ref, ref_ref, ref_cur;
    
    sign *= mask_dist;

    for (i = 0; i < h; i++) {
        start = i - ref_line_size;
        stop = i + ref_line_size + 1;
        if (start < 0)
            start = 0;
        if (stop >= h)
            stop = h - 1;

        for (k = 0; k < stop - start; k++) {
            j = k + start;
            if (imaskp[j * imaskstride] < 128 && imaskp[sign + imaskstride * j] < 128) {
                /* cur_cur = const_cur[i]; */
                /* cur_ref = const_cur[j]; */
                /* ref_cur = const_ref[i]; */
                /* ref_ref = const_ref[j]; */
                /* w_s = exp(-((j - i) * (j - i)) / (sigmaS * sigmaS)); */
                /* w_c = exp(-((ref_ref - ref_cur) * (ref_ref - ref_cur) + (cur_ref - cur_cur) * (cur_ref - cur_cur)) / (sigmaR * sigmaR)); */
                /* w_d = exp(-(ref_ref / cur_ref - ref_cur / cur_cur) * (ref_ref / cur_ref - ref_cur / cur_cur) / (sigmaD * sigmaD)); */
                w_s = exp(-((j - i) * (j - i)) / (sigmaS * sigmaS));
                w_c = exp(-((const_ref[j] - const_ref[i]) * (const_ref[j] - const_ref[i]) / (sigmaR * sigmaR)));
                w_d = exp(-(const_ref[j] / const_cur[j] - const_ref[i] / const_cur[i]) * (const_ref[j] / const_cur[j] - const_ref[i] / const_cur[i]) / (sigmaD * sigmaD));
                weights[k] = w_s * w_c * w_d;
            } else
                weights[k] = 0.0;
            /* if (i == 880) { */
            /*     dstp[j * stride] = weights[k]; */
            /*     dstp[i * stride] = 0.5; */
            /* } */
        }

        status = gsl_fit_wmul(const_cur + start, 1, weights, 1, const_ref + start, 1, stop - start, &c1, &cov11, &sumsq);

        if (!status && isfinite(c1)) 
            dstp[i * stride] *= c1;

    }

    free(weights);
    free(cur);
    free(ref);
}

static const VSFrame *VS_CC singlePlaneGetFrame(int n, int activationReason, void *instanceData, void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) {
    LinearRegressionData *d = (LinearRegressionData *)instanceData;

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
            ptrdiff_t imaskstride = vsapi->getStride(ignore_mask, 0);
            imaskp = vsapi->getReadPtr(ignore_mask, d->plane);
            if (d->top != 0) {
                for (int row = d->top - 1; row > -1; --row)
                    processRowSLRMasked(row, w, h, stride, dstp, imaskp, imaskstride, d->top - row);
            }
            if (d->bottom != 0) {
                for (int row = h - d->bottom; row < h; ++row)
                    processRowSLRMasked(row, w, h, stride, dstp, imaskp, imaskstride, d->bottom + row - h + 1);
            }
            if (d->left != 0) {
                for (int column = d->left - 1; column > -1; --column)
                    processColumnSLRMasked(column, w, h, stride, dstp, imaskp, imaskstride, d->left - column);
            }
            if (d->right != 0) {
                for (int column = w - d->right; column < w; ++column)
                    processColumnSLRMasked(column, w, h, stride, dstp, imaskp, imaskstride, d->right + column - w + 1);
            }
            vsapi->freeFrame(ignore_mask);
        } else {
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
            ptrdiff_t imaskstride = vsapi->getStride(ignore_mask, 0);
            imaskp = vsapi->getReadPtr(ignore_mask, d->plane);
            if (d->top != 0) {
                for (int row = d->top - 1; row > -1; --row)
                    processRowSLRRefMasked(row, w, h, stride, dstp, d->ref_line_size, imaskp, imaskstride, d->top - row);
            }
            if (d->bottom != 0) {
                for (int row = h - d->bottom; row < h; ++row)
                    processRowSLRRefMasked(row, w, h, stride, dstp, d->ref_line_size, imaskp, imaskstride, d->bottom + row - h + 1);
            }
            if (d->left != 0) {
                for (int column = d->left - 1; column > -1; --column)
                    processColumnSLRRefMasked(column, w, h, stride, dstp, d->ref_line_size, imaskp, imaskstride, d->left - column);
            }
            if (d->right != 0) {
                for (int column = w - d->right; column < w; ++column)
                    processColumnSLRRefMasked(column, w, h, stride, dstp, d->ref_line_size, imaskp, imaskstride, d->right + column - w + 1);
            }
            vsapi->freeFrame(ignore_mask);
        } else {
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
            ptrdiff_t imaskstride = vsapi->getStride(ignore_mask, 0);
            if (d->top != 0) {
                for (int row = d->top - 1; row > -1; --row)
                    processRowWSLRMasked(row, w, h, stride, dstp, d->ref_line_size, d->sigmaS, d->sigmaR, d->sigmaD, imaskp, imaskstride, d->top - row);
            }
            if (d->bottom != 0) {
                for (int row = h - d->bottom; row < h; ++row)
                    processRowWSLRMasked(row, w, h, stride, dstp, d->ref_line_size, d->sigmaS, d->sigmaR, d->sigmaD, imaskp, imaskstride, d->bottom + row - h + 1);
            }
            if (d->left != 0) {
                for (int column = d->left - 1; column > -1; --column)
                    processColumnWSLRMasked(column, w, h, stride, dstp, d->ref_line_size, d->sigmaS, d->sigmaR, d->sigmaD, imaskp, imaskstride, d->left - column);
            }
            if (d->right != 0) {
                for (int column = w - d->right; column < w; ++column)
                    processColumnWSLRMasked(column, w, h, stride, dstp, d->ref_line_size, d->sigmaS, d->sigmaR, d->sigmaD, imaskp, imaskstride, d->right + column - w + 1);
            }
            vsapi->freeFrame(ignore_mask);
        } else {
            if (d->top != 0) {
                for (int row = d->top - 1; row > -1; --row)
                    processRowWSLR(row, w, h, stride, dstp, d->ref_line_size, d->sigmaS, d->sigmaR, d->sigmaD);
            }
            if (d->bottom != 0) {
                for (int row = h - d->bottom; row < h; ++row)
                    processRowWSLR(row, w, h, stride, dstp, d->ref_line_size, d->sigmaS, d->sigmaR, d->sigmaD);
            }
            if (d->left != 0) {
                for (int column = d->left - 1; column > -1; --column)
                    processColumnWSLR(column, w, h, stride, dstp, d->ref_line_size, d->sigmaS, d->sigmaR, d->sigmaD);
            }
            if (d->right != 0) {
                for (int column = w - d->right; column < w; ++column)
                    processColumnWSLR(column, w, h, stride, dstp, d->ref_line_size, d->sigmaS, d->sigmaR, d->sigmaD);
            }
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
            ptrdiff_t imaskstride = vsapi->getStride(ignore_mask, 0);
            imaskp = vsapi->getReadPtr(ignore_mask, d->plane);
            if (d->top != 0) {
                for (int row = d->top - 1; row > -1; --row)
                    debugRowSLRMasked(row, w, h, stride, dstp, dst, vsapi, imaskp, imaskstride, d->top - row);
            }
            if (d->bottom != 0) {
                for (int row = h - d->bottom; row < h; ++row)
                    debugRowSLRMasked(row, w, h, stride, dstp, dst, vsapi, imaskp, imaskstride, d->bottom + row - h + 1);
            }
            if (d->left != 0) {
                for (int column = d->left - 1; column > -1; --column)
                    debugColumnSLRMasked(column, w, h, stride, dstp, dst, vsapi, imaskp, imaskstride, d->left - column);
            }
            if (d->right != 0) {
                for (int column = w - d->right; column < w; ++column)
                    debugColumnSLRMasked(column, w, h, stride, dstp, dst, vsapi, imaskp, imaskstride, d->right + column - w + 1);
            }
            vsapi->freeFrame(ignore_mask);
        } else {
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
        }
        
        vsapi->freeFrame(src);

        return dst;
    }

    return NULL;
}

static void VS_CC linearRegressionFree(void *instanceData, VSCore *core, const VSAPI *vsapi) {
    LinearRegressionData *d = (LinearRegressionData *)instanceData;
    vsapi->freeNode(d->node);
    vsapi->freeNode(d->ignore_mask);
    free(d);
}

static void VS_CC linearRegressionCreate(const VSMap *in, VSMap *out, void *userData, VSCore *core, const VSAPI *vsapi) {
    LinRegMode mode = (enum LinRegMode) userData;
    LinearRegressionData d;
    LinearRegressionData *data;
    int err;

    d.node = vsapi->mapGetNode(in, "clip", 0, 0);
    const VSVideoInfo *vi = vsapi->getVideoInfo(d.node);

    d.ignore_mask = vsapi->mapGetNode(in, "ignore_mask", 0, &err);
    if (err)
        d.ignore_mask = NULL;
    else {
        vsapi->freeNode(d.ignore_mask);
    }

    d.top = vsapi->mapGetInt(in, "top", 0, &err);
    if (err)
        d.top = 0;

    else if (d.top > vi->height / 2) {
        vsapi->mapSetError(out, "SinglePlaneDebug: top must be in [0, height / 2]");
        vsapi->freeNode(d.node);
        vsapi->freeNode(d.ignore_mask);
        return;
    }

    d.bottom = vsapi->mapGetInt(in, "bottom", 0, &err);
    if (err)
        d.bottom = 0;

    else if (d.bottom > vi->height / 2) {
        vsapi->mapSetError(out, "SinglePlaneDebug: bottom must be in [0, height / 2]");
        vsapi->freeNode(d.node);
        vsapi->freeNode(d.ignore_mask);
        return;
    }

    d.left = vsapi->mapGetInt(in, "left", 0, &err);
    if (err)
        d.left = 0;

    else if (d.left > vi->width / 2) {
        vsapi->mapSetError(out, "SinglePlaneDebug: left must be in [0, width / 2]");
        vsapi->freeNode(d.node);
        vsapi->freeNode(d.ignore_mask);
        return;
    }

    d.right = vsapi->mapGetInt(in, "right", 0, &err);
    if (err)
        d.right = 0;

    else if (d.right > vi->width / 2) {
        vsapi->mapSetError(out, "SinglePlaneDebug: right must be in [0, width / 2]");
        vsapi->freeNode(d.node);
        vsapi->freeNode(d.ignore_mask);
        return;
    }

    d.plane = vsapi->mapGetInt(in, "plane", 0, &err);
    if (err)
        d.plane = 0;

    if (!vsh_isConstantVideoFormat(vi) || vi->format.sampleType != stFloat) {
        vsapi->mapSetError(out, "SinglePlaneDebug: only constant format single float input");
        vsapi->freeNode(d.node);
        vsapi->freeNode(d.ignore_mask);
        return;
    }

    if (mode == LINREG_MODE_SINGLE_LIMITED || mode == LINREG_MODE_SINGLE_WEIGHTED) {
        d.ref_line_size = vsapi->mapGetInt(in, "ref_line_size", 0, &err);
        if (err)
            d.ref_line_size = 100;
    }

    if (mode == LINREG_MODE_SINGLE_WEIGHTED) {
        d.sigmaS = vsapi->mapGetFloat(in, "sigmaS", 0, &err);
        if (err)
            d.sigmaS = 50.0;

        d.sigmaR = vsapi->mapGetFloat(in, "sigmaR", 0, &err);
        if (err)
            d.sigmaR = 0.5;

        d.sigmaD = vsapi->mapGetFloat(in, "sigmaD", 0, &err);
        if (err)
            d.sigmaD = 1.5;
    }

    data = (LinearRegressionData *)malloc(sizeof(d));
    *data = d;

    VSFilterDependency deps[] = {{d.node, rpStrictSpatial}};
    if (mode == LINREG_MODE_SINGLE)
        vsapi->createVideoFilter(out, "SinglePlane", vi, singlePlaneGetFrame, linearRegressionFree, fmParallel, deps, 1, data, core);
    else if (mode == LINREG_MODE_MULTI)
        vsapi->createVideoFilter(out, "MultiPlane", vi, multiPlaneGetFrame, linearRegressionFree, fmParallel, deps, 1, data, core);
    else if (mode == LINREG_MODE_SINGLE_LIMITED)
        vsapi->createVideoFilter(out, "SinglePlaneLimited", vi, singlePlaneLimitedGetFrame, linearRegressionFree, fmParallel, deps, 1, data, core);
    else if (mode == LINREG_MODE_SINGLE_WEIGHTED)
        vsapi->createVideoFilter(out, "SinglePlaneWeighted", vi, singlePlaneWeightedGetFrame, linearRegressionFree, fmParallel, deps, 1, data, core);
    else // (mode == LINREG_MODE_SINGLE_DEBUG)
        vsapi->createVideoFilter(out, "SinglePlaneDebug", vi, singlePlaneDebugGetFrame, linearRegressionFree, fmParallel, deps, 1, data, core);
}

VS_EXTERNAL_API(void) VapourSynthPluginInit2(VSPlugin *plugin, const VSPLUGINAPI *vspapi) {
    vspapi->configPlugin("ng.opusga.bore", "bore", "bore plugin", VS_MAKE_VERSION(1, 0), VAPOURSYNTH_API_VERSION, 0, plugin);
    vspapi->registerFunction("SinglePlane", "clip:vnode;left:int:opt;right:int:opt;top:int:opt;bottom:int:opt;ignore_mask:vnode:opt;plane:int:opt;", "clip:vnode;", linearRegressionCreate, (void *)(LINREG_MODE_SINGLE), plugin);
    vspapi->registerFunction("MultiPlane", "clip:vnode;left:int:opt;right:int:opt;top:int:opt;bottom:int:opt;ignore_mask:vnode:opt;plane:int:opt;", "clip:vnode;", linearRegressionCreate, (void *)(LINREG_MODE_MULTI), plugin);
    vspapi->registerFunction("SinglePlaneLimited", "clip:vnode;left:int:opt;right:int:opt;top:int:opt;bottom:int:opt;ignore_mask:vnode:opt;ref_line_size:int:opt;plane:int:opt;", "clip:vnode;", linearRegressionCreate, (void *)(LINREG_MODE_SINGLE_LIMITED), plugin);
    vspapi->registerFunction("SinglePlaneWeighted", "clip:vnode;left:int:opt;right:int:opt;top:int:opt;bottom:int:opt;ignore_mask:vnode:opt;sigmaS:float:opt;sigmaR:float:opt;sigmaD:float:opt;ref_line_size:int:opt;plane:int:opt;", "clip:vnode;", linearRegressionCreate, (void *)(LINREG_MODE_SINGLE_WEIGHTED), plugin);
    vspapi->registerFunction("SinglePlaneDebug", "clip:vnode;left:int:opt;right:int:opt;top:int:opt;bottom:int:opt;ignore_mask:vnode:opt;plane:int:opt;", "clip:vnode;", linearRegressionCreate, (void *)(LINREG_MODE_SINGLE_DEBUG), plugin);
}

