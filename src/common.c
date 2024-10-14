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

#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>

void processRowSLR(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, ...)
{
    int sign = 1;
    if (row > h / 2)
        sign = -1;

    double* cur, * ref;
    int i;

    cur = malloc(sizeof(double) * w);
    ref = malloc(sizeof(double) * w);

    dstp += row * stride;
    for (i = 0; i < w; i++)
    {
        cur[i] = dstp[i];
        ref[i] = dstp[sign * stride + i];
    }

    double cov11, sumsq;

    double c1;

    int status = gsl_fit_mul(cur, 1, ref, 1, w, &c1, &cov11, &sumsq);

    if (!status && isfinite(c1))
    {
        // adjust each pixel
        for (i = 0; i < w; i++)
        {
            dstp[i] *= (float)c1;
        }
    }

    free(cur);
    free(ref);
}

void processColumnSLR(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, ...)
{
    int sign = 1;
    if (column > w / 2)
        sign = -1;

    double* cur, * ref;
    int i;

    cur = malloc(sizeof(double) * h);
    ref = malloc(sizeof(double) * h);

    dstp += column;
    for (i = 0; i < h; i++)
    {
        cur[i] = dstp[i * stride];
        ref[i] = dstp[sign + stride * i];
    }

    double cov11, sumsq;

    double c1;

    int status = gsl_fit_mul(cur, 1, ref, 1, h, &c1, &cov11, &sumsq);

    if (!status && isfinite(c1))
    {
        // adjust each pixel
        for (i = 0; i < h; i++)
        {
            dstp[i * stride] *= (float)c1;
        }
    }

    free(cur);
    free(ref);
}

void processRowSLRMasked(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, ...)
{
    va_list args;
    va_start(args, dstp);
    const unsigned char* imaskp = va_arg(args, const unsigned char*);
    const ptrdiff_t imaskstride = va_arg(args, ptrdiff_t);
    const int mask_dist = va_arg(args, int);
    va_end(args);

    int sign = 1;
    if (row > h / 2)
        sign = -1;

    double* cur, * ref;
    double* weights;
    int i;

    cur = malloc(sizeof(double) * w);
    ref = malloc(sizeof(double) * w);
    weights = malloc(sizeof(double) * w);

    dstp += row * stride;
    imaskp += row * imaskstride;
    for (i = 0; i < w; i++)
    {
        cur[i] = dstp[i];
        ref[i] = dstp[sign * stride + i];
        if (imaskp[i] < 128 && imaskp[sign * mask_dist * imaskstride + i] < 128)
            weights[i] = 1.0;
        else
            weights[i] = 0.0;
    }

    double cov11, sumsq;

    double c1;

    int status = gsl_fit_wmul(cur, 1, weights, 1, ref, 1, w, &c1, &cov11, &sumsq);

    if (!status && isfinite(c1))
    {
        // adjust each pixel
        for (i = 0; i < w; i++)
        {
            dstp[i] *= (float)c1;
        }
    }

    free(cur);
    free(ref);
    free(weights);
}

void processColumnSLRMasked(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, ...)
{
    va_list args;
    va_start(args, dstp);
    const unsigned char* imaskp = va_arg(args, const unsigned char*);
    const ptrdiff_t imaskstride = va_arg(args, ptrdiff_t);
    const int mask_dist = va_arg(args, int);
    va_end(args);

    int sign = 1;
    if (column > w / 2)
        sign = -1;

    double* cur, * ref;
    double* weights;
    int i;

    cur = malloc(sizeof(double) * h);
    ref = malloc(sizeof(double) * h);
    weights = malloc(sizeof(double) * h);

    imaskp += column;
    dstp += column;
    for (i = 0; i < h; i++)
    {
        cur[i] = dstp[i * stride];
        ref[i] = dstp[sign + stride * i];
        if (imaskp[i * imaskstride] < 128 && imaskp[sign * mask_dist + imaskstride * i] < 128)
            weights[i] = 1.0;
        else
            weights[i] = 0.0;
    }

    double cov11, sumsq;

    double c1;

    int status = gsl_fit_wmul(cur, 1, weights, 1, ref, 1, h, &c1, &cov11, &sumsq);

    if (!status && isfinite(c1))
    {
        // adjust each pixel
        for (i = 0; i < h; i++)
        {
            dstp[i * stride] *= (float)c1;
        }
    }

    free(cur);
    free(ref);
    free(weights);
}

void debugRowSLR(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, double** props)
{
    int sign = 1;
    if (row > h / 2)
        sign = -1;

    double* cur, * ref;
    int i;

    cur = malloc(sizeof(double) * w);
    ref = malloc(sizeof(double) * w);

    dstp += row * stride;
    for (i = 0; i < w; i++)
    {
        cur[i] = dstp[i];
        ref[i] = dstp[sign * stride + i];
    }

    double c1_cov11_sumsq[3];

    int status = gsl_fit_mul(cur, 1, ref, 1, w, &c1_cov11_sumsq[0], &c1_cov11_sumsq[1], &c1_cov11_sumsq[2]);

    if (status || !isfinite(c1_cov11_sumsq[0]))
    {
        c1_cov11_sumsq[0] = 1.0;
    }

    *props = c1_cov11_sumsq;

    free(cur);
    free(ref);
}

void debugColumnSLR(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, double** props)
{
    int sign = 1;
    if (column > w / 2)
        sign = -1;

    double* cur, * ref;
    int i;

    cur = malloc(sizeof(double) * h);
    ref = malloc(sizeof(double) * h);

    dstp += column;
    for (i = 0; i < h; i++)
    {
        cur[i] = dstp[i * stride];
        ref[i] = dstp[sign + stride * i];
    }

    double c1_cov11_sumsq[3];

    int status = gsl_fit_mul(cur, 1, ref, 1, h, &c1_cov11_sumsq[0], &c1_cov11_sumsq[1], &c1_cov11_sumsq[2]);

    if (status || !isfinite(c1_cov11_sumsq[0]))
    {
        c1_cov11_sumsq[0] = 1.0;
    }

    *props = c1_cov11_sumsq;

    free(cur);
    free(ref);
}

void debugRowSLRMasked(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, double** props,
    const unsigned char* __restrict imaskp, ptrdiff_t imaskstride, int mask_dist)
{
    int sign = 1;
    if (row > h / 2)
        sign = -1;

    double* cur, * ref;
    double* weights;
    int i;

    cur = malloc(sizeof(double) * w);
    ref = malloc(sizeof(double) * w);
    weights = malloc(sizeof(double) * w);

    dstp += row * stride;
    imaskp += row * imaskstride;
    for (i = 0; i < w; i++)
    {
        cur[i] = dstp[i];
        ref[i] = dstp[sign * stride + i];
        if (imaskp[i] < 128 && imaskp[sign * mask_dist * imaskstride + i] < 128)
            weights[i] = 1.0;
        else
            weights[i] = 0.0;
    }

    double c1_cov11_sumsq[3];

    int status = gsl_fit_wmul(cur, 1, weights, 1, ref, 1, w, &c1_cov11_sumsq[0], &c1_cov11_sumsq[1], &c1_cov11_sumsq[2]);

    if (status || !isfinite(c1_cov11_sumsq[0]))
    {
        c1_cov11_sumsq[0] = 1.0;
    }

    *props = c1_cov11_sumsq;

    free(cur);
    free(ref);
    free(weights);
}

void debugColumnSLRMasked(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, double** props,
    const unsigned char* __restrict imaskp, ptrdiff_t imaskstride, int mask_dist)
{
    int sign = 1;
    if (column > w / 2)
        sign = -1;

    double* cur, * ref;
    double* weights;
    int i;

    cur = malloc(sizeof(double) * h);
    ref = malloc(sizeof(double) * h);
    weights = malloc(sizeof(double) * h);

    imaskp += column;
    dstp += column;
    for (i = 0; i < h; i++)
    {
        cur[i] = dstp[i * stride];
        ref[i] = dstp[sign + stride * i];
        if (imaskp[i * imaskstride] < 128 && imaskp[sign * mask_dist + imaskstride * i] < 128)
            weights[i] = 1.0;
        else
            weights[i] = 0.0;
    }

    double c1_cov11_sumsq[3];

    int status = gsl_fit_wmul(cur, 1, weights, 1, ref, 1, h, &c1_cov11_sumsq[0], &c1_cov11_sumsq[1], &c1_cov11_sumsq[2]);

    if (status || !isfinite(c1_cov11_sumsq[0]))
    {
        c1_cov11_sumsq[0] = 1.0;
    }

    *props = c1_cov11_sumsq;

    free(cur);
    free(ref);
    free(weights);
}

void processRowMLR(int row, int w, int h, ptrdiff_t stride, float* dstp, float* dstp1, float* dstp2, float* dstp3)
{
    int i;
    int sign = 1;
    if (row > h / 2)
        sign = -1;

    dstp += stride * row;
    dstp1 += stride * row;
    dstp2 += stride * row;
    dstp3 += stride * row;
    gsl_vector* y = gsl_vector_alloc(w);
    gsl_matrix* x = gsl_matrix_alloc(w, 3);
    for (i = 0; i < w; i++)
    {
        gsl_vector_set(y, i, dstp[i + stride * sign]);
        gsl_matrix_set(x, i, 0, dstp1[i]);
        gsl_matrix_set(x, i, 1, dstp2[i]);
        gsl_matrix_set(x, i, 2, dstp3[i]);
    }

    gsl_multifit_linear_workspace* ws = gsl_multifit_linear_alloc(w, 3);
    double chisq;
    gsl_matrix* cov = gsl_matrix_alloc(3, 3);
    gsl_vector* b = gsl_vector_alloc(3);
    int status = gsl_multifit_linear(x, y, b, cov, &chisq, ws);

    if (!status)
    {
        // adjust each pixel
        for (i = 0; i < w; i++)
        {
            dstp[i] = (float)(gsl_vector_get(b, 0) * dstp1[i] + gsl_vector_get(b, 1) * dstp2[i] + gsl_vector_get(b, 2) * dstp3[i]);
        }
    }

    gsl_vector_free(b);
    gsl_matrix_free(cov);
    gsl_multifit_linear_free(ws);
    gsl_matrix_free(x);
    gsl_vector_free(y);
}

void processColumnMLR(int column, int w, int h, ptrdiff_t stride, float* dstp, float* dstp1, float* dstp2, float* dstp3)
{
    int i;
    int j;
    int sign = 1;
    if (column > w / 2)
        sign = -1;


    gsl_vector* y = gsl_vector_alloc(h);
    gsl_matrix* x = gsl_matrix_alloc(h, 3);
    for (i = 0; i < h; i++)
    {
        gsl_vector_set(y, i, dstp[i * stride + column + sign]);
        gsl_matrix_set(x, i, 0, dstp1[i * stride + column]);
        gsl_matrix_set(x, i, 1, dstp2[i * stride + column]);
        gsl_matrix_set(x, i, 2, dstp3[i * stride + column]);
    }
    gsl_multifit_linear_workspace* ws = gsl_multifit_linear_alloc(h, 3);
    double chisq;
    gsl_matrix* cov = gsl_matrix_alloc(3, 3);
    gsl_vector* b = gsl_vector_alloc(3);
    int status = gsl_multifit_linear(x, y, b, cov, &chisq, ws);

    if (!status)
    {
        // adjust each pixel
        for (i = 0; i < h; i++)
        {
            j = (int)(i * stride + column);
            dstp[j] = (float)(gsl_vector_get(b, 0) * dstp1[j] + gsl_vector_get(b, 1) * dstp2[j] + gsl_vector_get(b, 2) * dstp3[j]);
        }
    }

    gsl_vector_free(b);
    gsl_matrix_free(cov);
    gsl_multifit_linear_free(ws);
    gsl_matrix_free(x);
    gsl_vector_free(y);
}

void processRowSLRRef(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, ...)
{
    va_list args;
    va_start(args, dstp);
    const int ref_line_size = va_arg(args, int);
    va_end(args);

    int sign = 1;
    if (row > h / 2)
        sign = -1;

    int i;
    int start, stop;
    int status;

    double c1;
    double cov11, sumsq;

    double* cur, * ref;
    cur = malloc(sizeof(double) * w);
    ref = malloc(sizeof(double) * w);

    dstp += row * stride;
    for (i = 0; i < w; i++)
    {
        cur[i] = dstp[i];
        ref[i] = dstp[sign * stride + i];
    }

    for (i = 0; i < w; i++)
    {
        start = i - ref_line_size;
        stop = i + ref_line_size + 1;
        if (start < 0)
            start = 0;
        if (stop >= w)
            stop = w - 1;

        status = gsl_fit_mul(cur + start, 1, ref + start, 1, stop - start, &c1, &cov11, &sumsq);

        if (!status && isfinite(c1))
            dstp[i] *= (float)c1;
    }

    free(cur);
    free(ref);
}

void processColumnSLRRef(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, ...)
{
    va_list args;
    va_start(args, dstp);
    const int ref_line_size = va_arg(args, int);
    va_end(args);

    int sign = 1;
    if (column > w / 2)
        sign = -1;

    int i;
    int start, stop;
    int status;

    double c1;
    double cov11, sumsq;

    double* cur, * ref;
    cur = malloc(sizeof(double) * h);
    ref = malloc(sizeof(double) * h);

    dstp += column;
    for (i = 0; i < h; i++)
    {
        cur[i] = dstp[i * stride];
        ref[i] = dstp[sign + stride * i];
    }

    for (i = 0; i < h; i++)
    {
        start = i - ref_line_size;
        stop = i + ref_line_size + 1;
        if (start < 0)
            start = 0;
        if (stop >= h)
            stop = h - 1;

        status = gsl_fit_mul(cur + start, 1, ref + start, 1, stop - start, &c1, &cov11, &sumsq);

        if (!status && isfinite(c1))
            dstp[i * stride] *= (float)c1;
    }

    free(cur);
    free(ref);
}

void processRowSLRRefMasked(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, ...)
{
    va_list args;
    va_start(args, dstp);
    const int ref_line_size = va_arg(args, int);
    const unsigned char* imaskp = va_arg(args, const unsigned char*);
    const ptrdiff_t imaskstride = va_arg(args, ptrdiff_t);
    const int mask_dist = va_arg(args, int);
    va_end(args);

    int sign = 1;
    if (row > h / 2)
        sign = -1;

    int i;
    int start, stop;
    int status;

    double c1;
    double cov11, sumsq;

    double* cur, * ref;
    double* weights;

    cur = malloc(sizeof(double) * w);
    ref = malloc(sizeof(double) * w);
    weights = malloc(sizeof(double) * w);

    dstp += row * stride;
    imaskp += row * imaskstride;
    for (i = 0; i < w; i++)
    {
        cur[i] = dstp[i];
        ref[i] = dstp[sign * stride + i];
        if (imaskp[i] < 128 && imaskp[sign * mask_dist * imaskstride + i] < 128)
            weights[i] = 1.0;
        else
            weights[i] = 0.0;
    }

    for (i = 0; i < w; i++)
    {
        start = i - ref_line_size;
        stop = i + ref_line_size + 1;
        if (start < 0)
            start = 0;
        if (stop >= w)
            stop = w - 1;

        status = gsl_fit_wmul(cur + start, 1, weights + start, 1, ref + start, 1, stop - start, &c1, &cov11, &sumsq);

        if (!status && isfinite(c1))
            dstp[i] *= (float)c1;
    }

    free(cur);
    free(ref);
    free(weights);
}

void processColumnSLRRefMasked(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, ...)
{
    va_list args;
    va_start(args, dstp);
    const int ref_line_size = va_arg(args, int);
    const unsigned char* imaskp = va_arg(args, const unsigned char*);
    const ptrdiff_t imaskstride = va_arg(args, ptrdiff_t);
    const int mask_dist = va_arg(args, int);
    va_end(args);

    int sign = 1;
    if (column > w / 2)
        sign = -1;

    int i;
    int start, stop;
    int status;

    double c1;
    double cov11, sumsq;

    double* cur, * ref;
    double* weights;

    cur = malloc(sizeof(double) * h);
    ref = malloc(sizeof(double) * h);
    weights = malloc(sizeof(double) * h);

    imaskp += column;
    dstp += column;
    for (i = 0; i < h; i++)
    {
        cur[i] = dstp[i * stride];
        ref[i] = dstp[sign + stride * i];
        if (imaskp[i * imaskstride] < 128 && imaskp[sign * mask_dist + imaskstride * i] < 128)
            weights[i] = 1.0;
        else
            weights[i] = 1.0;
    }

    for (i = 0; i < h; i++)
    {
        start = i - ref_line_size;
        stop = i + ref_line_size + 1;
        if (start < 0)
            start = 0;
        if (stop >= h)
            stop = h - 1;

        status = gsl_fit_wmul(cur + start, 1, weights + start, 1, ref + start, 1, stop - start, &c1, &cov11, &sumsq);

        if (!status && isfinite(c1))
            dstp[i * stride] *= (float)c1;
    }

    free(cur);
    free(ref);
    free(weights);
}

void processRowWSLR(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, ...)
{
    va_list args;
    va_start(args, dstp);
    const int ref_line_size = va_arg(args, int);
    const double sigmaS = va_arg(args, double);
    const double sigmaR = va_arg(args, double);
    const double sigmaD = va_arg(args, double);
    va_end(args);

    int sign = 1;
    if (row > h / 2)
        sign = -1;

    int i;
    int k;
    int start, stop;
    int status;

    double c1;
    double cov11, sumsq;

    double* cur, * ref;
    cur = malloc(sizeof(double) * w);
    ref = malloc(sizeof(double) * w);

    dstp += row * stride;
    for (i = 0; i < w; i++)
    {
        cur[i] = dstp[i];
        ref[i] = dstp[sign * stride + i];
    }

    double w_s;
    double w_c;
    double w_d;

    for (i = 0; i < w; i++)
    {
        start = i - ref_line_size;
        stop = i + ref_line_size + 1;
        if (start < 0)
            start = 0;
        if (stop >= w)
            stop = w - 1;

        double* weights;
        weights = malloc(sizeof(double) * (stop - start));

        for (k = start; k < stop; k++)
        {
            w_s = exp(-((k - i) * (k - i)) / (sigmaS * sigmaS));
            w_c = exp(-((cur[k] - cur[i]) * (cur[k] - cur[i]) / (sigmaR * sigmaR)));
            w_d = exp(-(ref[k] / cur[k] - ref[i] / cur[i]) * (ref[k] / cur[k] - ref[i] / cur[i]) / (sigmaD * sigmaD));
            weights[k - start] = w_s * w_c * w_d;
            /* weights[k] = exp(-((j - i) * (j - i)) / (sigmaS * sigmaS) - ((cur[k] - cur[ref_line_size]) * (cur[k] - cur[ref_line_size])) / (sigmaR * sigmaR) - (ref[k] / cur[k] - cur[ref_line_size] / dstp[sign * stride + i]) * (ref[k] / cur[k] - cur[ref_line_size] / dstp[sign * stride + i]) / (sigmaD * sigmaD)); */
            /* weights[k] = exp(-((j - i) * (j - i)) / (sigmaS * sigmaS) - ((cur[k] - cur[ref_line_size]) * (cur[k] - cur[ref_line_size])) / (sigmaR * sigmaR)); */
            /* if (i == 1370) { */
            /*     dstp[j] = w_d;//(dstp[j] - dstp[i]) * (dstp[j] - dstp[i]); */
            /*     dstp[i] = 1.0; */
            /* } */
        }


        status = gsl_fit_wmul(cur + start, 1, weights, 1, ref + start, 1, stop - start, &c1, &cov11, &sumsq);

        if (!status && isfinite(c1))
            dstp[i] *= (float)c1;

        free(weights);
    }

    free(cur);
    free(ref);
}

void processColumnWSLR(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, ...)
{
    va_list args;
    va_start(args, dstp);
    const int ref_line_size = va_arg(args, int);
    const double sigmaS = va_arg(args, double);
    const double sigmaR = va_arg(args, double);
    const double sigmaD = va_arg(args, double);
    va_end(args);

    int sign = 1;
    if (column > w / 2)
        sign = -1;

    int i;
    int k;
    int start, stop;
    int status;

    double c1;
    double cov11, sumsq;

    double* cur, * ref;
    cur = malloc(sizeof(double) * h);
    ref = malloc(sizeof(double) * h);

    dstp += column;
    for (i = 0; i < h; i++)
    {
        cur[i] = dstp[i * stride];
        ref[i] = dstp[sign + stride * i];
    }

    for (i = 0; i < h; i++)
    {
        start = i - ref_line_size;
        stop = i + ref_line_size + 1;
        if (start < 0)
            start = 0;
        if (stop >= h)
            stop = h - 1;

        double* weights;
        weights = malloc(sizeof(double) * (stop - start));

        for (k = start; k < stop; k++)
        {
            weights[k - start] = exp(-((k - i) * (k - i)) / (sigmaS * sigmaS) - ((cur[k] - dstp[i * stride]) * (cur[k] - dstp[i * stride])) / (sigmaR * sigmaR) - (ref[k] / cur[k] - dstp[i] / dstp[sign * stride + i]) * (ref[k] / cur[k] - dstp[i] / dstp[sign * stride + i]) / (sigmaD * sigmaD));
        }

        status = gsl_fit_wmul(cur + start, 1, weights, 1, ref + start, 1, stop - start, &c1, &cov11, &sumsq);

        if (!status && isfinite(c1))
            dstp[i * stride] *= (float)c1;

        free(weights);
    }

    free(cur);
    free(ref);
}


void processRowWSLRMasked(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, ...)
{
    va_list args;
    va_start(args, dstp);
    const int ref_line_size = va_arg(args, int);
    const double sigmaS = va_arg(args, double);
    const double sigmaR = va_arg(args, double);
    const double sigmaD = va_arg(args, double);
    const unsigned char* imaskp = va_arg(args, const unsigned char*);
    const ptrdiff_t imaskstride = va_arg(args, ptrdiff_t);
    const int mask_dist = va_arg(args, int);
    va_end(args);

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

    double* cur, * ref;
    cur = malloc(sizeof(double) * w);
    ref = malloc(sizeof(double) * w);

    dstp += row * stride;
    imaskp += row * imaskstride;
    for (i = 0; i < w; i++)
    {
        cur[i] = dstp[i];
        ref[i] = dstp[sign * stride + i];
    }

    double w_s;
    double w_c;
    double w_d;

    double* weights;
    weights = malloc(sizeof(double) * (2 * ref_line_size));

    sign *= mask_dist;

    for (i = 0; i < w; i++)
    {
        start = i - ref_line_size;
        stop = i + ref_line_size + 1;
        if (start < 0)
            start = 0;
        if (stop >= w)
            stop = w - 1;

        for (k = 0; k < stop - start; k++)
        {
            j = k + start;
            if (imaskp[j] < 128 && imaskp[sign * imaskstride + j] < 128)
            {
                w_s = exp(-((j - i) * (j - i)) / (sigmaS * sigmaS));
                w_c = exp(-((ref[j] - ref[i]) * (ref[j] - ref[i]) / (sigmaR * sigmaR)));
                w_d = exp(-(ref[j] / cur[j] - ref[i] / cur[i]) * (ref[j] / cur[j] - ref[i] / cur[i]) / (sigmaD * sigmaD));
                weights[k] = w_s * w_c * w_d;
            }
            else
                weights[k] = 0.0;
            /* if (i == 1370) { */
            /*     dstp[j] = w_d;//(dstp[j] - dstp[i]) * (dstp[j] - dstp[i]); */
            /*     dstp[i] = 1.0; */
            /* } */
        }

        status = gsl_fit_wmul(cur + start, 1, weights, 1, ref + start, 1, stop - start, &c1, &cov11, &sumsq);

        if (!status && isfinite(c1))
            dstp[i] *= (float)c1;

    }

    free(weights);
    free(cur);
    free(ref);
}

void processColumnWSLRMasked(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, ...)
{
    va_list args;
    va_start(args, dstp);
    const int ref_line_size = va_arg(args, int);
    const double sigmaS = va_arg(args, double);
    const double sigmaR = va_arg(args, double);
    const double sigmaD = va_arg(args, double);
    const unsigned char* imaskp = va_arg(args, const unsigned char*);
    const ptrdiff_t imaskstride = va_arg(args, ptrdiff_t);
    const int mask_dist = va_arg(args, int);
    va_end(args);

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

    double* cur, * ref;
    cur = malloc(sizeof(double) * h);
    ref = malloc(sizeof(double) * h);

    double* weights;
    weights = malloc(sizeof(double) * (2 * ref_line_size));

    imaskp += column;
    dstp += column;
    for (i = 0; i < h; i++)
    {
        cur[i] = dstp[i * stride];
        ref[i] = dstp[sign + stride * i];
    }

    double w_s, w_c, w_d;

    sign *= mask_dist;

    for (i = 0; i < h; i++)
    {
        start = i - ref_line_size;
        stop = i + ref_line_size + 1;
        if (start < 0)
            start = 0;
        if (stop >= h)
            stop = h - 1;

        for (k = 0; k < stop - start; k++)
        {
            j = k + start;
            if (imaskp[j * imaskstride] < 128 && imaskp[sign + imaskstride * j] < 128)
            {
                /* cur_cur = cur[i]; */
                /* cur_ref = cur[j]; */
                /* ref_cur = ref[i]; */
                /* ref_ref = ref[j]; */
                /* w_s = exp(-((j - i) * (j - i)) / (sigmaS * sigmaS)); */
                /* w_c = exp(-((ref_ref - ref_cur) * (ref_ref - ref_cur) + (cur_ref - cur_cur) * (cur_ref - cur_cur)) / (sigmaR * sigmaR)); */
                /* w_d = exp(-(ref_ref / cur_ref - ref_cur / cur_cur) * (ref_ref / cur_ref - ref_cur / cur_cur) / (sigmaD * sigmaD)); */
                w_s = exp(-((j - i) * (j - i)) / (sigmaS * sigmaS));
                w_c = exp(-((ref[j] - ref[i]) * (ref[j] - ref[i]) / (sigmaR * sigmaR)));
                w_d = exp(-(ref[j] / cur[j] - ref[i] / cur[i]) * (ref[j] / cur[j] - ref[i] / cur[i]) / (sigmaD * sigmaD));
                weights[k] = w_s * w_c * w_d;
            }
            else
                weights[k] = 0.0;
            /* if (i == 880) { */
            /*     dstp[j * stride] = weights[k]; */
            /*     dstp[i * stride] = 0.5; */
            /* } */
        }

        status = gsl_fit_wmul(cur + start, 1, weights, 1, ref + start, 1, stop - start, &c1, &cov11, &sumsq);

        if (!status && isfinite(c1))
            dstp[i * stride] *= (float)c1;

    }

    free(weights);
    free(cur);
    free(ref);
}
