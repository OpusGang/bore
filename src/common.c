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

#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>

void boreInit(void)
{
    // GSL's default error handler calls abort(), which would take down the host
    // application (e.g. on an allocation failure inside gsl_multifit_linear).
    // All GSL return codes are checked, so report errors through them instead.
    gsl_set_error_handler_off();
}

struct BoreScratch
{
    double* cur;
    double* ref;
    double* weights;

    // SinglePlaneWeighted only
    double* ratio;
    int spatial_radius;
    double* spatial_exponent; // -(d * d) / (sigmaS * sigmaS), d in [-spatial_radius, spatial_radius]

    // MultiPlane only
    gsl_matrix* x;
    gsl_vector* y;
    gsl_vector* b;
    gsl_matrix* cov;
    gsl_multifit_linear_workspace* ws;
};

BoreScratch* boreScratchAlloc(const LinearRegressionData* d, int w, int h)
{
    const size_t n = (size_t)(w > h ? w : h);
    BoreScratch* s = calloc(1, sizeof(BoreScratch));
    if (!s)
        return NULL;

    if (d->mode == LINREG_MODE_MULTI)
    {
        s->x = gsl_matrix_alloc(n, 3);
        s->y = gsl_vector_alloc(n);
        s->b = gsl_vector_alloc(3);
        s->cov = gsl_matrix_alloc(3, 3);
        s->ws = gsl_multifit_linear_alloc(n, 3);
        if (!s->x || !s->y || !s->b || !s->cov || !s->ws)
        {
            boreScratchFree(s);
            return NULL;
        }
    }
    else
    {
        s->cur = malloc(sizeof(double) * n);
        s->ref = malloc(sizeof(double) * n);
        s->weights = malloc(sizeof(double) * n);
        if (!s->cur || !s->ref || !s->weights)
        {
            boreScratchFree(s);
            return NULL;
        }
    }

    if (d->mode == LINREG_MODE_SINGLE_WEIGHTED)
    {
        // a window never spans more than n samples, so |k - i| <= n - 1
        const int radius = d->ref_line_size < (int)n - 1 ? d->ref_line_size : (int)n - 1;
        const double sigmaS = d->sigmaS;

        s->spatial_radius = radius;
        s->ratio = malloc(sizeof(double) * n);
        s->spatial_exponent = malloc(sizeof(double) * (2 * (size_t)radius + 1));
        if (!s->ratio || !s->spatial_exponent)
        {
            boreScratchFree(s);
            return NULL;
        }

        for (int k = -radius; k <= radius; k++)
            s->spatial_exponent[k + radius] = -(k * k) / (sigmaS * sigmaS);
    }

    return s;
}

void boreScratchFree(BoreScratch* s)
{
    if (!s)
        return;

    free(s->cur);
    free(s->ref);
    free(s->weights);
    free(s->ratio);
    free(s->spatial_exponent);
    if (s->x)
        gsl_matrix_free(s->x);
    if (s->y)
        gsl_vector_free(s->y);
    if (s->b)
        gsl_vector_free(s->b);
    if (s->cov)
        gsl_matrix_free(s->cov);
    if (s->ws)
        gsl_multifit_linear_free(s->ws);
    free(s);
}

void processRowSLR(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p)
{
    int sign = 1;
    if (row >= h / 2)
        sign = -1;

    double* const cur = p->scratch->cur;
    double* const ref = p->scratch->ref;
    int i;


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
}

void processColumnSLR(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p)
{
    int sign = 1;
    if (column >= w / 2)
        sign = -1;

    double* const cur = p->scratch->cur;
    double* const ref = p->scratch->ref;
    int i;


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
}

void processRowSLRMasked(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p)
{
    const float* wmaskp = p->wmaskp;
    const ptrdiff_t wmaskstride = p->wmaskstride;
    const int mask_dist = p->mask_dist;

    int sign = 1;
    if (row >= h / 2)
        sign = -1;

    double* const cur = p->scratch->cur;
    double* const ref = p->scratch->ref;
    double* const weights = p->scratch->weights;
    int i;


    dstp += row * stride;
    wmaskp += row * wmaskstride;
    for (i = 0; i < w; i++)
    {
        cur[i] = dstp[i];
        ref[i] = dstp[sign * stride + i];
        weights[i] = wmaskp[i] * wmaskp[sign * mask_dist * wmaskstride + i];
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
}

void processColumnSLRMasked(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p)
{
    const float* wmaskp = p->wmaskp;
    const ptrdiff_t wmaskstride = p->wmaskstride;
    const int mask_dist = p->mask_dist;

    int sign = 1;
    if (column >= w / 2)
        sign = -1;

    double* const cur = p->scratch->cur;
    double* const ref = p->scratch->ref;
    double* const weights = p->scratch->weights;
    int i;


    wmaskp += column;
    dstp += column;
    for (i = 0; i < h; i++)
    {
        cur[i] = dstp[i * stride];
        ref[i] = dstp[sign + stride * i];
        weights[i] = wmaskp[i * wmaskstride] * wmaskp[sign * mask_dist + wmaskstride * i];
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
}

void debugRowSLR(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p)
{
    double* c1_cov11_sumsq = p->c1_cov11_sumsq;
    int sign = 1;
    if (row >= h / 2)
        sign = -1;

    double* const cur = p->scratch->cur;
    double* const ref = p->scratch->ref;
    int i;


    dstp += row * stride;
    for (i = 0; i < w; i++)
    {
        cur[i] = dstp[i];
        ref[i] = dstp[sign * stride + i];
    }

    int status = gsl_fit_mul(cur, 1, ref, 1, w, &c1_cov11_sumsq[0], &c1_cov11_sumsq[1], &c1_cov11_sumsq[2]);

    if (status || !isfinite(c1_cov11_sumsq[0]))
    {
        c1_cov11_sumsq[0] = 1.0;
    }
}

void debugColumnSLR(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p)
{
    double* c1_cov11_sumsq = p->c1_cov11_sumsq;
    int sign = 1;
    if (column >= w / 2)
        sign = -1;

    double* const cur = p->scratch->cur;
    double* const ref = p->scratch->ref;
    int i;


    dstp += column;
    for (i = 0; i < h; i++)
    {
        cur[i] = dstp[i * stride];
        ref[i] = dstp[sign + stride * i];
    }

    int status = gsl_fit_mul(cur, 1, ref, 1, h, &c1_cov11_sumsq[0], &c1_cov11_sumsq[1], &c1_cov11_sumsq[2]);

    if (status || !isfinite(c1_cov11_sumsq[0]))
    {
        c1_cov11_sumsq[0] = 1.0;
    }
}

void debugRowSLRMasked(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p)
{
    double* c1_cov11_sumsq = p->c1_cov11_sumsq;
    const float* wmaskp = p->wmaskp;
    const ptrdiff_t wmaskstride = p->wmaskstride;
    const int mask_dist = p->mask_dist;
    int sign = 1;
    if (row >= h / 2)
        sign = -1;

    double* const cur = p->scratch->cur;
    double* const ref = p->scratch->ref;
    double* const weights = p->scratch->weights;
    int i;


    dstp += row * stride;
    wmaskp += row * wmaskstride;
    for (i = 0; i < w; i++)
    {
        cur[i] = dstp[i];
        ref[i] = dstp[sign * stride + i];
        weights[i] = wmaskp[i] * wmaskp[sign * mask_dist * wmaskstride + i];
    }

    int status = gsl_fit_wmul(cur, 1, weights, 1, ref, 1, w, &c1_cov11_sumsq[0], &c1_cov11_sumsq[1], &c1_cov11_sumsq[2]);

    if (status || !isfinite(c1_cov11_sumsq[0]))
    {
        c1_cov11_sumsq[0] = 1.0;
    }
}

void debugColumnSLRMasked(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p)
{
    double* c1_cov11_sumsq = p->c1_cov11_sumsq;
    const float* wmaskp = p->wmaskp;
    const ptrdiff_t wmaskstride = p->wmaskstride;
    const int mask_dist = p->mask_dist;
    int sign = 1;
    if (column >= w / 2)
        sign = -1;

    double* const cur = p->scratch->cur;
    double* const ref = p->scratch->ref;
    double* const weights = p->scratch->weights;
    int i;


    wmaskp += column;
    dstp += column;
    for (i = 0; i < h; i++)
    {
        cur[i] = dstp[i * stride];
        ref[i] = dstp[sign + stride * i];
        weights[i] = wmaskp[i * wmaskstride] * wmaskp[sign * mask_dist + wmaskstride * i];
    }

    int status = gsl_fit_wmul(cur, 1, weights, 1, ref, 1, h, &c1_cov11_sumsq[0], &c1_cov11_sumsq[1], &c1_cov11_sumsq[2]);

    if (status || !isfinite(c1_cov11_sumsq[0]))
    {
        c1_cov11_sumsq[0] = 1.0;
    }
}

void processRowMLR(int row, int w, int h, ptrdiff_t stride, float* dstp, float* dstp1, float* dstp2, float* dstp3, BoreScratch* scratch)
{
    int i;
    int sign = 1;
    if (row >= h / 2)
        sign = -1;

    dstp += stride * row;
    dstp1 += stride * row;
    dstp2 += stride * row;
    dstp3 += stride * row;
    gsl_vector_view y = gsl_vector_subvector(scratch->y, 0, w);
    gsl_matrix_view x = gsl_matrix_submatrix(scratch->x, 0, 0, w, 3);
    for (i = 0; i < w; i++)
    {
        gsl_vector_set(&y.vector, i, dstp[i + stride * sign]);
        gsl_matrix_set(&x.matrix, i, 0, dstp1[i]);
        gsl_matrix_set(&x.matrix, i, 1, dstp2[i]);
        gsl_matrix_set(&x.matrix, i, 2, dstp3[i]);
    }

    double chisq;
    gsl_vector* b = scratch->b;
    int status = gsl_multifit_linear(&x.matrix, &y.vector, b, scratch->cov, &chisq, scratch->ws);

    if (!status)
    {
        // adjust each pixel
        for (i = 0; i < w; i++)
        {
            dstp[i] = (float)(gsl_vector_get(b, 0) * dstp1[i] + gsl_vector_get(b, 1) * dstp2[i] + gsl_vector_get(b, 2) * dstp3[i]);
        }
    }
}

void processColumnMLR(int column, int w, int h, ptrdiff_t stride, float* dstp, float* dstp1, float* dstp2, float* dstp3, BoreScratch* scratch)
{
    int i;
    int j;
    int sign = 1;
    if (column >= w / 2)
        sign = -1;

    gsl_vector_view y = gsl_vector_subvector(scratch->y, 0, h);
    gsl_matrix_view x = gsl_matrix_submatrix(scratch->x, 0, 0, h, 3);
    for (i = 0; i < h; i++)
    {
        gsl_vector_set(&y.vector, i, dstp[i * stride + column + sign]);
        gsl_matrix_set(&x.matrix, i, 0, dstp1[i * stride + column]);
        gsl_matrix_set(&x.matrix, i, 1, dstp2[i * stride + column]);
        gsl_matrix_set(&x.matrix, i, 2, dstp3[i * stride + column]);
    }

    double chisq;
    gsl_vector* b = scratch->b;
    int status = gsl_multifit_linear(&x.matrix, &y.vector, b, scratch->cov, &chisq, scratch->ws);

    if (!status)
    {
        // adjust each pixel
        for (i = 0; i < h; i++)
        {
            j = (int)(i * stride + column);
            dstp[j] = (float)(gsl_vector_get(b, 0) * dstp1[j] + gsl_vector_get(b, 1) * dstp2[j] + gsl_vector_get(b, 2) * dstp3[j]);
        }
    }
}

// Per-pixel fit of ref ~= c1 * cur over the window [i - r, i + r] of a line, i.e.
// c1 = sum(w * x * y) / sum(w * x * x) (w = 1 if weights is NULL), the same quantity
// gsl_fit_mul / gsl_fit_wmul compute. The window sums are updated incrementally, so a
// line costs O(n) instead of O(n * r); results can differ from GSL in the last bits.
// Windows containing non-finite terms fall back to GSL to keep its exact behaviour.
static void slidingWindowFit(const double* __restrict cur, const double* __restrict ref, const double* __restrict weights,
    int n, int r, float* __restrict dstp, ptrdiff_t dst_step)
{
    double sxy = 0.0, sxx = 0.0;
    int nonzero = 0;   // terms with w * x * x != 0 in the window (all zero -> fit undefined)
    int nonfinite = 0; // terms that are not finite (not included in the sums)
    int lo = 0, hi = 0;

    for (int i = 0; i < n; i++)
    {
        int start = i - r;
        int stop = i + r + 1;
        if (start < 0)
            start = 0;
        if (stop > n)
            stop = n;

        for (; hi < stop; hi++)
        {
            const double wx = weights ? weights[hi] * cur[hi] : cur[hi];
            const double txy = wx * ref[hi];
            const double txx = wx * cur[hi];
            if (isfinite(txy) && isfinite(txx))
            {
                sxy += txy;
                sxx += txx;
                nonzero += txx != 0.0;
            }
            else
                nonfinite++;
        }
        for (; lo < start; lo++)
        {
            const double wx = weights ? weights[lo] * cur[lo] : cur[lo];
            const double txy = wx * ref[lo];
            const double txx = wx * cur[lo];
            if (isfinite(txy) && isfinite(txx))
            {
                sxy -= txy;
                sxx -= txx;
                nonzero -= txx != 0.0;
            }
            else
                nonfinite--;
        }

        double c1;
        if (nonfinite)
        {
            double cov11, sumsq;
            int status = weights
                ? gsl_fit_wmul(cur + start, 1, weights + start, 1, ref + start, 1, stop - start, &c1, &cov11, &sumsq)
                : gsl_fit_mul(cur + start, 1, ref + start, 1, stop - start, &c1, &cov11, &sumsq);
            if (status)
                continue;
        }
        else if (nonzero)
            c1 = sxy / sxx;
        else
            continue;

        if (isfinite(c1))
            dstp[i * dst_step] *= (float)c1;
    }
}

void processRowSLRRef(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p)
{
    const int ref_line_size = p->ref_line_size;

    int sign = 1;
    if (row >= h / 2)
        sign = -1;

    int i;

    double* const cur = p->scratch->cur;
    double* const ref = p->scratch->ref;

    dstp += row * stride;
    for (i = 0; i < w; i++)
    {
        cur[i] = dstp[i];
        ref[i] = dstp[sign * stride + i];
    }

    slidingWindowFit(cur, ref, NULL, w, ref_line_size, dstp, 1);
}

void processColumnSLRRef(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p)
{
    const int ref_line_size = p->ref_line_size;

    int sign = 1;
    if (column >= w / 2)
        sign = -1;

    int i;

    double* const cur = p->scratch->cur;
    double* const ref = p->scratch->ref;

    dstp += column;
    for (i = 0; i < h; i++)
    {
        cur[i] = dstp[i * stride];
        ref[i] = dstp[sign + stride * i];
    }

    slidingWindowFit(cur, ref, NULL, h, ref_line_size, dstp, stride);
}

void processRowSLRRefMasked(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p)
{
    const int ref_line_size = p->ref_line_size;
    const float* wmaskp = p->wmaskp;
    const ptrdiff_t wmaskstride = p->wmaskstride;
    const int mask_dist = p->mask_dist;

    int sign = 1;
    if (row >= h / 2)
        sign = -1;

    int i;

    double* const cur = p->scratch->cur;
    double* const ref = p->scratch->ref;
    double* const weights = p->scratch->weights;

    dstp += row * stride;
    wmaskp += row * wmaskstride;
    for (i = 0; i < w; i++)
    {
        cur[i] = dstp[i];
        ref[i] = dstp[sign * stride + i];
        weights[i] = wmaskp[i] * wmaskp[sign * mask_dist * wmaskstride + i];
    }

    slidingWindowFit(cur, ref, weights, w, ref_line_size, dstp, 1);
}

void processColumnSLRRefMasked(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p)
{
    const int ref_line_size = p->ref_line_size;
    const float* wmaskp = p->wmaskp;
    const ptrdiff_t wmaskstride = p->wmaskstride;
    const int mask_dist = p->mask_dist;

    int sign = 1;
    if (column >= w / 2)
        sign = -1;

    int i;

    double* const cur = p->scratch->cur;
    double* const ref = p->scratch->ref;
    double* const weights = p->scratch->weights;

    wmaskp += column;
    dstp += column;
    for (i = 0; i < h; i++)
    {
        cur[i] = dstp[i * stride];
        ref[i] = dstp[sign + stride * i];
        weights[i] = wmaskp[i * wmaskstride] * wmaskp[sign * mask_dist + wmaskstride * i];
    }

    slidingWindowFit(cur, ref, weights, h, ref_line_size, dstp, stride);
}

// Fills the per-line ratio buffer ref[i] / cur[i] used by the weighted fits.
static void computeRatios(const double* __restrict cur, const double* __restrict ref, double* __restrict ratio, int n)
{
    for (int i = 0; i < n; i++)
        ratio[i] = ref[i] / cur[i];
}

void processRowWSLR(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p)
{
    const int ref_line_size = p->ref_line_size;
    const double sigmaR2 = p->sigmaR * p->sigmaR;
    const double sigmaD2 = p->sigmaD * p->sigmaD;
    // spatial_exponent[d] = -(d * d) / (sigmaS * sigmaS) for d in [-spatial_radius, spatial_radius]
    const double* const spatial_exponent = p->scratch->spatial_exponent + p->scratch->spatial_radius;

    int sign = 1;
    if (row >= h / 2)
        sign = -1;

    int i;
    int k;
    int start, stop;
    double c1;

    double* const cur = p->scratch->cur;
    double* const ref = p->scratch->ref;
    double* const ratio = p->scratch->ratio;

    dstp += row * stride;
    for (i = 0; i < w; i++)
    {
        cur[i] = dstp[i];
        ref[i] = dstp[sign * stride + i];
    }
    computeRatios(cur, ref, ratio, w);

    for (i = 0; i < w; i++)
    {
        start = i - ref_line_size;
        stop = i + ref_line_size + 1;
        if (start < 0)
            start = 0;
        if (stop > w)
            stop = w;

        double sxy = 0.0, sxx = 0.0;
        for (k = start; k < stop; k++)
        {
            // exp(a) * exp(b) * exp(c) evaluated as exp(a + b + c)
            const double wt = exp(spatial_exponent[k - i] - ((cur[k] - cur[i]) * (cur[k] - cur[i])) / sigmaR2 - (ratio[k] - ratio[i]) * (ratio[k] - ratio[i]) / sigmaD2);
            sxy += wt * cur[k] * ref[k];
            sxx += wt * cur[k] * cur[k];
        }

        // weighted least squares through the origin (what gsl_fit_wmul computes)
        c1 = sxy / sxx;

        if (isfinite(c1))
            dstp[i] *= (float)c1;
    }
}

void processColumnWSLR(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p)
{
    const int ref_line_size = p->ref_line_size;
    const double sigmaR2 = p->sigmaR * p->sigmaR;
    const double sigmaD2 = p->sigmaD * p->sigmaD;
    // spatial_exponent[d] = -(d * d) / (sigmaS * sigmaS) for d in [-spatial_radius, spatial_radius]
    const double* const spatial_exponent = p->scratch->spatial_exponent + p->scratch->spatial_radius;

    int sign = 1;
    if (column >= w / 2)
        sign = -1;

    int i;
    int k;
    int start, stop;
    double c1;

    double* const cur = p->scratch->cur;
    double* const ref = p->scratch->ref;
    double* const ratio = p->scratch->ratio;

    dstp += column;
    for (i = 0; i < h; i++)
    {
        cur[i] = dstp[i * stride];
        ref[i] = dstp[sign + stride * i];
    }
    computeRatios(cur, ref, ratio, h);

    for (i = 0; i < h; i++)
    {
        start = i - ref_line_size;
        stop = i + ref_line_size + 1;
        if (start < 0)
            start = 0;
        if (stop > h)
            stop = h;

        double sxy = 0.0, sxx = 0.0;
        for (k = start; k < stop; k++)
        {
            const double wt = exp(spatial_exponent[k - i] - ((cur[k] - cur[i]) * (cur[k] - cur[i])) / sigmaR2 - (ratio[k] - ratio[i]) * (ratio[k] - ratio[i]) / sigmaD2);
            sxy += wt * cur[k] * ref[k];
            sxx += wt * cur[k] * cur[k];
        }

        // weighted least squares through the origin (what gsl_fit_wmul computes)
        c1 = sxy / sxx;

        if (isfinite(c1))
            dstp[i * stride] *= (float)c1;
    }
}

void processRowWSLRMasked(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p)
{
    const int ref_line_size = p->ref_line_size;
    const double sigmaR2 = p->sigmaR * p->sigmaR;
    const double sigmaD2 = p->sigmaD * p->sigmaD;
    const double* const spatial_exponent = p->scratch->spatial_exponent + p->scratch->spatial_radius;
    const float* wmaskp = p->wmaskp;
    const ptrdiff_t wmaskstride = p->wmaskstride;
    const int mask_dist = p->mask_dist;

    int sign = 1;
    if (row >= h / 2)
        sign = -1;

    int i;
    int j;
    int k;
    int start, stop;
    double c1;

    double* const cur = p->scratch->cur;
    double* const ref = p->scratch->ref;
    double* const ratio = p->scratch->ratio;

    dstp += row * stride;
    wmaskp += row * wmaskstride;
    for (i = 0; i < w; i++)
    {
        cur[i] = dstp[i];
        ref[i] = dstp[sign * stride + i];
    }
    computeRatios(cur, ref, ratio, w);

    sign *= mask_dist;

    for (i = 0; i < w; i++)
    {
        start = i - ref_line_size;
        stop = i + ref_line_size + 1;
        if (start < 0)
            start = 0;
        if (stop > w)
            stop = w;

        double sxy = 0.0, sxx = 0.0;
        for (k = 0; k < stop - start; k++)
        {
            j = k + start;

            // exp(a) * exp(b) * exp(c) evaluated as exp(a + b + c)
            const double wt = exp(spatial_exponent[j - i] - ((ref[j] - ref[i]) * (ref[j] - ref[i])) / sigmaR2 - (ratio[j] - ratio[i]) * (ratio[j] - ratio[i]) / sigmaD2) * wmaskp[j] * wmaskp[sign * wmaskstride + j];
            sxy += wt * cur[j] * ref[j];
            sxx += wt * cur[j] * cur[j];
        }

        // weighted least squares through the origin (what gsl_fit_wmul computes)
        c1 = sxy / sxx;

        if (isfinite(c1))
            dstp[i] *= (float)c1;
    }
}

void processColumnWSLRMasked(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p)
{
    const int ref_line_size = p->ref_line_size;
    const double sigmaR2 = p->sigmaR * p->sigmaR;
    const double sigmaD2 = p->sigmaD * p->sigmaD;
    const double* const spatial_exponent = p->scratch->spatial_exponent + p->scratch->spatial_radius;
    const float* wmaskp = p->wmaskp;
    const ptrdiff_t wmaskstride = p->wmaskstride;
    const int mask_dist = p->mask_dist;

    int sign = 1;
    if (column >= w / 2)
        sign = -1;

    int i;
    int j;
    int k;
    int start, stop;
    double c1;

    double* const cur = p->scratch->cur;
    double* const ref = p->scratch->ref;
    double* const ratio = p->scratch->ratio;

    wmaskp += column;
    dstp += column;
    for (i = 0; i < h; i++)
    {
        cur[i] = dstp[i * stride];
        ref[i] = dstp[sign + stride * i];
    }
    computeRatios(cur, ref, ratio, h);

    sign *= mask_dist;

    for (i = 0; i < h; i++)
    {
        start = i - ref_line_size;
        stop = i + ref_line_size + 1;
        if (start < 0)
            start = 0;
        if (stop > h)
            stop = h;

        double sxy = 0.0, sxx = 0.0;
        for (k = 0; k < stop - start; k++)
        {
            j = k + start;

            // exp(a) * exp(b) * exp(c) evaluated as exp(a + b + c)
            const double wt = exp(spatial_exponent[j - i] - ((ref[j] - ref[i]) * (ref[j] - ref[i])) / sigmaR2 - (ratio[j] - ratio[i]) * (ratio[j] - ratio[i]) / sigmaD2) * wmaskp[j * wmaskstride] * wmaskp[sign + wmaskstride * j];
            sxy += wt * cur[j] * ref[j];
            sxx += wt * cur[j] * cur[j];
        }

        // weighted least squares through the origin (what gsl_fit_wmul computes)
        c1 = sxy / sxx;

        if (isfinite(c1))
            dstp[i * stride] *= (float)c1;
    }
}
