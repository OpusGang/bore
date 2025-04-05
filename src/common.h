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

#ifndef BORE_COMMON_H
#define BORE_COMMON_H

#include <stdarg.h>
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

    typedef enum
    {
        LINREG_MODE_SINGLE = 1,
        LINREG_MODE_MULTI = 2,
        LINREG_MODE_SINGLE_LIMITED = 3,
        LINREG_MODE_SINGLE_WEIGHTED = 4,
        LINREG_MODE_SINGLE_DEBUG = 5
    } LinRegMode;

    typedef struct
    {
        int plane;
        int top;
        int bottom;
        int left;
        int right;
        int ref_line_size;
        double sigmaS;
        double sigmaR;
        double sigmaD;
        void (*processRow)(int, int w, int h, ptrdiff_t stride, float* __restrict dstp, ...);
        void (*processColumn)(int, int w, int h, ptrdiff_t stride, float* __restrict dstp, ...);
    } LinearRegressionData;

    void processRowSLR(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, ...);
    void processColumnSLR(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, ...);
    void processRowSLRMasked(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, ...);
    void processColumnSLRMasked(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, ...);
    void debugRowSLR(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, double* c1_cov11_sumsq);
    void debugColumnSLR(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, double* c1_cov11_sumsq);
    void debugRowSLRMasked(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, double* c1_cov11_sumsq,
        const float* __restrict wmaskp, ptrdiff_t wmaskstride, int mask_dist);
    void debugColumnSLRMasked(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, double* c1_cov11_sumsq,
        const float* __restrict wmaskp, ptrdiff_t wmaskstride, int mask_dist);
    void processRowMLR(int row, int w, int h, ptrdiff_t stride, float* dstp, float* dstp1, float* dstp2, float* dstp3);
    void processColumnMLR(int column, int w, int h, ptrdiff_t stride, float* dstp, float* dstp1, float* dstp2, float* dstp3);
    void processRowSLRRef(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, ...);
    void processColumnSLRRef(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, ...);
    void processRowSLRRefMasked(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, ...);
    void processColumnSLRRefMasked(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, ...);
    void processRowWSLR(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, ...);
    void processColumnWSLR(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, ...);
    void processRowWSLRMasked(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, ...);
    void processColumnWSLRMasked(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, ...);

#ifdef __cplusplus
}
#endif

#endif
