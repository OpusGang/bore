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

    // Per-frame scratch buffers, sized for the longest line of the plane.
    typedef struct BoreScratch BoreScratch;

    BoreScratch* boreScratchAlloc(int w, int h, int multi); // multi: MultiPlane buffers
    void boreScratchFree(BoreScratch* s);

    // Per-line parameters passed to the line processing functions.
    typedef struct
    {
        int ref_line_size;
        double sigmaS;
        double sigmaR;
        double sigmaD;
        const float* wmaskp; // weight mask plane, NULL if unused
        ptrdiff_t wmaskstride;
        int mask_dist; // distance from the current line to the first unprocessed (clean) line
        double* c1_cov11_sumsq; // debug functions: fit results of the line
        BoreScratch* scratch;
    } BoreLineParams;

    typedef void (*BoreLineFunc)(int line, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p);

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
        BoreLineFunc processRow;
        BoreLineFunc processColumn;
    } LinearRegressionData;

    void boreInit(void);

    void processRowSLR(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p);
    void processColumnSLR(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p);
    void processRowSLRMasked(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p);
    void processColumnSLRMasked(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p);
    void debugRowSLR(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p);
    void debugColumnSLR(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p);
    void debugRowSLRMasked(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p);
    void debugColumnSLRMasked(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p);
    void processRowMLR(int row, int w, int h, ptrdiff_t stride, float* dstp, float* dstp1, float* dstp2, float* dstp3, BoreScratch* scratch);
    void processColumnMLR(int column, int w, int h, ptrdiff_t stride, float* dstp, float* dstp1, float* dstp2, float* dstp3, BoreScratch* scratch);
    void processRowSLRRef(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p);
    void processColumnSLRRef(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p);
    void processRowSLRRefMasked(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p);
    void processColumnSLRRefMasked(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p);
    void processRowWSLR(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p);
    void processColumnWSLR(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p);
    void processRowWSLRMasked(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p);
    void processColumnWSLRMasked(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, const BoreLineParams* p);

#ifdef __cplusplus
}
#endif

#endif
