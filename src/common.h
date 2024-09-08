#ifndef BORE_COMMON_H
#define BORE_COMMON_H

#include <stdlib.h>

typedef enum LinRegMode
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
    void (*processRow)(int, int, int, ptrdiff_t, float* __restrict, int, double, double, double, const unsigned char* __restrict, ptrdiff_t, int);
    void (*processColumn)(int, int, int, ptrdiff_t, float* __restrict, int, double, double, double, const unsigned char* __restrict, ptrdiff_t, int);
} LinearRegressionData;

void processRowSLR(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, int ref_line_size, double sigmaS, double sigmaR, double sigmaD, const unsigned char* __restrict imaskp, ptrdiff_t imaskstride, int mask_dist);
void processColumnSLR(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, int ref_line_size, double sigmaS, double sigmaR, double sigmaD, const unsigned char* __restrict imaskp, ptrdiff_t imaskstride, int mask_dist);
void processRowSLRMasked(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, int ref_line_size, double sigmaS, double sigmaR, double sigmaD, const unsigned char* __restrict imaskp, ptrdiff_t imaskstride, int mask_dist);
void processColumnSLRMasked(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, int ref_line_size, double sigmaS, double sigmaR, double sigmaD, const unsigned char* __restrict imaskp, ptrdiff_t imaskstride, int mask_dist);
void debugRowSLR(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, double** props);
void debugColumnSLR(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, double** props);
void debugRowSLRMasked(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, double** props, const unsigned char* __restrict imaskp, ptrdiff_t imaskstride, int mask_dist);
void debugColumnSLRMasked(int column, int w, int h, ptrdiff_t stride, float* __restrict dstp, double** props, const unsigned char* __restrict imaskp, ptrdiff_t imaskstride, int mask_dist);
void processRowMLR(int row, int w, int h, ptrdiff_t stride, float* __restrict dstp, float* __restrict dstp1, float* __restrict dstp2, float* __restrict dstp3);
void processColumnMLR(int column, int w, int h, ptrdiff_t stride, float*__restrict dstp, float*__restrict dstp1, float*__restrict dstp2, float*__restrict dstp3);
void processRowSLRRef(int row, int w, int h, ptrdiff_t stride, float*__restrict dstp, int ref_line_size, double sigmaS, double sigmaR, double sigmaD, const unsigned char* __restrict imaskp, ptrdiff_t imaskstride, int mask_dist);
void processColumnSLRRef(int column, int w, int h, ptrdiff_t stride, float*__restrict dstp, int ref_line_size, double sigmaS, double sigmaR, double sigmaD, const unsigned char* __restrict imaskp, ptrdiff_t imaskstride, int mask_dist);
void processRowSLRRefMasked(int row, int w, int h, ptrdiff_t stride, float*__restrict dstp, int ref_line_size, double sigmaS, double sigmaR, double sigmaD, const unsigned char* __restrict imaskp, ptrdiff_t imaskstride, int mask_dist);
void processColumnSLRRefMasked(int column, int w, int h, ptrdiff_t stride, float*__restrict dstp, int ref_line_size, double sigmaS, double sigmaR, double sigmaD, const unsigned char* __restrict imaskp, ptrdiff_t imaskstride, int mask_dist);
void processRowWSLR(int row, int w, int h, ptrdiff_t stride, float*__restrict dstp, int ref_line_size, double sigmaS, double sigmaR, double sigmaD, const unsigned char* __restrict imaskp, ptrdiff_t imaskstride, int mask_dist);
void processColumnWSLR(int column, int w, int h, ptrdiff_t stride, float*__restrict dstp, int ref_line_size, double sigmaS, double sigmaR, double sigmaD, const unsigned char* __restrict imaskp, ptrdiff_t imaskstride, int mask_dist);
void processRowWSLRMasked(int row, int w, int h, ptrdiff_t stride, float*__restrict dstp, int ref_line_size, double sigmaS, double sigmaR, double sigmaD, const unsigned char* __restrict imaskp, ptrdiff_t imaskstride, int mask_dist);
void processColumnWSLRMasked(int column, int w, int h, ptrdiff_t stride, float*__restrict dstp, int ref_line_size, double sigmaS, double sigmaR, double sigmaD, const unsigned char* __restrict imaskp, ptrdiff_t imaskstride, int mask_dist);

#endif
