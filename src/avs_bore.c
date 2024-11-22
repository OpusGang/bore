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

#include "avisynth/avisynth_c.h"
#include "common.h"

typedef struct
{
    union
    {
        LinearRegressionData data;
    } shared;
    AVS_Clip* weight_mask;
} AVS_LinearRegressionData;

static AVS_FORCEINLINE int avs_plane(const int plane_in, const AVS_VideoInfo* vi)
{
    const int planes_rgb[3] = { AVS_PLANAR_R, AVS_PLANAR_G, AVS_PLANAR_B };
    const int planes_yuv[3] = { AVS_PLANAR_Y, AVS_PLANAR_U, AVS_PLANAR_V };
    const int* planes = avs_is_rgb(vi) ? planes_rgb : planes_yuv;

    return planes[plane_in];
}

static AVS_VideoFrame* AVSC_CC singlePlaneGetFrame(AVS_FilterInfo* fi, int n)
{
    AVS_LinearRegressionData* d = (AVS_LinearRegressionData*)fi->user_data;

    AVS_VideoFrame* frame = avs_get_frame(fi->child, n);
    if (!frame)
        return NULL;

    avs_make_writable(fi->env, &frame);

    const int plane = avs_plane(d->shared.data.plane, &fi->vi);
    const int top = d->shared.data.top;
    const int bottom = d->shared.data.bottom;
    const int left = d->shared.data.left;
    const int right = d->shared.data.right;
    const int ref_line_size = d->shared.data.ref_line_size;
    const double sigmaS = d->shared.data.sigmaS;
    const double sigmaR = d->shared.data.sigmaR;
    const double sigmaD = d->shared.data.sigmaD;
    ptrdiff_t stride = avs_get_pitch_p(frame, plane) / 4;
    int w = avs_get_row_size_p(frame, plane) / 4;
    int h = avs_get_height_p(frame, plane);
    float* __restrict dstp = (float*)avs_get_write_ptr_p(frame, plane);

    if (d->weight_mask)
    {
        AVS_VideoFrame* weight_mask = avs_get_frame(d->weight_mask, n);
        ptrdiff_t wmaskstride = avs_get_pitch(weight_mask) / 4;
        const float* wmaskp = (float*) avs_get_read_ptr(weight_mask);
        if (top != 0)
        {
            for (int row = top - 1; row > -1; --row)
                d->shared.data.processRow(row, w, h, stride, dstp, ref_line_size, sigmaS, sigmaR, sigmaD, wmaskp, wmaskstride, top - row);
        }
        if (bottom != 0)
        {
            for (int row = h - bottom; row < h; ++row)
                d->shared.data.processRow(row, w, h, stride, dstp, ref_line_size, sigmaS, sigmaR, sigmaD, wmaskp, wmaskstride, bottom + row - h + 1);
        }
        if (left != 0)
        {
            for (int column = left - 1; column > -1; --column)
                d->shared.data.processColumn(column, w, h, stride, dstp, ref_line_size, sigmaS, sigmaR, sigmaD, wmaskp, wmaskstride, left - column);
        }
        if (right != 0)
        {
            for (int column = w - right; column < w; ++column)
                d->shared.data.processColumn(column, w, h, stride, dstp, ref_line_size, sigmaS, sigmaR, sigmaD, wmaskp, wmaskstride, right + column - w + 1);
        }

        avs_release_video_frame(weight_mask);
    }
    else
    {
        if (top != 0)
        {
            for (int row = top - 1; row > -1; --row)
                d->shared.data.processRow(row, w, h, stride, dstp, ref_line_size, sigmaS, sigmaR, sigmaD);
        }
        if (bottom != 0)
        {
            for (int row = h - bottom; row < h; ++row)
                d->shared.data.processRow(row, w, h, stride, dstp, ref_line_size, sigmaS, sigmaR, sigmaD);
        }
        if (left != 0)
        {
            for (int column = left - 1; column > -1; --column)
                d->shared.data.processColumn(column, w, h, stride, dstp, ref_line_size, sigmaS, sigmaR, sigmaD);
        }
        if (right != 0)
        {
            for (int column = w - right; column < w; ++column)
                d->shared.data.processColumn(column, w, h, stride, dstp, ref_line_size, sigmaS, sigmaR, sigmaD);
        }
    }

    return frame;
}

static AVS_VideoFrame* AVSC_CC multiPlaneGetFrame(AVS_FilterInfo* fi, int n)
{
    AVS_LinearRegressionData* d = (AVS_LinearRegressionData*)fi->user_data;

    AVS_VideoFrame* frame = avs_get_frame(fi->child, n);
    if (!frame)
        return NULL;

    avs_make_writable(fi->env, &frame);

    const int plane = avs_plane(d->shared.data.plane, &fi->vi);
    const int top = d->shared.data.top;
    const int bottom = d->shared.data.bottom;
    const int left = d->shared.data.left;
    const int right = d->shared.data.right;
    ptrdiff_t stride = avs_get_pitch_p(frame, plane) / 4;
    int w = avs_get_row_size_p(frame, plane) / 4;
    int h = avs_get_height_p(frame, plane);
    float* dstp = (float*)avs_get_write_ptr_p(frame, plane);
    float* dstp1 = (float*)avs_get_write_ptr_p(frame, avs_plane(0, &fi->vi));
    float* dstp2 = (float*)avs_get_write_ptr_p(frame, avs_plane(1, &fi->vi));
    float* dstp3 = (float*)avs_get_write_ptr_p(frame, avs_plane(2, &fi->vi));

    if (top != 0)
    {
        for (int row = top - 1; row > -1; --row)
            processRowMLR(row, w, h, stride, dstp, dstp1, dstp2, dstp3);
    }
    if (bottom != 0)
    {
        for (int row = h - bottom; row < h; ++row)
            processRowMLR(row, w, h, stride, dstp, dstp1, dstp2, dstp3);
    }
    if (left != 0)
    {
        for (int column = left - 1; column > -1; --column)
            processColumnMLR(column, w, h, stride, dstp, dstp1, dstp2, dstp3);
    }
    if (right != 0)
    {
        for (int column = w - right; column < w; ++column)
            processColumnMLR(column, w, h, stride, dstp, dstp1, dstp2, dstp3);
    }

    return frame;
}

static AVS_FORCEINLINE void set_frame_props(AVS_VideoFrame* frame, double* props, AVS_ScriptEnvironment* env)
{
    AVS_Map* __restrict dst_props = avs_get_frame_props_rw(env, frame);
    avs_prop_set_float(env, dst_props, "BoreAdjustment", props[0], AVS_PROPAPPENDMODE_APPEND);
    avs_prop_set_float(env, dst_props, "BoreCovariance", props[1], AVS_PROPAPPENDMODE_APPEND);
    avs_prop_set_float(env, dst_props, "BoreSumSquares", props[2], AVS_PROPAPPENDMODE_APPEND);
}

static AVS_VideoFrame* AVSC_CC singlePlaneDebugGetFrame(AVS_FilterInfo* fi, int n)
{
    AVS_LinearRegressionData* d = (AVS_LinearRegressionData*)fi->user_data;

    AVS_VideoFrame* frame = avs_get_frame(fi->child, n);
    if (!frame)
        return NULL;

    avs_make_writable(fi->env, &frame);

    const int plane = avs_plane(d->shared.data.plane, &fi->vi);
    const int top = d->shared.data.top;
    const int bottom = d->shared.data.bottom;
    const int left = d->shared.data.left;
    const int right = d->shared.data.right;
    ptrdiff_t stride = avs_get_pitch_p(frame, plane) / 4;
    int w = avs_get_row_size_p(frame, plane) / 4;
    int h = avs_get_height_p(frame, plane);
    float* __restrict dstp = (float*)avs_get_write_ptr_p(frame, plane);
    double c1_cov11_sumsq[3] = { 0.0 };
    double* props = c1_cov11_sumsq;

    if (d->weight_mask)
    {
        AVS_VideoFrame* weight_mask = avs_get_frame(d->weight_mask, n);
        ptrdiff_t wmaskstride = avs_get_pitch(weight_mask) / 4;
        const float* wmaskp = (float*) avs_get_read_ptr(weight_mask);
        if (top != 0)
        {
            for (int row = top - 1; row > -1; --row)
            {
                debugRowSLRMasked(row, w, h, stride, dstp, &props, wmaskp, wmaskstride, top - row);
                set_frame_props(frame, props, fi->env);
            }
        }
        if (bottom != 0)
        {
            for (int row = h - bottom; row < h; ++row)
            {
                debugRowSLRMasked(row, w, h, stride, dstp, &props, wmaskp, wmaskstride, bottom + row - h + 1);
                set_frame_props(frame, props, fi->env);
            }
        }
        if (left != 0)
        {
            for (int column = left - 1; column > -1; --column)
            {
                debugColumnSLRMasked(column, w, h, stride, dstp, &props, wmaskp, wmaskstride, left - column);
                set_frame_props(frame, props, fi->env);
            }
        }
        if (right != 0)
        {
            for (int column = w - right; column < w; ++column)
            {
                debugColumnSLRMasked(column, w, h, stride, dstp, &props, wmaskp, wmaskstride, right + column - w + 1);
                set_frame_props(frame, props, fi->env);
            }
        }

        avs_release_video_frame(weight_mask);
    }
    else
    {
        if (top != 0)
        {
            for (int row = top - 1; row > -1; --row)
            {
                debugRowSLR(row, w, h, stride, dstp, &props);
                set_frame_props(frame, props, fi->env);
            }
        }
        if (bottom != 0)
        {
            for (int row = h - bottom; row < h; ++row)
            {
                debugRowSLR(row, w, h, stride, dstp, &props);
                set_frame_props(frame, props, fi->env);
            }
        }
        if (left != 0)
        {
            for (int column = left - 1; column > -1; --column)
            {
                debugColumnSLR(column, w, h, stride, dstp, &props);
                set_frame_props(frame, props, fi->env);
            }
        }
        if (right != 0)
        {
            for (int column = w - right; column < w; ++column)
            {
                debugColumnSLR(column, w, h, stride, dstp, &props);
                set_frame_props(frame, props, fi->env);
            }
        }
    }

    return frame;
}

static void AVSC_CC free_bore(AVS_FilterInfo* fi)
{
    AVS_LinearRegressionData* d = (AVS_LinearRegressionData*)fi->user_data;

    if (d->weight_mask)
        avs_release_clip(d->weight_mask);

    free(d);
}


static int AVSC_CC set_cache_hints_bore(AVS_FilterInfo* fi, int cachehints, int frame_range)
{
    return cachehints == AVS_CACHE_GET_MTMODE ? 1 : 0;
}

static AVS_FORCEINLINE AVS_Value set_error(AVS_Clip* clip, AVS_Clip* weight_mask, const char* msg)
{
    if (weight_mask)
        avs_release_clip(weight_mask);

    avs_release_clip(clip);

    return avs_new_value_error(msg);
}

static AVS_Value AVSC_CC linearRegressionCreate(AVS_ScriptEnvironment* env, AVS_Value args, void* param)
{
    typedef enum
    {
        Clip,
        Left,
        Right,
        Top,
        Bottom,
        Weight_mask,
        Plane
    } args_name;

    const int mode = (LinRegMode)param;

    AVS_LinearRegressionData* d = (AVS_LinearRegressionData*)malloc(sizeof(AVS_LinearRegressionData));

    AVS_FilterInfo* fi;
    AVS_Clip* clip = avs_new_c_filter(env, &fi, avs_array_elt(args, Clip), 1);

    if (avs_check_version(env, 10))
        return set_error(clip, NULL, "bore: AviSynth+ version must be r3928 or later.");

    const int plane = avs_defined(avs_array_elt(args, Plane)) ? avs_as_int(avs_array_elt(args, Plane)) : 0;

    const int num_planes = avs_num_components(&fi->vi);

    if (!avs_is_planar(&fi->vi) || avs_component_size(&fi->vi) < 4)
        return set_error(clip, NULL, "bore: clip must be in 32-bit planar format.");
    if (num_planes < 3 && plane > 0)
        return set_error(clip, NULL, "bore: plane cannot be bigger than 0.");
    if (num_planes == 4)
        return set_error(clip, NULL, "bore: clip must have less than 4 planes.");

    AVS_Clip* weight_mask = avs_defined(avs_array_elt(args, Weight_mask)) ? (avs_take_clip(avs_array_elt(args, Weight_mask), env)) : NULL;
    if (weight_mask)
    {
        const AVS_VideoInfo* ivi = avs_get_video_info(weight_mask);

        if (!avs_is_planar(ivi) || avs_num_components(ivi) > 1)
            return set_error(clip, weight_mask, "bore: weight_mask must be in gray planar format.");
        if (avs_component_size(ivi) < 4)
            return set_error(clip, weight_mask, "bore: weight_mask bit depth must be 32-bit.");
        if (ivi->width != fi->vi.width || ivi->height != fi->vi.height)
            return set_error(clip, weight_mask, "bore: clip and weight_mask must have matching dimensions.");
        if (ivi->num_frames != fi->vi.num_frames)
            return set_error(clip, weight_mask, "bore: clip and weight_mask must have same frames number.");
    }

    const int is_rgb = avs_is_rgb(&fi->vi);
    const int half_h = ((plane == 0 || is_rgb) ? fi->vi.height
        : (fi->vi.height >> avs_get_plane_height_subsampling(&fi->vi, AVS_PLANAR_U))) >> 1;
    const int half_w = ((plane == 0 || is_rgb) ? fi->vi.width
        : (fi->vi.width >> avs_get_plane_width_subsampling(&fi->vi, AVS_PLANAR_U))) >> 1;

    const int top = avs_defined(avs_array_elt(args, Top)) ? avs_as_int(avs_array_elt(args, Top)) : 0;
    if (top > half_h)
        return set_error(clip, weight_mask, "bore: top must be in [0, height / 2].");

    const int bottom = avs_defined(avs_array_elt(args, Bottom)) ? avs_as_int(avs_array_elt(args, Bottom)) : 0;
    if (bottom > half_h)
        return set_error(clip, weight_mask, "bore: bottom must be in [0, height / 2].");

    const int left = avs_defined(avs_array_elt(args, Left)) ? avs_as_int(avs_array_elt(args, Left)) : 0;
    if (left > half_w)
        return set_error(clip, weight_mask, "bore: left must be in [0, width / 2].");

    const int right = avs_defined(avs_array_elt(args, Right)) ? avs_as_int(avs_array_elt(args, Right)) : 0;
    if (right > half_w)
        return set_error(clip, weight_mask, "bore: right must be in [0, width / 2].");

    d->weight_mask = weight_mask;
    d->shared.data.plane = plane;
    d->shared.data.top = top;
    d->shared.data.bottom = bottom;
    d->shared.data.left = left;
    d->shared.data.right = right;

    if (mode == LINREG_MODE_SINGLE_LIMITED || mode == LINREG_MODE_SINGLE_WEIGHTED)
    {
        typedef enum
        {
            Ref_line_size = 7
        } ref_line_arg_name;

        d->shared.data.ref_line_size = avs_defined(avs_array_elt(args, Ref_line_size)) ?
            avs_as_int(avs_array_elt(args, Ref_line_size)) : 100;
    }
    else
        d->shared.data.ref_line_size = 0;

    switch (mode)
    {
        case LINREG_MODE_SINGLE:
            if (weight_mask)
            {
                d->shared.data.processRow = &processRowSLRMasked;
                d->shared.data.processColumn = &processColumnSLRMasked;
            }
            else
            {
                d->shared.data.processRow = &processRowSLR;
                d->shared.data.processColumn = &processColumnSLR;
            }
            fi->get_frame = singlePlaneGetFrame;
            break;
        case LINREG_MODE_MULTI:
            if (num_planes == 1)
                return set_error(clip, weight_mask, "bore: clip must have 3 planes.");
            if (!is_rgb && !avs_is_444(&fi->vi))
                return set_error(clip, weight_mask, "bore: only 444 format is supported.");

            fi->get_frame = multiPlaneGetFrame;
            break;
        case LINREG_MODE_SINGLE_LIMITED:
            if (weight_mask)
            {
                d->shared.data.processRow = &processRowSLRRefMasked;
                d->shared.data.processColumn = &processColumnSLRRefMasked;
            }
            else
            {
                d->shared.data.processRow = &processRowSLRRef;
                d->shared.data.processColumn = &processColumnSLRRef;
            }
            fi->get_frame = singlePlaneGetFrame;
            break;
        case LINREG_MODE_SINGLE_WEIGHTED:
        {
            typedef enum
            {
                SigmaS = 8,
                SigmaR,
                SigmaD
            } sigma_arg_name;

            d->shared.data.sigmaS = avs_defined(avs_array_elt(args, SigmaS)) ? avs_as_float(avs_array_elt(args, SigmaS)) : 50.0;
            d->shared.data.sigmaR = avs_defined(avs_array_elt(args, SigmaR)) ? avs_as_float(avs_array_elt(args, SigmaR)) : 0.5;
            d->shared.data.sigmaD = avs_defined(avs_array_elt(args, SigmaD)) ? avs_as_float(avs_array_elt(args, SigmaD)) : 1.5;

            if (weight_mask)
            {
                d->shared.data.processRow = &processRowWSLRMasked;
                d->shared.data.processColumn = &processColumnWSLRMasked;
            }
            else
            {
                d->shared.data.processRow = &processRowWSLR;
                d->shared.data.processColumn = &processColumnWSLR;
            }
            fi->get_frame = singlePlaneGetFrame;
            break;
        }
        case LINREG_MODE_SINGLE_DEBUG:
            fi->get_frame = singlePlaneDebugGetFrame;
            break;
    }

    AVS_Value v = avs_new_value_clip(clip);

    fi->user_data = (void*)d;
    fi->set_cache_hints = set_cache_hints_bore;
    fi->free_filter = free_bore;

    avs_release_clip(clip);

    return v;
}

const char* AVSC_CC avisynth_c_plugin_init(AVS_ScriptEnvironment* env)
{
    avs_add_function(env, "bore_SinglePlane",
        "c"
        "[left]i"
        "[right]i"
        "[top]i"
        "[bottom]i"
        "[weight_mask]c"
        "[plane]i", linearRegressionCreate, (void*)LINREG_MODE_SINGLE);
    avs_add_function(env, "bore_MultiPlane",
        "c"
        "[left]i"
        "[right]i"
        "[top]i"
        "[bottom]i"
        "[weight_mask]c"
        "[plane]i", linearRegressionCreate, (void*)LINREG_MODE_MULTI);
    avs_add_function(env, "bore_SinglePlaneLimited",
        "c"
        "[left]i"
        "[right]i"
        "[top]i"
        "[bottom]i"
        "[weight_mask]c"
        "[plane]i"
        "[ref_line_size]i", linearRegressionCreate, (void*)LINREG_MODE_SINGLE_LIMITED);
    avs_add_function(env, "bore_SinglePlaneWeighted",
        "c"
        "[left]i"
        "[right]i"
        "[top]i"
        "[bottom]i"
        "[weight_mask]c"
        "[plane]i"
        "[ref_line_size]i"
        "[sigmaS]f"
        "[sigmaR]f"
        "[sigmaD]f", linearRegressionCreate, (void*)LINREG_MODE_SINGLE_WEIGHTED);
    avs_add_function(env, "bore_SinglePlaneDebug",
        "c"
        "[left]i"
        "[right]i"
        "[top]i"
        "[bottom]i"
        "[weight_mask]c"
        "[plane]i", linearRegressionCreate, (void*)LINREG_MODE_SINGLE_DEBUG);

    return "bore";
}
