#include <stdlib.h>
#include <vapoursynth/VapourSynth4.h>
#include <vapoursynth/VSHelper4.h>

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

VS_EXTERNAL_API(void) VapourSynthPluginInit2(VSPlugin *plugin, const VSPLUGINAPI *vspapi) {
    vspapi->configPlugin("ng.opusga.bore", "bore", "bore plugin", VS_MAKE_VERSION(1, 0), VAPOURSYNTH_API_VERSION, 0, plugin);
    vspapi->registerFunction("FixBrightness", "clip:vnode;top:int:opt;bottom:int:opt;left:int:opt;right:int:opt;lower:float:opt;upper:float:opt;thrlo:float:opt;thrhi:float:opt;step:int:opt;plane:int:opt;", "clip:vnode;", fixBrightnessCreate, NULL, plugin);
}

