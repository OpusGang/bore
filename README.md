# bore

bore fixes dirty borders: lines near the frame edge that are brighter or darker than they should be. Each border line is scaled to match its inner neighbour, with the scale factor found by [linear least squares](https://www.gnu.org/software/gsl/doc/html/lls.html). Lines are processed from the innermost to the outermost, so each one is compared with a line that has already been fixed. The plugin works with VapourSynth and AviSynth+.

## Filters

* **SinglePlane**: simple linear regression between each border line and its inner neighbour; one factor per line.
* **MultiPlane**: the same with multiple linear regression, using all three planes of the line as predictors for the target plane.
* **SinglePlaneLimited**: a separate factor for each pixel, fitted only on the nearest `2 * ref_line_size + 1` pixels of the line.
* **SinglePlaneWeighted**: SinglePlaneLimited with bilateral-filter-style weights for distance and pixel difference.
* **SinglePlaneDebug**: computes the same factors as SinglePlane but leaves the pixels unchanged and stores the values in frame properties instead.

SinglePlaneLimited and SinglePlaneWeighted are more experimental and should be used with care.

SinglePlaneDebug appends one value per processed line to these frame properties:
* `BoreAdjustment`
* `BoreCovariance`
* `BoreSumSquares`

The values are in order of operation: top, bottom, left, right, each from the innermost to the outermost line.

## Parameters

| Name | Type | Default | Valid range | Used by | Meaning |
|---|---|---|---|---|---|
| `clip` | clip | required | 32-bit float, planar, constant format | all | Input clip. |
| `left`, `right` | int | 0 | 0 to plane width / 2 | all | Number of columns to fix at the left/right border. |
| `top`, `bottom` | int | 0 | 0 to plane height / 2 | all | Number of rows to fix at the top/bottom border. |
| `weight_mask` | clip | none | 32-bit float, one plane | all except MultiPlane | Per-pixel weights for the fit, useful for excluding edges and clipped areas. A pixel's weight is its mask value times the mask value at the same position on the first clean line. |
| `plane` | int | 0 | 0 to number of planes - 1 | all | Plane to fix. |
| `ref_line_size` | int | 100 | >= 0 | Limited, Weighted | Half-width of the window around each pixel; the window is `2 * ref_line_size + 1` pixels. |
| `sigmaS` | float | 50.0 | > 0 | Weighted | Spatial sigma: weight by distance from the current pixel. |
| `sigmaR` | float | 0.5 | > 0 | Weighted | Range sigma: weight by the difference between pixels on the same line. |
| `sigmaD` | float | 1.5 | > 0 | Weighted | Weight by the difference between a neighbour's adjustment and the current pixel's adjustment. |

Notes:
* Only 32-bit float input is accepted (no integer or 16-bit half float). Convert first, e.g. with `core.resize.Point(clip, format=vs.YUV444PS)`.
* MultiPlane needs a three-plane clip in RGB or YUV 4:4:4. It accepts `weight_mask` but ignores it.
* `weight_mask` dimensions: in VapourSynth it must match the processed plane (so a chroma plane of a subsampled clip needs a chroma-sized mask); in AviSynth+ it must match the clip's luma size and have the same number of frames.
* AviSynth+ clips with an alpha plane are rejected.
* The positional order of the optional arguments differs between VapourSynth and AviSynth+ (see the signatures below). Use keyword arguments.

## VapourSynth usage

```
core.bore.SinglePlane(clip clip, int left=0, int right=0, int top=0, int bottom=0, clip weight_mask=None, int plane=0)
core.bore.MultiPlane(clip clip, int left=0, int right=0, int top=0, int bottom=0, clip weight_mask=None, int plane=0)
core.bore.SinglePlaneLimited(clip clip, int left=0, int right=0, int top=0, int bottom=0, clip weight_mask=None, int ref_line_size=100, int plane=0)
core.bore.SinglePlaneWeighted(clip clip, int left=0, int right=0, int top=0, int bottom=0, clip weight_mask=None, float sigmaS=50.0, float sigmaR=0.5, float sigmaD=1.5, int ref_line_size=100, int plane=0)
core.bore.SinglePlaneDebug(clip clip, int left=0, int right=0, int top=0, int bottom=0, clip weight_mask=None, int plane=0)
```

Example:

```python
import vapoursynth as vs
core = vs.core

clip = core.resize.Point(src, format=vs.YUV420PS)
clip = core.bore.SinglePlane(clip, top=2, bottom=1, plane=0)
clip = core.bore.SinglePlaneLimited(clip, left=1, ref_line_size=50, plane=1)
```

## AviSynth+ usage

```
bore_SinglePlane(clip input, int "left", int "right", int "top", int "bottom", clip "weight_mask", int "plane")
bore_MultiPlane(clip input, int "left", int "right", int "top", int "bottom", clip "weight_mask", int "plane")
bore_SinglePlaneLimited(clip input, int "left", int "right", int "top", int "bottom", clip "weight_mask", int "plane", int "ref_line_size")
bore_SinglePlaneWeighted(clip input, int "left", int "right", int "top", int "bottom", clip "weight_mask", int "plane", int "ref_line_size", float "sigmaS", float "sigmaR", float "sigmaD")
bore_SinglePlaneDebug(clip input, int "left", int "right", int "top", int "bottom", clip "weight_mask", int "plane")
```

Example:

```
ConvertBits(32)
bore_SinglePlane(top=2, bottom=1)
bore_SinglePlaneWeighted(left=1, sigmaS=30.0, plane=0)
```

## Installation

The VapourSynth plugin is available as a wheel:

```sh
pip install vapoursynth-bore
```

## Building from source

Requirements:
* [GSL](https://www.gnu.org/software/gsl/)
* VapourSynth headers (API R55 or newer) and/or AviSynth+ headers
* meson and ninja

```sh
meson setup build -Dlibtype=vapoursynth
ninja -C build
```

Options:
* `-Dlibtype=vapoursynth|avisynth|both`: which frame server to build for (default `vapoursynth`). The AviSynth+ build uses the `avs_c_api_loader` submodule (`git submodule update --init`).
* `-Dstatic_gsl=true`: link GSL statically.

If the headers are not found through pkg-config, add their location to `CFLAGS`/`CXXFLAGS`. `meson install -C build` copies the plugin into the frame server's plugin directory. For a Windows build from Linux, `cross-mingw-x86_64.txt` is a MinGW cross file (`meson setup build --cross-file cross-mingw-x86_64.txt`).

## Tests

The test suite uses the VapourSynth plugin and needs `vapoursynth`, `numpy` and `pytest` in the Python environment. Build the plugin, then point `BORE_PLUGIN` at it:

```sh
BORE_PLUGIN=build/libbore.so pytest tests
```

Without `BORE_PLUGIN` (or without VapourSynth) the tests are skipped.
