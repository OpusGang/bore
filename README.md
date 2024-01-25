# FixBrightness

This is a bbmod-style border deringer, which works by multiplying every pixel by the average of dividing each pixel in a line by its nearest pixel in the next line. Unlike bbmod, this doesn't attempt to localize the result and just assumes the same adjustment to be made for the entire row. It's effectively an automated FixBrightness/rektlvls, hence the name.

## Usage

```
core.bore.FixBrightness(clip clip, int top=0, int bottom=0, int left=0, int right=0, int step=1, float upper=1.0, float lower=0.0, int plane=0)
```

* `clip`: 32-bit float clip.
* `top = 0`, `bottom = 0`, `left = 0`, `right = 0`: number of lines from each border to adjust. Negative bottom/right values just count from the other side.
* `step = 1`: Speed up processing ever so slightly by lowering the number of pixels used to find the adjustment. Might be an unnecessary parameter.
* `upper = 1`: Upper limit of range, this allows excluding pixels.
* `lower = 0`: Lower limit of range, this allows exculding pixels.
* `plane = 0`: Plane to process.

## Compilation

```sh 
meson build
ninja -C build 
```

