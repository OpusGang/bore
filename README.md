# FixBrightness

This is a bbmod-style border deringer, which works by multiplying every pixel by the average of dividing each pixel in a line by its nearest pixel in the next line. Unlike bbmod, this doesn't attempt to localize the result and just assumes the same adjustment to be made for the entire row. It's effectively an automated FixBrightness/rektlvls, hence the name.

## Usage

```
core.bore.FixBrightness(clip clip, int top=0, int bottom=0, int left=0, int right=0, float lower=0.0, float upper=1.0, float thrlo=0.1, float thrhi=8.0, int step=1, int plane=0)
```

* `clip`: 32-bit float clip.
* `top = 0`, `bottom = 0`, `left = 0`, `right = 0`: number of lines from each border to adjust.
* `lower = 0`: Lower limit of range, this allows excluding pixels.
* `upper = 1`: Upper limit of range, this allows excluding pixels.
* `thrlo = 0.1`: Lower limit of adjustment. Any quotient below this will be ignored.
* `thrhi = 8.0`: Upper limit of adjustment. Any quotient above this will be ignored.
* `step = 1`: Speed up processing ever so slightly by lowering the number of pixels used to find the adjustment. Might be an unnecessary parameter.
* `plane = 0`: Plane to process.

## Compilation

```sh 
meson build
ninja -C build 
```

In Windows you can configure MinGW and run the meson build, but the easiest way is to [download LLVM](https://github.com/llvm/llvm-project/releases) (*-win64.exe file) and install with the **path option ticked**, then run:

```ps
clang -O2 -march=native -shared -o bore.dll -I"C:\Program Files\VapourSynth\sdk\include" src\bore.c
```
