# FixBrightness

This is a bbmod-style border deringer, which works by multiplying every pixel by the average of dividing each pixel in a line by its nearest pixel in the next line. Unlike bbmod, this doesn't attempt to localize the result and just assumes the same adjustment to be made for the entire row. It's effectively an automated FixBrightness/rektlvls, hence the name.

## Usage

```
core.bore.FixBrightness(clip clip, int top=0, int bottom=0, int left=0, int right=0, int step=1, float upper=1.0, float lower=0.0, int plane=0)
```

* `clip`: 32-bit float clip.
* `top = 0`, `bottom = 0`, `left = 0`, `right = 0`: number of lines from each border to adjust.
* `step = 1`: Speed up processing ever so slightly by lowering the number of pixels used to find the adjustment. Might be an unnecessary parameter.
* `upper = 1`: Upper limit of range, this allows excluding pixels.
* `lower = 0`: Lower limit of range, this allows excluding pixels.
* `plane = 0`: Plane to process.

# Balance

This approach to border deringing uses [linear least squares](https://www.gnu.org/software/gsl/doc/html/lls.html) to find a proper adjustment. 

In simple (0) `mode`, it functions similarly to `FixBrightness`, using simple linear regression instead of a simple mean. This is more robust than `FixBrightness`, although a bit slower.
In multiple (1) `mode`, it uses multiple linear regression with all three planes as input parameters, which may help with ringing that was added in a different color space. This should be considered a last resort effort. 

## Requirements
* [GSL](https://www.gnu.org/software/gsl/)

## Usage

```
core.bore.Balance(clip clip, int top=0, int bottom=0, int left=0, int right=0, float upper=1.0, float lower=0.0, int plane=0)
```

* `clip`: 32-bit float clip.
* `top = 0`, `bottom = 0`, `left = 0`, `right = 0`: number of lines from each border to adjust.
* `plane = 0`: Plane to adjust.
* `mode = 0`: Linear regression mode.

# Compilation

```sh 
meson build
ninja -C build 
```

In Windows you can configure MinGW and run the meson build, but the easiest way is to [download LLVM](https://github.com/llvm/llvm-project/releases) (*-win64.exe file) and install with the **path option ticked**, then run:

```ps
clang -O2 -march=native -shared -o bore.dll -I"C:\Program Files\VapourSynth\sdk\include" src\bore.c
```
