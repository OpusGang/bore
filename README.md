# FixBrightness

This is a bbmod-style border deringer, which works by multiplying every pixel by the average of dividing each pixel in a line by its nearest pixel in the next line. Unlike bbmod, this doesn't attempt to localize the result and just assumes the same adjustment to be made for the entire row. It's effectively an automated FixBrightness/rektlvls, hence the name.

## Usage

```
core.bore.FixBrightness(clip clip, int top=0, int bottom=0, int left=0, int right=0, clip ignore_mask=None, float thrlo=0.1, float thrhi=8.0, int step=1, int plane=0)
```

* `clip`: 32-bit float clip.
* `top = 0`, `bottom = 0`, `left = 0`, `right = 0`: number of lines from each border to adjust.
* `ignore_mask = None`: Ignore mask, needs to be 8-bit. Anything below 128 will be ignored during adjustment calculation.
* `thrlo = 0.1`: Lower limit of adjustment. Any quotient below this will be ignored.
* `thrhi = 8.0`: Upper limit of adjustment. Any quotient above this will be ignored.
* `step = 1`: Speed up processing ever so slightly by lowering the number of pixels used to find the adjustment. Might be an unnecessary parameter.
* `plane = 0`: Plane to process.

# Balance

This approach to border deringing uses [linear least squares](https://www.gnu.org/software/gsl/doc/html/lls.html) to find a proper adjustment. 

In simple (0) `mode`, it functions similarly to `FixBrightness`, using simple linear regression instead of a simple mean. This is more robust than `FixBrightness`, although a bit slower.
In multiple (1) `mode`, it uses multiple linear regression with all three planes as input parameters, which may help with ringing that was added in a different color space. This should be considered a last resort effort. 

## Requirements
* [GSL](https://www.gnu.org/software/gsl/)

## Usage

```
core.bore.Balance(clip clip, int top=0, int bottom=0, int left=0, int right=0, int plane=0, int mode=0)
```

* `clip`: 32-bit float clip.
* `top = 0`, `bottom = 0`, `left = 0`, `right = 0`: number of lines from each border to adjust.
* `plane = 0`: Plane to adjust.
* `mode = 0`: Linear regression mode, 0 for simple and 1 for multiple linear regression using all three planes as parameters (requires no subsampling).

# Compilation

```sh 
meson build
ninja -C build 
```

In Windows you can configure MinGW and run the meson build, but the easiest way is to [download LLVM](https://github.com/llvm/llvm-project/releases) (*-win64.exe file) and install with the **path option ticked**, then run:

```ps
clang -O2 -march=native -shared -o bore.dll -I"C:\Program Files\VapourSynth\sdk\include" src\bore.c
```
