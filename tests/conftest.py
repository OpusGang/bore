"""Shared helpers for the bore VapourSynth plugin tests.

The plugin under test is taken from the BORE_PLUGIN environment variable
(path to the built libbore.so / bore.dll).  Tests are skipped if it is unset
or if vapoursynth / numpy cannot be imported.
"""

import os

import pytest

PLUGIN_ENV = "BORE_PLUGIN"


def plugin_path():
    return os.environ.get(PLUGIN_ENV)


@pytest.fixture(scope="session")
def core():
    vs = pytest.importorskip("vapoursynth")
    path = plugin_path()
    if not path:
        pytest.skip(f"{PLUGIN_ENV} is not set")
    c = vs.core
    if not hasattr(c, "bore"):
        c.std.LoadPlugin(path)
    return c


# --------------------------------------------------------------------------
# Synthetic data
# --------------------------------------------------------------------------

def smooth_image(w, h, seed=0, offset=0.0):
    """Deterministic smooth, strictly positive float32 image with mild texture."""
    import numpy as np

    rng = np.random.default_rng(seed)
    y, x = np.mgrid[0:h, 0:w].astype(np.float64)
    img = (
        0.45
        + offset
        + 0.15 * np.sin(2 * np.pi * x / w * 1.3 + 0.4 + seed) * np.cos(2 * np.pi * y / h * 0.9)
        + 0.10 * x / w
        + 0.05 * y / h
        + 0.004 * rng.standard_normal((h, w))
    )
    return img.astype(np.float32)


def degrade(img, top=0, bottom=0, left=0, right=0, factors=(0.6, 0.75, 0.88, 0.95)):
    """Multiply the outermost lines of each side by factors[k] (k=0 is the edge)."""
    out = img.copy()
    h, w = out.shape
    for k in range(top):
        out[k, :] *= factors[k]
    for k in range(bottom):
        out[h - 1 - k, :] *= factors[k]
    for k in range(left):
        out[:, k] *= factors[k]
    for k in range(right):
        out[:, w - 1 - k] *= factors[k]
    return out


def make_clip(core, planes_per_frame, fmt):
    """Build a clip of format `fmt` from a list (frames) of lists (planes) of 2D arrays."""
    import numpy as np

    first = planes_per_frame[0]
    h, w = first[0].shape
    blank = core.std.BlankClip(format=fmt, width=w, height=h, length=len(planes_per_frame))
    data = [[np.ascontiguousarray(p, dtype=np.float32) for p in fr] for fr in planes_per_frame]

    def fill(n, f):
        fout = f.copy()
        for p, arr in enumerate(data[n]):
            np.asarray(fout[p])[...] = arr
        return fout

    return core.std.ModifyFrame(clip=blank, clips=blank, selector=fill)


def gray_clip(core, img, frames=1):
    import vapoursynth as vs

    return make_clip(core, [[img]] * frames, vs.GRAYS)


def get_planes(clip, n=0):
    import numpy as np

    f = clip.get_frame(n)
    return [np.array(f[p], copy=True) for p in range(f.format.num_planes)]
