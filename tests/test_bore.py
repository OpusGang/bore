"""Functional tests for the bore VapourSynth plugin.

Run with:  BORE_PLUGIN=/path/to/libbore.so python -m pytest tests
"""

import os
import subprocess
import sys
import textwrap

import pytest

np = pytest.importorskip("numpy")
vs = pytest.importorskip("vapoursynth", reason="vapoursynth python module not available")
if not os.environ.get("BORE_PLUGIN"):
    pytest.skip("BORE_PLUGIN env var not set (path to built bore plugin)", allow_module_level=True)

from conftest import degrade, get_planes, gray_clip, make_clip, smooth_image  # noqa: E402

W, H = 64, 48
SIDES = ("top", "bottom", "left", "right")
SINGLE_FILTERS = ("SinglePlane", "SinglePlaneLimited", "SinglePlaneWeighted", "SinglePlaneDebug")
ALL_FILTERS = SINGLE_FILTERS + ("MultiPlane",)
SUPPORTED_FORMATS = ("GRAYS", "YUV444PS", "RGBS")


def border_mask(top=0, bottom=0, left=0, right=0, w=W, h=H):
    m = np.zeros((h, w), dtype=bool)
    if top:
        m[:top, :] = True
    if bottom:
        m[h - bottom:, :] = True
    if left:
        m[:, :left] = True
    if right:
        m[:, w - right:] = True
    return m


def three_plane_clip(core, fmt, degraded_side=None, n=3, frames=1):
    """Correlated 3-plane clip; plane 0 = 0.8*plane1 + 0.2*plane2 (+ border degradation)."""
    p1 = smooth_image(W, H, seed=1)
    p2 = smooth_image(W, H, seed=2, offset=0.1)
    p0 = (0.8 * p1 + 0.2 * p2).astype(np.float32)
    clean = [p0, p1, p2]
    deg = [degrade(p0, **{degraded_side: n}) if degraded_side else p0, p1, p2]
    return clean, deg, make_clip(core, [deg] * frames, getattr(vs, fmt))


# --------------------------------------------------------------------------
# Basic operation on supported formats
# --------------------------------------------------------------------------

@pytest.mark.parametrize("fmt", SUPPORTED_FORMATS)
@pytest.mark.parametrize("name", ALL_FILTERS)
def test_runs_on_supported_formats(core, name, fmt):
    if name == "MultiPlane" and fmt == "GRAYS":
        pytest.skip("MultiPlane requires 3 planes")
    planes = [smooth_image(W, H, seed=s) for s in range(1 if fmt == "GRAYS" else 3)]
    src = make_clip(core, [planes] * 3, getattr(vs, fmt))
    out = getattr(core.bore, name)(src, top=2, bottom=2, left=2, right=2)
    assert out.format.id == src.format.id
    assert (out.width, out.height, out.num_frames) == (W, H, 3)
    for n in range(out.num_frames):
        res = get_planes(out, n)
        for p, arr in enumerate(res):
            assert np.isfinite(arr).all()
            if p != 0:  # only plane 0 is processed by default
                np.testing.assert_array_equal(arr, planes[p])


@pytest.mark.parametrize("name", ALL_FILTERS)
def test_no_borders_is_identity(core, name):
    planes = [smooth_image(W, H, seed=s) for s in range(3)]
    src = make_clip(core, [planes], vs.RGBS)
    out = get_planes(getattr(core.bore, name)(src))
    for a, b in zip(out, planes):
        np.testing.assert_array_equal(a, b)


@pytest.mark.parametrize("name", SINGLE_FILTERS)
def test_subsampled_chroma_plane(core, name):
    y = smooth_image(W, H, seed=3)
    u = degrade(smooth_image(W // 2, H // 2, seed=4), top=2, right=2)
    v = smooth_image(W // 2, H // 2, seed=5)
    src = make_clip(core, [[y, u, v]], vs.YUV420PS)
    mask = gray_clip(core, np.ones((H // 2, W // 2), np.float32))
    for kw in ({}, {"weight_mask": mask}):
        out = get_planes(getattr(core.bore, name)(src, top=2, right=2, plane=1, **kw))
        np.testing.assert_array_equal(out[0], y)
        np.testing.assert_array_equal(out[2], v)
        inner = ~border_mask(top=2, right=2, w=W // 2, h=H // 2)
        np.testing.assert_array_equal(out[1][inner], u[inner])
        assert np.isfinite(out[1]).all()


# --------------------------------------------------------------------------
# Correction quality
# --------------------------------------------------------------------------

CORRECTING = [
    ("SinglePlane", {}),
    ("SinglePlaneLimited", {"ref_line_size": 8}),
    ("SinglePlaneLimited", {}),
    ("SinglePlaneWeighted", {}),
    ("SinglePlaneWeighted", {"ref_line_size": 8, "sigmaS": 5.0, "sigmaR": 0.2, "sigmaD": 0.5}),
]


@pytest.mark.parametrize("side", SIDES)
@pytest.mark.parametrize("name,extra", CORRECTING, ids=lambda v: v if isinstance(v, str) else repr(v))
def test_single_plane_corrects_border(core, name, extra, side):
    clean = smooth_image(W, H)
    deg = degrade(clean, **{side: 3})
    out = get_planes(getattr(core.bore, name)(gray_clip(core, deg), **{side: 3}, **extra))[0]
    m = border_mask(**{side: 3})
    before = np.abs(deg - clean)[m]
    after = np.abs(out - clean)[m]
    assert after.mean() < before.mean() / 6
    assert after.max() < before.max() / 4
    # interior must be bit-exact
    np.testing.assert_array_equal(out[~m], deg[~m])


def test_all_sides_at_once(core):
    clean = smooth_image(W, H, seed=7)
    kw = dict(top=2, bottom=3, left=1, right=4)
    deg = degrade(clean, **kw)
    out = get_planes(core.bore.SinglePlane(gray_clip(core, deg), **kw))[0]
    m = border_mask(**kw)
    assert np.abs(out - clean)[m].mean() < np.abs(deg - clean)[m].mean() / 5
    np.testing.assert_array_equal(out[~m], deg[~m])


def test_exact_scale_recovery(core):
    """A row that is an exact multiple of its neighbour is restored to the neighbour."""
    row = np.linspace(0.2, 0.8, W, dtype=np.float32)
    img = np.tile(row, (H, 1))
    deg = img.copy()
    deg[0] *= 0.5
    deg[-1] *= 2.0
    out = get_planes(core.bore.SinglePlane(gray_clip(core, deg), top=1, bottom=1))[0]
    np.testing.assert_allclose(out, img, rtol=1e-5)


@pytest.mark.parametrize("side", SIDES)
def test_multiplane_corrects_border(core, side):
    clean, deg, src = three_plane_clip(core, "RGBS", side, n=2)
    out = get_planes(core.bore.MultiPlane(src, **{side: 2}))
    m = border_mask(**{side: 2})
    before = np.abs(deg[0] - clean[0])[m]
    after = np.abs(out[0] - clean[0])[m]
    assert after.mean() < before.mean() / 10
    np.testing.assert_array_equal(out[0][~m], deg[0][~m])
    np.testing.assert_array_equal(out[1], deg[1])
    np.testing.assert_array_equal(out[2], deg[2])


@pytest.mark.parametrize("plane", [1, 2])
def test_multiplane_other_plane(core, plane):
    p = [smooth_image(W, H, seed=s) for s in range(3)]
    deg = list(p)
    deg[plane] = degrade(p[plane], top=2)
    out = get_planes(core.bore.MultiPlane(make_clip(core, [deg], vs.YUV444PS), top=2, plane=plane))
    m = border_mask(top=2)
    assert np.abs(out[plane] - p[plane])[m].mean() < np.abs(deg[plane] - p[plane])[m].mean() / 5
    for q in range(3):
        if q != plane:
            np.testing.assert_array_equal(out[q], deg[q])


def test_limited_handles_varying_gain(core):
    """A gain that varies across the line is fixed better with a small ref_line_size."""
    clean = smooth_image(W, H, seed=9)
    deg = clean.copy()
    deg[0] *= np.linspace(0.5, 1.0, W, dtype=np.float32)
    src = gray_clip(core, deg)
    glob = get_planes(core.bore.SinglePlane(src, top=1))[0]
    loc = get_planes(core.bore.SinglePlaneLimited(src, top=1, ref_line_size=4))[0]
    e_glob = np.abs(glob[0] - clean[0]).mean()
    e_loc = np.abs(loc[0] - clean[0]).mean()
    assert e_loc < e_glob / 3


def test_ref_line_size_changes_output(core):
    src = gray_clip(core, degrade(smooth_image(W, H), top=2, left=2))
    a = get_planes(core.bore.SinglePlaneLimited(src, top=2, left=2, ref_line_size=3))[0]
    b = get_planes(core.bore.SinglePlaneLimited(src, top=2, left=2, ref_line_size=20))[0]
    assert not np.array_equal(a, b)


def test_limited_default_ref_line_size_is_100(core):
    src = gray_clip(core, degrade(smooth_image(W, H), top=2, left=2))
    a = get_planes(core.bore.SinglePlaneLimited(src, top=2, left=2))[0]
    b = get_planes(core.bore.SinglePlaneLimited(src, top=2, left=2, ref_line_size=100))[0]
    np.testing.assert_array_equal(a, b)


def test_weighted_defaults(core):
    src = gray_clip(core, degrade(smooth_image(W, H), top=2, right=2))
    a = get_planes(core.bore.SinglePlaneWeighted(src, top=2, right=2))[0]
    b = get_planes(core.bore.SinglePlaneWeighted(
        src, top=2, right=2, sigmaS=50.0, sigmaR=0.5, sigmaD=1.5, ref_line_size=100))[0]
    np.testing.assert_array_equal(a, b)


@pytest.mark.parametrize("param,value", [("sigmaS", 2.0), ("sigmaR", 0.02), ("sigmaD", 0.05), ("ref_line_size", 5)])
def test_weighted_params_change_output(core, param, value):
    clean = smooth_image(W, H, seed=11)
    deg = clean.copy()
    deg[0] *= np.linspace(0.6, 1.0, W, dtype=np.float32)
    deg[:, -1] *= np.linspace(0.7, 0.9, H, dtype=np.float32)
    src = gray_clip(core, deg)
    base = get_planes(core.bore.SinglePlaneWeighted(src, top=1, right=1))[0]
    alt = get_planes(core.bore.SinglePlaneWeighted(src, top=1, right=1, **{param: value}))[0]
    assert not np.array_equal(base, alt)
    assert np.isfinite(alt).all()
    m = border_mask(top=1, right=1)
    np.testing.assert_array_equal(alt[~m], deg[~m])


def test_weighted_huge_sigmas_matches_limited(core):
    """With flat weights the weighted fit degenerates to the limited fit."""
    src = gray_clip(core, degrade(smooth_image(W, H), top=2, bottom=1, left=2, right=1))
    kw = dict(top=2, bottom=1, left=2, right=1, ref_line_size=6)
    a = get_planes(core.bore.SinglePlaneWeighted(src, sigmaS=1e6, sigmaR=1e6, sigmaD=1e6, **kw))[0]
    b = get_planes(core.bore.SinglePlaneLimited(src, **kw))[0]
    np.testing.assert_allclose(a, b, rtol=1e-5, atol=1e-6)


@pytest.mark.xfail(reason="window end clamp `stop = w - 1` drops the last sample of each line, "
                          "so a window covering the whole line differs from SinglePlane", strict=False)
def test_limited_full_window_equals_single(core):
    src = gray_clip(core, degrade(smooth_image(W, H), top=2, left=2))
    a = get_planes(core.bore.SinglePlaneLimited(src, top=2, left=2, ref_line_size=1000))[0]
    b = get_planes(core.bore.SinglePlane(src, top=2, left=2))[0]
    np.testing.assert_allclose(a, b, rtol=1e-6)


# --------------------------------------------------------------------------
# weight_mask
# --------------------------------------------------------------------------

def outlier_case(core):
    clean = smooth_image(W, H, seed=13)
    deg = degrade(clean, top=1)
    deg[0, 20:30] = 5.0  # bright garbage in the border line (e.g. burnt-in logo)
    mask = np.ones((H, W), np.float32)
    mask[:, 20:30] = 0.0
    good = np.ones(W, bool)
    good[20:30] = False
    return clean, deg, mask, good


@pytest.mark.parametrize("name,extra", [
    ("SinglePlane", {}),
    ("SinglePlaneLimited", {"ref_line_size": 100}),
    ("SinglePlaneWeighted", {"sigmaR": 100.0, "sigmaD": 100.0}),
], ids=["SinglePlane", "SinglePlaneLimited", "SinglePlaneWeighted"])
def test_weight_mask_excludes_outliers(core, name, extra):
    clean, deg, mask, good = outlier_case(core)
    src = gray_clip(core, deg)
    f = getattr(core.bore, name)
    plain = get_planes(f(src, top=1, **extra))[0]
    masked = get_planes(f(src, top=1, weight_mask=gray_clip(core, mask), **extra))[0]
    e_plain = np.abs(plain[0, good] - clean[0, good]).mean()
    e_masked = np.abs(masked[0, good] - clean[0, good]).mean()
    assert e_masked < e_plain / 5
    np.testing.assert_array_equal(masked[1:], deg[1:])


@pytest.mark.parametrize("name", ("SinglePlane", "SinglePlaneLimited", "SinglePlaneWeighted"))
def test_weight_mask_all_ones_close_to_unmasked(core, name):
    src = gray_clip(core, degrade(smooth_image(W, H), top=2, bottom=2, left=2, right=2))
    kw = dict(top=2, bottom=2, left=2, right=2)
    f = getattr(core.bore, name)
    a = get_planes(f(src, **kw))[0]
    b = get_planes(f(src, weight_mask=gray_clip(core, np.ones((H, W), np.float32)), **kw))[0]
    if name == "SinglePlaneWeighted":
        # masked variant uses ref (not cur) for the range term, so only approximately equal
        np.testing.assert_allclose(a, b, atol=5e-3)
    else:
        np.testing.assert_allclose(a, b, rtol=1e-5, atol=1e-6)


@pytest.mark.parametrize("name", ("SinglePlane", "SinglePlaneLimited", "SinglePlaneWeighted", "SinglePlaneDebug"))
def test_weight_mask_all_zero_leaves_input(core, name):
    deg = degrade(smooth_image(W, H), top=2, left=2)
    src = gray_clip(core, deg)
    out = get_planes(getattr(core.bore, name)(
        src, top=2, left=2, weight_mask=gray_clip(core, np.zeros((H, W), np.float32))))[0]
    np.testing.assert_array_equal(out, deg)


def test_multiplane_ignores_weight_mask(core):
    _, _, src = three_plane_clip(core, "RGBS", "top", n=2)
    a = get_planes(core.bore.MultiPlane(src, top=2))
    b = get_planes(core.bore.MultiPlane(src, top=2, weight_mask=gray_clip(core, np.zeros((H, W), np.float32))))
    for x, y in zip(a, b):
        np.testing.assert_array_equal(x, y)


# --------------------------------------------------------------------------
# Debug
# --------------------------------------------------------------------------

def test_debug_props(core):
    clean = smooth_image(W, H)
    factors = (0.6, 0.75, 0.88, 0.95)
    kw = dict(top=2, bottom=1, left=3, right=1)
    deg = degrade(clean, factors=factors, **kw)
    out = core.bore.SinglePlaneDebug(gray_clip(core, deg, frames=2), **kw)
    for n in range(2):
        f = out.get_frame(n)
        np.testing.assert_array_equal(np.asarray(f[0]), deg)  # pixels untouched
        adj = np.asarray(f.props["BoreAdjustment"], dtype=float)
        cov = np.asarray(f.props["BoreCovariance"], dtype=float)
        ssq = np.asarray(f.props["BoreSumSquares"], dtype=float)
        total = sum(kw.values())
        assert adj.shape == cov.shape == ssq.shape == (total,)
        assert np.isfinite(adj).all() and np.isfinite(cov).all() and np.isfinite(ssq).all()
        assert (cov > 0).all() and (ssq >= 0).all()
        # Order: top (inner->outer), bottom, left (inner->outer), right; each entry is the
        # least-squares gain cur->neighbour.  Debug applies no corrections, so every line is
        # compared against the *uncorrected* neighbour.
        exp = []
        for row in range(kw["top"] - 1, -1, -1):
            exp.append(np.dot(deg[row], deg[row + 1]) / np.dot(deg[row], deg[row]))
        for row in range(H - kw["bottom"], H):
            exp.append(np.dot(deg[row], deg[row - 1]) / np.dot(deg[row], deg[row]))
        for col in range(kw["left"] - 1, -1, -1):
            exp.append(np.dot(deg[:, col], deg[:, col + 1]) / np.dot(deg[:, col], deg[:, col]))
        for col in range(W - kw["right"], W):
            exp.append(np.dot(deg[:, col], deg[:, col - 1]) / np.dot(deg[:, col], deg[:, col]))
        np.testing.assert_allclose(adj, exp, rtol=1e-5)
        # Sanity: inner top line (row 1 vs clean row 2) ~ 1/f1; outer (row 0 vs row 1) ~ f1/f0.
        assert abs(adj[0] - 1 / factors[1]) < 0.05
        assert abs(adj[1] - factors[1] / factors[0]) < 0.05


def test_debug_matches_singleplane_adjustment(core):
    """Applying the Debug adjustment for a single border line reproduces SinglePlane."""
    deg = degrade(smooth_image(W, H, seed=21), bottom=1)
    src = gray_clip(core, deg)
    adj = core.bore.SinglePlaneDebug(src, bottom=1).get_frame(0).props["BoreAdjustment"]
    out = get_planes(core.bore.SinglePlane(src, bottom=1))[0]
    np.testing.assert_allclose(out[-1], deg[-1] * np.float32(adj), rtol=1e-6)


def test_debug_no_borders_has_no_props(core):
    f = core.bore.SinglePlaneDebug(gray_clip(core, smooth_image(W, H))).get_frame(0)
    assert "BoreAdjustment" not in f.props


def test_debug_weight_mask(core):
    clean, deg, mask, good = outlier_case(core)
    src = gray_clip(core, deg)
    a = core.bore.SinglePlaneDebug(src, top=1).get_frame(0).props["BoreAdjustment"]
    b = core.bore.SinglePlaneDebug(src, top=1, weight_mask=gray_clip(core, mask)).get_frame(0).props["BoreAdjustment"]
    assert a < 0.5  # outlier drags the fit down
    assert abs(b - 1 / 0.6) < 0.05
    z = core.bore.SinglePlaneDebug(src, top=1, weight_mask=gray_clip(core, np.zeros((H, W), np.float32)))
    assert z.get_frame(0).props["BoreAdjustment"] == 1.0  # failed fit reported as 1.0


# --------------------------------------------------------------------------
# Argument / format validation
# --------------------------------------------------------------------------

@pytest.mark.parametrize("name", ALL_FILTERS)
@pytest.mark.parametrize("fmt", ("GRAY8", "GRAY16", "YUV420P8", "YUV444P16", "RGB24"))
def test_rejects_integer_formats(core, name, fmt):
    src = core.std.BlankClip(format=getattr(vs, fmt), width=W, height=H)
    with pytest.raises(vs.Error, match="float"):
        getattr(core.bore, name)(src, top=1)


@pytest.mark.parametrize("name", ALL_FILTERS)
def test_rejects_variable_format(core, name):
    a = core.std.BlankClip(format=vs.YUV444PS, width=W, height=H, length=1)
    b = core.std.BlankClip(format=vs.RGBS, width=W, height=H, length=1)
    src = core.std.Splice([a, b], mismatch=True)
    with pytest.raises(vs.Error):
        getattr(core.bore, name)(src, top=1)


@pytest.mark.parametrize("name", ALL_FILTERS)
def test_rejects_variable_resolution(core, name):
    a = core.std.BlankClip(format=vs.RGBS, width=W, height=H, length=1)
    b = core.std.BlankClip(format=vs.RGBS, width=W * 2, height=H, length=1)
    src = core.std.Splice([a, b], mismatch=True)
    # vsh_isConstantVideoFormat also checks dimensions, so the format error fires first
    with pytest.raises(vs.Error):
        getattr(core.bore, name)(src, top=1)


@pytest.mark.parametrize("name", ALL_FILTERS)
@pytest.mark.parametrize("side,dim", [("top", H), ("bottom", H), ("left", W), ("right", W)])
def test_border_size_limits(core, name, side, dim):
    src = core.std.BlankClip(format=vs.RGBS, width=W, height=H, color=[0.5, 0.4, 0.3])
    f = getattr(core.bore, name)
    f(src, **{side: dim // 2}).get_frame(0)  # max allowed
    with pytest.raises(vs.Error, match=side):
        f(src, **{side: dim // 2 + 1})


@pytest.mark.parametrize("name", SINGLE_FILTERS)
def test_border_limit_uses_plane_dimensions(core, name):
    src = core.std.BlankClip(format=vs.YUV420PS, width=W, height=H)
    f = getattr(core.bore, name)
    f(src, top=H // 4, plane=1).get_frame(0)
    with pytest.raises(vs.Error, match="top"):
        f(src, top=H // 4 + 1, plane=1)
    f(src, top=H // 2, plane=0).get_frame(0)


@pytest.mark.parametrize("name", ALL_FILTERS)
def test_rejects_plane_on_gray(core, name):
    src = core.std.BlankClip(format=vs.GRAYS, width=W, height=H)
    with pytest.raises(vs.Error):
        getattr(core.bore, name)(src, top=1, plane=1)


@pytest.mark.parametrize("ss", [(1, 1), (1, 0), (0, 1), (2, 2), (2, 0)])
def test_multiplane_rejects_subsampled(core, ss):
    fmt = core.query_video_format(vs.YUV, vs.FLOAT, 32, ss[0], ss[1])
    src = core.std.BlankClip(format=fmt, width=W, height=H)
    with pytest.raises(vs.Error, match="RGB or YUV444"):
        core.bore.MultiPlane(src, top=1)


def test_multiplane_rejects_gray(core):
    src = core.std.BlankClip(format=vs.GRAYS, width=W, height=H)
    with pytest.raises(vs.Error, match="3 planes"):
        core.bore.MultiPlane(src, top=1)


@pytest.mark.parametrize("name", ("SinglePlane", "SinglePlaneLimited", "SinglePlaneWeighted", "SinglePlaneDebug"))
@pytest.mark.parametrize("mask_kind,msg", [
    ("int", "float"),
    ("yuv", "gray"),
    ("size", "dimensions"),
])
def test_rejects_bad_weight_mask(core, name, mask_kind, msg):
    src = core.std.BlankClip(format=vs.GRAYS, width=W, height=H)
    mask = {
        "int": core.std.BlankClip(format=vs.GRAY8, width=W, height=H),
        "yuv": core.std.BlankClip(format=vs.YUV444PS, width=W, height=H),
        "size": core.std.BlankClip(format=vs.GRAYS, width=W // 2, height=H),
    }[mask_kind]
    with pytest.raises(vs.Error, match=msg):
        getattr(core.bore, name)(src, top=1, weight_mask=mask)


def test_rejects_full_size_mask_for_chroma_plane(core):
    src = core.std.BlankClip(format=vs.YUV420PS, width=W, height=H)
    with pytest.raises(vs.Error, match="dimensions"):
        core.bore.SinglePlane(src, top=1, plane=1, weight_mask=core.std.BlankClip(format=vs.GRAYS, width=W, height=H))


@pytest.mark.parametrize("name", SINGLE_FILTERS)
def test_rejects_wrong_arg_types(core, name):
    src = core.std.BlankClip(format=vs.GRAYS, width=W, height=H)
    with pytest.raises((vs.Error, ValueError, TypeError)):
        getattr(core.bore, name)(src, top="x")
    with pytest.raises(vs.Error):
        getattr(core.bore, name)(src, nonexistent=1)


@pytest.mark.parametrize("name", ("SinglePlane", "MultiPlane", "SinglePlaneDebug"))
def test_unknown_filter_specific_args_rejected(core, name):
    src = core.std.BlankClip(format=vs.RGBS, width=W, height=H)
    with pytest.raises(vs.Error):
        getattr(core.bore, name)(src, top=1, ref_line_size=5)
    with pytest.raises(vs.Error):
        getattr(core.bore, name)(src, top=1, sigmaS=5.0)


# --------------------------------------------------------------------------
# Known bugs: invalid args that currently crash or are silently accepted.
# Run in a subprocess so a segfault cannot take down the test session.
# --------------------------------------------------------------------------

def run_isolated(expr):
    code = textwrap.dedent(f"""
        import os, vapoursynth as vs
        core = vs.core
        if not hasattr(core, 'bore'):
            core.std.LoadPlugin(os.environ['BORE_PLUGIN'])
        try:
            clip = {expr}
            clip.get_frame(0)
        except vs.Error as e:
            print('VSERROR', e)
        else:
            print('ACCEPTED')
    """)
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=60)
    return r.returncode, r.stdout


INVALID_ARGS_BUGS = {
    "plane_too_big": "core.bore.SinglePlane(core.std.BlankClip(format=vs.YUV444PS, width=64, height=48), top=1, plane=3)",
    "plane_negative": "core.bore.SinglePlane(core.std.BlankClip(format=vs.GRAYS, width=64, height=48), top=1, plane=-1)",
    "multiplane_plane_too_big": "core.bore.MultiPlane(core.std.BlankClip(format=vs.RGBS, width=64, height=48), top=1, plane=3)",
    "ref_line_size_negative": "core.bore.SinglePlaneLimited(core.std.BlankClip(format=vs.GRAYS, width=64, height=48, color=0.5), top=1, ref_line_size=-5)",
    "half_float_input": "core.bore.SinglePlane(core.std.BlankClip(format=vs.GRAYH, width=64, height=48), top=1)",
    "negative_border": "core.bore.SinglePlane(core.std.BlankClip(format=vs.GRAYS, width=64, height=48), top=-1)",
}


@pytest.mark.xfail(reason="plugin does not validate these args: out-of-range plane / negative "
                          "ref_line_size crash (SIGSEGV), 16-bit float input and negative borders "
                          "are silently accepted", strict=False)
@pytest.mark.parametrize("case", sorted(INVALID_ARGS_BUGS))
def test_invalid_args_raise_cleanly(core, case):
    rc, out = run_isolated(INVALID_ARGS_BUGS[case])
    assert rc == 0, f"process crashed (rc={rc})"
    assert out.startswith("VSERROR"), out


@pytest.mark.xfail(reason="sign uses `row > h / 2` / `column > w / 2`: with bottom == h/2 (right == w/2) and even "
                          "size, the innermost border line references the line *outward* instead of inward",
                   strict=False)
@pytest.mark.parametrize("side", ("bottom", "right"))
def test_max_border_references_inner_line(core, side):
    deg = degrade(smooth_image(W, H, seed=31), **{side: 2})
    src = gray_clip(core, deg)
    if side == "bottom":
        line = H // 2
        adj = core.bore.SinglePlaneDebug(src, bottom=H // 2).get_frame(0).props["BoreAdjustment"][0]
        cur, ref = deg[line], deg[line - 1]
    else:
        line = W // 2
        adj = core.bore.SinglePlaneDebug(src, right=W // 2).get_frame(0).props["BoreAdjustment"][0]
        cur, ref = deg[:, line], deg[:, line - 1]
    assert adj == pytest.approx(np.dot(cur, ref) / np.dot(cur, cur), rel=1e-6)
