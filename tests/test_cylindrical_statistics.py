import numpy as np
import pytest

from pyR2D2.util.cylindrical_statistics import CylindricalStatistics


@pytest.fixture
def grid():
    y = np.linspace(-10.0, 10.0, 41)
    z = np.linspace(-8.0, 12.0, 37)
    return y, z


def _mean(cyl_stats, q, **kwargs):
    return cyl_stats.compute(q, stats=["mean"], **kwargs)["mean"]


def _rms(cyl_stats, q, **kwargs):
    return cyl_stats.compute(q, stats=["rms"], **kwargs)["rms"]


# ---------------------------------------------------------------------------
# mean
# ---------------------------------------------------------------------------


def test_uniform_ring_value_matches_bin_index_exactly(grid):
    """A-1: assign each grid point the value of its own radial bin index.

    Every grid point inside a given bin shares the exact same value, so the
    averaged result must reproduce that bin index with no discretization
    error introduced by the binning itself.
    """
    y, z = grid
    yc, zc = 1.5, -2.0
    dr = 0.5

    cyl_stats = CylindricalStatistics(y, z, yc=yc, zc=zc, dr=dr)

    radius = np.hypot(y[:, None] - yc, z[None, :] - zc)
    bin_id = np.floor(radius / dr)
    q = bin_id[None, :, :]

    result = _mean(cyl_stats, q)

    expected = np.where(
        cyl_stats.count > 0,
        np.arange(cyl_stats.nbin, dtype=np.float64),
        np.nan,
    )

    np.testing.assert_allclose(result[0], expected, equal_nan=True)


def test_known_radial_function_matches_independent_reference(grid):
    """A-2: average of a known function of r should match a reference
    computed independently (plain Python accumulation), without relying
    on the class's internal bincount-based implementation.

    We do not compare directly against f(bin center) here. compute()
    computes the exact average of f(r) over whichever grid points happen to
    fall in a bin, and those grid points' radii are spread across the bin
    width rather than sitting exactly at its nominal center. For a
    curved function like this one, that true average differs from
    f(bin center) by an amount that depends on dr and grid resolution,
    not on correctness. A tight tolerance against f(bin center) would
    fail even for a correct implementation; a loose one risks hiding a
    real bug (e.g. an off-by-one in bin assignment) behind the expected
    discretization error. Comparing against an independently computed
    per-bin average sidesteps that confound and lets us use a tight
    tolerance. See A-3 for a complementary, deliberately loose check
    against the function itself.
    """
    y, z = grid
    yc, zc = 0.7, 0.3
    dr = 0.4

    cyl_stats = CylindricalStatistics(y, z, yc=yc, zc=zc, dr=dr)

    radius = np.hypot(y[:, None] - yc, z[None, :] - zc)

    def f(r):
        return np.sin(r) + 2.0 * r

    q = f(radius)[None, :, :]

    result = _mean(cyl_stats, q)

    sums = np.zeros(cyl_stats.nbin)
    counts = np.zeros(cyl_stats.nbin, dtype=np.int64)
    rmax = cyl_stats.rmax
    for r in radius.ravel():
        if r <= rmax:
            b = int(np.floor(r / dr))
            sums[b] += f(r)
            counts[b] += 1

    expected = np.full(cyl_stats.nbin, np.nan)
    valid = counts > 0
    expected[valid] = sums[valid] / counts[valid]

    np.testing.assert_allclose(result[0], expected, rtol=1e-10, equal_nan=True)


def test_mean_of_r_stays_within_bin_width_of_its_center(grid):
    """A-3: sanity check directly against the analytic function f(r) = r.

    By definition of the binning, every grid point assigned to bin b has a
    radius inside [b*dr, (b+1)*dr). So the true average radius within
    that bin must also lie in that range, i.e. within dr/2 of the bin's
    nominal center (`cyl_stats.bin_radius`). This bound holds regardless of how
    grid points are distributed inside the bin, so it is a coarse but always-
    valid check -- unlike a tight tolerance, it needs no assumption about
    grid fineness. It complements A-2 by comparing directly against the
    known function instead of an independently recomputed reference.
    """
    y, z = grid
    dr = 0.5

    cyl_stats = CylindricalStatistics(y, z, dr=dr)

    radius = np.hypot(y[:, None], z[None, :])
    q = radius[None, :, :]

    result = _mean(cyl_stats, q)[0]

    valid = cyl_stats.count > 0
    diff = np.abs(result[valid] - cyl_stats.bin_radius[valid])

    assert np.all(diff <= dr / 2)


def test_2d_input_matches_3d_single_plane_equivalent(grid):
    """A 2D q (shape (ny, nz)) should be treated as a single x-plane and
    return a 1D (nbin,) result, identical to calling compute() with the
    equivalent 3D input q[None, :, :] and taking its only row.
    """
    y, z = grid
    yc, zc = -0.4, 1.1
    dr = 0.6

    cyl_stats = CylindricalStatistics(y, z, yc=yc, zc=zc, dr=dr)

    radius = np.hypot(y[:, None] - yc, z[None, :] - zc)
    q_2d = np.sin(radius) + radius

    result_2d = _mean(cyl_stats, q_2d)
    result_3d = _mean(cyl_stats, q_2d[None, :, :])

    assert result_2d.shape == (cyl_stats.nbin,)
    assert result_3d.shape == (1, cyl_stats.nbin)
    np.testing.assert_array_equal(result_2d, result_3d[0])


def test_invalid_ndim_raises_value_error(grid):
    """q with ndim other than 2 or 3 (e.g. 1D or 4D) must be rejected."""
    y, z = grid
    cyl_stats = CylindricalStatistics(y, z, dr=1.0)

    with pytest.raises(ValueError):
        cyl_stats.compute(np.zeros(cyl_stats.jx))

    with pytest.raises(ValueError):
        cyl_stats.compute(np.zeros((2, cyl_stats.jx, cyl_stats.kx, 3)))


def test_mismatched_yz_shape_raises_value_error(grid):
    """C-2: q whose (ny, nz) doesn't match the coordinates given to
    __init__ must be rejected, for both 2D and 3D input.
    """
    y, z = grid
    cyl_stats = CylindricalStatistics(y, z, dr=1.0)

    with pytest.raises(ValueError):
        cyl_stats.compute(np.zeros((cyl_stats.jx + 1, cyl_stats.kx)))

    with pytest.raises(ValueError):
        cyl_stats.compute(np.zeros((3, cyl_stats.jx, cyl_stats.kx - 1)))


def test_partial_nan_plane_averages_only_finite_grid_points(grid):
    """B-1: grid points marked NaN must be excluded from their bin's average.

    Exercises the branch in compute() where some but not all values in
    the chunk are finite, and checks against an independently computed
    reference that also skips NaNs.
    """
    y, z = grid
    yc, zc = 0.2, -0.6
    dr = 0.5

    cyl_stats = CylindricalStatistics(y, z, yc=yc, zc=zc, dr=dr)

    radius = np.hypot(y[:, None] - yc, z[None, :] - zc)
    q = radius.copy()
    nan_mask = (np.arange(q.size) % 3 == 0).reshape(q.shape)
    q[nan_mask] = np.nan

    result = _mean(cyl_stats, q)

    sums = np.zeros(cyl_stats.nbin)
    counts = np.zeros(cyl_stats.nbin, dtype=np.int64)
    for r, v in zip(radius.ravel(), q.ravel()):
        if r <= cyl_stats.rmax and np.isfinite(v):
            b = int(np.floor(r / dr))
            sums[b] += v
            counts[b] += 1

    expected = np.full(cyl_stats.nbin, np.nan)
    valid = counts > 0
    expected[valid] = sums[valid] / counts[valid]

    np.testing.assert_allclose(result, expected, rtol=1e-10, equal_nan=True)


def test_all_nan_plane_returns_all_nan(grid):
    """B-2: a plane where every grid point is NaN must produce an all-NaN
    result for every bin, including bins that would otherwise have data.
    """
    y, z = grid
    cyl_stats = CylindricalStatistics(y, z, dr=0.7)

    q = np.full((cyl_stats.jx, cyl_stats.kx), np.nan)
    result = _mean(cyl_stats, q)

    assert np.all(np.isnan(result))


def test_mixed_finite_and_nan_planes_in_same_chunk(grid):
    """B-3: within a single chunk, some x-planes are fully finite and
    one contains NaNs. The chunk-wide `isfinite().all()` branch falls
    back to per-grid-point counting for the whole chunk, so this checks that
    the finite planes still come out correct (not just the NaN plane).
    """
    y, z = grid
    yc, zc = -0.3, 0.5
    dr = 0.5
    cyl_stats = CylindricalStatistics(y, z, yc=yc, zc=zc, dr=dr)

    radius = np.hypot(y[:, None] - yc, z[None, :] - zc)

    n_planes = 5
    q = np.stack([radius * (i + 1) for i in range(n_planes)])

    nan_mask = (np.arange(radius.size) % 4 == 0).reshape(radius.shape)
    q[2][nan_mask] = np.nan

    result = _mean(cyl_stats, q, chunk_size=n_planes)

    for i in range(n_planes):
        sums = np.zeros(cyl_stats.nbin)
        counts = np.zeros(cyl_stats.nbin, dtype=np.int64)
        for r, v in zip(radius.ravel(), q[i].ravel()):
            if r <= cyl_stats.rmax and np.isfinite(v):
                b = int(np.floor(r / dr))
                sums[b] += v
                counts[b] += 1

        expected = np.full(cyl_stats.nbin, np.nan)
        valid = counts > 0
        expected[valid] = sums[valid] / counts[valid]

        np.testing.assert_allclose(
            result[i],
            expected,
            rtol=1e-10,
            equal_nan=True,
            err_msg=f"plane {i}",
        )


def test_dr_variation_recomputes_nbin_and_radius(grid):
    """D-1: nbin and the bin-center radius array must match the given dr,
    for a fixed rmax.
    """
    y, z = grid
    rmax = 5.0

    for dr in (0.3, 0.5, 1.0, 1.7):
        cyl_stats = CylindricalStatistics(y, z, dr=dr, rmax=rmax)

        expected_nbin = int(np.floor(rmax / dr)) + 1
        assert cyl_stats.nbin == expected_nbin

        expected_radius = (np.arange(expected_nbin, dtype=np.float64) + 0.5) * dr
        np.testing.assert_allclose(cyl_stats.bin_radius, expected_radius)


def test_rmax_defaults_to_grid_max_radius(grid):
    """D-2: rmax=None must default to the grid's actual max radius, while
    an explicit rmax must be kept exactly as given.
    """
    y, z = grid
    yc, zc = 0.3, -0.2
    dr = 0.5

    radius = np.hypot(y[:, None] - yc, z[None, :] - zc)

    cyl_stats_auto = CylindricalStatistics(y, z, yc=yc, zc=zc, dr=dr)
    assert cyl_stats_auto.rmax == pytest.approx(radius.max())

    explicit_rmax = 3.0
    cyl_stats_explicit = CylindricalStatistics(
        y, z, yc=yc, zc=zc, dr=dr, rmax=explicit_rmax
    )
    assert cyl_stats_explicit.rmax == explicit_rmax


def test_rmax_smaller_than_data_range_excludes_outer_grid_points(grid):
    """D-3: grid points farther than rmax must be excluded, so
    valid_grid_index and count only account for grid points within rmax.
    """
    y, z = grid
    dr = 0.5

    full_radius = np.hypot(y[:, None], z[None, :])
    rmax = full_radius.max() / 2  # well inside the data range

    cyl_stats = CylindricalStatistics(y, z, dr=dr, rmax=rmax)

    expected_included = np.count_nonzero(full_radius <= rmax)
    assert cyl_stats.valid_grid_index.size == expected_included
    assert cyl_stats.count.sum() == expected_included
    assert cyl_stats.nbin == int(np.floor(rmax / dr)) + 1


def test_chunk_size_does_not_affect_result(grid):
    """E: compute() must produce identical results regardless of chunk_size.

    Chunking is purely a memory/performance knob. Splitting the same
    input across different chunk boundaries must not change any value,
    including which planes end up grouped with a NaN-containing plane
    (which flips the finite-only vs NaN-aware branch for the whole
    chunk).
    """
    y, z = grid
    yc, zc = 0.1, -0.3
    dr = 0.4
    cyl_stats = CylindricalStatistics(y, z, yc=yc, zc=zc, dr=dr)

    radius = np.hypot(y[:, None] - yc, z[None, :] - zc)

    n_planes = 11
    q = np.stack([np.sin(radius) * (i + 1) + radius for i in range(n_planes)])

    # NaNs in some (not all) planes, so different chunk_size values mix
    # NaN and finite-only planes into chunks differently.
    nan_mask = (np.arange(radius.size) % 5 == 0).reshape(radius.shape)
    for i in (2, 3, 7):
        q[i][nan_mask] = np.nan

    baseline = _mean(cyl_stats, q, chunk_size=n_planes)  # everything in one chunk

    for chunk_size in (1, 2, 3, 4, 6, 100):
        result = _mean(cyl_stats, q, chunk_size=chunk_size)
        np.testing.assert_array_equal(
            result, baseline, err_msg=f"chunk_size={chunk_size}"
        )


def test_non_positive_dr_raises_value_error(grid):
    """dr must be strictly positive."""
    y, z = grid

    with pytest.raises(ValueError):
        CylindricalStatistics(y, z, dr=0.0)

    with pytest.raises(ValueError):
        CylindricalStatistics(y, z, dr=-1.0)


def test_dr_defaults_to_y_grid_spacing(grid):
    """dr, when omitted, must default to the mean spacing of `y` rather
    than an arbitrary fixed value.
    """
    y, z = grid

    cyl_stats = CylindricalStatistics(y, z)

    assert cyl_stats.dr == pytest.approx(np.mean(np.abs(np.diff(y))))


def test_non_positive_rmax_raises_value_error(grid):
    """An explicit rmax must be strictly positive."""
    y, z = grid

    with pytest.raises(ValueError):
        CylindricalStatistics(y, z, dr=0.5, rmax=0.0)

    with pytest.raises(ValueError):
        CylindricalStatistics(y, z, dr=0.5, rmax=-2.0)


def test_non_positive_chunk_size_raises_value_error(grid):
    """chunk_size must be strictly positive; previously a non-positive
    value silently produced an all-NaN result instead of erroring.
    """
    y, z = grid
    cyl_stats = CylindricalStatistics(y, z, dr=0.5)
    q = np.zeros((cyl_stats.jx, cyl_stats.kx))

    with pytest.raises(ValueError):
        cyl_stats.compute(q, chunk_size=0)

    with pytest.raises(ValueError):
        cyl_stats.compute(q, chunk_size=-5)


# ---------------------------------------------------------------------------
# stats= validation
# ---------------------------------------------------------------------------


def test_unknown_stat_name_raises_value_error(grid):
    y, z = grid
    cyl_stats = CylindricalStatistics(y, z, dr=0.5)
    q = np.zeros((cyl_stats.jx, cyl_stats.kx))

    with pytest.raises(ValueError):
        cyl_stats.compute(q, stats=["mean", "median"])


def test_empty_stats_raises_value_error(grid):
    y, z = grid
    cyl_stats = CylindricalStatistics(y, z, dr=0.5)
    q = np.zeros((cyl_stats.jx, cyl_stats.kx))

    with pytest.raises(ValueError):
        cyl_stats.compute(q, stats=[])


# ---------------------------------------------------------------------------
# rms (population standard deviation of qq around its own bin mean)
# ---------------------------------------------------------------------------


def test_rms_uniform_value_per_bin_is_zero(grid):
    """R-1: every grid point inside a given bin shares the exact same
    value, so each bin's own mean equals that value exactly and the rms
    (deviation from the bin's own mean) must be exactly 0.
    """
    y, z = grid
    yc, zc = 1.5, -2.0
    dr = 0.5

    cyl_stats = CylindricalStatistics(y, z, yc=yc, zc=zc, dr=dr)

    radius = np.hypot(y[:, None] - yc, z[None, :] - zc)
    bin_id = np.floor(radius / dr)
    q = bin_id[None, :, :]

    result = _rms(cyl_stats, q)

    expected = np.where(cyl_stats.count > 0, 0.0, np.nan)

    np.testing.assert_allclose(result[0], expected, atol=1e-9, equal_nan=True)


def _bin_included(r, v, rmax, cond_v):
    """Mirrors compute()'s grid-point inclusion rule: within rmax, finite,
    and (when a condition value is given) finite and positive.
    """
    if r > rmax or not np.isfinite(v):
        return False
    if cond_v is not None and not (np.isfinite(cond_v) and cond_v > 0):
        return False
    return True


def _reference_bin_mean(radius, q, rmax, dr, nbin, cond=None):
    sums = np.zeros(nbin)
    counts = np.zeros(nbin, dtype=np.int64)
    cond_flat = cond.ravel() if cond is not None else None
    for i, (r, v) in enumerate(zip(radius.ravel(), q.ravel())):
        cond_v = cond_flat[i] if cond_flat is not None else None
        if _bin_included(r, v, rmax, cond_v):
            b = int(np.floor(r / dr))
            sums[b] += v
            counts[b] += 1
    bin_mean = np.full(nbin, np.nan)
    valid = counts > 0
    bin_mean[valid] = sums[valid] / counts[valid]
    return bin_mean, counts


def _reference_bin_rms(radius, q, rmax, dr, nbin, cond=None):
    bin_mean, counts = _reference_bin_mean(radius, q, rmax, dr, nbin, cond=cond)
    sq_sums = np.zeros(nbin)
    cond_flat = cond.ravel() if cond is not None else None
    for i, (r, v) in enumerate(zip(radius.ravel(), q.ravel())):
        cond_v = cond_flat[i] if cond_flat is not None else None
        if _bin_included(r, v, rmax, cond_v):
            b = int(np.floor(r / dr))
            sq_sums[b] += (v - bin_mean[b]) ** 2
    expected = np.full(nbin, np.nan)
    valid = counts > 0
    expected[valid] = np.sqrt(sq_sums[valid] / counts[valid])
    return expected


def test_rms_known_radial_function_matches_independent_reference(grid):
    """R-2: rms of a known function of r should match a reference computed
    independently (plain Python accumulation over each bin's own mean),
    analogous to A-2 for mean().
    """
    y, z = grid
    yc, zc = 0.7, 0.3
    dr = 0.4

    cyl_stats = CylindricalStatistics(y, z, yc=yc, zc=zc, dr=dr)

    radius = np.hypot(y[:, None] - yc, z[None, :] - zc)

    def f(r):
        return np.sin(r) + 2.0 * r

    q = f(radius)
    result = _rms(cyl_stats, q)

    expected = _reference_bin_rms(radius, q, cyl_stats.rmax, dr, cyl_stats.nbin)

    np.testing.assert_allclose(result, expected, rtol=1e-9, equal_nan=True)


def test_rms_of_r_stays_within_half_bin_width(grid):
    """R-3: the population std of values confined to an interval of width
    dr is bounded above by dr/2 (the bound is tight only for a two-point
    mass at the interval's endpoints), so rms(r) must satisfy this bound
    in every populated bin. Complements R-2 with a coarse, always-valid
    check that needs no independent reference computation.
    """
    y, z = grid
    dr = 0.5

    cyl_stats = CylindricalStatistics(y, z, dr=dr)

    radius = np.hypot(y[:, None], z[None, :])
    q = radius[None, :, :]

    result = _rms(cyl_stats, q)[0]

    valid = cyl_stats.count > 0
    assert np.all(result[valid] <= dr / 2 + 1e-9)


def test_rms_2d_input_matches_3d_single_plane_equivalent(grid):
    y, z = grid
    yc, zc = -0.4, 1.1
    dr = 0.6

    cyl_stats = CylindricalStatistics(y, z, yc=yc, zc=zc, dr=dr)

    radius = np.hypot(y[:, None] - yc, z[None, :] - zc)
    q_2d = np.sin(radius) + radius

    result_2d = _rms(cyl_stats, q_2d)
    result_3d = _rms(cyl_stats, q_2d[None, :, :])

    assert result_2d.shape == (cyl_stats.nbin,)
    assert result_3d.shape == (1, cyl_stats.nbin)
    np.testing.assert_array_equal(result_2d, result_3d[0])


def test_rms_partial_nan_plane_averages_only_finite_grid_points(grid):
    y, z = grid
    yc, zc = 0.2, -0.6
    dr = 0.5

    cyl_stats = CylindricalStatistics(y, z, yc=yc, zc=zc, dr=dr)

    radius = np.hypot(y[:, None] - yc, z[None, :] - zc)
    q = radius.copy()
    nan_mask = (np.arange(q.size) % 3 == 0).reshape(q.shape)
    q[nan_mask] = np.nan

    result = _rms(cyl_stats, q)
    expected = _reference_bin_rms(radius, q, cyl_stats.rmax, dr, cyl_stats.nbin)

    np.testing.assert_allclose(result, expected, rtol=1e-9, equal_nan=True)


def test_rms_all_nan_plane_returns_all_nan(grid):
    y, z = grid
    cyl_stats = CylindricalStatistics(y, z, dr=0.7)

    q = np.full((cyl_stats.jx, cyl_stats.kx), np.nan)
    result = _rms(cyl_stats, q)

    assert np.all(np.isnan(result))


def test_rms_mixed_finite_and_nan_planes_in_same_chunk(grid):
    y, z = grid
    yc, zc = -0.3, 0.5
    dr = 0.5
    cyl_stats = CylindricalStatistics(y, z, yc=yc, zc=zc, dr=dr)

    radius = np.hypot(y[:, None] - yc, z[None, :] - zc)

    n_planes = 5
    q = np.stack([radius * (i + 1) for i in range(n_planes)])

    nan_mask = (np.arange(radius.size) % 4 == 0).reshape(radius.shape)
    q[2][nan_mask] = np.nan

    result = _rms(cyl_stats, q, chunk_size=n_planes)

    for i in range(n_planes):
        expected = _reference_bin_rms(radius, q[i], cyl_stats.rmax, dr, cyl_stats.nbin)
        np.testing.assert_allclose(
            result[i],
            expected,
            rtol=1e-9,
            equal_nan=True,
            err_msg=f"plane {i}",
        )


def test_rms_chunk_size_does_not_affect_result(grid):
    y, z = grid
    yc, zc = 0.1, -0.3
    dr = 0.4
    cyl_stats = CylindricalStatistics(y, z, yc=yc, zc=zc, dr=dr)

    radius = np.hypot(y[:, None] - yc, z[None, :] - zc)

    n_planes = 11
    q = np.stack([np.sin(radius) * (i + 1) + radius for i in range(n_planes)])

    nan_mask = (np.arange(radius.size) % 5 == 0).reshape(radius.shape)
    for i in (2, 3, 7):
        q[i][nan_mask] = np.nan

    baseline = _rms(cyl_stats, q, chunk_size=n_planes)

    for chunk_size in (1, 2, 3, 4, 6, 100):
        result = _rms(cyl_stats, q, chunk_size=chunk_size)
        np.testing.assert_array_equal(
            result, baseline, err_msg=f"chunk_size={chunk_size}"
        )


# ---------------------------------------------------------------------------
# combined mean+rms
# ---------------------------------------------------------------------------


def test_compute_mean_and_rms_together_matches_separate_calls(grid):
    """Requesting both stats in one call must give exactly the same
    values as requesting each separately -- the shared sums/counts
    bincounts are purely an internal reuse and must not change results.
    """
    y, z = grid
    yc, zc = 0.4, -0.7
    dr = 0.5
    cyl_stats = CylindricalStatistics(y, z, yc=yc, zc=zc, dr=dr)

    radius = np.hypot(y[:, None] - yc, z[None, :] - zc)
    n_planes = 6
    q = np.stack([np.sin(radius) * (i + 1) + radius for i in range(n_planes)])
    nan_mask = (np.arange(radius.size) % 4 == 0).reshape(radius.shape)
    q[1][nan_mask] = np.nan

    combined = cyl_stats.compute(q, stats=["mean", "rms"], chunk_size=4)
    separate_mean = _mean(cyl_stats, q, chunk_size=4)
    separate_rms = _rms(cyl_stats, q, chunk_size=4)

    np.testing.assert_array_equal(combined["mean"], separate_mean)
    np.testing.assert_array_equal(combined["rms"], separate_rms)


# ---------------------------------------------------------------------------
# conditional averaging (condition / condition_field)
# ---------------------------------------------------------------------------


def test_condition_defaults_to_qq_itself(grid):
    """condition_field defaults to qq: only grid points with qq > 0
    contribute to the bin's mean/rms.
    """
    y, z = grid
    yc, zc = 0.3, -0.5
    dr = 0.4
    cyl_stats = CylindricalStatistics(y, z, yc=yc, zc=zc, dr=dr)

    radius = np.hypot(y[:, None] - yc, z[None, :] - zc)
    q = np.sin(radius) * 3.0 - 1.0  # mix of positive and negative values

    result = cyl_stats.compute(q, stats=["mean", "rms"], condition=lambda v: v > 0)

    expected_mean, _ = _reference_bin_mean(
        radius, q, cyl_stats.rmax, dr, cyl_stats.nbin, cond=q
    )
    expected_rms = _reference_bin_rms(
        radius, q, cyl_stats.rmax, dr, cyl_stats.nbin, cond=q
    )

    np.testing.assert_allclose(result["mean"], expected_mean, rtol=1e-9, equal_nan=True)
    np.testing.assert_allclose(result["rms"], expected_rms, rtol=1e-9, equal_nan=True)


def test_condition_on_separate_field(grid):
    """condition_field distinct from qq: average qq only where a second
    field (e.g. a velocity component) is positive.
    """
    y, z = grid
    yc, zc = -0.2, 0.6
    dr = 0.5
    cyl_stats = CylindricalStatistics(y, z, yc=yc, zc=zc, dr=dr)

    radius = np.hypot(y[:, None] - yc, z[None, :] - zc)
    q = np.cos(radius) + radius
    vz = np.sin(radius * 2.0)  # unrelated sign pattern from q

    result = cyl_stats.compute(
        q, stats=["mean", "rms"], condition_field=vz, condition=lambda v: v > 0
    )

    expected_mean, _ = _reference_bin_mean(
        radius, q, cyl_stats.rmax, dr, cyl_stats.nbin, cond=vz
    )
    expected_rms = _reference_bin_rms(
        radius, q, cyl_stats.rmax, dr, cyl_stats.nbin, cond=vz
    )

    np.testing.assert_allclose(result["mean"], expected_mean, rtol=1e-9, equal_nan=True)
    np.testing.assert_allclose(result["rms"], expected_rms, rtol=1e-9, equal_nan=True)


def test_condition_field_nan_is_always_excluded(grid):
    """NaN in condition_field must be excluded even if `condition` itself
    would return True for it (e.g. a condition that tolerates NaN).
    """
    y, z = grid
    dr = 0.5
    cyl_stats = CylindricalStatistics(y, z, dr=dr)

    radius = np.hypot(y[:, None], z[None, :])
    q = radius.copy()

    vz = np.ones_like(radius)
    nan_mask = (np.arange(vz.size) % 5 == 0).reshape(vz.shape)
    vz[nan_mask] = np.nan

    # A condition that would (incorrectly, if not overridden) accept NaN.
    def tolerant_condition(v):
        return np.isnan(v) | (v > 0)

    result = cyl_stats.compute(
        q, stats=["mean"], condition_field=vz, condition=tolerant_condition
    )["mean"]

    expected, _ = _reference_bin_mean(
        radius, q, cyl_stats.rmax, dr, cyl_stats.nbin, cond=vz
    )

    np.testing.assert_allclose(result, expected, rtol=1e-9, equal_nan=True)


def test_condition_field_without_condition_raises_value_error(grid):
    y, z = grid
    cyl_stats = CylindricalStatistics(y, z, dr=0.5)
    q = np.zeros((cyl_stats.jx, cyl_stats.kx))

    with pytest.raises(ValueError):
        cyl_stats.compute(q, condition_field=q)


def test_condition_field_shape_mismatch_raises_value_error(grid):
    y, z = grid
    cyl_stats = CylindricalStatistics(y, z, dr=0.5)
    q = np.zeros((cyl_stats.jx, cyl_stats.kx))
    wrong_shape_field = np.zeros((cyl_stats.jx, cyl_stats.kx - 1))

    with pytest.raises(ValueError):
        cyl_stats.compute(
            q, condition_field=wrong_shape_field, condition=lambda v: v > 0
        )


def test_condition_chunk_size_does_not_affect_result(grid):
    y, z = grid
    yc, zc = 0.1, -0.3
    dr = 0.4
    cyl_stats = CylindricalStatistics(y, z, yc=yc, zc=zc, dr=dr)

    radius = np.hypot(y[:, None] - yc, z[None, :] - zc)
    n_planes = 11
    q = np.stack([np.sin(radius) * (i + 1) + radius for i in range(n_planes)])
    vz = np.stack([np.cos(radius) * (i + 1) - radius for i in range(n_planes)])

    baseline = cyl_stats.compute(
        q,
        stats=["mean", "rms"],
        condition_field=vz,
        condition=lambda v: v > 0,
        chunk_size=n_planes,
    )

    for chunk_size in (1, 2, 3, 4, 6, 100):
        result = cyl_stats.compute(
            q,
            stats=["mean", "rms"],
            condition_field=vz,
            condition=lambda v: v > 0,
            chunk_size=chunk_size,
        )
        np.testing.assert_array_equal(
            result["mean"], baseline["mean"], err_msg=f"chunk_size={chunk_size}"
        )
        np.testing.assert_array_equal(
            result["rms"], baseline["rms"], err_msg=f"chunk_size={chunk_size}"
        )


# ---------------------------------------------------------------------------
# cov / corr (population covariance and Pearson correlation between qq, qq2)
# ---------------------------------------------------------------------------


def _reference_bin_cov_corr(radius, q, q2, rmax, dr, nbin, cond=None):
    """Independent reference for cov/corr: a grid point contributes only
    where q, q2 (and cond, if given) are all finite and r <= rmax -- the
    "pairwise valid" set compute() uses internally for cov/corr, which
    may differ from the set mean()/rms() would use for q alone.
    """
    r_flat = radius.ravel()
    q_flat = q.ravel()
    q2_flat = q2.ravel()
    cond_flat = cond.ravel() if cond is not None else None

    def included(i):
        if (
            r_flat[i] > rmax
            or not np.isfinite(q_flat[i])
            or not np.isfinite(q2_flat[i])
        ):
            return False
        if cond_flat is not None:
            c = cond_flat[i]
            if not (np.isfinite(c) and c > 0):
                return False
        return True

    sums_q = np.zeros(nbin)
    sums_q2 = np.zeros(nbin)
    counts = np.zeros(nbin, dtype=np.int64)
    for i in range(r_flat.size):
        if included(i):
            b = int(np.floor(r_flat[i] / dr))
            sums_q[b] += q_flat[i]
            sums_q2[b] += q2_flat[i]
            counts[b] += 1

    valid = counts > 0
    mean_q = np.full(nbin, np.nan)
    mean_q[valid] = sums_q[valid] / counts[valid]
    mean_q2 = np.full(nbin, np.nan)
    mean_q2[valid] = sums_q2[valid] / counts[valid]

    sq_cov = np.zeros(nbin)
    sq_var_q = np.zeros(nbin)
    sq_var_q2 = np.zeros(nbin)
    for i in range(r_flat.size):
        if included(i):
            b = int(np.floor(r_flat[i] / dr))
            dq = q_flat[i] - mean_q[b]
            dq2 = q2_flat[i] - mean_q2[b]
            sq_cov[b] += dq * dq2
            sq_var_q[b] += dq**2
            sq_var_q2[b] += dq2**2

    cov = np.full(nbin, np.nan)
    cov[valid] = sq_cov[valid] / counts[valid]

    var_q = np.full(nbin, np.nan)
    var_q[valid] = sq_var_q[valid] / counts[valid]
    var_q2 = np.full(nbin, np.nan)
    var_q2[valid] = sq_var_q2[valid] / counts[valid]

    denom = np.sqrt(var_q * var_q2)
    corr = np.full(nbin, np.nan)
    nonzero = valid & (denom > 0)
    corr[nonzero] = cov[nonzero] / denom[nonzero]

    return cov, corr


def test_cov_known_fields_matches_independent_reference(grid):
    y, z = grid
    yc, zc = 0.3, -0.4
    dr = 0.5
    cyl_stats = CylindricalStatistics(y, z, yc=yc, zc=zc, dr=dr)

    radius = np.hypot(y[:, None] - yc, z[None, :] - zc)
    q = np.sin(radius) + radius
    q2 = np.cos(radius) - 0.5 * radius

    result = cyl_stats.compute(q, qq2=q2, stats=["cov"])["cov"]
    expected_cov, _ = _reference_bin_cov_corr(
        radius, q, q2, cyl_stats.rmax, dr, cyl_stats.nbin
    )

    np.testing.assert_allclose(result, expected_cov, rtol=1e-9, equal_nan=True)


def test_corr_known_fields_matches_independent_reference(grid):
    y, z = grid
    yc, zc = -0.6, 0.2
    dr = 0.4
    cyl_stats = CylindricalStatistics(y, z, yc=yc, zc=zc, dr=dr)

    radius = np.hypot(y[:, None] - yc, z[None, :] - zc)
    q = np.sin(radius) + radius
    q2 = np.cos(radius * 1.3) - 0.5 * radius

    result = cyl_stats.compute(q, qq2=q2, stats=["corr"])["corr"]
    _, expected_corr = _reference_bin_cov_corr(
        radius, q, q2, cyl_stats.rmax, dr, cyl_stats.nbin
    )

    np.testing.assert_allclose(result, expected_corr, rtol=1e-9, equal_nan=True)


def test_corr_of_field_with_itself_is_one(grid):
    """corr(q, q) must be exactly 1 in every bin with nonzero variance."""
    y, z = grid
    dr = 0.5
    cyl_stats = CylindricalStatistics(y, z, dr=dr)

    radius = np.hypot(y[:, None], z[None, :])
    q = np.sin(radius) * 2.0 + radius

    result = cyl_stats.compute(q, qq2=q, stats=["corr"])["corr"]

    finite = np.isfinite(result)
    assert finite.any()
    np.testing.assert_allclose(result[finite], 1.0, atol=1e-9)


def test_mean_unaffected_by_cov_when_qq2_has_different_nan_pattern(grid):
    """Requesting stats=["mean", "cov"] together must give the same
    "mean" as requesting ["mean"] alone, even though qq2 is NaN at
    different grid points than qq -- mean uses qq's own valid set, cov
    uses the narrower pairwise (qq and qq2 both finite) set.
    """
    y, z = grid
    yc, zc = 0.2, -0.3
    dr = 0.5
    cyl_stats = CylindricalStatistics(y, z, yc=yc, zc=zc, dr=dr)

    radius = np.hypot(y[:, None] - yc, z[None, :] - zc)
    q = radius.copy()
    q2 = radius * 2.0

    q_nan_mask = (np.arange(q.size) % 5 == 0).reshape(q.shape)
    q2_nan_mask = (np.arange(q.size) % 7 == 0).reshape(q.shape)
    q[q_nan_mask] = np.nan
    q2[q2_nan_mask] = np.nan

    mean_alone = _mean(cyl_stats, q)
    combined = cyl_stats.compute(q, qq2=q2, stats=["mean", "cov"])

    np.testing.assert_array_equal(combined["mean"], mean_alone)

    expected_cov, _ = _reference_bin_cov_corr(
        radius, q, q2, cyl_stats.rmax, dr, cyl_stats.nbin
    )
    np.testing.assert_allclose(combined["cov"], expected_cov, rtol=1e-9, equal_nan=True)


def test_cov_corr_without_qq2_raises_value_error(grid):
    y, z = grid
    cyl_stats = CylindricalStatistics(y, z, dr=0.5)
    q = np.zeros((cyl_stats.jx, cyl_stats.kx))

    with pytest.raises(ValueError):
        cyl_stats.compute(q, stats=["cov"])

    with pytest.raises(ValueError):
        cyl_stats.compute(q, stats=["corr"])


def test_qq2_shape_mismatch_raises_value_error(grid):
    y, z = grid
    cyl_stats = CylindricalStatistics(y, z, dr=0.5)
    q = np.zeros((cyl_stats.jx, cyl_stats.kx))
    wrong_shape_field = np.zeros((cyl_stats.jx, cyl_stats.kx - 1))

    with pytest.raises(ValueError):
        cyl_stats.compute(q, qq2=wrong_shape_field, stats=["cov"])


def test_cov_corr_chunk_size_does_not_affect_result(grid):
    y, z = grid
    yc, zc = 0.1, -0.3
    dr = 0.4
    cyl_stats = CylindricalStatistics(y, z, yc=yc, zc=zc, dr=dr)

    radius = np.hypot(y[:, None] - yc, z[None, :] - zc)
    n_planes = 11
    q = np.stack([np.sin(radius) * (i + 1) + radius for i in range(n_planes)])
    q2 = np.stack([np.cos(radius) * (i + 1) - radius for i in range(n_planes)])

    baseline = cyl_stats.compute(q, qq2=q2, stats=["cov", "corr"], chunk_size=n_planes)

    for chunk_size in (1, 2, 3, 4, 6, 100):
        result = cyl_stats.compute(
            q, qq2=q2, stats=["cov", "corr"], chunk_size=chunk_size
        )
        np.testing.assert_array_equal(
            result["cov"], baseline["cov"], err_msg=f"chunk_size={chunk_size}"
        )
        np.testing.assert_array_equal(
            result["corr"], baseline["corr"], err_msg=f"chunk_size={chunk_size}"
        )
