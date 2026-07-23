from collections.abc import Callable, Sequence

import numpy as np

_VALID_STATS = frozenset({"mean", "rms", "cov", "corr"})
_PAIR_STATS = frozenset({"cov", "corr"})


class CylindricalStatistics:
    """
    Computes per-x-plane statistics of qq[x, y, z] over concentric radial
    bins around a horizontal (y, z) center, producing e.g. qq_mean[x, r],
    qq_rms[x, r], and (given a second field) their covariance/correlation.
    """

    def __init__(
        self,
        y: np.ndarray,
        z: np.ndarray,
        *,
        yc: float = 0.0,
        zc: float = 0.0,
        dr: float = 1.0,
        rmax: float | None = None,
    ) -> None:
        """
        Parameters
        ----------
        y : ndarray
            1D coordinate array along the horizontal y-axis.
        z : ndarray
            1D coordinate array along the horizontal z-axis.
        yc : float, optional
            Center of the polar coordinate system along y.
        zc : float, optional
            Center of the polar coordinate system along z.
        dr : float, optional
            Width of each radial bin, in the same units as y/z.
            Determines the resolution of the r-axis in compute()'s
            output (see `bin_radius`).
        rmax : float or None, optional
            Maximum radius r to include; grid points farther than this are
            dropped, and it sets the outer edge of the last radial
            bin. Defaults to the largest radius present in the (y, z)
            grid.
        """
        if y.ndim != 1 or z.ndim != 1:
            raise ValueError("y and z must be 1D coordinate arrays")

        self.jx = y.size
        self.kx = z.size
        self.dr = float(dr)

        if self.dr <= 0:
            raise ValueError("dr must be a positive value")

        # shape = (jx, kx)
        grid_radius = np.hypot(
            y[:, None] - yc,
            z[None, :] - zc,
        )

        if rmax is None:
            rmax = grid_radius.max()

        self.rmax = float(rmax)

        if self.rmax <= 0:
            raise ValueError("rmax must be a positive value")

        self.nbin = int(np.floor(self.rmax / self.dr)) + 1

        # Grid points to include in the averaging
        mask = grid_radius <= self.rmax

        # Flat-array positions of grid points within rmax
        self.valid_grid_index = np.flatnonzero(mask.ravel())

        # Radial bin index for each valid grid point
        self.valid_radial_bin = np.floor(
            grid_radius.ravel()[self.valid_grid_index] / self.dr
        ).astype(np.int32)

        self.bin_radius = (
            np.arange(self.nbin, dtype=np.float64) + 0.5
        ) * self.dr

        # Grid point count per bin, shared across all x-planes when there are no NaNs
        self.count = np.bincount(
            self.valid_radial_bin,
            minlength=self.nbin,
        )

        # Cache for compute()'s per-chunk_size combined_bin array (see
        # compute()), keyed by chunk_size since that is the only compute()
        # input it depends on.
        self._combined_bin_cache_chunk_size: int | None = None
        self._combined_bin_cache: np.ndarray | None = None

    def _bin_average(
        self,
        combined_bin_valid: np.ndarray,
        weights: np.ndarray,
        counts: np.ndarray,
        chunk_len: int,
        out: np.ndarray | None = None,
    ) -> np.ndarray:
        """
        Sums `weights` per (plane, bin) via `combined_bin_valid` (see
        compute()) and divides by `counts`, writing into `out` if given
        (else a fresh NaN-filled buffer). Bins with counts == 0 keep
        whatever `out` already held there (NaN, by this class's
        NaN-default convention), since np.divide's `where` skips them.
        """
        if out is None:
            out = np.full((chunk_len, self.nbin), np.nan, dtype=np.float64)

        sums = np.bincount(
            combined_bin_valid,
            weights=weights,
            minlength=chunk_len * self.nbin,
        ).reshape(chunk_len, self.nbin)

        np.divide(sums, counts, out=out, where=counts > 0)
        return out

    def compute(
        self,
        qq: np.ndarray,
        *,
        chunk_size: int = 32,
        stats: Sequence[str] = ("mean",),
        qq2: np.ndarray | None = None,
        condition_field: np.ndarray | None = None,
        condition: Callable[[np.ndarray], np.ndarray] | None = None,
    ) -> dict[str, np.ndarray]:
        """
        Parameters
        ----------
        qq : ndarray, shape (jx, kx) or (ix, jx, kx)
            If 2D, treated as a single x-plane.
        chunk_size : int
            Number of x-planes processed at once. Tune based on
            available memory.
        stats : sequence of {"mean", "rms", "cov", "corr"}
            Which statistics to compute. "rms" is the population standard
            deviation (divisor N, not N-1) of qq's deviation from its own
            bin mean. "cov" and "corr" are the population covariance and
            Pearson correlation coefficient (also divisor N) between qq
            and `qq2`, and require `qq2` to be given. Requesting multiple
            stats in one call is cheaper than separate calls, since the
            gather and per-bin sums/counts are shared between them.
        qq2 : ndarray, shape matching qq, or None
            Second field for "cov"/"corr". A grid point contributes to
            "cov"/"corr" only where both qq and qq2 are finite (and pass
            `condition`, if given) -- a narrower set than "mean"/"rms"
            use for qq alone, so requesting e.g. ["mean", "cov"] together
            does not change what "mean" would have been on its own.
        condition_field : ndarray, shape matching qq, or None
            Field to evaluate `condition` against for conditional
            averaging (e.g. average qq only where vertical velocity is
            positive). Gathered and evaluated one chunk at a time, same
            as qq, so no full-size mask array is ever materialized.
            Defaults to qq itself when `condition` is given. Ignored if
            `condition` is None.
        condition : callable or None
            Called once per chunk as `condition(condition_field_chunk)`,
            returning a boolean array of the same shape; grid points
            where it is False are excluded from that call's stats, same
            as NaN grid points are. NaN in `condition_field` is always
            excluded regardless of what `condition` returns for it.

        Returns
        -------
        result : dict[str, ndarray]
            One entry per name in `stats`, each shape (nbin,) if the
            input was 2D, (ix, nbin) if 3D.

        Examples
        --------
        Average qq only over grid points where qq itself is positive
        (`condition_field` defaults to qq)::

            cyl_stats.compute(qq, condition=lambda v: v > 0)

        Average temperature only where vertical velocity vz is positive
        (conditional averaging on a different field)::

            cyl_stats.compute(temperature, condition_field=vz, condition=lambda v: v > 0)

        Covariance and correlation between temperature and vz::

            cyl_stats.compute(temperature, qq2=vz, stats=["cov", "corr"])
        """
        stats_set = set(stats)

        if not stats_set:
            raise ValueError("stats must not be empty")

        unknown = stats_set - _VALID_STATS
        if unknown:
            raise ValueError(
                f"Unknown stats requested: {sorted(unknown)}; "
                f"expected a subset of {sorted(_VALID_STATS)}"
            )

        need_mean = "mean" in stats_set
        need_rms = "rms" in stats_set
        need_cov = "cov" in stats_set
        need_corr = "corr" in stats_set
        need_pair = need_cov or need_corr

        if need_pair and qq2 is None:
            raise ValueError(
                f"stats {sorted(stats_set & _PAIR_STATS)} require `qq2` "
                f"to also be given"
            )

        if chunk_size <= 0:
            raise ValueError("chunk_size must be a positive integer")

        if qq.ndim not in (2, 3):
            raise ValueError(
                "qq.shape must be (jx, kx) or (ix, jx, kx)"
            )

        if condition_field is not None and condition is None:
            raise ValueError(
                "condition_field requires `condition` to also be given"
            )

        if condition is not None and condition_field is None:
            condition_field = qq

        squeeze_output = qq.ndim == 2

        if squeeze_output:
            qq = qq[None, :, :]

        ix, jx, kx = qq.shape

        if (jx, kx) != (self.jx, self.kx):
            raise ValueError(
                f"qq's yz shape {(jx, kx)} does not match "
                f"the expected {(self.jx, self.kx)}"
            )

        if condition is not None:
            if condition_field.ndim == 2:
                condition_field = condition_field[None, :, :]

            if condition_field.shape != qq.shape:
                raise ValueError(
                    f"condition_field shape {condition_field.shape} does "
                    f"not match qq's shape {qq.shape}"
                )

            condition_field_flat = condition_field.reshape(ix, jx * kx)

        if need_pair:
            if qq2.ndim == 2:
                qq2 = qq2[None, :, :]

            if qq2.shape != qq.shape:
                raise ValueError(
                    f"qq2 shape {qq2.shape} does not match qq's shape {qq.shape}"
                )

            qq2_flat = qq2.reshape(ix, jx * kx)

        # NaN default: bins with no valid grid points (count == 0) are never
        # written below, since np.divide's `where` skips them. Only the
        # requested output(s) are allocated.
        qq_mean = (
            np.full((ix, self.nbin), np.nan, dtype=np.float64)
            if need_mean else None
        )
        qq_rms = (
            np.full((ix, self.nbin), np.nan, dtype=np.float64)
            if need_rms else None
        )
        qq_cov = (
            np.full((ix, self.nbin), np.nan, dtype=np.float64)
            if need_cov else None
        )
        qq_corr = (
            np.full((ix, self.nbin), np.nan, dtype=np.float64)
            if need_corr else None
        )

        qq_flat = qq.reshape(ix, jx * kx)

        if self._combined_bin_cache_chunk_size != chunk_size:
            # np.bincount sums into one flat set of bins, so to average all
            # chunk_size planes with a single bincount call below (instead of
            # looping per plane), shift each plane's bin indices into its own
            # non-overlapping slice of a combined index space:
            #   plane 0 -> bins [0, nbin), plane 1 -> bins [nbin, 2*nbin), ...
            offset = (
                np.arange(chunk_size, dtype=np.int64)[:, None] * self.nbin
            )

            # valid_radial_bin (per grid point, shape (n_valid,)) broadcasts
            # against offset (per plane, shape (chunk_size, 1)) to
            # (chunk_size, n_valid): row p holds every grid point's bin index
            # shifted into plane p's slice. This depends only on chunk_size
            # and self.valid_radial_bin (fixed at construction), not on qq,
            # so it is cached here and reused across compute() calls that
            # share the same chunk_size instead of being rebuilt every call.
            self._combined_bin_cache = self.valid_radial_bin[None, :] + offset
            self._combined_bin_cache_chunk_size = chunk_size

        for chunk_start in range(0, ix, chunk_size):
            chunk_end = min(chunk_start + chunk_size, ix)
            chunk_len = chunk_end - chunk_start

            # shape = (chunk_len, number of valid yz grid points)
            # .take(axis=1) instead of qq_flat[chunk_start:chunk_end,
            # self.valid_grid_index]: that fancy-indexing form returns an
            # F-contiguous array for this slice+index-array combination,
            # which forces the .ravel() below to make a second full copy.
            # .take() returns a C-contiguous array directly, so .ravel()
            # becomes a free view.
            values = qq_flat[chunk_start:chunk_end].take(
                self.valid_grid_index, axis=1
            )

            # Slice of the cached combined_bin (built above); slicing rows
            # off the front of a C-contiguous array keeps it C-contiguous,
            # so .ravel() here stays a free view, not a copy. The summed
            # bincount result below is reshaped back to (chunk_len, nbin).
            combined_bin = self._combined_bin_cache[:chunk_len].ravel()

            values_flat = values.ravel()
            valid = np.isfinite(values_flat)

            if condition is not None:
                # Gathered/raveled the same way as qq itself (see above),
                # so no full-size condition array or mask is ever built.
                cond_values_flat = (
                    condition_field_flat[chunk_start:chunk_end]
                    .take(self.valid_grid_index, axis=1)
                    .ravel()
                )
                # NaN in condition_field is always excluded, regardless of
                # what `condition` itself returns for it. Kept as its own
                # mask (not folded straight into `valid`) so it can also
                # gate valid_pair below for cov/corr.
                cond_mask = np.isfinite(cond_values_flat) & condition(cond_values_flat)
                valid = valid & cond_mask
            else:
                cond_mask = None

            if condition is None and valid.all():
                # No NaNs and no extra filtering: reuse the grid point
                # counts computed once in __init__ instead of recomputing
                # them per chunk.
                counts = self.count[None, :]
            else:
                counts = np.bincount(
                    combined_bin[valid],
                    minlength=chunk_len * self.nbin,
                ).reshape(chunk_len, self.nbin)

            # rms needs the per-bin mean as an intermediate (it averages
            # squared deviations from it), so it is computed whenever
            # either output is requested, and the sums/counts bincounts
            # above are shared between both rather than redone per output.
            if need_mean or need_rms:
                chunk_mean = self._bin_average(
                    combined_bin[valid],
                    values_flat[valid],
                    counts,
                    chunk_len,
                    out=qq_mean[chunk_start:chunk_end] if need_mean else None,
                )

            if need_rms:
                # Broadcast each valid grid point's own bin mean back out via
                # the same combined_bin flat index space used for binning
                # (chunk_mean.ravel() is a free view: chunk_mean is always
                # C-contiguous, either as a fresh array or a row-slice of one).
                deviation = (
                    values_flat[valid] - chunk_mean.ravel()[combined_bin[valid]]
                )

                mean_sq = self._bin_average(
                    combined_bin[valid],
                    deviation ** 2,
                    counts,
                    chunk_len,
                    out=qq_rms[chunk_start:chunk_end],
                )
                np.sqrt(mean_sq, out=mean_sq)

            if need_pair:
                values2 = qq2_flat[chunk_start:chunk_end].take(
                    self.valid_grid_index, axis=1
                )
                values2_flat = values2.ravel()

                # A point contributes to cov/corr only where qq *and* qq2
                # are both finite (plus `condition`, if given) -- generally
                # a narrower set than `valid` above, so its own counts must
                # be computed separately rather than reusing `counts`.
                valid_pair = np.isfinite(values_flat) & np.isfinite(values2_flat)
                if cond_mask is not None:
                    valid_pair = valid_pair & cond_mask

                if condition is None and valid_pair.all():
                    counts_pair = self.count[None, :]
                else:
                    counts_pair = np.bincount(
                        combined_bin[valid_pair],
                        minlength=chunk_len * self.nbin,
                    ).reshape(chunk_len, self.nbin)

                # Own bin means over valid_pair (not the possibly-different
                # `chunk_mean` above), so cov/corr are internally consistent
                # with the sample they're actually computed over.
                mean_qq_pair = self._bin_average(
                    combined_bin[valid_pair],
                    values_flat[valid_pair],
                    counts_pair,
                    chunk_len,
                )
                mean_qq2_pair = self._bin_average(
                    combined_bin[valid_pair],
                    values2_flat[valid_pair],
                    counts_pair,
                    chunk_len,
                )

                deviation_qq = (
                    values_flat[valid_pair]
                    - mean_qq_pair.ravel()[combined_bin[valid_pair]]
                )
                deviation_qq2 = (
                    values2_flat[valid_pair]
                    - mean_qq2_pair.ravel()[combined_bin[valid_pair]]
                )

                chunk_cov = self._bin_average(
                    combined_bin[valid_pair],
                    deviation_qq * deviation_qq2,
                    counts_pair,
                    chunk_len,
                    out=qq_cov[chunk_start:chunk_end] if need_cov else None,
                )

                if need_corr:
                    var_qq_pair = self._bin_average(
                        combined_bin[valid_pair],
                        deviation_qq ** 2,
                        counts_pair,
                        chunk_len,
                    )
                    var_qq2_pair = self._bin_average(
                        combined_bin[valid_pair],
                        deviation_qq2 ** 2,
                        counts_pair,
                        chunk_len,
                    )

                    denom = np.sqrt(var_qq_pair * var_qq2_pair)
                    chunk_corr = qq_corr[chunk_start:chunk_end]
                    np.divide(chunk_cov, denom, out=chunk_corr, where=denom > 0)

        if squeeze_output:
            if need_mean:
                qq_mean = qq_mean[0]
            if need_rms:
                qq_rms = qq_rms[0]
            if need_cov:
                qq_cov = qq_cov[0]
            if need_corr:
                qq_corr = qq_corr[0]

        result: dict[str, np.ndarray] = {}
        if need_mean:
            result["mean"] = qq_mean
        if need_rms:
            result["rms"] = qq_rms
        if need_cov:
            result["cov"] = qq_cov
        if need_corr:
            result["corr"] = qq_corr

        return result
