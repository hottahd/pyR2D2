import os
import zipfile
from pathlib import Path

import numpy as np

import pyR2D2


class Data:
    """
    Class for managing R2D2 data

    Attributes
    ----------
    p : pyR2D2.Parameters
        Run parameters and grids (``x, y, z, nx, ...``); also reachable directly as
        ``d.x``, ``d.nx``, ... (attribute lookup falls through to ``d.p``).
    qf : pyR2D2.FullData
        Full 3D data (``remap/qq``).
        ``d.qf.read(n)`` -> ``d.qf.ro, vx, vy, vz, bx, by, bz, se, pr, te, op``
    qx : pyR2D2.XSelect
        2D plane at a given x.
        ``d.qx.read(xs, n)`` -> ``d.qx.ro, vx, vy, vz, bx, by, bz, se, ...``
    qz : pyR2D2.ZSelect
        2D plane at a given z.
        ``d.qz.read(zs, n)`` -> ``d.qz.ro, vx, ...``
    qm : pyR2D2.MPIRegion
        3D data of one remap MPI region.
        ``d.qm.read(ixrt, n)`` -> ``d.qm.ro, ...``
    qr : pyR2D2.RestrictedData
        3D data inside a sub-volume.
        ``d.qr.read(n, keys, x0=..., x1=..., ...)`` -> ``d.qr.ro, ...``
    qt : pyR2D2.OpticalDepth
        2D data on optical-depth surfaces (tau = 1, 0.1, 0.01).
        ``d.qt.read(n)`` -> ``d.qt.rt, he, fr, ro, se, vx, ...``
        (``rt`` at tau=1, ``rt01`` at 0.1, ``rt001`` at 0.01; multigroup ``rt1, rt2, ...``)
    vc : pyR2D2.OnTheFly
        On-the-fly means, RMS, correlations and xy/xz cuts (``remap/vl``).
        ``d.vc.read(n)`` -> ``d.vc.rom, sem, serms, se_xz, bx_xy, ...``;
        ``d.vc.read_et(n)`` -> ``d.vc.etm`` (R2D2plus only)
    qs : pyR2D2.Slice
        2D slice output (``slice/``).
        ``d.qs.read(n_slice, direc, n)`` -> ``d.qs.ro, vx, ..., pr, te``
    q2 : pyR2D2.TwoDimension
        Full data of a 2D run (dict access).
        ``d.q2.read(n)`` -> ``d.q2['ro'], d.q2['tu'], ...``
    ms : pyR2D2.ModelS
        Model S based stratification.  ``d.ms.read()``
    qp, qa : pyR2D2.Previous, pyR2D2.After
        Checkpoint before / after step n.
        ``d.qp.read(n, n_prev)``, ``d.qa.read(n, n_aftr)``
    time : float
        Time at a selected step. See :meth:`pyR2D2.Data.time_read`
    qc : numpy.ndarray, float
        3D checkpoint data. See :meth:`pyR2D2.Data.qc_read`
    sync : pyR2D2.Sync
    eos : pyR2D2.cpp_util.EOS
    yinyang : pyR2D2.cpp_util.YinYang

    ``help(d.qt)`` (or ``d.qt?`` in IPython) shows how to read and which variables appear.

    Output format
    -------------
    The format is detected automatically: when ``<datadir>/param/format.toml`` says
    ``format = "compressed"`` (R2D2plus DEC-576) the readers above read the new files
    with the same arguments, attribute names and array shapes, and compute ``pr``,
    ``te``, ``op`` from ``ro``, ``se`` with ``param/eos_table.tbl`` (not stored in that
    format). ``d.p.output_format`` is ``"legacy"`` or ``"compressed"``.
    See :py:mod:`pyR2D2.data_io.compressed`.
    """

    def __init__(self, datadir, verbose=False, self_old=None):
        """
        Initialize pyR2D2.Data
        """
        self.datadir = Path(datadir)
        self.p = pyR2D2.Parameters(self)
        self.qx = pyR2D2.XSelect(self)
        self.qz = pyR2D2.ZSelect(self)
        self.qm = pyR2D2.MPIRegion(self)
        self.qf = pyR2D2.FullData(self)
        self.qr = pyR2D2.RestrictedData(self)
        self.qt = pyR2D2.OpticalDepth(self)
        self.vc = pyR2D2.OnTheFly(self)
        self.qs = pyR2D2.Slice(self)
        self.q2 = pyR2D2.TwoDimension(self)
        self.ms = pyR2D2.ModelS(self)
        self.time = None
        self.qc = None
        self.qp = pyR2D2.Previous(self)
        self.qa = pyR2D2.After(self)
        self.sync = pyR2D2.Sync(self)

        if verbose:
            self.summary()

        # 台帳 (r2d2plus-ledger) へ未送信の変更があれば 1 行だけ知らせる。
        # R2D2plus のラン (run_summary.toml がある) だけ。失敗しても黙る。
        try:
            from .ledger import unsent_hint

            hint = unsent_hint(self.datadir)
            if hint:
                print(hint)
        except Exception:
            pass

        eosdir = self.datadir.parent / "input_data"
        # 新出力形式 (R2D2plus DEC-576) は param/eos_table.tbl を Parameters が読んでいる
        compressed_with_table = self.p.output_format == "compressed" and hasattr(
            self.p, "log_ro_e"
        )
        if (eosdir / "eos_table_sero.npz").exists() or compressed_with_table:
            self.eos = pyR2D2.cpp_util.EOS(
                self.log_ro_e.astype(np.float32),
                self.se_e.astype(np.float32),
                self.log_pr_e.astype(np.float32),
                self.log_en_e.astype(np.float32),
                self.log_te_e.astype(np.float32),
                self.log_op_e.astype(np.float32),
            )

        if self.p.geometry == "YinYang":
            self.yinyang = pyR2D2.cpp_util.YinYang(
                self.p.yg_yy,
                self.p.zg_yy,
                self.p.y,
                self.p.z,
            )

    def __getattr__(self, name):
        """
        When an attribute is not found in pyR2D2.Data, it is searched in pyR2D2.Data.p
        """
        if hasattr(self.p, name):
            attr = getattr(self.p, name)
            return attr

        raise AttributeError(
            f"'{type(self).__name__}' object has no attribute '{name}'"
        )

    def zip_time(self, remove_original=False):
        """
        Archive datadir/time into datadir/time.zip.

        Existing entries in time.zip are skipped.
        Newly added files under time/ are appended.

        Parameters
        ----------
        remove_original : bool, optional
            If True, remove the original files after archiving. Default is False.
        """
        src = self.datadir / "time"
        dst = self.datadir / "time.zip"

        if not src.exists():
            raise FileNotFoundError(f"{src} does not exist")

        # Files are only deleted after the zip archive below is fully closed
        # (its central directory flushed to disk), so a crash or kill
        # partway through archiving cannot delete originals that were never
        # actually persisted in a readable zip.
        pending_removal = []

        # zipfile.ZIP_STORED is used to avoid compression, which can be slow for many small files.
        # allowZip64=True is used to support large zip files.
        with zipfile.ZipFile(
            dst, "a", compression=zipfile.ZIP_STORED, allowZip64=True
        ) as zf:
            # set is used to avoid duplicate entries in the zip file
            existing = set(zf.namelist())

            for p in src.rglob("*"):
                # Skip directories and non-files
                if not p.is_file():
                    continue

                # Create a relative path for the archive name
                arcname = str(p.relative_to(self.datadir))

                # Skip if the file already exists in the zip archive
                if arcname not in existing:
                    zf.write(p, arcname)
                    existing.add(arcname)

                if remove_original:
                    pending_removal.append(p)

        for p in pending_removal:
            os.remove(p)

    def time_read(self, n, tau=False, verbose=True, use_zip=False):
        """
        Reads time at a selected time step
        The data is stored in self.t

        Parameters
        ----------
        n : int
            selected time step for data
        tau : bool
            if True time for optical depth (high cadence)
        verbose: bool
            if True, print a message indicating that the time variable has been stored in self.time
        use_zip: bool
            if True, read time from time.zip instead of the original files
            if False, read from the original file if it exists; otherwise, fall back to time.zip.

        Returns
        -------
        time : float
            time at a selected time step
        """

        subdir = "tau" if tau else "mhd"
        filename = f"t.dac.{n:08d}"

        filepath = self.datadir / "time" / subdir / filename

        if not filepath.exists() or use_zip:
            zippath = self.datadir / "time.zip"
            arcname = f"time/{subdir}/{filename}"

            with zipfile.ZipFile(zippath, "r") as zf:
                data = zf.read(arcname)
                self.time = np.frombuffer(
                    data,
                    dtype=self.endian + "d",
                    count=1,
                )[0]
        else:
            with open(filepath, "rb") as f:
                self.time = np.fromfile(f, self.endian + "d", 1)[0]

        if verbose:
            print("### time is stored in self.time ###")

        return self.time

    def qc_read(self, n, end_step=False):
        """
        Reads 3D full data for checkpoint
        The data is stored in self.qc dictionary

        Parameters
        ----------
        n : int
            A selected time step for data
        end_step : bool
            If true, checkpoint of end step is read.
        """

        step = f"{n:08d}"
        if end_step:
            if np.mod(self.nd, 2) == 0:
                step = "e"
            if np.mod(self.nd, 2) == 1:
                step = "o"

        with open(self.datadir / "qq" / f"qq.dac.{step}", "rb") as f:
            self.qc = np.fromfile(
                f, self.endian + "d", self.mtype * self.ixg * self.jxg * self.kxg
            ).reshape((self.ixg, self.jxg, self.kxg, self.mtype), order="F")
