import numpy as np

import pyR2D2


class Data:
    """
    Class for managing R2D2 data

    Attributes
    ----------
    p : pyR2D2.Parameters
        Instance of pyR2D2.Parameters
    qx : pyR2D2.XSelect
        Instance of pyR2D2.XSelect
    qz : pyR2D2.ZSelect
        Instance of pyR2D2.ZSelect
    qm : pyR2D2.MPIRegion
        Instance of pyR2D2.MPIRegion
    qf : pyR2D2.FullData
        Instance of pyR2D2.FullData
    qr : pyR2D2.RestrictedData
        Instance of pyR2D2.RestrictedData
    qt : pyR2D2.OpticalDepth
        Instance of pyR2D2.OpticalDepth
    vc : pyR2D2.OnTheFly
        Instance of pyR2D2.OnTheFly
    qs : pyR2D2.Slice
        Instance of pyR2D2.Slice
    ms : pyR2D2.ModelS
        Instance of pyR2D2.ModelS
    time : float
        Time at a selected time step. See :meth:`pyR2D2.Data.time_read`
    qc : numpy.ndarray, float
        3D full data for checkpoint. See :meth:`pyR2D2.Data.qc_read`

    sync : pyR2D2.Sync
        Instance of pyR2D2.Sync
    eos : pyR2D2.EOS
        Instance of pyR2D2.EOS

    """

    def __init__(self, datadir, verbose=False, self_old=None):
        """
        Initialize pyR2D2.Data
        """
        self.datadir = datadir
        self.p = pyR2D2.Parameters(self)
        self.qx = pyR2D2.XSelect(self)
        self.qz = pyR2D2.ZSelect(self)
        self.qm = pyR2D2.MPIRegion(self)
        self.qf = pyR2D2.FullData(self)
        self.qr = pyR2D2.RestrictedData(self)
        self.qt = pyR2D2.OpticalDepth(self)
        self.vc = pyR2D2.OnTheFly(self)
        self.qs = pyR2D2.Slice(self)
        self.ms = pyR2D2.ModelS(self)
        self.time = None
        self.qc = None
        self.sync = pyR2D2.Sync(self)
        self.eos = pyR2D2.cpp_util.EOS(
            self.log_ro_e,
            self.se_e,
            self.log_pr_e,
            self.log_en_e,
            self.log_te_e,
            self.log_op_e,
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

    def time_read(self, n, tau=False, verbose=True):
        """
        Reads time at a selected time step
        The data is stored in self.t

        Parameters
        ----------
        n : int
            selected time step for data
        tau : bool
            if True time for optical depth (high cadence)

        Returns
        -------
        time : float
            time at a selected time step
        """

        if tau:
            with open(
                self.datadir + "time/tau/t.dac." + "{0:08d}".format(n), "rb"
            ) as f:
                self.time = np.fromfile(f, self.endian + "d", 1).reshape(
                    (1), order="F"
                )[0]
        else:
            with open(
                self.datadir + "time/mhd/t.dac." + "{0:08d}".format(n), "rb"
            ) as f:
                self.time = np.fromfile(f, self.endian + "d", 1).reshape(
                    (1), order="F"
                )[0]

        if verbose:
            print("### variales are stored in self.time ###")

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

        step = "{0:08d}".format(n)
        if end_step:
            if np.mod(self.nd, 2) == 0:
                step = "e"
            if np.mod(self.nd, 2) == 1:
                step = "o"

        with open(self.datadir + "qq/qq.dac." + step, "rb") as f:
            self.qc = np.fromfile(
                f, self.endian + "d", self.mtype * self.ixg * self.jxg * self.kxg
            ).reshape((self.ixg, self.jxg, self.kxg, self.mtype), order="F")
