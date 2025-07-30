import os

import numpy as np


class _BaseReader:
    """
    Base class for data readers
    """

    def __getattr__(self, name):
        """
        When an attribute is not found in pyR2D2._BaseReader, it is searched in pyR2D2._BaseReader.data
        """
        if hasattr(self.data, name):
            attr = getattr(self.data, name)
            if not callable(attr):
                return attr


class _BaseRemapReader(_BaseReader):
    """
    Base class for remap data readers
    """

    def __init__(self, data):
        self.data = data

        for key in self.data.remap_kind + self.data.remap_kind_add:
            self.__dict__[key] = None

        self._add_docstring()

    def _allocate_remap_qq(self, ijk: tuple):
        """
        Paramters
        ---------
        qq_type : str
            type of data.

        ijk : tuple
            size of data
        """
        memflag = True
        if self.ro is not None:
            memflag = not self.ro.shape == ijk
        if self.ro is None or memflag:
            for key in self.remap_kind + self.remap_kind_add:
                self.__dict__[key] = np.zeros(ijk)

    def _dtype_remap_qq(self, np0):
        """
        returns dtype of remap/qq/

        Parameters
        ----------
        np0 : int
            MPI process number

        Returns
        -------
            dtype : np.dtype
                data type of remap/qq/
        """
        dtype = np.dtype(
            [
                (
                    "qq",
                    self.endian
                    + str(self.mtype * self.iixl[np0] * self.jjxl[np0] * self.kx)
                    + "f",
                ),
                (
                    "pr",
                    self.endian + str(self.iixl[np0] * self.jjxl[np0] * self.kx) + "f",
                ),
                (
                    "te",
                    self.endian + str(self.iixl[np0] * self.jjxl[np0] * self.kx) + "f",
                ),
                (
                    "op",
                    self.endian + str(self.iixl[np0] * self.jjxl[np0] * self.kx) + "f",
                ),
            ]
        )

        return dtype

    def _get_filepath_remap_qq(self, n, np0):
        """
        get file path of remap/qq/

        Parameters
        ----------
        n : int
            time step
        np0 : int
            MPI process number

        Returns
        -------
        filepath : str
            file path of remap/qq/
        """
        if os.path.isdir(self.datadir + "remap/qq/00000/"):
            cnou = "{0:05d}".format(np0 // 1000)
            cno = "{0:08d}".format(np0)
            filepath = (
                self.datadir
                + "remap/qq/"
                + cnou
                + "/"
                + cno
                + "/qq.dac."
                + "{0:08d}".format(n)
                + "."
                + "{0:08d}".format(np0)
            )
        else:
            # directoryを分けない古いバージョン対応
            filepath = (
                self.datadir
                + "remap/qq/qq.dac."
                + "{0:08d}".format(n)
                + "."
                + "{0:08d}".format(np0)
            )
        return filepath

    def _add_docstring(self):
        docstring = "\n"
        docstring += "    Attributes\n"
        docstring += "    ----------\n"
        docstring += "    ro: numpy.ndarray, float\n"
        docstring += "        density\n"
        docstring += "    vx: numpy.ndarray, float\n"
        docstring += "        x velocity\n"
        docstring += "    vy: numpy.ndarray, float\n"
        docstring += "        y velocity\n"
        docstring += "    vz: numpy.ndarray, float\n"
        docstring += "        z velocity\n"
        docstring += "    bx: numpy.ndarray, float\n"
        docstring += "        x magnetic field\n"
        docstring += "    by: numpy.ndarray, float\n"
        docstring += "        y magnetic field\n"
        docstring += "    bz: numpy.ndarray, float\n"
        docstring += "        z magnetic field\n"
        docstring += "    se: numpy.ndarray, float\n"
        docstring += "        specific entropy\n"
        docstring += "    ph: numpy.ndarray, float\n"
        docstring += "        div B cleaning variable\n"
        docstring += "    pr: numpy.ndarray, float\n"
        docstring += "        pressure\n"
        docstring += "    te: numpy.ndarray, float\n"
        docstring += "        temperature\n"
        docstring += "    op: numpy.ndarray, float\n"
        docstring += "        opacity\n"

        self.__class__.__doc__ = self.__class__.__doc__ + docstring


class XSelect(_BaseRemapReader):
    """
    Class for 2D selected data at a certain x

    Important
    ---------
    pyR2D2.Data class can access this class as :code:`pyR2D2.Data.qx`

    """

    def read(self, xs: float, n: int):
        """
        Reads 2D selected data at a certain x

        Parameters
        ----------
        xs : float
            A selected height for data
        n : int
            A selected time step for data
        """

        i0 = np.argmin(np.abs(self.x - xs))
        ir0 = self.i2ir[i0]

        self._allocate_remap_qq(ijk=[self.jx, self.kx])

        for jr0 in range(1, self.jxr + 1):
            np0 = self.np_ijr[ir0 - 1, jr0 - 1]

            if jr0 == self.jr[np0]:
                dtype = self._dtype_remap_qq(np0)
                filepath = self._get_filepath_remap_qq(n, np0)

                qqq_mem = np.memmap(filepath, dtype=dtype, mode="r", shape=(1,))
                qqq = qqq_mem[0]

                for key, m in zip(self.remap_kind, range(self.mtype)):
                    self.__dict__[key][self.jss[np0] : self.jee[np0] + 1, :] = qqq[
                        "qq"
                    ].reshape(
                        (self.iixl[np0], self.jjxl[np0], self.kx, self.mtype), order="F"
                    )[
                        i0 - self.iss[np0], :, :, m
                    ]

                for key in self.remap_kind_add:
                    self.__dict__[key][self.jss[np0] : self.jee[np0] + 1, :] = qqq[
                        key
                    ].reshape((self.iixl[np0], self.jjxl[np0], self.kx), order="F")[
                        i0 - self.iss[np0], :, :
                    ]

                self.info = {}
                self.info["xs"] = self.x[i0]


class ZSelect(_BaseRemapReader):
    """
    Class for 2D selected data at a certain z

    Important
    ---------
    pyR2D2.Data class can access this class as :code:`pyR2D2.Data.qz`

    """

    def read(self, zs: float, n: int):
        """
        Reads 2D slice data at constant z
        The data is stored in self.qz dictionary

        Parameters
        ----------
        zs : float
            A selected z for data
        n : int
            A selected time step for data

        """
        k0 = np.argmin(np.abs(self.z - zs))

        self._allocate_remap_qq(ijk=[self.ix, self.jx])

        for ir0 in range(1, self.ixr + 1):
            for jr0 in range(1, self.jxr + 1):
                np0 = self.np_ijr[ir0 - 1, jr0 - 1]

                dtype = self._dtype_remap_qq(np0)
                filepath = self._get_filepath_remap_qq(n, np0)

                qqq_mem = np.memmap(filepath, dtype=dtype, mode="r", shape=(1,))
                qqq = qqq_mem[0]

                for key, m in zip(self.remap_kind, range(self.mtype)):
                    self.__dict__[key][
                        self.iss[np0] : self.iee[np0] + 1,
                        self.jss[np0] : self.jee[np0] + 1,
                    ] = qqq["qq"].reshape(
                        (self.iixl[np0], self.jjxl[np0], self.kx, self.mtype), order="F"
                    )[
                        :, :, k0, m
                    ]

                for key in self.remap_kind_add:
                    self.__dict__[key][
                        self.iss[np0] : self.iee[np0] + 1,
                        self.jss[np0] : self.jee[np0] + 1,
                    ] = qqq[key].reshape(
                        (self.iixl[np0], self.jjxl[np0], self.kx), order="F"
                    )[
                        :, :, k0
                    ]

            self.info = {}
            self.info["zs"] = self.z[k0]


class MPIRegion(_BaseRemapReader):
    """
    Class for 3D data at a selected MPI process

    Important
    ---------
    pyR2D2.Data class can access this class as :code:`pyR2D2.Data.qm`

    """

    def read(self, ixrt, n):
        """
        Reads 3D data at a selected a MPI process in x-direction.
        corrensponding to ixr-th region in remap coordinate

        Parameters
        ----------
        ixrt : int
            A selected x-region (0<=ixrt<=ixr-1)
            Note
            ----
            ixr is the number of MPI process in x-direction

        n : int
            A selected time step for data
        """
        # corresponding MPI process
        nps = np.where(self.ir - 1 == ixrt)[0]
        # correnponding i range
        self.i_ixrt = np.where(self.i2ir - 1 == ixrt)[0]

        self._allocate_remap_qq(ijk=[len(self.i_ixrt), self.jx, self.kx])

        for np0 in nps:
            if self.iixl[np0] != 0:
                dtype = self._dtype_remap_qq(np0)
                filepath = self._get_filepath_remap_qq(n, np0)

                with open(filepath, "rb") as f:
                    qqq = np.fromfile(f, dtype=dtype, count=1)
                    for key, m in zip(self.remap_kind, range(self.mtype)):
                        self.__dict__[key][
                            :, self.jss[np0] : self.jee[np0] + 1, :
                        ] = qqq["qq"].reshape(
                            (self.iixl[np0], self.jjxl[np0], self.kx, self.mtype),
                            order="F",
                        )[
                            :, :, :, m
                        ]

                    for key in self.remap_kind_add:
                        self.__dict__[key][
                            :, self.jss[np0] : self.jee[np0] + 1, :
                        ] = qqq[key].reshape(
                            (self.iixl[np0], self.jjxl[np0], self.kx), order="F"
                        )[
                            :, :, :
                        ]


class FullData(_BaseRemapReader):
    """
    Class for 3D full data

    Important
    ---------
    pyR2D2.Data class can access this class as :code:`pyR2D2.Data.qf`

    """

    def read(self, n: int, value):
        """
        Reads 3D full data
        The data is stored in self.qq dictionary

        Parameters
        ----------
        n : int
            A selected time step for data
        value : str
            Kind of value. Options are:
                - "ro" : density
                - "vx", "vy", "vz": velocity
                - "bx", "by", "bz": magnetic field
                - "se": entropy
                - "ph": div B cleaning
                - "te": temperature
                - "op": Opacity
        Notes
        -----
        If value is 'all', all values are read
        """

        if type(value) == str:
            if value == "all":
                values_input = self.remap_kind + self.remap_kind_add
            else:
                values_input = [value]
        if type(value) == list:
            values_input = value

        for value in values_input:
            if value not in self.p.remap_kind + self.p.remap_kind_add:
                print("######")
                print("value =", value)
                print(
                    "value should be one of ", self.p.remap_kind + self.p.remap_kind_add
                )
                print("return")
                return

        self._allocate_remap_qq(ijk=[self.ix, self.jx, self.kx])

        for ir0 in range(1, self.ixr + 1):
            for jr0 in range(1, self.jxr + 1):
                np0 = self.np_ijr[ir0 - 1, jr0 - 1]
                if ir0 == self.ir[np0] and jr0 == self.jr[np0]:

                    dtype = self._dtype_remap_qq(np0)
                    filepath = self._get_filepath_remap_qq(n, np0)

                    with open(filepath, "rb") as f:
                        qqq = np.fromfile(f, dtype=dtype, count=1)

                        for value in values_input:
                            if value in self.p.remap_kind:
                                m = self.p.remap_kind.index(value)
                                self.__dict__[value][
                                    self.iss[np0] : self.iee[np0] + 1,
                                    self.jss[np0] : self.jee[np0] + 1,
                                    :,
                                ] = qqq["qq"].reshape(
                                    (
                                        self.iixl[np0],
                                        self.jjxl[np0],
                                        self.kx,
                                        self.mtype,
                                    ),
                                    order="F",
                                )[
                                    :, :, :, m
                                ]
                            else:
                                self.__dict__[value][
                                    self.iss[np0] : self.iee[np0] + 1,
                                    self.jss[np0] : self.jee[np0] + 1,
                                    :,
                                ] = qqq[value].reshape(
                                    (self.iixl[np0], self.jjxl[np0], self.kx), order="F"
                                )[
                                    :, :, :
                                ]


class RestrictedData(_BaseRemapReader):
    """
    Class for 3D restricted-volume data

    Important
    ---------
    pyR2D2.Data class can access this class as :code:`pyR2D2.Data.qr`

    """

    def read(
        self,
        n: int,
        value,
        x0: float,
        x1: float,
        y0: float,
        y1: float,
        z0: float,
        z1: float,
    ):
        """
        Reads 3D restricted-area data
        The data is stored in self.qr dictionary

        Parameters
        ----------
        n : int
            A selected time step for data
        value : str
            See :meth:`R2D2.Read.qq_3d` for options

        x0, y0, z0 : float
            Minimum x, y, z
        x1, y1, z1 : float
            Maximum x, y, z
        """

        i0, i1 = np.argmin(abs(self.x - x0)), np.argmin(abs(self.x - x1))
        j0, j1 = np.argmin(abs(self.y - y0)), np.argmin(abs(self.y - y1))
        k0, k1 = np.argmin(abs(self.z - z0)), np.argmin(abs(self.z - z1))
        ixr = i1 - i0 + 1
        jxr = j1 - j0 + 1
        kxr = k1 - k0 + 1

        if type(value) == str:
            if value == "all":
                values_input = self.remap_kind + self.remap_kind_add
            else:
                values_input = [value]
        if type(value) == list:
            values_input = value

        for value in values_input:
            self.__dict__[value] = np.zeros((ixr, jxr, kxr), dtype=np.float32)

        for ir0 in range(1, self.ixr + 1):
            for jr0 in range(1, self.jxr + 1):
                np0 = self.np_ijr[ir0 - 1, jr0 - 1]

                if not (
                    self.iss[np0] > i1
                    or self.iee[np0] < i0
                    or self.jss[np0] > j1
                    or self.jee[np0] < j0
                ):

                    dtype = self._dtype_remap_qq(np0)
                    filepath = self._get_filepath_remap_qq(n, np0)

                    with open(filepath, "rb") as f:
                        qqq = np.fromfile(f, dtype=dtype, count=1)

                        for value in values_input:
                            isrt_rcv = max([0  ,self.iss[np0]-i0  ])
                            iend_rcv = min([ixr,self.iee[np0]-i0+1])
                            jsrt_rcv = max([0  ,self.jss[np0]-j0  ])
                            jend_rcv = min([jxr,self.jee[np0]-j0+1])
                                
                            isrt_snd = isrt_rcv - (self.iss[np0]-i0)
                            iend_snd = isrt_snd + (iend_rcv - isrt_rcv)
                            jsrt_snd = jsrt_rcv - (self.jss[np0]-j0)
                            jend_snd = jsrt_snd + (jend_rcv - jsrt_rcv)
                            
                            if value in self.p.remap_kind:
                                m = self.p.remap_kind.index(value)
<<<<<<< HEAD
                                isrt_rcv = max([0, self.iss[np0] - i0])
                                iend_rcv = min([ixr, self.iee[np0] - i0 + 1])
                                jsrt_rcv = max([0, self.jss[np0] - j0])
                                jend_rcv = min([jxr, self.jee[np0] - j0 + 1])

                                isrt_snd = isrt_rcv - (self.iss[np0] - i0)
                                iend_snd = isrt_snd + (iend_rcv - isrt_rcv)
                                jsrt_snd = jsrt_rcv - (self.jss[np0] - j0)
                                jend_snd = jsrt_snd + (jend_rcv - jsrt_rcv)
=======
>>>>>>> f2c20fd (RestrictedDataクラスのreadメソッド内での受信および送信インデックス計算を整理)

                                self.__dict__[value][
                                    isrt_rcv:iend_rcv, jsrt_rcv:jend_rcv, :
                                ] = qqq["qq"].reshape(
                                    (
                                        self.iixl[np0],
                                        self.jjxl[np0],
                                        self.kx,
                                        self.mtype,
                                    ),
                                    order="F",
                                )[
                                    isrt_snd:iend_snd, jsrt_snd:jend_snd, k0 : k1 + 1, m
                                ]
                            else:
                                self.__dict__[value][
                                    isrt_rcv:iend_rcv, jsrt_rcv:jend_rcv, :
                                ] = qqq[value].reshape(
                                    (self.iixl[np0], self.jjxl[np0], self.kx), order="F"
                                )[
                                    isrt_snd:iend_snd, jsrt_snd:jend_snd, k0 : k1 + 1
                                ]


class OpticalDepth(_BaseReader):
    """
    Class for 2D data at certain optical depths

    Important
    ---------
    pyR2D2.Data class can access this class as :code:`pyR2D2.Data.qt`

    Attributes
    ----------
    rt : numpy.ndarray, float
        raditive intensity
    ro : numpy.ndarray, float
        density
    se : numpy.ndarray, float
        specific entropy
    pr : numpy.ndarray, float
        pressure
    te : numpy.ndarray, float
        temperature
    vx : numpy.ndarray, float
        x velocity
    vy : numpy.ndarray, float
        y velocity
    vz : numpy.ndarray, float
        z velocity
    bx : numpy.ndarray, float
        x magnetic field
    by : numpy.ndarray, float
        y magnetic field
    bz : numpy.ndarray, float
        z magnetic field
    he : numpy.ndarray, float
        height of the optical depth
    fr : numpy.ndarray, float
        raditveit flux at the optical depth

    Important
    ---------
    Each attibutes have three values at optical depths 1, 0.1, and 0.01
    These correspond to, for example, rt, rt01, and rt001


    """

    def __init__(self, data):
        self.data = data
        self.value_keys = [
            "rt",
            "ro",
            "se",
            "pr",
            "te",
            "vx",
            "vy",
            "vz",
            "bx",
            "by",
            "bz",
            "he",
            "fr",
        ]

        for key in self.value_keys:
            for tau in ["", "01", "001"]:
                self.__dict__[key + tau] = None

    def read(self, n: int):
        """
        Reads 2D data at certain optical depths.
        The data is stored in self.qt dictionary.
        In this version the selected optical depth is 1, 0.1, and 0.01

        Parameters
        ----------
        n : int
            A selected time step for data
        """
        with open(self.datadir + "tau/qq.dac." + "{0:08d}".format(n), "rb") as f:
            qq = np.fromfile(
                f, self.endian + "f", self.m_tu * self.m_in * self.jx * self.kx
            )

        for key, mk in zip(self.value_keys, range(len(self.value_keys))):
            for tau, mt in zip(["", "01", "001"], range(3)):
                self.__dict__[key + tau] = qq.reshape(
                    (self.m_tu, self.m_in, self.jx, self.kx), order="F"
                )[mt, mk, :, :]


class OnTheFly(_BaseReader):
    """
    Class for on-the-fly analysis data
    Mean, RMS, and correlation are done in longitudinal or z directions.


    Important
    ---------
    pyR2D2.Data class can access this class as :code:`pyR2D2.Data.vc`

    """

    def __init__(self, data):
        self.data = data

        m_total = self.m2d_xy + self.m2d_xz + self.m2d_flux
        if self.geometry == "YinYang":
            m_total += self.m2d_spex

        for m in range(m_total):
            self.__dict__[self.cl[m]] = None

        self._generate_docstring()

    def read(self, n):
        """
        Reads on the fly analysis data from fortran.
        The data is stored in self.vc dictionary

        Parameters
        ----------
        n : int
            A selected time step for data
        """

        # read xy plane data
        with open(
            self.datadir + "remap/vl/vl_xy.dac." + "{0:08d}".format(n), "rb"
        ) as f:
            vl = np.fromfile(
                f, self.endian + "f", self.m2d_xy * self.ix * self.jx
            ).reshape((self.ix, self.jx, self.m2d_xy), order="F")

        for m in range(self.m2d_xy):
            self.__dict__[self.cl[m]] = vl[:, :, m]

        # read xz plane data
        with open(
            self.datadir + "remap/vl/vl_xz.dac." + "{0:08d}".format(n), "rb"
        ) as f:
            vl = np.fromfile(
                f, self.endian + "f", self.m2d_xz * self.ix * self.kx
            ).reshape((self.ix, self.kx, self.m2d_xz), order="F")

        for m in range(self.m2d_xz):
            self.__dict__[self.cl[m + self.m2d_xy]] = vl[:, :, m]

        # read flux related value
        with open(
            self.datadir + "remap/vl/vl_flux.dac." + "{0:08d}".format(n), "rb"
        ) as f:
            vl = np.fromfile(
                f, self.endian + "f", self.m2d_flux * (self.ix + 1) * self.jx
            ).reshape((self.ix + 1, self.jx, self.m2d_flux), order="F")

        for m in range(self.m2d_flux):
            self.__dict__[self.cl[m + self.m2d_xy + self.m2d_xz]] = vl[:, :, m]

        # read spectra
        if self.geometry == "YinYang":
            with open(
                self.datadir + "remap/vl/vl_spex.dac." + "{0:08d}".format(n), "rb"
            ) as f:
                vl = np.fromfile(
                    f, self.endian + "f", self.m2d_spex * self.ix * self.kx // 4
                ).reshape((self.ix, self.kx // 4, self.m2d_spex), order="F")

            for m in range(self.m2d_spex):
                self.__dict__[
                    self.cl[m + self.m2d_xy + self.m2d_xz + self.m2d_flux]
                ] = vl[:, :, m]

    def _update_json_template(
        self,
        output_file=os.path.dirname(os.path.abspath(__file__)) + "/" + "OnTheFly.json",
    ):
        """
        Generate JSON template for pyR2D2.Parameters
        """
        import json

        # Extract attributes from the instance
        current_attributes = {k: "" for k in self.__dict__ if not k.startswith("_")}

        try:
            with open(output_file, "r") as f:
                existing_data = json.load(f)
        except FileNotFoundError:
            existing_data = {}

        updated_data = {**current_attributes, **existing_data}

        # Save the template as a JSON file
        with open(output_file, "w") as f:
            json.dump(updated_data, f, indent=4)

    def _generate_docstring(
        self,
        json_file=os.path.dirname(os.path.abspath(__file__)) + "/" + "OnTheFly.json",
        update_json=False,
    ):
        """
        Generate docstring for pyR2D2.Parameters
        """
        import json

        if update_json:
            self._update_json_template(json_file)

        gene_space = ""
        desc_space = " " * 4

        docstring = "\n"
        docstring += gene_space + "Attributes\n"
        docstring += gene_space + "----------\n"
        try:
            with open(json_file, "r") as f:
                saved_attributes = json.load(f)
        except FileNotFoundError:
            saved_attributes = {}
        for key in self.__dict__.keys():
            type_name = "numpy.ndarray, float"

            docstring += gene_space + f"{key} : {type_name}\n"
            if key in saved_attributes:
                docstring += gene_space + desc_space + saved_attributes[key] + "\n"
            else:
                docstring += "\n"

        self.__class__.__doc__ = self.__class__.__doc__ + docstring

    ##############################


class Slice(_BaseReader):
    """
    Class for 2D slice data

    Important
    ---------
    pyR2D2.Data class can access this class as :code:`pyR2D2.Data.qs`

    Attributes
    ----------
    ro : numpy.ndarray, float
        density
    vx : numpy.ndarray, float
        x velocity
    vy : numpy.ndarray, float
        y velocity
    vz : numpy.ndarray, float
        z velocity
    bx : numpy.ndarray, float
        x magnetic field
    by : numpy.ndarray, float
        y magnetic field
    bz : numpy.ndarray, float
        z magnetic field
    se : numpy.ndarray, float
        specific entropy
    ph : numpy.ndarray, float
        div B cleaning variable
    pr : numpy.ndarray, float
        pressure
    te : numpy.ndarray, float
        temperature

    """

    def __init__(self, data):
        self.data = data

        for key in self.data.remap_kind + self.data.remap_kind_add[:-1]:
            self.__dict__[key] = None

    def read(self, n_slice, direc, n):
        """
        Reads 2D data of slice.
        The data is stored in self.ql dictionary

        Parameters
        ----------
        n_slice : int
            index of slice
        direc : str
            slice direction. 'x', 'y', or 'z'
        n : int
            a selected time step for data
        """

        if self.geometry == "YinYang":
            postfixes = ["_yin", "_yan"]
        else:
            postfixes = [""]

        for postfix in postfixes:
            with open(
                self.datadir
                + "slice/qq"
                + direc
                + postfix
                + ".dac."
                + "{0:08d}".format(n)
                + "."
                + "{0:08d}".format(n_slice + 1),
                "rb",
            ) as f:
                if direc == "x":
                    if self.geometry == "YinYang":
                        n1, n2 = (
                            self.jx_yy + 2 * self.margin,
                            self.kx_yy + 2 * self.margin,
                        )
                    else:
                        n1, n2 = self.jx, self.kx
                if direc == "y":
                    n1, n2 = self.ix, self.kx
                if direc == "z":
                    n1, n2 = self.ix, self.jx
                qq = np.fromfile(f, self.endian + "f", (self.mtype + 2) * n1 * n2)

            for key, m in zip(
                self.remap_kind + self.remap_kind_add[:-1], range(self.mtype + 2)
            ):
                self.__dict__[key + postfix] = qq.reshape(
                    (n1, n2, self.mtype + 2), order="F"
                )[:, :, m]
            self.info = {}
            self.info["direc"] = direc
            if direc == "x":
                slice = self.x_slice
            if direc == "y":
                slice = self.y_slice
            if direc == "z":
                slice = self.z_slice
            self.info["slice"] = slice[n_slice]
            self.info["n_slice"] = n_slice


class TwoDimension(_BaseReader):
    """
    Class for 2D data
    """

    def __init__(self, data):
        self.data = data
        self.value_keys = [
            "ro",
            "vx",
            "vy",
            "vz",
            "bx",
            "by",
            "bz",
            "se",
            "pr",
            "te",
            "op",
            "tu",
            "fr",
        ]

        for key in self.value_keys:
            self.__dict__[key] = None

    def read(self, n):
        """
        Reads full data of 2D calculation
        The data is stored in self.q2 dictionary

        Parameters
        ----------
        n : int
            A selected time step for data
        """

        ### Only when memory is not allocated
        ### and the size of array is different
        ### memory is allocated
        memflag = True
        if hasattr(self, "ro"):
            memflag = not self.ro.shape == (self.ix, self.jx)
        if not hasattr(self, "ro") or memflag:
            for key in self.value_keys:
                self.__dict__[key] = np.zeros((self.ix, self.jx))

        dtype = np.dtype(
            [("qq", self.endian + str((self.mtype + 5) * self.ix * self.jx) + "f")]
        )
        with open(self.datadir + "remap/qq/qq.dac." + "{0:08d}".format(n), "rb") as f:
            qq = np.fromfile(f, dtype=dtype, count=1)

        for key, m in zip(value_keys, range(self.mtype)):
            self.__dict__[key] = qq["qq"].reshape(
                (self.mtype + 5, self.ix, self.jx), order="F"
            )[m, :, :]


class ModelS(_BaseReader):
    """
    Class for Model S based stratification data
    """

    def __init__(self, data):
        self.data = data

    def read(self):
        """
        Reads Model S based stratification data
        """

        with open(self.datadir + "../input_data/params.txt", "r") as f:
            self.__dict__["ix"] = int(f.read())

        endian = ">"
        dtype = np.dtype(
            [
                ("head", endian + "i"),
                ("x", endian + str(self.ix) + "d"),
                ("pr0", endian + str(self.ix) + "d"),
                ("ro0", endian + str(self.ix) + "d"),
                ("te0", endian + str(self.ix) + "d"),
                ("se0", endian + str(self.ix) + "d"),
                ("en0", endian + str(self.ix) + "d"),
                ("op0", endian + str(self.ix) + "d"),
                ("tu0", endian + str(self.ix) + "d"),
                ("dsedr0", endian + str(self.ix) + "d"),
                ("dtedr0", endian + str(self.ix) + "d"),
                ("dprdro", endian + str(self.ix) + "d"),
                ("dprdse", endian + str(self.ix) + "d"),
                ("dtedro", endian + str(self.ix) + "d"),
                ("dtedse", endian + str(self.ix) + "d"),
                ("dendro", endian + str(self.ix) + "d"),
                ("dendse", endian + str(self.ix) + "d"),
                ("gx", endian + str(self.ix) + "d"),
                ("kp", endian + str(self.ix) + "d"),
                ("cp", endian + str(self.ix) + "d"),
                ("tail", endian + "i"),
            ]
        )
        with open(self.datadir + "../input_data/value_cart.dac", "rb") as f:
            qq = np.fromfile(f, dtype=dtype, count=1)

        for key in qq.dtype.names:
            if qq[key].size == self.ix:
<<<<<<< HEAD
                self.__dict__[key] = qq[key].reshape((self.ix), order="F")
=======
                self.__dict__[key] = qq[key].reshape((self.ix),order='F')
>>>>>>> f2c20fd (RestrictedDataクラスのreadメソッド内での受信および送信インデックス計算を整理)
