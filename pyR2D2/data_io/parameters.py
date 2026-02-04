import os
import pydoc

import numpy as np

import pyR2D2


class Parameters:
    """
    Class for managing R2D2 basic parameters
    """

    def __init__(self, data):
        """
        Initialize pyR2D2.Parameters

        Parameters
        ----------
        data : pyR2D2.Data
            Instance of pyR2D2.Data
        """
        from scipy.io import FortranFile

        self.data = data
        self.datadir = data.datadir

        # Time control parameters
        with open(self.datadir + "param/nd.dac", "r") as f:
            nn = f.read().split()
            self.nd = int(nn[0])  # current maximum time step
            self.nd_tau = int(nn[1])  # current maximum time step for tau

        # data/time/tau のファイルの数を数え、nd_tauと比較して大きい方をnd_tauとする
        if os.path.isdir(self.datadir + "time/tau"):
            self.nd_tau = max(
                len(os.listdir(self.datadir + "time/tau")) - 1, self.nd_tau
            )

        # Read basic parameters
        with open(self.datadir + "param/params.dac", "r") as f:
            for line in f:
                parts = line.split()
                key = parts[1]
                value = parts[0]
                value_type = parts[2]

                if value_type == "i":  # integer
                    self.__dict__[key] = int(value)
                elif value_type == "d":  # double
                    self.__dict__[key] = float(value)
                elif value_type == "c":  # character
                    self.__dict__[key] = value
                elif value_type == "l":
                    if value == "F":
                        self.__dict__[key] = False
                    else:
                        self.__dict__[key] = True

        # self.rstarのない場合はself.rstar = pyR2D2.constant.RSUN
        if "rstar" not in self.__dict__:
            self.rstar = pyR2D2.constant.RSUN
        if "lstar" not in self.__dict__:
            self.lstar = pyR2D2.constant.LSUN

        # remapで出力している変数の種類
        self.remap_kind = ["ro", "vx", "vy", "vz", "bx", "by", "bz", "se", "ph"]
        self.remap_kind_add = ["pr", "te", "op"]

        # transform endiant for python
        if self.swap == 0:  # little
            self.endian = "<"
        else:
            self.endian = ">"

        # MPI information
        f = FortranFile(self.datadir + "param/xyz.dac", "r")
        shape = (self.ix0 * self.jx0 * self.kx0, 3)
        self.xyz = f.read_reals(dtype=np.int32).reshape(shape, order="F")

        ## number of grid in each direction
        self.ix = self.ix0 * self.nx
        self.jx = self.jx0 * self.ny
        self.kx = self.kx0 * self.nz

        self.ixg = self.ix + 2 * self.margin * (self.xdcheck - 1)
        self.jxg = self.jx + 2 * self.margin * (self.ydcheck - 1)
        self.kxg = self.kx + 2 * self.margin * (self.zdcheck - 1)

        # read background variables
        dtype = np.dtype(
            [
                (
                    "head",
                    self.endian + "i",
                ),
                ("x", self.endian + "d", (self.ixg,)),
                ("y", self.endian + "d", (self.jxg,)),
                ("z", self.endian + "d", (self.kxg,)),
                ("pr0", self.endian + "d", (self.ixg,)),
                ("te0", self.endian + "d", (self.ixg,)),
                ("ro0", self.endian + "d", (self.ixg,)),
                ("se0", self.endian + "d", (self.ixg,)),
                ("en0", self.endian + "d", (self.ixg,)),
                ("op0", self.endian + "d", (self.ixg,)),
                ("tu0", self.endian + "d", (self.ixg,)),
                ("dsedr0", self.endian + "d", (self.ixg,)),
                ("dtedr0", self.endian + "d", (self.ixg,)),
                ("dprdro", self.endian + "d", (self.ixg,)),
                ("dprdse", self.endian + "d", (self.ixg,)),
                ("dtedro", self.endian + "d", (self.ixg,)),
                ("dtedse", self.endian + "d", (self.ixg,)),
                ("dendro", self.endian + "d", (self.ixg,)),
                ("dendse", self.endian + "d", (self.ixg,)),
                ("gx", self.endian + "d", (self.ixg,)),
                ("cp", self.endian + "d", (self.ixg,)),
                ("fa", self.endian + "d", (self.ixg,)),
                ("sa", self.endian + "d", (self.ixg,)),
                ("xi", self.endian + "d", (self.ixg,)),
                ("tail", self.endian + "i"),
            ]
        )
        with open(self.datadir + "param/back.dac", "rb") as f:
            back = np.fromfile(f, dtype=dtype, count=1)

        # marginも含んだ座標
        self.xg = back["x"].reshape((self.ixg), order="F")
        self.yg = back["y"].reshape((self.jxg), order="F")
        self.zg = back["z"].reshape((self.kxg), order="F")

        # fluxはi+1/2で定義するので、そのための座標を定義
        self.x_flux = np.zeros(self.ix + 1)
        for i in range(0, self.ix + 1):
            self.x_flux[i] = 0.5 * (
                self.xg[i + self.margin] + self.xg[i + self.margin - 1]
            )

        self.ro0g = back["ro0"].reshape((self.ixg), order="F")
        self.se0g = back["se0"].reshape((self.ixg), order="F")

        for key in back.dtype.names:
            if back[key].size == self.ixg:
                self.__dict__[key] = back[key].reshape((self.ixg), order="F")[
                    self.margin : self.ixg - self.margin
                ]
            elif back[key].size == self.jxg:
                self.__dict__[key] = back[key].reshape((self.jxg), order="F")[
                    self.margin : self.jxg - self.margin
                ]
            elif back[key].size == self.kxg:
                self.__dict__[key] = back[key].reshape((self.kxg), order="F")[
                    self.margin : self.kxg - self.margin
                ]

        self.xr = self.x / self.rstar

        if self.geometry == "YinYang":
            # original geometry in Yin-Yang grid
            self.jx_yy = self.jx
            self.kx_yy = self.kx

            self.jxg_yy = self.jx + 2 * self.margin
            self.kxg_yy = self.kx + 2 * self.margin

            self.y_yy = self.y - 0.5 * (self.y[1] - self.y[0])
            self.z_yy = self.z - 0.5 * (self.z[1] - self.z[0])

            self.yg_yy = self.yg
            self.zg_yy = self.zg

            # converted spherical geometry
            self.jx = self.jx * 2
            self.kx = self.jx * 2

            dy = np.pi / self.jx
            dz = 2 * np.pi / self.kx

            y = np.zeros(self.jx)
            z = np.zeros(self.kx)

            y[0] = 0.5 * dy
            z[0] = 0.5 * dz - np.pi

            for j in range(1, self.jx):
                y[j] = y[j - 1] + dy

            for k in range(1, self.kx):
                z[k] = z[k - 1] + dz

            self.y = y
            self.z = z

        if self.zdcheck == 2:
            dimension = "3d"
        else:
            dimension = "2d"

        ##############################
        # read value information
        if dimension == "3d":
            with open(self.datadir + "remap/vl/c.dac", "r") as f:
                value = f.read().split()
                self.m2d_xy = int(value[0])
                self.m2d_xz = int(value[1])
                self.m2d_flux = int(value[2])
                if self.geometry == "YinYang":
                    self.m2d_spex = int(value[3])
                del value[0:3]
                if self.geometry == "YinYang":
                    del value[0]
                self.cl = list(map(str.strip, value))  ## strip space from character

            ##############################
            # read mpi information for remap
            dtyp = np.dtype(
                [
                    ("iss", self.endian + str(self.npe) + "i4"),
                    ("iee", self.endian + str(self.npe) + "i4"),
                    ("jss", self.endian + str(self.npe) + "i4"),
                    ("jee", self.endian + str(self.npe) + "i4"),
                    ("iixl", self.endian + str(self.npe) + "i4"),
                    ("jjxl", self.endian + str(self.npe) + "i4"),
                    ("np_ijr", self.endian + str(self.ixr * self.jxr) + "i4"),
                    ("ir", self.endian + str(self.npe) + "i4"),
                    ("jr", self.endian + str(self.npe) + "i4"),
                    ("i2ir", self.endian + str(self.ixg) + "i4"),
                    ("j2jr", self.endian + str(self.jxg) + "i4"),
                ]
            )

            with open(self.datadir + "remap/remap_info.dac", "rb") as f:
                mpi = np.fromfile(f, dtype=dtyp, count=1)

            for key in mpi.dtype.names:
                if key == "np_ijr":
                    self.__dict__[key] = mpi[key].reshape(
                        (self.ixr, self.jxr), order="F"
                    )
                else:
                    self.__dict__[key] = mpi[key].reshape((mpi[key].size), order="F")

            self.i2ir = self.i2ir[self.margin : self.ixg - self.margin]
            self.j2jr = self.j2jr[self.margin : self.jxg - self.margin]

            self.iss = self.iss - 1
            self.iee = self.iee - 1
            self.jss = self.jss - 1
            self.jee = self.jee - 1

            if os.path.isdir(self.datadir + "slice"):
                with open(self.datadir + "slice/params.dac", "r") as f:
                    line = f.readline()
                    while line:
                        self.__dict__[line.split()[1]] = int(line.split()[0])
                        line = f.readline()

                dtype = np.dtype(
                    [
                        ("x_slice", (self.endian + "d", (self.nx_slice,))),
                        ("y_slice", (self.endian + "d", (self.ny_slice,))),
                        ("z_slice", (self.endian + "d", (self.nz_slice,))),
                    ]
                )
                with open(self.datadir + "slice/slice.dac", "rb") as f:
                    slice = np.fromfile(f, dtype=dtype)

                    self.x_slice = slice["x_slice"].reshape(self.nx_slice, order="F")
                    self.y_slice = slice["y_slice"].reshape(self.ny_slice, order="F")
                    self.z_slice = slice["z_slice"].reshape(self.nz_slice, order="F")

        # read equation of state
        eosdir = self.datadir[:-5] + "input_data/"
        if os.path.exists(eosdir + "eos_table_sero.npz"):
            eos_d = np.load(eosdir + "eos_table_sero.npz")
            self.log_ro_e = eos_d["ro"]  # density is defined in logarithmic scale
            self.se_e = eos_d["se"]
            self.ix_e = len(self.log_ro_e)
            self.jx_e = len(self.se_e)

            self.log_pr_e = np.log(eos_d["pr"] + 1.0e-200)
            self.log_en_e = np.log(eos_d["en"] + 1.0e-200)
            self.log_te_e = np.log(eos_d["te"] + 1.0e-200)
            self.log_op_e = np.log(eos_d["op"] + 1.0e-200)
            self.log_dprdro_e = np.log(eos_d["dprdro"] + 1.0e-200)

            self.dlogro_e = self.log_ro_e[1] - self.log_ro_e[0]
            self.dse_e = self.se_e[1] - self.se_e[0]

        # read original data
        if os.path.exists(self.datadir + "cont_log.txt"):
            with open(self.datadir + "cont_log.txt") as f:
                self.origin = f.readlines()[6][-11:-7]
        else:
            self.origin = "N/A"

        self._generate_docstring()

    def yinyang_setup(self):
        """
        YinYangSet function sets up the YinYang geometry for
        YinYang direct plot.
        """
        if self.geometry == "YinYang":
            if not "Z_yy" in self.__dict__:
                self.Z_yy, self.Y_yy = np.meshgrid(self.z_yy, self.y_yy)
                self.Zg_yy, self.Yg_yy = np.meshgrid(self.zg_yy, self.yg_yy)

                # Geometry in Yang grid
                self.Yo_yy = np.arccos(np.sin(self.Y_yy) * np.sin(self.Z_yy))
                self.Zo_yy = np.arcsin(np.cos(self.Y_yy) / np.sin(self.Yo_yy))

                self.Yog_yy = np.arccos(np.sin(self.Yg_yy) * np.sin(self.Zg_yy))
                self.Zog_yy = np.arcsin(np.cos(self.Yg_yy) / np.sin(self.Yog_yy))

                sct = np.sin(self.Y_yy) * np.cos(self.Z_yy)
                sco = -np.sin(self.Yo_yy) * np.cos(self.Zo_yy)

                sctg = np.sin(self.Yg_yy) * np.cos(self.Zg_yy)
                scog = -np.sin(self.Yog_yy) * np.cos(self.Zog_yy)

                tmp = sct * sco < 0
                self.Zo_yy[tmp] = np.sign(self.Zo_yy[tmp]) * np.pi - self.Zo_yy[tmp]

                tmp = sctg * scog < 0
                self.Zog_yy[tmp] = np.sign(self.Zog_yy[tmp]) * np.pi - self.Zog_yy[tmp]

    def summary(self):
        """
        Show pyR2D2.Parameter summary.
        """

        RED = "\033[31m"
        END = "\033[0m"

        print(RED + "### Star information ###" + END)
        if "mstar" in self.__dict__:
            print(
                "mstar = ", "{:.2f}".format(self.mstar / pyR2D2.constant.MSUN) + " msun"
            )
        if "astar" in self.__dict__:
            print("astar = ", "{:.2e}".format(self.astar) + " yr")

        print("")
        if self.geometry == "Cartesian":
            print(RED + "### calculation domain ###" + END)
            print(
                "xmax - rstar = ",
                "{:6.2f}".format((self.xmax - self.rstar) * 1.0e-8),
                "[Mm], xmin - rstar = ",
                "{:6.2f}".format((self.xmin - self.rstar) * 1.0e-8),
                "[Mm]",
            )
            print(
                "ymax         = ",
                "{:6.2f}".format(self.ymax * 1.0e-8),
                "[Mm], ymin         = ",
                "{:6.2f}".format(self.ymin * 1.0e-8),
                "[Mm]",
            )
            print(
                "zmax         = ",
                "{:6.2f}".format(self.zmax * 1.0e-8),
                "[Mm], zmin         = ",
                "{:6.2f}".format(self.zmin * 1.0e-8),
                "[Mm]",
            )

        if self.geometry == "Spherical":
            pi2rad = 180 / np.pi
            print("### calculation domain ###")
            print(
                "xmax - rstar = ",
                "{:6.2f}".format((self.xmax - self.rstar) * 1.0e-8),
                "[Mm],  xmin - rstar = ",
                "{:6.2f}".format((self.xmin - self.rstar) * 1.0e-8),
                "[Mm]",
            )
            print(
                "ymax        = ",
                "{:6.2f}".format(self.ymax * pi2rad),
                "[rad], ymin        = ",
                "{:6.2f}".format(self.ymin * pi2rad),
                "[rad]",
            )
            print(
                "zmax        = ",
                "{:6.2f}".format(self.zmax * pi2rad),
                "[rad], zmin        = ",
                "{:6.2f}".format(self.zmin * pi2rad),
                "[rad]",
            )

        if self.geometry == "YinYang":
            pi2rad = 180 / np.pi
            print("### calculation domain ###")
            print(
                "xmax - rstar = ",
                "{:6.2f}".format((self.xmax - self.rstar) * 1.0e-8),
                "[Mm],  xmin - rstar = ",
                "{:6.2f}".format((self.xmin - self.rstar) * 1.0e-8),
                "[Mm]",
            )
            print("Yin-Yang grid is used to cover the whole sphere")

        print("")
        print(RED + "### number of grid ###" + END)
        print("(ix,jx,kx)=(", self.ix, ",", self.jx, ",", self.kx, ")")

        print("")
        print(RED + "### calculation time ###" + END)
        print("time step (nd) =", self.nd)
        print("time step (nd_tau) =", self.nd_tau)
        t = self.data.time_read(self.nd, verbose=False)
        print("time =", "{:.2f}".format(t / 3600), " [hour]")

    def _update_json_template(
        self,
        output_file=os.path.dirname(os.path.abspath(__file__))
        + "/"
        + "parameters.json",
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
        json_file=os.path.dirname(os.path.abspath(__file__)) + "/" + "parameters.json",
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

        docstring = ""
        docstring += "\n"
        docstring += gene_space + "Attributes\n"
        docstring += gene_space + "----------\n"
        try:
            with open(json_file, "r") as f:
                saved_attributes = json.load(f)
        except FileNotFoundError:
            saved_attributes = {}
        for key in self.__dict__.keys():
            value_type = type(self.__dict__[key])
            type_name = (
                value_type.__name__
                if value_type.__module__ == "builtins"
                else f"{value_type.__module__}.{value_type.__name__}"
            )

            docstring += gene_space + f"{key} : {type_name}\n"
            if key in saved_attributes:
                docstring += gene_space + desc_space + saved_attributes[key] + "\n"
            else:
                docstring += "\n"
        self.__class__.__doc__ = self.__class__.__doc__ + docstring
