import gc
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np

import pyR2D2


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
        raise AttributeError(
            f"'{self.__class__.__name__}' object has no attribute '{name}'"
        )


class _BaseRemapReader(_BaseReader):
    """
    Base class for remap data readers
    """

    def __init__(self, data):
        """
        Initialize the _BaseRemapReader

        Parameters
        ----------
        data : pyR2D2.Data
                pyR2D2.Data instance. This is used to access the attributes of pyR2D2.Data, such as x, y, z, etc.
        """

        self.data = data

        for key in self.remap_keys + self.remap_keys_add:
            self.__dict__[key] = None

        self._add_docstring()

    def _allocate_remap_qq(self, ijk: tuple, keys: list = None):
        """
        Allocates memory for remap/qq/ data

        Paramters
        ---------
        ijk : tuple
            size of data

        keys : list, optional
            list of keys to allocate. If None, all keys are allocated. If "all", all keys are allocated. By default None.
        """

        if keys is None:
            keys = self.remap_keys + self.remap_keys_add
        else:
            if isinstance(keys, str):
                if keys == "all":
                    keys = self.remap_keys + self.remap_keys_add
                else:
                    keys = [keys]

        for key in keys:
            if key not in self.remap_keys + self.remap_keys_add:
                raise ValueError(
                    f"key should be one of {self.remap_keys + self.remap_keys_add}, but got {key}"
                )
            memflag = True
            # check if the shape of allocated array is the same as ijk.
            # If not, allocate a new array. This is to avoid memory error when the shape of data changes.
            if self.__dict__.get(key, None) is not None:
                memflag = not self.__dict__[key].shape == ijk
            if self.__dict__.get(key, None) is None or memflag:
                self.__dict__[key] = np.empty(ijk, dtype=np.float32)

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
        filepath : pathlib.Path
            file path of remap/qq/
        """
        if (self.datadir / "remap" / "qq" / "00000").is_dir():
            cnou = f"{np0//1000:05d}"
            cno = f"{np0:08d}"
            filepath = (
                self.datadir / "remap" / "qq" / cnou / cno / f"qq.dac.{n:08d}.{np0:08d}"
            )
        else:
            # directoryを分けない古いバージョン対応
            filepath = self.datadir / "remap" / "qq" / f"qq.dac.{n:08d}.{np0:08d}"
        return filepath

    def _get_filepath_remap_zarr(self, n: int):
        """
        Filepath for remap date in zarr format

        Parameters
        ----------
        n : int
            time step

        Return
        -------
        filepath : pathlib.Path
            file path of remap/qq/ in zarr format
        """

        return self.datadir / "remap" / "qq" / "zarr" / f"qq.{n:08d}.zarr"

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
        docstring += "        div B cleaning variable (not included in zarr file)\n"
        docstring += "    pr: numpy.ndarray, float\n"
        docstring += "        pressure (not included in zarr file)\n"
        docstring += "    te: numpy.ndarray, float\n"
        docstring += "        temperature (not included in zarr file)\n"
        docstring += "    op: numpy.ndarray, float\n"
        docstring += "        opacity (not included in zarr file)\n"

        self.__class__.__doc__ = self.__class__.__doc__ + docstring


class XSelect(_BaseRemapReader):
    """
    Class for 2D selected data at a certain x

    Important
    ---------
    pyR2D2.Data class can access this class as :code:`pyR2D2.Data.qx`

    """

    def read(
        self,
        xs: float,
        n: int,
        keys="all",
        zarr_flag: bool = False,
        max_workers: int = 1,
        verbose: bool = False,
    ):
        """
        Reads 2D selected data at a certain x

        Parameters
        ----------
        xs : float
            A selected height for data
        n : int
            A selected time step for data
        keys : str, list, or tuple
            Kind of variable.
        zarr_flag : bool
            If True, read from zarr format instead of binary format
        max_workers : int
            The number of worker threads to use for reading files
        verbose : bool
            If True, print verbose output
        """

        i0 = np.argmin(np.abs(self.p.x - xs))
        ir0 = self.i2ir[i0]

        missing = False
        # check if original binary files exists
        for jr0 in range(1, self.jxr + 1):
            np0 = self.np_ijr[ir0 - 1, jr0 - 1]
            if jr0 == self.jr[np0]:
                filepath = self._get_filepath_remap_qq(n, np0)
                if not filepath.exists():
                    missing = True
                    break

        if missing and not zarr_flag:
            if verbose:
                print("Some original binary files are missing.")
                print("Trying to read from zarr file.")
            zarr_flag = True

        if zarr_flag:
            if keys == "all":
                names = keys
            else:
                if isinstance(keys, str):
                    names = [keys, "x", "y", "z"]
                else:
                    names = list(keys) + ["x", "y", "z"]

            # If zarr_flag is True, try to read from zarr file
            zarr_path = self._get_filepath_remap_zarr(n)
            if not zarr_path.exists():
                raise FileNotFoundError(f"Zarr file does not exist: {zarr_path}")
            qq = pyR2D2.zarr_util.load(zarr_path, names=names, i0=i0, i1=i0 + 1)
            for key in qq.keys():
                if key in ["x", "y", "z"]:
                    self.__dict__[key] = qq[key]
                else:
                    self.__dict__[key] = qq[key][0, :, :].squeeze()
        else:
            if isinstance(keys, str):
                if keys == "all":
                    keys_input = self.remap_keys + self.remap_keys_add
                else:
                    keys_input = [keys]
            elif isinstance(keys, (list, tuple)):
                keys_input = list(keys)
            else:
                raise TypeError("keys must be str, list, or tuple")

            for key in keys_input:
                if key not in self.remap_keys + self.remap_keys_add:
                    print("######")
                    print("key =", key)
                    print(
                        "key should be one of ", self.remap_keys + self.remap_keys_add
                    )
                    print("return")
                    return

            self._allocate_remap_qq(ijk=[self.jx, self.kx], keys=keys_input)
            dtype = np.dtype(self.endian + "f4")

            target_nps = []
            for jr0 in range(1, self.jxr + 1):
                np0 = self.np_ijr[ir0 - 1, jr0 - 1]

                if jr0 == self.jr[np0]:
                    target_nps.append(np0)

            def _read_one(np0: int):
                n_ijk = self.iixl[np0] * self.jjxl[np0] * self.kx
                byte_n_ijk = n_ijk * 4  # for float32
                byte_n_ijk_mtype = byte_n_ijk * self.mtype

                offset = {
                    "pr": byte_n_ijk_mtype,
                    "te": byte_n_ijk_mtype + byte_n_ijk,
                    "op": byte_n_ijk_mtype + 2 * byte_n_ijk,
                }
                filepath = self._get_filepath_remap_qq(n, np0)
                out = {}
                with open(filepath, "rb") as f:
                    for key in keys_input:
                        if key in self.remap_keys:
                            m_idx = self.remap_keys.index(key)
                            f.seek(m_idx * byte_n_ijk, os.SEEK_SET)
                            buf = np.fromfile(f, dtype=dtype, count=n_ijk)
                            out[key] = buf.reshape(
                                (self.iixl[np0], self.jjxl[np0], self.kx), order="F"
                            )[i0 - self.iss[np0], :, :]
                        elif key in self.remap_keys_add:
                            f.seek(offset[key], os.SEEK_SET)
                            buf = np.fromfile(f, dtype=dtype, count=n_ijk)
                            out[key] = buf.reshape(
                                (self.iixl[np0], self.jjxl[np0], self.kx), order="F"
                            )[i0 - self.iss[np0], :, :]
                    j0, j1 = self.jss[np0], self.jee[np0] + 1
                    return (np0, j0, j1, out)

            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = [executor.submit(_read_one, np0) for np0 in target_nps]
                for future in as_completed(futures):
                    np0, j0, j1, out = future.result()
                    for key, arr in out.items():
                        self.__dict__[key][j0:j1, :] = arr

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

                for key, m in zip(self.remap_keys, range(self.mtype)):
                    self.__dict__[key][
                        self.iss[np0] : self.iee[np0] + 1,
                        self.jss[np0] : self.jee[np0] + 1,
                    ] = qqq["qq"].reshape(
                        (self.iixl[np0], self.jjxl[np0], self.kx, self.mtype), order="F"
                    )[
                        :, :, k0, m
                    ]

                for key in self.remap_keys_add:
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
                    for key, m in zip(self.remap_keys, range(self.mtype)):
                        self.__dict__[key][
                            :, self.jss[np0] : self.jee[np0] + 1, :
                        ] = qqq["qq"].reshape(
                            (self.iixl[np0], self.jjxl[np0], self.kx, self.mtype),
                            order="F",
                        )[
                            :, :, :, m
                        ]
                    for key in self.remap_keys_add:
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

    zarr_keys = ["ro", "vx", "vy", "vz", "bx", "by", "bz", "se"]

    def read(self, n: int, keys="all", max_workers: int = 1, zarr_flag: bool = False, verbose: bool = False):
        """
        Reads 3D full data
        The data is stored in self.qq dictionary

        Parameters
        ----------
        n : int
            A selected time step for data
        keys : str
            Kind of variable. Options are:
                - "ro" : density
                - "vx", "vy", "vz": velocity
                - "bx", "by", "bz": magnetic field
                - "se": entropy
                - "ph": div B cleaning (not included in zarr file)
                - "te": temperature (not included in zarr file)
                - "op": Opacity (not included in zarr file)

        max_workers : int
            Number of workers for parallel reading of binary files

        zarr_flag : bool
            If True, read from zarr format instead of binary format

        verbose : bool
            If True, print verbose output

        Notes
        -----
        If keys is 'all', all variables are read
        """

        # check if original binary files exists
        missing = False
        for ir0 in range(1, self.ixr + 1):
            for jr0 in range(1, self.jxr + 1):
                np0 = self.np_ijr[ir0 - 1, jr0 - 1]
                if ir0 == self.ir[np0] and jr0 == self.jr[np0]:

                    filepath = self._get_filepath_remap_qq(n, np0)
                    if not filepath.exists():
                        missing = True
                        break
            if missing:
                break

        if missing and not zarr_flag:
            if verbose:
                print("Some original binary files are missing.")
                print("Trying to read from zarr file.")
            zarr_flag = True

        if zarr_flag:
            if keys == "all":
                names = keys
            else:
                if isinstance(keys, str):
                    names = [keys, "x", "y", "z"]
                else:
                    names = list(keys) + ["x", "y", "z"]

            zarr_filepath = self._get_filepath_remap_zarr(n)
            if not zarr_filepath.exists():
                print(f"Zarr file does not exist at n={n}.")
                return

            (qq, params) = pyR2D2.zarr_util.load(
                zarr_filepath, with_attrs=True, names=names
            )

            for key in qq.keys():
                self.__dict__[key] = qq[key]
            self.params = params

        else:

            if isinstance(keys, str):
                if keys == "all":
                    keys_input = self.remap_keys + self.remap_keys_add
                else:
                    keys_input = [keys]
            elif isinstance(keys, (list, tuple)):
                keys_input = list(keys)
            else:
                raise TypeError("keys must be str, list, or tuple")

            for key in keys_input:
                if key not in self.remap_keys + self.remap_keys_add:
                    print("######")
                    print("key =", key)
                    print(
                        "key should be one of ", self.remap_keys + self.remap_keys_add
                    )
                    print("return")
                    return

            self._allocate_remap_qq(ijk=[self.ix, self.jx, self.kx], keys=keys_input)
            dtype = np.dtype(self.endian + "f4")

            target_nps = []
            for ir0 in range(1, self.ixr + 1):
                for jr0 in range(1, self.jxr + 1):
                    np0 = self.np_ijr[ir0 - 1, jr0 - 1]
                    if ir0 == self.ir[np0] and jr0 == self.jr[np0]:
                        target_nps.append(np0)

            def _read_one(np0: int):
                n_ijk = self.iixl[np0] * self.jjxl[np0] * self.kx
                byte_n_ijk = n_ijk * 4  # for float32
                byte_n_ijk_mtype = byte_n_ijk * self.mtype

                offset = {
                    "pr": byte_n_ijk_mtype,
                    "te": byte_n_ijk_mtype + byte_n_ijk,
                    "op": byte_n_ijk_mtype + 2 * byte_n_ijk,
                }
                filepath = self._get_filepath_remap_qq(n, np0)
                out = {}
                with open(filepath, "rb") as f:
                    for key in keys_input:
                        if key in self.remap_keys:
                            m_idx = self.remap_keys.index(key)
                            f.seek(m_idx * byte_n_ijk, os.SEEK_SET)
                            buf = np.fromfile(f, dtype=dtype, count=n_ijk)
                            out[key] = buf.reshape(
                                (self.iixl[np0], self.jjxl[np0], self.kx), order="F"
                            )
                        elif key in self.remap_keys_add:
                            f.seek(offset[key], os.SEEK_SET)
                            buf = np.fromfile(f, dtype=dtype, count=n_ijk)
                            out[key] = buf.reshape(
                                (self.iixl[np0], self.jjxl[np0], self.kx), order="F"
                            )
                i0, i1 = self.iss[np0], self.iee[np0] + 1
                j0, j1 = self.jss[np0], self.jee[np0] + 1
                return (np0, i0, i1, j0, j1, out)

            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = [executor.submit(_read_one, np0) for np0 in target_nps]
                for future in as_completed(futures):
                    np0, i0, i1, j0, j1, out = future.result()
                    for key, arr in out.items():
                        self.__dict__[key][i0:i1, j0:j1, :] = arr

    def compress(
        self,
        n: int,
        zarr_filepath: str = None,
        i_start: int = None,
        i_size: int = None,
        j_start: int = None,
        j_size: int = None,
        k_start: int = None,
        k_size: int = None,
        chunks3d: tuple = None,
        keys: list = zarr_keys,
        overwrite: bool = False,
        max_workers: int = 1,
        lightweight: bool = False,
    ):
        """
        Compresses the remap/qq/ data for a given time step n into zarr format

        Parameters
        ----------
        n : int
            A selected time step for data
        zarr_filepath : str, optional
            File path for the output zarr file. If None, the default file path is used. By default None.
        i_start : int, optional
            Starting index in i direction for the data to be compressed. If None, it starts from 0. By default None.
        i_size : int, optional
            Size in i direction for the data to be compressed. If None, it uses the full size in i direction. By default None.
        j_start : int, optional
            Starting index in j direction for the data to be compressed. If None, it starts from 0. By default None.
        j_size : int, optional
            Size in j direction for the data to be compressed. If None, it uses the full size in j direction. By default None.
        k_start : int, optional
            Starting index in k direction for the data to be compressed. If None, it starts from 0. By default None.
        k_size : int, optional
            Size in k direction for the data to be compressed. If None, it uses the full size in k direction. By default None.
        chunks3d : tuple, optional
            Chunk size for the 3D data in zarr format. If None, the chunk size is automatically determined. By default None.
        keys : list, optional
            List of variables to be compressed. By default, it includes all variables in zarr_keys.
        overwrite : bool, optional
            If True, overwrite the existing zarr file. If False and the zarr file already exists, the function will print a message and return without doing anything. By default False.
        max_workers : int, optional
            Number of workers for parallel reading of binary files when lightweight is True. By default 1.
        lightweight : bool, optional
            If True, the function reads and writes each variable separately to save memory. This is useful when the data is too large to fit in memory. By default False (all variables are read and written together).

        """

        if zarr_filepath is None:
            zarr_filepath = self._get_filepath_remap_zarr(n)
        else:
            zarr_filepath = Path(zarr_filepath)

        if not overwrite and zarr_filepath.exists():
            print(
                f"File {zarr_filepath} already exists. Set overwrite=True to overwrite it."
            )
            return

        i_start = 0 if i_start is None else i_start
        i_size = self.ix if i_size is None else i_size
        j_start = 0 if j_start is None else j_start
        j_size = self.jx if j_size is None else j_size
        k_start = 0 if k_start is None else k_start
        k_size = self.kx if k_size is None else k_size

        if chunks3d is None:
            i_size_chunk = i_size
            j_size_chunk = j_size
            k_size_chunk = k_size

            target_size = 2**24  # 64 MB for float32

            i_count = 1
            jk_count = 1
            jk_priority = 3  # priority for j and k to be halved compared to i
            while i_size_chunk * j_size_chunk * k_size_chunk > target_size:
                if i_count <= jk_count * jk_priority:
                    i_size_chunk = max(i_size_chunk // 2, 1)
                    i_count += 1
                else:
                    j_size_chunk = max(j_size_chunk // 2, 1)
                    k_size_chunk = max(k_size_chunk // 2, 1)
                    jk_count += 1

            chunks3d = (i_size_chunk, j_size_chunk, k_size_chunk)

        params_dict = {
            "i_start": i_start,
            "i_size": i_size,
            "j_start": j_start,
            "j_size": j_size,
            "k_start": k_start,
            "k_size": k_size,
        }

        vars_dict = {}
        vars_dict["x"] = self.x[i_start : i_start + i_size]
        vars_dict["y"] = self.y[j_start : j_start + j_size]
        vars_dict["z"] = self.z[k_start : k_start + k_size]
        pyR2D2.zarr_util.save(
            zarr_filepath, vars_dict, params_dict, chunks3d=chunks3d, mode="w"
        )

        if lightweight:
            for key in keys:
                print(key)
                self.read(n=n, keys=key, max_workers=max_workers, zarr_flag=False)
                arr = self.__dict__[key][
                    i_start : i_start + i_size,
                    j_start : j_start + j_size,
                    k_start : k_start + k_size,
                ]
                pyR2D2.zarr_util.save(
                    zarr_filepath, {key: arr}, chunks3d=chunks3d, mode="a"
                )
                del self.__dict__[key]
                del arr
                gc.collect()
        else:

            self.read(n=n, keys=keys, max_workers=max_workers, zarr_flag=False)

            vars_dict = {}
            for key in keys:
                vars_dict[key] = self.__dict__[key][
                    i_start : i_start + i_size,
                    j_start : j_start + j_size,
                    k_start : k_start + k_size,
                ]
            pyR2D2.zarr_util.save(
                zarr_filepath, vars_dict, params_dict, chunks3d=chunks3d, mode="a"
            )

    @staticmethod
    def _check_core(qq1, qq2, key):
        if (abs(qq1 - qq2).sum() / (abs(qq1).sum() + 1e-12)) > 1e-6:
            print(f"Data mismatch for {key}.")
            return False
        return True

    def check(
        self,
        n: int,
        keys: list = zarr_keys,
        max_workers: int = 1,
        lightweight: bool = False,
    ):
        """
        Check if the remap/qq/ file exists for all MPI processes for a given time step n

        Parameters
        ----------
        n : int
            A selected time step for data
        keys : list or str
            List of values to check. If the check fails for any of the values, the function returns False. By default, all values are checked.
        lightweight : bool
            If True, check is done for each value separately to save memory. This is useful when the data is too large to fit in memory. By default, False (all values are checked together).
        """

        zarr_filepath = self._get_filepath_remap_zarr(n)
        if not zarr_filepath.exists():
            print(
                f"Zarr file does not exist at n={n}. Please run FullData.compress() to create it."
            )
            return False

        for key in keys:
            if not key in pyR2D2.zarr_util.list_vars(zarr_filepath):
                print(
                    f"Variable {key} does not exist in zarr file at n={n}. Please run FullData.compress() to create it."
                )
                return False

        for np0 in range(self.npe):
            filepath = self._get_filepath_remap_qq(n, np0)
            if self.iixl[np0] * self.jjxl[np0] != 0:
                if not filepath.exists():
                    print(f"File {filepath} does not exist. Anyway you can delete it.")
                    return True

        if lightweight:
            for i in range(0, self.ix, self.ix // 3):
                self.qx.read(self.x[i], n=n, zarr_flag=True)
                qq_copy = {}
                for key in keys:
                    if key in self.remap_keys + self.remap_keys_add:
                        qq_copy[key] = self.qx.__dict__[key].copy()
                self.qx.read(self.x[i], n=n, zarr_flag=False)
                for key in keys:
                    if key in self.remap_keys + self.remap_keys_add:
                        if not self._check_core(
                            qq_copy[key], self.qx.__dict__[key], key
                        ):
                            return False
        else:
            self.read(n=n, keys=keys, zarr_flag=True, max_workers=max_workers)
            qq_copy = {}
            for key in keys:
                if key in self.remap_keys + self.remap_keys_add:
                    qq_copy[key] = self.__dict__[key].copy()

            self.read(n=n, keys=keys, zarr_flag=False)
            for key in keys:
                if not self._check_core(qq_copy[key], self.__dict__[key], key):
                    return False

        print(f"Check passed at n={n}")
        return True

    def delete(self, n: int, force: bool = False, lightweight: bool = False):
        """
        Delete the remap/qq/ file for a given time step n

        Parameters
        ----------
        n : int
            A selected time step for data
        force : bool
            If True, delete the files even if the check fails
        lightweight : bool
            If True, delete the zarr file for each value separately to save memory. This is useful when the data is too large to fit in memory. By default, False (all values are deleted together).
        """

        if force:
            check_flag = True
        else:
            check_flag = self.check(n=n, lightweight=lightweight)

        if check_flag:
            if (self.datadir / "remap" / "qq" / "00000").is_dir():
                base = self.datadir / "remap" / "qq"
                print(
                    "Deleting files in",
                    f"{base}/*/*/qq.dac.{n:08d}.*",
                )
                for np0 in range(self.npe):

                    filepath = self._get_filepath_remap_qq(n, np0)
                    if filepath.exists():
                        filepath.unlink()
            else:
                print(
                    "Deleting",
                    f"{self.datadir / 'remap' / 'qq'}/qq.dac.{n:08d}.*",
                )
                pattern = f"qq.dac.{n:08d}.*"
                subprocess.run(
                    [
                        "find",
                        str(self.datadir / "remap" / "qq"),
                        "-maxdepth",
                        "1",
                        "-name",
                        pattern,
                        "-delete",
                    ],
                    check=True,
                )

    def clear(self, keys="all"):
        """
        Free the memory used for remap/qq/ data
        """

        if isinstance(keys, str):
            if keys == "all":
                keys = self.remap_keys + self.remap_keys_add
            else:
                keys = [keys]

        for key in keys:
            if hasattr(self, key):
                self.__dict__[key] = None

        gc.collect()


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
        keys: list,
        x0: float = None,
        x1: float = None,
        y0: float = None,
        y1: float = None,
        z0: float = None,
        z1: float = None,
        zarr_flag: bool = False,
        verbose: bool = False,
    ):
        """
        Reads 3D restricted-area data
        The data is stored in self.qr dictionary

        Parameters
        ----------
        n : int
            A selected time step for data
        keys : list of str
            See :meth:`R2D2.FullData` for options

        x0, y0, z0 : float
            Minimum x, y, z by default, the minimum of the domain
        x1, y1, z1 : float
            Maximum x, y, z by default, the maximum of the domain
        zarr_flag : bool
            If True, read from zarr format instead of binary format
        verbose : bool
            If True, print verbose output
        """

        if x0 is None:
            x0 = self.x[0]
        if x1 is None:
            x1 = self.x[-1]
        if y0 is None:
            y0 = self.y[0]
        if y1 is None:
            y1 = self.y[-1]
        if z0 is None:
            z0 = self.z[0]
        if z1 is None:
            z1 = self.z[-1]

        i0, i1 = np.argmin(abs(self.x - x0)), np.argmin(abs(self.x - x1))
        j0, j1 = np.argmin(abs(self.y - y0)), np.argmin(abs(self.y - y1))
        k0, k1 = np.argmin(abs(self.z - z0)), np.argmin(abs(self.z - z1))
        ixr = i1 - i0 + 1
        jxr = j1 - j0 + 1
        kxr = k1 - k0 + 1

        # check if original binary files exists
        missing = False
        for ir0 in range(1, self.ixr + 1):
            for jr0 in range(1, self.jxr + 1):
                np0 = self.np_ijr[ir0 - 1, jr0 - 1]
                if ir0 == self.ir[np0] and jr0 == self.jr[np0]:

                    filepath = self._get_filepath_remap_qq(n, np0)
                    if not filepath.exists():
                        missing = True
                        break
            if missing:
                break

        if missing and not zarr_flag:
            if verbose:
                print("Some original binary files are missing.")
                print("Trying to read from zarr file.")
            zarr_flag = True

        if zarr_flag:
            if keys == "all":
                names = keys
            else:
                if isinstance(keys, str):
                    names = [keys, "x", "y", "z"]
                else:
                    names = list(keys) + ["x", "y", "z"]

            zarr_filepath = self._get_filepath_remap_zarr(n)
            if not zarr_filepath.exists():
                print(f"Zarr file does not exist at n={n}.")
                return

            (qq, params) = pyR2D2.zarr_util.load(
                zarr_filepath, with_attrs=True, names=names, i0=i0, i1=i1 + 1, j0=j0, j1=j1 + 1, k0=k0, k1=k1 + 1
            )

            for key in qq.keys():
                if key in ["x", "y", "z"]:
                    self.__dict__[key] = qq[key]
                else:
                    self.__dict__[key] = qq[key][:, j0 : j1 + 1, :].squeeze()

        else:            
            if isinstance(keys, str):
                if keys == "all":
                    keys_input = self.remap_keys + self.remap_keys_add
                else:
                    keys_input = [keys]
            elif isinstance(keys, (list, tuple)):
                keys_input = list(keys)
            else:
                raise TypeError("keys must be str, list, or tuple")

            for key in keys_input:
                self.__dict__[key] = np.zeros((ixr, jxr, kxr), dtype=np.float32)

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

                            for key in keys_input:
                                isrt_rcv = max([0, self.iss[np0] - i0])
                                iend_rcv = min([ixr, self.iee[np0] - i0 + 1])
                                jsrt_rcv = max([0, self.jss[np0] - j0])
                                jend_rcv = min([jxr, self.jee[np0] - j0 + 1])

                                isrt_snd = isrt_rcv - (self.iss[np0] - i0)
                                iend_snd = isrt_snd + (iend_rcv - isrt_rcv)
                                jsrt_snd = jsrt_rcv - (self.jss[np0] - j0)
                                jend_snd = jsrt_snd + (jend_rcv - jsrt_rcv)

                                if key in self.remap_keys:
                                    m = self.remap_keys.index(key)

                                    self.__dict__[key][
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
                                    self.__dict__[key][
                                        isrt_rcv:iend_rcv, jsrt_rcv:jend_rcv, :
                                    ] = qqq[key].reshape(
                                        (self.iixl[np0], self.jjxl[np0], self.kx), order="F"
                                    )[
                                        isrt_snd:iend_snd, jsrt_snd:jend_snd, k0 : k1 + 1
                                    ]

        self.info = {"x0": self.x[i0], "x1": self.x[i1], "y0": self.y[j0], "y1": self.y[j1], "z0": self.z[k0], "z1": self.z[k1]}

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
        pressure (not included in zarr file)
    te : numpy.ndarray, float
        temperature (not included in zarr file)
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
        radiative flux at the optical depth

    Important
    ---------
    Each attibutes have three values at optical depths 1, 0.1, and 0.01
    These correspond to, for example, rt, rt01, and rt001


    """

    zarr_keys = [
        "rt",
        "ro",
        "ro01",
        "ro001",
        "se",
        "se01",
        "se001",
        "vx",
        "vx01",
        "vx001",
        "vy",
        "vy01",
        "vy001",
        "vz",
        "vz01",
        "vz001",
        "bx",
        "bx01",
        "bx001",
        "by",
        "by01",
        "by001",
        "bz",
        "bz01",
        "bz001",
        "he",
        "he01",
        "he001",
        "fr",
        "fr01",
        "fr001",
    ]

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

    def _get_filepath_optical_depth(self, n: int):
        return self.datadir / "tau" / f"qq.dac.{n:08d}"

    def _get_filepath_optical_depth_zarr(self, n: int):
        return self.datadir / "tau" / "zarr" / f"qq.{n:08d}.zarr"

    def read(
        self,
        n: int,
        zarr_flag: bool = False,
        verbose: bool = True,
        keys: str | list | tuple = "all",
    ):
        """
        Reads 2D data at certain optical depths.
        The data is stored in self.qt dictionary.
        In this version the selected optical depth is 1, 0.1, and 0.01

        Parameters
        ----------
        n : int
            A selected time step for data
        zarr_flag : bool, optional
            Whether to read from zarr file, by default True
        verbose : bool, optional
            Whether to print verbose output, by default True
        keys : str or list or tuple, optional
            List of values to read. If 'all', all values are read. By default 'all'.
            It works only when zarr_flag is True. If zarr_flag is False, all values are read regardless of keys.
        """

        if not zarr_flag:
            if not self._get_filepath_optical_depth(n).exists():
                if verbose:
                    print(
                        f"File {self._get_filepath_optical_depth(n)} does not exist. Trying to read from zarr file."
                    )
                zarr_flag = True

        if zarr_flag:
            if isinstance(keys, str):
                if keys == "all":
                    names = "all"
                else:
                    names = [keys]
            elif isinstance(keys, (list, tuple)):
                names = list(keys)
            else:
                raise TypeError("keys must be str, list, or tuple")

            zarr_filepath = self._get_filepath_optical_depth_zarr(n)
            if not zarr_filepath.exists():
                print(f"Zarr file does not exist at n={n}.")
                return

            qq = pyR2D2.zarr_util.load(zarr_filepath, names=names)

            for key in qq.keys():
                self.__dict__[key] = qq[key]
            return

        else:

            with open(self.datadir / "tau" / f"qq.dac.{n:08d}", "rb") as f:
                qq = np.fromfile(
                    f, self.endian + "f", self.m_tu * self.m_in * self.jx * self.kx
                )

            for key, mk in zip(self.value_keys, range(len(self.value_keys))):
                for tau, mt in zip(["", "01", "001"], range(3)):
                    self.__dict__[key + tau] = qq.reshape(
                        (self.m_tu, self.m_in, self.jx, self.kx), order="F"
                    )[mt, mk, :, :]

    def compress(
        self,
        n: int,
        zarr_filepath: str = None,
        keys: list = zarr_keys,
        overwrite: bool = False,
    ):

        if zarr_filepath is None:
            zarr_filepath = self._get_filepath_optical_depth_zarr(n)
        else:
            zarr_filepath = Path(zarr_filepath)

        if not overwrite and zarr_filepath.exists():
            print(
                f"File {zarr_filepath} already exists. Set overwrite=True to overwrite it."
            )
            return

        self.read(n=n, zarr_flag=False)
        chunk_max = 4096
        chunks3d = (1, min(self.jx, chunk_max), min(self.kx, chunk_max))
        vars_dict = {}
        for key in keys:
            vars_dict[key] = self.__dict__[key]

        params_dict = {}
        pyR2D2.zarr_util.save(
            zarr_filepath, vars_dict, params_dict, chunks3d=chunks3d, mode="w"
        )

    def check(self, n: int, keys: list = zarr_keys):
        """
        Check if the tau/qq/ file exists for a given time step n

        Parameters
        ----------
        n : int
            A selected time step for data
        keys : list or str
            List of values to check. If the check fails for any of the values, the function returns False. By default, all values are checked.
        """

        zarr_filepath = self._get_filepath_optical_depth_zarr(n)
        if not zarr_filepath.exists():
            print(
                f"Zarr file does not exist at n={n}. Please run OpticalDepth.compress() to create it."
            )
            return False

        for key in keys:
            if key not in pyR2D2.zarr_util.list_vars(zarr_filepath):
                print(
                    f"Variable {key} does not exist in zarr file at n={n}. Please run OpticalDepth.compress() to create it."
                )
                return False

        self.read(n=n, zarr_flag=True)
        qq_copy = {}
        for key in keys:
            qq_copy[key] = self.__dict__[key].copy()

        self.read(n=n, zarr_flag=False)
        for key in keys:
            if not np.allclose(qq_copy[key], self.__dict__[key], rtol=1e-6, atol=1e-12):
                print(f"Data mismatch for {key}.")
                return False

        print(f"Check passed at n={n}")
        return True

    def delete(self, n: int, force: bool = False):
        """
        Delete the tau/qq/ file for a given time step n

        Parameters
        ----------
        n : int
            A selected time step for data
        force : bool
            If True, delete the files even if the check fails
        """

        if force:
            check_flag = True
        else:
            check_flag = self.check(n=n)

        if check_flag:
            filepath = self._get_filepath_optical_depth(n)
            if filepath.exists():
                print(f"Deleting {filepath}.")
                filepath.unlink()


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

        if self.nx*self.ny*self.nz == 8:

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
        with open(self.datadir / "remap" / "vl" / f"vl_xy.dac.{n:08d}", "rb") as f:
            vl = np.fromfile(
                f, self.endian + "f", self.m2d_xy * self.ix * self.jx
            ).reshape((self.ix, self.jx, self.m2d_xy), order="F")

        for m in range(self.m2d_xy):
            self.__dict__[self.cl[m]] = vl[:, :, m]

        # read xz plane data
        with open(self.datadir / "remap" / "vl" / f"vl_xz.dac.{n:08d}", "rb") as f:
            vl = np.fromfile(
                f, self.endian + "f", self.m2d_xz * self.ix * self.kx
            ).reshape((self.ix, self.kx, self.m2d_xz), order="F")

        for m in range(self.m2d_xz):
            self.__dict__[self.cl[m + self.m2d_xy]] = vl[:, :, m]

        # read flux related value
        with open(self.datadir / "remap" / "vl" / f"vl_flux.dac.{n:08d}", "rb") as f:
            vl = np.fromfile(
                f, self.endian + "f", self.m2d_flux * (self.ix + 1) * self.jx
            ).reshape((self.ix + 1, self.jx, self.m2d_flux), order="F")

        for m in range(self.m2d_flux):
            self.__dict__[self.cl[m + self.m2d_xy + self.m2d_xz]] = vl[:, :, m]

        # read spectra
        if self.geometry == "YinYang":
            with open(
                self.datadir / "remap" / "vl" / f"vl_spex.dac.{n:08d}", "rb"
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
        output_file=Path(__file__).resolve().parent / "OnTheFly.json",
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
        json_file=Path(__file__).resolve().parent / "OnTheFly.json",
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

    zarr_keys = ["ro", "vx", "vy", "vz", "bx", "by", "bz", "se"]

    def __init__(self, data):
        self.data = data

        for key in self.data.remap_keys + self.data.remap_keys_add[:-1]:
            self.__dict__[key] = None

    def _get_filepath_slice(self, n_slice, direc, n, postfix):
        return (
            self.datadir / "slice" / f"qq{direc}{postfix}.dac.{n:08d}.{n_slice+1:08d}"
        )

    def _get_filepath_slice_zarr(self, n, direc):
        return self.datadir / "slice" / "zarr" / direc / f"qq.{n:08d}.zarr"

    def _get_postfixes(self):
        if self.geometry == "YinYang":
            return ["_yin", "_yan"]
        else:
            return [""]

    def xyz_slice_select(self, direc):
        if direc == "x":
            return self.x_slice
        elif direc == "y":
            return self.y_slice
        elif direc == "z":
            return self.z_slice
        else:
            raise ValueError("direc should be 'x', 'y', or 'z'")

    def read(self, n_slice, direc, n, zarr_flag: bool = False):
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
        zarr_flag : bool
            If True, read from zarr format instead of binary format. By default, False (read from binary format).
        """

        postfixes = self._get_postfixes()

        # check if the files exist
        if not zarr_flag:
            for postfix in postfixes:
                filepath = self._get_filepath_slice(n_slice, direc, n, postfix)
                if not filepath.exists():
                    zarr_flag = True
                    print(
                        f"File does not exist at n={n} for slice {direc} with postfix {postfix}."
                    )
                    print("Trying to read from zarr file.")
                    break

        if zarr_flag:
            zarr_filepath = self._get_filepath_slice_zarr(n, direc)
            if not zarr_filepath.exists():
                print(f"Zarr file does not exist at n={n} for slice {direc}.")
                return

            i0, j0, k0 = 0, 0, 0
            i1, j1, k1 = self.ix, self.jx, self.kx
            if direc == "x":
                i0, i1 = n_slice, n_slice + 1
            elif direc == "y":
                j0, j1 = n_slice, n_slice + 1
            elif direc == "z":
                k0, k1 = n_slice, n_slice + 1

            names = []
            for key in self.zarr_keys:
                for postfix in postfixes:
                    names.append(key + postfix)

            qq = pyR2D2.zarr_util.load(
                zarr_filepath,
                names=names,
                i0=i0,
                i1=i1,
                j0=j0,
                j1=j1,
                k0=k0,
                k1=k1,
            )
            for key, value in qq.items():
                self.__dict__[key] = value.squeeze()
        else:
            for postfix in postfixes:
                with open(
                    self._get_filepath_slice(n_slice, direc, n, postfix),
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
                    self.remap_keys + self.remap_keys_add[:-1], range(self.mtype + 2)
                ):
                    self.__dict__[key + postfix] = qq.reshape(
                        (n1, n2, self.mtype + 2), order="F"
                    )[:, :, m]

            # self.info = {}
            # self.info["direc"] = direc
            # if direc == "x":
            #     slice = self.x_slice
            # if direc == "y":
            #     slice = self.y_slice
            # if direc == "z":
            #     slice = self.z_slice
            # self.info["slice"] = slice[n_slice]
            # self.info["n_slice"] = n_slice

    def compress(
        self,
        n: int,
        direc: str,
        zarr_filepath: str = None,
        i_start: int = None,
        i_size: int = None,
        j_start: int = None,
        j_size: int = None,
        k_start: int = None,
        k_size: int = None,
        keys: list = zarr_keys,
        overwrite: bool = False,
    ):
        """
        Compress 2D slice data into zarr format

        Parameters
        ----------
        n : int
            A selected time step for data
        zarr_filepath : str, optional
            Filepath for the output zarr file, by default None
        direc : str, optional
            Slice direction. 'x', 'y', or 'z', by default "x"
        i_start : int, optional
            Start index for the i dimension, by default None
        i_size : int, optional
            Size of the i dimension, by default None
        j_start : int, optional
            Start index for the j dimension, by default None
        j_size : int, optional
            Size of the j dimension, by default None
        k_start : int, optional
            Start index for the k dimension, by default None
        k_size : int, optional
            Size of the k dimension, by default None
        keys : list, optional
            Keys for the variables to be compressed, by default zarr_keys
        overwrite : bool, optional
            Whether to overwrite the existing zarr file, by default False
        """

        postfixes = self._get_postfixes()

        if zarr_filepath is None:
            zarr_filepath = self._get_filepath_slice_zarr(n, direc)
        else:
            zarr_filepath = Path(zarr_filepath)

        if not overwrite and zarr_filepath.exists():
            print(
                f"File {zarr_filepath} already exists. Set overwrite=True to overwrite it."
            )
            return

        if not direc in ["x", "y", "z"]:
            print(f"Invalid slice direction: {direc}. Must be 'x', 'y', or 'z'.")
            return

        chunk_max = 4096
        if direc == "x":
            i_start, i_size = 0, 1
            if self.geometry == "YinYang":
                chunks3d = (1, min(self.jxg_yy, chunk_max), min(self.kxg_yy, chunk_max))
            else:
                chunks3d = (1, min(self.jx, chunk_max), min(self.kx, chunk_max))
        else:
            i_start = 0 if i_start is None else i_start
            i_size = self.ix if i_size is None else i_size

        if direc == "y":
            j_start, j_size = 0, 1
            chunks3d = (min(self.ix, chunk_max), 1, min(self.kx, chunk_max))
        else:
            j_start = 0 if j_start is None else j_start
            if self.geometry == "YinYang":
                j_size = self.jxg_yy if j_size is None else j_size
            else:
                j_size = self.jx if j_size is None else j_size

        if direc == "z":
            k_start, k_size = 0, 1
            chunks3d = (min(self.ix, chunk_max), min(self.jx, chunk_max), 1)
        else:
            k_start = 0 if k_start is None else k_start
            if self.geometry == "YinYang":
                k_size = self.kxg_yy if k_size is None else k_size
            else:
                k_size = self.kx if k_size is None else k_size

        params_dict = {
            "i_start": i_start,
            "i_size": i_size,
            "j_start": j_start,
            "j_size": j_size,
            "k_start": k_start,
            "k_size": k_size,
        }

        vars_dict = {}
        for key in keys:
            for postfix in postfixes:
                if direc == "x":
                    vars_dict[key + postfix] = np.empty(
                        (len(self.x_slice), j_size, k_size), dtype=np.float32
                    )
                elif direc == "y":
                    vars_dict[key + postfix] = np.empty(
                        (i_size, len(self.y_slice), k_size), dtype=np.float32
                    )
                elif direc == "z":
                    vars_dict[key + postfix] = np.empty(
                        (i_size, j_size, len(self.z_slice)), dtype=np.float32
                    )

        xyz_slice = self.xyz_slice_select(direc)
        for n_slice in range(len(xyz_slice)):
            self.read(n_slice=n_slice, direc=direc, n=n)
            for key in keys:
                for postfix in postfixes:
                    if direc == "x":
                        vars_dict[key + postfix][n_slice, :, :] = self.__dict__[
                            key + postfix
                        ]
                    elif direc == "y":
                        vars_dict[key + postfix][:, n_slice, :] = self.__dict__[
                            key + postfix
                        ]
                    elif direc == "z":
                        vars_dict[key + postfix][:, :, n_slice] = self.__dict__[
                            key + postfix
                        ]

        pyR2D2.zarr_util.save(
            zarr_filepath, vars_dict, params_dict, chunks3d=chunks3d, mode="w"
        )

    def check(
        self,
        n: int,
        direc: str,
        keys: list = zarr_keys,
    ):

        postfixes = self._get_postfixes()
        zarr_filepath = self._get_filepath_slice_zarr(n, direc)
        if not zarr_filepath.exists():
            print(
                f"Zarr file does not exist at n={n} and direc={direc}. Please run Slice.compress() to create it."
            )
            return False

        for key in keys:
            for postfix in postfixes:
                if not key + postfix in pyR2D2.zarr_util.list_vars(zarr_filepath):
                    print(
                        f"Variable {key + postfix} does not exist in zarr file at n={n} and direc={direc}. Please run Slice.compress() to create it."
                    )
                    return False

        xyz_slice = self.xyz_slice_select(direc)
        for postfix in postfixes:
            for n_slice in range(len(xyz_slice)):
                filepath = self._get_filepath_slice(n_slice, direc, n, postfix)
                if not filepath.exists():
                    print(f"File {filepath} does not exist. Anyway you can delete it.")
                    return True

        for n_slice in range(len(xyz_slice)):
            self.read(n_slice=n_slice, direc=direc, n=n, zarr_flag=True)
            qq_copy = {}
            for key in keys:
                for postfix in postfixes:
                    qq_copy[key + postfix] = self.__dict__[key + postfix].copy()

            self.read(n_slice=n_slice, direc=direc, n=n, zarr_flag=False)
            for key in keys:
                for postfix in postfixes:
                    if not np.allclose(
                        qq_copy[key + postfix],
                        self.__dict__[key + postfix],
                        rtol=1e-6,
                        atol=1e-12,
                    ):
                        print(
                            f"Data mismatch for {key + postfix} at slice {n_slice} and direc {direc}."
                        )
                        return False
        print(f"Check passed at n={n} and direc={direc}.")
        return True

    def delete(
        self,
        n: int,
        direc: str,
        force: bool = False,
    ):

        postfixes = self._get_postfixes()
        if force:
            check_flag = True
        else:
            check_flag = self.check(n=n, direc=direc)

        if check_flag:
            xyz_slice = self.xyz_slice_select(direc)
            for n_slice in range(len(xyz_slice)):
                for postfix in postfixes:
                    filepath = self._get_filepath_slice(
                        n_slice, direc, n, postfix=postfix
                    )
                    if filepath.exists():
                        print(f"Deleting {filepath}")
                        filepath.unlink()

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
        if self.ro is not None:
            memflag = self.ro.shape != (self.ix, self.jx)
        if self.ro is None or memflag:
            for key in self.value_keys:
                self.__dict__[key] = np.zeros((self.ix, self.jx))

        dtype = np.dtype(
            [("qq", self.endian + str((self.mtype + 5) * self.ix * self.jx) + "f")]
        )
        with open(self.datadir / "remap" / "qq" / f"qq.dac.{n:08d}", "rb") as f:
            qq = np.fromfile(f, dtype=dtype, count=1)

        for key, m in zip(self.value_keys, range(self.mtype)):
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

        with open(self.datadir.parent / "input_data" / "params.txt", "r") as f:
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
        with open(self.datadir.parent / "input_data" / "value_cart.dac", "rb") as f:
            qq = np.fromfile(f, dtype=dtype, count=1)

        for key in qq.dtype.names:
            if qq[key].size == self.ix:
                self.__dict__[key] = qq[key].reshape((self.ix), order="F")


class _BasePrevAftr(_BaseReader):
    """
    Base class for previous and after time step data
    """

    prev_aftr_kind = ["ro", "vx", "vy", "vz", "bx", "by", "bz", "se"]

    def __init__(self, data):
        self.data = data
        if not hasattr(self.data, "ib_rte_bot"):
            self.ib_rte_bot = 0
        self.ix_prev_aftr = (self.ix0 - self.ib_rte_bot) * self.nx

    def _allocate_prev_aftr_qq(self, dtype):
        """
        Allocate memory for previous and after time step data

        Parameters
        ----------
        dtype : data type
            Data type for allocation
        """
        memflag = True
        if self.ro is not None:
            memflag = not self.ro.shape == (self.ix_prev_aftr, self.jx, self.kx)
        if self.ro is None or memflag:
            for key in self.prev_aftr_kind:
                self.__dict__[key] = np.zeros(
                    (self.ix_prev_aftr, self.jx, self.kx), dtype=dtype
                )

    def _dtype_prev_aftr_qq(self, kind):
        """
        Data type for previous and after time step data

        Parameters
        ----------
        kind : str
            Data type kind for dtype
        """
        return np.dtype(
            [
                (
                    "qq",
                    self.endian
                    + str(self.nx * self.ny * self.nz * len(self.prev_aftr_kind))
                    + kind,
                )
            ]
        )

    def _get_filepath_prev_aftr_qq(self, n, n_prev_aftr, np0, prev_aftr):
        """
        Filepath for previous and after time step data

        Parameters
        ----------
        n : int
            A selected time step for data
        n_prev_aftr : int
            A selected previous or after time step for data
        np0 : int
            A selected MPI process number
        """
        cnou = "{0:05d}".format(np0 // 1000)
        cno = "{0:08d}".format(np0)

        return (
            self.datadir
            / prev_aftr
            / cnou
            / cno
            / f"qq.dac.{n:08d}.{n_prev_aftr:08d}.{cno}"
        )

    def _get_filepath_prev_aftr_zarr(self, n, n_prev_aftr, prev_aftr):
        """
        Filepath for previous and after time step data in zarr format

        Parameters
        ----------
        n : int
            A selected time step for data
        n_prev_aftr : int
            A selected previous or after time step for data
        """
        return self.datadir / prev_aftr / "zarr" / f"qq.{n:08d}.{n_prev_aftr:08d}.zarr"

    def _ijk_start_size(
        self,
        i_start: int = None,
        i_size: int = None,
        j_start: int = None,
        j_size: int = None,
        k_start: int = None,
        k_size: int = None,
    ):
        """
        Get start index and size in each direction for previous and after time step data

        Parameters
        ----------
        i_start : int, optional
            Starting index in x-direction
        i_size : int, optional
            Size in x-direction
        j_start : int, optional
            Starting index in y-direction
        j_size : int, optional
            Size in y-direction
        k_start : int, optional
            Starting index in z-direction
        k_size : int, optional
            Size in z-direction

        Returns
        -------
        i_start, i_size, j_start, j_size, k_start, k_size : int
            Start index and size in each direction
        """

        if i_start is None:
            i_start = 0
        if i_size is None:
            i_size = self.ix_prev_aftr - i_start
        if j_start is None:
            j_start = 0
        if j_size is None:
            j_size = self.jx - j_start
        if k_start is None:
            k_start = 0
        if k_size is None:
            k_size = self.kx - k_start

        return i_start, i_size, j_start, j_size, k_start, k_size

    def _read(
        self,
        n: int,
        n_prev_aftr: int,
        zarr_flag: bool = False,
    ):
        """
        Core function to read previous and after time step data

        Parameters
        ----------
        n : int
            A selected time step for data
        n_prev_aftr : int
            A selected previous or after time step for data
        Core function to read previous and after time step data

        Parameters
        ----------
        n : int
            A selected time step for data
        n_prev_aftr : int
            A selected previous or after time step for data
        """

        for np0 in range(self.npe):
            ib, jb, kb = self.xyz[np0]
            if ib == self.ib_rte_bot:
                filepath = self._get_filepath_prev_aftr_qq(
                    n, n_prev_aftr, np0, prev_aftr=self.prev_aftr
                )

                if not filepath.exists():
                    zarr_flag = True
                    break

        if zarr_flag:
            zarr_filepath = self._get_filepath_prev_aftr_zarr(
                n, n_prev_aftr, prev_aftr=self.prev_aftr
            )
            (qq, params) = pyR2D2.zarr_util.load(zarr_filepath, with_attrs=True)

            for key in self.prev_aftr_kind:
                self.__dict__[key] = qq[key]

            self.x_prev_aftr = qq["x"]
            self.y_prev_aftr = qq["y"]
            self.z_prev_aftr = qq["z"]

            for key in params.keys():
                self.__dict__[key] = params[key]
            return

        else:

            self.x_prev_aftr = self.data.x[-self.ix_prev_aftr :]
            self.y_prev_aftr = self.data.y
            self.z_prev_aftr = self.data.z

            self._allocate_prev_aftr_qq(dtype=np.float64)
            dtype = self._dtype_prev_aftr_qq(kind="d")

            for np0 in range(self.npe):
                ib, jb, kb = self.xyz[np0]
                if ib >= self.ib_rte_bot:
                    ibt = ib - self.ib_rte_bot

                    filepath = self._get_filepath_prev_aftr_qq(
                        n, n_prev_aftr, np0, prev_aftr=self.prev_aftr
                    )
                    qqq_mem = np.memmap(
                        filename=filepath,
                        dtype=dtype,
                        mode="r",
                        shape=(1,),
                    )
                    qqq = qqq_mem[0]
                    for key, m in zip(
                        self.prev_aftr_kind, range(len(self.prev_aftr_kind))
                    ):
                        self.__dict__[key][
                            ibt * self.nx : (ibt + 1) * self.nx,
                            jb * self.ny : (jb + 1) * self.ny,
                            kb * self.nz : (kb + 1) * self.nz,
                        ] = qqq["qq"].reshape(
                            (self.nx, self.ny, self.nz, len(self.prev_aftr_kind)),
                            order="F",
                        )[
                            :, :, :, m
                        ]

    def _compress(
        self,
        n: int,
        n_prev_aftr: int,
        zarr_filepath: str = None,
        i_start: int = None,
        i_size: int = None,
        j_start: int = None,
        j_size: int = None,
        k_start: int = None,
        k_size: int = None,
        overwrite: bool = False,
    ):
        """
        Core function to compress previous and after time step data into zarr format

        Parameters
        ----------
        n : int
            A selected time step for data
        n_prev_aftr : int
            A selected previous or after time step for data
        """

        if zarr_filepath is None:
            zarr_filepath = self._get_filepath_prev_aftr_zarr(
                n, n_prev_aftr, prev_aftr=self.prev_aftr
            )
        else:
            zarr_filepath = Path(zarr_filepath)

        if not overwrite and zarr_filepath.exists():
            print(
                "zarr file ",
                zarr_filepath,
                " already exists. Set overwrite=True to overwrite it.",
            )
            return
        else:

            i_start, i_size, j_start, j_size, k_start, k_size = self._ijk_start_size(
                i_start=i_start,
                i_size=i_size,
                j_start=j_start,
                j_size=j_size,
                k_start=k_start,
                k_size=k_size,
            )

            params_dict = {
                "i_start": i_start,
                "i_size": i_size,
                "j_start": j_start,
                "j_size": j_size,
                "k_start": k_start,
                "k_size": k_size,
            }

            # check if original files exist
            for np0 in range(self.npe):
                ib, jb, kb = self.xyz[np0]
                if ib == self.ib_rte_bot:
                    filepath = self._get_filepath_prev_aftr_qq(
                        n, n_prev_aftr, np0, prev_aftr=self.prev_aftr
                    )

                    if not filepath.exists():
                        print(
                            "Original file does not exist at n=",
                            n,
                            ", ",
                            self.prev_aftr,
                            "=",
                            n_prev_aftr,
                            ". Compression failed.",
                        )
                        return
            self._read(n, n_prev_aftr, zarr_flag=False)

            vars_dict = {}
            for key in self.prev_aftr_kind:
                vars_dict[key] = self.__dict__[key][
                    i_start : i_start + i_size,
                    j_start : j_start + j_size,
                    k_start : k_start + k_size,
                ]

            vars_dict["x"] = self.x_prev_aftr[i_start : i_start + i_size]
            vars_dict["y"] = self.y_prev_aftr[j_start : j_start + j_size]
            vars_dict["z"] = self.z_prev_aftr[k_start : k_start + k_size]

            pyR2D2.zarr_util.save(zarr_filepath, vars_dict, params_dict)

    def check(self, n: int, n_prev_aftr: int):
        """
        Core function to check previous and after time step data files

        Parameters
        ----------
        n : int
            A selected time step for data
        n_prev_aftr : int
            A selected previous or after time step for data
        """

        zarr_filepath = self._get_filepath_prev_aftr_zarr(
            n, n_prev_aftr, prev_aftr=self.prev_aftr
        )
        if not zarr_filepath.exists():
            print(f"Zarr file does not exist at n={n}, {self.prev_aftr}={n_prev_aftr}")
            return False

        for np0 in range(self.npe):
            ib, jb, kb = self.xyz[np0]
            if ib == self.ib_rte_bot:
                filepath = self._get_filepath_prev_aftr_qq(
                    n, n_prev_aftr, np0, prev_aftr=self.prev_aftr
                )

                if not filepath.exists():
                    return True

        self._read(n, n_prev_aftr, zarr_flag=False)

        qq_copy = {}
        for key in self.prev_aftr_kind:
            qq_copy[key] = np.copy(self.__dict__[key])

        self._read(n, n_prev_aftr, zarr_flag=True)
        for key in self.prev_aftr_kind:
            qq_copy[key] = qq_copy[key][
                self.i_start : self.i_start + self.i_size,
                self.j_start : self.j_start + self.j_size,
                self.k_start : self.k_start + self.k_size,
            ]

        for key in self.prev_aftr_kind:
            if (abs(qq_copy[key] - self.__dict__[key])).sum() / (
                abs(qq_copy[key]).sum() + 1e-12
            ) > 1.0e-6:
                print(
                    f"Check failed for {key} at n={n}, {self.prev_aftr}={n_prev_aftr}"
                )
                return False
        print(f"Check passed at n={n}, {self.prev_aftr}={n_prev_aftr}")
        return True

    def _delete(self, n: int, n_prev_aftr: int):
        """
        Core function to delete previous and after time step data files

        Parameters
        ----------
        n : int
            A selected time step for data
        n_prev_aftr : int
            A selected previous or after time step for data
        """

        for np0 in range(self.npe):
            ib, jb, kb = self.xyz[np0]
            if ib >= self.ib_rte_bot:
                filepath = self._get_filepath_prev_aftr_qq(
                    n, n_prev_aftr, np0, prev_aftr=self.prev_aftr
                )
                if filepath.exists():
                    filepath.unlink()


class Previous(_BasePrevAftr):
    """
    Class for previous and after time step data
    """

    prev_aftr = "prev"

    def read(self, n: int, n_prev: int, zarr_flag: bool = False):
        """
        Reads previous time step data

        Parameters
        ----------
        n : int
            A selected time step for data
        n_prev : int
            A selected previous time step for data
        """

        self._read(n, n_prev, zarr_flag=zarr_flag)

    def compress(
        self,
        n: int,
        n_prev: int,
        zarr_filepath: str = None,
        i_start: int = None,
        i_size: int = None,
        j_start: int = None,
        j_size: int = None,
        k_start: int = None,
        k_size: int = None,
        overwrite: bool = False,
    ):
        """
        Compress previous time step data into zarr format

        Parameters
        ----------
        n : int
            A selected time step for data
        n_prev : int
            A selected previous time step for data
        zarr_filepath : str
            Filepath for zarr file
        i_start : int, optional
            Starting index in x-direction for compression
        i_size : int, optional
            Size in x-direction for compression
        j_start : int, optional
            Starting index in y-direction for compression
        j_size : int, optional
            Size in y-direction for compression
        k_start : int, optional
            Starting index in z-direction for compression
        """

        self._compress(
            n,
            n_prev,
            zarr_filepath=zarr_filepath,
            i_start=i_start,
            i_size=i_size,
            j_start=j_start,
            j_size=j_size,
            k_start=k_start,
            k_size=k_size,
            overwrite=overwrite,
        )

    def delete(self, n: int, n_prev: int):
        """
        Delete previous time step data files

        Parameters
        ----------
        n : int
            A selected time step for data
        n_prev : int
            A selected previous time step for data
        """

        self._delete(n, n_prev)


class After(_BasePrevAftr):
    prev_aftr = "aftr"

    def read(self, n: int, n_aftr: int, zarr_flag: bool = False):
        """
        Reads after time step data

        Parameters
        ----------
        n : int
            A selected time step for data
        n_aftr : int
            A selected after time step for data
        """

        self._read(n, n_aftr, zarr_flag=zarr_flag)

    def compress(
        self,
        n: int,
        n_aftr: int,
        zarr_filepath: str = None,
        i_start: int = None,
        i_size: int = None,
        j_start: int = None,
        j_size: int = None,
        k_start: int = None,
        k_size: int = None,
        overwrite: bool = False,
    ):
        """
        Compress after time step data into zarr format

        Parameters
        ----------
        n : int
            A selected time step for data
        n_aftr : int
            A selected after time step for data
        zarr_filepath : str
            Filepath for zarr file
        i_start : int, optional
            Starting index in x-direction for compression
        i_size : int, optional
            Size in x-direction for compression
        j_start : int, optional
            Starting index in y-direction for compression
        j_size : int, optional
            Size in y-direction for compression
        k_start : int, optional
            Starting index in z-direction for compression
        """

        self._compress(
            n,
            n_aftr,
            zarr_filepath=zarr_filepath,
            i_start=i_start,
            i_size=i_size,
            j_start=j_start,
            j_size=j_size,
            k_start=k_start,
            k_size=k_size,
            overwrite=overwrite,
        )

    def delete(self, n: int, n_aftr: int):
        """
        Delete after time step data files

        Parameters
        ----------
        n : int
            A selected time step for data
        n_aftr : int
            A selected after time step for data
        """

        self._delete(n, n_aftr)
