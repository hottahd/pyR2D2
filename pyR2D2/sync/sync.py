import os
import subprocess

import numpy as np


class Sync:
    """
    Class for downloading data from remote server
    """

    def __init__(self, data):
        """
        Initialize pyR2D2.Sync

        Parameters
        ----------
        read : pyR2D2.Read
            Instance of pyR2D2.read
        """
        self.data = data

    def __getattr__(self, name):
        if hasattr(self.data, name):
            return getattr(self.data, name)

        raise AttributeError(
            f"'{type(self).__name__}' object has no attribute '{name}'"
        )

    def rsync_subprocess_wrapper(args):
        command = ["rsync", "-avP"] + args
        result = subprocess.run(
            command,
            # stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True,
            shell=False,
        )

        return result

    @staticmethod
    def setup(
        server, caseid, ssh="ssh", project=os.getcwd().split("/")[-2], dist="../run/"
    ):
        """
        Downloads setting data from remote server

        Parameters
        ----------
        server : str
            Name of remote server
        caseid : str
            caseid format of 'd001'
        ssh : str
            Type of ssh command
        project : str
            Name of project such as 'R2D2'
        dist :str
            Destination of data directory
        """

        args = [
            "--exclude=a.out",
            "--exclude=.git",
            "--exclude=make/*.o",
            "--exclude=make/*.lst",
            "--exclude=make/*.mod",
            "--exclude=data/qq",
            "--exclude=data/remap/qq",
            "--exclude=data/remap/vl/vl*",
            "--exclude=data/slice/qq*",
            "--exclude=data/tau/qq*",
            "--exclude=data/aftr",
            "--exclude=data/prev",
            "--exclude=output.*",
            "-e",
            ssh,
            server + ":work/" + project + "/run/" + caseid + "/",
            dist + caseid + "/",
        ]

        result = Sync.rsync_subprocess_wrapper(args)

    def tau(self, server, n: int = None, ssh="ssh", project=os.getcwd().split("/")[-2]):
        """
        Downloads data at constant optical depth

        Parameters
        ----------
        server : str
            Name of remote server
        ssh : str
            Type of ssh command
        project : str
            Name of project such as 'R2D2'
        """

        if n is None:
            filename = "*"
        else:
            filename = "tau/qq.dac." + str(n).zfill(8)

        caseid = self.datadir.split("/")[-3]
        args = [
            "--exclude=param",
            "--exclude=qq",
            "--exclude=remap",
            "--exclude=slice",
            "--exclude=time/mhd",
            "-e",
            ssh,
            server + ":work/" + project + "/run/" + caseid + "/data/" + filename,
            self.datadir + filename,
        ]

        result = Sync.rsync_subprocess_wrapper(args)

    def remap_qq(self, server, n, ssh="ssh", project=os.getcwd().split("/")[-2]):
        """
        Downloads full 3D remap data

        Parameters
        ----------
        n : int
            Target time step
        server : str
            Name of remote server
        project : str
            Name of project such as 'R2D2'
        ssh : str
            Type of ssh command
        """

        caseid = self.datadir.split("/")[-3]

        # check if file exists
        ssh_result = subprocess.run(
            [
                ssh,
                server,
                "ls work/"
                + project
                + "/run/"
                + caseid
                + "/data/remap/qq/00000/00000000/qq.dac."
                + str(n).zfill(8)
                + ".00000000",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        # remapを行ったMPIランクの洗い出し
        if ssh_result.returncode == 0:
            nps = np.char.zfill(self.np_ijr.flatten().astype(str), 8)
            for ns in nps:
                par_dir = str(int(ns) // 1000).zfill(5) + "/"
                chi_dir = str(int(ns)).zfill(8) + "/"

                os.makedirs(
                    self.datadir + "remap/qq/" + par_dir + chi_dir, exist_ok=True
                )
                args = [
                    "-e",
                    ssh,
                    server
                    + ":work/"
                    + project
                    + "/run/"
                    + caseid
                    + "/data/remap/qq/"
                    + par_dir
                    + chi_dir
                    + "qq.dac."
                    + str(n).zfill(8)
                    + "."
                    + ns,
                    self.datadir + "remap/qq/" + par_dir + chi_dir,
                ]
                Sync.rsync_subprocess_wrapper(args)
        else:
            print("File does not exist in " + server)

    def xselect(
        self, xs, server, n: int = None, ssh="ssh", project=os.getcwd().split("/")[-2]
    ):
        """
        Downloads data at certain height

        Parameters
        ----------
            xs : float
                Height to be downloaded
            server : str
                Name of remote server
            ssh : str
                Type of ssh command
            project : str
                Name of project such as 'R2D2'
        """

        i0 = np.argmin(np.abs(self.x - xs))
        ir0 = self.i2ir[i0]

        nps = np.char.zfill(self.np_ijr[ir0 - 1, :].astype(str), 8)
        caseid = self.datadir.split("/")[-3]

        if n is None:
            filename_part = "qq.dac.*."
        else:
            filename_part = "qq.dac." + str(n).zfill(8) + "."

        # check if file exists
        ssh_result = subprocess.run(
            [
                ssh,
                server,
                "ls work/"
                + project
                + "/run/"
                + caseid
                + "/data/remap/qq/00000/00000000/qq.dac."
                + str(n).zfill(8)
                + ".00000000",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        for ns in nps:
            if ssh_result.returncode == 0:
                par_dir = str(int(ns) // 1000).zfill(5) + "/"
                chi_dir = str(int(ns)).zfill(8) + "/"
            else:
                par_dir = ""
                chi_dir = ""

            os.makedirs(self.datadir + "remap/qq/" + par_dir + chi_dir, exist_ok=True)
            args = [
                "-e",
                ssh,
                server
                + ":work/"
                + project
                + "/run/"
                + caseid
                + "/data/remap/qq/"
                + par_dir
                + chi_dir
                + filename_part
                + ns,
                self.datadir + "remap/qq/" + par_dir + chi_dir,
            ]
            Sync.rsync_subprocess_wrapper(args)

    def vc(self, server, ssh="ssh", project=os.getcwd().split("/")[-2]):
        """
        Downloads pre analyzed data

        Parameters
        ----------
        server : str
            Name of remote server
        ssh : str
            Type of ssh command
        project : str
            Name of project such as 'R2D2'
        """

        caseid = self.datadir.split("/")[-3]
        Sync.setup(server, caseid, ssh=ssh, project=project)
        args = [
            "--exclude=time/mhd",
            "-e",
            ssh,
            server + ":work/" + project + "/run/" + caseid + "/data/remap/vl",
            self.datadir + "remap/",
        ]
        Sync.rsync_subprocess_wrapper(args)

    def check(
        self, server, n, ssh="ssh", project=os.getcwd().split("/")[-2], end_step=False
    ):
        """
        Downloads checkpoint data

        Parameters
        ----------
        n : int
            Step to be downloaded
        server : str
            Name of remote server
        ssh : str
            Type of ssh command
        project : str
            Name of project such as 'R2D2'
        end_step : bool
            If true, checkpoint of end step is read
        """

        step = str(n).zfill(8)

        if end_step:
            if np.mod(self.nd, 2) == 0:
                step = "e"
            if np.mod(self.nd, 2) == 1:
                step = "o"

        caseid = self.datadir.split("/")[-3]
        ssh_result = subprocess.run(
            [
                ssh,
                server,
                "ls work/" + project + "/run/" + caseid + "/data/qq/00000",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        if ssh_result.returncode == 0:
            io_type = "posixio"
        else:
            io_type = "mpiio"

        if io_type == "posixio":
            for ns in range(self.npe):
                par_dir = str(int(ns) // 1000).zfill(5) + "/"
                chi_dir = str(int(ns)).zfill(8) + "/"

                os.makedirs(self.datadir + "qq/" + par_dir + chi_dir, exist_ok=True)
                args = [
                    "-e",
                    ssh,
                    server
                    + ":work/"
                    + project
                    + "/run/"
                    + caseid
                    + "/data/qq/"
                    + par_dir
                    + chi_dir
                    + "qq.dac."
                    + step
                    + "."
                    + str(int(ns)).zfill(8),
                    self.datadir + "qq/" + par_dir + chi_dir,
                ]
                Sync.rsync_subprocess_wrapper(args)
        elif io_type == "mpiio":
            args = [
                "-e",
                ssh,
                server
                + ":work/"
                + project
                + "/run/"
                + caseid
                + "/data/qq/qq.dac."
                + step,
                self.datadir + "qq/",
            ]
            Sync.rsync_subprocess_wrapper(args)

    def slice(
        self,
        server,
        n=None,
        direc="*",
        n_slice=None,
        ssh="ssh",
        project=os.getcwd().split("/")[-2],
    ):
        """
        Downloads slice data

        Parameters
        ----------
        server : str
            Name of remote server
        n : int
            Step to be downloaded
        direc : str
            Direction of slice data, 'x', 'y', 'z' or '*'
        n_slice : int
            Slice number
        ssh : str
            Type of ssh command
        project : str
            Name of project such as 'R2D2'
        """

        if n is None:
            step = "*"
        else:
            step = str(n).zfill(8)

        if n_slice is None:
            n_slice_str = "*"
        else:
            n_slice_str = str(n_slice).zfill(8)

        caseid = self.datadir.split("/")[-3]
        args = [
            "-e",
            ssh,
            server + ":work/" + project + "/run/" + caseid + "/data/slice/slice.dac",
            self.datadir + "/slice",
        ]
        Sync.rsync_subprocess_wrapper(args)
        args = [
            "-e",
            ssh,
            server
            + ":work/"
            + project
            + "/run/"
            + caseid
            + "/data/slice/qq"
            + direc
            + ".dac."
            + step
            + "."
            + n_slice_str,
            self.datadir + "/slice",
        ]
        Sync.rsync_subprocess_wrapper(args)

    def all(
        self, server, ssh="ssh", project=os.getcwd().split("/")[-2], dist="../run/"
    ):
        """
        This method downloads all the data

        Parameters
        ----------
        server : str
            Name of remote server
        ssh : str
            Type of ssh command
        project :str
            Name of project such as 'R2D2'
        dist :str
            Destination of data directory
        """

        caseid = self.datadir.split("/")[-3]
        args = [
            "-e",
            ssh,
            server + ":work/" + project + "/run/" + caseid + "/",
            dist + caseid,
        ]
        Sync.rsync_subprocess_wrapper(args)
