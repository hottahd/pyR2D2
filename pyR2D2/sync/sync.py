import os
import re
import subprocess
from pathlib import Path

import numpy as np

try:  # Python 3.11+
    import tomllib
except ImportError:  # pragma: no cover
    tomllib = None

# R2D2plus がランの記録として data/param/ に写すソース一式と実行ファイル (DEC-638)。
# 解析には使わないので、どのメソッドでも持ってこない。
_PROVENANCE_EXCLUDES = [
    "--exclude=data/param/source",
    "--exclude=data/param/bin",
]


class Sync:
    """
    Class for downloading data from remote server

    Both the legacy output (Fortran R2D2 and R2D2plus ``format = "legacy"``) and
    the R2D2plus compressed output (``data/param/format.toml``, DEC-576) are
    supported. For the compressed output, the file names are taken from the
    local ``param/format.toml``, which :py:meth:`setup` downloads.
    The copies of the source and the executable that R2D2plus keeps in
    ``data/param/source`` and ``data/param/bin`` are never downloaded.
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

    def _caseid(self):
        return self.data.datadir.parts[-2]

    def __getattr__(self, name):
        if hasattr(self.data, name):
            return getattr(self.data, name)

        raise AttributeError(
            f"'{type(self).__name__}' object has no attribute '{name}'"
        )

    def _is_compressed(self):
        """R2D2plus の新出力形式 (``param/format.toml`` が ``format = "compressed"``) か。"""
        return getattr(self.data.p, "output_format", "legacy") == "compressed"

    def _format_pattern(self, kind, **fields):
        """
        ``format.toml`` のパスの雛形 ``kind`` を埋めて ``data/`` からの相対パスにする。
        値が ``None`` の欄は rsync (リモートのシェル) のワイルドカード ``*`` にする。
        """
        pattern = self.format_info[kind]
        if "rank" in fields and fields["rank"] is not None:
            fields.setdefault("group", fields["rank"] // 1000)
        for key in [k for k, v in fields.items() if v is None]:
            pattern = re.sub(r"\{" + key + r"(:[^}]*)?\}", "*", pattern)
            del fields[key]
        return pattern.format(**fields)

    @staticmethod
    def _remote_data(project, caseid):
        return Path("work") / project / "run" / caseid / "data"

    @staticmethod
    def _remote_exists(server, path, ssh="ssh"):
        ssh_result = subprocess.run(
            [ssh, server, "ls " + str(path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        return ssh_result.returncode == 0

    @staticmethod
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
        server, caseid, ssh="ssh", project=None, dist="../run/"
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

        if project is None:
            project = os.getcwd().split("/")[-2]

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
            # R2D2plus: restart と光球スラブ (DEC-641) は大きいので個別に持ってくる
            "--exclude=data/restart",
            "--exclude=data/photosphere",
            *_PROVENANCE_EXCLUDES,
            "--exclude=output.*",
            "-e",
            ssh,
            server + ":work/" + project + "/run/" + caseid + "/",
            dist + caseid + "/",
        ]

        result = Sync.rsync_subprocess_wrapper(args)

    def tau(self, server, n: int = None, ssh="ssh", project=None):
        """
        Downloads data at constant optical depth

        Parameters
        ----------
        server : str
            Name of remote server
        n : int
            Step to be downloaded. If None, all the steps (and data/time/) are downloaded
        ssh : str
            Type of ssh command
        project : str
            Name of project such as 'R2D2'
        """

        if project is None:
            project = os.getcwd().split("/")[-2]

        caseid = self._caseid()
        remote_data = self._remote_data(project, caseid)

        if n is None:
            # data/ 以下から tau/ と time/ だけを持ってくる
            args = [
                "--include=/tau/",
                "--include=/tau/**",
                "--exclude=/time/mhd",
                "--include=/time/",
                "--include=/time/**",
                "--exclude=*",
                "-e",
                ssh,
                server + ":" + str(remote_data) + "/",
                str(self.datadir) + "/",
            ]
        else:
            if self._is_compressed():
                filename = self._format_pattern("tau", nd_tau=n)
            else:
                filename = f"tau/qq.dac.{n:08d}"
            (self.datadir / filename).parent.mkdir(parents=True, exist_ok=True)
            args = [
                "-e",
                ssh,
                server + ":" + str(remote_data / filename),
                str(self.datadir / filename),
            ]

        result = Sync.rsync_subprocess_wrapper(args)

    def _remap_rank_file(self, remote_data, n, ns):
        """ランク ``ns`` の 3D remap のファイル (リモートのパス) とローカルの置き場所。"""
        if self._is_compressed():
            rel = self._format_pattern("remap_qq", rank=int(ns), nd=n)
            return remote_data / rel, self.datadir / Path(rel).parent

        par_dir = f"{int(ns) // 1000:05d}"
        chi_dir = f"{int(ns):08d}"
        nd = "*" if n is None else f"{n:08d}"
        remote_file = (
            remote_data / "remap" / "qq" / par_dir / chi_dir / f"qq.dac.{nd}.{ns}"
        )
        return remote_file, self.datadir / "remap" / "qq" / par_dir / chi_dir

    def _remap_nested(self, server, remote_data, ssh):
        """従来形式の remap/qq がランクごとの小ディレクトリに入っているか。新形式は常に入れ子。"""
        if self._is_compressed():
            return True
        return self._remote_exists(
            server,
            remote_data / "remap/qq/00000/00000000/qq.dac.00000000.00000000",
            ssh,
        )

    def remap_qq(self, server, n, ssh="ssh", project=None):
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

        if project is None:
            project = os.getcwd().split("/")[-2]

        caseid = self._caseid()
        remote_data = self._remote_data(project, caseid)

        # remapを行ったMPIランクの洗い出し
        if self._remap_nested(server, remote_data, ssh):
            nps = np.char.zfill(self.np_ijr.flatten().astype(str), 8)
            for ns in nps:
                remote_file, local_dir = self._remap_rank_file(remote_data, n, ns)
                local_dir.mkdir(parents=True, exist_ok=True)

                args = [
                    "-e",
                    ssh,
                    server + ":" + str(remote_file),
                    str(local_dir),
                ]
                Sync.rsync_subprocess_wrapper(args)
        else:
            print("File does not exist in " + server)

    def xselect(
        self, xs, server, n: int = None, ssh="ssh", project=None
    ):
        """
        Downloads data at certain height

        Parameters
        ----------
            xs : float
                Height to be downloaded
            server : str
                Name of remote server
            n : int
                Step to be downloaded. If None, all the steps are downloaded
            ssh : str
                Type of ssh command
            project : str
                Name of project such as 'R2D2'
        """

        if project is None:
            project = os.getcwd().split("/")[-2]

        i0 = np.argmin(np.abs(self.x - xs))
        ir0 = self.i2ir[i0]

        nps = np.char.zfill(self.np_ijr[ir0 - 1, :].astype(str), 8)
        caseid = self._caseid()
        remote_data = self._remote_data(project, caseid)

        nested = self._remap_nested(server, remote_data, ssh)

        for ns in nps:
            if nested:
                remote_file, local_dir = self._remap_rank_file(remote_data, n, ns)
            else:
                nd = "*" if n is None else f"{n:08d}"
                remote_file = remote_data / "remap" / "qq" / f"qq.dac.{nd}.{ns}"
                local_dir = self.datadir / "remap" / "qq"

            local_dir.mkdir(parents=True, exist_ok=True)
            args = [
                "-e",
                ssh,
                server + ":" + str(remote_file),
                str(local_dir),
            ]
            Sync.rsync_subprocess_wrapper(args)

    def vc(self, server, ssh="ssh", project=None):
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

        if project is None:
            project = os.getcwd().split("/")[-2]

        caseid = self._caseid()
        Sync.setup(server, caseid, ssh=ssh, project=project)
        args = [
            "--exclude=time/mhd",
            "-e",
            ssh,
            server + ":work/" + project + "/run/" + caseid + "/data/remap/vl",
            str(self.datadir / "remap") + "/",
        ]
        Sync.rsync_subprocess_wrapper(args)

    def check(
        self, server, n, ssh="ssh", project=None, end_step=False
    ):
        """
        Downloads checkpoint data

        For R2D2plus runs (``data/restart/`` exists on the remote server), the
        restart slot ``data/restart/<n:08d>/`` is downloaded (only the generation
        that ``meta.toml`` points to). With ``end_step=True``, the newer one of
        the alternating slots ``e`` and ``o`` is downloaded.

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

        if project is None:
            project = os.getcwd().split("/")[-2]

        caseid = self._caseid()
        remote_data = self._remote_data(project, caseid)

        if self._remote_exists(server, remote_data / "restart", ssh):
            self._check_restart(server, n, ssh, remote_data, end_step)
            return

        step = str(n).zfill(8)

        if end_step:
            if np.mod(self.nd, 2) == 0:
                step = "e"
            if np.mod(self.nd, 2) == 1:
                step = "o"

        if self._remote_exists(server, remote_data / "qq" / "00000", ssh):
            io_type = "posixio"
        else:
            io_type = "mpiio"

        if io_type == "posixio":
            for ns in range(self.npe):
                par_dir = f"{int(ns) // 1000:05d}"
                chi_dir = f"{int(ns):08d}"

                remote_base = remote_data / "qq" / par_dir / chi_dir
                remote_file = remote_base / f"qq.dac.{step}.{ns:08d}"
                local_dir = self.datadir / "qq" / par_dir / chi_dir

                local_dir.mkdir(parents=True, exist_ok=True)
                args = [
                    "-e",
                    ssh,
                    server + ":" + str(remote_file),
                    str(self.datadir / "qq" / par_dir / chi_dir),
                ]
                Sync.rsync_subprocess_wrapper(args)
        elif io_type == "mpiio":
            args = [
                "-e",
                ssh,
                server + ":" + str(remote_data / "qq" / ("qq.dac." + step)),
                str(self.datadir / "qq") + "/",
            ]
            Sync.rsync_subprocess_wrapper(args)

    def _fetch_restart_meta(self, server, slot, ssh, remote_data):
        """リモートの ``restart/<slot>/meta.toml`` をローカルへ写して読む。無ければ ``None``。"""
        local_dir = self.datadir / "restart" / slot
        local_dir.mkdir(parents=True, exist_ok=True)
        args = [
            "-e",
            ssh,
            server + ":" + str(remote_data / "restart" / slot / "meta.toml"),
            str(local_dir) + "/",
        ]
        Sync.rsync_subprocess_wrapper(args)
        meta = local_dir / "meta.toml"
        if not meta.exists():
            return None
        if tomllib is None:  # pragma: no cover
            raise ImportError("downloading R2D2plus restart needs Python >= 3.11 (tomllib)")
        with open(meta, "rb") as f:
            return tomllib.load(f)

    def _check_restart(self, server, n, ssh, remote_data, end_step):
        """R2D2plus の restart (``data/restart/<slot>/``, DEC-076) を持ってくる。"""
        if end_step:
            # e と o は交代で書かれるので、step の大きい方が最新
            metas = {}
            for slot in ["e", "o"]:
                meta = self._fetch_restart_meta(server, slot, ssh, remote_data)
                if meta is not None:
                    metas[slot] = meta
            if not metas:
                print("Restart (e/o) does not exist in " + server)
                return
            slot = max(metas, key=lambda s: int(metas[s]["restart"]["step"]))
            meta = metas[slot]
        else:
            slot = f"{n:08d}"
            meta = self._fetch_restart_meta(server, slot, ssh, remote_data)
            if meta is None:
                print(f"Restart {slot} does not exist in " + server)
                return

        remote_slot = remote_data / "restart" / slot
        local_slot = self.datadir / "restart" / slot
        generation = meta["restart"].get("generation")
        if generation:
            # v2: meta.toml が指す世代のディレクトリだけ (古い世代は持ってこない)
            args = [
                "-e",
                ssh,
                server + ":" + str(remote_slot / generation) + "/",
                str(local_slot / generation) + "/",
            ]
        else:
            # v1: slot の直下の state.<rank>.bin
            args = [
                "--exclude=.generation-*",
                "-e",
                ssh,
                server + ":" + str(remote_slot) + "/",
                str(local_slot) + "/",
            ]
        Sync.rsync_subprocess_wrapper(args)

    def slice(
        self,
        server,
        n=None,
        direc="*",
        n_slice=None,
        ssh="ssh",
        project=None,
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

        if project is None:
            project = os.getcwd().split("/")[-2]

        caseid = self._caseid()
        remote_data = self._remote_data(project, caseid)
        (self.datadir / "slice").mkdir(parents=True, exist_ok=True)

        # 断面の位置 (slice.dac) と個数 (params.dac)
        args = [
            "-e",
            ssh,
            server + ":" + str(remote_data / "slice" / "slice.dac"),
            server + ":" + str(remote_data / "slice" / "params.dac"),
            str(self.datadir / "slice") + "/",
        ]
        Sync.rsync_subprocess_wrapper(args)

        if self._is_compressed():
            # YinYang の _yin/_yan も拾うよう、方向の後ろにも * を付ける
            filename = self._format_pattern(
                "slice", name="qq" + direc + "*", nd_tau=n, n=n_slice
            )
        else:
            step = "*" if n is None else str(n).zfill(8)
            n_slice_str = "*" if n_slice is None else str(n_slice).zfill(8)
            filename = "slice/qq" + direc + ".dac." + step + "." + n_slice_str

        args = [
            "-e",
            ssh,
            server + ":" + str(remote_data / filename),
            str(self.datadir / "slice") + "/",
        ]
        Sync.rsync_subprocess_wrapper(args)

    def all(
        self, server, ssh="ssh", project=None, dist="../run/"
    ):
        """
        This method downloads all the data

        The copies of the source and the executable that R2D2plus keeps in
        ``data/param/source`` and ``data/param/bin`` are not downloaded.

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

        if project is None:
            project = os.getcwd().split("/")[-2]

        caseid = self._caseid()
        args = [
            *_PROVENANCE_EXCLUDES,
            "-e",
            ssh,
            server + ":work/" + project + "/run/" + caseid + "/",
            dist + caseid,
        ]
        Sync.rsync_subprocess_wrapper(args)
