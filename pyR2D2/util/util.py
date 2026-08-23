import pyR2D2


def init(
    main_locals,
    instance_name="d",
    verbose=True,
    datadir=None,
    google=False,
    copy_to_local=False,
):
    """
    Initialize pyR2D2.Data instance in main program.

    Parameters
    ----------
    main_locals : dict
        locals() in main program, which include local variables
    instance_name : str
        name of instance of pyR2D2.Data object
    verbose : bool
        If True, print self.summary()
    datadir : str
        directory path to pyR2D2 data


    Notes
    -----
    1. Select caseid from existing caseid in main_locals or from user input.
    2. Initialize instance of pyR2D2.Data object in main program.
    3. Assign the pyR2D2.Data.p dictionary keys to main_locals.

    Examples
    --------
    .. code-block:: python

        import pyR2D2
        pyR2D2.util.init(locals())

    """
    import pyR2D2

    if datadir is None:
        # if google is True, force to read caseid from user input.
        main_locals["caseid"] = caseid_select(main_locals, force_read=google)
        datadir = "../run/" + main_locals["caseid"] + "/data/"

    initialize_instance(main_locals, instance_name)

    print(
        "### An instance of pyR2D2.Data class is initialized as "
        + instance_name
        + " ###"
    )
    main_locals[instance_name] = pyR2D2.Data(datadir)

    if copy_to_local:
        print("### Copying pyR2D2.Data.p to locals() ###")
        locals_define(main_locals[instance_name], main_locals)

    print("")
    if google:
        import pyR2D2.write.google

        pyR2D2.write.google.set_top_line()
        pyR2D2.write.google.set_cells_gspread(main_locals[instance_name])

    if verbose:
        main_locals[instance_name].summary()


def initialize_instance(main_locals, instance_name):
    """
    Initializes arbitrary instance of pyR2D2.Data object in main program.

    Parameters
    ----------
    main_locals : dict
        locals() in main program, which include local variables
    instance_name : str
        name of instance of pyR2D2.Data object
    """
    if not instance_name in main_locals:
        main_locals[instance_name] = None


def caseid_select(main_locals, force_read=False):
    """
    Choose caseid from input or from user input.

    Parameters
    ----------
    main_locals : dict
        locals() in main program, which include local variables
    force_read : bool
        If True, force to read caseid from user input.
        If False, read caseid from main_locals if exists

    Returns
    -------
    caseid : str
        caseid
    """

    RED = "\033[31m"
    END = "\033[0m"

    if force_read or not "caseid" in main_locals:
        caseid = "d" + str(input(RED + "input caseid id (3 digit): " + END)).zfill(3)
    else:
        caseid = main_locals["caseid"]

    return caseid


def locals_define(data, main_locals):
    """
    Substitute selp.p to main_locals in main program.

    Parameters
    ----------
    data : pyR2D2.Data, or, pyR2D2.Read
        Instance of pyR2D2.Data or pyR2D2.Read classes
    main_locals : dict
        locals() in main program, which include local variables
    """
    for key, value in data.p.__dict__.items():
        main_locals[key] = value
    return


def define_n0(data, main_locals, nd_type="nd"):
    """
    Define n0 in main_locals if not exists.

    Parameters
    ----------
    data : pyR2D2.Data, or, pyR2D2.Read
        instance of pyR2D2.Data or pyR2D2.Read classes
    main_locals : dict
        locals() in main program, which include local variables
    nd_type : str
        type of nd. 'nd' for MHD output, 'nd_tau' for high cadence output.
    """
    if "n0" not in main_locals:
        main_locals["n0"] = 0
    if main_locals["n0"] > data.p.__dict__[nd_type]:
        main_locals["n0"] = data.p.__dict__[nd_type]

    return main_locals["n0"]


def get_best_unit(size, unit_multipliers):
    """
    Get the best unit for displaying the size.

    Parameters
    ----------
    size : int
        size of the file
    unit_multipliers : dict
        dictionary of unit multipliers
    """
    for unit in reversed(["B", "kB", "MB", "GB", "TB", "PB"]):
        if size >= unit_multipliers[unit]:
            return unit
    return "B"


def get_total_file_size(directory, unit=None):
    """
    Evaluate total size of files in directory.

    Parameters:
       directory (str): directory path
       unit (str): unit of file size. Choose from 'B', 'kB', 'MB', 'GB', 'TB', 'PB'.
    Returns:
       total_size (int): total size of files in directory in bytes
    """
    import os

    unit_multipliers = {
        "B": 1,
        "kB": 1024,
        "MB": 1024**2,
        "GB": 1024**3,
        "TB": 1024**4,
        "PB": 1024**5,
    }

    if unit is not None and unit not in unit_multipliers:
        raise ValueError(
            f"Invalid unit: {unit}. Choose from 'B', 'kB', 'MB', 'GB', 'TB', 'PB'."
        )

    total_size = 0
    for dirpath, dirnames, filenames in os.walk(directory):
        for filename in filenames:
            filepath = os.path.join(dirpath, filename)
            try:
                total_size += os.path.getsize(filepath)

                if unit is None:
                    display_unit = get_best_unit(total_size, unit_multipliers)
                else:
                    display_unit = unit
                converted_size = total_size / unit_multipliers[display_unit]
                print(
                    f"\r{' '*50}\rCurrent total size: {converted_size:.2f} {display_unit}",
                    end="",
                    flush=True,
                )
                # print(f"\rCurrent total size: {converted_size:.2f} {display_unit} ", end='', flush=True)
            except FileNotFoundError:
                # ファイルが見つからない場合は無視
                continue

    if unit is None:
        unit = get_best_unit(total_size, unit_multipliers)

    final_size = total_size / unit_multipliers[unit]
    print(f"\nFinal total size: {final_size:.2f} {unit}")

    return final_size, unit


def update_results_file(file_path, total_size, unit, caseid, dir_path):
    """
    Updates the results file with the size of a directory.

    Parameters
    ----------
    file_path : str
        Path to the results file
    total_size : float
        Total size of the directory
    unit : str
        Unit of the file size
    caseid : str
        The case ID of the directory
    dir_path : str
        The directory path

    """
    import os
    from datetime import datetime

    # 現在の日時を取得
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # シンボリックリンクの場合の実体パスを取得
    if os.path.islink(dir_path):
        real_path = os.path.abspath(os.readlink(dir_path))
    else:
        real_path = os.path.abspath(dir_path)

    # 既存のデータを読み込む
    if os.path.exists(file_path):
        with open(file_path, "r") as f:
            existing_data = f.readlines()
    else:
        existing_data = []

    # データを辞書に変換
    existing_dict = {}
    for line in existing_data:
        parts = line.strip().rsplit(maxsplit=3)
        if len(parts) == 4:
            existing_caseid = parts[0]
            size = parts[1]
            timestamp = parts[2]
            path = parts[3]
            existing_dict[existing_caseid] = f"{size} {timestamp} {path}"

    # ディレクトリサイズの情報を更新
    directory_size_str = f"{total_size:6.2f} {unit}"
    existing_dict[caseid] = f"{directory_size_str} {now} {real_path}"

    # 更新されたデータを書き込む
    with open(file_path, "w") as f:
        for caseid in sorted(existing_dict):
            size_timestamp_realpath = existing_dict[caseid]
            # フォーマット: ディレクトリ名、右揃えで6桁、小数点2桁、単位
            f.write(f"{caseid:<10} {size_timestamp_realpath}\n")


def eos(data, ro, se, var):
    """
    Returns the Table equation of state.

    Parameters
    ----------
    data : pyR2D2.Data, or, pyR2D2.Read
        Instance of pyR2D2.Data or pyR2D2.Read classes
    ro : float
        Density
    se : float)
        Entropy
    var : str
        Variable name, pr, te, en, op, dprodro

    Returns
    -------
    qq : float
        Corresponding variable

    Warning
    -------
    This method is very slow for large numpy array. Please consider using fortraun code.
    It also only accepts scalars; use :func:`pyR2D2.util.eos_table` for arrays.

    Warning
    -------
    **This reads the table unconditionally, which is not what R2D2 does.**
    The code switches to the linearized EOS whenever the relative amplitude
    ``ct = max(|ro1|/ro0, |se1|/se0)`` stays below ``ct0 = 3e-3``, and in a
    deep-convection run that is everywhere.  Comparing a plain table lookup
    against the ``pr``/``te`` in the output files can be off by four orders of
    magnitude *of the signal amplitude*, not four decimal places.  Use
    :func:`pyR2D2.util.eos_switch` to reproduce what the code actually did.
    """

    import numpy as np

    iro = int((np.log(ro) - data.log_ro_e[0]) // data.dlogro_e)
    ise = int((se - data.se_e[0]) // data.dse_e)

    dlogro = np.log(ro) - data.log_ro_e[iro]
    dse = se - data.se_e[ise]

    qq = np.exp(
        (
            +data.p.__dict__["log_" + var + "_e"][iro, ise]
            * (data.dlogro_e - dlogro)
            * (data.dse_e - dse)
            + data.p.__dict__["log_" + var + "_e"][iro + 1, ise]
            * (dlogro)
            * (data.dse_e - dse)
            + data.p.__dict__["log_" + var + "_e"][iro, ise + 1]
            * (data.dlogro_e - dlogro)
            * (dse)
            + data.p.__dict__["log_" + var + "_e"][iro + 1, ise + 1] * (dlogro) * (dse)
        )
        / data.dlogro_e
        / data.dse_e
    )
    return qq


def eos_table(data, ro, se, var):
    """
    Returns the table equation of state, vectorized over numpy arrays.

    This is the array-capable counterpart of :func:`pyR2D2.util.eos`, which
    only accepts scalars (it calls ``int()`` on the index).  The interpolation
    itself is identical: bilinear in ``(log rho, s)`` on the logarithm of the
    tabulated quantity, exponentiated afterwards.  Everything is done in
    float64, whereas ``pyR2D2.Data.eos.eval`` (the C++ helper) works in
    float32.

    Parameters
    ----------
    data : pyR2D2.Data, or, pyR2D2.Read
        Instance of pyR2D2.Data or pyR2D2.Read classes
    ro : float or numpy.ndarray
        Total density (background + perturbation)
    se : float or numpy.ndarray
        Total entropy (background + perturbation)
    var : str
        Variable name; pr, te, en, op, or dprdro

    Returns
    -------
    qq : numpy.ndarray
        Corresponding variable, same shape as ``ro``
    """

    import numpy as np

    log_table = np.asarray(data.p.__dict__["log_" + var + "_e"], dtype=np.float64)
    log_ro_e = np.asarray(data.log_ro_e, dtype=np.float64)
    se_e = np.asarray(data.se_e, dtype=np.float64)
    dlogro_e = float(data.dlogro_e)
    dse_e = float(data.dse_e)

    log_ro = np.log(np.asarray(ro, dtype=np.float64))
    se = np.asarray(se, dtype=np.float64)

    # R2D2 (eos_proc.F95 の eos_calc_init) と同じく、両端で 2x2 ステンシルが
    # 表からはみ出さないところまでクランプする。範囲外は外挿になる。
    iro = np.clip(((log_ro - log_ro_e[0]) / dlogro_e).astype(int), 0, log_ro_e.size - 2)
    ise = np.clip(((se - se_e[0]) / dse_e).astype(int), 0, se_e.size - 2)

    dlogro = log_ro - log_ro_e[iro]
    dse = se - se_e[ise]

    # fmt: off
    return np.exp(
        (
            + log_table[iro    , ise    ] * (dlogro_e - dlogro) * (dse_e - dse)
            + log_table[iro + 1, ise    ] * (           dlogro) * (dse_e - dse)
            + log_table[iro    , ise + 1] * (dlogro_e - dlogro) * (          dse)
            + log_table[iro + 1, ise + 1] * (           dlogro) * (          dse)
        )
        / dlogro_e
        / dse_e
    )
    # fmt: on


def eos_switch(data, ro1, se1, var, ct0=3.0e-3):
    """
    Returns the equation of state **as R2D2 itself evaluates it**, i.e. with
    the switch between the table EOS and the linearized EOS.

    Why this exists
    ---------------
    :func:`pyR2D2.util.eos`, :func:`pyR2D2.util.eos_table` and
    ``pyR2D2.Data.eos.eval`` all read the table unconditionally.  R2D2 does
    not.  Wherever the code needs pressure, temperature or internal energy it
    computes the relative amplitude

    .. math:: c_t = \\max(|\\rho_1|/\\rho_0, |s_1|/s_0)

    and uses the **linearized** EOS when :math:`c_t \\le c_{t0}`
    (``ct0 = 3e-3``, ``eos_def.F90``), falling back to the table only above
    it.  The switch is hard, not smooth::

        feos = sign(0.5, ct - ct0) + 0.5

    and appears identically in ``runge_kutta_cartesian.F90`` /
    ``runge_kutta_spherical.F90`` (pressure), ``artdif_cartesian.F90`` /
    ``artdif_spherical.F90`` (temperature and enthalpy), ``cfl.F90``
    (``dprdro``, hence the sound speed) and ``remap_calc.F90`` (the pr/te/en
    written into ``remap/vl``).

    In a deep-convection run :math:`c_t` never reaches :math:`c_{t0}`, so the
    code runs entirely on the linearized EOS while a plain table lookup
    returns something else.  That is the discrepancy this function removes.

    Parameters
    ----------
    data : pyR2D2.Data, or, pyR2D2.Read
        Instance of pyR2D2.Data or pyR2D2.Read classes
    ro1 : numpy.ndarray
        Density **perturbation**, i.e. ``d.qq['ro']`` as stored in the output.
        The first axis must be the vertical one (size ``ix``).
    se1 : numpy.ndarray
        Entropy **perturbation**, same shape as ``ro1``
    var : str
        Variable name; pr, te, en, op, or dprdro
    ct0 : float
        Switching threshold, ``eos_def.F90``'s ``ct0``.  Change it only to
        experiment; the code itself compiles it in as 3e-3.

    Returns
    -------
    qq : numpy.ndarray
        For ``pr``, ``te`` and ``en``: the **perturbation**, matching what
        ``remap_calc.F90`` writes into the output files (add ``d.p.pr0`` etc.
        for the total).  For ``op`` and ``dprdro``: the value itself, since
        neither has a linearized perturbation form.  ``op`` is never switched
        (the table is the only source); ``dprdro`` switches to the background
        profile ``d.p.dprdro``.
    """

    import numpy as np

    ro1 = np.asarray(ro1, dtype=np.float64)
    se1 = np.asarray(se1, dtype=np.float64)
    if ro1.shape != se1.shape:
        raise ValueError(
            f"ro1 and se1 must have the same shape, got {ro1.shape} and {se1.shape}"
        )

    def background(name):
        """鉛直1D の背景量を ro1 の形へ broadcast する。"""
        a = np.asarray(data.p.__dict__[name], dtype=np.float64)
        if ro1.ndim == 0:
            raise ValueError(
                "eos_switch needs arrays whose first axis is the vertical one; "
                "for a single point, index the background yourself and use eos_table"
            )
        if a.size != ro1.shape[0]:
            raise ValueError(
                f"the first axis of ro1 (size {ro1.shape[0]}) must be the vertical one "
                f"(background '{name}' has size {a.size})"
            )
        return a.reshape((a.size,) + (1,) * (ro1.ndim - 1))

    ro0 = background("ro0")
    se0 = background("se0")
    table = eos_table(data, ro1 + ro0, se1 + se0, var)

    # 不透明度には線形化した対応物が無く、R2D2 も常に表を引く
    # (rte_multiray.F90, rte_io.F90, remap_calc.F90)。
    if var == "op":
        return table

    # ct = max(|ro1|/ro0, |se1|/se0)、feos は 0 か 1 のどちらか。
    ct = np.maximum(np.abs(ro1) / ro0, np.abs(se1) / se0)
    feos = np.where(ct > ct0, 1.0, 0.0)

    # dprdro の線形側は背景プロファイルそのものである (cfl.F90)。
    if var == "dprdro":
        return table * feos + background("dprdro") * (1.0 - feos)

    linear = background("d" + var + "dro") * ro1 + background("d" + var + "dse") * se1
    return (table - background(var + "0")) * feos + linear * (1.0 - feos)
