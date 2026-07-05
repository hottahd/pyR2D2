import shutil
import zipfile
from pathlib import Path

import numpy as np
import zarr


def open_zarr_group(path: str, use_zip: bool = False):
    path = Path(path)

    # normal .zarr directory
    if path.is_dir() and not use_zip:
        return zarr.open_group(path, mode="r")

    # explicit .zarr.zip file
    if path.is_file() and path.suffix == ".zip":
        store = zarr.storage.ZipStore(path, mode="r")
        return zarr.open_group(store=store, mode="r")

    zip_path = Path(str(path) + ".zip")
    if zip_path.is_file():
        store = zarr.storage.ZipStore(zip_path, mode="r")
        return zarr.open_group(store=store, mode="r")

    raise FileNotFoundError(f"Zarr directory or zip file not found: {path}")


def _check_zarr_zip_equivalent(path, zip_path):
    root_org = open_zarr_group(path)
    root_zip = open_zarr_group(zip_path)

    if root_org.attrs.asdict() != root_zip.attrs.asdict():
        raise RuntimeError(
            "Attributes in the original and zipped zarr do not match. "
            "Aborting removal of original directory."
        )

    if set(root_org.array_keys()) != set(root_zip.array_keys()):
        raise RuntimeError(
            "Array names in the original and zipped zarr do not match. "
            "Aborting removal of original directory."
        )

    for name in root_org.array_keys():
        if root_org[name].shape != root_zip[name].shape:
            raise RuntimeError(
                f"Array shapes for '{name}' do not match. "
                "Aborting removal of original directory."
            )
        if root_org[name].dtype != root_zip[name].dtype:
            raise RuntimeError(
                f"Array dtypes for '{name}' do not match. "
                "Aborting removal of original directory."
            )


def zip_zarr(
    path: str,
    zip_path: str = None,
    overwrite: bool = False,
    remove_original: bool = False,
):
    """
    Convert a .zarr directory to a .zarr.zip file without compression.

    The zip file can be opened by zarr.storage.ZipStore.

    Parameters
    ----------
    path : str
        Path to the .zarr directory.
    zip_path : str, optional
        Path to the output .zarr.zip file. If None, defaults to path + ".zip".
    overwrite : bool, optional
        If True, overwrite the existing zip file. Default is False.
    remove_original : bool, optional
        If True, remove the original .zarr directory after zipping. Default is False.
    """

    path = Path(path)

    if not path.is_dir():
        raise FileNotFoundError(f"Zarr directory not found: {path}")

    if zip_path is None:
        zip_path = Path(str(path) + ".zip")
    else:
        zip_path = Path(zip_path)

    if zip_path.exists():
        if overwrite:
            zip_path.unlink()
        else:
            if remove_original:
                _check_zarr_zip_equivalent(path, zip_path)
                shutil.rmtree(path)
                return zip_path

            raise FileExistsError(f"Zip file already exists: {zip_path}")

    zip_path.parent.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(
        zip_path,
        mode="w",
        compression=zipfile.ZIP_STORED,
        allowZip64=True,
    ) as zf:
        for f in path.rglob("*"):
            if f.is_file():
                zf.write(f, arcname=f.relative_to(path))

    if remove_original:
        _check_zarr_zip_equivalent(path, zip_path)
        shutil.rmtree(path)
    return zip_path


def _check_chunks(chunks, shape_len):
    if chunks is None:
        return

    if not isinstance(chunks, (tuple, list)):
        raise TypeError(f"Chunks must be a tuple or list, got {type(chunks)}")

    if len(chunks) != shape_len:
        raise ValueError(
            f"Chunks length {len(chunks)} does not match array dimensions {shape_len}"
        )
    for chunk in chunks:
        if not isinstance(chunk, int):
            raise TypeError(f"Chunk sizes must be integers, got {type(chunk)}")
        if chunk <= 0:
            raise ValueError(f"Invalid chunk size {chunk}")


def save(
    path: str,
    vars_dict: dict,
    params: dict = None,
    max_chunk_size: int = 512,
    chunks1d: tuple = None,
    chunks2d: tuple = None,
    chunks3d: tuple = None,
    clevel: int = 5,
    mode: str = "w",
):
    """


    Parameters
    ----------
    path : str
        path to save zarr data
    vars_dict : dict
        dictionary of variables to save
    params : dict, optional
        additional parameters to save, by default None
    max_chunk_size : int, optional
        maximum chunk size for zarr arrays, by default 512
    chunks1d : tuple, optional
        chunk size for 1D arrays, by default None (calculated from max_chunk_size)
    chunks2d : tuple, optional
        chunk size for 2D arrays, by default None (calculated from max_chunk_size)
    chunks3d : tuple, optional
        chunk size for 3D arrays, by default None (calculated from max_chunk_size
    clevel : int, optional
        compression level for zarr, by default 5
    mode : str, optional
        file mode for zarr (e.g., "w" for write, "a" for append), by default "w"

    Raises
    ------
    ValueError
        If an array with unsupported dimensions is encountered.
    """
    root = zarr.open_group(path, mode=mode)

    codec_f32 = zarr.codecs.BloscCodec(
        cname="zstd",
        clevel=clevel,
        shuffle=zarr.codecs.BloscShuffle.bitshuffle,
        typesize=4,
    )

    if params is not None:
        root.attrs["params"] = params

    for name, array in vars_dict.items():
        array = np.asarray(array, dtype=np.float32, order="C")
        if array.ndim == 3:
            if chunks3d is None:
                chunks = tuple(min(max_chunk_size, array.shape[i]) for i in range(3))
            else:
                _check_chunks(chunks3d, 3)
                chunks = chunks3d
            root.create_array(
                name,
                data=array,
                chunks=chunks,
                compressors=codec_f32,
                overwrite=True,
            )
        elif array.ndim == 2:
            if chunks2d is None:
                chunks = tuple(min(max_chunk_size, array.shape[i]) for i in range(2))
            else:
                _check_chunks(chunks2d, 2)
                chunks = chunks2d
            root.create_array(
                name,
                data=array,
                chunks=chunks,
                compressors=codec_f32,
                overwrite=True,
            )
        elif array.ndim == 1:
            if chunks1d is None:
                chunks = (min(max_chunk_size, array.shape[0]),)
            else:
                _check_chunks(chunks1d, 1)
                chunks = chunks1d
            root.create_array(
                name,
                data=array,
                chunks=chunks,
                compressors=codec_f32,
                overwrite=True,
            )
        else:
            raise ValueError(
                f"Unsupported array dimension: {array.ndim} for variable {name}"
            )


def load(
    path: str,
    names="all",
    with_attrs: bool = False,
    i0: int = None,
    i1: int = None,
    j0: int = None,
    j1: int = None,
    k0: int = None,
    k1: int = None,
    use_zip: bool = False,
):
    """
    Load zarr data

    Parameters
    ----------
    path : str
        Path to the Zarr directory.
    names : str or list of str
        Variable name(s) to load.
        If "all", all arrays in the Zarr group are loaded.
    with_attrs : bool, optional
        If True, also return Zarr attributes (metadata). By default, False.
    i0, i1, j0, j1, k0, k1 : int, optional
        Index ranges for slicing the arrays. If None, the full range is used.
    use_zip : bool, optional
        If True, load from a .zarr.zip file instead of a directory. By default, False.
    Returns
    -------
    dict
        loaded array from zarr data
    """

    root = open_zarr_group(path, use_zip=use_zip)

    if i0 is None:
        i0 = 0
    if j0 is None:
        j0 = 0
    if k0 is None:
        k0 = 0

    if isinstance(names, str):
        if names == "all":
            names = list(root.array_keys())
        else:
            names = [names]

    for key in list(root.array_keys()):
        if root[key].ndim == 3:
            i1_tmp, j1_tmp, k1_tmp = root[key].shape
            break
        if root[key].ndim == 2:
            i1_tmp, j1_tmp = root[key].shape
            k1_tmp = 1
            break

    if i1 is None:
        i1 = i1_tmp
    if j1 is None:
        j1 = j1_tmp
    if k1 is None:
        k1 = k1_tmp

    if i0 >= i1:
        raise ValueError(f"Invalid i range: i0={i0} must be less than i1={i1}")
    if j0 >= j1:
        raise ValueError(f"Invalid j range: j0={j0} must be less than j1={j1}")
    if k0 >= k1:
        raise ValueError(f"Invalid k range: k0={k0} must be less than k1={k1}")

    data = {}
    for name in names:
        if name in ["x", "y", "z"]:
            data[name] = root[name]
        else:
            if root[name].ndim == 3:
                data[name] = root[name][i0:i1, j0:j1, k0:k1]
            elif root[name].ndim == 2:
                data[name] = root[name][i0:i1, j0:j1]

    if with_attrs:
        attrs = root.attrs.asdict()
        return data, attrs["params"]

    return data


def list_vars(path: str, use_zip: bool = False):
    """
    List variables in zarr data

    Parameters
    ----------
    path : str
        path to zarr data
    use_zip : bool, optional
        if True, read from a .zarr.zip file instead of a directory. By default, False.

    Returns
    -------
    list
        list of variable names in the zarr data
    """
    root = open_zarr_group(path, use_zip=use_zip)
    return list(root.array_keys())
