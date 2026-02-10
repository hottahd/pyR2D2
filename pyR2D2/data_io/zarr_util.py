import numpy as np
import zarr


def save(
    path: str,
    vars_dict: dict,
    params: dict = None,
    max_chunk_size: int = 512,
    clevel: int = 5,
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
    clevel : int, optional
        compression level for zarr, by default 5

    Raises
    ------
    ValueError
        If an array with unsupported dimensions is encountered.
    """
    root = zarr.open_group(path, mode="w")

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
            chunks = tuple(min(max_chunk_size, array.shape[i]) for i in range(3))
            root.create_array(
                name,
                data=array,
                chunks=chunks,
                compressors=codec_f32,
                overwrite=True,
            )
        elif array.ndim == 2:
            chunks = tuple(min(max_chunk_size, array.shape[i]) for i in range(2))
            root.create_array(
                name,
                data=array,
                chunks=chunks,
                compressors=codec_f32,
                overwrite=True,
            )
        elif array.ndim == 1:
            chunks = (min(max_chunk_size, array.shape[0]),)
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


def load(path: str, names="all", with_attrs: bool = False):
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

    Returns
    -------
    np.ndarray or dict
        loaded array from zarr data
    """
    root = zarr.open_group(path, mode="r")

    if isinstance(names, str):
        if names == "all":
            names = list(root.array_keys())
        else:
            names = [names]

    data = {}
    for name in names:
        data[name] = root[name][...]

    if len(data) == 1:
        data = data[names[0]]

    if with_attrs:
        attrs = root.attrs.asdict()
        return data, attrs["params"]

    return data


def list_vars(path: str):
    """
    List variables in zarr data

    Parameters
    ----------
    path : str
        path to zarr data

    Returns
    -------
    list
        list of variable names in the zarr data
    """
    root = zarr.open_group(path, mode="r")
    return list(root.array_keys())
