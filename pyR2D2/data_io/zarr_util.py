import numpy as np
import zarr


def save(
    path: str,
    vars_dict: dict,
    chunks3d: tuple = (64, 64, 64),
    chunks1d: tuple = (32,),
    clevel: int = 5,
):
    """


    Parameters
    ----------
    path : str
        path to save zarr data
    vars_dict : dict
        dictionary of variables to save
    chunks3d : tuple
        chunk size for 3D arrays
    chunks1d : tuple
        chunk size for 1D arrays
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

    for name, array in vars_dict.items():
        array = np.asarray(array, dtype=np.float32, order="C")
        if array.ndim == 3:
            root.create_array(
                name,
                data=array,
                chunks=chunks3d,
                compressors=codec_f32,
                overwrite=True,
            )
        elif array.ndim == 1:
            root.create_array(
                name,
                data=array,
                chunks=chunks1d,
                compressors=codec_f32,
                overwrite=True,
            )
        else:
            raise ValueError(
                f"Unsupported array dimension: {array.ndim} for variable {name}"
            )


def load(path: str, name: str):
    """
    Load zarr data

    Parameters
    ----------
    path : str
        path to zarr data
    name : str
        name of the variable to load

    Returns
    -------
    np.ndarray
        loaded array from zarr data
    """
    root = zarr.open_group(path, mode="r")
    return root[name][...]


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
