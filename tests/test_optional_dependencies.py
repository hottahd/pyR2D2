"""Optional dependencies are lazy and fail with actionable messages."""

import builtins

import pytest

from pyR2D2.data_io import zarr_util


def test_zarr_is_not_imported_eagerly():
    """Importing pyR2D2 must not initialize the optional Zarr package."""
    assert "zarr" not in zarr_util.__dict__


def test_missing_zarr_error_mentions_install_extra(monkeypatch):
    """A requested Zarr operation explains how to install the extra."""
    original_import = builtins.__import__

    def import_without_zarr(name, *args, **kwargs):
        if name == "zarr":
            raise ModuleNotFoundError("blocked by test", name="zarr")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", import_without_zarr)

    with pytest.raises(ModuleNotFoundError, match=r"pyR2D2\[zarr\]"):
        zarr_util._require_zarr()


def test_missing_zarr_transitive_dependency_is_not_misdiagnosed(monkeypatch):
    """A broken Zarr installation keeps its original missing-dependency error."""
    original_import = builtins.__import__

    def import_with_broken_zarr(name, *args, **kwargs):
        if name == "zarr":
            raise ModuleNotFoundError("missing numcodecs", name="numcodecs")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", import_with_broken_zarr)

    with pytest.raises(ModuleNotFoundError) as error:
        zarr_util._require_zarr()
    assert error.value.name == "numcodecs"
