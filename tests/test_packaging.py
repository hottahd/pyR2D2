"""Packaging metadata and runtime-data policy regression tests."""

import importlib.resources
import tomllib
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def load_pyproject():
    """Return the parsed PEP 621/build configuration."""
    with (PROJECT_ROOT / "pyproject.toml").open("rb") as stream:
        return tomllib.load(stream)


def test_supported_python_matches_dependency_floor():
    """NumPy 2.1 and Zarr 3 require a newer Python than the old 3.9 claim."""
    project = load_pyproject()["project"]
    assert project["requires-python"] == ">=3.11"


def test_optional_dependencies_are_not_core_requirements():
    """Zarr and plotting packages are installed only for requested features."""
    project = load_pyproject()["project"]
    dependencies = project["dependencies"]
    extras = project["optional-dependencies"]

    assert any(requirement.startswith("numpy") for requirement in dependencies)
    assert any(requirement.startswith("scipy") for requirement in dependencies)
    assert not any(requirement.startswith("zarr") for requirement in dependencies)
    assert not any(
        requirement.startswith("matplotlib") for requirement in dependencies
    )
    assert "zarr>=3.0.8,<4" in extras["zarr"]
    assert any(
        requirement.startswith("matplotlib") for requirement in extras["plot"]
    )


def test_distribution_is_explicitly_private():
    """Prevent accidental publication while no redistributable license exists."""
    project = load_pyproject()["project"]
    assert "Private :: Do Not Upload" in project["classifiers"]


def test_runtime_json_files_are_declared_and_available():
    """Explicit package data includes both JSON files used by data_io."""
    configuration = load_pyproject()
    declared = configuration["tool"]["setuptools"]["package-data"][
        "pyR2D2.data_io"
    ]
    assert declared == ["*.json"]

    data_io = importlib.resources.files("pyR2D2.data_io")
    assert data_io.joinpath("OnTheFly.json").is_file()
    assert data_io.joinpath("parameters.json").is_file()


def test_setuptools_is_the_only_extension_build_path():
    """Do not reintroduce a second CMake/Makefile configuration."""
    assert not (PROJECT_ROOT / "pyR2D2" / "Makefile").exists()
    assert not (PROJECT_ROOT / "pyR2D2" / "cpp_util" / "CMakeLists.txt").exists()
