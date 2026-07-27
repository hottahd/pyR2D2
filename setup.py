import sys

from pybind11.setup_helpers import Pybind11Extension
from setuptools import find_packages, setup

if sys.platform == "linux":
    extra_compile_args = ["-O3", "-DNDEBUG", "-march=native", "-fopenmp"]
    extra_link_args = ["-fopenmp"]
else:
    extra_compile_args = ["-O3", "-DNDEBUG"]
    extra_link_args = []


ext_modules = [
    Pybind11Extension(
        "pyR2D2.cpp_util.cpp_util",
        ["pyR2D2/cpp_util/bindings.cpp"],
        cxx_std=20,
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
    ),
]

setup(
    name="pyR2D2",
    use_scm_version={
        "root": ".",
        "relative_to": __file__,
    },
    packages=find_packages(include=["pyR2D2", "pyR2D2.*"]),  # 対象パッケージを指定
    include_package_data=True,
    zip_safe=False,
    ext_modules=ext_modules,
)
