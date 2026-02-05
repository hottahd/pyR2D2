import os
import subprocess
import sys

from pybind11.setup_helpers import Pybind11Extension
from setuptools import find_packages, setup
from setuptools.command.build_py import build_py
from setuptools_scm import get_version

if sys.platform == "linux":
    extra_compile_args = ["-O3", "-fopenmp"]
    extra_link_args = ["-fopenmp"]
else:
    extra_compile_args = ["-O3"]
    extra_link_args = []


class CustomBuildPy(build_py):
    def run(self):
        # Fortranコードが存在するディレクトリ
        fortran_dir = os.path.join(
            os.path.dirname(__file__), "pyR2D2/fortran_util/fortran_src"
        )

        # makeコマンドを実行
        print("Running 'make' to build Fortran library...")
        subprocess.check_call(["make"], cwd=fortran_dir)

        # ビルドが成功したら、親クラスのrun()を呼び出す
        super().run()


ext_modules = [
    Pybind11Extension(
        "pyR2D2.cpp_util.cpp_util",
        ["pyR2D2/cpp_util/bindings.cpp"],
        cxx_std=14,
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
    setup_requires=["setuptools_scm"],
    packages=find_packages(include=["pyR2D2", "pyR2D2.*"]),  # 対象パッケージを指定
    include_package_data=True,
    cmdclass={"build_py": CustomBuildPy},  # build_pyをカスタムコマンドに置き換え
    zip_safe=False,
    ext_modules=ext_modules,
)
