from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
from setuptools_scm import get_version
import os
import subprocess


class CustomBuildPy(build_py):
    def run(self):
        # Fortranコードが存在するディレクトリ
        fortran_dir = os.path.join(os.path.dirname(__file__), "pyR2D2/fortran_util/fortran_src")

        # makeコマンドを実行
        print("Running 'make' to build Fortran library...")
        subprocess.check_call(["make"], cwd=fortran_dir)

        # ビルドが成功したら、親クラスのrun()を呼び出す
        super().run()


setup(
    name="pyR2D2",
    use_scm_version={
        "root": ".",
        "relative_to": __file__,
    },
    setup_requires=["setuptools_scm"],
    packages=find_packages(),
    cmdclass={"build_py": CustomBuildPy},  # build_pyをカスタムコマンドに置き換え
    zip_safe=False,
)