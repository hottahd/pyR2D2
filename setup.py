import os
import shutil
import sys
from importlib.machinery import EXTENSION_SUFFIXES
from pathlib import Path

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import Command, setup


def _native_optimization_enabled():
    """Return whether this explicitly local build requested -march=native."""
    value = os.environ.get("PYR2D2_NATIVE", "").strip().lower()
    return value in {"1", "true", "yes", "on"}


def _compiler_flags():
    """Return compiler and linker flags for the current platform."""
    if sys.platform == "win32":
        return ["/O2", "/DNDEBUG", "/openmp"], []

    compile_args = ["-O3", "-DNDEBUG"]
    link_args = []
    if sys.platform == "linux":
        compile_args.append("-fopenmp")
        link_args.append("-fopenmp")
        if _native_optimization_enabled():
            compile_args.append("-march=native")
    return compile_args, link_args


def _is_absolute_rpath_flag(flag):
    """Identify an absolute GNU-style RPATH inherited from the build Python."""
    prefix = "-Wl,-rpath,"
    return flag.startswith(prefix) and os.path.isabs(flag[len(prefix) :])


class BuildExt(build_ext):
    """Build the extension without leaking the build environment's RPATH."""

    def _sanitize_compiler_linker(self):
        # setuptools 83 uses linker_so_cxx for C++ shared libraries, while
        # older versions use linker_so. Handle both without changing other
        # compiler or library-search options.
        for attribute in ("linker_so", "linker_so_cxx"):
            linker = getattr(self.compiler, attribute, None)
            if linker:
                setattr(
                    self.compiler,
                    attribute,
                    [
                        flag
                        for flag in linker
                        if not _is_absolute_rpath_flag(flag)
                    ],
                )

    def build_extensions(self):
        self._sanitize_compiler_linker()
        super().build_extensions()

    def build_extension(self, extension):
        # pybind11 may update compiler options inside build_extensions(), so
        # sanitize again at the last point before setuptools links the module.
        self._sanitize_compiler_linker()
        super().build_extension(extension)


class Clean(Command):
    """Remove build directories and extension binaries produced in place."""

    description = "remove pyR2D2 build artifacts"
    user_options = [("all", "a", "remove all build output")]
    boolean_options = ["all"]

    def initialize_options(self):
        self.all = False

    def finalize_options(self):
        pass

    def run(self):
        project_root = Path(__file__).parent
        build_dir = project_root / "build"
        if build_dir.exists():
            shutil.rmtree(build_dir)
            self.announce(f"removed '{build_dir}'", level=2)

        extension_dir = project_root / "pyR2D2" / "cpp_util"
        suffixes = set(EXTENSION_SUFFIXES) | {".so", ".pyd", ".dylib"}
        for suffix in suffixes:
            for artifact in extension_dir.glob(f"cpp_util*{suffix}"):
                artifact.unlink()
                self.announce(f"removed '{artifact}'", level=2)


extra_compile_args, extra_link_args = _compiler_flags()

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
    ext_modules=ext_modules,
    cmdclass={"build_ext": BuildExt, "clean": Clean},
)
