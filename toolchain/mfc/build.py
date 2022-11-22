import os, typing, dataclasses

from .state   import ARG
from .printer import cons
from .        import common


@dataclasses.dataclass
class MFCTarget:
    name:         str
    flags:        typing.List[str]
    requires:     typing.List[str]
    isDependency: bool
    isDefault:    bool


TARGETS: typing.List[MFCTarget] = [
    MFCTarget(name='fftw', flags=['-DMFC_FFTW=ON'],
              isDependency=True, isDefault=False, requires=[]),
    MFCTarget(name='hdf5', flags=['-DMFC_HDF5=ON'],
              isDependency=True, isDefault=False, requires=[]),
    MFCTarget(name='silo', flags=['-DMFC_SILO=ON'],
              isDependency=True, isDefault=False, requires=["hdf5"]),
    MFCTarget(name='pre_process', flags=['-DMFC_PRE_PROCESS=ON'],
              isDependency=False, isDefault=True, requires=[]),
    MFCTarget(name='simulation', flags=['-DMFC_SIMULATION=ON'],
              isDependency=False, isDefault=True, requires=["fftw"]),
    MFCTarget(name='post_process', flags=['-DMFC_POST_PROCESS=ON'],
              isDependency=False, isDefault=True, requires=['fftw', 'silo']),
    MFCTarget(name="documentation", flags=['-DMFC_DOCUMENTATION=ON'],
              isDependency=False, isDefault=False, requires=[])
]


def get_mfc_target_names() -> typing.List[str]:
    return [ target.name for target in TARGETS if target.isDefault ]


def get_dependencies_names() -> typing.List[str]:
    return [ target.name for target in TARGETS if target.isDependency ]


def get_target_names() -> typing.List[str]:
    return [ target.name for target in TARGETS ]


def get_target(name: str) -> MFCTarget:
    for target in TARGETS:
        if target.name == name:
            return target

    raise common.MFCException(f"Target '{name}' does not exist.")


# Get path to directory that will store the build files
def get_build_dirpath(target: MFCTarget) -> str:
    return os.sep.join([os.getcwd(), "build", target.name])


# Get the directory that contains the target's CMakeLists.txt
def get_cmake_dirpath(target: MFCTarget) -> str:
    subdir = ["", os.sep.join(["toolchain", "dependencies"])][int(target.isDependency)]
    
    return os.sep.join([os.getcwd(), subdir])


def get_install_dirpath() -> str:
    return os.sep.join([os.getcwd(), "build", "install"])


def is_target_configured(target: MFCTarget) -> bool:
    build_dirpath = get_build_dirpath(target)
    return os.path.isdir(build_dirpath)


def clean_target(name: str):
    cons.print(f"[bold]Cleaning [magenta]{name}[/magenta]:[/bold]")
    cons.indent()

    target = get_target(name)

    build_dirpath = get_build_dirpath(target)

    if not os.path.isdir(build_dirpath):
        cons.print("Target not configured. Nothing to clean.")
        cons.unindent()
        return

    clean = ["cmake", "--build",  build_dirpath, "--target", "clean",
                      "--config", "Debug" if ARG("debug") else "Release" ]

    if ARG("verbose"):
        clean.append("--verbose")

    common.system(clean, exception_text=f"Failed to clean the [bold magenta]{name}[/bold magenta] target.")

    cons.unindent()


def clean_targets(targets: typing.List[str]):
    for target in targets:
        clean_target(target)


def clean():
    clean_targets(ARG("targets"))


def build_target(name: str, history: typing.List[str] = None):
    cons.print(f"[bold]Building [magenta]{name}[/magenta]:[/bold]")
    cons.indent()

    if history is None:
        history = []

    if ARG("no_build"):
        cons.print("--no-build specified, skipping...")
        cons.unindent()
        return

    if name in history:
        cons.print("Already built, skipping...")
        cons.unindent()
        return

    history.append(name)

    target = get_target(name)

    if target.isDependency and ARG(f"no_{target.name}"):
        cons.print(f"--no-{target.name} given, skipping...")
        cons.unindent()
        return

    build_dirpath   = get_build_dirpath(target)
    cmake_dirpath   = get_cmake_dirpath(target)
    install_dirpath = get_install_dirpath()

    flags: list = target.flags.copy() + [
        # Disable CMake warnings intended for developers (us).
        # See: https://cmake.org/cmake/help/latest/manual/cmake.1.html.
        f"-Wno-dev",
        # Save a compile_commands.json file with the compile commands used to
        # build the configured targets. This is mostly useful for debugging.
        # See: https://cmake.org/cmake/help/latest/variable/CMAKE_EXPORT_COMPILE_COMMANDS.html.
        f"-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
        # Set build type (e.g Debug, Release, etc.).
        # See: https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html
        f"-DCMAKE_BUILD_TYPE={'Debug' if ARG('debug') else 'Release'}",
        # Used by FIND_PACKAGE (/FindXXX) to search for packages, with the
        # second heighest level of priority, still letting users manually
        # specify <PackageName>_ROOT, which has precedence over CMAKE_PREFIX_PATH.
        # See: https://cmake.org/cmake/help/latest/command/find_package.html.
        f"-DCMAKE_PREFIX_PATH={install_dirpath}",
        # First directory that FIND_LIBRARY searches.
        # See: https://cmake.org/cmake/help/latest/command/find_library.html.
        f"-DCMAKE_FIND_ROOT_PATH={install_dirpath}",
        # Location prefix to install bin/, lib/, include/, etc.
        # See: https://cmake.org/cmake/help/latest/command/install.html.
        f"-DCMAKE_INSTALL_PREFIX={install_dirpath}",
    ]

    if not target.isDependency:
        flags.append(f"-DMFC_MPI={    'ON' if ARG('mpi') else 'OFF'}")
        flags.append(f"-DMFC_OpenACC={'ON' if ARG('gpu') else 'OFF'}")

    configure = ["cmake"] + flags + ["-S", cmake_dirpath, "-B", build_dirpath]

    if common.does_command_exist("ninja"):
        configure.append('-GNinja')

    build     = ["cmake", "--build",  build_dirpath,
                          "--target", name,
                          "-j",       ARG("jobs"), 
                          "--config", 'Debug' if ARG('debug') else 'Release']
    if ARG('verbose'):
        build.append("--verbose")

    install   = ["cmake", "--install", build_dirpath]

    # Only configure the first time
    if not is_target_configured(target):
        for dependency_name in target.requires:
            build_target(dependency_name, history)

        common.create_directory(build_dirpath)

        if common.system(configure, no_exception=True) != 0:
            common.delete_directory(build_dirpath)

            raise common.MFCException(f"Failed to configure the [bold magenta]{name}[/bold magenta] target.")

    common.system(build,   exception_text=f"Failed to build the [bold magenta]{name}[/bold magenta] target.")
    common.system(install, exception_text=f"Failed to install the [bold magenta]{name}[/bold magenta] target.")
    cons.print(no_indent=True)

    cons.unindent()


def build_targets(targets):
    for target in targets:
        build_target(target)


def build():
    build_targets(ARG("targets"))
