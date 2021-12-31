#!/usr/bin/env python3

import os
import re
import sys
import shutil
import tarfile
import argparse
import dataclasses
import urllib.request


@dataclasses.dataclass
class MFCGlobalState:
    conf: dict = dataclasses.field(default_factory=dict)
    lock: dict = dataclasses.field(default_factory=dict)
    args: dict = dataclasses.field(default_factory=dict)

    root_path: str = ""


class MFCException(Exception):
    pass


def execute_shell_command_safe(command: str):
    if os.system(command) != 0:
        raise MFCException(f'Failed to execute command "{command}".')


def import_module_safe(import_name: str, pip_name: str):
    while True:
        try:
            globals()[import_name] = __import__(import_name)

            break
        except ImportError as exc:
            prompt = f"The Python package {exc.name} needs to be installed. Would you like to install it now? (yY/nN): "

            if input(prompt).upper().strip() != "Y":
                exit(1)

            if 0 != execute_shell_command_safe("python3 -m pip install --user {}".format(pip_name)):
                exit(1)


def clear_line():
    sys.stdout.write("\033[K")


def file_load_yaml(filepath: str):
    try:
        with open(filepath, "r") as f:
            return yaml.safe_load(f)
    except (IOError, yaml.YAMLError) as exc:
        raise MFCException(f"Failed to load YAML from \"{filepath}\": {exc}")


def file_dump_yaml(filepath: str, data):
    try:
        with open(filepath, "w") as f:
            yaml.dump(data, f)
    except (IOError, yaml.YAMLError) as exc:
        raise MFCException(f"Failed to dump YAML to \"{filepath}\": {exc}.")


# TODO: Better solution
def uncompress_archive_to(archive_filepath, destination):
    archive = tarfile.open(archive_filepath, "r")
    archive.extractall("/".join(destination.rstrip("/").split("/")[:-1]))
    archive.close()

    filename_wo_ext = archive_filepath.replace(".tar.gz", "")

    os.rename(filename_wo_ext, destination)


def delete_directory_recursive_safe(dirpath):
    if os.path.isdir(dirpath):
        shutil.rmtree(dirpath)


def center_ansi_escaped_text(message: str):
    nCols = shutil.get_terminal_size((80, 20)).columns

    to_escape = [re.escape(colorama.Style.RESET_ALL)]
    for key in dir(colorama.Fore):
        if not callable(getattr(colorama.Fore, key)) and not key.startswith("__"):
            to_escape.append(re.escape(getattr(colorama.Fore, key)))

    longest_string_len = max([len(re.compile("|".join(to_escape), flags=re.DOTALL).sub("", line)) for line in message.splitlines()])

    padding = " "*((nCols - longest_string_len) // 2)

    return "\n".join([f'{padding}{line}{padding}' for line in message.splitlines()])


def mfc_read_configuration_file():
    return file_load_yaml("mfc.cfg.yaml")


def mfc_read_lock_file():
    if not os.path.exists("mfc.lock.yaml"):
        with open("mfc.lock.yaml", 'w') as f:
            f.write("targets: []")

    return file_load_yaml("mfc.lock.yaml")


def mfc_write_lock_file(mfc_state: MFCGlobalState):
    file_dump_yaml("mfc.lock.yaml", mfc_state.lock)


def mfc_parse_arguments(mfc_conf: dict):
    parser = argparse.ArgumentParser(description="Wecome to the MFC master script.", )

    general = parser.add_argument_group(title="General")
    compiling = parser.add_argument_group(title="Compiling")

    compiler_target_names = [e["name"] for e in mfc_conf["targets"]]
    compiling.add_argument("-t", "--targets", nargs="+", type=str,
                           choices=compiler_target_names, default="",
                           help="The space-separated targets you wish to have built.")

    compiler_configuration_names = [e["name"] for e in mfc_conf["compilers"]["configurations"]]
    compiling.add_argument("-cc", "--compiler-configuration", type=str,
                           choices=compiler_configuration_names, default=compiler_configuration_names[1])

    compiling.add_argument("-j", "--jobs", metavar="N", type=int,
                           help="Allows for N concurrent jobs.", default=1)

    compiling.add_argument("-c", "--clean", action="store_true",
                           help="Cleans all binaries and build artificats.")

    general.add_argument("-f", "--force", action="store_true",
                         help="Force the operation.")

    return vars(parser.parse_args())


def mfc_print_header():
    print(center_ansi_escaped_text("""\

                                      {0}      ___            ___          ___      {2}  
                                      {0}     /__/\          /  /\        /  /\     {2}  
                                   |  {0}    |  |::\        /  /:/_      /  /:/     {2}  |
+----------------------------------+  {0}    |  |:|:\      /  /:/ /\    /  /:/      {2}  +----------------------------------+
|__/\________/\__/\__/\________/\__|  {0}  __|__|:|\:\    /  /:/ /:/   /  /:/  ___  {2}  |    Multi-component Flow Code     |
|_/__\__/\__/__\/__\/__\__/\__/__\_|  {0} /__/::::| \:\  /__/:/ /:/   /__/:/  /  /\ {2}  |                                  |
|/____\/__\/____________\/__\/____\|  {0} \  \:\~~\__\/  \  \:\/:/    \  \:\ /  /:/ {2}  | https://github.com/MFlowCode/MFC |
+----------------------------------+  {0}  \  \:\         \  \::/      \  \:\  /:/  {2}  +----------------------------------+
                                   |  {0}   \  \:\         \  \:\       \  \:\/:/   {2}  |
                                      {0}    \  \:\         \  \:\       \  \::/    {2}  
                                      {0}     \__\/          \__\/        \__\/     {2}  

""".format(colorama.Fore.BLUE, colorama.Fore.MAGENTA, colorama.Style.RESET_ALL)))


def mfc_get_dependency(obj: dict, name: str):
    matches = list(filter(lambda x: x["name"] == name, obj["targets"]))

    if len(matches) != 1:
        raise MFCException(f"Failed to find dependency \"{name}\" in \"{obj}\".")

    return matches[0]


def mfc_string_replace(mfc_state: MFCGlobalState, obj: dict, dependency_name: str, string: str, recursive=True):
    dep = mfc_get_dependency(obj, dependency_name)

    if dep["type"] in ["clone", "download"]:
        install_path = mfc_state.root_path + "/dependencies/build"
        source_path  = mfc_state.root_path + "/dependencies/src/" + dependency_name
    elif dep["type"] == "source":
        install_path = "ERR_INSTALL_PATH_IS_UNDEFINED"
        source_path  = dep["source"]["source"]
    else:
        raise MFCException(f'Unknown type "{dep["type"]}".')

    regulars, mpis = mfc_state.conf["compilers"]["regular"], mfc_state.conf["compilers"]["mpi"]

    compiler_cfg = None
    for c_cfg in mfc_state.conf["compilers"]["configurations"]:
        if c_cfg["name"] == mfc_state.args["compiler_configuration"]:
            compiler_cfg = c_cfg
            break

    if compiler_cfg is None:
        raise MFCException(
            "Failed to locate the compiler configuration \"{}\".".format(mfc_state.args["compiler_configuration"]))

    flags = compiler_cfg["flags"]

    replace_list = [
        ("${MFC_ROOT_PATH}",     mfc_state.root_path),
        ("${CONFIGURE_OPTIONS}", f'--prefix="{install_path}"'),
        ("${SOURCE_PATH}",       source_path),
        ("${INSTALL_PATH}",      install_path),
        ("${INSTALL_PATH}",      install_path),
        ("${MAKE_OPTIONS}",      f'-j {mfc_state.args["jobs"]}'),
        ("${COMPILER_FLAGS}",    f'CFLAGS="{flags["c"]}" CPPFLAGS="{flags["c++"]}" FFLAGS="{flags["fortran"]}"'),
        ("${REGULAR_COMPILERS}", f'CC="{regulars["c"]}"  CXX="{regulars["c++"]}"   FC="{regulars["fortran"]}"'),
        ("${MPI_COMPILERS}",     f'CC="{mpis["c"]}"      CXX="{mpis["c++"]}"       FC="{mpis["fortran"]}"')
    ]

    for e in replace_list:
        string = string.replace(*e)

    if recursive:
        for dep2_info in obj["targets"]:
            string = string.replace("${" + dep2_info["name"] + ":", "${")
            string = mfc_string_replace(mfc_state, obj, dep2_info["name"], string, recursive=False)

    return string


def mfc_was_dependency_updated(cfg_info, cur_info):
    if cfg_info["type"] != cur_info["type"]:
        return True

    if cfg_info["type"] == "source":
        return True

    t = cfg_info["type"]
    if cfg_info[t] != cur_info[t]:
        return True

    return False


def mfc_build_dependency(mfc_state: MFCGlobalState, name: str):
    cfg = mfc_get_dependency(mfc_state.conf, name)

    # Create Basic Directories
    for e in ["src", "build", "log"]:
        os.makedirs(mfc_state.root_path + "/dependencies/" + e, exist_ok=True)

    # Build the dependencies it depends on
    bHadToBuildADependency = False
    for depends_name in cfg.get("depends", []):
        bHadToBuildADependency |= mfc_build_dependency(mfc_state, depends_name)

    cur = None
    cur_matches = list(filter(lambda x: x["name"] == name, mfc_state.lock["targets"]))
    bWasNeverBuilt = len(cur_matches) != 1

    # TODO: desciption
    bDifferentCFGandCUR = False
    if len(cur_matches) == 1:
        cur = cur_matches[0]

        bDifferentCFGandCUR = mfc_was_dependency_updated(cfg, cur)

        # Remove old directories from previous build/download type
        if cfg["type"] != cur["type"]:
            if cur["type"] in ["clone", "download"]:
                delete_directory_recursive_safe(mfc_state.root_path + "/dependencies/src/" + name)

    if not (bWasNeverBuilt or bHadToBuildADependency or bDifferentCFGandCUR):
        return False

    # Create Log File
    log_filepath = mfc_state.root_path + "/dependencies/log/" + name + ".log"
    if os.path.isfile(log_filepath):
        os.remove(log_filepath)

    # 
    if cfg["type"] == "clone" or cfg["type"] == "download":
        source_path = mfc_state.root_path + "/dependencies/src/" + name

        if cfg["type"] == "clone":
            if cur != None and os.path.isdir(source_path):
                if cfg["clone"]["git"] != cur["clone"]["git"]:
                    delete_directory_recursive_safe(source_path)

            if not os.path.isdir(source_path):
                execute_shell_command_safe(
                    "git clone --recursive \"{}\" \"{}\" >> \"{}\" 2>&1".format(cfg["clone"]["git"], source_path,
                                                                                log_filepath))

            execute_shell_command_safe(
                "cd \"{}\" && git checkout \"{}\" >> \"{}\" 2>&1".format(source_path, cfg["clone"]["hash"],
                                                                         log_filepath))
        elif cfg["type"] == "download":
            delete_directory_recursive_safe(mfc_state.root_path + "/dependencies/src/" + name)

            download_link = cfg[cfg["type"]]["link"].replace("${version}", cfg[cfg["type"]]["version"])
            filename = download_link.split("/")[-1]

            urllib.request.urlretrieve(download_link, mfc_state.root_path + "/dependencies/src/" + filename)

            uncompress_archive_to(mfc_state.root_path + "/dependencies/src/" + filename,
                                  mfc_state.root_path + "/dependencies/src/" + name)

            os.remove(mfc_state.root_path + "/dependencies/src/" + filename)
    elif cfg["type"] == "source":
        pass
    else:
        raise MFCException("Dependency type \"{}\" is unsupported.", cfg["type"])

    # Build
    if cfg["type"] in ["clone", "download", "source"]:
        for command in cfg["build"]:
            if cfg["type"] in ["clone", "download", "source"]:
                command = mfc_string_replace(mfc_state, mfc_state.conf, name, f"""\
cd "${{SOURCE_PATH}}" && \
PYTHON="python3" PYTHON_CPPFLAGS="$PYTHON_CPPFLAGS $(python3-config --includes) $(python3-config --libs)" \
bash -c '{command}' >> "{log_filepath}" 2>&1""")

                print(command)
                execute_shell_command_safe(command)
            else:
                raise MFCException(f'Dependency type "{cfg["type"]}" is unsupported.')
    else:
        raise MFCException(f'Unknown type "{cfg["type"]}".')

    if cfg["type"] == "source":
        return True

    # Update CUR
    if len(cur_matches) == 0:
        mfc_state.lock["targets"].append(cfg)
    else:
        for index, dep in enumerate(mfc_state.lock["targets"]):
            if dep["name"] == name:
                mfc_state.lock["targets"][index] = cfg

    mfc_write_lock_file(mfc_state)

    return True


def mfc_environment_checks(mfc_state: MFCGlobalState):
    print("|--> Checking for the presence of required command-line utilities...", end='\r')

    required = ["python3", "python3-config", "make", "git"]
    required += mfc_state.conf["compilers"]["regular"].values()
    required += mfc_state.conf["compilers"]["mpi"].values()

    for index, utility in enumerate(required):
        clear_line()
        print(f"|--> {index+1}/{len(required)} Checking for {utility}...", end='\r')
        if shutil.which(utility) is None:
            raise MFCException(
                "Failed to find the command line utility \"{}\". Please install it or make it visible.".format(utility))

    # TODO: MacOS CHECKS
    if sys.platform == "darwin":  # MacOS
        pass
        #mfc_state.conf["compilers"]["mpi"]["fortran"]


def mfc_main():
    for module in [("yaml", "pyyaml"), ("colorama", "colorama")]:
        import_module_safe(*module)

    colorama.init()

    try:
        conf = mfc_read_configuration_file()
        lock = mfc_read_lock_file()
        args = mfc_parse_arguments(conf)
        root_path = os.path.dirname(os.path.realpath(__file__))
        mfc_state = MFCGlobalState(conf=conf, lock=lock, args=args, root_path=root_path)

        mfc_print_header()
        mfc_environment_checks(mfc_state)

        if mfc_state.args["clean"]:
            delete_directory_recursive_safe("./dependencies/")

            mfc_state.lock["targets"] = []

            mfc_write_lock_file(mfc_state)

        for target_name in mfc_state.args["targets"]:
            mfc_build_dependency(mfc_state, target_name)

        mfc_write_lock_file(mfc_state)
    except MFCException as exc:
        print(f"{colorama.Fore.RED}|--> {str(exc)}{colorama.Style.RESET_ALL}")
        exit(1)


if __name__ == "__main__":
    mfc_main()

