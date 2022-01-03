#!/usr/bin/env python3

import os
import re
import sys
import shutil
import tarfile
import argparse
import dataclasses
import urllib.request


MFC_USE_SUBDIR    = ".mfc"
MFC_CONF_FILEPATH = f"mfc.conf.yaml"
MFC_LOCK_FILEPATH = f"{MFC_USE_SUBDIR}/mfc.lock.yaml"


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
            prompt = f"The Python package {pip_name} needs to be installed. Would you like to install it now? (yY/nN): "

            if input(prompt).upper().strip() != "Y":
                raise MFCException(f'You requested not to download the Python package {pip_name}. Aborting...')

            if 0 != execute_shell_command_safe(f'python3 -m pip install --user {pip_name}'):
                raise MFCException(f'Failed to install Python package {pip_name} automatically. Please install it manually using Pip. Aborting...')


def clear_line():
    sys.stdout.write("\033[K")


def file_load_yaml(filepath: str):
    try:
        with open(filepath, "r") as f:
            return yaml.safe_load(f)
    except (IOError, yaml.YAMLError) as exc:
        raise MFCException(f'Failed to load YAML from "{filepath}": {exc}')


def file_dump_yaml(filepath: str, data):
    try:
        with open(filepath, "w") as f:
            yaml.dump(data, f)
    except (IOError, yaml.YAMLError) as exc:
        raise MFCException(f'Failed to dump YAML to "{filepath}": {exc}.')


# TODO: Better solution
def uncompress_archive_to(archive_filepath, destination):
    archive = tarfile.open(archive_filepath, "r")
    archive.extractall("/".join(destination.rstrip("/").split("/")[:-1]))
    archive.close()

    filename_wo_ext = archive_filepath.replace(".tar.gz", "")

    os.rename(filename_wo_ext, destination)


def delete_file_safe(filepath: str):
    if os.path.exists(filepath):
        os.remove(filepath)


def create_file_safe(filepath: str):
    if not os.path.exists(filepath):
        open(filepath, "w").close()


def delete_directory_recursive_safe(directory_path: str):
    if os.path.isdir(directory_path):
        shutil.rmtree(directory_path)


def create_directory_safe(directory_path: str):
    os.makedirs(directory_path, exist_ok=True)


def center_ansi_escaped_text(message: str):
    nCols = shutil.get_terminal_size((80, 20)).columns

    to_escape = [re.escape(colorama.Style.RESET_ALL)]
    for key in dir(colorama.Fore):
        if not callable(getattr(colorama.Fore, key)) and not key.startswith("__"):
            to_escape.append(re.escape(getattr(colorama.Fore, key)))

    longest_string_len = max([len(re.compile("|".join(to_escape), flags=re.DOTALL).sub("", line)) for line in message.splitlines()])

    padding = " "*((nCols - longest_string_len) // 2)

    return "\n".join([f'{padding}{line}{padding}' for line in message.splitlines()])


def clear_print(message, end='\n'):
    clear_line()
    print(message, end=end)


def mfc_read_lock_file():
    if not os.path.exists(f"{MFC_LOCK_FILEPATH}"):
        with open(f"{MFC_LOCK_FILEPATH}", 'w') as f:
            f.write("targets: []")

    return file_load_yaml(f"{MFC_LOCK_FILEPATH}")


def mfc_dump_lock_file(mfc_state: MFCGlobalState):
    file_dump_yaml(f"{MFC_LOCK_FILEPATH}", mfc_state.lock)


def mfc_read_conf_file():
    return file_load_yaml(MFC_CONF_FILEPATH)


def mfc_parse_arguments(mfc_conf: dict):
    parser = argparse.ArgumentParser(description="Wecome to the MFC master script.", )

    general   = parser.add_argument_group(title="General")
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
                           help="Cleans all binaries and build artifacts.")

    compiling.add_argument("-s", "--scratch", action="store_true",
                         help="Build all targets from scratch.")

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


def mfc_peek_target(obj: dict, name: str, restrict_cc: str = None):
    def peek_filter(e: dict):
        if e["name"] != name:
            return False
        
        if restrict_cc is None:
            return True
        
        return e.get("compiler_configuration", None) == restrict_cc

    return list(filter(peek_filter, obj["targets"]))


def mfc_unique_target_exists(obj: dict, name: str, restrict_cc: str = None):
    return len(mfc_peek_target(obj, name, restrict_cc=restrict_cc)) == 1


def mfc_get_target(obj: dict, name: str, restrict_cc: str = None):
    matches = mfc_peek_target(obj, name, restrict_cc=restrict_cc)

    if len(matches) == 0:
        raise MFCException(f'Failed to retrieve dependency "{name}" in "{obj}" with restrict_cc="{restrict_cc}".')

    if len(matches) > 1:
        raise MFCException(f'More than one dependency to choose from for "{name}", defined in "{obj}" with restrict_cc="{restrict_cc}".')

    return matches[0]


def mfc_string_replace(mfc_state: MFCGlobalState, obj: dict, dependency_name: str, string: str, restrict_cc: str = None, recursive=True):
    dep = mfc_get_target(obj, dependency_name, restrict_cc=restrict_cc)

    regulars, mpis = mfc_state.conf["compilers"]["regular"], mfc_state.conf["compilers"]["mpi"]

    compiler_cfg = None
    for c_cfg in mfc_state.conf["compilers"]["configurations"]:
        if c_cfg["name"] == mfc_state.args["compiler_configuration"]:
            compiler_cfg = c_cfg
            break

    if compiler_cfg is None:
        raise MFCException(
            f'Failed to locate the compiler configuration "{mfc_state.args["compiler_configuration"]}".')

    if dep["type"] in ["clone", "download"]:
        install_path = f'{mfc_state.root_path}/{MFC_USE_SUBDIR}/{compiler_cfg["name"]}/build'
        source_path  = f'{mfc_state.root_path}/{MFC_USE_SUBDIR}/{compiler_cfg["name"]}/src/{dependency_name}'
    elif dep["type"] == "source":
        install_path = "ERR_INSTALL_PATH_IS_UNDEFINED"
        source_path  = dep["source"]["source"]
    else:
        raise MFCException(f'Unknown type "{dep["type"]}".')

    flags = compiler_cfg["flags"]

    replace_list = [
        ("${MFC_ROOT_PATH}",     mfc_state.root_path),
        ("${CONFIGURE_OPTIONS}", f'--prefix="{install_path}"'),
        ("${SOURCE_PATH}",       source_path),
        ("${INSTALL_PATH}",      install_path),
        ("${INSTALL_PATH}",      install_path),
        ("${MAKE_OPTIONS}",      f'-j {mfc_state.args["jobs"]}'),
        ("${COMPILER_FLAGS}",    f'CFLAGS="{flags.get("c", "")}" CPPFLAGS="{flags.get("c++", "")}" FFLAGS="{flags.get("fortran", "")}"'),
        ("${REGULAR_COMPILERS}", f'CC="{regulars.get ("c", "")}" CXX="{regulars.get  ("c++", "")}" FC="{regulars.get ("fortran", "")}"'),
        ("${MPI_COMPILERS}",     f'CC="{mpis.get     ("c", "")}" CXX="{mpis.get      ("c++", "")}" FC="{mpis.get     ("fortran", "")}"')
    ]

    for e in replace_list: string = string.replace(*e)

    if recursive:
        for dep2_info in obj["targets"]:
            string = string.replace("${" + dep2_info["name"] + ":", "${")
            string = mfc_string_replace(mfc_state, obj, dep2_info["name"], string, recursive=False)

    return string


def mfc_get_dependency_names(obj: dict, name: str, restrict_cc: str = None, recursive=False, visited: list = None):
    result: list = []

    if visited == None:
        visited = []

    if name not in visited:
        visited.append(name)

        desc = mfc_get_target(obj, name, restrict_cc=restrict_cc)

        for dependency_name in desc.get("depends", []):
            result.append(dependency_name)

            if recursive:
                result  += mfc_get_dependency_names(obj, dependency_name, recursive=recursive, visited=visited)
                visited += result

    return list(set(result))


def mfc_is_target_build_satisfied(mfc_state: MFCGlobalState, name: str):
    # Check if it hasn't been built before
    if len(mfc_peek_target(mfc_state.lock, name, mfc_state.args["compiler_configuration"])) == 0:
        return False

    # Retrive CONF & LOCK descriptors
    conf_desc = mfc_get_target(mfc_state.conf, name)
    lock_desc = mfc_get_target(mfc_state.lock, name, mfc_state.args["compiler_configuration"])

    # Check if it needs updating (LOCK & CONFIG descriptions don't match)
    if conf_desc["type"] != lock_desc["type"] or \
       conf_desc["type"] == "source"          or \
       conf_desc[conf_desc["type"]] != conf_desc[conf_desc["type"]]:
       return False

    # Check if any of its dependencies needs updating
    for dependency_name in mfc_get_dependency_names(mfc_state.conf, name, recursive=True):
        if not mfc_is_target_build_satisfied(mfc_state, dependency_name):
            return False

    # Check for "scratch" flag
    if mfc_state.args["scratch"]:
        return False

    return True


def mfc_build_target__clean_previous(mfc_state: MFCGlobalState, name: str):
    if not mfc_unique_target_exists(mfc_state.lock, name, mfc_state.args["compiler_configuration"]):
        return
    
    conf_desc = mfc_get_target(mfc_state.conf, name)
    lock_desc = mfc_get_target(mfc_state.lock, name, mfc_state.args["compiler_configuration"])
    
    if ((    conf_desc["type"] != lock_desc["type"]
         and lock_desc["type"] in ["clone", "download"]
        ) or (mfc_state.args["scratch"])):
        delete_directory_recursive_safe(f'{mfc_state.root_path}/{MFC_USE_SUBDIR}/{lock_desc["compiler_configuration"]}/src/{name}')


def mfc_build_target__fetch(mfc_state: MFCGlobalState, name: str, logfile):
    conf = mfc_get_target(mfc_state.conf, name)
    
    base_path = f'{mfc_state.root_path}/{MFC_USE_SUBDIR}/{mfc_state.args["compiler_configuration"]}'

    if conf["type"] in ["clone", "download"]:
        source_path = f'{base_path}/src/{name}'

        if conf["type"] == "clone":
            lock_matches = mfc_peek_target(mfc_state.lock, name, mfc_state.args["compiler_configuration"])
            
            if ((   len(lock_matches)    == 1
                and conf["clone"]["git"] != mfc_get_target(mfc_state.lock, name)["clone"]["git"], mfc_state.args["compiler_configuration"])
                or (mfc_state.args["scratch"])):
                clear_print(f'|--> Package {name}: GIT repository changed. Updating...', end='\r')

                delete_directory_recursive_safe(source_path)

            if not os.path.isdir(source_path):
                clear_print(f'|--> Package {name}: Cloning repository...', end='\r')
                
                execute_shell_command_safe(
                    f'git clone --recursive "{conf["clone"]["git"]}" "{source_path}" >> "{logfile.name}" 2>&1')

            clear_print(f'|--> Package {name}: Checking out {conf["clone"]["hash"]}...', end='\r')

            execute_shell_command_safe(
                f'cd "{source_path}" && git checkout "{conf["clone"]["hash"]}" >> "{logfile.name}" 2>&1')
        elif conf["type"] == "download":
            clear_print(f'|--> Package {name}: Removing previously downloaded version...', end='\r')

            delete_directory_recursive_safe(f'{base_path}/src/{name}')

            download_link = conf[conf["type"]]["link"].replace("${version}", conf[conf["type"]]["version"])
            filename = download_link.split("/")[-1]

            clear_print(f'|--> Package {name}: Downloading source...', end='\r')

            urllib.request.urlretrieve(download_link, f'{base_path}/src/{filename}')

            clear_print(f'|--> Package {name}: Uncompressing archive...', end='\r')

            uncompress_archive_to(f'{base_path}/src/{filename}',
                                  f'{base_path}/src/{name}')

            os.remove(f'{base_path}/src/{filename}')
    elif conf["type"] == "source":
        pass
    else:
        raise MFCException(f'Dependency type "{conf["type"]}" is unsupported.')


def mfc_build_target__build(mfc_state: MFCGlobalState, name: str, logfile):
    conf = mfc_get_target(mfc_state.conf, name)

    if conf["type"] not in ["clone", "download", "source"]:
        raise MFCException(f'Unknown target type "{conf["type"]}".')
    
    for command in conf["build"]:
        command = mfc_string_replace(mfc_state, mfc_state.conf, name, f"""\
cd "${{SOURCE_PATH}}" && \
PYTHON="python3" PYTHON_CPPFLAGS="$PYTHON_CPPFLAGS $(python3-config --includes) $(python3-config --libs)" \
bash -c '{command}' >> "{logfile.name}" 2>&1""")

        logfile.write(f'\n--- ./mfc.py ---\n{command}\n--- ./mfc.py ---\n\n')

        clear_print(f'|--> Package {name}: Building (Logging to {logfile.name})...', end='\r')

        execute_shell_command_safe(command)                 


def mfc_build_target__update_lock(mfc_state: MFCGlobalState, name: str):
    conf = mfc_get_target(mfc_state.conf, name)

    clear_print(f'|--> Package {name}: Updating lock file...', end='\r')

    new_entry = dict(conf, **{"compiler_configuration": mfc_state.args["compiler_configuration"]})

    if len(mfc_peek_target(mfc_state.lock, name, mfc_state.args["compiler_configuration"])) == 0:
        mfc_state.lock["targets"].append(new_entry)
    else:
        for index, dep in enumerate(mfc_state.lock["targets"]):
            if dep["name"] == name:
                mfc_state.lock["targets"][index] = new_entry

    mfc_dump_lock_file(mfc_state)


def mfc_build_target(mfc_state: MFCGlobalState, name: str):
    # Check if it needs to be (re)built
    if mfc_is_target_build_satisfied(mfc_state, name):
        clear_print(f'|--> Package {name}: Nothing to do ({colorama.Fore.GREEN}Success{colorama.Style.RESET_ALL})')
        return False
    
    # Build its dependencies
    for dependency_names in mfc_get_dependency_names(mfc_state.conf, name, recursive=False):
        mfc_build_target(mfc_state, dependency_names)

    clear_print(f'|--> Package {name}: Preparing build...', end='\r')

    with open(f'{mfc_state.root_path}/{MFC_USE_SUBDIR}/{mfc_state.args["compiler_configuration"]}/log/{name}.log', "w") as logfile:
        mfc_build_target__clean_previous(mfc_state, name)          # Clean any old build artifacts
        mfc_build_target__fetch         (mfc_state, name, logfile) # Fetch Source Code
        mfc_build_target__build         (mfc_state, name, logfile) # Build
        mfc_build_target__update_lock   (mfc_state, name)          # Update LOCK

    clear_print(f'|--> Package {name}: Done. ({colorama.Fore.GREEN}Success{colorama.Style.RESET_ALL})')

    return True


def mfc_environment_checks(mfc_state: MFCGlobalState):
    print("|--> Checking for the presence of required command-line utilities...", end='\r')

    required = ["python3", "python3-config", "make", "git"]
    required += mfc_state.conf["compilers"]["regular"].values()
    required += mfc_state.conf["compilers"]["mpi"].values()

    for index, utility in enumerate(required):
        clear_print(f"|--> {index+1}/{len(required)} Checking for {utility}...", end='\r')

        if shutil.which(utility) is None:
            raise MFCException(
                f'Failed to find the command line utility "{utility}". Please install it or make it visible.')

    # TODO: MacOS Checks
    if sys.platform == "darwin": # MacOS
        pass
        #mfc_state.conf["compilers"]["mpi"]["fortran"]
    
    clear_print(f"|--> Build environment: Passing. ({colorama.Fore.GREEN}Success{colorama.Style.RESET_ALL})")


def mfc_setup_directories(conf: dict):
    create_directory_safe(MFC_USE_SUBDIR)

    for d in ["src", "build", "log"]:
        for cc in [ cc["name"] for cc in conf["compilers"]["configurations"] ]:
            create_directory_safe(f"{MFC_USE_SUBDIR}/{cc}/{d}")


def mfc_main():
    try:
        for module in [("yaml", "pyyaml"), ("colorama", "colorama")]:
            import_module_safe(*module)

        colorama.init()

        conf = mfc_read_conf_file()
        mfc_setup_directories(conf)

        lock = mfc_read_lock_file()
        args = mfc_parse_arguments(conf)
        root_path = os.path.dirname(os.path.realpath(__file__))
        mfc_state = MFCGlobalState(conf=conf, lock=lock, args=args, root_path=root_path)

        mfc_print_header()
        mfc_environment_checks(mfc_state)

        if mfc_state.args["clean"]:
            delete_directory_recursive_safe(f"./{MFC_USE_SUBDIR}/")

            mfc_state.lock["targets"] = []

            mfc_dump_lock_file(mfc_state)

        for target_name in mfc_state.args["targets"]:
            mfc_build_target(mfc_state, target_name)

        mfc_dump_lock_file(mfc_state)
    except MFCException as exc:
        print(f"{colorama.Fore.RED}|--> {str(exc)}{colorama.Style.RESET_ALL}")
        exit(1)


if __name__ == "__main__":
    mfc_main()

