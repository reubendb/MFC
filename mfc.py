#!/usr/bin/env python3

# Imports.

import os                 # 
import re                 # Regular expressions.
import sys
import shutil             # Archive file handling.
import tarfile
import argparse           # Parsing command-line arguments.
import urllib.request     # 

class MFCException(Exception):
    pass

def ExecuteShellCommandSafe(cmd):
    if os.system(cmd) != 0:
        raise MFCException("Failed to execute command \"{}\".".format(cmd))

MFC_ROOT_PATH = os.path.dirname(os.path.realpath(__file__))

for module in [{ "import": "yaml",     "pip": "pyyaml"   },
               { "import": "colorama", "pip": "colorama" }]:
    while True:
        try:
            globals()[module["import"]] = __import__(module["import"])

            break
        except ImportError as exc:
            prompt = f"The Python package {exc.name} needs to be installed. Would you like to install it now? (yY/nN): "

            if input(prompt).upper().strip() == "Y":
                status = ExecuteShellCommandSafe("python3 -m pip install --user {}".format(module["pip"]))
                
                if status != 0:
                    exit(status)
            else:
                exit(1)

# Utility functions.

def FileLoadYAML(filepath: str):
    try:
        f = open(filepath, "r")
        d = yaml.safe_load(f)
        f.close()

        return d
    except (IOError, yaml.YAMLError) as exc:
        raise MFCException(f"Failed to load YAML from \"{filepath}\": {exc}")

def FileDumpYAML(filepath: str, data):
    try:
        f = open(filepath, "w")
        yaml.dump(data, f)
        f.close()
    except (IOError, yaml.YAMLError) as exc:
        raise MFCException(f"Failed to dump YAML from \"{filepath}\": {exc}")

def SimplifyDependencyName(name: str):
    return re.sub(r"[0-9]+", "", name.strip())

def UncompressArchiveTo(archive_filepath, destination):
    archive = tarfile.open(archive_filepath, "r")
    archive.extractall("/".join(destination.rstrip("/").split("/")[:-1]))
    archive.close()

    # WARNING: Technicall doesn't always work
    # TODO: Better solution
    filename_wo_ext = archive_filepath.replace(".tar.gz", "")

    os.rename(filename_wo_ext, destination)

def DeleteDictoryRecursiveSafe(dirpath):
    if os.path.isdir(dirpath):
        shutil.rmtree(dirpath)

# Good Luck. © Henry Le Berre :)
CenterBlockOfText=lambda text: "\n".join(("{0}"+line+"{0}").format(" "*((shutil.get_terminal_size((80, 20)).columns - max([len(re.compile("|".join([re.escape(k) for k in [colorama.Style.RESET_ALL] + [ getattr(colorama.Fore, key) for key in dir(colorama.Fore) if not callable(getattr(colorama.Fore, key)) and not key.startswith("__") ]]), flags=re.DOTALL).sub("", line)) for line in text.splitlines()])) // 2)) for line in text.splitlines())

# Main MFC functions.

def MFC_ReadConfigurationFile():
    return FileLoadYAML("mfc.cfg.yaml")

def MFC_ReadLockFile():
    if not os.path.exists("mfc.lock.yaml"):
        f = open("mfc.lock.yaml", 'w')
        f.write("""\
targets: []
""")
        f.close()

    return FileLoadYAML("mfc.lock.yaml")

def MFC_WriteLockFile():
    FileDumpYAML("mfc.lock.yaml", MFC_CUR)

def MFC_ParseArguments():    
    parser = argparse.ArgumentParser(description="Wecome to the MFC master script.",)

    general   = parser.add_argument_group(title="General")
    compiling = parser.add_argument_group(title="Compiling")

    compiler_target_names = [ e["name"] for e in MFC_CFG["targets"] ]
    compiling.add_argument("-t", "--targets", nargs="+", type=str,
                           choices=compiler_target_names, default="",
                           help="The space-separated targets you wish to have built.")

    compiler_configuration_names = [ e["name"] for e in MFC_CFG["compilers"]["configurations"] ]
    compiling.add_argument("-cc", "--compiler-configuration", type=str,
                           choices=compiler_configuration_names, default=compiler_configuration_names[1])

    compiling.add_argument("-j", "--jobs", metavar="N", type=int,
                           help="Allows for N concurrent jobs.", default=1)

    compiling.add_argument("-c", "--clean", action="store_true",
                           help="Cleans all binaries and build artificats.")

    general.add_argument("-f", "--force", action="store_true",
                         help="Force the operation.")

    return vars(parser.parse_args())

def MFC_PrintHeader():
    print(CenterBlockOfText("""\

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

def MFC_Get_Dependency(name, obj):
    matches = list(filter(lambda x: x["name"]==name, obj["targets"]))

    if len(matches) != 1:
        raise MFCException(f"Failed to find dependency \"{name}\" in \"{obj}\".")

    return matches[0]

def MFC_String_Replace(string, dependency_name, obj, recursive=True):
    dep = MFC_Get_Dependency(dependency_name, obj)

    install_path = ""
    source_path  = ""
    if dep["type"] in ["clone", "download"]:
        install_path = MFC_ROOT_PATH + "/dependencies/build"
        source_path  = MFC_ROOT_PATH + "/dependencies/src/" + dependency_name
    elif dep["type"] == "source":
        install_path = "ERR_INSTALL_PATH_IS_UNDEFINED"
        source_path  = dep["source"]["source"]
    else:
        raise MFCException("hmmm....")

    regulars, mpis = MFC_CFG["compilers"]["regular"], MFC_CFG["compilers"]["mpi"]

    compiler_cfg = None
    for c_cfg in MFC_CFG["compilers"]["configurations"]:
        if c_cfg["name"] == MFC_ARGS["compiler_configuration"]:
            compiler_cfg = c_cfg
            break
    
    if compiler_cfg == None:
        raise MFCException("Failed to locate the compiler configuration \"{}\".".format(MFC_ARGS["compiler_configuration"]))

    compiler_flags = compiler_cfg["flags"]

    repl_l = [
        ("${MFC_ROOT_PATH}", MFC_ROOT_PATH), ("${CONFIGURE_OPTIONS}", f"--prefix=\"{install_path}\""),
        ("${SOURCE_PATH}",   source_path),   ("${INSTALL_PATH}",      install_path),
        ("${INSTALL_PATH}",  install_path),  ("${MAKE_OPTIONS}",      "-j {}".format(MFC_ARGS["jobs"])),
        ("${COMPILER_FLAGS}",                "CFLAGS=\"{}\" CPPFLAGS=\"{}\" FFLAGS=\"{}\"".format(compiler_flags["c"], compiler_flags["c++"], compiler_flags["fortran"])),
        ("${REGULAR_COMPILERS}", "CC=\"{}\" CXX=\"{}\" FC=\"{}\"".format(regulars["c"], regulars["c++"], regulars["fortran"])),
        ("${MPI_COMPILERS}",     "CC=\"{}\" CXX=\"{}\" FC=\"{}\"".format(mpis["c"],     mpis["c++"],     mpis["fortran"]))
    ]

    for a,b in repl_l:
        string = string.replace(a, b)

    if recursive:
        for dep2_info in obj["targets"]:
            string = string.replace("${"+dep2_info["name"]+":", "${")
            string = MFC_String_Replace(string, dep2_info["name"], obj, recursive=False)

    return string

def MFC_Was_Dependency_Updated(cfg_info, cur_info):
    if cfg_info["type"] != cur_info["type"]:
        return True
    
    if cfg_info["type"] == "source":
        return True

    t = cfg_info["type"]
    if cfg_info[t] != cur_info[t]:
        return True

    return False

def MFC_Build_Dependency(name):
    cfg = MFC_Get_Dependency(name, MFC_CFG)

    # Create Basic Directories
    for e in ["src", "build", "log"]:
        os.makedirs(MFC_ROOT_PATH + "/dependencies/"+e, exist_ok=True)

    # Build the dependencies it depends on
    bHadToBuildADependency=False
    for depends_name in cfg.get("depends", []):
        bHadToBuildADependency |= MFC_Build_Dependency(depends_name)

    cur            = None
    cur_matches    = list(filter(lambda x: x["name"]==name, MFC_CUR["targets"]))
    bWasNeverBuilt = len(cur_matches) != 1

    # TODO: desciption
    bDifferentCFGandCUR = False
    if len(cur_matches) == 1:
        cur = cur_matches[0]
        
        bDifferentCFGandCUR = MFC_Was_Dependency_Updated(cfg, cur)

        # Remove old directories from previous build/download type
        if cfg["type"] != cur["type"]:
            if cur["type"] in ["clone", "download"]:
                DeleteDictoryRecursiveSafe(MFC_ROOT_PATH + "/dependencies/src/" + name)

    if not (bWasNeverBuilt or bHadToBuildADependency or bDifferentCFGandCUR):
        return False

    # Create Log File
    log_filepath = MFC_ROOT_PATH + "/dependencies/log/" + name + ".log"
    if os.path.isfile(log_filepath):
        os.remove(log_filepath)

    # 
    if cfg["type"] == "clone" or cfg["type"] == "download":
        source_path = MFC_ROOT_PATH + "/dependencies/src/" + name

        if cfg["type"] == "clone":
            if cur != None and os.path.isdir(source_path):
                if cfg["clone"]["git"] != cur["clone"]["git"]:
                    DeleteDictoryRecursiveSafe(source_path)

            if not os.path.isdir(source_path):
                ExecuteShellCommandSafe("git clone --recursive \"{}\" \"{}\" >> \"{}\" 2>&1".format(cfg["clone"]["git"], source_path, log_filepath))
            
            ExecuteShellCommandSafe("cd \"{}\" && git checkout \"{}\" >> \"{}\" 2>&1".format(source_path, cfg["clone"]["hash"], log_filepath))   
        elif cfg["type"] == "download":
            DeleteDictoryRecursiveSafe(MFC_ROOT_PATH + "/dependencies/src/" + name)
            
            download_link = cfg[cfg["type"]]["link"].replace("${version}", cfg[cfg["type"]]["version"])
            filename      = download_link.split("/")[-1]

            urllib.request.urlretrieve(download_link, MFC_ROOT_PATH + "/dependencies/src/" + filename)

            UncompressArchiveTo(MFC_ROOT_PATH + "/dependencies/src/" + filename, MFC_ROOT_PATH + "/dependencies/src/" + name)

            os.remove(MFC_ROOT_PATH + "/dependencies/src/" + filename)
    elif cfg["type"] == "source":
        pass
    else:
        raise MFCException("Dependency type \"{}\" is unsupported.", cfg["type"])

    # Build
    if cfg["type"] in ["clone", "download", "source"]:
        for command in cfg["build"]:
            if cfg["type"] in ["clone", "download", "source"]:
                command = MFC_String_Replace(f"""\
cd "${{SOURCE_PATH}}" && \
PYTHON="python3" PYTHON_CPPFLAGS="$PYTHON_CPPFLAGS $(python3-config --includes) $(python3-config --libs)" \
bash -c '{command}' >> "{log_filepath}" 2>&1
""", name, MFC_CFG)

                print(command)
                ExecuteShellCommandSafe(command)
            else:
                raise MFCException("Dependency type \"{}\" is unsupported.", cfg["type"])
    else:
        raise MFCException("LKKLIAEIOAIAZEPIOKLQ...")

    if cfg["type"] == "source":
        return True

    # Update CUR
    if len(cur_matches) == 0:
        MFC_CUR["targets"].append(cfg)
    else:
        for dep in MFC_CUR["targets"]:
            if dep["name"] == name:
                dep = cfg
    
    MFC_WriteLockFile()
    
    return True

def MFC_PreflightChecks():
    REQUIRED_COMMANDS = ["python3", "python3-config", "make", "git"]
    REQUIRED_COMMANDS += MFC_CFG["compilers"]["regular"].values()
    REQUIRED_COMMANDS += MFC_CFG["compilers"]["mpi"].values()

    for utility in REQUIRED_COMMANDS:
        if shutil.which(utility) == None:
            raise MFCException("Failed to find the command line utility \"{}\". Please install it or make it visible.".format(utility))

    #TODO: MacOS CHECKS
    if sys.platform == "darwin": # MacOS
        MFC_CFG["compilers"]["mpi"]["fortran"]

# Code.

colorama.init()

try:
    MFC_CFG  = MFC_ReadConfigurationFile()
    MFC_CUR  = MFC_ReadLockFile()
    MFC_ARGS = MFC_ParseArguments()

    MFC_PrintHeader()
    MFC_PreflightChecks()

    if MFC_ARGS["clean"]:
        DeleteDictoryRecursiveSafe("./dependencies/")

        MFC_CUR["targets"] = []

        MFC_WriteLockFile()

    for target_name in MFC_ARGS["targets"]:
        MFC_Build_Dependency(target_name)

    MFC_WriteLockFile()
except MFCException as exc:
    print(colorama.Fore.RED + str(exc) + colorama.Style.RESET_ALL)
    exit(1)
