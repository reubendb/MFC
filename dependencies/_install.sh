#!/bin/bash

# Color ANSI Escape Sequences
FG_RED='\033[0;31m'
FG_GREEN='\033[0;32m'
FG_ORANGE='\033[0;33m'
FG_NONE='\033[0m'

# Exit with an error code if any command herein fails
set -e
set -o pipefail
set -o errtrace

echo 
echo -e $FG_ORANGE"Please ensure you have an appropriate build environment."$FG_NONE
echo -e $FG_ORANGE"Note: This file is meant to be called from within dependencies/."$FG_NONE

# Create some useful constants
dependencies_dir=$(pwd)
src_dir=$dependencies_dir"/src"
build_dir=$dependencies_dir"/build"
log_dir=$dependencies_dir"/log"

# Create required directories
mkdir -p src build log

# 0) Pre-Flight Safety
echo 
echo "-------------------------------------------------------"
echo "----------- Checking your build environment -----------"
echo "-------------------------------------------------------"
echo 

declare -a required_commands

required_commands[0]="make"
required_commands[1]="wget"
required_commands[2]="tar"
required_commands[3]="python3"
required_commands[4]="python3-config"

echo "+ Command Existance"

i=1
for required_command in "${required_commands[@]}"; do
    echo -e -n "+--+> ($i/${#required_commands[@]}) Checking existance of \"$required_command\"... "
    if ! command -v $required_command &> /dev/null; then
        echo -e $FG_RED"Not Found"$FG_NONE
        exit 1
    fi

    echo -e $FG_GREEN"Found"$FG_NONE

    i=$((i+1))
done

# If running on the Expanse supercomputer
if [[ $(env | grep -i 'expanse' | wc -c) -ne 0 ]]; then
    module load numactl
fi

# 1) Fetch

declare -a dependencies

dependencies[0]="FFTW3|http://www.fftw.org/fftw-3.3.10.tar.gz"
dependencies[1]="HDF5|https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.0/src/hdf5-1.12.0.tar.gz"
dependencies[2]="SILO|https://wci.llnl.gov/sites/wci/files/2021-09/silo-4.11.tgz"

echo 
echo "-------------------------------------------------------"
echo "---------------- Fetching Dependencies ----------------"
echo "-------------------------------------------------------"
echo 

cd $src_dir
    # For each dependency
    i=1
    for dependency in "${dependencies[@]}"; do    
        # Extract name and download link
        IFS="|" read -r -a arr <<< "${dependency}"

        name="${arr[0]}"
        link="${arr[1]}"

        echo -n "+--+> ($i/${#dependencies[@]}) "$name"... "

        # If we haven't downloaded it before (the directory named $name doesn't exist)
        if [ ! -d $(pwd)"/"$name ]; then
            archive_filename=$name".tar.gz"
            
            echo -e -n "\n|  |--> Downloading Source Code Archive... "
            wget -O $archive_filename -q $link
            echo -e $FG_GREEN"Success"$FG_NONE

            mkdir -p $name

            cd $name
                echo -n "|  |--> Uncompressing Downloaded Archive... "
                tar --strip-components 1 -xf "../"$archive_filename
                echo -e $FG_GREEN"Success"$FG_NONE
            cd ..
            
            rm $archive_filename
        else
            echo -e $FG_GREEN"Nothing to do"$FG_NONE
        fi

        i=$((i+1))
    done
cd ..

# 2) Add "build/<lib, bin, include, ...>" to the system's search path
#
# We append the export commands to every dotfile we can find because it's
# tricky to know which ones will be executed and when, especially if the
# user has multiple shells installed at the same time. There are also
# some complexities with login and non-login shell sessions.
#
# The code we add to the dotfiles has include guards to prevent
# errors and redundancy. Also, we only append to the dotfiles if
# the header guard isn't set to prevent duplicates.

echo 
echo -------------------------------------------------------
echo --------------- Exporting Library Paths ---------------
echo -------------------------------------------------------
echo 

found_dotfile_count=0
found_dotfile_list_string=""
if [[ -z "${MFC_ENV_SH_HEADER_GUARD}" ]]; then 
    export_cmds_0="export MFC_ENV_SH_HEADER_GUARD=\"SET\""
    export_cmds_1="export LD_LIBRARY_PATH=\"\$LD_LIBRARY_PATH:$build_dir/lib\""
    full_dotfile_string="\n# --- [Added by MFC | Start Section] --- #\nif [[ -z \"\${MFC_ENV_SH_HEADER_GUARD}\" ]]; then \n\t$export_cmds_0 \n\t$export_cmds_1 \nfi \n# --- [Added by MFC | End Section]   --- #\n"

    declare -a dotfile_paths

    dotfile_names[0]=".bashrc"
    dotfile_names[1]=".bash_profile"
    dotfile_names[2]=".bash_login"
    dotfile_names[3]=".profile"
    dotfile_names[4]=".zshrc"
    dotfile_names[5]=".cshrc"

    i=1

    for dotfile_name in "${dotfile_names[@]}"; do
        dotfile_path=$HOME"/"$dotfile_name

        echo -n "+--+> ($i/${#dotfile_names[@]}) $dotfile_path: "

        if [[ -a "$dotfile_path" ]]; then
            echo -e $FG_GREEN"Present."$FG_NONE
            echo -e $full_dotfile_string >> "$dotfile_path"

            found_dotfile_count=$((found_dotfile_count+1))
            
            if [ $i -ne "1" ]; then
                if [ $i -ne "${#dotfile_names[@]}" ]; then
                    found_dotfile_list_string="$found_dotfile_list_string, "
                else
                    found_dotfile_list_string="$found_dotfile_list_string, and "
                fi
            fi

            found_dotfile_list_string="$found_dotfile_list_string$dotfile_path"
        else
            echo -e $FG_ORANGE"Absent."$FG_NONE
        fi

        i=$((i+1))
    done

    if [ "$found_dotfile_count" -eq "0" ]; then
        echo -e "$FG_RED\n"
        echo "=================================================================================================="
        echo "| [ERROR] Could not find any dotfiles where we could export the path to the installed libraries. |"
        echo "=================================================================================================="
        echo -e "$FG_NONE"
        exit 1
    fi

    eval "$export_cmds_0"
    eval "$export_cmds_1"
else
    echo -e "+--+> "$FG_GREEN"The MFC header guard is already present."$FG_NONE$FG_ORANGE" No library path will be exported."$FG_NONE
fi

# 3) Build

echo 
echo -------------------------------------------------------
echo ---------------- Building Dependencies ----------------
echo -------------------------------------------------------
echo
echo -e $FG_ORANGE"Note: If any error occurs, please visit $log_dir/."$FG_NONE
echo -e $FG_ORANGE"Note: Building to $build_dir/."$FG_NONE
echo 

# FFTW3
echo "+--+> (1/3) FFTW3..."
log_filepath=$log_dir"/FFTW3.log"

cd $src_dir"/FFTW3"
    echo -n "|  |--> Configuring... "
    ./configure --prefix=$build_dir --enable-threads --enable-mpi > $log_filepath 2>&1
    echo -e $FG_GREEN"Success"$FG_NONE

    echo -n "|  |--> Building... " 
    make "$@" >> $log_filepath 2>&1
    echo -e $FG_GREEN"Success"$FG_NONE

    echo -n "|  |--> Installing... "
    make install >> $log_filepath 2>&1
    echo -e $FG_GREEN"Success"$FG_NONE

# HDF5
echo "+--+> (2/3) HDF5..."
log_filepath=$log_dir"/HDF5.log"

cd $src_dir"/HDF5"
    echo -n "|  |--> Configuring... "
    ./configure CC=mpicc CXX=mpicxx --prefix=$build_dir --enable-parallel --enable-deprecated-symbols > $log_filepath 2>&1
    echo -e $FG_GREEN"Success"$FG_NONE

    echo -n "|  |--> Building... "
    make "$@" >> $log_filepath 2>&1
    echo -e $FG_GREEN"Success"$FG_NONE

    echo -n "|  |--> Installing... "
    make install prefix=$build_dir >> $log_filepath 2>&1
    echo -e $FG_GREEN"Success"$FG_NONE

# SILO
echo "+--+> (3/3) SILO..."
log_filepath=$log_dir"/SILO.log"

cd $src_dir"/SILO"
    export PYTHON=python3
    export PYTHON_CPPFLAGS="$PYTHON_CPPFLAGS $(python3-config --includes) $(python3-config --libs)"

    # We use the following parameters
    # CC=mpicc CXX=mpicxx
    # So that "--with-hdf5" works...
    echo -n "|  |--> Configuring... "
    ./configure --prefix=$build_dir --enable-pythonmodule --enable-optimization \
                --disable-hzip      --disable-fpzip                             \
                FC=mpif90 F77=mpif77 CC=mpicc CXX=mpicxx                        \
                --with-hdf5=$build_dir"/include",$build_dir"/lib"               \
                --disable-silex > $log_filepath 2>&1
    echo -e $FG_GREEN"Success"$FG_NONE

    echo -n "|  |--> Building... "
    make "$@" >> $log_filepath 2>&1
    echo -e $FG_GREEN"Success"$FG_NONE

    echo -n "|  |--> Installing... "
    make install prefix=$build_dir >> $log_filepath 2>&1
    echo -e $FG_GREEN"Success"$FG_NONE

cd $src_dir

echo -e "\n$FG_GREEN"
echo "|-----------------------------------------------------|"
echo "|-------------- Completed Successfully ---------------|"
echo "|-----------------------------------------------------|"
echo -e "\n$FG_NONE"

if [ "$found_dotfile_count" -ne "0" ]; then
    echo -e "$FG_ORANGE\n[WARNING] MFC's dependency install script added code to $found_dotfile_count dotfiles ($found_dotfile_list_string) in order to correctly configure your environement variables (such as LD_LIBRARY_PATH). \n$FG_NONE"
    echo -e "$FG_GREEN""You are now in a new instance of your default shell ("$SHELL") and ready to build & run MFC! \n\n$FG_NONE"
    exec "$SHELL"
fi
