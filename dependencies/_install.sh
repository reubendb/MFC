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

# 1) Fetch

declare -a dependencies

dependencies[0]="FFTW3|http://www.fftw.org/fftw-3.3.10.tar.gz"
dependencies[1]="LAPACK|https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.10.0.tar.gz"
dependencies[2]="HDF5|https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.0/src/hdf5-1.12.0.tar.gz"
dependencies[3]="SILO|https://wci.llnl.gov/sites/wci/files/2021-09/silo-4.11.tgz"

dependencies_dir=$(pwd)
src_dir=$dependencies_dir"/src"
build_dir=$dependencies_dir"/build"
log_dir=$dependencies_dir"/log"

mkdir -p src build log

echo Please ensure you have an appropriate build environment
echo Note: This file is meant to be called from within dependencies/
echo 
echo "-------------------------------------------------------"
echo "---------------- Fetching Dependencies ----------------"
echo "-------------------------------------------------------"
echo 

cd $src_dir
    # For each dependency
    for dependency in "${dependencies[@]}"; do    
        # Extract name and download link
        IFS="|" read -r -a arr <<< "${dependency}"

        name="${arr[0]}"
        link="${arr[1]}"

        echo "+--+> "$name

        # If we haven't downloaded it before (the directory named $name doesn't exist)
        if [ ! -d $(pwd)"/"$name ]; then
            archive_filename=$name".tar.gz"
            
            echo "|  |--> Fetching "$name" @ "$link"."

            # Download Archive
            wget -O $archive_filename -q --show-progress $link

            mkdir -p $name

            cd $name
                echo "|  |--> Uncompressing "$archive_filename"."

                tar --strip-components 1 -xf "../"$archive_filename
            cd ..

            echo "|  |--> Removing "$archive_filename"."
            
            rm $archive_filename
        else
            echo "|  |--> Skipping "$name" because it has already been downloaded."
        fi
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
            echo "Present."
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
            echo "Absent."
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
    echo "+--+> The MFC header guard already present, no library path will be exported."
fi

# 3) Build

echo 
echo -------------------------------------------------------
echo ---------------- Building Dependencies ----------------
echo -------------------------------------------------------
echo Note: If any error occurs, please visit $log_dir/.
echo Note: Building to $build_dir/.
echo 

# FFTW3
echo "+--+> FFTW3"
log_filepath=$log_dir"/FFTW3.log"

cd $src_dir"/FFTW3"
    echo "|  |--> Configuring (w/ ./configure)..."
    ./configure --prefix=$build_dir --enable-threads --enable-mpi > $log_filepath 2>&1

    echo "|  |--> Building (w/ Make)..." 
    make "$@" >> $log_filepath 2>&1

    echo "|  |--> Installing (w/ Make)..."
    make install >> $log_filepath 2>&1

# LAPACK
echo "+--+> LAPACK"
log_filepath=$log_dir"/LAPACK.log"

lapack_staging_dir_absolute=$build_dir"/temp-LAPACK"

mkdir -p $lapack_staging_dir_absolute
cd $lapack_staging_dir_absolute
    echo "|  |--> Generating Build Files (w/ Cmake)..."
    cmake -DCMAKE_INSTALL_PREFIX=$build_dir               \
          -DCMAKE_INSTALL_LIBDIR=$build_dir"/lib"         \
          -DCMAKE_INSTALL_FULL_LIBDIR=$build_dir"/lib"    \
          -DCMAKE_INSTALL_INCLUDEDIR=$build_dir"/include" \
          -DCMAKE_INSTALL_BINDIR=$build_dir"/bin"         \
          $src_dir"/LAPACK" > $log_filepath 2>&1

    echo "|  |--> Building (w/ Cmake)..."
    cmake --build . --target install >> $log_filepath 2>&1

rm -rf $lapack_staging_dir_absolute

# HDF5
echo "+--+> HDF5"
log_filepath=$log_dir"/HDF5.log"

cd $src_dir"/HDF5"
    echo "|  |--> Configuring (./configure)..."
    ./configure CC=mpicc CXX=mpicxx --prefix=$build_dir --enable-parallel --enable-deprecated-symbols > $log_filepath 2>&1

    echo "|  |--> Building (w/ Make)..."
    make "$@" >> $log_filepath 2>&1

    echo "|  |--> Installing (w/ Make)..."
    make install prefix=$build_dir >> $log_filepath 2>&1

# SILO
echo "+--+> SILO"
log_filepath=$log_dir"/SILO.log"

cd $src_dir"/SILO"
    export PYTHON=python3
    export PYTHON_CPPFLAGS="$PYTHON_CPPFLAGS $(python3-config --cflags)"

    # We use the following parameters
    # CC=mpicc CXX=mpicxx
    # So that "--with-hdf5" works...
    echo "|  |--> Configuring (./configure)..."
    ./configure --prefix=$build_dir --enable-pythonmodule --enable-optimization \
                --disable-hzip      --disable-fpzip                             \
                FC=mpif90 F77=mpif77 CC=mpicc CXX=mpicxx                        \
                --with-hdf5=$build_dir"/include",$build_dir"/lib"               \
                --disable-silex > $log_filepath 2>&1

    echo "|  |--> Building (w/ Make)..."
    make "$@" >> $log_filepath 2>&1

    echo "|  |--> Installing (w/ Make)..."
    make install prefix=$build_dir >> $log_filepath 2>&1

echo -e "\n$FG_GREEN"
echo "|-----------------------------------------------------|"
echo "|-------------- Completed Successfully ---------------|"
echo "|-----------------------------------------------------|"
echo -e "\n$FG_NONE"

if [ "$found_dotfile_count" -ne "0" ]; then
    echo -e "$FG_ORANGE\n\n[WARNING] MFC's dependency install script added code to $found_dotfile_count dotfiles ($found_dotfile_list_string) in order to correctly configure your environement variables (such as LD_LIBRARY_PATH). Please start a new shell session (e.g exec \$SHELL) before running MFC or source one of the listed dotfiles (e.g source ~/.bashrc). \n\n$FG_NONE"
fi
