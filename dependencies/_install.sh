#!/bin/bash

# Exit with an error code if any command herein fails
set -e
set -o pipefail
set -o errtrace

# 1) Fetch

declare -a dependencies

dependencies[0]="FFTW3|http://www.fftw.org/fftw-3.3.10.tar.gz"
dependencies[1]="LAPACK|https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.10.0.tar.gz"
dependencies[2]="HDF5|https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.1/src/hdf5-1.12.1.tar.gz"
dependencies[3]="OPENMPI|https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.2.tar.gz"
dependencies[4]="SILO|https://wci.llnl.gov/sites/wci/files/2021-09/silo-4.11.tgz"

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
    cmake --build . -j --target install >> $log_filepath 2>&1

rm -rf $lapack_staging_dir_absolute

# HDF5
echo "+--+> HDF5"
log_filepath=$log_dir"/HDF5.log"

cd $src_dir"/HDF5"
    echo "|  |--> Configuring (./configure)..."
    ./configure --prefix=$build_dir --enable-parallel --enable-deprecated-symbols > $log_filepath 2>&1

    echo "|  |--> Building (w/ Make)..."
    make "$@" >> $log_filepath 2>&1

    echo "|  |--> Installing (w/ Make)..."
    make install prefix=$build_dir >> $log_filepath 2>&1

# OPENMPI
echo "+--+> OPENMPI"
log_filepath=$log_dir"/OPENMPI.log"

cd $src_dir"/OPENMPI"
    echo "|  |--> Configuring (./configure)..."
    ./configure --prefix=$build_dir > $log_filepath 2>&1

    echo "|  |--> Building (w/ Make)..."
    make "$@" >> $log_filepath 2>&1

    echo "|  |--> Installing (w/ Make)..."
    make install prefix=$build_dir >> $log_filepath 2>&1

# SILO
echo "+--+> SILO"
log_filepath=$log_dir"/SILO.log"

cd $src_dir"/SILO"
    export PYTHON_CPPFLAGS="$PYTHON_CPPFLAGS $(python3-config --cflags)"

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

echo 
echo "|-----------------------------------------------------|"
echo "|-------------- Completed Successfully ---------------|"
echo "|-----------------------------------------------------|"
echo 
