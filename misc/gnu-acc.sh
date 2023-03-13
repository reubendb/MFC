#!/bin/bash

# Authors:     Henry Le Berre (GaTech) & Alex Woods (GaTech).
# Description: Build a the GNU compiler suite with amdgcn-amdhsa offloading
#              from the devel/omp/gcc-12 branch.

set -e
set -x

work_dir="${TMPDIR:=/tmp}/temp-gcc-acc-work/build-gcc-amdgpu/"
build_dir="$work_dir/builds"
install_dir="$HOME/gcc-acc/"

mkdir -p "$work_dir" && cd "$work_dir"
mkdir -p builds

wget https://github.com/llvm/llvm-project/archive/refs/tags/llvmorg-13.0.1.tar.gz

mkdir -p llvm-project

tar -vxf llvmorg-13.0.1.tar.gz -C llvm-project --strip-components 1

rm llvmorg-13.0.1.tar.gz
 
llvmsrcdir="$work_dir/llvm-project"

mkdir -p "$build_dir/llvm" && cd "$build_dir/llvm"

cmake -D 'LLVM_TARGETS_TO_BUILD=X86;AMDGPU' \
      -D LLVM_ENABLE_PROJECTS=lld           \
      -G Ninja $llvmsrcdir/llvm
ninja

mkdir -p "$install_dir/amdgcn-amdhsa/bin"
cp -a bin/llvm-ar "$install_dir/amdgcn-amdhsa/bin/ar"
cp -a bin/llvm-ar "$install_dir/amdgcn-amdhsa/bin/ranlib"
cp -a bin/llvm-mc "$install_dir/amdgcn-amdhsa/bin/as"
cp -a bin/llvm-nm "$install_dir/amdgcn-amdhsa/bin/nm"
cp -a bin/lld     "$install_dir/amdgcn-amdhsa/bin/ld"

# Clone gcc, newlib and link newlib in
cd $work_dir
git clone --depth=1 --branch master           git://sourceware.org/git/newlib-cygwin.git newlib
git clone --depth=1 --branch devel/omp/gcc-12 git://gcc.gnu.org/git/gcc.git              gcc
cd gcc
contrib/download_prerequisites
ln -s ../newlib/newlib newlib
cd ..
target=$(gcc/config.guess)

# Build offloading GCC
mkdir -p "$build_dir/build-offload-gcc"
cd       "$build_dir/build-offload-gcc"

$work_dir/gcc/configure \
  --target=amdgcn-amdhsa    --enable-languages=c,lto,fortran \
  --disable-sjlj-exceptions --with-newlib                    \
  --disable-libquadmath     --prefix="$install_dir"          \
  --enable-as-accelerator-for=x86_64-pc-linux-gnu            \
  --with-build-time-tools="$install_dir/amdgcn-amdhsa/bin"

make -j $(nproc) && make install

cd .. && rm "$work_dir/gcc/newlib"

# Build host GCC
mkdir build-host-gcc
cd    build-host-gcc
$work_dir/gcc/configure \
    --build=x86_64-pc-linux-gnu --host=x86_64-pc-linux-gnu              \
    --target=x86_64-pc-linux-gnu                                        \
    --enable-offload-targets=amdgcn-amdhsa="$install_dir/amdgcn-amdhsa" \
    --enable-languages="c,c++,fortran,lto" \
    --disable-multilib --prefix="$install_dir"

make -j $(nproc) && make install

cd ..
