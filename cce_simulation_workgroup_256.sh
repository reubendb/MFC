jkKK#! /bin/bash

CCE_LLVM_PATH=${CRAY_CCE_CLANGSHARE}/../

## This was errantly 1024 by default in the generated binary
WGSIZE=256

## Turning bitcode into human-readable IR
echo "Disassembling"
${CCE_LLVM_PATH}/bin/llvm-dis build/simulation/simulation-cce-openmp-pre-llc.bc

## Find/replace the workgroup size to what it _should be_ in the (human readable) IR
echo "Globally setting amdgpu-flat-work-group-size size to 1,$WGSIZE"
sed "s/\"amdgpu-flat-work-group-size\"\=\"1,1024\"/\"amdgpu-flat-work-group-size\"\=\"1,${WGSIZE}\"/g" build/simulation/simulation-cce-openmp-pre-llc.ll > build/simulation/simulation-cce-openmp-pre-llc-wg${WGSIZE}.ll

## This is building to an AMD GPU code object, using internal copy/pasted Cray Fortran   flags.
## The flags may need to be adjusted for future compiler versions
echo "Invoking LLC to compile"
${CCE_LLVM_PATH}/bin/llc -mtriple=amdgcn-amd-amdhsa -disable-promote-alloca-to-lds -mcpu=gfx90a -amdgpu-dump-hsa-metadata build/simulation/simulation-cce-openmp-pre-llc-wg${WGSIZE}.ll -filetype=obj -o build/simulation/simulation-cce-openmp__llc_wg${WGSIZE}.amdgpu

## Wrapping AMD GPU code object in an object that CCE OpenACC runtime understands
echo "Linking to a CCE Offload module"
${CCE_LLVM_PATH}/bin/lld  -flavor gnu --no-undefined -shared -o build/simulation/simulation-wg${WGSIZE}.lld.exe build/simulation/simulation-cce-openmp__llc_wg${WGSIZE}.amdgpu

## Backend/hidden env var that tells the runtime where to find the offload object
echo "Now "
echo "export CRAY_ACC_MODULE=${PWD}/build/simulation/simulation-wg${WGSIZE}.lld.exe"
echo "to use the new GPU offload code."
echo "To use the original build"
echo "unset CRAY_ACC_MODULE"

# # This goes inside an sbatch script for multi-node submissions (> 1000 nodes important)
# # Requet an nvme via
# #SBATCH -C nvme
# # Put this in the sbatch script before execution via srun
# sbcast -pf ${PWD}/build/simulation/simulation-wg${WGSIZE}.lld.exe /mnt/bb/$USER/simulation-wg${WGSIZE}.lld.exe
# export CRAY_ACC_MODULE=/mnt/bb/$USER/simulation-wg${WGSIZE}.lld.exe"
