/*! \file   c_zfp_mpi.c
 *! \brief  Exposes a C interface for asynchronous MPI communication with ZFP compression.
 *! \author Henry Le Berre <hberre3@gatech.edu>
 */


#ifdef MFC_ZFP
#include <zfp.h>
#endif

#ifdef MFC_MPI
#include <mpi.h>
#endif

#if defined(__CUDACC__) || defined(__NVCOMPILER)
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#endif

#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#define FROM_MASK (0b001)
#define TO_MASK   (0b010)
#define ZFP_MASK  (0b100)

typedef struct {
    uint8_t bitfield;
    double* pDoubles;
    size_t  nDoubles;
    double  rate;
} config_t;

typedef struct {
    uint8_t* pBytes;
    size_t   nBytes;
} staging_t;

typedef struct {
    zfp_field*  field;  /* ZFP     field. */
    zfp_stream* stream; /* ZFP    stream. */
    bitstream*  bits;   /* ZFP bitstream. */
} impl_t;

typedef struct {
    config_t  config;
    staging_t staging;
    impl_t*   pImpl;
} worker_t;

void c_zfp_mpi_init(const config_t* const pConfig, worker_t* const pWorker) {
    if (!pConfig || !pWorker) return;
    if (!pConfig->pDoubles || !pConfig->nDoubles) return;

    pWorker->config = *pConfig;

    if (pConfig->bitfield & TO_MASK) {
        cuMalloc((void**)&pWorker->staging.pBytes, pWorker->staging.nBytes);
    } else {
        if (pConfig->bitfield & FROM_MASK) {
            cudaHostAlloc((void**)&pWorker->staging.pBytes, pWorker->staging.nBytes, cudaHostAllocMapped);
        } else {
            pWorker->staging.pBytes = (uint8_t*)malloc(pWorker->staging.nBytes);
        }

        if (!pWorker->staging.pBytes) return;
    }

    if (pConfig->bitfield & ZFP_MASK) {
        pWorker->pImpl = (impl_t*)malloc(sizeof(impl_t));
        if (!pWorker->pImpl) return;

        pWorker->pImpl->stream = stream_open(NULL);
        if (!pWorker->pImpl->stream) return;

        if (!zfp_stream_set_execution(pWorker->pImpl->stream, zfp_exec_cuda)) return;
    } else {
        pWorker->pImpl = NULL;
    }

    //if (!(pConfig->bitfield & ZFP_MASK)) {
    //    pWorker->staging.nBytes = pConfig->nDoubles * sizeof(double);
    //}
}
