/*! \file   m_zfp_mpi.c
 *! \brief  Exposes a C interface for asynchronous MPI communication with ZFP compression.
 *! \author Henry Le Berre <hberre3@gatech.edu>
 */


#include <zfp.h>
#include <mpi.h>

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

void zfp_mpi_init(const config_t* const pConfig, worker_t* const pWorker) {
    if (!pConfig || !pWorker) return;
    if (!pConfig->pDoubles || !pConfig->nDoubles) return;

    pWorker->config = *pConfig;

    if (pConfig->bitfield & 1) {
        cuMalloc((void**)&pWorker->staging.pBytes, pWorker->staging.nBytes);
    } else {
        if (pConfig->bitfield & 2) {
            cudaHostAlloc((void**)&pWorker->staging.pBytes, pWorker->staging.nBytes, cudaHostAllocMapped);
        } else {
            pWorker->staging.pBytes = (uint8_t*)malloc(pWorker->staging.nBytes);
        }

        if (!pWorker->staging.pBytes) return;
    }

    //if (pConfig->bitfield & 4) {
    //    pWorker->pImpl = (impl_t*)malloc(sizeof(impl_t));
    //    if (!pWorker->pImpl) return;
//
    //    pWorker->pImpl->stream = zfp_stream_open(NULL);
    //    if (!pWorker->pImpl->stream) return;
//
    //    if (!zfp_stream_set_execution(pWorker->pImpl->stream, zfp_exec_cuda)) return;
    //} else {
    //    pWorker->pImpl = NULL;
    //}

    printf("Hello from zfp_mpi_init()!");

    //if (!(pConfig->bitfield & ZFP_MASK)) {
    //    pWorker->staging.nBytes = pConfig->nDoubles * sizeof(double);
    //}
}

void zfp_mpi_sendrecv(const worker_t* const pWorker) {
    MPI_Sendrecv(pWorker->config.pDoubles, pWorker->config.nDoubles*sizeof(double), MPI_BYTE, 0, 0,
                 pWorker->staging.pBytes,  pWorker->staging.nBytes,                 MPI_BYTE, 0, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}
