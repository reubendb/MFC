/*! \file   c_compress.cpp
 *! \brief  Exposes an implementation-agnostic C interface for compression.
 *! \author Henry Le Berre <hberre3@gatech.edu>
 */

#ifdef MFC_ZFP
#include <zfp.h>
#endif // MFC_ZFP

#ifdef MFC_MPI
#include <mpi.h>
#endif // MFC_MPI

#if defined(__CUDACC__) || defined(__NVCOMPILER)
#define MFC_CUDA

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#endif // defined(__CUDACC__) || defined(__NVCOMPILER)

#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#define C_TYPE double

#ifdef MFC_ZFP
#define ZFP_TYPE zfp_type_double
#endif // MFC_ZFP

/* Enumeration for memory location. */
typedef enum {
    location_cpu = 0, /* Host   */
    location_gpu = 1  /* Device */
} e_location;


/* Struct holding internal pState-> It is implementation defined and
   may vary upon ifdefs and compression backend. */
typedef struct {
#ifdef MFC_ZFP
    zfp_field*  field;  /* ZFP     field. */
    zfp_stream* stream; /* ZFP    stream. */
    bitstream*  bits;   /* ZFP bitstream. */
#endif // MFC_ZFP
} impl_t;


typedef struct {
    bool       bZFP;  /* True to enable compression. */
    C_TYPE*    pBuff; /* Buffer of C_TYPEs to (de)compress. */
    size_t     count; /* Number of C_TYPEs to (de)compress. */    
    double     rate;  /* (Actual) ZFP Rate. */
    e_location from;  /* Location of the uncompressed buffer. */
    e_location to;    /* Location of the   compressed buffer. */
} config_t;


/* Struct holding compression state and information. */
typedef struct {
    /* Bookkeeping items. */
    config_t config; // (De)Compression configuration.
    impl_t*  pImpl;  // Pointer to an opaque object holding (internal) state.

    /* Staging buffer information. */
    struct {
        size_t   bytes; // (Maximum) Size of the buffer in bytes.
        uint8_t* pBuff; // Buffer holding the staging data.
    } staging;
} state_t;


bool c_compress_init(
    const config_t* const pConfig,
          state_t * const pState
) {
    pState->config = *pConfig;

    pState->staging.pBuff = NULL;
    pState->staging.bytes = 0;

    pState->pImpl = (impl_t*)malloc(sizeof(impl_t));

    if (!pState->pImpl) {
        printf("[c_compress_init] Error: malloc failed.");
        return false;
    }

#ifdef MFC_ZFP

    pState->pImpl->stream = zfp_stream_open(NULL);

#ifdef MFC_CUDA

    if (!zfp_stream_set_execution(pState->pImpl->stream, zfp_exec_cuda)) {
        printf("[c_compress_init] Error: zfp_exec_cuda not available.");

        free(pState->pImpl);

        return false;
    }

    printf("[c_compress_init] Cuda is enabled.\n");

#endif // MFC_CUDA

    pState->config.rate = zfp_stream_set_rate(pState->pImpl->stream, pConfig->rate, ZFP_TYPE, 1, false);

    printf("[zfp] Rate is %f (requested %f).\n", pConfig->rate, pConfig->rate);

    pState->pImpl->field  = zfp_field_1d(pConfig->pBuff, ZFP_TYPE, pConfig->count);
    pState->staging.bytes = zfp_stream_maximum_size(pState->pImpl->stream, pState->pImpl->field);
#else
    pState->staging.bytes = pState->config.count * sizeof(C_TYPE);
#endif // MFC_ZFP

    // Set pBytes
#ifdef MFC_CUDA
    if (to == location_device) {
        const cudaError_t cu_err = cudaMalloc((void**)&pState->staging.pBuff, pState->staging.bytes);

        if (cu_err != cudaSuccess) {
            printf("[c_compress_init] Error: cudaMalloc returned %d.", (int)cu_err);

            zfp_field_free(pState->pImpl->field);
            free(pState->pImpl);

            return false;
        }
    } else { 
#endif

        pState->staging.pBuff = (uint8_t*)malloc(pState->staging.bytes);
#ifdef MFC_CUDA
    }
#endif

#ifdef MFC_ZFP
    pState->pImpl->bits = stream_open(pState->staging.pBuff, pState->staging.bytes);

    zfp_stream_set_bit_stream(pState->pImpl->stream, pState->pImpl->bits);
#endif // MFC_ZFP

    return true;
}


size_t c_compress(state_t* const pState) {
#ifdef MFC_ZFP
    zfp_stream_rewind(pState->pImpl->stream);

    return zfp_compress(pState->pImpl->stream, pState->pImpl->field);
#else
    printf("WARNING: c_compress is not implemented.");

    return 0;
#endif
}


size_t c_decompress(state_t* const pState) {
#ifdef MFC_ZFP
    zfp_stream_rewind(pState->pImpl->stream);

    return zfp_decompress(pState->pImpl->stream, pState->pImpl->field);
#else
    printf("WARNING: c_compress is not implemented.");
    
    return 0;
#endif
}


void c_compress_finalize(state_t* const pState) {
#ifdef MFC_CUDA
    if (pState->to == location_device) {
        cudaFree(pState->pBytes);
    } else {
#else
        free(pState->staging.pBuff);
#endif
#ifdef MFC_CUDA
    }
#endif

#ifdef MFC_ZFP
    stream_close(pState->pImpl->bits);
    zfp_stream_close(pState->pImpl->stream);
    zfp_field_free(pState->pImpl->field);
#endif // MFC_ZFP

    free(pState->pImpl);
}




#define BENCH_N        10000
#define ZFP_BENCH_RATE (8.f*sizeof(C_TYPE)) / 2.f




void c_compress_bench() {
    /*
    C_TYPE* buffer        = (C_TYPE*)malloc(sizeof(C_TYPE)*BENCH_N);
    C_TYPE* buffer_origin = (C_TYPE*)malloc(sizeof(C_TYPE)*BENCH_N);

    srand(1);
    for (int i = 0; i < BENCH_N; ++i) {
        buffer[i] = (2 * (rand()/(C_TYPE)RAND_MAX)) - 1;
    }
    printf("\n");

    memcpy(buffer_origin, buffer, BENCH_N*sizeof(C_TYPE));

    state_t zfp = c_compress_init(buffer, BENCH_N, ZFP_BENCH_RATE, 0, 0);

    printf("c_compress:   %d\n", (int)c_compress(&zfp));
    memset(buffer, 0, BENCH_N*sizeof(C_TYPE));
    printf("c_decompress: %d\n", (int)c_decompress(&zfp));


    C_TYPE avg_abs = 0; C_TYPE min_abs = -1; C_TYPE max_abs = 0;
    C_TYPE avg_rel = 0; C_TYPE min_rel = -1; C_TYPE max_rel = 0;

    for (int i = 0; i < BENCH_N; ++i) {
        C_TYPE abs = buffer_origin[i] - buffer[i];
        C_TYPE rel = 0; 

        if (buffer_origin[i] != buffer[i]) {
            // assumes buffer_origin[i] != 0
            rel = abs / buffer_origin[i];
        }

        if (rel < 0) { rel = -rel; }
        if (abs < 0) { abs = -abs; }

        avg_abs += abs;
        avg_rel += rel;

        if (rel > max_rel) {
            max_rel = rel;
        }

        if (abs > max_abs) {
            max_abs = abs;
        }

        if (abs < min_abs || min_abs < 0) {
            min_abs = abs;
        }

        if (rel < min_rel || min_rel < 0) {
            min_rel = rel;
        }
    }

    for (int i = 0; i < BENCH_N; ++i) {
        //printf("%.15f -> %.15f\n", buffer_origin[i], buffer[i]);
    }

    printf("Absolute Error:\n");
    printf("  - Average: %.15e\n", avg_abs / BENCH_N);
    printf("  - Minimum: %.15e\n", min_abs);
    printf("  - Maximum: %.15e\n", max_abs);

    printf("Relative Error:\n");
    printf("  - Average: %.15e\n", avg_rel / BENCH_N);
    printf("  - Minimum: %.15e\n", min_rel);
    printf("  - Maximum: %.15e\n", max_rel);

    c_compress_finalize(&zfp);
    */
}
