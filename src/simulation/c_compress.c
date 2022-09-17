/*! \file   c_compress.c
 *! \brief  Exposes an implementation-agnostic C interface for compression.
 *! \author H. Le Berre
 */

#include <zfp.h>


#if defined(__CUDACC__) || defined(__NVCOMPILER)

#   include <cuda.h>
#   include <cuda_runtime.h>
#   include <device_launch_parameters.h>

#   define MFC_CUDA

#endif // defined(__CUDACC__) || defined(__NVCOMPILER)

#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#define C_TYPE   double
#define ZFP_TYPE zfp_type_double

#define BENCH_N        10000
#define ZFP_BENCH_RATE (8.f*sizeof(C_TYPE)) / 2.f

/* Enumeration for memory location. */
typedef enum c_compress_loc {
    c_compress_loc_host   = 0, /* Host   (e.g CPU). */
    c_compress_loc_device = 1  /* Device (e.g GPU). */
} c_compress_loc_t; // c_compress_loc_t


/* Struct holding internal state. It is implementation defined and
   may vary upon ifdefs and compression backend. */
typedef struct {
    /* ZFP internal types. */
    zfp_field*  field;  /* ZFP     field.  */
    zfp_stream* stream; /* ZFP    stream. */
    bitstream*  bits;   /* ZFP bitstream. */
} c_compress_internal_t;


/* Struct holding compression state and information. */
typedef struct {
    bool bActive; /* True if c_compress_init was successful. */

    /* User Configuration. */
    C_TYPE*          pDoubles; /* Buffer of doubles to compress from & decompress to. */
    size_t           nDoubles; /* Number of doubles to compress from & decompress to. */    
    double           rate;     /* (Actual) ZFP Rate. */
    c_compress_loc_t from;     /* Location of the uncompressed buffer. */
    c_compress_loc_t to;       /* Location of the   compressed buffer. */

    /* Compressed buffer information. */
    uint8_t* pBytes; /* Buffer holding the compressed data. */
    size_t   nBytes; /* (Maximum) Size of pBytes (in bytes). */

    /* Internal State. */
    c_compress_internal_t* pInternal;
} c_compress_state_t;


c_compress_state_t c_compress_init(
    const C_TYPE*          pDoubles, /* Buffer of doubles to compress from & decompress to. */
    const size_t           nDoubles, /* Number of doubles to compress from & decompress to. */
    const double           rate,     /* ZFP (de)compression rate. */
    const c_compress_loc_t from,     /* Location of the uncompressed buffer. */
    const c_compress_loc_t to        /* Location of the   compressed buffer. */
) {
    c_compress_state_t state;

    state.bActive = false;

    state.pDoubles = pDoubles;
    state.nDoubles = nDoubles;
    state.from     = from;
    state.to       = to;

    state.pBytes = NULL;
    state.nBytes = 0;

    state.pInternal = (c_compress_internal_t*)malloc(sizeof(c_compress_internal_t));
    state.pInternal->stream = zfp_stream_open(NULL);

#ifdef MFC_CUDA

    if (!zfp_stream_set_execution(state.pInternal->stream, zfp_exec_cuda)) {
        printf("[c_compress_init] Error: zfp_exec_cuda not available.");
        return state;
    }

    printf("[c_compress_init] Cuda is enabled.\n");

#endif // MFC_CUDA

    state.rate = zfp_stream_set_rate(state.pInternal->stream, rate, ZFP_TYPE, 1, false);

    printf("[zfp] Rate is %f (requested %f).\n", state.rate, rate);

    state.pInternal->field = zfp_field_1d(pDoubles, ZFP_TYPE, nDoubles);
    state.nBytes = zfp_stream_maximum_size(state.pInternal->stream, state.pInternal->field);

    // Set pBytes
#ifdef MFC_CUDA
    if (to == c_compress_loc_device) {
        const cudaError_t cu_err = cudaMalloc((void**)&state.pBytes, state.nBytes);

        if (cu_err != cudaSuccess) {
            zfp_field_free(state.pInternal->field);
            printf("[c_compress_init] Error: cudaMalloc returned %d.", (int)cu_err);
            return state;
        }
    } else { 
#endif
        state.pBytes = (uint8_t*)malloc(state.nBytes);
#ifdef MFC_CUDA
    }
#endif

    state.pInternal->bits = stream_open(state.pBytes, state.nBytes);

    zfp_stream_set_bit_stream(state.pInternal->stream, state.pInternal->bits);

    state.bActive = true;

    return state;
}


size_t c_compress(c_compress_state_t* const pState) {
    zfp_stream_rewind(pState->pInternal->stream);
    const size_t offset = zfp_compress(pState->pInternal->stream, pState->pInternal->field);

    return offset;
}


size_t c_decompress(c_compress_state_t* const pState) {
    zfp_stream_rewind(pState->pInternal->stream);

    const size_t offset = zfp_decompress(pState->pInternal->stream, pState->pInternal->field);

    printf("c_decompress = %d\n", (int)offset);

    return offset;
}


void c_compress_finalize(c_compress_state_t* const pState) {
    if (pState->bActive) {
#ifdef MFC_CUDA
        if (pState->to == c_compress_loc_device) {
            cudaFree(pState->pBytes);
        } else {
#else
            free(pState->pBytes);
#endif
#ifdef MFC_CUDA
        }
#endif

        stream_close(pState->pInternal->bits);
        zfp_stream_close(pState->pInternal->stream);
        zfp_field_free(pState->pInternal->field);

        free(pState->pInternal);

        pState->bActive = false;
    }
}


void c_bench() {
    C_TYPE* buffer        = (C_TYPE*)malloc(sizeof(C_TYPE)*BENCH_N);
    C_TYPE* buffer_origin = (C_TYPE*)malloc(sizeof(C_TYPE)*BENCH_N);

    srand(1);
    for (int i = 0; i < BENCH_N; ++i) {
        buffer[i] = (2 * (rand()/(C_TYPE)RAND_MAX)) - 1;
    }
    printf("\n");

    memcpy(buffer_origin, buffer, BENCH_N*sizeof(C_TYPE));

    c_compress_state_t zfp = c_compress_init(buffer, BENCH_N, ZFP_BENCH_RATE, 0, 0);

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
}
