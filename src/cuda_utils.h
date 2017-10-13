#ifndef om_cuda_tools_
#define om_cuda_tools_

#include <stdio.h>
#include <stdlib.h>

#ifndef __USE_POSIX199309
#define __USE_POSIX199309
#include <time.h>
#undef __USE_POSIX199309
#else
#include <time.h>
#endif

typedef struct timespec Timer;

static inline void init_timer(Timer* timer) {

  clock_gettime(CLOCK_MONOTONIC,timer);
}

static inline double get_timer(Timer* timer) {

  struct timespec end_timespec;
  clock_gettime(CLOCK_MONOTONIC,&end_timespec);
  return (end_timespec.tv_sec-timer->tv_sec)+
    (end_timespec.tv_nsec-timer->tv_nsec)*1e-9;
}

static inline double get_and_reset_timer(Timer* timer) {

  struct timespec end_timespec;
  clock_gettime(CLOCK_MONOTONIC,&end_timespec);
  double result=(end_timespec.tv_sec-timer->tv_sec)+
    (end_timespec.tv_nsec-timer->tv_nsec)*1e-9;
  *timer=end_timespec;
  return result;
}

static inline double get_timer_difference(Timer* timer_start,Timer* timer_end) {
  return (timer_end->tv_sec-timer_start->tv_sec)+
    (timer_end->tv_nsec-timer_start->tv_nsec)*1e-9;
}

// Goes through all available devices and prints their properties to the screen.
static inline void device_query() {

    int n_devices;

    cudaGetDeviceCount(&n_devices);

    for (int i = 0; i < n_devices; i++) {

      cudaDeviceProp prop;
      cudaGetDeviceProperties(&prop, i);

      printf("  Device Number:                          %d\n", i);
      printf("  Device Name:                            %s\n", prop.name);
      printf("  Memory Clock Rate (KHz):                %d\n",
             prop.memoryClockRate);
      printf("  Memory Bus Width (bits):                %d\n",
             prop.memoryBusWidth);
      printf("  Peak Memory Bandwidth (GB/s):           %f\n",
             2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e06);
      printf("  Maximum number of threads per block:    %d\n",
             prop.maxThreadsPerBlock);
      printf("  Max block dimensions:                   [ %d, %d, %d ]\n",
             prop.maxThreadsDim[0], prop.maxThreadsDim[1], prop.maxThreadsDim[2]);
      printf("  Max grid dimensions:                    [ %d, %d, %d ]\n\n",
             prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2]);
    }
}

#ifndef NDEBUG
#define cu_verify(x) do {                                                \
    cudaError_t result = x;                                              \
    if (result!=cudaSuccess) {                                           \
      fprintf(stderr,"%s:%i: error: cuda function call failed:\n"        \
              "  %s;\nmessage: %s\n",                                    \
              __FILE__,__LINE__,#x,cudaGetErrorString(result));          \
      exit(1);                                                           \
    }                                                                    \
  } while(0)
#define cu_verify_void(x) do {                                           \
    x;                                                                   \
    cudaError_t result = cudaGetLastError();                             \
    if (result!=cudaSuccess) {                                           \
      fprintf(stderr,"%s:%i: error: cuda function call failed:\n"        \
              "  %s;\nmessage: %s\n",                                    \
              __FILE__,__LINE__,#x,cudaGetErrorString(result));          \
      exit(1);                                                           \
    }                                                                    \
  } while(0)
#else
#define cu_verify(x) do {                                                \
    x;                                                                   \
    } while(0)
#define cu_verify_void(x) do {                                           \
    x;                                                                   \
    } while(0)
#endif

#endif
