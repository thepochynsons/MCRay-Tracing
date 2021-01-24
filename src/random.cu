//#ifndef RANDOM_H
//#define RANDOM_H

//#endif // RANDOM_H

//#include <stdio.h>
//#include <stdlib.h>
//#include <cuda.h>
//#include <curand_kernel.h>
//#include "device_launch_parameters.h"

//#define CUDA_CALL(x) do { if((x) != cudaSuccess) { \
//    printf("Error at %s:%d\n",__FILE__,__LINE__); \
//    }} while(0)

//__global__ void setup_kernel(curandState *state)
//{
//    int id = threadIdx.x + blockIdx.x * 64;
//    /* Each thread gets same seed, a different sequence number,
//       no offset */
//    curand_init(1234, id, 0, &state[id]);
//}
//__global__ float generate_kernel(curandState *state, int *result)
//{
//    int id = threadIdx.x + blockIdx.x * 64;
//    int count = 0;
//    float x;
//    /* Copy state to local memory for efficiency */
//            curandState localState = state[id];
//    /* Generate pseudo-random unsigned ints */
//    for(int n = 0; n < 100000; n++) {
//        x = curand_uniform(&localState);
//        /* Check if low bit set */

//    }
//    /* Copy state back to global memory */
//    state[id] = localState;
//    /* Store results */
//    result[id] += count;
//}



//__host__ void generate_random_numbers(float randoms[], unsigned int size)
//{
//    curandState *devStates;
//    int *devResults, *hostResults;


//    /* Allocate space for results on host */

//    //hostResults = static_cast<int*>(calloc(size, sizeof(int)));

//    /* Allocate space for results on device */

//    CUDA_CALL(cudaMalloc((void **)&devResults, size * sizeof(float)));

//    /* Set results to 0 */

//    CUDA_CALL(cudaMemset(devResults, 0, size * sizeof(float)));

//    /* Allocate space for prng states on device */

//    CUDA_CALL(cudaMalloc((void **)&devStates, size *
//                         sizeof(curandState)));

//    /* Setup prng states */
//    setup_kernel<<<1024, 1024>>>(devStates);
//    /* Generate and use pseudo-random */
//    generate_kernel<<<1024, 1024>>>(devStates, devResults);

//    /* Copy device memory to host */
//    CUDA_CALL(cudaMemcpy(randoms, devResults, size *
//                         sizeof(float), cudaMemcpyDeviceToHost));
//    /* Show result */



//    /* Cleanup device*/
//    CUDA_CALL(cudaFree(devStates));
//    CUDA_CALL(cudaFree(devResults));
//}
