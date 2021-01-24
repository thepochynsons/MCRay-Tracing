#include "rfimage.h"

using rf_image_ = rf_image<512, 100, 322>;

__global__ __device__
void rf_image_::cuda_convolve(float* axial_kernel, float* lateral_kernel){

}

void rf_image_::cuda_convolve_wrapper(float* axial_kernel, float* lateral_kernel, int axial_size, int lateral_size)
{
    float* dev_axial_kernel;
    float* dev_lateral_kernel;
    //cv::Mat dev_intensities;
    cudaMalloc( (void**)&dev_axial_kernel, sizeof(float)*axial_size );
    cudaMalloc( (void**)&dev_lateral_kernel, sizeof(float)*lateral_size);

    cudaMemcpy( dev_axial_kernel, axial_kernel, sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy( dev_lateral_kernel, lateral_kernel, sizeof(float), cudaMemcpyHostToDevice);
    //COPIAR DEV_INTENSITIES

    cuda_convolve<<<1,1>>>(axial_kernel, lateral_kernel);

    cudaFree(dev_axial_kernel);
    cudaFree(dev_lateral_kernel);
}
