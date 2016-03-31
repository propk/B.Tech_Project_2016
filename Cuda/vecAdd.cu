// MP 1
#include <stdio.h>
#include <stdlib.h>

__global__ void vecAdd(float *in1, float *in2, float *out, int len) {
  //@@ Insert code to implement vector addition here
	int i = threadIdx.x +blockDim.x*blockIdx.x;
	if (i < len) out[i] = in1[i] + in2[i];
}

int main() {
  int inputLength = 4;
  float hostInput1[] = {1,2,3,4};
  float hostInput2[] = {5,6,7,8};
  float hostOutput[] = {0,0,0,0};
  float *deviceInput1;
  float *deviceInput2;
  float *deviceOutput;

  //@@ Allocate GPU memory here
	
	cudaMalloc((void**) & deviceInput1, inputLength * sizeof(float));
	cudaMalloc((void**) & deviceInput2, inputLength * sizeof(float));
	cudaMalloc((void**) & deviceOutput, inputLength * sizeof(float));
	
  //@@ Copy memory to the GPU here
	
	cudaMemcpy(deviceInput1, hostInput1, inputLength * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(deviceInput2, hostInput2, inputLength * sizeof(float), cudaMemcpyHostToDevice);

  //@@ Initialize the grid and block dimensions here
	
	dim3 DimGrid((inputLength-1)/256+1, 1, 1);
	dim3 DimBlock(256, 1, 1);

  //@@ Launch the GPU Kernel here
	
	vecAdd<<< DimGrid, DimBlock >>>(deviceInput1, deviceInput2, deviceOutput, inputLength);

  cudaDeviceSynchronize();
  //@@ Copy the GPU memory back to the CPU here
	
	cudaMemcpy(hostOutput, deviceOutput, inputLength * sizeof(float), cudaMemcpyDeviceToHost);

  //@@ Free the GPU memory here
	
	cudaFree(deviceInput1);
	cudaFree(deviceInput2);
	cudaFree(deviceOutput);
	
	printf("%lf, %lf, %lf, %lf", hostOutput[0], hostOutput[1], hostOutput[2], hostOutput[3]);

  return 0;
}
