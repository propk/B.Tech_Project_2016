#include<stdio.h>
#include <stdlib.h>
__device__ void T2x2H(int *iCoeff, int valRound)
{
    int valT1, valT2;
    iCoeff[0] += iCoeff[3];
    iCoeff[1] -= iCoeff[2];
    valT1 = ((iCoeff[0] - iCoeff[1] + valRound) >> 1);
    valT2 = iCoeff[2];
    iCoeff[2] = valT1 - iCoeff[3];
    iCoeff[3] = valT1 - valT2;
    iCoeff[0] -= iCoeff[3];
    iCoeff[1] += iCoeff[2];
}

__device__ void InvTOdd(int *iCoeff, int dummy)
{
    iCoeff[1] += iCoeff[3];
    iCoeff[0] -= iCoeff[2];
    iCoeff[3] -= (iCoeff[1] >> 1);
    iCoeff[2] += ((iCoeff[0] + 1) >> 1);
    iCoeff[0] -= ((3* iCoeff[1] + 4) >> 3);
    iCoeff[1] += ((3* iCoeff[0] + 4) >> 3);
    iCoeff[2] -= ((3* iCoeff[3] + 4) >> 3);
    iCoeff[3] += (3* iCoeff[2] + 4) >> 3;
    iCoeff[2] -= ((iCoeff[1] + 1) >> 1);
    iCoeff[3] = ((iCoeff[0] + 1) >> 1) - iCoeff[3];
    iCoeff[1] += iCoeff[2];
    iCoeff[0] -= iCoeff[3];
}

__device__ void InvTOddOdd(int *iCoeff, int dummy)
{
    int valT1, valT2;
    iCoeff[3] += iCoeff[0];
    iCoeff[2] -= iCoeff[1];
    valT1 = iCoeff[3] >> 1;
    valT2 = iCoeff[2] >> 1;
    iCoeff[0] -= valT1;
    iCoeff[1] += valT2;
    iCoeff[0] -= ((iCoeff[1] * 3 + 3) >> 3);
    iCoeff[1] += ((iCoeff[0] * 3 + 3) >> 2);
    iCoeff[0] -= ((iCoeff[1] * 3 + 4) >> 3);
    iCoeff[1] -= valT2;
    iCoeff[0] += valT1;
    iCoeff[2] += iCoeff[1];
    iCoeff[3] -= iCoeff[0];
    iCoeff[1] = -iCoeff[1];
    iCoeff[2] = -iCoeff[2];
}

__device__ void InvPermute(int *arrayInput)
{
    int i;
    int arrayTemp[16];
    int InvPermArr[16] = {
        0, 8, 4, 13, 2, 15, 3, 14,
        1, 12, 5, 9, 7, 11, 6, 10
    };
    for (i = 0; i <= 15; i++)
        arrayTemp[InvPermArr[i]] = arrayInput[i];
    for (i = 0; i <= 15; i++)
        arrayInput[i] = arrayTemp[i];
}

__device__ void InvPermute2pt(int *arrayInput)
{
    int arrayTemp[2], i;
    arrayTemp[0] = arrayInput[1];
    arrayTemp[1] = arrayInput[0];
    for (i = 0; i <= 1; i++)
        arrayInput[i] = arrayTemp[i];
}

__device__ void T2pt(int *iCoeff)
{
    iCoeff[0] -= (iCoeff[1]+1) >> 1;
    iCoeff[1] += iCoeff[0];
}

__device__ void (*pointerFunct[8]) (int *arg1, int arg2) = {
    T2x2H, InvTOdd, InvTOdd, InvTOddOdd,
    T2x2H, T2x2H, T2x2H, T2x2H
};

__device__ void ICT4x4(int *iCoeff)
{
    InvPermute(iCoeff);
    int arrayLocal[4];
    int arrayTemp[8][4] = {
        { 0, 1, 4, 5},
        { 2, 3, 6, 7},
        { 8, 12, 9, 13},
        { 10, 11, 14, 15},
        { 0, 3, 12, 15},
        { 5, 6, 9, 10},
        { 1, 2, 13, 14},
        { 4, 7, 8, 11}
    };


    int arg2Array[8] = { 1, 0, 0, 0, 0, 0, 0, 0};
    int i, j;

    for(i = 0; i < 8; i++)
    {
        for(j = 0; j < 4; j++)
            arrayLocal[j] = iCoeff[arrayTemp[i][j]];

        (*pointerFunct[i]) (arrayLocal, arg2Array[i]);

        for(j = 0; j < 4; j++)
            iCoeff[arrayTemp[i][j]] = arrayLocal[j];
    }
}

__device__ void InvTOddOddPOST(int* iCoeff, int dummy)
{
	int valT1, valT2;
	iCoeff[3] += iCoeff[0];
	iCoeff[2] -= iCoeff[1];
	valT1 = iCoeff[3] >> 1;
	valT2 = iCoeff[2] >> 1;
	iCoeff[0] -= valT1;
	iCoeff[1] += valT2;
	iCoeff[0] -= (iCoeff[1] * 3 + 6) >> 3;
	iCoeff[1] += (iCoeff[0] * 3 + 2) >> 2;
	iCoeff[0] -= (iCoeff[1] * 3 + 4) >> 3;
	iCoeff[1] -= valT2;
	iCoeff[0] += valT1;
	iCoeff[2] += iCoeff[1];
	iCoeff[3] -= iCoeff[0];
}

__device__ void T2x2HPOST(int* iCoeff, int dummy)
{
	int valT1;
	iCoeff[1] -= iCoeff[2];
	iCoeff[0] += (iCoeff[3] * 3 + 4) >> 3;
	iCoeff[3] -= (iCoeff[1] >> 1);
	iCoeff[2] = ((iCoeff[0] - iCoeff[1]) >> 1) - iCoeff[2];
	valT1 = iCoeff[2];
	iCoeff[2] = iCoeff[3];
	iCoeff[3] = valT1;
	iCoeff[0] -= iCoeff[3];
	iCoeff[1] += iCoeff[2];
}

__device__ void InvScale (int* iCoeff, int dummy)
{
	iCoeff[0] += iCoeff[1];
	iCoeff[1] = (iCoeff[0] >> 1) - iCoeff[1];
	iCoeff[0] += (iCoeff[1] * 3 + 0) >> 3;
	iCoeff[1] += (iCoeff[0] * 3 + 0) >> 4;
	iCoeff[1] += (iCoeff[0] >> 7);
	iCoeff[1] -= (iCoeff[0] >> 10);
}

__device__ void InvRotate(int* iCoeff, int dummy)
{
	iCoeff[0] -= ((iCoeff[1] + 1) >> 1);
	iCoeff[1] += ((iCoeff[0] + 1) >> 1);
}

__device__ void OverlapPostFilter2(int* iCoeff)
{
	iCoeff[1] += ((iCoeff[0] + 2) >> 2);
	iCoeff[0] += ((iCoeff[1] + 1) >> 1);
	iCoeff[0] += (iCoeff[1] >> 5);
	iCoeff[0] += (iCoeff[1] >> 9);
	iCoeff[0] += (iCoeff[1] >> 13);
	iCoeff[1] += ((iCoeff[0] + 2) >> 2);
}

__device__ void OverlapPostFilter2x2(int* iCoeff)
{
	iCoeff[0] += iCoeff[3];
	iCoeff[1] += iCoeff[2];
	iCoeff[3] -= ((iCoeff[0] + 1) >> 1);
	iCoeff[2] -= ((iCoeff[1] + 1) >> 1);
	iCoeff[1] += ((iCoeff[0] + 2) >> 2);
	iCoeff[0] += ((iCoeff[1] + 1) >> 1);
	iCoeff[0] += (iCoeff[1] >> 5);
	iCoeff[0] += (iCoeff[1] >> 9);
	iCoeff[0] += (iCoeff[1] >> 13);
	iCoeff[1] += ((iCoeff[0] + 2) >> 2);
	iCoeff[3] += ((iCoeff[0] + 1) >> 1);
	iCoeff[2] += ((iCoeff[1] + 1) >> 1);
	iCoeff[0] -= iCoeff[3];
	iCoeff[1] -= iCoeff[2];
}

__device__ void OverlapPostFilter4(int *iCoeff)
{
    int arrayLocal[2];
	iCoeff[0] += iCoeff[3];
	iCoeff[1] += iCoeff[2];
	iCoeff[3] -= ((iCoeff[0] + 1) >> 1);
	iCoeff[2] -= ((iCoeff[1] + 1) >> 1);
	arrayLocal[0] = iCoeff[0], arrayLocal[1] = iCoeff[3];
	InvScale(arrayLocal, 0);
	iCoeff[0] = arrayLocal[0], iCoeff[3] = arrayLocal[1];

	arrayLocal[0] = iCoeff[1], arrayLocal[1] = iCoeff[2];
	InvScale(arrayLocal, 0);
	iCoeff[1] = arrayLocal[0], iCoeff[2] = arrayLocal[1];
	iCoeff[0] += ((iCoeff[3] * 3+ 4) >> 3);
	iCoeff[1] += ((iCoeff[2] * 3 + 4) >> 3);
	iCoeff[3] -= ( iCoeff[0] >> 1);
	iCoeff[2] -= ( iCoeff[1] >> 1);
	iCoeff[0] += iCoeff[3];
	iCoeff[1] += iCoeff[2];
	iCoeff[3] = -iCoeff[3];
	iCoeff[2] = -iCoeff[2];
	arrayLocal[0] = iCoeff[2], arrayLocal[1] = iCoeff[3];
	InvRotate(arrayLocal, 0);
	iCoeff[2] = arrayLocal[0], iCoeff[3] = arrayLocal[1];
	iCoeff[3] += ((iCoeff[0] + 1) >> 1);
	iCoeff[2] += ((iCoeff[1] + 1) >> 1);
	iCoeff[0] -= iCoeff[3];
	iCoeff[1] -= iCoeff[2];
}

__device__ void (*pointerFunct[17]) (int *arg1, int arg2) = {
    T2x2H, T2x2H, T2x2H, T2x2H,
    InvRotate, InvRotate, InvRotate, InvRotate,
    InvTOddOddPOST, InvScale, InvScale, InvScale, InvScale,
    T2x2HPOST, T2x2HPOST, T2x2HPOST, T2x2HPOST
};

__device__ void OverlapPostFilter4x4(int *iCoeff)
{
	int arrayLocal[4];
    int arrayTemp[17][4] = {
        { 0, 3, 12, 15},
        { 1, 2, 13, 14},
        { 4, 7, 8, 11},
        { 5, 6, 9, 10},

        { 13, 12, -1, -1},
        { 9, 8, -1, -1},
        { 7, 3, -1, -1},
        { 6, 2, -1, -1},

        { 10, 11, 14, 15},

        { 0, 15, -1, -1},
        { 1, 14, -1, -1},
        { 4, 11, -1, -1},
        { 5, 10, -1, -1},


        { 0, 3, 12, 15},
        { 1, 2, 13, 14},
        { 4, 7, 8, 11},
        { 5, 6, 9, 10}
    };


    int i, j;

    for(i = 0; i < 17; i++)
    {
        for(j = 0; j < 4; j++)
            if(arrayTemp[i][j] >= 0)
                arrayLocal[j] = iCoeff[arrayTemp[i][j]];

        (*pointerFunct[i]) (arrayLocal, 0);

        for(j = 0; j < 4; j++)
            if(arrayTemp[i][j] >= 0)
                iCoeff[arrayTemp[i][j]] = arrayLocal[j];
    }
}

__global__ void DecFirstStagePostFiltering(int *image, int numRows, int numCols)
{
    int i, j;
    int arrayLocal[16];
    int block_i = threadIdx.x, block_j = threadIdx.y;
    int macro_i = blockIdx.x, macro_j = blockIdx.y;

    for (i = 0; i < 4; ++i)
    {
        for (j = 0; j < 4; ++j)
        {
            arrayLocal[i*4 + j] = image[ (macro_i*16 + block_i*4 + i) * numCols + macro_j*16 + block_j*4 + j];
        }
    }

    ICT4x4(arrayLocal);

    for (i = 0; i < 4; ++i)
    {
        for (j = 0; j < 4; ++j)
        {
            image[(macro_i*16 + block_i*4 + i) * numCols + macro_j*16 + block_j*4 + j] = arrayLocal[i*4 + j];
        }
    }
}

__global__ void DecSecondStagePostFiltering(int * image, int numRows, int numCols)
{
    int i, j;
    int arrayLocal[16];
    int macro_i = blockIdx.x, macro_j = blockIdx.y;

    for (i = 0; i < 4; ++i)
    {
        for (j = 0; j < 4; ++j)
        {
            arrayLocal[i*4 + j] = image[ (macro_i*16 + i*4) * numCols + macro_j*16 + j*4];
        }
    }
    ICT4x4(arrayLocal);
    for (i = 0; i < 4; ++i)
    {
        for (j = 0; j < 4; ++j)
        {
            image[(macro_i*16 + i*4) * numCols + macro_j*16 + j*4] = arrayLocal[i*4 + j];
        }
    }
}

__global__ void DecFirstStageOverlapFilter(int* image, int numRows, int numCols)
{
    int arrayLocal_16[16], arrayLocal_4[4];
    int block_i = blockIdx.x, block_j = blockIdx.y, i, j;

    //numRows /= 4;
    //numCols /= 4;

    // 4x4 blocks
    for( i = 0; i < 4; i++)
        for( j = 0; j < 4; j++)
            arrayLocal_16[i*4+j] = image[((block_i*4 + i + 2)*4)*numCols + (block_j*4 + j + 2)*4];

    OverlapPostFilter4x4(arrayLocal_16);
    for( i = 0; i < 4; i++)
        for( j = 0; j < 4; j++)
            image[((block_i*4 + i + 2)*4)*numCols + (block_j*4 + j + 2)*4] = arrayLocal_16[i*4+j];
    //4x4 block end

    if(block_j == 0)
    {
        //left edge
        for(i = 0; i < 2; i++)
        {
            for(j = 0; j < 4; j++)
                arrayLocal_4[j] = image[((block_i*4+j+2)*4)*numCols + i*4];
            OverlapPostFilter4(arrayLocal_4);
            for(j = 0; j < 4; j++)
                image[((block_i*4+j+2)*4)*numCols + i*4] = arrayLocal_4[j];
        }
        // right edge
        for(i = numCols/4-2; i < numCols/4; i++)
        {
            for(j = 0; j < 4; j++)
                arrayLocal_4[j] = image[((block_i*4+j+2)*4)*numCols + i*4];
            OverlapPostFilter4(arrayLocal_4);
            for(j = 0; j < 4; j++)
                image[((block_i*4+j+2)*4)*numCols + i*4] = arrayLocal_4[j];
        }
    }


    if(block_i == 0)
    {
        //top edge
        for(i = 0; i < 2; i++)
        {
            for(j = 0; j < 4; j++)
                arrayLocal_4[j] = image[(i*4)*numCols + (block_j*4+j+2)*4];
            OverlapPostFilter4(arrayLocal_4);
            for(j = 0; j < 4; j++)
                image[(i*4)*numCols + (block_j*4+j+2)*4] = arrayLocal_4[j];
        }
        //bottom edge
        for(i = numRows/4-2; i < numRows/4; i++)
        {
            for(j = 0; j < 4; j++)
                arrayLocal_4[j] = image[(i*4)*numCols + (block_j*4+2+j)*4];
            OverlapPostFilter4(arrayLocal_4);
            for(j = 0; j < 4; j++)
                image[(i*4)*numCols + (block_j*4+j+2)*4] = arrayLocal_4[j];
        }
    }

    if(block_j == 0 && block_i == 0)
    {
        // top left
        arrayLocal_4[0] = image[(0)*numCols + 0], arrayLocal_4[1] = image[(0)*numCols + 4];
        arrayLocal_4[2] = image[(4)*numCols + 0], arrayLocal_4[3] = image[(4)*numCols + 4];
        OverlapPostFilter4(arrayLocal_4);
        image[(0)*numCols + 0] = arrayLocal_4[0], image[(0)*numCols + 4] = arrayLocal_4[1];
        image[(4)*numCols + 0] = arrayLocal_4[2], image[(4)*numCols + 4] = arrayLocal_4[3];

        // top right
        arrayLocal_4[0] = image[(0)*numCols + (numCols/4-2)*4], arrayLocal_4[1] = image[(0)*numCols + (numCols/4-1)*4];
        arrayLocal_4[2] = image[(4)*numCols + (numCols/4-2)*4], arrayLocal_4[3] = image[(4)*numCols + (numCols/4-1)*4];
        OverlapPostFilter4(arrayLocal_4);
        image[(0)*numCols + (numCols/4-2)*4] = arrayLocal_4[0], image[(0)*numCols + (numCols/4-1)*4] = arrayLocal_4[1];
        image[(4)*numCols + (numCols/4-2)*4] = arrayLocal_4[2], image[(4)*numCols + (numCols/4-1)*4] = arrayLocal_4[3];

        // bottom left
        arrayLocal_4[0] = image[((numRows/4-2)*4)*numCols + 0], arrayLocal_4[1] = image[((numRows/4-2)*4)*numCols + 4];
        arrayLocal_4[2] = image[((numRows/4-1)*4)*numCols + 0], arrayLocal_4[3] = image[((numRows/4-1)*4)*numCols + 4];
        OverlapPostFilter4(arrayLocal_4);
        image[((numRows/4-2)*4)*numCols + 0] = arrayLocal_4[0], image[((numRows/4-2)*4)*numCols + 4] = arrayLocal_4[1];
        image[((numRows/4-1)*4)*numCols + 0] = arrayLocal_4[2], image[((numRows/4-1)*4)*numCols + 4] = arrayLocal_4[3];

        // bottom right
        arrayLocal_4[0] = image[((numRows/4-2)*4)*numCols + (numCols/4-2)*4], arrayLocal_4[1] = image[((numRows/4-2)*4)*numCols + (numCols/4-1)*4];
        arrayLocal_4[2] = image[((numRows/4-1)*4)*numCols + (numCols/4-2)*4], arrayLocal_4[3] = image[((numRows/4-1)*4)*numCols + (numCols/4-1)*4];
        OverlapPostFilter4(arrayLocal_4);
        image[((numRows/4-2)*4)*numCols + (numCols/4-2)*4] = arrayLocal_4[0], image[((numRows/4-2)*4)*numCols + (numCols/4-1)*4] = arrayLocal_4[1];
        image[((numRows/4-1)*4)*numCols + (numCols/4-2)*4] = arrayLocal_4[2], image[((numRows/4-1)*4)*numCols + (numCols/4-1)*4] = arrayLocal_4[3];
    }
}

__global__ void DecSecondStageOverlapFilter(int *image, int numRows, int numCols)
{
    int arrayLocal_16[16], arrayLocal_4[4];
    int block_i = blockIdx.x, block_j = blockIdx.y, i, j;

    // 4x4 blocks
    for( i = 0; i < 4; i++)
    {
        for( j = 0; j < 4; j++)
            arrayLocal_16[i*4+j] = image[ (block_i*4 + 2 + i)*numCols + block_j*4 + 2 + j];
    }
    OverlapPostFilter4x4(arrayLocal_16);
    for( i = 0; i < 4; i++)
    {
        for( j = 0; j < 4; j++)
            image[(block_i*4 + 2 + i)*numCols + block_j*4 + 2 + j] = arrayLocal_16[i*4+j];
    }
    //4x4 block end

    if(block_j == 0)
    {
        //left edge
        for(i = 0; i < 2; i++)
        {
            for(j = 0; j < 4; j++)
                arrayLocal_4[j] = image[(block_i*4 + 2+j)*numCols + i];
            OverlapPostFilter4(arrayLocal_4);
            for(j = 0; j < 4; j++)
                image[(block_i*4 + 2+j)*numCols + i] = arrayLocal_4[j];
        }
        // right edge
        for(i = numCols-2; i < numCols; i++)
        {
            for(j = 0; j < 4; j++)
                arrayLocal_4[j] = image[(block_i*4 + 2 + j)*numCols + i];
            OverlapPostFilter4(arrayLocal_4);
            for(j = 0; j < 4; j++)
                image[(block_i*4 + 2+j)*numCols + i] = arrayLocal_4[j];
        }
    }


    if(block_i == 0)
    {
        //top edge
        for(i = 0; i < 2; i++)
        {
            for(j = 0; j < 4; j++)
                arrayLocal_4[j] = image[(i)*numCols + block_j*4 + 2+j];
            OverlapPostFilter4(arrayLocal_4);
            for(j = 0; j < 4; j++)
                image[(i)*numCols + block_j*4 + 2+j] = arrayLocal_4[j];
        }
        //bottom edge
        for(i = numRows-2; i < numRows; i++)
        {
            for(j = 0; j < 4; j++)
                arrayLocal_4[j] = image[(i)*numCols + block_j*4 + 2 + j];
            OverlapPostFilter4(arrayLocal_4);
            for(j = 0; j < 4; j++)
                image[(i)*numCols + block_j*4 + 2 + j] = arrayLocal_4[j];
        }
    }

    if(block_i == 0 && block_j == 0)
    {
        // top left
        arrayLocal_4[0] = image[(0)*numCols + 0], arrayLocal_4[1] = image[(0)*numCols + 1];
        arrayLocal_4[2] = image[(1)*numCols + 0], arrayLocal_4[3] = image[(1)*numCols + 1];
        OverlapPostFilter4(arrayLocal_4);
        image[(0)*numCols + 0] = arrayLocal_4[0], image[(0)*numCols + 1] = arrayLocal_4[1];
        image[(1)*numCols + 0] = arrayLocal_4[2], image[(1)*numCols + 1] = arrayLocal_4[3];

        // top right
        arrayLocal_4[0] = image[(0)*numCols + numCols-2], arrayLocal_4[1] = image[(0)*numCols + numCols-1];
        arrayLocal_4[2] = image[(1)*numCols + numCols-2], arrayLocal_4[3] = image[(1)*numCols + numCols-1];
        OverlapPostFilter4(arrayLocal_4);
        image[(0)*numCols + numCols-2] = arrayLocal_4[0], image[(0)*numCols + numCols-1] = arrayLocal_4[1];
        image[(1)*numCols + numCols-2] = arrayLocal_4[2], image[(1)*numCols + numCols-1] = arrayLocal_4[3];

        // bottom left
        arrayLocal_4[0] = image[(numRows-2)*numCols + 0], arrayLocal_4[1] = image[(numRows-2)*numCols + 1];
        arrayLocal_4[2] = image[(numRows-1)*numCols + 0], arrayLocal_4[3] = image[(numRows-1)*numCols + 1];
        OverlapPostFilter4(arrayLocal_4);
        image[(numRows-2)*numCols + 0] = arrayLocal_4[0], image[(numRows-2)*numCols + 1] = arrayLocal_4[1];
        image[(numRows-1)*numCols + 0] = arrayLocal_4[2], image[(numRows-1)*numCols + 1] = arrayLocal_4[3];

        // bottom right
        arrayLocal_4[0] = image[(numRows-2)*numCols + numCols-2], arrayLocal_4[1] = image[(numRows-2)*numCols + numCols-1];
        arrayLocal_4[2] = image[(numRows-1)*numCols + numCols-2], arrayLocal_4[3] = image[(numRows-1)*numCols + numCols-1];
        OverlapPostFilter4(arrayLocal_4);
        image[(numRows-2)*numCols + numCols-2] = arrayLocal_4[0], image[(numRows-2)*numCols + numCols-1] = arrayLocal_4[1];
        image[(numRows-1)*numCols + numCols-2] = arrayLocal_4[2], image[(numRows-1)*numCols + numCols-1] = arrayLocal_4[3];
    }
}


int main()
{
    // read image in host
    int imageWidth = 112, imageHeight = 128;
    //scanf("%d %d", &imageHeight, &imageWidth);
    int image[128][112]; // = (int**) malloc(imageHeight * sizeof(int*) );
    int i, j;
    for(i = 0; i < imageHeight; i++){
        //image[i] = (int*) malloc(imageWidth * sizeof(int) );

        for(j = 0; j < imageWidth; j++)
            scanf( "%d", &image[i][j]);
    }

    // allocate & copy image memory in device
    int *imageDevice;
    int size = imageWidth * imageHeight * sizeof(int);
    cudaMalloc((void**) &imageDevice, size );
    cudaMemcpy(imageDevice, image, size, cudaMemcpyHostToDevice);

    /* kernel invocation start */
    dim3 DimGrid(imageHeight/16, imageWidth/16);
    dim3 DimBlock(4, 4);
    dim3 DimGrid2(imageHeight/4-1, imageWidth/4-1);
    dim3 DimGrid3(imageHeight/16-1, imageWidth/16-1);
    // second stage Post-filtering
    DecSecondStagePostFiltering<<< DimGrid, 1>>>(imageDevice, imageHeight, imageWidth);
    // first stage frequency transform
    DecFirstStageOverlapFilter<<< DimGrid3, 1>>>(imageDevice, imageHeight, imageWidth);
    // first stage Post-filtering
    DecFirstStagePostFiltering<<< DimGrid, DimBlock>>>(imageDevice, imageHeight, imageWidth);
    // second stage frequency transform
    DecSecondStageOverlapFilter<<< DimGrid2, 1>>>(imageDevice, imageHeight, imageWidth);

    /* kernel function invocation end */

    // copy from device to host
    cudaMemcpy(image, imageDevice, size, cudaMemcpyDeviceToHost);

    //free device memory
    cudaFree(imageDevice);

    //store processed image in file
    for( i = 0; i < imageHeight; i++)
    {
        for( j = 0; j < imageWidth; j++)
            fprintf(out, "%d ", image[i][j] );
        fprintf(out, "\n");
    }

    return 0;
}
