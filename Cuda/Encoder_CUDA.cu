#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
//__global__ int *arrayLocal_16, *arrayLocal_4;

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

__device__ void TOdd(int *iCoeff, int dummy)
{
    iCoeff[1] -= iCoeff[2];
    iCoeff[0] += iCoeff[3];
    iCoeff[2] += ((iCoeff[1] + 1) >> 1);
    iCoeff[3] = ((iCoeff[0] + 1) >> 1) - iCoeff[3];
    iCoeff[1] -= ((3* iCoeff[0] + 4) >> 3);
    iCoeff[0] += ((3* iCoeff[1] + 4) >> 3);
    iCoeff[3] -= ((3* iCoeff[2] + 4) >> 3);
    iCoeff[2] += ((3* iCoeff[3] + 4) >> 3);
    iCoeff[3] += (iCoeff[1] >> 1);
    iCoeff[2] -= ((iCoeff[0] + 1) >> 1);
    iCoeff[1] -= iCoeff[3];
    iCoeff[0] += iCoeff[2];
}

__device__ void TOddOdd(int *iCoeff, int dummy)
{
    int valT1, valT2;
    iCoeff[1] = -iCoeff[1];
    iCoeff[2] = -iCoeff[2];
    iCoeff[3] += iCoeff[0];
    iCoeff[2] -= iCoeff[1];
    iCoeff[0] -= (valT1 = iCoeff[3] >> 1);
    iCoeff[1] += (valT2 = iCoeff[2] >> 1);
    iCoeff[0] += ((iCoeff[1] * 3 + 4) >> 3);
    iCoeff[1] -= ((iCoeff[0] * 3 + 3) >> 2);
    iCoeff[0] += ((iCoeff[1] * 3 + 3) >> 3);
    iCoeff[1] -= valT2;
    iCoeff[0] += valT1;
    iCoeff[2] += iCoeff[1];
    iCoeff[3] -= iCoeff[0];
}

__device__ void FwdPermute(int *arrayInput)
{
    int i;
    int arrayTemp[16];
    int FwdPermArr[16] = {
        0, 8, 4, 6, 2, 10, 14, 12,
        1, 11, 15, 13, 9, 3, 7, 5
    };
    for (i = 0; i <= 15; i++)
        arrayTemp[FwdPermArr[i]] = arrayInput[i];
    for (i = 0; i <= 15; i++)
        arrayInput[i] = arrayTemp[i];
}

__device__ void (*pointerFunct_FCT[8]) (int *arg1, int arg2) = {
    T2x2H, T2x2H, T2x2H, T2x2H,
    T2x2H, TOdd, TOdd, TOddOdd
};
__device__ void FCT4x4(int *iCoeff)
{
    int arrayLocal[4];
    int arrayTemp[8][4] = {
        { 0, 3, 12, 15},
        { 5, 6, 9, 10},
        { 1, 2, 13, 14},
        { 4, 7, 8, 11},
        { 0, 1, 4, 5},
        { 2, 3, 6, 7},
        { 8, 12, 9, 13},
        { 10, 11, 14, 15}
    };


    int arg2Array[8] = { 0, 0, 0, 0, 1, 0, 0, 0};
    int i, j;

    for(i = 0; i < 8; i++)
    {
        for(j = 0; j < 4; j++)
            arrayLocal[j] = iCoeff[arrayTemp[i][j]];

        (*pointerFunct_FCT[i]) (arrayLocal, arg2Array[i]);

        for(j = 0; j < 4; j++)
            iCoeff[arrayTemp[i][j]] = arrayLocal[j];
    }

    FwdPermute(iCoeff);
}

__device__ void T2x2hEnc(int *iCoeff, int dummy)
{
    int valT1, valT2;
    iCoeff[0] += iCoeff[3];
    iCoeff[1] -= iCoeff[2];
    valT1 = iCoeff[3];
    valT2 = iCoeff[2];
    iCoeff[2] = ((iCoeff[0] - iCoeff[1]) >> 1) - valT1;
    iCoeff[3] = valT2 + (iCoeff[1] >> 1);
    iCoeff[1] += iCoeff[2];
    iCoeff[0] -= (iCoeff[3] * 3 + 4) >> 3;
}

__device__ void FwdRotate(int *iCoeff, int dummy)
{
    iCoeff[1] -= (iCoeff[0] + 1) >> 1;
    iCoeff[0] += (iCoeff[1] + 1) >> 1;
}

__device__ void FwdScale(int *iCoeff, int dummy)
{
    iCoeff[1] -= (iCoeff[0] * 3 + 0) >> 4;
    iCoeff[1] -= (iCoeff[0] >> 7);
    iCoeff[1] += (iCoeff[0] >> 10);
    iCoeff[0] -= (iCoeff[1] * 3 + 0) >> 3;
    iCoeff[1] = (iCoeff[0] >> 1) - iCoeff[1];
    iCoeff[0] -= iCoeff[1];
}

__device__ void FwdTOddOdd(int *iCoeff, int dummy)
{
    int valT1, valT2;
    iCoeff[3] += iCoeff[0];
    iCoeff[2] -= iCoeff[1];
    valT1 = iCoeff[3] >> 1;
    valT2 = iCoeff[2] >> 1;
    iCoeff[0] -= valT1;
    iCoeff[1] += valT2;
    iCoeff[0] += (iCoeff[1] * 3 + 4) >>3;
    iCoeff[1] -= (iCoeff[0] * 3 + 2) >> 2;
    iCoeff[0] += (iCoeff[1] * 3 + 6) >> 3;
    iCoeff[1] -= valT2;
    iCoeff[0] += valT1;
    iCoeff[2] += iCoeff[1];
    iCoeff[3] -= iCoeff[0];
}

__device__ void (*pointerFunct[17]) (int *arg1, int arg2) = {
    T2x2hEnc, T2x2hEnc, T2x2hEnc, T2x2hEnc,
    FwdScale, FwdScale, FwdScale, FwdScale,
    FwdRotate, FwdRotate, FwdRotate, FwdRotate,
    FwdTOddOdd, T2x2H, T2x2H, T2x2H, T2x2H
};

__device__ void OverlapPreFilter4x4(int *iCoeff)
{
    int arrayLocal[4];
    int arrayTemp[17][4] = {
        { 0, 3, 12, 15},
        { 1, 2, 13, 14},
        { 4, 7, 8, 11},
        { 5, 6, 9, 10},

        { 0, 15, -1, -1},
        { 1, 14, -1, -1},
        { 4, 11, -1, -1},
        { 5, 10, -1, -1},

        { 13, 12, -1, -1},
        { 9, 8, -1, -1},
        { 7, 3, -1, -1},
        { 6, 2, -1, -1},

        { 10, 11, 14, 15},
        { 0, 12, 3, 15},
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

__device__ void OverlapPreFilter4(int *iCoeff)
{
    int arrayLocal[2];
    iCoeff[0] += iCoeff[3];
    iCoeff[1] += iCoeff[2];
    iCoeff[3] -= ((iCoeff[0] + 1) >> 1);
    iCoeff[2] -= ((iCoeff[1] + 1) >> 1);
    arrayLocal[0] = iCoeff[2], arrayLocal[1] = iCoeff[3];
    FwdRotate(arrayLocal, 0);
    iCoeff[2] = arrayLocal[0], iCoeff[3] = arrayLocal[1];
    iCoeff[3] = -iCoeff[3];
    iCoeff[2] = -iCoeff[2];
    iCoeff[0] -= iCoeff[3];
    iCoeff[1] -= iCoeff[2];
    iCoeff[3] += ( iCoeff[0] >> 1);
    iCoeff[2] += ( iCoeff[1] >> 1);
    iCoeff[0] -= ((iCoeff[3] * 3 + 4) >> 3);
    iCoeff[1] -= ((iCoeff[2] * 3 + 4) >> 3);
    arrayLocal[0] = iCoeff[0], arrayLocal[1] = iCoeff[3];
    FwdScale(arrayLocal, 0);
    iCoeff[0] = arrayLocal[0], iCoeff[3] = arrayLocal[1];

    arrayLocal[0] = iCoeff[1], arrayLocal[1] = iCoeff[2];
    FwdScale(arrayLocal, 0);
    iCoeff[1] = arrayLocal[0], iCoeff[2] = arrayLocal[1];
    iCoeff[3] += ((iCoeff[0] + 1) >> 1);
    iCoeff[2] += ((iCoeff[1] + 1) >> 1);
    iCoeff[0] -= iCoeff[3];
    iCoeff[1] -= iCoeff[2];
}

__device__ void OverlapPreFilter2x2(int *iCoeff)
{
    iCoeff[0] += iCoeff[3];
    iCoeff[1] += iCoeff[2];
    iCoeff[3] -= ((iCoeff[0] + 1) >> 1);
    iCoeff[2] -= ((iCoeff[1] + 1) >> 1);
    iCoeff[1] -= ((iCoeff[0] + 2) >> 2);
    iCoeff[0] -= (iCoeff[1] >> 5);
    iCoeff[0] -= (iCoeff[1] >> 9);
    iCoeff[0] -= (iCoeff[1] >> 13);
    iCoeff[0] -= ((iCoeff[1] + 1) >> 1);
    iCoeff[1] -= ((iCoeff[0] + 2) >> 2);
    iCoeff[3] += ((iCoeff[0] + 1) >> 1);
    iCoeff[2] += ((iCoeff[1] + 1) >> 1);
    iCoeff[0] -= iCoeff[3];
    iCoeff[1] -= iCoeff[2];
}

__device__ void OverlapPreFilter2(int *iCoeff)
{
    iCoeff[1] -= ((iCoeff[0] + 2) >> 2);
    iCoeff[0] -= (iCoeff[1] >> 13);
    iCoeff[0] -= (iCoeff[1] >> 9);
    iCoeff[0] -= (iCoeff[1] >> 5);
    iCoeff[0] -= ((iCoeff[1] + 1) >> 1);
    iCoeff[1] -= ((iCoeff[0] + 2) >> 2);
}

__global__ void EncFirstStagePreFiltering(int *image, int numRows, int numCols)
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

    FCT4x4(arrayLocal);

    for (i = 0; i < 4; ++i)
    {
        for (j = 0; j < 4; ++j)
        {
            image[(macro_i*16 + block_i*4 + i) * numCols + macro_j*16 + block_j*4 + j] = arrayLocal[i*4 + j];
        }
    }
}

__global__ void EncSecondStagePreFiltering(int* image, int numRows, int numCols)
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
    FCT4x4(arrayLocal);
    for (i = 0; i < 4; ++i)
    {
        for (j = 0; j < 4; ++j)
        {
            image[(macro_i*16 + i*4) * numCols + macro_j*16 + j*4] = arrayLocal[i*4 + j];
        }
    }
}

__global__ void EncFirstStageOverlapFilter(int* image, int numRows, int numCols)
{
    int arrayLocal_16[16], arrayLocal_4[4];
    int block_i = blockIdx.x, block_j = blockIdx.y, i, j;

    //numRows /= 4;
    //numCols /= 4;

    // 4x4 blocks
    for( i = 0; i < 4; i++)
        for( j = 0; j < 4; j++)
            arrayLocal_16[i*4+j] = image[((block_i*4 + i + 2)*4)*numCols + (block_j*4 + j + 2)*4];

    OverlapPreFilter4x4(arrayLocal_16);
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
            OverlapPreFilter4(arrayLocal_4);
            for(j = 0; j < 4; j++)
                image[((block_i*4+j+2)*4)*numCols + i*4] = arrayLocal_4[j];
        }
        // right edge
        for(i = numCols/4-2; i < numCols/4; i++)
        {
            for(j = 0; j < 4; j++)
                arrayLocal_4[j] = image[((block_i*4+j+2)*4)*numCols + i*4];
            OverlapPreFilter4(arrayLocal_4);
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
            OverlapPreFilter4(arrayLocal_4);
            for(j = 0; j < 4; j++)
                image[(i*4)*numCols + (block_j*4+j+2)*4] = arrayLocal_4[j];
        }
        //bottom edge
        for(i = numRows/4-2; i < numRows/4; i++)
        {
            for(j = 0; j < 4; j++)
                arrayLocal_4[j] = image[(i*4)*numCols + (block_j*4+2+j)*4];
            OverlapPreFilter4(arrayLocal_4);
            for(j = 0; j < 4; j++)
                image[(i*4)*numCols + (block_j*4+j+2)*4] = arrayLocal_4[j];
        }
    }

    if(block_j == 0 && block_i == 0)
    {
        // top left
        arrayLocal_4[0] = image[(0)*numCols + 0], arrayLocal_4[1] = image[(0)*numCols + 4];
        arrayLocal_4[2] = image[(4)*numCols + 0], arrayLocal_4[3] = image[(4)*numCols + 4];
        OverlapPreFilter4(arrayLocal_4);
        image[(0)*numCols + 0] = arrayLocal_4[0], image[(0)*numCols + 4] = arrayLocal_4[1];
        image[(4)*numCols + 0] = arrayLocal_4[2], image[(4)*numCols + 4] = arrayLocal_4[3];

        // top right
        arrayLocal_4[0] = image[(0)*numCols + (numCols/4-2)*4], arrayLocal_4[1] = image[(0)*numCols + (numCols/4-1)*4];
        arrayLocal_4[2] = image[(4)*numCols + (numCols/4-2)*4], arrayLocal_4[3] = image[(4)*numCols + (numCols/4-1)*4];
        OverlapPreFilter4(arrayLocal_4);
        image[(0)*numCols + (numCols/4-2)*4] = arrayLocal_4[0], image[(0)*numCols + (numCols/4-1)*4] = arrayLocal_4[1];
        image[(4)*numCols + (numCols/4-2)*4] = arrayLocal_4[2], image[(4)*numCols + (numCols/4-1)*4] = arrayLocal_4[3];

        // bottom left
        arrayLocal_4[0] = image[((numRows/4-2)*4)*numCols + 0], arrayLocal_4[1] = image[((numRows/4-2)*4)*numCols + 4];
        arrayLocal_4[2] = image[((numRows/4-1)*4)*numCols + 0], arrayLocal_4[3] = image[((numRows/4-1)*4)*numCols + 4];
        OverlapPreFilter4(arrayLocal_4);
        image[((numRows/4-2)*4)*numCols + 0] = arrayLocal_4[0], image[((numRows/4-2)*4)*numCols + 4] = arrayLocal_4[1];
        image[((numRows/4-1)*4)*numCols + 0] = arrayLocal_4[2], image[((numRows/4-1)*4)*numCols + 4] = arrayLocal_4[3];

        // bottom right
        arrayLocal_4[0] = image[((numRows/4-2)*4)*numCols + (numCols/4-2)*4], arrayLocal_4[1] = image[((numRows/4-2)*4)*numCols + (numCols/4-1)*4];
        arrayLocal_4[2] = image[((numRows/4-1)*4)*numCols + (numCols/4-2)*4], arrayLocal_4[3] = image[((numRows/4-1)*4)*numCols + (numCols/4-1)*4];
        OverlapPreFilter4(arrayLocal_4);
        image[((numRows/4-2)*4)*numCols + (numCols/4-2)*4] = arrayLocal_4[0], image[((numRows/4-2)*4)*numCols + (numCols/4-1)*4] = arrayLocal_4[1];
        image[((numRows/4-1)*4)*numCols + (numCols/4-2)*4] = arrayLocal_4[2], image[((numRows/4-1)*4)*numCols + (numCols/4-1)*4] = arrayLocal_4[3];
    }
}

__global__ void EncSecondStageOverlapFilter(int* image, int numRows, int numCols)
{
    int arrayLocal_16[16], arrayLocal_4[4];
    int block_i = blockIdx.x, block_j = blockIdx.y, i, j;

    // 4x4 blocks
    for( i = 0; i < 4; i++)
    {
        for( j = 0; j < 4; j++)
            arrayLocal_16[i*4+j] = image[ (block_i*4 + 2 + i)*numCols + block_j*4 + 2 + j];
    }
    OverlapPreFilter4x4(arrayLocal_16);
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
            OverlapPreFilter4(arrayLocal_4);
            for(j = 0; j < 4; j++)
                image[(block_i*4 + 2+j)*numCols + i] = arrayLocal_4[j];
        }
        // right edge
        for(i = numCols-2; i < numCols; i++)
        {
            for(j = 0; j < 4; j++)
                arrayLocal_4[j] = image[(block_i*4 + 2 + j)*numCols + i];
            OverlapPreFilter4(arrayLocal_4);
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
            OverlapPreFilter4(arrayLocal_4);
            for(j = 0; j < 4; j++)
                image[(i)*numCols + block_j*4 + 2+j] = arrayLocal_4[j];
        }
        //bottom edge
        for(i = numRows-2; i < numRows; i++)
        {
            for(j = 0; j < 4; j++)
                arrayLocal_4[j] = image[(i)*numCols + block_j*4 + 2 + j];
            OverlapPreFilter4(arrayLocal_4);
            for(j = 0; j < 4; j++)
                image[(i)*numCols + block_j*4 + 2 + j] = arrayLocal_4[j];
        }
    }

    if(block_i == 0 && block_j == 0)
    {
        // top left
        arrayLocal_4[0] = image[(0)*numCols + 0], arrayLocal_4[1] = image[(0)*numCols + 1];
        arrayLocal_4[2] = image[(1)*numCols + 0], arrayLocal_4[3] = image[(1)*numCols + 1];
        OverlapPreFilter4(arrayLocal_4);
        image[(0)*numCols + 0] = arrayLocal_4[0], image[(0)*numCols + 1] = arrayLocal_4[1];
        image[(1)*numCols + 0] = arrayLocal_4[2], image[(1)*numCols + 1] = arrayLocal_4[3];

        // top right
        arrayLocal_4[0] = image[(0)*numCols + numCols-2], arrayLocal_4[1] = image[(0)*numCols + numCols-1];
        arrayLocal_4[2] = image[(1)*numCols + numCols-2], arrayLocal_4[3] = image[(1)*numCols + numCols-1];
        OverlapPreFilter4(arrayLocal_4);
        image[(0)*numCols + numCols-2] = arrayLocal_4[0], image[(0)*numCols + numCols-1] = arrayLocal_4[1];
        image[(1)*numCols + numCols-2] = arrayLocal_4[2], image[(1)*numCols + numCols-1] = arrayLocal_4[3];

        // bottom left
        arrayLocal_4[0] = image[(numRows-2)*numCols + 0], arrayLocal_4[1] = image[(numRows-2)*numCols + 1];
        arrayLocal_4[2] = image[(numRows-1)*numCols + 0], arrayLocal_4[3] = image[(numRows-1)*numCols + 1];
        OverlapPreFilter4(arrayLocal_4);
        image[(numRows-2)*numCols + 0] = arrayLocal_4[0], image[(numRows-2)*numCols + 1] = arrayLocal_4[1];
        image[(numRows-1)*numCols + 0] = arrayLocal_4[2], image[(numRows-1)*numCols + 1] = arrayLocal_4[3];

        // bottom right
        arrayLocal_4[0] = image[(numRows-2)*numCols + numCols-2], arrayLocal_4[1] = image[(numRows-2)*numCols + numCols-1];
        arrayLocal_4[2] = image[(numRows-1)*numCols + numCols-2], arrayLocal_4[3] = image[(numRows-1)*numCols + numCols-1];
        OverlapPreFilter4(arrayLocal_4);
        image[(numRows-2)*numCols + numCols-2] = arrayLocal_4[0], image[(numRows-2)*numCols + numCols-1] = arrayLocal_4[1];
        image[(numRows-1)*numCols + numCols-2] = arrayLocal_4[2], image[(numRows-1)*numCols + numCols-1] = arrayLocal_4[3];
    }
}


int main()
{
    struct timeval tim[50];
    int t = 0;

    FILE *time_log = fopen("timing_log.txt", "w");

    gettimeofday(&tim[t], NULL); t++;
    FILE *ip = fopen("BImage", "r");
    FILE *op = fopen("encoded.txt", "w");

    // read image in host
    int imageWidth = 16000, imageHeight=8000;
    //scanf("%d %d", &imageHeight, &imageWidth);
    //int image[128][112]; // = 
    int *image = (int*) malloc(imageHeight * imageWidth * sizeof(int) );
    int i, j;

    for(i = 0; i < imageHeight; i++){
        //image[i] = (int*) malloc(imageWidth * sizeof(int) );

        for(j = 0; j < imageWidth; j++)
            fscanf(ip, "%d", &image[i * imageWidth + j]);
    }
    gettimeofday(&tim[t], NULL); t++;

    // allocate & copy image memory in device
    int *imageDevice;
    size_t size = imageWidth * imageHeight * sizeof(int);
    cudaMalloc((void**) &imageDevice, size );
    cudaMemcpy(imageDevice, image, size, cudaMemcpyHostToDevice);

    /* kernel invocation start*/
    dim3 DimGrid(imageHeight/16, imageWidth/16);
    dim3 DimBlock(4, 4);
    dim3 DimGrid2(imageHeight/4-1, imageWidth/4-1);
    dim3 DimGrid3(imageHeight/16-1, imageWidth/16-1);

    gettimeofday(&tim[t], NULL); t++;

    // second stage frequency transform
    EncSecondStageOverlapFilter<<< DimGrid2, 1>>>(imageDevice, imageHeight, imageWidth);
    gettimeofday(&tim[t], NULL); t++;

    // first stage pre-filtering
    EncFirstStagePreFiltering<<< DimGrid, DimBlock>>>(imageDevice, imageHeight, imageWidth);
    gettimeofday(&tim[t], NULL); t++;

    // first stage frequency transform
    EncFirstStageOverlapFilter<<< DimGrid3, 1>>>(imageDevice, imageHeight, imageWidth);
    gettimeofday(&tim[t], NULL); t++;

    // second stage pre-filtering
    EncSecondStagePreFiltering<<< DimGrid, 1>>>(imageDevice, imageHeight, imageWidth);
    gettimeofday(&tim[t], NULL); t++;

    /* kernel function invocation end*/
    cudaDeviceSynchronize();
    gettimeofday(&tim[t], NULL); t++;

    // copy from device to host
    cudaMemcpy(image, imageDevice, size, cudaMemcpyDeviceToHost);

    //free device memory
    cudaFree(imageDevice);
    gettimeofday(&tim[t], NULL); t++;

    //store processed image in file
    for( i = 0; i < imageHeight; i++)
    {
        for( j = 0; j < imageWidth; j++)
            fprintf(op, "%d ", image[i * imageWidth + j] );
        fprintf(op, "\n");
    }
    gettimeofday(&tim[t], NULL); t++;

    //time logs
    i = 0;

    fprintf(time_log, "Image Read from file: %lu seconds %lu microseconds\n", tim[i+1].tv_sec - tim[i].tv_sec, tim[i+1].tv_usec - tim[i].tv_usec);
    i++;
    fprintf(time_log, "Image Copy To Device: %lu seconds %lu microseconds\n", tim[i+1].tv_sec - tim[i].tv_sec, tim[i+1].tv_usec - tim[i].tv_usec);
    i++;
    fprintf(time_log, "Stage 1: %lu seconds %lu microseconds\n", tim[i+1].tv_sec - tim[i].tv_sec, tim[i+1].tv_usec - tim[i].tv_usec);
    i++;
    fprintf(time_log, "Stage 2: %lu seconds %lu microseconds\n", tim[i+1].tv_sec - tim[i].tv_sec, tim[i+1].tv_usec - tim[i].tv_usec);
    i++;
    fprintf(time_log, "Stage 3: %lu seconds %lu microseconds\n", tim[i+1].tv_sec - tim[i].tv_sec, tim[i+1].tv_usec - tim[i].tv_usec);
    i++;
    fprintf(time_log, "Stage 4: %lu seconds %lu microseconds\n", tim[i+1].tv_sec - tim[i].tv_sec, tim[i+1].tv_usec - tim[i].tv_usec);
    i++;
    fprintf(time_log, "Cuda Device Synchronise: %lu seconds %lu microseconds\n", tim[i+1].tv_sec - tim[i].tv_sec, tim[i+1].tv_usec - tim[i].tv_usec);
    i++;
    fprintf(time_log, "Image COpy back to Host: %lu seconds %lu microseconds\n", tim[i+1].tv_sec - tim[i].tv_sec, tim[i+1].tv_usec - tim[i].tv_usec);
    i++;
    fprintf(time_log, "Writing Encoded to File: %lu seconds %lu microseconds\n", tim[i+1].tv_sec - tim[i].tv_sec, tim[i+1].tv_usec - tim[i].tv_usec);

    fclose(ip);
    fclose(op);
    fclose(time_log);
    return 0;
}
