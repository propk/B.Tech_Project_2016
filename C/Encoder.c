#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

void T2x2H(int *iCoeff, int valRound)
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

void TOdd(int *iCoeff, int dummy)
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

void TOddOdd(int *iCoeff, int dummy)
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

void FwdPermute(int *arrayInput)
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

void FCT4x4(int *iCoeff)
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

    void (*pointerFunct[8]) (int *arg1, int arg2) = {
        T2x2H, T2x2H, T2x2H, T2x2H,
        T2x2H, TOdd, TOdd, TOddOdd
    };
    int arg2Array[8] = { 0, 0, 0, 0, 1, 0, 0, 0};
    int i, j;

    for(i = 0; i < 8; i++)
    {
        for(j = 0; j < 4; j++)
            arrayLocal[j] = iCoeff[arrayTemp[i][j]];

        (*pointerFunct[i]) (arrayLocal, arg2Array[i]);

        for(j = 0; j < 4; j++)
            iCoeff[arrayTemp[i][j]] = arrayLocal[j];
    }

    FwdPermute(iCoeff);
}

void T2x2hEnc(int *iCoeff, int dummy)
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

void FwdRotate(int *iCoeff, int dummy)
{
    iCoeff[1] -= (iCoeff[0] + 1) >> 1;
    iCoeff[0] += (iCoeff[1] + 1) >> 1;
}

void FwdScale(int *iCoeff, int dummy)
{
    iCoeff[1] -= (iCoeff[0] * 3 + 0) >> 4;
    iCoeff[1] -= (iCoeff[0] >> 7);
    iCoeff[1] += (iCoeff[0] >> 10);
    iCoeff[0] -= (iCoeff[1] * 3 + 0) >> 3;
    iCoeff[1] = (iCoeff[0] >> 1) - iCoeff[1];
    iCoeff[0] -= iCoeff[1];
}

void FwdTOddOdd(int *iCoeff, int dummy)
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

void OverlapPreFilter4x4(int *iCoeff)
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

    void (*pointerFunct[17]) (int *arg1, int arg2) = {
        T2x2hEnc, T2x2hEnc, T2x2hEnc, T2x2hEnc,
        FwdScale, FwdScale, FwdScale, FwdScale,
        FwdRotate, FwdRotate, FwdRotate, FwdRotate,
        FwdTOddOdd, T2x2H, T2x2H, T2x2H, T2x2H
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

void OverlapPreFilter4(int *iCoeff)
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

void OverlapPreFilter2x2(int *iCoeff)
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

void OverlapPreFilter2(int *iCoeff)
{
    iCoeff[1] -= ((iCoeff[0] + 2) >> 2);
    iCoeff[0] -= (iCoeff[1] >> 13);
    iCoeff[0] -= (iCoeff[1] >> 9);
    iCoeff[0] -= (iCoeff[1] >> 5);
    iCoeff[0] -= ((iCoeff[1] + 1) >> 1);
    iCoeff[1] -= ((iCoeff[0] + 2) >> 2);
}

void EncFirstStagePreFiltering(int ** image, int numRows, int numCols)
{
    int macro_i, macro_j, block_i, block_j, i, j;
    int arrayLocal[16];

    for(macro_i = 0; macro_i < numRows/16; macro_i++)
    {
        for(macro_j = 0; macro_j < numCols/16; macro_j++)
        {
            for(block_i = 0; block_i < 4; block_i++)
            {
                for(block_j = 0; block_j < 4; block_j++)
                {
                    for (i = 0; i < 4; ++i)
                    {
                        for (j = 0; j < 4; ++j)
                        {
                            arrayLocal[i*4 + j] = image[macro_i*16 + block_i*4 + i][macro_j*16 + block_j*4 + j];
                        }
                    }
                    FCT4x4(arrayLocal);
                    for (i = 0; i < 4; ++i)
                    {
                        for (j = 0; j < 4; ++j)
                        {
                            image[macro_i*16 + block_i*4 + i][macro_j*16 + block_j*4 + j] = arrayLocal[i*4 + j];
                        }
                    }
                }
            } // end of block
        }
    }// end of macro block
}

void EncSecondStagePreFiltering(int ** image, int numRows, int numCols)
{
    int arrayLocal[16], i, j, macro_i, macro_j;
    for(macro_i = 0; macro_i < numRows/16; macro_i++)
    {
        for(macro_j = 0; macro_j < numCols/16; macro_j++)
        {
            for (i = 0; i < 4; ++i)
            {
                for (j = 0; j < 4; ++j)
                {
                    arrayLocal[i*4 + j] = image[macro_i*16 + i*4][macro_j*16 + j*4];
                }
            }
            FCT4x4(arrayLocal);
            for (i = 0; i < 4; ++i)
            {
                for (j = 0; j < 4; ++j)
                {
                    image[macro_i*16 + i*4][macro_j*16 + j*4] = arrayLocal[i*4 + j];
                }
            }
        }
    }
}

void EncFirstStageOverlapFilter(int** image, int numRows, int numCols)
{
    int arrayLocal_16[16], arrayLocal_4[4];
    int block_i, block_j, i, j;

    numRows /= 4;
    numCols /= 4;

    // 4x4 blocks
    for(block_i = 2; block_i < numRows-2; block_i += 4)
    {
        for(block_j = 2; block_j < numCols-2; block_j += 4)
        {
            for( i = 0; i < 4; i++)
            {
                for( j = 0; j < 4; j++)
                    arrayLocal_16[i*4+j] = image[(block_i + i)*4][(block_j + j)*4];
            }
            OverlapPreFilter4x4(arrayLocal_16);
            for( i = 0; i < 4; i++)
            {
                for( j = 0; j < 4; j++)
                    image[(block_i + i)*4][(block_j + j)*4] = arrayLocal_16[i*4+j];
            }
        }
    }
    //4x4 block end

    for(block_i = 2; block_i < numRows-2; block_i += 4)
    {
        //left edge
        for(i = 0; i < 2; i++)
        {
            for(j = 0; j < 4; j++)
                arrayLocal_4[j] = image[(block_i+j)*4][i*4];
            OverlapPreFilter4(arrayLocal_4);
            for(j = 0; j < 4; j++)
                image[(block_i+j)*4][i*4] = arrayLocal_4[j];
        }
        // right edge
        for(i = numCols-2; i < numCols-2; i++)
        {
            for(j = 0; j < 4; j++)
                arrayLocal_4[j] = image[(block_i+j)*4][i*4];
            OverlapPreFilter4(arrayLocal_4);
            for(j = 0; j < 4; j++)
                image[(block_i+j)*4][i*4] = arrayLocal_4[j];
        }
    }


    for(block_j = 2; block_j < numCols-2; block_j += 4)
    {
        //top edge
        for(i = 0; i < 2; i++)
        {
            for(j = 0; j < 4; j++)
                arrayLocal_4[j] = image[i*4][(block_j+j)*4];
            OverlapPreFilter4(arrayLocal_4);
            for(j = 0; j < 4; j++)
                image[i*4][(block_j+j)*4] = arrayLocal_4[j];
        }
        //bottom edge
        for(i = numRows-2; i < numRows; i++)
        {
            for(j = 0; j < 4; j++)
                arrayLocal_4[j] = image[i*4][(block_j+j)*4];
            OverlapPreFilter4(arrayLocal_4);
            for(j = 0; j < 4; j++)
                image[i*4][(block_j+j)*4] = arrayLocal_4[j];
        }
    }

    // top left
    arrayLocal_4[0] = image[0][0], arrayLocal_4[1] = image[0][4];
    arrayLocal_4[2] = image[4][0], arrayLocal_4[3] = image[4][4];
    OverlapPreFilter4(arrayLocal_4);
    image[0][0] = arrayLocal_4[0], image[0][4] = arrayLocal_4[1];
    image[4][0] = arrayLocal_4[2], image[4][4] = arrayLocal_4[3];

    // top right
    arrayLocal_4[0] = image[0][(numCols-2)*4], arrayLocal_4[1] = image[0][(numCols-1)*4];
    arrayLocal_4[2] = image[4][(numCols-2)*4], arrayLocal_4[3] = image[4][(numCols-1)*4];
    OverlapPreFilter4(arrayLocal_4);
    image[0][(numCols-2)*4] = arrayLocal_4[0], image[0][(numCols-1)*4] = arrayLocal_4[1];
    image[4][(numCols-2)*4] = arrayLocal_4[2], image[4][(numCols-1)*4] = arrayLocal_4[3];

    // bottom left
    arrayLocal_4[0] = image[(numRows-2)*4][0], arrayLocal_4[1] = image[(numRows-2)*4][4];
    arrayLocal_4[2] = image[(numRows-1)*4][0], arrayLocal_4[3] = image[(numRows-1)*4][4];
    OverlapPreFilter4(arrayLocal_4);
    image[(numRows-2)*4][0] = arrayLocal_4[0], image[(numRows-2)*4][4] = arrayLocal_4[1];
    image[(numRows-1)*4][0] = arrayLocal_4[2], image[(numRows-1)*4][4] = arrayLocal_4[3];

    // bottom right
    arrayLocal_4[0] = image[(numRows-2)*4][(numCols-2)*4], arrayLocal_4[1] = image[(numRows-2)*4][(numCols-1)*4];
    arrayLocal_4[2] = image[(numRows-1)*4][(numCols-2)*4], arrayLocal_4[3] = image[(numRows-1)*4][(numCols-1)*4];
    OverlapPreFilter4(arrayLocal_4);
    image[(numRows-2)*4][(numCols-2)*4] = arrayLocal_4[0], image[(numRows-2)*4][(numCols-1)*4] = arrayLocal_4[1];
    image[(numRows-1)*4][(numCols-2)*4] = arrayLocal_4[2], image[(numRows-1)*4][(numCols-1)*4] = arrayLocal_4[3];
}

void EncSecondStageOverlapFilter(int **image, int numRows, int numCols)
{
    int arrayLocal_16[16], arrayLocal_4[4];
    int block_i, block_j, i, j;
    // 4x4 blocks
    for(block_i = 2; block_i < numRows-2; block_i += 4)
    {
        for(block_j = 2; block_j < numCols-2; block_j += 4)
        {
            for( i = 0; i < 4; i++)
            {
                for( j = 0; j < 4; j++)
                    arrayLocal_16[i*4+j] = image[block_i + i][block_j + j];
            }
            OverlapPreFilter4x4(arrayLocal_16);
            for( i = 0; i < 4; i++)
            {
                for( j = 0; j < 4; j++)
                    image[block_i + i][block_j + j] = arrayLocal_16[i*4+j];
            }
        }
    }
    //4x4 block end

    for(block_i = 2; block_i < numRows-2; block_i += 4)
    {
        //left edge
        for(i = 0; i < 2; i++)
        {
            for(j = 0; j < 4; j++)
                arrayLocal_4[j] = image[block_i+j][i];
            OverlapPreFilter4(arrayLocal_4);
            for(j = 0; j < 4; j++)
                image[block_i+j][i] = arrayLocal_4[j];
        }
        // right edge
        for(i = numCols-2; i < numCols; i++)
        {
            for(j = 0; j < 4; j++)
                arrayLocal_4[j] = image[block_i+j][i];
            OverlapPreFilter4(arrayLocal_4);
            for(j = 0; j < 4; j++)
                image[block_i+j][i] = arrayLocal_4[j];
        }
    }


    for(block_j = 2; block_j < numCols-2; block_j += 4)
    {
        //top edge
        for(i = 0; i < 2; i++)
        {
            for(j = 0; j < 4; j++)
                arrayLocal_4[j] = image[i][block_j+j];
            OverlapPreFilter4(arrayLocal_4);
            for(j = 0; j < 4; j++)
                image[i][block_j+j] = arrayLocal_4[j];
        }
        //bottom edge
        for(i = numRows-2; i < numRows; i++)
        {
            for(j = 0; j < 4; j++)
                arrayLocal_4[j] = image[i][block_j+j];
            OverlapPreFilter4(arrayLocal_4);
            for(j = 0; j < 4; j++)
                image[i][block_j+j] = arrayLocal_4[j];
        }
    }

    // top left
    arrayLocal_4[0] = image[0][0], arrayLocal_4[1] = image[0][1];
    arrayLocal_4[2] = image[1][0], arrayLocal_4[3] = image[1][1];
    OverlapPreFilter4(arrayLocal_4);
    image[0][0] = arrayLocal_4[0], image[0][1] = arrayLocal_4[1];
    image[1][0] = arrayLocal_4[2], image[1][1] = arrayLocal_4[3];

    // top right
    arrayLocal_4[0] = image[0][numCols-2], arrayLocal_4[1] = image[0][numCols-1];
    arrayLocal_4[2] = image[1][numCols-2], arrayLocal_4[3] = image[1][numCols-1];
    OverlapPreFilter4(arrayLocal_4);
    image[0][numCols-2] = arrayLocal_4[0], image[0][numCols-1] = arrayLocal_4[1];
    image[1][numCols-2] = arrayLocal_4[2], image[1][numCols-1] = arrayLocal_4[3];

    // bottom left
    arrayLocal_4[0] = image[numRows-2][0], arrayLocal_4[1] = image[numRows-2][1];
    arrayLocal_4[2] = image[numRows-1][0], arrayLocal_4[3] = image[numRows-1][1];
    OverlapPreFilter4(arrayLocal_4);
    image[numRows-2][0] = arrayLocal_4[0], image[numRows-2][1] = arrayLocal_4[1];
    image[numRows-1][0] = arrayLocal_4[2], image[numRows-1][1] = arrayLocal_4[3];

    // bottom right
    arrayLocal_4[0] = image[numRows-2][numCols-2], arrayLocal_4[1] = image[numRows-2][numCols-1];
    arrayLocal_4[2] = image[numRows-1][numCols-2], arrayLocal_4[3] = image[numRows-1][numCols-1];
    OverlapPreFilter4(arrayLocal_4);
    image[numRows-2][numCols-2] = arrayLocal_4[0], image[numRows-2][numCols-1] = arrayLocal_4[1];
    image[numRows-1][numCols-2] = arrayLocal_4[2], image[numRows-1][numCols-1] = arrayLocal_4[3];
}

int main()
{
  struct timeval tim[50];
  int t = 0;
        
  FILE *stg1 = fopen("EncoderStage1.log", "w");
  FILE *stg2 = fopen("EncoderStage2.log", "w");
  FILE *stg3 = fopen("EncoderStage3.log", "w");
  FILE *stg4 = fopen("EncoderStage4.log", "w");
  
  int imageWidth,imageHeight;
  //scanf("%d %d", &imageHeight, &imageWidth);
  for(imageHeight = 256, imageWidth = 512; imageHeight <= 8000 && imageWidth <= 16000; imageHeight += 256, imageWidth += 512 )
  {
      t = 0;
      // read image in host
      FILE *ip = fopen("BImage", "r");
      FILE *op = fopen("encoded.txt", "w");
      
      int **image = (int**) malloc(imageHeight * sizeof(int*) );
      int i, j;
      for(i = 0; i < imageHeight; i++){
          image[i] = (int*) malloc(imageWidth * sizeof(int) );
    
          for(j = 0; j < imageWidth; j++)
              fscanf(ip, "%d", &image[i][j]);
      }
    
      /* kernel invocation start */
        gettimeofday(&tim[t], NULL); t++;
      // second stage frequency transform
      EncSecondStageOverlapFilter(image, imageHeight, imageWidth);
        gettimeofday(&tim[t], NULL); t++;
      // first stage pre-filtering
      EncFirstStagePreFiltering(image, imageHeight, imageWidth);
        gettimeofday(&tim[t], NULL); t++;
      // first stage frequency transform
      EncFirstStageOverlapFilter(image, imageHeight, imageWidth);
        gettimeofday(&tim[t], NULL); t++;
      // second stage pre-filtering
      EncSecondStagePreFiltering(image, imageHeight, imageWidth);
        gettimeofday(&tim[t], NULL); t++;
      /* kernel function invocation end */
    
      //store processed image in file
      for( i = 0; i < imageHeight; i++)
      {
          for( j = 0; j < imageWidth; j++)
              fprintf(op, "%d ", image[i][j] );
          fprintf(op, "\n");
      }
      
      i = 0;
      fprintf(stg1 ,"%d %lu %lu\n", imageHeight*imageWidth, tim[i+1].tv_sec - tim[i].tv_sec, tim[i+1].tv_usec - tim[i].tv_usec );
      i++;
      fprintf(stg2 ,"%d %lu %lu\n", imageHeight*imageWidth, tim[i+1].tv_sec - tim[i].tv_sec, tim[i+1].tv_usec - tim[i].tv_usec );
      i++;
      fprintf(stg3 ,"%d %lu %lu\n", imageHeight*imageWidth, tim[i+1].tv_sec - tim[i].tv_sec, tim[i+1].tv_usec - tim[i].tv_usec );
      i++;
      fprintf(stg4 ,"%d %lu %lu\n", imageHeight*imageWidth, tim[i+1].tv_sec - tim[i].tv_sec, tim[i+1].tv_usec - tim[i].tv_usec );
      
      fclose(ip);
      fclose(op);
      free(image);
      
    }
    fclose(stg1); fclose(stg2); fclose(stg3); fclose(stg4);
    return 0;
}
