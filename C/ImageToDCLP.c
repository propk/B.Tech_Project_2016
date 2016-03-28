#include <stdio.h>

int main()
{
  FILE *ip = fopen("thirdstg.txt", "r+");
  FILE *op = fopen("thirdstgDCLP.txt", "w+");
  int i, j, a;
  for(i = 0; i < 128; i+= 1)
  {
    for(j = 0; j < 112; j+= 1)
    {
      fscanf(ip, "%d", &a);
      if(i%4 == 0 && j%4 == 0)
        fprintf(op, "%d ", a);
    }
    if(i%4 == 0 && j%4 == 0)
    fprintf(op, "\n" );
  }
}
