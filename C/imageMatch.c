#include <stdio.h>

int main()
{
    FILE *ip = fopen("image.txt", "r");
    FILE *ip2 = fopen("test.txt", "r");
    FILE *op = fopen("checkFinal.txt", "w");

    int i, j, a,b,yo=0;
    for(i = 0; i < 128; i++)
    {
        for(j = 0; j < 112; j++)
        {
            fscanf(ip, "%d", &a);
            fscanf(ip2, "%d", &b);

            if(a == b)
              fprintf(op, "%d ", 0 );
            else{
              fprintf(op, "%d ", 1 );
              yo = 1;
            }
        }
        fprintf(op, "\n" );
    }
    fclose(op);
    printf("%d\n",yo );
    scanf("%d\n", &yo);
}
