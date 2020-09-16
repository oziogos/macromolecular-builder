
#include"builder.h"

void permute(int N,int **res)
{
    int M,L,i,j,k,f,E,E1,E2,p;
    int *array,d,c;
    
    M=N;
    
    L=pow(M,N);
    array=(int*)malloc(N*sizeof(int));
    //printf("%d\n",L);
    c=0;
    for(i=0;i<L;++i)
    {
        k=-1;
        for(j=1;j<=N;++j)
        {
            f=(int)pow(M,N-j);
            // 1 + Floor[i/f[j, NN, MM]]
            E1=1+(int)floor(i/f);
            // MM - MM Floor[1 + Floor[i/f[j, NN, MM]]/MM]
            p=(int)floor(1+(int)floor(i/f)/M);
            E2=M-M*p;
            //
            E=E1+E2;
            k=k+1;
            array[k]=E;
            //printf("%d\t",E);
        }
        //printf("\n");

        d=0;
        for(j=0;j<N-1;++j)
        {
            for(k=j+1;k<N;++k)
            {
                if(array[j]==array[k])
                {
                    d=1;
                }
            }
        }
        if(d==0)
        {
            c=c+1;
            for(j=0;j<N;++j)
                res[c-1][j]=array[j];
                //printf("%d\t",array[j]);
            //printf("\n");
        }
    }
    
    //printf("%d\n",c);
    
    free(array);
}