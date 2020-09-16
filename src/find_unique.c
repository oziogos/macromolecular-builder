
#include"builder.h"

void find_unique(int N, int M, int **array,int **out)
{
    int i,j,k,sum;
    int *init,*final;
    
    int equal;
    
    init=(int*)malloc(N*sizeof(int));
    final=(int*)malloc(N*sizeof(int));
    
    for(i=0;i<N;++i)
    {
        init[i]=i+1;
        final[i]=i+1;
    }
    
    for(i=0;i<N-1;++i)
    {
        for(j=i+1;j<N;++j)
        {
            sum=0;
            for(k=0;k<M;++k){equal=0;if(array[i][k]==array[j][k])equal=1;sum=sum+equal;}
            if(sum==M)
            {
                final[j]=-1;
            }
        }
    }
    
    j=0;
    for(i=0;i<N;++i)
    {
        if(final[i]!=-1)
        {
            j=j+1;
            final[i]=j;
        }
    }
    
    for(i=0;i<N;++i)
        out[i][2]=final[i];
    
    for(i=0;i<N;++i)
    {
        if(final[i]==-1)
        {
            for(j=0;j<N;++j)
            {
                if(final[j]!=-1)
                {
                    sum=0;
                    for(k=0;k<M;++k){equal=0;if(array[i][k]==array[j][k])equal=1;sum=sum+equal;}
                    if(sum==M)
                    {
                        final[i]=final[j];
                        break;
                    }
                }
            }
        }
    }
    
    for(i=0;i<N;++i)
    {
        out[i][0]=init[i];
        out[i][1]=final[i];
    }
    
    free(init);free(final);
    
}