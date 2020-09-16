
#include"builder.h"

int one_two_init(int *B1, int *B2, int atoms, int bonds);
void one_two_build(int **one_two, int *B1, int *B2, int atoms, int bonds, int max_1_2);

void topo_neighbors(int general_atoms, int general_bonds, int *master_B1_core_array, int *master_B2_core_array,
                    
                    int ***one_two, int ***one_three, int ***one_four)
{
    int max_1_2,i,j,k,l,m,n,o;
    int one_two_temp_list[4];
    int one_three_temp_list[12];
    int one_four_temp_list[36];
    
    //max_1_2=one_two_init(master_B1_core_array,master_B2_core_array,general_atoms,general_bonds);
    max_1_2=4;
    
    *one_two=(int**)malloc(general_atoms*sizeof(int*));for(i=0;i<general_atoms;++i)(*one_two)[i]=(int*)malloc(max_1_2*sizeof(int));
    one_two_build(*one_two,master_B1_core_array,master_B2_core_array,general_atoms,general_bonds,max_1_2);
    *one_three=(int**)malloc(general_atoms*sizeof(int*));for(i=0;i<general_atoms;++i)(*one_three)[i]=(int*)malloc(12*sizeof(int));
    *one_four=(int**)malloc(general_atoms*sizeof(int*));for(i=0;i<general_atoms;++i)(*one_four)[i]=(int*)malloc(36*sizeof(int));
    for(i=0;i<general_atoms;++i){for(j=0;j<12;++j)(*one_three)[i][j]=0;for(j=0;j<36;++j)(*one_four)[i][j]=0;}
    
    for(i=0;i<general_atoms;++i)
    {
        
        //
        
        for(l=0;l<4;++l)one_two_temp_list[l]=0;
        for(l=0;l<12;++l)one_three_temp_list[l]=0;
        for(l=0;l<36;++l)one_four_temp_list[l]=0;
        
        for(j=0;j<max_1_2;++j)
        {
            one_two_temp_list[j]=(*one_two)[i][j];
        }
        /*
         printf("%d 1-2 neighbors:\n",i+1);
         for(j=0;j<4;++j)
         //if(one_two_temp_list[j]!=0)
         printf("%d\n",one_two_temp_list[j]);
         */
        l=-1;
        for(j=0;j<4;++j)
        {
            if(one_two_temp_list[j]!=0)
            {
                for(k=0;k<max_1_2;++k)
                {
                    if((*one_two)[one_two_temp_list[j]-1][k]!=0 && (*one_two)[one_two_temp_list[j]-1][k]!=i+1)
                    {
                        l=l+1;
                        one_three_temp_list[l]=(*one_two)[one_two_temp_list[j]-1][k];
                    }
                }
            }
        }
        
        for(j=0;j<12-1;++j)
        {
            for(k=j+1;k<12;++k)
            {
                if(one_three_temp_list[k]==one_three_temp_list[j])
                {
                    one_three_temp_list[k]=0;
                }
            }
        }

        /*
         printf("%d 1-3 neighbors:\n",i+1);
         for(j=0;j<12;++j)
         //if(one_three_temp_list[j]!=0)
         printf("%d\n",one_three_temp_list[j]);
         */
        k=-1;
        for(j=0;j<12;++j){
            if(one_three_temp_list[j]!=0){
                k=k+1;(*one_three)[i][k]=one_three_temp_list[j];}}
        o=-1;
        for(j=0;j<12;++j)
        {
            if(one_three_temp_list[j]!=0)
            {
                for(k=0;k<max_1_2;++k)
                {
                    l=(*one_two)[one_three_temp_list[j]-1][k];
                    n=0;
                    for(m=0;m<4;++m)
                    {
                        if(l==one_two_temp_list[m])
                        {
                            n=1;
                            break;
                        }
                    }
                    if(n==0)
                    {
                        o=o+1;
                        one_four_temp_list[o]=l;
                    }
                }
            }
        }
        
        for(j=0;j<36-1;++j)
        {
            for(k=j+1;k<36;++k)
            {
                if(one_four_temp_list[k]==one_four_temp_list[j])
                {
                    one_four_temp_list[k]=0;
                }
            }
        }
        
        /*
         printf("%d 1-4 neighbors:\n",i+1);
         for(j=0;j<36;++j)
         //if(one_four_temp_list[j]!=0)
         printf("%d\n",one_four_temp_list[j]);
         */
        k=-1;
        for(j=0;j<36;++j){
            if(one_four_temp_list[j]!=0){
                k=k+1;(*one_four)[i][k]=one_four_temp_list[j];}}
        
        
        
        
        
        
    }
    /*
    for(i=0;i<general_atoms;++i){
        printf("[%d]\t",i+1);
        for(j=0;j<4;++j)if((*one_two)[i][j]!=0)printf("%d\t",(*one_two)[i][j]);
        printf("|\t");
        for(j=0;j<12;++j)if((*one_three)[i][j]!=0)printf("%d\t",(*one_three)[i][j]);
        printf("|\t");
        for(j=0;j<36;++j)if((*one_four)[i][j]!=0)printf("%d\t",(*one_four)[i][j]);
        printf("\n");}
     */




}
