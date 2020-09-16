
#include"builder.h"

void initialize_Cn(int molecules, int *max_growth_step_array, char **quad_matrix, int *mol_type_array,
                   int **bb_growth_step_array,int **added_backbone_flag, double **l_array, double **l2_array, int **Re_start, int **Re_stop, int **prev_link_atom_array)
{
    char word[cmax_length];
    int i,int_buffer;
    
    
    *bb_growth_step_array=(int*)malloc(molecules*sizeof(int));
    *added_backbone_flag=(int*)malloc(molecules*sizeof(int));
    *l_array=(double*)malloc(molecules*sizeof(double));
    *l2_array=(double*)malloc(molecules*sizeof(double));
    *Re_start=(int*)malloc(molecules*sizeof(int));
    *Re_stop=(int*)malloc(molecules*sizeof(int));
    *prev_link_atom_array=(int*)malloc(molecules*sizeof(int));
    for(i=0;i<molecules;++i)
    {
        (*bb_growth_step_array)[i]=0;
        (*Re_start)[i]=-1;(*Re_stop)[i]=-1;(*prev_link_atom_array)[i]=-1;
        if(max_growth_step_array[i]!=0)
        {
            sscanf(quad_matrix[mol_type_array[i]-1],"%d\t%s\t%d\t%d",&int_buffer,word,&int_buffer,&(*Re_start)[i]);
            //printf("%d\t%s\n",(*Re_start)[i],quad_matrix[mol_type_array[i]-1]);
        }
    }

    //for(i=0;i<molecules;++i)printf("[%d]\t%d\t|\t%d\t%d\n",i+1,mol_type_array[i],(*Re_start)[i],(*Re_stop)[i]);
    //getchar();
}
