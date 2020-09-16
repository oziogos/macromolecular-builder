
#include"builder.h"

void initialize_tacticity(int master_species_rows, char **master_atom_name_chain_array, int mol_types, int *MA_start, int *MA_stop, int tacticity, int molecules, char **tacticity_id,
                          int *mol_type_array, char **tacticity_type,
                          
                          int **yet_one_more_t_array, int **tacticity_counter_array, int **tacticity_counter_array_v2, char ***tacticity_search, int *tacticity_array_length,
                          double **tacticity_array,
                          int **t_tracker_1,int **t_tracker_2,int **t_tracker_3,int **t_tracker_4,int **t_tracker_mol
                          )
{
    int i,j,k,l,sum,int_buffer;
    char word[cmax_length];

    *yet_one_more_t_array=(int*)malloc(molecules*sizeof(int));
    *tacticity_counter_array=(int*)malloc(molecules*sizeof(int));
    
    *tacticity_counter_array_v2=(int*)malloc(molecules*sizeof(int));
    
    for(i=0;i<molecules;++i)(*tacticity_counter_array)[i]=0;
    
    //for(i=0;i<tacticity;++i)printf("[%d]\t%s\n",i+1,tacticity_id[i]);
    // search for tacticity_id in master_atom_name_chain_array
    *tacticity_search=(char**)malloc(mol_types*sizeof(char*));
    for(i=0;i<mol_types;++i)(*tacticity_search)[i]=(char*)malloc(cmax_length*sizeof(char));
    for(i=0;i<mol_types;++i)sprintf((*tacticity_search)[i],"%s","X 0 X");
    for(i=0;i<mol_types;++i)
    {
        l=0;
        for(j=MA_start[i]-1;j<MA_stop[i];++j)
        {
            for(k=0;k<tacticity;++k)
            {
                if(strcmp(master_atom_name_chain_array[j],tacticity_id[k])==0)
                {
                    l=l+1;
                    sprintf((*tacticity_search)[i],"%s %d %s",master_atom_name_chain_array[j],l,tacticity_type[k]);
                }
            }
        }
    }
    //for(i=0;i<mol_types;++i)printf("%s\n",(*tacticity_search)[i]);
    //for(i=0;i<molecules;++i)printf("%d\n",(*tacticity_counter_array)[i]);
    
    //
    *tacticity_array_length=0;
    
    for(i=0;i<molecules;++i)(*yet_one_more_t_array)[i]=-1;
    
    sum=0;
    
    for(i=0;i<molecules;++i)
    {
        //printf("molecule %d --tacticity--> %s\n",i+1,(*tacticity_search)[mol_type_array[i]-1]);
        sscanf((*tacticity_search)[mol_type_array[i]-1],"%s\t%d",word,&int_buffer);
        
        sum=sum+int_buffer;
        //printf("%d\n",sum);
        
        if(int_buffer!=0){
            *tacticity_array_length=*tacticity_array_length+int_buffer;
            (*yet_one_more_t_array)[i]=sum-int_buffer;
            
            
        }
    }
    //printf("tacticity array length = %d\n",*tacticity_array_length);
    
    *tacticity_array=(double*)malloc((*tacticity_array_length)*sizeof(double));
    
    *t_tracker_1=(int*)malloc(*tacticity_array_length*sizeof(int));
    *t_tracker_2=(int*)malloc(*tacticity_array_length*sizeof(int));
    *t_tracker_3=(int*)malloc(*tacticity_array_length*sizeof(int));
    *t_tracker_4=(int*)malloc(*tacticity_array_length*sizeof(int));
    *t_tracker_mol=(int*)malloc(*tacticity_array_length*sizeof(int));
    
    for(i=0;i<*tacticity_array_length;++i){
        (*tacticity_array)[i]=0.0;
        (*t_tracker_1)[i]=0;
        (*t_tracker_2)[i]=0;
        (*t_tracker_3)[i]=0;
        (*t_tracker_4)[i]=0;
        (*t_tracker_mol)[i]=0;
    }
    
    //for(i=0;i<molecules;++i)printf("%d\n",(*yet_one_more_t_array)[i]);
    
    //for(i=0;i<molecules;++i)printf("molecule %d --tacticity--> %s\n",i+1,(*tacticity_search)[mol_type_array[i]-1]);
    //getchar();

    
}
