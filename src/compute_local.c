
#include"builder.h"

void compute_local(char *current_folder,double *segment_x_array,double *segment_y_array,double *segment_z_array,
                   int atoms_block_size,int bonds_block_size,int *local_topo_array,char **segment_species_array,
                   int *segment_B1_internal_core_array,int *segment_B2_internal_core_array,char **segment_B_type_core_array,
                   
                   int *local_atoms,int *local_bonds,char ***local_species,int **local_B1,int **local_B2,char ***local_B_type,int **local_map
                   )
{
    FILE *fp;
    char file_path[cmax_length];
    int j,k;
    
    //
    /*
    // xyz preview: open an appending file pointer
    sprintf(file_path,"%s/%s",current_folder,"local.dat");
    fp=fopen(file_path,"a");
    */
    // count the number of local atoms
    *local_atoms=0;
    //for(j=0;j<atoms_block_size;++j)if(local_topo_array[j]==1)*local_atoms=*local_atoms+1;
    for(j=0;j<atoms_block_size;++j)if(local_topo_array[j]!=0)*local_atoms=*local_atoms+1;
    // here we preallocate two arrays:
    // int* local_map has length local_atoms and will store the internal atomic (original) ids of the local atoms
    // char** local_species will hold the species information
    *local_map=(int*)malloc(*local_atoms*sizeof(int));
    *local_species=(char**)malloc(*local_atoms*sizeof(char*));
    for(j=0;j<*local_atoms;++j)(*local_species)[j]=(char*)malloc(sub_length*sizeof(char));
    // populate local_map and local_species
    k=-1;
    for(j=0;j<atoms_block_size;++j)
    {
        //if(local_topo_array[j]==1)
        if(local_topo_array[j]!=0)
        {
            k=k+1;
            sprintf((*local_species)[k],"%s",segment_species_array[j]);
            (*local_map)[k]=j+1;
        }
    }
    // we now need to find out how many bonds connect the local atoms
    *local_bonds=0;
    for(j=0;j<bonds_block_size;++j)
    {
        //if(local_topo_array[segment_B1_internal_core_array[j]-1]==1 && local_topo_array[segment_B2_internal_core_array[j]-1]==1)
        if(local_topo_array[segment_B1_internal_core_array[j]-1]!=0 && local_topo_array[segment_B2_internal_core_array[j]-1]!=0)
        {
            *local_bonds=*local_bonds+1;
        }
    }
    // prealloc and populate local bonds arrays: int* local_B1, local_B2
    *local_B1=(int*)malloc(*local_bonds*sizeof(int));
    *local_B2=(int*)malloc(*local_bonds*sizeof(int));
    *local_B_type=(char**)malloc(*local_bonds*sizeof(char*));for(j=0;j<*local_bonds;++j)(*local_B_type)[j]=(char*)malloc(sub_length*sizeof(char));
    k=-1;
    for(j=0;j<bonds_block_size;++j)
    {
        //if(local_topo_array[segment_B1_internal_core_array[j]-1]==1 && local_topo_array[segment_B2_internal_core_array[j]-1]==1)
        if(local_topo_array[segment_B1_internal_core_array[j]-1]!=0 && local_topo_array[segment_B2_internal_core_array[j]-1]!=0)
        {
            k=k+1;
            (*local_B1)[k]=segment_B1_internal_core_array[j];
            (*local_B2)[k]=segment_B2_internal_core_array[j];
            sprintf((*local_B_type)[k],"%s",segment_B_type_core_array[j]);
        }
    }
    /*
    // xyz preview: write to file
    fprintf(fp,"%d\n\n",*local_atoms);
    for(j=0;j<atoms_block_size;++j)
    {
        if(local_topo_array[j]==1)
        {
            fprintf(fp,"%c\t%lf\t%lf\t%lf\n",segment_species_array[j][0],segment_x_array[j],segment_y_array[j],segment_z_array[j]);
        }
    }
    fclose(fp);
    */
    /*
    // console output
    printf("** local_map:\n");
    for(j=0;j<*local_atoms;++j)printf("[%d\t%d]\t%s\n",j+1,(*local_map)[j],(*local_species)[j]);
    printf("** local bonds BEFORE mapping:\n");
    for(j=0;j<*local_bonds;++j)printf("[%d]\t%d\t%d\t%s\n",j+1,(*local_B1)[j],(*local_B2)[j],(*local_B_type)[j]);
    */
    // apply mapping: from original internal ids to consistent internal ids (starting from 1)
    for(j=0;j<*local_bonds;++j)
    {
        for(k=0;k<*local_atoms;++k)if((*local_map)[k]==(*local_B1)[j])break;(*local_B1)[j]=k+1;
        for(k=0;k<*local_atoms;++k)if((*local_map)[k]==(*local_B2)[j])break;(*local_B2)[j]=k+1;
        
    }
    // console output
    //printf("** local bonds AFTER mapping:\n");
    //for(j=0;j<*local_bonds;++j)printf("[%d]\t%d\t%d\t%s\n",j+1,(*local_B1)[j],(*local_B2)[j],(*local_B_type)[j]);
    //getchar();
    
}
