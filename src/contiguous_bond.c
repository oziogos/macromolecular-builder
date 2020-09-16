
#include"builder.h"

void contiguous_bonds(int atoms,int bonds,int molecules,int mol_types,int *mol_type_array,int *type_array,char *current_folder,char **core_file,
                      int *master_id_array,int *master_id_internal_array,int *master_mol_array,
                      char ***master_B_type_core_array,int **master_B_ID_core_array,int **master_B1_core_array,int **master_B2_core_array,int **mol_id_array_BONDS,
                      int **master_B1_internal_core_array,int **master_B2_internal_core_array)
{
    FILE *fp;
    char file_path[cmax_length],buffer[cmax_length];
    int i,j,k,l,int_buffer;
    int atoms_core,bonds_core;
    
    // preallocate arrays for bonding data
    *master_B_type_core_array=(char**)malloc(bonds*sizeof(char*));
    for(i=0;i<bonds;++i)(*master_B_type_core_array)[i]=(char*)malloc(sub_length*sizeof(char));
    *master_B1_core_array=(int*)malloc(bonds*sizeof(int));
    *master_B2_core_array=(int*)malloc(bonds*sizeof(int));
    *mol_id_array_BONDS=(int*)malloc(bonds*sizeof(int));
    *master_B1_internal_core_array=(int*)malloc(bonds*sizeof(int));
    *master_B2_internal_core_array=(int*)malloc(bonds*sizeof(int));
    *master_B_ID_core_array=(int*)malloc(bonds*sizeof(int));
    // read bonding data
    l=-1;
    for(i=0;i<molecules;++i)
    {
        for(j=0;j<mol_types;++j)
        {
            if(mol_type_array[i]==type_array[j])
                break;
        }
        sprintf(file_path,"%s/%s",current_folder,core_file[j]);
        fp=fopen(file_path,"r");
        if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
        while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"@<TRIPOS>MOLECULE\n")==0)break;          // locate the MOLECULE section
        fgets(buffer,cmax_length,fp);fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d",&atoms_core,&bonds_core);    // skip comment line, read atoms and bonds from the next one
        while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"@<TRIPOS>BOND\n")==0)break;          // locate the BOND section
        for(k=0;k<bonds_core;++k)
        {
            l=l+1;
            fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d\t%d\t%s",&int_buffer,&(*master_B1_core_array)[l],&(*master_B2_core_array)[l],(*master_B_type_core_array)[l]);
            (*mol_id_array_BONDS)[l]=i+1;
        }
        fclose(fp);
    }
    // backup bonds
    for(i=0;i<bonds;++i){(*master_B1_internal_core_array)[i]=(*master_B1_core_array)[i];(*master_B2_internal_core_array)[i]=(*master_B2_core_array)[i];}
    // normalize bond info indices
    for(i=0;i<bonds;++i)
    {
        for(j=0;j<atoms;++j)
        {
            if((*master_B1_core_array)[i]==master_id_internal_array[j] && (*mol_id_array_BONDS)[i]==master_mol_array[j])
            {
                (*master_B1_core_array)[i]=master_id_array[j];
                break;
            }
        }
        
        for(j=0;j<atoms;++j)
        {
            if((*master_B2_core_array)[i]==master_id_internal_array[j] && (*mol_id_array_BONDS)[i]==master_mol_array[j])
            {
                (*master_B2_core_array)[i]=master_id_array[j];
                break;
            }
        }
        
        //
        (*master_B_ID_core_array)[i]=i+1;
        
    }
}
