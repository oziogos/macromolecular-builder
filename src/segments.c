
#include"builder.h"

void segments(int mol_id,
              int atoms,
              int *master_id_array,int *master_mol_array,char **master_species_array,double *master_q_array,double *master_x_array,double *master_y_array,double *master_z_array,
              int *host_id_array,int *chain_id_array,int *master_id_internal_array,int *master_mol_type_array,
              int *atoms_block_size,
              int **segment_id_array,int **segment_mol_array,char ***segment_species_array,double **segment_q_array,double **segment_x_array,double **segment_y_array,double **segment_z_array,
              int **segment_host_id_array,int **segment_chain_id_array,int **segment_id_internal_array,int **segment_mol_type_array,
              int bonds,
              int *master_B_ID_core_array,int *master_B1_core_array,int *master_B2_core_array,int *mol_id_array_BONDS,int *master_B1_internal_core_array,int *master_B2_internal_core_array,char **master_B_type_core_array,
              int *bonds_block_size,
              int **segment_B_ID_core_array,int **segment_B1_core_array,int **segment_B2_core_array,int **segment_mol_id_array_BONDS,int **segment_B1_internal_core_array,int **segment_B2_internal_core_array,char ***segment_B_type_core_array,
              char ***segment_atom_name_core,char ***segment_subst_name_core,int **segment_B_ID_core,
              int *atoms_first_entry,int *atoms_last_entry,int *bonds_first_entry,int *bonds_last_entry,
              char **master_atom_name_array
              )
{
    int j,k;
    int first_entry,last_entry;
    
    // atoms
    // find first entry
    first_entry=0;
    for(j=0;j<atoms;++j)
    {
        if(master_mol_array[j]==mol_id)
        {
            first_entry=j;
            break;
        }
    }
    // find last entry
    last_entry=0;
    for(j=first_entry;j<atoms;++j)
    {
        if(master_mol_array[j]!=mol_id)
        {
            last_entry=j-1;
            break;
        }
        if(j==atoms-1)
            last_entry=j;
    }
    // size
    *atoms_block_size=last_entry-first_entry+1;
    //
    //printf("%d\t%d\t%d\t%d\n",mol_id,first_entry,last_entry,*atoms_block_size);
    // return entries to main
    *atoms_first_entry=first_entry;
    *atoms_last_entry=last_entry;
    // preallocate arrays
    *segment_species_array=(char**)malloc(*atoms_block_size*sizeof(char*));for(j=0;j<*atoms_block_size;++j)(*segment_species_array)[j]=(char*)malloc(sub_length*sizeof(char));
    *segment_id_array=(int*)malloc(*atoms_block_size*sizeof(int));
    *segment_mol_array=(int*)malloc(*atoms_block_size*sizeof(int));
    *segment_host_id_array=(int*)malloc(*atoms_block_size*sizeof(int));
    *segment_chain_id_array=(int*)malloc(*atoms_block_size*sizeof(int));
    *segment_id_internal_array=(int*)malloc(*atoms_block_size*sizeof(int));
    *segment_mol_type_array=(int*)malloc(*atoms_block_size*sizeof(int));
    *segment_q_array=(double*)malloc(*atoms_block_size*sizeof(double));
    *segment_x_array=(double*)malloc(*atoms_block_size*sizeof(double));
    *segment_y_array=(double*)malloc(*atoms_block_size*sizeof(double));
    *segment_z_array=(double*)malloc(*atoms_block_size*sizeof(double));
    //
    *segment_atom_name_core=(char**)malloc(*atoms_block_size*sizeof(char*));for(j=0;j<*atoms_block_size;++j)(*segment_atom_name_core)[j]=(char*)malloc(sub_length*sizeof(char));
    *segment_subst_name_core=(char**)malloc(*atoms_block_size*sizeof(char*));for(j=0;j<*atoms_block_size;++j)(*segment_subst_name_core)[j]=(char*)malloc(sub_length*sizeof(char));

    // populate arrays
    k=-1;
    for(j=first_entry;j<=last_entry;++j)
    {
        k=k+1;
        (*segment_id_array)[k]=master_id_array[j];
        (*segment_mol_array)[k]=master_mol_array[j];
        (*segment_host_id_array)[k]=host_id_array[j];
        (*segment_chain_id_array)[k]=chain_id_array[j];
        (*segment_id_internal_array)[k]=master_id_internal_array[j];
        (*segment_mol_type_array)[k]=master_mol_type_array[j];
        (*segment_q_array)[k]=master_q_array[j];
        (*segment_x_array)[k]=master_x_array[j];
        (*segment_y_array)[k]=master_y_array[j];
        (*segment_z_array)[k]=master_z_array[j];
        sprintf((*segment_species_array)[k],"%s",master_species_array[j]);
        sprintf((*segment_atom_name_core)[k],"%s",master_atom_name_array[j]);
        sprintf((*segment_subst_name_core)[k],"LIG%d",master_mol_array[j]);
    }
    
    // bonds
    // find first entry
    first_entry=0;
    for(j=0;j<bonds;++j)
    {
        if(mol_id_array_BONDS[j]==mol_id)
        {
            first_entry=j;
            break;
        }
    }
    // find last entry
    last_entry=0;
    for(j=first_entry;j<bonds;++j)
    {
        if(mol_id_array_BONDS[j]!=mol_id)
        {
            last_entry=j-1;
            break;
        }
        if(j==bonds-1)
            last_entry=j;
    }
    // size
    *bonds_block_size=last_entry-first_entry+1;
    //
    //printf("%d\t%d\t%d\t%d\n",mol_id,first_entry,last_entry,*bonds_block_size);
    // return entries to main
    *bonds_first_entry=first_entry;
    *bonds_last_entry=last_entry;
    // preallocate arrays
    *segment_B_type_core_array=(char**)malloc(*bonds_block_size*sizeof(char*));for(j=0;j<*bonds_block_size;++j)(*segment_B_type_core_array)[j]=(char*)malloc(sub_length*sizeof(char));
    *segment_B_ID_core_array=(int*)malloc(*bonds_block_size*sizeof(int));
    *segment_B1_core_array=(int*)malloc(*bonds_block_size*sizeof(int));
    *segment_B2_core_array=(int*)malloc(*bonds_block_size*sizeof(int));
    *segment_mol_id_array_BONDS=(int*)malloc(*bonds_block_size*sizeof(int));
    *segment_B1_internal_core_array=(int*)malloc(*bonds_block_size*sizeof(int));
    *segment_B2_internal_core_array=(int*)malloc(*bonds_block_size*sizeof(int));
    //
    *segment_B_ID_core=(int*)malloc(*bonds_block_size*sizeof(int));
    // populate
    // populate arrays
    k=-1;
    for(j=first_entry;j<=last_entry;++j)
    {
        k=k+1;
        (*segment_B_ID_core_array)[k]=master_B_ID_core_array[j];
        (*segment_B1_core_array)[k]=master_B1_core_array[j];
        (*segment_B2_core_array)[k]=master_B2_core_array[j];
        (*segment_mol_id_array_BONDS)[k]=mol_id_array_BONDS[j];
        (*segment_B1_internal_core_array)[k]=master_B1_internal_core_array[j];
        (*segment_B2_internal_core_array)[k]=master_B2_internal_core_array[j];
        sprintf((*segment_B_type_core_array)[k],"%s",master_B_type_core_array[j]);
        //
        (*segment_B_ID_core)[k]=k+1;
    }


}