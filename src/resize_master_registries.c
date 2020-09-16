
#include"builder.h"

void shrink_2D_char(char ***input,int rows,int cols,int remove);

void resize_master_registries(int *atoms,int new_entries_atom_total,int atoms_last_entry,
                              int **master_id_array,int **master_mol_array,int **host_id_array,int **chain_id_array,int **master_id_internal_array,int **master_mol_type_array,
                              double **master_q_array,double **master_x_array,double **master_y_array,double **master_z_array,
                              char ***master_species_array,
                              int *bonds,int new_entries_bond_total,int bonds_last_entry,
                              int **master_B_ID_core_array,int **master_B1_core_array,int **master_B2_core_array,int **mol_id_array_BONDS,
                              int **master_B1_internal_core_array,int **master_B2_internal_core_array,
                              char ***master_B_type_core_array,
                              char ***master_atom_name_array)
{

    int j;
    
    // push down entries
    if(new_entries_atom_total>0){
        
        /*
        printf("%d\t%d\n",*atoms,new_entries_atom_total);
        for(j=0;j<*atoms;++j)printf("%d\t$\n",j);
        for(j=1;j<=new_entries_atom_total;++j)printf("%d\n",*atoms-1+j);
        printf("%d\t<--\n",*atoms+new_entries_atom_total-1);
        getchar();
        */
        
        // realloc: atoms
        *master_id_array=(int*)realloc(*master_id_array,(*atoms+new_entries_atom_total)*sizeof(int));
        *master_mol_array=(int*)realloc(*master_mol_array,(*atoms+new_entries_atom_total)*sizeof(int));
        *host_id_array=(int*)realloc(*host_id_array,(*atoms+new_entries_atom_total)*sizeof(int));
        *chain_id_array=(int*)realloc(*chain_id_array,(*atoms+new_entries_atom_total)*sizeof(int));
        *master_id_internal_array=(int*)realloc(*master_id_internal_array,(*atoms+new_entries_atom_total)*sizeof(int));
        *master_mol_type_array=(int*)realloc(*master_mol_type_array,(*atoms+new_entries_atom_total)*sizeof(int));
        *master_q_array=(double*)realloc(*master_q_array,(*atoms+new_entries_atom_total)*sizeof(double));
        *master_x_array=(double*)realloc(*master_x_array,(*atoms+new_entries_atom_total)*sizeof(double));
        *master_y_array=(double*)realloc(*master_y_array,(*atoms+new_entries_atom_total)*sizeof(double));
        *master_z_array=(double*)realloc(*master_z_array,(*atoms+new_entries_atom_total)*sizeof(double));
        *master_species_array=(char**)realloc(*master_species_array,(*atoms+new_entries_atom_total)*sizeof(char*));
        for(j=1;j<=new_entries_atom_total;++j)(*master_species_array)[*atoms-1+j]=(char*)malloc(sub_length*sizeof(char));
        *master_atom_name_array=(char**)realloc(*master_atom_name_array,(*atoms+new_entries_atom_total)*sizeof(char*));
        for(j=1;j<=new_entries_atom_total;++j)(*master_atom_name_array)[*atoms-1+j]=(char*)malloc(sub_length*sizeof(char));
        if( atoms_last_entry != (*atoms-1) )
        {
            for(j=0;j<*atoms-atoms_last_entry-1;++j)
            {
                (*master_id_array)[*atoms+new_entries_atom_total-1-j]=(*master_id_array)[*atoms-1-j];
                (*master_mol_array)[*atoms+new_entries_atom_total-1-j]=(*master_mol_array)[*atoms-1-j];
                (*host_id_array)[*atoms+new_entries_atom_total-1-j]=(*host_id_array)[*atoms-1-j];
                (*chain_id_array)[*atoms+new_entries_atom_total-1-j]=(*chain_id_array)[*atoms-1-j];
                (*master_id_internal_array)[*atoms+new_entries_atom_total-1-j]=(*master_id_internal_array)[*atoms-1-j];
                (*master_mol_type_array)[*atoms+new_entries_atom_total-1-j]=(*master_mol_type_array)[*atoms-1-j];
                (*master_q_array)[*atoms+new_entries_atom_total-1-j]=(*master_q_array)[*atoms-1-j];
                (*master_x_array)[*atoms+new_entries_atom_total-1-j]=(*master_x_array)[*atoms-1-j];
                (*master_y_array)[*atoms+new_entries_atom_total-1-j]=(*master_y_array)[*atoms-1-j];
                (*master_z_array)[*atoms+new_entries_atom_total-1-j]=(*master_z_array)[*atoms-1-j];
                sprintf((*master_species_array)[*atoms+new_entries_atom_total-1-j],"%s",(*master_species_array)[*atoms-1-j]);
                sprintf((*master_atom_name_array)[*atoms+new_entries_atom_total-1-j],"%s",(*master_atom_name_array)[*atoms-1-j]);
            }
        }
    }
    else if(new_entries_atom_total<0)
    {
        if( atoms_last_entry != (*atoms-1) )
        {
            for(j=atoms_last_entry+1;j<*atoms;++j)
            {
                (*master_id_array)[j+new_entries_atom_total]=(*master_id_array)[j];
                (*master_mol_array)[j+new_entries_atom_total]=(*master_mol_array)[j];
                (*host_id_array)[j+new_entries_atom_total]=(*host_id_array)[j];
                (*chain_id_array)[j+new_entries_atom_total]=(*chain_id_array)[j];
                (*master_id_internal_array)[j+new_entries_atom_total]=(*master_id_internal_array)[j];
                (*master_mol_type_array)[j+new_entries_atom_total]=(*master_mol_type_array)[j];
                (*master_q_array)[j+new_entries_atom_total]=(*master_q_array)[j];
                (*master_x_array)[j+new_entries_atom_total]=(*master_x_array)[j];
                (*master_y_array)[j+new_entries_atom_total]=(*master_y_array)[j];
                (*master_z_array)[j+new_entries_atom_total]=(*master_z_array)[j];
                sprintf((*master_species_array)[j+new_entries_atom_total],"%s",(*master_species_array)[j]);
                sprintf((*master_atom_name_array)[j+new_entries_atom_total],"%s",(*master_atom_name_array)[j]);
            }
        }
        *master_id_array=(int*)realloc(*master_id_array,(*atoms+new_entries_atom_total)*sizeof(int));
        *master_mol_array=(int*)realloc(*master_mol_array,(*atoms+new_entries_atom_total)*sizeof(int));
        *host_id_array=(int*)realloc(*host_id_array,(*atoms+new_entries_atom_total)*sizeof(int));
        *chain_id_array=(int*)realloc(*chain_id_array,(*atoms+new_entries_atom_total)*sizeof(int));
        *master_id_internal_array=(int*)realloc(*master_id_internal_array,(*atoms+new_entries_atom_total)*sizeof(int));
        *master_mol_type_array=(int*)realloc(*master_mol_type_array,(*atoms+new_entries_atom_total)*sizeof(int));
        *master_q_array=(double*)realloc(*master_q_array,(*atoms+new_entries_atom_total)*sizeof(double));
        *master_x_array=(double*)realloc(*master_x_array,(*atoms+new_entries_atom_total)*sizeof(double));
        *master_y_array=(double*)realloc(*master_y_array,(*atoms+new_entries_atom_total)*sizeof(double));
        *master_z_array=(double*)realloc(*master_z_array,(*atoms+new_entries_atom_total)*sizeof(double));
        shrink_2D_char(master_species_array,*atoms,sub_length,-new_entries_atom_total);
        shrink_2D_char(master_atom_name_array,*atoms,sub_length,-new_entries_atom_total);
        
    }
    
    // augment
    *atoms=*atoms+new_entries_atom_total;
    // write ids
    for(j=0;j<*atoms;++j)(*master_id_array)[j]=j+1;
    
    // push down entries
    if(new_entries_bond_total>0)
    {
        // realloc: bonds
        // master_B_ID_core_array,master_B1_core_array,master_B2_core_array,master_B_type_core_array,
        // mol_id_array_BONDS,master_B1_internal_core_array,master_B2_internal_core_array
        *master_B_ID_core_array=(int*)realloc(*master_B_ID_core_array,(*bonds+new_entries_bond_total)*sizeof(int));
        *master_B1_core_array=(int*)realloc(*master_B1_core_array,(*bonds+new_entries_bond_total)*sizeof(int));
        *master_B2_core_array=(int*)realloc(*master_B2_core_array,(*bonds+new_entries_bond_total)*sizeof(int));
        *mol_id_array_BONDS=(int*)realloc(*mol_id_array_BONDS,(*bonds+new_entries_bond_total)*sizeof(int));
        *master_B1_internal_core_array=(int*)realloc(*master_B1_internal_core_array,(*bonds+new_entries_bond_total)*sizeof(int));
        *master_B2_internal_core_array=(int*)realloc(*master_B2_internal_core_array,(*bonds+new_entries_bond_total)*sizeof(int));
        *master_B_type_core_array=(char**)realloc(*master_B_type_core_array,(*bonds+new_entries_bond_total)*sizeof(char*));
        for(j=1;j<=new_entries_bond_total;++j)(*master_B_type_core_array)[*bonds-1+j]=(char*)malloc(sub_length*sizeof(char));
        if( bonds_last_entry != (*bonds-1) )
        {
            for(j=0;j<*bonds-bonds_last_entry-1;++j)
            {
                (*master_B_ID_core_array)[*bonds+new_entries_bond_total-1-j]=(*master_B_ID_core_array)[*bonds-1-j];
                (*master_B1_core_array)[*bonds+new_entries_bond_total-1-j]=(*master_B1_core_array)[*bonds-1-j];
                (*master_B2_core_array)[*bonds+new_entries_bond_total-1-j]=(*master_B2_core_array)[*bonds-1-j];
                (*mol_id_array_BONDS)[*bonds+new_entries_bond_total-1-j]=(*mol_id_array_BONDS)[*bonds-1-j];
                (*master_B1_internal_core_array)[*bonds+new_entries_bond_total-1-j]=(*master_B1_internal_core_array)[*bonds-1-j];
                (*master_B2_internal_core_array)[*bonds+new_entries_bond_total-1-j]=(*master_B2_internal_core_array)[*bonds-1-j];
                sprintf((*master_B_type_core_array)[*bonds+new_entries_bond_total-1-j],"%s",(*master_B_type_core_array)[*bonds-1-j]);
            }
        }
    }
    else if(new_entries_bond_total<0)
    {
        if( bonds_last_entry != (*bonds-1) )
        {
            for(j=bonds_last_entry+1;j<*bonds;++j)
            {
                (*master_B_ID_core_array)[j+new_entries_bond_total]=(*master_B_ID_core_array)[j];
                (*master_B1_core_array)[j+new_entries_bond_total]=(*master_B1_core_array)[j];
                (*master_B2_core_array)[j+new_entries_bond_total]=(*master_B2_core_array)[j];
                (*mol_id_array_BONDS)[j+new_entries_bond_total]=(*mol_id_array_BONDS)[j];
                (*master_B1_internal_core_array)[j+new_entries_bond_total]=(*master_B1_internal_core_array)[j];
                (*master_B2_internal_core_array)[j+new_entries_bond_total]=(*master_B2_internal_core_array)[j];
                sprintf((*master_B_type_core_array)[j+new_entries_bond_total],"%s",(*master_B_type_core_array)[j]);
            }
        }
        *master_B_ID_core_array=(int*)realloc(*master_B_ID_core_array,(*bonds+new_entries_bond_total)*sizeof(int));
        *master_B1_core_array=(int*)realloc(*master_B1_core_array,(*bonds+new_entries_bond_total)*sizeof(int));
        *master_B2_core_array=(int*)realloc(*master_B2_core_array,(*bonds+new_entries_bond_total)*sizeof(int));
        *mol_id_array_BONDS=(int*)realloc(*mol_id_array_BONDS,(*bonds+new_entries_bond_total)*sizeof(int));
        *master_B1_internal_core_array=(int*)realloc(*master_B1_internal_core_array,(*bonds+new_entries_bond_total)*sizeof(int));
        *master_B2_internal_core_array=(int*)realloc(*master_B2_internal_core_array,(*bonds+new_entries_bond_total)*sizeof(int));
        shrink_2D_char(master_B_type_core_array,*bonds,sub_length,-new_entries_bond_total);
        
    }
    // augment
    //*bonds=*bonds+new_entries_bond_total;
    // write ids
    //for(j=0;j<*bonds;++j)(*master_B_ID_core_array)[j]=j+1;
    for(j=0;j<(*bonds+new_entries_bond_total);++j)(*master_B_ID_core_array)[j]=j+1;

}

void shrink_2D_char(char ***input,int rows,int cols,int remove);
void shrink_2D_char(char ***input,int rows,int cols,int remove)
{
    char **buffer;
    int i,new_rows;
    
    buffer=(char**)malloc(rows*sizeof(char*));
    for(i=0;i<rows;++i)buffer[i]=(char*)malloc(cols*sizeof(char));
    
    for(i=0;i<rows;++i)sscanf((*input)[i],"%s",buffer[i]);
    
    for(i=0;i<rows;++i)free((*input)[i]);free(*input);
    
    new_rows=rows-remove;
    *input=(char**)malloc(new_rows*sizeof(char*));
    for(i=0;i<new_rows;++i)(*input)[i]=(char*)malloc(cols*sizeof(char));
    
    for(i=0;i<new_rows;++i)sprintf((*input)[i],"%s",buffer[i]);
    
    for(i=0;i<rows;++i)free(buffer[i]);free(buffer);
}
