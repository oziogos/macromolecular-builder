
#include"builder.h"

void update_master(int i, int atoms, int bonds, int atoms_block_size, int bonds_block_size, int atoms_first_entry, int bonds_first_entry, int *master_id_array,
                   int *master_mol_array, int *host_id_array, int *chain_id_array, int *master_id_internal_array, int *master_mol_type_array,
                   double *master_q_array, double *master_x_array, double *master_y_array, double *master_z_array, char **master_species_array,
                   char **master_B_type_core_array, int *mol_id_array_BONDS, int *master_B1_core_array, int *master_B2_core_array, int *master_B1_internal_core_array, int *master_B2_internal_core_array,
                   int *segment_mol_array, int *segment_host_id_array, int *segment_chain_id_array, int *segment_id_internal_array, int *mol_type_array,
                   double *segment_q_array, double *segment_x_array, double *segment_y_array, double *segment_z_array, char **segment_species_array,
                   char **segment_B_type_core_array, int *segment_B1_internal_core_array, int *segment_B2_internal_core_array,
                   int *atom_scaling_array,
                   char **master_atom_name_array,char **segment_atom_name_core)
{
    int j,k;
    
    // write atomic data in reserved space
    for(j=0;j<atoms_block_size;++j)
    {
        master_mol_array[atoms_first_entry+j]=segment_mol_array[j];
        host_id_array[atoms_first_entry+j]=segment_host_id_array[j];
        chain_id_array[atoms_first_entry+j]=segment_chain_id_array[j];
        master_id_internal_array[atoms_first_entry+j]=segment_id_internal_array[j];
        master_mol_type_array[atoms_first_entry+j]=mol_type_array[i];
        master_q_array[atoms_first_entry+j]=segment_q_array[j];
        master_x_array[atoms_first_entry+j]=segment_x_array[j];
        master_y_array[atoms_first_entry+j]=segment_y_array[j];
        master_z_array[atoms_first_entry+j]=segment_z_array[j];
        sprintf(master_species_array[atoms_first_entry+j],"%s",segment_species_array[j]);
        sprintf(master_atom_name_array[atoms_first_entry+j],"%s",segment_atom_name_core[j]);
    }
    // update bond data
    for(j=0;j<bonds_block_size;++j)
    {
        sprintf(master_B_type_core_array[j+bonds_first_entry],"%s",segment_B_type_core_array[j]);
        mol_id_array_BONDS[j+bonds_first_entry]=segment_mol_array[0];
        master_B1_internal_core_array[j+bonds_first_entry]=segment_B1_internal_core_array[j];
        master_B2_internal_core_array[j+bonds_first_entry]=segment_B2_internal_core_array[j];
        //master_B1_core_array[j+bonds_first_entry]=segment_B1_internal_core_array[j]+atom_scaling_array[i];
        //master_B2_core_array[j+bonds_first_entry]=segment_B2_internal_core_array[j]+atom_scaling_array[i];
    }
    // all
    for(j=0;j<bonds;++j)
    {
        master_B1_core_array[j]=master_B1_internal_core_array[j]+atom_scaling_array[mol_id_array_BONDS[j]-1];
        master_B2_core_array[j]=master_B2_internal_core_array[j]+atom_scaling_array[mol_id_array_BONDS[j]-1];
    }
    
    /*
    for(j=0;j<bonds;++j)
    {
        for(k=0;k<atoms;++k)
        {
            if(master_mol_array[k]==mol_id_array_BONDS[j] && master_id_internal_array[k]==master_B1_internal_core_array[j])
            {
                master_B1_core_array[j]=master_id_array[k];
                break;
            }
            
        }
        for(k=0;k<atoms;++k)
        {
            
            if(master_mol_array[k]==mol_id_array_BONDS[j] && master_id_internal_array[k]==master_B2_internal_core_array[j])
            {
                master_B2_core_array[j]=master_id_array[k];
                break;
            }
        }
    }
    */
}
