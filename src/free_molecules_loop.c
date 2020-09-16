
#include"builder.h"

int global_topo_flag;
int general_bonds,angles_core,dihedrals_core,impropers_core;
int **global_bonds_array,**global_angles_array,**global_dihedrals_array,**global_impropers_array;
int atom_types,bond_types,angle_types,dihedral_types,improper_types;
char **global_species,***global_bond_type_array,***global_angle_type_array,***global_dihedral_type_array,***global_improper_type_array;
int *global_neighbors;

int molecule_bond_types,molecule_angle_types,molecule_dihedral_types,molecule_improper_types;
int **btmapping_matrix,**atmapping_matrix,**dtmapping_matrix,**itmapping_matrix;

int *added_atoms_array;
int remove_N,*rescale_array;

int *local_map;

int *init_bo_array;
int grows,**growth_registry;
int *chain_atoms_array,*chain_bonds_array;
int total_chain_atoms,total_chain_bonds;
char **species_chain_array,**B_type_chain_array;
int *B1_chain_array,*B2_chain_array;
int atoms_block_size,bonds_block_size;
char **segment_species_array;
int *segment_id_array,*segment_mol_array,*segment_host_id_array,*segment_chain_id_array,*segment_id_internal_array,*segment_mol_type_array;
double *segment_q_array,*segment_x_array,*segment_y_array,*segment_z_array;
char **segment_B_type_core_array;
int *segment_B_ID_core_array,*segment_B1_core_array,*segment_B2_core_array;
int *segment_mol_id_array_BONDS,*segment_B1_internal_core_array,*segment_B2_internal_core_array;
char **segment_atom_name_core,**segment_subst_name_core;
int *segment_B_ID_core;
double *score;
int *atom_id_new,*B1_new,*B2_new;
double *x_new,*y_new,*z_new;
int GSs;
char **current_GS_reg,**GS_reg;
int nup;
int **unique_permutations;
char **atom_name_core_B,**atom_type_core_B,**subst_name_core_B,**B_type_core_B;
int atoms_core_B,bonds_core_B;
int *atom_id_core_B,*subst_id_core_B,*B_ID_core_B,*B1_core_B,*B2_core_B;
double *x_core_B,*y_core_B,*z_core_B,*charge_core_B;
int *host_id_array_B,*chain_id_array_B;
int *one_two_core_vector;
int *one_two_list;
int *one_three_list;
int *one_four_list;
int *local_topo_array;
int *local_topo_array_B;
char **atom_name_chain_array;
int molecule_atoms,molecule_bonds;
char **molecule_species_registry, ***molecule_bond_types_registry, ***molecule_angle_types_registry, ***molecule_dihedral_types_registry, ***molecule_improper_types_registry;
int **molecule_bonds_registry, **molecule_angles_registry, **molecule_dihedrals_registry, **molecule_impropers_registry;
int molecule_angles, molecule_dihedrals, molecule_impropers;
int molecule_atom_types, molecule_bond_types, molecule_angle_types, molecule_dihedral_types, molecule_improper_types;
int molecule_atoms_B,molecule_bonds_B;
char **molecule_species_registry_B, ***molecule_bond_types_registry_B, ***molecule_angle_types_registry_B, ***molecule_dihedral_types_registry_B, ***molecule_improper_types_registry_B;
int **molecule_bonds_registry_B, **molecule_angles_registry_B, **molecule_dihedrals_registry_B, **molecule_impropers_registry_B;
int molecule_angles_B, molecule_dihedrals_B, molecule_impropers_B;
int molecule_atom_types_B, molecule_bond_types_B, molecule_angle_types_B, molecule_dihedral_types_B, molecule_improper_types_B;
int general_bonds_B, general_angles_B, general_dihedrals_B, general_impropers_B;
int general_atom_types_B, general_bond_types_B, general_angle_types_B, general_dihedral_types_B, general_improper_types_B;
char **general_species_registry_B, ***general_bond_types_registry_B, ***general_angle_types_registry_B, ***general_dihedral_types_registry_B, ***general_improper_types_registry_B;

int *general_neighbors_registry_B;

void free_molecules_loop()
{
    int j,k;
    
    if(global_topo_flag==1)
    {
        for(j=0;j<general_bonds;++j)free(global_bonds_array[j]);free(global_bonds_array);
        for(j=0;j<angles_core;++j)free(global_angles_array[j]);free(global_angles_array);
        for(j=0;j<dihedrals_core;++j)free(global_dihedrals_array[j]);free(global_dihedrals_array);
        if(impropers_core>0){for(j=0;j<impropers_core;++j)free(global_impropers_array[j]);free(global_impropers_array);}
        for(k=0;k<atom_types;++k)free(global_species[k]);free(global_species);
        for(k=0;k<bond_types;++k)
            for(j=0;j<3;++j)
                free(global_bond_type_array[k][j]);
        for(k=0;k<bond_types;++k)free(global_bond_type_array[k]);
        free(global_bond_type_array);
        for(k=0;k<angle_types;++k)
            for(j=0;j<5;++j)
                free(global_angle_type_array[k][j]);
        for(k=0;k<angle_types;++k)free(global_angle_type_array[k]);
        free(global_angle_type_array);
        for(k=0;k<dihedral_types;++k)
            for(j=0;j<7;++j)
                free(global_dihedral_type_array[k][j]);
        for(k=0;k<dihedral_types;++k)free(global_dihedral_type_array[k]);
        free(global_dihedral_type_array);
        free(global_neighbors);
        if(impropers_core>0)
        {
            for(k=0;k<improper_types;++k)
                for(j=0;j<4;++j)
                    free(global_improper_type_array[k][j]);
            for(k=0;k<improper_types;++k)free(global_improper_type_array[k]);
            free(global_improper_type_array);
        }
    }
    
    for(j=0;j<molecule_bond_types;++j)free(btmapping_matrix[j]);free(btmapping_matrix);
    for(j=0;j<molecule_angle_types;++j)free(atmapping_matrix[j]);free(atmapping_matrix);
    for(j=0;j<molecule_dihedral_types;++j)free(dtmapping_matrix[j]);free(dtmapping_matrix);
    for(j=0;j<molecule_improper_types;++j)free(itmapping_matrix[j]);free(itmapping_matrix);
    
    free(added_atoms_array);
    if(remove_N>0)free(rescale_array);
    
    free(local_map);
    
    free(init_bo_array);
    for(j=0;j<grows;++j)free(growth_registry[j]);free(growth_registry);
    free(chain_atoms_array);free(chain_bonds_array);
    for(j=0;j<total_chain_atoms;++j)free(species_chain_array[j]);free(species_chain_array);
    for(j=0;j<total_chain_bonds;++j)free(B_type_chain_array[j]);free(B_type_chain_array);
    free(B1_chain_array);free(B2_chain_array);
    for(j=0;j<atoms_block_size;++j)free(segment_species_array[j]);free(segment_species_array);
    free(segment_id_array);free(segment_mol_array);free(segment_host_id_array);free(segment_chain_id_array);free(segment_id_internal_array);free(segment_mol_type_array);
    free(segment_q_array);free(segment_x_array);free(segment_y_array);free(segment_z_array);
    for(j=0;j<bonds_block_size;++j)free(segment_B_type_core_array[j]);free(segment_B_type_core_array);
    free(segment_B_ID_core_array);free(segment_B1_core_array);free(segment_B2_core_array);
    free(segment_mol_id_array_BONDS);free(segment_B1_internal_core_array);free(segment_B2_internal_core_array);
    for(j=0;j<atoms_block_size;++j)free(segment_atom_name_core[j]);free(segment_atom_name_core);
    for(j=0;j<atoms_block_size;++j)free(segment_subst_name_core[j]);free(segment_subst_name_core);
    free(segment_B_ID_core);
    
    free(score);
    
    free(x_new);free(y_new);free(z_new);free(atom_id_new);free(B1_new);free(B2_new);
    for(j=0;j<GSs;++j)free(current_GS_reg[j]);free(current_GS_reg);
    for(j=0;j<nup;++j)free(unique_permutations[j]);free(unique_permutations);
    
    free(atom_id_core_B);free(subst_id_core_B);free(x_core_B);free(y_core_B);free(z_core_B);free(charge_core_B);
    for(j=0;j<atoms_core_B;++j){free(atom_name_core_B[j]);free(atom_type_core_B[j]);free(subst_name_core_B[j]);}
    free(host_id_array_B);free(chain_id_array_B);
    free(atom_name_core_B);free(atom_type_core_B);free(subst_name_core_B);
    free(B_ID_core_B);free(B1_core_B);free(B2_core_B);
    for(j=0;j<bonds_core_B;++j)free(B_type_core_B[j]);free(B_type_core_B);
    
    for(j=0;j<GSs;++j)free(GS_reg[j]);free(GS_reg);
    free(one_two_core_vector);
    free(one_two_list);free(one_three_list);free(one_four_list);
    free(local_topo_array);free(local_topo_array_B);
    for(j=0;j<total_chain_atoms;++j)free(atom_name_chain_array[j]);free(atom_name_chain_array);
    
    for(j=0;j<molecule_bonds;++j)free(molecule_bonds_registry[j]);free(molecule_bonds_registry);
    for(j=0;j<molecule_angles;++j)free(molecule_angles_registry[j]);free(molecule_angles_registry);
    if(molecule_dihedrals>0){for(j=0;j<molecule_dihedrals;++j)free(molecule_dihedrals_registry[j]);free(molecule_dihedrals_registry);}
    if(molecule_impropers>0){for(j=0;j<molecule_impropers;++j)free(molecule_impropers_registry[j]);free(molecule_impropers_registry);}
    for(j=0;j<molecule_atom_types;++j)free(molecule_species_registry[j]);free(molecule_species_registry);
    for(k=0;k<molecule_bond_types;++k)
        for(j=0;j<3;++j)
            free(molecule_bond_types_registry[k][j]);
    for(k=0;k<molecule_bond_types;++k)free(molecule_bond_types_registry[k]);
    free(molecule_bond_types_registry);
    for(k=0;k<molecule_angle_types;++k)
        for(j=0;j<5;++j)
            free(molecule_angle_types_registry[k][j]);
    for(k=0;k<molecule_angle_types;++k)free(molecule_angle_types_registry[k]);
    free(molecule_angle_types_registry);
    if(molecule_dihedrals>0)
    {
        for(k=0;k<molecule_dihedral_types;++k)
            for(j=0;j<7;++j)
                free(molecule_dihedral_types_registry[k][j]);
        for(k=0;k<molecule_dihedral_types;++k)free(molecule_dihedral_types_registry[k]);
        free(molecule_dihedral_types_registry);
    }
    if(molecule_impropers>0)
    {
        for(k=0;k<molecule_improper_types;++k)
            for(j=0;j<4;++j)
                free(molecule_improper_types_registry[k][j]);
        for(k=0;k<molecule_improper_types;++k)free(molecule_improper_types_registry[k]);
        free(molecule_improper_types_registry);
    }
    //-- TMB: free (in the footer of the master loop)
    for(j=0;j<molecule_bonds_B;++j)free(molecule_bonds_registry_B[j]);free(molecule_bonds_registry_B);
    for(j=0;j<molecule_angles_B;++j)free(molecule_angles_registry_B[j]);free(molecule_angles_registry_B);
    if(molecule_dihedrals_B>0){for(j=0;j<molecule_dihedrals_B;++j)free(molecule_dihedrals_registry_B[j]);free(molecule_dihedrals_registry_B);}
    if(molecule_impropers_B>0){for(j=0;j<molecule_impropers_B;++j)free(molecule_impropers_registry_B[j]);free(molecule_impropers_registry_B);}
    for(j=0;j<molecule_atom_types_B;++j)free(molecule_species_registry_B[j]);free(molecule_species_registry_B);
    for(k=0;k<molecule_bond_types_B;++k)
        for(j=0;j<3;++j)
            free(molecule_bond_types_registry_B[k][j]);
    for(k=0;k<molecule_bond_types_B;++k)free(molecule_bond_types_registry_B[k]);
    free(molecule_bond_types_registry_B);
    for(k=0;k<molecule_angle_types_B;++k)
        for(j=0;j<5;++j)
            free(molecule_angle_types_registry_B[k][j]);
    for(k=0;k<molecule_angle_types_B;++k)free(molecule_angle_types_registry_B[k]);
    free(molecule_angle_types_registry_B);
    if(molecule_dihedrals_B>0)
    {
        for(k=0;k<molecule_dihedral_types_B;++k)
            for(j=0;j<7;++j)
                free(molecule_dihedral_types_registry_B[k][j]);
        for(k=0;k<molecule_dihedral_types_B;++k)free(molecule_dihedral_types_registry_B[k]);
        free(molecule_dihedral_types_registry_B);
    }
    if(molecule_impropers_B>0)
    {
        for(k=0;k<molecule_improper_types_B;++k)
            for(j=0;j<4;++j)
                free(molecule_improper_types_registry_B[k][j]);
        for(k=0;k<molecule_improper_types_B;++k)free(molecule_improper_types_registry_B[k]);
        free(molecule_improper_types_registry_B);
    }
    //-- TGBp: free (in the footer of the master loop)
    for(j=0;j<general_atom_types_B;++j)free(general_species_registry_B[j]);free(general_species_registry_B);
    for(k=0;k<general_bond_types_B;++k)
        for(j=0;j<3;++j)
            free(general_bond_types_registry_B[k][j]);
    for(k=0;k<general_bond_types_B;++k)free(general_bond_types_registry_B[k]);
    free(general_bond_types_registry_B);
    for(k=0;k<general_angle_types_B;++k)
        for(j=0;j<5;++j)
            free(general_angle_types_registry_B[k][j]);
    for(k=0;k<general_angle_types_B;++k)free(general_angle_types_registry_B[k]);
    free(general_angle_types_registry_B);
    if(general_dihedrals_B>0)
    {
        for(k=0;k<general_dihedral_types_B;++k)
            for(j=0;j<7;++j)
                free(general_dihedral_types_registry_B[k][j]);
        for(k=0;k<general_dihedral_types_B;++k)free(general_dihedral_types_registry_B[k]);
        free(general_dihedral_types_registry_B);
    }
    if(general_impropers_B>0)
    {
        for(k=0;k<general_improper_types_B;++k)
            for(j=0;j<4;++j)
                free(general_improper_types_registry_B[k][j]);
        for(k=0;k<general_improper_types_B;++k)free(general_improper_types_registry_B[k]);
        free(general_improper_types_registry_B);
    }
    //
    free(general_neighbors_registry_B);


}
