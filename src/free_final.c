
#include"builder.h"

char **core_file;
char **quad_matrix;
int molecules,mol_types;
int total_chains;
int *type_array,*chains_array;
int *MGR_start,*MGR_stop;
int *max_growth_step_array,*growth_step_array;
int *mol_id_array,*mol_type_array,*com_flag_array;
double *X,*Y,*Z,*alpha,*beta,*my_gamma;
char **master_species_array;
int *master_id_internal_array;
int *master_mol_array;
int *master_id_array,*master_mol_type_array;
double *master_x_array,*master_y_array,*master_z_array;
int *MA_start,*MA_stop,*MB_start,*MB_stop;
char **master_B_type_core_array;
int *master_B_ID_core_array;
int *master_B1_core_array,*master_B2_core_array;
int *mol_id_array_BONDS;
int *master_B1_internal_core_array;
int *master_B2_internal_core_array;
int *host_id_array,*chain_id_array;
int master_growth_registry_rows,master_species_rows,master_topo_rows;
char **master_species_chain_array,**master_B_type_chain_array;
int **master_growth_registry;
int *master_B1_chain_array,*master_B2_chain_array;
double *master_q_array;
int entries_UFF;
char **species_UFF_array;
double *r0_UFF_array,*theta0_UFF_array,*x_UFF_array,*D_UFF_array,*zeta_UFF_array,*Zstar_UFF_array,*chi_UFF_array,*Vtor_UFF_array,*Utor_UFF_array;
int general_atoms,general_bonds;
char **general_species_registry, ***general_bond_types_registry, ***general_angle_types_registry, ***general_dihedral_types_registry, ***general_improper_types_registry;
int *general_neighbors_registry, **general_bonds_registry, **general_angles_registry, **general_dihedrals_registry, **general_impropers_registry;
int general_angles, general_dihedrals, general_impropers;
int general_atom_types, general_bond_types, general_angle_types, general_dihedral_types, general_improper_types;
int **topo_boundaries;
int *atoms_per_molecule;
int *atom_scaling_array;
char **stmap;
int **scm;
int stmap_rows;
char ***btmap;
int **bcm;
int btmap_rows;
char ***atmap;
int **acm;
int atmap_rows;
char ***dtmap;
int **dcm;
int dtmap_rows;
char ***itmap;
int **icm;
int itmap_rows;
char **master_atom_name_chain_array,**master_atom_name_array;
int backbones;
char **backbone_id;
double *bb_dx,**bb_entries;
int bb_lines_max;
int tacticity;
char **tacticity_id,**tacticity_type;
char **tacticity_search;
int *tacticity_counter_array;
double *tacticity_array;
int *yet_one_more_t_array;
int *t_tracker_1,*t_tracker_2,*t_tracker_3,*t_tracker_4;
int *t_tracker_mol;
int *tacticity_counter_array_v2;
//int *Re_start,*Re_stop,*bb_growth_step_array;
//int *prev_link_atom_array;
//int *added_backbone_flag;
//double *l_array,*l2_array;

int backbone_type_counter;
char **backbone_type_array;

int *head,**neighbors;

double *bb_pdf_dx,**bb_pdf_entries;
int bb_pdf_lines_max;

int CBMCG_flag;
double *pe_array;

int entries_DRE;
char **species_DRE_array;
double *R0_DRE_array,*theta0_DRE_array,*Rvdw0_DRE_array,*D0_DRE_array;

void free_final()
{
    int i,j,k;
        
    for(j=0;j<general_bonds;++j)free(general_bonds_registry[j]);free(general_bonds_registry);
    for(j=0;j<general_angles;++j)free(general_angles_registry[j]);free(general_angles_registry);
    for(j=0;j<general_dihedrals;++j)free(general_dihedrals_registry[j]);free(general_dihedrals_registry);
    if(general_impropers>0){for(j=0;j<general_impropers;++j)free(general_impropers_registry[j]);free(general_impropers_registry);}
    for(k=0;k<general_atom_types;++k)free(general_species_registry[k]);free(general_species_registry);
    for(k=0;k<general_bond_types;++k)
        for(j=0;j<3;++j)
            free(general_bond_types_registry[k][j]);
    for(k=0;k<general_bond_types;++k)free(general_bond_types_registry[k]);
    free(general_bond_types_registry);
    for(k=0;k<general_angle_types;++k)
        for(j=0;j<5;++j)
            free(general_angle_types_registry[k][j]);
    for(k=0;k<general_angle_types;++k)free(general_angle_types_registry[k]);
    free(general_angle_types_registry);
    for(k=0;k<general_dihedral_types;++k)
        for(j=0;j<7;++j)
            free(general_dihedral_types_registry[k][j]);
    for(k=0;k<general_dihedral_types;++k)free(general_dihedral_types_registry[k]);
    free(general_dihedral_types_registry);
    free(general_neighbors_registry);
    if(general_impropers>0)
    {
        for(k=0;k<general_improper_types;++k)
            for(j=0;j<4;++j)
                free(general_improper_types_registry[k][j]);
        for(k=0;k<general_improper_types;++k)free(general_improper_types_registry[k]);
        free(general_improper_types_registry);
    }
    for(j=0;j<stmap_rows;++j)free(stmap[j]);free(stmap);
    for(j=0;j<molecules;++j)free(scm[j]);free(scm);
    
    for(k=0;k<btmap_rows;++k)
        for(j=0;j<3;++j)
            free(btmap[k][j]);
    for(k=0;k<btmap_rows;++k)free(btmap[k]);
    free(btmap);
    for(j=0;j<molecules;++j)free(bcm[j]);free(bcm);
    for(k=0;k<atmap_rows;++k)
        for(j=0;j<5;++j)
            free(atmap[k][j]);
    for(k=0;k<atmap_rows;++k)free(atmap[k]);
    free(atmap);
    for(j=0;j<molecules;++j)free(acm[j]);free(acm);
    for(k=0;k<dtmap_rows;++k)
        for(j=0;j<7;++j)
            free(dtmap[k][j]);
    for(k=0;k<dtmap_rows;++k)free(dtmap[k]);
    free(dtmap);
    for(j=0;j<molecules;++j)free(dcm[j]);free(dcm);
    for(k=0;k<itmap_rows;++k)
        for(j=0;j<4;++j)
            free(itmap[k][j]);
    for(k=0;k<itmap_rows;++k)free(itmap[k]);
    free(itmap);
    for(j=0;j<molecules;++j)free(icm[j]);free(icm);
    
    for(i=0;i<master_growth_registry_rows;++i)free(master_growth_registry[i]);free(master_growth_registry);
    for(i=0;i<mol_types;++i)free(core_file[i]);free(core_file);
    free(type_array);free(chains_array);
    for(i=0;i<total_chains;++i)free(quad_matrix[i]);free(quad_matrix);
    free(MGR_start);free(MGR_stop);free(max_growth_step_array);free(growth_step_array);
    free(mol_id_array);free(mol_type_array);free(com_flag_array);
    free(X);free(Y);free(Z);free(alpha);free(beta);free(my_gamma);
    free(master_id_internal_array);free(master_mol_array);free(master_x_array);free(master_y_array);free(master_z_array);
    for(i=0;i<general_atoms;++i)free(master_species_array[i]);free(master_species_array);
    for(i=0;i<general_atoms;++i)free(master_atom_name_array[i]);free(master_atom_name_array);
    for(i=0;i<master_species_rows;++i)free(master_species_chain_array[i]);free(master_species_chain_array);
    for(i=0;i<master_topo_rows;++i)free(master_B_type_chain_array[i]);free(master_B_type_chain_array);
    free(master_B1_chain_array);free(master_B2_chain_array);
    free(MA_start);free(MA_stop);free(MB_start);free(MB_stop);
    for(i=0;i<general_bonds;++i)free(master_B_type_core_array[i]);free(master_B_type_core_array);
    free(master_B1_core_array);free(master_B2_core_array);free(mol_id_array_BONDS);
    free(master_id_array);free(master_mol_type_array);
    free(master_B1_internal_core_array);free(master_B2_internal_core_array);
    free(host_id_array);free(chain_id_array);
    free(master_q_array);free(master_B_ID_core_array);
    for(i=0;i<entries_UFF;++i)free(species_UFF_array[i]);free(species_UFF_array);
    free(r0_UFF_array);free(theta0_UFF_array);free(x_UFF_array);free(D_UFF_array);
    free(zeta_UFF_array);free(Zstar_UFF_array);free(chi_UFF_array);free(Vtor_UFF_array);free(Utor_UFF_array);
    for(i=0;i<molecules;++i)free(topo_boundaries[i]);free(topo_boundaries);
    free(atoms_per_molecule);free(atom_scaling_array);
    for(i=0;i<master_species_rows;++i)free(master_atom_name_chain_array[i]);free(master_atom_name_chain_array);
    if(backbones!=0){for(i=0;i<backbones;++i)free(backbone_id[i]);free(backbone_id);}
    if(backbones>0){
        free(bb_dx);
        for(i=0;i<bb_lines_max;++i)free(bb_entries[i]);free(bb_entries);
        free(bb_pdf_dx);
        for(i=0;i<bb_pdf_lines_max;++i)free(bb_pdf_entries[i]);free(bb_pdf_entries);
    }
    if(tacticity>0)
    {
        for(i=0;i<tacticity;++i)free(tacticity_id[i]);free(tacticity_id);
        for(i=0;i<tacticity;++i)free(tacticity_type[i]);free(tacticity_type);
        for(i=0;i<mol_types;++i)free(tacticity_search[i]);free(tacticity_search);
        free(tacticity_counter_array);free(yet_one_more_t_array);free(tacticity_array);
        free(tacticity_counter_array_v2);
        free(t_tracker_1);
        free(t_tracker_2);
        free(t_tracker_3);
        free(t_tracker_4);
        
        free(t_tracker_mol);
    }
    if(backbones>0)
    {
        //free(Re_start);free(Re_stop);free(prev_link_atom_array);
        //free(bb_growth_step_array);free(added_backbone_flag);
        //free(l_array);free(l2_array);
        for(i=0;i<backbone_type_counter;++i)free(backbone_type_array[i]);free(backbone_type_array);
    }
    
    if(CBMCG_flag==1){
    free(head);
    for(i=0;i<26;++i)free(neighbors[i]);free(neighbors);
    }
    free(pe_array);
    
    for(i=0;i<entries_DRE;++i)free(species_DRE_array[i]);free(species_DRE_array);
    free(R0_DRE_array);free(theta0_DRE_array);free(Rvdw0_DRE_array);free(D0_DRE_array);
    
}
