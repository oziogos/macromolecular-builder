
void growth_reg_init(int *type, char **core_mol2, int *chains, int mol_types, int total_chains, char **quad_matrix, int *rows_out, int *cols_out, int *species_rows, int *topo_rows,
                     int *master_growth_registry_rows,int *master_growth_registry_cols,int *master_species_rows,int *master_topo_rows,
                     char ***master_species_chain_array,char ***master_B_type_chain_array,int ***master_growth_registry,
                     int **master_B1_chain_array,int **master_B2_chain_array,
                     char ***master_atom_name_chain_array,
                     int *MGR_start, int *MGR_stop, int *MA_start, int *MA_stop, int *MB_start, int *MB_stop,
                     int molecules, int *mol_type_array, int **max_growth_step_array, int **growth_step_array);

void builder2_initialize(char *current_folder,char *argv1,int *molecules,int *mol_types,double *xlo,double *xhi,double *ylo,double *yhi,double *zlo,double *zhi,
                         char ***core_file,int **type_array,int **chains_array,int **mol_id_array,int **mol_type_array,int **com_flag_array,
                         double **X,double **Y,double **Z,double **alpha,double **beta,double **gamma,
                         int *total_chains,char ***quad_matrix,
                         int *master_growth_registry_rows,int *master_growth_registry_cols,int *master_species_rows,int *master_topo_rows,
                         int **MGR_start,int **MGR_stop,int **MA_start,int **MA_stop,int **MB_start,int **MB_stop,
                         int *entries_UFF,
                         char ***species_UFF_array,
                         double **D_UFF_array,double **x_UFF_array,double **r0_UFF_array,double **chi_UFF_array,double **Zstar_UFF_array,
                         double **theta0_UFF_array,double **Vtor_UFF_array,double **Utor_UFF_array,double **zeta_UFF_array,
                         int *entries_DRE,
                         char ***species_DRE_array,
                         double **R0_DRE_array,double **theta0_DRE_array,double **Rvdw0_DRE_array,double **D0_DRE_array,
                         int *backbones,char ***backbone_id,double **bb_dx,double ***bb_entries,int *bb_lines_max,
                         int *tacticity,char ***tacticity_id,char ***tacticity_type,
                         int *backbone_type_counter, char ***backbone_type_array,
                         char *mpi_from_config, char *lmp_from_config,
                         double **bb_pdf_dx,double ***bb_pdf_entries,int *bb_pdf_lines_max);

void builder2_euler(char *current_folder,char *argv1,int atoms,int molecules,
                    int *com_flag_array,int *master_mol_array,int *mol_id_array,int *master_id_internal_array,char **master_species_array,
                    double *X,double *Y,double *Z,double *alpha,double *beta,double *gamma,
                    double *master_x_array,double *master_y_array,double *master_z_array,int reshuffle_euler,int reshuffle_coms,double rand_disp,
                    double xlo,double ylo,double zlo,double xhi,double yhi,double zhi);

void contiguous_atom(int *atoms,int *bonds,int molecules,int mol_types,int *mol_id_array,int *mol_type_array,int *type_array,char **core_file,
                     char *current_folder,int **master_mol_array,int **master_id_internal_array,double **master_x_array,double **master_y_array,double **master_z_array,
                     double **master_q_array,
                     char ***master_species_array,int **master_id_array,int **master_mol_type_array,int total_chains,char **quad_matrix,int **host_id_array,int **chain_id_array,
                     int **atoms_per_molecule,int **atom_scaling_array,
                     char ***master_atom_name_array);

void contiguous_bonds(int atoms,int bonds,int molecules,int mol_types,int *mol_type_array,int *type_array,char *current_folder,char **core_file,
                      int *master_id_array,int *master_id_internal_array,int *master_mol_array,
                      char ***master_B_type_core_array,int **master_B_ID_core_array,int **master_B1_core_array,int **master_B2_core_array,int **mol_id_array_BONDS,
                      int **master_B1_internal_core_array,int **master_B2_internal_core_array);

void retrieve_mol_type_data(int current_mol_type,int mol_types,int *type_array,int *MGR_start,int *MGR_stop,int *MA_start,int *MA_stop,int *MB_start,int *MB_stop,
                            int total_chains,char **quad_matrix,int *chains_array,
                            char **master_species_chain_array,int *master_B1_chain_array,int *master_B2_chain_array,char **master_B_type_chain_array,
                            int gcols,int **master_growth_registry,
                            int *grows,int **init_bo_array,int ***growth_registry,
                            int *total_chain_atoms,int *total_chain_bonds,int **chain_atoms_array,int **chain_bonds_array,char ***species_chain_array,
                            int **B1_chain_array,int **B2_chain_array,char ***B_type_chain_array,
                            char **master_atom_name_chain_array,char ***atom_name_chain_array);

void topo(int atoms,char **species,int bonds,int *init_B1,int *init_B2,char **init_B_type,
          char ***global_species,char ****global_bond_type_array,char ****global_angle_type_array,char ****global_dihedral_type_array,char ****global_improper_type_array,
          int **global_neighbors,int ***global_bonds_array,int ***global_angles_array,int ***global_dihedrals_array,int ***global_impropers_array,
          int *angles_core,int *dihedrals_core,int *impropers_core,
          int *atom_types,int *bond_types,int *angle_types,int *dihedral_types,int *improper_types,
          int output_flag,char *current_folder,char *file_name);

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
              );

void builder_loop_preamble(int growth_step,int *chain_index,int **growth_registry,int atoms_core,int *host_id_array,int *chain_id_array,int *host_index,
                           int *max_1_2_core,int *B1_core,int *B2_core,int bonds_core,int *GS_array,char **atom_type_core,int *GSs,
                           int **one_two_core_vector,char ***GS_reg,
                           int *chain_atoms_array,int *init_bo_array,int *chain_bonds_array,int *B1_chain_array,int *B2_chain_array,char **B_type_chain_array,
                           char **species_chain_array,int *nup,int ***unique_permutations,
                           int *new_entries_atom,int *new_entries_bond,int **atom_id_new,double **x_new,double **y_new,double **z_new,
                           int **B1_new,int **B2_new,
                           int *atoms_core_B,int *bonds_core_B,int **atom_id_core_B,int **subst_id_core_B,double **x_core_B,double **y_core_B,double **z_core_B,double **charge_core_B,
                           char ***atom_name_core_B,char ***atom_type_core_B,char ***subst_name_core_B,int **host_id_array_B,int **chain_id_array_B,
                           int **B_ID_core_B,int **B1_core_B,int **B2_core_B,char ***B_type_core_B,
                           int *atom_id_core,int *subst_id_core,double *x_core,double *y_core,double *z_core,double *charge_core,
                           char **atom_name_core,char **subst_name_core,
                           int *B_ID_core,char **B_type_core,
                           char ***current_GS_reg,double **score,
                           int *new_entries_atom_total,int *new_entries_bond_total,
                           int *one_two_counter,int **one_two_list,int *one_three_counter,int **one_three_list,int *one_four_counter,int **one_four_list,
                           char **atom_name_chain_array,
                           int **local_topo_array,int **local_topo_array_B);

void builder_core_funct(int chain_index,int nup_index,int nup,
                        int GSs,int host_index,int max_1_2_core,
                        int new_entries_atom,int new_entries_bond,
                        int *atoms_core,int *bonds_core,
                        int atoms_core_B,int bonds_core_B,
                        int *one_two_core_vector,int *GS_array,
                        char **GS_reg,char **current_GS_reg,int **unique_permutations,
                        int **atom_id_core,int **subst_id_core,int **B_ID_core,int **B1_core,int **B2_core,int **host_id_array,int **chain_id_array,
                        double **x_core,double **y_core,double **z_core,double **charge_core,
                        char ***atom_name_core,char ***atom_type_core,char ***subst_name_core,char ***B_type_core,
                        int *atom_id_core_B,int *subst_id_core_B,int *B_ID_core_B,int *B1_core_B,int *B2_core_B,int *host_id_array_B,int *chain_id_array_B,
                        double *x_core_B,double *y_core_B,double *z_core_B,double *charge_core_B,
                        char **atom_name_core_B,char **atom_type_core_B,char **subst_name_core_B,char **B_type_core_B,
                        int *atom_id_new,double *x_new,double *y_new,double *z_new,int *B1_new,int *B2_new,
                        int **local_topo_array,int *local_topo_array_B,
                        int *remove_N,int *added_atoms,int **added_atoms_array,
                        int **rescale_array,int perc_flag
                        );

void resize_master_registries(int *atoms,int new_entries_atom_total,int atoms_last_entry,
                              int **master_id_array,int **master_mol_array,int **host_id_array,int **chain_id_array,int **master_id_internal_array,int **master_mol_type_array,
                              double **master_q_array,double **master_x_array,double **master_y_array,double **master_z_array,
                              char ***master_species_array,
                              int *bonds,int new_entries_bond_total,int bonds_last_entry,
                              int **master_B_ID_core_array,int **master_B1_core_array,int **master_B2_core_array,int **mol_id_array_BONDS,
                              int **master_B1_internal_core_array,int **master_B2_internal_core_array,
                              char ***master_B_type_core_array,
                              char ***master_atom_name_array);

void write_xyz_preview(char *current_folder,int atoms,char **master_species_array,double *master_x_array,double *master_y_array,double *master_z_array,
                       double xlo,double xhi,double ylo,double yhi,double zlo,double zhi,int append_flag);

void topo_molecule_alloc_populate(int i,int atoms_block_size,int **topo_boundaries,
                                  int general_bond_types,int general_angle_types,int general_dihedral_types,int general_improper_types,
                                  int **general_bonds_registry,int **general_angles_registry,int **general_dihedrals_registry,int **general_impropers_registry,
                                  char ***general_bond_types_registry,char ***general_angle_types_registry,char ***general_dihedral_types_registry,char ***general_improper_types_registry,
                                  int *molecule_atoms,int *molecule_bonds,int *molecule_angles,int *molecule_dihedrals,int *molecule_impropers,
                                  int *molecule_atom_types,int *molecule_bond_types,int *molecule_angle_types,int *molecule_dihedral_types,int *molecule_improper_types,
                                  int ***molecule_bonds_registry,int ***molecule_angles_registry,int ***molecule_dihedrals_registry,int ***molecule_impropers_registry,
                                  char ***molecule_species_registry,
                                  char ****molecule_bond_types_registry,char ****molecule_angle_types_registry,char ****molecule_dihedral_types_registry,char ****molecule_improper_types_registry,
                                  int *molecule_atoms_B,int *molecule_bonds_B,int *molecule_angles_B,int *molecule_dihedrals_B,int *molecule_impropers_B,
                                  int *molecule_atom_types_B,int *molecule_bond_types_B,int *molecule_angle_types_B,int *molecule_dihedral_types_B,int *molecule_improper_types_B,
                                  int ***molecule_bonds_registry_B,int ***molecule_angles_registry_B,int ***molecule_dihedrals_registry_B,int ***molecule_impropers_registry_B,
                                  char ***molecule_species_registry_B,
                                  char ****molecule_bond_types_registry_B,char ****molecule_angle_types_registry_B,char ****molecule_dihedral_types_registry_B,char ****molecule_improper_types_registry_B,
                                  
                                  int verb
                                  );

void topo_general_backup(int general_bonds,int general_angles,int general_dihedrals,int general_impropers,
                         int general_atom_types,int general_bond_types,int general_angle_types,int general_dihedral_types,int general_improper_types,
                         char **general_species_registry,
                         char ***general_bond_types_registry,char ***general_angle_types_registry,char ***general_dihedral_types_registry,char ***general_improper_types_registry,
                         int *general_bonds_B,int *general_angles_B,int *general_dihedrals_B,int *general_impropers_B,
                         int *general_atom_types_B,int *general_bond_types_B,int *general_angle_types_B,int *general_dihedral_types_B,int *general_improper_types_B,
                         char ***general_species_registry_B,
                         char ****general_bond_types_registry_B,char ****general_angle_types_registry_B,char ****general_dihedral_types_registry_B,char ****general_improper_types_registry_B,
                         int general_atoms,int *general_neighbors_registry,int **general_neighbors_registry_B);

void alter_molecular_topo(int local_bonds,int local_angles,int local_dihedrals,int local_impropers,
                          int local_atom_types,
                          
                          int *local_map,
                          
                          char **local_species_registry,
                          int **local_bonds_registry,int **local_angles_registry,int **local_dihedrals_registry,int **local_impropers_registry,
                          char ***local_bond_types_registry,char ***local_angle_types_registry,char ***local_dihedral_types_registry,char ***local_improper_types_registry,
                          
                          // output
                          
                          int *molecule_bonds,int *molecule_angles,int *molecule_dihedrals,int *molecule_impropers,
                          int *molecule_atom_types,int *molecule_bond_types,int *molecule_angle_types,int *molecule_dihedral_types,int *molecule_improper_types,
                          
                          char ***molecule_species_registry,
                          int ***molecule_bonds_registry,int ***molecule_angles_registry,int ***molecule_dihedrals_registry,int ***molecule_impropers_registry,
                          char ****molecule_bond_types_registry,char ****molecule_angle_types_registry,char ****molecule_dihedral_types_registry,char ****molecule_improper_types_registry,
                          
                          
                          int verb
                          );

void compute_local(char *current_folder,double *segment_x_array,double *segment_y_array,double *segment_z_array,
                   int atoms_block_size,int bonds_block_size,int *local_topo_array,char **segment_species_array,
                   int *segment_B1_internal_core_array,int *segment_B2_internal_core_array,char **segment_B_type_core_array,
                   
                   int *local_atoms,int *local_bonds,char ***local_species,int **local_B1,int **local_B2,char ***local_B_type,int **local_map
                   );

void alter_general_topo(int i, int molecules,
                        int molecule_atoms, int molecule_atoms_B, int delta_atoms, int molecule_atom_types,
                        int molecule_bonds, int molecule_bonds_B, int molecule_bond_types,
                        int molecule_angles, int molecule_angles_B, int molecule_angle_types,
                        int molecule_dihedrals, int molecule_dihedrals_B, int molecule_dihedral_types,
                        int molecule_impropers, int molecule_impropers_B, int molecule_improper_types,
                        int *atoms_per_molecule, int *atom_scaling_array,
                        char **molecule_species_registry,
                        char ***molecule_bond_types_registry, int **molecule_bonds_registry,
                        char ***molecule_angle_types_registry, int **molecule_angles_registry,
                        char ***molecule_dihedral_types_registry, int **molecule_dihedrals_registry,
                        char ***molecule_improper_types_registry, int **molecule_impropers_registry,
                        int *general_bonds, int *general_angles, int *general_dihedrals, int *general_impropers,
                        int **topo_boundaries,
                        //
                        int *stmap_rows, char ***stmap, int ***scm,
                        int *general_atom_types,
                        char ***general_species_registry,
                        
                        int *btmap_rows, char ****btmap, int ***bcm,
                        int *general_bond_types, int ***general_bonds_registry,
                        char ****general_bond_types_registry,
                        int ***btmapping_matrix,
                        
                        int *atmap_rows, char ****atmap, int ***acm,
                        int *general_angle_types, int ***general_angles_registry,
                        char ****general_angle_types_registry,
                        int ***atmapping_matrix,
                        
                        int *dtmap_rows, char ****dtmap, int ***dcm,
                        int *general_dihedral_types, int ***general_dihedrals_registry,
                        char ****general_dihedral_types_registry,
                        int ***dtmapping_matrix,
                        
                        int *itmap_rows, char ****itmap, int ***icm,
                        int *general_improper_types, int ***general_impropers_registry,
                        char ****general_improper_types_registry,
                        int ***itmapping_matrix,
                        
                        //
                        
                        int verb
                        
                        );

void tmap_topo_boundaries_init(int molecules, int *master_mol_array, char **master_species_array,
                               int atoms, int general_atom_types, char **general_species_registry,
                               int general_bonds, int general_bond_types, char ***general_bond_types_registry, int **general_bonds_registry,
                               int general_angles, int general_angle_types, char ***general_angle_types_registry, int **general_angles_registry,
                               int general_dihedrals, int general_dihedral_types, char ***general_dihedral_types_registry, int **general_dihedrals_registry,
                               int general_impropers, int general_improper_types, char ***general_improper_types_registry, int **general_impropers_registry,
                               //
                               int *stmap_rows, char ***stmap, int ***scm,
                               int *btmap_rows, char ****btmap, int ***bcm, int *atmap_rows, char ****atmap, int ***acm,
                               int *dtmap_rows, char ****dtmap, int ***dcm, int *itmap_rows, char ****itmap, int ***icm,
                               int ***topo_boundaries,
                               int verb);

void update_master(int i, int atoms, int bonds, int atoms_block_size, int bonds_block_size, int atoms_first_entry, int bonds_first_entry, int *master_id_array,
                   int *master_mol_array, int *host_id_array, int *chain_id_array, int *master_id_internal_array, int *master_mol_type_array,
                   double *master_q_array, double *master_x_array, double *master_y_array, double *master_z_array, char **master_species_array,
                   char **master_B_type_core_array, int *mol_id_array_BONDS, int *master_B1_core_array, int *master_B2_core_array, int *master_B1_internal_core_array, int *master_B2_internal_core_array,
                   int *segment_mol_array, int *segment_host_id_array, int *segment_chain_id_array, int *segment_id_internal_array, int *mol_type_array,
                   double *segment_q_array, double *segment_x_array, double *segment_y_array, double *segment_z_array, char **segment_species_array,
                   char **segment_B_type_core_array, int *segment_B1_internal_core_array, int *segment_B2_internal_core_array,
                   int *atom_scaling_array,
                   char **master_atom_name_array,char **segment_atom_name_core);

void topo_diagnostic(int general_bonds, int **general_bonds_registry, char ***general_bond_types_registry,
                     int general_angles, int **general_angles_registry, char ***general_angle_types_registry,
                     int general_dihedrals, int **general_dihedrals_registry, char ***general_dihedral_types_registry,
                     int general_impropers, int **general_impropers_registry, char ***general_improper_types_registry,
                     int bonds, int **global_bonds_array, char ***global_bond_type_array,
                     int angles_core, int **global_angles_array, char ***global_angle_type_array,
                     int dihedrals_core, int **global_dihedrals_array, char ***global_dihedral_type_array,
                     int impropers_core, int **global_impropers_array, char ***global_improper_type_array);

void rotate(double ux, double uy, double uz, double x_center, double y_center, double z_center, double theta, double x_old, double y_old, double z_old, double *x_new, double *y_new, double *z_new);

void remove_deleted(int molecule_atoms,int *rescale_array,
                    char **segment_species_array,
                    int atoms_block_size,int bonds_block_size,
                    int *segment_B1_internal_core_array,int *segment_B2_internal_core_array,char **segment_B_type_core_array,
                    
                    int *molecule_atom_types,
                    int *molecule_bond_types,int *molecule_angle_types,int *molecule_dihedral_types,int *molecule_improper_types,
                    int *molecule_bonds,int *molecule_angles,int *molecule_dihedrals,int molecule_impropers,
                    char ***molecule_species_registry,
                    char ****molecule_bond_types_registry,char ****molecule_angle_types_registry,char ****molecule_dihedral_types_registry,char ****molecule_improper_types_registry,
                    int ***molecule_bonds_registry,int ***molecule_angles_registry,int ***molecule_dihedrals_registry,int **molecule_impropers_registry
                    );

void initialize_tacticity(int master_species_rows, char **master_atom_name_chain_array, int mol_types, int *MA_start, int *MA_stop, int tacticity, int molecules, char **tacticity_id,
                          int *mol_type_array, char **tacticity_type,
                          
                          int **yet_one_more_t_array, int **tacticity_counter_array, int **tacticity_counter_array_v2, char ***tacticity_search, int *tacticity_array_length,
                          double **tacticity_array,
                          int **t_tracker_1,int **t_tracker_2,int **t_tracker_3,int **t_tracker_4,int **t_tracker_mol
                          );

void apply_backbone_dihedral(int i, int backbones, int local_atoms,int *local_topo_array, int *local_map, char **segment_atom_name_core, char **backbone_id, int *atom_scaling_array,
                             double **bb_entries,double *bb_dx,
                             double *master_x_array,double *master_y_array,double *master_z_array,
                             //int *Re_stop,int *bb_growth_step_array, int *prev_link_atom_array,int *added_backbone_flag,
                             int backbone_type_counter, char **backbone_type_array
                             );

void force_tacticity(int i, int tacticity, int *mol_type_array, char **tacticity_search, int local_atoms, char **segment_atom_name_core, int *local_map, int *local_topo_array,
                     int *atom_scaling_array, char **backbone_id, int *yet_one_more_t_array, int backbones,
                     int backbone_type_counter, char **backbone_type_array,
                     
                     int molecules,int tacticity_array_length,
                     
                     int *tacticity_counter_array,
                     int *general_t_counter, int *t_tracker_1, int *t_tracker_2, int *t_tracker_3, int *t_tracker_4,int *t_tracker_mol,
                     double *tacticity_array, double *master_x_array, double *master_y_array, double *master_z_array,
                     int *ft);
void free_final();
void free_locals();
void free_molecules_loop();

void minimize(char *current_folder,int counter,int l,
              double xlo,double xhi,double ylo,double yhi,double zlo,double zhi,
              int *master_id_array,int *master_mol_array,double *master_q_array,double *master_x_array,double *master_y_array,double *master_z_array,
              int general_atoms, int general_bonds, int general_angles, int general_dihedrals, int general_impropers,
              int general_atom_types, int general_bond_types, int general_angle_types, int general_dihedral_types, int general_improper_types,
              char **master_species_array,
              char **general_species_registry,
              int **general_bonds_registry, int **general_angles_registry, int **general_dihedrals_registry, int **general_impropers_registry,
              char ***general_bond_types_registry,char ***general_angle_types_registry,char ***general_dihedral_types_registry,char ***general_improper_types_registry,
              int entries_UFF, char **species_UFF_array,
              double *D_UFF_array,double *x_UFF_array,double *r0_UFF_array,double *chi_UFF_array,double *Zstar_UFF_array,
              double *theta0_UFF_array,double *Vtor_UFF_array,double *Utor_UFF_array,double *zeta_UFF_array,
              int entries_DRE,
              char **species_DRE_array,
              double *R0_DRE_array,double *theta0_DRE_array,double *Rvdw0_DRE_array,double *D0_DRE_array,
              int *general_neighbors_registry,
              int nproc,
              char *mpi_from_config,char *lmp_from_config,
              int FF);

void write_mol2_file(char *current_folder,int atoms,int bonds,
                     int *master_id_array,char **master_atom_name_array,
                     double *master_x_array,double *master_y_array,double *master_z_array,
                     char **master_species_array,int *master_mol_array,
                     int *master_B1_core_array,int *master_B2_core_array,char **master_B_type_core_array,char *mol2_file_name);

void calculate_tacticity(int molecules, int general_t_counter, int tacticity_array_length,
                         int *atom_scaling_array, double *master_x_array, double *master_y_array, double *master_z_array,
                         int *tacticity_counter_array, int *tacticity_counter_array_v2, int *yet_one_more_t_array,
                         int *t_tracker_1, int *t_tracker_2, int *t_tracker_3, int *t_tracker_4, int *t_tracker_mol,
                         double *tacticity_array,int counter,char *current_folder);

void update_neighbors_registry(int i,int general_atoms,int delta_atoms,int molecules,int molecule_atoms,int atoms_first_entry,int atoms_last_entry,
                               char **master_species_array,int **general_neighbors_registry);
/*
void initialize_Cn(int molecules, int *max_growth_step_array, char **quad_matrix, int *mol_type_array,
                   int **bb_growth_step_array,int **added_backbone_flag, double **l_array, double **l2_array, int **Re_start, int **Re_stop, int **prev_link_atom_array);

void calculate_Cn(char *current_folder, int counter, int molecules, int general_bonds, int backbones,
                  int *max_growth_step_array, int *atom_scaling_array, int *master_B1_core_array, int *master_B2_core_array, int *mol_id_array_BONDS,
                  double *master_x_array, double *master_y_array, double *master_z_array,char **master_atom_name_array, char **backbone_id,
                  double *l_array, double *l2_array, int *Re_start, int *Re_stop, int *bb_growth_step_array,
                  int backbone_type_counter, char **backbone_type_array);
*/

void post_exec(int myseed);

void topo_neighbors(int general_atoms, int general_bonds, int *master_B1_core_array, int *master_B2_core_array,
                    
                    int ***one_two, int ***one_three, int ***one_four);

void linked_list_cell_initialize(int cutoff,double xlo,double xhi,double ylo,double yhi,double zlo,double zhi,int general_atoms,
                                 int *Mx,int *My,int *Mz,int *M3,double *lx,double *ly,double *lz,double *lcx,double *lcy,double *lcz,int ***neighbors,int **head);

double lj_pe_brute_force(int general_atoms, double lx, double ly, double lz, double cutoff, int **one_two, int **one_three, int **one_four,double **nb_FF,int *nb_type_array,
                         double *master_x_array, double *master_y_array,double *master_z_array,int *master_nx_array, int *master_ny_array,int *master_nz_array,
                         int *pairs_bf,int *host_id_array);

double lj_pe_linked_list_cell(int general_atoms, double lx, double ly, double lz, double cutoff, int **one_two, int **one_three, int **one_four,double **nb_FF,int *nb_type_array,
                              double *master_x_array, double *master_y_array,double *master_z_array,int *master_nx_array, int *master_ny_array,int *master_nz_array,
                              
                              int M3, int *head, int *list, int **neighbors,
                              
                              int *pairs_llc,int *host_id_array
                              );

void calc_head_list(int general_atoms, int Mx, int M2, int M3,
                    double xlo, double ylo, double zlo, double xhi, double yhi, double zhi, double lx, double ly, double lz, double lcx, double lcy, double lcz,
                    double *master_x_array,double *master_y_array, double *master_z_array, int *master_nx_array, int *master_ny_array, int *master_nz_array,
                    int *head, int **list);

void UFF_lj_params(int general_atoms, int general_atom_types,char **master_species_array, char **general_species_registry,
              int entries_UFF, char **species_UFF_array, double *D_UFF_array, double *x_UFF_array,
              double ***nb_FF, int **nb_type_array);

void calc_wrapped_coords_n_arrays(int general_atoms,double xlo,double ylo,double zlo,double xhi,double yhi,double zhi,double lx,double ly,double lz,
                                  double *master_x_array,double *master_y_array,double *master_z_array,int *master_nx_array,int *master_ny_array,int *master_nz_array);

int equal(double A,double B);
int equal(double A,double B)
{
    int value;
    if(fabs(A-B)<equal_diff){value=1;}else{value=0;}
    return value;
}

int one_two_init(int *B1, int *B2, int atoms, int bonds);
int one_two_init(int *B1, int *B2, int atoms, int bonds)
{
    int i,max_1_2,*one_two_counter_array;
    one_two_counter_array=(int*)malloc(atoms*sizeof(int));
    for(i=0;i<atoms;++i)one_two_counter_array[i]=0;
    for(i=0;i<bonds;++i)
    {
        one_two_counter_array[B1[i]-1]=one_two_counter_array[B1[i]-1]+1;
        one_two_counter_array[B2[i]-1]=one_two_counter_array[B2[i]-1]+1;
    }
    max_1_2=1;
    for(i=0;i<atoms;++i)
        if(one_two_counter_array[i]>max_1_2)max_1_2=one_two_counter_array[i];
    free(one_two_counter_array);
    return max_1_2;
}

void one_two_build(int **one_two, int *B1, int *B2, int atoms, int bonds, int max_1_2);
void one_two_build(int **one_two, int *B1, int *B2, int atoms, int bonds, int max_1_2)
{
    int i,j;
    //int c=0;
    for(i=0;i<atoms;++i)
        for(j=0;j<max_1_2;++j)
            one_two[i][j]=0;
    for(i=0;i<bonds;++i)
    {
        j=0;
        while(one_two[B1[i]-1][j]!=0)
        {
            //c=c+1;
            j=j+1;
        }
        one_two[B1[i]-1][j]=B2[i];
        j=0;
        while(one_two[B2[i]-1][j]!=0)
        {
            //c=c+1;
            j=j+1;
        }
        one_two[B2[i]-1][j]=B1[i];
    }
    //printf("@@##\t%d\t%d\n",bonds,c);//getchar();
}

//void rotate(double ux, double uy, double uz, double x_center, double y_center, double z_center, double theta, double x_old, double y_old, double z_old, double *x_new, double *y_new, double *z_new);

void CBMCG(int i,int local_atoms,int backbone_type_counter,int backbones,int CBMCG_N_trial,char *current_folder,
           int *local_topo_array,int *local_map,char **segment_atom_name_core,char **backbone_type_array,int *atom_scaling_array,char **backbone_id,
           double *bb_dx,double **bb_entries,double *bb_pdf_dx,double **bb_pdf_entries,
           
           int tacticity,int *mol_type_array,char **tacticity_search,int *yet_one_more_t_array,
           int molecules,int tacticity_array_length,
           int *tacticity_counter_array,
           int *general_t_counter, int *t_tracker_1, int *t_tracker_2, int *t_tracker_3, int *t_tracker_4,int *t_tracker_mol,
           double *tacticity_array,
           
           char **master_species_array,double *master_x_array,double *master_y_array,double *master_z_array,
           
           double cutoff,double xlo,double ylo,double zlo,double xhi,double yhi,double zhi,int general_atoms,int general_bonds,int general_atom_types,
           int Mx,int My,int Mz,int M3,double lx,double ly,double lz,double lcx,double lcy,double lcz,int **neighbors,
           int *master_B1_core_array,int *master_B2_core_array,char **general_species_registry,
           
           int entries_UFF,char **species_UFF_array,double *D_UFF_array,double *x_UFF_array,
           int entries_DRE,char **species_DRE_array,double *D0_DRE_array,double *Rvdw0_DRE_array,
           
           int *host_id_array,
           double *pe_array,int counter,
           
           double T,int check_with_cells,
           
           int FF
           
           );

void pe_calc_module_v1(int general_atoms,int general_bonds,int general_atom_types,int entries_UFF,int entries_DRE,double cutoff,
                       double xlo,double ylo,double zlo,double xhi,double yhi,double zhi,double lx,double ly,double lz,double lcx,double lcy,double lcz,
                       int Mx,int My,int M3,
                       int *head,int **neighbors,double *pe_array,int counter,
                       double *master_x_array,double *master_y_array,double *master_z_array,int *host_id_array,
                       int *master_B1_core_array,int *master_B2_core_array,char **master_species_array,char **general_species_registry,
                       char **species_UFF_array,double *D_UFF_array,double *x_UFF_array,char **species_DRE_array,double *D0_DRE_array,double *Rvdw0_DRE_array,
                       int FF);

void calculate_dihedrals(char *current_folder,int general_dihedrals,int backbones,double *master_x_array,double *master_y_array,double *master_z_array,
                         int **general_dihedrals_registry,char **backbone_id,char **master_atom_name_array,double **bb_pdf_entries,double *bb_pdf_dx);

