char current_folder[cmax_length],command[cmax_length];// strings
char **core_file;                                                           // array to store core mol2 filenames (size: mol_types)
char **quad_matrix;             // matrix to store contiguous side chain info (molecular type index, side chain mol2, graft_id, host_id, init_bo, side chain atoms, side chain bonds)
int i,j,l,m;       // indices and buffer
int molecules,mol_types;        // number of molecules and molecule types
int total_chains;               // chains counter over all molecular types
int *type_array,*chains_array;  // arrays to store molecular type index and number of side chains per molecular type (size: mol_types)
double xlo,xhi,ylo,yhi,zlo,zhi; // supercell parameters
int rows,cols;                  // stores output from growth_reg_init
int *MGR_start,*MGR_stop;   // arrays to store the boundaries inside **master_growth_registry (size: mol_types) and summation variable
int *max_growth_step_array,*growth_step_array;      // arrays to store the maximum and current growth index per molecule (size: molecules)
int *mol_id_array,*mol_type_array,*com_flag_array;  // molecular data (size: molecules)
double *X,*Y,*Z,*alpha,*beta,*my_gamma;                // molecular data (size: molecules)
// contiguous atomic data (size: atoms)
char **master_species_array;                            // species
int *master_id_internal_array;                          // unaltered atomic id
int *master_mol_array;                                  // molecule id
int *master_id_array,*master_mol_type_array;            // contiguous atomic id and molecular type arrays
double *master_x_array,*master_y_array,*master_z_array; // coordinates
//
int species_rows,topo_rows;                             // stores output from growth_reg_init
int *MA_start,*MA_stop,*MB_start,*MB_stop;    // arrays to store the boundaries inside atomic and bond contiguous registries (size: mol_types) and summation variables
// contiguous bonding data (size: bonds)
char **master_B_type_core_array;                        // contiguous bond topology info (type)
int *master_B_ID_core_array;
int *master_B1_core_array,*master_B2_core_array;        // contiguous bond topology info (B1,B2)
int *mol_id_array_BONDS;                                // array to hold the molecular id of the molecule to which the bond belongs
int *master_B1_internal_core_array;                     // unaltered bond topology info
int *master_B2_internal_core_array;                     // unaltered bond topology info
// mol_type isolated variables
int *chain_atoms_array,*chain_bonds_array;              // number atoms and bonds per side chain
int *init_bo_array;                                     // array to hold init bo per side chain
int grows,**growth_registry;                      // isolated growth_registry variables
int total_chain_atoms,total_chain_bonds;                // counter of atoms and bonds per molecular type
char **species_chain_array;                             // isolated atomic data: species
char **B_type_chain_array;                              // isolated bonding data: type
int *B1_chain_array,*B2_chain_array;                    // isolated bonding data: B1, B2
// building related variables
int *host_id_array,*chain_id_array;                     // key building vectors (size: atoms)
int master_growth_registry_rows,master_growth_registry_cols,master_species_rows,master_topo_rows;
char **master_species_chain_array,**master_B_type_chain_array;
int **master_growth_registry;
int *master_B1_chain_array,*master_B2_chain_array;
double *master_q_array;
// topo
char **global_species,***global_bond_type_array,***global_angle_type_array,***global_dihedral_type_array,***global_improper_type_array;
int *global_neighbors;
int **global_bonds_array,**global_angles_array,**global_dihedrals_array,**global_impropers_array;
int angles_core,dihedrals_core,impropers_core,atom_types,bond_types,angle_types,dihedral_types,improper_types;
// UFF data
int entries_UFF;                                                                // number of entries in the prm file
char **species_UFF_array;
double *r0_UFF_array,*theta0_UFF_array,*x_UFF_array,*D_UFF_array,*zeta_UFF_array,*Zstar_UFF_array,*chi_UFF_array,*Vtor_UFF_array,*Utor_UFF_array;
//
int atoms_block_size,bonds_block_size;
char **segment_species_array;
int *segment_id_array,*segment_mol_array,*segment_host_id_array,*segment_chain_id_array,*segment_id_internal_array,*segment_mol_type_array;
double *segment_q_array,*segment_x_array,*segment_y_array,*segment_z_array;
char **segment_B_type_core_array;
int *segment_B_ID_core_array,*segment_B1_core_array,*segment_B2_core_array,*segment_mol_id_array_BONDS,*segment_B1_internal_core_array,*segment_B2_internal_core_array;
char **segment_atom_name_core,**segment_subst_name_core;
int *segment_B_ID_core;
//
int *one_two_core_vector;                                                       // array to store the 1-2 neighbors of the host atom
char **GS_reg;
double *score;
int chain_index,host_index,max_1_2_core,GS_array[MAX_GS],GSs,nup,**unique_permutations,new_entries_atom,new_entries_bond;
int *atom_id_new,*B1_new,*B2_new;                                               // arrays to hold atomic ids and bonding info for the new atoms after a successful building step
double *x_new,*y_new,*z_new;                                                    // arrays to hold the coordinates of the new atoms after a successful building step
char **current_GS_reg;                                                          // matrix to hold GS_reg and extra info in order to build in the nup loop
// backup mol2 variables for the core molecu;e
char **atom_name_core_B,**atom_type_core_B,**subst_name_core_B,**B_type_core_B;
int atoms_core_B,bonds_core_B;
int *atom_id_core_B,*subst_id_core_B,*B_ID_core_B,*B1_core_B,*B2_core_B;
double *x_core_B,*y_core_B,*z_core_B,*charge_core_B;
int *host_id_array_B,*chain_id_array_B;
int atoms_first_entry,bonds_first_entry;
int atoms_last_entry,bonds_last_entry;
int new_entries_atom_total,new_entries_bond_total;
int counter;
int MAXMAX;
// lammps minimization variables
int lmp,nproc;
int one_two_counter;
int one_three_counter;
int one_four_counter;
int *one_two_list;
int *one_three_list;
int *one_four_list;
int *local_topo_array;
int *local_topo_array_B;
int local_atoms,local_bonds,*local_B1,*local_B2;
char **local_B_type;
char **local_species;
int *local_map;
char **local_species_registry, ***local_bond_types_registry, ***local_angle_types_registry, ***local_dihedral_types_registry, ***local_improper_types_registry;
int *local_neighbors_registry, **local_bonds_registry, **local_angles_registry, **local_dihedrals_registry, **local_impropers_registry;
int local_angles, local_dihedrals, local_impropers;
int local_atom_types, local_bond_types, local_angle_types, local_dihedral_types, local_improper_types;
int general_atoms,general_bonds;
char **general_species_registry, ***general_bond_types_registry, ***general_angle_types_registry, ***general_dihedral_types_registry, ***general_improper_types_registry;
int *general_neighbors_registry, **general_bonds_registry, **general_angles_registry, **general_dihedrals_registry, **general_impropers_registry;
int general_angles, general_dihedrals, general_impropers;
int general_atom_types, general_bond_types, general_angle_types, general_dihedral_types, general_improper_types;
int **topo_boundaries;
int molecule_atoms,molecule_bonds;
char **molecule_species_registry, ***molecule_bond_types_registry, ***molecule_angle_types_registry, ***molecule_dihedral_types_registry, ***molecule_improper_types_registry;
int **molecule_bonds_registry, **molecule_angles_registry, **molecule_dihedrals_registry, **molecule_impropers_registry;
int molecule_angles, molecule_dihedrals, molecule_impropers;
int molecule_atom_types, molecule_bond_types, molecule_angle_types, molecule_dihedral_types, molecule_improper_types;
int remove_N;
int molecule_atoms_B,molecule_bonds_B;
char **molecule_species_registry_B, ***molecule_bond_types_registry_B, ***molecule_angle_types_registry_B, ***molecule_dihedral_types_registry_B, ***molecule_improper_types_registry_B;
int **molecule_bonds_registry_B, **molecule_angles_registry_B, **molecule_dihedrals_registry_B, **molecule_impropers_registry_B;
int molecule_angles_B, molecule_dihedrals_B, molecule_impropers_B;
int molecule_atom_types_B, molecule_bond_types_B, molecule_angle_types_B, molecule_dihedral_types_B, molecule_improper_types_B;
int general_bonds_B, general_angles_B, general_dihedrals_B, general_impropers_B;
int general_atom_types_B, general_bond_types_B, general_angle_types_B, general_dihedral_types_B, general_improper_types_B;
char **general_species_registry_B, ***general_bond_types_registry_B, ***general_angle_types_registry_B, ***general_dihedral_types_registry_B, ***general_improper_types_registry_B;
int *atoms_per_molecule;
int *atom_scaling_array;
//
char **stmap;   // species type map array
int **scm;      // species count matrix
int stmap_rows;
char ***btmap;
int **bcm;
int btmap_rows;
int **btmapping_matrix;
char ***atmap;
int **acm;
int atmap_rows;
int **atmapping_matrix;
char ***dtmap;
int **dcm;
int dtmap_rows;
int **dtmapping_matrix;
char ***itmap;
int **icm;
int itmap_rows;
int **itmapping_matrix;
int delta_atoms;
int added_atoms,*added_atoms_array;
int *general_neighbors_registry_B;
char **master_atom_name_chain_array,**atom_name_chain_array,**master_atom_name_array;
int backbones;
char **backbone_id;
//
int write_interval;
int global_topo_flag;
int timer;
int nup_flag;
//
double *bb_dx,**bb_entries;
int bb_lines_max;
int myseed;
//
int *rescale_array;
//
int tacticity;
char **tacticity_id,**tacticity_type;
char **tacticity_search;
int *tacticity_counter_array;
int tacticity_array_length;
double *tacticity_array;
int *yet_one_more_t_array;
int *t_tracker_1,*t_tracker_2,*t_tracker_3,*t_tracker_4,general_t_counter=0;
int *t_tracker_mol;
int *tacticity_counter_array_v2;

//int *Re_start,*Re_stop,*bb_growth_step_array,*prev_link_atom_array,*added_backbone_flag;
//double *l_array,*l2_array;

//##########################################################################
int found,k;
char **central_bo;
//##########################################################################

int backbone_type_counter;
char **backbone_type_array;

//
int **one_two,**one_three,**one_four;

double cutoff=12.0; // send cutoff in minimize() as well...
double pe_bf,pe_llc;

double lx,ly,lz,lcx,lcy,lcz;
int M3,Mx,My,Mz;
int *head,*list,**neighbors;

int *master_nx_array,*master_ny_array,*master_nz_array;

double **nb_FF_for_lj;
int *nb_type_array;

char mpi_from_config[cmax_length],lmp_from_config[cmax_length];

double *bb_pdf_dx,**bb_pdf_entries;
int bb_pdf_lines_max;

int CBMCG_flag,CBMCG_N_trial;

int ft;

int pairs_bf,pairs_llc;

int total_steps;
double *pe_array;

int line,ii,diagnose;
int **one_two_local,**one_three_local,**one_four_local;
int **one_two_general,**one_three_general,**one_four_general;
double T;
int perc_flag;

int check_with_cells;
int reshuffle_euler,reshuffle_coms;
double rand_disp;
char init_time[cmax_length],final_time[cmax_length];

int entries_DRE;
char **species_DRE_array;
double *R0_DRE_array,*theta0_DRE_array,*Rvdw0_DRE_array,*D0_DRE_array;

int FF;


