
#include "builder.h"
#include "builder3.h"
#include "builder3_vars.h"
#include<time.h>

int main(int argc,char **argv)
{
    //--- Initialization -------------------------------------------------------
    clock_t c0f,c1f;time_t mytime;                              // timing variables
    mytime = time(NULL);sprintf(init_time,"%s",ctime(&mytime)); // initialize total execution time variable and store entry to char
    c0f=clock();                                                // header timing initialization
    //... Flags and command line arguments .....................................
    write_interval=1;              // defines the frequency for file dumps
    global_topo_flag=0;             // switches on total topology derivation for debugging purposes
    timer=0;                        // switches on module timing
    nup_flag=0;                     // switches on the nup loop    <-- DEPRECATED -- leave equal to 0 --!!
    lmp=atoi(argv[2]);              // switches on LAMMPS minimization
    nproc=atoi(argv[3]);            // defines number of processes for LAMMPS MPI calls
    myseed=atoi(argv[4]);           // seed for rand()
    CBMCG_flag=atoi(argv[5]);       // switches on MC growth
    CBMCG_N_trial=40;               // number of MC trial dihedrals
    T=atof(argv[6]);                // MC temperature
    check_with_cells=0;             // switches on lj pe calculation with linked list cells for debugging
    perc_flag=1;                    // progress percentage flag
    reshuffle_euler=atoi(argv[7]);  // randomize Euler angles
    reshuffle_coms=atoi(argv[8]);   // randomize CoM positions
    rand_disp=atof(argv[9]);        // max radius for CoM displacement
    FF=1;
    //..........................................................................
    // Random number generator initialization
    srand(myseed);
    // I.1
    builder2_initialize(current_folder,argv[1],&molecules,&mol_types,&xlo,&xhi,&ylo,&yhi,&zlo,&zhi,
                        &core_file,&type_array,&chains_array,&mol_id_array,&mol_type_array,&com_flag_array,
                        &X,&Y,&Z,&alpha,&beta,&my_gamma,
                        &total_chains,&quad_matrix,
                        &master_growth_registry_rows,&master_growth_registry_cols,&master_species_rows,&master_topo_rows,
                        &MGR_start,&MGR_stop,&MA_start,&MA_stop,&MB_start,&MB_stop,
                        &entries_UFF,
                        &species_UFF_array,
                        &D_UFF_array,&x_UFF_array,&r0_UFF_array,&chi_UFF_array,&Zstar_UFF_array,
                        &theta0_UFF_array,&Vtor_UFF_array,&Utor_UFF_array,&zeta_UFF_array,
                        &entries_DRE,
                        &species_DRE_array,
                        &R0_DRE_array,&theta0_DRE_array,&Rvdw0_DRE_array,&D0_DRE_array,
                        &backbones,&backbone_id,&bb_dx,&bb_entries,&bb_lines_max,
                        &tacticity,&tacticity_id,&tacticity_type,
                        &backbone_type_counter,&backbone_type_array,
                        mpi_from_config,lmp_from_config,
                        &bb_pdf_dx,&bb_pdf_entries,&bb_pdf_lines_max);
    counter=0;  // counts eligible molecule loop steps                                      // move to builder2_initialize()
    //... Console out ..........................................................
    printf("\n$- parameters ------------------------------------------\n");
    printf("            lammps:\t%d\n",lmp);
    printf("             nproc:\t%d\n",nproc);
    printf("             CBMCG:\t%d\tCBMCG_N_trial:\t%d\n",CBMCG_flag,CBMCG_N_trial);
    printf("        write step:\t%d\n",write_interval);
    printf("             gtopo:\t%d\n",global_topo_flag);
    printf("       cells check:\t%d\n",check_with_cells);
    printf("         perc flag:\t%d\n",perc_flag);
    printf("   reshuffle Euler:\t%d\n",reshuffle_euler);
    printf("    reshuffle CoMs:\t%d\trand_disp:\t%lf\n",reshuffle_coms,rand_disp);
    printf("             timer:\t%d\n",timer);
    printf("              seed:\t%d\n",myseed);
    printf("               mpi:\t%s\n",mpi_from_config);
    printf("               lmp:\t%s\n",lmp_from_config);
    printf("                 T:\t%lf\t(K)\n",T);
    printf("                FF:\t%d\n",FF);
    //..........................................................................
    // I.2
    // Create a contiguous registry for the building process: side chain info
    growth_reg_init(type_array, core_file, chains_array, mol_types, total_chains, quad_matrix, &rows, &cols, &species_rows, &topo_rows,
                    &master_growth_registry_rows,&master_growth_registry_cols,&master_species_rows,&master_topo_rows,
                    &master_species_chain_array,&master_B_type_chain_array,&master_growth_registry,
                    &master_B1_chain_array,&master_B2_chain_array,
                    &master_atom_name_chain_array,
                    MGR_start, MGR_stop, MA_start, MA_stop, MB_start, MB_stop,
                    molecules,mol_type_array,&max_growth_step_array,&growth_step_array);
    total_steps=0;for(i=0;i<molecules;++i)total_steps=total_steps+max_growth_step_array[i]; // move to growth_reg_init()
    pe_array=(double*)malloc((total_steps+1)*sizeof(double));                               // move to growth_reg_init()
    // I.3
    if(tacticity>0)
        initialize_tacticity(master_species_rows, master_atom_name_chain_array, mol_types, MA_start, MA_stop, tacticity, molecules, tacticity_id, mol_type_array, tacticity_type,
                             &yet_one_more_t_array, &tacticity_counter_array, &tacticity_counter_array_v2, &tacticity_search, &tacticity_array_length, &tacticity_array,
                             &t_tracker_1,&t_tracker_2,&t_tracker_3,&t_tracker_4,&t_tracker_mol);
    // I.4
    // Create a contiguous registry for the building process: core molecules info
    //
    // Here we introduce two extra int* arrays: atoms_per_molecule and atom_scaling_array
    // Atoms_per_molecule clearly holds the number of atoms per molecule...
    // Atom_scaling_array is derived from atoms_per_molecule and is used to map internal ids to general ids:
    // for a given molecule i, adding atom_scaling_array[i] to its internal id, results to the general id
    //
    contiguous_atom(&general_atoms,&general_bonds,molecules,mol_types,mol_id_array,mol_type_array,type_array,core_file,
                    current_folder,&master_mol_array,&master_id_internal_array,&master_x_array,&master_y_array,&master_z_array,
                    &master_q_array,
                    &master_species_array,&master_id_array,&master_mol_type_array,total_chains,quad_matrix,&host_id_array,&chain_id_array,
                    &atoms_per_molecule,&atom_scaling_array,
                    &master_atom_name_array);
    contiguous_bonds(general_atoms,general_bonds,molecules,mol_types,mol_type_array,type_array,current_folder,core_file,
                     master_id_array,master_id_internal_array,master_mol_array,
                     &master_B_type_core_array,&master_B_ID_core_array,&master_B1_core_array,&master_B2_core_array,&mol_id_array_BONDS,
                     &master_B1_internal_core_array,&master_B2_internal_core_array);
    // I.5
    /*
    // initialize_Cn.
    if(backbones>0){
        initialize_Cn(molecules,max_growth_step_array,quad_matrix,mol_type_array,
                      &bb_growth_step_array,&added_backbone_flag,&l_array,&l2_array,&Re_start,&Re_stop,&prev_link_atom_array);
    }
    */
    // I.6
    // Initial placement and rotations. -- !!!!! without PBC !!!!!
    builder2_euler(current_folder,argv[1],general_atoms,molecules,
                   com_flag_array,master_mol_array,mol_id_array,master_id_internal_array,master_species_array,
                   X,Y,Z,alpha,beta,my_gamma,
                   master_x_array,master_y_array,master_z_array,reshuffle_euler,reshuffle_coms,rand_disp,xlo,ylo,zlo,xhi,yhi,zhi);
    // I.7
    write_xyz_preview(current_folder,general_atoms,master_species_array,master_x_array,master_y_array,master_z_array,xlo,xhi,ylo,yhi,zlo,zhi,0);
    // I.8
    // General topology
    topo(general_atoms,master_species_array,general_bonds,master_B1_core_array,master_B2_core_array,master_B_type_core_array,
         &general_species_registry,&general_bond_types_registry,&general_angle_types_registry,&general_dihedral_types_registry,&general_improper_types_registry,
         &general_neighbors_registry,&general_bonds_registry,&general_angles_registry,&general_dihedrals_registry,&general_impropers_registry,
         &general_angles,&general_dihedrals,&general_impropers,
         &general_atom_types,&general_bond_types,&general_angle_types,&general_dihedral_types,&general_improper_types,
         1,current_folder,"topo_0_0.txt");
    // I.9
    tmap_topo_boundaries_init(molecules,master_mol_array,master_species_array,general_atoms,general_atom_types,general_species_registry,
                              general_bonds,general_bond_types,general_bond_types_registry,general_bonds_registry,
                              general_angles,general_angle_types,general_angle_types_registry,general_angles_registry,
                              general_dihedrals,general_dihedral_types,general_dihedral_types_registry,general_dihedrals_registry,
                              general_impropers,general_improper_types,general_improper_types_registry,general_impropers_registry,
                              &stmap_rows,&stmap,&scm,&btmap_rows,&btmap,&bcm,&atmap_rows,&atmap,&acm,&dtmap_rows,&dtmap,&dcm,&itmap_rows,&itmap,&icm,&topo_boundaries,
                              0);
    //
    if(CBMCG_flag==1){
        // only need to call once (supercell does not change)
        linked_list_cell_initialize(cutoff,xlo,xhi,ylo,yhi,zlo,zhi,general_atoms,&Mx,&My,&Mz,&M3,&lx,&ly,&lz,&lcx,&lcy,&lcz,&neighbors,&head);
        // console out
        printf("             cells:\t[%dx%dx%d][M3:%d]\n",Mx,My,Mz,M3);
        printf("               lcx:\t%lf\t(Ang)\n",lcx);
        printf("               lcy:\t%lf\t(Ang)\n",lcy);
        printf("               lcz:\t%lf\t(Ang)\n",lcz);
        //
        if(check_with_cells==1)
            pe_calc_module_v1(general_atoms,general_bonds,general_atom_types,entries_UFF,entries_DRE,cutoff,xlo,ylo,zlo,xhi,yhi,zhi,lx,ly,lz,lcx,lcy,lcz,Mx,My,M3,
                              head,neighbors,pe_array,counter,master_x_array,master_y_array,master_z_array,host_id_array,
                              master_B1_core_array,master_B2_core_array,master_species_array,general_species_registry,
                              species_UFF_array,D_UFF_array,x_UFF_array,species_DRE_array,D0_DRE_array,Rvdw0_DRE_array,FF);
    }
    // console out
    printf("$-------------------------------------------------------\n");
    // write mol2 file for init state
    write_mol2_file(current_folder,general_atoms,general_bonds,master_id_array,master_atom_name_array,master_x_array,master_y_array,master_z_array,master_species_array,
                    master_mol_array,master_B1_core_array,master_B2_core_array,master_B_type_core_array,"preview_init.mol2");
    // console out
    printf("\n$- initialization done! --------------------------------\n");
    // Timing: header
    c1f=clock();printf("\n$ header execution time: %lf\n\n",(double)(c1f-c0f)/CLOCKS_PER_SEC);
    //
    // MASTER LOOP -------------------------------------------------------------
    //
    // @@@@
    // Master building loop: init max step and counter
    MAXMAX=max_growth_step_array[0];
    for(i=0;i<molecules;++i)if(max_growth_step_array[i]>MAXMAX)MAXMAX=max_growth_step_array[i];
    //counter=0;  // counts eligible molecule loop steps; commented because its definition has been moved up top!
    // @@@@
    // Master loop
    for(m=0;m<MAXMAX;++m)
    {
        // Molecules loop
        for(i=0;i<molecules;++i)
        {
            // @@@@
            //if(backbones>0)added_backbone_flag[i]=0;
            // @@@@
            // Eligible molecule for growth
            if(growth_step_array[i]<max_growth_step_array[i])
            {
                // @@@@
                counter=counter+1;
                // @@@@
                // M.1 ---------------------------------------------------------
                // locate and store segment
                if(timer==1)c0f=clock();
                segments(i+1,general_atoms,master_id_array,master_mol_array,master_species_array,master_q_array,master_x_array,master_y_array,master_z_array,
                         host_id_array,chain_id_array,master_id_internal_array,master_mol_type_array,
                         &atoms_block_size,&segment_id_array,&segment_mol_array,&segment_species_array,&segment_q_array,&segment_x_array,&segment_y_array,&segment_z_array,
                         &segment_host_id_array,&segment_chain_id_array,&segment_id_internal_array,&segment_mol_type_array,
                         general_bonds,master_B_ID_core_array,master_B1_core_array,master_B2_core_array,mol_id_array_BONDS,master_B1_internal_core_array,master_B2_internal_core_array,master_B_type_core_array,
                         &bonds_block_size,&segment_B_ID_core_array,&segment_B1_core_array,&segment_B2_core_array,&segment_mol_id_array_BONDS,&segment_B1_internal_core_array,&segment_B2_internal_core_array,&segment_B_type_core_array,
                         &segment_atom_name_core,&segment_subst_name_core,&segment_B_ID_core,
                         &atoms_first_entry,&atoms_last_entry,&bonds_first_entry,&bonds_last_entry,
                         master_atom_name_array);
                if(timer==1){c1f=clock();printf("$ segments execution time: %lf\n",(double)(c1f-c0f)/CLOCKS_PER_SEC);}
                // M.2 ---------------------------------------------------------
                if(timer==1)c0f=clock();
                retrieve_mol_type_data(mol_type_array[i],mol_types,type_array,MGR_start,MGR_stop,MA_start,MA_stop,MB_start,MB_stop,
                                       total_chains,quad_matrix,chains_array,master_species_chain_array,master_B1_chain_array,master_B2_chain_array,master_B_type_chain_array,
                                       master_growth_registry_cols,master_growth_registry,
                                       &grows,&init_bo_array,&growth_registry,
                                       &total_chain_atoms,&total_chain_bonds,&chain_atoms_array,&chain_bonds_array,&species_chain_array,
                                       &B1_chain_array,&B2_chain_array,&B_type_chain_array,
                                       master_atom_name_chain_array,&atom_name_chain_array);
                if(timer==1){c1f=clock();printf("$ retrieve_mol_type_data execution time: %lf\n",(double)(c1f-c0f)/CLOCKS_PER_SEC);}
                // M.3 ---------------------------------------------------------
                //
                //-- TM: alloc, populate (in the master loop, before the nup loop)
                //-- TMB: alloc, populate (in the master loop, before the nup loop)
                //
                if(timer==1)c0f=clock();
                // remove _B!!
                topo_molecule_alloc_populate(i,atoms_block_size,topo_boundaries,
                                             general_bond_types,general_angle_types,general_dihedral_types,general_improper_types,
                                             general_bonds_registry,general_angles_registry,general_dihedrals_registry,general_impropers_registry,
                                             general_bond_types_registry,general_angle_types_registry,general_dihedral_types_registry,general_improper_types_registry,
                                             &molecule_atoms,&molecule_bonds,&molecule_angles,&molecule_dihedrals,&molecule_impropers,
                                             &molecule_atom_types,&molecule_bond_types,&molecule_angle_types,&molecule_dihedral_types,&molecule_improper_types,
                                             &molecule_bonds_registry,&molecule_angles_registry,&molecule_dihedrals_registry,&molecule_impropers_registry,
                                             &molecule_species_registry,
                                             &molecule_bond_types_registry,&molecule_angle_types_registry,&molecule_dihedral_types_registry,&molecule_improper_types_registry,
                                             &molecule_atoms_B,&molecule_bonds_B,&molecule_angles_B,&molecule_dihedrals_B,&molecule_impropers_B,
                                             &molecule_atom_types_B,&molecule_bond_types_B,&molecule_angle_types_B,&molecule_dihedral_types_B,&molecule_improper_types_B,
                                             &molecule_bonds_registry_B,&molecule_angles_registry_B,&molecule_dihedrals_registry_B,&molecule_impropers_registry_B,
                                             &molecule_species_registry_B,
                                             &molecule_bond_types_registry_B,&molecule_angle_types_registry_B,&molecule_dihedral_types_registry_B,&molecule_improper_types_registry_B,
                                             0);
                
               
                
                if(timer==1){c1f=clock();printf("$ topo_molecule_alloc_populate execution time: %lf\n",(double)(c1f-c0f)/CLOCKS_PER_SEC);}
                // M.4 ---------------------------------------------------------
                //-- TGBp: alloc, populate (in the master loop, before the nup loop); p: PARTIAL backup!
                if(timer==1)c0f=clock();
                // function is to be removed!!
                topo_general_backup(general_bonds,general_angles,general_dihedrals,general_impropers,
                                    general_atom_types,general_bond_types,general_angle_types,general_dihedral_types,general_improper_types,
                                    general_species_registry,
                                    general_bond_types_registry,general_angle_types_registry,general_dihedral_types_registry,general_improper_types_registry,
                                    &general_bonds_B,&general_angles_B,&general_dihedrals_B,&general_impropers_B,
                                    &general_atom_types_B,&general_bond_types_B,&general_angle_types_B,&general_dihedral_types_B,&general_improper_types_B,
                                    &general_species_registry_B,
                                    &general_bond_types_registry_B,&general_angle_types_registry_B,&general_dihedral_types_registry_B,&general_improper_types_registry_B,
                                    general_atoms,general_neighbors_registry,&general_neighbors_registry_B);
                if(timer==1){c1f=clock();printf("$ topo_general_backup execution time: %lf\n",(double)(c1f-c0f)/CLOCKS_PER_SEC);}
                // M.5 ---------------------------------------------------------
                if(timer==1)c0f=clock();
                builder_loop_preamble(growth_step_array[i],&chain_index,growth_registry,atoms_block_size,segment_host_id_array,segment_chain_id_array,&host_index,
                                      &max_1_2_core,segment_B1_internal_core_array,segment_B2_internal_core_array,bonds_block_size,GS_array,segment_species_array,&GSs,
                                      &one_two_core_vector,&GS_reg,
                                      chain_atoms_array,init_bo_array,chain_bonds_array,B1_chain_array,B2_chain_array,B_type_chain_array,
                                      species_chain_array,&nup,&unique_permutations,
                                      &new_entries_atom,&new_entries_bond,&atom_id_new,&x_new,&y_new,&z_new,
                                      &B1_new,&B2_new,
                                      &atoms_core_B,&bonds_core_B,&atom_id_core_B,&subst_id_core_B,&x_core_B,&y_core_B,&z_core_B,&charge_core_B,
                                      &atom_name_core_B,&atom_type_core_B,&subst_name_core_B,&host_id_array_B,&chain_id_array_B,
                                      &B_ID_core_B,&B1_core_B,&B2_core_B,&B_type_core_B,
                                      segment_id_internal_array,segment_mol_array,segment_x_array,segment_y_array,segment_z_array,segment_q_array,
                                      segment_atom_name_core,segment_subst_name_core,
                                      segment_B_ID_core,segment_B_type_core_array,
                                      &current_GS_reg,&score,
                                      &new_entries_atom_total,&new_entries_bond_total,
                                      &one_two_counter,&one_two_list,&one_three_counter,&one_three_list,&one_four_counter,&one_four_list,
                                      atom_name_chain_array,
                                      &local_topo_array,&local_topo_array_B);
                //
                if(timer==1){c1f=clock();printf("$ builder_loop_preamble execution time: %lf\n",(double)(c1f-c0f)/CLOCKS_PER_SEC);}
                // M.6 ---------------------------------------------------------
                if(timer==1)c0f=clock();
                resize_master_registries(&general_atoms,new_entries_atom_total,atoms_last_entry,
                                         &master_id_array,&master_mol_array,&host_id_array,&chain_id_array,&master_id_internal_array,&master_mol_type_array,
                                         &master_q_array,&master_x_array,&master_y_array,&master_z_array,
                                         &master_species_array,
                                         &general_bonds,new_entries_bond_total,bonds_last_entry,
                                         &master_B_ID_core_array,&master_B1_core_array,&master_B2_core_array,&mol_id_array_BONDS,
                                         &master_B1_internal_core_array,&master_B2_internal_core_array,
                                         &master_B_type_core_array,
                                         &master_atom_name_array);
                if(timer==1){c1f=clock();printf("$ resize_master_registries execution time: %lf\n",(double)(c1f-c0f)/CLOCKS_PER_SEC);}
                // @@@@
                if(nup_flag==1)
                {
                    void;
                }
                else
                {
                    // this is the else clause to bypass the nup loop
                    //l=nup-1;
                    l=0;
                }
                // @@@@
                // M.7 ---------------------------------------------------------
                if(timer==1)c0f=clock();
                builder_core_funct(chain_index,l,nup,
                                   GSs,host_index,max_1_2_core,
                                   new_entries_atom,new_entries_bond,
                                   &atoms_block_size,&bonds_block_size,
                                   atoms_core_B,bonds_core_B,
                                   one_two_core_vector,GS_array,
                                   GS_reg,current_GS_reg,unique_permutations,
                                   &segment_id_internal_array,&segment_mol_array,&segment_B_ID_core,&segment_B1_internal_core_array,&segment_B2_internal_core_array,
                                   &segment_host_id_array,&segment_chain_id_array,
                                   &segment_x_array,&segment_y_array,&segment_z_array,&segment_q_array,
                                   &segment_atom_name_core,&segment_species_array,&segment_subst_name_core,&segment_B_type_core_array,
                                   atom_id_core_B,subst_id_core_B,B_ID_core_B,B1_core_B,B2_core_B,host_id_array_B,chain_id_array_B,
                                   x_core_B,y_core_B,z_core_B,charge_core_B,
                                   atom_name_core_B,atom_type_core_B,subst_name_core_B,B_type_core_B,
                                   atom_id_new,x_new,y_new,z_new,B1_new,B2_new,
                                   &local_topo_array,local_topo_array_B,
                                   &remove_N,&added_atoms,&added_atoms_array,
                                   &rescale_array,perc_flag);
                if(timer==1){c1f=clock();printf("$ builder_core_funct execution time: %lf\n",(double)(c1f-c0f)/CLOCKS_PER_SEC);}
                //
                // !!!!! come back later !!!!!
                //
                // M.8 ---------------------------------------------------------
                if(remove_N>0){
                    if(timer==1)c0f=clock();
                    remove_deleted(molecule_atoms,rescale_array,
                                   segment_species_array,
                                   atoms_block_size,bonds_block_size,
                                   segment_B1_internal_core_array,segment_B2_internal_core_array,segment_B_type_core_array,
                                   &molecule_atom_types,
                                   &molecule_bond_types,&molecule_angle_types,&molecule_dihedral_types,&molecule_improper_types,
                                   &molecule_bonds,&molecule_angles,&molecule_dihedrals,molecule_impropers,
                                   &molecule_species_registry,
                                   &molecule_bond_types_registry,&molecule_angle_types_registry,&molecule_dihedral_types_registry,&molecule_improper_types_registry,
                                   &molecule_bonds_registry,&molecule_angles_registry,&molecule_dihedrals_registry,molecule_impropers_registry);
                    if(timer==1){c1f=clock();printf("$ remove_deleted execution time: %lf\n",(double)(c1f-c0f)/CLOCKS_PER_SEC);}
                }
                // @@@@
                molecule_atoms=atoms_block_size;
                delta_atoms=molecule_atoms-molecule_atoms_B;
                // @@@@
                // M.9 ---------------------------------------------------------
                if(timer==1)c0f=clock();
                compute_local(current_folder,segment_x_array,segment_y_array,segment_z_array,atoms_block_size,bonds_block_size,local_topo_array,
                              segment_species_array,segment_B1_internal_core_array,segment_B2_internal_core_array,segment_B_type_core_array,
                              &local_atoms,&local_bonds,&local_species,&local_B1,&local_B2,&local_B_type,&local_map);
                if(timer==1){c1f=clock();printf("$ compute_local execution time: %lf\n",(double)(c1f-c0f)/CLOCKS_PER_SEC);}
                // M.10 --------------------------------------------------------
                if(timer==1)c0f=clock();
                topo(local_atoms,local_species,local_bonds,local_B1,local_B2,local_B_type,
                     &local_species_registry,&local_bond_types_registry,&local_angle_types_registry,&local_dihedral_types_registry,&local_improper_types_registry,
                     &local_neighbors_registry,&local_bonds_registry,&local_angles_registry,&local_dihedrals_registry,&local_impropers_registry,
                     &local_angles,&local_dihedrals,&local_impropers,
                     &local_atom_types,&local_bond_types,&local_angle_types,&local_dihedral_types,&local_improper_types,
                     0,current_folder,"topo_local.txt");
                if(timer==1){c1f=clock();printf("$ topo_local execution time: %lf\n",(double)(c1f-c0f)/CLOCKS_PER_SEC);}
                // @@@@
                for(j=0;j<local_atoms;++j)free(local_species[j]);free(local_species);free(local_B1);free(local_B2);for(j=0;j<local_bonds;++j)free(local_B_type[j]);free(local_B_type);
                // @@@@
                // M.11 --------------------------------------------------------
                if(timer==1)c0f=clock();
                /*
                printf("molecule: bond types:\n");
                for(j=0;j<molecule_bond_types;++j)printf("[%d]\t(%s)--%s--(%s)\n",j+1,molecule_bond_types_registry[j][0],molecule_bond_types_registry[j][1],molecule_bond_types_registry[j][2]);
                printf("molecule: bonds:\n");
                for(j=0;j<molecule_bonds;++j)printf("[%d]\t%d\t%d\t%d\t|\t(%s)--%s--(%s)\n",j+1,
                                                    molecule_bonds_registry[j][0],
                                                    molecule_bonds_registry[j][1],
                                                    molecule_bonds_registry[j][2],
                                                    molecule_bond_types_registry[molecule_bonds_registry[j][0]-1][0],
                                                    molecule_bond_types_registry[molecule_bonds_registry[j][0]-1][1],
                                                    molecule_bond_types_registry[molecule_bonds_registry[j][0]-1][2]);
                */
                alter_molecular_topo(local_bonds,local_angles,local_dihedrals,local_impropers,
                                     local_atom_types,
                                     local_map,
                                     local_species_registry,
                                     local_bonds_registry,local_angles_registry,local_dihedrals_registry,local_impropers_registry,
                                     local_bond_types_registry,local_angle_types_registry,local_dihedral_types_registry,local_improper_types_registry,
                                     &molecule_bonds,&molecule_angles,&molecule_dihedrals,&molecule_impropers,
                                     &molecule_atom_types,&molecule_bond_types,&molecule_angle_types,&molecule_dihedral_types,&molecule_improper_types,
                                     &molecule_species_registry,
                                     &molecule_bonds_registry,&molecule_angles_registry,&molecule_dihedrals_registry,&molecule_impropers_registry,
                                     &molecule_bond_types_registry,&molecule_angle_types_registry,&molecule_dihedral_types_registry,&molecule_improper_types_registry,
                                     0);
                /*
                printf("AFTER - molecule: bond types:\n");
                for(j=0;j<molecule_bond_types;++j)printf("[%d]\t(%s)--%s--(%s)\n",j+1,molecule_bond_types_registry[j][0],molecule_bond_types_registry[j][1],molecule_bond_types_registry[j][2]);
                printf("AFTER - molecule: bonds:\n");
                for(j=0;j<molecule_bonds;++j)printf("[%d]\t%d\t%d\t%d\t|\t(%s)--%s--(%s)\n",j+1,
                                                    molecule_bonds_registry[j][0],
                                                    molecule_bonds_registry[j][1],
                                                    molecule_bonds_registry[j][2],
                                                    molecule_bond_types_registry[molecule_bonds_registry[j][0]-1][0],
                                                    molecule_bond_types_registry[molecule_bonds_registry[j][0]-1][1],
                                                    molecule_bond_types_registry[molecule_bonds_registry[j][0]-1][2]);
                */
                if(timer==1){c1f=clock();printf("$ alter_molecular_topo execution time: %lf\n",(double)(c1f-c0f)/CLOCKS_PER_SEC);}
                // M.12 --------------------------------------------------------
                if(timer==1)c0f=clock();
                alter_general_topo(i,molecules,
                                   molecule_atoms, molecule_atoms_B, delta_atoms, molecule_atom_types,
                                   molecule_bonds, molecule_bonds_B, molecule_bond_types,
                                   molecule_angles, molecule_angles_B, molecule_angle_types,
                                   molecule_dihedrals, molecule_dihedrals_B, molecule_dihedral_types,
                                   molecule_impropers, molecule_impropers_B, molecule_improper_types,
                                   atoms_per_molecule, atom_scaling_array,
                                   molecule_species_registry,
                                   molecule_bond_types_registry, molecule_bonds_registry,
                                   molecule_angle_types_registry,molecule_angles_registry,
                                   molecule_dihedral_types_registry, molecule_dihedrals_registry,
                                   molecule_improper_types_registry, molecule_impropers_registry,
                                   &general_bonds, &general_angles, &general_dihedrals, &general_impropers,
                                   topo_boundaries,
                                   &stmap_rows, &stmap, &scm, &general_atom_types, &general_species_registry,
                                   &btmap_rows, &btmap, &bcm,&general_bond_types, &general_bonds_registry, &general_bond_types_registry, &btmapping_matrix,
                                   &atmap_rows, &atmap, &acm, &general_angle_types, &general_angles_registry, &general_angle_types_registry, &atmapping_matrix,
                                   &dtmap_rows, &dtmap, &dcm, &general_dihedral_types, &general_dihedrals_registry, &general_dihedral_types_registry, &dtmapping_matrix,
                                   &itmap_rows, &itmap, &icm, &general_improper_types, &general_impropers_registry, &general_improper_types_registry, &itmapping_matrix,
                                   0);
                if(timer==1){c1f=clock();printf("$ alter_general_topo execution time: %lf\n",(double)(c1f-c0f)/CLOCKS_PER_SEC);}
                // M.13 --------------------------------------------------------
                if(timer==1)c0f=clock();
                free_locals();
                if(timer==1){c1f=clock();printf("$ free_locals execution time: %lf\n",(double)(c1f-c0f)/CLOCKS_PER_SEC);}
                // M.14 --------------------------------------------------------
                if(timer==1)c0f=clock();
                update_master(i,general_atoms,general_bonds,atoms_block_size,bonds_block_size,atoms_first_entry,bonds_first_entry,master_id_array,
                              master_mol_array,host_id_array,chain_id_array,master_id_internal_array,master_mol_type_array,
                              master_q_array,master_x_array,master_y_array,master_z_array,master_species_array,
                              master_B_type_core_array,mol_id_array_BONDS,master_B1_core_array,master_B2_core_array,master_B1_internal_core_array,master_B2_internal_core_array,
                              segment_mol_array,segment_host_id_array,segment_chain_id_array,segment_id_internal_array,mol_type_array,
                              segment_q_array,segment_x_array,segment_y_array,segment_z_array,segment_species_array,
                              segment_B_type_core_array,segment_B1_internal_core_array,segment_B2_internal_core_array,
                              atom_scaling_array,
                              master_atom_name_array,segment_atom_name_core);
                if(timer==1){c1f=clock();printf("$ update_master execution time: %lf\n",(double)(c1f-c0f)/CLOCKS_PER_SEC);}
                // M.15 --------------------------------------------------------
                //
                // conditional code fragment A !!!!
                //
                if(timer==1)c0f=clock();
                update_neighbors_registry(i,general_atoms,delta_atoms,molecules,molecule_atoms,atoms_first_entry,atoms_last_entry,master_species_array,&general_neighbors_registry);
                if(timer==1){c1f=clock();printf("$ update_neighbors_registry execution time: %lf\n",(double)(c1f-c0f)/CLOCKS_PER_SEC);}
                // M.16 --------------------------------------------------------
                //
                // diagnostic for debug only!
                //
                if(global_topo_flag==1)
                {
                    topo(general_atoms,master_species_array,general_bonds,master_B1_core_array,master_B2_core_array,master_B_type_core_array,
                         &global_species,&global_bond_type_array,&global_angle_type_array,&global_dihedral_type_array,&global_improper_type_array,
                         &global_neighbors,&global_bonds_array,&global_angles_array,&global_dihedrals_array,&global_impropers_array,
                         &angles_core,&dihedrals_core,&impropers_core,
                         &atom_types,&bond_types,&angle_types,&dihedral_types,&improper_types,
                         //0,NULL,NULL);
                         1,current_folder,"topo.txt");
                    topo_diagnostic(general_bonds,general_bonds_registry,general_bond_types_registry,general_angles,general_angles_registry,general_angle_types_registry,
                                    general_dihedrals,general_dihedrals_registry,general_dihedral_types_registry,general_impropers,general_impropers_registry,general_improper_types_registry,
                                    general_bonds,global_bonds_array,global_bond_type_array,angles_core,global_angles_array,global_angle_type_array,
                                    dihedrals_core,global_dihedrals_array,global_dihedral_type_array,impropers_core,global_impropers_array,global_improper_type_array);
                    for(j=0;j<general_atoms;++j)
                        if(global_neighbors[j]!=general_neighbors_registry[j]){
                            printf("[%d] %s %d %d -- general neigh inconsistency!!\n",j+1,master_species_array[j],global_neighbors[j],general_neighbors_registry[j]);exit(-2);}
                    if(perc_flag==0)printf("$$$ topo diagnostic passed!!!\n");
                }
                //==============================================================
                // this is the part where conformation and tacticity is enforced ONCE!
                // - place it in an if-clause controlled by an external flag
                // - check mallocs in force_tacticity()!!
                if(CBMCG_flag==0){
                    // M.17 --------------------------------------------------------
                    if(backbones!=0){
                        if(timer==1)c0f=clock();
                        apply_backbone_dihedral(i,backbones,local_atoms,local_topo_array,local_map,segment_atom_name_core,backbone_id,atom_scaling_array,
                                                bb_entries,bb_dx,master_x_array,master_y_array,master_z_array,
                                                //Re_stop,bb_growth_step_array,prev_link_atom_array,added_backbone_flag,
                                                backbone_type_counter,backbone_type_array);
                        if(timer==1){c1f=clock();printf("$ apply_backbone_dihedral execution time: %lf\n",(double)(c1f-c0f)/CLOCKS_PER_SEC);}
                    }
                    // M.18 --------------------------------------------------------
                    //
                    // isotactic polymer: all dyads must be meso
                    // syndiotactic polymer: all dyads must be racemic
                    // atactic polymer: dyads take random meso or racemic states
                    //
                    if(tacticity!=0){
                        if(timer==1)c0f=clock();
                        force_tacticity(i,tacticity,mol_type_array,tacticity_search,local_atoms,segment_atom_name_core,local_map,local_topo_array,
                                        atom_scaling_array,backbone_id,yet_one_more_t_array,backbones,backbone_type_counter,backbone_type_array,
                                        
                                        molecules,tacticity_array_length,
                                        
                                        tacticity_counter_array,
                                        &general_t_counter,t_tracker_1,t_tracker_2,t_tracker_3,t_tracker_4,t_tracker_mol,
                                        tacticity_array,master_x_array,master_y_array,master_z_array,
                                        &ft);
                        if(timer==1){c1f=clock();printf("$ force_tacticity execution time: %lf\n",(double)(c1f-c0f)/CLOCKS_PER_SEC);}
                    }
                }
                //==============================================================
                // when M.17, M.18 are not triggered, configurational bias Monte Carlo growth (CBMCG) is employed - access via else{}
                // - the key here is to force tacticity - if applicable - for every backbone configurational step
                // - course of action: do not reuse M.17, M.18 routines; combine them into a new one!
                //
                // CBMCG
                //
                else
                {
                    if(timer==1)c0f=clock();
                    CBMCG(i,local_atoms,backbone_type_counter,backbones,CBMCG_N_trial,current_folder,
                          local_topo_array,local_map,segment_atom_name_core,backbone_type_array,atom_scaling_array,backbone_id,
                          bb_dx,bb_entries,bb_pdf_dx,bb_pdf_entries,
                          tacticity,mol_type_array,tacticity_search,yet_one_more_t_array,
                          molecules,tacticity_array_length,
                          tacticity_counter_array,
                          &general_t_counter,t_tracker_1,t_tracker_2,t_tracker_3,t_tracker_4,t_tracker_mol,
                          tacticity_array,
                          master_species_array,master_x_array,master_y_array,master_z_array,
                          cutoff,xlo,ylo,zlo,xhi,yhi,zhi,general_atoms,general_bonds,general_atom_types,
                          Mx,My,Mz,M3,lx,ly,lz,lcx,lcy,lcz,neighbors,
                          master_B1_core_array,master_B2_core_array,general_species_registry,
                          entries_UFF,species_UFF_array,D_UFF_array,x_UFF_array,
                          entries_DRE,species_DRE_array,D0_DRE_array,Rvdw0_DRE_array,
                          host_id_array,pe_array,counter,T,check_with_cells,FF);
                    if(timer==1){c1f=clock();printf("$ CBMCG execution time: %lf\n",(double)(c1f-c0f)/CLOCKS_PER_SEC);}
                } // CBMCG else
                //
                //==============================================================
                
                
                
                
                // M.19 --------------------------------------------------------
                if(timer==1)c0f=clock();
                free_molecules_loop();
                if(timer==1){c1f=clock();printf("$ free_molecules_loop execution time: %lf\n",(double)(c1f-c0f)/CLOCKS_PER_SEC);}
                // @@@@
                growth_step_array[i]=growth_step_array[i]+1;
                // @@@@
                // Console out:
                /*
                printf("State: %d | atoms: %d\n",counter,general_atoms);
                printf("Growth steps:\n");
                printf("mol\tm_type\tcurrent\tmax\n");
                printf("%d\t%d\t%d\t%d\n",i+1,mol_type_array[i],growth_step_array[i],max_growth_step_array[i]);
                */
                if(perc_flag==1){
                    printf("\rIn progress: %.0lf%%", ((double)m/(MAXMAX-1))*100.0);
                }
                else
                {
                    printf("State: %d | atoms: %d\n",counter,general_atoms);
                    printf("Growth steps:\n");
                    printf("mol\tm_type\tcurrent\tmax\n");
                    printf("%d\t%d\t%d\t%d\n",i+1,mol_type_array[i],growth_step_array[i],max_growth_step_array[i]);
                }
                //==============================================================
                if(CBMCG_flag==1){
                    if(timer==1)c0f=clock();
                    if(check_with_cells==1)
                        pe_calc_module_v1(general_atoms,general_bonds,general_atom_types,entries_UFF,entries_DRE,cutoff,xlo,ylo,zlo,xhi,yhi,zhi,lx,ly,lz,lcx,lcy,lcz,Mx,My,M3,
                                          head,neighbors,pe_array,counter,master_x_array,master_y_array,master_z_array,host_id_array,
                                          master_B1_core_array,master_B2_core_array,master_species_array,general_species_registry,
                                          species_UFF_array,D_UFF_array,x_UFF_array,species_DRE_array,D0_DRE_array,Rvdw0_DRE_array,FF);
                    if(timer==1){c1f=clock();printf("$ pe_calc execution time: %lf\n",(double)(c1f-c0f)/CLOCKS_PER_SEC);}
                }
                //==============================================================
            }
        }   // end of molecules loop
        // M.20 ----------------------------------------------------------------
        if(lmp==1)
        {
            minimize(current_folder,counter,l,
                     xlo,xhi,ylo,yhi,zlo,zhi,
                     master_id_array,master_mol_array,master_q_array,master_x_array,master_y_array,master_z_array,
                     general_atoms,general_bonds,general_angles,general_dihedrals,general_impropers,
                     general_atom_types,general_bond_types,general_angle_types,general_dihedral_types,general_improper_types,
                     master_species_array,
                     general_species_registry,
                     general_bonds_registry,general_angles_registry,general_dihedrals_registry,general_impropers_registry,
                     general_bond_types_registry,general_angle_types_registry,general_dihedral_types_registry,general_improper_types_registry,
                     entries_UFF,species_UFF_array,
                     D_UFF_array,x_UFF_array,r0_UFF_array,chi_UFF_array,Zstar_UFF_array,
                     theta0_UFF_array,Vtor_UFF_array,Utor_UFF_array,zeta_UFF_array,
                     entries_DRE,
                     species_DRE_array,
                     R0_DRE_array,theta0_DRE_array,Rvdw0_DRE_array,D0_DRE_array,
                     general_neighbors_registry,
                     nproc,
                     mpi_from_config,lmp_from_config,FF);
            //==============================================================
            if(CBMCG_flag==1){
                if(timer==1)c0f=clock();
                //printf("    $ post minimize:\n");
                if(check_with_cells==1)
                    pe_calc_module_v1(general_atoms,general_bonds,general_atom_types,entries_UFF,entries_DRE,cutoff,xlo,ylo,zlo,xhi,yhi,zhi,lx,ly,lz,lcx,lcy,lcz,Mx,My,M3,
                                      head,neighbors,pe_array,counter,master_x_array,master_y_array,master_z_array,host_id_array,
                                      master_B1_core_array,master_B2_core_array,master_species_array,general_species_registry,
                                      species_UFF_array,D_UFF_array,x_UFF_array,species_DRE_array,D0_DRE_array,Rvdw0_DRE_array,FF);
                if(timer==1){c1f=clock();printf("$ pe_calc execution time: %lf\n",(double)(c1f-c0f)/CLOCKS_PER_SEC);}
            }
            //==============================================================
        }
        // M.21 ----------------------------------------------------------------
        if(timer==1)c0f=clock();
        if((m+1)%write_interval==0)
        //if(counter%write_interval==0)
        {
            // xyz
            write_xyz_preview(current_folder,general_atoms,master_species_array,master_x_array,master_y_array,master_z_array,xlo,xhi,ylo,yhi,zlo,zhi,1);
            // mol2
            write_mol2_file(current_folder,general_atoms,general_bonds,master_id_array,master_atom_name_array,master_x_array,master_y_array,master_z_array,master_species_array,
                            master_mol_array,master_B1_core_array,master_B2_core_array,master_B_type_core_array,"preview.mol2");
        }
        if(timer==1){c1f=clock();printf("$ write execution time: %lf\n",(double)(c1f-c0f)/CLOCKS_PER_SEC);}
        // M.22 ----------------------------------------------------------------
        if(tacticity!=0){
            if(timer==1)c0f=clock();
            calculate_tacticity(molecules,general_t_counter,tacticity_array_length,atom_scaling_array,master_x_array,master_y_array,master_z_array,tacticity_counter_array,
                                tacticity_counter_array_v2,yet_one_more_t_array,t_tracker_1,t_tracker_2,t_tracker_3,t_tracker_4,t_tracker_mol,tacticity_array,
                                counter,current_folder);
            if(timer==1){c1f=clock();printf("$$$ calculate_tacticity execution time: %lf\n",(double)(c1f-c0f)/CLOCKS_PER_SEC);}
        }
        // M.23 ----------------------------------------------------------------
        /*
        if(backbones>0){
            if(timer==1)c0f=clock();
            calculate_Cn(current_folder,counter,molecules,general_bonds,backbones,max_growth_step_array,atom_scaling_array,master_B1_core_array,master_B2_core_array,
                         mol_id_array_BONDS,master_x_array,master_y_array,master_z_array,master_atom_name_array,backbone_id,l_array,l2_array,Re_start,Re_stop,bb_growth_step_array,
                         backbone_type_counter,backbone_type_array);
            if(timer==1){c1f=clock();printf("$$$ calculate_Cn execution time: %lf\n",(double)(c1f-c0f)/CLOCKS_PER_SEC);}
        }
        */
        
        //getchar();
        
    }   // end of master loop
    printf("\n\n");
    //
    //--------------------------------------------------------------------------
    //
    write_xyz_preview(current_folder,general_atoms,master_species_array,master_x_array,master_y_array,master_z_array,xlo,xhi,ylo,yhi,zlo,zhi,1);
    // mol2
    write_mol2_file(current_folder,general_atoms,general_bonds,master_id_array,master_atom_name_array,master_x_array,master_y_array,master_z_array,master_species_array,
                    master_mol_array,master_B1_core_array,master_B2_core_array,master_B_type_core_array,"preview.mol2");
    //
    //for(j=0;j<tacticity_array_length;++j)printf("[%d]\t%lf\t{%d}\t%d\t%d\t%d\t%d\t-- %d --\n",j+1,tacticity_array[j],t_tracker_mol[j],t_tracker_1[j],t_tracker_2[j],t_tracker_3[j],t_tracker_4[j],yet_one_more_t_array[t_tracker_mol[j]-1]);
    for(j=0;j<tacticity_array_length;++j)printf("[%d]\t%lf\n",j+1,tacticity_array[j]);
    //
    if(nproc==0 && MAXMAX==0)
    {
        minimize(current_folder,counter,l,
                 xlo,xhi,ylo,yhi,zlo,zhi,
                 master_id_array,master_mol_array,master_q_array,master_x_array,master_y_array,master_z_array,
                 general_atoms,general_bonds,general_angles,general_dihedrals,general_impropers,
                 general_atom_types,general_bond_types,general_angle_types,general_dihedral_types,general_improper_types,
                 master_species_array,
                 general_species_registry,
                 general_bonds_registry,general_angles_registry,general_dihedrals_registry,general_impropers_registry,
                 general_bond_types_registry,general_angle_types_registry,general_dihedral_types_registry,general_improper_types_registry,
                 entries_UFF,species_UFF_array,
                 D_UFF_array,x_UFF_array,r0_UFF_array,chi_UFF_array,Zstar_UFF_array,
                 theta0_UFF_array,Vtor_UFF_array,Utor_UFF_array,zeta_UFF_array,
                 entries_DRE,
                 species_DRE_array,
                 R0_DRE_array,theta0_DRE_array,Rvdw0_DRE_array,D0_DRE_array,
                 general_neighbors_registry,
                 nproc,
                 mpi_from_config,lmp_from_config,FF);
    }
    // calculate backbone dihedral angle distributions
    calculate_dihedrals(current_folder,general_dihedrals,backbones,master_x_array,master_y_array,master_z_array,
                        general_dihedrals_registry,backbone_id,master_atom_name_array,bb_pdf_entries,bb_pdf_dx);
    //==========================================================================
    // free memory
    free_final();
    //
    post_exec(myseed);
    //
    mytime = time(NULL);
    sprintf(final_time,"%s",ctime(&mytime));
    printf("$ timing info:\n");
    printf("$  started at: %s",init_time);
    printf("$    ended at: %s",final_time);
    printf("\n");
    //
    return 0;
}
