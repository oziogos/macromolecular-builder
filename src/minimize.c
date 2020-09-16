
#include"builder.h"

void SYBYLtoDRE(int atoms_core,char **atom_type_core,char **atom_type_core_FF);

void DRE_topo(int general_atoms,int general_atom_types,
              int general_bonds,int *intermed_bond_types,
              int general_angles,int *intermed_angle_types,
              int general_dihedrals,int *intermed_dihedral_types,
              int general_impropers,int *intermed_improper_types,
              char **master_species_array,char **atom_type_core_FF,
              char **species_intermed_array,
              char ****global_bond_type_intermed_array,char ****global_angle_type_intermed_array,char ****global_dihedral_type_intermed_array,char ****global_improper_type_intermed_array,
              int **general_bonds_registry,int **general_angles_registry,int **general_dihedrals_registry,int **general_impropers_registry,
              int entries_DRE,char **species_DRE_array,
              double *R0_DRE_array,double *theta0_DRE_array,double *Rvdw0_DRE_array,double *D0_DRE_array,
              double ***nb_FF,double ***bonds_FF,double ***angles_FF,double ***dihedrals_FF,double ***impropers_FF,
              int **b_refine_array,int **a_refine_array,int **d_refine_array,int **i_refine_array,
              int *equiv_B_N,char ***equiv_B,int *equiv_A_N,char ***equiv_A,int *equiv_D_N,char ***equiv_D);

void DRE_write_lammps(char *current_folder,int counter,int l,
                      int general_atoms,int general_bonds,int general_angles,int general_dihedrals,int general_impropers,
                      int general_atom_types,int intermed_bond_types,int intermed_angle_types,int intermed_dihedral_types,int intermed_improper_types,
                      char **atom_type_core_FF,char **species_intermed_array,
                      double xlo,double xhi,double ylo,double yhi,double zlo,double zhi,
                      int *master_id_array,int *master_mol_array,double *master_q_array,double *master_x_array,double *master_y_array,double *master_z_array,
                      int **general_bonds_registry,int **general_angles_registry,int **general_dihedrals_registry,int **general_impropers_registry,
                      char ***global_bond_type_intermed_array,char ***global_angle_type_intermed_array,char ***global_dihedral_type_intermed_array,char ***global_improper_type_intermed_array,
                      double **nb_FF,double **bonds_FF,double **angles_FF,double **dihedrals_FF,double **impropers_FF,
                      int **b_refine_array,int **a_refine_array,int **d_refine_array,int **i_refine_array,
                      int *D_multi,
                      int equiv_B_N,char **equiv_B,int equiv_A_N,char **equiv_A,int equiv_D_N,char **equiv_D,
                      char ***general_bond_types_registry,char ***general_angle_types_registry,char ***general_dihedral_types_registry);

void SYBYLtoUFF(int atoms_core,char **atom_type_core,char **atom_type_core_FF);

void UFF_topo(int general_atoms,int general_atom_types,
              int general_bonds,int *intermed_bond_types,
              int general_angles,int *intermed_angle_types,
              int general_dihedrals,int *intermed_dihedral_types,
              int general_impropers,int *intermed_improper_types,
              char **master_species_array,char **atom_type_core_FF,
              char **species_intermed_array,
              char ****global_bond_type_intermed_array,char ****global_angle_type_intermed_array,char ****global_dihedral_type_intermed_array,char ****global_improper_type_intermed_array,
              int **general_bonds_registry,int **general_angles_registry,int **general_dihedrals_registry,int **general_impropers_registry,
              int entries_UFF,char **species_UFF_array,
              double *D_UFF_array,double *x_UFF_array,double *r0_UFF_array,double *chi_UFF_array,double *Zstar_UFF_array,double *theta0_UFF_array,double *Utor_UFF_array,double *Vtor_UFF_array,
              double ***nb_FF,double ***bonds_FF,double ***angles_FF,double ***dihedrals_FF,double ***impropers_FF,
              int **b_refine_array, int **a_refine_array, int **d_refine_array, int **i_refine_array);

void UFF_write_lammps(char *current_folder,int counter,int l,
                      int general_atoms,int general_bonds,int general_angles,int general_dihedrals,int general_impropers,
                      int general_atom_types,int intermed_bond_types,int intermed_angle_types,int intermed_improper_types,
                      char **atom_type_core_FF,char **species_intermed_array,int *D_multi,
                      double xlo,double xhi,double ylo,double yhi,double zlo,double zhi,
                      int *master_id_array,int *master_mol_array,double *master_q_array,double *master_x_array,double *master_y_array,double *master_z_array,
                      int **general_bonds_registry,int **general_angles_registry,int **general_dihedrals_registry,int **general_impropers_registry,
                      char ***global_bond_type_intermed_array,char ***global_angle_type_intermed_array,char ***global_dihedral_type_intermed_array,char ***global_improper_type_intermed_array,
                      double **nb_FF,double **bonds_FF,double **angles_FF,double **dihedrals_FF,double **impropers_FF,
                      int **b_refine_array, int **a_refine_array, int **d_refine_array, int **i_refine_array);

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
              int FF)
{
    FILE *fp;
    
    char **atom_type_core_FF;
    int j,k;
    char **species_intermed_array;
    int intermed_bond_types,intermed_angle_types,intermed_dihedral_types,intermed_improper_types;
    char ***global_bond_type_intermed_array,***global_angle_type_intermed_array,***global_dihedral_type_intermed_array,***global_improper_type_intermed_array;
    double **bonds_FF,**angles_FF,**dihedrals_FF,**impropers_FF,**nb_FF;            // matrices to hold force field parameters
    int *D_multi;
    int *b_refine_array,*a_refine_array,*d_refine_array,*i_refine_array;
    
    char command[cmax_length],file_path[cmax_length],buffer[cmax_length],word[cmax_length];
    double d_buffer;
    int int_buffer;
    int nx,ny,nz;
    double time_elapsed;
    
    char **equiv_B,**equiv_A,**equiv_D;
    int equiv_B_N,equiv_A_N,equiv_D_N;
    
    // UFF
    if(FF==0)
    {
        // the program invokes the appropriate atom type identification and conversion
        // routine that reads char **master_species_array and stores in char **atom_type_core_FF
        // the force field compatible atom types
        atom_type_core_FF=(char**)malloc(general_atoms*sizeof(char*));
        for(j=0;j<general_atoms;++j)atom_type_core_FF[j]=(char*)malloc(6*sizeof(char));
        //
        SYBYLtoUFF(general_atoms,master_species_array,atom_type_core_FF);
        
        // here we create a family of intermediate arrays to store atom, bond, angle, dihedral and improper types
        // the reason we use these intermediate arrays is because the force field topology routine OVERWRITES and
        // alters the SIZE of the type arrays
        species_intermed_array=(char**)malloc(general_atom_types*sizeof(char*));
        for(j=0;j<general_atom_types;++j)species_intermed_array[j]=(char*)malloc(sub_length*sizeof(char));
        for(j=0;j<general_atom_types;++j)sprintf(species_intermed_array[j],"%s",general_species_registry[j]);
        intermed_bond_types=general_bond_types;
        intermed_angle_types=general_angle_types;
        intermed_dihedral_types=general_dihedral_types;
        intermed_improper_types=general_improper_types;
        global_bond_type_intermed_array=(char***)malloc(intermed_bond_types*sizeof(char**));
        for(j=0;j<intermed_bond_types;++j)global_bond_type_intermed_array[j]=(char**)malloc(3*sizeof(char*));
        for(j=0;j<intermed_bond_types;++j)for(k=0;k<3;++k)global_bond_type_intermed_array[j][k]=(char*)malloc(sub_length*sizeof(char));
        for(j=0;j<intermed_bond_types;++j)for(k=0;k<3;++k)sprintf(global_bond_type_intermed_array[j][k],"%s",general_bond_types_registry[j][k]);
        global_angle_type_intermed_array=(char***)malloc(intermed_angle_types*sizeof(char**));
        for(j=0;j<intermed_angle_types;++j)global_angle_type_intermed_array[j]=(char**)malloc(5*sizeof(char*));
        for(j=0;j<intermed_angle_types;++j)for(k=0;k<5;++k)global_angle_type_intermed_array[j][k]=(char*)malloc(sub_length*sizeof(char));
        for(j=0;j<intermed_angle_types;++j)for(k=0;k<5;++k)sprintf(global_angle_type_intermed_array[j][k],"%s",general_angle_types_registry[j][k]);
        global_dihedral_type_intermed_array=(char***)malloc(intermed_dihedral_types*sizeof(char**));
        for(j=0;j<intermed_dihedral_types;++j)global_dihedral_type_intermed_array[j]=(char**)malloc(7*sizeof(char*));
        for(j=0;j<intermed_dihedral_types;++j)for(k=0;k<7;++k)global_dihedral_type_intermed_array[j][k]=(char*)malloc(sub_length*sizeof(char));
        for(j=0;j<intermed_dihedral_types;++j)for(k=0;k<7;++k)sprintf(global_dihedral_type_intermed_array[j][k],"%s",general_dihedral_types_registry[j][k]);
        global_improper_type_intermed_array=(char***)malloc(intermed_improper_types*sizeof(char**));
        for(j=0;j<intermed_improper_types;++j)global_improper_type_intermed_array[j]=(char**)malloc(4*sizeof(char*));
        for(j=0;j<intermed_improper_types;++j)for(k=0;k<4;++k)global_improper_type_intermed_array[j][k]=(char*)malloc(sub_length*sizeof(char));
        for(j=0;j<intermed_improper_types;++j)for(k=0;k<4;++k)sprintf(global_improper_type_intermed_array[j][k],"%s",general_improper_types_registry[j][k]);

        //
        UFF_topo(general_atoms,general_atom_types,
                 general_bonds,&intermed_bond_types,
                 general_angles,&intermed_angle_types,
                 general_dihedrals,&intermed_dihedral_types,
                 general_impropers,&intermed_improper_types,
                 master_species_array,atom_type_core_FF,
                 species_intermed_array,
                 &global_bond_type_intermed_array,&global_angle_type_intermed_array,&global_dihedral_type_intermed_array,&global_improper_type_intermed_array,
                 general_bonds_registry,general_angles_registry,general_dihedrals_registry,general_impropers_registry,
                 entries_UFF,species_UFF_array,
                 D_UFF_array,x_UFF_array,r0_UFF_array,chi_UFF_array,Zstar_UFF_array,theta0_UFF_array,Utor_UFF_array,Vtor_UFF_array,
                 &nb_FF,&bonds_FF,&angles_FF,&dihedrals_FF,&impropers_FF,
                 &b_refine_array,&a_refine_array,&d_refine_array,&i_refine_array);

        D_multi=(int*)malloc(general_dihedrals*sizeof(int));
        for(j=0;j<general_dihedrals;++j)
        {
            D_multi[j]=(general_neighbors_registry[general_dihedrals_registry[j][2]-1]-1)*(general_neighbors_registry[general_dihedrals_registry[j][3]-1]-1);
        }
        UFF_write_lammps(current_folder,counter,l,
                         general_atoms,general_bonds,general_angles,general_dihedrals,general_impropers,
                         general_atom_types,intermed_bond_types,intermed_angle_types,intermed_improper_types,
                         atom_type_core_FF,species_intermed_array,D_multi,
                         xlo,xhi,ylo,yhi,zlo,zhi,
                         master_id_array,master_mol_array,master_q_array,master_x_array,master_y_array,master_z_array,
                         general_bonds_registry,general_angles_registry,general_dihedrals_registry,general_impropers_registry,
                         global_bond_type_intermed_array,global_angle_type_intermed_array,global_dihedral_type_intermed_array,global_improper_type_intermed_array,
                         nb_FF,bonds_FF,angles_FF,dihedrals_FF,impropers_FF,
                         &b_refine_array,&a_refine_array,&d_refine_array,&i_refine_array);
        free(b_refine_array);free(a_refine_array);free(d_refine_array);if(general_improper_types>0)free(i_refine_array);

        for(j=0;j<general_atoms;++j)free(atom_type_core_FF[j]);free(atom_type_core_FF);
        for(j=0;j<intermed_bond_types;++j)free(bonds_FF[j]);free(bonds_FF);
        for(j=0;j<intermed_angle_types;++j)free(angles_FF[j]);free(angles_FF);
        for(j=0;j<intermed_dihedral_types;++j)free(dihedrals_FF[j]);free(dihedrals_FF);
        if(general_impropers>0){
            for(j=0;j<intermed_improper_types;++j)free(impropers_FF[j]);free(impropers_FF);}
        for(j=0;j<general_atom_types;++j)free(nb_FF[j]);free(nb_FF);
        free(D_multi);
        for(j=0;j<general_atom_types;++j)free(species_intermed_array[j]);free(species_intermed_array);
        for(k=0;k<intermed_bond_types;++k)
            for(j=0;j<3;++j)
                free(global_bond_type_intermed_array[k][j]);
        for(k=0;k<intermed_bond_types;++k)free(global_bond_type_intermed_array[k]);
        free(global_bond_type_intermed_array);
        for(k=0;k<intermed_angle_types;++k)
            for(j=0;j<5;++j)
                free(global_angle_type_intermed_array[k][j]);
        for(k=0;k<intermed_angle_types;++k)free(global_angle_type_intermed_array[k]);
        free(global_angle_type_intermed_array);
        for(k=0;k<intermed_dihedral_types;++k)
            for(j=0;j<7;++j)
                free(global_dihedral_type_intermed_array[k][j]);
        for(k=0;k<intermed_dihedral_types;++k)free(global_dihedral_type_intermed_array[k]);
        free(global_dihedral_type_intermed_array);
        for(k=0;k<intermed_improper_types;++k)
            for(j=0;j<4;++j)
                free(global_improper_type_intermed_array[k][j]);
        for(k=0;k<intermed_improper_types;++k)free(global_improper_type_intermed_array[k]);
        free(global_improper_type_intermed_array);
        
    }
    else if(FF==1)
    {
    // DREIDING
    
        atom_type_core_FF=(char**)malloc(general_atoms*sizeof(char*));
        for(j=0;j<general_atoms;++j)atom_type_core_FF[j]=(char*)malloc(6*sizeof(char));
        //
        SYBYLtoDRE(general_atoms,master_species_array,atom_type_core_FF);
        //
        species_intermed_array=(char**)malloc(general_atom_types*sizeof(char*));
        for(j=0;j<general_atom_types;++j)species_intermed_array[j]=(char*)malloc(sub_length*sizeof(char));
        for(j=0;j<general_atom_types;++j)sprintf(species_intermed_array[j],"%s",general_species_registry[j]);
        intermed_bond_types=general_bond_types;
        intermed_angle_types=general_angle_types;
        intermed_dihedral_types=general_dihedral_types;
        intermed_improper_types=general_improper_types;
        global_bond_type_intermed_array=(char***)malloc(intermed_bond_types*sizeof(char**));
        for(j=0;j<intermed_bond_types;++j)global_bond_type_intermed_array[j]=(char**)malloc(3*sizeof(char*));
        for(j=0;j<intermed_bond_types;++j)for(k=0;k<3;++k)global_bond_type_intermed_array[j][k]=(char*)malloc(sub_length*sizeof(char));
        for(j=0;j<intermed_bond_types;++j)for(k=0;k<3;++k)sprintf(global_bond_type_intermed_array[j][k],"%s",general_bond_types_registry[j][k]);
        global_angle_type_intermed_array=(char***)malloc(intermed_angle_types*sizeof(char**));
        for(j=0;j<intermed_angle_types;++j)global_angle_type_intermed_array[j]=(char**)malloc(5*sizeof(char*));
        for(j=0;j<intermed_angle_types;++j)for(k=0;k<5;++k)global_angle_type_intermed_array[j][k]=(char*)malloc(sub_length*sizeof(char));
        for(j=0;j<intermed_angle_types;++j)for(k=0;k<5;++k)sprintf(global_angle_type_intermed_array[j][k],"%s",general_angle_types_registry[j][k]);
        global_dihedral_type_intermed_array=(char***)malloc(intermed_dihedral_types*sizeof(char**));
        for(j=0;j<intermed_dihedral_types;++j)global_dihedral_type_intermed_array[j]=(char**)malloc(7*sizeof(char*));
        for(j=0;j<intermed_dihedral_types;++j)for(k=0;k<7;++k)global_dihedral_type_intermed_array[j][k]=(char*)malloc(sub_length*sizeof(char));
        for(j=0;j<intermed_dihedral_types;++j)for(k=0;k<7;++k)sprintf(global_dihedral_type_intermed_array[j][k],"%s",general_dihedral_types_registry[j][k]);
        global_improper_type_intermed_array=(char***)malloc(intermed_improper_types*sizeof(char**));
        for(j=0;j<intermed_improper_types;++j)global_improper_type_intermed_array[j]=(char**)malloc(4*sizeof(char*));
        for(j=0;j<intermed_improper_types;++j)for(k=0;k<4;++k)global_improper_type_intermed_array[j][k]=(char*)malloc(sub_length*sizeof(char));
        for(j=0;j<intermed_improper_types;++j)for(k=0;k<4;++k)sprintf(global_improper_type_intermed_array[j][k],"%s",general_improper_types_registry[j][k]);
        //
        DRE_topo(general_atoms,general_atom_types,
                 general_bonds,&intermed_bond_types,
                 general_angles,&intermed_angle_types,
                 general_dihedrals,&intermed_dihedral_types,
                 general_impropers,&intermed_improper_types,
                 master_species_array,atom_type_core_FF,
                 species_intermed_array,
                 &global_bond_type_intermed_array,&global_angle_type_intermed_array,&global_dihedral_type_intermed_array,&global_improper_type_intermed_array,
                 general_bonds_registry,general_angles_registry,general_dihedrals_registry,general_impropers_registry,
                 entries_DRE,species_DRE_array,
                 R0_DRE_array,theta0_DRE_array,Rvdw0_DRE_array,D0_DRE_array,
                 &nb_FF,&bonds_FF,&angles_FF,&dihedrals_FF,&impropers_FF,
                 &b_refine_array,&a_refine_array,&d_refine_array,&i_refine_array,
                 &equiv_B_N,&equiv_B,&equiv_A_N,&equiv_A,&equiv_D_N,&equiv_D);
       
        //--------------------------------------------------------------------------
        // multiplicity calculation
        //
        // Mayo1990:
        // The VjK is taken as the total barrier after adding all possible I and L terms to the energy expression,
        // but the energy is renormalized by the total number of terms having a common J and K.
        //                   -----------------------------------------------------------------
        // Rappe1992:
        // For a given central J-K bond, all torsions about this bond are considered, with each torsional barrier
        // being divided by the number of torsions present about this J-K bond.
        //       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        //

        D_multi=(int*)malloc(general_dihedrals*sizeof(int));
        for(j=0;j<general_dihedrals;++j)
        {
            D_multi[j]=(general_neighbors_registry[general_dihedrals_registry[j][2]-1]-1)*(general_neighbors_registry[general_dihedrals_registry[j][3]-1]-1);
        }

        //
        DRE_write_lammps(current_folder,counter,l,
                         general_atoms,general_bonds,general_angles,general_dihedrals,general_impropers,
                         general_atom_types,intermed_bond_types,intermed_angle_types,intermed_dihedral_types,intermed_improper_types,
                         atom_type_core_FF,species_intermed_array,
                         xlo,xhi,ylo,yhi,zlo,zhi,
                         master_id_array,master_mol_array,master_q_array,master_x_array,master_y_array,master_z_array,
                         general_bonds_registry,general_angles_registry,general_dihedrals_registry,general_impropers_registry,
                         global_bond_type_intermed_array,global_angle_type_intermed_array,global_dihedral_type_intermed_array,global_improper_type_intermed_array,
                         nb_FF,bonds_FF,angles_FF,dihedrals_FF,impropers_FF,
                         &b_refine_array,&a_refine_array,&d_refine_array,&i_refine_array,
                         D_multi,
                         equiv_B_N,equiv_B,equiv_A_N,equiv_A,equiv_D_N,equiv_D,
                         general_bond_types_registry,general_angle_types_registry,general_dihedral_types_registry);
        //
        free(b_refine_array);free(a_refine_array);free(d_refine_array);if(general_improper_types>0)free(i_refine_array);
        for(j=0;j<equiv_B_N;++j)free(equiv_B[j]);free(equiv_B);
        for(j=0;j<equiv_A_N;++j)free(equiv_A[j]);free(equiv_A);
        for(j=0;j<equiv_D_N;++j)free(equiv_D[j]);free(equiv_D);
        //
        
        for(j=0;j<general_atoms;++j)free(atom_type_core_FF[j]);free(atom_type_core_FF);
        for(j=0;j<intermed_bond_types;++j)free(bonds_FF[j]);free(bonds_FF);
        for(j=0;j<intermed_angle_types;++j)free(angles_FF[j]);free(angles_FF);
        for(j=0;j<intermed_dihedral_types;++j)free(dihedrals_FF[j]);free(dihedrals_FF);
        if(general_impropers>0){
            for(j=0;j<intermed_improper_types;++j)free(impropers_FF[j]);free(impropers_FF);}
        for(j=0;j<general_atom_types;++j)free(nb_FF[j]);free(nb_FF);
        free(D_multi);
        for(j=0;j<general_atom_types;++j)free(species_intermed_array[j]);free(species_intermed_array);
        for(k=0;k<intermed_bond_types;++k)
            for(j=0;j<3;++j)
                free(global_bond_type_intermed_array[k][j]);
        for(k=0;k<intermed_bond_types;++k)free(global_bond_type_intermed_array[k]);
        free(global_bond_type_intermed_array);
        for(k=0;k<intermed_angle_types;++k)
            for(j=0;j<5;++j)
                free(global_angle_type_intermed_array[k][j]);
        for(k=0;k<intermed_angle_types;++k)free(global_angle_type_intermed_array[k]);
        free(global_angle_type_intermed_array);
        for(k=0;k<intermed_dihedral_types;++k)
            for(j=0;j<7;++j)
                free(global_dihedral_type_intermed_array[k][j]);
        for(k=0;k<intermed_dihedral_types;++k)free(global_dihedral_type_intermed_array[k]);
        free(global_dihedral_type_intermed_array);
        for(k=0;k<intermed_improper_types;++k)
            for(j=0;j<4;++j)
                free(global_improper_type_intermed_array[k][j]);
        for(k=0;k<intermed_improper_types;++k)free(global_improper_type_intermed_array[k]);
        free(global_improper_type_intermed_array);
    }
    //getchar();
    
    //
    
    if(nproc>0)
    {
        //
        // minimize (+write minimized coords)
        // lammps
        // call lammps to perform the minimization
        //sprintf(command,"mpirun -np %d ./lmp_micro < mol_%d.in > silence.dat",nproc,1);
        //sprintf(command,"/usr/bin/mpirun.mpich2 -np %d ./lmp_micro < mol_%d.in > silence.dat",nproc,1);
        //sprintf(command,"mpirun -np %d ./lmp_openmpi_gnu < mol_%d.in > silence.dat",nproc,1);
        if(FF==0)sprintf(command,"%s -np %d ./%s < mol_%d.in > silence.dat",mpi_from_config,nproc,lmp_from_config,1);
        if(FF==1)sprintf(command,"%s -np %d ./%s < mol_DREIDING_%d.in > silence.dat",mpi_from_config,nproc,lmp_from_config,1);
        system(command);
        // read energy from log file
        sprintf(file_path,"%s/log.lammps",current_folder);
        fp=fopen(file_path,"r");
        if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
        while(fgets(buffer,cmax_length,fp)!=NULL)
        {
            sscanf(buffer,"%s",word);
            if(strcmp(word,"Step")==0)break;
        }
        while(fgets(buffer,cmax_length,fp)!=NULL)
        {
            sscanf(buffer,"%s",word);
            if(strcmp(word,"Loop")==0){sscanf(buffer,"Loop time of %lf",&time_elapsed);break;}
            sscanf(buffer,"%d\t%lf",&int_buffer,&d_buffer);
        }
        fclose(fp);
        // console out
        //printf("$$$ minimization loop time: %lf\t|\tenergy: %lf\n",time_elapsed,d_buffer);
        
        // read coords from lammps snapshot file!!!
        sprintf(file_path,"%s/snapshot_%d.dat",current_folder,1);
        fp=fopen(file_path,"r");
        if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
        fgets(buffer,cmax_length,fp);
        fgets(buffer,cmax_length,fp);
        fgets(buffer,cmax_length,fp);
        fgets(buffer,cmax_length,fp);
        fgets(buffer,cmax_length,fp);
        fgets(buffer,cmax_length,fp);
        fgets(buffer,cmax_length,fp);
        fgets(buffer,cmax_length,fp);
        fgets(buffer,cmax_length,fp);
        for(j=0;j<general_atoms;++j)
        {
            fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%d\t%d\t%d",&int_buffer,&int_buffer,&int_buffer,&d_buffer,&master_x_array[j],&master_y_array[j],&master_z_array[j],&nx,&ny,&nz);
            master_x_array[j]=master_x_array[j]+nx*(xhi-xlo);
            master_y_array[j]=master_y_array[j]+ny*(yhi-ylo);
            master_z_array[j]=master_z_array[j]+nz*(zhi-zlo);
        }
        fclose(fp);
    }
}

