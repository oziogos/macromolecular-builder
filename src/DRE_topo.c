
#define kappa_bonds_ref_DRE     700.0   // kcal/(mol*Ang^2)
#define delta_DRE               0.01    // Ang
#define kappa_angles_ref_DRE    100.0   // kcal/(mol*rad^2)

#include"builder.h"

int equal(double A,double B);

void find_unique_d(int N, int M, double **array,int **out);

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
              int *equiv_B_N,char ***equiv_B,int *equiv_A_N,char ***equiv_A,int *equiv_D_N,char ***equiv_D)
{
    char word[cmax_length];
    int j,n,igl,jgl,kgl,lgl,int_buffer;
    int M,N;
    double BO;
    int **bonds_out,**angles_out,**dihedrals_out,**impropers_out;
    char **refine_buffer;
    //
    double kappa,w,n0,phi0,theta_n,C0;
    int col16,exocyclic;
    double dih_d,dih_n,omega;
    //
    //--------------------------------------------------------------------------
    // apply conversion to species
    // loop on species_intermed_array (holds SYBYL atom types) using j and on
    // master_species_array using n in order to overwrite species_intermed_array
    // with DREIDING compatible types
    for(j=0;j<general_atom_types;++j)
    {
        for(n=0;n<general_atoms;++n)
        {
            if(strcmp(species_intermed_array[j],master_species_array[n])==0)
            {
                sprintf(species_intermed_array[j],"%s",atom_type_core_FF[n]);
                break;
            }
        }
    }
    //--------------------------------------------------------------------------
    // resolve nb parameters
    // we use general_atom_types unaltered because DREIDING utilizes a "1-1" mapping
    // for the atoms and hybridization states that are used in the building code
    // nb_FF will hold DREIDING D0 and Rvdw0 parameters
    *nb_FF=(double**)malloc(general_atom_types*sizeof(double*));
    for(j=0;j<general_atom_types;++j)(*nb_FF)[j]=(double*)malloc(2*sizeof(double));
    for(j=0;j<general_atom_types;++j)
    {
        for(n=0;n<entries_DRE;++n)
        {
            // identify FF parameters
            if(strcmp(species_intermed_array[j],species_DRE_array[n])==0)
            {
                (*nb_FF)[j][0]=D0_DRE_array[n];
                (*nb_FF)[j][1]=Rvdw0_DRE_array[n];
                break;
            }
        }
    }
    //--------------------------------------------------------------------------
    
    //
    // 01. Overwrite global_bond_type_intermed_array according to the FF type mapping.
    // 02. Preallocate bonds_FF: intermed_bond_types x 2; will hold bond FF parameters.
    // 03. Loop on converted global_bond_type_intermed_array and seek out matching positions
    //     igl, jgl in species_DRE_array. Position indices igl, jgl point to appropriate FF
    //     parameters array (R0_DRE_array).
    // 04. Resolve bond order.
    // 05. Calculate and store FF parameters.
    // 06. Use find_unique_d() to resolve multiple occurrencies in bonds_FF. Creates bonds_out.
    // 07. b_refine_array is introduced. The use of this array is crucial! The first column
    //     of general_bonds_registry holds the TOPOLOGICAL unique type of each bond. This
    //     information is also kept in the first column of bonds_out. When we will strip
    //     multiple occurrencies down below in order to keep only unique FF entries, we end
    //     up loosing entries from the first row that will be needed for the type mapping!!
    //     This is why we store this information in b_refine_array. Inside it we place the
    //     entries of the second column of bonds out in order to have a "1-1" mapping from
    //     topological types (probably with force field multiplicities) to unique FF types!
    // 08. Store in refine_buffer all essential FF information. bonds_out is kept for debugging
    //     purposes.
    // 09. Recalculate intermed_bond_types after refinement in order to count unique entries.
    // 10. Flush global_bond_type_intermed_array, bonds_FF and bonds_out.
    // 11. Allocate global_bond_type_intermed_array, bonds_FF and bonds_out using the updated
    //     value of intermed_bond_types.
    // 12. Read final data from refine_buffer.
    //
    
    // apply conversion to bonds
    for(j=0;j<*intermed_bond_types;++j)
    {
        for(n=0;n<general_atoms;++n)
        {
            if(strcmp((*global_bond_type_intermed_array)[j][0],master_species_array[n])==0)
            {
                sprintf((*global_bond_type_intermed_array)[j][0],"%s",atom_type_core_FF[n]);
                break;
            }
        }
        for(n=0;n<general_atoms;++n)
        {
            if(strcmp((*global_bond_type_intermed_array)[j][2],master_species_array[n])==0)
            {
                sprintf((*global_bond_type_intermed_array)[j][2],"%s",atom_type_core_FF[n]);
                break;
            }
        }
    }
    // now the bond type array global_bond_type_intermed_array is populated with
    // DREIDING compatible types
    // work on bonds:
    N=*intermed_bond_types;
    M=2;
    *bonds_FF=(double**)malloc(N*sizeof(double*));
    for(j=0;j<N;++j)(*bonds_FF)[j]=(double*)malloc(M*sizeof(double));
    //
    for(j=0;j<N;++j)
    {
        for(igl=0;igl<entries_DRE;++igl)
        {
            if(strcmp((*global_bond_type_intermed_array)[j][0],species_DRE_array[igl])==0)
                break;
        }
        for(jgl=0;jgl<entries_DRE;++jgl)
        {
            if(strcmp((*global_bond_type_intermed_array)[j][2],species_DRE_array[jgl])==0)
                break;
        }
        // resolve bond order
        //BO=1.0;//printf("--------- %s %s\n",species_DRE_array[igl],species_DRE_array[jgl]);
        kappa=kappa_bonds_ref_DRE;
        //if((strcmp(species_UFF_array[igl],"C_R")==0 && strcmp(species_UFF_array[jgl],"N_R")==0)||(strcmp(species_UFF_array[jgl],"C_R")==0 && strcmp(species_UFF_array[igl],"N_R")==0))
        //    BO=1.41;
        //else if(species_UFF_array[igl][2]=='R' && species_UFF_array[jgl][2]=='R')
        /*
         if(species_DRE_array[igl][2]=='R' && species_DRE_array[jgl][2]=='R'){
         BO=1.50;kappa=BO*kappa_bonds_ref_DRE;}
         else if(species_DRE_array[igl][2]=='1' && species_DRE_array[jgl][2]=='1'){
         BO=3.00;kappa=BO*kappa_bonds_ref_DRE;}
         else if(species_DRE_array[igl][2]=='2' && species_DRE_array[jgl][2]=='2'){
         BO=2.0;kappa=BO*kappa_bonds_ref_DRE;}
         */
        BO=1.0;
        if(strcmp((*global_bond_type_intermed_array)[j][1],"ar")==0)
            BO=1.5;
        else if(strcmp((*global_bond_type_intermed_array)[j][1],"2")==0)
            BO=2.0;
        else if(strcmp((*global_bond_type_intermed_array)[j][1],"3")==0)
            BO=3.0;
        // Kij
        (*bonds_FF)[j][0]=kappa*BO;
        // rij
        (*bonds_FF)[j][1]=R0_DRE_array[igl]+R0_DRE_array[jgl]-delta_DRE;
    }
    // the matrix bonds_FF holds the parameters for the 1-2 bonding interactions
    // refine bonds
    // here we chech for duplicates
    // the family of find_unique functions receive as input:
    // 1. the number of lines of the matrix to be refined
    // 2. the number of columns
    // 3. the matrix itselt
    // and returns the output always in a ( number_of_lines x 3) matrix
    // the first and second columns of the output matrix are the matching rules;
    // -1 in the third column indicates duplicates
    bonds_out=(int**)malloc(N*sizeof(int*));
    for(j=0;j<N;++j)bonds_out[j]=(int*)malloc(3*sizeof(int));
    find_unique_d(N,M,*bonds_FF,bonds_out);
    // apply type refinement
    /*
     for(j=0;j<general_bonds;++j)
     general_bonds_registry[j][0]=bonds_out[general_bonds_registry[j][0]-1][1];
     */
    // b_refine_array stores the second column of the output of find_unique_d, i.e.
    // the new unique enumeration
    *b_refine_array=(int*)malloc((*intermed_bond_types)*sizeof(int));
    for(j=0;j<*intermed_bond_types;++j)(*b_refine_array)[j]=bonds_out[j][1];
    // the refine_buffer matrix stores the following data in tab separated format:
    // 1. the bond type, e.g. C.3  C3
    // 2. the associated results from find_unique_d (init type id, scaled type id
    //    and -1 multiple occurence flag)
    // 3. the actual force field parameters
    // the use of this buffer array is paramount since the type array, the output
    // from find_unique_d() and the force field parameter array are about to
    // be flushed!
    refine_buffer=(char**)malloc(N*sizeof(char*));
    for(j=0;j<N;++j)refine_buffer[j]=(char*)malloc(cmax_length*sizeof(char));
    for(j=0;j<N;++j)sprintf(refine_buffer[j],"%s\t%s\t%s\t%d\t%d\t%d\t%lf\t%lf\t",
                            (*global_bond_type_intermed_array)[j][0],
                            (*global_bond_type_intermed_array)[j][1],
                            (*global_bond_type_intermed_array)[j][2],
                            bonds_out[j][0],
                            bonds_out[j][1],
                            bonds_out[j][2],
                            (*bonds_FF)[j][0],
                            (*bonds_FF)[j][1]
                            );
    // count unique bonds after refinement
    *intermed_bond_types=0;for(j=0;j<N;++j)if(bonds_out[j][2]!=-1)*intermed_bond_types=*intermed_bond_types+1;
    // flush! (we could have used conditional resizes, but flushing works just as well!)
    for(igl=0;igl<N;++igl)
        for(jgl=0;jgl<3;++jgl)
            free((*global_bond_type_intermed_array)[igl][jgl]);
    for(igl=0;igl<N;++igl)free((*global_bond_type_intermed_array)[igl]);
    free(*global_bond_type_intermed_array);
    for(j=0;j<N;++j)free((*bonds_FF)[j]);free((*bonds_FF));
    for(j=0;j<N;++j)free(bonds_out[j]);free(bonds_out);
    // preallocate using the unique number of bonds
    *global_bond_type_intermed_array=(char***)malloc(*intermed_bond_types*sizeof(char**));
    for(igl=0;igl<*intermed_bond_types;++igl)(*global_bond_type_intermed_array)[igl]=(char**)malloc(3*sizeof(char*));
    for(igl=0;igl<*intermed_bond_types;++igl)
        for(jgl=0;jgl<3;++jgl)
            (*global_bond_type_intermed_array)[igl][jgl]=(char*)malloc(sub_length*sizeof(char));
    bonds_out=(int**)malloc(*intermed_bond_types*sizeof(int*));
    for(j=0;j<*intermed_bond_types;++j)bonds_out[j]=(int*)malloc(3*sizeof(int));
    *bonds_FF=(double**)malloc(*intermed_bond_types*sizeof(double*));
    for(j=0;j<*intermed_bond_types;++j)(*bonds_FF)[j]=(double*)malloc(M*sizeof(double));
    // read the buffer
    // when you find a non -1 value in the occurence flag column, augment the jgl counter and
    // write the appropriate data to memory
    jgl=-1;
    for(igl=0;igl<N;++igl)
    {
        sscanf(refine_buffer[igl],"%s\t%s\t%s\t%d\t%d\t%d",word,word,word,&int_buffer,&int_buffer,&kgl);
        if(kgl!=-1)
        {
            jgl=jgl+1;
            sscanf(refine_buffer[igl],"%s\t%s\t%s\t%d\t%d\t%d\t%lf\t%lf",
                   (*global_bond_type_intermed_array)[jgl][0],
                   (*global_bond_type_intermed_array)[jgl][1],
                   (*global_bond_type_intermed_array)[jgl][2],
                   &bonds_out[jgl][0],
                   &bonds_out[jgl][1],
                   &bonds_out[jgl][2],
                   &(*bonds_FF)[jgl][0],
                   &(*bonds_FF)[jgl][1]
                   );
        }
    }
    
    *equiv_B_N=N;
    *equiv_B=(char**)malloc(N*sizeof(char*));for(j=0;j<N;++j)(*equiv_B)[j]=(char*)malloc(cmax_length*sizeof(char));
    for(j=0;j<N;++j)sprintf((*equiv_B)[j],"%s",refine_buffer[j]);
    
    //printf("DRE:\n");for(igl=0;igl<N;++igl)printf("[%d]\t%s\n",igl+1,refine_buffer[igl]);printf("----------------------\n");//getchar();
    
    // free no longer needed arrays
    for(j=0;j<N;++j)free(refine_buffer[j]);free(refine_buffer);
    for(j=0;j<*intermed_bond_types;++j)free(bonds_out[j]);free(bonds_out);
    
    //--------------------------------------------------------------------------
    
    /// apply conversion to angles
    for(j=0;j<*intermed_angle_types;++j)
    {
        for(n=0;n<general_atoms;++n)
        {
            if(strcmp((*global_angle_type_intermed_array)[j][0],master_species_array[n])==0)
            {
                sprintf((*global_angle_type_intermed_array)[j][0],"%s",atom_type_core_FF[n]);
                break;
            }
        }
        for(n=0;n<general_atoms;++n)
        {
            if(strcmp((*global_angle_type_intermed_array)[j][2],master_species_array[n])==0)
            {
                sprintf((*global_angle_type_intermed_array)[j][2],"%s",atom_type_core_FF[n]);
                break;
            }
        }
        for(n=0;n<general_atoms;++n)
        {
            if(strcmp((*global_angle_type_intermed_array)[j][4],master_species_array[n])==0)
            {
                sprintf((*global_angle_type_intermed_array)[j][4],"%s",atom_type_core_FF[n]);
                break;
            }
        }
        
    }
    // work on angles:
    N=*intermed_angle_types;
    M=2;
    *angles_FF=(double**)malloc(N*sizeof(double*));
    for(j=0;j<N;++j)(*angles_FF)[j]=(double*)malloc(M*sizeof(double));
    //
    for(j=0;j<N;++j)
    {
        for(jgl=0;jgl<entries_DRE;++jgl)
        {
            if(strcmp((*global_angle_type_intermed_array)[j][2],species_DRE_array[jgl])==0)
                break;
        }
        if(equal(theta0_DRE_array[jgl],180.0)==1)
        {
            kappa=kappa_angles_ref_DRE;
            theta_n=180.0;
        }
        else
        {
            C0=sin(theta0_DRE_array[jgl]*pi/180.0);
            C0=C0*C0;
            kappa=kappa_angles_ref_DRE/C0;
            theta_n=theta0_DRE_array[jgl];
        }
        // assign
        (*angles_FF)[j][0]=kappa;
        (*angles_FF)[j][1]=theta_n;
    }
    // refine angles
    angles_out=(int**)malloc(N*sizeof(int*));
    for(j=0;j<N;++j)angles_out[j]=(int*)malloc(3*sizeof(int));
    find_unique_d(N,M,*angles_FF,angles_out);
    /*
     // apply type refinement
     for(j=0;j<general_angles;++j)
     general_angles_registry[j][0]=angles_out[general_angles_registry[j][0]-1][1];
     */
    *a_refine_array=(int*)malloc((*intermed_angle_types)*sizeof(int));
    for(j=0;j<*intermed_angle_types;++j)(*a_refine_array)[j]=angles_out[j][1];
    
    //
    refine_buffer=(char**)malloc(N*sizeof(char*));
    for(j=0;j<N;++j)refine_buffer[j]=(char*)malloc(cmax_length*sizeof(char));
    for(j=0;j<N;++j)sprintf(refine_buffer[j],"%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%lf\t%lf\t",
                            (*global_angle_type_intermed_array)[j][0],
                            (*global_angle_type_intermed_array)[j][1],
                            (*global_angle_type_intermed_array)[j][2],
                            (*global_angle_type_intermed_array)[j][3],
                            (*global_angle_type_intermed_array)[j][4],
                            angles_out[j][0],
                            angles_out[j][1],
                            angles_out[j][2],
                            (*angles_FF)[j][0],
                            (*angles_FF)[j][1]
                            );
    //
    *intermed_angle_types=0;for(j=0;j<N;++j)if(angles_out[j][2]!=-1)*intermed_angle_types=*intermed_angle_types+1;
    for(igl=0;igl<N;++igl)
        for(jgl=0;jgl<5;++jgl)
            free((*global_angle_type_intermed_array)[igl][jgl]);
    for(igl=0;igl<N;++igl)free((*global_angle_type_intermed_array)[igl]);
    free(*global_angle_type_intermed_array);
    for(j=0;j<N;++j)free((*angles_FF)[j]);free(*angles_FF);
    for(j=0;j<N;++j)free(angles_out[j]);free(angles_out);
    *global_angle_type_intermed_array=(char***)malloc(*intermed_angle_types*sizeof(char**));
    for(igl=0;igl<*intermed_angle_types;++igl)(*global_angle_type_intermed_array)[igl]=(char**)malloc(5*sizeof(char*));
    for(igl=0;igl<*intermed_angle_types;++igl)
        for(jgl=0;jgl<5;++jgl)
            (*global_angle_type_intermed_array)[igl][jgl]=(char*)malloc(sub_length*sizeof(char));
    angles_out=(int**)malloc(*intermed_angle_types*sizeof(int*));
    for(j=0;j<*intermed_angle_types;++j)angles_out[j]=(int*)malloc(3*sizeof(int));
    *angles_FF=(double**)malloc(*intermed_angle_types*sizeof(double*));
    for(j=0;j<*intermed_angle_types;++j)(*angles_FF)[j]=(double*)malloc(M*sizeof(double));
    jgl=-1;
    for(igl=0;igl<N;++igl)
    {
        sscanf(refine_buffer[igl],"%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d",word,word,word,word,word,&int_buffer,&int_buffer,&kgl);
        if(kgl!=-1)
        {
            jgl=jgl+1;
            sscanf(refine_buffer[igl],"%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%lf\t%lf",
                   (*global_angle_type_intermed_array)[jgl][0],
                   (*global_angle_type_intermed_array)[jgl][1],
                   (*global_angle_type_intermed_array)[jgl][2],
                   (*global_angle_type_intermed_array)[jgl][3],
                   (*global_angle_type_intermed_array)[jgl][4],
                   &angles_out[jgl][0],
                   &angles_out[jgl][1],
                   &angles_out[jgl][2],
                   &(*angles_FF)[jgl][0],
                   &(*angles_FF)[jgl][1]
                   );
        }
    }

    *equiv_A_N=N;
    *equiv_A=(char**)malloc(N*sizeof(char*));for(j=0;j<N;++j)(*equiv_A)[j]=(char*)malloc(cmax_length*sizeof(char));
    for(j=0;j<N;++j)sprintf((*equiv_A)[j],"%s",refine_buffer[j]);
    
    //printf("DRE:\n");for(igl=0;igl<N;++igl)printf("[%d]\t%s\n",igl+1,refine_buffer[igl]);printf("----------------------\n");//getchar();
    
    //
    for(j=0;j<N;++j)free(refine_buffer[j]);free(refine_buffer);
    //
    N=*intermed_angle_types;
    for(j=0;j<N;++j)free(angles_out[j]);free(angles_out);
    //
    
    //==========================================================================
    // Proper dihedral angles in DREIDING
    //--------------------------------------------------------------------------
    // (a) central atoms: _3; BO=1.0
    //--------------------------------------------------------------------------
    // (b) central atoms: _3 + ( _R || _2 ); BO=1.0
    //     e.g.: acetic acid:
    //      H
    //      |
    //   H--C--C==O
    //      |  |
    //      H  O--H
    //--------------------------------------------------------------------------
    // (c) central atoms: _2; BO=2.0
    //--------------------------------------------------------------------------
    // (d) central atoms: _R; BO=1.5
    //--------------------------------------------------------------------------
    // (e) central atoms: _2 || _R; BO=1.0
    //     e.g.: butadiene
    //      H  H  H  H
    //      |  |  |  |
    //   H--C==C--C==C--H
    //      |  |  |  |
    //      H  H  H  H
    //--------------------------------------------------------------------------
    // (f) exception to (e): exocyclic dihedral single bond between two sp2 atoms
    //     central atoms: _R + ( _R || _2 ); BO=1.0 && ( _R || _2 ):exocyclic
    //     e.g.: phenyl ester or biphenyl
    //
    //      H      H
    //       \    /
    //        C--C    O--H
    //       /    \   |
    //   H--C (Ar) C--C==O
    //       \    /
    //        C--C
    //       /    \
    //      H      H
    //--------------------------------------------------------------------------
    // (g) BO=3.0
    //--------------------------------------------------------------------------
    // (h) central atoms: _3 of col16; BO=1.0
    //--------------------------------------------------------------------------
    // (i) central atoms: _3 of col16 + ( _2 || _R ) !col16; BO=1.0
    //--------------------------------------------------------------------------
    // ## triggered as alternative to (b) ##
    // (j) central atoms: ( _2 || _R ) + _3; BO=1.0
    //     AND: atom linked to sp2 !( _2 || _R )
    //     e.g.: propene
    //   H      H
    //    \    /
    //     C==C   H
    //    /    \ /
    //   H      C
    //         / \
    //        H   H
    //==========================================================================
    //
    // Selection tree:
    /*
     if (bo==1.0)
     {
     if (h_state == (_3,_3))
     {
     [h]--else-->[a]
     }
     else if (h_state == (_3,( _R || _2 )))
     {
     [i]--else-->[j]--else-->[b]
     }
     else if (h_state == (( _R || _2),( _R || _2 )))
     {
     [f]--else-->[e]
     }
     }
     else if (bo==1.5)
     {
     [d]
     }
     else if (bo==2.0)
     {
     [c]
     }
     else if (bo==3.0)
     {
     [g]
     }
     */
    //
    // PMMA dihedral types:
    // [1]	C_3	C_3	C_3	C_3     (a)                 3
    // [2]	C_3	C_3	C_3	H_      (a)                 1,4
    // [3]	C_3	C_3	C_2	O_2     (b)                 12
    // [4]	C_3	C_3	C_2	O_3     (j)                 7
    // [5]	C_3	C_3	C_3	C_2     (a)                 5,11
    // [6]	C_3	C_2	O_3	C_3     (i)                 8
    // [7]	C_2	C_3	C_3	H_      (a)                 2,6
    // [8]	C_2	O_3	C_3	H_      (a)                 10
    // [9]	O_2	C_2	O_3	C_3     (i)                 9
    
    //--------------------------------------------------------------------------
    
    /// apply conversion to dihedrals
    for(j=0;j<*intermed_dihedral_types;++j)
    {
        for(n=0;n<general_atoms;++n)
        {
            if(strcmp((*global_dihedral_type_intermed_array)[j][0],master_species_array[n])==0)
            {
                sprintf((*global_dihedral_type_intermed_array)[j][0],"%s",atom_type_core_FF[n]);
                break;
            }
        }
        for(n=0;n<general_atoms;++n)
        {
            if(strcmp((*global_dihedral_type_intermed_array)[j][2],master_species_array[n])==0)
            {
                sprintf((*global_dihedral_type_intermed_array)[j][2],"%s",atom_type_core_FF[n]);
                break;
            }
        }
        for(n=0;n<general_atoms;++n)
        {
            if(strcmp((*global_dihedral_type_intermed_array)[j][4],master_species_array[n])==0)
            {
                sprintf((*global_dihedral_type_intermed_array)[j][4],"%s",atom_type_core_FF[n]);
                break;
            }
        }
        for(n=0;n<general_atoms;++n)
        {
            if(strcmp((*global_dihedral_type_intermed_array)[j][6],master_species_array[n])==0)
            {
                sprintf((*global_dihedral_type_intermed_array)[j][6],"%s",atom_type_core_FF[n]);
                break;
            }
        }
    }
    
    // work on dihedrals:
    N=*intermed_dihedral_types;
    M=4;    // K,n0,phi0,w
    *dihedrals_FF=(double**)malloc(N*sizeof(double*));
    for(j=0;j<N;++j)(*dihedrals_FF)[j]=(double*)malloc(M*sizeof(double));
    //
    for(j=0;j<N;++j)
    {
        
        kappa=0.0;
        //n0=1.0;
        //phi0=0.0;
        dih_d=1.0;
        dih_n=0.0;
        w=1.0;
        
        for(igl=0;igl<entries_DRE;++igl)
        {
            if(strcmp((*global_dihedral_type_intermed_array)[j][0],species_DRE_array[igl])==0)
                break;
        }
        for(jgl=0;jgl<entries_DRE;++jgl)
        {
            if(strcmp((*global_dihedral_type_intermed_array)[j][2],species_DRE_array[jgl])==0)
                break;
        }
        for(kgl=0;kgl<entries_DRE;++kgl)
        {
            if(strcmp((*global_dihedral_type_intermed_array)[j][4],species_DRE_array[kgl])==0)
                break;
        }
        for(lgl=0;lgl<entries_DRE;++lgl)
        {
            if(strcmp((*global_dihedral_type_intermed_array)[j][6],species_DRE_array[lgl])==0)
                break;
        }
        
        
        //printf("$ (%s)--%s--(%s)--%s--(%s)--%s--(%s)\n",(*global_dihedral_type_intermed_array)[j][0],(*global_dihedral_type_intermed_array)[j][1],(*global_dihedral_type_intermed_array)[j][2],(*global_dihedral_type_intermed_array)[j][3],(*global_dihedral_type_intermed_array)[j][4],(*global_dihedral_type_intermed_array)[j][5],(*global_dihedral_type_intermed_array)[j][6]);
        
        //printf("** [%d]\t%s-%s-%s-%s\n",j+1,species_DRE_array[igl],species_DRE_array[jgl],species_DRE_array[kgl],species_DRE_array[lgl]);
        //printf("[%c][%c]\n",species_DRE_array[jgl][2],species_DRE_array[kgl][2]);
        
        // bond order refinement
        // 1.0
        if(strcmp((*global_dihedral_type_intermed_array)[j][3],"1")==0)
        {
            // look for '_3'-'_3' central bond
            if(species_DRE_array[jgl][2]=='3' && species_DRE_array[kgl][2]=='3')
            {
                // check for col16 first!!
                col16=0;
                if(((species_DRE_array[jgl][0]=='O' && species_DRE_array[jgl][1]=='_') || (species_DRE_array[jgl][0]=='S' && species_DRE_array[jgl][1]=='_'))
                   &&
                   ((species_DRE_array[kgl][0]=='O' && species_DRE_array[kgl][1]=='_') || (species_DRE_array[kgl][0]=='S' && species_DRE_array[kgl][1]=='_'))
                   )
                {
                    col16=1;
                }
                if(col16==1)
                {
                    // [h]
                    kappa=1.0;
                    //n0=2.0;
                    //phi0=0.0;
                    dih_d=1.0;
                    dih_n=2.0;
                    w=1.0;
                    
                    //printf("$ [%d->h] (%s)--%s--(%s)--%s--(%s)--%s--(%s)\n",j+1,(*global_dihedral_type_intermed_array)[j][0],(*global_dihedral_type_intermed_array)[j][1],(*global_dihedral_type_intermed_array)[j][2],(*global_dihedral_type_intermed_array)[j][3],(*global_dihedral_type_intermed_array)[j][4],(*global_dihedral_type_intermed_array)[j][5],(*global_dihedral_type_intermed_array)[j][6]);
                    
                    
                }
                else
                {
                    // [a]
                    kappa=1.0;
                    //n0=3.0;
                    //phi0=0.0;
                    dih_d=1.0;
                    dih_n=3.0;
                    w=1.0;
                    
                    //printf("$ [%d->a] (%s)--%s--(%s)--%s--(%s)--%s--(%s)\n",j+1,(*global_dihedral_type_intermed_array)[j][0],(*global_dihedral_type_intermed_array)[j][1],(*global_dihedral_type_intermed_array)[j][2],(*global_dihedral_type_intermed_array)[j][3],(*global_dihedral_type_intermed_array)[j][4],(*global_dihedral_type_intermed_array)[j][5],(*global_dihedral_type_intermed_array)[j][6]);
                    
                }
            }
            //
            else if((species_DRE_array[jgl][2]=='R' || species_DRE_array[jgl][2]=='2') && (species_DRE_array[kgl][2]=='R' || species_DRE_array[kgl][2]=='2'))
            {
                exocyclic=0;
                // ADD EXOCYCLIC!!!
                // ...
                
                if(exocyclic==1)
                {
                    void;
                }
                else
                {
                    // [e]
                    kappa=2.5;
                    //n0=2.0;
                    //phi0=180.0;
                    dih_d=-1.0;
                    dih_n=2.0;
                    w=1.0;
                    
                    //printf("$ [%d->e] (%s)--%s--(%s)--%s--(%s)--%s--(%s)\n",j+1,(*global_dihedral_type_intermed_array)[j][0],(*global_dihedral_type_intermed_array)[j][1],(*global_dihedral_type_intermed_array)[j][2],(*global_dihedral_type_intermed_array)[j][3],(*global_dihedral_type_intermed_array)[j][4],(*global_dihedral_type_intermed_array)[j][5],(*global_dihedral_type_intermed_array)[j][6]);
                    
                    
                }
            }
            //
            // plus reverse!!!
            else if(species_DRE_array[jgl][2]=='3' && (species_DRE_array[kgl][2]=='R' || species_DRE_array[kgl][2]=='2'))
            {
                
                
                col16=0;
                if(((species_DRE_array[jgl][0]=='O' && species_DRE_array[jgl][1]=='_') || (species_DRE_array[jgl][0]=='S' && species_DRE_array[jgl][1]=='_'))
                   &&
                   !((species_DRE_array[kgl][0]=='O' && species_DRE_array[kgl][1]=='_') || (species_DRE_array[kgl][0]=='S' && species_DRE_array[kgl][1]=='_')))
                {
                    col16=1;
                }
                if(col16==1)
                {
                    // [i]
                    kappa=1.0;
                    //n0=2.0;
                    //phi0=180.0;
                    dih_d=-1.0;
                    dih_n=2.0;
                    w=1.0;
                    
                    //printf("$ [%d->i] (%s)--%s--(%s)--%s--(%s)--%s--(%s)\n",j+1,(*global_dihedral_type_intermed_array)[j][0],(*global_dihedral_type_intermed_array)[j][1],(*global_dihedral_type_intermed_array)[j][2],(*global_dihedral_type_intermed_array)[j][3],(*global_dihedral_type_intermed_array)[j][4],(*global_dihedral_type_intermed_array)[j][5],(*global_dihedral_type_intermed_array)[j][6]);
                    
                }
                else
                {
                    if(!(species_DRE_array[lgl][2]=='R' || species_DRE_array[lgl][2]=='2'))
                    {
                        // [j]
                        kappa=1.0;
                        //n0=3.0;
                        //phi0=0.0;
                        dih_d=1.0;
                        dih_n=3.0;
                        w=1.0;
                        
                        //printf("$ [%d->j] (%s)--%s--(%s)--%s--(%s)--%s--(%s)\n",j+1,(*global_dihedral_type_intermed_array)[j][0],(*global_dihedral_type_intermed_array)[j][1],(*global_dihedral_type_intermed_array)[j][2],(*global_dihedral_type_intermed_array)[j][3],(*global_dihedral_type_intermed_array)[j][4],(*global_dihedral_type_intermed_array)[j][5],(*global_dihedral_type_intermed_array)[j][6]);
                        
                    }
                    else
                    {
                        // [b]
                        kappa=0.5;
                        //n0=6.0;
                        //phi0=180.0;
                        dih_d=-1.0;
                        dih_n=6.0;
                        w=1.0;
                        
                        //printf("$ [%d->b] (%s)--%s--(%s)--%s--(%s)--%s--(%s)\n",j+1,(*global_dihedral_type_intermed_array)[j][0],(*global_dihedral_type_intermed_array)[j][1],(*global_dihedral_type_intermed_array)[j][2],(*global_dihedral_type_intermed_array)[j][3],(*global_dihedral_type_intermed_array)[j][4],(*global_dihedral_type_intermed_array)[j][5],(*global_dihedral_type_intermed_array)[j][6]);
                        
                        
                    }
                }
            }
            // reverse!
            else if(species_DRE_array[kgl][2]=='3' && (species_DRE_array[jgl][2]=='R' || species_DRE_array[jgl][2]=='2'))
            {
                
                
                
                col16=0;
                if(((species_DRE_array[kgl][0]=='O' && species_DRE_array[kgl][1]=='_') || (species_DRE_array[kgl][0]=='S' && species_DRE_array[kgl][1]=='_'))
                   &&
                   !((species_DRE_array[jgl][0]=='O' && species_DRE_array[jgl][1]=='_') || (species_DRE_array[jgl][0]=='S' && species_DRE_array[jgl][1]=='_')))
                {
                    col16=1;
                }
                if(col16==1)
                {
                    // [i]
                    kappa=1.0;
                    //n0=2.0;
                    //phi0=180.0;
                    dih_d=-1.0;
                    dih_n=2.0;
                    w=1.0;
                    
                    //printf("$ [%d->i]REV (%s)--%s--(%s)--%s--(%s)--%s--(%s)\n",j+1,(*global_dihedral_type_intermed_array)[j][0],(*global_dihedral_type_intermed_array)[j][1],(*global_dihedral_type_intermed_array)[j][2],(*global_dihedral_type_intermed_array)[j][3],(*global_dihedral_type_intermed_array)[j][4],(*global_dihedral_type_intermed_array)[j][5],(*global_dihedral_type_intermed_array)[j][6]);
                    
                    
                }
                else
                {
                    if(!(species_DRE_array[igl][2]=='R' || species_DRE_array[igl][2]=='2'))
                    {
                        // [j]
                        kappa=1.0;
                        //n0=3.0;
                        //phi0=0.0;
                        dih_d=1.0;
                        dih_n=3.0;
                        w=1.0;
                        
                        //printf("$ [%d->j]REV (%s)--%s--(%s)--%s--(%s)--%s--(%s)\n",j+1,(*global_dihedral_type_intermed_array)[j][0],(*global_dihedral_type_intermed_array)[j][1],(*global_dihedral_type_intermed_array)[j][2],(*global_dihedral_type_intermed_array)[j][3],(*global_dihedral_type_intermed_array)[j][4],(*global_dihedral_type_intermed_array)[j][5],(*global_dihedral_type_intermed_array)[j][6]);
                        
                        
                    }
                    else
                    {
                        // [b]
                        kappa=0.5;
                        //n0=6.0;
                        //phi0=180.0;
                        dih_d=-1.0;
                        dih_n=6.0;
                        w=1.0;
                        
                        //printf("$ [%d->b]REV (%s)--%s--(%s)--%s--(%s)--%s--(%s)\n",j+1,(*global_dihedral_type_intermed_array)[j][0],(*global_dihedral_type_intermed_array)[j][1],(*global_dihedral_type_intermed_array)[j][2],(*global_dihedral_type_intermed_array)[j][3],(*global_dihedral_type_intermed_array)[j][4],(*global_dihedral_type_intermed_array)[j][5],(*global_dihedral_type_intermed_array)[j][6]);
                        
                        
                    }
                }
            }
        }
        // ar
        else if(strcmp((*global_dihedral_type_intermed_array)[j][3],"ar")==0)
        {
            // [d]
            // no need to check atom types; 'ar' always connects '_R' types!
            kappa=12.5;
            //n0=2.0;
            //phi0=180.0;
            dih_d=-1.0;
            dih_n=2.0;
            w=1.0;
            
            //printf("$ [%d->d] (%s)--%s--(%s)--%s--(%s)--%s--(%s)\n",j+1,(*global_dihedral_type_intermed_array)[j][0],(*global_dihedral_type_intermed_array)[j][1],(*global_dihedral_type_intermed_array)[j][2],(*global_dihedral_type_intermed_array)[j][3],(*global_dihedral_type_intermed_array)[j][4],(*global_dihedral_type_intermed_array)[j][5],(*global_dihedral_type_intermed_array)[j][6]);
            
            
        }
        // 2.0
        else if(strcmp((*global_dihedral_type_intermed_array)[j][3],"2")==0)
        {
            // [c]
            // no need to check atom types; '2' always connects '_2' types!
            kappa=22.5;
            //n0=2.0;
            //phi0=180.0;
            dih_d=-1.0;
            dih_n=2.0;
            w=1.0;
            
            //printf("$ [%d->c] (%s)--%s--(%s)--%s--(%s)--%s--(%s)\n",j+1,(*global_dihedral_type_intermed_array)[j][0],(*global_dihedral_type_intermed_array)[j][1],(*global_dihedral_type_intermed_array)[j][2],(*global_dihedral_type_intermed_array)[j][3],(*global_dihedral_type_intermed_array)[j][4],(*global_dihedral_type_intermed_array)[j][5],(*global_dihedral_type_intermed_array)[j][6]);
            
            
        }
        // 3.0
        else
        {
            // [g]
            kappa=0.0;
            //n0=1.0;
            //phi0=0.0;
            dih_d=1.0;
            dih_n=0.0;
            w=1.0;
            
            //printf("$ [%d->g] (%s)--%s--(%s)--%s--(%s)--%s--(%s)\n",j+1,(*global_dihedral_type_intermed_array)[j][0],(*global_dihedral_type_intermed_array)[j][1],(*global_dihedral_type_intermed_array)[j][2],(*global_dihedral_type_intermed_array)[j][3],(*global_dihedral_type_intermed_array)[j][4],(*global_dihedral_type_intermed_array)[j][5],(*global_dihedral_type_intermed_array)[j][6]);
            
            
        }
        
        // assign
        (*dihedrals_FF)[j][0]=kappa;
        //(*dihedrals_FF)[j][1]=n0;
        //(*dihedrals_FF)[j][2]=phi0;
        (*dihedrals_FF)[j][1]=dih_d;
        (*dihedrals_FF)[j][2]=dih_n;
        (*dihedrals_FF)[j][3]=w;
    }
    
    // refine dihedrals
    dihedrals_out=(int**)malloc(N*sizeof(int*));
    for(j=0;j<N;++j)dihedrals_out[j]=(int*)malloc(3*sizeof(int));
    find_unique_d(N,M,*dihedrals_FF,dihedrals_out);
    
    *d_refine_array=(int*)malloc((*intermed_dihedral_types)*sizeof(int));
    for(j=0;j<*intermed_dihedral_types;++j)(*d_refine_array)[j]=dihedrals_out[j][1];
    
    //
    refine_buffer=(char**)malloc(N*sizeof(char*));
    for(j=0;j<N;++j)refine_buffer[j]=(char*)malloc(cmax_length*sizeof(char));
    for(j=0;j<N;++j)sprintf(refine_buffer[j],"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t",
                            (*global_dihedral_type_intermed_array)[j][0],
                            (*global_dihedral_type_intermed_array)[j][1],
                            (*global_dihedral_type_intermed_array)[j][2],
                            (*global_dihedral_type_intermed_array)[j][3],
                            (*global_dihedral_type_intermed_array)[j][4],
                            (*global_dihedral_type_intermed_array)[j][5],
                            (*global_dihedral_type_intermed_array)[j][6],
                            dihedrals_out[j][0],
                            dihedrals_out[j][1],
                            dihedrals_out[j][2],
                            (*dihedrals_FF)[j][0],
                            (*dihedrals_FF)[j][1],
                            (*dihedrals_FF)[j][2],
                            (*dihedrals_FF)[j][3]
                            );
    
    //
    *intermed_dihedral_types=0;for(j=0;j<N;++j)if(dihedrals_out[j][2]!=-1)*intermed_dihedral_types=*intermed_dihedral_types+1;
    for(igl=0;igl<N;++igl)
        for(jgl=0;jgl<7;++jgl)
            free((*global_dihedral_type_intermed_array)[igl][jgl]);
    for(igl=0;igl<N;++igl)free((*global_dihedral_type_intermed_array)[igl]);
    free(*global_dihedral_type_intermed_array);
    for(j=0;j<N;++j)free((*dihedrals_FF)[j]);free(*dihedrals_FF);
    for(j=0;j<N;++j)free(dihedrals_out[j]);free(dihedrals_out);
    *global_dihedral_type_intermed_array=(char***)malloc(*intermed_dihedral_types*sizeof(char**));
    for(igl=0;igl<*intermed_dihedral_types;++igl)(*global_dihedral_type_intermed_array)[igl]=(char**)malloc(7*sizeof(char*));
    for(igl=0;igl<*intermed_dihedral_types;++igl)
        for(jgl=0;jgl<7;++jgl)
            (*global_dihedral_type_intermed_array)[igl][jgl]=(char*)malloc(sub_length*sizeof(char));
    dihedrals_out=(int**)malloc(*intermed_dihedral_types*sizeof(int*));
    for(j=0;j<*intermed_dihedral_types;++j)dihedrals_out[j]=(int*)malloc(3*sizeof(int));
    *dihedrals_FF=(double**)malloc(*intermed_dihedral_types*sizeof(double*));
    for(j=0;j<*intermed_dihedral_types;++j)(*dihedrals_FF)[j]=(double*)malloc(M*sizeof(double));
    jgl=-1;
    for(igl=0;igl<N;++igl)
    {
        sscanf(refine_buffer[igl],"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d",word,word,word,word,word,word,word,&int_buffer,&int_buffer,&kgl);
        if(kgl!=-1)
        {
            jgl=jgl+1;
            sscanf(refine_buffer[igl],"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf",
                   (*global_dihedral_type_intermed_array)[jgl][0],
                   (*global_dihedral_type_intermed_array)[jgl][1],
                   (*global_dihedral_type_intermed_array)[jgl][2],
                   (*global_dihedral_type_intermed_array)[jgl][3],
                   (*global_dihedral_type_intermed_array)[jgl][4],
                   (*global_dihedral_type_intermed_array)[jgl][5],
                   (*global_dihedral_type_intermed_array)[jgl][6],
                   &dihedrals_out[jgl][0],
                   &dihedrals_out[jgl][1],
                   &dihedrals_out[jgl][2],
                   &(*dihedrals_FF)[jgl][0],
                   &(*dihedrals_FF)[jgl][1],
                   &(*dihedrals_FF)[jgl][2],
                   &(*dihedrals_FF)[jgl][3]
                   );
        }
    }
    
    *equiv_D_N=N;
    *equiv_D=(char**)malloc(N*sizeof(char*));for(j=0;j<N;++j)(*equiv_D)[j]=(char*)malloc(cmax_length*sizeof(char));
    for(j=0;j<N;++j)sprintf((*equiv_D)[j],"%s",refine_buffer[j]);
    
    //printf("DRE:\n");for(igl=0;igl<N;++igl)printf("[%d]\t%s\n",igl+1,refine_buffer[igl]);printf("----------------------\n");//getchar();
    
    
    //
    for(j=0;j<N;++j)free(refine_buffer[j]);free(refine_buffer);
    //
    N=*intermed_dihedral_types;
    for(j=0;j<N;++j)free(dihedrals_out[j]);free(dihedrals_out);
    //
    
    //--------------------------------------------------------------------------
    // impropers
    if(general_impropers>0)
    {
        /// apply conversion to dihedrals
        for(j=0;j<*intermed_improper_types;++j)
        {
            for(n=0;n<general_atoms;++n)
            {
                if(strcmp((*global_improper_type_intermed_array)[j][0],master_species_array[n])==0)
                {
                    sprintf((*global_improper_type_intermed_array)[j][0],"%s",atom_type_core_FF[n]);
                    break;
                }
            }
            for(n=0;n<general_atoms;++n)
            {
                if(strcmp((*global_improper_type_intermed_array)[j][1],master_species_array[n])==0)
                {
                    sprintf((*global_improper_type_intermed_array)[j][1],"%s",atom_type_core_FF[n]);
                    break;
                }
            }
            for(n=0;n<general_atoms;++n)
            {
                if(strcmp((*global_improper_type_intermed_array)[j][2],master_species_array[n])==0)
                {
                    sprintf((*global_improper_type_intermed_array)[j][2],"%s",atom_type_core_FF[n]);
                    break;
                }
            }
            for(n=0;n<general_atoms;++n)
            {
                if(strcmp((*global_improper_type_intermed_array)[j][3],master_species_array[n])==0)
                {
                    sprintf((*global_improper_type_intermed_array)[j][3],"%s",atom_type_core_FF[n]);
                    break;
                }
            }
            
        }
        //
        // work on impropers
        N=*intermed_improper_types;
        M=2;
        *impropers_FF=(double**)malloc(N*sizeof(double*));
        for(j=0;j<N;++j)(*impropers_FF)[j]=(double*)malloc(M*sizeof(double));
        for(j=0;j<N;++j)
        {
            // initialize to zero
            kappa=0.0;omega=0.0;
            for(igl=0;igl<entries_DRE;++igl)
            {
                if(strcmp((*global_improper_type_intermed_array)[j][0],species_DRE_array[igl])==0)
                    break;
            }
            for(jgl=0;jgl<entries_DRE;++jgl)
            {
                if(strcmp((*global_improper_type_intermed_array)[j][1],species_DRE_array[jgl])==0)
                    break;
            }
            for(kgl=0;kgl<entries_DRE;++kgl)
            {
                if(strcmp((*global_improper_type_intermed_array)[j][2],species_DRE_array[kgl])==0)
                    break;
            }
            for(lgl=0;lgl<entries_DRE;++lgl)
            {
                if(strcmp((*global_improper_type_intermed_array)[j][3],species_DRE_array[lgl])==0)
                    break;
            }
            // resolve parameters
            if((*global_improper_type_intermed_array)[j][0][2]=='2' || (*global_improper_type_intermed_array)[j][0][2]=='R')
            {
                kappa=40.0;
                omega=0.0;
            }
            
            // assign
            (*impropers_FF)[j][0]=kappa;
            (*impropers_FF)[j][1]=omega;
        }
        // refine impropers
        impropers_out=(int**)malloc(N*sizeof(int*));
        for(j=0;j<N;++j)impropers_out[j]=(int*)malloc(3*sizeof(int));
        find_unique_d(N,M,*impropers_FF,impropers_out);
        
        *i_refine_array=(int*)malloc((*intermed_improper_types)*sizeof(int));
        for(j=0;j<*intermed_improper_types;++j)(*i_refine_array)[j]=impropers_out[j][1];
        
        //
        refine_buffer=(char**)malloc(N*sizeof(char*));
        for(j=0;j<N;++j)refine_buffer[j]=(char*)malloc(cmax_length*sizeof(char));
        for(j=0;j<N;++j)sprintf(refine_buffer[j],"%s\t%s\t%s\t%s\t%d\t%d\t%d\t%lf\t%lf\t",
                                (*global_improper_type_intermed_array)[j][0],
                                (*global_improper_type_intermed_array)[j][1],
                                (*global_improper_type_intermed_array)[j][2],
                                (*global_improper_type_intermed_array)[j][3],
                                impropers_out[j][0],
                                impropers_out[j][1],
                                impropers_out[j][2],
                                (*impropers_FF)[j][0],
                                (*impropers_FF)[j][1]
                                );
        //
        *intermed_improper_types=0;for(j=0;j<N;++j)if(impropers_out[j][2]!=-1)*intermed_improper_types=*intermed_improper_types+1;
        for(igl=0;igl<N;++igl)
            for(jgl=0;jgl<4;++jgl)
                free((*global_improper_type_intermed_array)[igl][jgl]);
        for(igl=0;igl<N;++igl)free((*global_improper_type_intermed_array)[igl]);
        free(*global_improper_type_intermed_array);
        for(j=0;j<N;++j)free((*impropers_FF)[j]);free(*impropers_FF);
        for(j=0;j<N;++j)free(impropers_out[j]);free(impropers_out);
        *global_improper_type_intermed_array=(char***)malloc(*intermed_improper_types*sizeof(char**));
        for(igl=0;igl<*intermed_improper_types;++igl)(*global_improper_type_intermed_array)[igl]=(char**)malloc(4*sizeof(char*));
        for(igl=0;igl<*intermed_improper_types;++igl)
            for(jgl=0;jgl<4;++jgl)
                (*global_improper_type_intermed_array)[igl][jgl]=(char*)malloc(sub_length*sizeof(char));
        impropers_out=(int**)malloc(*intermed_improper_types*sizeof(int*));
        for(j=0;j<*intermed_improper_types;++j)impropers_out[j]=(int*)malloc(3*sizeof(int));
        *impropers_FF=(double**)malloc(*intermed_improper_types*sizeof(double*));
        for(j=0;j<*intermed_improper_types;++j)(*impropers_FF)[j]=(double*)malloc(M*sizeof(double));
        jgl=-1;
        for(igl=0;igl<N;++igl)
        {
            sscanf(refine_buffer[igl],"%s\t%s\t%s\t%s\t%d\t%d\t%d",word,word,word,word,&int_buffer,&int_buffer,&kgl);
            if(kgl!=-1)
            {
                jgl=jgl+1;
                sscanf(refine_buffer[igl],"%s\t%s\t%s\t%s\t%d\t%d\t%d\t%lf\t%lf\t",
                       (*global_improper_type_intermed_array)[jgl][0],
                       (*global_improper_type_intermed_array)[jgl][1],
                       (*global_improper_type_intermed_array)[jgl][2],
                       (*global_improper_type_intermed_array)[jgl][3],
                       &impropers_out[jgl][0],
                       &impropers_out[jgl][1],
                       &impropers_out[jgl][2],
                       &(*impropers_FF)[jgl][0],
                       &(*impropers_FF)[jgl][1]
                       );
            }
        }
        
        //printf("DRE:\n");for(igl=0;igl<N;++igl)printf("[%d]\t%s\n",igl+1,refine_buffer[igl]);printf("----------------------\n");//getchar();
        
        //
        for(j=0;j<N;++j)free(refine_buffer[j]);free(refine_buffer);
        //
        N=*intermed_improper_types;
        for(j=0;j<N;++j)free(impropers_out[j]);free(impropers_out);
    }
    
    
}
