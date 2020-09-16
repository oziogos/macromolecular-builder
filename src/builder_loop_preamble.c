
// --- start the building here ---------------------------------------------
//
// variables used:
//
// general:         int     chains                      the number of chains; read from the driver file
//                                                      -- used in the *host_id_array and *chain_id_array initialization loop
//                          grows                       the number of rows of the growth_registry matrix
//                                                      -- controls the master building loop
//                          *init_bo_array              array containing the initial bond order per bond; read from the driver file
//                                                      -- used only at the first building step
//                          **growth_registry           the MASTER martix that describes the building sequence
//                                                      -- its entries define: host_id, graft_id and the chain index chain_index
//
// core related:    char    **atom_type_core            mol2 core array: species (e.g. 'O.2') - SYBYL atom type
//                          **atom_name_core            mol2 core array: label (e.g. 'O')
//                          **B_type_core               mol2 core array: bond type (e.g. '1', '2', 'ar', ...)
//                          **subst_name_core           mol2 core array: the name of the substructure containing the atom
//                  int     atoms_core                  number of core atoms
//                          bonds_core                  number of core bonds
//                          *B1_core                    mol2 core array: 1-2
//                          *B2_core                    mol2 core array: 1-2
//                          *atom_id_core               mol2 core array: atomic id
//                          *subst_id_core              mol2 core array: the ID number of the substructure containing the atom
//                          *B_ID_core                  mol2 core array: bond id
//                  double  *x_core                     mol2 core array: x coord
//                          *y_core                     mol2 core array: y coord
//                          *z_core                     mol2 core array: z coord
//                          *charge_core                mol2 core array: charge
//
// chain related:   char    **species_chain_array       cumulative side chain species array
//                          **B_type_chain_array        cumulative side chain bond type array
//                  int     *chain_atoms_array          array to hold the number of atoms per side chain
//                          *chain_bonds_array          array to hold the number of bonds per side chain
//                          *B1_chain_array             cumulative side chain 1-2 array
//                          *B2_chain_array             cumulative side chain 1-2 array
//
// -------------------------------------------------------------------------

#include"builder.h"

int one_two_init(int *B1, int *B2, int atoms, int bonds);
void one_two_build(int **one_two, int *B1, int *B2, int atoms, int bonds, int max_1_2);

void permute(int N,int **res);

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
                           int **local_topo_array,int **local_topo_array_B)
{
    char word[cmax_length];
    int j,l,m,n;
    int graft_id,sum,internal_index,bo,GS,diff;
    int counter;
    char ***up_matrix;
    
    int host_id;
    int **one_two_core;
    int factorial;
    int **permute_res;
    
    int k,i,o;
    int one_two_temp_list[4];
    int one_three_temp_list[12];
    int one_four_temp_list[36];
    /*
    int one_two_counter;
    int one_three_counter;
    int one_four_counter;
    int *one_two_list;
    int *one_three_list;
    int *one_four_list;
    */
    
    *one_two_counter=0;
    *one_three_counter=0;
    *one_four_counter=0;
    
    // console out
    //printf("GROWTH INDEX: %d\n",growth_step+1);
    //printf("chain index: %d\n",growth_registry[growth_step][0]);
    // chain_index: store the chain index
    *chain_index=growth_registry[growth_step][0];
    // resolve GS_array
    j=0;
    // resolve host and graft ids
    host_id=growth_registry[growth_step][2*(j+1)+1];
    //printf("host_id: %d\n",host_id);
    // locate the host atom
    for(l=0;l<atoms_core;++l)if(host_id_array[l]==host_id && chain_id_array[l]==*chain_index){*host_index=l;break;}

    // resolve 1-2 neighbors of the core molecule
    *max_1_2_core=one_two_init(B1_core,B2_core,atoms_core,bonds_core);

    one_two_core=(int**)malloc(atoms_core*sizeof(int*));for(l=0;l<atoms_core;++l)one_two_core[l]=(int*)malloc(*max_1_2_core*sizeof(int));
    
    one_two_build(one_two_core,B1_core,B2_core,atoms_core,bonds_core,*max_1_2_core);
    
    //--------------------------------------------------------------------------
    /*
    for(l=0;l<atoms_core;++l)
    {
        printf("[%d]\t",l+1);
        for(m=0;m<*max_1_2_core;++m)if(one_two_core[l][m]!=0)printf("%d\t",one_two_core[l][m]);
        printf("\n");
    }
    */
    for(l=0;l<4;++l)one_two_temp_list[l]=0;
    for(l=0;l<12;++l)one_three_temp_list[l]=0;
    for(l=0;l<36;++l)one_four_temp_list[l]=0;
    
    for(i=0;i<*max_1_2_core;++i)
    {
        one_two_temp_list[i]=one_two_core[*host_index][i];
    }
    /*
    printf("%d 1-2 neighbors:\n",*host_index+1);
    for(i=0;i<4;++i)
        if(one_two_temp_list[i]!=0)printf("%d\n",one_two_temp_list[i]);
    */
    l=-1;
    for(i=0;i<4;++i)
    {
        if(one_two_temp_list[i]!=0)
        {
            for(k=0;k<*max_1_2_core;++k)
            {
                if(one_two_core[one_two_temp_list[i]-1][k]!=0 && one_two_core[one_two_temp_list[i]-1][k]!=*host_index+1)
                {
                    l=l+1;
                    one_three_temp_list[l]=one_two_core[one_two_temp_list[i]-1][k];
                }
            }
        }
    }
    
    for(i=0;i<12-1;++i)
    {
        for(k=i+1;k<12;++k)
        {
            if(one_three_temp_list[k]==one_three_temp_list[i])
            {
                one_three_temp_list[k]=0;
            }
        }
    }
    /*
    printf("%d 1-3 neighbors:\n",*host_index+1);
    for(i=0;i<12;++i)
        if(one_three_temp_list[i]!=0)
            printf("%d\n",one_three_temp_list[i]);
    */
    o=-1;
    for(i=0;i<12;++i)
    {
        if(one_three_temp_list[i]!=0)
        {
            for(k=0;k<*max_1_2_core;++k)
            {
                l=one_two_core[one_three_temp_list[i]-1][k];
                n=0;
                for(m=0;m<4;++m)
                {
                    if(l==one_two_temp_list[m])
                    {
                        n=1;
                        break;
                    }
                }
                if(n==0)
                {
                    o=o+1;
                    one_four_temp_list[o]=l;
                }
            }
        }
    }
    
    for(i=0;i<36-1;++i)
    {
        for(k=i+1;k<36;++k)
        {
            if(one_four_temp_list[k]==one_four_temp_list[i])
            {
                one_four_temp_list[k]=0;
            }
        }
    }
    /*
    printf("%d 1-4 neighbors:\n",*host_index+1);
    for(i=0;i<36;++i)
        if(one_four_temp_list[i]!=0)
            printf("%d\n",one_four_temp_list[i]);
    */
    
    for(i=0;i<4;++i)if(one_two_temp_list[i]!=0)*one_two_counter=*one_two_counter+1;
    for(i=0;i<12;++i)if(one_three_temp_list[i]!=0)*one_three_counter=*one_three_counter+1;
    for(i=0;i<36;++i)if(one_four_temp_list[i]!=0)*one_four_counter=*one_four_counter+1;
    *one_two_list=(int*)malloc(*one_two_counter*sizeof(int));
    *one_three_list=(int*)malloc(*one_three_counter*sizeof(int));
    *one_four_list=(int*)malloc(*one_four_counter*sizeof(int));
    k=-1;
    for(i=0;i<4;++i)if(one_two_temp_list[i]!=0){k=k+1;(*one_two_list)[k]=one_two_temp_list[i];}
    k=-1;
    for(i=0;i<12;++i)if(one_three_temp_list[i]!=0){k=k+1;(*one_three_list)[k]=one_three_temp_list[i];}
    k=-1;
    for(i=0;i<36;++i)if(one_four_temp_list[i]!=0){k=k+1;(*one_four_list)[k]=one_four_temp_list[i];}
    /*
    printf("atom %d\n",*host_index+1);
    printf("1-2:\n");
    for(i=0;i<one_two_counter;++i)printf("%d\n",one_two_list[i]);
    printf("1-3:\n");
    for(i=0;i<one_three_counter;++i)printf("%d\n",one_three_list[i]);
    printf("1-4:\n");
    for(i=0;i<one_four_counter;++i)printf("%d\n",one_four_list[i]);
    */
    //free(one_two_list);free(one_three_list);free(one_four_list);
    //getchar();
    //--------------------------------------------------------------------------
    
    // populate GS_array
    for(l=0;l<MAX_GS;++l)GS_array[l]=0;
    
    
    // populate array with H atomic ids bonded to the host atom
    m=-1;
    for(l=0;l<*max_1_2_core;++l)
    {
        
        if(one_two_core[*host_index][l]!=0 && strcmp(atom_type_core[one_two_core[*host_index][l]-1],"H")==0){m=m+1;GS_array[m]=one_two_core[*host_index][l];}
    }
    // resolve GSs
    *GSs=m+1;
    
    // backup vector
    *one_two_core_vector=(int*)malloc(*max_1_2_core*sizeof(int));
    for(l=0;l<*max_1_2_core;++l)(*one_two_core_vector)[l]=one_two_core[*host_index][l];
    for(l=0;l<atoms_core;++l)free(one_two_core[l]);free(one_two_core);
    //
    // *** at this point, GS_array holds all GSs (their core atomic ids) and GSs holds the number of them
    //
    // we now create GS_reg, a char matrix to hold the following information:
    // # the type of the graft atom (e.g. C.3, C.2, ...)
    // # the bond order (no resonance!!)
    // # the core atomic id of the GS
    // # the graft_id, i.e. the side chain atomic id
    //
    // preallocate and initialize to passivating hydrogens
    *GS_reg=(char**)malloc(*GSs*sizeof(char*));for(l=0;l<*GSs;++l)(*GS_reg)[l]=(char*)malloc(cmax_length*sizeof(char));
    //for(l=0;l<*GSs;++l)sprintf((*GS_reg)[l],"H_1\t%d\t-1",GS_array[l]);
    for(l=0;l<*GSs;++l)sprintf((*GS_reg)[l],"H_1\t%d\t-1\tEND",GS_array[l]);
    
    //
    // populate GS_reg
    m=-1;
    for(j=0;j<growth_registry[growth_step][1];j=j+1)
    {
        // resolve host and graft ids
        host_id=growth_registry[growth_step][2*(j+1)+1];
        graft_id=growth_registry[growth_step][2*(j+1)];
        // resolve the internal index
        sum=0;
        for(l=0;l<=*chain_index-1;++l)sum=sum+chain_atoms_array[l];
        internal_index=sum+graft_id-1;
        // resolve bond order
        if(host_id==0)
        {
            bo=init_bo_array[*chain_index];
        }
        else
        {
            sum=0;
            for(l=0;l<*chain_index;++l)sum=sum+chain_bonds_array[l];
            
            for(l=sum;l<sum+chain_bonds_array[*chain_index];++l)
            {
                if((B1_chain_array[l]==graft_id && B2_chain_array[l]==host_id)||(B2_chain_array[l]==graft_id && B1_chain_array[l]==host_id))
                {
                    // WARNING!! supports only numerical entries!
                    sscanf(B_type_chain_array[l],"%d",&bo);
                    break;
                }
            }
        }
        //
        m=m+1;
        sscanf((*GS_reg)[m],"%s\t%d",word,&GS);
        //sprintf((*GS_reg)[m],"%s_%d\t%d\t%d",species_chain_array[internal_index],bo,GS,graft_id);
        sprintf((*GS_reg)[m],"%s_%d\t%d\t%d\t%s",species_chain_array[internal_index],bo,GS,graft_id,atom_name_chain_array[internal_index]);
        diff=bo-1;
        if(diff>0)
        {
            for(l=0;l<diff;++l)
            {
                m=m+1;
                sscanf((*GS_reg)[m],"%s\t%d",word,&GS);
                //sprintf((*GS_reg)[m],"%s_%d\t%d\t-1","X",0,GS);
                sprintf((*GS_reg)[m],"%s_%d\t%d\t-1\t%s","X",0,GS,"X");
            }
        }
    }
    
    //
    // *** at this point we have a fully populated GS_reg matrix that holds
    //     crucial information regarding the building process
    //     In particular, this registry associates graft atoms with particular
    //     growth sites
    //
    // permutations
    factorial=1;
    for(l=0;l<*GSs;++l)
        factorial=factorial*(l+1);
    permute_res=(int**)malloc(factorial*sizeof(int*));
    for(l=0;l<factorial;++l)permute_res[l]=(int*)malloc(*GSs*sizeof(int));
    permute(*GSs,permute_res);
    
    //
    // *** permute_res is the permutation matrix: factorial is the number of rows
    //     and GSs is the number of columns
    //
    // find unique permutations
    up_matrix=(char***)malloc(factorial*sizeof(char**));
    for(l=0;l<factorial;++l)up_matrix[l]=(char**)malloc(*GSs*sizeof(char*));
    for(l=0;l<factorial;++l)
        for(m=0;m<*GSs;++m)
            up_matrix[l][m]=(char*)malloc(cmax_length*sizeof(char));
    
    for(l=0;l<factorial;++l)
    {
        for(m=0;m<*GSs;++m)
        {
            sscanf((*GS_reg)[permute_res[l][m]-1],"%s",word);
            sprintf(up_matrix[l][m],"%s",word);
        }
    }
    for(l=0;l<factorial-1;++l)
    {
        for(m=l+1;m<factorial;++m)
        {
            counter=0;
            for(n=0;n<*GSs;++n)
            {
                if(strcmp(up_matrix[l][n],up_matrix[m][n])==0)
                    counter=counter+1;
            }
            if(counter==*GSs)
                for(n=0;n<*GSs;++n)
                {
                    sprintf(up_matrix[m][n],"%s","_____");
                    permute_res[m][n]=-1;
                }
        }
    }
    //
    *nup=0;
    for(l=0;l<factorial;++l)if(permute_res[l][0]!=-1)*nup=*nup+1;
    *unique_permutations=(int**)malloc(*nup*sizeof(int*));
    for(l=0;l<*nup;++l)(*unique_permutations)[l]=(int*)malloc(*GSs*sizeof(int));
    n=-1;
    for(l=0;l<factorial;++l)
    {
        if(permute_res[l][0]!=-1)
        {
            n=n+1;
            for(m=0;m<*GSs;++m)
                (*unique_permutations)[n][m]=permute_res[l][m];
        }
    }
    
    //
    
    for(l=0;l<factorial;++l)
        for(m=0;m<*GSs;++m)
            free(up_matrix[l][m]);
    for(l=0;l<factorial;++l)free(up_matrix[l]);
    free(up_matrix);
    
    for(l=0;l<factorial;++l)free(permute_res[l]);free(permute_res);
    
    //
    
    //
    *new_entries_atom_total=0;
    *new_entries_bond_total=0;
    
    // resolve number of new entries
    *new_entries_atom=0;
    *new_entries_bond=0;
    for(l=0;l<*GSs;++l)
    {
        sscanf((*GS_reg)[l],"%s",word);
        if(strcmp(word,"C.3_1")==0 || strcmp(word,"C.2_1")==0){
            *new_entries_atom=*new_entries_atom+3;*new_entries_bond=*new_entries_bond+3;
            *new_entries_atom_total=*new_entries_atom_total+3;*new_entries_bond_total=*new_entries_bond_total+3;
        }
        else if(strcmp(word,"C.3_2")==0 || strcmp(word,"C.2_2")==0){
            *new_entries_atom=*new_entries_atom+2;*new_entries_bond=*new_entries_bond+2;
            *new_entries_atom_total=*new_entries_atom_total+2;*new_entries_bond_total=*new_entries_bond_total+2;
        }
        else if(strcmp(word,"O.3_1")==0){
            *new_entries_atom=*new_entries_atom+1;*new_entries_bond=*new_entries_bond+1;
            *new_entries_atom_total=*new_entries_atom_total+1;*new_entries_bond_total=*new_entries_bond_total+1;
        }
        else if(strcmp(word,"O.2_2")==0){void;}
        else if(strcmp(word,"H_1")==0){void;}
        else if(strcmp(word,"X_0")==0){*new_entries_atom_total=*new_entries_atom_total-1;*new_entries_bond_total=*new_entries_bond_total-1;}
        else if(strcmp(word,"N.3_1")==0){
            *new_entries_atom=*new_entries_atom+2;*new_entries_bond=*new_entries_bond+2;
            *new_entries_atom_total=*new_entries_atom_total+2;*new_entries_bond_total=*new_entries_bond_total+2;
        }
        else{printf("$$\t UNSUPPORTED MOVE!");exit(-1);}
    }
    // atom_id_new, x_new, y_new, z_new, B1_new, B2_new
    // int *atom_id_new,double *x_new,double *y_new,double *z_new,int *B1_new,int *B2_new
    *atom_id_new=(int*)malloc(*new_entries_atom*sizeof(int));
    *x_new=(double*)malloc(*new_entries_atom*sizeof(double));
    *y_new=(double*)malloc(*new_entries_atom*sizeof(double));
    *z_new=(double*)malloc(*new_entries_atom*sizeof(double));
    *B1_new=(int*)malloc(*new_entries_bond*sizeof(int));
    *B2_new=(int*)malloc(*new_entries_bond*sizeof(int));
    
    // backup!!
    *atoms_core_B=atoms_core;
    *bonds_core_B=bonds_core;
    // atoms
    *atom_id_core_B=(int*)malloc(*atoms_core_B*sizeof(int));
    *subst_id_core_B=(int*)malloc(*atoms_core_B*sizeof(int));
    *x_core_B=(double*)malloc(*atoms_core_B*sizeof(double));
    *y_core_B=(double*)malloc(*atoms_core_B*sizeof(double));
    *z_core_B=(double*)malloc(*atoms_core_B*sizeof(double));
    *charge_core_B=(double*)malloc(*atoms_core_B*sizeof(double));
    *atom_name_core_B=(char**)malloc(*atoms_core_B*sizeof(char*));
    for(j=0;j<*atoms_core_B;++j)(*atom_name_core_B)[j]=(char*)malloc(sub_length*sizeof(char));
    *atom_type_core_B=(char**)malloc(*atoms_core_B*sizeof(char*));
    for(j=0;j<*atoms_core_B;++j)(*atom_type_core_B)[j]=(char*)malloc(sub_length*sizeof(char));
    *subst_name_core_B=(char**)malloc(*atoms_core_B*sizeof(char*));
    for(j=0;j<*atoms_core_B;++j)(*subst_name_core_B)[j]=(char*)malloc(sub_length*sizeof(char));
    *host_id_array_B=(int*)malloc(*atoms_core_B*sizeof(int));
    *chain_id_array_B=(int*)malloc(*atoms_core_B*sizeof(int));
    // bonds
    *B_ID_core_B=(int*)malloc(*bonds_core_B*sizeof(int));
    *B1_core_B=(int*)malloc(*bonds_core_B*sizeof(int));
    *B2_core_B=(int*)malloc(*bonds_core_B*sizeof(int));
    *B_type_core_B=(char**)malloc(*bonds_core_B*sizeof(char*));
    for(j=0;j<*bonds_core_B;++j)(*B_type_core_B)[j]=(char*)malloc(sub_length*sizeof(char));
    //
    for(j=0;j<atoms_core;++j)
    {
        (*atom_id_core_B)[j]=atom_id_core[j];
        (*subst_id_core_B)[j]=subst_id_core[j];
        (*x_core_B)[j]=x_core[j];
        (*y_core_B)[j]=y_core[j];
        (*z_core_B)[j]=z_core[j];
        (*charge_core_B)[j]=charge_core[j];
        sprintf((*atom_name_core_B)[j],"%s",atom_name_core[j]);
        sprintf((*atom_type_core_B)[j],"%s",atom_type_core[j]);
        sprintf((*subst_name_core_B)[j],"%s",subst_name_core[j]);
        (*host_id_array_B)[j]=host_id_array[j];
        (*chain_id_array_B)[j]=chain_id_array[j];
    }
    for(j=0;j<bonds_core;++j)
    {
        (*B_ID_core_B)[j]=B_ID_core[j];
        (*B1_core_B)[j]=B1_core[j];
        (*B2_core_B)[j]=B2_core[j];
        sprintf((*B_type_core_B)[j],"%s",B_type_core[j]);
    }
    // end of backup

    //
    // the nup loop
    *current_GS_reg=(char**)malloc(*GSs*sizeof(char*));for(l=0;l<*GSs;++l)(*current_GS_reg)[l]=(char*)malloc(cmax_length*sizeof(char));
    
    //
    *score=(double*)malloc(*nup*sizeof(double));
    
    //
    
    //
    //-- local -------------------------------------------------------------------------------------------------------##
    //
    // As a new addition to builder_loop_preamble(), the number and atomic ids of 1-2, 1-3 and 1-4 neighbors
    // of the host site are identified and stored to memory:
    // one_two_counter,one_two_list,one_three_counter,one_three_list,one_four_counter,one_four_list
    // Note: the ids are internal!
    //
    // Then int* local_topo_array is prealloced, initialized and populated. This array's length is equal to
    // the number of atoms in the molecule (atoms_block_size) and hold the information whether an atom is
    // 1-2, 1-3 or 1-4 neighbor of the host site (value 1) or not (value 0). The host atoms is also given
    // the value 1.
    // Note: the array is backed up ('_B')!
    //
    // New mapping:
    // 1: host
    // 2: 1-2
    // 3: 1-3
    // 4: 1-4
    // -1: GS
    // -2: newH
    //
    // Later on, all newly added atoms by builder_core_funct() are also given the value 1.
    //
    
    // atoms_block_size --> atoms_core
    
    *local_topo_array=(int*)malloc(atoms_core*sizeof(int));
    *local_topo_array_B=(int*)malloc(atoms_core*sizeof(int));
    for(j=0;j<atoms_core;++j)(*local_topo_array)[j]=0;
    (*local_topo_array)[*host_index]=1;
    for(j=0;j<*one_two_counter;++j)(*local_topo_array)[(*one_two_list)[j]-1]=2;
    for(j=0;j<*one_three_counter;++j)(*local_topo_array)[(*one_three_list)[j]-1]=3;
    for(j=0;j<*one_four_counter;++j)(*local_topo_array)[(*one_four_list)[j]-1]=4;
    for(j=0;j<atoms_core;++j)(*local_topo_array_B)[j]=(*local_topo_array)[j];
    //----------------------------------------------------------------------------------------------------------------##
    
}
