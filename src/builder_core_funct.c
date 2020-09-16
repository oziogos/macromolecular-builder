
#include"builder.h"

void grow_ch3(int host_index,int selected_GS,int bonds,int max_1_2,int *atom_id,double *x,double *y,double *z,char **atom_type,char **atom_name,int *B1,int *B2,char **B_type,int *one_two,double *new_x,double *new_y,double *new_z,char *back);
void grow_oh(int host_index,int selected_GS,int bonds,int max_1_2,int *atom_id,double *x,double *y,double *z,char **atom_type,char **atom_name,int *B1,int *B2,char **B_type,int *one_two,double *new_x,double *new_y,double *new_z,char *back);
void grow_ch2(int host_index,int selected_GS,int bonds,int max_1_2,int *atom_id,double *x,double *y,double *z,char **atom_type,char **atom_name,int *B1,int *B2,char **B_type,int *one_two,double *new_x,double *new_y,double *new_z,char *back);
void grow_nh2(int host_index,int selected_GS,int bonds,int max_1_2,int *atom_id,double *x,double *y,double *z,char **atom_type,char **atom_name,int *B1,int *B2,char **B_type,int *one_two,double *new_x,double *new_y,double *new_z,char *back);

// new_entries_atom and new_entries_bond are resolved by builder_loop_preamble() and correspond to the number of newly added atoms and bonds per building step
// atoms_core and bonds_core receive atoms_block_size and bonds_block_size respectively which hold the number of atoms and bonds PRIOR to the building step (resolved by segments())

// atom_id_core stores segment_id_internal_array

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
                        )
//,int **remove_list
{
    
    // local
    char word[cmax_length],species[cmax_length];
    //int j,m,v1,v2,new_index,selected_GS,graft_id,bo,*rescale_array,temp_atoms_core,temp_bonds_core;
    int j,m,v1,v2,new_index,selected_GS,graft_id,bo,temp_atoms_core,temp_bonds_core;
    double *x_from_funct,*y_from_funct,*z_from_funct;
    double d,d0,lamda;//,d_buffer;
    int *temp_atom_id_core,*temp_subst_id_core,*temp_B_ID_core,*temp_B1_core,*temp_B2_core,*temp_host_id_array,*temp_chain_id_array;
    double *temp_x_core,*temp_y_core,*temp_z_core,*temp_charge_core;
    char **temp_atom_name_core,**temp_atom_type_core,**temp_subst_name_core,**temp_B_type_core;
    int *temp_local_topo_array;
    //
    //int added_atoms,*added_atoms_array,n;
    int n;
    //
    char back[cmax_length];
    //
    
    // remove: int *added_atoms,int **added_atoms_array

    // count host_index
    *added_atoms=1;
    
    // read from GS_reg and store in current_GS_reg
    for(m=0;m<GSs;++m)
    {
        //sscanf(GS_reg[unique_permutations[nup_index][m]-1],"%s\t%d\t%d",word,&v1,&v2);
        sscanf(GS_reg[unique_permutations[nup_index][m]-1],"%s\t%d\t%d\t%s",word,&v1,&v2,back);
        //sprintf(current_GS_reg[m],"%s\t%d\t%d",word,GS_array[m],v2);
        sprintf(current_GS_reg[m],"%s\t%d\t%d\t%s",word,GS_array[m],v2,back);
        //
        if( strcmp(word,"C.3_1")==0 || strcmp(word,"C.2_1")==0 )
            *added_atoms=*added_atoms+4;
        if( strcmp(word,"O.3_1")==0 )
            *added_atoms=*added_atoms+2;
        if( strcmp(word,"C.2_2")==0 || strcmp(word,"C.3_2")==0 )
            *added_atoms=*added_atoms+3;
        if( strcmp(word,"O.2_2")==0 )
            *added_atoms=*added_atoms+1;
        // N.3
        if( strcmp(word,"N.3_1")==0 )
            *added_atoms=*added_atoms+3;
    }
    *added_atoms_array=(int*)malloc((*added_atoms)*sizeof(int));
    n=0;
    (*added_atoms_array)[n]=host_index+1;
    
    // console output
    //__console_out__
    if(perc_flag==0){
        printf("\n### building... ########################### [%d/%d]\n",nup_index+1,nup);
        for(m=0;m<GSs;++m){printf("# %s",current_GS_reg[m]);printf("\n");}
        printf("###########################################\n");
    }
    // flush original
    free(*atom_id_core);free(*subst_id_core);free(*x_core);free(*y_core);free(*z_core);free(*charge_core);
    for(j=0;j<*atoms_core;++j){free((*atom_name_core)[j]);free((*atom_type_core)[j]);free((*subst_name_core)[j]);}
    free(*atom_name_core);free(*atom_type_core);free(*subst_name_core);free(*B_ID_core);free(*B1_core);free(*B2_core);
    for(j=0;j<*bonds_core;++j)free((*B_type_core)[j]);free(*B_type_core);free(*host_id_array);free(*chain_id_array);
    free(*local_topo_array);
    
    // prealloc original
    // atoms
    *atom_id_core=(int*)malloc(atoms_core_B*sizeof(int));*subst_id_core=(int*)malloc(atoms_core_B*sizeof(int));
    *x_core=(double*)malloc(atoms_core_B*sizeof(double));*y_core=(double*)malloc(atoms_core_B*sizeof(double));*z_core=(double*)malloc(atoms_core_B*sizeof(double));
    *charge_core=(double*)malloc(atoms_core_B*sizeof(double));
    *atom_name_core=(char**)malloc(atoms_core_B*sizeof(char*));for(j=0;j<atoms_core_B;++j)(*atom_name_core)[j]=(char*)malloc(sub_length*sizeof(char));
    *atom_type_core=(char**)malloc(atoms_core_B*sizeof(char*));for(j=0;j<atoms_core_B;++j)(*atom_type_core)[j]=(char*)malloc(sub_length*sizeof(char));
    *subst_name_core=(char**)malloc(atoms_core_B*sizeof(char*));for(j=0;j<atoms_core_B;++j)(*subst_name_core)[j]=(char*)malloc(sub_length*sizeof(char));
    *host_id_array=(int*)malloc(atoms_core_B*sizeof(int));*chain_id_array=(int*)malloc(atoms_core_B*sizeof(int));
    *local_topo_array=(int*)malloc(atoms_core_B*sizeof(int));
    // bonds
    *B_ID_core=(int*)malloc(bonds_core_B*sizeof(int));*B1_core=(int*)malloc(bonds_core_B*sizeof(int));*B2_core=(int*)malloc(bonds_core_B*sizeof(int));
    *B_type_core=(char**)malloc(bonds_core_B*sizeof(char*));for(j=0;j<bonds_core_B;++j)(*B_type_core)[j]=(char*)malloc(sub_length*sizeof(char));

    // read from backup
    for(j=0;j<atoms_core_B;++j)
    {
        (*atom_id_core)[j]=atom_id_core_B[j];(*subst_id_core)[j]=subst_id_core_B[j];
        (*x_core)[j]=x_core_B[j];(*y_core)[j]=y_core_B[j];(*z_core)[j]=z_core_B[j];(*charge_core)[j]=charge_core_B[j];
        sprintf((*atom_name_core)[j],"%s",atom_name_core_B[j]);
        sprintf((*atom_type_core)[j],"%s",atom_type_core_B[j]);
        sprintf((*subst_name_core)[j],"%s",subst_name_core_B[j]);
        (*host_id_array)[j]=host_id_array_B[j];(*chain_id_array)[j]=chain_id_array_B[j];
        (*local_topo_array)[j]=local_topo_array_B[j];
    }
    for(j=0;j<bonds_core_B;++j)
    {
        (*B_ID_core)[j]=B_ID_core_B[j];(*B1_core)[j]=B1_core_B[j];(*B2_core)[j]=B2_core_B[j];
        sprintf((*B_type_core)[j],"%s",B_type_core_B[j]);
    }
    
    //
    *atoms_core=atoms_core_B;
    *bonds_core=bonds_core_B;
    new_index=-1;
    // read info from current_GS_reg
    for(m=0;m<GSs;++m)
    {
        //
        //sscanf(current_GS_reg[m],"%s\t%d\t%d",word,&selected_GS,&graft_id);
        sscanf(current_GS_reg[m],"%s\t%d\t%d\t%s",word,&selected_GS,&graft_id,back);
        for(j=0;j<cmax_length;++j){if(word[j]=='_')break;species[j]=word[j];}species[j]='\0';
        bo=word[j+1]-'0';
        // grow -ch3
        if((strcmp(species,"C.3")==0 && bo==1)||(strcmp(species,"C.2")==0 && bo==1))
        {
            // info I need to keep:
            // - host_index
            // - selected_GS-1
            // - *atoms_core+new_index+1 x 3
            //
            // added atoms: 1+1+3-1=4
            
            n=n+1;
            (*added_atoms_array)[n]=selected_GS;
            
            //(*local_topo_array)[selected_GS-1]=1;
            (*local_topo_array)[selected_GS-1]=-1;
            (*host_id_array)[selected_GS-1]=graft_id;
            (*chain_id_array)[selected_GS-1]=chain_index;
            //
            x_from_funct=(double*)malloc(3*sizeof(double));y_from_funct=(double*)malloc(3*sizeof(double));z_from_funct=(double*)malloc(3*sizeof(double));
            grow_ch3(host_index,selected_GS,*bonds_core,max_1_2_core,*atom_id_core,*x_core,*y_core,*z_core,*atom_type_core,*atom_name_core,*B1_core,*B2_core,*B_type_core,one_two_core_vector,x_from_funct,y_from_funct,z_from_funct,back);
            //
            for(j=0;j<3;++j)
            {
                new_index=new_index+1;
                atom_id_new[new_index]=*atoms_core+new_index+1;
                x_new[new_index]=x_from_funct[j];y_new[new_index]=y_from_funct[j];z_new[new_index]=z_from_funct[j];
                B1_new[new_index]=selected_GS;B2_new[new_index]=atom_id_new[new_index];
                
                n=n+1;
                (*added_atoms_array)[n]=*atoms_core+new_index+1;
                
            }
            //
            free(x_from_funct);free(y_from_funct);free(z_from_funct);
        }
        // grow -oh
        else if(strcmp(species,"O.3")==0 && bo==1)
        {
            // info I need to keep:
            // - host_index
            // - selected_GS-1
            // - *atoms_core+new_index+1 x 1
            //
            // added atoms: 1+1+1-1=2
            
            n=n+1;
            (*added_atoms_array)[n]=selected_GS;
            
            //(*local_topo_array)[selected_GS-1]=1;
            (*local_topo_array)[selected_GS-1]=-1;
            (*host_id_array)[selected_GS-1]=graft_id;
            (*chain_id_array)[selected_GS-1]=chain_index;
            //
            x_from_funct=(double*)malloc(1*sizeof(double));y_from_funct=(double*)malloc(1*sizeof(double));z_from_funct=(double*)malloc(1*sizeof(double));
            grow_oh(host_index,selected_GS,*bonds_core,max_1_2_core,*atom_id_core,*x_core,*y_core,*z_core,*atom_type_core,*atom_name_core,*B1_core,*B2_core,*B_type_core,one_two_core_vector,x_from_funct,y_from_funct,z_from_funct,back);
            //
            j=0;
            
            new_index=new_index+1;
            atom_id_new[new_index]=*atoms_core+new_index+1;
            x_new[new_index]=x_from_funct[j];y_new[new_index]=y_from_funct[j];z_new[new_index]=z_from_funct[j];
            B1_new[new_index]=selected_GS;B2_new[new_index]=atom_id_new[new_index];
            
            n=n+1;
            (*added_atoms_array)[n]=*atoms_core+new_index+1;
            
            //
            free(x_from_funct);free(y_from_funct);free(z_from_funct);
        }
        // grow =ch2
        else if((strcmp(species,"C.2")==0 && bo==2)||(strcmp(species,"C.3")==0 && bo==2))
        {
            // info I need to keep:
            // - host_index
            // - selected_GS-1
            // - *atoms_core+new_index+1 x 2
            //
            // added atoms: 1+1+2-1=3
            
            n=n+1;
            (*added_atoms_array)[n]=selected_GS;
            
            //(*local_topo_array)[selected_GS-1]=1;
            (*local_topo_array)[selected_GS-1]=-1;
            (*host_id_array)[selected_GS-1]=graft_id;
            (*chain_id_array)[selected_GS-1]=chain_index;
            //
            x_from_funct=(double*)malloc(2*sizeof(double));y_from_funct=(double*)malloc(2*sizeof(double));z_from_funct=(double*)malloc(2*sizeof(double));
            grow_ch2(host_index,selected_GS,*bonds_core,max_1_2_core,*atom_id_core,*x_core,*y_core,*z_core,*atom_type_core,*atom_name_core,*B1_core,*B2_core,*B_type_core,one_two_core_vector,x_from_funct,y_from_funct,z_from_funct,back);
            // *** check the type of the host atom; if it's not X.2 change it!
            sprintf((*atom_type_core)[host_index],"%s","C.2");
            //
            for(j=0;j<2;++j)
            {
                new_index=new_index+1;
                atom_id_new[new_index]=*atoms_core+new_index+1;
                x_new[new_index]=x_from_funct[j];y_new[new_index]=y_from_funct[j];z_new[new_index]=z_from_funct[j];
                B1_new[new_index]=selected_GS;B2_new[new_index]=atom_id_new[new_index];
                
                n=n+1;
                (*added_atoms_array)[n]=*atoms_core+new_index+1;
                
            }
            //
            free(x_from_funct);free(y_from_funct);free(z_from_funct);
        }
        // delete h
        else if(strcmp(species,"X")==0)
        {
            
            //printf("### -- %d -- will be deleted!!\n",(*atom_id_core)[selected_GS-1]);
            //getchar();
            
            (*atom_id_core)[selected_GS-1]=-1;
            for(j=0;j<*bonds_core;++j)
            {
                if((*B1_core)[j]==selected_GS || (*B2_core)[j]==selected_GS)
                {
                    (*B_ID_core)[j]=-1;
                    break;
                }
            }
        }
        // grow =o
        else if(strcmp(species,"O.2")==0 && bo==2)
        {
            // info I need to keep:
            // - host_index
            // - selected_GS-1
            //
            // added atoms: 1+1-1=1
            
            n=n+1;
            (*added_atoms_array)[n]=selected_GS;
            
            // 2. displace the GS along the 1-2 bond in order to achieve the X=O.2 bond length
            // 2.1 resolve target bond length
            if(strcmp((*atom_type_core)[host_index],"C.3")==0)d=d_C_2_O_2;
            else
            {
                printf("@@ ERROR!! unsupported 1-2 pair!!!\n\n");exit(-2);
            }
            // 2.2 displace
            d0=sqrt(((*x_core)[selected_GS-1]-(*x_core)[host_index])*((*x_core)[selected_GS-1]-(*x_core)[host_index])+((*y_core)[selected_GS-1]-(*y_core)[host_index])*((*y_core)[selected_GS-1]-(*y_core)[host_index])+((*z_core)[selected_GS-1]-(*z_core)[host_index])*((*z_core)[selected_GS-1]-(*z_core)[host_index]));
            lamda=d/d0;
            (*x_core)[selected_GS-1]=(*x_core)[host_index]+lamda*((*x_core)[selected_GS-1]-(*x_core)[host_index]);
            (*y_core)[selected_GS-1]=(*y_core)[host_index]+lamda*((*y_core)[selected_GS-1]-(*y_core)[host_index]);
            (*z_core)[selected_GS-1]=(*z_core)[host_index]+lamda*((*z_core)[selected_GS-1]-(*z_core)[host_index]);
            // 3. overwrite atom type: from H to O.2
            sprintf((*atom_type_core)[selected_GS-1],"%s","O.2");
            sprintf((*atom_name_core)[selected_GS-1],"%s","O");
            //
            //(*local_topo_array)[selected_GS-1]=1;
            (*local_topo_array)[selected_GS-1]=-1;
            (*host_id_array)[selected_GS-1]=graft_id;
            (*chain_id_array)[selected_GS-1]=chain_index;
            // 4. overwrite bond order: store '2'
            for(j=0;j<*bonds_core;++j)
            {
                if(((*B1_core)[j]==(*atom_id_core)[host_index] && (*B2_core)[j]==selected_GS)||((*B2_core)[j]==(*atom_id_core)[host_index] && (*B1_core)[j]==selected_GS))
                {
                    sprintf((*B_type_core)[j],"%s","2");
                    break;
                }
            }
            // 5. check the type of the host atom; if it's not X.2 change it!
            sprintf((*atom_type_core)[host_index],"%s","C.2");
        }
        // grow -nh2
        if((strcmp(species,"N.3")==0 && bo==1))
        {
            
            n=n+1;
            (*added_atoms_array)[n]=selected_GS;
            
            (*local_topo_array)[selected_GS-1]=-1;
            (*host_id_array)[selected_GS-1]=graft_id;
            (*chain_id_array)[selected_GS-1]=chain_index;
            //
            x_from_funct=(double*)malloc(2*sizeof(double));y_from_funct=(double*)malloc(2*sizeof(double));z_from_funct=(double*)malloc(2*sizeof(double));
            grow_nh2(host_index,selected_GS,*bonds_core,max_1_2_core,*atom_id_core,*x_core,*y_core,*z_core,*atom_type_core,*atom_name_core,*B1_core,*B2_core,*B_type_core,one_two_core_vector,x_from_funct,y_from_funct,z_from_funct,back);
            //
            for(j=0;j<2;++j)
            {
                new_index=new_index+1;
                atom_id_new[new_index]=*atoms_core+new_index+1;
                x_new[new_index]=x_from_funct[j];y_new[new_index]=y_from_funct[j];z_new[new_index]=z_from_funct[j];
                B1_new[new_index]=selected_GS;B2_new[new_index]=atom_id_new[new_index];
                
                n=n+1;
                (*added_atoms_array)[n]=*atoms_core+new_index+1;
                
            }
            //
            free(x_from_funct);free(y_from_funct);free(z_from_funct);
        }
    }
    
    //for(j=0;j<*added_atoms;++j)printf("[%d]\t%d\n",j+1,(*added_atoms_array)[j]);//getchar();
    
    // rescale_array v2
    *remove_N=0;
    for(j=0;j<*atoms_core;++j)
        if((*atom_id_core)[j]==-1)
            *remove_N=*remove_N+1;
    
    if(*remove_N>0){
    *rescale_array=(int*)malloc(*atoms_core*sizeof(int));
    for(j=0;j<*atoms_core;++j)(*rescale_array)[j]=0;
    m=0;
    for(j=0;j<*atoms_core;++j)
        if((*atom_id_core)[j]!=-1)
        {
            m=m+1;
            (*rescale_array)[j]=m;
        }
    }
    //

    // realloc arrays
    *atom_id_core=(int*)realloc(*atom_id_core,(*atoms_core+new_entries_atom)*sizeof(int));
    *subst_id_core=(int*)realloc(*subst_id_core,(*atoms_core+new_entries_atom)*sizeof(int));
    *x_core=(double*)realloc(*x_core,(*atoms_core+new_entries_atom)*sizeof(double));
    *y_core=(double*)realloc(*y_core,(*atoms_core+new_entries_atom)*sizeof(double));
    *z_core=(double*)realloc(*z_core,(*atoms_core+new_entries_atom)*sizeof(double));
    *charge_core=(double*)realloc(*charge_core,(*atoms_core+new_entries_atom)*sizeof(double));
    *atom_name_core=(char**)realloc(*atom_name_core,(*atoms_core+new_entries_atom)*sizeof(char*));
    for(j=1;j<=new_entries_atom;++j)(*atom_name_core)[*atoms_core-1+j]=(char*)malloc(sub_length*sizeof(char));
    *atom_type_core=(char**)realloc(*atom_type_core,(*atoms_core+new_entries_atom)*sizeof(char*));
    for(j=1;j<=new_entries_atom;++j)(*atom_type_core)[*atoms_core-1+j]=(char*)malloc(sub_length*sizeof(char));
    *subst_name_core=(char**)realloc(*subst_name_core,(*atoms_core+new_entries_atom)*sizeof(char*));
    for(j=1;j<=new_entries_atom;++j)(*subst_name_core)[*atoms_core-1+j]=(char*)malloc(sub_length*sizeof(char));
    *host_id_array=(int*)realloc(*host_id_array,(*atoms_core+new_entries_atom)*sizeof(int));
    *chain_id_array=(int*)realloc(*chain_id_array,(*atoms_core+new_entries_atom)*sizeof(int));
    *local_topo_array=(int*)realloc(*local_topo_array,(*atoms_core+new_entries_atom)*sizeof(int));
    
    // rescale_array v2
    if(*remove_N>0){
        *rescale_array=(int*)realloc(*rescale_array,(*atoms_core+new_entries_atom)*sizeof(int));}
    //
    
    // populate arrays
    for(j=1;j<=new_entries_atom;++j)
    {
        sprintf((*atom_type_core)[*atoms_core-1+j],"%s","H");
        sprintf((*atom_name_core)[*atoms_core-1+j],"%s","H");
        (*atom_id_core)[*atoms_core-1+j]=*atoms_core+j;
        (*charge_core)[*atoms_core-1+j]=0.0;
        sprintf((*subst_name_core)[(*atoms_core)-1+j],"%s",(*subst_name_core)[0]);
        (*subst_id_core)[*atoms_core-1+j]=(*subst_id_core)[0];
        //(*host_id_array)[*atoms_core-1+j]=-1;
        (*host_id_array)[*atoms_core-1+j]=-2;
        (*chain_id_array)[*atoms_core-1+j]=chain_index;
        (*x_core)[*atoms_core-1+j]=x_new[j-1];
        (*y_core)[*atoms_core-1+j]=y_new[j-1];
        (*z_core)[*atoms_core-1+j]=z_new[j-1];
        //(*local_topo_array)[*atoms_core-1+j]=1;
        (*local_topo_array)[*atoms_core-1+j]=-2;
        // rescale_array v2
        if(*remove_N>0){
        m=m+1;
            (*rescale_array)[*atoms_core-1+j]=m;}
        //
    }
    // add new bonds
    // realloc
    *B_ID_core=(int*)realloc(*B_ID_core,(*bonds_core+new_entries_bond)*sizeof(int));
    *B1_core=(int*)realloc(*B1_core,(*bonds_core+new_entries_bond)*sizeof(int));
    *B2_core=(int*)realloc(*B2_core,(*bonds_core+new_entries_bond)*sizeof(int));
    *B_type_core=(char**)realloc(*B_type_core,(*bonds_core+new_entries_bond)*sizeof(char*));
    for(j=1;j<=new_entries_bond;++j)(*B_type_core)[*bonds_core-1+j]=(char*)malloc(sub_length*sizeof(char));
    // populate
    for(j=1;j<=new_entries_bond;++j)
    {
        (*B_ID_core)[*bonds_core-1+j]=*bonds_core+j;
        sprintf((*B_type_core)[*bonds_core-1+j],"%s","1");
        (*B1_core)[*bonds_core-1+j]=B1_new[j-1];
        (*B2_core)[*bonds_core-1+j]=B2_new[j-1];
    }
    *atoms_core=*atoms_core+new_entries_atom;
    *bonds_core=*bonds_core+new_entries_bond;
    
    /*
    *remove_N=0;
    for(j=0;j<*atoms_core;++j)
        if((*atom_id_core)[j]==-1)
            *remove_N=*remove_N+1;
    */
    
    if(*remove_N>0)
    {
        /*
        *remove_list=(int*)malloc(*remove_N*sizeof(int));
        m=-1;
        for(j=0;j<*atoms_core;++j)
            if((*atom_id_core)[j]==-1)
            {
                m=m+1;
                (*remove_list)[m]=j+1;
            }
        */
        
        /*
        rescale_array=(int*)malloc(*atoms_core*sizeof(int));
        for(j=0;j<*atoms_core;++j)rescale_array[j]=0;
        m=0;
        for(j=0;j<*atoms_core;++j)
            if((*atom_id_core)[j]!=-1)
            {
                m=m+1;
                rescale_array[j]=m;
            }
        */
        //
        
        // console out
        //for(j=0;j<*atoms_core;++j)printf("{%d\t%d}\n",(*atom_id_core)[j],(*rescale_array)[j]);
        
        //
        
        temp_atoms_core=*atoms_core-*remove_N;
        temp_bonds_core=*bonds_core-*remove_N;
        // preallocate arrays
        // atoms
        temp_atom_id_core=(int*)malloc(temp_atoms_core*sizeof(int));
        temp_subst_id_core=(int*)malloc(temp_atoms_core*sizeof(int));
        temp_x_core=(double*)malloc(temp_atoms_core*sizeof(double));
        temp_y_core=(double*)malloc(temp_atoms_core*sizeof(double));
        temp_z_core=(double*)malloc(temp_atoms_core*sizeof(double));
        temp_charge_core=(double*)malloc(temp_atoms_core*sizeof(double));
        temp_atom_name_core=(char**)malloc(temp_atoms_core*sizeof(char*));
        for(j=0;j<temp_atoms_core;++j)temp_atom_name_core[j]=(char*)malloc(sub_length*sizeof(char));
        temp_atom_type_core=(char**)malloc(temp_atoms_core*sizeof(char*));
        for(j=0;j<temp_atoms_core;++j)temp_atom_type_core[j]=(char*)malloc(sub_length*sizeof(char));
        temp_subst_name_core=(char**)malloc(temp_atoms_core*sizeof(char*));
        for(j=0;j<temp_atoms_core;++j)temp_subst_name_core[j]=(char*)malloc(sub_length*sizeof(char));
        temp_host_id_array=(int*)malloc(temp_atoms_core*sizeof(int));
        temp_chain_id_array=(int*)malloc(temp_atoms_core*sizeof(int));
        temp_local_topo_array=(int*)malloc(temp_atoms_core*sizeof(int));
        // bonds
        temp_B_ID_core=(int*)malloc(temp_bonds_core*sizeof(int));
        temp_B1_core=(int*)malloc(temp_bonds_core*sizeof(int));
        temp_B2_core=(int*)malloc(temp_bonds_core*sizeof(int));
        temp_B_type_core=(char**)malloc(temp_bonds_core*sizeof(char*));
        for(j=0;j<temp_bonds_core;++j)temp_B_type_core[j]=(char*)malloc(sub_length*sizeof(char));
        //
        m=-1;
        for(j=0;j<*atoms_core;++j)
        {
            if((*atom_id_core)[j]!=-1)
            {
                m=m+1;
                temp_atom_id_core[m]=(*rescale_array)[j];temp_subst_id_core[m]=(*subst_id_core)[j];
                temp_x_core[m]=(*x_core)[j];temp_y_core[m]=(*y_core)[j];temp_z_core[m]=(*z_core)[j];
                temp_charge_core[m]=(*charge_core)[j];
                sprintf(temp_atom_name_core[m],"%s",(*atom_name_core)[j]);
                sprintf(temp_atom_type_core[m],"%s",(*atom_type_core)[j]);
                sprintf(temp_subst_name_core[m],"%s",(*subst_name_core)[j]);
                temp_host_id_array[m]=(*host_id_array)[j];temp_chain_id_array[m]=(*chain_id_array)[j];
                temp_local_topo_array[m]=(*local_topo_array)[j];
            }
        }
        m=-1;
        for(j=0;j<*bonds_core;++j)
        {
            if((*B_ID_core)[j]!=-1)
            {
                m=m+1;
                temp_B_ID_core[m]=m+1;sprintf(temp_B_type_core[m],"%s",(*B_type_core)[j]);
                temp_B1_core[m]=(*rescale_array)[(*B1_core)[j]-1];temp_B2_core[m]=(*rescale_array)[(*B2_core)[j]-1];
            }
        }
        
        //
        free(*atom_id_core);free(*subst_id_core);free(*x_core);free(*y_core);free(*z_core);free(*charge_core);free(*host_id_array);free(*chain_id_array);free(*local_topo_array);
        for(j=0;j<*atoms_core;++j)
        {
            free((*atom_name_core)[j]);free((*atom_type_core)[j]);free((*subst_name_core)[j]);
        }
        free(*atom_name_core);free(*atom_type_core);free(*subst_name_core);
        free(*B_ID_core);free(*B1_core);free(*B2_core);
        for(j=0;j<*bonds_core;++j)free((*B_type_core)[j]);free(*B_type_core);
        *atoms_core=temp_atoms_core;
        *bonds_core=temp_bonds_core;
        // preallocate arrays
        // atoms
        *atom_id_core=(int*)malloc(*atoms_core*sizeof(int));
        *subst_id_core=(int*)malloc(*atoms_core*sizeof(int));
        *x_core=(double*)malloc(*atoms_core*sizeof(double));
        *y_core=(double*)malloc(*atoms_core*sizeof(double));
        *z_core=(double*)malloc(*atoms_core*sizeof(double));
        *charge_core=(double*)malloc(*atoms_core*sizeof(double));
        *atom_name_core=(char**)malloc(*atoms_core*sizeof(char*));
        for(j=0;j<*atoms_core;++j)(*atom_name_core)[j]=(char*)malloc(sub_length*sizeof(char));
        *atom_type_core=(char**)malloc(*atoms_core*sizeof(char*));
        for(j=0;j<*atoms_core;++j)(*atom_type_core)[j]=(char*)malloc(sub_length*sizeof(char));
        *subst_name_core=(char**)malloc(*atoms_core*sizeof(char*));
        for(j=0;j<*atoms_core;++j)(*subst_name_core)[j]=(char*)malloc(sub_length*sizeof(char));
        *host_id_array=(int*)malloc(*atoms_core*sizeof(int));
        *chain_id_array=(int*)malloc(*atoms_core*sizeof(int));
        *local_topo_array=(int*)malloc(*atoms_core*sizeof(int));
        // bonds
        *B_ID_core=(int*)malloc(*bonds_core*sizeof(int));
        *B1_core=(int*)malloc(*bonds_core*sizeof(int));
        *B2_core=(int*)malloc(*bonds_core*sizeof(int));
        *B_type_core=(char**)malloc(*bonds_core*sizeof(char*));
        for(j=0;j<*bonds_core;++j)(*B_type_core)[j]=(char*)malloc(sub_length*sizeof(char));
        for(j=0;j<*atoms_core;++j)
        {
            (*atom_id_core)[j]=temp_atom_id_core[j];(*subst_id_core)[j]=temp_subst_id_core[j];
            (*x_core)[j]=temp_x_core[j];(*y_core)[j]=temp_y_core[j];(*z_core)[j]=temp_z_core[j];
            (*charge_core)[j]=temp_charge_core[j];
            sprintf((*atom_name_core)[j],"%s",temp_atom_name_core[j]);
            sprintf((*atom_type_core)[j],"%s",temp_atom_type_core[j]);
            sprintf((*subst_name_core)[j],"%s",temp_subst_name_core[j]);
            (*host_id_array)[j]=temp_host_id_array[j];(*chain_id_array)[j]=temp_chain_id_array[j];
            (*local_topo_array)[j]=temp_local_topo_array[j];
        }
        for(j=0;j<*bonds_core;++j)
        {
            (*B_ID_core)[j]=temp_B_ID_core[j];(*B1_core)[j]=temp_B1_core[j];(*B2_core)[j]=temp_B2_core[j];sprintf((*B_type_core)[j],"%s",temp_B_type_core[j]);
        }
        //
        //free(rescale_array);
        free(temp_atom_id_core);free(temp_subst_id_core);free(temp_x_core);free(temp_y_core);free(temp_z_core);free(temp_charge_core);free(temp_host_id_array);free(temp_chain_id_array);free(temp_local_topo_array);
        for(j=0;j<temp_atoms_core;++j)
        {
            free(temp_atom_name_core[j]);free(temp_atom_type_core[j]);free(temp_subst_name_core[j]);
        }
        free(temp_atom_name_core);free(temp_atom_type_core);free(temp_subst_name_core);
        free(temp_B_ID_core);free(temp_B1_core);free(temp_B2_core);
        for(j=0;j<temp_bonds_core;++j)free(temp_B_type_core[j]);free(temp_B_type_core);
        
    }
    
    //free(added_atoms_array);
    /*
    if(*remove_N>0){
        printf("[%d]\n",(*atoms_core+new_entries_atom));
        for(j=0;j<(*atoms_core+new_entries_atom);++j)printf("%d\t%d\n",j+1,(*rescale_array)[j]);}
    */
    //for(j=0;j<*bonds_core;++j)printf("[%d]\t%d\t%d\t%s\n",(*B_ID_core)[j],(*B1_core)[j],(*B2_core)[j],(*B_type_core)[j]);
}

