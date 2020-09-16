
#include "builder.h"

int equal(double A,double B);
void lb(FILE *fp,int start,int *branches_out,int *depth_out,int ***lb_output_matrix);
void find_unique(int N, int M, int **array,int **out);

int one_two_init(int *B1, int *B2, int atoms, int bonds);
void one_two_build(int **one_two, int *B1, int *B2, int atoms, int bonds, int max_1_2);

void growth_reg_init(int *type, char **core_mol2, int *chains, int mol_types, int total_chains, char **quad_matrix, int *rows_out, int *cols_out, int *species_rows, int *topo_rows,
                     int *master_growth_registry_rows,int *master_growth_registry_cols,int *master_species_rows,int *master_topo_rows,
                     char ***master_species_chain_array,char ***master_B_type_chain_array,int ***master_growth_registry,
                     int **master_B1_chain_array,int **master_B2_chain_array,
                     char ***master_atom_name_chain_array,
                     int *MGR_start, int *MGR_stop, int *MA_start, int *MA_stop, int *MB_start, int *MB_stop,
                     int molecules, int *mol_type_array, int **max_growth_step_array, int **growth_step_array)
{
    FILE *fp;
    char current_folder[cmax_length],file_path[cmax_length],buffer[cmax_length],word[cmax_length];
    int i,j,k,l,m,int_buffer;
    char **atom_name_core,**atom_type_core,**subst_name_core,**B_type_core;
    int atoms_core,bonds_core;
    int *atom_id_core,*subst_id_core,*B_ID_core,*B1_core,*B2_core;
    double *x_core,*y_core,*z_core,*charge_core;
    char **chain_file_name;
    int final_rows,final_columns;
    int host_id,graft_id,bo;
    int *init_host_array,*init_graft_array,*init_bo_array;
    int *chain_atoms_array,*chain_bonds_array;
    int atoms,bonds,branches,depth,unique,columns,rows,max_columns,target,master_counter;
    int **fu_in,**fu_out;
    int **fu2_in,**fu2_out;
    int **sequence_matrix;
    int length,*length_array;
    int **partition_matrix;
    int **growth_registry,grows,gcols;
    int **BACKUP;
    int total_chain_atoms,total_chain_bonds;
    char **species_chain_array,**B_type_chain_array;
    int *B1_chain_array,*B2_chain_array;
    double d_buffer;
    int **lb_output_matrix;
    char **atom_name_chain_array;
    
    int sum,sum1,sum2,n;
    
    //==========================================================================
    //
    // The function is invoked inside a mol_types i-loop. The main function
    // provides:
    // - type_array[i]   -> type (int)
    // - core_file[i]    -> core_mol2 (char*)
    // - chains_array[i] -> chains (int)
    // - total_chains    -> total_chains (int)
    // - quad_matrix     -> quad_matrix (char**)
    //
    // The function opens the core mol2 file and stores all data in memory. Then
    // the quad_matrix array is read and side chain mol2 filename(s), growth and
    // host id(s) and initial bond order(s) corresponding to the defined molecular
    // type are stored.
    //
    // Inside a k-loop over the number of side chains:
    // - The associated side chain mol2 file is opened
    // - The number of atoms and bonds are stored in int *chain_atoms_array and
    //   *chain_bonds_array respectively
    // - Function lb() is invoked. The function processes the side chain mol2
    //   file, having the information of the graft atom id via init_graft_array[k]
    //   and returns int **lb_output_matrix that has int branches rows and int
    //   depth columns. Each row of lb_output_matrix contains the atomic id
    //   sequence that defines a unique path from the graft atom to a side chain
    //   terminal atom.
    
    getcwd(current_folder,cmax_length);
    
    sum=0;sum1=0;sum2=0;
    for(n=0;n<mol_types;++n)
    {
        
        sprintf(file_path,"%s/%s",current_folder,core_mol2[n]);
        
        fp=fopen(file_path,"r");
        if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
        while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"@<TRIPOS>MOLECULE\n")==0)break;          // locate the MOLECULE section
        fgets(buffer,cmax_length,fp);fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d",&atoms_core,&bonds_core);    // skip comment line, read atoms and bonds from the next one
        // preallocate arrays
        // atoms
        atom_id_core=(int*)malloc(atoms_core*sizeof(int));
        subst_id_core=(int*)malloc(atoms_core*sizeof(int));
        x_core=(double*)malloc(atoms_core*sizeof(double));
        y_core=(double*)malloc(atoms_core*sizeof(double));
        z_core=(double*)malloc(atoms_core*sizeof(double));
        charge_core=(double*)malloc(atoms_core*sizeof(double));
        atom_name_core=(char**)malloc(atoms_core*sizeof(char*));
        for(i=0;i<atoms_core;++i)atom_name_core[i]=(char*)malloc(sub_length*sizeof(char));
        atom_type_core=(char**)malloc(atoms_core*sizeof(char*));
        for(i=0;i<atoms_core;++i)atom_type_core[i]=(char*)malloc(sub_length*sizeof(char));
        subst_name_core=(char**)malloc(atoms_core*sizeof(char*));
        for(i=0;i<atoms_core;++i)subst_name_core[i]=(char*)malloc(sub_length*sizeof(char));
        // bonds
        B_ID_core=(int*)malloc(bonds_core*sizeof(int));
        B1_core=(int*)malloc(bonds_core*sizeof(int));
        B2_core=(int*)malloc(bonds_core*sizeof(int));
        B_type_core=(char**)malloc(bonds_core*sizeof(char*));
        for(i=0;i<bonds_core;++i)B_type_core[i]=(char*)malloc(sub_length*sizeof(char));
        // read data
        while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"@<TRIPOS>ATOM\n")==0)break;
        for(i=0;i<atoms_core;++i)
        {
            fgets(buffer,cmax_length,fp);
            sscanf(buffer,"%d\t%s\t%lf\t%lf\t%lf\t%s\t%d\t%s\t%lf",&atom_id_core[i],atom_name_core[i],&x_core[i],&y_core[i],&z_core[i],atom_type_core[i],&subst_id_core[i],subst_name_core[i],&charge_core[i]);
        }
        while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"@<TRIPOS>BOND\n")==0)break;
        for(i=0;i<bonds_core;++i)
        {
            fgets(buffer,cmax_length,fp);
            sscanf(buffer,"%d\t%d\t%d\t%s",&B_ID_core[i],&B1_core[i],&B2_core[i],B_type_core[i]);
        }
        fclose(fp);
        
        //
        
        chain_file_name=(char**)malloc(chains[n]*sizeof(char*));
        for(i=0;i<chains[n];++i)chain_file_name[i]=(char*)malloc(sub_length*sizeof(char));
        init_host_array=(int*)malloc(chains[n]*sizeof(int));
        init_graft_array=(int*)malloc(chains[n]*sizeof(int));
        init_bo_array=(int*)malloc(chains[n]*sizeof(int));
        
        j=-1;
        for(i=0;i<total_chains;++i)
        {
            sscanf(quad_matrix[i],"%d\t%s\t%d\t%d\t%d",&int_buffer,word,&graft_id,&host_id,&bo);
            if(int_buffer==type[n])
            {
                j=j+1;
                sprintf(chain_file_name[j],"%s",word);
                init_graft_array[j]=graft_id;
                init_host_array[j]=host_id;
                init_bo_array[j]=bo;
            }
        }
        
        //
        
        // this code segment stores the number of atoms and bonds in every side chain and determines the number of rows and columns for the sequence matrix
        
        chain_atoms_array=(int*)malloc(chains[n]*sizeof(int));
        chain_bonds_array=(int*)malloc(chains[n]*sizeof(int));
        final_rows=0;final_columns=1;
        
        for(k=0;k<chains[n];++k)
        {
            sprintf(file_path,"%s/%s",current_folder,chain_file_name[k]);
            fp=fopen(file_path,"r");
            if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
            // read atoms and bonds
            while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"@<TRIPOS>MOLECULE\n")==0)break;          // locate the MOLECULE section
            fgets(buffer,cmax_length,fp);fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d",&atoms,&bonds);    // skip comment line, read atoms and bonds from the next one
            chain_atoms_array[k]=atoms;
            chain_bonds_array[k]=bonds;
            rewind(fp);
            lb(fp,init_graft_array[k],&branches,&depth,&lb_output_matrix);
            
            /*
             printf("lb_output_matrix:\n");
             for(j=0;j<branches;++j)
             {
             for(l=0;l<depth;++l)printf("%d\t",lb_output_matrix[j][l]);
             printf("\n");
             }
             printf("*****************\n\n\n\n");
             */
            
            fclose(fp);
            fu_in =(int**)malloc(branches*sizeof(int*));for(i=0;i<branches;++i)fu_in[i] =(int*)malloc(2*sizeof(int));
            fu_out=(int**)malloc(branches*sizeof(int*));for(i=0;i<branches;++i)fu_out[i]=(int*)malloc(3*sizeof(int));
            // count...
            rows=0;max_columns=1;
            for(l=0;l<depth-1;++l)
            {
                for(i=0;i<branches;++i)
                {
                    fu_in[i][0]=lb_output_matrix[i][l];
                    fu_in[i][1]=lb_output_matrix[i][l+1];
                }
                find_unique(branches,2,fu_in,fu_out);
                
                //
                //printf("store dyad of columns and invoke find_unique() -- fu_in | fu_out\n");
                //
                
                unique=0;
                for(i=0;i<branches;++i)
                {
                    if(fu_in[i][0]!=0 && fu_in[i][1]!=0)
                    {
                        if(fu_out[i][2]!=-1)
                            unique=unique+1;
                    }
                    
                    //
                    //printf("$$\t%d\t%d\t|\t%d\t%d\t%d\n",fu_in[i][0],fu_in[i][1],fu_out[i][0],fu_out[i][1],fu_out[i][2]);
                    //
                    
                }
                
                //
                //printf("unique = %d\n",unique);
                //
                
                fu2_in =(int**)malloc(unique*sizeof(int*));for(i=0;i<unique;++i)fu2_in[i] =(int*)malloc(1*sizeof(int));
                fu2_out =(int**)malloc(unique*sizeof(int*));for(i=0;i<unique;++i)fu2_out[i] =(int*)malloc(3*sizeof(int));
                
                //
                //printf("transfer unique entries first column to memory and invoke find_unique() -- fu2_in | fu2_out\n");
                //
                
                j=0;
                for(i=0;i<branches;++i)
                {
                    if(fu_in[i][0]!=0 && fu_in[i][1]!=0)
                    {
                        if(fu_out[i][2]!=-1)
                        {
                            j=j+1;
                            fu2_in[j-1][0]=fu_in[i][0];
                        }
                    }
                }
                find_unique(unique,1,fu2_in,fu2_out);
                for(i=0;i<unique;++i)
                {
                    if(fu2_out[i][2]!=-1)
                    {
                        rows=rows+1;
                        target=fu2_in[fu2_out[i][0]-1][0];
                        columns=0;
                        for(j=0;j<branches;++j)
                        {
                            if(fu_in[j][0]!=0 && fu_in[j][1]!=0)
                            {
                                if(fu_out[j][2]!=-1)
                                {
                                    if(fu_in[j][0]==target)
                                    {
                                        columns=columns+1;
                                    }
                                }
                            }
                        }
                        
                        //
                        //printf("$$\t%d\t|\t%d\t%d\t%d\t",fu2_in[i][0],fu2_out[i][0],fu2_out[i][1],fu2_out[i][2]);
                        //printf("##\t[%d]\tcolumns = %d\n",i+1,columns);
                        //
                        
                        if(columns>max_columns)max_columns=columns;
                        
                    }
                    
                }
                //printf("*****************\n\n\n\n");
                for(i=0;i<unique;++i){free(fu2_in[i]);free(fu2_out[i]);}free(fu2_in);free(fu2_out);
                //if(columns>max_columns)max_columns=columns;
            }
            
            
            // free mem
            for(i=0;i<branches;++i)free(lb_output_matrix[i]);free(lb_output_matrix);    // global!
            for(i=0;i<branches;++i){free(fu_in[i]);free(fu_out[i]);}free(fu_in);free(fu_out);
            // end of counting
            final_rows=final_rows+rows;
            if(max_columns>final_columns)final_columns=max_columns;
            
            
        }
        
        //printf("*****************\tfinal_columns = %d\n\n\n\n",final_columns);
        
        // generate the sequence matrix
        
        sequence_matrix=(int**)malloc(final_rows*sizeof(int*));
        for(i=0;i<final_rows;++i)sequence_matrix[i]=(int*)malloc((final_columns*2+2)*sizeof(int));
        for(i=0;i<final_rows;++i)
            for(j=0;j<(final_columns*2+2);++j)
                sequence_matrix[i][j]=0;
        master_counter=0;
        // k index runs over the side chains
        for(k=0;k<chains[n];++k)
        {
            // read side chain mol2 file and invoke the lb routine
            sprintf(file_path,"%s/%s",current_folder,chain_file_name[k]);
            fp=fopen(file_path,"r");
            if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
            lb(fp,init_graft_array[k],&branches,&depth,&lb_output_matrix);
            fclose(fp);
            // preallocate matrices for the find_unique routine
            fu_in =(int**)malloc(branches*sizeof(int*));for(i=0;i<branches;++i)fu_in[i] =(int*)malloc(2*sizeof(int));
            fu_out=(int**)malloc(branches*sizeof(int*));for(i=0;i<branches;++i)fu_out[i]=(int*)malloc(3*sizeof(int));
            // traverse the lb matrix
            // the l index loops over the columns
            for(l=0;l<depth-1;++l)
            {
                // populate the fu_in matrix and call find_unique
                for(i=0;i<branches;++i)
                {
                    fu_in[i][0]=lb_output_matrix[i][l];
                    fu_in[i][1]=lb_output_matrix[i][l+1];
                }
                find_unique(branches,2,fu_in,fu_out);
                // count unique entries
                unique=0;
                for(i=0;i<branches;++i)
                {
                    if(fu_in[i][0]!=0 && fu_in[i][1]!=0)
                    {
                        if(fu_out[i][2]!=-1)
                            unique=unique+1;
                    }
                }
                // preallocate arrays for the find_unique routine to refine the branching multiplicity
                fu2_in =(int**)malloc(unique*sizeof(int*));for(i=0;i<unique;++i)fu2_in[i] =(int*)malloc(1*sizeof(int));
                fu2_out =(int**)malloc(unique*sizeof(int*));for(i=0;i<unique;++i)fu2_out[i] =(int*)malloc(3*sizeof(int));
                // populate the fu2_in matrix
                j=0;
                for(i=0;i<branches;++i)
                {
                    if(fu_in[i][0]!=0 && fu_in[i][1]!=0)
                    {
                        if(fu_out[i][2]!=-1)
                        {
                            j=j+1;
                            fu2_in[j-1][0]=fu_in[i][0];
                        }
                    }
                }
                // refine fu2_in
                find_unique(unique,1,fu2_in,fu2_out);
                //
                for(i=0;i<unique;++i)
                {
                    m=0;
                    if(fu2_out[i][2]!=-1)
                    {
                        master_counter=master_counter+1;
                        sequence_matrix[master_counter-1][0]=k;
                        target=fu2_in[fu2_out[i][0]-1][0];
                        for(j=0;j<branches;++j)
                        {
                            if(fu_in[j][0]!=0 && fu_in[j][1]!=0)
                            {
                                if(fu_out[j][2]!=-1)
                                {
                                    if(fu_in[j][0]==target)
                                    {
                                        m=m+2;
                                        sequence_matrix[master_counter-1][m]=fu_in[j][1];
                                        sequence_matrix[master_counter-1][m+1]=fu_in[j][0];
                                    }
                                }
                            }
                        }
                    }
                }
                // free mem
                for(i=0;i<unique;++i){free(fu2_in[i]);free(fu2_out[i]);}free(fu2_in);free(fu2_out);
            }
            // free mem
            for(i=0;i<branches;++i)free(lb_output_matrix[i]);free(lb_output_matrix);    // global!
            for(i=0;i<branches;++i){free(fu_in[i]);free(fu_out[i]);}free(fu_in);free(fu_out);
        }
        //printf("** Initial sequence matrix:\n");for(i=0;i<final_rows;++i){for(j=0;j<(final_columns*2+2);++j)printf("%d\t",sequence_matrix[i][j]);printf("\n");}printf("\n");
        for(i=0;i<final_rows;++i)
        {
            for(j=2;j<(final_columns*2+2);++j)
            {
                if(sequence_matrix[i][j]==0)
                {
                    sequence_matrix[i][1]=j/2-1;
                    break;
                }
                if(j==(final_columns*2+2-1))
                {
                    sequence_matrix[i][1]=(j+1)/2-1;
                }
            }
        }
        //printf("** Final sequence matrix:\n");for(i=0;i<final_rows;++i){for(j=0;j<(final_columns*2+2);++j)printf("%d\t",sequence_matrix[i][j]);printf("\n");}printf("\n");
        
        // using the partition matrix to define the shuffling, define the growth_registry matrix
        
        length_array=(int*)malloc(chains[n]*sizeof(int));
        for(i=0;i<chains[n];++i)length_array[i]=0;
        for(i=0;i<final_rows;++i)
        {
            length_array[sequence_matrix[i][0]]=length_array[sequence_matrix[i][0]]+1;
        }
        // mem workaround in order to include molecules that don't have anything growing out of them...
        if(chains[n]>0)
            length=length_array[0];
        else
            length=0;
        //
        for(i=0;i<chains[n];++i)if(length_array[i]>length)length=length_array[i];
        partition_matrix=(int**)malloc(length*sizeof(int*));
        for(i=0;i<length;++i)partition_matrix[i]=(int*)malloc(chains[n]*sizeof(int));
        for(i=0;i<length;++i)for(j=0;j<chains[n];++j)partition_matrix[i][j]=-1;
        for(k=0;k<chains[n];++k)
        {
            j=-1;
            for(i=0;i<final_rows;++i)
            {
                if(sequence_matrix[i][0]==k)
                {
                    j=j+1;
                    partition_matrix[j][k]=i;
                }
            }
        }
        //for(i=0;i<length;++i){for(j=0;j<chains;++j)printf("%d\t",partition_matrix[i][j]);printf("\n");}
        grows=0;
        for(i=0;i<length;++i)
        {
            for(j=0;j<chains[n];++j)
            {
                if(partition_matrix[i][j]!=-1)
                {
                    //for(k=0;k<(final_columns*2+2);++k)
                    //    printf("%d\t",sequence_matrix[partition_matrix[i][j]][k]);
                    //printf("\n");
                    grows=grows+1;
                }
            }
        }
        //.. incorporate bootstrap
        grows=grows+chains[n];
        gcols=final_columns*2+2;
        growth_registry=(int**)malloc(grows*sizeof(int*));
        for(i=0;i<grows;++i)growth_registry[i]=(int*)malloc(gcols*sizeof(int));
        //.. incorporate bootstrap
        l=chains[n]-1;
        for(i=0;i<chains[n];++i)for(j=0;j<gcols;++j){growth_registry[i][j]=0;}
        for(i=0;i<chains[n];++i){growth_registry[i][0]=i;growth_registry[i][1]=1;growth_registry[i][2]=init_graft_array[i];}
        for(i=0;i<length;++i)
        {
            for(j=0;j<chains[n];++j)
            {
                if(partition_matrix[i][j]!=-1)
                {
                    l=l+1;
                    for(k=0;k<(final_columns*2+2);++k)
                        growth_registry[l][k]=sequence_matrix[partition_matrix[i][j]][k];
                }
            }
        }
        //zprintf("\n** Growth registry:\n\n");for(i=0;i<grows;++i){for(j=0;j<gcols;++j)printf("%d\t",growth_registry[i][j]);printf("\n");}printf("\n");
        
        //getchar();
        
        *rows_out=grows;
        *cols_out=gcols;
        
        if(*master_growth_registry_rows==0)
        {
            *master_growth_registry=(int**)malloc(grows*sizeof(int*));
            for(i=0;i<grows;++i)(*master_growth_registry)[i]=(int*)malloc(gcols*sizeof(int));
            *master_growth_registry_rows=grows;
            *master_growth_registry_cols=gcols;
            for(i=0;i<grows;++i)for(j=0;j<gcols;++j)(*master_growth_registry)[i][j]=growth_registry[i][j];
        }
        else
        {
            BACKUP=(int**)malloc(*master_growth_registry_rows*sizeof(int*));
            for(i=0;i<*master_growth_registry_rows;++i)BACKUP[i]=(int*)malloc(*master_growth_registry_cols*sizeof(int));
            for(i=0;i<*master_growth_registry_rows;++i)for(j=0;j<*master_growth_registry_cols;++j)BACKUP[i][j]=(*master_growth_registry)[i][j];
            for(i=0;i<*master_growth_registry_rows;++i)free((*master_growth_registry)[i]);free(*master_growth_registry);
            *master_growth_registry=(int**)malloc((*master_growth_registry_rows+grows)*sizeof(int*));
            if(gcols>*master_growth_registry_cols)
            {
                for(i=0;i<(*master_growth_registry_rows+grows);++i)(*master_growth_registry)[i]=(int*)malloc(gcols*sizeof(int));
                for(i=0;i<(*master_growth_registry_rows+grows);++i)for(j=0;j<gcols;++j)(*master_growth_registry)[i][j]=0;
            }
            else
            {
                for(i=0;i<(*master_growth_registry_rows+grows);++i)(*master_growth_registry)[i]=(int*)malloc(*master_growth_registry_cols*sizeof(int));
                for(i=0;i<(*master_growth_registry_rows+grows);++i)for(j=0;j<*master_growth_registry_cols;++j)(*master_growth_registry)[i][j]=0;
            }
            k=-1;
            for(i=0;i<*master_growth_registry_rows;++i){k=k+1;for(j=0;j<*master_growth_registry_cols;++j){(*master_growth_registry)[k][j]=BACKUP[i][j];}}
            for(i=0;i<grows;++i){k=k+1;for(j=0;j<gcols;++j){(*master_growth_registry)[k][j]=growth_registry[i][j];}}
            for(i=0;i<*master_growth_registry_rows;++i)free(BACKUP[i]);free(BACKUP);
            *master_growth_registry_rows=*master_growth_registry_rows+grows;
            if(gcols>*master_growth_registry_cols)
                *master_growth_registry_cols=gcols;
            
        }
        
        //
        
        // store the atomic species and 1-2 topology and bonding info in contiguous format
        
        total_chain_atoms=0;total_chain_bonds=0;
        for(k=0;k<chains[n];++k)
        {
            total_chain_atoms=total_chain_atoms+chain_atoms_array[k];
            total_chain_bonds=total_chain_bonds+chain_bonds_array[k];
        }
        
        //printf("%d\t%d\n",total_chain_atoms,total_chain_bonds);
        
        species_chain_array=(char**)malloc(total_chain_atoms*sizeof(char*));
        for(i=0;i<total_chain_atoms;++i)species_chain_array[i]=(char*)malloc(sub_length*sizeof(char));
        atom_name_chain_array=(char**)malloc(total_chain_atoms*sizeof(char*));
        for(i=0;i<total_chain_atoms;++i)atom_name_chain_array[i]=(char*)malloc(sub_length*sizeof(char));
        B1_chain_array=(int*)malloc(total_chain_bonds*sizeof(int));
        B2_chain_array=(int*)malloc(total_chain_bonds*sizeof(int));
        B_type_chain_array=(char**)malloc(total_chain_bonds*sizeof(char*));
        for(i=0;i<total_chain_bonds;++i)B_type_chain_array[i]=(char*)malloc(sub_length*sizeof(char));
        
        j=-1;m=-1;
        for(k=0;k<chains[n];++k)
        {
            sprintf(file_path,"%s/%s",current_folder,chain_file_name[k]);
            fp=fopen(file_path,"r");
            if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
            while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"@<TRIPOS>ATOM\n")==0)break;
            for(i=0;i<chain_atoms_array[k];++i)
            {
                j=j+1;
                fgets(buffer,cmax_length,fp);
                sscanf(buffer,"%d\t%s\t%lf\t%lf\t%lf\t%s\t%d\t%s\t%lf",&int_buffer,atom_name_chain_array[j],&d_buffer,&d_buffer,&d_buffer,species_chain_array[j],&int_buffer,word,&d_buffer);
            }
            while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"@<TRIPOS>BOND\n")==0)break;
            for(i=0;i<chain_bonds_array[k];++i)
            {
                m=m+1;
                fgets(buffer,cmax_length,fp);
                sscanf(buffer,"%d\t%d\t%d\t%s",&int_buffer,&B1_chain_array[m],&B2_chain_array[m],B_type_chain_array[m]);
            }
            fclose(fp);
        }
        
        //for(i=0;i<total_chain_atoms;++i)printf("%s\n",species_chain_array[i]);
        //for(i=0;i<total_chain_bonds;++i)printf("%d\t%d\t%s\n",B1_chain_array[i],B2_chain_array[i],B_type_chain_array[i]);
        
        *species_rows=total_chain_atoms;
        *topo_rows=total_chain_bonds;
        
        //    char **master_species_chain_array,**master_B_type_chain_array;
        //    int *master_B1_chain_array,*master_B2_chain_array;
        //    int *master_species_rows,*master_topo_rows;
        
        if(*master_species_rows==0)
        {
            
            *master_species_chain_array=(char**)malloc(total_chain_atoms*sizeof(char*));
            for(i=0;i<total_chain_atoms;++i)(*master_species_chain_array)[i]=(char*)malloc(sub_length*sizeof(char));
            // atom name
            *master_atom_name_chain_array=(char**)malloc(total_chain_atoms*sizeof(char*));
            for(i=0;i<total_chain_atoms;++i)(*master_atom_name_chain_array)[i]=(char*)malloc(sub_length*sizeof(char));
            
            *master_species_rows=total_chain_atoms;
            
            for(i=0;i<total_chain_atoms;++i)sprintf((*master_species_chain_array)[i],"%s",species_chain_array[i]);
            // atom name
            for(i=0;i<total_chain_atoms;++i)sprintf((*master_atom_name_chain_array)[i],"%s",atom_name_chain_array[i]);
            
            //
            
            *master_B_type_chain_array=(char**)malloc(total_chain_bonds*sizeof(char*));
            for(i=0;i<total_chain_bonds;++i)(*master_B_type_chain_array)[i]=(char*)malloc(sub_length*sizeof(char));
            *master_B1_chain_array=(int*)malloc(total_chain_atoms*sizeof(int));
            *master_B2_chain_array=(int*)malloc(total_chain_atoms*sizeof(int));
            
            *master_topo_rows=total_chain_bonds;
            
            for(i=0;i<total_chain_bonds;++i)
            {
                sprintf((*master_B_type_chain_array)[i],"%s",B_type_chain_array[i]);
                (*master_B1_chain_array)[i]=B1_chain_array[i];
                (*master_B2_chain_array)[i]=B2_chain_array[i];
            }
            
        }
        else
        {
            *master_species_chain_array=(char**)realloc(*master_species_chain_array,(*master_species_rows+total_chain_atoms)*sizeof(char*));
            for(i=1;i<=total_chain_atoms;++i)(*master_species_chain_array)[*master_species_rows-1+i]=(char*)malloc(sub_length*sizeof(char));
            // atom name
            *master_atom_name_chain_array=(char**)realloc(*master_atom_name_chain_array,(*master_species_rows+total_chain_atoms)*sizeof(char*));
            for(i=1;i<=total_chain_atoms;++i)(*master_atom_name_chain_array)[*master_species_rows-1+i]=(char*)malloc(sub_length*sizeof(char));
            
            j=-1;
            for(i=1;i<=total_chain_atoms;++i)
            {
                j=j+1;
                sprintf((*master_species_chain_array)[*master_species_rows-1+i],"%s",species_chain_array[j]);
                // atom name
                sprintf((*master_atom_name_chain_array)[*master_species_rows-1+i],"%s",atom_name_chain_array[j]);
            }
            
            *master_species_rows=*master_species_rows+total_chain_atoms;
            
            //
            
            *master_B_type_chain_array=(char**)realloc(*master_B_type_chain_array,(*master_topo_rows+total_chain_bonds)*sizeof(char*));
            for(i=1;i<=total_chain_bonds;++i)(*master_B_type_chain_array)[*master_topo_rows-1+i]=(char*)malloc(sub_length*sizeof(char));
            *master_B1_chain_array=(int*)realloc(*master_B1_chain_array,(*master_topo_rows+total_chain_bonds)*sizeof(int));
            *master_B2_chain_array=(int*)realloc(*master_B2_chain_array,(*master_topo_rows+total_chain_bonds)*sizeof(int));
            
            j=-1;
            for(i=1;i<=total_chain_bonds;++i)
            {
                j=j+1;
                sprintf((*master_B_type_chain_array)[*master_topo_rows-1+i],"%s",B_type_chain_array[j]);
                (*master_B1_chain_array)[*master_topo_rows-1+i]=B1_chain_array[j];
                (*master_B2_chain_array)[*master_topo_rows-1+i]=B2_chain_array[j];
            }
            
            *master_topo_rows=*master_topo_rows+total_chain_bonds;
            
        }
        
        //
        
        // mol2 - core
        free(atom_id_core);free(subst_id_core);free(x_core);free(y_core);free(z_core);free(charge_core);
        for(i=0;i<atoms_core;++i)
        {
            free(atom_name_core[i]);free(atom_type_core[i]);free(subst_name_core[i]);
        }
        free(atom_name_core);free(atom_type_core);free(subst_name_core);
        free(B_ID_core);free(B1_core);free(B2_core);
        for(i=0;i<bonds_core;++i)free(B_type_core[i]);free(B_type_core);
        
        
        free(chain_atoms_array);free(chain_bonds_array);
        for(i=0;i<chains[n];++i)free(chain_file_name[i]);free(chain_file_name);
        free(init_host_array);free(init_graft_array);free(init_bo_array);
        
        for(i=0;i<final_rows;++i)free(sequence_matrix[i]);free(sequence_matrix);
        free(length_array);
        for(i=0;i<length;++i)free(partition_matrix[i]);free(partition_matrix);
        for(i=0;i<grows;++i)free(growth_registry[i]);free(growth_registry);
        
        for(i=0;i<total_chain_atoms;++i)free(species_chain_array[i]);free(species_chain_array);
        free(B1_chain_array);free(B2_chain_array);
        for(i=0;i<total_chain_bonds;++i)free(B_type_chain_array[i]);free(B_type_chain_array);
        
        for(i=0;i<total_chain_atoms;++i)free(atom_name_chain_array[i]);free(atom_name_chain_array);
        
        MGR_start[n]=1+sum;sum=sum+*rows_out;MGR_stop[n]=sum;
        MA_start[n]=1+sum1;sum1=sum1+*species_rows;MA_stop[n]=sum1;
        MB_start[n]=1+sum2;sum2=sum2+*topo_rows;MB_stop[n]=sum2;
        
    }
    // Here we count the number of iterations.
    // Every molecule in assigned with an entry inside the following two arrays:
    *max_growth_step_array=(int*)malloc(molecules*sizeof(int));
    *growth_step_array=(int*)malloc(molecules*sizeof(int));
    //
    // The first array is the number of growth steps that every molecule must complete while
    // the second array is a counter reflecting the current building step per molecule.
    for(i=0;i<molecules;++i)
    {
        // Initialize counter.
        (*growth_step_array)[i]=0;
        // Populate max array using MGR registries per type.
        for(j=0;j<mol_types;++j)
        {
            if(mol_type_array[i]==type[j])
            {
                (*max_growth_step_array)[i]=MGR_stop[j]-MGR_start[j]+1;
                break;
            }
        }
    }
}
