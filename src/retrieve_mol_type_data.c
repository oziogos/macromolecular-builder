
#include"builder.h"

void retrieve_mol_type_data(int current_mol_type,int mol_types,int *type_array,int *MGR_start,int *MGR_stop,int *MA_start,int *MA_stop,int *MB_start,int *MB_stop,
                            int total_chains,char **quad_matrix,int *chains_array,
                            char **master_species_chain_array,int *master_B1_chain_array,int *master_B2_chain_array,char **master_B_type_chain_array,
                            int gcols,int **master_growth_registry,
                            int *grows,int **init_bo_array,int ***growth_registry,
                            int *total_chain_atoms,int *total_chain_bonds,int **chain_atoms_array,int **chain_bonds_array,char ***species_chain_array,
                            int **B1_chain_array,int **B2_chain_array,char ***B_type_chain_array,
                            char **master_atom_name_chain_array,char ***atom_name_chain_array)
{
    char word[cmax_length];
    int i,j,k,l,m,int_buffer;
    int current_mol_type_index;
    int chains;
    
    for(i=0;i<mol_types;++i)
    {
        if(current_mol_type==type_array[i])
        {
            current_mol_type_index=i;
            break;
        }
    }
    
    // chains
    chains=chains_array[current_mol_type_index];
    
    // grows
    *grows=MGR_stop[current_mol_type_index]-MGR_start[current_mol_type_index]+1;

    // *init_bo_array
    *init_bo_array=(int*)malloc(chains*sizeof(int));

    
    k=-1;
    for(j=0;j<total_chains;++j)
    {
        sscanf(quad_matrix[j],"%d",&int_buffer);
        if(int_buffer==current_mol_type)
        {
            k=k+1;
            sscanf(quad_matrix[j],"%d\t%s\t%d\t%d\t%d",&int_buffer,word,&int_buffer,&int_buffer,&(*init_bo_array)[k]);
        }
    }
    
    // **growth_registry
    *growth_registry=(int**)malloc(*grows*sizeof(int*));
    for(j=0;j<*grows;++j)(*growth_registry)[j]=(int*)malloc(gcols*sizeof(int));
    l=-1;
    // @ loop: nested_2
    for(j=MGR_start[i]-1;j<=MGR_stop[i]-1;++j)
    {
        l=l+1;
        for(k=0;k<gcols;++k)
            (*growth_registry)[l][k]=master_growth_registry[j][k];
    }
    
    // *chain_atoms_array, *chain_bonds_array, total_chain_atoms, total_chain_bonds
    *chain_atoms_array=(int*)malloc(chains*sizeof(int));
    *chain_bonds_array=(int*)malloc(chains*sizeof(int));
    *total_chain_atoms=0;
    *total_chain_bonds=0;
    k=-1;
    for(j=0;j<total_chains;++j)
    {
        sscanf(quad_matrix[j],"%d",&int_buffer);
        if(int_buffer==current_mol_type)
        {
            k=k+1;
            sscanf(quad_matrix[j],"%d\t%s\t%d\t%d\t%d\t%d\t%d",&int_buffer,word,&int_buffer,&int_buffer,&int_buffer,&(*chain_atoms_array)[k],&(*chain_bonds_array)[k]);
            *total_chain_atoms=*total_chain_atoms+(*chain_atoms_array)[k];
            *total_chain_bonds=*total_chain_bonds+(*chain_bonds_array)[k];
        }
    }
    
    // **species_chain_array
    *species_chain_array=(char**)malloc(*total_chain_atoms*sizeof(char*));
    for(j=0;j<*total_chain_atoms;++j)(*species_chain_array)[j]=(char*)malloc(sub_length*sizeof(char));
    //
    *atom_name_chain_array=(char**)malloc(*total_chain_atoms*sizeof(char*));
    for(j=0;j<*total_chain_atoms;++j)(*atom_name_chain_array)[j]=(char*)malloc(sub_length*sizeof(char));
    l=MA_start[i]-1;m=-1;
    // @ loop: nested_2
    for(j=0;j<chains;++j)
    {
        for(k=0;k<(*chain_atoms_array)[j];++k)
        {
            l=l+1;
            m=m+1;
            sprintf((*species_chain_array)[m],"%s",master_species_chain_array[l-1]);
            sprintf((*atom_name_chain_array)[m],"%s",master_atom_name_chain_array[l-1]);
        }
    }
    
    // *B1_chain_array, *B2_chain_array, **B_type_chain_array
    *B_type_chain_array=(char**)malloc(*total_chain_bonds*sizeof(char*));
    for(j=0;j<*total_chain_bonds;++j)(*B_type_chain_array)[j]=(char*)malloc(sub_length*sizeof(char));
    *B1_chain_array=(int*)malloc(*total_chain_bonds*sizeof(int));
    *B2_chain_array=(int*)malloc(*total_chain_bonds*sizeof(int));
    l=MB_start[i]-1;m=-1;
    // @ loop: nested_2
    for(j=0;j<chains;++j)
    {
        for(k=0;k<(*chain_bonds_array)[j];++k)
        {
            l=l+1;
            m=m+1;
            (*B1_chain_array)[m]=master_B1_chain_array[l-1];
            (*B2_chain_array)[m]=master_B2_chain_array[l-1];
            sprintf((*B_type_chain_array)[m],"%s",master_B_type_chain_array[l-1]);
        }
    }
    
    /*
     //
     printf("type=%d\n",type);
     printf("chains=%d\n",chains);
     printf("grows=%d\tgcols=%d\n",grows,gcols);
     printf("init_bo_array:\n");
     for(j=0;j<chains;++j)printf("%d\t",init_bo_array[j]);printf("\n");
     printf("growth_matrix:\n");
     for(j=0;j<grows;++j)
     {
     for(k=0;k<gcols;++k)
     printf("%d\t",growth_registry[j][k]);
     printf("\n");
     }
     printf("atoms_core=%d\n",atoms_core[type-1]);
     printf("bonds_core=%d\n",bonds_core[type-1]);
     printf("total_chain_atoms=%d\n",total_chain_atoms);
     printf("total_chain_bonds=%d\n",total_chain_bonds);
     printf("chain_atoms_array:\t");for(j=0;j<chains;++j)printf("%d\t",chain_atoms_array[j]);printf("\n");
     printf("chain_bonds_array:\t");for(j=0;j<chains;++j)printf("%d\t",chain_bonds_array[j]);printf("\n");
     printf("species_chain_array:\n");for(j=0;j<total_chain_atoms;++j)printf("%d\t%s\n",j+1,species_chain_array[j]);
     printf("B1_chain_array, B2_chain_array, B_type_chain_array:\n");for(j=0;j<total_chain_bonds;++j)printf("%d\t%d\t%d\t%s\n",j+1,B1_chain_array[j],B2_chain_array[j],B_type_chain_array[j]);
     */
    
}
