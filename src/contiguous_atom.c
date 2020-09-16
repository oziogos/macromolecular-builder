
#include"builder.h"

void contiguous_atom(int *atoms,int *bonds,int molecules,int mol_types,int *mol_id_array,int *mol_type_array,int *type_array,char **core_file,
                     char *current_folder,int **master_mol_array,int **master_id_internal_array,double **master_x_array,double **master_y_array,double **master_z_array,
                     double **master_q_array,
                     char ***master_species_array,int **master_id_array,int **master_mol_type_array,int total_chains,char **quad_matrix,int **host_id_array,int **chain_id_array,
                     int **atoms_per_molecule,int **atom_scaling_array,
                     char ***master_atom_name_array)
{
    FILE *fp;
    char file_path[cmax_length],buffer[cmax_length],word[cmax_length];
    int i,j,k,l,int_buffer,graft_id,host_id;
    int atoms_core,bonds_core;
    
    //
    *atoms_per_molecule=(int*)malloc(molecules*sizeof(int));
    *atom_scaling_array=(int*)malloc(molecules*sizeof(int));
    
    *atoms=0;*bonds=0;                                    // total number of atoms and bonds
    
    // open core mol2 files and store atoms and bonds (per mol type and total)
    for(i=0;i<molecules;++i)
    {
        for(j=0;j<mol_types;++j)
        {
            if(mol_type_array[i]==type_array[j])
                break;
        }
        sprintf(file_path,"%s/%s",current_folder,core_file[j]);
        fp=fopen(file_path,"r");
        if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
        while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"@<TRIPOS>MOLECULE\n")==0)break;          // locate the MOLECULE section
        fgets(buffer,cmax_length,fp);fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d",&atoms_core,&bonds_core);    // skip comment line, read atoms and bonds from the next one
        fclose(fp);
        *atoms=*atoms+atoms_core;
        *bonds=*bonds+bonds_core;
        //
        (*atoms_per_molecule)[i]=atoms_core;
    }
    
    //
    (*atom_scaling_array)[0]=0;
    for(i=1;i<molecules;++i)
    {
        (*atom_scaling_array)[i]=(*atom_scaling_array)[i-1]+(*atoms_per_molecule)[i-1];
    }
    
    // preallocate arrays to store atomic data in a contiguous manner
    *master_mol_array=(int*)malloc(*atoms*sizeof(int));
    *master_id_internal_array=(int*)malloc(*atoms*sizeof(int));
    *master_x_array=(double*)malloc(*atoms*sizeof(double));
    *master_y_array=(double*)malloc(*atoms*sizeof(double));
    *master_z_array=(double*)malloc(*atoms*sizeof(double));
    *master_q_array=(double*)malloc(*atoms*sizeof(double));
    *master_species_array=(char**)malloc(*atoms*sizeof(char*));
    for(i=0;i<*atoms;++i)(*master_species_array)[i]=(char*)malloc(sub_length*sizeof(char));
    *master_atom_name_array=(char**)malloc(*atoms*sizeof(char*));
    for(i=0;i<*atoms;++i)(*master_atom_name_array)[i]=(char*)malloc(sub_length*sizeof(char));
    *master_id_array=(int*)malloc(*atoms*sizeof(int));
    *master_mol_type_array=(int*)malloc(*atoms*sizeof(int));
    
    // store contiguous data
    l=-1;
    for(i=0;i<molecules;++i)
    {
        for(j=0;j<mol_types;++j)
        {
            if(mol_type_array[i]==type_array[j])
                break;
        }
        sprintf(file_path,"%s/%s",current_folder,core_file[j]);
        fp=fopen(file_path,"r");
        if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
        while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"@<TRIPOS>MOLECULE\n")==0)break;          // locate the MOLECULE section
        fgets(buffer,cmax_length,fp);fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d",&atoms_core,&bonds_core);    // skip comment line, read atoms and bonds from the next one
        
        while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"@<TRIPOS>ATOM\n")==0)break;
        for(k=0;k<atoms_core;++k)
        {
            l=l+1;
            fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%s\t%lf\t%lf\t%lf\t%s",&(*master_id_internal_array)[l],(*master_atom_name_array)[l],&(*master_x_array)[l],&(*master_y_array)[l],&(*master_z_array)[l],(*master_species_array)[l]);
            (*master_mol_array)[l]=mol_id_array[i];
        }
        fclose(fp);
    }
    
    // running atomic index and molecular type arrays
    for(i=0;i<*atoms;++i)(*master_id_array)[i]=i+1;
    for(i=0;i<*atoms;++i)(*master_mol_type_array)[i]=mol_type_array[(*master_mol_array)[i]-1];
    
    // preallocate and initialize *host_id_array,*chain_id_array
    *host_id_array=(int*)malloc(*atoms*sizeof(int));
    *chain_id_array=(int*)malloc(*atoms*sizeof(int));
    // initialize
    for(i=0;i<*atoms;++i)
    {
        if(strcmp((*master_species_array)[i],"H")==0){
            (*host_id_array)[i]=-2;}
        else{
            (*host_id_array)[i]=-1;}
        (*chain_id_array)[i]=-1;
    }
    // define the host atoms; label '0' (host_id)
    for(i=0;i<molecules;++i)
    {
        k=-1;
        for(j=0;j<total_chains;++j)
        {
            sscanf(quad_matrix[j],"%d",&int_buffer);
            if(int_buffer==mol_type_array[i])
            {
                k=k+1;
                sscanf(quad_matrix[j],"%d\t%s\t%d\t%d",&int_buffer,word,&graft_id,&host_id);
                
                for(l=0;l<*atoms;++l)
                {
                    if((*master_mol_array)[l]==mol_id_array[i] && (*master_id_internal_array)[l]==host_id)
                    {
                        (*host_id_array)[l]=0;
                        (*chain_id_array)[l]=k;
                        break;
                    }
                }
                
            }
        }
    }

    //
    for(i=0;i<*atoms;++i)(*master_q_array)[i]=0.0;
}
