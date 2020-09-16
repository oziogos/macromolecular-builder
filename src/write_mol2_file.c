
#include"builder.h"

void write_mol2_file(char *current_folder,int atoms,int bonds,
                     int *master_id_array,char **master_atom_name_array,
                     double *master_x_array,double *master_y_array,double *master_z_array,
                     char **master_species_array,int *master_mol_array,
                     int *master_B1_core_array,int *master_B2_core_array,char **master_B_type_core_array,
                     char *mol2_file_name)
{
    
    FILE *fp;
    char file_path[cmax_length];
    int j;
    
    sprintf(file_path,"%s/%s",current_folder,mol2_file_name);
    fp=fopen(file_path,"w+");
    fprintf(fp,"@<TRIPOS>MOLECULE\n");
    fprintf(fp,"*****\n");
    fprintf(fp,"%d\t%d\t%d\t%d\t%d\n",atoms,bonds,0,0,0);
    fprintf(fp,"SMALL\n");
    fprintf(fp,"NO_CHARGES\n");
    fprintf(fp,"\n");
    fprintf(fp,"@<TRIPOS>ATOM\n");
    for(j=0;j<atoms;++j)fprintf(fp,"%d\t%s\t%lf\t%lf\t%lf\t%s\t%d\tLIG%d\n",master_id_array[j],master_atom_name_array[j],master_x_array[j],master_y_array[j],master_z_array[j],master_species_array[j],master_mol_array[j],master_mol_array[j]);
    fprintf(fp,"@<TRIPOS>BOND\n");
    for(j=0;j<bonds;++j)fprintf(fp,"%d\t%d\t%d\t%s\n",j+1,master_B1_core_array[j],master_B2_core_array[j],master_B_type_core_array[j]);
    fclose(fp);
}
