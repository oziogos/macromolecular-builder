
#include"builder.h"

void topo_general_backup(int general_bonds,int general_angles,int general_dihedrals,int general_impropers,
                         int general_atom_types,int general_bond_types,int general_angle_types,int general_dihedral_types,int general_improper_types,
                         char **general_species_registry,
                         char ***general_bond_types_registry,char ***general_angle_types_registry,char ***general_dihedral_types_registry,char ***general_improper_types_registry,
                         int *general_bonds_B,int *general_angles_B,int *general_dihedrals_B,int *general_impropers_B,
                         int *general_atom_types_B,int *general_bond_types_B,int *general_angle_types_B,int *general_dihedral_types_B,int *general_improper_types_B,
                         char ***general_species_registry_B,
                         char ****general_bond_types_registry_B,char ****general_angle_types_registry_B,char ****general_dihedral_types_registry_B,char ****general_improper_types_registry_B,
                         int general_atoms,int *general_neighbors_registry,int **general_neighbors_registry_B)
{
    
    int j,k;

    // backup general topology elements
    *general_bonds_B=general_bonds;
    *general_angles_B=general_angles;
    *general_dihedrals_B=general_dihedrals;
    *general_impropers_B=general_impropers;
    *general_atom_types_B=general_atom_types;
    *general_bond_types_B=general_bond_types;
    *general_angle_types_B=general_angle_types;
    *general_dihedral_types_B=general_dihedral_types;
    *general_improper_types_B=general_improper_types;
    //
    *general_species_registry_B=(char**)malloc(*general_atom_types_B*sizeof(char*));
    for(j=0;j<*general_atom_types_B;++j)(*general_species_registry_B)[j]=(char*)malloc(sub_length*sizeof(char));
    *general_bond_types_registry_B=(char***)malloc(*general_bond_types_B*sizeof(char**));
    for(j=0;j<*general_bond_types_B;++j)(*general_bond_types_registry_B)[j]=(char**)malloc(3*sizeof(char*));
    for(j=0;j<*general_bond_types_B;++j)for(k=0;k<3;++k)(*general_bond_types_registry_B)[j][k]=(char*)malloc(sub_length*sizeof(char));
    *general_angle_types_registry_B=(char***)malloc(*general_angle_types_B*sizeof(char**));
    for(j=0;j<*general_angle_types_B;++j)(*general_angle_types_registry_B)[j]=(char**)malloc(5*sizeof(char*));
    for(j=0;j<*general_angle_types_B;++j)for(k=0;k<5;++k)(*general_angle_types_registry_B)[j][k]=(char*)malloc(sub_length*sizeof(char));
    if(*general_dihedrals_B>0)
    {
        *general_dihedral_types_registry_B=(char***)malloc(*general_dihedral_types_B*sizeof(char**));
        for(j=0;j<*general_dihedral_types_B;++j)(*general_dihedral_types_registry_B)[j]=(char**)malloc(7*sizeof(char*));
        for(j=0;j<*general_dihedral_types_B;++j)for(k=0;k<7;++k)(*general_dihedral_types_registry_B)[j][k]=(char*)malloc(sub_length*sizeof(char));
    }
    if(*general_impropers_B>0)
    {
        *general_improper_types_registry_B=(char***)malloc(*general_improper_types_B*sizeof(char**));
        for(j=0;j<*general_improper_types_B;++j)(*general_improper_types_registry_B)[j]=(char**)malloc(4*sizeof(char*));
        for(j=0;j<*general_improper_types_B;++j)for(k=0;k<4;++k)(*general_improper_types_registry_B)[j][k]=(char*)malloc(sub_length*sizeof(char));
    }
    //
    for(j=0;j<general_atom_types;++j)
        sprintf((*general_species_registry_B)[j],"%s",general_species_registry[j]);
    for(j=0;j<general_bond_types;++j){
        sprintf((*general_bond_types_registry_B)[j][0],"%s",general_bond_types_registry[j][0]);
        sprintf((*general_bond_types_registry_B)[j][1],"%s",general_bond_types_registry[j][1]); // BO
        sprintf((*general_bond_types_registry_B)[j][2],"%s",general_bond_types_registry[j][2]);
    }
    for(j=0;j<general_angle_types;++j)
    {
        sprintf((*general_angle_types_registry_B)[j][0],"%s",general_angle_types_registry[j][0]);
        sprintf((*general_angle_types_registry_B)[j][1],"%s",general_angle_types_registry[j][1]);   // BO
        sprintf((*general_angle_types_registry_B)[j][2],"%s",general_angle_types_registry[j][2]);
        sprintf((*general_angle_types_registry_B)[j][3],"%s",general_angle_types_registry[j][3]);   // BO
        sprintf((*general_angle_types_registry_B)[j][4],"%s",general_angle_types_registry[j][4]);
    }
    if(*general_dihedrals_B>0)
    {
        for(j=0;j<general_dihedral_types;++j)
        {
            sprintf((*general_dihedral_types_registry_B)[j][0],"%s",general_dihedral_types_registry[j][0]);
            sprintf((*general_dihedral_types_registry_B)[j][1],"%s",general_dihedral_types_registry[j][1]); // BO
            sprintf((*general_dihedral_types_registry_B)[j][2],"%s",general_dihedral_types_registry[j][2]);
            sprintf((*general_dihedral_types_registry_B)[j][3],"%s",general_dihedral_types_registry[j][3]); // BO
            sprintf((*general_dihedral_types_registry_B)[j][4],"%s",general_dihedral_types_registry[j][4]);
            sprintf((*general_dihedral_types_registry_B)[j][5],"%s",general_dihedral_types_registry[j][5]); // BO
            sprintf((*general_dihedral_types_registry_B)[j][6],"%s",general_dihedral_types_registry[j][6]);
        }
        
    }
    if(*general_impropers_B>0)
    {
        for(j=0;j<general_improper_types;++j)
        {
            sprintf((*general_improper_types_registry_B)[j][0],"%s",general_improper_types_registry[j][0]);
            sprintf((*general_improper_types_registry_B)[j][1],"%s",general_improper_types_registry[j][1]);
            sprintf((*general_improper_types_registry_B)[j][2],"%s",general_improper_types_registry[j][2]);
            sprintf((*general_improper_types_registry_B)[j][3],"%s",general_improper_types_registry[j][3]);
        }
    }
    //

    //
    *general_neighbors_registry_B=(int*)malloc(general_atoms*sizeof(int));
    for(j=0;j<general_atoms;++j)(*general_neighbors_registry_B)[j]=general_neighbors_registry[j];


}
