
#include"builder.h"

int local_bonds,local_angles,local_dihedrals,local_impropers;
int local_atom_types,local_bond_types,local_angle_types,local_dihedral_types,local_improper_types;
int **local_bonds_registry,**local_angles_registry,**local_dihedrals_registry,**local_impropers_registry;
char **local_species_registry,***local_bond_types_registry,***local_angle_types_registry,***local_dihedral_types_registry,***local_improper_types_registry;
int *local_neighbors_registry;

void free_locals()
{
    int j,k;

    for(j=0;j<local_bonds;++j)free(local_bonds_registry[j]);free(local_bonds_registry);
    for(j=0;j<local_angles;++j)free(local_angles_registry[j]);free(local_angles_registry);
    for(j=0;j<local_dihedrals;++j)free(local_dihedrals_registry[j]);free(local_dihedrals_registry);
    if(local_impropers>0){for(j=0;j<local_impropers;++j)free(local_impropers_registry[j]);free(local_impropers_registry);}
    for(k=0;k<local_atom_types;++k)free(local_species_registry[k]);free(local_species_registry);
    for(k=0;k<local_bond_types;++k)
        for(j=0;j<3;++j)
            free(local_bond_types_registry[k][j]);
    for(k=0;k<local_bond_types;++k)free(local_bond_types_registry[k]);
    free(local_bond_types_registry);
    for(k=0;k<local_angle_types;++k)
        for(j=0;j<5;++j)
            free(local_angle_types_registry[k][j]);
    for(k=0;k<local_angle_types;++k)free(local_angle_types_registry[k]);
    free(local_angle_types_registry);
    for(k=0;k<local_dihedral_types;++k)
        for(j=0;j<7;++j)
            free(local_dihedral_types_registry[k][j]);
    for(k=0;k<local_dihedral_types;++k)free(local_dihedral_types_registry[k]);
    free(local_dihedral_types_registry);
    free(local_neighbors_registry);
    if(local_impropers>0)
    {
        for(k=0;k<local_improper_types;++k)
            for(j=0;j<4;++j)
                free(local_improper_types_registry[k][j]);
        for(k=0;k<local_improper_types;++k)free(local_improper_types_registry[k]);
        free(local_improper_types_registry);
    }

}
