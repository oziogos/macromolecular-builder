
#include"builder.h"

void topo_diagnostic(int general_bonds, int **general_bonds_registry, char ***general_bond_types_registry,
                     int general_angles, int **general_angles_registry, char ***general_angle_types_registry,
                     int general_dihedrals, int **general_dihedrals_registry, char ***general_dihedral_types_registry,
                     int general_impropers, int **general_impropers_registry, char ***general_improper_types_registry,
                     int bonds, int **global_bonds_array, char ***global_bond_type_array,
                     int angles_core, int **global_angles_array, char ***global_angle_type_array,
                     int dihedrals_core, int **global_dihedrals_array, char ***global_dihedral_type_array,
                     int impropers_core, int **global_impropers_array, char ***global_improper_type_array)
{
    int j,k,*checked;

    // bonds
    checked=(int*)malloc(general_bonds*sizeof(int));
    for(j=0;j<general_bonds;++j)checked[j]=0;
    for(j=0;j<general_bonds;++j)
    {
        for(k=0;k<bonds;++k)
        {
            if((general_bonds_registry[j][1]==global_bonds_array[k][1] &&
                general_bonds_registry[j][2]==global_bonds_array[k][2] &&
                strcmp(general_bond_types_registry[general_bonds_registry[j][0]-1][0],global_bond_type_array[global_bonds_array[k][0]-1][0])==0 &&
                strcmp(general_bond_types_registry[general_bonds_registry[j][0]-1][1],global_bond_type_array[global_bonds_array[k][0]-1][1])==0 &&
                strcmp(general_bond_types_registry[general_bonds_registry[j][0]-1][2],global_bond_type_array[global_bonds_array[k][0]-1][2])==0
                )
               ||
               (general_bonds_registry[j][1]==global_bonds_array[k][2] &&
                general_bonds_registry[j][2]==global_bonds_array[k][1] &&
                strcmp(general_bond_types_registry[general_bonds_registry[j][0]-1][0],global_bond_type_array[global_bonds_array[k][0]-1][0])==0 &&
                strcmp(general_bond_types_registry[general_bonds_registry[j][0]-1][1],global_bond_type_array[global_bonds_array[k][0]-1][1])==0 &&
                strcmp(general_bond_types_registry[general_bonds_registry[j][0]-1][2],global_bond_type_array[global_bonds_array[k][0]-1][2])==0
                )
               ||
               (general_bonds_registry[j][1]==global_bonds_array[k][1] &&
                general_bonds_registry[j][2]==global_bonds_array[k][2] &&
                strcmp(general_bond_types_registry[general_bonds_registry[j][0]-1][0],global_bond_type_array[global_bonds_array[k][0]-1][2])==0 &&
                strcmp(general_bond_types_registry[general_bonds_registry[j][0]-1][1],global_bond_type_array[global_bonds_array[k][0]-1][1])==0 &&
                strcmp(general_bond_types_registry[general_bonds_registry[j][0]-1][2],global_bond_type_array[global_bonds_array[k][0]-1][0])==0
                )
               ||
               (general_bonds_registry[j][1]==global_bonds_array[k][2] &&
                general_bonds_registry[j][2]==global_bonds_array[k][1] &&
                strcmp(general_bond_types_registry[general_bonds_registry[j][0]-1][0],global_bond_type_array[global_bonds_array[k][0]-1][2])==0 &&
                strcmp(general_bond_types_registry[general_bonds_registry[j][0]-1][1],global_bond_type_array[global_bonds_array[k][0]-1][1])==0 &&
                strcmp(general_bond_types_registry[general_bonds_registry[j][0]-1][2],global_bond_type_array[global_bonds_array[k][0]-1][0])==0
                )
               )
            {
                checked[j]=checked[j]+1;
                break;
            }
        }
    }
    for(j=0;j<general_bonds;++j)
    {
        //printf("{%d}\t%d\t%d\t%s\t%s\n",checked[j],general_bonds_registry[j][1],general_bonds_registry[j][2],general_bond_types_registry[general_bonds_registry[j][0]-1][0],general_bond_types_registry[general_bonds_registry[j][0]-1][1]);
        if(checked[j]!=1){printf("B comparison failed!!\n\n");exit(-4);}
    }
    free(checked);
    // angles
    checked=(int*)malloc(general_angles*sizeof(int));
    for(j=0;j<general_angles;++j)checked[j]=0;
    for(j=0;j<general_angles;++j)
    {
        for(k=0;k<angles_core;++k)
        {
            if((general_angles_registry[j][1]==global_angles_array[k][1] &&
                general_angles_registry[j][2]==global_angles_array[k][2] &&
                general_angles_registry[j][3]==global_angles_array[k][3] &&
                strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][0],global_angle_type_array[global_angles_array[k][0]-1][0])==0 &&
                strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][1],global_angle_type_array[global_angles_array[k][0]-1][1])==0 &&
                strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][2],global_angle_type_array[global_angles_array[k][0]-1][2])==0 &&
                strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][3],global_angle_type_array[global_angles_array[k][0]-1][3])==0 &&
                strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][4],global_angle_type_array[global_angles_array[k][0]-1][4])==0
                )
               ||
               (general_angles_registry[j][1]==global_angles_array[k][3] &&
                general_angles_registry[j][2]==global_angles_array[k][2] &&
                general_angles_registry[j][3]==global_angles_array[k][1] &&
                strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][0],global_angle_type_array[global_angles_array[k][0]-1][0])==0 &&
                strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][1],global_angle_type_array[global_angles_array[k][0]-1][1])==0 &&
                strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][2],global_angle_type_array[global_angles_array[k][0]-1][2])==0 &&
                strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][3],global_angle_type_array[global_angles_array[k][0]-1][3])==0 &&
                strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][4],global_angle_type_array[global_angles_array[k][0]-1][4])==0
                )
               ||
               (general_angles_registry[j][1]==global_angles_array[k][1] &&
                general_angles_registry[j][2]==global_angles_array[k][2] &&
                general_angles_registry[j][3]==global_angles_array[k][3] &&
                strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][0],global_angle_type_array[global_angles_array[k][0]-1][4])==0 &&
                strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][1],global_angle_type_array[global_angles_array[k][0]-1][3])==0 &&
                strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][2],global_angle_type_array[global_angles_array[k][0]-1][2])==0 &&
                strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][3],global_angle_type_array[global_angles_array[k][0]-1][1])==0 &&
                strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][4],global_angle_type_array[global_angles_array[k][0]-1][0])==0
                )
               ||
               (general_angles_registry[j][1]==global_angles_array[k][3] &&
                general_angles_registry[j][2]==global_angles_array[k][2] &&
                general_angles_registry[j][3]==global_angles_array[k][1] &&
                strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][0],global_angle_type_array[global_angles_array[k][0]-1][4])==0 &&
                strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][1],global_angle_type_array[global_angles_array[k][0]-1][3])==0 &&
                strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][2],global_angle_type_array[global_angles_array[k][0]-1][2])==0 &&
                strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][3],global_angle_type_array[global_angles_array[k][0]-1][1])==0 &&
                strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][4],global_angle_type_array[global_angles_array[k][0]-1][0])==0
                )
               )
            {
                checked[j]=checked[j]+1;
                break;
            }
        }
    }
    for(j=0;j<general_angles;++j)
    {
        //printf("{%d}\t%d\t%d\t%d\t%s\t%s\t%s\n",checked[j],general_angles_registry[j][1],general_angles_registry[j][2],general_angles_registry[j][3],general_angle_types_registry[general_angles_registry[j][0]-1][0],general_angle_types_registry[general_angles_registry[j][0]-1][1],general_angle_types_registry[general_angles_registry[j][0]-1][2]);
        if(checked[j]!=1){printf("A comparison failed!!\n\n");exit(-4);}
    }
    free(checked);
    // dihedrals
    checked=(int*)malloc(general_dihedrals*sizeof(int));
    for(j=0;j<general_dihedrals;++j)checked[j]=0;
    for(j=0;j<general_dihedrals;++j)
    {
        for(k=0;k<dihedrals_core;++k)
        {
            if((general_dihedrals_registry[j][1]==global_dihedrals_array[k][1] &&
                general_dihedrals_registry[j][2]==global_dihedrals_array[k][2] &&
                general_dihedrals_registry[j][3]==global_dihedrals_array[k][3] &&
                general_dihedrals_registry[j][4]==global_dihedrals_array[k][4] &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][0],global_dihedral_type_array[global_dihedrals_array[k][0]-1][0])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][1],global_dihedral_type_array[global_dihedrals_array[k][0]-1][1])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][2],global_dihedral_type_array[global_dihedrals_array[k][0]-1][2])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][3],global_dihedral_type_array[global_dihedrals_array[k][0]-1][3])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][4],global_dihedral_type_array[global_dihedrals_array[k][0]-1][4])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][5],global_dihedral_type_array[global_dihedrals_array[k][0]-1][5])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][6],global_dihedral_type_array[global_dihedrals_array[k][0]-1][6])==0
                )
               ||
               (general_dihedrals_registry[j][1]==global_dihedrals_array[k][4] &&
                general_dihedrals_registry[j][2]==global_dihedrals_array[k][3] &&
                general_dihedrals_registry[j][3]==global_dihedrals_array[k][2] &&
                general_dihedrals_registry[j][4]==global_dihedrals_array[k][1] &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][0],global_dihedral_type_array[global_dihedrals_array[k][0]-1][0])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][1],global_dihedral_type_array[global_dihedrals_array[k][0]-1][1])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][2],global_dihedral_type_array[global_dihedrals_array[k][0]-1][2])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][3],global_dihedral_type_array[global_dihedrals_array[k][0]-1][3])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][4],global_dihedral_type_array[global_dihedrals_array[k][0]-1][4])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][5],global_dihedral_type_array[global_dihedrals_array[k][0]-1][5])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][6],global_dihedral_type_array[global_dihedrals_array[k][0]-1][6])==0
                )
               ||
               (general_dihedrals_registry[j][1]==global_dihedrals_array[k][1] &&
                general_dihedrals_registry[j][2]==global_dihedrals_array[k][2] &&
                general_dihedrals_registry[j][3]==global_dihedrals_array[k][3] &&
                general_dihedrals_registry[j][4]==global_dihedrals_array[k][4] &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][0],global_dihedral_type_array[global_dihedrals_array[k][0]-1][6])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][1],global_dihedral_type_array[global_dihedrals_array[k][0]-1][5])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][2],global_dihedral_type_array[global_dihedrals_array[k][0]-1][4])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][3],global_dihedral_type_array[global_dihedrals_array[k][0]-1][3])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][4],global_dihedral_type_array[global_dihedrals_array[k][0]-1][2])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][5],global_dihedral_type_array[global_dihedrals_array[k][0]-1][1])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][6],global_dihedral_type_array[global_dihedrals_array[k][0]-1][0])==0
                )
               ||
               (general_dihedrals_registry[j][1]==global_dihedrals_array[k][4] &&
                general_dihedrals_registry[j][2]==global_dihedrals_array[k][3] &&
                general_dihedrals_registry[j][3]==global_dihedrals_array[k][2] &&
                general_dihedrals_registry[j][4]==global_dihedrals_array[k][1] &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][0],global_dihedral_type_array[global_dihedrals_array[k][0]-1][6])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][1],global_dihedral_type_array[global_dihedrals_array[k][0]-1][5])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][2],global_dihedral_type_array[global_dihedrals_array[k][0]-1][4])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][3],global_dihedral_type_array[global_dihedrals_array[k][0]-1][3])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][4],global_dihedral_type_array[global_dihedrals_array[k][0]-1][2])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][5],global_dihedral_type_array[global_dihedrals_array[k][0]-1][1])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][6],global_dihedral_type_array[global_dihedrals_array[k][0]-1][0])==0
                )
               )
            {
                checked[j]=checked[j]+1;
                break;
            }
        }
    }
    for(j=0;j<general_dihedrals;++j)
    {
        //printf("{%d}\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n",checked[j],general_dihedrals_registry[j][1],general_dihedrals_registry[j][2],general_dihedrals_registry[j][3],general_dihedrals_registry[j][4],general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][0],general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][1],general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][2],general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][3]);
        if(checked[j]!=1){printf("D comparison failed!!\n\n");exit(-4);}
    }
    free(checked);
    // impropers
    checked=(int*)malloc(general_impropers*sizeof(int));
    for(j=0;j<general_impropers;++j)checked[j]=0;
    for(j=0;j<general_impropers;++j)
    {
        for(k=0;k<impropers_core;++k)
        {
            if((general_impropers_registry[j][1]==global_impropers_array[k][1] &&
                general_impropers_registry[j][2]==global_impropers_array[k][2] &&
                general_impropers_registry[j][3]==global_impropers_array[k][3] &&
                general_impropers_registry[j][4]==global_impropers_array[k][4]
                )
               ||
               (general_impropers_registry[j][1]==global_impropers_array[k][1] &&
                general_impropers_registry[j][2]==global_impropers_array[k][2] &&
                general_impropers_registry[j][3]==global_impropers_array[k][4] &&
                general_impropers_registry[j][4]==global_impropers_array[k][3]
                )
               ||
               (general_impropers_registry[j][1]==global_impropers_array[k][1] &&
                general_impropers_registry[j][2]==global_impropers_array[k][3] &&
                general_impropers_registry[j][3]==global_impropers_array[k][2] &&
                general_impropers_registry[j][4]==global_impropers_array[k][4]
                )
               ||
               (general_impropers_registry[j][1]==global_impropers_array[k][1] &&
                general_impropers_registry[j][2]==global_impropers_array[k][3] &&
                general_impropers_registry[j][3]==global_impropers_array[k][4] &&
                general_impropers_registry[j][4]==global_impropers_array[k][2]
                )
               ||
               (general_impropers_registry[j][1]==global_impropers_array[k][1] &&
                general_impropers_registry[j][2]==global_impropers_array[k][4] &&
                general_impropers_registry[j][3]==global_impropers_array[k][2] &&
                general_impropers_registry[j][4]==global_impropers_array[k][3]
                )
               ||
               (general_impropers_registry[j][1]==global_impropers_array[k][1] &&
                general_impropers_registry[j][2]==global_impropers_array[k][4] &&
                general_impropers_registry[j][3]==global_impropers_array[k][3] &&
                general_impropers_registry[j][4]==global_impropers_array[k][2]
                )
               )
            {
                
                if((strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][0],global_improper_type_array[global_impropers_array[k][0]-1][0])==0 &&
                    strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][1],global_improper_type_array[global_impropers_array[k][0]-1][1])==0 &&
                    strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][2],global_improper_type_array[global_impropers_array[k][0]-1][2])==0 &&
                    strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][3],global_improper_type_array[global_impropers_array[k][0]-1][3])==0
                    )
                   ||
                   (strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][0],global_improper_type_array[global_impropers_array[k][0]-1][0])==0 &&
                    strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][1],global_improper_type_array[global_impropers_array[k][0]-1][1])==0 &&
                    strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][2],global_improper_type_array[global_impropers_array[k][0]-1][3])==0 &&
                    strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][3],global_improper_type_array[global_impropers_array[k][0]-1][2])==0
                    )
                   ||
                   (strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][0],global_improper_type_array[global_impropers_array[k][0]-1][0])==0 &&
                    strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][1],global_improper_type_array[global_impropers_array[k][0]-1][2])==0 &&
                    strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][2],global_improper_type_array[global_impropers_array[k][0]-1][1])==0 &&
                    strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][3],global_improper_type_array[global_impropers_array[k][0]-1][3])==0
                    )
                   ||
                   (strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][0],global_improper_type_array[global_impropers_array[k][0]-1][0])==0 &&
                    strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][1],global_improper_type_array[global_impropers_array[k][0]-1][2])==0 &&
                    strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][2],global_improper_type_array[global_impropers_array[k][0]-1][3])==0 &&
                    strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][3],global_improper_type_array[global_impropers_array[k][0]-1][1])==0
                    )
                   ||
                   (strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][0],global_improper_type_array[global_impropers_array[k][0]-1][0])==0 &&
                    strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][1],global_improper_type_array[global_impropers_array[k][0]-1][3])==0 &&
                    strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][2],global_improper_type_array[global_impropers_array[k][0]-1][1])==0 &&
                    strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][3],global_improper_type_array[global_impropers_array[k][0]-1][2])==0
                    )
                   ||
                   (strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][0],global_improper_type_array[global_impropers_array[k][0]-1][0])==0 &&
                    strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][1],global_improper_type_array[global_impropers_array[k][0]-1][3])==0 &&
                    strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][2],global_improper_type_array[global_impropers_array[k][0]-1][2])==0 &&
                    strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][3],global_improper_type_array[global_impropers_array[k][0]-1][1])==0
                    )
                   )
                {
                    checked[j]=checked[j]+1;
                    break;
                }
            }
        }
    }
    for(j=0;j<general_impropers;++j)
    {
        //printf("{%d}\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n",checked[j],general_impropers_registry[j][1],general_impropers_registry[j][2],general_impropers_registry[j][3],general_impropers_registry[j][4],general_improper_types_registry[general_impropers_registry[j][0]-1][0],general_improper_types_registry[general_impropers_registry[j][0]-1][1],general_improper_types_registry[general_impropers_registry[j][0]-1][2],general_improper_types_registry[general_impropers_registry[j][0]-1][3]);
        if(checked[j]!=1){printf("I comparison failed!!\n\n");exit(-4);}
    }
    free(checked);
    // end of diagnostic
    
}
