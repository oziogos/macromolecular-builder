
#include"builder.h"

void tmap_topo_boundaries_init(int molecules, int *master_mol_array, char **master_species_array,
                               int atoms, int general_atom_types, char **general_species_registry,
                               int general_bonds, int general_bond_types, char ***general_bond_types_registry, int **general_bonds_registry,
                               int general_angles, int general_angle_types, char ***general_angle_types_registry, int **general_angles_registry,
                               int general_dihedrals, int general_dihedral_types, char ***general_dihedral_types_registry, int **general_dihedrals_registry,
                               int general_impropers, int general_improper_types, char ***general_improper_types_registry, int **general_impropers_registry,
                               //
                               int *stmap_rows, char ***stmap, int ***scm,
                               int *btmap_rows, char ****btmap, int ***bcm, int *atmap_rows, char ****atmap, int ***acm,
                               int *dtmap_rows, char ****dtmap, int ***dcm, int *itmap_rows, char ****itmap, int ***icm,
                               int ***topo_boundaries,
                               
                               int verb)
{
    
    int i,j,k,l;

    // __tmap__________
    // species
    
    *stmap_rows=general_atom_types;
    *stmap=(char**)malloc((*stmap_rows)*sizeof(char*)); // prealloc and populate stmap
    for(j=0;j<*stmap_rows;++j)(*stmap)[j]=(char*)malloc(sub_length*sizeof(char));
    for(j=0;j<*stmap_rows;++j)sprintf((*stmap)[j],"%s",general_species_registry[j]);
    *scm=(int**)malloc(molecules*sizeof(int*));  // prealloc and initialize scm
    for(j=0;j<molecules;++j)(*scm)[j]=(int*)malloc((*stmap_rows)*sizeof(int));
    for(j=0;j<molecules;++j)
        for(k=0;k<*stmap_rows;++k)
            (*scm)[j][k]=0;
    // populate scm
    for(j=0;j<atoms;++j)
    {
        // locate species in stmap
        for(k=0;k<*stmap_rows;++k)
            if(strcmp(master_species_array[j],(*stmap)[k])==0)
                break;
        (*scm)[master_mol_array[j]-1][k]=1;
    }
    /*
     // console out
     for(j=0;j<*stmap_rows;++j)printf("%s\t",(*stmap)[j]);printf("\n");
     for(j=0;j<molecules;++j){
     for(k=0;k<*stmap_rows;++k)
     printf("%d\t",(*scm)[j][k]);
     printf("\n");
     }
     */
    // __tmap__________
    // bonds
    
    *btmap_rows=general_bond_types;
    *btmap=(char***)malloc((*btmap_rows)*sizeof(char**)); // prealloc and populate tmap
    for(j=0;j<*btmap_rows;++j)(*btmap)[j]=(char**)malloc(3*sizeof(char*));
    for(j=0;j<*btmap_rows;++j)for(k=0;k<3;++k)(*btmap)[j][k]=(char*)malloc(cmax_length*sizeof(char));
    for(j=0;j<*btmap_rows;++j){
        sprintf((*btmap)[j][0],"%s",general_bond_types_registry[j][0]);
        sprintf((*btmap)[j][1],"%s",general_bond_types_registry[j][1]); // BO
        sprintf((*btmap)[j][2],"%s",general_bond_types_registry[j][2]);
    }
    *bcm=(int**)malloc(molecules*sizeof(int*));  // prealloc and initialize cm
    for(j=0;j<molecules;++j)(*bcm)[j]=(int*)malloc((*btmap_rows)*sizeof(int));
    for(j=0;j<molecules;++j)
        for(k=0;k<*btmap_rows;++k)
            (*bcm)[j][k]=0;
    // populate cm
    for(j=0;j<general_bonds;++j)
    {
        // locate in tmap
        for(k=0;k<*btmap_rows;++k)
            if(
               (strcmp(general_bond_types_registry[general_bonds_registry[j][0]-1][0],(*btmap)[k][0])==0 &&
                strcmp(general_bond_types_registry[general_bonds_registry[j][0]-1][1],(*btmap)[k][1])==0 &&
                strcmp(general_bond_types_registry[general_bonds_registry[j][0]-1][2],(*btmap)[k][2])==0
                )||
               (strcmp(general_bond_types_registry[general_bonds_registry[j][0]-1][0],(*btmap)[k][2])==0 &&
                strcmp(general_bond_types_registry[general_bonds_registry[j][0]-1][1],(*btmap)[k][1])==0 &&
                strcmp(general_bond_types_registry[general_bonds_registry[j][0]-1][2],(*btmap)[k][0])==0
                )
               )
                break;
        (*bcm)[master_mol_array[general_bonds_registry[j][1]-1]-1][k]=1;
    }
    /*
     // console out
     for(j=0;j<*btmap_rows;++j)printf("%s-%s\t",(*btmap)[j][0],(*btmap)[j][1]);printf("\n");
     for(j=0;j<molecules;++j){
     for(k=0;k<*btmap_rows;++k)
     printf("%d\t",(*bcm)[j][k]);
     printf("\n");
     }
     */
    // __tmap__________
    // angles
    
    *atmap_rows=general_angle_types;
    *atmap=(char***)malloc((*atmap_rows)*sizeof(char**)); // prealloc and populate tmap
    for(j=0;j<*atmap_rows;++j)(*atmap)[j]=(char**)malloc(5*sizeof(char*));
    for(j=0;j<*atmap_rows;++j)for(k=0;k<5;++k)(*atmap)[j][k]=(char*)malloc(cmax_length*sizeof(char));
    for(j=0;j<*atmap_rows;++j){
        sprintf((*atmap)[j][0],"%s",general_angle_types_registry[j][0]);
        sprintf((*atmap)[j][1],"%s",general_angle_types_registry[j][1]);    // BO
        sprintf((*atmap)[j][2],"%s",general_angle_types_registry[j][2]);
        sprintf((*atmap)[j][3],"%s",general_angle_types_registry[j][3]);    // BO
        sprintf((*atmap)[j][4],"%s",general_angle_types_registry[j][4]);
    }
    *acm=(int**)malloc(molecules*sizeof(int*));  // prealloc and initialize cm
    for(j=0;j<molecules;++j)(*acm)[j]=(int*)malloc((*atmap_rows)*sizeof(int));
    for(j=0;j<molecules;++j)
        for(k=0;k<*atmap_rows;++k)
            (*acm)[j][k]=0;
    // populate cm
    for(j=0;j<general_angles;++j)
    {
        // locate in tmap
        for(k=0;k<*atmap_rows;++k)
            if(
               (strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][0],(*atmap)[k][0])==0 &&
                strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][1],(*atmap)[k][1])==0 &&
                strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][2],(*atmap)[k][2])==0 &&
                strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][3],(*atmap)[k][3])==0 &&
                strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][4],(*atmap)[k][4])==0
                )||
               (strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][0],(*atmap)[k][4])==0 &&
                strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][1],(*atmap)[k][3])==0 &&
                strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][2],(*atmap)[k][2])==0 &&
                strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][3],(*atmap)[k][1])==0 &&
                strcmp(general_angle_types_registry[general_angles_registry[j][0]-1][4],(*atmap)[k][0])==0
                )
               )
                break;
        (*acm)[master_mol_array[general_angles_registry[j][1]-1]-1][k]=1;
    }
    /*
     // console out
     for(j=0;j<*atmap_rows;++j)printf("%s-%s-%s\t",(*atmap)[j][0],(*atmap)[j][1],(*atmap)[j][2]);printf("\n");
     for(j=0;j<molecules;++j){
     for(k=0;k<*atmap_rows;++k)
     printf("%d\t",(*acm)[j][k]);
     printf("\n");
     }
     */
    // __tmap__________
    // dihedrals
    
    *dtmap_rows=general_dihedral_types;
    *dtmap=(char***)malloc((*dtmap_rows)*sizeof(char**)); // prealloc and populate tmap
    for(j=0;j<*dtmap_rows;++j)(*dtmap)[j]=(char**)malloc(7*sizeof(char*));
    for(j=0;j<*dtmap_rows;++j)for(k=0;k<7;++k)(*dtmap)[j][k]=(char*)malloc(cmax_length*sizeof(char));
    for(j=0;j<*dtmap_rows;++j){
        sprintf((*dtmap)[j][0],"%s",general_dihedral_types_registry[j][0]);
        sprintf((*dtmap)[j][1],"%s",general_dihedral_types_registry[j][1]); // BO
        sprintf((*dtmap)[j][2],"%s",general_dihedral_types_registry[j][2]);
        sprintf((*dtmap)[j][3],"%s",general_dihedral_types_registry[j][3]); // BO
        sprintf((*dtmap)[j][4],"%s",general_dihedral_types_registry[j][4]);
        sprintf((*dtmap)[j][5],"%s",general_dihedral_types_registry[j][5]); // BO
        sprintf((*dtmap)[j][6],"%s",general_dihedral_types_registry[j][6]);
    
    }
    *dcm=(int**)malloc(molecules*sizeof(int*));  // prealloc and initialize cm
    for(j=0;j<molecules;++j)(*dcm)[j]=(int*)malloc((*dtmap_rows)*sizeof(int));
    for(j=0;j<molecules;++j)
        for(k=0;k<*dtmap_rows;++k)
            (*dcm)[j][k]=0;
    // populate cm
    for(j=0;j<general_dihedrals;++j)
    {
        // locate in tmap
        for(k=0;k<*dtmap_rows;++k)
            if(
               (strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][0],(*dtmap)[k][0])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][1],(*dtmap)[k][1])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][2],(*dtmap)[k][2])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][3],(*dtmap)[k][3])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][4],(*dtmap)[k][4])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][5],(*dtmap)[k][5])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][6],(*dtmap)[k][6])==0
                )||
               (strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][0],(*dtmap)[k][6])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][1],(*dtmap)[k][5])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][2],(*dtmap)[k][4])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][3],(*dtmap)[k][3])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][4],(*dtmap)[k][2])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][5],(*dtmap)[k][1])==0 &&
                strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][6],(*dtmap)[k][0])==0
                )
               )
                break;
        (*dcm)[master_mol_array[general_dihedrals_registry[j][1]-1]-1][k]=1;
    }
    /*
     // console out
     for(j=0;j<*dtmap_rows;++j)printf("%s-%s-%s-%s\t",(*dtmap)[j][0],(*dtmap)[j][1],(*dtmap)[j][2],(*dtmap)[j][3]);printf("\n");
     for(j=0;j<molecules;++j){
     for(k=0;k<*dtmap_rows;++k)
     printf("%d\t",(*dcm)[j][k]);
     printf("\n");
     }
     */
    // __tmap__________
    // impropers
    
    *itmap_rows=general_improper_types;
    *itmap=(char***)malloc((*itmap_rows)*sizeof(char**)); // prealloc and populate tmap
    for(j=0;j<*itmap_rows;++j)(*itmap)[j]=(char**)malloc(4*sizeof(char*));
    for(j=0;j<*itmap_rows;++j)for(k=0;k<4;++k)(*itmap)[j][k]=(char*)malloc(cmax_length*sizeof(char));
    for(j=0;j<*itmap_rows;++j){
        sprintf((*itmap)[j][0],"%s",general_improper_types_registry[j][0]);
        sprintf((*itmap)[j][1],"%s",general_improper_types_registry[j][1]);
        sprintf((*itmap)[j][2],"%s",general_improper_types_registry[j][2]);
        sprintf((*itmap)[j][3],"%s",general_improper_types_registry[j][3]);}
    *icm=(int**)malloc(molecules*sizeof(int*));  // prealloc and initialize cm
    for(j=0;j<molecules;++j)(*icm)[j]=(int*)malloc((*itmap_rows)*sizeof(int));
    for(j=0;j<molecules;++j)
        for(k=0;k<*itmap_rows;++k)
            (*icm)[j][k]=0;
    // populate cm
    for(j=0;j<general_impropers;++j)
    {
        // locate in tmap
        for(k=0;k<*itmap_rows;++k)
            if(
               (strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][0],(*itmap)[k][0])==0 &&
                strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][1],(*itmap)[k][1])==0 &&
                strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][2],(*itmap)[k][2])==0 &&
                strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][3],(*itmap)[k][3])==0
                )||
               (strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][0],(*itmap)[k][0])==0 &&
                strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][1],(*itmap)[k][1])==0 &&
                strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][2],(*itmap)[k][3])==0 &&
                strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][3],(*itmap)[k][2])==0
                )||
               (strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][0],(*itmap)[k][0])==0 &&
                strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][1],(*itmap)[k][2])==0 &&
                strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][2],(*itmap)[k][1])==0 &&
                strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][3],(*itmap)[k][3])==0
                )||
               (strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][0],(*itmap)[k][0])==0 &&
                strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][1],(*itmap)[k][2])==0 &&
                strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][2],(*itmap)[k][3])==0 &&
                strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][3],(*itmap)[k][1])==0
                )||
               (strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][0],(*itmap)[k][0])==0 &&
                strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][1],(*itmap)[k][3])==0 &&
                strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][2],(*itmap)[k][2])==0 &&
                strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][3],(*itmap)[k][1])==0
                )||
               (strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][0],(*itmap)[k][0])==0 &&
                strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][1],(*itmap)[k][3])==0 &&
                strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][2],(*itmap)[k][1])==0 &&
                strcmp(general_improper_types_registry[general_impropers_registry[j][0]-1][3],(*itmap)[k][2])==0
                )
               
               
               )
                break;
        (*icm)[master_mol_array[general_impropers_registry[j][1]-1]-1][k]=1;
    }
    /*
     // console out
     for(j=0;j<*itmap_rows;++j)printf("%s-%s-%s-%s\t",(*itmap)[j][0],(*itmap)[j][1],(*itmap)[j][2],(*itmap)[j][3]);printf("\n");
     for(j=0;j<molecules;++j){
     for(k=0;k<*itmap_rows;++k)
     printf("%d\t",(*icm)[j][k]);
     printf("\n");
     }
     */
    // topo_boundaries
    // preallocate and initialize
    *topo_boundaries=(int**)malloc(molecules*sizeof(int*));
    for(i=0;i<molecules;++i)(*topo_boundaries)[i]=(int*)malloc(8*sizeof(int));
    for(i=0;i<molecules;++i)
        for(j=0;j<8;++j)
            (*topo_boundaries)[i][j]=-1;
    // populate
    for(i=0;i<molecules;++i)
    {
        k=0;for(j=0;j<general_bonds;++j)if(master_mol_array[general_bonds_registry[j][1]-1]==i+1){k=1;l=j;break;}
        if(k==1){
            (*topo_boundaries)[i][0]=l;
            for(j=l;j<general_bonds;++j)if(master_mol_array[general_bonds_registry[j][1]-1]!=i+1)break;(*topo_boundaries)[i][1]=j-1;
        }
        k=0;for(j=0;j<general_angles;++j)if(master_mol_array[general_angles_registry[j][1]-1]==i+1){k=1;l=j;break;}
        if(k==1){
            (*topo_boundaries)[i][2]=l;
            for(j=l;j<general_angles;++j)if(master_mol_array[general_angles_registry[j][1]-1]!=i+1)break;(*topo_boundaries)[i][3]=j-1;
        }
        k=0;for(j=0;j<general_dihedrals;++j)if(master_mol_array[general_dihedrals_registry[j][1]-1]==i+1){k=1;l=j;break;}
        if(k==1){
            (*topo_boundaries)[i][4]=l;
            for(j=l;j<general_dihedrals;++j)if(master_mol_array[general_dihedrals_registry[j][1]-1]!=i+1)break;(*topo_boundaries)[i][5]=j-1;
        }
        if(general_impropers>0)
        {
            k=0;for(j=0;j<general_impropers;++j)if(master_mol_array[general_impropers_registry[j][1]-1]==i+1){k=1;l=j;break;}
            if(k==1){
                (*topo_boundaries)[i][6]=l;
                for(j=l;j<general_impropers;++j)if(master_mol_array[general_impropers_registry[j][1]-1]!=i+1)break;(*topo_boundaries)[i][7]=j-1;
            }
        }
    }
    
     // console out
	if(verb==1){
     printf("\n$ TOPOLOGY (initial state)\n");

    
     printf("\n$ topo_boundaries (initial state):\n");
     for(i=0;i<molecules;++i)
     {
     printf("[%d]\t",i+1);
     for(j=0;j<8;++j)printf("%d\t",(*topo_boundaries)[i][j]);
     printf("\n");
     }
    
    //
    
    printf("\n$ general bonds registry (gbr):\n");
    for(j=0;j<general_bonds;++j)printf("[%d]\t%d\t%d\t%d\n",j+1,general_bonds_registry[j][0],general_bonds_registry[j][1],general_bonds_registry[j][2]);
    printf("$ general bond types registry (gbtr):\n");
    for(j=0;j<general_bond_types;++j)printf("[%d]\t(%s)--%s--(%s)\n",j+1,general_bond_types_registry[j][0],general_bond_types_registry[j][1],general_bond_types_registry[j][2]);
    printf("$ bond type mapping (btmap):\n");
    for(j=0;j<*btmap_rows;++j)printf("{%d}\t(%s)--%s--(%s)\n",j+1,(*btmap)[j][0],(*btmap)[j][1],(*btmap)[j][2]);
    printf("$ bond type counting matrix (bcm):\n");
    for(j=0;j<molecules;++j)
    {
        printf("{%d}\t",j+1);
        for(k=0;k<*btmap_rows;++k)printf("%d\t",(*bcm)[j][k]);
        printf("\n");
    }

    //
    
    printf("\n$ general angles registry (gar):\n");
    for(j=0;j<general_angles;++j)printf("[%d]\t%d\t%d\t%d\t%d\n",j+1,general_angles_registry[j][0],general_angles_registry[j][1],general_angles_registry[j][2],general_angles_registry[j][3]);
    printf("$ general angle types registry (gatr):\n");
    for(j=0;j<general_angle_types;++j)printf("[%d]\t(%s)--%s--(%s)--%s--(%s)\n",j+1,general_angle_types_registry[j][0],general_angle_types_registry[j][1],general_angle_types_registry[j][2],general_angle_types_registry[j][3],general_angle_types_registry[j][4]);
    printf("$ angle type mapping (atmap):\n");
    for(j=0;j<*atmap_rows;++j)printf("{%d}\t(%s)--%s--(%s)--%s--(%s)\n",j+1,(*atmap)[j][0],(*atmap)[j][1],(*atmap)[j][2],(*atmap)[j][3],(*atmap)[j][4]);
    printf("$ angle type counting matrix (acm):\n");
    for(j=0;j<molecules;++j)
    {
        printf("{%d}\t",j+1);
        for(k=0;k<*atmap_rows;++k)printf("%d\t",(*acm)[j][k]);
        printf("\n");
    }

    //
    
    printf("\n$ general dihedrals registry (gdr):\n");
    for(j=0;j<general_dihedrals;++j)printf("[%d]\t%d\t%d\t%d\t%d\t%d\n",j+1,general_dihedrals_registry[j][0],general_dihedrals_registry[j][1],general_dihedrals_registry[j][2],general_dihedrals_registry[j][3],general_dihedrals_registry[j][4]);
    printf("$ general dihedral types registry (gdtr):\n");
    for(j=0;j<general_dihedral_types;++j)printf("[%d]\t(%s)--%s--(%s)--%s--(%s)--%s--(%s)\n",j+1,general_dihedral_types_registry[j][0],general_dihedral_types_registry[j][1],general_dihedral_types_registry[j][2],general_dihedral_types_registry[j][3],general_dihedral_types_registry[j][4],general_dihedral_types_registry[j][5],general_dihedral_types_registry[j][6]);
    printf("$ dihedral type mapping (dtmap):\n");
    for(j=0;j<*dtmap_rows;++j)printf("{%d}\t(%s)--%s--(%s)--%s--(%s)--%s--(%s)\n",j+1,(*dtmap)[j][0],(*dtmap)[j][1],(*dtmap)[j][2],(*dtmap)[j][3],(*dtmap)[j][4],(*dtmap)[j][5],(*dtmap)[j][6]);
    printf("$ dihedral type counting matrix (dcm):\n");
    for(j=0;j<molecules;++j)
    {
        printf("{%d}\t",j+1);
        for(k=0;k<*dtmap_rows;++k)printf("%d\t",(*dcm)[j][k]);
        printf("\n");
    }

    //
    
    printf("\n$ general impropers registry (gdr):\n");
    for(j=0;j<general_impropers;++j)printf("[%d]\t%d\t%d\t%d\t%d\t%d\n",j+1,general_impropers_registry[j][0],general_impropers_registry[j][1],general_impropers_registry[j][2],general_impropers_registry[j][3],general_impropers_registry[j][4]);
    printf("$ general improper types registry (gdtr):\n");
    for(j=0;j<general_improper_types;++j)printf("[%d]\t%s\t%s\t%s\t%s\n",j+1,general_improper_types_registry[j][0],general_improper_types_registry[j][1],general_improper_types_registry[j][2],general_improper_types_registry[j][3]);
    printf("$ improper type mapping (dtmap):\n");
    for(j=0;j<*itmap_rows;++j)printf("{%d}\t%s\t%s\t%s\t%s\n",j+1,(*itmap)[j][0],(*itmap)[j][1],(*itmap)[j][2],(*itmap)[j][3]);
    printf("$ improper type counting matrix (dcm):\n");
    for(j=0;j<molecules;++j)
    {
        printf("{%d}\t",j+1);
        for(k=0;k<*itmap_rows;++k)printf("%d\t",(*icm)[j][k]);
        printf("\n");
    }
	}
    //getchar();
    
}
