
#include"builder.h"

void topo_molecule_alloc_populate(int i,int atoms_block_size,int **topo_boundaries,
                                  int general_bond_types,int general_angle_types,int general_dihedral_types,int general_improper_types,
                                  int **general_bonds_registry,int **general_angles_registry,int **general_dihedrals_registry,int **general_impropers_registry,
                                  char ***general_bond_types_registry,char ***general_angle_types_registry,char ***general_dihedral_types_registry,char ***general_improper_types_registry,
                                  int *molecule_atoms,int *molecule_bonds,int *molecule_angles,int *molecule_dihedrals,int *molecule_impropers,
                                  int *molecule_atom_types,int *molecule_bond_types,int *molecule_angle_types,int *molecule_dihedral_types,int *molecule_improper_types,
                                  int ***molecule_bonds_registry,int ***molecule_angles_registry,int ***molecule_dihedrals_registry,int ***molecule_impropers_registry,
                                  char ***molecule_species_registry,
                                  char ****molecule_bond_types_registry,char ****molecule_angle_types_registry,char ****molecule_dihedral_types_registry,char ****molecule_improper_types_registry,
                                  int *molecule_atoms_B,int *molecule_bonds_B,int *molecule_angles_B,int *molecule_dihedrals_B,int *molecule_impropers_B,
                                  int *molecule_atom_types_B,int *molecule_bond_types_B,int *molecule_angle_types_B,int *molecule_dihedral_types_B,int *molecule_improper_types_B,
                                  int ***molecule_bonds_registry_B,int ***molecule_angles_registry_B,int ***molecule_dihedrals_registry_B,int ***molecule_impropers_registry_B,
                                  char ***molecule_species_registry_B,
                                  char ****molecule_bond_types_registry_B,char ****molecule_angle_types_registry_B,char ****molecule_dihedral_types_registry_B,char ****molecule_improper_types_registry_B,
                                  
                                  int verb
                                  )
{
    char **dummy;
    int j,k,rescale;
    int *type_map,*type_map_index;
    
    // read topo_boundaries: set topology.molecule.counters.N.*
    *molecule_atoms=atoms_block_size;
    *molecule_bonds=topo_boundaries[i][1]-topo_boundaries[i][0]+1;
    *molecule_angles=topo_boundaries[i][3]-topo_boundaries[i][2]+1;
    if(topo_boundaries[i][4]==-1)
        *molecule_dihedrals=0;
    else
        *molecule_dihedrals=topo_boundaries[i][5]-topo_boundaries[i][4]+1;
    if(topo_boundaries[i][6]==-1)
        *molecule_impropers=0;
    else
        *molecule_impropers=topo_boundaries[i][7]-topo_boundaries[i][6]+1;
    
    // preallocate topology.molecule.arrays.*
    *molecule_bonds_registry=(int**)malloc(*molecule_bonds*sizeof(int*));for(j=0;j<*molecule_bonds;++j)(*molecule_bonds_registry)[j]=(int*)malloc(3*sizeof(int));
    *molecule_angles_registry=(int**)malloc(*molecule_angles*sizeof(int*));for(j=0;j<*molecule_angles;++j)(*molecule_angles_registry)[j]=(int*)malloc(4*sizeof(int));
    if(*molecule_dihedrals>0){
        *molecule_dihedrals_registry=(int**)malloc(*molecule_dihedrals*sizeof(int*));
        for(j=0;j<*molecule_dihedrals;++j)(*molecule_dihedrals_registry)[j]=(int*)malloc(5*sizeof(int));
    }
    if(*molecule_impropers>0){
        *molecule_impropers_registry=(int**)malloc(*molecule_impropers*sizeof(int*));
        for(j=0;j<*molecule_impropers;++j)(*molecule_impropers_registry)[j]=(int*)malloc(5*sizeof(int));
    }
    
    // populate topology.molecule.arrays.*
    
    // BONDS
    // these are two mapping arrays to normalize the type indices
    type_map=(int*)malloc(general_bond_types*sizeof(int));
    type_map_index=(int*)malloc(general_bond_types*sizeof(int));
    for(j=0;j<general_bond_types;++j)type_map[j]=0;
    // resolve bonds
    k=-1;
    for(j=topo_boundaries[i][0];j<=topo_boundaries[i][1];++j)
    {
        k=k+1;
        (*molecule_bonds_registry)[k][0]=general_bonds_registry[j][0];
        (*molecule_bonds_registry)[k][1]=general_bonds_registry[j][1];
        (*molecule_bonds_registry)[k][2]=general_bonds_registry[j][2];
        // mapping
        type_map[general_bonds_registry[j][0]-1]=1;
    }
    // rescale atomic ids to unity
    // find rescaling constant
    rescale=(*molecule_bonds_registry)[0][1];
    for(k=0;k<*molecule_bonds;++k)
    {
        if((*molecule_bonds_registry)[k][1]<rescale)rescale=(*molecule_bonds_registry)[k][1];
        if((*molecule_bonds_registry)[k][2]<rescale)rescale=(*molecule_bonds_registry)[k][2];
    }
    // apply rescaling
    for(k=0;k<*molecule_bonds;++k)
    {
        (*molecule_bonds_registry)[k][1]=(*molecule_bonds_registry)[k][1]-rescale+1;
        (*molecule_bonds_registry)[k][2]=(*molecule_bonds_registry)[k][2]-rescale+1;
    }
    // resolve bond types
    *molecule_bond_types=0;
    for(j=0;j<general_bond_types;++j)if(type_map[j]==1)*molecule_bond_types=*molecule_bond_types+1;
    // resolve mapping
    k=0;
    for(j=0;j<general_bond_types;++j)if(type_map[j]!=0){k=k+1;type_map_index[j]=k;}
    //for(j=0;j<general_bond_types;++j)if(type_map[j]!=0)printf("%d --> %d\n",j+1,type_map_index[j]);
    // console out
    //for(j=0;j<general_bond_types;++j)if(type_map[j]==1)printf("%d\t%s\t%s\n",j+1,general_bond_types_registry[j][0],general_bond_types_registry[j][1]);
    // preallocate bond types array
    *molecule_bond_types_registry=(char***)malloc(*molecule_bond_types*sizeof(char**));
    for(j=0;j<*molecule_bond_types;++j)(*molecule_bond_types_registry)[j]=(char**)malloc(3*sizeof(char*));
    for(j=0;j<*molecule_bond_types;++j)for(k=0;k<3;++k)(*molecule_bond_types_registry)[j][k]=(char*)malloc(sub_length*sizeof(char));
    // populate bond types array
    for(j=0;j<general_bond_types;++j)
        if(type_map[j]!=0)
        {
            sprintf((*molecule_bond_types_registry)[type_map_index[j]-1][0],"%s",general_bond_types_registry[j][0]);
            sprintf((*molecule_bond_types_registry)[type_map_index[j]-1][1],"%s",general_bond_types_registry[j][1]);    // BO
            sprintf((*molecule_bond_types_registry)[type_map_index[j]-1][2],"%s",general_bond_types_registry[j][2]);
        }
    // apply type mapping
    for(j=0;j<*molecule_bonds;++j)(*molecule_bonds_registry)[j][0]=type_map_index[(*molecule_bonds_registry)[j][0]-1];
    // free mapping arrays
    free(type_map);free(type_map_index);
    // console out
    //printf("bond types = %d\n",*molecule_bond_types);
    //for(j=0;j<*molecule_bond_types;++j)printf("[%d]\t%s\t%s\n",j+1,(*molecule_bond_types_registry)[j][0],(*molecule_bond_types_registry)[j][1]);
    
    //
    // SPECIES
    // resolve species from the bonding info
    //
    
    dummy=(char**)malloc((2*(*molecule_bond_types))*sizeof(char*));
    for(j=0;j<2*(*molecule_bond_types);++j)dummy[j]=(char*)malloc(sub_length*sizeof(char));
    k=-1;
    for(j=0;j<*molecule_bond_types;++j)
    {
        k=k+1;
        sprintf(dummy[k],"%s",(*molecule_bond_types_registry)[j][0]);
        k=k+1;
        sprintf(dummy[k],"%s",(*molecule_bond_types_registry)[j][2]);       // index changed!
    }
    
    for(j=0;j<2*(*molecule_bond_types)-1;++j)
    {
        for(k=j+1;k<2*(*molecule_bond_types);++k)
        {
            if(strcmp(dummy[j],dummy[k])==0)sprintf(dummy[k],"%s","___");
        }
    }
    
    //for(j=0;j<2*(*molecule_bond_types);++j)printf("%s\n",dummy[j]);
    
    *molecule_atom_types=2*(*molecule_bond_types);
    for(j=0;j<2*(*molecule_bond_types);++j)
        if(strcmp(dummy[j],"___")==0)*molecule_atom_types=*molecule_atom_types-1;
    
    *molecule_species_registry=(char**)malloc(*molecule_atom_types*sizeof(char*));
    for(j=0;j<*molecule_atom_types;++j)(*molecule_species_registry)[j]=(char*)malloc(sub_length*sizeof(char));
    
    k=-1;
    for(j=0;j<2*(*molecule_bond_types);++j)
        if(strcmp(dummy[j],"___")!=0)
        {
            k=k+1;
            sprintf((*molecule_species_registry)[k],"%s",dummy[j]);
        }
    
    for(j=0;j<2*(*molecule_bond_types);++j)free(dummy[j]);free(dummy);
    
    //printf("atom types = %d\n",*molecule_atom_types);
    //for(j=0;j<*molecule_atom_types;++j)printf("%s\n",(*molecule_species_registry)[j]);
    
    // ANGLES
    // these are two mapping arrays to normalize the type indices
    type_map=(int*)malloc(general_angle_types*sizeof(int));
    type_map_index=(int*)malloc(general_angle_types*sizeof(int));
    for(j=0;j<general_angle_types;++j)type_map[j]=0;
    // resolve angles
    k=-1;
    for(j=topo_boundaries[i][2];j<=topo_boundaries[i][3];++j)
    {
        k=k+1;
        (*molecule_angles_registry)[k][0]=general_angles_registry[j][0];
        (*molecule_angles_registry)[k][1]=general_angles_registry[j][1];
        (*molecule_angles_registry)[k][2]=general_angles_registry[j][2];
        (*molecule_angles_registry)[k][3]=general_angles_registry[j][3];
        // mapping
        type_map[general_angles_registry[j][0]-1]=1;
    }
    // rescale atomic ids to unity
    // apply rescaling
    for(k=0;k<*molecule_angles;++k)
    {
        (*molecule_angles_registry)[k][1]=(*molecule_angles_registry)[k][1]-rescale+1;
        (*molecule_angles_registry)[k][2]=(*molecule_angles_registry)[k][2]-rescale+1;
        (*molecule_angles_registry)[k][3]=(*molecule_angles_registry)[k][3]-rescale+1;
    }
    // resolve angle types
    *molecule_angle_types=0;
    for(j=0;j<general_angle_types;++j)if(type_map[j]==1)*molecule_angle_types=*molecule_angle_types+1;
    // resolve mapping
    k=0;
    for(j=0;j<general_angle_types;++j)if(type_map[j]!=0){k=k+1;type_map_index[j]=k;}
    // console out
    //for(j=0;j<general_angle_types;++j)if(type_map[j]==1)printf("%d\t%s\t%s\t%s\n",j+1,general_angle_types_registry[j][0],general_angle_types_registry[j][1],general_angle_types_registry[j][2]);
    // preallocate angle types array
    *molecule_angle_types_registry=(char***)malloc(*molecule_angle_types*sizeof(char**));
    for(j=0;j<*molecule_angle_types;++j)(*molecule_angle_types_registry)[j]=(char**)malloc(5*sizeof(char*));
    for(j=0;j<*molecule_angle_types;++j)for(k=0;k<5;++k)(*molecule_angle_types_registry)[j][k]=(char*)malloc(sub_length*sizeof(char));
    // populate angle types array
    for(j=0;j<general_angle_types;++j)
        if(type_map[j]!=0)
        {
            sprintf((*molecule_angle_types_registry)[type_map_index[j]-1][0],"%s",general_angle_types_registry[j][0]);
            sprintf((*molecule_angle_types_registry)[type_map_index[j]-1][1],"%s",general_angle_types_registry[j][1]);  // BO
            sprintf((*molecule_angle_types_registry)[type_map_index[j]-1][2],"%s",general_angle_types_registry[j][2]);
            sprintf((*molecule_angle_types_registry)[type_map_index[j]-1][3],"%s",general_angle_types_registry[j][3]);  // BO
            sprintf((*molecule_angle_types_registry)[type_map_index[j]-1][4],"%s",general_angle_types_registry[j][4]);
        }
    // apply type mapping
    for(j=0;j<*molecule_angles;++j)(*molecule_angles_registry)[j][0]=type_map_index[(*molecule_angles_registry)[j][0]-1];
    // free mapping arrays
    free(type_map);free(type_map_index);
    // console out
    //printf("angle types = %d\n",*molecule_angle_types);
    //for(j=0;j<*molecule_angle_types;++j)printf("[%d]\t%s\t%s\t%s\n",j+1,(*molecule_angle_types_registry)[j][0],(*molecule_angle_types_registry)[j][1],(*molecule_angle_types_registry)[j][2]);
    
    // DIHEDRALS
    
    // added...
    *molecule_dihedral_types=0;
    
    if(*molecule_dihedrals>0)
    {
        // these are two mapping arrays to normalize the type indices
        type_map=(int*)malloc(general_dihedral_types*sizeof(int));
        type_map_index=(int*)malloc(general_dihedral_types*sizeof(int));
        for(j=0;j<general_dihedral_types;++j)type_map[j]=0;
        // resolve dihedrals
        k=-1;
        for(j=topo_boundaries[i][4];j<=topo_boundaries[i][5];++j)
        {
            k=k+1;
            (*molecule_dihedrals_registry)[k][0]=general_dihedrals_registry[j][0];
            (*molecule_dihedrals_registry)[k][1]=general_dihedrals_registry[j][1];
            (*molecule_dihedrals_registry)[k][2]=general_dihedrals_registry[j][2];
            (*molecule_dihedrals_registry)[k][3]=general_dihedrals_registry[j][3];
            (*molecule_dihedrals_registry)[k][4]=general_dihedrals_registry[j][4];
            // mapping
            type_map[general_dihedrals_registry[j][0]-1]=1;
        }
        // rescale atomic ids to unity
        // apply rescaling
        for(k=0;k<*molecule_dihedrals;++k)
        {
            (*molecule_dihedrals_registry)[k][1]=(*molecule_dihedrals_registry)[k][1]-rescale+1;
            (*molecule_dihedrals_registry)[k][2]=(*molecule_dihedrals_registry)[k][2]-rescale+1;
            (*molecule_dihedrals_registry)[k][3]=(*molecule_dihedrals_registry)[k][3]-rescale+1;
            (*molecule_dihedrals_registry)[k][4]=(*molecule_dihedrals_registry)[k][4]-rescale+1;
        }
        // resolve dihedral types
        *molecule_dihedral_types=0;
        for(j=0;j<general_dihedral_types;++j)if(type_map[j]==1)*molecule_dihedral_types=*molecule_dihedral_types+1;
        // resolve mapping
        k=0;
        for(j=0;j<general_dihedral_types;++j)if(type_map[j]!=0){k=k+1;type_map_index[j]=k;}
        // console out
        //for(j=0;j<general_dihedral_types;++j)if(type_map[j]==1)printf("%d\t%s\t%s\t%s\t%s\n",j+1,general_dihedral_types_registry[j][0],general_dihedral_types_registry[j][1],general_dihedral_types_registry[j][2],general_dihedral_types_registry[j][3]);
        // preallocate dihedral types array
        *molecule_dihedral_types_registry=(char***)malloc(*molecule_dihedral_types*sizeof(char**));
        for(j=0;j<*molecule_dihedral_types;++j)(*molecule_dihedral_types_registry)[j]=(char**)malloc(7*sizeof(char*));
        for(j=0;j<*molecule_dihedral_types;++j)for(k=0;k<7;++k)(*molecule_dihedral_types_registry)[j][k]=(char*)malloc(sub_length*sizeof(char));
        // populate dihedral types array
        for(j=0;j<general_dihedral_types;++j)
            if(type_map[j]!=0)
            {
                sprintf((*molecule_dihedral_types_registry)[type_map_index[j]-1][0],"%s",general_dihedral_types_registry[j][0]);
                sprintf((*molecule_dihedral_types_registry)[type_map_index[j]-1][1],"%s",general_dihedral_types_registry[j][1]);    // BO
                sprintf((*molecule_dihedral_types_registry)[type_map_index[j]-1][2],"%s",general_dihedral_types_registry[j][2]);
                sprintf((*molecule_dihedral_types_registry)[type_map_index[j]-1][3],"%s",general_dihedral_types_registry[j][3]);    // BO
                sprintf((*molecule_dihedral_types_registry)[type_map_index[j]-1][4],"%s",general_dihedral_types_registry[j][4]);
                sprintf((*molecule_dihedral_types_registry)[type_map_index[j]-1][5],"%s",general_dihedral_types_registry[j][5]);    // BO
                sprintf((*molecule_dihedral_types_registry)[type_map_index[j]-1][6],"%s",general_dihedral_types_registry[j][6]);
            }
        // apply type mapping
        for(j=0;j<*molecule_dihedrals;++j)(*molecule_dihedrals_registry)[j][0]=type_map_index[(*molecule_dihedrals_registry)[j][0]-1];
        // free mapping arrays
        free(type_map);free(type_map_index);
        // console out
        //printf("dihedral types = %d\n",*molecule_dihedral_types);
        //for(j=0;j<*molecule_dihedral_types;++j)printf("[%d]\t%s\t%s\t%s\t%s\n",j+1,(*molecule_dihedral_types_registry)[j][0],(*molecule_dihedral_types_registry)[j][1],(*molecule_dihedral_types_registry)[j][2],(*molecule_dihedral_types_registry)[j][3]);
    }
    
    // IMPROPERS
    
    // added...
    *molecule_improper_types=0;
    
    if(*molecule_impropers>0)
    {
        // these are two mapping arrays to normalize the type indices
        type_map=(int*)malloc(general_improper_types*sizeof(int));
        type_map_index=(int*)malloc(general_improper_types*sizeof(int));
        for(j=0;j<general_improper_types;++j)type_map[j]=0;
        // resolve impropers
        k=-1;
        for(j=topo_boundaries[i][6];j<=topo_boundaries[i][7];++j)
        {
            k=k+1;
            (*molecule_impropers_registry)[k][0]=general_impropers_registry[j][0];
            (*molecule_impropers_registry)[k][1]=general_impropers_registry[j][1];
            (*molecule_impropers_registry)[k][2]=general_impropers_registry[j][2];
            (*molecule_impropers_registry)[k][3]=general_impropers_registry[j][3];
            (*molecule_impropers_registry)[k][4]=general_impropers_registry[j][4];
            // mapping
            type_map[general_impropers_registry[j][0]-1]=1;
        }
        // rescale atomic ids to unity
        // apply rescaling
        for(k=0;k<*molecule_impropers;++k)
        {
            (*molecule_impropers_registry)[k][1]=(*molecule_impropers_registry)[k][1]-rescale+1;
            (*molecule_impropers_registry)[k][2]=(*molecule_impropers_registry)[k][2]-rescale+1;
            (*molecule_impropers_registry)[k][3]=(*molecule_impropers_registry)[k][3]-rescale+1;
            (*molecule_impropers_registry)[k][4]=(*molecule_impropers_registry)[k][4]-rescale+1;
        }
        // resolve improper types
        *molecule_improper_types=0;
        for(j=0;j<general_improper_types;++j)if(type_map[j]==1)*molecule_improper_types=*molecule_improper_types+1;
        // resolve mapping
        k=0;
        for(j=0;j<general_improper_types;++j)if(type_map[j]!=0){k=k+1;type_map_index[j]=k;}
        // console out
        //for(j=0;j<general_improper_types;++j)if(type_map[j]==1)printf("%d\t%s\t%s\t%s\t%s\n",j+1,general_improper_types_registry[j][0],general_improper_types_registry[j][1],general_improper_types_registry[j][2],general_improper_types_registry[j][3]);
        // preallocate improper types array
        *molecule_improper_types_registry=(char***)malloc(*molecule_improper_types*sizeof(char**));
        for(j=0;j<*molecule_improper_types;++j)(*molecule_improper_types_registry)[j]=(char**)malloc(4*sizeof(char*));
        for(j=0;j<*molecule_improper_types;++j)for(k=0;k<4;++k)(*molecule_improper_types_registry)[j][k]=(char*)malloc(sub_length*sizeof(char));
        // populate improper types array
        for(j=0;j<general_improper_types;++j)
            if(type_map[j]!=0)
            {
                sprintf((*molecule_improper_types_registry)[type_map_index[j]-1][0],"%s",general_improper_types_registry[j][0]);
                sprintf((*molecule_improper_types_registry)[type_map_index[j]-1][1],"%s",general_improper_types_registry[j][1]);
                sprintf((*molecule_improper_types_registry)[type_map_index[j]-1][2],"%s",general_improper_types_registry[j][2]);
                sprintf((*molecule_improper_types_registry)[type_map_index[j]-1][3],"%s",general_improper_types_registry[j][3]);
            }
        // apply type mapping
        for(j=0;j<*molecule_impropers;++j)(*molecule_impropers_registry)[j][0]=type_map_index[(*molecule_impropers_registry)[j][0]-1];
        // free mapping arrays
        free(type_map);free(type_map_index);
        // console out
        //printf("improper types = %d\n",*molecule_improper_types);
        //for(j=0;j<*molecule_improper_types;++j)printf("[%d]\t%s\t%s\t%s\t%s\n",j+1,(*molecule_improper_types_registry)[j][0],(*molecule_improper_types_registry)[j][1],(*molecule_improper_types_registry)[j][2],(*molecule_improper_types_registry)[j][3]);
    }
    
    // print
    
    if(verb==1){
    
    printf("\n$ This is output from topo_molecule_alloc_populate() --> topo:molecule\n");
    
    printf("\n");
    
    printf("Number of atoms = %d\n",*molecule_atoms);
    printf("Number of bonds = %d\n",*molecule_bonds);
    printf("Number of angles = %d\n",*molecule_angles);
    printf("Number of dihedrals = %d\n",*molecule_dihedrals);
    printf("Number of impropers = %d\n",*molecule_impropers);
    
    printf("\n");
    
    printf("Number of atom types = %d\n",*molecule_atom_types);
    printf("Number of bond types = %d\n",*molecule_bond_types);
    printf("Number of angle types = %d\n",*molecule_angle_types);
    //if(*molecule_dihedrals>0)
    printf("Number of dihedral types = %d\n",*molecule_dihedral_types);
    //if(*molecule_impropers>0)
    printf("Number of improper types = %d\n",*molecule_improper_types);
    
    printf("\nAtom types\n\n");
    
    for(j=0;j<*molecule_atom_types;++j)printf("[%d]\t%s\n",j+1,(*molecule_species_registry)[j]);
    
    printf("\nBond types\n\n");
    
    for(j=0;j<*molecule_bond_types;++j)printf("[%d]\t(%s)--%s--(%s)\n",j+1,(*molecule_bond_types_registry)[j][0],(*molecule_bond_types_registry)[j][1],(*molecule_bond_types_registry)[j][2]);
    
    printf("\nAngle types\n\n");
    
    for(j=0;j<*molecule_angle_types;++j)printf("[%d]\t(%s)--%s--(%s)--%s--(%s)\n",j+1,(*molecule_angle_types_registry)[j][0],(*molecule_angle_types_registry)[j][1],(*molecule_angle_types_registry)[j][2],(*molecule_angle_types_registry)[j][3],(*molecule_angle_types_registry)[j][4]);
    //if(*molecule_dihedrals>0)
    //{
        
        printf("\nDihedral types\n\n");
        
        for(j=0;j<*molecule_dihedral_types;++j)printf("[%d]\t(%s)--%s--(%s)--%s--(%s)--%s--(%s)\n",j+1,(*molecule_dihedral_types_registry)[j][0],(*molecule_dihedral_types_registry)[j][1],(*molecule_dihedral_types_registry)[j][2],(*molecule_dihedral_types_registry)[j][3],(*molecule_dihedral_types_registry)[j][4],(*molecule_dihedral_types_registry)[j][5],(*molecule_dihedral_types_registry)[j][6]);
    //}
    //if(*molecule_impropers>0)
    //{
        
        printf("\nImproper types\n\n");
        
        for(j=0;j<*molecule_improper_types;++j)printf("[%d]\t%s\t%s\t%s\t%s\n",j+1,(*molecule_improper_types_registry)[j][0],(*molecule_improper_types_registry)[j][1],(*molecule_improper_types_registry)[j][2],(*molecule_improper_types_registry)[j][3]);
    //}
    
    printf("\nBonds\n\n");
    
    for(j=0;j<*molecule_bonds;++j)printf("[%d]\t%d\t%d\t%d\n",j+1,(*molecule_bonds_registry)[j][0],(*molecule_bonds_registry)[j][1],(*molecule_bonds_registry)[j][2]);
    
    printf("\nAngles\n\n");
    
    for(j=0;j<*molecule_angles;++j)printf("[%d]\t%d\t%d\t%d\t%d\n",j+1,(*molecule_angles_registry)[j][0],(*molecule_angles_registry)[j][1],(*molecule_angles_registry)[j][2],(*molecule_angles_registry)[j][3]);
    //if(*molecule_dihedrals>0)
    //{
        
        printf("\nDihedrals\n\n");
        
        for(j=0;j<*molecule_dihedrals;++j)printf("[%d]\t%d\t%d\t%d\t%d\t%d\n",j+1,(*molecule_dihedrals_registry)[j][0],(*molecule_dihedrals_registry)[j][1],(*molecule_dihedrals_registry)[j][2],(*molecule_dihedrals_registry)[j][3],(*molecule_dihedrals_registry)[j][4]);
    //}
    //if(*molecule_impropers>0)
    //{
        
        printf("\nImpropers\n\n");
        
        for(j=0;j<*molecule_impropers;++j)printf("[%d]\t%d\t%d\t%d\t%d\t%d\n",j+1,(*molecule_impropers_registry)[j][0],(*molecule_impropers_registry)[j][1],(*molecule_impropers_registry)[j][2],(*molecule_impropers_registry)[j][3],(*molecule_impropers_registry)[j][4]);
    //}
	}
    
    // backup: topology.molecule_B.*
    *molecule_atoms_B=atoms_block_size;
    *molecule_bonds_B=*molecule_bonds;
    *molecule_angles_B=*molecule_angles;
    *molecule_dihedrals_B=*molecule_dihedrals;
    *molecule_impropers_B=*molecule_impropers;
    *molecule_atom_types_B=*molecule_atom_types;
    *molecule_bond_types_B=*molecule_bond_types;
    *molecule_angle_types_B=*molecule_angle_types;
    *molecule_dihedral_types_B=*molecule_dihedral_types;
    *molecule_improper_types_B=*molecule_improper_types;
    // prealloc
    *molecule_bonds_registry_B=(int**)malloc(*molecule_bonds_B*sizeof(int*));for(j=0;j<*molecule_bonds_B;++j)(*molecule_bonds_registry_B)[j]=(int*)malloc(3*sizeof(int));
    *molecule_angles_registry_B=(int**)malloc(*molecule_angles_B*sizeof(int*));for(j=0;j<*molecule_angles_B;++j)(*molecule_angles_registry_B)[j]=(int*)malloc(4*sizeof(int));
    if(*molecule_dihedrals_B>0){
        *molecule_dihedrals_registry_B=(int**)malloc(*molecule_dihedrals_B*sizeof(int*));
        for(j=0;j<*molecule_dihedrals_B;++j)(*molecule_dihedrals_registry_B)[j]=(int*)malloc(5*sizeof(int));
    }
    if(*molecule_impropers_B>0){
        *molecule_impropers_registry_B=(int**)malloc(*molecule_impropers_B*sizeof(int*));
        for(j=0;j<*molecule_impropers_B;++j)(*molecule_impropers_registry_B)[j]=(int*)malloc(5*sizeof(int));
    }
    *molecule_species_registry_B=(char**)malloc(*molecule_atom_types_B*sizeof(char*));
    for(j=0;j<*molecule_atom_types_B;++j)(*molecule_species_registry_B)[j]=(char*)malloc(sub_length*sizeof(char));
    *molecule_bond_types_registry_B=(char***)malloc(*molecule_bond_types_B*sizeof(char**));
    for(j=0;j<*molecule_bond_types_B;++j)(*molecule_bond_types_registry_B)[j]=(char**)malloc(3*sizeof(char*));
    for(j=0;j<*molecule_bond_types_B;++j)for(k=0;k<3;++k)(*molecule_bond_types_registry_B)[j][k]=(char*)malloc(sub_length*sizeof(char));
    *molecule_angle_types_registry_B=(char***)malloc(*molecule_angle_types_B*sizeof(char**));
    for(j=0;j<*molecule_angle_types_B;++j)(*molecule_angle_types_registry_B)[j]=(char**)malloc(5*sizeof(char*));
    for(j=0;j<*molecule_angle_types_B;++j)for(k=0;k<5;++k)(*molecule_angle_types_registry_B)[j][k]=(char*)malloc(sub_length*sizeof(char));
    if(*molecule_dihedrals_B>0)
    {
        *molecule_dihedral_types_registry_B=(char***)malloc(*molecule_dihedral_types_B*sizeof(char**));
        for(j=0;j<*molecule_dihedral_types_B;++j)(*molecule_dihedral_types_registry_B)[j]=(char**)malloc(7*sizeof(char*));
        for(j=0;j<*molecule_dihedral_types_B;++j)for(k=0;k<7;++k)(*molecule_dihedral_types_registry_B)[j][k]=(char*)malloc(sub_length*sizeof(char));
    }
    if(*molecule_impropers_B>0)
    {
        *molecule_improper_types_registry_B=(char***)malloc(*molecule_improper_types_B*sizeof(char**));
        for(j=0;j<*molecule_improper_types_B;++j)(*molecule_improper_types_registry_B)[j]=(char**)malloc(4*sizeof(char*));
        for(j=0;j<*molecule_improper_types_B;++j)for(k=0;k<4;++k)(*molecule_improper_types_registry_B)[j][k]=(char*)malloc(sub_length*sizeof(char));
    }
    
    // populate
    for(j=0;j<*molecule_atom_types;++j)
        sprintf((*molecule_species_registry_B)[j],"%s",(*molecule_species_registry)[j]);
    for(j=0;j<*molecule_bond_types;++j){
        sprintf((*molecule_bond_types_registry_B)[j][0],"%s",(*molecule_bond_types_registry)[j][0]);
        sprintf((*molecule_bond_types_registry_B)[j][1],"%s",(*molecule_bond_types_registry)[j][1]);    // BO
        sprintf((*molecule_bond_types_registry_B)[j][2],"%s",(*molecule_bond_types_registry)[j][2]);
    }
    for(j=0;j<*molecule_angle_types;++j)
    {
        sprintf((*molecule_angle_types_registry_B)[j][0],"%s",(*molecule_angle_types_registry)[j][0]);
        sprintf((*molecule_angle_types_registry_B)[j][1],"%s",(*molecule_angle_types_registry)[j][1]);  // BO
        sprintf((*molecule_angle_types_registry_B)[j][2],"%s",(*molecule_angle_types_registry)[j][2]);
        sprintf((*molecule_angle_types_registry_B)[j][3],"%s",(*molecule_angle_types_registry)[j][3]);  // BO
        sprintf((*molecule_angle_types_registry_B)[j][4],"%s",(*molecule_angle_types_registry)[j][4]);
    }
    if(*molecule_dihedrals>0)
    {
        for(j=0;j<*molecule_dihedral_types;++j)
        {
            sprintf((*molecule_dihedral_types_registry_B)[j][0],"%s",(*molecule_dihedral_types_registry)[j][0]);
            sprintf((*molecule_dihedral_types_registry_B)[j][1],"%s",(*molecule_dihedral_types_registry)[j][1]);    // BO
            sprintf((*molecule_dihedral_types_registry_B)[j][2],"%s",(*molecule_dihedral_types_registry)[j][2]);
            sprintf((*molecule_dihedral_types_registry_B)[j][3],"%s",(*molecule_dihedral_types_registry)[j][3]);    // BO
            sprintf((*molecule_dihedral_types_registry_B)[j][4],"%s",(*molecule_dihedral_types_registry)[j][4]);
            sprintf((*molecule_dihedral_types_registry_B)[j][5],"%s",(*molecule_dihedral_types_registry)[j][5]);    // BO
            sprintf((*molecule_dihedral_types_registry_B)[j][6],"%s",(*molecule_dihedral_types_registry)[j][6]);
        }
        
    }
    if(*molecule_impropers>0)
    {
        for(j=0;j<*molecule_improper_types;++j)
        {
            sprintf((*molecule_improper_types_registry_B)[j][0],"%s",(*molecule_improper_types_registry)[j][0]);
            sprintf((*molecule_improper_types_registry_B)[j][1],"%s",(*molecule_improper_types_registry)[j][1]);
            sprintf((*molecule_improper_types_registry_B)[j][2],"%s",(*molecule_improper_types_registry)[j][2]);
            sprintf((*molecule_improper_types_registry_B)[j][3],"%s",(*molecule_improper_types_registry)[j][3]);
        }
    }
    for(j=0;j<*molecule_bonds;++j)
    {
        (*molecule_bonds_registry_B)[j][0]=(*molecule_bonds_registry)[j][0];
        (*molecule_bonds_registry_B)[j][1]=(*molecule_bonds_registry)[j][1];
        (*molecule_bonds_registry_B)[j][2]=(*molecule_bonds_registry)[j][2];
    }
    for(j=0;j<*molecule_angles;++j)
    {
        (*molecule_angles_registry_B)[j][0]=(*molecule_angles_registry)[j][0];
        (*molecule_angles_registry_B)[j][1]=(*molecule_angles_registry)[j][1];
        (*molecule_angles_registry_B)[j][2]=(*molecule_angles_registry)[j][2];
        (*molecule_angles_registry_B)[j][3]=(*molecule_angles_registry)[j][3];
    }
    if(*molecule_dihedrals>0)
    {
        for(j=0;j<*molecule_dihedrals;++j)
        {
            (*molecule_dihedrals_registry_B)[j][0]=(*molecule_dihedrals_registry)[j][0];
            (*molecule_dihedrals_registry_B)[j][1]=(*molecule_dihedrals_registry)[j][1];
            (*molecule_dihedrals_registry_B)[j][2]=(*molecule_dihedrals_registry)[j][2];
            (*molecule_dihedrals_registry_B)[j][3]=(*molecule_dihedrals_registry)[j][3];
            (*molecule_dihedrals_registry_B)[j][4]=(*molecule_dihedrals_registry)[j][4];
        }
    }
    if(*molecule_impropers>0)
    {
        for(j=0;j<*molecule_impropers;++j)
        {
            (*molecule_impropers_registry_B)[j][0]=(*molecule_impropers_registry)[j][0];
            (*molecule_impropers_registry_B)[j][1]=(*molecule_impropers_registry)[j][1];
            (*molecule_impropers_registry_B)[j][2]=(*molecule_impropers_registry)[j][2];
            (*molecule_impropers_registry_B)[j][3]=(*molecule_impropers_registry)[j][3];
            (*molecule_impropers_registry_B)[j][4]=(*molecule_impropers_registry)[j][4];
        }
    }
    
}
