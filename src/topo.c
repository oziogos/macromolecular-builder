
#include"builder.h"

//#include"time.h"

void topo(int atoms,char **species,int bonds,int *init_B1,int *init_B2,char **init_B_type,
          char ***global_species,char ****global_bond_type_array,char ****global_angle_type_array,char ****global_dihedral_type_array,char ****global_improper_type_array,
          int **global_neighbors,int ***global_bonds_array,int ***global_angles_array,int ***global_dihedrals_array,int ***global_impropers_array,
          int *angles_core,int *dihedrals_core,int *impropers_core,
          int *atom_types,int *bond_types,int *angle_types,int *dihedral_types,int *improper_types,
          int output_flag,char *current_folder,char *file_name)
{
    //clock_t c0,c1;
    
    FILE *fp;
    char file_path[cmax_length];
    
    //
    int atom_types_internal;
    int i,j,k,l,m,n;
    int angles;
    int *atomic_ID;
    char **unique_species;
    
    int flag,temp;
    int element,T1,T2,T3,dihedrals,init_dihedrals;
    int max_neighbors,*neighbors;
    
    int *B1,*B2;
    int **bonds_array;
    
    int *A1,*A2,*A3;
    
    int *init_D1,*init_D2,*init_D3,*init_D4;
    int *D1,*D2,*D3,*D4;
    
    int **bonds_ID;
    int **bonds_ID_final,bonds_ID_counter;
    int *bond_ID_export;
    
    int **angles_ID;
    int **angles_ID_final,angles_ID_counter;
    int *angles_ID_export;
    
    int **dihedrals_ID;
    int **dihedrals_ID_final,dihedrals_ID_counter;
    int *dihedrals_ID_export;
    //
    
    int imp_matrix[3][4]={{1,3,5,7},{1,5,7,3},{1,7,3,5}};
    int impropers;
    int *I1,*I2,*I3,*I4;
    int *I_type;
    int I_types;
    
    //
    
    char ***connectivity_matrix;
    int global_max_neighbors;
    
    char **bonds_ID_final_type;
    
    char **BO_array,word[cmax_length];
    
    char ***BO_record_matrix;
    
    int pos1,pos2,pos3;
    int pos1_i,pos2_i,pos1_j,pos2_j;
    
    char **BO_array_A1_A2,**BO_array_A2_A3;
    char **BO_A1_A2_record_matrix,**BO_A2_A3_record_matrix;
    
    char **BO_array_D1_D2,**BO_array_D2_D3,**BO_array_D3_D4;
    char **BO_D1_D2_record_matrix,**BO_D2_D3_record_matrix,**BO_D3_D4_record_matrix;
    
    // !!!
    *improper_types=0;
    
    //
    //c0=clock();
    
    // preallocate the arrays for the full bonding info
    B1=(int*)malloc(2*bonds*sizeof(int));
    B2=(int*)malloc(2*bonds*sizeof(int));
    
    BO_array=(char**)malloc(2*bonds*sizeof(char*));
    for(i=0;i<2*bonds;++i)BO_array[i]=(char*)malloc(sub_length*sizeof(char));
    
    // populate
    for(i=0;i<bonds;++i)
    {
        B1[i]=init_B1[i];
        B2[i]=init_B2[i];
        sprintf(BO_array[i],"%s",init_B_type[i]);
        B1[i+bonds]=init_B2[i];
        B2[i+bonds]=init_B1[i];
        sprintf(BO_array[i+bonds],"%s",init_B_type[i]);
    }
    // bubble sort with respect to the B1 array
    flag=1;
    while(flag==1)
    {
        flag=0;
        for(i=1;i<2*bonds;++i)
        {
            if(B1[i]<B1[i-1])
            {
                temp=B1[i-1];
                B1[i-1]=B1[i];
                B1[i]=temp;
                
                temp=B2[i-1];
                B2[i-1]=B2[i];
                B2[i]=temp;
                
                sprintf(word,"%s",BO_array[i-1]);
                sprintf(BO_array[i-1],"%s",BO_array[i]);
                sprintf(BO_array[i],"%s",word);
                
                flag=1;
            }
        }
    }
    
    // array to store the number of atoms bonded to every single atom
    neighbors=(int*)malloc(atoms*sizeof(int));
    // initialize
    for(i=0;i<atoms;++i){neighbors[i]=0;}
    // populate array
    for(i=0;i<2*bonds;++i){neighbors[B1[i]-1]=neighbors[B1[i]-1]+1;}
    // calculate the maximum number of bonded atoms - needed for the following preallocation
    max_neighbors=0;
    for(i=0;i<atoms;++i){if(neighbors[i]>=max_neighbors){max_neighbors=neighbors[i];}}
    
    // bonding info array preallocation and initialization
    bonds_array=(int**)malloc(atoms*sizeof(int*));
    for(i=0;i<atoms;++i){bonds_array[i]=(int*)malloc(max_neighbors*sizeof(int));}
    for(i=0;i<atoms;++i){for(j=0;j<max_neighbors;++j){bonds_array[i][j]=0;}}
    
    //
    BO_record_matrix=(char***)malloc(atoms*sizeof(char**));
    for(i=0;i<atoms;++i)BO_record_matrix[i]=(char**)malloc(max_neighbors*sizeof(char*));
    for(i=0;i<atoms;++i)for(j=0;j<max_neighbors;++j)BO_record_matrix[i][j]=(char*)malloc(sub_length*sizeof(char));
    for(i=0;i<atoms;++i)for(j=0;j<max_neighbors;++j)sprintf(BO_record_matrix[i][j],"%s","__");
    
    // populate bonds array -- the following routine works correctly provided that the B2 array is sorted with respect to the B1 array!
    k=-1;
    for(i=0;i<atoms;++i)
    {
        for(j=0;j<neighbors[i];++j)
        {
            k=k+1;
            bonds_array[i][j]=B2[k];
            sprintf(BO_record_matrix[i][j],"%s",BO_array[k]);
        }
    }
    //c1=clock();printf("%d\t%lf\t",atoms,(double)(c1-c0)/CLOCKS_PER_SEC);
    /*
    printf("1-2 neighbors:\n");
    for(i=0;i<atoms;++i){printf("[%d]\t",i+1);for(j=0;j<max_neighbors;++j)printf("%d\t",bonds_array[i][j]);printf("\n");}
    printf("bond orders:\n");
    for(i=0;i<atoms;++i){printf("[%d]\t",i+1);for(j=0;j<max_neighbors;++j)printf("%s\t",BO_record_matrix[i][j]);printf("\n");}
    //getchar();
    */
    
    //c0=clock();
    // calculate the total number of angles
    angles=0;
    for(i=0;i<atoms;++i){angles=angles+(neighbors[i]*(neighbors[i]-1))/2;}
    
    A1=(int*)malloc(angles*sizeof(int));
    A2=(int*)malloc(angles*sizeof(int));
    A3=(int*)malloc(angles*sizeof(int));
    
    // populate angles array
    l=-1;
    for(i=0;i<atoms;++i)
    {
        // the combinations double loop
        for(j=0;j<neighbors[i]-1;++j)
        {
            for(k=j+1;k<neighbors[i];++k)
            {
                l=l+1;
                A1[l]=bonds_array[i][j];
                A2[l]=i+1;
                A3[l]=bonds_array[i][k];
            }
        }
    }
    //c1=clock();printf("%lf\t",(double)(c1-c0)/CLOCKS_PER_SEC);
    
    //c0=clock();
    // calculate the initial number of dihedrals examining the branching info
    init_dihedrals=0;
    for(i=0;i<angles;++i)
    {
        // check the branching starting from the left side
        element=A1[i];
        T1=A2[i];
        T2=A3[i];
        for(j=0;j<neighbors[element-1];++j)
        {
            T3=bonds_array[element-1][j];
            if(T3!=T1 && T3!=T2)
            {
                init_dihedrals=init_dihedrals+1;
            }
        }
        // check the branching starting from the right side
        element=A3[i];
        T1=A1[i];
        T2=A2[i];
        for(j=0;j<neighbors[element-1];++j)
        {
            T3=bonds_array[element-1][j];
            if(T3!=T1 && T3!=T2)
            {
                init_dihedrals=init_dihedrals+1;
            }
        }
    }
    
    // preallocate the array for all dihedrals
    init_D1=(int*)malloc(init_dihedrals*sizeof(int));
    init_D2=(int*)malloc(init_dihedrals*sizeof(int));
    init_D3=(int*)malloc(init_dihedrals*sizeof(int));
    init_D4=(int*)malloc(init_dihedrals*sizeof(int));
    
    // populate
    k=0;
    for(i=0;i<angles;++i)
    {
        // left
        element=A1[i];
        T1=A2[i];
        T2=A3[i];
        for(j=0;j<neighbors[element-1];++j)
        {
            T3=bonds_array[element-1][j];
            if(T3!=T1 && T3!=T2)
            {
                if(T3<T2)	// (1)<(4)
                {
                    init_D1[k]=T3;
                    init_D2[k]=element;
                    init_D3[k]=T1;
                    init_D4[k]=T2;
                }
                else
                {
                    init_D1[k]=T2;
                    init_D2[k]=T1;
                    init_D3[k]=element;
                    init_D4[k]=T3;
                }
                k=k+1;
            }
        }
        // right
        element=A3[i];
        T1=A1[i];
        T2=A2[i];
        for(j=0;j<neighbors[element-1];++j)
        {
            T3=bonds_array[element-1][j];
            if(T3!=T1 && T3!=T2)
            {
                if(T1<T3)	// (1)<(4)
                {
                    init_D1[k]=T1;
                    init_D2[k]=T2;
                    init_D3[k]=element;
                    init_D4[k]=T3;
                }
                else
                {
                    init_D1[k]=T3;
                    init_D2[k]=element;
                    init_D3[k]=T2;
                    init_D4[k]=T1;
                }
                k=k+1;
            }
        }
    }
    
    // erase multiple entries
    for(i=0;i<init_dihedrals;++i)
    {
        for(j=0;j<i-1;++j)
        {
            if(init_D1[j]==init_D1[i] && init_D2[j]==init_D2[i] && init_D3[j]==init_D3[i] && init_D4[j]==init_D4[i])
            {
                init_D1[j]=0;
                init_D2[j]=0;
                init_D3[j]=0;
                init_D4[j]=0;
            }
        }
        for(j=i+1;j<init_dihedrals;++j)
        {
            if(init_D1[j]==init_D1[i] && init_D2[j]==init_D2[i] && init_D3[j]==init_D3[i] && init_D4[j]==init_D4[i])
            {
                init_D1[j]=0;
                init_D2[j]=0;
                init_D3[j]=0;
                init_D4[j]=0;
            }
        }
    }
    
    // bubble sort
    flag=1;
    while(flag==1)
    {
        flag=0;
        for(i=1;i<init_dihedrals;++i)
        {
            if(init_D1[i]<init_D1[i-1])
            {
                temp=init_D1[i-1];init_D1[i-1]=init_D1[i];init_D1[i]=temp;
                temp=init_D2[i-1];init_D2[i-1]=init_D2[i];init_D2[i]=temp;
                temp=init_D3[i-1];init_D3[i-1]=init_D3[i];init_D3[i]=temp;
                temp=init_D4[i-1];init_D4[i-1]=init_D4[i];init_D4[i]=temp;
                
                flag=1;
            }
        }
    }
    
    // count the final number of dihedrals based on non-zero entries
    dihedrals=0;
    for(i=0;i<init_dihedrals;++i){if(init_D1[i]!=0){dihedrals=dihedrals+1;}}
    
    // preallocate and populate final dihedrals array
    D1=(int*)malloc(dihedrals*sizeof(int));
    D2=(int*)malloc(dihedrals*sizeof(int));
    D3=(int*)malloc(dihedrals*sizeof(int));
    D4=(int*)malloc(dihedrals*sizeof(int));
    
    k=-1;
    for(i=0;i<init_dihedrals;++i)
    {
        if(init_D1[i]!=0)
        {
            k=k+1;
            D1[k]=init_D1[i];
            D2[k]=init_D2[i];
            D3[k]=init_D3[i];
            D4[k]=init_D4[i];
        }
    }
    //c1=clock();printf("%lf\t",(double)(c1-c0)/CLOCKS_PER_SEC);
    //------------------------------------------------------------------
    //c0=clock();
    // assign an integer counter starting from 1 for every different species present in the dat file
    atomic_ID=(int*)malloc(atoms*sizeof(int));			// the array for the new species representation based on integers
    unique_species=(char**)malloc(atoms*sizeof(char*));	// the array to store the unique species
    for (i=0;i<atoms;++i){unique_species[i]=(char*)malloc(sub_length*sizeof(char));}
    
    // identify unique secies
    sprintf(unique_species[0],"%s",species[0]);
    k=1;
    // check the rest species
    for(i=1;i<atoms;++i)
    {
        flag=0;
        for(j=0;j<k;++j)
        {
            if(strcmp(species[i],unique_species[j])==0){flag=1;}
        }
        if(flag==0)
        {
            sprintf(unique_species[k],"%s",species[i]);
            k=k+1;
        }
    }
    // save atom types
    atom_types_internal=k;
    
    // populate the atomic_ID array
    for(i=0;i<atoms;++i)
    {
        for(j=0;j<k;++j)
        {
            if(strcmp(species[i],unique_species[j])==0)
            {
                atomic_ID[i]=j+1;
            }
        }
    }
    //c1=clock();printf("%lf\t",(double)(c1-c0)/CLOCKS_PER_SEC);
    
    //..........................................................................
    //......................... bo refinements below ...........................
    //..........................................................................
    
    //--------------------------------------------------------------------------
    
    // Bonds
    
    //c0=clock();
    bonds_ID=(int**)malloc(bonds*sizeof(int*));		// 2D array to store the species counters
    for(i=0;i<bonds;++i){bonds_ID[i]=(int*)malloc(2*sizeof(int));}
    // populate
    for(i=0;i<bonds;++i)
    {
        bonds_ID[i][0]=atomic_ID[init_B1[i]-1];
        bonds_ID[i][1]=atomic_ID[init_B2[i]-1];
    }
    /*
    printf("atom type mapping:\n'1-1' equivalence between SYBYL species and incrementing integers starting from 1\n");
    for(i=0;i<atoms;++i)printf("%d\t%s\n",atomic_ID[i],species[i]);
    printf("---\n");
    printf("bonds with numerical atom types prior to type refinement:\n");
    for(i=0;i<bonds;++i)printf("%d\t\%d\t|\t\n",bonds_ID[i][0],bonds_ID[i][1]);
    */
    // erase multiple entries
    for(i=0;i<bonds;++i)
    {
        for(j=0;j<i-1;++j)
        {
            if(bonds_ID[j][0]==bonds_ID[i][0]&&bonds_ID[j][1]==bonds_ID[i][1]&&strcmp(init_B_type[j],init_B_type[i])==0)
            {
                bonds_ID[j][0]=0;
                bonds_ID[j][1]=0;
            }
            else if(bonds_ID[j][1]==bonds_ID[i][0]&&bonds_ID[j][0]==bonds_ID[i][1]&&strcmp(init_B_type[j],init_B_type[i])==0)
            {
                bonds_ID[j][0]=0;
                bonds_ID[j][1]=0;
            }
        }
        for(j=i+1;j<bonds;++j)
        {
            if(bonds_ID[j][0]==bonds_ID[i][0]&&bonds_ID[j][1]==bonds_ID[i][1]&&strcmp(init_B_type[j],init_B_type[i])==0)
            {
                bonds_ID[j][0]=0;
                bonds_ID[j][1]=0;
            }
            else if(bonds_ID[j][1]==bonds_ID[i][0]&&bonds_ID[j][0]==bonds_ID[i][1]&&strcmp(init_B_type[j],init_B_type[i])==0)
            {
                bonds_ID[j][0]=0;
                bonds_ID[j][1]=0;
            }
        }
    }
    /*
    printf("---\n");
    printf("bonds with numerical atom types after type refinement:\n");
    for(i=0;i<bonds;++i)printf("%d\t\%d\t|\t(%d)--%s--(%d)\n",bonds_ID[i][0],bonds_ID[i][1],init_B1[i],init_B_type[i],init_B2[i]);
    //getchar();
    */
    // calculate number of unique bonds
    bonds_ID_counter=0;
    for(i=0;i<bonds;++i)
    {
        if(bonds_ID[i][0]!=0)
        {
            bonds_ID_counter=bonds_ID_counter+1;
        }
    }
    
    // save the unique bonding info
    bonds_ID_final=(int**)malloc(bonds_ID_counter*sizeof(int*));	// preallocate
    for(i=0;i<bonds_ID_counter;++i){bonds_ID_final[i]=(int*)malloc(2*sizeof(int));}

    bonds_ID_final_type=(char**)malloc(bonds_ID_counter*sizeof(char*));
    for(i=0;i<bonds_ID_counter;++i)bonds_ID_final_type[i]=(char*)malloc(sub_length*sizeof(char));
    
    // populate
    j=-1;
    for(i=0;i<bonds;++i)
    {
        if(bonds_ID[i][0]!=0)
        {
            j=j+1;
            bonds_ID_final[j][0]=bonds_ID[i][0];
            bonds_ID_final[j][1]=bonds_ID[i][1];
            sprintf(bonds_ID_final_type[j],"%s",init_B_type[i]);
        }
    }
    /*
    printf("---\n");
    printf("resolved bond types using atom types and bond order:\n");
    printf("bond types:\n");
    for(i=0;i<bonds_ID_counter;++i)printf("[%d]\t(%s)--%s--(%s)\n",i+1,unique_species[bonds_ID_final[i][0]-1],bonds_ID_final_type[i],unique_species[bonds_ID_final[i][1]-1]);
    */
    // bond type counter
    bond_ID_export=(int*)malloc(bonds*sizeof(int));
    // populate based on comparison
    for(i=0;i<bonds;++i)
    {
        for(j=0;j<bonds_ID_counter;++j)
        {
            if(atomic_ID[init_B1[i]-1]==bonds_ID_final[j][0] && atomic_ID[init_B2[i]-1]==bonds_ID_final[j][1] && strcmp(init_B_type[i],bonds_ID_final_type[j])==0)
            {
                bond_ID_export[i]=j+1;
            }
            else if(atomic_ID[init_B1[i]-1]==bonds_ID_final[j][1] && atomic_ID[init_B2[i]-1]==bonds_ID_final[j][0] && strcmp(init_B_type[i],bonds_ID_final_type[j])==0)
            {
                bond_ID_export[i]=j+1;
            }
        }
    }
    /*
    printf("bonds:\n");
    for(i=0;i<bonds;++i)printf("[%d]\t%d\t%d\t%d\t|\t(%s)--%s--(%s)\n",i+1,bond_ID_export[i],init_B1[i],init_B2[i],unique_species[atomic_ID[init_B1[i]-1]-1],init_B_type[i],unique_species[atomic_ID[init_B2[i]-1]-1]);
    */
    //c1=clock();printf("%lf\t",(double)(c1-c0)/CLOCKS_PER_SEC);
    
    //--------------------------------------------------------------------------

    // Angles
    
    //c0=clock();
    angles_ID=(int**)malloc(angles*sizeof(int*));	// array to store the species ID
    for(i=0;i<angles;++i){angles_ID[i]=(int*)malloc(3*sizeof(int));}
    
    BO_array_A1_A2=(char**)malloc(angles*sizeof(char*));for(i=0;i<angles;++i)BO_array_A1_A2[i]=(char*)malloc(sub_length*sizeof(char));
    BO_array_A2_A3=(char**)malloc(angles*sizeof(char*));for(i=0;i<angles;++i)BO_array_A2_A3[i]=(char*)malloc(sub_length*sizeof(char));
    
    // populate
    for(i=0;i<angles;++i)
    {
        angles_ID[i][0]=atomic_ID[A1[i]-1];
        angles_ID[i][1]=atomic_ID[A2[i]-1];
        angles_ID[i][2]=atomic_ID[A3[i]-1];
    }
    /*
    printf("atom type mapping:\n'1-1' equivalence between SYBYL species and incrementing integers starting from 1\n");
    for(i=0;i<atoms;++i)printf("%d\t%s\n",atomic_ID[i],species[i]);
    printf("---\n");
    printf("angles with numerical atom types prior to type refinement:\n");
    for(i=0;i<angles;++i)printf("%d\t\%d\t%d\t|\t\n",angles_ID[i][0],angles_ID[i][1],angles_ID[i][2]);
    */
    // store bond orders to memory
    BO_A1_A2_record_matrix=(char**)malloc(angles*sizeof(char*));for(i=0;i<angles;++i)BO_A1_A2_record_matrix[i]=(char*)malloc(sub_length*sizeof(char));
    BO_A2_A3_record_matrix=(char**)malloc(angles*sizeof(char*));for(i=0;i<angles;++i)BO_A2_A3_record_matrix[i]=(char*)malloc(sub_length*sizeof(char));

    for(i=0;i<angles;++i)
    {
        for(k=0;k<max_neighbors;++k)
        {
            if(bonds_array[A1[i]-1][k]==A2[i]){pos1=k;break;}
        }
        for(k=0;k<max_neighbors;++k)
        {
            if(bonds_array[A2[i]-1][k]==A3[i]){pos2=k;break;}
        }
        sprintf(BO_A1_A2_record_matrix[i],"%s",BO_record_matrix[A1[i]-1][pos1]);
        sprintf(BO_A2_A3_record_matrix[i],"%s",BO_record_matrix[A2[i]-1][pos2]);
    }
    
    // erase multiple entries
    for(i=0;i<angles;++i)
    {
        for(j=0;j<i-1;++j)
        {
            if(angles_ID[j][0]==angles_ID[i][0] &&
               strcmp(BO_A1_A2_record_matrix[i],BO_A1_A2_record_matrix[j])==0 &&
               angles_ID[j][1]==angles_ID[i][1] &&
               strcmp(BO_A2_A3_record_matrix[i],BO_A2_A3_record_matrix[j])==0 &&
               angles_ID[j][2]==angles_ID[i][2]
               )
            {
                angles_ID[j][0]=0;
                angles_ID[j][1]=0;
                angles_ID[j][2]=0;
            }
            else if(angles_ID[j][2]==angles_ID[i][0] &&
                    strcmp(BO_A2_A3_record_matrix[i],BO_A1_A2_record_matrix[j])==0 &&
                    angles_ID[j][1]==angles_ID[i][1] &&
                    strcmp(BO_A1_A2_record_matrix[i],BO_A2_A3_record_matrix[j])==0 &&
                    angles_ID[j][0]==angles_ID[i][2]
                    )
            {
                angles_ID[j][0]=0;
                angles_ID[j][1]=0;
                angles_ID[j][2]=0;
            }
        }
        for(j=i+1;j<angles;++j)
        {
            if(angles_ID[j][0]==angles_ID[i][0] &&
               strcmp(BO_A1_A2_record_matrix[i],BO_A1_A2_record_matrix[j])==0 &&
               angles_ID[j][1]==angles_ID[i][1] &&
               strcmp(BO_A2_A3_record_matrix[i],BO_A2_A3_record_matrix[j])==0 &&
               angles_ID[j][2]==angles_ID[i][2]
               )
            {
                angles_ID[j][0]=0;
                angles_ID[j][1]=0;
                angles_ID[j][2]=0;
            }
            else if(angles_ID[j][2]==angles_ID[i][0] &&
                    strcmp(BO_A2_A3_record_matrix[i],BO_A1_A2_record_matrix[j])==0 &&
                    angles_ID[j][1]==angles_ID[i][1] &&
                    strcmp(BO_A1_A2_record_matrix[i],BO_A2_A3_record_matrix[j])==0 &&
                    angles_ID[j][0]==angles_ID[i][2]
                    )
            {
                angles_ID[j][0]=0;
                angles_ID[j][1]=0;
                angles_ID[j][2]=0;
            }
        }
    }
    /*
    printf("---\n");
    printf("angles with numerical atom types after type refinement:\n");
    for(i=0;i<angles;++i)
    {
        printf("%d\t\%d\t%d\t|\t(%d)--%s--(%d)--%s--(%d)\n",angles_ID[i][0],angles_ID[i][1],angles_ID[i][2],A1[i],BO_A1_A2_record_matrix[i],A2[i],BO_A2_A3_record_matrix[i],A3[i]);
    }
    */
    // unique angles counter
    angles_ID_counter=0;
    for(i=0;i<angles;++i)
    {
        if(angles_ID[i][0]!=0)
        {
            angles_ID_counter=angles_ID_counter+1;
        }
    }
    
    // store unique angles info
    angles_ID_final=(int**)malloc(angles_ID_counter*sizeof(int*));	// preallocate
    for(i=0;i<angles_ID_counter;++i){angles_ID_final[i]=(int*)malloc(3*sizeof(int));}
    
    // populate
    j=-1;
    for(i=0;i<angles;++i)
    {
        if(angles_ID[i][0]!=0)
        {
            j=j+1;
            angles_ID_final[j][0]=angles_ID[i][0];
            angles_ID_final[j][1]=angles_ID[i][1];
            angles_ID_final[j][2]=angles_ID[i][2];
            
            sprintf(BO_array_A1_A2[j],"%s",BO_A1_A2_record_matrix[i]);
            sprintf(BO_array_A2_A3[j],"%s",BO_A2_A3_record_matrix[i]);
            
        }
    }
    /*
    printf("---\n");
    printf("resolved angle types using atom types and bond order:\n");
    printf("angle types:\n");
    for(i=0;i<angles_ID_counter;++i)printf("[%d]\t(%s)--%s--(%s)--%s--(%s)\n",
                                           i+1,
                                           unique_species[angles_ID_final[i][0]-1],
                                           BO_array_A1_A2[i],
                                           unique_species[angles_ID_final[i][1]-1],
                                           BO_array_A2_A3[i],
                                           unique_species[angles_ID_final[i][2]-1]);
    */
    // angle type counter
    angles_ID_export=(int*)malloc(angles*sizeof(int));
    
    // populate
    for(i=0;i<angles;++i)
    {
        for(j=0;j<angles_ID_counter;++j)
        {
            if(atomic_ID[A1[i]-1]==angles_ID_final[j][0] &&
               strcmp(BO_A1_A2_record_matrix[i],BO_array_A1_A2[j])==0 &&
               atomic_ID[A2[i]-1]==angles_ID_final[j][1] &&
               strcmp(BO_A2_A3_record_matrix[i],BO_array_A2_A3[j])==0 &&
               atomic_ID[A3[i]-1]==angles_ID_final[j][2]
               )
            {
                angles_ID_export[i]=j+1;
            }
            else if(atomic_ID[A1[i]-1]==angles_ID_final[j][2] &&
                    strcmp(BO_A1_A2_record_matrix[i],BO_array_A2_A3[j])==0 &&
                    atomic_ID[A2[i]-1]==angles_ID_final[j][1] &&
                    strcmp(BO_A2_A3_record_matrix[i],BO_array_A1_A2[j])==0 &&
                    atomic_ID[A3[i]-1]==angles_ID_final[j][0]
                    )
            {
                angles_ID_export[i]=j+1;
            }
        }
    }
    /*
    printf("angles:\n");
    for(i=0;i<angles;++i){printf("[%d]\t%d\t%d\t%d\t%d\t|\t(%s)--%s--(%s)--%s--(%s)\n",
                                i+1,
                                angles_ID_export[i],
                                A1[i],
                                A2[i],
                                A3[i],
                                unique_species[atomic_ID[A1[i]-1]-1],
                                BO_A1_A2_record_matrix[i],
                                unique_species[atomic_ID[A2[i]-1]-1],
                                BO_A2_A3_record_matrix[i],
                                unique_species[atomic_ID[A3[i]-1]-1]
                                );
    }
    //getchar();
    */
    //c1=clock();printf("%lf\t",(double)(c1-c0)/CLOCKS_PER_SEC);
    
    //------------------------------------------------------------------
    
    // Dihedrals
    //c0=clock();
    dihedrals_ID=(int**)malloc(dihedrals*sizeof(int*));	// array to store the species ID
    for(i=0;i<dihedrals;++i){dihedrals_ID[i]=(int*)malloc(4*sizeof(int));}
    
    BO_array_D1_D2=(char**)malloc(dihedrals*sizeof(char*));for(i=0;i<dihedrals;++i)BO_array_D1_D2[i]=(char*)malloc(sub_length*sizeof(char));
    BO_array_D2_D3=(char**)malloc(dihedrals*sizeof(char*));for(i=0;i<dihedrals;++i)BO_array_D2_D3[i]=(char*)malloc(sub_length*sizeof(char));
    BO_array_D3_D4=(char**)malloc(dihedrals*sizeof(char*));for(i=0;i<dihedrals;++i)BO_array_D3_D4[i]=(char*)malloc(sub_length*sizeof(char));
    
    // populate
    for(i=0;i<dihedrals;++i)
    {
        dihedrals_ID[i][0]=atomic_ID[D1[i]-1];
        dihedrals_ID[i][1]=atomic_ID[D2[i]-1];
        dihedrals_ID[i][2]=atomic_ID[D3[i]-1];
        dihedrals_ID[i][3]=atomic_ID[D4[i]-1];
    }
    /*
    printf("atom type mapping:\n'1-1' equivalence between SYBYL species and incrementing integers starting from 1\n");
    for(i=0;i<atoms;++i)printf("%d\t%s\n",atomic_ID[i],species[i]);
    printf("---\n");
    printf("dihedrals with numerical atom types prior to type refinement:\n");
    for(i=0;i<dihedrals;++i)printf("%d\t\%d\t%d\t%d\t|\t\n",dihedrals_ID[i][0],dihedrals_ID[i][1],dihedrals_ID[i][2],dihedrals_ID[i][3]);
    */
    // store bond orders to memory
    BO_D1_D2_record_matrix=(char**)malloc(dihedrals*sizeof(char*));for(i=0;i<dihedrals;++i)BO_D1_D2_record_matrix[i]=(char*)malloc(sub_length*sizeof(char));
    BO_D2_D3_record_matrix=(char**)malloc(dihedrals*sizeof(char*));for(i=0;i<dihedrals;++i)BO_D2_D3_record_matrix[i]=(char*)malloc(sub_length*sizeof(char));
    BO_D3_D4_record_matrix=(char**)malloc(dihedrals*sizeof(char*));for(i=0;i<dihedrals;++i)BO_D3_D4_record_matrix[i]=(char*)malloc(sub_length*sizeof(char));
    
    for(i=0;i<dihedrals;++i)
    {
        for(k=0;k<max_neighbors;++k)
        {
            if(bonds_array[D1[i]-1][k]==D2[i]){pos1=k;break;}
        }
        for(k=0;k<max_neighbors;++k)
        {
            if(bonds_array[D2[i]-1][k]==D3[i]){pos2=k;break;}
        }
        for(k=0;k<max_neighbors;++k)
        {
            if(bonds_array[D3[i]-1][k]==D4[i]){pos3=k;break;}
        }
        sprintf(BO_D1_D2_record_matrix[i],"%s",BO_record_matrix[D1[i]-1][pos1]);
        sprintf(BO_D2_D3_record_matrix[i],"%s",BO_record_matrix[D2[i]-1][pos2]);
        sprintf(BO_D3_D4_record_matrix[i],"%s",BO_record_matrix[D3[i]-1][pos3]);
    }
    
    // erase multiple entries
    for(i=0;i<dihedrals;++i)
    {
        for(j=0;j<i-1;++j)
        {
            if(dihedrals_ID[j][0]==dihedrals_ID[i][0] &&
               strcmp(BO_D1_D2_record_matrix[i],BO_D1_D2_record_matrix[j])==0 &&
               dihedrals_ID[j][1]==dihedrals_ID[i][1] &&
               strcmp(BO_D2_D3_record_matrix[i],BO_D2_D3_record_matrix[j])==0 &&
               dihedrals_ID[j][2]==dihedrals_ID[i][2] &&
               strcmp(BO_D3_D4_record_matrix[i],BO_D3_D4_record_matrix[j])==0 &&
               dihedrals_ID[j][3]==dihedrals_ID[i][3]
               )
            {
                dihedrals_ID[j][0]=0;
                dihedrals_ID[j][1]=0;
                dihedrals_ID[j][2]=0;
                dihedrals_ID[j][3]=0;
            }
            
            else if(dihedrals_ID[j][3]==dihedrals_ID[i][0] &&
                    strcmp(BO_D3_D4_record_matrix[i],BO_D1_D2_record_matrix[j])==0 &&
                    dihedrals_ID[j][2]==dihedrals_ID[i][1] &&
                    strcmp(BO_D2_D3_record_matrix[i],BO_D2_D3_record_matrix[j])==0 &&
                    dihedrals_ID[j][1]==dihedrals_ID[i][2] &&
                    strcmp(BO_D1_D2_record_matrix[i],BO_D3_D4_record_matrix[j])==0 &&
                    dihedrals_ID[j][0]==dihedrals_ID[i][3]
                    )
            {
                dihedrals_ID[j][0]=0;
                dihedrals_ID[j][1]=0;
                dihedrals_ID[j][2]=0;
                dihedrals_ID[j][3]=0;
            }
            
        }
        for(j=i+1;j<dihedrals;++j)
        {
            if(dihedrals_ID[j][0]==dihedrals_ID[i][0] &&
               strcmp(BO_D1_D2_record_matrix[i],BO_D1_D2_record_matrix[j])==0 &&
               dihedrals_ID[j][1]==dihedrals_ID[i][1] &&
               strcmp(BO_D2_D3_record_matrix[i],BO_D2_D3_record_matrix[j])==0 &&
               dihedrals_ID[j][2]==dihedrals_ID[i][2] &&
               strcmp(BO_D3_D4_record_matrix[i],BO_D3_D4_record_matrix[j])==0 &&
               dihedrals_ID[j][3]==dihedrals_ID[i][3]
               )
            {
                dihedrals_ID[j][0]=0;
                dihedrals_ID[j][1]=0;
                dihedrals_ID[j][2]=0;
                dihedrals_ID[j][3]=0;
            }
            
            else if(dihedrals_ID[j][3]==dihedrals_ID[i][0] &&
                    strcmp(BO_D3_D4_record_matrix[i],BO_D1_D2_record_matrix[j])==0 &&
                    dihedrals_ID[j][2]==dihedrals_ID[i][1] &&
                    strcmp(BO_D2_D3_record_matrix[i],BO_D2_D3_record_matrix[j])==0 &&
                    dihedrals_ID[j][1]==dihedrals_ID[i][2] &&
                    strcmp(BO_D1_D2_record_matrix[i],BO_D3_D4_record_matrix[j])==0 &&
                    dihedrals_ID[j][0]==dihedrals_ID[i][3]
                    )
            {
                dihedrals_ID[j][0]=0;
                dihedrals_ID[j][1]=0;
                dihedrals_ID[j][2]=0;
                dihedrals_ID[j][3]=0;
            }
            
        }
    }
    /*
    printf("---\n");
    printf("dihedrals with numerical atom types after type refinement:\n");
    for(i=0;i<dihedrals;++i)
    {
        printf("%d\t\%d\t%d\t%d\t|\t(%d)--%s--(%d)--%s--(%d)--%s--(%d)\n",dihedrals_ID[i][0],dihedrals_ID[i][1],dihedrals_ID[i][2],dihedrals_ID[i][3],D1[i],BO_D1_D2_record_matrix[i],D2[i],BO_D2_D3_record_matrix[i],D3[i],BO_D3_D4_record_matrix[i],D4[i]);
    }
    */
    // unique dihedrals counter
    dihedrals_ID_counter=0;
    for(i=0;i<dihedrals;++i)
    {
        if(dihedrals_ID[i][0]!=0)
        {
            dihedrals_ID_counter=dihedrals_ID_counter+1;
        }
    }
    
    // store unique dihedrals info
    dihedrals_ID_final=(int**)malloc(dihedrals_ID_counter*sizeof(int*));	// preallocate
    for(i=0;i<dihedrals_ID_counter;++i){dihedrals_ID_final[i]=(int*)malloc(4*sizeof(int));}
    
    // populate
    j=-1;
    for(i=0;i<dihedrals;++i)
    {
        if(dihedrals_ID[i][0]!=0)
        {
            j=j+1;
            dihedrals_ID_final[j][0]=dihedrals_ID[i][0];
            dihedrals_ID_final[j][1]=dihedrals_ID[i][1];
            dihedrals_ID_final[j][2]=dihedrals_ID[i][2];
            dihedrals_ID_final[j][3]=dihedrals_ID[i][3];
            
            sprintf(BO_array_D1_D2[j],"%s",BO_D1_D2_record_matrix[i]);
            sprintf(BO_array_D2_D3[j],"%s",BO_D2_D3_record_matrix[i]);
            sprintf(BO_array_D3_D4[j],"%s",BO_D3_D4_record_matrix[i]);

        }
    }
    /*
    printf("---\n");
    printf("resolved dihedral types using atom types and bond order:\n");
    printf("dihedral types:\n");
    for(i=0;i<dihedrals_ID_counter;++i)printf("[%d]\t(%s)--%s--(%s)--%s--(%s)--%s--(%s)\n",
                                              i+1,
                                              unique_species[dihedrals_ID_final[i][0]-1],
                                              BO_array_D1_D2[i],
                                              unique_species[dihedrals_ID_final[i][1]-1],
                                              BO_array_D2_D3[i],
                                              unique_species[dihedrals_ID_final[i][2]-1],
                                              BO_array_D3_D4[i],
                                              unique_species[dihedrals_ID_final[i][3]-1]);
    */
    // dihedrals type counter
    dihedrals_ID_export=(int*)malloc(dihedrals*sizeof(int));
    
    // populate
    for(i=0;i<dihedrals;++i)
    {
        for(j=0;j<dihedrals_ID_counter;++j)
        {
            if(atomic_ID[D1[i]-1]==dihedrals_ID_final[j][0] &&
               strcmp(BO_D1_D2_record_matrix[i],BO_array_D1_D2[j])==0 &&
               atomic_ID[D2[i]-1]==dihedrals_ID_final[j][1] &&
               strcmp(BO_D2_D3_record_matrix[i],BO_array_D2_D3[j])==0 &&
               atomic_ID[D3[i]-1]==dihedrals_ID_final[j][2] &&
               strcmp(BO_D3_D4_record_matrix[i],BO_array_D3_D4[j])==0 &&
               atomic_ID[D4[i]-1]==dihedrals_ID_final[j][3]
               )
            {
                dihedrals_ID_export[i]=j+1;
            }
            if(atomic_ID[D1[i]-1]==dihedrals_ID_final[j][3] &&
               strcmp(BO_D1_D2_record_matrix[i],BO_array_D3_D4[j])==0 &&
               atomic_ID[D2[i]-1]==dihedrals_ID_final[j][2] &&
               strcmp(BO_D2_D3_record_matrix[i],BO_array_D2_D3[j])==0 &&
               atomic_ID[D3[i]-1]==dihedrals_ID_final[j][1] &&
               strcmp(BO_D3_D4_record_matrix[i],BO_array_D1_D2[j])==0 &&
               atomic_ID[D4[i]-1]==dihedrals_ID_final[j][0]
               )
            {
                dihedrals_ID_export[i]=j+1;
            }
        }
        
    }
    /*
    printf("dihedrals:\n");
    for(i=0;i<dihedrals;++i){printf("[%d]\t%d\t%d\t%d\t%d\t%d\t|\t(%s)--%s--(%s)--%s--(%s)--%s--(%s)\n",
                                    i+1,
                                    dihedrals_ID_export[i],
                                    D1[i],
                                    D2[i],
                                    D3[i],
                                    D4[i],
                                    unique_species[atomic_ID[D1[i]-1]-1],
                                    BO_D1_D2_record_matrix[i],
                                    unique_species[atomic_ID[D2[i]-1]-1],
                                    BO_D2_D3_record_matrix[i],
                                    unique_species[atomic_ID[D3[i]-1]-1],
                                    BO_D3_D4_record_matrix[i],
                                    unique_species[atomic_ID[D4[i]-1]-1]
                                    );}
    //getchar();
    */
    //c1=clock();printf("%lf\t",(double)(c1-c0)/CLOCKS_PER_SEC);
    //------------------------------------------------------------------
    // save to global variables
    //c0=clock();
    //global_bonds=bonds;
    *angles_core=angles;
    *dihedrals_core=dihedrals;
    //
    *atom_types=atom_types_internal;
    *global_species=(char**)malloc(atom_types_internal*sizeof(char*));
    for(i=0;i<atom_types_internal;++i)(*global_species)[i]=(char*)malloc(sub_length*sizeof(char));
    for(i=0;i<atom_types_internal;++i)sprintf((*global_species)[i],"%s",unique_species[i]);
    //##########################################################################
    // variables involved:
    // int bonds_ID_counter, int bonds, int **bonds_ID, int *bond_ID_export, char **unique_species
    //
    // bonds_ID_counter holds the number of unique bond types; resolved via int **bonds_ID
    // int **bonds_ID is a (bonds x 2) matrix.
    // int bonds is a read-only variable from the input that corresponds to the number of bonds.
    //
    // this is the final number of unique bond types; returned as output
    *bond_types=bonds_ID_counter;
    // allocate bond types output matrix
    *global_bond_type_array=(char***)malloc(bonds_ID_counter*sizeof(char**));
    for(i=0;i<bonds_ID_counter;++i)(*global_bond_type_array)[i]=(char**)malloc(3*sizeof(char*));
    for(i=0;i<bonds_ID_counter;++i)
        for(j=0;j<3;++j)
            (*global_bond_type_array)[i][j]=(char*)malloc(sub_length*sizeof(char));
    // populate the bond types output matrix
    // the uniqueness criterium is a non-zero value in the first column of bonds_ID
    for(i=0;i<bonds;++i)
    {
        if(bonds_ID[i][0]!=0)
        {
            sprintf((*global_bond_type_array)[bond_ID_export[i]-1][0],"%s",unique_species[bonds_ID[i][0]-1]);
            sprintf((*global_bond_type_array)[bond_ID_export[i]-1][1],"%s",init_B_type[i]);
            sprintf((*global_bond_type_array)[bond_ID_export[i]-1][2],"%s",unique_species[bonds_ID[i][1]-1]);
        }
    }
    //##########################################################################
    //
    *angle_types=angles_ID_counter;
    *global_angle_type_array=(char***)malloc(angles_ID_counter*sizeof(char**));
    for(i=0;i<angles_ID_counter;++i)(*global_angle_type_array)[i]=(char**)malloc(5*sizeof(char*));
    for(i=0;i<angles_ID_counter;++i)
        for(j=0;j<5;++j)
            (*global_angle_type_array)[i][j]=(char*)malloc(sub_length*sizeof(char));
    for(i=0;i<angles;++i)
    {
        if(angles_ID[i][0]!=0)
        {
            sprintf((*global_angle_type_array)[angles_ID_export[i]-1][0],"%s",unique_species[angles_ID[i][0]-1]);
            sprintf((*global_angle_type_array)[angles_ID_export[i]-1][1],"%s",BO_A1_A2_record_matrix[i]);
            sprintf((*global_angle_type_array)[angles_ID_export[i]-1][2],"%s",unique_species[angles_ID[i][1]-1]);
            sprintf((*global_angle_type_array)[angles_ID_export[i]-1][3],"%s",BO_A2_A3_record_matrix[i]);
            sprintf((*global_angle_type_array)[angles_ID_export[i]-1][4],"%s",unique_species[angles_ID[i][2]-1]);
        }
    }
    //
    *dihedral_types=dihedrals_ID_counter;
    *global_dihedral_type_array=(char***)malloc(dihedrals_ID_counter*sizeof(char**));
    for(i=0;i<dihedrals_ID_counter;++i)(*global_dihedral_type_array)[i]=(char**)malloc(7*sizeof(char*));
    for(i=0;i<dihedrals_ID_counter;++i)
        for(j=0;j<7;++j)
            (*global_dihedral_type_array)[i][j]=(char*)malloc(sub_length*sizeof(char));
    for(i=0;i<dihedrals;++i)
    {
        if(dihedrals_ID[i][0]!=0)
        {
            sprintf((*global_dihedral_type_array)[dihedrals_ID_export[i]-1][0],"%s",unique_species[dihedrals_ID[i][0]-1]);
            sprintf((*global_dihedral_type_array)[dihedrals_ID_export[i]-1][1],"%s",BO_D1_D2_record_matrix[i]);
            sprintf((*global_dihedral_type_array)[dihedrals_ID_export[i]-1][2],"%s",unique_species[dihedrals_ID[i][1]-1]);
            sprintf((*global_dihedral_type_array)[dihedrals_ID_export[i]-1][3],"%s",BO_D2_D3_record_matrix[i]);
            sprintf((*global_dihedral_type_array)[dihedrals_ID_export[i]-1][4],"%s",unique_species[dihedrals_ID[i][2]-1]);
            sprintf((*global_dihedral_type_array)[dihedrals_ID_export[i]-1][5],"%s",BO_D3_D4_record_matrix[i]);
            sprintf((*global_dihedral_type_array)[dihedrals_ID_export[i]-1][6],"%s",unique_species[dihedrals_ID[i][3]-1]);
        }
    }
    //
    *global_bonds_array=(int**)malloc(bonds*sizeof(int*));
    for(i=0;i<bonds;++i)(*global_bonds_array)[i]=(int*)malloc(3*sizeof(int));
    for(i=0;i<bonds;++i)
    {
        (*global_bonds_array)[i][0]=bond_ID_export[i];
        (*global_bonds_array)[i][1]=init_B1[i];
        (*global_bonds_array)[i][2]=init_B2[i];
    }
    //
    *global_angles_array=(int**)malloc(angles*sizeof(int*));
    for(i=0;i<angles;++i)(*global_angles_array)[i]=(int*)malloc(4*sizeof(int));
    for(i=0;i<angles;++i)
    {
        (*global_angles_array)[i][0]=angles_ID_export[i];
        (*global_angles_array)[i][1]=A1[i];
        (*global_angles_array)[i][2]=A2[i];
        (*global_angles_array)[i][3]=A3[i];
    }
    //
    *global_dihedrals_array=(int**)malloc(dihedrals*sizeof(int*));
    for(i=0;i<dihedrals;++i)(*global_dihedrals_array)[i]=(int*)malloc(5*sizeof(int));
    for(i=0;i<dihedrals;++i)
    {
        (*global_dihedrals_array)[i][0]=dihedrals_ID_export[i];
        (*global_dihedrals_array)[i][1]=D1[i];
        (*global_dihedrals_array)[i][2]=D2[i];
        (*global_dihedrals_array)[i][3]=D3[i];
        (*global_dihedrals_array)[i][4]=D4[i];
    }
    
    //
    global_max_neighbors=max_neighbors;
    connectivity_matrix=(char***)malloc(atoms*sizeof(char**));
    for(i=0;i<atoms;++i)connectivity_matrix[i]=(char**)malloc(2*(max_neighbors+1)*sizeof(char*));
    for(i=0;i<atoms;++i)
        for(j=0;j<2*(max_neighbors+1);++j)
            connectivity_matrix[i][j]=(char*)malloc(sub_length*sizeof(char));
    for(i=0;i<atoms;++i)
        for(j=0;j<2*(max_neighbors+1);++j)
        {
            connectivity_matrix[i][j][0]=' ';
            connectivity_matrix[i][j][0]='\0';
        }
    
    for(i=0;i<atoms;++i)
    {
        k=0;
        sprintf(connectivity_matrix[i][k],"%d",i+1);
        k=1;
        sprintf(connectivity_matrix[i][k],"%s",unique_species[atomic_ID[i]-1]);
        
        for(j=0;j<max_neighbors;++j)
        {
            
            if(bonds_array[i][j]!=0)
            {
                k=k+1;
                sprintf(connectivity_matrix[i][k],"%d",bonds_array[i][j]);
                k=k+1;
                sprintf(connectivity_matrix[i][k],"%s",unique_species[atomic_ID[bonds_array[i][j]-1]-1]);
            }
        }
    }
    //c1=clock();printf("%lf\t",(double)(c1-c0)/CLOCKS_PER_SEC);
    // ***
    *global_neighbors=(int*)malloc(atoms*sizeof(int));
    for(i=0;i<atoms;++i)(*global_neighbors)[i]=neighbors[i];
    
    //c0=clock();
    // impropers
    // count
    /*
     impropers=0;
     for(i=0;i<atoms;++i)
     {
     // exact 3 neighbors
     if((*global_neighbors)[i]==3)
     {
     for(n=0;n<3;++n)
     {
     impropers=impropers+1;
     }
     }
     }
     */
    impropers=0;
    for(i=0;i<atoms;++i)
    {
        // exact 3 neighbors
        if((*global_neighbors)[i]==3)
        {
            impropers=impropers+1;
        }
    }
    
    *impropers_core=impropers;
    if(impropers>0)
    {
        // prealloc
        I1=(int*)malloc(impropers*sizeof(int));
        I2=(int*)malloc(impropers*sizeof(int));
        I3=(int*)malloc(impropers*sizeof(int));
        I4=(int*)malloc(impropers*sizeof(int));
        // populate
        k=0;
        for(i=0;i<atoms;++i)
        {
            // exact 3 neighbors
            if((*global_neighbors)[i]==3)
            {
                /*
                 for(n=0;n<3;++n)
                 {
                 k=k+1;
                 m=0;j=imp_matrix[n][m]-1;sscanf(connectivity_matrix[i][j],"%d",&I1[k-1]);
                 m=1;j=imp_matrix[n][m]-1;sscanf(connectivity_matrix[i][j],"%d",&I2[k-1]);
                 m=2;j=imp_matrix[n][m]-1;sscanf(connectivity_matrix[i][j],"%d",&I3[k-1]);
                 m=3;j=imp_matrix[n][m]-1;sscanf(connectivity_matrix[i][j],"%d",&I4[k-1]);
                 }
                 */
                n=0;
                k=k+1;
                m=0;j=imp_matrix[n][m]-1;sscanf(connectivity_matrix[i][j],"%d",&I1[k-1]);
                m=1;j=imp_matrix[n][m]-1;sscanf(connectivity_matrix[i][j],"%d",&I2[k-1]);
                m=2;j=imp_matrix[n][m]-1;sscanf(connectivity_matrix[i][j],"%d",&I3[k-1]);
                m=3;j=imp_matrix[n][m]-1;sscanf(connectivity_matrix[i][j],"%d",&I4[k-1]);
                
            }
        }
        
        //for(i=0;i<impropers;++i)printf("[%d]\t%d\t%d\t%d\t%d\n",i+1,I1[i],I2[i],I3[i],I4[i]);
        
        //for(i=0;i<impropers;++i)printf("[%d]\t%s\t%s\t%s\t%s\n",i+1,species[I1[i]-1],species[I2[i]-1],species[I3[i]-1],species[I4[i]-1]);
        
        I_type=(int*)malloc(impropers*sizeof(int));
        for(i=0;i<impropers;++i)I_type[i]=0;
        for(i=0;i<impropers-1;++i)
        {
            for(j=i+1;j<impropers;++j)
            {
                if((strcmp(species[I1[i]-1],species[I1[j]-1])==0 && strcmp(species[I2[i]-1],species[I2[j]-1])==0 && strcmp(species[I3[i]-1],species[I3[j]-1])==0 && strcmp(species[I4[i]-1],species[I4[j]-1])==0) ||
                   (strcmp(species[I1[i]-1],species[I1[j]-1])==0 && strcmp(species[I2[i]-1],species[I2[j]-1])==0 && strcmp(species[I3[i]-1],species[I4[j]-1])==0 && strcmp(species[I4[i]-1],species[I3[j]-1])==0) ||
                   (strcmp(species[I1[i]-1],species[I1[j]-1])==0 && strcmp(species[I2[i]-1],species[I3[j]-1])==0 && strcmp(species[I3[i]-1],species[I2[j]-1])==0 && strcmp(species[I4[i]-1],species[I4[j]-1])==0) ||
                   (strcmp(species[I1[i]-1],species[I1[j]-1])==0 && strcmp(species[I2[i]-1],species[I3[j]-1])==0 && strcmp(species[I3[i]-1],species[I4[j]-1])==0 && strcmp(species[I4[i]-1],species[I2[j]-1])==0) ||
                   (strcmp(species[I1[i]-1],species[I1[j]-1])==0 && strcmp(species[I2[i]-1],species[I4[j]-1])==0 && strcmp(species[I3[i]-1],species[I2[j]-1])==0 && strcmp(species[I4[i]-1],species[I3[j]-1])==0) ||
                   (strcmp(species[I1[i]-1],species[I1[j]-1])==0 && strcmp(species[I2[i]-1],species[I4[j]-1])==0 && strcmp(species[I3[i]-1],species[I3[j]-1])==0 && strcmp(species[I4[i]-1],species[I2[j]-1])==0)
                   )
                {
                    I_type[j]=-1;
                }
            }
        }
        
        I_types=0;
        for(i=0;i<impropers;++i)if(I_type[i]!=-1)I_types=I_types+1;
        //for(i=0;i<impropers;++i)if(I_type[i]==0)printf("[%d]\t%s\t%s\t%s\t%s\n",i+1,species[I1[i]-1],species[I2[i]-1],species[I3[i]-1],species[I4[i]-1]);
        
        *improper_types=I_types;
        *global_improper_type_array=(char***)malloc(*improper_types*sizeof(char**));
        for(i=0;i<*improper_types;++i)(*global_improper_type_array)[i]=(char**)malloc(4*sizeof(char*));
        for(i=0;i<*improper_types;++i)
            for(j=0;j<4;++j)
                (*global_improper_type_array)[i][j]=(char*)malloc(sub_length*sizeof(char));
        
        j=-1;
        for(i=0;i<impropers;++i)
        {
            if(I_type[i]==0)
            {
                j=j+1;
                sprintf((*global_improper_type_array)[j][0],"%s",species[I1[i]-1]);
                sprintf((*global_improper_type_array)[j][1],"%s",species[I2[i]-1]);
                sprintf((*global_improper_type_array)[j][2],"%s",species[I3[i]-1]);
                sprintf((*global_improper_type_array)[j][3],"%s",species[I4[i]-1]);
            }
        }
        
        //reset
        for(i=0;i<impropers;++i)I_type[i]=0;
        for(i=0;i<impropers;++i)
        {
            for(j=0;j<*improper_types;++j)
            {
                if((strcmp(species[I1[i]-1],(*global_improper_type_array)[j][0])==0 &&
                    strcmp(species[I2[i]-1],(*global_improper_type_array)[j][1])==0 &&
                    strcmp(species[I3[i]-1],(*global_improper_type_array)[j][2])==0 &&
                    strcmp(species[I4[i]-1],(*global_improper_type_array)[j][3])==0) ||
                   
                   (strcmp(species[I1[i]-1],(*global_improper_type_array)[j][0])==0 &&
                    strcmp(species[I2[i]-1],(*global_improper_type_array)[j][1])==0 &&
                    strcmp(species[I3[i]-1],(*global_improper_type_array)[j][3])==0 &&
                    strcmp(species[I4[i]-1],(*global_improper_type_array)[j][2])==0) ||
                   
                   (strcmp(species[I1[i]-1],(*global_improper_type_array)[j][0])==0 &&
                    strcmp(species[I2[i]-1],(*global_improper_type_array)[j][2])==0 &&
                    strcmp(species[I3[i]-1],(*global_improper_type_array)[j][1])==0 &&
                    strcmp(species[I4[i]-1],(*global_improper_type_array)[j][3])==0) ||
                   
                   (strcmp(species[I1[i]-1],(*global_improper_type_array)[j][0])==0 &&
                    strcmp(species[I2[i]-1],(*global_improper_type_array)[j][2])==0 &&
                    strcmp(species[I3[i]-1],(*global_improper_type_array)[j][3])==0 &&
                    strcmp(species[I4[i]-1],(*global_improper_type_array)[j][1])==0) ||
                   
                   (strcmp(species[I1[i]-1],(*global_improper_type_array)[j][0])==0 &&
                    strcmp(species[I2[i]-1],(*global_improper_type_array)[j][3])==0 &&
                    strcmp(species[I3[i]-1],(*global_improper_type_array)[j][1])==0 &&
                    strcmp(species[I4[i]-1],(*global_improper_type_array)[j][2])==0) ||
                   
                   (strcmp(species[I1[i]-1],(*global_improper_type_array)[j][0])==0 &&
                    strcmp(species[I2[i]-1],(*global_improper_type_array)[j][3])==0 &&
                    strcmp(species[I3[i]-1],(*global_improper_type_array)[j][2])==0 &&
                    strcmp(species[I4[i]-1],(*global_improper_type_array)[j][1])==0)
                   )
                {
                    I_type[i]=j+1;
                }
            }
        }
        
        // save to global
        //*impropers_core=impropers;
        *global_impropers_array=(int**)malloc(*impropers_core*sizeof(int*));
        for(i=0;i<*impropers_core;++i)(*global_impropers_array)[i]=(int*)malloc(5*sizeof(int));
        for(i=0;i<*impropers_core;++i)
        {
            (*global_impropers_array)[i][0]=I_type[i];
            (*global_impropers_array)[i][1]=I1[i];
            (*global_impropers_array)[i][2]=I2[i];
            (*global_impropers_array)[i][3]=I3[i];
            (*global_impropers_array)[i][4]=I4[i];
        }
        
        
        //for(i=0;i<impropers;++i)printf("[%d]\t%d\t%d\t%d\t%d\t%d\n",i+1,I_type[i],I1[i],I2[i],I3[i],I4[i]);
    }
    /*else
     {
     *global_improper_type_array=(char***)malloc(*improper_types*sizeof(char**));
     for(i=0;i<*improper_types;++i)(*global_improper_type_array)[i]=(char**)malloc(4*sizeof(char*));
     for(i=0;i<*improper_types;++i)
     for(j=0;j<4;++j)
     (*global_improper_type_array)[i][j]=(char*)malloc(sub_length*sizeof(char));
     }*/
    //c1=clock();printf("%lf\n",(double)(c1-c0)/CLOCKS_PER_SEC);
    
    
    
    // output
    if(output_flag==1)
    {
        sprintf(file_path,"%s/%s",current_folder,file_name);
        fp=fopen(file_path,"w+");
        
        fprintf(fp,"\n# This is output from topo()\n");
        
        fprintf(fp,"\n");
        
        fprintf(fp,"Number of atoms = %d\n",atoms);
        fprintf(fp,"Number of bonds = %d\n",bonds);
        fprintf(fp,"Number of angles = %d\n",*angles_core);
        fprintf(fp,"Number of dihedrals = %d\n",*dihedrals_core);
        fprintf(fp,"Number of impropers = %d\n",*impropers_core);
        
        fprintf(fp,"\n");
        
        fprintf(fp,"Number of atom types = %d\n",*atom_types);
        fprintf(fp,"Number of bond types = %d\n",*bond_types);
        fprintf(fp,"Number of angle types = %d\n",*angle_types);
        fprintf(fp,"Number of dihedral types = %d\n",*dihedral_types);
        fprintf(fp,"Number of improper types = %d\n",*improper_types);
        
        fprintf(fp,"\nAtom types\n\n");
        
        for(i=0;i<*atom_types;++i)
            fprintf(fp,"[%d]\t%s\n",i+1,(*global_species)[i]);
        
        fprintf(fp,"\nBond types\n\n");
        
        for(i=0;i<*bond_types;++i)
            fprintf(fp,"[%d]\t(%s)--%s--(%s)\n",i+1,(*global_bond_type_array)[i][0],(*global_bond_type_array)[i][1],(*global_bond_type_array)[i][2]);
        
        fprintf(fp,"\nAngle types\n\n");
        
        for(i=0;i<*angle_types;++i)
            fprintf(fp,"[%d]\t(%s)--%s--(%s)--%s--(%s)\n",i+1,(*global_angle_type_array)[i][0],(*global_angle_type_array)[i][1],(*global_angle_type_array)[i][2],(*global_angle_type_array)[i][3],(*global_angle_type_array)[i][4]);
        
        fprintf(fp,"\nDihedral types\n\n");
        
        for(i=0;i<*dihedral_types;++i)
            fprintf(fp,"[%d]\t(%s)--%s--(%s)--%s--(%s)--%s--(%s)\n",i+1,(*global_dihedral_type_array)[i][0],(*global_dihedral_type_array)[i][1],(*global_dihedral_type_array)[i][2],(*global_dihedral_type_array)[i][3],(*global_dihedral_type_array)[i][4],(*global_dihedral_type_array)[i][5],(*global_dihedral_type_array)[i][6]);
        
        fprintf(fp,"\nImproper types\n\n");
        
        for(i=0;i<*improper_types;++i)
            fprintf(fp,"[%d]\t%s\t%s\t%s\t%s\n",i+1,(*global_improper_type_array)[i][0],(*global_improper_type_array)[i][1],(*global_improper_type_array)[i][2],(*global_improper_type_array)[i][3]);
        
        fprintf(fp,"\nBonds\n\n");
        
        for(i=0;i<bonds;++i)
            fprintf(fp,"[%d]\t%d\t%d\t%d\n",i+1,(*global_bonds_array)[i][0],(*global_bonds_array)[i][1],(*global_bonds_array)[i][2]);
        
        fprintf(fp,"\nAngles\n\n");
        
        for(i=0;i<*angles_core;++i)
            fprintf(fp,"[%d]\t%d\t%d\t%d\t%d\n",i+1,(*global_angles_array)[i][0],(*global_angles_array)[i][1],(*global_angles_array)[i][2],(*global_angles_array)[i][3]);
        
        fprintf(fp,"\nDihedrals\n\n");
        
        for(i=0;i<*dihedrals_core;++i)
            fprintf(fp,"[%d]\t%d\t%d\t%d\t%d\t%d\n",i+1,(*global_dihedrals_array)[i][0],(*global_dihedrals_array)[i][1],(*global_dihedrals_array)[i][2],(*global_dihedrals_array)[i][3],(*global_dihedrals_array)[i][4]);
        
        fprintf(fp,"\nImpropers\n\n");
        
        for(i=0;i<*impropers_core;++i)
            fprintf(fp,"[%d]\t%d\t%d\t%d\t%d\t%d\n",i+1,(*global_impropers_array)[i][0],(*global_impropers_array)[i][1],(*global_impropers_array)[i][2],(*global_impropers_array)[i][3],(*global_impropers_array)[i][4]);
        
        fclose(fp);
    }
    else if(output_flag==2)
    {
        
        printf("\n# This is output from topo()\n");
        
        printf("\n");
        
        printf("Number of atoms = %d\n",atoms);
        printf("Number of bonds = %d\n",bonds);
        printf("Number of angles = %d\n",*angles_core);
        printf("Number of dihedrals = %d\n",*dihedrals_core);
        printf("Number of impropers = %d\n",*impropers_core);
        
        printf("\n");
        
        printf("Number of atom types = %d\n",*atom_types);
        printf("Number of bond types = %d\n",*bond_types);
        printf("Number of angle types = %d\n",*angle_types);
        printf("Number of dihedral types = %d\n",*dihedral_types);
        printf("Number of improper types = %d\n",*improper_types);
        
        printf("\nAtom types\n\n");
        
        for(i=0;i<*atom_types;++i)
            printf("[%d]\t%s\n",i+1,(*global_species)[i]);
        
        printf("\nBond types\n\n");
        
        for(i=0;i<*bond_types;++i)
            printf("[%d]\t(%s)--%s--(%s)\n",i+1,(*global_bond_type_array)[i][0],(*global_bond_type_array)[i][1],(*global_bond_type_array)[i][2]);
        
        printf("\nAngle types\n\n");
        
        for(i=0;i<*angle_types;++i)
            printf("[%d]\t(%s)--%s--(%s)--%s--(%s)\n",i+1,(*global_angle_type_array)[i][0],(*global_angle_type_array)[i][1],(*global_angle_type_array)[i][2],(*global_angle_type_array)[i][3],(*global_angle_type_array)[i][4]);
        
        printf("\nDihedral types\n\n");
        
        for(i=0;i<*dihedral_types;++i)
            printf("[%d]\t(%s)--%s--(%s)--%s--(%s)--%s--(%s)\n",i+1,(*global_dihedral_type_array)[i][0],(*global_dihedral_type_array)[i][1],(*global_dihedral_type_array)[i][2],(*global_dihedral_type_array)[i][3],(*global_dihedral_type_array)[i][4],(*global_dihedral_type_array)[i][5],(*global_dihedral_type_array)[i][6]);
        
        printf("\nImproper types\n\n");
        
        for(i=0;i<*improper_types;++i)
            printf("[%d]\t%s\t%s\t%s\t%s\n",i+1,(*global_improper_type_array)[i][0],(*global_improper_type_array)[i][1],(*global_improper_type_array)[i][2],(*global_improper_type_array)[i][3]);
        
        printf("\nBonds\n\n");
        
        for(i=0;i<bonds;++i)
            printf("[%d]\t%d\t%d\t%d\n",i+1,(*global_bonds_array)[i][0],(*global_bonds_array)[i][1],(*global_bonds_array)[i][2]);
        
        printf("\nAngles\n\n");
        
        for(i=0;i<*angles_core;++i)
            printf("[%d]\t%d\t%d\t%d\t%d\n",i+1,(*global_angles_array)[i][0],(*global_angles_array)[i][1],(*global_angles_array)[i][2],(*global_angles_array)[i][3]);
        
        printf("\nDihedrals\n\n");
        
        for(i=0;i<*dihedrals_core;++i)
            printf("[%d]\t%d\t%d\t%d\t%d\t%d\n",i+1,(*global_dihedrals_array)[i][0],(*global_dihedrals_array)[i][1],(*global_dihedrals_array)[i][2],(*global_dihedrals_array)[i][3],(*global_dihedrals_array)[i][4]);
        
        printf("\nImpropers\n\n");
        
        for(i=0;i<*impropers_core;++i)
            printf("[%d]\t%d\t%d\t%d\t%d\t%d\n",i+1,(*global_impropers_array)[i][0],(*global_impropers_array)[i][1],(*global_impropers_array)[i][2],(*global_impropers_array)[i][3],(*global_impropers_array)[i][4]);
        
    }
    
    
    
    //------------------------------------------------------------------
    // free memory

    for (i=0;i<dihedrals_ID_counter;++i){free(dihedrals_ID_final[i]);}free(dihedrals_ID_final);
    for (i=0;i<dihedrals;++i){free(dihedrals_ID[i]);}free(dihedrals_ID);
    for (i=0;i<angles_ID_counter;++i){free(angles_ID_final[i]);}free(angles_ID_final);
    for (i=0;i<angles;++i){free(angles_ID[i]);}free(angles_ID);
    for (i=0;i<bonds_ID_counter;++i){free(bonds_ID_final[i]);}free(bonds_ID_final);
    for (i=0;i<bonds_ID_counter;++i){free(bonds_ID_final_type[i]);}free(bonds_ID_final_type);
    for (i=0;i<bonds;++i){free(bonds_ID[i]);}free(bonds_ID);
    for (i=0;i<atoms;++i){free(bonds_array[i]);}free(bonds_array);
    for (i=0;i<atoms;++i){free(unique_species[i]);}free(unique_species);
    free(neighbors);
    free(bond_ID_export);free(angles_ID_export);free(dihedrals_ID_export);
    free(atomic_ID);
    free(B1);free(B2);
    

    
    for(i=0;i<2*bonds;++i)free(BO_array[i]);free(BO_array);
    
    
    free(A1);free(A2);free(A3);
    free(init_D1);free(init_D2);free(init_D3);free(init_D4);free(D1);free(D2);free(D3);free(D4);
    
    if(impropers>0){
        free(I1);free(I2);free(I3);free(I4);free(I_type);
    }
    //
    //
    //
    
    /*
     for(i=0;i<atoms;++i)
     {
     printf("[%d]\t",i+1);
     for(j=0;j<2*(global_max_neighbors+1);++j)printf("%s\t",connectivity_matrix[i][j]);
     printf("\n");
     }
     getchar();
     */
    

    
    for(i=0;i<atoms;++i)
        for(j=0;j<2*(global_max_neighbors+1);++j)
            free(connectivity_matrix[i][j]);
    for(i=0;i<atoms;++i)free(connectivity_matrix[i]);
    free(connectivity_matrix);

    for(i=0;i<atoms;++i)for(j=0;j<max_neighbors;++j)free(BO_record_matrix[i][j]);
    for(i=0;i<atoms;++i)free(BO_record_matrix[i]);
    free(BO_record_matrix);
    
    for(i=0;i<angles;++i)free(BO_array_A1_A2[i]);free(BO_array_A1_A2);
    for(i=0;i<angles;++i)free(BO_array_A2_A3[i]);free(BO_array_A2_A3);
    
    for(i=0;i<angles;++i)free(BO_A1_A2_record_matrix[i]);free(BO_A1_A2_record_matrix);
    for(i=0;i<angles;++i)free(BO_A2_A3_record_matrix[i]);free(BO_A2_A3_record_matrix);
    
    for(i=0;i<dihedrals;++i)free(BO_array_D1_D2[i]);free(BO_array_D1_D2);
    for(i=0;i<dihedrals;++i)free(BO_array_D2_D3[i]);free(BO_array_D2_D3);
    for(i=0;i<dihedrals;++i)free(BO_array_D3_D4[i]);free(BO_array_D3_D4);

    for(i=0;i<dihedrals;++i)free(BO_D1_D2_record_matrix[i]);free(BO_D1_D2_record_matrix);
    for(i=0;i<dihedrals;++i)free(BO_D2_D3_record_matrix[i]);free(BO_D2_D3_record_matrix);
    for(i=0;i<dihedrals;++i)free(BO_D3_D4_record_matrix[i]);free(BO_D3_D4_record_matrix);

    //printf("topo executed successfully!\n");

}
