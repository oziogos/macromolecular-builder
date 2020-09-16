
#include"builder.h"

void one_two_build(int **one_two, int *B1, int *B2, int atoms, int bonds, int max_1_2);

void remove_deleted(int molecule_atoms,int *rescale_array,
                    char **segment_species_array,
                    int atoms_block_size,int bonds_block_size,
                    int *segment_B1_internal_core_array,int *segment_B2_internal_core_array,char **segment_B_type_core_array,
                    
                    int *molecule_atom_types,
                    int *molecule_bond_types,int *molecule_angle_types,int *molecule_dihedral_types,int *molecule_improper_types,
                    int *molecule_bonds,int *molecule_angles,int *molecule_dihedrals,int molecule_impropers,
                    char ***molecule_species_registry,
                    char ****molecule_bond_types_registry,char ****molecule_angle_types_registry,char ****molecule_dihedral_types_registry,char ****molecule_improper_types_registry,
                    int ***molecule_bonds_registry,int ***molecule_angles_registry,int ***molecule_dihedrals_registry,int **molecule_impropers_registry
                    )
{

    int j,k,row,found,erase;
    int **itemp;
    
    int max_1_2=4,**one_two,subset;
    
    char ***BO_record_matrix;
    char w1[sub_length],w2[sub_length],w3[sub_length],w4[sub_length],w5[sub_length],w6[sub_length],w7[sub_length];
    
    // this routine is triggered when we have H atom deletions in the active molecule
    
    // print info
    //printf("\n$ rescale array:\n");
    //for(j=0;j<atoms_block_size;++j)printf("{%d\t%d}\n",j+1,rescale_array[j]);
    
    // species
    // erase molecular species registry and reconstruct from char **segment_species_array; safe array: deleted atoms have been already excluded!
    // reminder: the building process acts on the 'segment_' arrays while the 'molecule_' arrays act as
    //           intermediate memory elements to be used for 'general_' topology updates
    for(j=0;j<*molecule_atom_types;++j)free((*molecule_species_registry)[j]);free(*molecule_species_registry);
    *molecule_atom_types=1;
    *molecule_species_registry=(char**)malloc(1*sizeof(char*));
    (*molecule_species_registry)[0]=(char*)malloc(sub_length*sizeof(char));
    sprintf((*molecule_species_registry)[0],"%s",segment_species_array[0]);
    for(j=1;j<molecule_atoms;++j)
    {
        found=0;
        for(k=0;k<*molecule_atom_types;++k)
        {
            if(strcmp(segment_species_array[j],(*molecule_species_registry)[k])==0)
            {
                found=1;break;
            }
        }
        if(found==0)
        {
            *molecule_atom_types=*molecule_atom_types+1;
            *molecule_species_registry=(char**)realloc(*molecule_species_registry,(*molecule_atom_types)*sizeof(char*));
            (*molecule_species_registry)[*molecule_atom_types-1]=(char*)malloc(sub_length*sizeof(char));
            sprintf((*molecule_species_registry)[*molecule_atom_types-1],"%s",segment_species_array[j]);
        }
    }
    // print
    //printf("\n$ molecular species:\n");
    //for(j=0;j<*molecule_atom_types;++j)printf("[%d]\t%s\n",j+1,(*molecule_species_registry)[j]);
    //--------------------------------------------------------------------------
    // work on bonds
    // mark entries for deletion with -1
    erase=0;                        // count number of bonds to be deleted
    for(j=0;j<*molecule_bonds;++j)  // loop on bonds registry; int *rescale_array holds the new atom id mapping after the deletions - zeroes correspond to deleted H atoms
    {
        if(rescale_array[(*molecule_bonds_registry)[j][1]-1]==0 || rescale_array[(*molecule_bonds_registry)[j][2]-1]==0)
        {
            erase=erase+1;
            (*molecule_bonds_registry)[j][0]=-1;    // mark deletion in the bonds registry
        }
        else
        {
            (*molecule_bonds_registry)[j][1]=rescale_array[(*molecule_bonds_registry)[j][1]-1]; // apply id rescaling
            (*molecule_bonds_registry)[j][2]=rescale_array[(*molecule_bonds_registry)[j][2]-1];
        }
    }
    // quick and dirty registry resize...
    // store valid, rescaled registry in itemp
    itemp=(int**)malloc((*molecule_bonds-erase)*sizeof(int*));
    for(j=0;j<*molecule_bonds-erase;++j)itemp[j]=(int*)malloc(3*sizeof(int));
    k=-1;
    for(j=0;j<*molecule_bonds;++j)
    {
        if((*molecule_bonds_registry)[j][0]!=-1)
        {
            k=k+1;
            itemp[k][0]=(*molecule_bonds_registry)[j][0];
            itemp[k][1]=(*molecule_bonds_registry)[j][1];
            itemp[k][2]=(*molecule_bonds_registry)[j][2];
        }
    }
    // trim current registry from 'below'
    for(j=*molecule_bonds-1;j>*molecule_bonds-1-erase;--j)free((*molecule_bonds_registry)[j]);
    // paste entries from itemp
    for(j=0;j<*molecule_bonds-erase;++j)
    {
        (*molecule_bonds_registry)[j][0]=itemp[j][0];
        (*molecule_bonds_registry)[j][1]=itemp[j][1];
        (*molecule_bonds_registry)[j][2]=itemp[j][2];
    }
    // free itemp
    for(j=0;j<*molecule_bonds-erase;++j)free(itemp[j]);free(itemp);
    
    
    // overwrite number of bonds
    *molecule_bonds=*molecule_bonds-erase;
    
    subset=1;
    for(j=0;j<*molecule_bonds;++j)
    {
        if(!(((*molecule_bonds_registry)[j][1]==segment_B1_internal_core_array[j] && (*molecule_bonds_registry)[j][2]==segment_B2_internal_core_array[j])
           ||
           ((*molecule_bonds_registry)[j][2]==segment_B1_internal_core_array[j] && (*molecule_bonds_registry)[j][1]==segment_B2_internal_core_array[j])))
        {
            subset=0;printf("$ remove_deleted() subset issue!\n");exit(-1);
        }
        //else{printf("%d\t%d\t%d\t%d\t%s\n",(*molecule_bonds_registry)[j][1],(*molecule_bonds_registry)[j][2],segment_B1_internal_core_array[j],segment_B2_internal_core_array[j],segment_B_type_core_array[j]);}
    }
    
    //
    
    one_two=(int**)malloc(molecule_atoms*sizeof(int*));for(j=0;j<molecule_atoms;++j)one_two[j]=(int*)malloc(max_1_2*sizeof(int));
    one_two_build(one_two,segment_B1_internal_core_array,segment_B2_internal_core_array,atoms_block_size,bonds_block_size,max_1_2);

    BO_record_matrix=(char***)malloc(atoms_block_size*sizeof(char**));
    for(j=0;j<atoms_block_size;++j)BO_record_matrix[j]=(char**)malloc(max_1_2*sizeof(char*));
    for(j=0;j<atoms_block_size;++j)for(k=0;k<max_1_2;++k)BO_record_matrix[j][k]=(char*)malloc(sub_length*sizeof(char));
    for(j=0;j<atoms_block_size;++j)for(k=0;k<max_1_2;++k)sprintf(BO_record_matrix[j][k],"%s","__");
    for(j=0;j<*molecule_bonds;++j)
    {
        row=(*molecule_bonds_registry)[j][1]-1;
        for(k=0;k<max_1_2;++k)if(one_two[row][k]==(*molecule_bonds_registry)[j][2])break;
        sprintf(BO_record_matrix[row][k],"%s",segment_B_type_core_array[j]);
        row=(*molecule_bonds_registry)[j][2]-1;
        for(k=0;k<max_1_2;++k)if(one_two[row][k]==(*molecule_bonds_registry)[j][1])break;
        sprintf(BO_record_matrix[row][k],"%s",segment_B_type_core_array[j]);
    }
    /*
    for(j=0;j<atoms_block_size;++j){printf("[%d]\t",j+1);for(k=0;k<max_1_2;++k)printf("%d\t",one_two[j][k]);printf("\n");}
    for(j=0;j<atoms_block_size;++j){printf("[%d]\t",j+1);for(k=0;k<max_1_2;++k)printf("%s\t",BO_record_matrix[j][k]);printf("\n");}
    */
    //
    
    /*
    for(j=0;j<*molecule_bond_types;++j)printf("[%d]\t(%s)--%s--(%s)\n",j+1,(*molecule_bond_types_registry)[j][0],(*molecule_bond_types_registry)[j][1],(*molecule_bond_types_registry)[j][2]);
    for(j=0;j<*molecule_bonds;++j)printf("[%d]\t%d\t%d\t%d\t|\t|\t%s\t%s\n",
                                        j+1,(*molecule_bonds_registry)[j][0],(*molecule_bonds_registry)[j][1],(*molecule_bonds_registry)[j][2],
                                        segment_species_array[(*molecule_bonds_registry)[j][1]-1],segment_species_array[(*molecule_bonds_registry)[j][2]-1]);
    */
    
    //
    
    
    
    // the following code segment is utilized when H atoms are deleted from a molecule
    // after we apply the topological connectivity rescaling to the molecules' bond registry, discarding this way the deleted bonds, we turn our attention to the bond types registry of the molecule
    // what the following segment does:
    // - flush char ***molecule_bond_types_registry
    // - reset molecule_bond_types to 1
    // - take the first bond, resolve its type and write it in the first entry of molecule_bond_types_registry
    // - loop on the rest molecule bonds and conditionally add their type to molecule_bond_types_registry (augmenting molecule_bond_types) via appropriate memory reallocations
    // - once this is done, loop over all bonds and update type mapping in molecule_bonds_registry
    
    // resolve types and match
    for(j=1;j<*molecule_bond_types;++j)
    {
        free((*molecule_bond_types_registry)[j][0]);
        free((*molecule_bond_types_registry)[j][1]);
        free((*molecule_bond_types_registry)[j][2]);
        free((*molecule_bond_types_registry)[j]);
    }
    *molecule_bond_types=1;
    
    j=0;
    sprintf((*molecule_bond_types_registry)[j][0],"%s",segment_species_array[(*molecule_bonds_registry)[j][1]-1]);
    for(k=0;k<max_1_2;++k)if(one_two[(*molecule_bonds_registry)[j][1]-1][k]==(*molecule_bonds_registry)[j][2])break;
    sprintf((*molecule_bond_types_registry)[j][1],"%s",BO_record_matrix[(*molecule_bonds_registry)[j][1]-1][k]);
    sprintf((*molecule_bond_types_registry)[j][2],"%s",segment_species_array[(*molecule_bonds_registry)[j][2]-1]);
    //printf("first type: (%s)--%s--(%s)\n",(*molecule_bond_types_registry)[j][0],(*molecule_bond_types_registry)[j][1],(*molecule_bond_types_registry)[j][2]);
    (*molecule_bonds_registry)[j][0]=1;
    
    for(j=1;j<*molecule_bonds;++j)
    {
        for(k=0;k<max_1_2;++k)if(one_two[(*molecule_bonds_registry)[j][1]-1][k]==(*molecule_bonds_registry)[j][2])break;
        //printf("[%d]\t(%s)--%s--(%s)\n",j+1,
        //       segment_species_array[(*molecule_bonds_registry)[j][1]-1],
         //      BO_record_matrix[(*molecule_bonds_registry)[j][1]-1][k],
          //     segment_species_array[(*molecule_bonds_registry)[j][2]-1]);
        sprintf(w1,"%s",segment_species_array[(*molecule_bonds_registry)[j][1]-1]);
        sprintf(w2,"%s",BO_record_matrix[(*molecule_bonds_registry)[j][1]-1][k]);
        sprintf(w3,"%s",segment_species_array[(*molecule_bonds_registry)[j][2]-1]);
        
        found=0;
        for(k=0;k<*molecule_bond_types;++k)
        {
            if((strcmp(w1,(*molecule_bond_types_registry)[k][0])==0 &&
                strcmp(w2,(*molecule_bond_types_registry)[k][1])==0 &&
                strcmp(w3,(*molecule_bond_types_registry)[k][2])==0)
               ||
               (strcmp(w3,(*molecule_bond_types_registry)[k][0])==0 &&
                strcmp(w2,(*molecule_bond_types_registry)[k][1])==0 &&
                strcmp(w1,(*molecule_bond_types_registry)[k][2])==0))
            {
                found=1;
                (*molecule_bonds_registry)[j][0]=k+1;
                break;
            }
        }
        if(found==0)
        {
            *molecule_bond_types=*molecule_bond_types+1;
            *molecule_bond_types_registry=(char***)realloc(*molecule_bond_types_registry,(*molecule_bond_types)*sizeof(char**));
            (*molecule_bond_types_registry)[*molecule_bond_types-1]=(char**)malloc(3*sizeof(char*));
            (*molecule_bond_types_registry)[*molecule_bond_types-1][0]=(char*)malloc(sub_length*sizeof(char));
            (*molecule_bond_types_registry)[*molecule_bond_types-1][1]=(char*)malloc(sub_length*sizeof(char));
            (*molecule_bond_types_registry)[*molecule_bond_types-1][2]=(char*)malloc(sub_length*sizeof(char));
            sprintf((*molecule_bond_types_registry)[*molecule_bond_types-1][0],"%s",w1);
            sprintf((*molecule_bond_types_registry)[*molecule_bond_types-1][1],"%s",w2);
            sprintf((*molecule_bond_types_registry)[*molecule_bond_types-1][2],"%s",w3);
            (*molecule_bonds_registry)[j][0]=*molecule_bond_types;
        }
        
    }
    /*
    printf("** remove_deleted()\n");
    for(j=0;j<*molecule_bond_types;++j)printf("[%d]\t(%s)--%s--(%s)\n",j+1,(*molecule_bond_types_registry)[j][0],(*molecule_bond_types_registry)[j][1],(*molecule_bond_types_registry)[j][2]);
    for(j=0;j<*molecule_bonds;++j)printf("[%d]\t%d\t%d\t%d\t|\t|\t%s\t%s\n",
                                         j+1,(*molecule_bonds_registry)[j][0],(*molecule_bonds_registry)[j][1],(*molecule_bonds_registry)[j][2],
                                         segment_species_array[(*molecule_bonds_registry)[j][1]-1],segment_species_array[(*molecule_bonds_registry)[j][2]-1]);
    */
    //--------------------------------------------------------------------------
    
    // work on angles
    // mark entries for deletion with -1
    erase=0;
    for(j=0;j<*molecule_angles;++j)
    {
        if(rescale_array[(*molecule_angles_registry)[j][1]-1]==0 || rescale_array[(*molecule_angles_registry)[j][2]-1]==0 || rescale_array[(*molecule_angles_registry)[j][3]-1]==0)
        {
            erase=erase+1;
            (*molecule_angles_registry)[j][0]=-1;
        }
        else
        {
            (*molecule_angles_registry)[j][1]=rescale_array[(*molecule_angles_registry)[j][1]-1];
            (*molecule_angles_registry)[j][2]=rescale_array[(*molecule_angles_registry)[j][2]-1];
            (*molecule_angles_registry)[j][3]=rescale_array[(*molecule_angles_registry)[j][3]-1];
        }
    }
    // quick and dirty registry resize...
    itemp=(int**)malloc((*molecule_angles-erase)*sizeof(int*));
    for(j=0;j<*molecule_angles-erase;++j)itemp[j]=(int*)malloc(4*sizeof(int));
    k=-1;
    for(j=0;j<*molecule_angles;++j)
    {
        if((*molecule_angles_registry)[j][0]!=-1)
        {
            k=k+1;
            itemp[k][0]=(*molecule_angles_registry)[j][0];
            itemp[k][1]=(*molecule_angles_registry)[j][1];
            itemp[k][2]=(*molecule_angles_registry)[j][2];
            itemp[k][3]=(*molecule_angles_registry)[j][3];
        }
    }
    for(j=*molecule_angles-1;j>*molecule_angles-1-erase;--j)free((*molecule_angles_registry)[j]);
    for(j=0;j<*molecule_angles-erase;++j)
    {
        (*molecule_angles_registry)[j][0]=itemp[j][0];
        (*molecule_angles_registry)[j][1]=itemp[j][1];
        (*molecule_angles_registry)[j][2]=itemp[j][2];
        (*molecule_angles_registry)[j][3]=itemp[j][3];
    }
    for(j=0;j<*molecule_angles-erase;++j)free(itemp[j]);free(itemp);
    *molecule_angles=*molecule_angles-erase;
    /*
    for(j=0;j<*molecule_angle_types;++j)printf("[%d]\t(%s)--%s--(%s)--%s--(%s)\n",j+1,(*molecule_angle_types_registry)[j][0],(*molecule_angle_types_registry)[j][1],(*molecule_angle_types_registry)[j][2],(*molecule_angle_types_registry)[j][3],(*molecule_angle_types_registry)[j][4]);
    for(j=0;j<*molecule_angles;++j)printf("[%d]\t%d\t%d\t%d\t%d\t|\t|\t%s\t%s\t%s\n",
                                         j+1,(*molecule_angles_registry)[j][0],(*molecule_angles_registry)[j][1],(*molecule_angles_registry)[j][2],(*molecule_angles_registry)[j][3],
                                         segment_species_array[(*molecule_angles_registry)[j][1]-1],segment_species_array[(*molecule_angles_registry)[j][2]-1],segment_species_array[(*molecule_angles_registry)[j][3]-1]);
    */
    // resolve types and match
    for(j=1;j<*molecule_angle_types;++j)
    {
        free((*molecule_angle_types_registry)[j][0]);
        free((*molecule_angle_types_registry)[j][1]);
        free((*molecule_angle_types_registry)[j][2]);
        free((*molecule_angle_types_registry)[j][3]);
        free((*molecule_angle_types_registry)[j][4]);
        free((*molecule_angle_types_registry)[j]);
    }
    
    *molecule_angle_types=1;

    j=0;
    sprintf((*molecule_angle_types_registry)[j][0],"%s",segment_species_array[(*molecule_angles_registry)[j][1]-1]);
    
    for(k=0;k<max_1_2;++k)if(one_two[(*molecule_angles_registry)[j][1]-1][k]==(*molecule_angles_registry)[j][2])break;
    sprintf((*molecule_angle_types_registry)[j][1],"%s",BO_record_matrix[(*molecule_angles_registry)[j][1]-1][k]);
    
    sprintf((*molecule_angle_types_registry)[j][2],"%s",segment_species_array[(*molecule_angles_registry)[j][2]-1]);
    
    for(k=0;k<max_1_2;++k)if(one_two[(*molecule_angles_registry)[j][2]-1][k]==(*molecule_angles_registry)[j][3])break;
    sprintf((*molecule_angle_types_registry)[j][3],"%s",BO_record_matrix[(*molecule_angles_registry)[j][2]-1][k]);
    
    sprintf((*molecule_angle_types_registry)[j][4],"%s",segment_species_array[(*molecule_angles_registry)[j][3]-1]);

    //printf("first type: (%s)--%s--(%s)--%s--(%s)\n",(*molecule_angle_types_registry)[j][0],(*molecule_angle_types_registry)[j][1],(*molecule_angle_types_registry)[j][2],(*molecule_angle_types_registry)[j][3],(*molecule_angle_types_registry)[j][4]);
    (*molecule_angles_registry)[j][0]=1;

    for(j=1;j<*molecule_angles;++j)
    {
        sprintf(w1,"%s",segment_species_array[(*molecule_angles_registry)[j][1]-1]);

        for(k=0;k<max_1_2;++k)if(one_two[(*molecule_angles_registry)[j][1]-1][k]==(*molecule_angles_registry)[j][2])break;
        sprintf(w2,"%s",BO_record_matrix[(*molecule_angles_registry)[j][1]-1][k]);
        
        sprintf(w3,"%s",segment_species_array[(*molecule_angles_registry)[j][2]-1]);
        
        for(k=0;k<max_1_2;++k)if(one_two[(*molecule_angles_registry)[j][2]-1][k]==(*molecule_angles_registry)[j][3])break;
        sprintf(w4,"%s",BO_record_matrix[(*molecule_angles_registry)[j][2]-1][k]);
        
        sprintf(w5,"%s",segment_species_array[(*molecule_angles_registry)[j][3]-1]);
        
        found=0;
        for(k=0;k<*molecule_angle_types;++k)
        {
            if((strcmp(w1,(*molecule_angle_types_registry)[k][0])==0 &&
                strcmp(w2,(*molecule_angle_types_registry)[k][1])==0 &&
                strcmp(w3,(*molecule_angle_types_registry)[k][2])==0 &&
                strcmp(w4,(*molecule_angle_types_registry)[k][3])==0 &&
                strcmp(w5,(*molecule_angle_types_registry)[k][4])==0
                )
               ||
               (strcmp(w1,(*molecule_angle_types_registry)[k][4])==0 &&
                strcmp(w2,(*molecule_angle_types_registry)[k][3])==0 &&
                strcmp(w3,(*molecule_angle_types_registry)[k][2])==0 &&
                strcmp(w4,(*molecule_angle_types_registry)[k][1])==0 &&
                strcmp(w5,(*molecule_angle_types_registry)[k][0])==0
                ))
            {
                found=1;
                (*molecule_angles_registry)[j][0]=k+1;
                break;
            }
        }
        if(found==0)
        {
            *molecule_angle_types=*molecule_angle_types+1;
            *molecule_angle_types_registry=(char***)realloc(*molecule_angle_types_registry,(*molecule_angle_types)*sizeof(char**));
            (*molecule_angle_types_registry)[*molecule_angle_types-1]=(char**)malloc(5*sizeof(char*));
            (*molecule_angle_types_registry)[*molecule_angle_types-1][0]=(char*)malloc(sub_length*sizeof(char));
            (*molecule_angle_types_registry)[*molecule_angle_types-1][1]=(char*)malloc(sub_length*sizeof(char));
            (*molecule_angle_types_registry)[*molecule_angle_types-1][2]=(char*)malloc(sub_length*sizeof(char));
            (*molecule_angle_types_registry)[*molecule_angle_types-1][3]=(char*)malloc(sub_length*sizeof(char));
            (*molecule_angle_types_registry)[*molecule_angle_types-1][4]=(char*)malloc(sub_length*sizeof(char));
            sprintf((*molecule_angle_types_registry)[*molecule_angle_types-1][0],"%s",w1);
            sprintf((*molecule_angle_types_registry)[*molecule_angle_types-1][1],"%s",w2);
            sprintf((*molecule_angle_types_registry)[*molecule_angle_types-1][2],"%s",w3);
            sprintf((*molecule_angle_types_registry)[*molecule_angle_types-1][3],"%s",w4);
            sprintf((*molecule_angle_types_registry)[*molecule_angle_types-1][4],"%s",w5);
            (*molecule_angles_registry)[j][0]=*molecule_angle_types;
        }
        
    }
    /*
    for(j=0;j<*molecule_angle_types;++j)printf("[%d]\t(%s)--%s--(%s)--%s--(%s)\n",j+1,(*molecule_angle_types_registry)[j][0],(*molecule_angle_types_registry)[j][1],(*molecule_angle_types_registry)[j][2],(*molecule_angle_types_registry)[j][3],(*molecule_angle_types_registry)[j][4]);
    for(j=0;j<*molecule_angles;++j)printf("[%d]\t%d\t%d\t%d\t%d\t|\t|\t%s\t%s\t%s\n",
                                          j+1,(*molecule_angles_registry)[j][0],(*molecule_angles_registry)[j][1],(*molecule_angles_registry)[j][2],(*molecule_angles_registry)[j][3],
                                         segment_species_array[(*molecule_angles_registry)[j][1]-1],segment_species_array[(*molecule_angles_registry)[j][2]-1],segment_species_array[(*molecule_angles_registry)[j][3]-1]);
    */
    //--------------------------------------------------------------------------
    
    //
    if(*molecule_dihedrals>0)
    {
        // work on dihedrals
        // mark entries for deletion with -1
        erase=0;
        for(j=0;j<*molecule_dihedrals;++j)
        {
            if(rescale_array[(*molecule_dihedrals_registry)[j][1]-1]==0 || rescale_array[(*molecule_dihedrals_registry)[j][2]-1]==0 || rescale_array[(*molecule_dihedrals_registry)[j][3]-1]==0 || rescale_array[(*molecule_dihedrals_registry)[j][4]-1]==0)
            {
                erase=erase+1;
                (*molecule_dihedrals_registry)[j][0]=-1;
            }
            else
            {
                (*molecule_dihedrals_registry)[j][1]=rescale_array[(*molecule_dihedrals_registry)[j][1]-1];
                (*molecule_dihedrals_registry)[j][2]=rescale_array[(*molecule_dihedrals_registry)[j][2]-1];
                (*molecule_dihedrals_registry)[j][3]=rescale_array[(*molecule_dihedrals_registry)[j][3]-1];
                (*molecule_dihedrals_registry)[j][4]=rescale_array[(*molecule_dihedrals_registry)[j][4]-1];
            }
        }
        // quick and dirty registry resize...
        itemp=(int**)malloc((*molecule_dihedrals-erase)*sizeof(int*));
        for(j=0;j<*molecule_dihedrals-erase;++j)itemp[j]=(int*)malloc(5*sizeof(int));
        k=-1;
        for(j=0;j<*molecule_dihedrals;++j)
        {
            if((*molecule_dihedrals_registry)[j][0]!=-1)
            {
                k=k+1;
                itemp[k][0]=(*molecule_dihedrals_registry)[j][0];
                itemp[k][1]=(*molecule_dihedrals_registry)[j][1];
                itemp[k][2]=(*molecule_dihedrals_registry)[j][2];
                itemp[k][3]=(*molecule_dihedrals_registry)[j][3];
                itemp[k][4]=(*molecule_dihedrals_registry)[j][4];
            }
        }
        for(j=*molecule_dihedrals-1;j>*molecule_dihedrals-1-erase;--j)free((*molecule_dihedrals_registry)[j]);
        for(j=0;j<*molecule_dihedrals-erase;++j)
        {
            (*molecule_dihedrals_registry)[j][0]=itemp[j][0];
            (*molecule_dihedrals_registry)[j][1]=itemp[j][1];
            (*molecule_dihedrals_registry)[j][2]=itemp[j][2];
            (*molecule_dihedrals_registry)[j][3]=itemp[j][3];
            (*molecule_dihedrals_registry)[j][4]=itemp[j][4];
        }
        for(j=0;j<*molecule_dihedrals-erase;++j)free(itemp[j]);free(itemp);
        *molecule_dihedrals=*molecule_dihedrals-erase;
        /*
        for(j=0;j<*molecule_dihedral_types;++j)printf("[%d]\t(%s)--%s--(%s)--%s--(%s)--%s--(%s)\n",j+1,(*molecule_dihedral_types_registry)[j][0],(*molecule_dihedral_types_registry)[j][1],(*molecule_dihedral_types_registry)[j][2],(*molecule_dihedral_types_registry)[j][3],(*molecule_dihedral_types_registry)[j][4],(*molecule_dihedral_types_registry)[j][5],(*molecule_dihedral_types_registry)[j][6]);
        for(j=0;j<*molecule_dihedrals;++j)printf("[%d]\t%d\t%d\t%d\t%d\t%d\t|\t|\t%s\t%s\t%s\t%s\n",
                                                 
                                                 j+1,
                                                 (*molecule_dihedrals_registry)[j][0],
                                                 (*molecule_dihedrals_registry)[j][1],
                                                 (*molecule_dihedrals_registry)[j][2],
                                                 (*molecule_dihedrals_registry)[j][3],
                                                 (*molecule_dihedrals_registry)[j][4],
                                                 
                                                 segment_species_array[(*molecule_dihedrals_registry)[j][1]-1],
                                                 segment_species_array[(*molecule_dihedrals_registry)[j][2]-1],
                                                 segment_species_array[(*molecule_dihedrals_registry)[j][3]-1],
                                                 segment_species_array[(*molecule_dihedrals_registry)[j][4]-1]);
        */        
        // resolve types and match
        for(j=1;j<*molecule_dihedral_types;++j)
        {
            free((*molecule_dihedral_types_registry)[j][0]);
            free((*molecule_dihedral_types_registry)[j][1]);
            free((*molecule_dihedral_types_registry)[j][2]);
            free((*molecule_dihedral_types_registry)[j][3]);
            free((*molecule_dihedral_types_registry)[j][4]);
            free((*molecule_dihedral_types_registry)[j][5]);
            free((*molecule_dihedral_types_registry)[j][6]);
            free((*molecule_dihedral_types_registry)[j]);
        }
        
        *molecule_dihedral_types=1;
        
        j=0;
        sprintf((*molecule_dihedral_types_registry)[j][0],"%s",segment_species_array[(*molecule_dihedrals_registry)[j][1]-1]);
        
        for(k=0;k<max_1_2;++k)if(one_two[(*molecule_dihedrals_registry)[j][1]-1][k]==(*molecule_dihedrals_registry)[j][2])break;
        sprintf((*molecule_dihedral_types_registry)[j][1],"%s",BO_record_matrix[(*molecule_dihedrals_registry)[j][1]-1][k]);
        
        sprintf((*molecule_dihedral_types_registry)[j][2],"%s",segment_species_array[(*molecule_dihedrals_registry)[j][2]-1]);
        
        for(k=0;k<max_1_2;++k)if(one_two[(*molecule_dihedrals_registry)[j][2]-1][k]==(*molecule_dihedrals_registry)[j][3])break;
        sprintf((*molecule_dihedral_types_registry)[j][3],"%s",BO_record_matrix[(*molecule_dihedrals_registry)[j][2]-1][k]);
        
        sprintf((*molecule_dihedral_types_registry)[j][4],"%s",segment_species_array[(*molecule_dihedrals_registry)[j][3]-1]);
        
        for(k=0;k<max_1_2;++k)if(one_two[(*molecule_dihedrals_registry)[j][3]-1][k]==(*molecule_dihedrals_registry)[j][4])break;
        sprintf((*molecule_dihedral_types_registry)[j][5],"%s",BO_record_matrix[(*molecule_dihedrals_registry)[j][3]-1][k]);
        
        sprintf((*molecule_dihedral_types_registry)[j][6],"%s",segment_species_array[(*molecule_dihedrals_registry)[j][4]-1]);
        
        //printf("first type: (%s)--%s--(%s)--%s--(%s)--%s--(%s)\n",(*molecule_dihedral_types_registry)[j][0],(*molecule_dihedral_types_registry)[j][1],(*molecule_dihedral_types_registry)[j][2],(*molecule_dihedral_types_registry)[j][3],(*molecule_dihedral_types_registry)[j][4],(*molecule_dihedral_types_registry)[j][5],(*molecule_dihedral_types_registry)[j][6]);
        (*molecule_dihedrals_registry)[j][0]=1;
        
        for(j=1;j<*molecule_dihedrals;++j)
        {
            sprintf(w1,"%s",segment_species_array[(*molecule_dihedrals_registry)[j][1]-1]);
            
            for(k=0;k<max_1_2;++k)if(one_two[(*molecule_dihedrals_registry)[j][1]-1][k]==(*molecule_dihedrals_registry)[j][2])break;
            sprintf(w2,"%s",BO_record_matrix[(*molecule_dihedrals_registry)[j][1]-1][k]);
            
            sprintf(w3,"%s",segment_species_array[(*molecule_dihedrals_registry)[j][2]-1]);
            
            for(k=0;k<max_1_2;++k)if(one_two[(*molecule_dihedrals_registry)[j][2]-1][k]==(*molecule_dihedrals_registry)[j][3])break;
            sprintf(w4,"%s",BO_record_matrix[(*molecule_dihedrals_registry)[j][2]-1][k]);
            
            sprintf(w5,"%s",segment_species_array[(*molecule_dihedrals_registry)[j][3]-1]);
            
            for(k=0;k<max_1_2;++k)if(one_two[(*molecule_dihedrals_registry)[j][3]-1][k]==(*molecule_dihedrals_registry)[j][4])break;
            sprintf(w6,"%s",BO_record_matrix[(*molecule_dihedrals_registry)[j][3]-1][k]);
            
            sprintf(w7,"%s",segment_species_array[(*molecule_dihedrals_registry)[j][4]-1]);
            
            found=0;
            for(k=0;k<*molecule_dihedral_types;++k)
            {
                if((strcmp(w1,(*molecule_dihedral_types_registry)[k][0])==0 &&
                    strcmp(w2,(*molecule_dihedral_types_registry)[k][1])==0 &&
                    strcmp(w3,(*molecule_dihedral_types_registry)[k][2])==0 &&
                    strcmp(w4,(*molecule_dihedral_types_registry)[k][3])==0 &&
                    strcmp(w5,(*molecule_dihedral_types_registry)[k][4])==0 &&
                    strcmp(w6,(*molecule_dihedral_types_registry)[k][5])==0 &&
                    strcmp(w7,(*molecule_dihedral_types_registry)[k][6])==0
                    )
                   ||
                   (strcmp(w1,(*molecule_dihedral_types_registry)[k][6])==0 &&
                    strcmp(w2,(*molecule_dihedral_types_registry)[k][5])==0 &&
                    strcmp(w3,(*molecule_dihedral_types_registry)[k][4])==0 &&
                    strcmp(w4,(*molecule_dihedral_types_registry)[k][3])==0 &&
                    strcmp(w5,(*molecule_dihedral_types_registry)[k][2])==0 &&
                    strcmp(w6,(*molecule_dihedral_types_registry)[k][1])==0 &&
                    strcmp(w7,(*molecule_dihedral_types_registry)[k][0])==0
                    ))
                {
                    found=1;
                    (*molecule_dihedrals_registry)[j][0]=k+1;
                    break;
                }
            }
            if(found==0)
            {
                *molecule_dihedral_types=*molecule_dihedral_types+1;
                *molecule_dihedral_types_registry=(char***)realloc(*molecule_dihedral_types_registry,(*molecule_dihedral_types)*sizeof(char**));
                (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1]=(char**)malloc(7*sizeof(char*));
                (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][0]=(char*)malloc(sub_length*sizeof(char));
                (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][1]=(char*)malloc(sub_length*sizeof(char));
                (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][2]=(char*)malloc(sub_length*sizeof(char));
                (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][3]=(char*)malloc(sub_length*sizeof(char));
                (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][4]=(char*)malloc(sub_length*sizeof(char));
                (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][5]=(char*)malloc(sub_length*sizeof(char));
                (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][6]=(char*)malloc(sub_length*sizeof(char));
                sprintf((*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][0],"%s",w1);
                sprintf((*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][1],"%s",w2);
                sprintf((*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][2],"%s",w3);
                sprintf((*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][3],"%s",w4);
                sprintf((*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][4],"%s",w5);
                sprintf((*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][5],"%s",w6);
                sprintf((*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][6],"%s",w7);
                (*molecule_dihedrals_registry)[j][0]=*molecule_dihedral_types;
            }
            
        }
        /*
        for(j=0;j<*molecule_dihedral_types;++j)printf("[%d]\t(%s)--%s--(%s)--%s--(%s)--%s--(%s)\n",j+1,(*molecule_dihedral_types_registry)[j][0],(*molecule_dihedral_types_registry)[j][1],(*molecule_dihedral_types_registry)[j][2],(*molecule_dihedral_types_registry)[j][3],(*molecule_dihedral_types_registry)[j][4],(*molecule_dihedral_types_registry)[j][5],(*molecule_dihedral_types_registry)[j][6]);
            for(j=0;j<*molecule_dihedrals;++j)printf("[%d]\t%d\t%d\t%d\t%d\t%d\t|\t|\t%s\t%s\t%s\t%s\n",
                                                  j+1,(*molecule_dihedrals_registry)[j][0],(*molecule_dihedrals_registry)[j][1],(*molecule_dihedrals_registry)[j][2],(*molecule_dihedrals_registry)[j][3],(*molecule_dihedrals_registry)[j][4],
                                                  segment_species_array[(*molecule_dihedrals_registry)[j][1]-1],segment_species_array[(*molecule_dihedrals_registry)[j][2]-1],segment_species_array[(*molecule_dihedrals_registry)[j][3]-1],segment_species_array[(*molecule_dihedrals_registry)[j][4]-1]);
            
        */
            
            
            
            
            
    }
    
    // impropers
    if(molecule_impropers>0)
    {
        erase=0;
        for(j=0;j<molecule_impropers;++j)
        {
            if(rescale_array[molecule_impropers_registry[j][1]-1]==0 || rescale_array[molecule_impropers_registry[j][2]-1]==0 || rescale_array[molecule_impropers_registry[j][3]-1]==0 || rescale_array[molecule_impropers_registry[j][4]-1]==0)
            {
                erase=erase+1;
            }
            else
            {
                molecule_impropers_registry[j][1]=rescale_array[molecule_impropers_registry[j][1]-1];
                molecule_impropers_registry[j][2]=rescale_array[molecule_impropers_registry[j][2]-1];
                molecule_impropers_registry[j][3]=rescale_array[molecule_impropers_registry[j][3]-1];
                molecule_impropers_registry[j][4]=rescale_array[molecule_impropers_registry[j][4]-1];
            }
        }
        if(erase>0)
        {
            printf("$ rescaling encountered improper hydrogen!!!\n\n");exit(-10);
        }
        
        // resolve types and match
        for(j=1;j<*molecule_improper_types;++j)
        {
            free((*molecule_improper_types_registry)[j][0]);
            free((*molecule_improper_types_registry)[j][1]);
            free((*molecule_improper_types_registry)[j][2]);
            free((*molecule_improper_types_registry)[j][3]);
            free((*molecule_improper_types_registry)[j]);
        }
        *molecule_improper_types=1;
        sprintf((*molecule_improper_types_registry)[0][0],"%s",segment_species_array[molecule_impropers_registry[0][1]-1]);
        sprintf((*molecule_improper_types_registry)[0][1],"%s",segment_species_array[molecule_impropers_registry[0][2]-1]);
        sprintf((*molecule_improper_types_registry)[0][2],"%s",segment_species_array[molecule_impropers_registry[0][3]-1]);
        sprintf((*molecule_improper_types_registry)[0][3],"%s",segment_species_array[molecule_impropers_registry[0][4]-1]);
        for(j=1;j<molecule_impropers;++j)
        {
            found=0;
            for(k=0;k<*molecule_improper_types;++k)
            {
                if((strcmp(segment_species_array[molecule_impropers_registry[j][1]-1],(*molecule_improper_types_registry)[k][0])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][2]-1],(*molecule_improper_types_registry)[k][1])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][3]-1],(*molecule_improper_types_registry)[k][2])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][4]-1],(*molecule_improper_types_registry)[k][3])==0)
                   ||
                   (strcmp(segment_species_array[molecule_impropers_registry[j][1]-1],(*molecule_improper_types_registry)[k][0])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][2]-1],(*molecule_improper_types_registry)[k][1])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][3]-1],(*molecule_improper_types_registry)[k][3])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][4]-1],(*molecule_improper_types_registry)[k][2])==0)
                   ||
                   (strcmp(segment_species_array[molecule_impropers_registry[j][1]-1],(*molecule_improper_types_registry)[k][0])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][2]-1],(*molecule_improper_types_registry)[k][2])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][3]-1],(*molecule_improper_types_registry)[k][1])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][4]-1],(*molecule_improper_types_registry)[k][3])==0)
                   ||
                   (strcmp(segment_species_array[molecule_impropers_registry[j][1]-1],(*molecule_improper_types_registry)[k][0])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][2]-1],(*molecule_improper_types_registry)[k][2])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][3]-1],(*molecule_improper_types_registry)[k][3])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][4]-1],(*molecule_improper_types_registry)[k][1])==0)
                   ||
                   (strcmp(segment_species_array[molecule_impropers_registry[j][1]-1],(*molecule_improper_types_registry)[k][0])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][2]-1],(*molecule_improper_types_registry)[k][3])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][3]-1],(*molecule_improper_types_registry)[k][1])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][4]-1],(*molecule_improper_types_registry)[k][2])==0)
                   ||
                   (strcmp(segment_species_array[molecule_impropers_registry[j][1]-1],(*molecule_improper_types_registry)[k][0])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][2]-1],(*molecule_improper_types_registry)[k][3])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][3]-1],(*molecule_improper_types_registry)[k][2])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][4]-1],(*molecule_improper_types_registry)[k][1])==0))
                {
                    found=1;break;
                }
            }
            if(found==0)
            {
                *molecule_improper_types=*molecule_improper_types+1;
                *molecule_improper_types_registry=(char***)realloc(*molecule_improper_types_registry,(*molecule_improper_types)*sizeof(char**));
                (*molecule_improper_types_registry)[*molecule_improper_types-1]=(char**)malloc(4*sizeof(char*));
                (*molecule_improper_types_registry)[*molecule_improper_types-1][0]=(char*)malloc(sub_length*sizeof(char));
                (*molecule_improper_types_registry)[*molecule_improper_types-1][1]=(char*)malloc(sub_length*sizeof(char));
                (*molecule_improper_types_registry)[*molecule_improper_types-1][2]=(char*)malloc(sub_length*sizeof(char));
                (*molecule_improper_types_registry)[*molecule_improper_types-1][3]=(char*)malloc(sub_length*sizeof(char));
                sprintf((*molecule_improper_types_registry)[*molecule_improper_types-1][0],"%s",segment_species_array[molecule_impropers_registry[j][1]-1]);
                sprintf((*molecule_improper_types_registry)[*molecule_improper_types-1][1],"%s",segment_species_array[molecule_impropers_registry[j][2]-1]);
                sprintf((*molecule_improper_types_registry)[*molecule_improper_types-1][2],"%s",segment_species_array[molecule_impropers_registry[j][3]-1]);
                sprintf((*molecule_improper_types_registry)[*molecule_improper_types-1][3],"%s",segment_species_array[molecule_impropers_registry[j][4]-1]);
            }
        }
        for(j=0;j<molecule_impropers;++j)
        {
            for(k=0;k<*molecule_improper_types;++k)
            {
                if((strcmp(segment_species_array[molecule_impropers_registry[j][1]-1],(*molecule_improper_types_registry)[k][0])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][2]-1],(*molecule_improper_types_registry)[k][1])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][3]-1],(*molecule_improper_types_registry)[k][2])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][4]-1],(*molecule_improper_types_registry)[k][3])==0)
                   ||
                   (strcmp(segment_species_array[molecule_impropers_registry[j][1]-1],(*molecule_improper_types_registry)[k][0])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][2]-1],(*molecule_improper_types_registry)[k][1])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][3]-1],(*molecule_improper_types_registry)[k][3])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][4]-1],(*molecule_improper_types_registry)[k][2])==0)
                   ||
                   (strcmp(segment_species_array[molecule_impropers_registry[j][1]-1],(*molecule_improper_types_registry)[k][0])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][2]-1],(*molecule_improper_types_registry)[k][2])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][3]-1],(*molecule_improper_types_registry)[k][1])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][4]-1],(*molecule_improper_types_registry)[k][3])==0)
                   ||
                   (strcmp(segment_species_array[molecule_impropers_registry[j][1]-1],(*molecule_improper_types_registry)[k][0])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][2]-1],(*molecule_improper_types_registry)[k][2])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][3]-1],(*molecule_improper_types_registry)[k][3])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][4]-1],(*molecule_improper_types_registry)[k][1])==0)
                   ||
                   (strcmp(segment_species_array[molecule_impropers_registry[j][1]-1],(*molecule_improper_types_registry)[k][0])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][2]-1],(*molecule_improper_types_registry)[k][3])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][3]-1],(*molecule_improper_types_registry)[k][1])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][4]-1],(*molecule_improper_types_registry)[k][2])==0)
                   ||
                   (strcmp(segment_species_array[molecule_impropers_registry[j][1]-1],(*molecule_improper_types_registry)[k][0])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][2]-1],(*molecule_improper_types_registry)[k][3])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][3]-1],(*molecule_improper_types_registry)[k][2])==0 &&
                    strcmp(segment_species_array[molecule_impropers_registry[j][4]-1],(*molecule_improper_types_registry)[k][1])==0))
                {
                    molecule_impropers_registry[j][0]=k+1;break;
                }
            }
        }
        
        
        // print
        /*
        printf("---------------------------\n");
        printf("\n$ molecular improper types:\n");
        for(j=0;j<*molecule_improper_types;++j)printf("[%d]\t%s\t%s\t%s\t%s\n",j+1,(*molecule_improper_types_registry)[j][0],(*molecule_improper_types_registry)[j][1],(*molecule_improper_types_registry)[j][2],(*molecule_improper_types_registry)[j][3]);
        printf("\n$ molecular impropers:\n");
        */
    }

    //
    
    for(j=0;j<molecule_atoms;++j)free(one_two[j]);free(one_two);

    for(j=0;j<atoms_block_size;++j)for(k=0;k<max_1_2;++k)free(BO_record_matrix[j][k]);
    for(j=0;j<atoms_block_size;++j)free(BO_record_matrix[j]);
    free(BO_record_matrix);

    //getchar();

}
