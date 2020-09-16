
#include"builder.h"

void alter_molecular_topo(int local_bonds,int local_angles,int local_dihedrals,int local_impropers,
                          int local_atom_types,
                          
                          int *local_map,
                          
                          char **local_species_registry,
                          int **local_bonds_registry,int **local_angles_registry,int **local_dihedrals_registry,int **local_impropers_registry,
                          char ***local_bond_types_registry,char ***local_angle_types_registry,char ***local_dihedral_types_registry,char ***local_improper_types_registry,
                          
                          // output
                          
                          int *molecule_bonds,int *molecule_angles,int *molecule_dihedrals,int *molecule_impropers,
                          int *molecule_atom_types,int *molecule_bond_types,int *molecule_angle_types,int *molecule_dihedral_types,int *molecule_improper_types,
                          
                          char ***molecule_species_registry,
                          int ***molecule_bonds_registry,int ***molecule_angles_registry,int ***molecule_dihedrals_registry,int ***molecule_impropers_registry,
                          char ****molecule_bond_types_registry,char ****molecule_angle_types_registry,char ****molecule_dihedral_types_registry,char ****molecule_improper_types_registry,
                          
                          int verb
                          )
{
    
    int j,k,n;
    int present,found;
        
    // ATOM TYPES
    
    // j-loop on topology.local.types.P
    for(j=0;j<local_atom_types;++j)
    {
        // present flag
        present=0;
        // k-loop on topology.molecule.types.P
        for(k=0;k<*molecule_atom_types;++k)
        {
            // check if the type exists; yes --> present=1 and break!
            if(strcmp(local_species_registry[j],(*molecule_species_registry)[k])==0)
            {
                present=1;
                break;
            }
        }
        // present == 0 means that the current local type is NOT PRESENT in the molecular registry!
        if(present==0)
        {
            // console out
            //printf("type %s is not present...\n",local_species_registry[j]);
            // augment
            *molecule_atom_types=*molecule_atom_types+1;
            // resize array
            *molecule_species_registry=(char**)realloc(*molecule_species_registry,*molecule_atom_types*sizeof(char*));
            (*molecule_species_registry)[*molecule_atom_types-1]=(char*)malloc(sub_length*sizeof(char));
            // store new type
            sprintf((*molecule_species_registry)[*molecule_atom_types-1],"%s",local_species_registry[j]);
        }
    }
    /*
    // console out
    for(j=0;j<*molecule_atom_types;++j)
    {
        printf("[%d]\t%s\n",j+1,(*molecule_species_registry)[j]);
    }
    */
    // BONDS!!
    
    
    // local bonds loop; j-loop on topology.local.arrays.B
    for(j=0;j<local_bonds;++j)
    {
        // present flag --> 0
        present=0;
        // molecule bonds loop; k-loop on topology.molecule.arrays.B
        for(k=0;k<*molecule_bonds;++k)
        {
            // search for matching IDs: found?
            if(((*molecule_bonds_registry)[k][1]==local_map[local_bonds_registry[j][1]-1] && (*molecule_bonds_registry)[k][2]==local_map[local_bonds_registry[j][2]-1])
               ||
               ((*molecule_bonds_registry)[k][1]==local_map[local_bonds_registry[j][2]-1] && (*molecule_bonds_registry)[k][2]==local_map[local_bonds_registry[j][1]-1])
               )
                // (search for matching IDs: found?) YES!
            {
                /*
                // console out
                printf("bond %d-%d is present!\n",local_map[local_bonds_registry[j][1]-1],local_map[local_bonds_registry[j][2]-1]);
                printf("... checking bond type: molecule --> (%s)--%s--(%s) | local --> (%s)--%s--(%s)\n",
                       (*molecule_bond_types_registry)[(*molecule_bonds_registry)[k][0]-1][0],
                       (*molecule_bond_types_registry)[(*molecule_bonds_registry)[k][0]-1][1],
                       (*molecule_bond_types_registry)[(*molecule_bonds_registry)[k][0]-1][2],
                       local_bond_types_registry[local_bonds_registry[j][0]-1][0],
                       local_bond_types_registry[local_bonds_registry[j][0]-1][1],
                       local_bond_types_registry[local_bonds_registry[j][0]-1][2]);
                */
                // check the type: is it the same?
                if((strcmp((*molecule_bond_types_registry)[(*molecule_bonds_registry)[k][0]-1][0],local_bond_types_registry[local_bonds_registry[j][0]-1][0])==0 &&
                    strcmp((*molecule_bond_types_registry)[(*molecule_bonds_registry)[k][0]-1][1],local_bond_types_registry[local_bonds_registry[j][0]-1][1])==0 &&
                    strcmp((*molecule_bond_types_registry)[(*molecule_bonds_registry)[k][0]-1][2],local_bond_types_registry[local_bonds_registry[j][0]-1][2])==0
                    )
                   ||
                   (strcmp((*molecule_bond_types_registry)[(*molecule_bonds_registry)[k][0]-1][0],local_bond_types_registry[local_bonds_registry[j][0]-1][2])==0 &&
                    strcmp((*molecule_bond_types_registry)[(*molecule_bonds_registry)[k][0]-1][1],local_bond_types_registry[local_bonds_registry[j][0]-1][1])==0 &&
                    strcmp((*molecule_bond_types_registry)[(*molecule_bonds_registry)[k][0]-1][2],local_bond_types_registry[local_bonds_registry[j][0]-1][0])==0
                    )
                   
                   )
                    // (check the type: is it the same?) YES!
                {
                    // console out
                    //printf("... (%d) same type - the bond is unaltered!!\n",(*molecule_bonds_registry)[k][0]);
                    // set present flag --> 1 and break!
                    present=1;
                    break; // break the k-loop
                }
                else
                    // (check the type: is it the same?) NO!!
                {
                    // search molecule bond types: found?
                    for(n=0;n<*molecule_bond_types;++n)
                    {
                        //printf("... comparing with [%d] %s-%s ...\n",n+1,(*molecule_bond_types_registry)[n][0],(*molecule_bond_types_registry)[n][1]);
                        if((strcmp((*molecule_bond_types_registry)[n][0],local_bond_types_registry[local_bonds_registry[j][0]-1][0])==0 &&
                            strcmp((*molecule_bond_types_registry)[n][1],local_bond_types_registry[local_bonds_registry[j][0]-1][1])==0 &&
                            strcmp((*molecule_bond_types_registry)[n][2],local_bond_types_registry[local_bonds_registry[j][0]-1][2])==0
                            )
                           ||
                           (strcmp((*molecule_bond_types_registry)[n][0],local_bond_types_registry[local_bonds_registry[j][0]-1][2])==0 &&
                            strcmp((*molecule_bond_types_registry)[n][1],local_bond_types_registry[local_bonds_registry[j][0]-1][1])==0 &&
                            strcmp((*molecule_bond_types_registry)[n][2],local_bond_types_registry[local_bonds_registry[j][0]-1][0])==0
                            ))
                            // (search molecule bond types: found?) YES!
                        {
                            // update type in the molecule bonds registry
                            (*molecule_bonds_registry)[k][0]=n+1;
                            // set present flag --> 1 and break!
                            present=1;
                            break; // break the n-loop
                        }
                    }
                    // (search molecule bond types: found?) NO!!
                    if(present==0)
                    {
                        /*
                        // add new type to molecule bond types
                        printf("+++ I will need to add type (%s)--%s--(%s) to the registry...\n",
                               local_bond_types_registry[local_bonds_registry[j][0]-1][0],
                               local_bond_types_registry[local_bonds_registry[j][0]-1][1],
                               local_bond_types_registry[local_bonds_registry[j][0]-1][2]);
                        */
                        // augment molecule bond types
                        *molecule_bond_types=*molecule_bond_types+1;
                        // realloc
                        *molecule_bond_types_registry=(char***)realloc(*molecule_bond_types_registry,*molecule_bond_types*sizeof(char**));
                        (*molecule_bond_types_registry)[*molecule_bond_types-1]=(char**)malloc(3*sizeof(char*));
                        (*molecule_bond_types_registry)[*molecule_bond_types-1][0]=(char*)malloc(sub_length*sizeof(char));
                        (*molecule_bond_types_registry)[*molecule_bond_types-1][1]=(char*)malloc(sub_length*sizeof(char));
                        (*molecule_bond_types_registry)[*molecule_bond_types-1][2]=(char*)malloc(sub_length*sizeof(char));
                        // write type
                        sprintf((*molecule_bond_types_registry)[*molecule_bond_types-1][0],"%s",local_bond_types_registry[local_bonds_registry[j][0]-1][0]);
                        sprintf((*molecule_bond_types_registry)[*molecule_bond_types-1][1],"%s",local_bond_types_registry[local_bonds_registry[j][0]-1][1]);
                        sprintf((*molecule_bond_types_registry)[*molecule_bond_types-1][2],"%s",local_bond_types_registry[local_bonds_registry[j][0]-1][2]);
                        // update type in the molecule bonds registry
                        (*molecule_bonds_registry)[k][0]=*molecule_bond_types;
                        // check
                        //for(n=0;n<*molecule_bond_types;++n)printf("$$$ %s-%s\n",(*molecule_bond_types_registry)[n][0],(*molecule_bond_types_registry)[n][1]);
                        // set present flag --> 1 and break!
                        present=1;
                        break; // break the k-loop
                    }
                    else
                    {
                        // present is already equal to 1!!
                        break; // break the k-loop
                    }
                }
            }
        }
        if(present==0)
        {
            /*
            printf("*** %d-%d\n",local_bonds_registry[j][1]-1,local_bonds_registry[j][2]-1);
            
            printf("bond %d-%d is NOT present!\n",local_map[local_bonds_registry[j][1]-1],local_map[local_bonds_registry[j][2]-1]);
            
            printf("+++ its type is (%s)--%s--(%s)\n",local_bond_types_registry[local_bonds_registry[j][0]-1][0],local_bond_types_registry[local_bonds_registry[j][0]-1][1],local_bond_types_registry[local_bonds_registry[j][0]-1][2]);
            */
            // add new bond to the molecule bonds registry
            // augment
            *molecule_bonds=*molecule_bonds+1;
            // realloc
            *molecule_bonds_registry=(int**)realloc(*molecule_bonds_registry,*molecule_bonds*sizeof(int*));
            (*molecule_bonds_registry)[*molecule_bonds-1]=(int*)malloc(3*sizeof(int));
            // populate: ids
            (*molecule_bonds_registry)[*molecule_bonds-1][1]=local_map[local_bonds_registry[j][1]-1];
            (*molecule_bonds_registry)[*molecule_bonds-1][2]=local_map[local_bonds_registry[j][2]-1];
            // search molecule bond types: found?
            // initialize flag
            found=0;
            for(n=0;n<*molecule_bond_types;++n)
            {
                //printf("... comparing with [%d] %s-%s ...\n",n+1,(*molecule_bond_types_registry)[n][0],(*molecule_bond_types_registry)[n][1]);
                if((strcmp((*molecule_bond_types_registry)[n][0],local_bond_types_registry[local_bonds_registry[j][0]-1][0])==0 &&
                    strcmp((*molecule_bond_types_registry)[n][1],local_bond_types_registry[local_bonds_registry[j][0]-1][1])==0 &&
                    strcmp((*molecule_bond_types_registry)[n][2],local_bond_types_registry[local_bonds_registry[j][0]-1][2])==0)
                   ||
                   (strcmp((*molecule_bond_types_registry)[n][0],local_bond_types_registry[local_bonds_registry[j][0]-1][2])==0 &&
                    strcmp((*molecule_bond_types_registry)[n][1],local_bond_types_registry[local_bonds_registry[j][0]-1][1])==0 &&
                    strcmp((*molecule_bond_types_registry)[n][2],local_bond_types_registry[local_bonds_registry[j][0]-1][0])==0))
                    // (search molecule bond types: found?) YES!
                {
                    // update type in the molecule bonds registry
                    (*molecule_bonds_registry)[*molecule_bonds-1][0]=n+1;
                    // set found flag --> 1 and break!
                    found=1;
                    break; // break the n-loop
                }
            }
            // (search molecule bond types: found?) NO!!
            if(found==0)
            {
                /*
                // add new type to molecule bond types
                printf("+++ I will need to add type (%s)--%s--(%s) to the registry...\n",
                       local_bond_types_registry[local_bonds_registry[j][0]-1][0],
                       local_bond_types_registry[local_bonds_registry[j][0]-1][1],
                       local_bond_types_registry[local_bonds_registry[j][0]-1][2]);
                */
                // augment molecule bond types
                *molecule_bond_types=*molecule_bond_types+1;
                // realloc
                *molecule_bond_types_registry=(char***)realloc(*molecule_bond_types_registry,*molecule_bond_types*sizeof(char**));
                (*molecule_bond_types_registry)[*molecule_bond_types-1]=(char**)malloc(3*sizeof(char*));
                (*molecule_bond_types_registry)[*molecule_bond_types-1][0]=(char*)malloc(sub_length*sizeof(char));
                (*molecule_bond_types_registry)[*molecule_bond_types-1][1]=(char*)malloc(sub_length*sizeof(char));
                (*molecule_bond_types_registry)[*molecule_bond_types-1][2]=(char*)malloc(sub_length*sizeof(char));
                // write type
                sprintf((*molecule_bond_types_registry)[*molecule_bond_types-1][0],"%s",local_bond_types_registry[local_bonds_registry[j][0]-1][0]);
                sprintf((*molecule_bond_types_registry)[*molecule_bond_types-1][1],"%s",local_bond_types_registry[local_bonds_registry[j][0]-1][1]);
                sprintf((*molecule_bond_types_registry)[*molecule_bond_types-1][2],"%s",local_bond_types_registry[local_bonds_registry[j][0]-1][2]);
                // update type in the molecule bonds registry
                (*molecule_bonds_registry)[*molecule_bonds-1][0]=*molecule_bond_types;
                // check
                //for(n=0;n<*molecule_bond_types;++n)printf("$$$ %s-%s\n",(*molecule_bond_types_registry)[n][0],(*molecule_bond_types_registry)[n][1]);
                
            }
        }
    }
    
    
    
    // ANGLES!!
    
    // local angles loop
    for(j=0;j<local_angles;++j)
    {
        // present flag --> 0
        present=0;
        // molecule angles loop
        for(k=0;k<*molecule_angles;++k)
        {
            // search for matching IDs: found?
            if(((*molecule_angles_registry)[k][1]==local_map[local_angles_registry[j][1]-1] &&
                (*molecule_angles_registry)[k][2]==local_map[local_angles_registry[j][2]-1] &&
                (*molecule_angles_registry)[k][3]==local_map[local_angles_registry[j][3]-1]
                ) ||
               ((*molecule_angles_registry)[k][1]==local_map[local_angles_registry[j][3]-1] &&
                (*molecule_angles_registry)[k][2]==local_map[local_angles_registry[j][2]-1] &&
                (*molecule_angles_registry)[k][3]==local_map[local_angles_registry[j][1]-1]
                )
               )
                // (search for matching IDs: found?) YES!
            {
                /*
                printf("angle %d-%d-%d is present!\n",
                       local_map[local_angles_registry[j][1]-1],
                       local_map[local_angles_registry[j][2]-1],
                       local_map[local_angles_registry[j][3]-1]);
                
                printf("... checking angle type: molecule --> %s-%s-%s | local --> %s-%s-%s\n",
                       (*molecule_angle_types_registry)[(*molecule_angles_registry)[k][0]-1][0],
                       (*molecule_angle_types_registry)[(*molecule_angles_registry)[k][0]-1][1],
                       (*molecule_angle_types_registry)[(*molecule_angles_registry)[k][0]-1][2],
                       local_angle_types_registry[local_angles_registry[j][0]-1][0],
                       local_angle_types_registry[local_angles_registry[j][0]-1][1],
                       local_angle_types_registry[local_angles_registry[j][0]-1][2]
                       );
                */
                // check the type: is it the same?
                if((strcmp((*molecule_angle_types_registry)[(*molecule_angles_registry)[k][0]-1][0],local_angle_types_registry[local_angles_registry[j][0]-1][0])==0 &&
                    strcmp((*molecule_angle_types_registry)[(*molecule_angles_registry)[k][0]-1][1],local_angle_types_registry[local_angles_registry[j][0]-1][1])==0 &&
                    strcmp((*molecule_angle_types_registry)[(*molecule_angles_registry)[k][0]-1][2],local_angle_types_registry[local_angles_registry[j][0]-1][2])==0 &&
                    strcmp((*molecule_angle_types_registry)[(*molecule_angles_registry)[k][0]-1][3],local_angle_types_registry[local_angles_registry[j][0]-1][3])==0 &&
                    strcmp((*molecule_angle_types_registry)[(*molecule_angles_registry)[k][0]-1][4],local_angle_types_registry[local_angles_registry[j][0]-1][4])==0
                    )
                   ||
                   (strcmp((*molecule_angle_types_registry)[(*molecule_angles_registry)[k][0]-1][0],local_angle_types_registry[local_angles_registry[j][0]-1][4])==0 &&
                    strcmp((*molecule_angle_types_registry)[(*molecule_angles_registry)[k][0]-1][1],local_angle_types_registry[local_angles_registry[j][0]-1][3])==0 &&
                    strcmp((*molecule_angle_types_registry)[(*molecule_angles_registry)[k][0]-1][2],local_angle_types_registry[local_angles_registry[j][0]-1][2])==0 &&
                    strcmp((*molecule_angle_types_registry)[(*molecule_angles_registry)[k][0]-1][3],local_angle_types_registry[local_angles_registry[j][0]-1][1])==0 &&
                    strcmp((*molecule_angle_types_registry)[(*molecule_angles_registry)[k][0]-1][4],local_angle_types_registry[local_angles_registry[j][0]-1][0])==0
                    ))
                    // (check the type: is it the same?) YES!
                {
                    //printf("... (%d) same type - the angle is unaltered!!\n",(*molecule_angles_registry)[k][0]);
                    // set present flag --> 1 and break!
                    present=1;
                    break; // break the k-loop
                }
                else
                    // (check the type: is it the same?) NO!!
                {
                    // search molecule angle types: found?
                    for(n=0;n<*molecule_angle_types;++n)
                    {
                        /*
                        printf("... comparing with [%d] %s-%s-%s ...\n",n+1,
                               (*molecule_angle_types_registry)[n][0],
                               (*molecule_angle_types_registry)[n][1],
                               (*molecule_angle_types_registry)[n][2]);
                        */
                        if((strcmp((*molecule_angle_types_registry)[n][0],local_angle_types_registry[local_angles_registry[j][0]-1][0])==0 &&
                            strcmp((*molecule_angle_types_registry)[n][1],local_angle_types_registry[local_angles_registry[j][0]-1][1])==0 &&
                            strcmp((*molecule_angle_types_registry)[n][2],local_angle_types_registry[local_angles_registry[j][0]-1][2])==0 &&
                            strcmp((*molecule_angle_types_registry)[n][3],local_angle_types_registry[local_angles_registry[j][0]-1][3])==0 &&
                            strcmp((*molecule_angle_types_registry)[n][4],local_angle_types_registry[local_angles_registry[j][0]-1][4])==0
                            )
                           ||
                           (strcmp((*molecule_angle_types_registry)[n][0],local_angle_types_registry[local_angles_registry[j][0]-1][4])==0 &&
                            strcmp((*molecule_angle_types_registry)[n][1],local_angle_types_registry[local_angles_registry[j][0]-1][3])==0 &&
                            strcmp((*molecule_angle_types_registry)[n][2],local_angle_types_registry[local_angles_registry[j][0]-1][2])==0 &&
                            strcmp((*molecule_angle_types_registry)[n][3],local_angle_types_registry[local_angles_registry[j][0]-1][1])==0 &&
                            strcmp((*molecule_angle_types_registry)[n][4],local_angle_types_registry[local_angles_registry[j][0]-1][0])==0
                            ))
                            // (search molecule angle types: found?) YES!
                        {
                            // update type in the molecule angles registry
                            (*molecule_angles_registry)[k][0]=n+1;
                            // set present flag --> 1 and break!
                            present=1;
                            break; // break the n-loop
                        }
                    }
                    // (search molecule angle types: found?) NO!!
                    if(present==0)
                    {
                        // add new type to molecule angle types
                        /*
                        printf("+++ I will need to add type %s-%s-%s to the registry...\n",
                               local_angle_types_registry[local_angles_registry[j][0]-1][0],
                               local_angle_types_registry[local_angles_registry[j][0]-1][1],
                               local_angle_types_registry[local_angles_registry[j][0]-1][2]);
                        */
                        // augment molecule angle types
                        *molecule_angle_types=*molecule_angle_types+1;
                        // realloc
                        *molecule_angle_types_registry=(char***)realloc(*molecule_angle_types_registry,*molecule_angle_types*sizeof(char**));
                        (*molecule_angle_types_registry)[*molecule_angle_types-1]=(char**)malloc(5*sizeof(char*));
                        (*molecule_angle_types_registry)[*molecule_angle_types-1][0]=(char*)malloc(sub_length*sizeof(char));
                        (*molecule_angle_types_registry)[*molecule_angle_types-1][1]=(char*)malloc(sub_length*sizeof(char));
                        (*molecule_angle_types_registry)[*molecule_angle_types-1][2]=(char*)malloc(sub_length*sizeof(char));
                        (*molecule_angle_types_registry)[*molecule_angle_types-1][3]=(char*)malloc(sub_length*sizeof(char));
                        (*molecule_angle_types_registry)[*molecule_angle_types-1][4]=(char*)malloc(sub_length*sizeof(char));
                        // write type
                        sprintf((*molecule_angle_types_registry)[*molecule_angle_types-1][0],"%s",local_angle_types_registry[local_angles_registry[j][0]-1][0]);
                        sprintf((*molecule_angle_types_registry)[*molecule_angle_types-1][1],"%s",local_angle_types_registry[local_angles_registry[j][0]-1][1]);
                        sprintf((*molecule_angle_types_registry)[*molecule_angle_types-1][2],"%s",local_angle_types_registry[local_angles_registry[j][0]-1][2]);
                        sprintf((*molecule_angle_types_registry)[*molecule_angle_types-1][3],"%s",local_angle_types_registry[local_angles_registry[j][0]-1][3]);
                        sprintf((*molecule_angle_types_registry)[*molecule_angle_types-1][4],"%s",local_angle_types_registry[local_angles_registry[j][0]-1][4]);
                        
                        // update type in the molecule angles registry
                        (*molecule_angles_registry)[k][0]=*molecule_angle_types;
                        // check
                        //for(n=0;n<*molecule_angle_types;++n)printf("$$$ %s-%s-%s\n",(*molecule_angle_types_registry)[n][0],(*molecule_angle_types_registry)[n][1],(*molecule_angle_types_registry)[n][2]);
                        // set present flag --> 1 and break!
                        present=1;
                        break; // break the k-loop
                    }
                    else
                    {
                        // present is already equal to 1!!
                        break; // break the k-loop
                    }
                }
            }
        }
        if(present==0)
        {
            /*
            printf("*** %d-%d-%d~~\n",local_angles_registry[j][1]-1,local_angles_registry[j][2]-1,local_angles_registry[j][3]-1);
            
            printf("angle %d-%d-%d is NOT present!\n",local_map[local_angles_registry[j][1]-1],local_map[local_angles_registry[j][2]-1],local_map[local_angles_registry[j][3]-1]);
            
            printf("+++ its type is %s-%s-%s\n",local_angle_types_registry[local_angles_registry[j][0]-1][0],local_angle_types_registry[local_angles_registry[j][0]-1][1],local_angle_types_registry[local_angles_registry[j][0]-1][2]);
            */
            // add new angle to the molecule angles registry
            // augment
            *molecule_angles=*molecule_angles+1;
            // realloc
            *molecule_angles_registry=(int**)realloc(*molecule_angles_registry,*molecule_angles*sizeof(int*));
            (*molecule_angles_registry)[*molecule_angles-1]=(int*)malloc(4*sizeof(int));
            // populate: ids
            (*molecule_angles_registry)[*molecule_angles-1][1]=local_map[local_angles_registry[j][1]-1];
            (*molecule_angles_registry)[*molecule_angles-1][2]=local_map[local_angles_registry[j][2]-1];
            (*molecule_angles_registry)[*molecule_angles-1][3]=local_map[local_angles_registry[j][3]-1];
            // search molecule angle types: found?
            // initialize flag
            found=0;
            for(n=0;n<*molecule_angle_types;++n)
            {
                //printf("... comparing with [%d] %s-%s-%s ...\n",n+1,(*molecule_angle_types_registry)[n][0],(*molecule_angle_types_registry)[n][1],(*molecule_angle_types_registry)[n][2]);
                if((strcmp((*molecule_angle_types_registry)[n][0],local_angle_types_registry[local_angles_registry[j][0]-1][0])==0 &&
                    strcmp((*molecule_angle_types_registry)[n][1],local_angle_types_registry[local_angles_registry[j][0]-1][1])==0 &&
                    strcmp((*molecule_angle_types_registry)[n][2],local_angle_types_registry[local_angles_registry[j][0]-1][2])==0 &&
                    strcmp((*molecule_angle_types_registry)[n][3],local_angle_types_registry[local_angles_registry[j][0]-1][3])==0 &&
                    strcmp((*molecule_angle_types_registry)[n][4],local_angle_types_registry[local_angles_registry[j][0]-1][4])==0
                    )
                   ||
                   (strcmp((*molecule_angle_types_registry)[n][0],local_angle_types_registry[local_angles_registry[j][0]-1][4])==0 &&
                    strcmp((*molecule_angle_types_registry)[n][1],local_angle_types_registry[local_angles_registry[j][0]-1][3])==0 &&
                    strcmp((*molecule_angle_types_registry)[n][2],local_angle_types_registry[local_angles_registry[j][0]-1][2])==0 &&
                    strcmp((*molecule_angle_types_registry)[n][3],local_angle_types_registry[local_angles_registry[j][0]-1][1])==0 &&
                    strcmp((*molecule_angle_types_registry)[n][4],local_angle_types_registry[local_angles_registry[j][0]-1][0])==0
                    ))
                    // (search molecule angle types: found?) YES!
                {
                    // update type in the molecule angles registry
                    (*molecule_angles_registry)[*molecule_angles-1][0]=n+1;
                    // set found flag --> 1 and break!
                    found=1;
                    break; // break the n-loop
                }
            }
            // (search molecule angle types: found?) NO!!
            if(found==0)
            {
                /*
                // add new type to molecule angle types
                printf("+++ I will need to add type %s-%s-%s to the registry...\n",
                       local_angle_types_registry[local_angles_registry[j][0]-1][0],
                       local_angle_types_registry[local_angles_registry[j][0]-1][1],
                       local_angle_types_registry[local_angles_registry[j][0]-1][2]);
                */
                // augment molecule angle types
                *molecule_angle_types=*molecule_angle_types+1;
                // realloc
                *molecule_angle_types_registry=(char***)realloc(*molecule_angle_types_registry,*molecule_angle_types*sizeof(char**));
                (*molecule_angle_types_registry)[*molecule_angle_types-1]=(char**)malloc(5*sizeof(char*));
                (*molecule_angle_types_registry)[*molecule_angle_types-1][0]=(char*)malloc(sub_length*sizeof(char));
                (*molecule_angle_types_registry)[*molecule_angle_types-1][1]=(char*)malloc(sub_length*sizeof(char));
                (*molecule_angle_types_registry)[*molecule_angle_types-1][2]=(char*)malloc(sub_length*sizeof(char));
                (*molecule_angle_types_registry)[*molecule_angle_types-1][3]=(char*)malloc(sub_length*sizeof(char));
                (*molecule_angle_types_registry)[*molecule_angle_types-1][4]=(char*)malloc(sub_length*sizeof(char));
                // write type
                sprintf((*molecule_angle_types_registry)[*molecule_angle_types-1][0],"%s",local_angle_types_registry[local_angles_registry[j][0]-1][0]);
                sprintf((*molecule_angle_types_registry)[*molecule_angle_types-1][1],"%s",local_angle_types_registry[local_angles_registry[j][0]-1][1]);
                sprintf((*molecule_angle_types_registry)[*molecule_angle_types-1][2],"%s",local_angle_types_registry[local_angles_registry[j][0]-1][2]);
                sprintf((*molecule_angle_types_registry)[*molecule_angle_types-1][3],"%s",local_angle_types_registry[local_angles_registry[j][0]-1][3]);
                sprintf((*molecule_angle_types_registry)[*molecule_angle_types-1][4],"%s",local_angle_types_registry[local_angles_registry[j][0]-1][4]);
                // update type in the molecule angles registry
                (*molecule_angles_registry)[*molecule_angles-1][0]=*molecule_angle_types;
                // check
                //for(n=0;n<*molecule_angle_types;++n)printf("$$$ %s-%s-%s\n",(*molecule_angle_types_registry)[n][0],(*molecule_angle_types_registry)[n][1],(*molecule_angle_types_registry)[n][2]);
                
            }
        }
    }
    
    // !! WARNING !!
    // a serious issue with dihedral and improper registries:
    // if all core molecules have no dihedrals or impropers (e.g. CH4 or HOH molecules) in their initial form,
    // the following code segments should be altered, since the appropriate arrays are not preallocated!!
    
    // DIHEDRALS!!
    
    //printf("%d\t%d\n",local_dihedrals,*molecule_dihedrals);
    
    
    
    // local dihedrals loop
    for(j=0;j<local_dihedrals;++j)
    {
        // present flag --> 0
        present=0;
        // molecule dihedrals loop
        for(k=0;k<*molecule_dihedrals;++k)
        {
            // search for matching IDs: found?
            if(((*molecule_dihedrals_registry)[k][1]==local_map[local_dihedrals_registry[j][1]-1] &&
                (*molecule_dihedrals_registry)[k][2]==local_map[local_dihedrals_registry[j][2]-1] &&
                (*molecule_dihedrals_registry)[k][3]==local_map[local_dihedrals_registry[j][3]-1] &&
                (*molecule_dihedrals_registry)[k][4]==local_map[local_dihedrals_registry[j][4]-1]
                ) ||
               ((*molecule_dihedrals_registry)[k][1]==local_map[local_dihedrals_registry[j][4]-1] &&
                (*molecule_dihedrals_registry)[k][2]==local_map[local_dihedrals_registry[j][3]-1] &&
                (*molecule_dihedrals_registry)[k][3]==local_map[local_dihedrals_registry[j][2]-1] &&
                (*molecule_dihedrals_registry)[k][4]==local_map[local_dihedrals_registry[j][1]-1]
                )
               )
                // (search for matching IDs: found?) YES!
            {
                /*
                printf("dihedral %d-%d-%d-%d is present!\n",
                       local_map[local_dihedrals_registry[j][1]-1],
                       local_map[local_dihedrals_registry[j][2]-1],
                       local_map[local_dihedrals_registry[j][3]-1],
                       local_map[local_dihedrals_registry[j][4]-1]);
                
                printf("... checking dihedral type: molecule --> %s-%s-%s-%s | local --> %s-%s-%s-%s\n",
                       (*molecule_dihedral_types_registry)[(*molecule_dihedrals_registry)[k][0]-1][0],
                       (*molecule_dihedral_types_registry)[(*molecule_dihedrals_registry)[k][0]-1][1],
                       (*molecule_dihedral_types_registry)[(*molecule_dihedrals_registry)[k][0]-1][2],
                       (*molecule_dihedral_types_registry)[(*molecule_dihedrals_registry)[k][0]-1][3],
                       local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][0],
                       local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][1],
                       local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][2],
                       local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][3]
                       );
                */
                // check the type: is it the same?
                if((strcmp((*molecule_dihedral_types_registry)[(*molecule_dihedrals_registry)[k][0]-1][0],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][0])==0 &&
                    strcmp((*molecule_dihedral_types_registry)[(*molecule_dihedrals_registry)[k][0]-1][1],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][1])==0 &&
                    strcmp((*molecule_dihedral_types_registry)[(*molecule_dihedrals_registry)[k][0]-1][2],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][2])==0 &&
                    strcmp((*molecule_dihedral_types_registry)[(*molecule_dihedrals_registry)[k][0]-1][3],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][3])==0 &&
                    strcmp((*molecule_dihedral_types_registry)[(*molecule_dihedrals_registry)[k][0]-1][4],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][4])==0 &&
                    strcmp((*molecule_dihedral_types_registry)[(*molecule_dihedrals_registry)[k][0]-1][5],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][5])==0 &&
                    strcmp((*molecule_dihedral_types_registry)[(*molecule_dihedrals_registry)[k][0]-1][6],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][6])==0
                    )
                   ||
                   (strcmp((*molecule_dihedral_types_registry)[(*molecule_dihedrals_registry)[k][0]-1][0],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][6])==0 &&
                    strcmp((*molecule_dihedral_types_registry)[(*molecule_dihedrals_registry)[k][0]-1][1],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][5])==0 &&
                    strcmp((*molecule_dihedral_types_registry)[(*molecule_dihedrals_registry)[k][0]-1][2],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][4])==0 &&
                    strcmp((*molecule_dihedral_types_registry)[(*molecule_dihedrals_registry)[k][0]-1][3],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][3])==0 &&
                    strcmp((*molecule_dihedral_types_registry)[(*molecule_dihedrals_registry)[k][0]-1][4],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][2])==0 &&
                    strcmp((*molecule_dihedral_types_registry)[(*molecule_dihedrals_registry)[k][0]-1][5],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][1])==0 &&
                    strcmp((*molecule_dihedral_types_registry)[(*molecule_dihedrals_registry)[k][0]-1][6],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][0])==0
                    ))
                    // (check the type: is it the same?) YES!
                {
                    //printf("... (%d) same type - the dihedral is unaltered!!\n",(*molecule_dihedrals_registry)[k][0]);
                    // set present flag --> 1 and break!
                    present=1;
                    break; // break the k-loop
                }
                else
                    // (check the type: is it the same?) NO!!
                {
                    // search molecule dihedral types: found?
                    for(n=0;n<*molecule_dihedral_types;++n)
                    {
                        /*
                        printf("... comparing with [%d] %s-%s-%s-%s ...\n",n+1,
                               (*molecule_dihedral_types_registry)[n][0],
                               (*molecule_dihedral_types_registry)[n][1],
                               (*molecule_dihedral_types_registry)[n][2],
                               (*molecule_dihedral_types_registry)[n][3]);
                        */
                        if((strcmp((*molecule_dihedral_types_registry)[n][0],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][0])==0 &&
                            strcmp((*molecule_dihedral_types_registry)[n][1],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][1])==0 &&
                            strcmp((*molecule_dihedral_types_registry)[n][2],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][2])==0 &&
                            strcmp((*molecule_dihedral_types_registry)[n][3],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][3])==0 &&
                            strcmp((*molecule_dihedral_types_registry)[n][4],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][4])==0 &&
                            strcmp((*molecule_dihedral_types_registry)[n][5],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][5])==0 &&
                            strcmp((*molecule_dihedral_types_registry)[n][6],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][6])==0
                            )
                           ||
                           (strcmp((*molecule_dihedral_types_registry)[n][0],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][6])==0 &&
                            strcmp((*molecule_dihedral_types_registry)[n][1],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][5])==0 &&
                            strcmp((*molecule_dihedral_types_registry)[n][2],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][4])==0 &&
                            strcmp((*molecule_dihedral_types_registry)[n][3],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][3])==0 &&
                            strcmp((*molecule_dihedral_types_registry)[n][4],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][2])==0 &&
                            strcmp((*molecule_dihedral_types_registry)[n][5],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][1])==0 &&
                            strcmp((*molecule_dihedral_types_registry)[n][6],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][0])==0
                            ))
                            // (search molecule dihedral types: found?) YES!
                        {
                            // update type in the molecule dihedrals registry
                            (*molecule_dihedrals_registry)[k][0]=n+1;
                            // set present flag --> 1 and break!
                            present=1;
                            break; // break the n-loop
                        }
                    }
                    // (search molecule dihedral types: found?) NO!!
                    if(present==0)
                    {
                        // add new type to molecule dihedral types
                        /*
                        printf("+++ I will need to add type %s-%s-%s-%s to the registry...\n",
                               local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][0],
                               local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][1],
                               local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][2],
                               local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][3]);
                        */
                        // augment molecule dihedral types
                        *molecule_dihedral_types=*molecule_dihedral_types+1;
                        // realloc
                        *molecule_dihedral_types_registry=(char***)realloc(*molecule_dihedral_types_registry,*molecule_dihedral_types*sizeof(char**));
                        (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1]=(char**)malloc(7*sizeof(char*));
                        (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][0]=(char*)malloc(sub_length*sizeof(char));
                        (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][1]=(char*)malloc(sub_length*sizeof(char));
                        (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][2]=(char*)malloc(sub_length*sizeof(char));
                        (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][3]=(char*)malloc(sub_length*sizeof(char));
                        (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][4]=(char*)malloc(sub_length*sizeof(char));
                        (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][5]=(char*)malloc(sub_length*sizeof(char));
                        (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][6]=(char*)malloc(sub_length*sizeof(char));
                        // write type
                        sprintf((*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][0],"%s",local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][0]);
                        sprintf((*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][1],"%s",local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][1]);
                        sprintf((*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][2],"%s",local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][2]);
                        sprintf((*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][3],"%s",local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][3]);
                        sprintf((*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][4],"%s",local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][4]);
                        sprintf((*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][5],"%s",local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][5]);
                        sprintf((*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][6],"%s",local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][6]);
                        
                        // update type in the molecule dihedrals registry
                        (*molecule_dihedrals_registry)[k][0]=*molecule_dihedral_types;
                        // check
                        //for(n=0;n<*molecule_dihedral_types;++n)printf("$$$ %s-%s-%s-%s\n",(*molecule_dihedral_types_registry)[n][0],(*molecule_dihedral_types_registry)[n][1],(*molecule_dihedral_types_registry)[n][2],(*molecule_dihedral_types_registry)[n][3]);
                        // set present flag --> 1 and break!
                        present=1;
                        break; // break the k-loop
                    }
                    else
                    {
                        // present is already equal to 1!!
                        break; // break the k-loop
                    }
                }
            }
        }
        if(present==0)
        {
            
            //printf("*** %d-%d-%d-%d\n",local_dihedrals_registry[j][1]-1,local_dihedrals_registry[j][2]-1,local_dihedrals_registry[j][3]-1,local_dihedrals_registry[j][4]-1);
            
            //printf("dihedral %d-%d-%d-%d is NOT present!\n",local_map[local_dihedrals_registry[j][1]-1],local_map[local_dihedrals_registry[j][2]-1],local_map[local_dihedrals_registry[j][3]-1],local_map[local_dihedrals_registry[j][4]-1]);
            
            //printf("+++ its type is %s-%s-%s-%s\n",local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][0],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][1],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][2],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][3]);
            
            // add new dihedral to the molecule dihedrals registry
            // augment
            
            //printf("%d\n",*molecule_dihedrals);
            
            if(*molecule_dihedrals==0)
            {
                *molecule_dihedrals=*molecule_dihedrals+1;
                // alloc!!
                *molecule_dihedrals_registry=(int**)malloc(*molecule_dihedrals*sizeof(int*));
                (*molecule_dihedrals_registry)[*molecule_dihedrals-1]=(int*)malloc(5*sizeof(int));
            }
            else
            {
                *molecule_dihedrals=*molecule_dihedrals+1;
                // realloc
                *molecule_dihedrals_registry=(int**)realloc(*molecule_dihedrals_registry,*molecule_dihedrals*sizeof(int*));
                (*molecule_dihedrals_registry)[*molecule_dihedrals-1]=(int*)malloc(5*sizeof(int));
            }
            
            
            
            // populate: ids
            (*molecule_dihedrals_registry)[*molecule_dihedrals-1][1]=local_map[local_dihedrals_registry[j][1]-1];
            (*molecule_dihedrals_registry)[*molecule_dihedrals-1][2]=local_map[local_dihedrals_registry[j][2]-1];
            (*molecule_dihedrals_registry)[*molecule_dihedrals-1][3]=local_map[local_dihedrals_registry[j][3]-1];
            (*molecule_dihedrals_registry)[*molecule_dihedrals-1][4]=local_map[local_dihedrals_registry[j][4]-1];
            // search molecule dihedral types: found?
            // initialize flag
            found=0;
            for(n=0;n<*molecule_dihedral_types;++n)
            {
                //printf("... comparing with [%d] %s-%s-%s-%s ...\n",n+1,(*molecule_dihedral_types_registry)[n][0],(*molecule_dihedral_types_registry)[n][1],(*molecule_dihedral_types_registry)[n][2],(*molecule_dihedral_types_registry)[n][3]);
                if((strcmp((*molecule_dihedral_types_registry)[n][0],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][0])==0 &&
                    strcmp((*molecule_dihedral_types_registry)[n][1],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][1])==0 &&
                    strcmp((*molecule_dihedral_types_registry)[n][2],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][2])==0 &&
                    strcmp((*molecule_dihedral_types_registry)[n][3],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][3])==0 &&
                    strcmp((*molecule_dihedral_types_registry)[n][4],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][4])==0 &&
                    strcmp((*molecule_dihedral_types_registry)[n][5],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][5])==0 &&
                    strcmp((*molecule_dihedral_types_registry)[n][6],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][6])==0
                    )
                   ||
                   (strcmp((*molecule_dihedral_types_registry)[n][0],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][6])==0 &&
                    strcmp((*molecule_dihedral_types_registry)[n][1],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][5])==0 &&
                    strcmp((*molecule_dihedral_types_registry)[n][2],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][4])==0 &&
                    strcmp((*molecule_dihedral_types_registry)[n][3],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][3])==0 &&
                    strcmp((*molecule_dihedral_types_registry)[n][4],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][2])==0 &&
                    strcmp((*molecule_dihedral_types_registry)[n][5],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][1])==0 &&
                    strcmp((*molecule_dihedral_types_registry)[n][6],local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][0])==0
                    ))
                    // (search molecule dihedral types: found?) YES!
                {
                    // update type in the molecule dihedrals registry
                    (*molecule_dihedrals_registry)[*molecule_dihedrals-1][0]=n+1;
                    // set found flag --> 1 and break!
                    found=1;
                    break; // break the n-loop
                }
            }
            // (search molecule dihedral types: found?) NO!!

            //printf("%d\n",*molecule_dihedral_types);
            
            if(found==0)
            {
                // add new type to molecule dihedral types
                /*
                printf("+++ I will need to add type %s-%s-%s-%s to the registry...\n",
                       local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][0],
                       local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][1],
                       local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][2],
                       local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][3]);
                */
                if(*molecule_dihedral_types==0)
                {
                    // augment molecule dihedral types
                    *molecule_dihedral_types=*molecule_dihedral_types+1;
                    // alloc!!
                    *molecule_dihedral_types_registry=(char***)malloc(*molecule_dihedral_types*sizeof(char**));
                    (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1]=(char**)malloc(7*sizeof(char*));
                    (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][0]=(char*)malloc(sub_length*sizeof(char));
                    (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][1]=(char*)malloc(sub_length*sizeof(char));
                    (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][2]=(char*)malloc(sub_length*sizeof(char));
                    (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][3]=(char*)malloc(sub_length*sizeof(char));
                    (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][4]=(char*)malloc(sub_length*sizeof(char));
                    (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][5]=(char*)malloc(sub_length*sizeof(char));
                    (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][6]=(char*)malloc(sub_length*sizeof(char));
                }
                else
                {
                    // augment molecule dihedral types
                    *molecule_dihedral_types=*molecule_dihedral_types+1;
                    // realloc
                    *molecule_dihedral_types_registry=(char***)realloc(*molecule_dihedral_types_registry,*molecule_dihedral_types*sizeof(char**));
                    (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1]=(char**)malloc(7*sizeof(char*));
                    (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][0]=(char*)malloc(sub_length*sizeof(char));
                    (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][1]=(char*)malloc(sub_length*sizeof(char));
                    (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][2]=(char*)malloc(sub_length*sizeof(char));
                    (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][3]=(char*)malloc(sub_length*sizeof(char));
                    (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][4]=(char*)malloc(sub_length*sizeof(char));
                    (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][5]=(char*)malloc(sub_length*sizeof(char));
                    (*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][6]=(char*)malloc(sub_length*sizeof(char));
                }

                
                
                
                // write type
                sprintf((*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][0],"%s",local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][0]);
                sprintf((*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][1],"%s",local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][1]);
                sprintf((*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][2],"%s",local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][2]);
                sprintf((*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][3],"%s",local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][3]);
                sprintf((*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][4],"%s",local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][4]);
                sprintf((*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][5],"%s",local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][5]);
                sprintf((*molecule_dihedral_types_registry)[*molecule_dihedral_types-1][6],"%s",local_dihedral_types_registry[local_dihedrals_registry[j][0]-1][6]);
                // update type in the molecule dihedrals registry
                (*molecule_dihedrals_registry)[*molecule_dihedrals-1][0]=*molecule_dihedral_types;
                // check
                //for(n=0;n<*molecule_dihedral_types;++n)printf("$$$ %s-%s-%s-%s\n",(*molecule_dihedral_types_registry)[n][0],(*molecule_dihedral_types_registry)[n][1],(*molecule_dihedral_types_registry)[n][2],(*molecule_dihedral_types_registry)[n][3]);
                
            }
        }
    }
    
    // IMPROPERS!!
    
    // local impropers loop
    for(j=0;j<local_impropers;++j)
    {
        // present flag --> 0
        present=0;
        // molecule impropers loop
        for(k=0;k<*molecule_impropers;++k)
        {
            // search for matching IDs: found?
            if(((*molecule_impropers_registry)[k][1]==local_map[local_impropers_registry[j][1]-1] &&
                (*molecule_impropers_registry)[k][2]==local_map[local_impropers_registry[j][2]-1] &&
                (*molecule_impropers_registry)[k][3]==local_map[local_impropers_registry[j][3]-1] &&
                (*molecule_impropers_registry)[k][4]==local_map[local_impropers_registry[j][4]-1]
                ) ||
               ((*molecule_impropers_registry)[k][1]==local_map[local_impropers_registry[j][4]-1] &&
                (*molecule_impropers_registry)[k][2]==local_map[local_impropers_registry[j][3]-1] &&
                (*molecule_impropers_registry)[k][3]==local_map[local_impropers_registry[j][2]-1] &&
                (*molecule_impropers_registry)[k][4]==local_map[local_impropers_registry[j][1]-1]
                )
               )
                // (search for matching IDs: found?) YES!
            {
                /*
                printf("improper %d-%d-%d-%d is present!\n",
                       local_map[local_impropers_registry[j][1]-1],
                       local_map[local_impropers_registry[j][2]-1],
                       local_map[local_impropers_registry[j][3]-1],
                       local_map[local_impropers_registry[j][4]-1]);
                
                printf("... checking improper type: molecule --> %s-%s-%s-%s | local --> %s-%s-%s-%s\n",
                       (*molecule_improper_types_registry)[(*molecule_impropers_registry)[k][0]-1][0],
                       (*molecule_improper_types_registry)[(*molecule_impropers_registry)[k][0]-1][1],
                       (*molecule_improper_types_registry)[(*molecule_impropers_registry)[k][0]-1][2],
                       (*molecule_improper_types_registry)[(*molecule_impropers_registry)[k][0]-1][3],
                       local_improper_types_registry[local_impropers_registry[j][0]-1][0],
                       local_improper_types_registry[local_impropers_registry[j][0]-1][1],
                       local_improper_types_registry[local_impropers_registry[j][0]-1][2],
                       local_improper_types_registry[local_impropers_registry[j][0]-1][3]
                       );
                */
                // check the type: is it the same?
                if((strcmp((*molecule_improper_types_registry)[(*molecule_impropers_registry)[k][0]-1][0],local_improper_types_registry[local_impropers_registry[j][0]-1][0])==0 &&
                    strcmp((*molecule_improper_types_registry)[(*molecule_impropers_registry)[k][0]-1][1],local_improper_types_registry[local_impropers_registry[j][0]-1][1])==0 &&
                    strcmp((*molecule_improper_types_registry)[(*molecule_impropers_registry)[k][0]-1][2],local_improper_types_registry[local_impropers_registry[j][0]-1][2])==0 &&
                    strcmp((*molecule_improper_types_registry)[(*molecule_impropers_registry)[k][0]-1][3],local_improper_types_registry[local_impropers_registry[j][0]-1][3])==0
                    )
                   ||
                   (strcmp((*molecule_improper_types_registry)[(*molecule_impropers_registry)[k][0]-1][0],local_improper_types_registry[local_impropers_registry[j][0]-1][3])==0 &&
                    strcmp((*molecule_improper_types_registry)[(*molecule_impropers_registry)[k][0]-1][1],local_improper_types_registry[local_impropers_registry[j][0]-1][2])==0 &&
                    strcmp((*molecule_improper_types_registry)[(*molecule_impropers_registry)[k][0]-1][2],local_improper_types_registry[local_impropers_registry[j][0]-1][1])==0 &&
                    strcmp((*molecule_improper_types_registry)[(*molecule_impropers_registry)[k][0]-1][3],local_improper_types_registry[local_impropers_registry[j][0]-1][0])==0
                    ))
                    // (check the type: is it the same?) YES!
                {
                    //printf("... (%d) same type - the improper is unaltered!!\n",(*molecule_impropers_registry)[k][0]);
                    // set present flag --> 1 and break!
                    present=1;
                    break; // break the k-loop
                }
                else
                    // (check the type: is it the same?) NO!!
                {
                    // search molecule improper types: found?
                    for(n=0;n<*molecule_improper_types;++n)
                    {
                        /*
                        printf("... comparing with [%d] %s-%s-%s-%s ...\n",n+1,
                               (*molecule_improper_types_registry)[n][0],
                               (*molecule_improper_types_registry)[n][1],
                               (*molecule_improper_types_registry)[n][2],
                               (*molecule_improper_types_registry)[n][3]);
                        */
                        if((strcmp((*molecule_improper_types_registry)[n][0],local_improper_types_registry[local_impropers_registry[j][0]-1][0])==0 &&
                            strcmp((*molecule_improper_types_registry)[n][1],local_improper_types_registry[local_impropers_registry[j][0]-1][1])==0 &&
                            strcmp((*molecule_improper_types_registry)[n][2],local_improper_types_registry[local_impropers_registry[j][0]-1][2])==0 &&
                            strcmp((*molecule_improper_types_registry)[n][3],local_improper_types_registry[local_impropers_registry[j][0]-1][3])==0
                            )
                           ||
                           (strcmp((*molecule_improper_types_registry)[n][0],local_improper_types_registry[local_impropers_registry[j][0]-1][3])==0 &&
                            strcmp((*molecule_improper_types_registry)[n][1],local_improper_types_registry[local_impropers_registry[j][0]-1][2])==0 &&
                            strcmp((*molecule_improper_types_registry)[n][2],local_improper_types_registry[local_impropers_registry[j][0]-1][1])==0 &&
                            strcmp((*molecule_improper_types_registry)[n][3],local_improper_types_registry[local_impropers_registry[j][0]-1][0])==0
                            ))
                            // (search molecule improper types: found?) YES!
                        {
                            // update type in the molecule impropers registry
                            (*molecule_impropers_registry)[k][0]=n+1;
                            // set present flag --> 1 and break!
                            present=1;
                            break; // break the n-loop
                        }
                    }
                    // (search molecule improper types: found?) NO!!
                    if(present==0)
                    {
                        // add new type to molecule improper types
                        /*
                        printf("+++ I will need to add type %s-%s-%s-%s to the registry...\n",
                               local_improper_types_registry[local_impropers_registry[j][0]-1][0],
                               local_improper_types_registry[local_impropers_registry[j][0]-1][1],
                               local_improper_types_registry[local_impropers_registry[j][0]-1][2],
                               local_improper_types_registry[local_impropers_registry[j][0]-1][3]);
                        */
                        // augment molecule improper types
                        *molecule_improper_types=*molecule_improper_types+1;
                        // realloc
                        *molecule_improper_types_registry=(char***)realloc(*molecule_improper_types_registry,*molecule_improper_types*sizeof(char**));
                        (*molecule_improper_types_registry)[*molecule_improper_types-1]=(char**)malloc(4*sizeof(char*));
                        (*molecule_improper_types_registry)[*molecule_improper_types-1][0]=(char*)malloc(sub_length*sizeof(char));
                        (*molecule_improper_types_registry)[*molecule_improper_types-1][1]=(char*)malloc(sub_length*sizeof(char));
                        (*molecule_improper_types_registry)[*molecule_improper_types-1][2]=(char*)malloc(sub_length*sizeof(char));
                        (*molecule_improper_types_registry)[*molecule_improper_types-1][3]=(char*)malloc(sub_length*sizeof(char));
                        // write type
                        sprintf((*molecule_improper_types_registry)[*molecule_improper_types-1][0],"%s",local_improper_types_registry[local_impropers_registry[j][0]-1][0]);
                        sprintf((*molecule_improper_types_registry)[*molecule_improper_types-1][1],"%s",local_improper_types_registry[local_impropers_registry[j][0]-1][1]);
                        sprintf((*molecule_improper_types_registry)[*molecule_improper_types-1][2],"%s",local_improper_types_registry[local_impropers_registry[j][0]-1][2]);
                        sprintf((*molecule_improper_types_registry)[*molecule_improper_types-1][3],"%s",local_improper_types_registry[local_impropers_registry[j][0]-1][3]);
                        
                        // update type in the molecule impropers registry
                        (*molecule_impropers_registry)[k][0]=*molecule_improper_types;
                        // check
                        //for(n=0;n<*molecule_improper_types;++n)printf("$$$ %s-%s-%s-%s\n",(*molecule_improper_types_registry)[n][0],(*molecule_improper_types_registry)[n][1],(*molecule_improper_types_registry)[n][2],(*molecule_improper_types_registry)[n][3]);
                        // set present flag --> 1 and break!
                        present=1;
                        break; // break the k-loop
                    }
                    else
                    {
                        // present is already equal to 1!!
                        break; // break the k-loop
                    }
                }
            }
        }
        if(present==0)
        {
            
            //printf("*** %d-%d-%d-%d\n",local_impropers_registry[j][1]-1,local_impropers_registry[j][2]-1,local_impropers_registry[j][3]-1,local_impropers_registry[j][4]-1);
            
            //printf("improper %d-%d-%d-%d is NOT present!\n",local_map[local_impropers_registry[j][1]-1],local_map[local_impropers_registry[j][2]-1],local_map[local_impropers_registry[j][3]-1],local_map[local_impropers_registry[j][4]-1]);
            
            //printf("+++ its type is %s-%s-%s-%s\n",local_improper_types_registry[local_impropers_registry[j][0]-1][0],local_improper_types_registry[local_impropers_registry[j][0]-1][1],local_improper_types_registry[local_impropers_registry[j][0]-1][2],local_improper_types_registry[local_impropers_registry[j][0]-1][3]);
            
            if(*molecule_impropers==0)
            {
                *molecule_impropers=*molecule_impropers+1;
                // alloc!!
                *molecule_impropers_registry=(int**)malloc(*molecule_impropers*sizeof(int*));
                (*molecule_impropers_registry)[*molecule_impropers-1]=(int*)malloc(5*sizeof(int));
            }
            else
            {
                *molecule_impropers=*molecule_impropers+1;
                // realloc
                *molecule_impropers_registry=(int**)realloc(*molecule_impropers_registry,*molecule_impropers*sizeof(int*));
                (*molecule_impropers_registry)[*molecule_impropers-1]=(int*)malloc(5*sizeof(int));
            }
            
            // add new improper to the molecule impropers registry
            // augment
            //*molecule_impropers=*molecule_impropers+1;
            // realloc
            //*molecule_impropers_registry=(int**)realloc(*molecule_impropers_registry,*molecule_impropers*sizeof(int*));
            //(*molecule_impropers_registry)[*molecule_impropers-1]=(int*)malloc(5*sizeof(int));
            
            
            // populate: ids
            (*molecule_impropers_registry)[*molecule_impropers-1][1]=local_map[local_impropers_registry[j][1]-1];
            (*molecule_impropers_registry)[*molecule_impropers-1][2]=local_map[local_impropers_registry[j][2]-1];
            (*molecule_impropers_registry)[*molecule_impropers-1][3]=local_map[local_impropers_registry[j][3]-1];
            (*molecule_impropers_registry)[*molecule_impropers-1][4]=local_map[local_impropers_registry[j][4]-1];
            // search molecule improper types: found?
            // initialize flag
            found=0;
            for(n=0;n<*molecule_improper_types;++n)
            {
                //printf("... comparing with [%d] %s-%s-%s-%s ...\n",n+1,(*molecule_improper_types_registry)[n][0],(*molecule_improper_types_registry)[n][1],(*molecule_improper_types_registry)[n][2],(*molecule_improper_types_registry)[n][3]);
                if((strcmp((*molecule_improper_types_registry)[n][0],local_improper_types_registry[local_impropers_registry[j][0]-1][0])==0 &&
                    strcmp((*molecule_improper_types_registry)[n][1],local_improper_types_registry[local_impropers_registry[j][0]-1][1])==0 &&
                    strcmp((*molecule_improper_types_registry)[n][2],local_improper_types_registry[local_impropers_registry[j][0]-1][2])==0 &&
                    strcmp((*molecule_improper_types_registry)[n][3],local_improper_types_registry[local_impropers_registry[j][0]-1][3])==0
                    )
                   ||
                   (strcmp((*molecule_improper_types_registry)[n][0],local_improper_types_registry[local_impropers_registry[j][0]-1][3])==0 &&
                    strcmp((*molecule_improper_types_registry)[n][1],local_improper_types_registry[local_impropers_registry[j][0]-1][2])==0 &&
                    strcmp((*molecule_improper_types_registry)[n][2],local_improper_types_registry[local_impropers_registry[j][0]-1][1])==0 &&
                    strcmp((*molecule_improper_types_registry)[n][3],local_improper_types_registry[local_impropers_registry[j][0]-1][0])==0
                    ))
                    // (search molecule improper types: found?) YES!
                {
                    // update type in the molecule impropers registry
                    (*molecule_impropers_registry)[*molecule_impropers-1][0]=n+1;
                    // set found flag --> 1 and break!
                    found=1;
                    break; // break the n-loop
                }
            }
            // (search molecule improper types: found?) NO!!
            if(found==0)
            {
                // add new type to molecule improper types
                /*
                printf("+++ I will need to add type %s-%s-%s-%s to the registry...\n",
                       local_improper_types_registry[local_impropers_registry[j][0]-1][0],
                       local_improper_types_registry[local_impropers_registry[j][0]-1][1],
                       local_improper_types_registry[local_impropers_registry[j][0]-1][2],
                       local_improper_types_registry[local_impropers_registry[j][0]-1][3]);
                */
                
                if(*molecule_improper_types==0)
                {
                    // augment molecule dihedral types
                    *molecule_improper_types=*molecule_improper_types+1;
                    // alloc!!
                    *molecule_improper_types_registry=(char***)malloc(*molecule_improper_types*sizeof(char**));
                    (*molecule_improper_types_registry)[*molecule_improper_types-1]=(char**)malloc(4*sizeof(char*));
                    (*molecule_improper_types_registry)[*molecule_improper_types-1][0]=(char*)malloc(sub_length*sizeof(char));
                    (*molecule_improper_types_registry)[*molecule_improper_types-1][1]=(char*)malloc(sub_length*sizeof(char));
                    (*molecule_improper_types_registry)[*molecule_improper_types-1][2]=(char*)malloc(sub_length*sizeof(char));
                    (*molecule_improper_types_registry)[*molecule_improper_types-1][3]=(char*)malloc(sub_length*sizeof(char));
                }
                else
                {
                    // augment molecule dihedral types
                    *molecule_improper_types=*molecule_improper_types+1;
                    // realloc
                    *molecule_improper_types_registry=(char***)realloc(*molecule_improper_types_registry,*molecule_improper_types*sizeof(char**));
                    (*molecule_improper_types_registry)[*molecule_improper_types-1]=(char**)malloc(4*sizeof(char*));
                    (*molecule_improper_types_registry)[*molecule_improper_types-1][0]=(char*)malloc(sub_length*sizeof(char));
                    (*molecule_improper_types_registry)[*molecule_improper_types-1][1]=(char*)malloc(sub_length*sizeof(char));
                    (*molecule_improper_types_registry)[*molecule_improper_types-1][2]=(char*)malloc(sub_length*sizeof(char));
                    (*molecule_improper_types_registry)[*molecule_improper_types-1][3]=(char*)malloc(sub_length*sizeof(char));
                }
                
                // augment molecule improper types
                //*molecule_improper_types=*molecule_improper_types+1;
                // realloc
                //*molecule_improper_types_registry=(char***)realloc(*molecule_improper_types_registry,*molecule_improper_types*sizeof(char**));
                //(*molecule_improper_types_registry)[*molecule_improper_types-1]=(char**)malloc(4*sizeof(char*));
                //(*molecule_improper_types_registry)[*molecule_improper_types-1][0]=(char*)malloc(sub_length*sizeof(char));
                //(*molecule_improper_types_registry)[*molecule_improper_types-1][1]=(char*)malloc(sub_length*sizeof(char));
                //(*molecule_improper_types_registry)[*molecule_improper_types-1][2]=(char*)malloc(sub_length*sizeof(char));
                //(*molecule_improper_types_registry)[*molecule_improper_types-1][3]=(char*)malloc(sub_length*sizeof(char));
                // write type
                sprintf((*molecule_improper_types_registry)[*molecule_improper_types-1][0],"%s",local_improper_types_registry[local_impropers_registry[j][0]-1][0]);
                sprintf((*molecule_improper_types_registry)[*molecule_improper_types-1][1],"%s",local_improper_types_registry[local_impropers_registry[j][0]-1][1]);
                sprintf((*molecule_improper_types_registry)[*molecule_improper_types-1][2],"%s",local_improper_types_registry[local_impropers_registry[j][0]-1][2]);
                sprintf((*molecule_improper_types_registry)[*molecule_improper_types-1][3],"%s",local_improper_types_registry[local_impropers_registry[j][0]-1][3]);
                // update type in the molecule impropers registry
                (*molecule_impropers_registry)[*molecule_impropers-1][0]=*molecule_improper_types;
                // check
                //for(n=0;n<*molecule_improper_types;++n)printf("$$$ %s-%s-%s-%s\n",(*molecule_improper_types_registry)[n][0],(*molecule_improper_types_registry)[n][1],(*molecule_improper_types_registry)[n][2],(*molecule_improper_types_registry)[n][3]);
                
            }
        }
    }
    
    //
    /*
    printf("\n# This is output from alter_molecular_topo()\n");
    
    printf("\n");
    
    //printf("Number of atoms = %d\n",*molecule_atoms);
    printf("Number of atoms = N/A\n");
    printf("Number of bonds = %d\n",*molecule_bonds);
    printf("Number of angles = %d\n",*molecule_angles);
    printf("Number of dihedrals = %d\n",*molecule_dihedrals);
    printf("Number of impropers = %d\n",*molecule_impropers);
    
    printf("\n");
    
    printf("Number of atom types = %d\n",*molecule_atom_types);
    printf("Number of bond types = %d\n",*molecule_bond_types);
    printf("Number of angle types = %d\n",*molecule_angle_types);
    if(*molecule_dihedrals>0)printf("Number of dihedral types = %d\n",*molecule_dihedral_types);
    if(*molecule_impropers>0)printf("Number of improper types = %d\n",*molecule_improper_types);
    
    printf("\nAtom types\n\n");
    
    for(j=0;j<*molecule_atom_types;++j)printf("[%d]\t%s\n",j+1,(*molecule_species_registry)[j]);
    
    printf("\nBond types\n\n");
    
    for(j=0;j<*molecule_bond_types;++j)printf("[%d]\t%s\t%s\n",j+1,(*molecule_bond_types_registry)[j][0],(*molecule_bond_types_registry)[j][1]);
    
    printf("\nAngle types\n\n");
    
    for(j=0;j<*molecule_angle_types;++j)printf("[%d]\t%s\t%s\t%s\n",j+1,(*molecule_angle_types_registry)[j][0],(*molecule_angle_types_registry)[j][1],(*molecule_angle_types_registry)[j][2]);
    if(*molecule_dihedrals>0)
    {
        
        printf("\nDihedral types\n\n");
        
        for(j=0;j<*molecule_dihedral_types;++j)printf("[%d]\t%s\t%s\t%s\t%s\n",j+1,(*molecule_dihedral_types_registry)[j][0],(*molecule_dihedral_types_registry)[j][1],(*molecule_dihedral_types_registry)[j][2],(*molecule_dihedral_types_registry)[j][3]);
    }
    if(*molecule_impropers>0)
    {
        
        printf("\nImproper types\n\n");
        
        for(j=0;j<*molecule_improper_types;++j)printf("[%d]\t%s\t%s\t%s\t%s\n",j+1,(*molecule_improper_types_registry)[j][0],(*molecule_improper_types_registry)[j][1],(*molecule_improper_types_registry)[j][2],(*molecule_improper_types_registry)[j][3]);
    }
    
    printf("\nBonds\n\n");
    
    for(j=0;j<*molecule_bonds;++j)printf("[%d]\t%d\t%d\t%d\n",j+1,(*molecule_bonds_registry)[j][0],(*molecule_bonds_registry)[j][1],(*molecule_bonds_registry)[j][2]);
    
    printf("\nAngles\n\n");
    
    for(j=0;j<*molecule_angles;++j)printf("[%d]\t%d\t%d\t%d\t%d\n",j+1,(*molecule_angles_registry)[j][0],(*molecule_angles_registry)[j][1],(*molecule_angles_registry)[j][2],(*molecule_angles_registry)[j][3]);
    if(*molecule_dihedrals>0)
    {
        
        printf("\nDihedrals\n\n");
        
        for(j=0;j<*molecule_dihedrals;++j)printf("[%d]\t%d\t%d\t%d\t%d\t%d\n",j+1,(*molecule_dihedrals_registry)[j][0],(*molecule_dihedrals_registry)[j][1],(*molecule_dihedrals_registry)[j][2],(*molecule_dihedrals_registry)[j][3],(*molecule_dihedrals_registry)[j][4]);
    }
    if(*molecule_impropers>0)
    {
        
        printf("\nImpropers\n\n");
        
        for(j=0;j<*molecule_impropers;++j)printf("[%d]\t%d\t%d\t%d\t%d\t%d\n",j+1,(*molecule_impropers_registry)[j][0],(*molecule_impropers_registry)[j][1],(*molecule_impropers_registry)[j][2],(*molecule_impropers_registry)[j][3],(*molecule_impropers_registry)[j][4]);
    }
    */
    //
    if(verb==1){
    printf("\n$ New species:\n");
    for(j=0;j<*molecule_atom_types;++j)printf("[%d]\t%s\n",j+1,(*molecule_species_registry)[j]);
    
    printf("\n$ New bonds topology:\n");
    printf("\n$ Bond types:\n");
    for(j=0;j<*molecule_bond_types;++j)printf("[%d]\t(%s)--%s--(%s)\n",j+1,(*molecule_bond_types_registry)[j][0],(*molecule_bond_types_registry)[j][1],(*molecule_bond_types_registry)[j][2]);
    printf("\n$ Bonds:\n");
    for(j=0;j<*molecule_bonds;++j)printf("[%d]\t%d\t%d\t%d\n",j+1,(*molecule_bonds_registry)[j][0],(*molecule_bonds_registry)[j][1],(*molecule_bonds_registry)[j][2]);
    
    printf("\n$ New angles topology:\n");
    printf("\n$ Angle types:\n");
    for(j=0;j<*molecule_angle_types;++j)printf("[%d]\t(%s)--%s--(%s)--%s--(%s)\n",j+1,(*molecule_angle_types_registry)[j][0],(*molecule_angle_types_registry)[j][1],(*molecule_angle_types_registry)[j][2],(*molecule_angle_types_registry)[j][3],(*molecule_angle_types_registry)[j][4]);
    printf("\n$ Angles:\n");
    for(j=0;j<*molecule_angles;++j)printf("[%d]\t%d\t%d\t%d\t%d\n",j+1,(*molecule_angles_registry)[j][0],(*molecule_angles_registry)[j][1],(*molecule_angles_registry)[j][2],(*molecule_angles_registry)[j][3]);
    
    printf("\n$ New dihedrals topology:\n");
    printf("\n$ Dihedral types:\n");
    for(j=0;j<*molecule_dihedral_types;++j)printf("[%d]\t(%s)--%s--(%s)--%s--(%s)--%s--(%s)\n",j+1,(*molecule_dihedral_types_registry)[j][0],(*molecule_dihedral_types_registry)[j][1],(*molecule_dihedral_types_registry)[j][2],(*molecule_dihedral_types_registry)[j][3],(*molecule_dihedral_types_registry)[j][4],(*molecule_dihedral_types_registry)[j][5],(*molecule_dihedral_types_registry)[j][6]);
    printf("\n$ Dihedrals:\n");
    for(j=0;j<*molecule_dihedrals;++j)printf("[%d]\t%d\t%d\t%d\t%d\t%d\n",j+1,(*molecule_dihedrals_registry)[j][0],(*molecule_dihedrals_registry)[j][1],(*molecule_dihedrals_registry)[j][2],(*molecule_dihedrals_registry)[j][3],(*molecule_dihedrals_registry)[j][4]);
    
    printf("\n$ New impropers topology:\n");
    printf("\n$ Improper types:\n");
    for(j=0;j<*molecule_improper_types;++j)printf("[%d]\t%s\t%s\t%s\t%s\n",j+1,(*molecule_improper_types_registry)[j][0],(*molecule_improper_types_registry)[j][1],(*molecule_improper_types_registry)[j][2],(*molecule_improper_types_registry)[j][3]);
    printf("\n$ Impropers:\n");
    for(j=0;j<*molecule_impropers;++j)printf("[%d]\t%d\t%d\t%d\t%d\t%d\n",j+1,(*molecule_impropers_registry)[j][0],(*molecule_impropers_registry)[j][1],(*molecule_impropers_registry)[j][2],(*molecule_impropers_registry)[j][3],(*molecule_impropers_registry)[j][4]);
	}
    
    
   
    
    //getchar();
    
}
