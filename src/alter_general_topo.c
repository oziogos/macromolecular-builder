
#include"builder.h"

void alter_general_topo(int i, int molecules,
                        int molecule_atoms, int molecule_atoms_B, int delta_atoms, int molecule_atom_types,
                        int molecule_bonds, int molecule_bonds_B, int molecule_bond_types,
                        int molecule_angles, int molecule_angles_B, int molecule_angle_types,
                        int molecule_dihedrals, int molecule_dihedrals_B, int molecule_dihedral_types,
                        int molecule_impropers, int molecule_impropers_B, int molecule_improper_types,
                        int *atoms_per_molecule, int *atom_scaling_array,
                        char **molecule_species_registry,
                        char ***molecule_bond_types_registry, int **molecule_bonds_registry,
                        char ***molecule_angle_types_registry, int **molecule_angles_registry,
                        char ***molecule_dihedral_types_registry, int **molecule_dihedrals_registry,
                        char ***molecule_improper_types_registry, int **molecule_impropers_registry,
                        int *general_bonds, int *general_angles, int *general_dihedrals, int *general_impropers,
                        int **topo_boundaries,
                        //
                        int *stmap_rows, char ***stmap, int ***scm,
                        int *general_atom_types,
                        char ***general_species_registry,
                        
                        int *btmap_rows, char ****btmap, int ***bcm,
                        int *general_bond_types, int ***general_bonds_registry,
                        char ****general_bond_types_registry,
                        int ***btmapping_matrix,
                        
                        int *atmap_rows, char ****atmap, int ***acm,
                        int *general_angle_types, int ***general_angles_registry,
                        char ****general_angle_types_registry,
                        int ***atmapping_matrix,
                        
                        int *dtmap_rows, char ****dtmap, int ***dcm,
                        int *general_dihedral_types, int ***general_dihedrals_registry,
                        char ****general_dihedral_types_registry,
                        int ***dtmapping_matrix,
                        
                        int *itmap_rows, char ****itmap, int ***icm,
                        int *general_improper_types, int ***general_impropers_registry,
                        char ****general_improper_types_registry,
                        int ***itmapping_matrix,
                        
                        //
                        
                        int verb
                        
                        )
{
    int j,k,n,found,new_stypes,delete_flag;
    int *scm_sum,*bcm_sum,*acm_sum,*dcm_sum,*icm_sum;    
    int non_minus_above,non_minus_below,non_minus_current,delta,upper_0,lower_0,in_array[2],tcv,tav;
    int *gtmapping_array,counter;
    char ***backup;
    
    //
    if(verb==1){
	printf("\n$ scaling and counting arrays (unaltered):\n");
	for(j=0;j<molecules;++j)printf("[%d]\t%d\t%d\n",j+1,atoms_per_molecule[j],atom_scaling_array[j]);
	}
    // update current molecule's entry in atoms_per_molecule using delta_atoms
    atoms_per_molecule[i]=atoms_per_molecule[i]+delta_atoms;
    // repopulate the scaling array - the first entry is ommited because it is always zero!
    for(j=1;j<molecules;++j)
    {
        atom_scaling_array[j]=atom_scaling_array[j-1]+atoms_per_molecule[j-1];
    }

	//
	if(verb==1){
	printf("\n$ scaling and counting arrays:\n");
	for(j=0;j<molecules;++j)printf("[%d]\t%d\t%d\n",j+1,atoms_per_molecule[j],atom_scaling_array[j]);
	}
    //----
    
    // alter types
    // !!
    // !! implement tmap/cm for atomic species !!
    // !!
    
    // update stmap and scm
    for(j=0;j<molecule_atom_types;++j)
    {
        found=0;
        for(k=0;k<*stmap_rows;++k)
        {
            if(strcmp(molecule_species_registry[j],(*stmap)[k])==0)
            {
                found=1;break;
            }
        }
        if(found==0)
        {
            // augment
            *stmap_rows=*stmap_rows+1;
            // stmap
            // resize
            *stmap=(char**)realloc(*stmap,(*stmap_rows)*sizeof(char*));
            (*stmap)[*stmap_rows-1]=(char*)malloc(sub_length*sizeof(char));
            // write
            // write type
            sprintf((*stmap)[*stmap_rows-1],"%s",molecule_species_registry[j]);
            // scm
            for(k=0;k<molecules;++k)
            {
                (*scm)[k]=(int*)realloc((*scm)[k],(*stmap_rows)*sizeof(int*));
                (*scm)[k][*stmap_rows-1]=0;
            }
        }
    }
    
    // reset scm(i,:)
    for(j=0;j<*stmap_rows;++j)(*scm)[i][j]=0;
    // loop on mstr (molecule species type registry) using stmap and update scm
    for(j=0;j<molecule_atom_types;++j)
    {
        for(k=0;k<*stmap_rows;++k)
            if(strcmp(molecule_species_registry[j],(*stmap)[k])==0)
                break;
        (*scm)[i][k]=1;
    }
    
     // console out
     if(verb==1){
     printf("\n$ stmap/scm:\n");
     for(j=0;j<*stmap_rows;++j)printf("%s\t",(*stmap)[j]);printf("\n");
     for(j=0;j<molecules;++j){
     for(k=0;k<*stmap_rows;++k)
     printf("%d\t",(*scm)[j][k]);
     printf("\n");
	}
     }
     
    // global species through stmap/scm:
    scm_sum=(int*)malloc((*stmap_rows)*sizeof(int));
    n=0;
    for(k=0;k<*stmap_rows;++k)
    {
        scm_sum[k]=0;
        for(j=0;j<molecules;++j)
        {
            scm_sum[k]=scm_sum[k]+(*scm)[j][k];
        }
        //printf("%d\t",scm_sum[k]);
        if(scm_sum[k]!=0){n=n+1;}//printf("[%d]\t%s\n",n,(*stmap)[k]);}
    }
    new_stypes=n;
    if(new_stypes==*general_atom_types)
    {
        n=-1;
        for(k=0;k<*stmap_rows;++k)
            if(scm_sum[k]!=0){n=n+1;sprintf((*general_species_registry)[n],"%s",(*stmap)[k]);}
    }
    else
    {
        for(k=0;k<*general_atom_types;++k)free((*general_species_registry)[k]);free(*general_species_registry);
        *general_atom_types=new_stypes;
        *general_species_registry=(char**)malloc((*general_atom_types)*sizeof(char*));
        for(k=0;k<*general_atom_types;++k)(*general_species_registry)[k]=(char*)malloc(sub_length*sizeof(char));
        n=-1;
        for(k=0;k<*stmap_rows;++k)
            if(scm_sum[k]!=0){n=n+1;sprintf((*general_species_registry)[n],"%s",(*stmap)[k]);}
    }
    if(verb==1){
    printf("\n$ general_species_registry:\n");
    for(j=0;j<*general_atom_types;++j)printf("[%d]\t%s\n",j+1,(*general_species_registry)[j]);
	}
    free(scm_sum);
    
    /*
     ---------------------------------------------------------------------------
     what happens during a building step:
     
     the number of bonds can increase
     the number of bonds can remain the same, e.g. CH3-CH3 (7) --> CH3-COOH (7)
     the number of bonds can decrease, e.g. CH3-CH3 (7) --> CH3-CH=O (6)
     
     the number of angles can increase
     the number of angles can remain the same, e.g. CH4 (6) --> C2H4 (6)
     the number of angles can decrease, e.g. CH3-CH3 (12) --> CH3-COOH (10)
     
     BUT every molecule has always a finite number of bonds and angles, since 
     the smallest molecules supported are CH4, HOH, CH2=O
     
     this means that the alteration steps must be invoked for bond and angle
     registries
     
     let us consider the case of proper and improper dihedral angles:
     steps from [3] to [8] are to be applied only if the updated molecular
     topology contains dihedrals, e.i. molecule_dihedrals>0, molecule_impropers>0
     
     if this is not the case, a conditional id rescaling loop should suffice
     
     for memory compatibility reasons, only steps [7] and [8] are bypassed
     
     currently the bypass is applied only to impropers...
     
     ---------------------------------------------------------------------------
    */
    
    // alter bonds
    
    // [3]
    
    
    // update tmap and cm
    for(j=0;j<molecule_bond_types;++j)
    {
        found=0;
        // check if mbtr(j) belongs to btmap
        for(k=0;k<*btmap_rows;++k)
        {
            if(
               (strcmp(molecule_bond_types_registry[j][0],(*btmap)[k][0])==0 &&
                strcmp(molecule_bond_types_registry[j][1],(*btmap)[k][1])==0 &&
                strcmp(molecule_bond_types_registry[j][2],(*btmap)[k][2])==0
                )||
               (strcmp(molecule_bond_types_registry[j][0],(*btmap)[k][2])==0 &&
                strcmp(molecule_bond_types_registry[j][1],(*btmap)[k][1])==0 &&
                strcmp(molecule_bond_types_registry[j][2],(*btmap)[k][0])==0
                )
               )
            {
                found=1;break;
            }
        }
        // if you found it, do nothing...
        // if you didn't find it, add it to btmap and update bcm
        if(found==0)
        {
            // augment
            *btmap_rows=*btmap_rows+1;
            // tmap
            // resize
            *btmap=(char***)realloc(*btmap,(*btmap_rows)*sizeof(char**));
            (*btmap)[*btmap_rows-1]=(char**)malloc(3*sizeof(char*));
            (*btmap)[*btmap_rows-1][0]=(char*)malloc(sub_length*sizeof(char));
            (*btmap)[*btmap_rows-1][1]=(char*)malloc(sub_length*sizeof(char));
            (*btmap)[*btmap_rows-1][2]=(char*)malloc(sub_length*sizeof(char));
            // write
            // write type
            sprintf((*btmap)[*btmap_rows-1][0],"%s",molecule_bond_types_registry[j][0]);
            sprintf((*btmap)[*btmap_rows-1][1],"%s",molecule_bond_types_registry[j][1]);
            sprintf((*btmap)[*btmap_rows-1][2],"%s",molecule_bond_types_registry[j][2]);
            // cm
            for(k=0;k<molecules;++k)
            {
                (*bcm)[k]=(int*)realloc((*bcm)[k],(*btmap_rows)*sizeof(int*));
                (*bcm)[k][*btmap_rows-1]=0;
            }

        }

    }

    // [4]
    
    // reset cm(i,:)
    for(j=0;j<*btmap_rows;++j)(*bcm)[i][j]=0;
    // loop on mbtr (molecule bond type registry) using btmap and update bcm
    for(j=0;j<molecule_bond_types;++j)
    {
        for(k=0;k<*btmap_rows;++k)
            if((
                strcmp(molecule_bond_types_registry[j][0],(*btmap)[k][0])==0 &&
                strcmp(molecule_bond_types_registry[j][1],(*btmap)[k][1])==0 &&
                strcmp(molecule_bond_types_registry[j][2],(*btmap)[k][2])==0
                )
               ||
               (
                strcmp(molecule_bond_types_registry[j][0],(*btmap)[k][2])==0 &&
                strcmp(molecule_bond_types_registry[j][1],(*btmap)[k][1])==0 &&
                strcmp(molecule_bond_types_registry[j][2],(*btmap)[k][0])==0
                )
               )
                break;
        (*bcm)[i][k]=1;
    }
    
     // console out
     if(verb==1){
     printf("\n$ btmap/bcm:\n");
     for(j=0;j<*btmap_rows;++j)printf("(%s)--%s--(%s)\t",(*btmap)[j][0],(*btmap)[j][1],(*btmap)[j][2]);printf("\n");
     for(j=0;j<molecules;++j){
     for(k=0;k<*btmap_rows;++k)
     printf("%d\t",(*bcm)[j][k]);
     printf("\n");
	}
     }
    
    // [5]
    
    // global bond types through btmap/bcm:
    bcm_sum=(int*)malloc((*btmap_rows)*sizeof(int));
    for(k=0;k<*btmap_rows;++k)
    {
        bcm_sum[k]=0;
        for(j=0;j<molecules;++j)
        {
            bcm_sum[k]=bcm_sum[k]+(*bcm)[j][k];
        }
    }
    //
    if(verb==1){
    printf("\n$ bcm_sum array:\n");for(j=0;j<*btmap_rows;++j)printf("%d\t",bcm_sum[j]);printf("\n");
	}
	// [6]
	*btmapping_matrix=(int**)malloc(molecule_bond_types*sizeof(int*));
	for(j=0;j<molecule_bond_types;++j)(*btmapping_matrix)[j]=(int*)malloc(2*sizeof(int));
	
	n=-1;
	for(j=0;j<molecule_bond_types;++j)
	{
		found=0;
		for(k=0;k<*general_bond_types;++k)
		{
			if(
			   (strcmp(molecule_bond_types_registry[j][0],(*general_bond_types_registry)[k][0])==0 &&
				strcmp(molecule_bond_types_registry[j][1],(*general_bond_types_registry)[k][1])==0 &&
                strcmp(molecule_bond_types_registry[j][2],(*general_bond_types_registry)[k][2])==0
				)
			   ||
               (strcmp(molecule_bond_types_registry[j][0],(*general_bond_types_registry)[k][2])==0 &&
                strcmp(molecule_bond_types_registry[j][1],(*general_bond_types_registry)[k][1])==0 &&
                strcmp(molecule_bond_types_registry[j][2],(*general_bond_types_registry)[k][0])==0
                )
			   )
			{
				n=n+1;
				//printf("* molecular %d is found at general %d\n",j+1,k+1);
				(*btmapping_matrix)[n][0]=j+1;
				(*btmapping_matrix)[n][1]=k+1;
				found=1;
				break;
			}
		}
		if(found==0)
		{
			// augment
			*general_bond_types=*general_bond_types+1;
			// resize
			*general_bond_types_registry=(char***)realloc(*general_bond_types_registry,(*general_bond_types)*sizeof(char**));
			(*general_bond_types_registry)[*general_bond_types-1]=(char**)malloc(3*sizeof(char*));
			(*general_bond_types_registry)[*general_bond_types-1][0]=(char*)malloc(sub_length*sizeof(char));
			(*general_bond_types_registry)[*general_bond_types-1][1]=(char*)malloc(sub_length*sizeof(char));
            (*general_bond_types_registry)[*general_bond_types-1][2]=(char*)malloc(sub_length*sizeof(char));
            // write
			// write type
			sprintf((*general_bond_types_registry)[*general_bond_types-1][0],"%s",molecule_bond_types_registry[j][0]);
			sprintf((*general_bond_types_registry)[*general_bond_types-1][1],"%s",molecule_bond_types_registry[j][1]);
            sprintf((*general_bond_types_registry)[*general_bond_types-1][2],"%s",molecule_bond_types_registry[j][2]);
            //
			n=n+1;
			//printf("~ molecular %d is placed at general %d\n",j+1,*general_bond_types);
			(*btmapping_matrix)[n][0]=j+1;
			(*btmapping_matrix)[n][1]=*general_bond_types;
		}
		
	}
	if(verb==1){
	printf("\n$ btmapping_matrix:\n");
	for(j=0;j<molecule_bond_types;++j)printf("~ map %d --> %d\n",(*btmapping_matrix)[j][0],(*btmapping_matrix)[j][1]);

	printf("\n$ general bond types registry:\n");
	for(j=0;j<*general_bond_types;++j)printf("[%d]\t(%s)--%s--(%s)\n",j+1,(*general_bond_types_registry)[j][0],(*general_bond_types_registry)[j][1],(*general_bond_types_registry)[j][2]);
	}
	// [7]
	
	// initialize
	non_minus_above=0;
	non_minus_below=0;
	non_minus_current=0;
	// check current molecule i even column registry
	if(topo_boundaries[i][0]!=-1)non_minus_current=1;
	// traverse even column from 0 to i-1
	for(j=0;j<=i-1;++j)if(topo_boundaries[j][0]!=-1){non_minus_above=1;break;}
	// traverse even column from i+1 to molecules
	for(j=i+1;j<molecules;++j)if(topo_boundaries[j][0]!=-1){non_minus_below=1;break;}
	//
	if(verb==1){
	printf("\n$ minus flags:\n");
	printf("non_minus_above = %d\n",non_minus_above);
	printf("non_minus_below = %d\n",non_minus_below);
	printf("non_minus_current = %d\n\n",non_minus_current);
	}
	// case 1
	if(non_minus_above==0 && non_minus_below==0)
	{
        // for single chain
        if(non_minus_current==1)
        {
            // delta holds the number of bonds prior to the building step
            delta=topo_boundaries[i][1]-topo_boundaries[i][0]+1;
            // in_array holds the new values for the current molecule due to the building step
            in_array[0]=0;in_array[1]=molecule_bonds-1;
            if(verb==1)printf("$ case 1 in_array = [%d %d]\n",in_array[0],in_array[1]);
            // update current molecule values
            topo_boundaries[i][0]=in_array[0];
            topo_boundaries[i][1]=in_array[1];
            // tcv (topo current value) holds the difference in elements (delta_bonds)
            tcv=topo_boundaries[i][1]-topo_boundaries[i][0]+1-delta;
            if(verb==1)printf("$ case 1 tcv = %d\n",tcv);
            //
            if(tcv>0)
            {
                // resize adding extra space
                if(verb==1)printf("$ case 1 init / final sizes: %d / %d\n",*general_bonds,*general_bonds+tcv);
                *general_bonds_registry=(int**)realloc(*general_bonds_registry,(*general_bonds+tcv)*sizeof(int*));
                for(j=0;j<tcv;++j)(*general_bonds_registry)[*general_bonds+j]=(int*)malloc(3*sizeof(int));
                // place current molecule entries
                k=-1;
                for(j=topo_boundaries[i][0];j<=topo_boundaries[i][1];++j)
                {
                    k=k+1;
                    if(verb==1)printf("$ placing %d at %d...\n",k,j);
                    (*general_bonds_registry)[j][0]=(*btmapping_matrix)[molecule_bonds_registry[k][0]-1][1];    // type mapping
                    (*general_bonds_registry)[j][1]=molecule_bonds_registry[k][1]+atom_scaling_array[i];        // rescale
                    (*general_bonds_registry)[j][2]=molecule_bonds_registry[k][2]+atom_scaling_array[i];        // rescale
                }
            }
            else if(tcv<0)
            {
                // resize removing space
                if(verb==1){
                    printf("$ case 1 init / final sizes: %d / %d\n",*general_bonds,*general_bonds+tcv);
                    printf("\n$ tcv<0\n\n");
                    printf("\n$ topo_boundaries:\t%d\t%d\n",topo_boundaries[i][0],topo_boundaries[i][1]);}
                // place current molecule entries
                k=-1;
                for(j=topo_boundaries[i][0];j<=topo_boundaries[i][1];++j)
                {
                    k=k+1;
                    if(verb==1)printf("$ placing %d at %d...\n",k,j);
                    (*general_bonds_registry)[j][0]=(*btmapping_matrix)[molecule_bonds_registry[k][0]-1][1];    // type mapping
                    (*general_bonds_registry)[j][1]=molecule_bonds_registry[k][1]+atom_scaling_array[i];        // rescale
                    (*general_bonds_registry)[j][2]=molecule_bonds_registry[k][2]+atom_scaling_array[i];        // rescale
                }
                // free extra memory
                for(j=(*general_bonds-1);j>(*general_bonds+tcv-1);--j)
                {
                    free((*general_bonds_registry)[j]);
                }
            }
            else
            {
                // size kept the same
                if(verb==1)printf("$ case 1 init / final sizes: %d / %d\n",*general_bonds,*general_bonds+tcv);
                // place current molecule entries
                k=-1;
                for(j=topo_boundaries[i][0];j<=topo_boundaries[i][1];++j)
                {
                    k=k+1;
                    if(verb==1)printf("$ placing %d at %d...\n",k,j);
                    (*general_bonds_registry)[j][0]=(*btmapping_matrix)[molecule_bonds_registry[k][0]-1][1];    // type mapping
                    (*general_bonds_registry)[j][1]=molecule_bonds_registry[k][1]+atom_scaling_array[i];        // rescale
                    (*general_bonds_registry)[j][2]=molecule_bonds_registry[k][2]+atom_scaling_array[i];        // rescale
                }
            }
            // alter total
            *general_bonds=*general_bonds+tcv;
            //
            if(verb==1)printf("$ general bonds registry:\n");
            if(verb==1)for(j=0;j<*general_bonds;++j)printf("[%d]\t%d\t%d\t%d\n",j+1,(*general_bonds_registry)[j][0],(*general_bonds_registry)[j][1],(*general_bonds_registry)[j][2]);
        }
        else
        {
            printf("case 1 for bonds N/A!! How did you get in here??\n");exit(-1);
        }
    }
	
	// case 2
    // this case is triggered for bonds only when the current molecule is the *LAST* molecule!!!!
	else if(non_minus_above==1 && non_minus_below==0)
	{
        // delta holds the number of bonds prior to the building step
		delta=0;if(non_minus_current==1)delta=topo_boundaries[i][1]-topo_boundaries[i][0]+1;
		if(verb==1)printf("$ case 2 delta = %d\n",delta);
        // tav (topo above value) holds the largest eligible position value above the current molecule
		tav=topo_boundaries[0][1];for(j=0;j<=i-1;++j)if(topo_boundaries[j][1]>tav)tav=topo_boundaries[j][1];
		if(verb==1)printf("$ case 2 tav = %d\n",tav);
        // in_array holds the new values for the current molecule due to the building step
        in_array[0]=0+tav+1;in_array[1]=molecule_bonds-1+tav+1;
        if(verb==1)printf("$ case 2 in_array = [%d %d]\n",in_array[0],in_array[1]);
        // update current molecule values
		topo_boundaries[i][0]=in_array[0];
		topo_boundaries[i][1]=in_array[1];
        // tcv (topo current value) holds the difference in elements (delta_bonds)
		tcv=topo_boundaries[i][1]-topo_boundaries[i][0]+1-delta;
		if(verb==1){printf("$ case 2 tcv = %d\n",tcv);
        //
		printf("$ case 2 updated topo_boundaries:\n");
		for(j=0;j<molecules;++j)printf("[%d]\t%d\t%d\n",j+1,topo_boundaries[j][0],topo_boundaries[j][1]);}
        //
		if(tcv>0)
		{
            // resize adding extra space
			if(verb==1)printf("$ case 2 init / final sizes: %d / %d\n",*general_bonds,*general_bonds+tcv);
			*general_bonds_registry=(int**)realloc(*general_bonds_registry,(*general_bonds+tcv)*sizeof(int*));
			for(j=0;j<tcv;++j)(*general_bonds_registry)[*general_bonds+j]=(int*)malloc(3*sizeof(int));
		}
		else if (tcv<0)
		{
            // resize removing space
			printf("$ case 2 init / final sizes: %d / %d\n",*general_bonds,*general_bonds+tcv);
			printf("\n$ tcv<0\n\n");
			exit(-1);
		}
		else
		{
            // size kept the same
			if(verb==1)printf("$ case 2 init / final sizes: %d / %d\n",*general_bonds,*general_bonds+tcv);
		}
        // place current molecule entries
		k=-1;
		for(j=topo_boundaries[i][0];j<=topo_boundaries[i][1];++j)
		{
			k=k+1;
			if(verb==1)printf("$ placing %d at %d...\n",k,j);
			(*general_bonds_registry)[j][0]=(*btmapping_matrix)[molecule_bonds_registry[k][0]-1][1];    // utilize type mapping
			(*general_bonds_registry)[j][1]=molecule_bonds_registry[k][1]+atom_scaling_array[i];        // rescale ids
			(*general_bonds_registry)[j][2]=molecule_bonds_registry[k][2]+atom_scaling_array[i];        // rescale ids
		}
        // alter total
		*general_bonds=*general_bonds+tcv;
        //
		if(verb==1)printf("$ general bonds registry:\n");
		if(verb==1)for(j=0;j<*general_bonds;++j)printf("[%d]\t%d\t%d\t%d\n",j+1,(*general_bonds_registry)[j][0],(*general_bonds_registry)[j][1],(*general_bonds_registry)[j][2]);		
	}

	// case 3
    // this case is triggered for bonds only when the current molecule is the *FIRST* molecule!!!!
	else if(non_minus_above==0 && non_minus_below==1)
	{
        // delta holds the number of bonds prior to the building step
		delta=0;if(non_minus_current==1)delta=topo_boundaries[i][1]-topo_boundaries[i][0]+1;
		if(verb==1)printf("$ case 3 delta = %d\n",delta);
        // upper_0 and lower_0 define the boundaries for the values below
		for(j=i+1;j<molecules;++j)if(topo_boundaries[j][0]!=-1){upper_0=topo_boundaries[j][0];break;}
		lower_0=topo_boundaries[i][1];
		for(j=i+1;j<molecules;++j)if(topo_boundaries[j][1]>lower_0){lower_0=topo_boundaries[j][1];}
		if(verb==1){printf("$ case 3 upper_0 = %d\n",upper_0);
		printf("$ case 3 lower_0 = %d\n",lower_0);}
        // in_array holds the new values for the current molecule due to the building step
		in_array[0]=0;in_array[1]=molecule_bonds-1;
		if(verb==1)printf("$ case 3 in_array = [%d %d]\n",in_array[0],in_array[1]);
        // update current molecule values
		topo_boundaries[i][0]=in_array[0];
		topo_boundaries[i][1]=in_array[1];
        // tcv (topo current value) holds the difference in elements (delta_bonds)
		tcv=topo_boundaries[i][1]-topo_boundaries[i][0]+1-delta;
		if(verb==1)printf("$ case 3 tcv = %d\n",tcv);
        // update non -1 values below current molecule
		for(j=i+1;j<molecules;++j){if(topo_boundaries[j][0]!=-1){topo_boundaries[j][0]=topo_boundaries[j][0]+tcv;topo_boundaries[j][1]=topo_boundaries[j][1]+tcv;}}
		//
        if(verb==1){printf("$ case 3 updated topo_boundaries:\n");
		for(j=0;j<molecules;++j)printf("[%d]\t%d\t%d\n",j+1,topo_boundaries[j][0],topo_boundaries[j][1]);}
		//
        if(tcv>0)
		{
            // resize adding extra space
			if(verb==1)printf("$ case 3 init / final sizes: %d / %d\n",*general_bonds,*general_bonds+tcv);
			*general_bonds_registry=(int**)realloc(*general_bonds_registry,(*general_bonds+tcv)*sizeof(int*));
			for(j=0;j<tcv;++j)(*general_bonds_registry)[*general_bonds+j]=(int*)malloc(3*sizeof(int));
            // push down entries
			if(verb==1)printf("$ case 3 push down: %d:%d --> %d:%d\n",upper_0,lower_0,upper_0+tcv,lower_0+tcv);
			for(j=lower_0;j>=upper_0;--j)
			{
				if(verb==1)printf("%d --> %d\n",j,j+tcv);
				(*general_bonds_registry)[j+tcv][0]=(*general_bonds_registry)[j][0];
				(*general_bonds_registry)[j+tcv][1]=(*general_bonds_registry)[j][1]+delta_atoms;    // rescale
				(*general_bonds_registry)[j+tcv][2]=(*general_bonds_registry)[j][2]+delta_atoms;    // rescale
			}
		}
		else if (tcv<0)
		{
            // resize removing space
			printf("$ case 3 init / final sizes: %d / %d\n",*general_bonds,*general_bonds+tcv);
			printf("\n$ tcv<0\n\n");
			exit(-1);
		}
		else
		{
            // size kept the same
			if(verb==1)printf("$ case 3 init / final sizes: %d / %d\n",*general_bonds,*general_bonds+tcv);
            // update for rescaling purposes
			if(verb==1)printf("$ case 3 no shift: %d:%d <--> %d:%d\n",upper_0,lower_0,upper_0+tcv,lower_0+tcv);
			for(j=lower_0;j>=upper_0;--j)
			{
				if(verb==1)printf("%d --> %d\n",j,j+tcv);
				(*general_bonds_registry)[j+tcv][0]=(*general_bonds_registry)[j][0];
				(*general_bonds_registry)[j+tcv][1]=(*general_bonds_registry)[j][1]+delta_atoms;    // rescale
				(*general_bonds_registry)[j+tcv][2]=(*general_bonds_registry)[j][2]+delta_atoms;    // rescale
			}
		}
        // place current molecule entries
		k=-1;
		for(j=topo_boundaries[i][0];j<=topo_boundaries[i][1];++j)
		{
			k=k+1;
			if(verb==1)printf("$ placing %d at %d...\n",k,j);
			(*general_bonds_registry)[j][0]=(*btmapping_matrix)[molecule_bonds_registry[k][0]-1][1];    // type mapping
			(*general_bonds_registry)[j][1]=molecule_bonds_registry[k][1]+atom_scaling_array[i];        // rescale
			(*general_bonds_registry)[j][2]=molecule_bonds_registry[k][2]+atom_scaling_array[i];        // rescale
		}
        // alter total
		*general_bonds=*general_bonds+tcv;
        //
		if(verb==1)printf("$ general bonds registry:\n");
		if(verb==1)for(j=0;j<*general_bonds;++j)printf("[%d]\t%d\t%d\t%d\n",j+1,(*general_bonds_registry)[j][0],(*general_bonds_registry)[j][1],(*general_bonds_registry)[j][2]);
	}
	
	// case 4
    // active for all intermediate cases...
	else if(non_minus_above==1 && non_minus_below==1)
	{
        // delta holds the number of bonds prior to the building step
		delta=0;if(non_minus_current==1)delta=topo_boundaries[i][1]-topo_boundaries[i][0]+1;
		if(verb==1)printf("$ case 4 delta = %d\n",delta);
        // upper_0 and lower_0 define the boundaries for the values below
		for(j=i+1;j<molecules;++j)if(topo_boundaries[j][0]!=-1){upper_0=topo_boundaries[j][0];break;}
		lower_0=topo_boundaries[i][1];
		for(j=i+1;j<molecules;++j)if(topo_boundaries[j][1]>lower_0){lower_0=topo_boundaries[j][1];}
		if(verb==1){printf("$ case 4 upper_0 = %d\n",upper_0);
		printf("$ case 4 lower_0 = %d\n",lower_0);}
        // tav (topo above value) holds the largest eligible position value above the current molecule
		tav=topo_boundaries[0][1];for(j=0;j<=i-1;++j)if(topo_boundaries[j][1]>tav)tav=topo_boundaries[j][1];
		if(verb==1)printf("$ case 4 tav = %d\n",tav);
        // in_array holds the new values for the current molecule due to the building step
        in_array[0]=0+tav+1;in_array[1]=molecule_bonds-1+tav+1;
		if(verb==1)printf("$ case 4 in_array = [%d %d]\n",in_array[0],in_array[1]);
        // update current molecule values
		topo_boundaries[i][0]=in_array[0];
		topo_boundaries[i][1]=in_array[1];
        // tcv (topo current value) holds the difference in elements (delta_bonds)
		tcv=topo_boundaries[i][1]-topo_boundaries[i][0]+1-delta;
		if(verb==1)printf("$ case 4 tcv = %d\n",tcv);
        // update non -1 values below current molecule
		for(j=i+1;j<molecules;++j){if(topo_boundaries[j][0]!=-1){topo_boundaries[j][0]=topo_boundaries[j][0]+tcv;topo_boundaries[j][1]=topo_boundaries[j][1]+tcv;}}
        //
		if(verb==1){printf("$ case 4 updated topo_boundaries:\n");
		for(j=0;j<molecules;++j)printf("[%d]\t%d\t%d\n",j+1,topo_boundaries[j][0],topo_boundaries[j][1]);}
        //
		if(tcv>0)
		{
            // resize adding extra space
			if(verb==1)printf("$ case 4 init / final sizes: %d / %d\n",*general_bonds,*general_bonds+tcv);
			*general_bonds_registry=(int**)realloc(*general_bonds_registry,(*general_bonds+tcv)*sizeof(int*));
			for(j=0;j<tcv;++j)(*general_bonds_registry)[*general_bonds+j]=(int*)malloc(3*sizeof(int));
            // push down entries
			if(verb==1)printf("$ case 4 push down: %d:%d --> %d:%d\n",upper_0,lower_0,upper_0+tcv,lower_0+tcv);
			for(j=lower_0;j>=upper_0;--j)
			{
				if(verb==1)printf("%d --> %d\n",j,j+tcv);
				(*general_bonds_registry)[j+tcv][0]=(*general_bonds_registry)[j][0];
				(*general_bonds_registry)[j+tcv][1]=(*general_bonds_registry)[j][1]+delta_atoms;    // rescale
				(*general_bonds_registry)[j+tcv][2]=(*general_bonds_registry)[j][2]+delta_atoms;    // rescale
			}
		}
		else if (tcv<0)
		{
            // resize removing space
			printf("$ case 4 init / final sizes: %d / %d\n",*general_bonds,*general_bonds+tcv);
			printf("\n$ tcv<0\n\n");
			exit(-1);
		}
		else
		{
            // size kept the same
			if(verb==1){printf("$ case 4 init / final sizes: %d / %d\n",*general_bonds,*general_bonds+tcv);
			printf("$ case 4 no shift: %d:%d <--> %d:%d\n",upper_0,lower_0,upper_0+tcv,lower_0+tcv);}
            // update for rescaling purposes
			for(j=lower_0;j>=upper_0;--j)
			{
				if(verb==1)printf("%d --> %d\n",j,j+tcv);
				(*general_bonds_registry)[j+tcv][0]=(*general_bonds_registry)[j][0];
				(*general_bonds_registry)[j+tcv][1]=(*general_bonds_registry)[j][1]+delta_atoms;    // rescale
				(*general_bonds_registry)[j+tcv][2]=(*general_bonds_registry)[j][2]+delta_atoms;    // rescale
			}
		}
        // place current molecule entries
		k=-1;
		for(j=topo_boundaries[i][0];j<=topo_boundaries[i][1];++j)
		{
			k=k+1;
			if(verb==1)printf("$ placing %d at %d...\n",k,j);
			(*general_bonds_registry)[j][0]=(*btmapping_matrix)[molecule_bonds_registry[k][0]-1][1];    // type mapping
			(*general_bonds_registry)[j][1]=molecule_bonds_registry[k][1]+atom_scaling_array[i];        // rescale
			(*general_bonds_registry)[j][2]=molecule_bonds_registry[k][2]+atom_scaling_array[i];        // rescale
		}
        // alter total
		*general_bonds=*general_bonds+tcv;
        //
		if(verb==1)printf("$ general bonds registry:\n");
		if(verb==1)for(j=0;j<*general_bonds;++j)printf("[%d]\t%d\t%d\t%d\n",j+1,(*general_bonds_registry)[j][0],(*general_bonds_registry)[j][1],(*general_bonds_registry)[j][2]);
	}
	else
	{
		printf("unsupported bonds combination!!\n\n");exit(-1);
	}

	delete_flag=0;
    for(k=0;k<*btmap_rows;++k)
        if(bcm_sum[k]==0)
            delete_flag=delete_flag+1;
    
    if(delete_flag>0)
    {
		if(verb==1)printf("\n$ delete flag>0!!!\n\n");
        gtmapping_array=(int*)malloc((*general_bond_types)*sizeof(int));
        backup=(char***)malloc((*general_bond_types)*sizeof(char**));
        for(j=0;j<*general_bond_types;++j)backup[j]=(char**)malloc(3*sizeof(char*));
        for(j=0;j<*general_bond_types;++j)for(k=0;k<3;++k)backup[j][k]=(char*)malloc(sub_length*sizeof(char));
        counter=0;
        for(j=0;j<*general_bond_types;++j)
        {
            for(k=0;k<*btmap_rows;++k)
            {
                if((strcmp((*general_bond_types_registry)[j][0],(*btmap)[k][0])==0 &&
                    strcmp((*general_bond_types_registry)[j][1],(*btmap)[k][1])==0 &&
                    strcmp((*general_bond_types_registry)[j][2],(*btmap)[k][2])==0)
                   ||
                   (strcmp((*general_bond_types_registry)[j][0],(*btmap)[k][2])==0 &&
                    strcmp((*general_bond_types_registry)[j][1],(*btmap)[k][1])==0 &&
                    strcmp((*general_bond_types_registry)[j][2],(*btmap)[k][0])==0)
                   )
                    break;
            }
            if(verb==1){printf("$ (%s)--%s--(%s) is at position {%d}\n",(*general_bond_types_registry)[j][0],(*general_bond_types_registry)[j][1],(*general_bond_types_registry)[j][2],k+1);
            printf("$ checking cm_sum: %d\n",bcm_sum[k]);}
            if(bcm_sum[k]==0)
            {
                gtmapping_array[j]=0;
            }
            else
            {
                counter=counter+1;
                gtmapping_array[j]=counter;
                sprintf(backup[counter-1][0],"%s",(*general_bond_types_registry)[j][0]);
                sprintf(backup[counter-1][1],"%s",(*general_bond_types_registry)[j][1]);
                sprintf(backup[counter-1][2],"%s",(*general_bond_types_registry)[j][2]);
            }
            
        }
        if(verb==1){
        printf("$ gtmapping_array:\n");
        for(j=0;j<*general_bond_types;++j)printf("[%d]\t%d\n",j+1,gtmapping_array[j]);
        printf("$ altered types:\n");
        for(j=0;j<counter;++j)printf("[%d]\t(%s)--%s--(%s)\n",j+1,backup[j][0],backup[j][1],backup[j][2]);}
        
        for(j=0;j<*general_bonds;++j)(*general_bonds_registry)[j][0]=gtmapping_array[(*general_bonds_registry)[j][0]-1];
        for(j=0;j<counter;++j)
        {
            sprintf((*general_bond_types_registry)[j][0],"%s",backup[j][0]);
            sprintf((*general_bond_types_registry)[j][1],"%s",backup[j][1]);
            sprintf((*general_bond_types_registry)[j][2],"%s",backup[j][2]);
        }
        for(j=counter;j<*general_bond_types;++j)
        {
            for(k=0;k<3;++k)
                free((*general_bond_types_registry)[j][k]);
            free((*general_bond_types_registry)[j]);
        }
        for(j=0;j<*general_bond_types;++j)
            for(k=0;k<3;++k)
                free(backup[j][k]);
        for(j=0;j<*general_bond_types;++j)free(backup[j]);
        free(backup);
        *general_bond_types=counter;
        free(gtmapping_array);
        
        //
		if(verb==1)printf("$ general bonds registry:\n");
		if(verb==1)for(j=0;j<*general_bonds;++j)printf("[%d]\t%d\t%d\t%d\n",j+1,(*general_bonds_registry)[j][0],(*general_bonds_registry)[j][1],(*general_bonds_registry)[j][2]);
		        
	}

    free(bcm_sum);
    
    //--------------------------------------------------------------------------
    
    // alter angles
    
    // [3]
    
    // update tmap and cm
    for(j=0;j<molecule_angle_types;++j)
    {
        found=0;
        for(k=0;k<*atmap_rows;++k)
        {
            if(
               (strcmp(molecule_angle_types_registry[j][0],(*atmap)[k][0])==0 &&
                strcmp(molecule_angle_types_registry[j][1],(*atmap)[k][1])==0 &&
                strcmp(molecule_angle_types_registry[j][2],(*atmap)[k][2])==0 &&
                strcmp(molecule_angle_types_registry[j][3],(*atmap)[k][3])==0 &&
                strcmp(molecule_angle_types_registry[j][4],(*atmap)[k][4])==0
                )||
               (strcmp(molecule_angle_types_registry[j][0],(*atmap)[k][4])==0 &&
                strcmp(molecule_angle_types_registry[j][1],(*atmap)[k][3])==0 &&
                strcmp(molecule_angle_types_registry[j][2],(*atmap)[k][2])==0 &&
                strcmp(molecule_angle_types_registry[j][3],(*atmap)[k][1])==0 &&
                strcmp(molecule_angle_types_registry[j][4],(*atmap)[k][0])==0
                )
               )
            {
                found=1;break;
            }
        }
        if(found==0)
        {
            // augment
            *atmap_rows=*atmap_rows+1;
            // tmap
            // resize
            *atmap=(char***)realloc(*atmap,(*atmap_rows)*sizeof(char**));
            (*atmap)[*atmap_rows-1]=(char**)malloc(5*sizeof(char*));
            (*atmap)[*atmap_rows-1][0]=(char*)malloc(sub_length*sizeof(char));
            (*atmap)[*atmap_rows-1][1]=(char*)malloc(sub_length*sizeof(char));
            (*atmap)[*atmap_rows-1][2]=(char*)malloc(sub_length*sizeof(char));
            (*atmap)[*atmap_rows-1][3]=(char*)malloc(sub_length*sizeof(char));
            (*atmap)[*atmap_rows-1][4]=(char*)malloc(sub_length*sizeof(char));
            // write
            // write type
            sprintf((*atmap)[*atmap_rows-1][0],"%s",molecule_angle_types_registry[j][0]);
            sprintf((*atmap)[*atmap_rows-1][1],"%s",molecule_angle_types_registry[j][1]);
            sprintf((*atmap)[*atmap_rows-1][2],"%s",molecule_angle_types_registry[j][2]);
            sprintf((*atmap)[*atmap_rows-1][3],"%s",molecule_angle_types_registry[j][3]);
            sprintf((*atmap)[*atmap_rows-1][4],"%s",molecule_angle_types_registry[j][4]);
            // cm
            for(k=0;k<molecules;++k)
            {
                (*acm)[k]=(int*)realloc((*acm)[k],(*atmap_rows)*sizeof(int*));
                (*acm)[k][*atmap_rows-1]=0;
            }
        }
    }
    
    // [4]
    
    // reset cm(i,:)
    for(j=0;j<*atmap_rows;++j)(*acm)[i][j]=0;
    // loop on matr (molecule angle type registry) using atmap and update acm
    for(j=0;j<molecule_angle_types;++j)
    {
        for(k=0;k<*atmap_rows;++k)
            if((
                strcmp(molecule_angle_types_registry[j][0],(*atmap)[k][0])==0 &&
                strcmp(molecule_angle_types_registry[j][1],(*atmap)[k][1])==0 &&
                strcmp(molecule_angle_types_registry[j][2],(*atmap)[k][2])==0 &&
                strcmp(molecule_angle_types_registry[j][3],(*atmap)[k][3])==0 &&
                strcmp(molecule_angle_types_registry[j][4],(*atmap)[k][4])==0
                )
               ||
               (
                strcmp(molecule_angle_types_registry[j][0],(*atmap)[k][4])==0 &&
                strcmp(molecule_angle_types_registry[j][1],(*atmap)[k][3])==0 &&
                strcmp(molecule_angle_types_registry[j][2],(*atmap)[k][2])==0 &&
                strcmp(molecule_angle_types_registry[j][3],(*atmap)[k][1])==0 &&
                strcmp(molecule_angle_types_registry[j][4],(*atmap)[k][0])==0
                )
               )
                break;
        (*acm)[i][k]=1;
    }
    
     // console out
     if(verb==1){
     printf("\n$ atmap/acm:\n");
     for(j=0;j<*atmap_rows;++j)printf("(%s)--%s--(%s)--%s--(%s)\t",(*atmap)[j][0],(*atmap)[j][1],(*atmap)[j][2],(*atmap)[j][3],(*atmap)[j][4]);printf("\n");
     for(j=0;j<molecules;++j){
     for(k=0;k<*atmap_rows;++k)
     printf("%d\t",(*acm)[j][k]);
     printf("\n");
     }
    }
    // [5]
     
    // global angle types through atmap/acm:
    acm_sum=(int*)malloc((*atmap_rows)*sizeof(int));
    for(k=0;k<*atmap_rows;++k)
    {
        acm_sum[k]=0;
        for(j=0;j<molecules;++j)
        {
            acm_sum[k]=acm_sum[k]+(*acm)[j][k];
        }
    }
    //
    if(verb==1){
    printf("\n$ acm_sum array:\n");for(j=0;j<*atmap_rows;++j)printf("%d\t",acm_sum[j]);printf("\n");}
    
    // [6]
	*atmapping_matrix=(int**)malloc(molecule_angle_types*sizeof(int*));
	for(j=0;j<molecule_angle_types;++j)(*atmapping_matrix)[j]=(int*)malloc(2*sizeof(int));
	
	n=-1;
	for(j=0;j<molecule_angle_types;++j)
	{
		found=0;
		for(k=0;k<*general_angle_types;++k)
		{
			if(
			   (strcmp(molecule_angle_types_registry[j][0],(*general_angle_types_registry)[k][0])==0 &&
				strcmp(molecule_angle_types_registry[j][1],(*general_angle_types_registry)[k][1])==0 &&
				strcmp(molecule_angle_types_registry[j][2],(*general_angle_types_registry)[k][2])==0 &&
                strcmp(molecule_angle_types_registry[j][3],(*general_angle_types_registry)[k][3])==0 &&
                strcmp(molecule_angle_types_registry[j][4],(*general_angle_types_registry)[k][4])==0
				)
			   ||
               (strcmp(molecule_angle_types_registry[j][0],(*general_angle_types_registry)[k][4])==0 &&
                strcmp(molecule_angle_types_registry[j][1],(*general_angle_types_registry)[k][3])==0 &&
                strcmp(molecule_angle_types_registry[j][2],(*general_angle_types_registry)[k][2])==0 &&
                strcmp(molecule_angle_types_registry[j][3],(*general_angle_types_registry)[k][1])==0 &&
                strcmp(molecule_angle_types_registry[j][4],(*general_angle_types_registry)[k][0])==0
                )
			   )
			{
				n=n+1;
				//printf("* molecular %d is found at general %d\n",j+1,k+1);
				(*atmapping_matrix)[n][0]=j+1;
				(*atmapping_matrix)[n][1]=k+1;
				found=1;
				break;
			}
		}
		if(found==0)
		{
			// augment
			*general_angle_types=*general_angle_types+1;
			// resize
			*general_angle_types_registry=(char***)realloc(*general_angle_types_registry,(*general_angle_types)*sizeof(char**));
			(*general_angle_types_registry)[*general_angle_types-1]=(char**)malloc(5*sizeof(char*));
			(*general_angle_types_registry)[*general_angle_types-1][0]=(char*)malloc(sub_length*sizeof(char));
			(*general_angle_types_registry)[*general_angle_types-1][1]=(char*)malloc(sub_length*sizeof(char));
			(*general_angle_types_registry)[*general_angle_types-1][2]=(char*)malloc(sub_length*sizeof(char));
            (*general_angle_types_registry)[*general_angle_types-1][3]=(char*)malloc(sub_length*sizeof(char));
            (*general_angle_types_registry)[*general_angle_types-1][4]=(char*)malloc(sub_length*sizeof(char));
			// write
			// write type
			sprintf((*general_angle_types_registry)[*general_angle_types-1][0],"%s",molecule_angle_types_registry[j][0]);
			sprintf((*general_angle_types_registry)[*general_angle_types-1][1],"%s",molecule_angle_types_registry[j][1]);
			sprintf((*general_angle_types_registry)[*general_angle_types-1][2],"%s",molecule_angle_types_registry[j][2]);
            sprintf((*general_angle_types_registry)[*general_angle_types-1][3],"%s",molecule_angle_types_registry[j][3]);
            sprintf((*general_angle_types_registry)[*general_angle_types-1][4],"%s",molecule_angle_types_registry[j][4]);
			//
			n=n+1;
			//printf("~ molecular %d is placed at general %d\n",j+1,*general_angle_types);
			(*atmapping_matrix)[n][0]=j+1;
			(*atmapping_matrix)[n][1]=*general_angle_types;
		}
		
	}
	if(verb==1){
	printf("\n$ atmapping_matrix:\n");
	for(j=0;j<molecule_angle_types;++j)printf("~ map %d --> %d\n",(*atmapping_matrix)[j][0],(*atmapping_matrix)[j][1]);

	printf("\n$ general angle types registry:\n");
	for(j=0;j<*general_angle_types;++j)printf("[%d]\t(%s)--%s--(%s)--%s--(%s)\n",j+1,(*general_angle_types_registry)[j][0],(*general_angle_types_registry)[j][1],(*general_angle_types_registry)[j][2],(*general_angle_types_registry)[j][3],(*general_angle_types_registry)[j][4]);
	}
	// [7]
	
	// initialize
	non_minus_above=0;
	non_minus_below=0;
	non_minus_current=0;
	// check current molecule i even column registry
	if(topo_boundaries[i][2]!=-1)non_minus_current=1;
	// traverse even column from 0 to i-1
	for(j=0;j<=i-1;++j)if(topo_boundaries[j][2]!=-1){non_minus_above=1;break;}
	// traverse even column from i+1 to molecules
	for(j=i+1;j<molecules;++j)if(topo_boundaries[j][2]!=-1){non_minus_below=1;break;}
	//
	if(verb==1){
	printf("\n$ minus flags:\n");
	printf("non_minus_above = %d\n",non_minus_above);
	printf("non_minus_below = %d\n",non_minus_below);
	printf("non_minus_current = %d\n\n",non_minus_current);}

	// case 1
	if(non_minus_above==0 && non_minus_below==0)
	{
        // for single chain
        if(non_minus_current==1)
        {
            // delta holds the number of angles prior to the building step
            delta=topo_boundaries[i][3]-topo_boundaries[i][2]+1;
            // in_array holds the new values for the current molecule due to the building step
            in_array[0]=0;in_array[1]=molecule_angles-1;
            if(verb==1)printf("$ case 1 in_array = [%d %d]\n",in_array[0],in_array[1]);
            // update current molecule values
            topo_boundaries[i][2]=in_array[0];
            topo_boundaries[i][3]=in_array[1];
            // tcv (topo current value) holds the difference in elements (delta_angles)
            tcv=topo_boundaries[i][3]-topo_boundaries[i][2]+1-delta;
            if(verb==1)printf("$ case 1 tcv = %d\n",tcv);
            //
            if(tcv>0)
            {
                // resize adding extra space
                if(verb==1)printf("$ case 1 init / final sizes: %d / %d\n",*general_angles,*general_angles+tcv);
                *general_angles_registry=(int**)realloc(*general_angles_registry,(*general_angles+tcv)*sizeof(int*));
                for(j=0;j<tcv;++j)(*general_angles_registry)[*general_angles+j]=(int*)malloc(4*sizeof(int));
                // place current molecule entries
                k=-1;
                for(j=topo_boundaries[i][2];j<=topo_boundaries[i][3];++j)
                {
                    k=k+1;
                    if(verb==1)printf("$ placing %d at %d...\n",k,j);
                    (*general_angles_registry)[j][0]=(*atmapping_matrix)[molecule_angles_registry[k][0]-1][1];    // type mapping
                    (*general_angles_registry)[j][1]=molecule_angles_registry[k][1]+atom_scaling_array[i];        // rescale
                    (*general_angles_registry)[j][2]=molecule_angles_registry[k][2]+atom_scaling_array[i];        // rescale
                    (*general_angles_registry)[j][3]=molecule_angles_registry[k][3]+atom_scaling_array[i];        // rescale
                }
            }
            else if(tcv<0)
            {
                // resize removing space
                if(verb==1){
                    printf("$ case 1 init / final sizes: %d / %d\n",*general_angles,*general_angles+tcv);
                    printf("\n$ tcv<0\n\n");
                    printf("\n$ topo_boundaries:\t%d\t%d\n",topo_boundaries[i][2],topo_boundaries[i][3]);}
                // place current molecule entries
                k=-1;
                for(j=topo_boundaries[i][2];j<=topo_boundaries[i][3];++j)
                {
                    k=k+1;
                    if(verb==1)printf("$ placing %d at %d...\n",k,j);
                    (*general_angles_registry)[j][0]=(*atmapping_matrix)[molecule_angles_registry[k][0]-1][1];    // type mapping
                    (*general_angles_registry)[j][1]=molecule_angles_registry[k][1]+atom_scaling_array[i];        // rescale
                    (*general_angles_registry)[j][2]=molecule_angles_registry[k][2]+atom_scaling_array[i];        // rescale
                    (*general_angles_registry)[j][3]=molecule_angles_registry[k][3]+atom_scaling_array[i];        // rescale
                }
                // free extra memory
                for(j=(*general_angles-1);j>(*general_angles+tcv-1);--j)
                {
                    free((*general_angles_registry)[j]);
                }
            }
            else
            {
                // size kept the same
                if(verb==1)printf("$ case 1 init / final sizes: %d / %d\n",*general_angles,*general_angles+tcv);
                // place current molecule entries
                k=-1;
                for(j=topo_boundaries[i][2];j<=topo_boundaries[i][3];++j)
                {
                    k=k+1;
                    if(verb==1)printf("$ placing %d at %d...\n",k,j);
                    (*general_angles_registry)[j][0]=(*btmapping_matrix)[molecule_angles_registry[k][0]-1][1];    // type mapping
                    (*general_angles_registry)[j][1]=molecule_angles_registry[k][1]+atom_scaling_array[i];        // rescale
                    (*general_angles_registry)[j][2]=molecule_angles_registry[k][2]+atom_scaling_array[i];        // rescale
                    (*general_angles_registry)[j][3]=molecule_angles_registry[k][3]+atom_scaling_array[i];        // rescale
                }
            }
            // alter total
            *general_angles=*general_angles+tcv;
            //
            if(verb==1)printf("$ general angles registry:\n");
            if(verb==1)for(j=0;j<*general_angles;++j)printf("[%d]\t%d\t%d\t%d\t%d\n",j+1,(*general_angles_registry)[j][0],(*general_angles_registry)[j][1],(*general_angles_registry)[j][2],(*general_angles_registry)[j][3]);
        }
        else
        {
            printf("case 1 for angles N/A!! How did you get in here??\n");exit(-1);
        }
        
        //getchar();
        
	}

	// case 2
    // this case is triggered for angles only when the current molecule is the *LAST* molecule!!!!
	else if(non_minus_above==1 && non_minus_below==0)
	{
        // delta holds the number of angles prior to the building step
		delta=0;if(non_minus_current==1)delta=topo_boundaries[i][3]-topo_boundaries[i][2]+1;
		if(verb==1)printf("$ case 2 delta = %d\n",delta);
        // tav (topo above value) holds the largest eligible position value above the current molecule
		tav=topo_boundaries[0][3];for(j=0;j<=i-1;++j)if(topo_boundaries[j][3]>tav)tav=topo_boundaries[j][3];
		if(verb==1)printf("$ case 2 tav = %d\n",tav);
        // in_array holds the new values for the current molecule due to the building step
        in_array[0]=0+tav+1;in_array[1]=molecule_angles-1+tav+1;
        if(verb==1)printf("$ case 2 in_array = [%d %d]\n",in_array[0],in_array[1]);
        // update current molecule values
		topo_boundaries[i][2]=in_array[0];
		topo_boundaries[i][3]=in_array[1];
        // tcv (topo current value) holds the difference in elements (delta_angles)
		tcv=topo_boundaries[i][3]-topo_boundaries[i][2]+1-delta;
		if(verb==1){printf("$ case 2 tcv = %d\n",tcv);
		printf("$ case 2 updated topo_boundaries:\n");
		for(j=0;j<molecules;++j)printf("[%d]\t%d\t%d\n",j+1,topo_boundaries[j][2],topo_boundaries[j][3]);}
        //
		if(tcv>0)
		{
            // resize adding extra space
			if(verb==1)printf("$ case 2 init / final sizes: %d / %d\n",*general_angles,*general_angles+tcv);
			*general_angles_registry=(int**)realloc(*general_angles_registry,(*general_angles+tcv)*sizeof(int*));
			for(j=0;j<tcv;++j)(*general_angles_registry)[*general_angles+j]=(int*)malloc(4*sizeof(int));
			// place current molecule entries
			k=-1;
			for(j=topo_boundaries[i][2];j<=topo_boundaries[i][3];++j)
			{
				k=k+1;
				if(verb==1)printf("$ placing %d at %d...\n",k,j);
				(*general_angles_registry)[j][0]=(*atmapping_matrix)[molecule_angles_registry[k][0]-1][1];    // type mapping
				(*general_angles_registry)[j][1]=molecule_angles_registry[k][1]+atom_scaling_array[i];        // rescale
				(*general_angles_registry)[j][2]=molecule_angles_registry[k][2]+atom_scaling_array[i];        // rescale
				(*general_angles_registry)[j][3]=molecule_angles_registry[k][3]+atom_scaling_array[i];        // rescale
			}
		}
		else if (tcv<0)
		{
            // resize removing space
			if(verb==1){printf("$ case 2 init / final sizes: %d / %d\n",*general_angles,*general_angles+tcv);
			printf("\n$ tcv<0\n\n");
			printf("\n$ topo_boundaries:\t%d\t%d\n",topo_boundaries[i][2],topo_boundaries[i][3]);}
			// place current molecule entries
			k=-1;
			for(j=topo_boundaries[i][2];j<=topo_boundaries[i][3];++j)
			{
				k=k+1;
				if(verb==1)printf("$ placing %d at %d...\n",k,j);
				(*general_angles_registry)[j][0]=(*atmapping_matrix)[molecule_angles_registry[k][0]-1][1];    // type mapping
				(*general_angles_registry)[j][1]=molecule_angles_registry[k][1]+atom_scaling_array[i];        // rescale
				(*general_angles_registry)[j][2]=molecule_angles_registry[k][2]+atom_scaling_array[i];        // rescale
				(*general_angles_registry)[j][3]=molecule_angles_registry[k][3]+atom_scaling_array[i];        // rescale
			}
			// free extra memory
			for(j=(*general_angles-1);j>(*general_angles+tcv-1);--j)
			{
				free((*general_angles_registry)[j]);
			}
		}
		else
		{
            // size kept the same
			if(verb==1)printf("$ case 2 init / final sizes: %d / %d\n",*general_angles,*general_angles+tcv);
			// place current molecule entries
			k=-1;
			for(j=topo_boundaries[i][2];j<=topo_boundaries[i][3];++j)
			{
				k=k+1;
				if(verb==1)printf("$ placing %d at %d...\n",k,j);
				(*general_angles_registry)[j][0]=(*atmapping_matrix)[molecule_angles_registry[k][0]-1][1];    // type mapping
				(*general_angles_registry)[j][1]=molecule_angles_registry[k][1]+atom_scaling_array[i];        // rescale
				(*general_angles_registry)[j][2]=molecule_angles_registry[k][2]+atom_scaling_array[i];        // rescale
				(*general_angles_registry)[j][3]=molecule_angles_registry[k][3]+atom_scaling_array[i];        // rescale
			}
		}
		/*
        // place current molecule entries
		k=-1;
		for(j=topo_boundaries[i][2];j<=topo_boundaries[i][3];++j)
		{
			k=k+1;
			printf("$ placing %d at %d...\n",k,j);
			(*general_angles_registry)[j][0]=(*atmapping_matrix)[molecule_angles_registry[k][0]-1][1];  // utilize type mapping
			(*general_angles_registry)[j][1]=molecule_angles_registry[k][1]+atom_scaling_array[i];      // rescale ids
			(*general_angles_registry)[j][2]=molecule_angles_registry[k][2]+atom_scaling_array[i];      // rescale ids
			(*general_angles_registry)[j][3]=molecule_angles_registry[k][3]+atom_scaling_array[i];      // rescale ids
		}
		*/
        // alter total
		*general_angles=*general_angles+tcv;
        //
		if(verb==1)printf("$ general angles registry:\n");
		if(verb==1)for(j=0;j<*general_angles;++j)printf("[%d]\t%d\t%d\t%d\t%d\n",j+1,(*general_angles_registry)[j][0],(*general_angles_registry)[j][1],(*general_angles_registry)[j][2],(*general_angles_registry)[j][3]);		
	}
	
	// case 3
    // this case is triggered for angles only when the current molecule is the *FIRST* molecule!!!!
	else if(non_minus_above==0 && non_minus_below==1)
	{
        // delta holds the number of angles prior to the building step
		delta=0;if(non_minus_current==1)delta=topo_boundaries[i][3]-topo_boundaries[i][2]+1;
		if(verb==1)printf("$ case 3 delta = %d\n",delta);
        // upper_0 and lower_0 define the boundaries for the values below
		for(j=i+1;j<molecules;++j)if(topo_boundaries[j][2]!=-1){upper_0=topo_boundaries[j][2];break;}
		lower_0=topo_boundaries[i][3];
		for(j=i+1;j<molecules;++j)if(topo_boundaries[j][3]>lower_0){lower_0=topo_boundaries[j][3];}
		if(verb==1){printf("$ case 3 upper_0 = %d\n",upper_0);
		printf("$ case 3 lower_0 = %d\n",lower_0);}
        // in_array holds the new values for the current molecule due to the building step
		in_array[0]=0;in_array[1]=molecule_angles-1;
		if(verb==1)printf("$ case 3 in_array = [%d %d]\n",in_array[0],in_array[1]);
        // update current molecule values
		topo_boundaries[i][2]=in_array[0];
		topo_boundaries[i][3]=in_array[1];
        // tcv (topo current value) holds the difference in elements (delta_angles)
		tcv=topo_boundaries[i][3]-topo_boundaries[i][2]+1-delta;
		if(verb==1)printf("$ case 3 tcv = %d\n",tcv);
        // update non -1 values below current molecule
		for(j=i+1;j<molecules;++j){if(topo_boundaries[j][2]!=-1){topo_boundaries[j][2]=topo_boundaries[j][2]+tcv;topo_boundaries[j][3]=topo_boundaries[j][3]+tcv;}}
        //
		if(verb==1){printf("$ case 3 updated topo_boundaries:\n");
		for(j=0;j<molecules;++j)printf("[%d]\t%d\t%d\n",j+1,topo_boundaries[j][2],topo_boundaries[j][3]);}
        //
		if(tcv>0)
		{
            // resize adding extra space
			if(verb==1)printf("$ case 3 init / final sizes: %d / %d\n",*general_angles,*general_angles+tcv);
			*general_angles_registry=(int**)realloc(*general_angles_registry,(*general_angles+tcv)*sizeof(int*));
			for(j=0;j<tcv;++j)(*general_angles_registry)[*general_angles+j]=(int*)malloc(4*sizeof(int));
            // push down entries
			if(verb==1)printf("$ case 3 push down: %d:%d --> %d:%d\n",upper_0,lower_0,upper_0+tcv,lower_0+tcv);
			for(j=lower_0;j>=upper_0;--j)
			{
				if(verb==1)printf("%d --> %d\n",j,j+tcv);
				(*general_angles_registry)[j+tcv][0]=(*general_angles_registry)[j][0];
				(*general_angles_registry)[j+tcv][1]=(*general_angles_registry)[j][1]+delta_atoms;  // rescale
				(*general_angles_registry)[j+tcv][2]=(*general_angles_registry)[j][2]+delta_atoms;  // rescale
				(*general_angles_registry)[j+tcv][3]=(*general_angles_registry)[j][3]+delta_atoms;  // rescale
			}
			// place current molecule entries
			k=-1;
			for(j=topo_boundaries[i][2];j<=topo_boundaries[i][3];++j)
			{
				k=k+1;
				if(verb==1)printf("$ placing %d at %d...\n",k,j);
				(*general_angles_registry)[j][0]=(*atmapping_matrix)[molecule_angles_registry[k][0]-1][1];    // type mapping
				(*general_angles_registry)[j][1]=molecule_angles_registry[k][1]+atom_scaling_array[i];        // rescale
				(*general_angles_registry)[j][2]=molecule_angles_registry[k][2]+atom_scaling_array[i];        // rescale
				(*general_angles_registry)[j][3]=molecule_angles_registry[k][3]+atom_scaling_array[i];        // rescale
			}
		}
		else if (tcv<0)
		{
            // resize removing space
			if(verb==1){printf("$ case 3 init / final sizes: %d / %d\n",*general_angles,*general_angles+tcv);
			printf("\n$ tcv<0\n\n");
			printf("\n$ topo_boundaries:\t%d\t%d\n",topo_boundaries[i][2],topo_boundaries[i][3]);}
            // pull up!!
            for(j=upper_0;j<=lower_0;++j)
            {
				//printf("%d --> %d\n",j,j+tcv);
				(*general_angles_registry)[j+tcv][0]=(*general_angles_registry)[j][0];
				(*general_angles_registry)[j+tcv][1]=(*general_angles_registry)[j][1]+delta_atoms;
				(*general_angles_registry)[j+tcv][2]=(*general_angles_registry)[j][2]+delta_atoms;
				(*general_angles_registry)[j+tcv][3]=(*general_angles_registry)[j][3]+delta_atoms;
			}
            // place current molecule entries
			k=-1;
			for(j=topo_boundaries[i][2];j<=topo_boundaries[i][3];++j)
			{
				k=k+1;
				if(verb==1)printf("$ placing %d at %d...\n",k,j);
				(*general_angles_registry)[j][0]=(*atmapping_matrix)[molecule_angles_registry[k][0]-1][1];    // type mapping
				(*general_angles_registry)[j][1]=molecule_angles_registry[k][1]+atom_scaling_array[i];        // rescale
				(*general_angles_registry)[j][2]=molecule_angles_registry[k][2]+atom_scaling_array[i];        // rescale
				(*general_angles_registry)[j][3]=molecule_angles_registry[k][3]+atom_scaling_array[i];        // rescale
			}
            // free extra memory
			for(j=(*general_angles-1);j>(*general_angles+tcv-1);--j)
			{
				free((*general_angles_registry)[j]);
			}
		}
		else
		{
            // size kept the same
			if(verb==1){printf("$ case 3 init / final sizes: %d / %d\n",*general_angles,*general_angles+tcv);
            // update for rescaling purposes
			printf("$ case 3 no shift: %d:%d <--> %d:%d\n",upper_0,lower_0,upper_0+tcv,lower_0+tcv);}
			for(j=lower_0;j>=upper_0;--j)
			{
				if(verb==1)printf("%d --> %d\n",j,j+tcv);
				(*general_angles_registry)[j+tcv][0]=(*general_angles_registry)[j][0];
				(*general_angles_registry)[j+tcv][1]=(*general_angles_registry)[j][1]+delta_atoms;  // rescale
				(*general_angles_registry)[j+tcv][2]=(*general_angles_registry)[j][2]+delta_atoms;  // rescale
				(*general_angles_registry)[j+tcv][3]=(*general_angles_registry)[j][3]+delta_atoms;  // rescale
			}
			// place current molecule entries
			k=-1;
			for(j=topo_boundaries[i][2];j<=topo_boundaries[i][3];++j)
			{
				k=k+1;
				if(verb==1)printf("$ placing %d at %d...\n",k,j);
				(*general_angles_registry)[j][0]=(*atmapping_matrix)[molecule_angles_registry[k][0]-1][1];    // type mapping
				(*general_angles_registry)[j][1]=molecule_angles_registry[k][1]+atom_scaling_array[i];        // rescale
				(*general_angles_registry)[j][2]=molecule_angles_registry[k][2]+atom_scaling_array[i];        // rescale
				(*general_angles_registry)[j][3]=molecule_angles_registry[k][3]+atom_scaling_array[i];        // rescale
			}
		}
		/*
        // place current molecule entries
		k=-1;
		for(j=topo_boundaries[i][2];j<=topo_boundaries[i][3];++j)
		{
			k=k+1;
			printf("$ placing %d at %d...\n",k,j);
			(*general_angles_registry)[j][0]=(*atmapping_matrix)[molecule_angles_registry[k][0]-1][1];  // type mapping
			(*general_angles_registry)[j][1]=molecule_angles_registry[k][1]+atom_scaling_array[i];      // rescale
			(*general_angles_registry)[j][2]=molecule_angles_registry[k][2]+atom_scaling_array[i];      // rescale
			(*general_angles_registry)[j][3]=molecule_angles_registry[k][3]+atom_scaling_array[i];      // rescale
		}
		*/ 
        // alter total
		*general_angles=*general_angles+tcv;
        //
		if(verb==1)printf("$ general angles registry:\n");
		if(verb==1)for(j=0;j<*general_angles;++j)printf("[%d]\t%d\t%d\t%d\t%d\n",j+1,(*general_angles_registry)[j][0],(*general_angles_registry)[j][1],(*general_angles_registry)[j][2],(*general_angles_registry)[j][3]);
	}
	
	// case 4
    // active for all intermediate cases...
	else if(non_minus_above==1 && non_minus_below==1)
	{
        // delta holds the number of angles prior to the building step
		delta=0;if(non_minus_current==1)delta=topo_boundaries[i][3]-topo_boundaries[i][2]+1;
		if(verb==1)printf("$ case 4 delta = %d\n",delta);
        // upper_0 and lower_0 define the boundaries for the values below
		for(j=i+1;j<molecules;++j)if(topo_boundaries[j][2]!=-1){upper_0=topo_boundaries[j][2];break;}
		lower_0=topo_boundaries[i][3];
		for(j=i+1;j<molecules;++j)if(topo_boundaries[j][3]>lower_0){lower_0=topo_boundaries[j][3];}
		if(verb==1){printf("$ case 4 upper_0 = %d\n",upper_0);
		printf("$ case 4 lower_0 = %d\n",lower_0);}
        // tav (topo above value) holds the largest eligible position value above the current molecule
		tav=topo_boundaries[0][3];for(j=0;j<=i-1;++j)if(topo_boundaries[j][3]>tav)tav=topo_boundaries[j][3];
		if(verb==1)printf("$ case 4 tav = %d\n",tav);
        // in_array holds the new values for the current molecule due to the building step
		in_array[0]=0+tav+1;in_array[1]=molecule_angles-1+tav+1;
		if(verb==1)printf("$ case 4 in_array = [%d %d]\n",in_array[0],in_array[1]);
        // update current molecule values
		topo_boundaries[i][2]=in_array[0];
		topo_boundaries[i][3]=in_array[1];
        // tcv (topo current value) holds the difference in elements (delta_angles)
		tcv=topo_boundaries[i][3]-topo_boundaries[i][2]+1-delta;
		if(verb==1)printf("$ case 4 tcv = %d\n",tcv);
        // update non -1 values below current molecule
		for(j=i+1;j<molecules;++j){if(topo_boundaries[j][2]!=-1){topo_boundaries[j][2]=topo_boundaries[j][2]+tcv;topo_boundaries[j][3]=topo_boundaries[j][3]+tcv;}}
		//
        if(verb==1){printf("$ case 4 updated topo_boundaries:\n");
		for(j=0;j<molecules;++j)printf("[%d]\t%d\t%d\n",j+1,topo_boundaries[j][2],topo_boundaries[j][3]);}
		//
        if(tcv>0)
		{
            // resize adding extra space
			if(verb==1)printf("$ case 4 init / final sizes: %d / %d\n",*general_angles,*general_angles+tcv);
			*general_angles_registry=(int**)realloc(*general_angles_registry,(*general_angles+tcv)*sizeof(int*));
			for(j=0;j<tcv;++j)(*general_angles_registry)[*general_angles+j]=(int*)malloc(4*sizeof(int));
            // push down entries
			if(verb==1)printf("$ case 4 push down: %d:%d --> %d:%d\n",upper_0,lower_0,upper_0+tcv,lower_0+tcv);
			for(j=lower_0;j>=upper_0;--j)
			{
				if(verb==1)printf("%d --> %d\n",j,j+tcv);
				(*general_angles_registry)[j+tcv][0]=(*general_angles_registry)[j][0];
				(*general_angles_registry)[j+tcv][1]=(*general_angles_registry)[j][1]+delta_atoms;  // rescale
				(*general_angles_registry)[j+tcv][2]=(*general_angles_registry)[j][2]+delta_atoms;  // rescale
				(*general_angles_registry)[j+tcv][3]=(*general_angles_registry)[j][3]+delta_atoms;  // rescale
			}
			// place current molecule entries
			k=-1;
			for(j=topo_boundaries[i][2];j<=topo_boundaries[i][3];++j)
			{
				k=k+1;
				if(verb==1)printf("$ placing %d at %d...\n",k,j);
				(*general_angles_registry)[j][0]=(*atmapping_matrix)[molecule_angles_registry[k][0]-1][1];    // type mapping
				(*general_angles_registry)[j][1]=molecule_angles_registry[k][1]+atom_scaling_array[i];        // rescale
				(*general_angles_registry)[j][2]=molecule_angles_registry[k][2]+atom_scaling_array[i];        // rescale
				(*general_angles_registry)[j][3]=molecule_angles_registry[k][3]+atom_scaling_array[i];        // rescale
			}
		}
		else if (tcv<0)
		{
            // resize removing space
			if(verb==1){printf("$ case 4 init / final sizes: %d / %d\n",*general_angles,*general_angles+tcv);
			printf("\n$ tcv<0\n\n");}
			// pull up!!
            for(j=upper_0;j<=lower_0;++j)
            {
				//printf("%d --> %d\n",j,j+tcv);
				(*general_angles_registry)[j+tcv][0]=(*general_angles_registry)[j][0];
				(*general_angles_registry)[j+tcv][1]=(*general_angles_registry)[j][1]+delta_atoms;
				(*general_angles_registry)[j+tcv][2]=(*general_angles_registry)[j][2]+delta_atoms;
				(*general_angles_registry)[j+tcv][3]=(*general_angles_registry)[j][3]+delta_atoms;
			}
            // place current molecule entries
			k=-1;
			for(j=topo_boundaries[i][2];j<=topo_boundaries[i][3];++j)
			{
				k=k+1;
				if(verb==1)printf("$ placing %d at %d...\n",k,j);
				(*general_angles_registry)[j][0]=(*atmapping_matrix)[molecule_angles_registry[k][0]-1][1];    // type mapping
				(*general_angles_registry)[j][1]=molecule_angles_registry[k][1]+atom_scaling_array[i];        // rescale
				(*general_angles_registry)[j][2]=molecule_angles_registry[k][2]+atom_scaling_array[i];        // rescale
				(*general_angles_registry)[j][3]=molecule_angles_registry[k][3]+atom_scaling_array[i];        // rescale
			}
			// free extra memory
			for(j=(*general_angles-1);j>(*general_angles+tcv-1);--j)
			{
				free((*general_angles_registry)[j]);
			}
		}
		else
		{
            // size kept the same
			if(verb==1){printf("$ case 4 init / final sizes: %d / %d\n",*general_angles,*general_angles+tcv);
			printf("$ case 4 no shift: %d:%d <--> %d:%d\n",upper_0,lower_0,upper_0+tcv,lower_0+tcv);
            // update for rescaling purposes
			printf("$ case 4 no shift: %d:%d <--> %d:%d\n",upper_0,lower_0,upper_0+tcv,lower_0+tcv);}
			for(j=lower_0;j>=upper_0;--j)
			{
				if(verb==1)printf("%d --> %d\n",j,j+tcv);
				(*general_angles_registry)[j+tcv][0]=(*general_angles_registry)[j][0];
				(*general_angles_registry)[j+tcv][1]=(*general_angles_registry)[j][1]+delta_atoms;  // rescale
				(*general_angles_registry)[j+tcv][2]=(*general_angles_registry)[j][2]+delta_atoms;  // rescale
				(*general_angles_registry)[j+tcv][3]=(*general_angles_registry)[j][3]+delta_atoms;  // rescale
			}
			// place current molecule entries
			k=-1;
			for(j=topo_boundaries[i][2];j<=topo_boundaries[i][3];++j)
			{
				k=k+1;
				if(verb==1)printf("$ placing %d at %d...\n",k,j);
				(*general_angles_registry)[j][0]=(*atmapping_matrix)[molecule_angles_registry[k][0]-1][1];    // type mapping
				(*general_angles_registry)[j][1]=molecule_angles_registry[k][1]+atom_scaling_array[i];        // rescale
				(*general_angles_registry)[j][2]=molecule_angles_registry[k][2]+atom_scaling_array[i];        // rescale
				(*general_angles_registry)[j][3]=molecule_angles_registry[k][3]+atom_scaling_array[i];        // rescale
			}
		}
		/*
        // place current molecule entries
		k=-1;
		for(j=topo_boundaries[i][2];j<=topo_boundaries[i][3];++j)
		{
			k=k+1;
			printf("$ placing %d at %d...\n",k,j);
			(*general_angles_registry)[j][0]=(*atmapping_matrix)[molecule_angles_registry[k][0]-1][1];  // type mapping
			(*general_angles_registry)[j][1]=molecule_angles_registry[k][1]+atom_scaling_array[i];      // rescale
			(*general_angles_registry)[j][2]=molecule_angles_registry[k][2]+atom_scaling_array[i];      // rescale
			(*general_angles_registry)[j][3]=molecule_angles_registry[k][3]+atom_scaling_array[i];      // rescale
		}
		*/
        // alter total
		*general_angles=*general_angles+tcv;
        //
		if(verb==1)printf("$ general angles registry:\n");
		if(verb==1)for(j=0;j<*general_angles;++j)printf("[%d]\t%d\t%d\t%d\t%d\n",j+1,(*general_angles_registry)[j][0],(*general_angles_registry)[j][1],(*general_angles_registry)[j][2],(*general_angles_registry)[j][3]);
	}
	else
	{
		printf("unsupported angles combination!!\n\n");exit(-1);
	}

	delete_flag=0;
    for(k=0;k<*atmap_rows;++k)
        if(acm_sum[k]==0)
            delete_flag=delete_flag+1;
    
    if(delete_flag>0)
    {
        if(verb==1)printf("\n$ delete flag>0!!!\n\n");
        gtmapping_array=(int*)malloc((*general_angle_types)*sizeof(int));
        backup=(char***)malloc((*general_angle_types)*sizeof(char**));
        for(j=0;j<*general_angle_types;++j)backup[j]=(char**)malloc(5*sizeof(char*));
        for(j=0;j<*general_angle_types;++j)for(k=0;k<5;++k)backup[j][k]=(char*)malloc(sub_length*sizeof(char));
        counter=0;
        for(j=0;j<*general_angle_types;++j)
        {
            for(k=0;k<*atmap_rows;++k)
            {
                if((strcmp((*general_angle_types_registry)[j][0],(*atmap)[k][0])==0 &&
                    strcmp((*general_angle_types_registry)[j][1],(*atmap)[k][1])==0 &&
                    strcmp((*general_angle_types_registry)[j][2],(*atmap)[k][2])==0 &&
                    strcmp((*general_angle_types_registry)[j][3],(*atmap)[k][3])==0 &&
                    strcmp((*general_angle_types_registry)[j][4],(*atmap)[k][4])==0)
                   ||
                   (strcmp((*general_angle_types_registry)[j][0],(*atmap)[k][4])==0 &&
                    strcmp((*general_angle_types_registry)[j][1],(*atmap)[k][3])==0 &&
                    strcmp((*general_angle_types_registry)[j][2],(*atmap)[k][2])==0 &&
                    strcmp((*general_angle_types_registry)[j][3],(*atmap)[k][1])==0 &&
                    strcmp((*general_angle_types_registry)[j][4],(*atmap)[k][0])==0))
                    break;
            }
            if(verb==1){printf("$ (%s)--%s--(%s)--%s--(%s) is at position {%d}\n",(*general_angle_types_registry)[j][0],(*general_angle_types_registry)[j][1],(*general_angle_types_registry)[j][2],(*general_angle_types_registry)[j][3],(*general_angle_types_registry)[j][4],k+1);
            printf("$ checking cm_sum: %d\n",acm_sum[k]);}
            if(acm_sum[k]==0)
            {
                gtmapping_array[j]=0;
            }
            else
            {
                counter=counter+1;
                gtmapping_array[j]=counter;
                sprintf(backup[counter-1][0],"%s",(*general_angle_types_registry)[j][0]);
                sprintf(backup[counter-1][1],"%s",(*general_angle_types_registry)[j][1]);
                sprintf(backup[counter-1][2],"%s",(*general_angle_types_registry)[j][2]);
                sprintf(backup[counter-1][3],"%s",(*general_angle_types_registry)[j][3]);
                sprintf(backup[counter-1][4],"%s",(*general_angle_types_registry)[j][4]);
            }
            
        }
        if(verb==1){
        printf("$ gtmapping_array:\n");
        for(j=0;j<*general_angle_types;++j)printf("[%d]\t%d\n",j+1,gtmapping_array[j]);
        printf("$ altered types:\n");
        for(j=0;j<counter;++j)printf("[%d]\t(%s)--%s--(%s)--%s--(%s)\n",j+1,backup[j][0],backup[j][1],backup[j][2],backup[j][3],backup[j][4]);}
        
        for(j=0;j<*general_angles;++j)(*general_angles_registry)[j][0]=gtmapping_array[(*general_angles_registry)[j][0]-1];
        for(j=0;j<counter;++j)
        {
            sprintf((*general_angle_types_registry)[j][0],"%s",backup[j][0]);
            sprintf((*general_angle_types_registry)[j][1],"%s",backup[j][1]);
            sprintf((*general_angle_types_registry)[j][2],"%s",backup[j][2]);
            sprintf((*general_angle_types_registry)[j][3],"%s",backup[j][3]);
            sprintf((*general_angle_types_registry)[j][4],"%s",backup[j][4]);
        }
        for(j=counter;j<*general_angle_types;++j)
        {
            for(k=0;k<5;++k)
                free((*general_angle_types_registry)[j][k]);
            free((*general_angle_types_registry)[j]);
        }
        for(j=0;j<*general_angle_types;++j)
            for(k=0;k<5;++k)
                free(backup[j][k]);
        for(j=0;j<*general_angle_types;++j)free(backup[j]);
        free(backup);
        *general_angle_types=counter;
        free(gtmapping_array);
    }
    
    free(acm_sum);

    //--------------------------------------------------------------------------
    
    // alter dihedrals
    
    // [3]
    
    // update tmap and cm
    for(j=0;j<molecule_dihedral_types;++j)
    {
        found=0;
        for(k=0;k<*dtmap_rows;++k)
        {
            if(
               (strcmp(molecule_dihedral_types_registry[j][0],(*dtmap)[k][0])==0 &&
                strcmp(molecule_dihedral_types_registry[j][1],(*dtmap)[k][1])==0 &&
                strcmp(molecule_dihedral_types_registry[j][2],(*dtmap)[k][2])==0 &&
                strcmp(molecule_dihedral_types_registry[j][3],(*dtmap)[k][3])==0 &&
                strcmp(molecule_dihedral_types_registry[j][4],(*dtmap)[k][4])==0 &&
                strcmp(molecule_dihedral_types_registry[j][5],(*dtmap)[k][5])==0 &&
                strcmp(molecule_dihedral_types_registry[j][6],(*dtmap)[k][6])==0
                )||
               (strcmp(molecule_dihedral_types_registry[j][0],(*dtmap)[k][6])==0 &&
                strcmp(molecule_dihedral_types_registry[j][1],(*dtmap)[k][5])==0 &&
                strcmp(molecule_dihedral_types_registry[j][2],(*dtmap)[k][4])==0 &&
                strcmp(molecule_dihedral_types_registry[j][3],(*dtmap)[k][3])==0 &&
                strcmp(molecule_dihedral_types_registry[j][4],(*dtmap)[k][2])==0 &&
                strcmp(molecule_dihedral_types_registry[j][5],(*dtmap)[k][1])==0 &&
                strcmp(molecule_dihedral_types_registry[j][6],(*dtmap)[k][0])==0
                )
               )
            {
                found=1;break;
            }
        }
        if(found==0)
        {
            // augment
            *dtmap_rows=*dtmap_rows+1;
            // tmap
            // resize
            *dtmap=(char***)realloc(*dtmap,(*dtmap_rows)*sizeof(char**));
            (*dtmap)[*dtmap_rows-1]=(char**)malloc(7*sizeof(char*));
            (*dtmap)[*dtmap_rows-1][0]=(char*)malloc(sub_length*sizeof(char));
            (*dtmap)[*dtmap_rows-1][1]=(char*)malloc(sub_length*sizeof(char));
            (*dtmap)[*dtmap_rows-1][2]=(char*)malloc(sub_length*sizeof(char));
            (*dtmap)[*dtmap_rows-1][3]=(char*)malloc(sub_length*sizeof(char));
            (*dtmap)[*dtmap_rows-1][4]=(char*)malloc(sub_length*sizeof(char));
            (*dtmap)[*dtmap_rows-1][5]=(char*)malloc(sub_length*sizeof(char));
            (*dtmap)[*dtmap_rows-1][6]=(char*)malloc(sub_length*sizeof(char));
            // write
            // write type
            sprintf((*dtmap)[*dtmap_rows-1][0],"%s",molecule_dihedral_types_registry[j][0]);
            sprintf((*dtmap)[*dtmap_rows-1][1],"%s",molecule_dihedral_types_registry[j][1]);
            sprintf((*dtmap)[*dtmap_rows-1][2],"%s",molecule_dihedral_types_registry[j][2]);
            sprintf((*dtmap)[*dtmap_rows-1][3],"%s",molecule_dihedral_types_registry[j][3]);
            sprintf((*dtmap)[*dtmap_rows-1][4],"%s",molecule_dihedral_types_registry[j][4]);
            sprintf((*dtmap)[*dtmap_rows-1][5],"%s",molecule_dihedral_types_registry[j][5]);
            sprintf((*dtmap)[*dtmap_rows-1][6],"%s",molecule_dihedral_types_registry[j][6]);
            // cm
            for(k=0;k<molecules;++k)
            {
                (*dcm)[k]=(int*)realloc((*dcm)[k],(*dtmap_rows)*sizeof(int*));
                (*dcm)[k][*dtmap_rows-1]=0;
            }
        }
    }
    
    // [4]
    
    // reset cm(i,:)
    for(j=0;j<*dtmap_rows;++j)(*dcm)[i][j]=0;
    // loop on mdtr (molecule dihedral type registry) using dtmap and update dcm
    for(j=0;j<molecule_dihedral_types;++j)
    {
        for(k=0;k<*dtmap_rows;++k)
            if((
                strcmp(molecule_dihedral_types_registry[j][0],(*dtmap)[k][0])==0 &&
                strcmp(molecule_dihedral_types_registry[j][1],(*dtmap)[k][1])==0 &&
                strcmp(molecule_dihedral_types_registry[j][2],(*dtmap)[k][2])==0 &&
                strcmp(molecule_dihedral_types_registry[j][3],(*dtmap)[k][3])==0 &&
                strcmp(molecule_dihedral_types_registry[j][4],(*dtmap)[k][4])==0 &&
                strcmp(molecule_dihedral_types_registry[j][5],(*dtmap)[k][5])==0 &&
                strcmp(molecule_dihedral_types_registry[j][6],(*dtmap)[k][6])==0
                )
               ||
               (
                strcmp(molecule_dihedral_types_registry[j][0],(*dtmap)[k][6])==0 &&
                strcmp(molecule_dihedral_types_registry[j][1],(*dtmap)[k][5])==0 &&
                strcmp(molecule_dihedral_types_registry[j][2],(*dtmap)[k][4])==0 &&
                strcmp(molecule_dihedral_types_registry[j][3],(*dtmap)[k][3])==0 &&
                strcmp(molecule_dihedral_types_registry[j][4],(*dtmap)[k][2])==0 &&
                strcmp(molecule_dihedral_types_registry[j][5],(*dtmap)[k][1])==0 &&
                strcmp(molecule_dihedral_types_registry[j][6],(*dtmap)[k][0])==0
                )
               )
                break;
        (*dcm)[i][k]=1;
    }
    
     // console out
     if(verb==1){
     printf("\n$ dtmap/dcm:\n");
     for(j=0;j<*dtmap_rows;++j)printf("(%s)--%s--(%s)--%s--(%s)--%s--(%s)\t",(*dtmap)[j][0],(*dtmap)[j][1],(*dtmap)[j][2],(*dtmap)[j][3],(*dtmap)[j][4],(*dtmap)[j][5],(*dtmap)[j][6]);printf("\n");
     for(j=0;j<molecules;++j){
     for(k=0;k<*dtmap_rows;++k)
     printf("%d\t",(*dcm)[j][k]);
     printf("\n");
     }
	}
    // [5]
     
    // global dihedral types through dtmap/dcm:
    dcm_sum=(int*)malloc((*dtmap_rows)*sizeof(int));
    for(k=0;k<*dtmap_rows;++k)
    {
        dcm_sum[k]=0;
        for(j=0;j<molecules;++j)
        {
            dcm_sum[k]=dcm_sum[k]+(*dcm)[j][k];
        }
    }
    
    if(verb==1){
    printf("\n$ dcm_sum array:\n");for(j=0;j<*dtmap_rows;++j)printf("%d\t",dcm_sum[j]);printf("\n");}
    
    // [6]
    
	*dtmapping_matrix=(int**)malloc(molecule_dihedral_types*sizeof(int*));
	for(j=0;j<molecule_dihedral_types;++j)(*dtmapping_matrix)[j]=(int*)malloc(2*sizeof(int));
	
	n=-1;
	for(j=0;j<molecule_dihedral_types;++j)
	{
		found=0;
		for(k=0;k<*general_dihedral_types;++k)
		{
			if(
			   (strcmp(molecule_dihedral_types_registry[j][0],(*general_dihedral_types_registry)[k][0])==0 &&
				strcmp(molecule_dihedral_types_registry[j][1],(*general_dihedral_types_registry)[k][1])==0 &&
				strcmp(molecule_dihedral_types_registry[j][2],(*general_dihedral_types_registry)[k][2])==0 &&
				strcmp(molecule_dihedral_types_registry[j][3],(*general_dihedral_types_registry)[k][3])==0 &&
                strcmp(molecule_dihedral_types_registry[j][4],(*general_dihedral_types_registry)[k][4])==0 &&
                strcmp(molecule_dihedral_types_registry[j][5],(*general_dihedral_types_registry)[k][5])==0 &&
                strcmp(molecule_dihedral_types_registry[j][6],(*general_dihedral_types_registry)[k][6])==0
				)
			   ||
               (strcmp(molecule_dihedral_types_registry[j][0],(*general_dihedral_types_registry)[k][6])==0 &&
                strcmp(molecule_dihedral_types_registry[j][1],(*general_dihedral_types_registry)[k][5])==0 &&
                strcmp(molecule_dihedral_types_registry[j][2],(*general_dihedral_types_registry)[k][4])==0 &&
                strcmp(molecule_dihedral_types_registry[j][3],(*general_dihedral_types_registry)[k][3])==0 &&
                strcmp(molecule_dihedral_types_registry[j][4],(*general_dihedral_types_registry)[k][2])==0 &&
                strcmp(molecule_dihedral_types_registry[j][5],(*general_dihedral_types_registry)[k][1])==0 &&
                strcmp(molecule_dihedral_types_registry[j][6],(*general_dihedral_types_registry)[k][0])==0
                )
			   )
			{
				n=n+1;
				//printf("* molecular %d is found at general %d\n",j+1,k+1);
				(*dtmapping_matrix)[n][0]=j+1;
				(*dtmapping_matrix)[n][1]=k+1;
				found=1;
				break;
			}
		}
		if(found==0)
		{
			// augment
			*general_dihedral_types=*general_dihedral_types+1;
			// resize
			*general_dihedral_types_registry=(char***)realloc(*general_dihedral_types_registry,(*general_dihedral_types)*sizeof(char**));
			(*general_dihedral_types_registry)[*general_dihedral_types-1]=(char**)malloc(7*sizeof(char*));
			(*general_dihedral_types_registry)[*general_dihedral_types-1][0]=(char*)malloc(sub_length*sizeof(char));
			(*general_dihedral_types_registry)[*general_dihedral_types-1][1]=(char*)malloc(sub_length*sizeof(char));
			(*general_dihedral_types_registry)[*general_dihedral_types-1][2]=(char*)malloc(sub_length*sizeof(char));
			(*general_dihedral_types_registry)[*general_dihedral_types-1][3]=(char*)malloc(sub_length*sizeof(char));
            (*general_dihedral_types_registry)[*general_dihedral_types-1][4]=(char*)malloc(sub_length*sizeof(char));
            (*general_dihedral_types_registry)[*general_dihedral_types-1][5]=(char*)malloc(sub_length*sizeof(char));
            (*general_dihedral_types_registry)[*general_dihedral_types-1][6]=(char*)malloc(sub_length*sizeof(char));
			// write
			// write type
			sprintf((*general_dihedral_types_registry)[*general_dihedral_types-1][0],"%s",molecule_dihedral_types_registry[j][0]);
			sprintf((*general_dihedral_types_registry)[*general_dihedral_types-1][1],"%s",molecule_dihedral_types_registry[j][1]);
			sprintf((*general_dihedral_types_registry)[*general_dihedral_types-1][2],"%s",molecule_dihedral_types_registry[j][2]);
			sprintf((*general_dihedral_types_registry)[*general_dihedral_types-1][3],"%s",molecule_dihedral_types_registry[j][3]);
            sprintf((*general_dihedral_types_registry)[*general_dihedral_types-1][4],"%s",molecule_dihedral_types_registry[j][4]);
            sprintf((*general_dihedral_types_registry)[*general_dihedral_types-1][5],"%s",molecule_dihedral_types_registry[j][5]);
            sprintf((*general_dihedral_types_registry)[*general_dihedral_types-1][6],"%s",molecule_dihedral_types_registry[j][6]);
			//
			n=n+1;
			//printf("~ molecular %d is placed at general %d\n",j+1,*general_angle_types);
			(*dtmapping_matrix)[n][0]=j+1;
			(*dtmapping_matrix)[n][1]=*general_dihedral_types;
		}
		
	}
	if(verb==1){
	printf("\n$ dtmapping_matrix:\n");
	for(j=0;j<molecule_dihedral_types;++j)printf("~ map %d --> %d\n",(*dtmapping_matrix)[j][0],(*dtmapping_matrix)[j][1]);

	printf("\n$ general dihedral types registry:\n");
	for(j=0;j<*general_dihedral_types;++j)printf("[%d]\t(%s)--%s--(%s)--%s--(%s)--%s--(%s)\n",j+1,(*general_dihedral_types_registry)[j][0],(*general_dihedral_types_registry)[j][1],(*general_dihedral_types_registry)[j][2],(*general_dihedral_types_registry)[j][3],(*general_dihedral_types_registry)[j][4],(*general_dihedral_types_registry)[j][5],(*general_dihedral_types_registry)[j][6]);
	}
	// [7]
	
    if(molecule_dihedrals==0){printf("\n$ molecule_dihedrals = 0 ??\n\n");exit(-1);}
    
	// initialize
	non_minus_above=0;
	non_minus_below=0;
	non_minus_current=0;
	// check current molecule i even column registry
	if(topo_boundaries[i][4]!=-1)non_minus_current=1;
	// traverse even column from 0 to i-1
	for(j=0;j<=i-1;++j)if(topo_boundaries[j][4]!=-1){non_minus_above=1;break;}
	// traverse even column from i+1 to molecules
	for(j=i+1;j<molecules;++j)if(topo_boundaries[j][4]!=-1){non_minus_below=1;break;}
	//
	if(verb==1){
	printf("\n$ minus flags:\n");
	printf("non_minus_above = %d\n",non_minus_above);
	printf("non_minus_below = %d\n",non_minus_below);
	printf("non_minus_current = %d\n\n",non_minus_current);}
    
    // case 1
    // this case is triggered when a dihedral is first introduced to the system, e.g. CH4 --> 2HC==CH2
    // or when there is only one dihedral in an intermediate slot
    if(non_minus_above==0 && non_minus_below==0)
    {
		if(non_minus_current==0)
		{
            // in_array holds the new values for the current molecule due to the building step
			in_array[0]=0;in_array[1]=molecule_dihedrals-1;
			if(verb==1)printf("$ case 1 in_array = [%d %d]\n",in_array[0],in_array[1]);
            // update current molecule values
			topo_boundaries[i][4]=in_array[0];
			topo_boundaries[i][5]=in_array[1];
            // tcv (topo current value) holds the difference in elements (delta_dihedrals)
			tcv=topo_boundaries[i][5]-topo_boundaries[i][4]+1;
			if(verb==1){printf("$ case 1 tcv = %d\n",tcv);
			printf("$ case 1 updated topo_boundaries:\n");
			for(j=0;j<molecules;++j)printf("[%d]\t%d\t%d\n",j+1,topo_boundaries[j][4],topo_boundaries[j][5]);}
            // allocate space
			free(*general_dihedrals_registry);*general_dihedrals_registry=(int**)malloc(tcv*sizeof(int*));
			for(j=0;j<tcv;++j)(*general_dihedrals_registry)[j]=(int*)malloc(5*sizeof(int));
            // place current molecule entries
			k=-1;
			for(j=topo_boundaries[i][4];j<=topo_boundaries[i][5];++j)
			{
				k=k+1;
				if(verb==1)printf("$ placing %d at %d...\n",k,j);
				(*general_dihedrals_registry)[j][0]=(*dtmapping_matrix)[molecule_dihedrals_registry[k][0]-1][1];    // type mapping
				(*general_dihedrals_registry)[j][1]=molecule_dihedrals_registry[k][1]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][2]=molecule_dihedrals_registry[k][2]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][3]=molecule_dihedrals_registry[k][3]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][4]=molecule_dihedrals_registry[k][4]+atom_scaling_array[i];        // rescale
			}
            // alter total
			*general_dihedrals=*general_dihedrals+tcv;
            //
			if(verb==1)printf("$ general dihedrals registry:\n");
			if(verb==1)for(j=0;j<*general_dihedrals;++j)printf("[%d]\t%d\t%d\t%d\t%d\t%d\n",j+1,(*general_dihedrals_registry)[j][0],(*general_dihedrals_registry)[j][1],(*general_dihedrals_registry)[j][2],(*general_dihedrals_registry)[j][3],(*general_dihedrals_registry)[j][4]);
		}
		else if(non_minus_current==1)
		{
            // delta holds the number of dihedrals prior to the building step
            delta=topo_boundaries[i][5]-topo_boundaries[i][4]+1;
            // in_array holds the new values for the current molecule due to the building step
            in_array[0]=0;in_array[1]=molecule_dihedrals-1;
            if(verb==1)printf("$ case 1 in_array = [%d %d]\n",in_array[0],in_array[1]);
            // update current molecule values
            topo_boundaries[i][4]=in_array[0];
            topo_boundaries[i][5]=in_array[1];
            // tcv (topo current value) holds the difference in elements (delta_dihedrals)
            tcv=topo_boundaries[i][5]-topo_boundaries[i][4]+1-delta;
            if(verb==1)printf("$ case 1 tcv = %d\n",tcv);
            //
            if(tcv>0)
            {
                // resize adding extra space
                if(verb==1)printf("$ case 1 init / final sizes: %d / %d\n",*general_dihedrals,*general_dihedrals+tcv);
                *general_dihedrals_registry=(int**)realloc(*general_dihedrals_registry,(*general_dihedrals+tcv)*sizeof(int*));
                for(j=0;j<tcv;++j)(*general_dihedrals_registry)[*general_dihedrals+j]=(int*)malloc(5*sizeof(int));
                // place current molecule entries
				k=-1;
				for(j=topo_boundaries[i][4];j<=topo_boundaries[i][5];++j)
				{
					k=k+1;
					if(verb==1)printf("$ placing %d at %d...\n",k,j);
					(*general_dihedrals_registry)[j][0]=(*dtmapping_matrix)[molecule_dihedrals_registry[k][0]-1][1];    // type mapping
					(*general_dihedrals_registry)[j][1]=molecule_dihedrals_registry[k][1]+atom_scaling_array[i];        // rescale
					(*general_dihedrals_registry)[j][2]=molecule_dihedrals_registry[k][2]+atom_scaling_array[i];        // rescale
					(*general_dihedrals_registry)[j][3]=molecule_dihedrals_registry[k][3]+atom_scaling_array[i];        // rescale
					(*general_dihedrals_registry)[j][4]=molecule_dihedrals_registry[k][4]+atom_scaling_array[i];        // rescale
				}
            }
            else if(tcv<0)
            {
                // resize removing space
                if(verb==1){
                printf("$ case 1 init / final sizes: %d / %d\n",*general_dihedrals,*general_dihedrals+tcv);
                printf("\n$ tcv<0\n\n");
                printf("\n$ topo_boundaries:\t%d\t%d\n",topo_boundaries[i][4],topo_boundaries[i][5]);}
                // place current molecule entries
				k=-1;
				for(j=topo_boundaries[i][4];j<=topo_boundaries[i][5];++j)
				{
					k=k+1;
					if(verb==1)printf("$ placing %d at %d...\n",k,j);
					(*general_dihedrals_registry)[j][0]=(*dtmapping_matrix)[molecule_dihedrals_registry[k][0]-1][1];    // type mapping
					(*general_dihedrals_registry)[j][1]=molecule_dihedrals_registry[k][1]+atom_scaling_array[i];        // rescale
					(*general_dihedrals_registry)[j][2]=molecule_dihedrals_registry[k][2]+atom_scaling_array[i];        // rescale
					(*general_dihedrals_registry)[j][3]=molecule_dihedrals_registry[k][3]+atom_scaling_array[i];        // rescale
					(*general_dihedrals_registry)[j][4]=molecule_dihedrals_registry[k][4]+atom_scaling_array[i];        // rescale
				}
                // free extra memory
                for(j=(*general_dihedrals-1);j>(*general_dihedrals+tcv-1);--j)
                {
					free((*general_dihedrals_registry)[j]);
				}
            }
            else
            {
                // size kept the same
                if(verb==1)printf("$ case 1 init / final sizes: %d / %d\n",*general_dihedrals,*general_dihedrals+tcv);
                // place current molecule entries
				k=-1;
				for(j=topo_boundaries[i][4];j<=topo_boundaries[i][5];++j)
				{
					k=k+1;
					if(verb==1)printf("$ placing %d at %d...\n",k,j);
					(*general_dihedrals_registry)[j][0]=(*dtmapping_matrix)[molecule_dihedrals_registry[k][0]-1][1];    // type mapping
					(*general_dihedrals_registry)[j][1]=molecule_dihedrals_registry[k][1]+atom_scaling_array[i];        // rescale
					(*general_dihedrals_registry)[j][2]=molecule_dihedrals_registry[k][2]+atom_scaling_array[i];        // rescale
					(*general_dihedrals_registry)[j][3]=molecule_dihedrals_registry[k][3]+atom_scaling_array[i];        // rescale
					(*general_dihedrals_registry)[j][4]=molecule_dihedrals_registry[k][4]+atom_scaling_array[i];        // rescale
				}
            }
            /*
            // place current molecule entries
            k=-1;
            for(j=topo_boundaries[i][4];j<=topo_boundaries[i][5];++j)
            {
                k=k+1;
                printf("$ placing %d at %d...\n",k,j);
                (*general_dihedrals_registry)[j][0]=(*dtmapping_matrix)[molecule_dihedrals_registry[k][0]-1][1];    // type mapping
                (*general_dihedrals_registry)[j][1]=molecule_dihedrals_registry[k][1]+atom_scaling_array[i];        // rescale
                (*general_dihedrals_registry)[j][2]=molecule_dihedrals_registry[k][2]+atom_scaling_array[i];        // rescale
                (*general_dihedrals_registry)[j][3]=molecule_dihedrals_registry[k][3]+atom_scaling_array[i];        // rescale
                (*general_dihedrals_registry)[j][4]=molecule_dihedrals_registry[k][4]+atom_scaling_array[i];        // rescale
            }
            */
            // alter total
            *general_dihedrals=*general_dihedrals+tcv;
            //
            if(verb==1)printf("$ general dihedrals registry:\n");
            if(verb==1)for(j=0;j<*general_dihedrals;++j)printf("[%d]\t%d\t%d\t%d\t%d\t%d\n",j+1,(*general_dihedrals_registry)[j][0],(*general_dihedrals_registry)[j][1],(*general_dihedrals_registry)[j][2],(*general_dihedrals_registry)[j][3],(*general_dihedrals_registry)[j][4]);
		}
	}
	
	// case 2
	else if(non_minus_above==1 && non_minus_below==0)
	{
        // delta holds the number of dihedrals prior to the building step
		delta=0;if(non_minus_current==1)delta=topo_boundaries[i][5]-topo_boundaries[i][4]+1;
		if(verb==1)printf("$ case 2 delta = %d\n",delta);
        // tav (topo above value) holds the largest eligible position value above the current molecule
		tav=topo_boundaries[0][5];for(j=0;j<=i-1;++j)if(topo_boundaries[j][5]>tav)tav=topo_boundaries[j][5];
		if(verb==1)printf("$ case 2 tav = %d\n",tav);
        // in_array holds the new values for the current molecule due to the building step
		in_array[0]=0+tav+1;in_array[1]=molecule_dihedrals-1+tav+1;
		if(verb==1)printf("$ case 2 in_array = [%d %d]\n",in_array[0],in_array[1]);
        // update current molecule values
		topo_boundaries[i][4]=in_array[0];
		topo_boundaries[i][5]=in_array[1];
        // tcv (topo current value) holds the difference in elements (delta_dihedrals)
		tcv=topo_boundaries[i][5]-topo_boundaries[i][4]+1-delta;
		if(verb==1){printf("$ case 2 tcv = %d\n",tcv);
		printf("$ case 2 updated topo_boundaries:\n");
		for(j=0;j<molecules;++j)printf("[%d]\t%d\t%d\n",j+1,topo_boundaries[j][4],topo_boundaries[j][5]);}
        //
		if(tcv>0)
		{
            // resize adding extra space
			if(verb==1)printf("$ case 2 init / final sizes: %d / %d\n",*general_dihedrals,*general_dihedrals+tcv);
			*general_dihedrals_registry=(int**)realloc(*general_dihedrals_registry,(*general_dihedrals+tcv)*sizeof(int*));
			for(j=0;j<tcv;++j)(*general_dihedrals_registry)[*general_dihedrals+j]=(int*)malloc(5*sizeof(int));
			// place current molecule entries
			k=-1;
			for(j=topo_boundaries[i][4];j<=topo_boundaries[i][5];++j)
			{
				k=k+1;
				if(verb==1)printf("$ placing %d at %d...\n",k,j);
				(*general_dihedrals_registry)[j][0]=(*dtmapping_matrix)[molecule_dihedrals_registry[k][0]-1][1];    // type mapping
				(*general_dihedrals_registry)[j][1]=molecule_dihedrals_registry[k][1]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][2]=molecule_dihedrals_registry[k][2]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][3]=molecule_dihedrals_registry[k][3]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][4]=molecule_dihedrals_registry[k][4]+atom_scaling_array[i];        // rescale
			}
		}
		else if (tcv<0)
		{
            // resize removing space
            if(verb==1){
			printf("$ case 2 init / final sizes: %d / %d\n",*general_dihedrals,*general_dihedrals+tcv);
			printf("\n$ tcv<0\n\n");
            printf("\n$ topo_boundaries:\t%d\t%d\n",topo_boundaries[i][4],topo_boundaries[i][5]);}
			// place current molecule entries
			k=-1;
			for(j=topo_boundaries[i][4];j<=topo_boundaries[i][5];++j)
			{
				k=k+1;
				if(verb==1)printf("$ placing %d at %d...\n",k,j);
				(*general_dihedrals_registry)[j][0]=(*dtmapping_matrix)[molecule_dihedrals_registry[k][0]-1][1];    // type mapping
				(*general_dihedrals_registry)[j][1]=molecule_dihedrals_registry[k][1]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][2]=molecule_dihedrals_registry[k][2]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][3]=molecule_dihedrals_registry[k][3]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][4]=molecule_dihedrals_registry[k][4]+atom_scaling_array[i];        // rescale
			}
			// free extra memory
			for(j=(*general_dihedrals-1);j>(*general_dihedrals+tcv-1);--j)
			{
				free((*general_dihedrals_registry)[j]);
			}
		}
		else
		{
            // size kept the same
			if(verb==1)printf("$ case 2 init / final sizes: %d / %d\n",*general_dihedrals,*general_dihedrals+tcv);
			// place current molecule entries
			k=-1;
			for(j=topo_boundaries[i][4];j<=topo_boundaries[i][5];++j)
			{
				k=k+1;
				if(verb==1)printf("$ placing %d at %d...\n",k,j);
				(*general_dihedrals_registry)[j][0]=(*dtmapping_matrix)[molecule_dihedrals_registry[k][0]-1][1];    // type mapping
				(*general_dihedrals_registry)[j][1]=molecule_dihedrals_registry[k][1]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][2]=molecule_dihedrals_registry[k][2]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][3]=molecule_dihedrals_registry[k][3]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][4]=molecule_dihedrals_registry[k][4]+atom_scaling_array[i];        // rescale
			}
		}
		/*
        // place current molecule entries
		k=-1;
		for(j=topo_boundaries[i][4];j<=topo_boundaries[i][5];++j)
		{
			k=k+1;
			printf("$ placing %d at %d...\n",k,j);
			(*general_dihedrals_registry)[j][0]=(*dtmapping_matrix)[molecule_dihedrals_registry[k][0]-1][1];    // type mapping
			(*general_dihedrals_registry)[j][1]=molecule_dihedrals_registry[k][1]+atom_scaling_array[i];        // rescale
			(*general_dihedrals_registry)[j][2]=molecule_dihedrals_registry[k][2]+atom_scaling_array[i];        // rescale
			(*general_dihedrals_registry)[j][3]=molecule_dihedrals_registry[k][3]+atom_scaling_array[i];        // rescale
			(*general_dihedrals_registry)[j][4]=molecule_dihedrals_registry[k][4]+atom_scaling_array[i];        // rescale
		}
		*/
        // alter total
		*general_dihedrals=*general_dihedrals+tcv;
        //
		if(verb==1)printf("$ general dihedrals registry:\n");
		if(verb==1)for(j=0;j<*general_dihedrals;++j)printf("[%d]\t%d\t%d\t%d\t%d\t%d\n",j+1,(*general_dihedrals_registry)[j][0],(*general_dihedrals_registry)[j][1],(*general_dihedrals_registry)[j][2],(*general_dihedrals_registry)[j][3],(*general_dihedrals_registry)[j][4]);		
	}
	
	// case 3
	else if(non_minus_above==0 && non_minus_below==1)
	{
        // delta holds the number of dihedrals prior to the building step
		delta=0;if(non_minus_current==1)delta=topo_boundaries[i][5]-topo_boundaries[i][4]+1;
		if(verb==1)printf("$ case 3 delta = %d\n",delta);
        // upper_0 and lower_0 define the boundaries for the values below
		for(j=i+1;j<molecules;++j)if(topo_boundaries[j][4]!=-1){upper_0=topo_boundaries[j][4];break;}
		lower_0=topo_boundaries[i][5];
		for(j=i+1;j<molecules;++j)if(topo_boundaries[j][5]>lower_0){lower_0=topo_boundaries[j][5];}
		if(verb==1){printf("$ case 3 upper_0 = %d\n",upper_0);
		printf("$ case 3 lower_0 = %d\n",lower_0);}
        // in_array holds the new values for the current molecule due to the building step
		in_array[0]=0;in_array[1]=molecule_dihedrals-1;
		if(verb==1)printf("$ case 3 in_array = [%d %d]\n",in_array[0],in_array[1]);
        // update current molecule values
		topo_boundaries[i][4]=in_array[0];
		topo_boundaries[i][5]=in_array[1];
        // tcv (topo current value) holds the difference in elements (delta_dihedrals)
		tcv=topo_boundaries[i][5]-topo_boundaries[i][4]+1-delta;
		if(verb==1)printf("$ case 3 tcv = %d\n",tcv);
        // update non -1 values below current molecule
		for(j=i+1;j<molecules;++j){if(topo_boundaries[j][4]!=-1){topo_boundaries[j][4]=topo_boundaries[j][4]+tcv;topo_boundaries[j][5]=topo_boundaries[j][5]+tcv;}}
        //
		if(verb==1){printf("$ case 3 updated topo_boundaries:\n");
		for(j=0;j<molecules;++j)printf("[%d]\t%d\t%d\n",j+1,topo_boundaries[j][4],topo_boundaries[j][5]);}
        //
		if(tcv>0)
		{
            // resize adding extra space
			if(verb==1)printf("$ case 3 init / final sizes: %d / %d\n",*general_dihedrals,*general_dihedrals+tcv);
			*general_dihedrals_registry=(int**)realloc(*general_dihedrals_registry,(*general_dihedrals+tcv)*sizeof(int*));
			for(j=0;j<tcv;++j)(*general_dihedrals_registry)[*general_dihedrals+j]=(int*)malloc(5*sizeof(int));
            // push down entries
			if(verb==1)printf("$ case 3 push down: %d:%d --> %d:%d\n",upper_0,lower_0,upper_0+tcv,lower_0+tcv);
			for(j=lower_0;j>=upper_0;--j)
			{
				if(verb==1)printf("%d --> %d\n",j,j+tcv);
				(*general_dihedrals_registry)[j+tcv][0]=(*general_dihedrals_registry)[j][0];
				(*general_dihedrals_registry)[j+tcv][1]=(*general_dihedrals_registry)[j][1]+delta_atoms;    // rescale
				(*general_dihedrals_registry)[j+tcv][2]=(*general_dihedrals_registry)[j][2]+delta_atoms;    // rescale
				(*general_dihedrals_registry)[j+tcv][3]=(*general_dihedrals_registry)[j][3]+delta_atoms;    // rescale
				(*general_dihedrals_registry)[j+tcv][4]=(*general_dihedrals_registry)[j][4]+delta_atoms;    // rescale
			}
			// place current molecule entries
			k=-1;
			for(j=topo_boundaries[i][4];j<=topo_boundaries[i][5];++j)
			{
				k=k+1;
				if(verb==1)printf("$ placing %d at %d...\n",k,j);
				(*general_dihedrals_registry)[j][0]=(*dtmapping_matrix)[molecule_dihedrals_registry[k][0]-1][1];    // type mapping
				(*general_dihedrals_registry)[j][1]=molecule_dihedrals_registry[k][1]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][2]=molecule_dihedrals_registry[k][2]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][3]=molecule_dihedrals_registry[k][3]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][4]=molecule_dihedrals_registry[k][4]+atom_scaling_array[i];        // rescale
			}
		}
		else if (tcv<0)
		{
            // resize removing space
            if(verb==1){
			printf("$ case 3 init / final sizes: %d / %d\n",*general_dihedrals,*general_dihedrals+tcv);
			printf("\n$ tcv<0\n\n");
            printf("\n$ topo_boundaries:\t%d\t%d\n",topo_boundaries[i][4],topo_boundaries[i][5]);}
            // pull up!!
            for(j=upper_0;j<=lower_0;++j)
            {
				//printf("%d --> %d\n",j,j+tcv);
				(*general_dihedrals_registry)[j+tcv][0]=(*general_dihedrals_registry)[j][0];
				(*general_dihedrals_registry)[j+tcv][1]=(*general_dihedrals_registry)[j][1]+delta_atoms;
				(*general_dihedrals_registry)[j+tcv][2]=(*general_dihedrals_registry)[j][2]+delta_atoms;
				(*general_dihedrals_registry)[j+tcv][3]=(*general_dihedrals_registry)[j][3]+delta_atoms;
				(*general_dihedrals_registry)[j+tcv][4]=(*general_dihedrals_registry)[j][4]+delta_atoms;
			}
            // place current molecule entries
			k=-1;
			for(j=topo_boundaries[i][4];j<=topo_boundaries[i][5];++j)
			{
				k=k+1;
				if(verb==1)printf("$ placing %d at %d...\n",k,j);
				(*general_dihedrals_registry)[j][0]=(*dtmapping_matrix)[molecule_dihedrals_registry[k][0]-1][1];    // type mapping
				(*general_dihedrals_registry)[j][1]=molecule_dihedrals_registry[k][1]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][2]=molecule_dihedrals_registry[k][2]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][3]=molecule_dihedrals_registry[k][3]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][4]=molecule_dihedrals_registry[k][4]+atom_scaling_array[i];        // rescale
			}
            // free extra memory
			for(j=(*general_dihedrals-1);j>(*general_dihedrals+tcv-1);--j)
			{
				free((*general_dihedrals_registry)[j]);
			}
		}
		else
		{
            // size kept the same
			if(verb==1){printf("$ case 3 init / final sizes: %d / %d\n",*general_dihedrals,*general_dihedrals+tcv);
            // update for rescaling purposes
			printf("$ case 3 no shift: %d:%d <--> %d:%d\n",upper_0,lower_0,upper_0+tcv,lower_0+tcv);}
			for(j=lower_0;j>=upper_0;--j)
			{
				if(verb==1)printf("%d --> %d\n",j,j+tcv);
				(*general_dihedrals_registry)[j+tcv][0]=(*general_dihedrals_registry)[j][0];
				(*general_dihedrals_registry)[j+tcv][1]=(*general_dihedrals_registry)[j][1]+delta_atoms;    // rescale
				(*general_dihedrals_registry)[j+tcv][2]=(*general_dihedrals_registry)[j][2]+delta_atoms;    // rescale
				(*general_dihedrals_registry)[j+tcv][3]=(*general_dihedrals_registry)[j][3]+delta_atoms;    // rescale
				(*general_dihedrals_registry)[j+tcv][4]=(*general_dihedrals_registry)[j][4]+delta_atoms;    // rescale
			}
			// place current molecule entries
			k=-1;
			for(j=topo_boundaries[i][4];j<=topo_boundaries[i][5];++j)
			{
				k=k+1;
				if(verb==1)printf("$ placing %d at %d...\n",k,j);
				(*general_dihedrals_registry)[j][0]=(*dtmapping_matrix)[molecule_dihedrals_registry[k][0]-1][1];    // type mapping
				(*general_dihedrals_registry)[j][1]=molecule_dihedrals_registry[k][1]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][2]=molecule_dihedrals_registry[k][2]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][3]=molecule_dihedrals_registry[k][3]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][4]=molecule_dihedrals_registry[k][4]+atom_scaling_array[i];        // rescale
			}
		}
		/*
        // place current molecule entries
		k=-1;
		for(j=topo_boundaries[i][4];j<=topo_boundaries[i][5];++j)
		{
			k=k+1;
			printf("$ placing %d at %d...\n",k,j);
			(*general_dihedrals_registry)[j][0]=(*dtmapping_matrix)[molecule_dihedrals_registry[k][0]-1][1];    // type mapping
			(*general_dihedrals_registry)[j][1]=molecule_dihedrals_registry[k][1]+atom_scaling_array[i];        // rescale
			(*general_dihedrals_registry)[j][2]=molecule_dihedrals_registry[k][2]+atom_scaling_array[i];        // rescale
			(*general_dihedrals_registry)[j][3]=molecule_dihedrals_registry[k][3]+atom_scaling_array[i];        // rescale
			(*general_dihedrals_registry)[j][4]=molecule_dihedrals_registry[k][4]+atom_scaling_array[i];        // rescale
		}
		*/ 
        // alter total
		*general_dihedrals=*general_dihedrals+tcv;
        //
		if(verb==1)printf("$ general dihedrals registry:\n");
		if(verb==1)for(j=0;j<*general_dihedrals;++j)printf("[%d]\t%d\t%d\t%d\t%d\t%d\n",j+1,(*general_dihedrals_registry)[j][0],(*general_dihedrals_registry)[j][1],(*general_dihedrals_registry)[j][2],(*general_dihedrals_registry)[j][3],(*general_dihedrals_registry)[j][4]);
	}
	
	// case 4
	else if(non_minus_above==1 && non_minus_below==1)
	{
        // delta holds the number of dihedrals prior to the building step
		delta=0;if(non_minus_current==1)delta=topo_boundaries[i][5]-topo_boundaries[i][4]+1;
		if(verb==1)printf("$ case 4 delta = %d\n",delta);
        // upper_0 and lower_0 define the boundaries for the values below
		for(j=i+1;j<molecules;++j)if(topo_boundaries[j][4]!=-1){upper_0=topo_boundaries[j][4];break;}
		lower_0=topo_boundaries[i][5];
		for(j=i+1;j<molecules;++j)if(topo_boundaries[j][5]>lower_0){lower_0=topo_boundaries[j][5];}
		if(verb==1){printf("$ case 4 upper_0 = %d\n",upper_0);
		printf("$ case 4 lower_0 = %d\n",lower_0);}
        // tav (topo above value) holds the largest eligible position value above the current molecule
		tav=topo_boundaries[0][5];for(j=0;j<=i-1;++j)if(topo_boundaries[j][5]>tav)tav=topo_boundaries[j][5];
		if(verb==1)printf("$ case 4 tav = %d\n",tav);
        // in_array holds the new values for the current molecule due to the building step
        in_array[0]=0+tav+1;in_array[1]=molecule_dihedrals-1+tav+1;
        if(verb==1)printf("$ case 4 in_array = [%d %d]\n",in_array[0],in_array[1]);
        // update current molecule values
		topo_boundaries[i][4]=in_array[0];
		topo_boundaries[i][5]=in_array[1];
        // tcv (topo current value) holds the difference in elements (delta_dihedrals)
		tcv=topo_boundaries[i][5]-topo_boundaries[i][4]+1-delta;
		if(verb==1)printf("$ case 4 tcv = %d\n",tcv);
        // update non -1 values below current molecule
		for(j=i+1;j<molecules;++j){if(topo_boundaries[j][4]!=-1){topo_boundaries[j][4]=topo_boundaries[j][4]+tcv;topo_boundaries[j][5]=topo_boundaries[j][5]+tcv;}}
        //
		if(verb==1){printf("$ case 4 updated topo_boundaries:\n");
		for(j=0;j<molecules;++j)printf("[%d]\t%d\t%d\n",j+1,topo_boundaries[j][4],topo_boundaries[j][5]);}
        //
		if(tcv>0)
		{
            // resize adding extra space
			if(verb==1)printf("$ case 4 init / final sizes: %d / %d\n",*general_dihedrals,*general_dihedrals+tcv);
			*general_dihedrals_registry=(int**)realloc(*general_dihedrals_registry,(*general_dihedrals+tcv)*sizeof(int*));
			for(j=0;j<tcv;++j)(*general_dihedrals_registry)[*general_dihedrals+j]=(int*)malloc(5*sizeof(int));
            // push down entries
			if(verb==1)printf("$ case 4 push down: %d:%d --> %d:%d\n",upper_0,lower_0,upper_0+tcv,lower_0+tcv);
			for(j=lower_0;j>=upper_0;--j)
			{
				if(verb==1)printf("%d --> %d\n",j,j+tcv);
				(*general_dihedrals_registry)[j+tcv][0]=(*general_dihedrals_registry)[j][0];
				(*general_dihedrals_registry)[j+tcv][1]=(*general_dihedrals_registry)[j][1]+delta_atoms;    // rescale
				(*general_dihedrals_registry)[j+tcv][2]=(*general_dihedrals_registry)[j][2]+delta_atoms;    // rescale
				(*general_dihedrals_registry)[j+tcv][3]=(*general_dihedrals_registry)[j][3]+delta_atoms;    // rescale
				(*general_dihedrals_registry)[j+tcv][4]=(*general_dihedrals_registry)[j][4]+delta_atoms;    // rescale
			}
			// place current molecule entries
			k=-1;
			for(j=topo_boundaries[i][4];j<=topo_boundaries[i][5];++j)
			{
				k=k+1;
				if(verb==1)printf("$ placing %d at %d...\n",k,j);
				(*general_dihedrals_registry)[j][0]=(*dtmapping_matrix)[molecule_dihedrals_registry[k][0]-1][1];    // type mapping
				(*general_dihedrals_registry)[j][1]=molecule_dihedrals_registry[k][1]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][2]=molecule_dihedrals_registry[k][2]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][3]=molecule_dihedrals_registry[k][3]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][4]=molecule_dihedrals_registry[k][4]+atom_scaling_array[i];        // rescale
			}
		}
		else if (tcv<0)
		{
            // resize removing space
            if(verb==1){
			printf("$ case 4 init / final sizes: %d / %d\n",*general_dihedrals,*general_dihedrals+tcv);
			printf("\n$ tcv<0\n\n");}
			// pull up!!
            for(j=upper_0;j<=lower_0;++j)
            {
				//printf("%d --> %d\n",j,j+tcv);
				(*general_dihedrals_registry)[j+tcv][0]=(*general_dihedrals_registry)[j][0];
				(*general_dihedrals_registry)[j+tcv][1]=(*general_dihedrals_registry)[j][1]+delta_atoms;
				(*general_dihedrals_registry)[j+tcv][2]=(*general_dihedrals_registry)[j][2]+delta_atoms;
				(*general_dihedrals_registry)[j+tcv][3]=(*general_dihedrals_registry)[j][3]+delta_atoms;
				(*general_dihedrals_registry)[j+tcv][4]=(*general_dihedrals_registry)[j][4]+delta_atoms;
			}
            // place current molecule entries
			k=-1;
			for(j=topo_boundaries[i][4];j<=topo_boundaries[i][5];++j)
			{
				k=k+1;
				if(verb==1)printf("$ placing %d at %d...\n",k,j);
				(*general_dihedrals_registry)[j][0]=(*dtmapping_matrix)[molecule_dihedrals_registry[k][0]-1][1];    // type mapping
				(*general_dihedrals_registry)[j][1]=molecule_dihedrals_registry[k][1]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][2]=molecule_dihedrals_registry[k][2]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][3]=molecule_dihedrals_registry[k][3]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][4]=molecule_dihedrals_registry[k][4]+atom_scaling_array[i];        // rescale
			}
			// free extra memory
			for(j=(*general_dihedrals-1);j>(*general_dihedrals+tcv-1);--j)
			{
				free((*general_dihedrals_registry)[j]);
			}
		}
		else
		{
            // size kept the same
            if(verb==1){
			printf("$ case 4 init / final sizes: %d / %d\n",*general_dihedrals,*general_dihedrals+tcv);
			printf("$ case 4 no shift: %d:%d <--> %d:%d\n",upper_0,lower_0,upper_0+tcv,lower_0+tcv);}
            // update for rescaling purposes
			for(j=lower_0;j>=upper_0;--j)
			{
				if(verb==1)printf("%d --> %d\n",j,j+tcv);
				(*general_dihedrals_registry)[j+tcv][0]=(*general_dihedrals_registry)[j][0];
				(*general_dihedrals_registry)[j+tcv][1]=(*general_dihedrals_registry)[j][1]+delta_atoms;    // rescale
				(*general_dihedrals_registry)[j+tcv][2]=(*general_dihedrals_registry)[j][2]+delta_atoms;    // rescale
				(*general_dihedrals_registry)[j+tcv][3]=(*general_dihedrals_registry)[j][3]+delta_atoms;    // rescale
				(*general_dihedrals_registry)[j+tcv][4]=(*general_dihedrals_registry)[j][4]+delta_atoms;    // rescale
			}
			// place current molecule entries
			k=-1;
			for(j=topo_boundaries[i][4];j<=topo_boundaries[i][5];++j)
			{
				k=k+1;
				if(verb==1)printf("$ placing %d at %d...\n",k,j);
				(*general_dihedrals_registry)[j][0]=(*dtmapping_matrix)[molecule_dihedrals_registry[k][0]-1][1];    // type mapping
				(*general_dihedrals_registry)[j][1]=molecule_dihedrals_registry[k][1]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][2]=molecule_dihedrals_registry[k][2]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][3]=molecule_dihedrals_registry[k][3]+atom_scaling_array[i];        // rescale
				(*general_dihedrals_registry)[j][4]=molecule_dihedrals_registry[k][4]+atom_scaling_array[i];        // rescale
			}
		}
		/*
        // place current molecule entries
		k=-1;
		for(j=topo_boundaries[i][4];j<=topo_boundaries[i][5];++j)
		{
			k=k+1;
			printf("$ placing %d at %d...\n",k,j);
			(*general_dihedrals_registry)[j][0]=(*dtmapping_matrix)[molecule_dihedrals_registry[k][0]-1][1];    // type mapping
			(*general_dihedrals_registry)[j][1]=molecule_dihedrals_registry[k][1]+atom_scaling_array[i];        // rescale
			(*general_dihedrals_registry)[j][2]=molecule_dihedrals_registry[k][2]+atom_scaling_array[i];        // rescale
			(*general_dihedrals_registry)[j][3]=molecule_dihedrals_registry[k][3]+atom_scaling_array[i];        // rescale
			(*general_dihedrals_registry)[j][4]=molecule_dihedrals_registry[k][4]+atom_scaling_array[i];        // rescale
		}
		*/
        // alter total
		*general_dihedrals=*general_dihedrals+tcv;
        //
		if(verb==1)printf("$ general dihedrals registry:\n");
		if(verb==1)for(j=0;j<*general_dihedrals;++j)printf("[%d]\t%d\t%d\t%d\t%d\t%d\n",j+1,(*general_dihedrals_registry)[j][0],(*general_dihedrals_registry)[j][1],(*general_dihedrals_registry)[j][2],(*general_dihedrals_registry)[j][3],(*general_dihedrals_registry)[j][4]);
	}
	
	else
	{
		printf("unsupported dihedrals combination!!\n\n");exit(-1);
	}
	
	delete_flag=0;
    for(k=0;k<*dtmap_rows;++k)
        if(dcm_sum[k]==0)
            delete_flag=delete_flag+1;
    
    if(delete_flag>0)
    {
        if(verb==1)printf("\n$ delete flag>0!!!\n\n");
        gtmapping_array=(int*)malloc((*general_dihedral_types)*sizeof(int));
        backup=(char***)malloc((*general_dihedral_types)*sizeof(char**));
        for(j=0;j<*general_dihedral_types;++j)backup[j]=(char**)malloc(7*sizeof(char*));
        for(j=0;j<*general_dihedral_types;++j)for(k=0;k<7;++k)backup[j][k]=(char*)malloc(sub_length*sizeof(char));
        counter=0;
        for(j=0;j<*general_dihedral_types;++j)
        {
            for(k=0;k<*dtmap_rows;++k)
            {
                if((strcmp((*general_dihedral_types_registry)[j][0],(*dtmap)[k][0])==0 &&
                    strcmp((*general_dihedral_types_registry)[j][1],(*dtmap)[k][1])==0 &&
                    strcmp((*general_dihedral_types_registry)[j][2],(*dtmap)[k][2])==0 &&
                    strcmp((*general_dihedral_types_registry)[j][3],(*dtmap)[k][3])==0 &&
                    strcmp((*general_dihedral_types_registry)[j][4],(*dtmap)[k][4])==0 &&
                    strcmp((*general_dihedral_types_registry)[j][5],(*dtmap)[k][5])==0 &&
                    strcmp((*general_dihedral_types_registry)[j][6],(*dtmap)[k][6])==0
                    )
                   ||
                   (strcmp((*general_dihedral_types_registry)[j][0],(*dtmap)[k][6])==0 &&
                    strcmp((*general_dihedral_types_registry)[j][1],(*dtmap)[k][5])==0 &&
                    strcmp((*general_dihedral_types_registry)[j][2],(*dtmap)[k][4])==0 &&
                    strcmp((*general_dihedral_types_registry)[j][3],(*dtmap)[k][3])==0 &&
                    strcmp((*general_dihedral_types_registry)[j][4],(*dtmap)[k][2])==0 &&
                    strcmp((*general_dihedral_types_registry)[j][5],(*dtmap)[k][1])==0 &&
                    strcmp((*general_dihedral_types_registry)[j][6],(*dtmap)[k][0])==0
                    ))
                    break;
            }
            if(verb==1){printf("$ (%s)--%s--(%s)--%s--(%s)--%s--(%s) is at position {%d}\n",(*general_dihedral_types_registry)[j][0],(*general_dihedral_types_registry)[j][1],(*general_dihedral_types_registry)[j][2],(*general_dihedral_types_registry)[j][3],(*general_dihedral_types_registry)[j][4],(*general_dihedral_types_registry)[j][5],(*general_dihedral_types_registry)[j][6],k+1);
            printf("$ checking cm_sum: %d\n",dcm_sum[k]);}
            if(dcm_sum[k]==0)
            {
                gtmapping_array[j]=0;
            }
            else
            {
                counter=counter+1;
                gtmapping_array[j]=counter;
                sprintf(backup[counter-1][0],"%s",(*general_dihedral_types_registry)[j][0]);
                sprintf(backup[counter-1][1],"%s",(*general_dihedral_types_registry)[j][1]);
                sprintf(backup[counter-1][2],"%s",(*general_dihedral_types_registry)[j][2]);
                sprintf(backup[counter-1][3],"%s",(*general_dihedral_types_registry)[j][3]);
                sprintf(backup[counter-1][4],"%s",(*general_dihedral_types_registry)[j][4]);
                sprintf(backup[counter-1][5],"%s",(*general_dihedral_types_registry)[j][5]);
                sprintf(backup[counter-1][6],"%s",(*general_dihedral_types_registry)[j][6]);
            }
            
        }
        if(verb==1){
        printf("$ gtmapping_array:\n");
        for(j=0;j<*general_dihedral_types;++j)printf("[%d]\t%d\n",j+1,gtmapping_array[j]);
        printf("$ altered types:\n");
        for(j=0;j<counter;++j)printf("[%d]\t(%s)--%s--(%s)--%s--(%s)--%s--(%s)\n",j+1,backup[j][0],backup[j][1],backup[j][2],backup[j][3],backup[j][4],backup[j][5],backup[j][6]);}
        
        for(j=0;j<*general_dihedrals;++j)(*general_dihedrals_registry)[j][0]=gtmapping_array[(*general_dihedrals_registry)[j][0]-1];
        for(j=0;j<counter;++j)
        {
            sprintf((*general_dihedral_types_registry)[j][0],"%s",backup[j][0]);
            sprintf((*general_dihedral_types_registry)[j][1],"%s",backup[j][1]);
            sprintf((*general_dihedral_types_registry)[j][2],"%s",backup[j][2]);
            sprintf((*general_dihedral_types_registry)[j][3],"%s",backup[j][3]);
            sprintf((*general_dihedral_types_registry)[j][4],"%s",backup[j][4]);
            sprintf((*general_dihedral_types_registry)[j][5],"%s",backup[j][5]);
            sprintf((*general_dihedral_types_registry)[j][6],"%s",backup[j][6]);
        }
        for(j=counter;j<*general_dihedral_types;++j)
        {
            for(k=0;k<7;++k)
                free((*general_dihedral_types_registry)[j][k]);
            free((*general_dihedral_types_registry)[j]);
        }
        for(j=0;j<*general_dihedral_types;++j)
            for(k=0;k<7;++k)
                free(backup[j][k]);
        for(j=0;j<*general_dihedral_types;++j)free(backup[j]);
        free(backup);
        *general_dihedral_types=counter;
        free(gtmapping_array);
	}
    
    free(dcm_sum);
    
    //--------------------------------------------------------------------------
    
    // alter impropers
    
    //if(molecule_impropers>0){
    
    // [3]
    
    // update tmap and cm
    for(j=0;j<molecule_improper_types;++j)
    {
        found=0;
        for(k=0;k<*itmap_rows;++k)
        {
            if(
               (strcmp(molecule_improper_types_registry[j][0],(*itmap)[k][0])==0 &&
                strcmp(molecule_improper_types_registry[j][1],(*itmap)[k][1])==0 &&
                strcmp(molecule_improper_types_registry[j][2],(*itmap)[k][2])==0 &&
                strcmp(molecule_improper_types_registry[j][3],(*itmap)[k][3])==0
                )||
               (strcmp(molecule_improper_types_registry[j][0],(*itmap)[k][0])==0 &&
                strcmp(molecule_improper_types_registry[j][1],(*itmap)[k][1])==0 &&
                strcmp(molecule_improper_types_registry[j][2],(*itmap)[k][3])==0 &&
                strcmp(molecule_improper_types_registry[j][3],(*itmap)[k][2])==0
                )||
               (strcmp(molecule_improper_types_registry[j][0],(*itmap)[k][0])==0 &&
                strcmp(molecule_improper_types_registry[j][1],(*itmap)[k][2])==0 &&
                strcmp(molecule_improper_types_registry[j][2],(*itmap)[k][1])==0 &&
                strcmp(molecule_improper_types_registry[j][3],(*itmap)[k][3])==0
                )||
               (strcmp(molecule_improper_types_registry[j][0],(*itmap)[k][0])==0 &&
                strcmp(molecule_improper_types_registry[j][1],(*itmap)[k][2])==0 &&
                strcmp(molecule_improper_types_registry[j][2],(*itmap)[k][3])==0 &&
                strcmp(molecule_improper_types_registry[j][3],(*itmap)[k][1])==0
                )||
               (strcmp(molecule_improper_types_registry[j][0],(*itmap)[k][0])==0 &&
                strcmp(molecule_improper_types_registry[j][1],(*itmap)[k][3])==0 &&
                strcmp(molecule_improper_types_registry[j][2],(*itmap)[k][1])==0 &&
                strcmp(molecule_improper_types_registry[j][3],(*itmap)[k][2])==0
                )||
               (strcmp(molecule_improper_types_registry[j][0],(*itmap)[k][0])==0 &&
                strcmp(molecule_improper_types_registry[j][1],(*itmap)[k][3])==0 &&
                strcmp(molecule_improper_types_registry[j][2],(*itmap)[k][2])==0 &&
                strcmp(molecule_improper_types_registry[j][3],(*itmap)[k][1])==0
                )
               )
            {
                found=1;break;
            }
        }
        if(found==0)
        {
            // augment
            *itmap_rows=*itmap_rows+1;
            // tmap
            // resize
            *itmap=(char***)realloc(*itmap,(*itmap_rows)*sizeof(char**));
            (*itmap)[*itmap_rows-1]=(char**)malloc(4*sizeof(char*));
            (*itmap)[*itmap_rows-1][0]=(char*)malloc(sub_length*sizeof(char));
            (*itmap)[*itmap_rows-1][1]=(char*)malloc(sub_length*sizeof(char));
            (*itmap)[*itmap_rows-1][2]=(char*)malloc(sub_length*sizeof(char));
            (*itmap)[*itmap_rows-1][3]=(char*)malloc(sub_length*sizeof(char));
            // write
            // write type
            sprintf((*itmap)[*itmap_rows-1][0],"%s",molecule_improper_types_registry[j][0]);
            sprintf((*itmap)[*itmap_rows-1][1],"%s",molecule_improper_types_registry[j][1]);
            sprintf((*itmap)[*itmap_rows-1][2],"%s",molecule_improper_types_registry[j][2]);
            sprintf((*itmap)[*itmap_rows-1][3],"%s",molecule_improper_types_registry[j][3]);
            // cm
            for(k=0;k<molecules;++k)
            {
                (*icm)[k]=(int*)realloc((*icm)[k],(*itmap_rows)*sizeof(int*));
                (*icm)[k][*itmap_rows-1]=0;
            }
        }
    }
    
    // [4]
    
    // reset cm(i,:)
    for(j=0;j<*itmap_rows;++j)(*icm)[i][j]=0;
    // loop on mitr (molecule improper type registry) using itmap and update icm
    for(j=0;j<molecule_improper_types;++j)
    {
        for(k=0;k<*itmap_rows;++k)
            if((
                strcmp(molecule_improper_types_registry[j][0],(*itmap)[k][0])==0 &&
                strcmp(molecule_improper_types_registry[j][1],(*itmap)[k][1])==0 &&
                strcmp(molecule_improper_types_registry[j][2],(*itmap)[k][2])==0 &&
                strcmp(molecule_improper_types_registry[j][3],(*itmap)[k][3])==0
                )
               ||
               (
                strcmp(molecule_improper_types_registry[j][0],(*itmap)[k][0])==0 &&
                strcmp(molecule_improper_types_registry[j][1],(*itmap)[k][1])==0 &&
                strcmp(molecule_improper_types_registry[j][2],(*itmap)[k][3])==0 &&
                strcmp(molecule_improper_types_registry[j][3],(*itmap)[k][2])==0
                )
               ||
               (
                strcmp(molecule_improper_types_registry[j][0],(*itmap)[k][0])==0 &&
                strcmp(molecule_improper_types_registry[j][1],(*itmap)[k][2])==0 &&
                strcmp(molecule_improper_types_registry[j][2],(*itmap)[k][1])==0 &&
                strcmp(molecule_improper_types_registry[j][3],(*itmap)[k][3])==0
                )
               ||
               (
                strcmp(molecule_improper_types_registry[j][0],(*itmap)[k][0])==0 &&
                strcmp(molecule_improper_types_registry[j][1],(*itmap)[k][2])==0 &&
                strcmp(molecule_improper_types_registry[j][2],(*itmap)[k][3])==0 &&
                strcmp(molecule_improper_types_registry[j][3],(*itmap)[k][1])==0
                )
               ||
               (
                strcmp(molecule_improper_types_registry[j][0],(*itmap)[k][0])==0 &&
                strcmp(molecule_improper_types_registry[j][1],(*itmap)[k][3])==0 &&
                strcmp(molecule_improper_types_registry[j][2],(*itmap)[k][1])==0 &&
                strcmp(molecule_improper_types_registry[j][3],(*itmap)[k][2])==0
                )
               ||
               (
                strcmp(molecule_improper_types_registry[j][0],(*itmap)[k][0])==0 &&
                strcmp(molecule_improper_types_registry[j][1],(*itmap)[k][3])==0 &&
                strcmp(molecule_improper_types_registry[j][2],(*itmap)[k][2])==0 &&
                strcmp(molecule_improper_types_registry[j][3],(*itmap)[k][1])==0
                )
               
               )
                break;
        (*icm)[i][k]=1;
    }
    
     // console out
     if(verb==1){
     printf("\n$ itmap/icm:\n");
     for(j=0;j<*itmap_rows;++j)printf("%s-%s-%s-%s\t",(*itmap)[j][0],(*itmap)[j][1],(*itmap)[j][2],(*itmap)[j][3]);printf("\n");
     for(j=0;j<molecules;++j){
     for(k=0;k<*itmap_rows;++k)
     printf("%d\t",(*icm)[j][k]);
     printf("\n");
     }
     }
     
    // [5]
     
    // global improper types through itmap/icm:
    icm_sum=(int*)malloc((*itmap_rows)*sizeof(int));
    for(k=0;k<*itmap_rows;++k)
    {
        icm_sum[k]=0;
        for(j=0;j<molecules;++j)
        {
            icm_sum[k]=icm_sum[k]+(*icm)[j][k];
        }
    }
    //
    if(verb==1){
    printf("\n$ icm_sum array:\n");for(j=0;j<*itmap_rows;++j)printf("%d\t",icm_sum[j]);printf("\n");}
    
    // [6]
    *itmapping_matrix=(int**)malloc(molecule_improper_types*sizeof(int*));
	for(j=0;j<molecule_improper_types;++j)(*itmapping_matrix)[j]=(int*)malloc(2*sizeof(int));
    if(*general_impropers==0 && molecule_impropers>0){*general_improper_types_registry=(char***)malloc(0);}
    n=-1;
	for(j=0;j<molecule_improper_types;++j)
	{
		found=0;
		for(k=0;k<*general_improper_types;++k)
		{
			if(
			   (strcmp(molecule_improper_types_registry[j][0],(*general_improper_types_registry)[k][0])==0 &&
				strcmp(molecule_improper_types_registry[j][1],(*general_improper_types_registry)[k][1])==0 &&
				strcmp(molecule_improper_types_registry[j][2],(*general_improper_types_registry)[k][2])==0 &&
				strcmp(molecule_improper_types_registry[j][3],(*general_improper_types_registry)[k][3])==0
				)
			   ||
			   (strcmp(molecule_improper_types_registry[j][0],(*general_improper_types_registry)[k][0])==0 &&
				strcmp(molecule_improper_types_registry[j][1],(*general_improper_types_registry)[k][1])==0 &&
				strcmp(molecule_improper_types_registry[j][2],(*general_improper_types_registry)[k][3])==0 &&
				strcmp(molecule_improper_types_registry[j][3],(*general_improper_types_registry)[k][2])==0
				)
			   ||
			   (strcmp(molecule_improper_types_registry[j][0],(*general_improper_types_registry)[k][0])==0 &&
				strcmp(molecule_improper_types_registry[j][1],(*general_improper_types_registry)[k][2])==0 &&
				strcmp(molecule_improper_types_registry[j][2],(*general_improper_types_registry)[k][1])==0 &&
				strcmp(molecule_improper_types_registry[j][3],(*general_improper_types_registry)[k][3])==0
				)
			   ||
			   (strcmp(molecule_improper_types_registry[j][0],(*general_improper_types_registry)[k][0])==0 &&
				strcmp(molecule_improper_types_registry[j][1],(*general_improper_types_registry)[k][2])==0 &&
				strcmp(molecule_improper_types_registry[j][2],(*general_improper_types_registry)[k][3])==0 &&
				strcmp(molecule_improper_types_registry[j][3],(*general_improper_types_registry)[k][1])==0
				)
			   ||
			   (strcmp(molecule_improper_types_registry[j][0],(*general_improper_types_registry)[k][0])==0 &&
				strcmp(molecule_improper_types_registry[j][1],(*general_improper_types_registry)[k][3])==0 &&
				strcmp(molecule_improper_types_registry[j][2],(*general_improper_types_registry)[k][1])==0 &&
				strcmp(molecule_improper_types_registry[j][3],(*general_improper_types_registry)[k][2])==0
				)
			   ||
			   (strcmp(molecule_improper_types_registry[j][0],(*general_improper_types_registry)[k][0])==0 &&
				strcmp(molecule_improper_types_registry[j][1],(*general_improper_types_registry)[k][3])==0 &&
				strcmp(molecule_improper_types_registry[j][2],(*general_improper_types_registry)[k][2])==0 &&
				strcmp(molecule_improper_types_registry[j][3],(*general_improper_types_registry)[k][1])==0
				)
			   
			   )
			{
				n=n+1;
				//printf("* molecular %d is found at general %d\n",j+1,k+1);
				(*itmapping_matrix)[n][0]=j+1;
				(*itmapping_matrix)[n][1]=k+1;
				found=1;
				break;
			}
		}
		if(found==0)
		{
			// augment
			*general_improper_types=*general_improper_types+1;
			// resize
			*general_improper_types_registry=(char***)realloc(*general_improper_types_registry,(*general_improper_types)*sizeof(char**));
			(*general_improper_types_registry)[*general_improper_types-1]=(char**)malloc(4*sizeof(char*));
			(*general_improper_types_registry)[*general_improper_types-1][0]=(char*)malloc(sub_length*sizeof(char));
			(*general_improper_types_registry)[*general_improper_types-1][1]=(char*)malloc(sub_length*sizeof(char));
			(*general_improper_types_registry)[*general_improper_types-1][2]=(char*)malloc(sub_length*sizeof(char));
			(*general_improper_types_registry)[*general_improper_types-1][3]=(char*)malloc(sub_length*sizeof(char));

			// write
			// write type
			sprintf((*general_improper_types_registry)[*general_improper_types-1][0],"%s",molecule_improper_types_registry[j][0]);
			sprintf((*general_improper_types_registry)[*general_improper_types-1][1],"%s",molecule_improper_types_registry[j][1]);
			sprintf((*general_improper_types_registry)[*general_improper_types-1][2],"%s",molecule_improper_types_registry[j][2]);
			sprintf((*general_improper_types_registry)[*general_improper_types-1][3],"%s",molecule_improper_types_registry[j][3]);
			//
			n=n+1;
			//printf("~ molecular %d is placed at general %d\n",j+1,*general_improper_types);
			(*itmapping_matrix)[n][0]=j+1;
			(*itmapping_matrix)[n][1]=*general_improper_types;
		}
		
	}
	if(verb==1){
	printf("\n$ itmapping_matrix:\n");
	for(j=0;j<molecule_improper_types;++j)printf("~ map %d --> %d\n",(*itmapping_matrix)[j][0],(*itmapping_matrix)[j][1]);

	printf("\n$ general improper types registry:\n");
	for(j=0;j<*general_improper_types;++j)printf("[%d]\t%s\t%s\t%s\t%s\n",j+1,(*general_improper_types_registry)[j][0],(*general_improper_types_registry)[j][1],(*general_improper_types_registry)[j][2],(*general_improper_types_registry)[j][3]);
	}
    if(molecule_impropers>0){
    
	// [7]
	
	// initialize
	non_minus_above=0;
	non_minus_below=0;
	non_minus_current=0;
	// check current molecule i even column registry
	if(topo_boundaries[i][6]!=-1)non_minus_current=1;
	// traverse even column from 0 to i-1
	for(j=0;j<=i-1;++j)if(topo_boundaries[j][6]!=-1){non_minus_above=1;break;}
	// traverse even column from i+1 to molecules
	for(j=i+1;j<molecules;++j)if(topo_boundaries[j][6]!=-1){non_minus_below=1;break;}
	//
	if(verb==1){
	printf("\n$ minus flags:\n");
	printf("non_minus_above = %d\n",non_minus_above);
	printf("non_minus_below = %d\n",non_minus_below);
	printf("non_minus_current = %d\n\n",non_minus_current);}
	
	// case 1
    // this case is triggered when an improper is first introduced to the system, e.g. CH4 --> 2HC==CH2
    // or when there is only one improper in an intermediate slot
    if(non_minus_above==0 && non_minus_below==0)
    {
		if(non_minus_current==0)
		{
            // in_array holds the new values for the current molecule due to the building step
			in_array[0]=0;in_array[1]=molecule_impropers-1;
			if(verb==1)printf("$ case 1 in_array = [%d %d]\n",in_array[0],in_array[1]);
            // update current molecule values
			topo_boundaries[i][6]=in_array[0];
			topo_boundaries[i][7]=in_array[1];
            // tcv (topo current value) holds the difference in elements (delta_impropers)
			tcv=topo_boundaries[i][7]-topo_boundaries[i][6]+1;
			if(verb==1){printf("$ case 1 tcv = %d\n",tcv);
            printf("$ case 1 updated topo_boundaries:\n");
			for(j=0;j<molecules;++j)printf("[%d]\t%d\t%d\n",j+1,topo_boundaries[j][6],topo_boundaries[j][7]);}
            // allocate space
			*general_impropers_registry=(int**)malloc(tcv*sizeof(int*));
			for(j=0;j<tcv;++j)(*general_impropers_registry)[j]=(int*)malloc(5*sizeof(int));
            // place current molecule entries
			k=-1;
			for(j=topo_boundaries[i][6];j<=topo_boundaries[i][7];++j)
			{
				k=k+1;
				if(verb==1)printf("$ placing %d at %d...\n",k,j);
				(*general_impropers_registry)[j][0]=(*itmapping_matrix)[molecule_impropers_registry[k][0]-1][1];    // type mapping
				(*general_impropers_registry)[j][1]=molecule_impropers_registry[k][1]+atom_scaling_array[i];        // rescale
				(*general_impropers_registry)[j][2]=molecule_impropers_registry[k][2]+atom_scaling_array[i];        // rescale
				(*general_impropers_registry)[j][3]=molecule_impropers_registry[k][3]+atom_scaling_array[i];        // rescale
				(*general_impropers_registry)[j][4]=molecule_impropers_registry[k][4]+atom_scaling_array[i];        // rescale
			}
            // alter total
			*general_impropers=*general_impropers+tcv;
            //
			if(verb==1)printf("$ general impropers registry:\n");
			if(verb==1)for(j=0;j<*general_impropers;++j)printf("[%d]\t%d\t%d\t%d\t%d\t%d\n",j+1,(*general_impropers_registry)[j][0],(*general_impropers_registry)[j][1],(*general_impropers_registry)[j][2],(*general_impropers_registry)[j][3],(*general_impropers_registry)[j][4]);
		}
        else if(non_minus_current==1)
        {
            // delta holds the number of impropers prior to the building step
            delta=topo_boundaries[i][7]-topo_boundaries[i][6]+1;
            // in_array holds the new values for the current molecule due to the building step
            in_array[0]=0;in_array[1]=molecule_impropers-1;
            if(verb==1)printf("$ case 1 in_array = [%d %d]\n",in_array[0],in_array[1]);
            // update current molecule values
            topo_boundaries[i][6]=in_array[0];
            topo_boundaries[i][7]=in_array[1];
            // tcv (topo current value) holds the difference in elements (delta_impropers)
            tcv=topo_boundaries[i][7]-topo_boundaries[i][6]+1-delta;
            if(verb==1)printf("$ case 1 tcv = %d\n",tcv);
            //
            if(tcv>0)
            {
                // resize adding extra space
                if(verb==1)printf("$ case 1 init / final sizes: %d / %d\n",*general_impropers,*general_impropers+tcv);
                *general_impropers_registry=(int**)realloc(*general_impropers_registry,(*general_impropers+tcv)*sizeof(int*));
                for(j=0;j<tcv;++j)(*general_impropers_registry)[*general_impropers+j]=(int*)malloc(5*sizeof(int));
            }
            else if(tcv<0)
            {
                // resize removing space
                printf("$ case 1 init / final sizes: %d / %d\n",*general_impropers,*general_impropers+tcv);
                printf("\n$ tcv<0\n\n");
                exit(-1);
            }
            else
            {
                // size kept the same
                if(verb==1)printf("$ case 1 init / final sizes: %d / %d\n",*general_impropers,*general_impropers+tcv);
            }
            // place current molecule entries
            k=-1;
            for(j=topo_boundaries[i][6];j<=topo_boundaries[i][7];++j)
            {
                k=k+1;
                if(verb==1)printf("$ placing %d at %d...\n",k,j);
                (*general_impropers_registry)[j][0]=(*itmapping_matrix)[molecule_impropers_registry[k][0]-1][1];    // type mapping
                (*general_impropers_registry)[j][1]=molecule_impropers_registry[k][1]+atom_scaling_array[i];        // rescale
                (*general_impropers_registry)[j][2]=molecule_impropers_registry[k][2]+atom_scaling_array[i];        // rescale
                (*general_impropers_registry)[j][3]=molecule_impropers_registry[k][3]+atom_scaling_array[i];        // rescale
                (*general_impropers_registry)[j][4]=molecule_impropers_registry[k][4]+atom_scaling_array[i];        // rescale
            }
            // alter total
            *general_impropers=*general_impropers+tcv;
            //
            if(verb==1)printf("$ general dihedrals registry:\n");
            if(verb==1)for(j=0;j<*general_impropers;++j)printf("[%d]\t%d\t%d\t%d\t%d\t%d\n",j+1,(*general_impropers_registry)[j][0],(*general_impropers_registry)[j][1],(*general_impropers_registry)[j][2],(*general_impropers_registry)[j][3],(*general_impropers_registry)[j][4]);
        }
	}
	
	// case 2
	else if(non_minus_above==1 && non_minus_below==0)
	{
        // delta holds the number of impropers prior to the building step
		delta=0;if(non_minus_current==1)delta=topo_boundaries[i][7]-topo_boundaries[i][6]+1;
		if(verb==1)printf("$ case 2 delta = %d\n",delta);
        // tav (topo above value) holds the largest eligible position value above the current molecule
		tav=topo_boundaries[0][7];for(j=0;j<=i-1;++j)if(topo_boundaries[j][7]>tav)tav=topo_boundaries[j][7];
		if(verb==1)printf("$ case 2 tav = %d\n",tav);
        // in_array holds the new values for the current molecule due to the building step
		in_array[0]=0+tav+1;in_array[1]=molecule_impropers-1+tav+1;
		if(verb==1)printf("$ case 2 in_array = [%d %d]\n",in_array[0],in_array[1]);
        // update current molecule values
		topo_boundaries[i][6]=in_array[0];
		topo_boundaries[i][7]=in_array[1];
        // tcv (topo current value) holds the difference in elements (delta_impropers)
		tcv=topo_boundaries[i][7]-topo_boundaries[i][6]+1-delta;
		if(verb==1){printf("$ case 2 tcv = %d\n",tcv);
        //
		printf("$ case 2 updated topo_boundaries:\n");
		for(j=0;j<molecules;++j)printf("[%d]\t%d\t%d\n",j+1,topo_boundaries[j][6],topo_boundaries[j][7]);}
        //
		if(tcv>0)
		{
            // resize adding extra space
			if(verb==1)printf("$ case 2 init / final sizes: %d / %d\n",*general_impropers,*general_impropers+tcv);
			*general_impropers_registry=(int**)realloc(*general_impropers_registry,(*general_impropers+tcv)*sizeof(int*));
			for(j=0;j<tcv;++j)(*general_impropers_registry)[*general_impropers+j]=(int*)malloc(5*sizeof(int));
		}
		else if (tcv<0)
		{
            // resize removing space
			printf("$ case 2 init / final sizes: %d / %d\n",*general_impropers,*general_impropers+tcv);
			printf("\n$ tcv<0\n\n");
			exit(-1);
		}
		else
		{
            // size kept the same
			if(verb==1)printf("$ case 2 init / final sizes: %d / %d\n",*general_impropers,*general_impropers+tcv);
		}
        // place current molecule entries
		k=-1;
		for(j=topo_boundaries[i][6];j<=topo_boundaries[i][7];++j)
		{
			k=k+1;
			if(verb==1)printf("$ placing %d at %d...\n",k,j);
			(*general_impropers_registry)[j][0]=(*itmapping_matrix)[molecule_impropers_registry[k][0]-1][1];    // type mapping
			(*general_impropers_registry)[j][1]=molecule_impropers_registry[k][1]+atom_scaling_array[i];        // rescale
			(*general_impropers_registry)[j][2]=molecule_impropers_registry[k][2]+atom_scaling_array[i];        // rescale
			(*general_impropers_registry)[j][3]=molecule_impropers_registry[k][3]+atom_scaling_array[i];        // rescale
			(*general_impropers_registry)[j][4]=molecule_impropers_registry[k][4]+atom_scaling_array[i];        // rescale
		}
        // alter total
		*general_impropers=*general_impropers+tcv;
        //
		if(verb==1)printf("$ general impropers registry:\n");
		if(verb==1)for(j=0;j<*general_impropers;++j)printf("[%d]\t%d\t%d\t%d\t%d\t%d\n",j+1,(*general_impropers_registry)[j][0],(*general_impropers_registry)[j][1],(*general_impropers_registry)[j][2],(*general_impropers_registry)[j][3],(*general_impropers_registry)[j][4]);		
	}
	
	// case 3
	else if(non_minus_above==0 && non_minus_below==1)
	{
        // delta holds the number of impropers prior to the building step
		delta=0;if(non_minus_current==1)delta=topo_boundaries[i][7]-topo_boundaries[i][6]+1;
		if(verb==1)printf("$ case 3 delta = %d\n",delta);
        // upper_0 and lower_0 define the boundaries for the values below
		for(j=i+1;j<molecules;++j)if(topo_boundaries[j][6]!=-1){upper_0=topo_boundaries[j][6];break;}
		lower_0=topo_boundaries[i][7];
		for(j=i+1;j<molecules;++j)if(topo_boundaries[j][7]>lower_0){lower_0=topo_boundaries[j][7];}
		if(verb==1){printf("$ case 3 upper_0 = %d\n",upper_0);
		printf("$ case 3 lower_0 = %d\n",lower_0);}
        // in_array holds the new values for the current molecule due to the building step
		in_array[0]=0;in_array[1]=molecule_impropers-1;
		if(verb==1)printf("$ case 3 in_array = [%d %d]\n",in_array[0],in_array[1]);
        // update current molecule values
		topo_boundaries[i][6]=in_array[0];
		topo_boundaries[i][7]=in_array[1];
        // tcv (topo current value) holds the difference in elements (delta_impropers)
		tcv=topo_boundaries[i][7]-topo_boundaries[i][6]+1-delta;
		if(verb==1)printf("$ case 3 tcv = %d\n",tcv);
        // update non -1 values below current molecule
		for(j=i+1;j<molecules;++j){if(topo_boundaries[j][6]!=-1){topo_boundaries[j][6]=topo_boundaries[j][6]+tcv;topo_boundaries[j][7]=topo_boundaries[j][7]+tcv;}}
        //
		if(verb==1){printf("$ case 3 updated topo_boundaries:\n");
		for(j=0;j<molecules;++j)printf("[%d]\t%d\t%d\n",j+1,topo_boundaries[j][6],topo_boundaries[j][7]);}
        //
		if(tcv>0)
		{
            // resize adding extra space
			if(verb==1)printf("$ case 3 init / final sizes: %d / %d\n",*general_impropers,*general_impropers+tcv);
			*general_impropers_registry=(int**)realloc(*general_impropers_registry,(*general_impropers+tcv)*sizeof(int*));
			for(j=0;j<tcv;++j)(*general_impropers_registry)[*general_impropers+j]=(int*)malloc(5*sizeof(int));
            // push down entries
			if(verb==1)printf("$ case 3 push down: %d:%d --> %d:%d\n",upper_0,lower_0,upper_0+tcv,lower_0+tcv);
			for(j=lower_0;j>=upper_0;--j)
			{
				if(verb==1)printf("%d --> %d\n",j,j+tcv);
				(*general_impropers_registry)[j+tcv][0]=(*general_impropers_registry)[j][0];
				(*general_impropers_registry)[j+tcv][1]=(*general_impropers_registry)[j][1]+delta_atoms;    // rescale
				(*general_impropers_registry)[j+tcv][2]=(*general_impropers_registry)[j][2]+delta_atoms;    // rescale
				(*general_impropers_registry)[j+tcv][3]=(*general_impropers_registry)[j][3]+delta_atoms;    // rescale
				(*general_impropers_registry)[j+tcv][4]=(*general_impropers_registry)[j][4]+delta_atoms;    // rescale
			}
		}
		else if (tcv<0)
		{
            // resize removing space
			printf("$ case 3 init / final sizes: %d / %d\n",*general_impropers,*general_impropers+tcv);
			printf("\n$ tcv<0\n\n");
			exit(-1);
		}
		else
		{
            // size kept the same
            if(verb==1){
			printf("$ case 3 init / final sizes: %d / %d\n",*general_impropers,*general_impropers+tcv);
            // update for rescaling purposes
			printf("$ case 3 no shift: %d:%d <--> %d:%d\n",upper_0,lower_0,upper_0+tcv,lower_0+tcv);}
			for(j=lower_0;j>=upper_0;--j)
			{
				if(verb==1)printf("%d --> %d\n",j,j+tcv);
				(*general_impropers_registry)[j+tcv][0]=(*general_impropers_registry)[j][0];
				(*general_impropers_registry)[j+tcv][1]=(*general_impropers_registry)[j][1]+delta_atoms;    // rescale
				(*general_impropers_registry)[j+tcv][2]=(*general_impropers_registry)[j][2]+delta_atoms;    // rescale
				(*general_impropers_registry)[j+tcv][3]=(*general_impropers_registry)[j][3]+delta_atoms;    // rescale
				(*general_impropers_registry)[j+tcv][4]=(*general_impropers_registry)[j][4]+delta_atoms;    // rescale
			}
		}
        // place current molecule entries
		k=-1;
		for(j=topo_boundaries[i][6];j<=topo_boundaries[i][7];++j)
		{
			k=k+1;
			if(verb==1)printf("$ placing %d at %d...\n",k,j);
			(*general_impropers_registry)[j][0]=(*itmapping_matrix)[molecule_impropers_registry[k][0]-1][1];    // type mapping
			(*general_impropers_registry)[j][1]=molecule_impropers_registry[k][1]+atom_scaling_array[i];        // rescale
			(*general_impropers_registry)[j][2]=molecule_impropers_registry[k][2]+atom_scaling_array[i];        // rescale
			(*general_impropers_registry)[j][3]=molecule_impropers_registry[k][3]+atom_scaling_array[i];        // rescale
			(*general_impropers_registry)[j][4]=molecule_impropers_registry[k][4]+atom_scaling_array[i];        // rescale
		}
        // alter total
		*general_impropers=*general_impropers+tcv;
        //
		if(verb==1)printf("$ general impropers registry:\n");
		if(verb==1)for(j=0;j<*general_impropers;++j)printf("[%d]\t%d\t%d\t%d\t%d\t%d\n",j+1,(*general_impropers_registry)[j][0],(*general_impropers_registry)[j][1],(*general_impropers_registry)[j][2],(*general_impropers_registry)[j][3],(*general_impropers_registry)[j][4]);
	}
	
	// case 4
	else if(non_minus_above==1 && non_minus_below==1)
	{
        // delta holds the number of impropers prior to the building step
		delta=0;if(non_minus_current==1)delta=topo_boundaries[i][7]-topo_boundaries[i][6]+1;
		if(verb==1)printf("$ case 4 delta = %d\n",delta);
        // upper_0 and lower_0 define the boundaries for the values below
		for(j=i+1;j<molecules;++j)if(topo_boundaries[j][6]!=-1){upper_0=topo_boundaries[j][6];break;}
		lower_0=topo_boundaries[i][7];
		for(j=i+1;j<molecules;++j)if(topo_boundaries[j][7]>lower_0){lower_0=topo_boundaries[j][7];}
		if(verb==1){printf("$ case 4 upper_0 = %d\n",upper_0);
		printf("$ case 4 lower_0 = %d\n",lower_0);}
        // tav (topo above value) holds the largest eligible position value above the current molecule
		tav=topo_boundaries[0][7];for(j=0;j<=i-1;++j)if(topo_boundaries[j][7]>tav)tav=topo_boundaries[j][7];
		if(verb==1)printf("$ case 4 tav = %d\n",tav);
        // in_array holds the new values for the current molecule due to the building step
		in_array[0]=0+tav+1;in_array[1]=molecule_impropers-1+tav+1;
		if(verb==1)printf("$ case 4 in_array = [%d %d]\n",in_array[0],in_array[1]);
        // update current molecule values
		topo_boundaries[i][6]=in_array[0];
		topo_boundaries[i][7]=in_array[1];
        // tcv (topo current value) holds the difference in elements (delta_impropers)
		tcv=topo_boundaries[i][7]-topo_boundaries[i][6]+1-delta;
		if(verb==1)printf("$ case 4 tcv = %d\n",tcv);
        // update non -1 values below current molecule
		for(j=i+1;j<molecules;++j){if(topo_boundaries[j][6]!=-1){topo_boundaries[j][6]=topo_boundaries[j][6]+tcv;topo_boundaries[j][7]=topo_boundaries[j][7]+tcv;}}
        //
		if(verb==1){printf("$ case 4 updated topo_boundaries:\n");
		for(j=0;j<molecules;++j)printf("[%d]\t%d\t%d\n",j+1,topo_boundaries[j][6],topo_boundaries[j][7]);}
        //
		if(tcv>0)
		{
            // resize adding extra space
			if(verb==1)printf("$ case 4 init / final sizes: %d / %d\n",*general_impropers,*general_impropers+tcv);
			*general_impropers_registry=(int**)realloc(*general_impropers_registry,(*general_impropers+tcv)*sizeof(int*));
			for(j=0;j<tcv;++j)(*general_impropers_registry)[*general_impropers+j]=(int*)malloc(5*sizeof(int));
            // push down entries
			if(verb==1)printf("$ case 4 push down: %d:%d --> %d:%d\n",upper_0,lower_0,upper_0+tcv,lower_0+tcv);
			for(j=lower_0;j>=upper_0;--j)
			{
				if(verb==1)printf("%d --> %d\n",j,j+tcv);
				(*general_impropers_registry)[j+tcv][0]=(*general_impropers_registry)[j][0];
				(*general_impropers_registry)[j+tcv][1]=(*general_impropers_registry)[j][1]+delta_atoms;    // rescale
				(*general_impropers_registry)[j+tcv][2]=(*general_impropers_registry)[j][2]+delta_atoms;    // rescale
				(*general_impropers_registry)[j+tcv][3]=(*general_impropers_registry)[j][3]+delta_atoms;    // rescale
				(*general_impropers_registry)[j+tcv][4]=(*general_impropers_registry)[j][4]+delta_atoms;    // rescale
			}
		}
		else if (tcv<0)
		{
            // resize removing space
			printf("$ case 4 init / final sizes: %d / %d\n",*general_impropers,*general_impropers+tcv);
			printf("\n$ tcv<0\n\n");
			exit(-1);
		}
		else
		{
            // size kept the same
			if(verb==1){printf("$ case 4 init / final sizes: %d / %d\n",*general_impropers,*general_impropers+tcv);
            // update for rescaling purposes
			printf("$ case 4 no shift: %d:%d <--> %d:%d\n",upper_0,lower_0,upper_0+tcv,lower_0+tcv);}
			for(j=lower_0;j>=upper_0;--j)
			{
				if(verb==1)printf("%d --> %d\n",j,j+tcv);
				(*general_impropers_registry)[j+tcv][0]=(*general_impropers_registry)[j][0];
				(*general_impropers_registry)[j+tcv][1]=(*general_impropers_registry)[j][1]+delta_atoms;    // rescale
				(*general_impropers_registry)[j+tcv][2]=(*general_impropers_registry)[j][2]+delta_atoms;    // rescale
				(*general_impropers_registry)[j+tcv][3]=(*general_impropers_registry)[j][3]+delta_atoms;    // rescale
				(*general_impropers_registry)[j+tcv][4]=(*general_impropers_registry)[j][4]+delta_atoms;    // rescale
			}
		}
        // place current molecule entries
		k=-1;
		for(j=topo_boundaries[i][6];j<=topo_boundaries[i][7];++j)
		{
			k=k+1;
			if(verb==1)printf("$ placing %d at %d...\n",k,j);
			(*general_impropers_registry)[j][0]=(*itmapping_matrix)[molecule_impropers_registry[k][0]-1][1];    // type mapping
			(*general_impropers_registry)[j][1]=molecule_impropers_registry[k][1]+atom_scaling_array[i];        // rescale
			(*general_impropers_registry)[j][2]=molecule_impropers_registry[k][2]+atom_scaling_array[i];        // rescale
			(*general_impropers_registry)[j][3]=molecule_impropers_registry[k][3]+atom_scaling_array[i];        // rescale
			(*general_impropers_registry)[j][4]=molecule_impropers_registry[k][4]+atom_scaling_array[i];        // rescale
		}
        // alter total
		*general_impropers=*general_impropers+tcv;
        //
		if(verb==1)printf("$ general impropers registry:\n");
		if(verb==1)for(j=0;j<*general_impropers;++j)printf("[%d]\t%d\t%d\t%d\t%d\t%d\n",j+1,(*general_impropers_registry)[j][0],(*general_impropers_registry)[j][1],(*general_impropers_registry)[j][2],(*general_impropers_registry)[j][3],(*general_impropers_registry)[j][4]);
	}
	
	else
	{
		printf("unsupported impropers combination!!\n\n");exit(-1);
	}
	
	delete_flag=0;
    for(k=0;k<*itmap_rows;++k)
        if(icm_sum[k]==0)
            delete_flag=delete_flag+1;
    
    if(delete_flag>0)
    {
        if(verb==1)printf("\n$ delete flag>01!!!\n\n");
        gtmapping_array=(int*)malloc((*general_improper_types)*sizeof(int));
        backup=(char***)malloc((*general_improper_types)*sizeof(char**));
        for(j=0;j<*general_improper_types;++j)backup[j]=(char**)malloc(4*sizeof(char*));
        for(j=0;j<*general_improper_types;++j)for(k=0;k<4;++k)backup[j][k]=(char*)malloc(sub_length*sizeof(char));
        counter=0;
        for(j=0;j<*general_improper_types;++j)
        {
            for(k=0;k<*itmap_rows;++k)
            {
                if((strcmp((*general_improper_types_registry)[j][0],(*itmap)[k][0])==0 &&
                    strcmp((*general_improper_types_registry)[j][1],(*itmap)[k][1])==0 &&
                    strcmp((*general_improper_types_registry)[j][2],(*itmap)[k][2])==0 &&
                    strcmp((*general_improper_types_registry)[j][3],(*itmap)[k][3])==0)
                   ||
                   (strcmp((*general_improper_types_registry)[j][0],(*itmap)[k][0])==0 &&
                    strcmp((*general_improper_types_registry)[j][1],(*itmap)[k][1])==0 &&
                    strcmp((*general_improper_types_registry)[j][2],(*itmap)[k][3])==0 &&
                    strcmp((*general_improper_types_registry)[j][3],(*itmap)[k][2])==0)
                   ||
                   (strcmp((*general_improper_types_registry)[j][0],(*itmap)[k][0])==0 &&
                    strcmp((*general_improper_types_registry)[j][1],(*itmap)[k][2])==0 &&
                    strcmp((*general_improper_types_registry)[j][2],(*itmap)[k][1])==0 &&
                    strcmp((*general_improper_types_registry)[j][3],(*itmap)[k][3])==0)
                   ||
                   (strcmp((*general_improper_types_registry)[j][0],(*itmap)[k][0])==0 &&
                    strcmp((*general_improper_types_registry)[j][1],(*itmap)[k][2])==0 &&
                    strcmp((*general_improper_types_registry)[j][2],(*itmap)[k][3])==0 &&
                    strcmp((*general_improper_types_registry)[j][3],(*itmap)[k][1])==0)
                   ||
                   (strcmp((*general_improper_types_registry)[j][0],(*itmap)[k][0])==0 &&
                    strcmp((*general_improper_types_registry)[j][1],(*itmap)[k][3])==0 &&
                    strcmp((*general_improper_types_registry)[j][2],(*itmap)[k][1])==0 &&
                    strcmp((*general_improper_types_registry)[j][3],(*itmap)[k][2])==0)
                   ||
                   (strcmp((*general_improper_types_registry)[j][0],(*itmap)[k][0])==0 &&
                    strcmp((*general_improper_types_registry)[j][1],(*itmap)[k][3])==0 &&
                    strcmp((*general_improper_types_registry)[j][2],(*itmap)[k][2])==0 &&
                    strcmp((*general_improper_types_registry)[j][3],(*itmap)[k][1])==0))
                    break;
            }
            if(verb==1){printf("$ %s-%s-%s-%s is at position {%d}\n",(*general_improper_types_registry)[j][0],(*general_improper_types_registry)[j][1],(*general_improper_types_registry)[j][2],(*general_improper_types_registry)[j][3],k+1);
            printf("$ checking cm_sum: %d\n",icm_sum[k]);}
            if(icm_sum[k]==0)
            {
                gtmapping_array[j]=0;
            }
            else
            {
                counter=counter+1;
                gtmapping_array[j]=counter;
                sprintf(backup[counter-1][0],"%s",(*general_improper_types_registry)[j][0]);
                sprintf(backup[counter-1][1],"%s",(*general_improper_types_registry)[j][1]);
                sprintf(backup[counter-1][2],"%s",(*general_improper_types_registry)[j][2]);
                sprintf(backup[counter-1][3],"%s",(*general_improper_types_registry)[j][3]);
            }
            
        }
        if(verb==1){
        printf("$ gtmapping_array:\n");
        for(j=0;j<*general_improper_types;++j)printf("[%d]\t%d\n",j+1,gtmapping_array[j]);
        printf("$ altered types:\n");
        for(j=0;j<counter;++j)printf("[%d]\t%s-%s-%s-%s\n",j+1,backup[j][0],backup[j][1],backup[j][2],backup[j][3]);}
        
        for(j=0;j<*general_impropers;++j)(*general_impropers_registry)[j][0]=gtmapping_array[(*general_impropers_registry)[j][0]-1];
        for(j=0;j<counter;++j)
        {
            sprintf((*general_improper_types_registry)[j][0],"%s",backup[j][0]);
            sprintf((*general_improper_types_registry)[j][1],"%s",backup[j][1]);
            sprintf((*general_improper_types_registry)[j][2],"%s",backup[j][2]);
            sprintf((*general_improper_types_registry)[j][3],"%s",backup[j][3]);
        }
        for(j=counter;j<*general_improper_types;++j)
        {
            for(k=0;k<4;++k)
                free((*general_improper_types_registry)[j][k]);
            free((*general_improper_types_registry)[j]);
        }
        for(j=0;j<*general_improper_types;++j)
            for(k=0;k<4;++k)
                free(backup[j][k]);
        for(j=0;j<*general_improper_types;++j)free(backup[j]);
        free(backup);
        *general_improper_types=counter;
        free(gtmapping_array);
    }
	
    free(icm_sum);
        
    }
    else
    {
        if(verb==1)printf("\n$ current molecule has no improper angles!!\n\n");
        // initialize
        non_minus_above=0;
        non_minus_below=0;
        non_minus_current=0;
        // check current molecule i even column registry
        if(topo_boundaries[i][6]!=-1)non_minus_current=1;
        // traverse even column from 0 to i-1
        for(j=0;j<=i-1;++j)if(topo_boundaries[j][6]!=-1){non_minus_above=1;break;}
        // traverse even column from i+1 to molecules
        for(j=i+1;j<molecules;++j)if(topo_boundaries[j][6]!=-1){non_minus_below=1;break;}
        //
        if(verb==1){
        printf("\n$ minus flags:\n");
        printf("non_minus_above = %d\n",non_minus_above);
        printf("non_minus_below = %d\n",non_minus_below);
        printf("non_minus_current = %d\n\n",non_minus_current);}
        
        if(non_minus_below==1)
        {
			if(verb==1){
            printf("\n$ improper registries exist below!!\n");
            printf("\n$ delta_atoms = %d\n",delta_atoms);}
            // upper_0 and lower_0 define the boundaries for the values below
            for(j=i+1;j<molecules;++j)if(topo_boundaries[j][6]!=-1){upper_0=topo_boundaries[j][6];break;}
            lower_0=topo_boundaries[i][7];
            for(j=i+1;j<molecules;++j)if(topo_boundaries[j][7]>lower_0){lower_0=topo_boundaries[j][7];}
            if(verb==1){printf("$ case 3 upper_0 = %d\n",upper_0);
            printf("$ case 3 lower_0 = %d\n",lower_0);}
            tcv=0;
            for(j=lower_0;j>=upper_0;--j)
            {
                if(verb==1)printf("%d --> %d\n",j,j+tcv);
                (*general_impropers_registry)[j+tcv][0]=(*general_impropers_registry)[j][0];
                (*general_impropers_registry)[j+tcv][1]=(*general_impropers_registry)[j][1]+delta_atoms;    // rescale
                (*general_impropers_registry)[j+tcv][2]=(*general_impropers_registry)[j][2]+delta_atoms;    // rescale
                (*general_impropers_registry)[j+tcv][3]=(*general_impropers_registry)[j][3]+delta_atoms;    // rescale
                (*general_impropers_registry)[j+tcv][4]=(*general_impropers_registry)[j][4]+delta_atoms;    // rescale
            }
        }
        free(icm_sum);
        //getchar();
    }

    //

}
