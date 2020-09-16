#include"builder.h"

void SYBYLtoUFF(int atoms_core,char **atom_type_core,char **atom_type_core_FF);

void SYBYLtoDRE(int atoms_core,char **atom_type_core,char **atom_type_core_FF);

void LJ_params(int general_atoms, int general_atom_types,char **master_species_array, char **general_species_registry,
               int entries_UFF, char **species_UFF_array, double *D_UFF_array, double *x_UFF_array,
               int entries_DRE, char **species_DRE_array, double *D0_DRE_array,double *Rvdw0_DRE_array,
               double ***nb_FF, int **nb_type_array, int FF)
{
    char **atom_type_core_FF;
    int i,j,n;
    char **species_intermed_array;
    double ljconvert;

    ljconvert=pow(2.0,-1.0/6.0);

    atom_type_core_FF=(char**)malloc(general_atoms*sizeof(char*));
    for(j=0;j<general_atoms;++j)atom_type_core_FF[j]=(char*)malloc(6*sizeof(char));
    *nb_type_array=(int*)malloc(general_atoms*sizeof(int));
    
    if(FF==0)
    {
        // UFF
        SYBYLtoUFF(general_atoms,master_species_array,atom_type_core_FF);

        species_intermed_array=(char**)malloc(general_atom_types*sizeof(char*));for(j=0;j<general_atom_types;++j)species_intermed_array[j]=(char*)malloc(sub_length*sizeof(char));
        for(j=0;j<general_atom_types;++j)sprintf(species_intermed_array[j],"%s",general_species_registry[j]);

        for(j=0;j<general_atom_types;++j)
        {
            for(n=0;n<general_atoms;++n)
            {
                if(strcmp(species_intermed_array[j],master_species_array[n])==0)
                {
                    sprintf(species_intermed_array[j],"%s",atom_type_core_FF[n]);
                    break;
                }
            }
        }

        for(i=0;i<general_atoms;++i)
        {
            for(j=0;j<general_atom_types;++j)if(strcmp(atom_type_core_FF[i],species_intermed_array[j])==0){(*nb_type_array)[i]=j;break;}
        }

        //--------------------------------------------------------------------------
        // resolve nb parameters
        // we use general_atom_types unaltered because UFF utilizes a "1-1" mapping
        // for the atoms and hybridization states that are used in the building code
        // nb_FF will hold UFF D and x parameters

        *nb_FF=(double**)malloc(general_atom_types*sizeof(double*));
        for(j=0;j<general_atom_types;++j)(*nb_FF)[j]=(double*)malloc(2*sizeof(double));
        for(j=0;j<general_atom_types;++j)
        {
            for(n=0;n<entries_UFF;++n)
            {
                // identify FF parameters
                if(strcmp(species_intermed_array[j],species_UFF_array[n])==0)
                {
                    (*nb_FF)[j][0]=D_UFF_array[n];
                    (*nb_FF)[j][1]=x_UFF_array[n];
                    // convert
                    (*nb_FF)[j][0]=sqrt((*nb_FF)[j][0]);
                    (*nb_FF)[j][1]=sqrt((*nb_FF)[j][1]*ljconvert);
                    break;
                }
            }
        }
        
        for(j=0;j<general_atoms;++j)free(atom_type_core_FF[j]);free(atom_type_core_FF);
        for(j=0;j<general_atom_types;++j)free(species_intermed_array[j]);free(species_intermed_array);
    }
    else if (FF==1)
    {
        // DRE
        SYBYLtoDRE(general_atoms,master_species_array,atom_type_core_FF);
        
        species_intermed_array=(char**)malloc(general_atom_types*sizeof(char*));for(j=0;j<general_atom_types;++j)species_intermed_array[j]=(char*)malloc(sub_length*sizeof(char));
        for(j=0;j<general_atom_types;++j)sprintf(species_intermed_array[j],"%s",general_species_registry[j]);
        
        for(j=0;j<general_atom_types;++j)
        {
            for(n=0;n<general_atoms;++n)
            {
                if(strcmp(species_intermed_array[j],master_species_array[n])==0)
                {
                    sprintf(species_intermed_array[j],"%s",atom_type_core_FF[n]);
                    break;
                }
            }
        }
        
        for(i=0;i<general_atoms;++i)
        {
            for(j=0;j<general_atom_types;++j)if(strcmp(atom_type_core_FF[i],species_intermed_array[j])==0){(*nb_type_array)[i]=j;break;}
        }
        
        //--------------------------------------------------------------------------
        // resolve nb parameters
        // we use general_atom_types unaltered because DRE utilizes a "1-1" mapping
        // for the atoms and hybridization states that are used in the building code
        // nb_FF will hold UFF D and x parameters
        //
        // - hb not included!
        //
        
        *nb_FF=(double**)malloc(general_atom_types*sizeof(double*));
        for(j=0;j<general_atom_types;++j)(*nb_FF)[j]=(double*)malloc(2*sizeof(double));
        for(j=0;j<general_atom_types;++j)
        {
            for(n=0;n<entries_DRE;++n)
            {
                // identify FF parameters
                if(strcmp(species_intermed_array[j],species_DRE_array[n])==0)
                {
                    (*nb_FF)[j][0]=D0_DRE_array[n];
                    (*nb_FF)[j][1]=Rvdw0_DRE_array[n];
                    // convert
                    (*nb_FF)[j][0]=sqrt((*nb_FF)[j][0]);
                    (*nb_FF)[j][1]=sqrt((*nb_FF)[j][1]*ljconvert);
                    break;
                }
            }
        }
        
        for(j=0;j<general_atoms;++j)free(atom_type_core_FF[j]);free(atom_type_core_FF);
        for(j=0;j<general_atom_types;++j)free(species_intermed_array[j]);free(species_intermed_array);
    }
    
}
