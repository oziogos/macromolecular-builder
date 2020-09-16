
#include"builder.h"

void SYBYLtoUFF(int atoms_core,char **atom_type_core,char **atom_type_core_FF)
{
    int i;
    
    for(i=0;i<atoms_core;++i)
    {
        if(strcmp(atom_type_core[i],"H")==0)sprintf(atom_type_core_FF[i],"%s","H_");
        else if(strcmp(atom_type_core[i],"C.3")==0)sprintf(atom_type_core_FF[i],"%s","C_3");
        else if(strcmp(atom_type_core[i],"C.2")==0)sprintf(atom_type_core_FF[i],"%s","C_2");
        else if(strcmp(atom_type_core[i],"C.ar")==0)sprintf(atom_type_core_FF[i],"%s","C_R");
        else if(strcmp(atom_type_core[i],"O.3")==0)sprintf(atom_type_core_FF[i],"%s","O_3");
        else if(strcmp(atom_type_core[i],"O.2")==0)sprintf(atom_type_core_FF[i],"%s","O_2");
        else if(strcmp(atom_type_core[i],"N.3")==0 || strcmp(atom_type_core[i],"N.am")==0)sprintf(atom_type_core_FF[i],"%s","N_3");
        else
        {
            printf("Not supported SYBYL type by UFF!\n");exit(-1);
        }
    }
}
