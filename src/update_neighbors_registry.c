
#include"builder.h"

void update_neighbors_registry(int i,int general_atoms,int delta_atoms,int molecules,int molecule_atoms,int atoms_first_entry,int atoms_last_entry,
                               char **master_species_array,int **general_neighbors_registry)
{
    int j;
    if(delta_atoms>0)
    {
        (*general_neighbors_registry)=(int*)realloc((*general_neighbors_registry),general_atoms*sizeof(int));
        if(i!=molecules-1)
        {
            for(j=general_atoms-delta_atoms-1;j>atoms_last_entry;--j)
            {
                (*general_neighbors_registry)[j+delta_atoms]=(*general_neighbors_registry)[j];
            }
        }
        
        for(j=atoms_first_entry;j<atoms_first_entry+molecule_atoms;++j)
        {
            if(strcmp(master_species_array[j],"C.3")==0)
                (*general_neighbors_registry)[j]=4;
            else if(strcmp(master_species_array[j],"C.2")==0)
                (*general_neighbors_registry)[j]=3;
            else if(strcmp(master_species_array[j],"C.ar")==0)
                (*general_neighbors_registry)[j]=3;
            else if(strcmp(master_species_array[j],"O.3")==0)
                (*general_neighbors_registry)[j]=2;
            else if(strcmp(master_species_array[j],"O.2")==0)
                (*general_neighbors_registry)[j]=1;
            else if(strcmp(master_species_array[j],"H")==0)
                (*general_neighbors_registry)[j]=1;
            else if(strcmp(master_species_array[j],"N.3")==0 || strcmp(master_species_array[j],"N.am")==0)
                (*general_neighbors_registry)[j]=3;
            else
            {printf("cond 1: %s\n",master_species_array[j]);exit(-1);}
            
        }
    }
    else if(delta_atoms==0)
    {
        for(j=atoms_first_entry;j<atoms_first_entry+molecule_atoms;++j)
        {
            if(strcmp(master_species_array[j],"C.3")==0)
                (*general_neighbors_registry)[j]=4;
            else if(strcmp(master_species_array[j],"C.2")==0)
                (*general_neighbors_registry)[j]=3;
            else if(strcmp(master_species_array[j],"C.ar")==0)
                (*general_neighbors_registry)[j]=3;
            else if(strcmp(master_species_array[j],"O.3")==0)
                (*general_neighbors_registry)[j]=2;
            else if(strcmp(master_species_array[j],"O.2")==0)
                (*general_neighbors_registry)[j]=1;
            else if(strcmp(master_species_array[j],"H")==0)
                (*general_neighbors_registry)[j]=1;
            else if(strcmp(master_species_array[j],"N.3")==0 || strcmp(master_species_array[j],"N.am")==0)
                (*general_neighbors_registry)[j]=3;
            else
            {printf("cond 1: %s\n",master_species_array[j]);exit(-1);}
            
        }
    }
    else
    {
        printf("delta atoms<0\n\n");
        exit(-4);
    }
}
