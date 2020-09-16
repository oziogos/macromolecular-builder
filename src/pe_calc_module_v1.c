
#include"builder.h"


void topo_neighbors(int general_atoms, int general_bonds, int *master_B1_core_array, int *master_B2_core_array,
                    
                    int ***one_two, int ***one_three, int ***one_four);

double lj_pe_brute_force(int general_atoms, double lx, double ly, double lz, double cutoff, int **one_two, int **one_three, int **one_four,double **nb_FF,int *nb_type_array,
                         double *master_x_array, double *master_y_array,double *master_z_array,int *master_nx_array, int *master_ny_array,int *master_nz_array,
                         int *pairs_bf,int *host_id_array);

double lj_pe_linked_list_cell(int general_atoms, double lx, double ly, double lz, double cutoff, int **one_two, int **one_three, int **one_four,double **nb_FF,int *nb_type_array,
                              double *master_x_array, double *master_y_array,double *master_z_array,int *master_nx_array, int *master_ny_array,int *master_nz_array,
                              
                              int M3, int *head, int *list, int **neighbors,
                              
                              int *pairs_llc,int *host_id_array
                              );

void calc_head_list(int general_atoms, int Mx, int M2, int M3,
                    double xlo, double ylo, double zlo, double xhi, double yhi, double zhi, double lx, double ly, double lz, double lcx, double lcy, double lcz,
                    double *master_x_array,double *master_y_array, double *master_z_array, int *master_nx_array, int *master_ny_array, int *master_nz_array,
                    int *head, int **list);

void LJ_params(int general_atoms, int general_atom_types,char **master_species_array, char **general_species_registry,
               int entries_UFF, char **species_UFF_array, double *D_UFF_array, double *x_UFF_array,
               int entries_DRE, char **species_DRE_array, double *D0_DRE_array,double *Rvdw0_DRE_array,
               double ***nb_FF, int **nb_type_array, int FF);

void calc_wrapped_coords_n_arrays(int general_atoms,double xlo,double ylo,double zlo,double xhi,double yhi,double zhi,double lx,double ly,double lz,
                                  double *master_x_array,double *master_y_array,double *master_z_array,int *master_nx_array,int *master_ny_array,int *master_nz_array);

void pe_calc_module_v1(int general_atoms,int general_bonds,int general_atom_types,int entries_UFF,int entries_DRE,double cutoff,
                       double xlo,double ylo,double zlo,double xhi,double yhi,double zhi,double lx,double ly,double lz,double lcx,double lcy,double lcz,
                       int Mx,int My,int M3,
                       int *head,int **neighbors,double *pe_array,int counter,
                       double *master_x_array,double *master_y_array,double *master_z_array,int *host_id_array,
                       int *master_B1_core_array,int *master_B2_core_array,char **master_species_array,char **general_species_registry,
                       char **species_UFF_array,double *D_UFF_array,double *x_UFF_array,char **species_DRE_array,double *D0_DRE_array,double *Rvdw0_DRE_array,
                       int FF)
{
    int j;
    int **one_two,**one_three,**one_four;
    double **nb_FF_for_lj;
    int *nb_type_array;
    int *master_nx_array,*master_ny_array,*master_nz_array;
    int *list;
    double pe_llc;
    int pairs_llc;

    //
    // pe calculation
    //
    //-------------------------------------------------------------------------- not affected by coords ---
    topo_neighbors(general_atoms, general_bonds, master_B1_core_array, master_B2_core_array,&one_two, &one_three, &one_four);
    //-------------------------------------------------------------------------- not affected by coords ---
    LJ_params(general_atoms,general_atom_types,master_species_array,general_species_registry,
              entries_UFF,species_UFF_array,D_UFF_array,x_UFF_array,
              entries_DRE,species_DRE_array,D0_DRE_array,Rvdw0_DRE_array,
              &nb_FF_for_lj,&nb_type_array,FF);
    //--------------------------------------------------------------------------
    master_nx_array=(int*)malloc(general_atoms*sizeof(int));
    master_ny_array=(int*)malloc(general_atoms*sizeof(int));
    master_nz_array=(int*)malloc(general_atoms*sizeof(int));
    //.......................................................................... can be placed in loop ...
    // the latter arrays must be resized if delta_atoms!=0
    calc_wrapped_coords_n_arrays(general_atoms,xlo,ylo,zlo,xhi,yhi,zhi,lx,ly,lz,master_x_array,master_y_array,master_z_array,master_nx_array,master_ny_array,master_nz_array);
    //--------------------------------------------------------------------------
    calc_head_list(general_atoms,Mx,Mx*My,M3,xlo,ylo,zlo,xhi,yhi,zhi,lx,ly,lz,lcx,lcy,lcz,
                   master_x_array,master_y_array,master_z_array,master_nx_array,master_ny_array,master_nz_array,head,&list);
    //--------------------------------------------------------------------------
    /*
     pe_bf=lj_pe_brute_force(general_atoms,lx,ly,lz,cutoff,one_two,one_three,one_four,nb_FF_for_lj,nb_type_array,
     master_x_array,master_y_array,master_z_array,master_nx_array,master_ny_array,master_nz_array,
     &pairs_bf,host_id_array);
     printf("    pe:\t%lf\tbrute\tpairs: %d\n",pe_bf,pairs_bf);
     */
    //--------------------------------------------------------------------------
    pe_llc=lj_pe_linked_list_cell(general_atoms,lx,ly,lz,cutoff,one_two,one_three,one_four,nb_FF_for_lj,nb_type_array,
                                  master_x_array,master_y_array,master_z_array,master_nx_array,master_ny_array,master_nz_array,M3,head,list,neighbors,
                                  &pairs_llc,host_id_array);
    pe_array[counter]=pe_llc;
    //printf("    pe:\t%lf (kcal/mol)\tpairs: %d\n",pe_llc,pairs_llc);
    //printf("    pe:\t%lf\tcells\tpairs: %d\tdiff: %e%%\n",pe_llc,pairs_llc,(pe_llc/pe_bf-1.0)*100.0);
    //--------------------------------------------------------------------------
    // free memory: due to malloc in calc_head_list()
    free(list);
    //.......................................................................... end of code that can be placed in loop ...
    // free memory
    for(j=0;j<general_atom_types;++j)free(nb_FF_for_lj[j]);free(nb_FF_for_lj);
    for(j=0;j<general_atoms;++j)free(one_two[j]);free(one_two);
    for(j=0;j<general_atoms;++j)free(one_three[j]);free(one_three);
    for(j=0;j<general_atoms;++j)free(one_four[j]);free(one_four);
    free(nb_type_array);
    free(master_nx_array);free(master_ny_array);free(master_nz_array);

}
