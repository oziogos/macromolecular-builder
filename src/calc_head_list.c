#include "builder.h"


void calc_head_list(int general_atoms, int Mx, int M2, int M3,
                    double xlo, double ylo, double zlo, double xhi, double yhi, double zhi, double lx, double ly, double lz, double lcx, double lcy, double lcz,
                    double *master_x_array,double *master_y_array, double *master_z_array, int *master_nx_array, int *master_ny_array, int *master_nz_array,
                    int *head, int **list)
{
    int i,cindex;
    *list=(int*)malloc(general_atoms*sizeof(int));// list
    for (i=0;i<M3;++i){head[i]=-1;}
    for (i=0;i<general_atoms;++i){(*list)[i]=-1;}
    for (i=0;i<general_atoms;++i)
    {
        
        // populate
        cindex=floor((master_x_array[i]-master_nx_array[i]*lx-xlo)/lcx)+floor((master_y_array[i]-master_ny_array[i]*ly-ylo)/lcy)*Mx+floor((master_z_array[i]-master_nz_array[i]*lz-zlo)/lcz)*M2;
        (*list)[i]=head[cindex];
        head[cindex]=i;
        //
        //printf("atom %d: %.24lf,%.24lf,%.24lf\t cell=%d\n",i+1,master_x_array[i]-master_nx_array[i]*lx,master_y_array[i]-master_ny_array[i]*ly,master_z_array[i]-master_nz_array[i]*lz,cindex+1);
        //printf("\t\t\t\t\t\t\t\tunwrapped %.24lf,%.24lf,%.24lf\n",master_x_array[i],master_y_array[i],master_z_array[i]);
        
    }
}
