
#include"builder.h"

#define bound_tol 1.0e-12

void calc_wrapped_coords_n_arrays(int general_atoms,double xlo,double ylo,double zlo,double xhi,double yhi,double zhi,double lx,double ly,double lz,
                                  double *master_x_array,double *master_y_array,double *master_z_array,int *master_nx_array,int *master_ny_array,int *master_nz_array)
{
    int i;
    double nx,ny,nz;
    double x_shift,y_shift,z_shift;
    double x_temp,y_temp,z_temp;
    for (i=0;i<general_atoms;++i)
    {
        
        // shift to place in supercell; fixes rounding issues...
        x_shift=0.0;y_shift=0.0;z_shift=0.0;
        if(fabs(master_x_array[i]-xlo)<bound_tol){x_shift=bound_tol;}//printf("@ supercell edge coords rounding: +x_shift: %d\n",i+1);}
        else if(fabs(master_x_array[i]-xhi)<bound_tol){x_shift=-bound_tol;}//printf("@ supercell edge coords rounding: -x_shift: %d\n",i+1);}
        if(fabs(master_y_array[i]-ylo)<bound_tol){y_shift=bound_tol;}//printf("@ supercell edge coords rounding: +y_shift: %d\n",i+1);}
        else if(fabs(master_y_array[i]-yhi)<bound_tol){y_shift=-bound_tol;}//printf("@ supercell edge coords rounding: -y_shift: %d\n",i+1);}
        if(fabs(master_z_array[i]-zlo)<bound_tol){z_shift=bound_tol;}//printf("@ supercell edge coords rounding: +z_shift: %d\n",i+1);}
        else if(fabs(master_z_array[i]-zhi)<bound_tol){z_shift=-bound_tol;}//printf("@ supercell edge coords rounding: -z_shift: %d\n",i+1);}
        x_temp=master_x_array[i]+x_shift;
        y_temp=master_y_array[i]+y_shift;
        z_temp=master_z_array[i]+z_shift;
        
        // unwrap
        nx=0.0;ny=0.0;nz=0.0;
        if(x_temp>xhi)
        {
            nx=(x_temp-xlo)/lx;
            nx=floor(nx);
        }
        else if(x_temp<xlo)
        {
            nx=(xhi-x_temp)/lx;
            nx=-floor(nx);
        }
        if(y_temp>yhi)
        {
            ny=(y_temp-ylo)/ly;
            ny=floor(ny);
        }
        else if(y_temp<ylo)
        {
            ny=(yhi-y_temp)/ly;
            ny=-floor(ny);
        }
        if(z_temp>zhi)
        {
            nz=(z_temp-zlo)/lz;
            nz=floor(nz);
        }
        else if(z_temp<zlo)
        {
            nz=(zhi-z_temp)/lz;
            nz=-floor(nz);
        }
        master_nx_array[i]=(int)nx;
        master_ny_array[i]=(int)ny;
        master_nz_array[i]=(int)nz;
        // update
        master_x_array[i]=x_temp;
        master_y_array[i]=y_temp;
        master_z_array[i]=z_temp;
        
    }
    
}
