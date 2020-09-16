
#include"builder.h"

void write_xyz_preview(char *current_folder,int atoms,char **master_species_array,double *master_x_array,double *master_y_array,double *master_z_array,
                       double xlo,double xhi,double ylo,double yhi,double zlo,double zhi,int append_flag)
{
    FILE *fp;
    char file_path[cmax_length];
    int j;
    double lx,ly,lz;
    int nx,ny,nz;
    
    // xyz preview
    sprintf(file_path,"%s/preview.dat",current_folder);
    if(append_flag==0)
        fp=fopen(file_path,"w+");
    else
        fp=fopen(file_path,"a");
    fprintf(fp,"%d\n\n",atoms);
    for(j=0;j<atoms;++j)fprintf(fp,"%c\t%lf\t%lf\t%lf\n",master_species_array[j][0],master_x_array[j],master_y_array[j],master_z_array[j]);
    fclose(fp);
    
    // xyz preview
    lx=xhi-xlo;
    ly=yhi-ylo;
    lz=zhi-zlo;
    sprintf(file_path,"%s/preview_periodic.dat",current_folder);
    if(append_flag==0)
        fp=fopen(file_path,"w+");
    else
        fp=fopen(file_path,"a");
    fprintf(fp,"%d\n\n",atoms);
    for(j=0;j<atoms;++j)
    {
        nx=0;ny=0;nz=0;
        
        if(master_x_array[j]>xhi)
        {
            nx=(int)floor((master_x_array[j]-xlo)/lx);
        }
        if(master_x_array[j]<xlo)
        {
            nx=-(int)floor((xhi-master_x_array[j])/lx);
        }
        
        if(master_y_array[j]>yhi)
        {
            ny=(int)floor((master_y_array[j]-ylo)/ly);
        }
        if(master_y_array[j]<ylo)
        {
            ny=-(int)floor((yhi-master_y_array[j])/ly);
        }
        
        if(master_z_array[j]>zhi)
        {
            nz=(int)floor((master_z_array[j]-zlo)/lz);
        }
        if(master_z_array[j]<zlo)
        {
            nz=-(int)floor((zhi-master_z_array[j])/lz);
        }
        
        fprintf(fp,"%c\t%lf\t%lf\t%lf\n",master_species_array[j][0],master_x_array[j]-(double)nx*lx,master_y_array[j]-(double)ny*ly,master_z_array[j]-(double)nz*lz);
    }
    fclose(fp);
}
