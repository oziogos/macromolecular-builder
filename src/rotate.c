
#include"builder.h"

void rotate(double ux, double uy, double uz, double x_center, double y_center, double z_center, double theta, double x_old, double y_old, double z_old, double *x_new, double *y_new, double *z_new)
{
    double norm2,norm;
    double rads;
    rads=theta*pi/180.0;
    norm2=ux*ux+uy*uy+uz*uz;norm=sqrt(norm2);
    ux=ux/norm;uy=uy/norm;uz=uz/norm;
    *x_new=(ux*ux*(1.0-cos(rads))+cos(rads))*(x_old-x_center)+(ux*uy*(1.0-cos(rads))-uz*sin(rads))*(y_old-y_center)+(ux*uz*(1.0-cos(rads))+uy*sin(rads))*(z_old-z_center)+x_center;
    *y_new=(ux*uy*(1.0-cos(rads))+uz*sin(rads))*(x_old-x_center)+(uy*uy*(1.0-cos(rads))+cos(rads))*(y_old-y_center)+(uy*uz*(1.0-cos(rads))-ux*sin(rads))*(z_old-z_center)+y_center;
    *z_new=(ux*uz*(1.0-cos(rads))-uy*sin(rads))*(x_old-x_center)+(uy*uz*(1.0-cos(rads))+ux*sin(rads))*(y_old-y_center)+(uz*uz*(1.0-cos(rads))+cos(rads))*(z_old-z_center)+z_center;

}
