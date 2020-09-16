
#define hoh_O_x_ref 0.0
#define hoh_O_y_ref 0.0
#define hoh_O_z_ref 0.0
#define hoh_H1_x_ref 0.0
#define hoh_H1_y_ref 0.0
#define hoh_H1_z_ref 0.99
#define hoh_H2_x_ref 0.958466
#define hoh_H2_y_ref 0.0
#define hoh_H2_z_ref -0.247876

#include"builder.h"

void rotate(double ux, double uy, double uz, double x_center, double y_center, double z_center, double theta, double x_old, double y_old, double z_old, double *x_new, double *y_new, double *z_new);
int equal(double A,double B);

void grow_oh(int host_index,int selected_GS,int bonds,int max_1_2,int *atom_id,double *x,double *y,double *z,char **atom_type,char **atom_name,int *B1,int *B2,char **B_type,int *one_two,double *new_x,double *new_y,double *new_z,char *back)
{
    int i,j;
    
    double lamda,d,d0;
    double hoh_O_x,hoh_O_y,hoh_O_z;
    double hoh_H1_x,hoh_H1_y,hoh_H1_z;
    double hoh_H2_x,hoh_H2_y,hoh_H2_z;
    
    double theta,dot,norm1,norm2,arg;
    double vx,vy,vz,wx,wy,wz;
    double ux,uy,uz,crossx,crossy,crossz;
    int align,D4_atom;
    double sign,coeff,phi,cross;
    double xi,yi,zi;
    double xj,yj,zj;
    double xk,yk,zk;
    double xl,yl,zl;
    
    //
    
    // 2. displace the GS along the 1-2 bond in order to achieve the X-O.3 bond length
    // 2.1 resolve target bond length
    if(strcmp(atom_type[host_index],"C.ar")==0)d=d_C_ar_O_3;
    else if(strcmp(atom_type[host_index],"C.3")==0)d=d_C_3_O_3;
    else if(strcmp(atom_type[host_index],"C.2")==0)d=d_C_2_O_3;
    else
    {
        printf("@@ ERROR!! unsupported 1-2 pair!!!\n\n");exit(-2);
    }
    // 2.2 displace
    d0=sqrt((x[selected_GS-1]-x[host_index])*(x[selected_GS-1]-x[host_index])+(y[selected_GS-1]-y[host_index])*(y[selected_GS-1]-y[host_index])+(z[selected_GS-1]-z[host_index])*(z[selected_GS-1]-z[host_index]));
    lamda=d/d0;
    x[selected_GS-1]=x[host_index]+lamda*(x[selected_GS-1]-x[host_index]);
    y[selected_GS-1]=y[host_index]+lamda*(y[selected_GS-1]-y[host_index]);
    z[selected_GS-1]=z[host_index]+lamda*(z[selected_GS-1]-z[host_index]);
    // 3. overwrite atom type: from H to O.3
    sprintf(atom_type[selected_GS-1],"%s","O.3");
    //sprintf(atom_name[selected_GS-1],"%s","O");
    sprintf(atom_name[selected_GS-1],"%s",back);
    // 4. overwrite bond order: store '1' (not needed in this case)
    for(i=0;i<bonds;++i)
    {
        if((B1[i]==atom_id[host_index] && B2[i]==selected_GS)||(B2[i]==atom_id[host_index] && B1[i]==selected_GS))
        {
            sprintf(B_type[i],"%s","1");
            break;
        }
    }
    // 5. place HOH
    // 5.1 move hoh
    hoh_O_x=hoh_O_x_ref+x[selected_GS-1];
    hoh_O_y=hoh_O_y_ref+y[selected_GS-1];
    hoh_O_z=hoh_O_z_ref+z[selected_GS-1];
    hoh_H1_x=hoh_H1_x_ref+x[selected_GS-1];
    hoh_H1_y=hoh_H1_y_ref+y[selected_GS-1];
    hoh_H1_z=hoh_H1_z_ref+z[selected_GS-1];
    hoh_H2_x=hoh_H2_x_ref+x[selected_GS-1];
    hoh_H2_y=hoh_H2_y_ref+y[selected_GS-1];
    hoh_H2_z=hoh_H2_z_ref+z[selected_GS-1];
    // 6. store the bond vector
    vx=x[host_index]-hoh_O_x;
    vy=y[host_index]-hoh_O_y;
    vz=z[host_index]-hoh_O_z;
    // 7. store the O.3-H alignment vector of HOH
    wx=hoh_H1_x-hoh_O_x;
    wy=hoh_H1_y-hoh_O_y;
    wz=hoh_H1_z-hoh_O_z;
    // 8. calculate the unit normal vector and 9. apply appropriate rotation (or inversion)
    // calculate angle
    dot=vx*wx+vy*wy+vz*wz;
    norm1=sqrt(vx*vx+vy*vy+vz*vz);
    norm2=sqrt(wx*wx+wy*wy+wz*wz);
    arg=dot/(norm1*norm2);
    // check the fragment's alignment status
    align=0;
    if(equal(arg,1.0)==1)align=1;
    if(equal(arg,-1.0)==1)align=-1;
    // rotate
    if(align==0)
    {
        theta=acos(arg)*180.0/pi;
        // calculate cross product
        crossx=-vz*wy + vy*wz;crossy=vz*wx - vx*wz;crossz=-vy*wx + vx*wy;
        norm1=sqrt(crossx*crossx+crossy*crossy+crossz*crossz);
        ux=crossx/norm1;uy=crossy/norm1;uz=crossz/norm1;
        // rotate
        rotate(ux,uy,uz,hoh_O_x,hoh_O_y,hoh_O_z,-theta,hoh_H2_x,hoh_H2_y,hoh_H2_z,&hoh_H2_x,&hoh_H2_y,&hoh_H2_z);
    }
    // inverse
    else if(align==-1)
    {
        void;
    }
    // 10. pick the H atom from -OH
    // no need...
    // 11. pick a 1-2 neighbor of the host atom (besides the GS)
    for(j=0;j<max_1_2;++j)
    {
        if(one_two[j]!=selected_GS)
        {
            D4_atom=one_two[j];
            break;
        }
    }
    // 12. calculate the proper dihedral angle H-O-X-X
    //
    xi=hoh_H2_x;yi=hoh_H2_y;zi=hoh_H2_z;
    xj=hoh_O_x;yj=hoh_O_y;zj=hoh_O_z;
    xk=x[host_index];yk=y[host_index];zk=z[host_index];
    xl=x[D4_atom-1];yl=y[D4_atom-1];zl=z[D4_atom-1];
    cross=(-xj*yi + xk*yi + xi*yj - xk*yj - xi*yk + xj*yk)*(-xk*yj + xl*yj + xj*yk - xl*yk - xj*yl + xk*yl) + (xj*zi - xk*zi - xi*zj + xk*zj + xi*zk - xj*zk)*(xk*zj - xl*zj - xj*zk + xl*zk + xj*zl - xk*zl) + (-yj*zi + yk*zi + yi*zj - yk*zj - yi*zk + yj*zk)*(-yk*zj + yl*zj + yj*zk - yl*zk - yj*zl + yk*zl);
    norm1=sqrt((-xj*yi + xk*yi + xi*yj - xk*yj - xi*yk + xj*yk)*(-xj*yi + xk*yi + xi*yj - xk*yj - xi*yk + xj*yk)+(xj*zi - xk*zi - xi*zj + xk*zj + xi*zk - xj*zk)*(xj*zi - xk*zi - xi*zj + xk*zj + xi*zk - xj*zk)+(-yj*zi + yk*zi + yi*zj - yk*zj - yi*zk + yj*zk)*(-yj*zi + yk*zi + yi*zj - yk*zj - yi*zk + yj*zk));
    norm2=sqrt((-xk*yj + xl*yj + xj*yk - xl*yk - xj*yl + xk*yl)*(-xk*yj + xl*yj + xj*yk - xl*yk - xj*yl + xk*yl)+(xk*zj - xl*zj - xj*zk + xl*zk + xj*zl - xk*zl)*(xk*zj - xl*zj - xj*zk + xl*zk + xj*zl - xk*zl)+(-yk*zj + yl*zj + yj*zk - yl*zk - yj*zl + yk*zl)*(-yk*zj + yl*zj + yj*zk - yl*zk - yj*zl + yk*zl));
    sign=(-xk*yj + xl*yj + xj*yk - xl*yk - xj*yl + xk*yl)*(-zi + zj) + (-yi + yj)*(xk*zj - xl*zj - xj*zk + xl*zk + xj*zl - xk*zl) + (-xi + xj)*(-yk*zj + yl*zj + yj*zk - yl*zk - yj*zl + yk*zl);
    coeff=cross/(norm1*norm2);
    if(coeff>1.0)coeff=1.0;
    if(coeff<-1.0)coeff=-1.0;
    // phi angle
    phi=0.0;
    if(sign>0.0){
        phi=acos(coeff)*180.0/pi;}
    else if (sign<0.0){phi=-acos(coeff)*180.0/pi;}
    //
    // 13. rotate about the bond vector to achieve a trans conformation
    norm1=sqrt(vx*vx+vy*vy+vz*vz);
    vx=vx/norm1;vy=vy/norm1;vz=vz/norm1;
    rotate(vx,vy,vz,hoh_O_x,hoh_O_y,hoh_O_z,phi+180.0,hoh_H2_x,hoh_H2_y,hoh_H2_z,&hoh_H2_x,&hoh_H2_y,&hoh_H2_z);
    
    new_x[0]=hoh_H2_x;new_y[0]=hoh_H2_y;new_z[0]=hoh_H2_z;

}