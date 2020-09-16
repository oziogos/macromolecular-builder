
#define ch2o_C_x_ref 0.0
#define ch2o_C_y_ref 0.0
#define ch2o_C_z_ref 0.0
#define ch2o_O_x_ref 0.0
#define ch2o_O_y_ref 0.0
#define ch2o_O_z_ref 1.21900
#define ch2o_H1_x_ref 0.938772
#define ch2o_H1_y_ref 0.0
#define ch2o_H1_z_ref -0.542
#define ch2o_H2_x_ref -0.938772
#define ch2o_H2_y_ref 0.0
#define ch2o_H2_z_ref -0.542

#include"builder.h"

void rotate(double ux, double uy, double uz, double x_center, double y_center, double z_center, double theta, double x_old, double y_old, double z_old, double *x_new, double *y_new, double *z_new);
int equal(double A,double B);

void grow_ch2(int host_index,int selected_GS,int bonds,int max_1_2,int *atom_id,double *x,double *y,double *z,char **atom_type,char **atom_name,int *B1,int *B2,char **B_type,int *one_two,double *new_x,double *new_y,double *new_z,char *back)
{

    int i,j;
    
    double lamda,d,d0;
    double ch2o_C_x,ch2o_C_y,ch2o_C_z;
    double ch2o_O_x,ch2o_O_y,ch2o_O_z;
    double ch2o_H1_x,ch2o_H1_y,ch2o_H1_z;
    double ch2o_H2_x,ch2o_H2_y,ch2o_H2_z;
    
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
        
    // 2. displace the GS along the 1-2 bond in order to achieve the X=C.2 bond length
    // 2.1 resolve target bond length
    if(strcmp(atom_type[host_index],"C.3")==0)d=d_C_2_C_3;
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
    // 3. overwrite atom type: from H to C.2
    sprintf(atom_type[selected_GS-1],"%s","C.2");
    //sprintf(atom_name[selected_GS-1],"%s","C");
    sprintf(atom_name[selected_GS-1],"%s",back);
    // 4. overwrite bond order: store '2'
    for(i=0;i<bonds;++i)
    {
        if((B1[i]==atom_id[host_index] && B2[i]==selected_GS)||(B2[i]==atom_id[host_index] && B1[i]==selected_GS))
        {
            sprintf(B_type[i],"%s","2");
            break;
        }
    }
    // 5. place O=CH2
    // 5.1 move O=CH2
    ch2o_C_x=ch2o_C_x_ref+x[selected_GS-1];
    ch2o_C_y=ch2o_C_y_ref+y[selected_GS-1];
    ch2o_C_z=ch2o_C_z_ref+z[selected_GS-1];
    ch2o_O_x=ch2o_O_x_ref+x[selected_GS-1];
    ch2o_O_y=ch2o_O_y_ref+y[selected_GS-1];
    ch2o_O_z=ch2o_O_z_ref+z[selected_GS-1];
    ch2o_H1_x=ch2o_H1_x_ref+x[selected_GS-1];
    ch2o_H1_y=ch2o_H1_y_ref+y[selected_GS-1];
    ch2o_H1_z=ch2o_H1_z_ref+z[selected_GS-1];
    ch2o_H2_x=ch2o_H2_x_ref+x[selected_GS-1];
    ch2o_H2_y=ch2o_H2_y_ref+y[selected_GS-1];
    ch2o_H2_z=ch2o_H2_z_ref+z[selected_GS-1];
    // 6. store the bond vector
    vx=x[host_index]-ch2o_C_x;
    vy=y[host_index]-ch2o_C_y;
    vz=z[host_index]-ch2o_C_z;
    // 7. store the O.2=C.2 alignment vector of O=CH2
    wx=ch2o_O_x-ch2o_C_x;
    wy=ch2o_O_y-ch2o_C_y;
    wz=ch2o_O_z-ch2o_C_z;
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
        rotate(ux,uy,uz,ch2o_C_x,ch2o_C_y,ch2o_C_z,-theta,ch2o_H2_x,ch2o_H2_y,ch2o_H2_z,&ch2o_H2_x,&ch2o_H2_y,&ch2o_H2_z);
        rotate(ux,uy,uz,ch2o_C_x,ch2o_C_y,ch2o_C_z,-theta,ch2o_H1_x,ch2o_H1_y,ch2o_H1_z,&ch2o_H1_x,&ch2o_H1_y,&ch2o_H1_z);
    }
    // inverse
    else if(align==-1)
    {
        void;
    }
    // 10. pick a H atom from =CH2
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
    // 12. calculate the proper dihedral angle H-C=X-X
    //
    xi=ch2o_H1_x;yi=ch2o_H1_y;zi=ch2o_H1_z;
    xj=ch2o_C_x;yj=ch2o_C_y;zj=ch2o_C_z;
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
    // 13. rotate about the bond vector to achieve a cis conformation
    norm1=sqrt(vx*vx+vy*vy+vz*vz);
    vx=vx/norm1;vy=vy/norm1;vz=vz/norm1;
    rotate(vx,vy,vz,ch2o_C_x,ch2o_C_y,ch2o_C_z,phi+180.0,ch2o_H2_x,ch2o_H2_y,ch2o_H2_z,&ch2o_H2_x,&ch2o_H2_y,&ch2o_H2_z);
    rotate(vx,vy,vz,ch2o_C_x,ch2o_C_y,ch2o_C_z,phi+180.0,ch2o_H1_x,ch2o_H1_y,ch2o_H1_z,&ch2o_H1_x,&ch2o_H1_y,&ch2o_H1_z);

    new_x[0]=ch2o_H2_x;new_y[0]=ch2o_H2_y;new_z[0]=ch2o_H2_z;
    new_x[1]=ch2o_H1_x;new_y[1]=ch2o_H1_y;new_z[1]=ch2o_H1_z;
    
}