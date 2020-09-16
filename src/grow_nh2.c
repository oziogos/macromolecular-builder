
#define nh3_N1_x_ref  0.00000
#define nh3_N1_y_ref  0.00000
#define nh3_N1_z_ref  0.00000
#define nh3_H1_x_ref  1.04311
#define nh3_H1_y_ref -0.00709
#define nh3_H1_z_ref  0.01926
#define nh3_H2_x_ref -0.33802
#define nh3_H2_y_ref  0.61789
#define nh3_H2_z_ref  0.76971
#define nh3_H3_x_ref -0.33803
#define nh3_H3_y_ref -0.96948
#define nh3_H3_z_ref  0.18528


#include"builder.h"

void rotate(double ux, double uy, double uz, double x_center, double y_center, double z_center, double theta, double x_old, double y_old, double z_old, double *x_new, double *y_new, double *z_new);
int equal(double A,double B);

void grow_nh2(int host_index,int selected_GS,int bonds,int max_1_2,int *atom_id,double *x,double *y,double *z,char **atom_type,char **atom_name,int *B1,int *B2,char **B_type,int *one_two,double *new_x,double *new_y,double *new_z,char *back)
{
    int i,j;
    
    double lamda,d,d0;
    double nh3_N1_x,nh3_N1_y,nh3_N1_z;
    double nh3_H1_x,nh3_H1_y,nh3_H1_z;
    double nh3_H2_x,nh3_H2_y,nh3_H2_z;
    double nh3_H3_x,nh3_H3_y,nh3_H3_z;
    
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
        
    // 2. displace the GS along the 1-2 bond in order to achieve the X-N.3 bond length
    // 2.1 resolve target bond length
    if(strcmp(atom_type[host_index],"C.ar")==0)d=d_C_ar_N_3;
    else if(strcmp(atom_type[host_index],"C.3")==0)d=d_C_3_N_3;
    else if(strcmp(atom_type[host_index],"O.3")==0)d=d_O_3_N_3;
    else if(strcmp(atom_type[host_index],"C.2")==0)d=d_C_2_N_3;
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
    // 3. overwrite atom type: from H to N.3
    sprintf(atom_type[selected_GS-1],"%s","N.3");
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
    // 5. place NH3
    // 5.1 move ch3
    nh3_N1_x=nh3_N1_x_ref+x[selected_GS-1];
    nh3_N1_y=nh3_N1_y_ref+y[selected_GS-1];
    nh3_N1_z=nh3_N1_z_ref+z[selected_GS-1];
    nh3_H1_x=nh3_H1_x_ref+x[selected_GS-1];
    nh3_H1_y=nh3_H1_y_ref+y[selected_GS-1];
    nh3_H1_z=nh3_H1_z_ref+z[selected_GS-1];
    nh3_H2_x=nh3_H2_x_ref+x[selected_GS-1];
    nh3_H2_y=nh3_H2_y_ref+y[selected_GS-1];
    nh3_H2_z=nh3_H2_z_ref+z[selected_GS-1];
    nh3_H3_x=nh3_H3_x_ref+x[selected_GS-1];
    nh3_H3_y=nh3_H3_y_ref+y[selected_GS-1];
    nh3_H3_z=nh3_H3_z_ref+z[selected_GS-1];
    // 6. store the bond vector
    vx=x[host_index]-nh3_N1_x;
    vy=y[host_index]-nh3_N1_y;
    vz=z[host_index]-nh3_N1_z;
    // 7. store the N.3-H alignment vector of NH3
    wx=nh3_H1_x-nh3_N1_x;
    wy=nh3_H1_y-nh3_N1_y;
    wz=nh3_H1_z-nh3_N1_z;
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
        rotate(ux,uy,uz,nh3_N1_x,nh3_N1_y,nh3_N1_z,-theta,nh3_H2_x,nh3_H2_y,nh3_H2_z,&nh3_H2_x,&nh3_H2_y,&nh3_H2_z);
        rotate(ux,uy,uz,nh3_N1_x,nh3_N1_y,nh3_N1_z,-theta,nh3_H3_x,nh3_H3_y,nh3_H3_z,&nh3_H3_x,&nh3_H3_y,&nh3_H3_z);
    }
    // inverse
    else if(align==-1)
    {
        void;
    }
    // 10. pick a H from the -NH2 fragment
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
    // 12. calculate the proper dihedral angle H-N-X-X
    //
    xi=nh3_H2_x;yi=nh3_H2_y;zi=nh3_H2_z;
    xj=nh3_N1_x;yj=nh3_N1_y;zj=nh3_N1_z;
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
    // 13. rotate about the bond vector to achieve the staggered conformation
    norm1=sqrt(vx*vx+vy*vy+vz*vz);
    vx=vx/norm1;vy=vy/norm1;vz=vz/norm1;
    rotate(vx,vy,vz,nh3_N1_x,nh3_N1_y,nh3_N1_z,phi+60.0,nh3_H2_x,nh3_H2_y,nh3_H2_z,&nh3_H2_x,&nh3_H2_y,&nh3_H2_z);
    rotate(vx,vy,vz,nh3_N1_x,nh3_N1_y,nh3_N1_z,phi+60.0,nh3_H3_x,nh3_H3_y,nh3_H3_z,&nh3_H3_x,&nh3_H3_y,&nh3_H3_z);
    
    new_x[0]=nh3_H2_x;new_y[0]=nh3_H2_y;new_z[0]=nh3_H2_z;
    new_x[1]=nh3_H3_x;new_y[1]=nh3_H3_y;new_z[1]=nh3_H3_z;
    
}


