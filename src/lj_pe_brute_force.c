
#include "builder.h"

double lj_pe_brute_force(int general_atoms, double lx, double ly, double lz, double cutoff, int **one_two, int **one_three, int **one_four,double **nb_FF,int *nb_type_array,
                         double *master_x_array, double *master_y_array,double *master_z_array,int *master_nx_array, int *master_ny_array,int *master_nz_array,
                         int *pairs_bf,int *host_id_array)
{
    //
    int i,j,k;
    double pe;
    double xi,yi,zi,xj,yj,zj,r2,r,epsilon,sigma,r6,r12,sigma2,sigma6,sigma12;
    int select,one_four_flag;
    
    *pairs_bf=0;
    pe=0.0;
    for(i=0;i<general_atoms;++i)
    {
        for(j=0;j<general_atoms;++j)
        {
            //if(i!=j)
            if(i!=j && (host_id_array[i]!=-2 && host_id_array[j]!=-2))
            {
                
                xi=master_x_array[i]-master_nx_array[i]*lx;
                yi=master_y_array[i]-master_ny_array[i]*ly;
                zi=master_z_array[i]-master_nz_array[i]*lz;
                
                xj=master_x_array[j]-master_nx_array[j]*lx;
                yj=master_y_array[j]-master_ny_array[j]*ly;
                zj=master_z_array[j]-master_nz_array[j]*lz;
                
                if(xj-xi>0.5*lx) xj=xj-lx;
                if(xj-xi<-0.5*lx)xj=xj+lx;
                if(yj-yi>0.5*ly) yj=yj-ly;
                if(yj-yi<-0.5*ly)yj=yj+ly;
                if(zj-zi>0.5*lz) zj=zj-lz;
                if(zj-zi<-0.5*lz)zj=zj+lz;
                
                r2 = (xj-xi) * (xj-xi) + (yj-yi) * (yj-yi) + (zj-zi) * (zj-zi);
                r=sqrt(r2);
                
                if(r<cutoff)
                {
                    
                    select=1;one_four_flag=0;
                    
                    for(k=0;k<4;++k)
                    {
                        if(one_two[i][k]==j+1)
                        {
                            select=0;break;
                        }
                    }
                    
                    if(select==1)
                    {
                        for(k=0;k<12;++k)
                        {
                            if(one_three[i][k]==j+1)
                            {
                                select=0;
                                break;
                            }
                        }
                    }
                    
                    if(select==1)
                    {
                        for(k=0;k<36;++k)
                        {
                            if(one_four[i][k]==j+1)
                            {
                                select=0;one_four_flag=1;
                                break;
                            }
                        }
                    }
                    
                    if(select==1)
                    {
                        
                        
                        *pairs_bf=*pairs_bf+1;
                        
                        epsilon=nb_FF[nb_type_array[i]][0]*nb_FF[nb_type_array[j]][0];
                        sigma=nb_FF[nb_type_array[i]][1]*nb_FF[nb_type_array[j]][1];
                        
                        r6=r2*r2*r2;
                        r12=r6*r6;
                        sigma2=sigma*sigma;
                        sigma6=sigma2*sigma2*sigma2;
                        sigma12=sigma6*sigma6;
                        pe=pe+0.5*4.0*epsilon*(sigma12/r12-sigma6/r6);
                        
                        //printf("%d-%d is a valid LJ pair\tr=%lf\tepsilon=%lf\tsigma=%lf\n",i+1,j+1,r,epsilon,sigma);
                    }
                    if(one_four_flag==1)
                    {
                        
                        
                        *pairs_bf=*pairs_bf+1;
                        
                        epsilon=nb_FF[nb_type_array[i]][0]*nb_FF[nb_type_array[j]][0];
                        sigma=nb_FF[nb_type_array[i]][1]*nb_FF[nb_type_array[j]][1];
                        
                        r6=r2*r2*r2;
                        r12=r6*r6;
                        sigma2=sigma*sigma;
                        sigma6=sigma2*sigma2*sigma2;
                        sigma12=sigma6*sigma6;
                        pe=pe+0.5*4.0*epsilon*(sigma12/r12-sigma6/r6);
                        
                        //printf("%d-%d is a valid LJ 1-4 pair\tr=%lf\tepsilon=%lf\tsigma=%lf\n",i+1,j+1,r,epsilon,sigma);
                    }
                    
                    
                }
                
                
                
            }
            
            
            
            
            
            
            
        }
        
        
        
        
        
    }
    return pe;
}
