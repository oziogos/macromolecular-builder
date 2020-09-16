
// what lj_pe_linked_list_cell() needs to calculate the LJ potential energy:
// updated:
// - general atoms
// - one_two, one_three, one_four
// - nb_FF_for_lj, nb_type_array
// - master_x_array, master_y_array, master_z_array, master_nx_array, master_ny_array, master_nz_array
// - list

#include "builder.h"

double lj_pe_linked_list_cell(int general_atoms, double lx, double ly, double lz, double cutoff, int **one_two, int **one_three, int **one_four,double **nb_FF,int *nb_type_array,
                              double *master_x_array, double *master_y_array,double *master_z_array,int *master_nx_array, int *master_ny_array,int *master_nz_array,
                              
                              int M3, int *head, int *list, int **neighbors,
                              
                              int *pairs_llc,int *host_id_array
                              )
{
    //
    int i,j,k;
    double pe;
    double xi,yi,zi,xj,yj,zj,r2,r,epsilon,sigma,r6,r12,sigma2,sigma6,sigma12;
    int select,one_four_flag;
    
    int cellindex,cellpart,headpointer,listpointer,locali,localj,neighborcellindex,neighborcellpart;
    int *cellpartarray, *neighborcellpartarray;
    
    cellpartarray=(int*)malloc(general_atoms*sizeof(int));
    neighborcellpartarray=(int*)malloc(general_atoms*sizeof(int));
    
    *pairs_llc=0;
    
    pe=0.0;
    
    // cellindex-loop: run over all cells
    for (cellindex=0;cellindex<M3;++cellindex)
    {
        // populate cellpartarray: array to hold atom ids belonging to the central cellindex cell
        cellpart=0; // cellpart is the atom counter per cell
        headpointer=head[cellindex];
        if (headpointer!=-1)
        {
            cellpartarray[cellpart]=headpointer;
            cellpart=cellpart+1;
            listpointer=list[headpointer];
            while (listpointer!=-1)
            {
                if (listpointer!=-1)
                {
                    cellpartarray[cellpart]=listpointer;
                    cellpart=cellpart+1;
                    listpointer=list[listpointer];
                }
            }
        }
        // trigger calculation only if there are present atoms in cellindex; if-clause could be removed due to cellpart i-loop...
        if (cellpart!=0)
        {
            //printf("cell L %d (atoms %d):\t",cellindex+1,cellpart);
            //for(i=0;i<cellpart;++i)printf("%d\t",cellpartarray[i]+1);printf("\n");
            //getchar();
            // double i,j-loop and self-exclusion and implicit H if-clause
            // this code segment calculates interactions between atoms in the central cell
            for (i=0;i<cellpart;++i)
            {
                for (j=0;j<cellpart;++j)
                {
                    
                    locali=cellpartarray[i];
                    localj=cellpartarray[j];
                    
                    //if (i!=j)
                    if(i!=j && (host_id_array[locali]!=-2 && host_id_array[localj]!=-2))
                    {
                        //locali=cellpartarray[i];
                        //localj=cellpartarray[j];
                        
                        //
                        xi=master_x_array[locali]-master_nx_array[locali]*lx;
                        yi=master_y_array[locali]-master_ny_array[locali]*ly;
                        zi=master_z_array[locali]-master_nz_array[locali]*lz;
                        
                        xj=master_x_array[localj]-master_nx_array[localj]*lx;
                        yj=master_y_array[localj]-master_ny_array[localj]*ly;
                        zj=master_z_array[localj]-master_nz_array[localj]*lz;
                        
                        if(xj-xi>0.5*lx) xj=xj-lx;
                        if(xj-xi<-0.5*lx)xj=xj+lx;
                        if(yj-yi>0.5*ly) yj=yj-ly;
                        if(yj-yi<-0.5*ly)yj=yj+ly;
                        if(zj-zi>0.5*lz) zj=zj-lz;
                        if(zj-zi<-0.5*lz)zj=zj+lz;
                        
                        r2 = (xj-xi) * (xj-xi) + (yj-yi) * (yj-yi) + (zj-zi) * (zj-zi);
                        r=sqrt(r2);
                        
                        //printf("$$ %d-%d\tr=%lf\n",locali+1,localj+1,r);

                        if(r<cutoff)
                        {
                            
                            select=1;one_four_flag=0;
                            
                            for(k=0;k<4;++k)
                            {
                                if(one_two[locali][k]==localj+1)
                                {
                                    select=0;break;
                                }
                            }
                            
                            if(select==1)
                            {
                                for(k=0;k<12;++k)
                                {
                                    if(one_three[locali][k]==localj+1)
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
                                    if(one_four[locali][k]==localj+1)
                                    {
                                        select=0;one_four_flag=1;
                                        break;
                                    }
                                }
                            }
                            
                            if(select==1)
                            {
                                
                                
                                *pairs_llc=*pairs_llc+1;
                                
                                epsilon=nb_FF[nb_type_array[locali]][0]*nb_FF[nb_type_array[localj]][0];
                                sigma=nb_FF[nb_type_array[locali]][1]*nb_FF[nb_type_array[localj]][1];
                                
                                r6=r2*r2*r2;
                                r12=r6*r6;
                                sigma2=sigma*sigma;
                                sigma6=sigma2*sigma2*sigma2;
                                sigma12=sigma6*sigma6;
                                pe=pe+0.5*4.0*epsilon*(sigma12/r12-sigma6/r6);
                                
                                //printf("%d-%d is a valid LJ pair\tr=%lf\tepsilon=%lf\tsigma=%lf\n",locali+1,localj+1,r,epsilon,sigma);
                                
                            }
                            if(one_four_flag==1)
                            {
                                
                                
                                *pairs_llc=*pairs_llc+1;
                                
                                epsilon=nb_FF[nb_type_array[locali]][0]*nb_FF[nb_type_array[localj]][0];
                                sigma=nb_FF[nb_type_array[locali]][1]*nb_FF[nb_type_array[localj]][1];
                                
                                r6=r2*r2*r2;
                                r12=r6*r6;
                                sigma2=sigma*sigma;
                                sigma6=sigma2*sigma2*sigma2;
                                sigma12=sigma6*sigma6;
                                pe=pe+0.5*4.0*epsilon*(sigma12/r12-sigma6/r6);
                                
                                //printf("%d-%d is a valid LJ 1-4 pair\tr=%lf\tepsilon=%lf\tsigma=%lf\n",locali+1,localj+1,r,epsilon,sigma);

                            }
                            
                            
                        }
                        //
                    }
                }
            }
        }
        // loop on the 26 first neighboring cells
        for (neighborcellindex=0;neighborcellindex<26;++neighborcellindex)
        {
            // populate neighborcellpartarray: array to hold atomic ids of neighborcellindex cell
            neighborcellpart=0; // counter
            headpointer=head[neighbors[neighborcellindex][cellindex]];
            if (headpointer!=-1)
            {
                neighborcellpartarray[neighborcellpart]=headpointer;
                neighborcellpart=neighborcellpart+1;
                listpointer=list[headpointer];
                while (listpointer!=-1)
                {
                    if (listpointer!=-1)
                    {
                        neighborcellpartarray[neighborcellpart]=listpointer;
                        neighborcellpart=neighborcellpart+1;
                        listpointer=list[listpointer];
                    }
                }
            }
            // apply only when there are atoms in the current neighborcellindex cell
            if (neighborcellpart!=0)
            {
                //printf("cell L %d (atoms %d) cell N %d (atoms %d):\t",cellindex+1,cellpart,neighbors[neighborcellindex][cellindex]+1,neighborcellpart);
                //for(i=0;i<neighborcellpart;++i)printf("%d\t",neighborcellpartarray[i]+1);printf("\n");
                //getchar();
                // double i,j-loop and implicit H if-clause
                // this code segment calculates interaction between atoms belonging to the central cellindex cell and the adjacent neighborcellindex cell
                for (i=0;i<cellpart;++i)
                {
                    for (j=0;j<neighborcellpart;++j)
                    {
                        locali=cellpartarray[i];
                        localj=neighborcellpartarray[j];
                        
                        if(host_id_array[locali]!=-2 && host_id_array[localj]!=-2)
                        {
                            xi=master_x_array[locali]-master_nx_array[locali]*lx;
                            yi=master_y_array[locali]-master_ny_array[locali]*ly;
                            zi=master_z_array[locali]-master_nz_array[locali]*lz;
                            
                            xj=master_x_array[localj]-master_nx_array[localj]*lx;
                            yj=master_y_array[localj]-master_ny_array[localj]*ly;
                            zj=master_z_array[localj]-master_nz_array[localj]*lz;
                            
                            if(xj-xi>0.5*lx) xj=xj-lx;
                            if(xj-xi<-0.5*lx)xj=xj+lx;
                            if(yj-yi>0.5*ly) yj=yj-ly;
                            if(yj-yi<-0.5*ly)yj=yj+ly;
                            if(zj-zi>0.5*lz) zj=zj-lz;
                            if(zj-zi<-0.5*lz)zj=zj+lz;
                            
                            r2 = (xj-xi) * (xj-xi) + (yj-yi) * (yj-yi) + (zj-zi) * (zj-zi);
                            r=sqrt(r2);
                            
                            //printf("$$ %d-%d\tr=%lf\n",locali+1,localj+1,r);

                            if(r<cutoff)
                            {
                                
                                select=1;one_four_flag=0;
                                
                                for(k=0;k<4;++k)
                                {
                                    if(one_two[locali][k]==localj+1)
                                    {
                                        select=0;break;
                                    }
                                }
                                
                                if(select==1)
                                {
                                    for(k=0;k<12;++k)
                                    {
                                        if(one_three[locali][k]==localj+1)
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
                                        if(one_four[locali][k]==localj+1)
                                        {
                                            select=0;one_four_flag=1;
                                            break;
                                        }
                                    }
                                }
                                
                                if(select==1)
                                {
                                    
                                    
                                    *pairs_llc=*pairs_llc+1;
                                    
                                    epsilon=nb_FF[nb_type_array[locali]][0]*nb_FF[nb_type_array[localj]][0];
                                    sigma=nb_FF[nb_type_array[locali]][1]*nb_FF[nb_type_array[localj]][1];
                                    
                                    r6=r2*r2*r2;
                                    r12=r6*r6;
                                    sigma2=sigma*sigma;
                                    sigma6=sigma2*sigma2*sigma2;
                                    sigma12=sigma6*sigma6;
                                    pe=pe+0.5*4.0*epsilon*(sigma12/r12-sigma6/r6);
                                    
                                    //printf("%d-%d is a valid LJ pair\tr=%lf\tepsilon=%lf\tsigma=%lf\n",locali+1,localj+1,r,epsilon,sigma);

                                }
                                if(one_four_flag==1)
                                {
                                    
                                    
                                    *pairs_llc=*pairs_llc+1;
                                    
                                    epsilon=nb_FF[nb_type_array[locali]][0]*nb_FF[nb_type_array[localj]][0];
                                    sigma=nb_FF[nb_type_array[locali]][1]*nb_FF[nb_type_array[localj]][1];
                                    
                                    r6=r2*r2*r2;
                                    r12=r6*r6;
                                    sigma2=sigma*sigma;
                                    sigma6=sigma2*sigma2*sigma2;
                                    sigma12=sigma6*sigma6;
                                    pe=pe+0.5*4.0*epsilon*(sigma12/r12-sigma6/r6);
                                    
                                    //printf("%d-%d is a valid LJ 1-4 pair\tr=%lf\tepsilon=%lf\tsigma=%lf\n",locali+1,localj+1,r,epsilon,sigma);

                                }
                                
                                
                            }
                            //
                        
                        }
                        //
                        
                    }
                }
            }
            //
        }
    }
    
    free(cellpartarray);free(neighborcellpartarray);
    return pe;
}
