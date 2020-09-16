
#include"builder.h"

void calculate_dihedrals(char *current_folder,int general_dihedrals,int backbones,double *master_x_array,double *master_y_array,double *master_z_array,
                         int **general_dihedrals_registry,char **backbone_id,char **master_atom_name_array,double **bb_pdf_entries,double *bb_pdf_dx)
{

    FILE *fp;
    char file_path[cmax_length];
    int j,k,found,current_k;
    char D1[cmax_length],D2[cmax_length],D3[cmax_length],D4[cmax_length];
    double min,max,dphi;
    int interval,bins,b;
    double **hist,**ref_hist;
    int conf1,conf2,conf3,conf4;
    double xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl,cross,norm1,norm2,sign,coeff,phi;
    double sum;
    
    min=-180.0;
    max= 180.0;
    dphi=  0.5;
    
    interval=max-min;
    bins=floor(interval/dphi+0.5);
    
    hist=(double**)malloc(bins*sizeof(double*));
    ref_hist=(double**)malloc(bins*sizeof(double*));
    for(j=0;j<bins;++j)hist[j]=(double*)malloc(backbones*sizeof(double));
    for(j=0;j<bins;++j)ref_hist[j]=(double*)malloc(backbones*sizeof(double));
    for(j=0;j<bins;++j)for(k=0;k<backbones;++k)hist[j][k]=0.0;
    /*
     for(j=0;j<backbones;++j)
     {
     printf("[%d]\tbackbone type:\t%s\n",j+1,backbone_id[j]);
     }
     */
    for(j=0;j<general_dihedrals;++j)
    {
        found=0;
        for(k=0;k<backbones;++k)
        {
            sscanf(backbone_id[k],"%s\t%s\t%s\t%s",D1,D2,D3,D4);
            if((strcmp(D1,master_atom_name_array[general_dihedrals_registry[j][1]-1])==0 &&
                strcmp(D2,master_atom_name_array[general_dihedrals_registry[j][2]-1])==0 &&
                strcmp(D3,master_atom_name_array[general_dihedrals_registry[j][3]-1])==0 &&
                strcmp(D4,master_atom_name_array[general_dihedrals_registry[j][4]-1])==0)
               ||
               (strcmp(D1,master_atom_name_array[general_dihedrals_registry[j][4]-1])==0 &&
                strcmp(D2,master_atom_name_array[general_dihedrals_registry[j][3]-1])==0 &&
                strcmp(D3,master_atom_name_array[general_dihedrals_registry[j][2]-1])==0 &&
                strcmp(D4,master_atom_name_array[general_dihedrals_registry[j][1]-1])==0))
            {
                found=1;current_k=k;break;
            }
        }
        if(found==1){
            
            conf1=general_dihedrals_registry[j][1];
            conf2=general_dihedrals_registry[j][2];
            conf3=general_dihedrals_registry[j][3];
            conf4=general_dihedrals_registry[j][4];
            
            xi=master_x_array[conf1-1];yi=master_y_array[conf1-1];zi=master_z_array[conf1-1];
            xj=master_x_array[conf2-1];yj=master_y_array[conf2-1];zj=master_z_array[conf2-1];
            xk=master_x_array[conf3-1];yk=master_y_array[conf3-1];zk=master_z_array[conf3-1];
            xl=master_x_array[conf4-1];yl=master_y_array[conf4-1];zl=master_z_array[conf4-1];
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
            
            b=floor((phi-min)/dphi);
            hist[b][current_k]=hist[b][current_k]+1.0;
            /*
             printf("[%d]\t{%d}\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%lf\n",
             j+1,current_k+1,
             general_dihedrals_registry[j][1],general_dihedrals_registry[j][2],general_dihedrals_registry[j][3],general_dihedrals_registry[j][4],
             master_atom_name_array[general_dihedrals_registry[j][1]-1],
             master_atom_name_array[general_dihedrals_registry[j][2]-1],
             master_atom_name_array[general_dihedrals_registry[j][3]-1],
             master_atom_name_array[general_dihedrals_registry[j][4]-1],
             phi);
             */
        }
    }
    
    for(j=0;j<backbones;++j)
    {
        sum=0.0;
        for(k=0;k<bins-1;++k)sum=sum+0.5*dphi*(hist[k][j]+hist[k+1][j]);
        for(k=0;k<bins;++k)hist[k][j]=hist[k][j]/sum;
    }
    
    for(j=0;j<backbones;++j)
    {
        for(k=0;k<bins;++k)
        {
            phi=min+((k+1)-0.5)*dphi;
            b=(int)floor((phi-min)/bb_pdf_dx[j]);
            //printf("%lf\t[%lf\t%lf]\n",phi,bb_pdf_entries[b][j],bb_pdf_entries[b+1][j]);
            ref_hist[k][j]=bb_pdf_entries[b][1+j*2]+(bb_pdf_entries[b+1][1+j*2]-bb_pdf_entries[b][1+j*2])*(phi-bb_pdf_entries[b][0+j*2])/bb_pdf_dx[j];
            
        }
    }
    
    for(j=0;j<backbones;++j)
    {
        sum=0.0;
        for(k=0;k<bins-1;++k)sum=sum+0.5*dphi*(ref_hist[k][j]+ref_hist[k+1][j]);
        for(k=0;k<bins;++k)ref_hist[k][j]=ref_hist[k][j]/sum;
    }
    
    sprintf(file_path,"%s/%s",current_folder,"backbone_statistics.dat");
    fp=fopen(file_path,"w+");
    for(j=0;j<bins;++j)
    {
        fprintf(fp,"%lf\t",min+((j+1)-0.5)*dphi);
        for(k=0;k<backbones;++k)fprintf(fp,"%lf\t%lf\t",hist[j][k],ref_hist[j][k]);fprintf(fp,"\n");
    }
    fclose(fp);
    
    for(j=0;j<bins;++j)free(hist[j]);free(hist);
    for(j=0;j<bins;++j)free(ref_hist[j]);free(ref_hist);

}
