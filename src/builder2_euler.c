
#include"builder.h"

void tokenize(char *buffer,int *ignore_out,int *elements_out);

void rotate(double ux, double uy, double uz, double x_center, double y_center, double z_center, double theta, double x_old, double y_old, double z_old, double *x_new, double *y_new, double *z_new);

void builder2_euler(char *current_folder,char *argv1,int atoms,int molecules,
                    int *com_flag_array,int *master_mol_array,int *mol_id_array,int *master_id_internal_array,char **master_species_array,
                    double *X,double *Y,double *Z,double *alpha,double *beta,double *gamma,
                    double *master_x_array,double *master_y_array,double *master_z_array,int reshuffle_euler,int reshuffle_coms,double rand_disp,
                    double xlo,double ylo,double zlo,double xhi,double yhi,double zhi)
{
    FILE *fp;
    char file_path[cmax_length],buffer[cmax_length],word[cmax_length];
    int i,j,k,l,m;
    double M,xcom,ycom,zcom,mass;                           // CoM related variables
    int *com_id_array;                                      // array to hold atomic ids for CoM calculation
    int ignore,elements;                                    // tokenize output variables
    int start,stop,*pos;                                    // tokenize discretization variables

    int i_rand;
    double alpha_rs,beta_rs,gamma_rs,d_rand,disp;
    double rand_alpha_min,rand_alpha_max,rand_beta_min,rand_beta_max,rand_gamma_min,rand_gamma_max;
    
    if(reshuffle_euler==1)
    {
        printf("$ reshuffled Euler angles:\n");
        printf("mol\talpha\tbeta\tgamma\n");

        
        
        rand_alpha_min=0.0;rand_alpha_max=360.0;
        rand_beta_min=0.0;rand_beta_max=360.0;
        rand_gamma_min=0.0;rand_gamma_max=360.0;
        
        for(j=0;j<molecules;++j){
        
            alpha_rs=0.0;
            beta_rs=0.0;
            gamma_rs=0.0;
            
            i_rand=rand();
            d_rand=(double)i_rand/RAND_MAX;
            disp=rand_alpha_min+d_rand*(rand_alpha_max-rand_alpha_min);
            alpha_rs=alpha_rs+disp;

            i_rand=rand();
            d_rand=(double)i_rand/RAND_MAX;
            disp=rand_beta_min+d_rand*(rand_beta_max-rand_beta_min);
            beta_rs=beta_rs+disp;
       
            i_rand=rand();
            d_rand=(double)i_rand/RAND_MAX;
            disp=rand_gamma_min+d_rand*(rand_gamma_max-rand_gamma_min);
            gamma_rs=gamma_rs+disp;
            
            alpha[j]=alpha_rs;
            beta[j]=beta_rs;
            gamma[j]=gamma_rs;
            
            printf("%d\t%lf\t%lf\t%lf\n",j+1,alpha[j],beta[j],gamma[j]);
            
        }
    }
    
    if(reshuffle_coms==1)
    {
        printf("$ reshuffled CoMs:\n");
        printf("mol\tCoM_X\tCoM_Y\tCoM_z\n");
        
        for(j=0;j<molecules;++j)
        {
            i_rand=rand();
            d_rand=(double)i_rand/RAND_MAX;
            disp=-rand_disp+d_rand*(2.0*rand_disp);
            X[j]=X[j]+disp;
            if(X[j]>xhi)X[j]=X[j]-(xhi-xlo);
            if(X[j]<xlo)X[j]=X[j]+(xhi-xlo);
            
            i_rand=rand();
            d_rand=(double)i_rand/RAND_MAX;
            disp=-rand_disp+d_rand*(2.0*rand_disp);
            Y[j]=Y[j]+disp;
            if(Y[j]>yhi)Y[j]=Y[j]-(yhi-ylo);
            if(Y[j]<ylo)Y[j]=Y[j]+(yhi-ylo);
            
            i_rand=rand();
            d_rand=(double)i_rand/RAND_MAX;
            disp=-rand_disp+d_rand*(2.0*rand_disp);
            Z[j]=Z[j]+disp;
            if(Z[j]>zhi)Z[j]=Z[j]-(zhi-zlo);
            if(Z[j]<zlo)Z[j]=Z[j]+(zhi-zlo);
            
            printf("%d\t%lf\t%lf\t%lf\n",j+1,X[j],Y[j],Z[j]);
            
        }
        
        
    }
    
    // calculate CoMs and apply Euler rotations
    for(i=0;i<molecules;++i)
    {
        if(com_flag_array[i]==0)
        {
            M=0.0;xcom=0.0;ycom=0.0;zcom=0.0;
            for(j=0;j<atoms;++j)
            {
                if(master_mol_array[j]==mol_id_array[i])
                {
                    if(master_species_array[j][0]=='H')mass=mass_H;
                    else if(master_species_array[j][0]=='C')mass=mass_C;
                    else if(master_species_array[j][0]=='O')mass=mass_O;
                    else if(master_species_array[j][0]=='N')mass=mass_N;
                    else{printf("error with masses!!\n");exit(-2);}
                    M=M+mass;
                    xcom=xcom+master_x_array[j]*mass;
                    ycom=ycom+master_y_array[j]*mass;
                    zcom=zcom+master_z_array[j]*mass;
                }
            }
            xcom=xcom/M;ycom=ycom/M;zcom=zcom/M;
        }
        else
        {
            com_id_array=(int*)malloc(com_flag_array[i]*sizeof(int));
            sprintf(file_path,"%s/%s",current_folder,argv1);
            fp=fopen(file_path,"r");
            if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
            while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"molecules")==0)break;}
            for(j=0;j<molecules;++j)
            {
                fgets(buffer,cmax_length,fp);
                if(i==j)
                {
                    tokenize(buffer,&ignore,&elements);
                    
                    pos=(int*)malloc(elements*sizeof(int));
                    
                    l=-1;
                    for(k=0;k<cmax_length;++k)
                    {
                        if(buffer[k]=='\0')break;
                        if(buffer[k]=='$')
                        {
                            l=l+1;
                            pos[l]=k;
                        }
                    }
                    for(k=9;k<9+com_flag_array[i];++k)
                    {
                        start=pos[k-1]+1;
                        stop=pos[k]-1;
                        sprintf(word,"%s","");
                        m=-1;
                        for(l=start;l<=stop;++l)
                        {
                            m=m+1;
                            
                            word[m]=buffer[l];
                        }
                        m=m+1;word[m]='\0';
                        com_id_array[k-9]=atoi(word);
                    }
                    free(pos);
                    break;
                }
            }
            fclose(fp);
            //
            M=0.0;xcom=0.0;ycom=0.0;zcom=0.0;
            for(j=0;j<atoms;++j)
            {
                if(master_mol_array[j]==mol_id_array[i])
                {
                    for(k=0;k<com_flag_array[i];++k)
                    {
                        if(com_id_array[k]==master_id_internal_array[j])
                        {
                            if(master_species_array[j][0]=='H')mass=mass_H;
                            else if(master_species_array[j][0]=='C')mass=mass_C;
                            else if(master_species_array[j][0]=='O')mass=mass_O;
                            else if(master_species_array[j][0]=='N')mass=mass_N;
                            else{printf("error with masses!!\n");exit(-2);}
                            M=M+mass;
                            xcom=xcom+master_x_array[j]*mass;
                            ycom=ycom+master_y_array[j]*mass;
                            zcom=zcom+master_z_array[j]*mass;
                            break;
                        }
                    }
                }
            }
            xcom=xcom/M;ycom=ycom/M;zcom=zcom/M;
            //
            free(com_id_array);
        }
        for(j=0;j<atoms;++j)
        {
            if(master_mol_array[j]==mol_id_array[i])
            {
                master_x_array[j]=master_x_array[j]-xcom+X[i];
                master_y_array[j]=master_y_array[j]-ycom+Y[i];
                master_z_array[j]=master_z_array[j]-zcom+Z[i];
            }
        }
        xcom=X[i];ycom=Y[i];zcom=Z[i];
        // rotate
        for(j=0;j<atoms;++j)
        {
            if(master_mol_array[j]==mol_id_array[i])
            {
                rotate(0.0,0.0,1.0,xcom,ycom,zcom,alpha[i],master_x_array[j],master_y_array[j],master_z_array[j],&master_x_array[j],&master_y_array[j],&master_z_array[j]);
                rotate(1.0,0.0,0.0,xcom,ycom,zcom,beta[i],master_x_array[j],master_y_array[j],master_z_array[j],&master_x_array[j],&master_y_array[j],&master_z_array[j]);
                rotate(0.0,0.0,1.0,xcom,ycom,zcom,gamma[i],master_x_array[j],master_y_array[j],master_z_array[j],&master_x_array[j],&master_y_array[j],&master_z_array[j]);
            }
        }
    }

}

