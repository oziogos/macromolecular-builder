
#include"builder.h"

void calculate_tacticity(int molecules, int general_t_counter, int tacticity_array_length,
                         int *atom_scaling_array, double *master_x_array, double *master_y_array, double *master_z_array,
                         int *tacticity_counter_array, int *tacticity_counter_array_v2, int *yet_one_more_t_array,
                         int *t_tracker_1, int *t_tracker_2, int *t_tracker_3, int *t_tracker_4, int *t_tracker_mol,
                         double *tacticity_array,int counter,char *current_folder)
{
    int j;
    double xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl,coeff;
    
    FILE *fp;
    char file_path[cmax_length];
    
    for(j=0;j<molecules;++j)tacticity_counter_array_v2[j]=-1;
    
    //for(j=0;j<molecules;++j)printf("[%d]\t%d\n",j+1,tacticity_counter_array[j]);
    
    //for(j=0;j<general_t_counter;++j)printf("[%d]\t%d belonging to mol %d\n",j+1,t_tracker_2[j]+atom_scaling_array[t_tracker_mol[j]-1],t_tracker_mol[j]);
    
    for(j=0;j<general_t_counter;++j)
    {
        tacticity_counter_array_v2[t_tracker_mol[j]-1]=tacticity_counter_array_v2[t_tracker_mol[j]-1]+1;
        /*
        printf("---------------- molecule %d\t\t%d\n",t_tracker_mol[j],yet_one_more_t_array[t_tracker_mol[j]-1]+tacticity_counter_array_v2[t_tracker_mol[j]-1]);
        printf("[%d]\n",t_tracker_1[j]+atom_scaling_array[t_tracker_mol[j]-1]);
        printf(" |\n");
        printf("[%d]-[%d]\n",t_tracker_2[j]+atom_scaling_array[t_tracker_mol[j]-1],t_tracker_4[j]+atom_scaling_array[t_tracker_mol[j]-1]);
        printf(" |\n");
        printf("[%d]\n",t_tracker_3[j]+atom_scaling_array[t_tracker_mol[j]-1]);
        */
        // recalc here
        
        //-- mixed product calculator --
        
        // update tacticity vector
        xi=master_x_array[t_tracker_1[j]+atom_scaling_array[t_tracker_mol[j]-1]-1];
        yi=master_y_array[t_tracker_1[j]+atom_scaling_array[t_tracker_mol[j]-1]-1];
        zi=master_z_array[t_tracker_1[j]+atom_scaling_array[t_tracker_mol[j]-1]-1];
        xj=master_x_array[t_tracker_2[j]+atom_scaling_array[t_tracker_mol[j]-1]-1];
        yj=master_y_array[t_tracker_2[j]+atom_scaling_array[t_tracker_mol[j]-1]-1];
        zj=master_z_array[t_tracker_2[j]+atom_scaling_array[t_tracker_mol[j]-1]-1];
        xk=master_x_array[t_tracker_3[j]+atom_scaling_array[t_tracker_mol[j]-1]-1];
        yk=master_y_array[t_tracker_3[j]+atom_scaling_array[t_tracker_mol[j]-1]-1];
        zk=master_z_array[t_tracker_3[j]+atom_scaling_array[t_tracker_mol[j]-1]-1];
        xl=master_x_array[t_tracker_4[j]+atom_scaling_array[t_tracker_mol[j]-1]-1];
        yl=master_y_array[t_tracker_4[j]+atom_scaling_array[t_tracker_mol[j]-1]-1];
        zl=master_z_array[t_tracker_4[j]+atom_scaling_array[t_tracker_mol[j]-1]-1];
        coeff=(yl-yj)*(xj*zi-xk*zi-xi*zj+xk*zj+xi*zk-xj*zk)+(xl-xj)*(-yj*zi+yk*zi+yi*zj-yk*zj-yi*zk+yj*zk)+(-xj*yi+xk*yi+xi*yj-xk*yj-xi*yk+xj*yk)*(zl-zj);
        
        tacticity_array[yet_one_more_t_array[t_tracker_mol[j]-1]+tacticity_counter_array_v2[t_tracker_mol[j]-1]]=coeff;
        
        //tacticity_array[tacticity_counter_array[i]-1+yet_one_more_t_array[i]]=coeff;
        
        
        //--
        
        //
        
    }

        
        sprintf(file_path,"%s/tacticity_history.dat",current_folder);
        fp=fopen(file_path,"a");
        fprintf(fp,"State: %d\n",counter);
        for(j=0;j<tacticity_array_length;++j)fprintf(fp,"[%d]\t%lf\n",j,tacticity_array[j]);
        fclose(fp);
    
    //for(j=0;j<tacticity_array_length;++j)printf("[%d]\t%lf\n",j,tacticity_array[j]);
}
