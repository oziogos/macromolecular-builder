
#include"builder.h"

void calculate_Cn(char *current_folder, int counter, int molecules, int general_bonds, int backbones,
                  int *max_growth_step_array, int *atom_scaling_array, int *master_B1_core_array, int *master_B2_core_array, int *mol_id_array_BONDS,
                  double *master_x_array, double *master_y_array, double *master_z_array,char **master_atom_name_array, char **backbone_id,
                  double *l_array, double *l2_array, int *Re_start, int *Re_stop, int *bb_growth_step_array,
                  int backbone_type_counter, char **backbone_type_array)
{
    FILE *fp;
    int i,j,k,found;
    double norm2,Re2,Cn;
    char file_path[cmax_length];
    double sum;int N;
    
    for(i=0;i<molecules;++i){l_array[i]=0.0;l2_array[i]=0.0;}
    
    for(j=0;j<general_bonds;++j)
    {
        if(strcmp(master_atom_name_array[master_B1_core_array[j]-1],master_atom_name_array[master_B2_core_array[j]-1])==0)
        {
            found=0;
            /*
            for(k=0;k<backbones;++k)
            {
                if(strcmp(master_atom_name_array[master_B1_core_array[j]-1],backbone_id[k])==0)
                {
                    found=1;
                    break;
                }
            }
            */
            for(k=0;k<backbone_type_counter;++k)
            {
                if(strcmp(master_atom_name_array[master_B1_core_array[j]-1],backbone_type_array[k])==0)
                {
                    found=1;
                    break;
                }
            }
            if(found==1)
            {
                norm2=
                (master_x_array[master_B1_core_array[j]-1]-master_x_array[master_B2_core_array[j]-1])*
                (master_x_array[master_B1_core_array[j]-1]-master_x_array[master_B2_core_array[j]-1])+
                (master_y_array[master_B1_core_array[j]-1]-master_y_array[master_B2_core_array[j]-1])*
                (master_y_array[master_B1_core_array[j]-1]-master_y_array[master_B2_core_array[j]-1])+
                (master_z_array[master_B1_core_array[j]-1]-master_z_array[master_B2_core_array[j]-1])*
                (master_z_array[master_B1_core_array[j]-1]-master_z_array[master_B2_core_array[j]-1]);
                l_array[mol_id_array_BONDS[j]-1]=l_array[mol_id_array_BONDS[j]-1]+sqrt(norm2);
                l2_array[mol_id_array_BONDS[j]-1]=l2_array[mol_id_array_BONDS[j]-1]+norm2;
            }
        }
    }
    
    sprintf(file_path,"%s/Cn_history.dat",current_folder);
    fp=fopen(file_path,"a");
    
    fprintf(fp,"State: %d\n",counter);
    
    sum=0.0;N=0;
    for(i=0;i<molecules;++i)
    {
        if(max_growth_step_array[i]!=0 && bb_growth_step_array[i]!=0)// && added_backbone_flag[i]==1)
        {
            
            Re2=
            (master_x_array[Re_start[i]+atom_scaling_array[i]-1]-master_x_array[Re_stop[i]+atom_scaling_array[i]-1])*
            (master_x_array[Re_start[i]+atom_scaling_array[i]-1]-master_x_array[Re_stop[i]+atom_scaling_array[i]-1])+
            (master_y_array[Re_start[i]+atom_scaling_array[i]-1]-master_y_array[Re_stop[i]+atom_scaling_array[i]-1])*
            (master_y_array[Re_start[i]+atom_scaling_array[i]-1]-master_y_array[Re_stop[i]+atom_scaling_array[i]-1])+
            (master_z_array[Re_start[i]+atom_scaling_array[i]-1]-master_z_array[Re_stop[i]+atom_scaling_array[i]-1])*
            (master_z_array[Re_start[i]+atom_scaling_array[i]-1]-master_z_array[Re_stop[i]+atom_scaling_array[i]-1]);
            
            l_array[i]=l_array[i]/bb_growth_step_array[i];
            l2_array[i]=l2_array[i]/bb_growth_step_array[i];
            
            Cn=Re2/(bb_growth_step_array[i]*l2_array[i]);
            
            printf("[molecule %d]\tR^2 = %lf A^2\tn = %d\t<l> = %lf A\t<l2> = %lf A^2\tCn = %lf\n",i+1,Re2,bb_growth_step_array[i],l_array[i],l2_array[i],Cn);

            //printf("%lf\t%d\t%lf\n",Re2,bb_growth_step_array[i],l2_array[i]);
            
            fprintf(fp,"[molecule %d]\tR^2 = %lf A^2\tn = %d\t<l> = %lf A\t<l2> = %lf A^2\tCn = %lf\n",i+1,Re2,bb_growth_step_array[i],l_array[i],l2_array[i],Cn);

            N=N+1;
            sum=sum+Cn;
            
        }
    }
    sum=sum/N;
    
    //getchar();
    
    fclose(fp);
    
    sprintf(file_path,"%s/Cn_vs_build.dat",current_folder);
    fp=fopen(file_path,"a");
    fprintf(fp,"%d\t%lf\n",counter,sum);
    fclose(fp);

    //getchar();

}
