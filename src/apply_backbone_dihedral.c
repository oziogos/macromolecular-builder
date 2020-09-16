//
// 10/12/2016: Change from atom name backbone identification scheme to dihedral type identification scheme
//
#include"builder.h"

void rotate(double ux, double uy, double uz, double x_center, double y_center, double z_center, double theta, double x_old, double y_old, double z_old, double *x_new, double *y_new, double *z_new);

void apply_backbone_dihedral(int i, int backbones, int local_atoms,int *local_topo_array, int *local_map, char **segment_atom_name_core, char **backbone_id, int *atom_scaling_array,
                             double **bb_entries,double *bb_dx,
                             double *master_x_array,double *master_y_array,double *master_z_array,
                             //int *Re_stop,int *bb_growth_step_array, int *prev_link_atom_array,int *added_backbone_flag,
                             int backbone_type_counter, char **backbone_type_array
                             )
{
    
    int found,j,k,conf1,conf2,conf3,conf4,f1,f2,f3,bb_i,bb_index;//,for_rotation;
    double xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl,cross,norm1,norm2,sign,coeff,phi,bb_z,bb_zmap,vx,vy,vz;
    char D1[sub_length],D2[sub_length],D3[sub_length],D4[sub_length];
    char D1_current[sub_length],D2_current[sub_length],D3_current[sub_length],D4_current[sub_length];
    int identified_D=0,is_bb;
    
    // new course of action
    // apply the backward from the GS search using backbone_type_counter and backbone_type_array
    // once we establish that we have a backbone proper dihedral angle, loop on backbone_id to resolve type
    
    // flag
    found=0;
    // loop on local atoms
    for(j=0;j<local_atoms;++j)
    {
        // loop on backbone types
        for(k=0;k<backbone_type_counter;++k)
        {
            // -1 corresponds to the added atom: growth site (GS)
            // found is set equal to 1 if the GS is a backbone atom
            if(local_topo_array[local_map[j]-1]==-1 && strcmp(segment_atom_name_core[local_map[j]-1],backbone_type_array[k])==0)
            {
                found=1;
                break;
            }
        }
        if(found==1)break;
    }
    // if the GS is a backbone atom, we examine the rest local atoms
    // -- the dihedral atoms search is done in a backwards fashion!!
    if(found==1)
    {
        // we store the id of the GS in conf1, before index j is overwriten
        conf1=local_map[j]+atom_scaling_array[i];
        
        // store the unscaled by local_map id for the end-to-end distance calculation for molecule i
        //Re_stop[i]=local_map[j];
        
        // f1 is the flag for the host atom
        f1=0;
        // loop on local atoms
        for(j=0;j<local_atoms;++j)
        {
            
            // loop on backbone types
            for(k=0;k<backbone_type_counter;++k)
            {
                // 1 corresponds to the host atom
                // if the type is consistent with the k backbone type index, set f1 equal to 1 and store the atomic id in conf2
                if(local_topo_array[local_map[j]-1]==1 && strcmp(segment_atom_name_core[local_map[j]-1],backbone_type_array[k])==0)
                {
                    f1=1;conf2=local_map[j]+atom_scaling_array[i];break;
                }
            }
            if(f1==1)break;
        }
        /*
        // Cn calculation arrays
        // I believe prev_link_atom_array and added_backbone_flag are obsolete...
        if(f1==1)
        {
            bb_growth_step_array[i]=bb_growth_step_array[i]+1;
            prev_link_atom_array[i]=local_map[j];
            added_backbone_flag[i]=1;
            //printf("$$$ added a backbone link: %d-%d\n",conf1,conf2);
            //getchar();
        }
        //
        */
        // f2 is the flag for the 1-2 neighbor
        f2=0;
        for(j=0;j<local_atoms;++j)
        {
            // loop on backbone types
            for(k=0;k<backbone_type_counter;++k)
            {
                // 2 is for the 1-2 neighbor
                // if the type is consistent with the k backbone type index, set f2 equal to 1 and store the atomic id in conf3
                if(local_topo_array[local_map[j]-1]==2 && strcmp(segment_atom_name_core[local_map[j]-1],backbone_type_array[k])==0)
                {
                    f2=1;conf3=local_map[j]+atom_scaling_array[i];break;
                }
            }
        }
        
        // f3 is the flag for the 1-3 neighbor
        f3=0;
        for(j=0;j<local_atoms;++j)
        {
            // loop on backbone types
            for(k=0;k<backbone_type_counter;++k)
            {
                // 3 is for the 1-3 neighbor
                // if the type is consistent with the k backbone type index, set f3 equal to 1 and store the atomic id in conf4
                if(local_topo_array[local_map[j]-1]==3 && strcmp(segment_atom_name_core[local_map[j]-1],backbone_type_array[k])==0)
                {
                    f3=1;conf4=local_map[j]+atom_scaling_array[i];break;
                }
            }
        }
        
        // we are already in a found==1 if block
        // so, if f1==f2==f3==1, we have located a backbone dihedral angle
        if(f1==1 && f2==1 && f3==1)
        {
            // we now need to resolve the backbone dihedral type!
            
            //printf("found a backbone dihedral angle!\n");
            
            sprintf(D1_current,"%s",segment_atom_name_core[conf1-atom_scaling_array[i]-1]);
            sprintf(D2_current,"%s",segment_atom_name_core[conf2-atom_scaling_array[i]-1]);
            sprintf(D3_current,"%s",segment_atom_name_core[conf3-atom_scaling_array[i]-1]);
            sprintf(D4_current,"%s",segment_atom_name_core[conf4-atom_scaling_array[i]-1]);
            
            //printf("%s\t%s\t%s\t%s\n",D1_current,D2_current,D3_current,D4_current);//getchar();
            
            for(k=0;k<backbones;++k)
            {
                sscanf(backbone_id[k],"%s\t%s\t%s\t%s",D1,D2,D3,D4);
                if((strcmp(D1,D1_current)==0 && strcmp(D2,D2_current)==0 && strcmp(D3,D3_current)==0 && strcmp(D4,D4_current)==0)
                   ||
                   (strcmp(D1,D4_current)==0 && strcmp(D2,D3_current)==0 && strcmp(D3,D2_current)==0 && strcmp(D4,D1_current)==0))
                {
                    identified_D=1;
                    break;
                
                }
            }
            
            //printf("%d\n",identified_D);getchar();
            
            if(identified_D==1)
            {
                // calculate the proper dihedral angle conf1-conf2-conf3-conf4
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
                //printf("dihedral: %lf\n",phi);
                
                // store k type index in bb_i
                bb_i=k;
                // draw random number between 0 and 1: bb_z
                bb_z=(double)rand()/RAND_MAX;
                // select dihedral bb_zmap from the type bb_i invCDF
                bb_index=(int)floor(bb_z/bb_dx[bb_i]);
                bb_zmap=bb_entries[bb_index][1+bb_i*2]+(bb_entries[bb_index+1][1+bb_i*2]-bb_entries[bb_index][1+bb_i*2])*(bb_z-bb_entries[bb_index][0+bb_i*2])/bb_dx[bb_i];
                //bb_zmap=180.0;
                //printf("new angle: %lf\n",bb_zmap);
                
                // rotation vector
                vx=master_x_array[conf2-1]-master_x_array[conf3-1];
                vy=master_y_array[conf2-1]-master_y_array[conf3-1];
                vz=master_z_array[conf2-1]-master_z_array[conf3-1];
                norm1=sqrt(vx*vx+vy*vy+vz*vz);
                vx=vx/norm1;vy=vy/norm1;vz=vz/norm1;
                //printf("v vector: %lf\t%lf\t%lf\n",vx,vy,vz);
                
                // counter for debug purposes only
                //for_rotation=0;
                // loop on local atoms, seek rotatable atoms and apply rotation
                for(j=0;j<local_atoms;++j)
                {
                    // rotate the GS
                    if(local_topo_array[local_map[j]-1]==-1)
                    {
                        //for_rotation=for_rotation+1;
                        //printf("### -1 rot %d %s\n",local_map[j]+atom_scaling_array[i],segment_atom_name_core[local_map[j]-1]);
                        rotate(vx,vy,vz,
                               master_x_array[conf3-1],master_y_array[conf3-1],master_z_array[conf3-1],
                               -phi+bb_zmap,
                               master_x_array[local_map[j]+atom_scaling_array[i]-1],
                               master_y_array[local_map[j]+atom_scaling_array[i]-1],
                               master_z_array[local_map[j]+atom_scaling_array[i]-1],
                               &master_x_array[local_map[j]+atom_scaling_array[i]-1],
                               &master_y_array[local_map[j]+atom_scaling_array[i]-1],
                               &master_z_array[local_map[j]+atom_scaling_array[i]-1]);
                    }
                    // rotate the newly added hydrogen atoms
                    if(local_topo_array[local_map[j]-1]==-2)
                    {
                        //for_rotation=for_rotation+1;
                        //printf("### -2 rot %d %s\n",local_map[j]+atom_scaling_array[i],segment_atom_name_core[local_map[j]-1]);
                        rotate(vx,vy,vz,
                               master_x_array[conf3-1],master_y_array[conf3-1],master_z_array[conf3-1],
                               -phi+bb_zmap,
                               master_x_array[local_map[j]+atom_scaling_array[i]-1],
                               master_y_array[local_map[j]+atom_scaling_array[i]-1],
                               master_z_array[local_map[j]+atom_scaling_array[i]-1],
                               &master_x_array[local_map[j]+atom_scaling_array[i]-1],
                               &master_y_array[local_map[j]+atom_scaling_array[i]-1],
                               &master_z_array[local_map[j]+atom_scaling_array[i]-1]);
                    }
                    // rotate the rest 1-2 atoms
                    //if(local_topo_array[local_map[j]-1]==2 && strcmp(segment_atom_name_core[local_map[j]-1],backbone_id[k])!=0)
                    if(local_topo_array[local_map[j]-1]==2)
                    {
                        is_bb=0;
                        for(k=0;k<backbone_type_counter;++k)
                        {
                            if(strcmp(segment_atom_name_core[local_map[j]-1],backbone_type_array[k])==0)
                            {
                                is_bb=1;
                                break;
                            }
                        }
                        if(is_bb==0)
                        {
                        //for_rotation=for_rotation+1;
                        //printf("### 2 rot %d %s\n",local_map[j]+atom_scaling_array[i],segment_atom_name_core[local_map[j]-1]);
                        rotate(vx,vy,vz,
                               master_x_array[conf3-1],master_y_array[conf3-1],master_z_array[conf3-1],
                               -phi+bb_zmap,
                               master_x_array[local_map[j]+atom_scaling_array[i]-1],
                               master_y_array[local_map[j]+atom_scaling_array[i]-1],
                               master_z_array[local_map[j]+atom_scaling_array[i]-1],
                               &master_x_array[local_map[j]+atom_scaling_array[i]-1],
                               &master_y_array[local_map[j]+atom_scaling_array[i]-1],
                               &master_z_array[local_map[j]+atom_scaling_array[i]-1]);
                        }
                    }
                }
            }
            //getchar();
            
        }
    }
}
