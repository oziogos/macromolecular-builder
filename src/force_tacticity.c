
#include"builder.h"

void force_tacticity(int i, int tacticity, int *mol_type_array, char **tacticity_search, int local_atoms, char **segment_atom_name_core, int *local_map, int *local_topo_array,
                     int *atom_scaling_array, char **backbone_id, int *yet_one_more_t_array, int backbones,
                     int backbone_type_counter, char **backbone_type_array,
                     
                     int molecules,int tacticity_array_length,
                     
                     int *tacticity_counter_array,
                     int *general_t_counter, int *t_tracker_1, int *t_tracker_2, int *t_tracker_3, int *t_tracker_4,int *t_tracker_mol,
                     double *tacticity_array, double *master_x_array, double *master_y_array, double *master_z_array,
                     int *ft
                     )
{
    
    char word[cmax_length],tact_type[sub_length];
    int int_buffer,found,j,k,conf1,conf2,conf3,conf4,meso;
    double xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl,coeff,xsub,ysub,zsub,vx,vy,vz,norm,reflect[3][3],xb,yb,zb,bb_z;
    
    int jj;
    
    // console out
    //printf("working on molecule type %d\n",mol_type_array[i]);
    // read tacticity info for current molecular type
    sscanf(tacticity_search[mol_type_array[i]-1],"%s\t%d\t%s",word,&int_buffer,tact_type);
    // !=0 means that there is a tacticity to be forced
    if(int_buffer!=0)
    {
        // console out
        //printf("I must look for tacticity type %s\n",word);
        // search for the pendant group atom
        *ft=0;
        //found=0;
        for(j=0;j<local_atoms;++j)
        {
            if(strcmp(word,segment_atom_name_core[local_map[j]-1])==0 && local_topo_array[local_map[j]-1]==-1)
            {
                *ft=1;
                break;
            }
            
        }
        // if you find a pendant group:
        if(*ft==1)
        {
            
            /*
            //
            printf("$ prior:\n");
            printf("tacticity_counter_array:\t");
            for(jj=0;jj<molecules;++jj)printf("%d\t",tacticity_counter_array[jj]);printf("\n");
            printf("general_t_counter:\t");
            printf("%d\n",*general_t_counter);
            if(*general_t_counter==0){printf("t_tracker_1\tN\\A\n");printf("t_tracker_2\tN\\A\n");printf("t_tracker_3\tN\\A\n");printf("t_tracker_4\tN\\A\n");printf("t_tracker_mol\tN\\A\n");}
            else
            {
                printf("t_tracker_1:\t");for(jj=0;jj<*general_t_counter;++jj)printf("%d\t",t_tracker_1[jj]);printf("\n");
                printf("t_tracker_2:\t");for(jj=0;jj<*general_t_counter;++jj)printf("%d\t",t_tracker_2[jj]);printf("\n");
                printf("t_tracker_3:\t");for(jj=0;jj<*general_t_counter;++jj)printf("%d\t",t_tracker_3[jj]);printf("\n");
                printf("t_tracker_4:\t");for(jj=0;jj<*general_t_counter;++jj)printf("%d\t",t_tracker_4[jj]);printf("\n");
                printf("t_tracker_mol:\t");for(jj=0;jj<*general_t_counter;++jj)printf("%d\t",t_tracker_mol[jj]);printf("\n");
            }
            printf("tacticity_array:\t");for(jj=0;jj<tacticity_array_length;++jj)printf("%lf\t",tacticity_array[jj]);printf("\n");
            //
            */
            
            // store local-j index of pendant group to conf4
            conf4=j;
            // augment the tacticity counter for the current molecule
            tacticity_counter_array[i]=tacticity_counter_array[i]+1;
            //
            *general_t_counter=*general_t_counter+1;
            /*
            if(*general_t_counter==1)
            {
                *t_tracker_1=(int*)malloc(*general_t_counter*sizeof(int));
                *t_tracker_2=(int*)malloc(*general_t_counter*sizeof(int));
                *t_tracker_3=(int*)malloc(*general_t_counter*sizeof(int));
                *t_tracker_4=(int*)malloc(*general_t_counter*sizeof(int));
                *t_tracker_mol=(int*)malloc(*general_t_counter*sizeof(int));
            }
            else
            {
                *t_tracker_1=(int*)realloc(*t_tracker_1,*general_t_counter*sizeof(int));
                *t_tracker_2=(int*)realloc(*t_tracker_2,*general_t_counter*sizeof(int));
                *t_tracker_3=(int*)realloc(*t_tracker_3,*general_t_counter*sizeof(int));
                *t_tracker_4=(int*)realloc(*t_tracker_4,*general_t_counter*sizeof(int));
                *t_tracker_mol=(int*)realloc(*t_tracker_mol,*general_t_counter*sizeof(int));
            }
            */
            //
            t_tracker_4[*general_t_counter-1]=local_map[conf4];
            // console out
            //printf("this is the pendant group atom:\n");
            //printf("**** found %d %s\n",local_map[conf4]+atom_scaling_array[i],segment_atom_name_core[local_map[conf4]-1]);
            //for(k=0;k<local_atoms;++k)
            //    printf("%d\t%s\t%d\n",local_topo_array[local_map[k]-1],segment_atom_name_core[local_map[k]-1],local_map[k]+atom_scaling_array[i]);
            //printf("I will have to look for the backbone atom to which the pendant group %s is directly linked...\n",word);
            // look for the backbone atom to which the pendant atom is directly linked
            for(k=0;k<local_atoms;++k)
            {
                if(local_topo_array[local_map[k]-1]==1)
                {
                    break;
                }
            }
            // store linking backbone atom local-j index to conf2
            conf2=k;
            t_tracker_2[*general_t_counter-1]=local_map[conf2];
            t_tracker_mol[*general_t_counter-1]=i+1;
            // console out
            //printf("**** found %d %s\n",local_map[conf2]+atom_scaling_array[i],segment_atom_name_core[local_map[conf2]-1]);
            //printf("now I must look for the inbound backbone atom before the linking backbone atom...\n");
            // look for the inbound backbone atom
            for(k=0;k<local_atoms;++k)
            {
                if(local_topo_array[local_map[k]-1]==2)
                {
                    break;
                }
            }
            // store inbound backbone atom local-j index to conf1
            conf1=k;
            t_tracker_1[*general_t_counter-1]=local_map[conf1];
            // console out
            //printf("**** found %d %s\n",local_map[conf1]+atom_scaling_array[i],segment_atom_name_core[local_map[conf1]-1]);
            //printf("now I must look for the outbound backbone atom after the linking backbone atom...\n");
            // look for the outbound backbone atom
            for(k=0;k<local_atoms;++k)
            {
                found=0;
                /*
                for(j=0;j<backbones;++j)
                {
                    if(local_topo_array[local_map[k]-1]==-1 && strcmp(segment_atom_name_core[local_map[k]-1],backbone_id[j])==0)
                    {
                        found=1;
                        break;}
                }
                */
                for(j=0;j<backbone_type_counter;++j)
                {
                    if(local_topo_array[local_map[k]-1]==-1 && strcmp(segment_atom_name_core[local_map[k]-1],backbone_type_array[j])==0)
                    {
                        found=1;
                        break;}
                }
                if(found==1)break;
            }
            // store outbound backbone local-j index to conf3
            conf3=k;
            t_tracker_3[*general_t_counter-1]=local_map[conf3];
            // console out
            //printf("**** found %d %s\n",local_map[conf3]+atom_scaling_array[i],segment_atom_name_core[local_map[conf3]-1]);
            // mixed product:
            // (yl-yj)*(xj*zi-xk*zi-xi*zj+xk*zj+xi*zk-xj*zk)+(xl-xj)*(-yj*zi+yk*zi+yi*zj-yk*zj-yi*zk+yj*zk)+(-xj*yi+xk*yi+xi*yj-xk*yj-xi*yk+xj*yk)*(zl-zj)
            //
            // (i)-->(j)-->(k)
            //        |
            //       (l)
            // mixed product coords
            xi=master_x_array[local_map[conf1]+atom_scaling_array[i]-1];
            yi=master_y_array[local_map[conf1]+atom_scaling_array[i]-1];
            zi=master_z_array[local_map[conf1]+atom_scaling_array[i]-1];
            xj=master_x_array[local_map[conf2]+atom_scaling_array[i]-1];
            yj=master_y_array[local_map[conf2]+atom_scaling_array[i]-1];
            zj=master_z_array[local_map[conf2]+atom_scaling_array[i]-1];
            xk=master_x_array[local_map[conf3]+atom_scaling_array[i]-1];
            yk=master_y_array[local_map[conf3]+atom_scaling_array[i]-1];
            zk=master_z_array[local_map[conf3]+atom_scaling_array[i]-1];
            xl=master_x_array[local_map[conf4]+atom_scaling_array[i]-1];
            yl=master_y_array[local_map[conf4]+atom_scaling_array[i]-1];
            zl=master_z_array[local_map[conf4]+atom_scaling_array[i]-1];
            // store mixed product in coeff
            coeff=(yl-yj)*(xj*zi-xk*zi-xi*zj+xk*zj+xi*zk-xj*zk)+(xl-xj)*(-yj*zi+yk*zi+yi*zj-yk*zj-yi*zk+yj*zk)+(-xj*yi+xk*yi+xi*yj-xk*yj-xi*yk+xj*yk)*(zl-zj);
            // console out
            //printf("### mixed product: %lf\n",coeff);
            //printf("### tacticity type: %s\n",tact_type);
            //printf("t counter [%d] = %d\n",i+1,tacticity_counter_array[i]);
            //printf("go write at %d\n",tacticity_counter_array[i]-1+yet_one_more_t_array[i]);
            // save mixed product in memory
            tacticity_array[tacticity_counter_array[i]-1+yet_one_more_t_array[i]]=coeff;
            //printf("@@@\t%d\n",tacticity_counter_array[i]-1+yet_one_more_t_array[i]);
            // dyad check: >1
            if(tacticity_counter_array[i]>1)
            {
                // console out
                /*
                printf("$$$ my mixed product is %lf\n",tacticity_array[tacticity_counter_array[i]-1+yet_one_more_t_array[i]]);
                printf("$$$ my neighbor's mixed product is %lf\n",tacticity_array[-1+tacticity_counter_array[i]-1+yet_one_more_t_array[i]]);
                printf("$$$ atoms to define the reflection plane:\n");
                printf("$ %d(%s)\n",local_map[conf1]+atom_scaling_array[i],segment_atom_name_core[local_map[conf1]-1]);
                printf("$ %d(%s)\n",local_map[conf2]+atom_scaling_array[i],segment_atom_name_core[local_map[conf2]-1]);
                printf("$ %d(%s)\n",local_map[conf3]+atom_scaling_array[i],segment_atom_name_core[local_map[conf3]-1]);
                printf("$ atoms to reflect:\n");
                for(j=0;j<local_atoms;++j)
                {
                    if(local_topo_array[local_map[j]-1]==-1 ||
                       local_topo_array[local_map[j]-1]==-2 ||
                       (local_topo_array[local_map[j]-1]==2 && j!=conf1))
                    {
                        printf("--> %d(%s)\n",local_map[j]+atom_scaling_array[i],segment_atom_name_core[local_map[j]-1]);
                        
                    }
                    
                }
                */
                // resolve if the dyad is meso
                meso=0;
                if(tacticity_array[tacticity_counter_array[i]-1+yet_one_more_t_array[i]]*tacticity_array[-1+tacticity_counter_array[i]-1+yet_one_more_t_array[i]]>0.0)
                {
                    // console out
                    //printf("@@@ this is a MESO dyad\n");
                    // meso flag
                    meso=1;
                }
                else
                {
                    // console out
                    //printf("@@@ this is a RACEMO dyad\n");
                }
                // chech syndiotactic case
                if(strcmp(tact_type,"s")==0)
                {
                    // apply racemo only if meso==1
                    if(meso==1)
                    {
                        // bring first atom to axis origin: apply to all relevant atoms
                        // inbound coords
                        xsub=master_x_array[local_map[conf1]+atom_scaling_array[i]-1];
                        ysub=master_y_array[local_map[conf1]+atom_scaling_array[i]-1];
                        zsub=master_z_array[local_map[conf1]+atom_scaling_array[i]-1];
                        // reset inbound
                        master_x_array[local_map[conf1]+atom_scaling_array[i]-1]=master_x_array[local_map[conf1]+atom_scaling_array[i]-1]-xsub;
                        master_y_array[local_map[conf1]+atom_scaling_array[i]-1]=master_y_array[local_map[conf1]+atom_scaling_array[i]-1]-ysub;
                        master_z_array[local_map[conf1]+atom_scaling_array[i]-1]=master_z_array[local_map[conf1]+atom_scaling_array[i]-1]-zsub;
                        // reset linking
                        master_x_array[local_map[conf2]+atom_scaling_array[i]-1]=master_x_array[local_map[conf2]+atom_scaling_array[i]-1]-xsub;
                        master_y_array[local_map[conf2]+atom_scaling_array[i]-1]=master_y_array[local_map[conf2]+atom_scaling_array[i]-1]-ysub;
                        master_z_array[local_map[conf2]+atom_scaling_array[i]-1]=master_z_array[local_map[conf2]+atom_scaling_array[i]-1]-zsub;
                        // reset outbound and atoms to be reflected
                        for(j=0;j<local_atoms;++j)
                        {
                            if(local_topo_array[local_map[j]-1]==-1 ||
                               local_topo_array[local_map[j]-1]==-2 ||
                               (local_topo_array[local_map[j]-1]==2 && j!=conf1))
                            {
                                master_x_array[local_map[j]+atom_scaling_array[i]-1]=master_x_array[local_map[j]+atom_scaling_array[i]-1]-xsub;
                                master_y_array[local_map[j]+atom_scaling_array[i]-1]=master_y_array[local_map[j]+atom_scaling_array[i]-1]-ysub;
                                master_z_array[local_map[j]+atom_scaling_array[i]-1]=master_z_array[local_map[j]+atom_scaling_array[i]-1]-zsub;
                            }
                        }
                        // console out
                        /*
                        printf("!!! reflection plane atoms:\n");
                        printf("$ %d(%s)\t%lf\t%lf\t%lf\n",local_map[conf1]+atom_scaling_array[i],segment_atom_name_core[local_map[conf1]-1],
                               master_x_array[local_map[conf1]+atom_scaling_array[i]-1],
                               master_y_array[local_map[conf1]+atom_scaling_array[i]-1],
                               master_z_array[local_map[conf1]+atom_scaling_array[i]-1]);
                        printf("$ %d(%s)\t%lf\t%lf\t%lf\n",local_map[conf2]+atom_scaling_array[i],segment_atom_name_core[local_map[conf2]-1],
                               master_x_array[local_map[conf2]+atom_scaling_array[i]-1],
                               master_y_array[local_map[conf2]+atom_scaling_array[i]-1],
                               master_z_array[local_map[conf2]+atom_scaling_array[i]-1]);
                        printf("$ %d(%s)\t%lf\t%lf\t%lf\n",local_map[conf3]+atom_scaling_array[i],segment_atom_name_core[local_map[conf3]-1],
                               master_x_array[local_map[conf3]+atom_scaling_array[i]-1],
                               master_y_array[local_map[conf3]+atom_scaling_array[i]-1],
                               master_z_array[local_map[conf3]+atom_scaling_array[i]-1]);
                        */
                        //-- mixed product calculator --
                        
                        // normal plane vector
                        // coords
                        xi=master_x_array[local_map[conf1]+atom_scaling_array[i]-1];
                        yi=master_y_array[local_map[conf1]+atom_scaling_array[i]-1];
                        zi=master_z_array[local_map[conf1]+atom_scaling_array[i]-1];
                        xj=master_x_array[local_map[conf2]+atom_scaling_array[i]-1];
                        yj=master_y_array[local_map[conf2]+atom_scaling_array[i]-1];
                        zj=master_z_array[local_map[conf2]+atom_scaling_array[i]-1];
                        xk=master_x_array[local_map[conf3]+atom_scaling_array[i]-1];
                        yk=master_y_array[local_map[conf3]+atom_scaling_array[i]-1];
                        zk=master_z_array[local_map[conf3]+atom_scaling_array[i]-1];
                        // vectors:
                        //    x       y       z
                        // (xi-xj) (yi-yj) (zi-zj)
                        // (xk-xj) (yk-yj) (zk-zj)
                        // cross product
                        vx=(yi-yj)*(zk-zj)-(zi-zj)*(yk-yj);
                        vy=(zi-zj)*(xk-xj)-(xi-xj)*(zk-zj);
                        vz=(xi-xj)*(yk-yj)-(yi-yj)*(xk-xj);
                        norm=sqrt(vx*vx+vy*vy+vz*vz);
                        vx=vx/norm;vy=vy/norm;vz=vz/norm;
                        // console out
                        //printf("!!! normalized plane vector:\n");
                        //printf("%lf\t%lf\t%lf\n",vx,vy,vz);
                        // define reflection matrix
                        reflect[0][0]=1.0-2.0*vx*vx;
                        reflect[0][1]=   -2.0*vx*vy;
                        reflect[0][2]=   -2.0*vx*vz;
                        reflect[1][0]=   -2.0*vy*vx;
                        reflect[1][1]=1.0-2.0*vy*vy;
                        reflect[1][2]=   -2.0*vy*vz;
                        reflect[2][0]=   -2.0*vz*vx;
                        reflect[2][1]=   -2.0*vz*vy;
                        reflect[2][2]=1.0-2.0*vz*vz;
                        // apply the reflection matrix
                        for(j=0;j<local_atoms;++j)
                        {
                            if(local_topo_array[local_map[j]-1]==-1 ||
                               local_topo_array[local_map[j]-1]==-2 ||
                               (local_topo_array[local_map[j]-1]==2 && j!=conf1))
                            {
                                xb=master_x_array[local_map[j]+atom_scaling_array[i]-1];
                                yb=master_y_array[local_map[j]+atom_scaling_array[i]-1];
                                zb=master_z_array[local_map[j]+atom_scaling_array[i]-1];
                                master_x_array[local_map[j]+atom_scaling_array[i]-1]=reflect[0][0]*xb + reflect[0][1]*yb + reflect[0][2]*zb;
                                master_y_array[local_map[j]+atom_scaling_array[i]-1]=reflect[1][0]*xb + reflect[1][1]*yb + reflect[1][2]*zb;
                                master_z_array[local_map[j]+atom_scaling_array[i]-1]=reflect[2][0]*xb + reflect[2][1]*yb + reflect[2][2]*zb;
                            }
                        }
                        // restore coords
                        // inbound
                        master_x_array[local_map[conf1]+atom_scaling_array[i]-1]=master_x_array[local_map[conf1]+atom_scaling_array[i]-1]+xsub;
                        master_y_array[local_map[conf1]+atom_scaling_array[i]-1]=master_y_array[local_map[conf1]+atom_scaling_array[i]-1]+ysub;
                        master_z_array[local_map[conf1]+atom_scaling_array[i]-1]=master_z_array[local_map[conf1]+atom_scaling_array[i]-1]+zsub;
                        // linking
                        master_x_array[local_map[conf2]+atom_scaling_array[i]-1]=master_x_array[local_map[conf2]+atom_scaling_array[i]-1]+xsub;
                        master_y_array[local_map[conf2]+atom_scaling_array[i]-1]=master_y_array[local_map[conf2]+atom_scaling_array[i]-1]+ysub;
                        master_z_array[local_map[conf2]+atom_scaling_array[i]-1]=master_z_array[local_map[conf2]+atom_scaling_array[i]-1]+zsub;
                        // outbound and rest
                        for(j=0;j<local_atoms;++j)
                        {
                            if(local_topo_array[local_map[j]-1]==-1 ||
                               local_topo_array[local_map[j]-1]==-2 ||
                               (local_topo_array[local_map[j]-1]==2 && j!=conf1))
                            {
                                master_x_array[local_map[j]+atom_scaling_array[i]-1]=master_x_array[local_map[j]+atom_scaling_array[i]-1]+xsub;
                                master_y_array[local_map[j]+atom_scaling_array[i]-1]=master_y_array[local_map[j]+atom_scaling_array[i]-1]+ysub;
                                master_z_array[local_map[j]+atom_scaling_array[i]-1]=master_z_array[local_map[j]+atom_scaling_array[i]-1]+zsub;
                            }
                        }
                        // update tacticity vector
                        xi=master_x_array[local_map[conf1]+atom_scaling_array[i]-1];
                        yi=master_y_array[local_map[conf1]+atom_scaling_array[i]-1];
                        zi=master_z_array[local_map[conf1]+atom_scaling_array[i]-1];
                        xj=master_x_array[local_map[conf2]+atom_scaling_array[i]-1];
                        yj=master_y_array[local_map[conf2]+atom_scaling_array[i]-1];
                        zj=master_z_array[local_map[conf2]+atom_scaling_array[i]-1];
                        xk=master_x_array[local_map[conf3]+atom_scaling_array[i]-1];
                        yk=master_y_array[local_map[conf3]+atom_scaling_array[i]-1];
                        zk=master_z_array[local_map[conf3]+atom_scaling_array[i]-1];
                        xl=master_x_array[local_map[conf4]+atom_scaling_array[i]-1];
                        yl=master_y_array[local_map[conf4]+atom_scaling_array[i]-1];
                        zl=master_z_array[local_map[conf4]+atom_scaling_array[i]-1];
                        coeff=(yl-yj)*(xj*zi-xk*zi-xi*zj+xk*zj+xi*zk-xj*zk)+(xl-xj)*(-yj*zi+yk*zi+yi*zj-yk*zj-yi*zk+yj*zk)+(-xj*yi+xk*yi+xi*yj-xk*yj-xi*yk+xj*yk)*(zl-zj);
                        tacticity_array[tacticity_counter_array[i]-1+yet_one_more_t_array[i]]=coeff;
                        //printf("@@@\t%d\n",tacticity_counter_array[i]-1+yet_one_more_t_array[i]);
                        
                        //--
                        
                        // console out
                        //printf("### mixed product: %lf\n",coeff);
                    }
                }
                else if (strcmp(tact_type,"i")==0)
                {
                    if(meso==0)
                    {
                        // bring first atom to axis origin: apply to all relevant atoms
                        // inbound coords
                        xsub=master_x_array[local_map[conf1]+atom_scaling_array[i]-1];
                        ysub=master_y_array[local_map[conf1]+atom_scaling_array[i]-1];
                        zsub=master_z_array[local_map[conf1]+atom_scaling_array[i]-1];
                        // reset inbound
                        master_x_array[local_map[conf1]+atom_scaling_array[i]-1]=master_x_array[local_map[conf1]+atom_scaling_array[i]-1]-xsub;
                        master_y_array[local_map[conf1]+atom_scaling_array[i]-1]=master_y_array[local_map[conf1]+atom_scaling_array[i]-1]-ysub;
                        master_z_array[local_map[conf1]+atom_scaling_array[i]-1]=master_z_array[local_map[conf1]+atom_scaling_array[i]-1]-zsub;
                        // reset linking
                        master_x_array[local_map[conf2]+atom_scaling_array[i]-1]=master_x_array[local_map[conf2]+atom_scaling_array[i]-1]-xsub;
                        master_y_array[local_map[conf2]+atom_scaling_array[i]-1]=master_y_array[local_map[conf2]+atom_scaling_array[i]-1]-ysub;
                        master_z_array[local_map[conf2]+atom_scaling_array[i]-1]=master_z_array[local_map[conf2]+atom_scaling_array[i]-1]-zsub;
                        // reset outbound and atoms to be reflected
                        for(j=0;j<local_atoms;++j)
                        {
                            if(local_topo_array[local_map[j]-1]==-1 ||
                               local_topo_array[local_map[j]-1]==-2 ||
                               (local_topo_array[local_map[j]-1]==2 && j!=conf1))
                            {
                                master_x_array[local_map[j]+atom_scaling_array[i]-1]=master_x_array[local_map[j]+atom_scaling_array[i]-1]-xsub;
                                master_y_array[local_map[j]+atom_scaling_array[i]-1]=master_y_array[local_map[j]+atom_scaling_array[i]-1]-ysub;
                                master_z_array[local_map[j]+atom_scaling_array[i]-1]=master_z_array[local_map[j]+atom_scaling_array[i]-1]-zsub;
                            }
                        }
                        // console out
                        /*
                        printf("!!! reflection plane atoms:\n");
                        printf("$ %d(%s)\t%lf\t%lf\t%lf\n",local_map[conf1]+atom_scaling_array[i],segment_atom_name_core[local_map[conf1]-1],
                               master_x_array[local_map[conf1]+atom_scaling_array[i]-1],
                               master_y_array[local_map[conf1]+atom_scaling_array[i]-1],
                               master_z_array[local_map[conf1]+atom_scaling_array[i]-1]);
                        printf("$ %d(%s)\t%lf\t%lf\t%lf\n",local_map[conf2]+atom_scaling_array[i],segment_atom_name_core[local_map[conf2]-1],
                               master_x_array[local_map[conf2]+atom_scaling_array[i]-1],
                               master_y_array[local_map[conf2]+atom_scaling_array[i]-1],
                               master_z_array[local_map[conf2]+atom_scaling_array[i]-1]);
                        printf("$ %d(%s)\t%lf\t%lf\t%lf\n",local_map[conf3]+atom_scaling_array[i],segment_atom_name_core[local_map[conf3]-1],
                               master_x_array[local_map[conf3]+atom_scaling_array[i]-1],
                               master_y_array[local_map[conf3]+atom_scaling_array[i]-1],
                               master_z_array[local_map[conf3]+atom_scaling_array[i]-1]);
                        */
                        //-- mixed product calculator --
                        
                        // normal plane vector
                        // coords
                        xi=master_x_array[local_map[conf1]+atom_scaling_array[i]-1];
                        yi=master_y_array[local_map[conf1]+atom_scaling_array[i]-1];
                        zi=master_z_array[local_map[conf1]+atom_scaling_array[i]-1];
                        xj=master_x_array[local_map[conf2]+atom_scaling_array[i]-1];
                        yj=master_y_array[local_map[conf2]+atom_scaling_array[i]-1];
                        zj=master_z_array[local_map[conf2]+atom_scaling_array[i]-1];
                        xk=master_x_array[local_map[conf3]+atom_scaling_array[i]-1];
                        yk=master_y_array[local_map[conf3]+atom_scaling_array[i]-1];
                        zk=master_z_array[local_map[conf3]+atom_scaling_array[i]-1];
                        // vectors:
                        //    x       y       z
                        // (xi-xj) (yi-yj) (zi-zj)
                        // (xk-xj) (yk-yj) (zk-zj)
                        // cross product
                        vx=(yi-yj)*(zk-zj)-(zi-zj)*(yk-yj);
                        vy=(zi-zj)*(xk-xj)-(xi-xj)*(zk-zj);
                        vz=(xi-xj)*(yk-yj)-(yi-yj)*(xk-xj);
                        norm=sqrt(vx*vx+vy*vy+vz*vz);
                        vx=vx/norm;vy=vy/norm;vz=vz/norm;
                        // console out
                        //printf("!!! normalized plane vector:\n");
                        //printf("%lf\t%lf\t%lf\n",vx,vy,vz);
                        // define reflection matrix
                        reflect[0][0]=1.0-2.0*vx*vx;
                        reflect[0][1]=   -2.0*vx*vy;
                        reflect[0][2]=   -2.0*vx*vz;
                        reflect[1][0]=   -2.0*vy*vx;
                        reflect[1][1]=1.0-2.0*vy*vy;
                        reflect[1][2]=   -2.0*vy*vz;
                        reflect[2][0]=   -2.0*vz*vx;
                        reflect[2][1]=   -2.0*vz*vy;
                        reflect[2][2]=1.0-2.0*vz*vz;
                        // apply the reflection matrix
                        for(j=0;j<local_atoms;++j)
                        {
                            if(local_topo_array[local_map[j]-1]==-1 ||
                               local_topo_array[local_map[j]-1]==-2 ||
                               (local_topo_array[local_map[j]-1]==2 && j!=conf1))
                            {
                                xb=master_x_array[local_map[j]+atom_scaling_array[i]-1];
                                yb=master_y_array[local_map[j]+atom_scaling_array[i]-1];
                                zb=master_z_array[local_map[j]+atom_scaling_array[i]-1];
                                master_x_array[local_map[j]+atom_scaling_array[i]-1]=reflect[0][0]*xb + reflect[0][1]*yb + reflect[0][2]*zb;
                                master_y_array[local_map[j]+atom_scaling_array[i]-1]=reflect[1][0]*xb + reflect[1][1]*yb + reflect[1][2]*zb;
                                master_z_array[local_map[j]+atom_scaling_array[i]-1]=reflect[2][0]*xb + reflect[2][1]*yb + reflect[2][2]*zb;
                            }
                        }
                        // restore coords
                        // inbound
                        master_x_array[local_map[conf1]+atom_scaling_array[i]-1]=master_x_array[local_map[conf1]+atom_scaling_array[i]-1]+xsub;
                        master_y_array[local_map[conf1]+atom_scaling_array[i]-1]=master_y_array[local_map[conf1]+atom_scaling_array[i]-1]+ysub;
                        master_z_array[local_map[conf1]+atom_scaling_array[i]-1]=master_z_array[local_map[conf1]+atom_scaling_array[i]-1]+zsub;
                        // linking
                        master_x_array[local_map[conf2]+atom_scaling_array[i]-1]=master_x_array[local_map[conf2]+atom_scaling_array[i]-1]+xsub;
                        master_y_array[local_map[conf2]+atom_scaling_array[i]-1]=master_y_array[local_map[conf2]+atom_scaling_array[i]-1]+ysub;
                        master_z_array[local_map[conf2]+atom_scaling_array[i]-1]=master_z_array[local_map[conf2]+atom_scaling_array[i]-1]+zsub;
                        // outbound and rest
                        for(j=0;j<local_atoms;++j)
                        {
                            if(local_topo_array[local_map[j]-1]==-1 ||
                               local_topo_array[local_map[j]-1]==-2 ||
                               (local_topo_array[local_map[j]-1]==2 && j!=conf1))
                            {
                                master_x_array[local_map[j]+atom_scaling_array[i]-1]=master_x_array[local_map[j]+atom_scaling_array[i]-1]+xsub;
                                master_y_array[local_map[j]+atom_scaling_array[i]-1]=master_y_array[local_map[j]+atom_scaling_array[i]-1]+ysub;
                                master_z_array[local_map[j]+atom_scaling_array[i]-1]=master_z_array[local_map[j]+atom_scaling_array[i]-1]+zsub;
                            }
                        }
                        // update tacticity vector
                        xi=master_x_array[local_map[conf1]+atom_scaling_array[i]-1];
                        yi=master_y_array[local_map[conf1]+atom_scaling_array[i]-1];
                        zi=master_z_array[local_map[conf1]+atom_scaling_array[i]-1];
                        xj=master_x_array[local_map[conf2]+atom_scaling_array[i]-1];
                        yj=master_y_array[local_map[conf2]+atom_scaling_array[i]-1];
                        zj=master_z_array[local_map[conf2]+atom_scaling_array[i]-1];
                        xk=master_x_array[local_map[conf3]+atom_scaling_array[i]-1];
                        yk=master_y_array[local_map[conf3]+atom_scaling_array[i]-1];
                        zk=master_z_array[local_map[conf3]+atom_scaling_array[i]-1];
                        xl=master_x_array[local_map[conf4]+atom_scaling_array[i]-1];
                        yl=master_y_array[local_map[conf4]+atom_scaling_array[i]-1];
                        zl=master_z_array[local_map[conf4]+atom_scaling_array[i]-1];
                        coeff=(yl-yj)*(xj*zi-xk*zi-xi*zj+xk*zj+xi*zk-xj*zk)+(xl-xj)*(-yj*zi+yk*zi+yi*zj-yk*zj-yi*zk+yj*zk)+(-xj*yi+xk*yi+xi*yj-xk*yj-xi*yk+xj*yk)*(zl-zj);
                        tacticity_array[tacticity_counter_array[i]-1+yet_one_more_t_array[i]]=coeff;
                        //printf("@@@\t%d\n",tacticity_counter_array[i]-1+yet_one_more_t_array[i]);
                        
                        //--
                        
                        // console out
                        //printf("### mixed product: %lf\n",coeff);
                    }
                }
                else if (strcmp(tact_type,"a")==0)
                {
                    bb_z=(double)rand()/RAND_MAX;
                    
                    //printf("%lf\n",bb_z);
                    //getchar();
                    
                    // force meso if <0.5
                    if(bb_z<0.5 && meso==0)
                    {
                        // bring first atom to axis origin: apply to all relevant atoms
                        // inbound coords
                        xsub=master_x_array[local_map[conf1]+atom_scaling_array[i]-1];
                        ysub=master_y_array[local_map[conf1]+atom_scaling_array[i]-1];
                        zsub=master_z_array[local_map[conf1]+atom_scaling_array[i]-1];
                        // reset inbound
                        master_x_array[local_map[conf1]+atom_scaling_array[i]-1]=master_x_array[local_map[conf1]+atom_scaling_array[i]-1]-xsub;
                        master_y_array[local_map[conf1]+atom_scaling_array[i]-1]=master_y_array[local_map[conf1]+atom_scaling_array[i]-1]-ysub;
                        master_z_array[local_map[conf1]+atom_scaling_array[i]-1]=master_z_array[local_map[conf1]+atom_scaling_array[i]-1]-zsub;
                        // reset linking
                        master_x_array[local_map[conf2]+atom_scaling_array[i]-1]=master_x_array[local_map[conf2]+atom_scaling_array[i]-1]-xsub;
                        master_y_array[local_map[conf2]+atom_scaling_array[i]-1]=master_y_array[local_map[conf2]+atom_scaling_array[i]-1]-ysub;
                        master_z_array[local_map[conf2]+atom_scaling_array[i]-1]=master_z_array[local_map[conf2]+atom_scaling_array[i]-1]-zsub;
                        // reset outbound and atoms to be reflected
                        for(j=0;j<local_atoms;++j)
                        {
                            if(local_topo_array[local_map[j]-1]==-1 ||
                               local_topo_array[local_map[j]-1]==-2 ||
                               (local_topo_array[local_map[j]-1]==2 && j!=conf1))
                            {
                                master_x_array[local_map[j]+atom_scaling_array[i]-1]=master_x_array[local_map[j]+atom_scaling_array[i]-1]-xsub;
                                master_y_array[local_map[j]+atom_scaling_array[i]-1]=master_y_array[local_map[j]+atom_scaling_array[i]-1]-ysub;
                                master_z_array[local_map[j]+atom_scaling_array[i]-1]=master_z_array[local_map[j]+atom_scaling_array[i]-1]-zsub;
                            }
                        }
                        // console out
                        /*
                        printf("!!! reflection plane atoms:\n");
                        printf("$ %d(%s)\t%lf\t%lf\t%lf\n",local_map[conf1]+atom_scaling_array[i],segment_atom_name_core[local_map[conf1]-1],
                               master_x_array[local_map[conf1]+atom_scaling_array[i]-1],
                               master_y_array[local_map[conf1]+atom_scaling_array[i]-1],
                               master_z_array[local_map[conf1]+atom_scaling_array[i]-1]);
                        printf("$ %d(%s)\t%lf\t%lf\t%lf\n",local_map[conf2]+atom_scaling_array[i],segment_atom_name_core[local_map[conf2]-1],
                               master_x_array[local_map[conf2]+atom_scaling_array[i]-1],
                               master_y_array[local_map[conf2]+atom_scaling_array[i]-1],
                               master_z_array[local_map[conf2]+atom_scaling_array[i]-1]);
                        printf("$ %d(%s)\t%lf\t%lf\t%lf\n",local_map[conf3]+atom_scaling_array[i],segment_atom_name_core[local_map[conf3]-1],
                               master_x_array[local_map[conf3]+atom_scaling_array[i]-1],
                               master_y_array[local_map[conf3]+atom_scaling_array[i]-1],
                               master_z_array[local_map[conf3]+atom_scaling_array[i]-1]);
                        */
                        //-- mixed product calculator --
                        
                        // normal plane vector
                        // coords
                        xi=master_x_array[local_map[conf1]+atom_scaling_array[i]-1];
                        yi=master_y_array[local_map[conf1]+atom_scaling_array[i]-1];
                        zi=master_z_array[local_map[conf1]+atom_scaling_array[i]-1];
                        xj=master_x_array[local_map[conf2]+atom_scaling_array[i]-1];
                        yj=master_y_array[local_map[conf2]+atom_scaling_array[i]-1];
                        zj=master_z_array[local_map[conf2]+atom_scaling_array[i]-1];
                        xk=master_x_array[local_map[conf3]+atom_scaling_array[i]-1];
                        yk=master_y_array[local_map[conf3]+atom_scaling_array[i]-1];
                        zk=master_z_array[local_map[conf3]+atom_scaling_array[i]-1];
                        // vectors:
                        //    x       y       z
                        // (xi-xj) (yi-yj) (zi-zj)
                        // (xk-xj) (yk-yj) (zk-zj)
                        // cross product
                        vx=(yi-yj)*(zk-zj)-(zi-zj)*(yk-yj);
                        vy=(zi-zj)*(xk-xj)-(xi-xj)*(zk-zj);
                        vz=(xi-xj)*(yk-yj)-(yi-yj)*(xk-xj);
                        norm=sqrt(vx*vx+vy*vy+vz*vz);
                        vx=vx/norm;vy=vy/norm;vz=vz/norm;
                        // console out
                        //printf("!!! normalized plane vector:\n");
                        //printf("%lf\t%lf\t%lf\n",vx,vy,vz);
                        // define reflection matrix
                        reflect[0][0]=1.0-2.0*vx*vx;
                        reflect[0][1]=   -2.0*vx*vy;
                        reflect[0][2]=   -2.0*vx*vz;
                        reflect[1][0]=   -2.0*vy*vx;
                        reflect[1][1]=1.0-2.0*vy*vy;
                        reflect[1][2]=   -2.0*vy*vz;
                        reflect[2][0]=   -2.0*vz*vx;
                        reflect[2][1]=   -2.0*vz*vy;
                        reflect[2][2]=1.0-2.0*vz*vz;
                        // apply the reflection matrix
                        for(j=0;j<local_atoms;++j)
                        {
                            if(local_topo_array[local_map[j]-1]==-1 ||
                               local_topo_array[local_map[j]-1]==-2 ||
                               (local_topo_array[local_map[j]-1]==2 && j!=conf1))
                            {
                                xb=master_x_array[local_map[j]+atom_scaling_array[i]-1];
                                yb=master_y_array[local_map[j]+atom_scaling_array[i]-1];
                                zb=master_z_array[local_map[j]+atom_scaling_array[i]-1];
                                master_x_array[local_map[j]+atom_scaling_array[i]-1]=reflect[0][0]*xb + reflect[0][1]*yb + reflect[0][2]*zb;
                                master_y_array[local_map[j]+atom_scaling_array[i]-1]=reflect[1][0]*xb + reflect[1][1]*yb + reflect[1][2]*zb;
                                master_z_array[local_map[j]+atom_scaling_array[i]-1]=reflect[2][0]*xb + reflect[2][1]*yb + reflect[2][2]*zb;
                            }
                        }
                        // restore coords
                        // inbound
                        master_x_array[local_map[conf1]+atom_scaling_array[i]-1]=master_x_array[local_map[conf1]+atom_scaling_array[i]-1]+xsub;
                        master_y_array[local_map[conf1]+atom_scaling_array[i]-1]=master_y_array[local_map[conf1]+atom_scaling_array[i]-1]+ysub;
                        master_z_array[local_map[conf1]+atom_scaling_array[i]-1]=master_z_array[local_map[conf1]+atom_scaling_array[i]-1]+zsub;
                        // linking
                        master_x_array[local_map[conf2]+atom_scaling_array[i]-1]=master_x_array[local_map[conf2]+atom_scaling_array[i]-1]+xsub;
                        master_y_array[local_map[conf2]+atom_scaling_array[i]-1]=master_y_array[local_map[conf2]+atom_scaling_array[i]-1]+ysub;
                        master_z_array[local_map[conf2]+atom_scaling_array[i]-1]=master_z_array[local_map[conf2]+atom_scaling_array[i]-1]+zsub;
                        // outbound and rest
                        for(j=0;j<local_atoms;++j)
                        {
                            if(local_topo_array[local_map[j]-1]==-1 ||
                               local_topo_array[local_map[j]-1]==-2 ||
                               (local_topo_array[local_map[j]-1]==2 && j!=conf1))
                            {
                                master_x_array[local_map[j]+atom_scaling_array[i]-1]=master_x_array[local_map[j]+atom_scaling_array[i]-1]+xsub;
                                master_y_array[local_map[j]+atom_scaling_array[i]-1]=master_y_array[local_map[j]+atom_scaling_array[i]-1]+ysub;
                                master_z_array[local_map[j]+atom_scaling_array[i]-1]=master_z_array[local_map[j]+atom_scaling_array[i]-1]+zsub;
                            }
                        }
                        // update tacticity vector
                        xi=master_x_array[local_map[conf1]+atom_scaling_array[i]-1];
                        yi=master_y_array[local_map[conf1]+atom_scaling_array[i]-1];
                        zi=master_z_array[local_map[conf1]+atom_scaling_array[i]-1];
                        xj=master_x_array[local_map[conf2]+atom_scaling_array[i]-1];
                        yj=master_y_array[local_map[conf2]+atom_scaling_array[i]-1];
                        zj=master_z_array[local_map[conf2]+atom_scaling_array[i]-1];
                        xk=master_x_array[local_map[conf3]+atom_scaling_array[i]-1];
                        yk=master_y_array[local_map[conf3]+atom_scaling_array[i]-1];
                        zk=master_z_array[local_map[conf3]+atom_scaling_array[i]-1];
                        xl=master_x_array[local_map[conf4]+atom_scaling_array[i]-1];
                        yl=master_y_array[local_map[conf4]+atom_scaling_array[i]-1];
                        zl=master_z_array[local_map[conf4]+atom_scaling_array[i]-1];
                        coeff=(yl-yj)*(xj*zi-xk*zi-xi*zj+xk*zj+xi*zk-xj*zk)+(xl-xj)*(-yj*zi+yk*zi+yi*zj-yk*zj-yi*zk+yj*zk)+(-xj*yi+xk*yi+xi*yj-xk*yj-xi*yk+xj*yk)*(zl-zj);
                        tacticity_array[tacticity_counter_array[i]-1+yet_one_more_t_array[i]]=coeff;
                        //printf("@@@\t%d\n",tacticity_counter_array[i]-1+yet_one_more_t_array[i]);
                        
                        meso=1;
                        
                        //--
                        
                        // console out
                        //printf("### mixed product: %lf\n",coeff);
                        
                    }
                    
                    // force racemo if >=0.5
                    if(bb_z>0.5 && meso==1)
                    {
                        // bring first atom to axis origin: apply to all relevant atoms
                        // inbound coords
                        xsub=master_x_array[local_map[conf1]+atom_scaling_array[i]-1];
                        ysub=master_y_array[local_map[conf1]+atom_scaling_array[i]-1];
                        zsub=master_z_array[local_map[conf1]+atom_scaling_array[i]-1];
                        // reset inbound
                        master_x_array[local_map[conf1]+atom_scaling_array[i]-1]=master_x_array[local_map[conf1]+atom_scaling_array[i]-1]-xsub;
                        master_y_array[local_map[conf1]+atom_scaling_array[i]-1]=master_y_array[local_map[conf1]+atom_scaling_array[i]-1]-ysub;
                        master_z_array[local_map[conf1]+atom_scaling_array[i]-1]=master_z_array[local_map[conf1]+atom_scaling_array[i]-1]-zsub;
                        // reset linking
                        master_x_array[local_map[conf2]+atom_scaling_array[i]-1]=master_x_array[local_map[conf2]+atom_scaling_array[i]-1]-xsub;
                        master_y_array[local_map[conf2]+atom_scaling_array[i]-1]=master_y_array[local_map[conf2]+atom_scaling_array[i]-1]-ysub;
                        master_z_array[local_map[conf2]+atom_scaling_array[i]-1]=master_z_array[local_map[conf2]+atom_scaling_array[i]-1]-zsub;
                        // reset outbound and atoms to be reflected
                        for(j=0;j<local_atoms;++j)
                        {
                            if(local_topo_array[local_map[j]-1]==-1 ||
                               local_topo_array[local_map[j]-1]==-2 ||
                               (local_topo_array[local_map[j]-1]==2 && j!=conf1))
                            {
                                master_x_array[local_map[j]+atom_scaling_array[i]-1]=master_x_array[local_map[j]+atom_scaling_array[i]-1]-xsub;
                                master_y_array[local_map[j]+atom_scaling_array[i]-1]=master_y_array[local_map[j]+atom_scaling_array[i]-1]-ysub;
                                master_z_array[local_map[j]+atom_scaling_array[i]-1]=master_z_array[local_map[j]+atom_scaling_array[i]-1]-zsub;
                            }
                        }
                        // console out
                        /*
                        printf("!!! reflection plane atoms:\n");
                        printf("$ %d(%s)\t%lf\t%lf\t%lf\n",local_map[conf1]+atom_scaling_array[i],segment_atom_name_core[local_map[conf1]-1],
                               master_x_array[local_map[conf1]+atom_scaling_array[i]-1],
                               master_y_array[local_map[conf1]+atom_scaling_array[i]-1],
                               master_z_array[local_map[conf1]+atom_scaling_array[i]-1]);
                        printf("$ %d(%s)\t%lf\t%lf\t%lf\n",local_map[conf2]+atom_scaling_array[i],segment_atom_name_core[local_map[conf2]-1],
                               master_x_array[local_map[conf2]+atom_scaling_array[i]-1],
                               master_y_array[local_map[conf2]+atom_scaling_array[i]-1],
                               master_z_array[local_map[conf2]+atom_scaling_array[i]-1]);
                        printf("$ %d(%s)\t%lf\t%lf\t%lf\n",local_map[conf3]+atom_scaling_array[i],segment_atom_name_core[local_map[conf3]-1],
                               master_x_array[local_map[conf3]+atom_scaling_array[i]-1],
                               master_y_array[local_map[conf3]+atom_scaling_array[i]-1],
                               master_z_array[local_map[conf3]+atom_scaling_array[i]-1]);
                        */
                        //-- mixed product calculator --
                        
                        // normal plane vector
                        // coords
                        xi=master_x_array[local_map[conf1]+atom_scaling_array[i]-1];
                        yi=master_y_array[local_map[conf1]+atom_scaling_array[i]-1];
                        zi=master_z_array[local_map[conf1]+atom_scaling_array[i]-1];
                        xj=master_x_array[local_map[conf2]+atom_scaling_array[i]-1];
                        yj=master_y_array[local_map[conf2]+atom_scaling_array[i]-1];
                        zj=master_z_array[local_map[conf2]+atom_scaling_array[i]-1];
                        xk=master_x_array[local_map[conf3]+atom_scaling_array[i]-1];
                        yk=master_y_array[local_map[conf3]+atom_scaling_array[i]-1];
                        zk=master_z_array[local_map[conf3]+atom_scaling_array[i]-1];
                        // vectors:
                        //    x       y       z
                        // (xi-xj) (yi-yj) (zi-zj)
                        // (xk-xj) (yk-yj) (zk-zj)
                        // cross product
                        vx=(yi-yj)*(zk-zj)-(zi-zj)*(yk-yj);
                        vy=(zi-zj)*(xk-xj)-(xi-xj)*(zk-zj);
                        vz=(xi-xj)*(yk-yj)-(yi-yj)*(xk-xj);
                        norm=sqrt(vx*vx+vy*vy+vz*vz);
                        vx=vx/norm;vy=vy/norm;vz=vz/norm;
                        // console out
                        //printf("!!! normalized plane vector:\n");
                        //printf("%lf\t%lf\t%lf\n",vx,vy,vz);
                        // define reflection matrix
                        reflect[0][0]=1.0-2.0*vx*vx;
                        reflect[0][1]=   -2.0*vx*vy;
                        reflect[0][2]=   -2.0*vx*vz;
                        reflect[1][0]=   -2.0*vy*vx;
                        reflect[1][1]=1.0-2.0*vy*vy;
                        reflect[1][2]=   -2.0*vy*vz;
                        reflect[2][0]=   -2.0*vz*vx;
                        reflect[2][1]=   -2.0*vz*vy;
                        reflect[2][2]=1.0-2.0*vz*vz;
                        // apply the reflection matrix
                        for(j=0;j<local_atoms;++j)
                        {
                            if(local_topo_array[local_map[j]-1]==-1 ||
                               local_topo_array[local_map[j]-1]==-2 ||
                               (local_topo_array[local_map[j]-1]==2 && j!=conf1))
                            {
                                xb=master_x_array[local_map[j]+atom_scaling_array[i]-1];
                                yb=master_y_array[local_map[j]+atom_scaling_array[i]-1];
                                zb=master_z_array[local_map[j]+atom_scaling_array[i]-1];
                                master_x_array[local_map[j]+atom_scaling_array[i]-1]=reflect[0][0]*xb + reflect[0][1]*yb + reflect[0][2]*zb;
                                master_y_array[local_map[j]+atom_scaling_array[i]-1]=reflect[1][0]*xb + reflect[1][1]*yb + reflect[1][2]*zb;
                                master_z_array[local_map[j]+atom_scaling_array[i]-1]=reflect[2][0]*xb + reflect[2][1]*yb + reflect[2][2]*zb;
                            }
                        }
                        // restore coords
                        // inbound
                        master_x_array[local_map[conf1]+atom_scaling_array[i]-1]=master_x_array[local_map[conf1]+atom_scaling_array[i]-1]+xsub;
                        master_y_array[local_map[conf1]+atom_scaling_array[i]-1]=master_y_array[local_map[conf1]+atom_scaling_array[i]-1]+ysub;
                        master_z_array[local_map[conf1]+atom_scaling_array[i]-1]=master_z_array[local_map[conf1]+atom_scaling_array[i]-1]+zsub;
                        // linking
                        master_x_array[local_map[conf2]+atom_scaling_array[i]-1]=master_x_array[local_map[conf2]+atom_scaling_array[i]-1]+xsub;
                        master_y_array[local_map[conf2]+atom_scaling_array[i]-1]=master_y_array[local_map[conf2]+atom_scaling_array[i]-1]+ysub;
                        master_z_array[local_map[conf2]+atom_scaling_array[i]-1]=master_z_array[local_map[conf2]+atom_scaling_array[i]-1]+zsub;
                        // outbound and rest
                        for(j=0;j<local_atoms;++j)
                        {
                            if(local_topo_array[local_map[j]-1]==-1 ||
                               local_topo_array[local_map[j]-1]==-2 ||
                               (local_topo_array[local_map[j]-1]==2 && j!=conf1))
                            {
                                master_x_array[local_map[j]+atom_scaling_array[i]-1]=master_x_array[local_map[j]+atom_scaling_array[i]-1]+xsub;
                                master_y_array[local_map[j]+atom_scaling_array[i]-1]=master_y_array[local_map[j]+atom_scaling_array[i]-1]+ysub;
                                master_z_array[local_map[j]+atom_scaling_array[i]-1]=master_z_array[local_map[j]+atom_scaling_array[i]-1]+zsub;
                            }
                        }
                        // update tacticity vector
                        xi=master_x_array[local_map[conf1]+atom_scaling_array[i]-1];
                        yi=master_y_array[local_map[conf1]+atom_scaling_array[i]-1];
                        zi=master_z_array[local_map[conf1]+atom_scaling_array[i]-1];
                        xj=master_x_array[local_map[conf2]+atom_scaling_array[i]-1];
                        yj=master_y_array[local_map[conf2]+atom_scaling_array[i]-1];
                        zj=master_z_array[local_map[conf2]+atom_scaling_array[i]-1];
                        xk=master_x_array[local_map[conf3]+atom_scaling_array[i]-1];
                        yk=master_y_array[local_map[conf3]+atom_scaling_array[i]-1];
                        zk=master_z_array[local_map[conf3]+atom_scaling_array[i]-1];
                        xl=master_x_array[local_map[conf4]+atom_scaling_array[i]-1];
                        yl=master_y_array[local_map[conf4]+atom_scaling_array[i]-1];
                        zl=master_z_array[local_map[conf4]+atom_scaling_array[i]-1];
                        coeff=(yl-yj)*(xj*zi-xk*zi-xi*zj+xk*zj+xi*zk-xj*zk)+(xl-xj)*(-yj*zi+yk*zi+yi*zj-yk*zj-yi*zk+yj*zk)+(-xj*yi+xk*yi+xi*yj-xk*yj-xi*yk+xj*yk)*(zl-zj);
                        tacticity_array[tacticity_counter_array[i]-1+yet_one_more_t_array[i]]=coeff;
                        //printf("@@@\t%d\n",tacticity_counter_array[i]-1+yet_one_more_t_array[i]);
                        
                        meso=0;
                        
                        //--
                        
                        // console out
                        //printf("### mixed product: %lf\n",coeff);
                        
                        
                        
                        
                    }
                    //printf("$$$ [%d] bb_z=%lf\tmeso=%d\n",*general_t_counter,bb_z,meso);
                    //getchar();
                    
                }
                
            }
            /*
            //
            printf("------------------------------------- - - - - --- - - -- \n$ after:\n");
            printf("tacticity_counter_array:\t");
            for(jj=0;jj<molecules;++jj)printf("%d\t",tacticity_counter_array[jj]);printf("\n");
            printf("general_t_counter:\t");
            printf("%d\n",*general_t_counter);
            if(*general_t_counter==0){printf("t_tracker_1\tN\\A\n");printf("t_tracker_2\tN\\A\n");printf("t_tracker_3\tN\\A\n");printf("t_tracker_4\tN\\A\n");printf("t_tracker_mol\tN\\A\n");}
            else
            {
                printf("t_tracker_1:\t");for(jj=0;jj<*general_t_counter;++jj)printf("%d\t",t_tracker_1[jj]);printf("\n");
                printf("t_tracker_2:\t");for(jj=0;jj<*general_t_counter;++jj)printf("%d\t",t_tracker_2[jj]);printf("\n");
                printf("t_tracker_3:\t");for(jj=0;jj<*general_t_counter;++jj)printf("%d\t",t_tracker_3[jj]);printf("\n");
                printf("t_tracker_4:\t");for(jj=0;jj<*general_t_counter;++jj)printf("%d\t",t_tracker_4[jj]);printf("\n");
                printf("t_tracker_mol:\t");for(jj=0;jj<*general_t_counter;++jj)printf("%d\t",t_tracker_mol[jj]);printf("\n");
            }
            printf("tacticity_array:\t");for(jj=0;jj<tacticity_array_length;++jj)printf("%lf\t",tacticity_array[jj]);printf("\n");
            //printf("yet_one_more_t_array:\t");for(jj=0;jj<molecules;++jj)printf("%d\t",yet_one_more_t_array[jj]);printf("\n");
            getchar();
            //
            */
        }
        //getchar();
    }

}
