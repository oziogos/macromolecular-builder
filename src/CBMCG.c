
#include"builder.h"

#define debug 0

void topo_neighbors(int general_atoms, int general_bonds, int *master_B1_core_array, int *master_B2_core_array,
                    
                    int ***one_two, int ***one_three, int ***one_four);

double lj_pe_brute_force(int general_atoms, double lx, double ly, double lz, double cutoff, int **one_two, int **one_three, int **one_four,double **nb_FF,int *nb_type_array,
                         double *master_x_array, double *master_y_array,double *master_z_array,int *master_nx_array, int *master_ny_array,int *master_nz_array,
                         int *pairs_bf,int *host_id_array);

double lj_pe_linked_list_cell(int general_atoms, double lx, double ly, double lz, double cutoff, int **one_two, int **one_three, int **one_four,double **nb_FF,int *nb_type_array,
                              double *master_x_array, double *master_y_array,double *master_z_array,int *master_nx_array, int *master_ny_array,int *master_nz_array,
                              
                              int M3, int *head, int *list, int **neighbors,
                              
                              int *pairs_llc,int *host_id_array
                              );

void calc_head_list(int general_atoms, int Mx, int M2, int M3,
                    double xlo, double ylo, double zlo, double xhi, double yhi, double zhi, double lx, double ly, double lz, double lcx, double lcy, double lcz,
                    double *master_x_array,double *master_y_array, double *master_z_array, int *master_nx_array, int *master_ny_array, int *master_nz_array,
                    int *head, int **list);

void LJ_params(int general_atoms, int general_atom_types,char **master_species_array, char **general_species_registry,
               int entries_UFF, char **species_UFF_array, double *D_UFF_array, double *x_UFF_array,
               int entries_DRE, char **species_DRE_array, double *D0_DRE_array,double *Rvdw0_DRE_array,
               double ***nb_FF, int **nb_type_array, int FF);

void calc_wrapped_coords_n_arrays(int general_atoms,double xlo,double ylo,double zlo,double xhi,double yhi,double zhi,double lx,double ly,double lz,
                                  double *master_x_array,double *master_y_array,double *master_z_array,int *master_nx_array,int *master_ny_array,int *master_nz_array);

void rotate(double ux, double uy, double uz, double x_center, double y_center, double z_center, double theta, double x_old, double y_old, double z_old, double *x_new, double *y_new, double *z_new);

void force_tacticity(int i, int tacticity, int *mol_type_array, char **tacticity_search, int local_atoms, char **segment_atom_name_core, int *local_map, int *local_topo_array,
                     int *atom_scaling_array, char **backbone_id, int *yet_one_more_t_array, int backbones,
                     int backbone_type_counter, char **backbone_type_array,
                     
                     int molecules,int tacticity_array_length,
                     
                     int *tacticity_counter_array,
                     int *general_t_counter, int *t_tracker_1, int *t_tracker_2, int *t_tracker_3, int *t_tracker_4,int *t_tracker_mol,
                     double *tacticity_array, double *master_x_array, double *master_y_array, double *master_z_array,
                     int *ft);

double lj_pe(int cellpart,int new_atoms_counter,double lx,double ly,double lz,double cutoff,
             int *cellpartarray,int *new_atoms_array,int *host_id_array,int **one_two,int **one_three,int **one_four,int *nb_type_array,double **nb_FF,
             double *master_x_array,double *master_y_array,double *master_z_array,int *master_nx_array,int *master_ny_array,int *master_nz_array,
             int *pairs_llc);

void CBMCG(int i,int local_atoms,int backbone_type_counter,int backbones,int CBMCG_N_trial,char *current_folder,
           int *local_topo_array,int *local_map,char **segment_atom_name_core,char **backbone_type_array,int *atom_scaling_array,char **backbone_id,
           double *bb_dx,double **bb_entries,double *bb_pdf_dx,double **bb_pdf_entries,
           
           int tacticity,int *mol_type_array,char **tacticity_search,int *yet_one_more_t_array,
           int molecules,int tacticity_array_length,
           int *tacticity_counter_array,
           int *general_t_counter, int *t_tracker_1, int *t_tracker_2, int *t_tracker_3, int *t_tracker_4,int *t_tracker_mol,
           double *tacticity_array,
           
           char **master_species_array,double *master_x_array,double *master_y_array,double *master_z_array,
           
           double cutoff,double xlo,double ylo,double zlo,double xhi,double yhi,double zhi,int general_atoms,int general_bonds,int general_atom_types,
           int Mx,int My,int Mz,int M3,double lx,double ly,double lz,double lcx,double lcy,double lcz,int **neighbors,
           int *master_B1_core_array,int *master_B2_core_array,char **general_species_registry,
           
           
           int entries_UFF,char **species_UFF_array,double *D_UFF_array,double *x_UFF_array,
           int entries_DRE,char **species_DRE_array,double *D0_DRE_array,double *Rvdw0_DRE_array,
           
           int *host_id_array,
           double *pe_array,int counter,
           
           double T,int check_with_cells,
           
           int FF
           
           )
{
    FILE *fp;
    char file_path[cmax_length];
    int check_t,found,j,k,conf1,conf2,conf3,conf4,f1,f2,f3,identified_D,bb_i,is_bb,o,ft;
    char D1_current[sub_length],D2_current[sub_length],D3_current[sub_length],D4_current[sub_length];
    char D1[sub_length],D2[sub_length],D3[sub_length],D4[sub_length];
    double local_x_backup[local_backup_length],local_y_backup[local_backup_length],local_z_backup[local_backup_length];
    double xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl,cross,norm1,norm2,sign,coeff,phi,vx,vy,vz;
    double bb_z,bb_zmap,bb_pdf;
    int bb_index,bb_pdf_index;
    double BB_d_angle_value[CBMCG_N_trial],BB_d_angle_pdf_value[CBMCG_N_trial];
    //
    
    int **one_two,**one_three,**one_four;
    double **nb_FF_for_lj;
    int *nb_type_array;
    int *master_nx_array,*master_ny_array,*master_nz_array;
    int pairs_bf,pairs_llc;
    int *list;
    double pe_bf,pe_llc;
    int head[M3];
    double local_pe_array[CBMCG_N_trial],boltzmann_w[CBMCG_N_trial],deltaU,N,sum,final[CBMCG_N_trial],ksi;
    int selected_k;
    
    //
    
    double one_master_x_array[1],one_master_y_array[1],one_master_z_array[1];
    int one_master_nx_array[1],one_master_ny_array[1],one_master_nz_array[1];
    int cindex,M2,h;
    int local_cells[27];
    int cells_registry[125];
    
    //
    
    int unique_cells,unique_cells_array[M3],cell_position_index,new_atoms_array[general_atoms],new_atoms_counter,l,unique_found,cells_registry_index,eligible[26],m;
    int cellpartarray[general_atoms],cellpart,headpointer,cellindex,listpointer;
    double local_delta_pe_array[CBMCG_N_trial];
    int pairs;
    
    double deltaU_test;
    
    //
    check_t=1;
    
    // check if we have a backbone D type
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
            
            //printf("%s\t%s\t%s\t%s\n",D1_current,D2_current,D3_current,D4_current);
            
            identified_D=0;
            for(k=0;k<backbones;++k)
            {
                sscanf(backbone_id[k],"%s\t%s\t%s\t%s",D1,D2,D3,D4);
                if((strcmp(D1,D1_current)==0 && strcmp(D2,D2_current)==0 && strcmp(D3,D3_current)==0 && strcmp(D4,D4_current)==0)
                   ||
                   (strcmp(D1,D4_current)==0 && strcmp(D2,D3_current)==0 && strcmp(D3,D2_current)==0 && strcmp(D4,D1_current)==0))
                {
                    identified_D=1;
                    bb_i=k;
                    break;
                    
                }
            }
            if(identified_D==1){
                // for backup
                /*
                printf("# local atoms: %d\n",local_atoms);
                printf("# candidate rotatable atoms:\n");
                for(j=0;j<local_atoms;++j)printf("[%d]\t%d\t%lf\t%lf\t%lf\n",j+1,local_map[j]+atom_scaling_array[i],master_x_array[local_map[j]+atom_scaling_array[i]-1],master_y_array[local_map[j]+atom_scaling_array[i]-1],master_z_array[local_map[j]+atom_scaling_array[i]-1]);
                */
                // backup coords
                if(local_atoms>local_backup_length){printf("$ local_atoms>local_backup_length\n\n");exit(-1);}
                for(j=0;j<local_atoms;++j)
                {
                    local_x_backup[j]=master_x_array[local_map[j]+atom_scaling_array[i]-1];
                    local_y_backup[j]=master_y_array[local_map[j]+atom_scaling_array[i]-1];
                    local_z_backup[j]=master_z_array[local_map[j]+atom_scaling_array[i]-1];
                }
                
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
                
                
                
                //sprintf(file_path,"%s/%s",current_folder,"rot_prev.xyz");
                //fp=fopen(file_path,"w+");
                /*
                 fprintf(fp,"%d\n\n",local_atoms);
                 for(j=0;j<local_atoms;++j)fprintf(fp,"%c\t%lf\t%lf\t%lf\n",master_species_array[local_map[j]+atom_scaling_array[i]-1][0],local_x_backup[j],local_y_backup[j],local_z_backup[j]);
                 */
                //
                //printf("# will generate %d trial conformations!\n",CBMCG_N_trial);
                //
                //-------------------------------------------------------------------------- not affected by coords ---
                topo_neighbors(general_atoms, general_bonds, master_B1_core_array, master_B2_core_array,&one_two, &one_three, &one_four);
                //-------------------------------------------------------------------------- not affected by coords ---
                LJ_params(general_atoms,general_atom_types,master_species_array,general_species_registry,
                          entries_UFF,species_UFF_array,D_UFF_array,x_UFF_array,
                          entries_DRE,species_DRE_array,D0_DRE_array,Rvdw0_DRE_array,
                          &nb_FF_for_lj,&nb_type_array,FF);
                //--------------------------------------------------------------------------
                master_nx_array=(int*)malloc(general_atoms*sizeof(int));
                master_ny_array=(int*)malloc(general_atoms*sizeof(int));
                master_nz_array=(int*)malloc(general_atoms*sizeof(int));
                //
                N=0.0;
                for(k=0;k<CBMCG_N_trial;++k)
                {
                    // draw random number between 0 and 1: bb_z
                    bb_z=(double)rand()/RAND_MAX;
                    // select dihedral bb_zmap from the type bb_i invCDF
                    bb_index=(int)floor(bb_z/bb_dx[bb_i]);
                    bb_zmap=bb_entries[bb_index][1+bb_i*2]+(bb_entries[bb_index+1][1+bb_i*2]-bb_entries[bb_index][1+bb_i*2])*(bb_z-bb_entries[bb_index][0+bb_i*2])/bb_dx[bb_i];
                    // resolve PDF
                    bb_pdf_index=(int)floor((bb_zmap+180.0)/bb_pdf_dx[bb_i]);
                    bb_pdf=bb_pdf_entries[bb_pdf_index][1+bb_i*2]+(bb_pdf_entries[bb_pdf_index+1][1+bb_i*2]-bb_pdf_entries[bb_pdf_index][1+bb_i*2])*(bb_zmap-bb_pdf_entries[bb_pdf_index][0+bb_i*2])/bb_pdf_dx[bb_i];
                    
                    
                    //printf("[%d]\t%lf\t%lf\n",k+1,bb_zmap,bb_pdf);
                    BB_d_angle_value[k]=bb_zmap;
                    BB_d_angle_pdf_value[k]=bb_pdf;
                    
                    
                    //
                    
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
                                   local_x_backup[j],
                                   local_y_backup[j],
                                   local_z_backup[j],
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
                                   local_x_backup[j],
                                   local_y_backup[j],
                                   local_z_backup[j],
                                   &master_x_array[local_map[j]+atom_scaling_array[i]-1],
                                   &master_y_array[local_map[j]+atom_scaling_array[i]-1],
                                   &master_z_array[local_map[j]+atom_scaling_array[i]-1]);
                        }
                        // rotate the rest 1-2 atoms
                        //if(local_topo_array[local_map[j]-1]==2 && strcmp(segment_atom_name_core[local_map[j]-1],backbone_id[k])!=0)
                        if(local_topo_array[local_map[j]-1]==2)
                        {
                            is_bb=0;
                            for(o=0;o<backbone_type_counter;++o)
                            {
                                if(strcmp(segment_atom_name_core[local_map[j]-1],backbone_type_array[o])==0)
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
                                       local_x_backup[j],
                                       local_y_backup[j],
                                       local_z_backup[j],
                                       &master_x_array[local_map[j]+atom_scaling_array[i]-1],
                                       &master_y_array[local_map[j]+atom_scaling_array[i]-1],
                                       &master_z_array[local_map[j]+atom_scaling_array[i]-1]);
                            }
                        }
                    }

                    //
                    
                    check_t=0;
                    
                    //
                    
                    if(tacticity!=0)
                    force_tacticity(i,tacticity,mol_type_array,tacticity_search,local_atoms,segment_atom_name_core,local_map,local_topo_array,
                                    atom_scaling_array,backbone_id,yet_one_more_t_array,backbones,backbone_type_counter,backbone_type_array,
                                    
                                    molecules,tacticity_array_length,
                                    
                                    tacticity_counter_array,
                                    //&general_t_counter,t_tracker_1,t_tracker_2,t_tracker_3,t_tracker_4,t_tracker_mol,
                                    general_t_counter,t_tracker_1,t_tracker_2,t_tracker_3,t_tracker_4,t_tracker_mol,
                                    tacticity_array,master_x_array,master_y_array,master_z_array,
                                    &ft);
                    
                    //
                    
                    
                    // CBMCG scales terribly with system size since it evaluates CBMCG_N_trial times the potential energy calculation function that currently
                    // loops on all cells!!
                    //
                    // a candidate solution:
                    // - loop on local atoms and locate -1s only: these are the newly added atoms. This approach should be valid since we carry out implicit hydrogen
                    //   energy calculations, so ignoring hydrogens on the host atom and the newly added atoms will pose no problem.
                    // - find the cell index of each -1 atom
                    
                    unique_cells=0;new_atoms_counter=0;
                    for(j=0;j<local_atoms;++j)
                    {
                        //printf("[%d]\t%d\t%d\n",j+1,local_map[j]+atom_scaling_array[i],local_topo_array[local_map[j]-1]);
                        
                        if(local_topo_array[local_map[j]-1]==-1)
                        {
                            new_atoms_counter=new_atoms_counter+1;
                            M2=Mx*My;
                            //h=segment_id_internal_array[host_index]+atom_scaling_array[i]-1;
                            h=local_map[j]-1+atom_scaling_array[i];
                            //printf("$ noH new atom: %d\n",h+1);
                            one_master_x_array[0]=master_x_array[h];
                            one_master_y_array[0]=master_y_array[h];
                            one_master_z_array[0]=master_z_array[h];
                            calc_wrapped_coords_n_arrays(1,xlo,ylo,zlo,xhi,yhi,zhi,lx,ly,lz,
                                                         one_master_x_array,one_master_y_array,one_master_z_array,one_master_nx_array,one_master_ny_array,one_master_nz_array);
                            cindex=floor((one_master_x_array[0]-one_master_nx_array[0]*lx-xlo)/lcx)+floor((one_master_y_array[0]-one_master_ny_array[0]*ly-ylo)/lcy)*Mx+floor((one_master_z_array[0]-one_master_nz_array[0]*lz-zlo)/lcz)*M2;
                            //printf("%lf\t%lf\t%lf\n",one_master_x_array[0],one_master_y_array[0],one_master_z_array[0]);
                            //printf("%d\t%d\t%d\n",one_master_nx_array[0],one_master_ny_array[0],one_master_nz_array[0]);
                            //printf("cindex: %d\n",cindex);
                            
                            if(unique_cells==0)
                            {
                                unique_cells=1;
                                unique_cells_array[0]=cindex;
                                new_atoms_array[0]=h;
                            }
                            else
                            {
                                new_atoms_array[new_atoms_counter-1]=h;
                                
                                unique_found=0;
                                for(l=0;l<unique_cells;++l)
                                {
                                    if(unique_cells_array[l]==cindex){unique_found=1;break;}
                                }
                                if(unique_found==0)
                                {
                                    unique_cells=unique_cells+1;
                                    unique_cells_array[unique_cells-1]=cindex;
                                }
                                
                            }
                            
                            //local_cells[0]=cindex;
                            //for(k=0;k<26;++k)local_cells[k+1]=neighbors[k][cindex];
                            //printf("local_cells:\n");
                            //for(k=0;k<27;++k)printf("%d\t",local_cells[k]);printf("\n");
                            
                            
                            
                        }
                        
                        
                    }
                    
                    //for(j=0;j<new_atoms_counter;++j)printf("[%d]\t%d\n",j+1,new_atoms_array[j]);
                    //for(j=0;j<unique_cells;++j)printf("cell {%d}\n",unique_cells_array[j]);
                    
                    cells_registry_index=0;
                    cells_registry[cells_registry_index]=unique_cells_array[0];
                    for(j=0;j<26;++j){
                        cells_registry_index=cells_registry_index+1;
                        cells_registry[cells_registry_index]=neighbors[j][unique_cells_array[0]];
                    }
                    for(j=1;j<unique_cells;++j)
                    {
                        unique_found=0;
                        for(l=0;l<cells_registry_index+1;++l)if(cells_registry[l]==unique_cells_array[j]){unique_found=1;break;}
                        if(unique_found==0){
                            cells_registry_index=cells_registry_index+1;
                            cells_registry[cells_registry_index]=unique_cells_array[j];}
                        for(l=0;l<26;++l){
                            eligible[l]=1;
                            for(m=0;m<cells_registry_index+1;++m)
                            {
                                if(cells_registry[m]==neighbors[l][unique_cells_array[j]])eligible[l]=0;
                            }
                            //printf("[%d]\t%d\n",neighbors[l][unique_cells_array[j]],eligible[l]);
                        }
                        
                        for(l=0;l<26;++l)
                        {
                            if(eligible[l]==1)
                            {
                                cells_registry_index=cells_registry_index+1;
                                cells_registry[cells_registry_index]=neighbors[l][unique_cells_array[j]];
                            }
                        }
                        
                    }
                    
                    //for(j=0;j<cells_registry_index+1;++j)printf("[%d]\t-- %d --\n",j+1,cells_registry[j]);
                    
                    //----
                    
                    // here we will carry out the pe calculation
                    
                    //==============================================================
                    //
                    // pe calculation
                    //
                    //.......................................................................... can be placed in loop ...
                    // the latter arrays must be resized if delta_atoms!=0
                    calc_wrapped_coords_n_arrays(general_atoms,xlo,ylo,zlo,xhi,yhi,zhi,lx,ly,lz,master_x_array,master_y_array,master_z_array,master_nx_array,master_ny_array,master_nz_array);
                    //--------------------------------------------------------------------------
                    calc_head_list(general_atoms,Mx,Mx*My,M3,xlo,ylo,zlo,xhi,yhi,zhi,lx,ly,lz,lcx,lcy,lcz,
                                   master_x_array,master_y_array,master_z_array,master_nx_array,master_ny_array,master_nz_array,head,&list);
                    
                    // find atoms
                    
                    cellpart=0; // cellpart is the atom counter per cell
                    
                    for(j=0;j<cells_registry_index+1;++j){
                    cellindex=cells_registry[j];
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
                    }
                    //printf("%d\n",cellpart);
                    //printf("%d\n",general_atoms);
                    
                    //for(j=0;j<cellpart;++j)printf("[%d]\t%d\n",j+1,cellpartarray[j]);
                    
                    //##########################################################
                    // cellpartarray (length: cellpart) holds all atoms that are eligible for interaction with the newly added atoms
                    // the newly added atoms are stored in new_atoms_array (length: new_atoms_counter)

                    
                    local_delta_pe_array[k]=lj_pe(cellpart,new_atoms_counter,lx,ly,lz,cutoff,
                                                  cellpartarray,new_atoms_array,host_id_array,one_two,one_three,one_four,nb_type_array,nb_FF_for_lj,
                                                  master_x_array,master_y_array,master_z_array,master_nx_array,master_ny_array,master_nz_array,
                                                  &pairs);
                    
                    
                    
                    
                    //##########################################################
                    
                    
                    //--------------------------------------------------------------------------
                    //pe_bf=lj_pe_brute_force(general_atoms,lx,ly,lz,cutoff,one_two,one_three,one_four,nb_FF_for_lj,nb_type_array,
                    //                        master_x_array,master_y_array,master_z_array,master_nx_array,master_ny_array,master_nz_array,
                    //                        &pairs_bf,host_id_array);
                    //printf("    pe:\t%lf\tbrute\tpairs: %d\n",pe_bf,pairs_bf);
                    //--------------------------------------------------------------------------
                    if(check_with_cells==1){
                        pe_llc=lj_pe_linked_list_cell(general_atoms,lx,ly,lz,cutoff,one_two,one_three,one_four,nb_FF_for_lj,nb_type_array,
                                                  master_x_array,master_y_array,master_z_array,master_nx_array,master_ny_array,master_nz_array,M3,head,list,neighbors,
                                                  &pairs_llc,host_id_array);
                        local_pe_array[k]=pe_llc;}
                    //pe_array[counter]=pe_llc;
                    //printf("    pe:\t%lf\tcells\tpairs: %d\tdiff: %e%%\n",pe_llc,pairs_llc,(pe_llc/pe_bf-1.0)*100.0);
                    //--------------------------------------------------------------------------
                    // free memory: due to malloc in calc_head_list()
                    free(list);
                    //.......................................................................... end of code that can be placed in loop ...
                    
                    
                    //
                    deltaU=local_delta_pe_array[k];
                    boltzmann_w[k]=exp(-deltaU/(kB*T))*BB_d_angle_pdf_value[k];
                    N=N+boltzmann_w[k];
                    //
                    
                    if(check_with_cells==1){deltaU_test=local_pe_array[k]-pe_array[counter-1];if(fabs(deltaU-deltaU_test)>1.0e-8){printf("%e\n",deltaU-deltaU_test);exit(-5);}}

                    
                    //
                    //
                    //==============================================================

                    // reset
                    if(tacticity!=0 && ft==1){
                    // this is a devel workaround to augment indices!!
                    //if(ft==1 && k<CBMCG_N_trial-1){
                        *general_t_counter=*general_t_counter-1;
                        tacticity_counter_array[i]=tacticity_counter_array[i]-1;
                        //getchar();
                    }
                    
                    
                    
                    

                    //fprintf(fp,"%d\n\n",local_atoms);
                    //for(j=0;j<local_atoms;++j)fprintf(fp,"%c\t%lf\t%lf\t%lf\n",master_species_array[local_map[j]+atom_scaling_array[i]-1][0],master_x_array[local_map[j]+atom_scaling_array[i]-1],master_y_array[local_map[j]+atom_scaling_array[i]-1],master_z_array[local_map[j]+atom_scaling_array[i]-1]);
                    
                    //
                    
                } // k-CBMCG_N_trial loop
                
                // free memory
                for(j=0;j<general_atom_types;++j)free(nb_FF_for_lj[j]);free(nb_FF_for_lj);
                for(j=0;j<general_atoms;++j)free(one_two[j]);free(one_two);
                for(j=0;j<general_atoms;++j)free(one_three[j]);free(one_three);
                for(j=0;j<general_atoms;++j)free(one_four[j]);free(one_four);
                free(nb_type_array);
                free(master_nx_array);free(master_ny_array);free(master_nz_array);
                
                // here we will evaluate all k-loop results
                //if(check_t==0)printf("$$ tacticity will be dealt in CBMC loop\n");
                
                
                N=1.0/N;
                if(debug==1){
                printf("N: %e\n",N);
                printf("counter: %d\n",counter);
                //printf("pe[%d]: %lf\n",counter-1,pe_array[counter-1]);
                printf("[k]\tphi\tPDF\tpe\tD_pe\tD_pe_pair\t|diff|\tw\tfinal\n");
                }
                sum=0.0;
                for(k=0;k<CBMCG_N_trial;++k)
                {
                    
                    sum=sum+boltzmann_w[k];
                    final[k]=N*sum;
                    if(debug==1)printf("[%d]\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%e\t%lf\n",
                                       k+1,BB_d_angle_value[k],
                                       BB_d_angle_pdf_value[k],
                                       
                                       0.0,
                                       //local_pe_array[k],
                                       0.0,
                                       //local_pe_array[k]-pe_array[counter-1],
                                       local_delta_pe_array[k],
                                       0.0,
                                       //fabs(local_pe_array[k]-pe_array[counter-1]-local_delta_pe_array[k]),
                                       boltzmann_w[k],
                                       final[k]);
                
                
                
                }
                
                //fclose(fp);
                //printf("preview ready!\n");
                
                ksi=(double)rand()/RAND_MAX;
                for(k=0;k<CBMCG_N_trial;++k)
                {
                    if(ksi<final[k])
                    {
                        selected_k=k;break;
                    }
                }
                
                k=selected_k;

                if(debug==1){
                    printf("ksi: %lf\n",ksi);
                    printf("selected state: %d\n",selected_k+1);
                    printf("pairs: %d\n",pairs);
                    printf("[%d]\t%lf\t%lf\t%lf\t%lf\t%e\t%lf\n",
                           k+1,
                           BB_d_angle_value[k],
                           BB_d_angle_pdf_value[k],
                           0.0,
                           //local_pe_array[k],
                           0.0,
                           //local_pe_array[k]-pe_array[counter-1],
                           boltzmann_w[k],
                           final[k]);
                    
                    //getchar();
                    
                }
                // alter coords: conformation and tacticity
                
                bb_zmap=BB_d_angle_value[k];
                
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
                               local_x_backup[j],
                               local_y_backup[j],
                               local_z_backup[j],
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
                               local_x_backup[j],
                               local_y_backup[j],
                               local_z_backup[j],
                               &master_x_array[local_map[j]+atom_scaling_array[i]-1],
                               &master_y_array[local_map[j]+atom_scaling_array[i]-1],
                               &master_z_array[local_map[j]+atom_scaling_array[i]-1]);
                    }
                    // rotate the rest 1-2 atoms
                    //if(local_topo_array[local_map[j]-1]==2 && strcmp(segment_atom_name_core[local_map[j]-1],backbone_id[k])!=0)
                    if(local_topo_array[local_map[j]-1]==2)
                    {
                        is_bb=0;
                        for(o=0;o<backbone_type_counter;++o)
                        {
                            if(strcmp(segment_atom_name_core[local_map[j]-1],backbone_type_array[o])==0)
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
                                   local_x_backup[j],
                                   local_y_backup[j],
                                   local_z_backup[j],
                                   &master_x_array[local_map[j]+atom_scaling_array[i]-1],
                                   &master_y_array[local_map[j]+atom_scaling_array[i]-1],
                                   &master_z_array[local_map[j]+atom_scaling_array[i]-1]);
                        }
                    }
                }
                
                if(tacticity!=0)
                force_tacticity(i,tacticity,mol_type_array,tacticity_search,local_atoms,segment_atom_name_core,local_map,local_topo_array,
                                atom_scaling_array,backbone_id,yet_one_more_t_array,backbones,backbone_type_counter,backbone_type_array,
                                
                                molecules,tacticity_array_length,
                                
                                tacticity_counter_array,
                                //&general_t_counter,t_tracker_1,t_tracker_2,t_tracker_3,t_tracker_4,t_tracker_mol,
                                general_t_counter,t_tracker_1,t_tracker_2,t_tracker_3,t_tracker_4,t_tracker_mol,
                                tacticity_array,master_x_array,master_y_array,master_z_array,
                                &ft);
                
            } // identified_D if
        } // f1,f2,f3 flags if
    } // found flag if
    if(check_t==1)
    {
        if(tacticity!=0){
        //printf("$$ tacticity will be dealt outside CBMC loop\n");//getchar();
        force_tacticity(i,tacticity,mol_type_array,tacticity_search,local_atoms,segment_atom_name_core,local_map,local_topo_array,
                        atom_scaling_array,backbone_id,yet_one_more_t_array,backbones,backbone_type_counter,backbone_type_array,
                        
                        molecules,tacticity_array_length,
                        
                        tacticity_counter_array,
                        //&general_t_counter,t_tracker_1,t_tracker_2,t_tracker_3,t_tracker_4,t_tracker_mol,
                        general_t_counter,t_tracker_1,t_tracker_2,t_tracker_3,t_tracker_4,t_tracker_mol,
                        tacticity_array,master_x_array,master_y_array,master_z_array,
                        &ft);
        }
    }
}

// cellpartarray (length: cellpart) holds all atoms that are eligible for interaction with the newly added atoms
// the newly added atoms are stored in new_atoms_array (length: new_atoms_counter)

double lj_pe(int cellpart,int new_atoms_counter,double lx,double ly,double lz,double cutoff,
             int *cellpartarray,int *new_atoms_array,int *host_id_array,int **one_two,int **one_three,int **one_four,int *nb_type_array,double **nb_FF,
             double *master_x_array,double *master_y_array,double *master_z_array,int *master_nx_array,int *master_ny_array,int *master_nz_array,
             int *pairs_llc)
{
    
    int i,j,k,locali,localj,select,one_four_flag;
    double xi,yi,zi,xj,yj,zj,r2,r,epsilon,sigma,r6,r12,sigma2,sigma6,sigma12,pe;
    
    pe=0.0;*pairs_llc=0;
    
    for (i=0;i<cellpart;++i)
    {
        for (j=0;j<new_atoms_counter;++j)
        {
            
            locali=cellpartarray[i];
            localj=new_atoms_array[j];
            
            if(locali!=localj && (host_id_array[locali]!=-2 && host_id_array[localj]!=-2))
            {
                
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
    
    return pe*2.0;
    
}


