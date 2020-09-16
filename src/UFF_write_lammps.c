
#include"builder.h"

int equal(double A,double B);

void UFF_write_lammps(char *current_folder,int counter,int l,
                      int general_atoms,int general_bonds,int general_angles,int general_dihedrals,int general_impropers,
                      int general_atom_types,int intermed_bond_types,int intermed_angle_types,int intermed_improper_types,
                      char **atom_type_core_FF,char **species_intermed_array,int *D_multi,
                      double xlo,double xhi,double ylo,double yhi,double zlo,double zhi,
                      int *master_id_array,int *master_mol_array,double *master_q_array,double *master_x_array,double *master_y_array,double *master_z_array,
                      int **general_bonds_registry,int **general_angles_registry,int **general_dihedrals_registry,int **general_impropers_registry,
                      char ***global_bond_type_intermed_array,char ***global_angle_type_intermed_array,char ***global_dihedral_type_intermed_array,char ***global_improper_type_intermed_array,
                      double **nb_FF,double **bonds_FF,double **angles_FF,double **dihedrals_FF,double **impropers_FF,
                      int **b_refine_array, int **a_refine_array, int **d_refine_array, int **i_refine_array)
{
    FILE *fp,*fp2;
    char file_path[cmax_length],buffer[cmax_length];
    int j,n,D_type_counter,found,igl,int_buffer;
    int *type_core;
    int *D_type,*reg1,*reg2;
    char **dihedrals_multi_text;
    double **dihedrals_FF_multi;
    double mass,ljconvert;
    double lx,ly,lz;
    double nx,ny,nz;
    
    // define supercell edges length; needed for periodic wrapping flags
    lx=xhi-xlo;
    ly=yhi-ylo;
    lz=zhi-zlo;
    
    // prepare the type array
    type_core=(int*)malloc(general_atoms*sizeof(int));
    for(j=0;j<general_atoms;++j)
    {
        for(n=0;n<general_atom_types;++n)
            if(strcmp(atom_type_core_FF[j],species_intermed_array[n])==0)
                type_core[j]=n+1;
    }

    //
    
    // dihedral multiplicity workaround
    D_type=(int*)malloc(general_dihedrals*sizeof(int));
    //for(j=0;j<general_dihedrals;++j)D_type[j]=general_dihedrals_registry[j][0];
    for(j=0;j<general_dihedrals;++j)D_type[j]=(*d_refine_array)[general_dihedrals_registry[j][0]-1];
    
    D_type_counter=1;
    reg1=(int*)malloc(D_type_counter*sizeof(int));
    reg2=(int*)malloc(D_type_counter*sizeof(int));
    reg1[0]=D_type[0];
    reg2[0]=D_multi[0];
    for(j=1;j<general_dihedrals;++j)
    {
        found=0;
        for(n=0;n<D_type_counter;++n)
        {
            if(D_type[j]==reg1[n] && D_multi[j]==reg2[n]){found=1;break;}
        }
        if(found==0)
        {
            D_type_counter=D_type_counter+1;
            reg1=(int*)realloc(reg1,D_type_counter*sizeof(int));
            reg2=(int*)realloc(reg2,D_type_counter*sizeof(int));
            reg1[D_type_counter-1]=D_type[j];
            reg2[D_type_counter-1]=D_multi[j];
        }
    }
    
    dihedrals_multi_text=(char**)malloc(D_type_counter*sizeof(char*));
    for(j=0;j<D_type_counter;++j)dihedrals_multi_text[j]=(char*)malloc(cmax_length*sizeof(char));
    
    dihedrals_FF_multi=(double**)malloc(D_type_counter*sizeof(double*));
    for(j=0;j<D_type_counter;++j)dihedrals_FF_multi[j]=(double*)malloc(3*sizeof(double));
    
    for(j=0;j<D_type_counter;++j)
    {
        sprintf(dihedrals_multi_text[j],"(%s)--%s--(%s)--%s--(%s)--%s--(%s) [No. %d - multi: %d]",
                global_dihedral_type_intermed_array[reg1[j]-1][0],
                global_dihedral_type_intermed_array[reg1[j]-1][1],
                global_dihedral_type_intermed_array[reg1[j]-1][2],
                global_dihedral_type_intermed_array[reg1[j]-1][3],
                global_dihedral_type_intermed_array[reg1[j]-1][4],
                global_dihedral_type_intermed_array[reg1[j]-1][5],
                global_dihedral_type_intermed_array[reg1[j]-1][6],
                reg1[j],reg2[j]);
        dihedrals_FF_multi[j][0]=dihedrals_FF[reg1[j]-1][0]/((double)reg2[j]);
        dihedrals_FF_multi[j][1]=dihedrals_FF[reg1[j]-1][1];
        dihedrals_FF_multi[j][2]=dihedrals_FF[reg1[j]-1][2];
    }
    for(j=0;j<general_dihedrals;++j)
    {
        for(n=0;n<D_type_counter;++n)
        {
            if(D_type[j]==reg1[n] && D_multi[j]==reg2[n])
            {
                D_type[j]=n+1;
                break;
            }
        }
    }
    
    //
    
    //
    sprintf(file_path,"%s/mol_%d.lammps",current_folder,l+1);
    fp=fopen(file_path,"w+");
    fprintf(fp,"mol_%d_%d\n",counter+1,l+1);
    fprintf(fp,"\n");
    fprintf(fp,"%d atoms\n",general_atoms);
    fprintf(fp,"%d bonds\n",general_bonds);
    fprintf(fp,"%d angles\n",general_angles);
    fprintf(fp,"%d dihedrals\n",general_dihedrals);
    if(general_impropers>0)fprintf(fp,"%d impropers\n",general_impropers);
    fprintf(fp,"\n");
    fprintf(fp,"%d atom types\n",general_atom_types);
    fprintf(fp,"%d bond types\n",intermed_bond_types);
    fprintf(fp,"%d angle types\n",intermed_angle_types);
    fprintf(fp,"%d dihedral types\n",D_type_counter);
    if(general_impropers>0)fprintf(fp,"%d improper types\n",intermed_improper_types);
    fprintf(fp,"\n");
    fprintf(fp,"%lf %lf xlo xhi\n",xlo,xhi);
    fprintf(fp,"%lf %lf ylo yhi\n",ylo,yhi);
    fprintf(fp,"%lf %lf zlo zhi\n",zlo,zhi);
    fprintf(fp,"\n");
    fprintf(fp,"Masses\n");
    fprintf(fp,"\n");
    for(j=0;j<general_atom_types;++j)
    {
        mass=0.0;
        if(species_intermed_array[j][0]=='H' && species_intermed_array[j][1]=='_')mass=mass_H;
        if(species_intermed_array[j][0]=='C' && species_intermed_array[j][1]=='_')mass=mass_C;
        if(species_intermed_array[j][0]=='O' && species_intermed_array[j][1]=='_')mass=mass_O;
        if(species_intermed_array[j][0]=='N' && species_intermed_array[j][1]=='_')mass=mass_N;
        fprintf(fp,"%d\t%lf\n",j+1,mass);
    }
    fprintf(fp,"\n");
    fprintf(fp,"Atoms\n");
    fprintf(fp,"\n");
    for(j=0;j<general_atoms;++j)
    {
        nx=0.0;ny=0.0;nz=0.0;
        
        if(master_x_array[j]>xhi)
        {
            nx=(master_x_array[j]-xlo)/lx;
            nx=floor(nx);
        }
        if(master_x_array[j]<xlo)
        {
            nx=(xhi-master_x_array[j])/lx;
            nx=-floor(nx);
        }
        
        if(master_y_array[j]>yhi)
        {
            ny=(master_y_array[j]-ylo)/ly;
            ny=floor(ny);
        }
        if(master_y_array[j]<ylo)
        {
            ny=(yhi-master_y_array[j])/ly;
            ny=-floor(ny);
        }
        
        if(master_z_array[j]>zhi)
        {
            nz=(master_z_array[j]-zlo)/lz;
            nz=floor(nz);
        }
        if(master_z_array[j]<zlo)
        {
            nz=(zhi-master_z_array[j])/lz;
            nz=-floor(nz);
        }
        
        fprintf(fp,"%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%d\t%d\t%d\n",master_id_array[j],master_mol_array[j],type_core[j],master_q_array[j],master_x_array[j]-nx*lx,master_y_array[j]-ny*ly,master_z_array[j]-nz*lz,(int)nx,(int)ny,(int)nz);
    }
    
    fprintf(fp,"\n");
    fprintf(fp,"Bonds\n");
    fprintf(fp,"\n");
    for(igl=0;igl<general_bonds;++igl)
        //fprintf(fp,"%d\t%d\t%d\t%d\n",igl+1,general_bonds_registry[igl][0],general_bonds_registry[igl][1],general_bonds_registry[igl][2]);
        fprintf(fp,"%d\t%d\t%d\t%d\n",igl+1,(*b_refine_array)[general_bonds_registry[igl][0]-1],general_bonds_registry[igl][1],general_bonds_registry[igl][2]);
    fprintf(fp,"\n");
    fprintf(fp,"Angles\n");
    fprintf(fp,"\n");
    for(igl=0;igl<general_angles;++igl)
        //fprintf(fp,"%d\t%d\t%d\t%d\t%d\n",igl+1,general_angles_registry[igl][0],general_angles_registry[igl][1],general_angles_registry[igl][2],general_angles_registry[igl][3]);
        fprintf(fp,"%d\t%d\t%d\t%d\t%d\n",igl+1,(*a_refine_array)[general_angles_registry[igl][0]-1],general_angles_registry[igl][1],general_angles_registry[igl][2],general_angles_registry[igl][3]);
    fprintf(fp,"\n");
    fprintf(fp,"Dihedrals\n");
    fprintf(fp,"\n");
    for(igl=0;igl<general_dihedrals;++igl)
        fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\n",igl+1,D_type[igl],general_dihedrals_registry[igl][1],general_dihedrals_registry[igl][2],general_dihedrals_registry[igl][3],general_dihedrals_registry[igl][4]);
    fprintf(fp,"\n");
    if(general_impropers>0)
    {
        fprintf(fp,"Impropers\n");
        fprintf(fp,"\n");
        for(igl=0;igl<general_impropers;++igl)
            fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\n",igl+1,(*i_refine_array)[general_impropers_registry[igl][0]-1],general_impropers_registry[igl][1],general_impropers_registry[igl][2],general_impropers_registry[igl][3],general_impropers_registry[igl][4]);
            //fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\n",igl+1,general_impropers_registry[igl][0],general_impropers_registry[igl][1],general_impropers_registry[igl][2],general_impropers_registry[igl][3],general_impropers_registry[igl][4]);
        fprintf(fp,"\n");
    }
    
    fclose(fp);
    
    //sprintf(buffer,"cp %s/mol_%d.lammps %s/%d_mol_%d.lammps",current_folder,l+1,current_folder,counter,l+1);
    //system(buffer);
    
    //
    
    sprintf(file_path,"%s/mol_%d.in",current_folder,l+1);
    fp=fopen(file_path,"w+");
    //fprintf(fp,"\nvariable E_tol equal 1.0e-5\nvariable F_tol equal 1.0e-4\nvariable N_iter equal 20000\nvariable N_eval equal 50000\n\n");
    fprintf(fp,"\nvariable E_tol equal 1.0e-4\nvariable F_tol equal 1.0e-3\nvariable N_iter equal 10000\nvariable N_eval equal 25000\n\n");
    //fprintf(fp,"\nvariable E_tol equal 1.0e-14\nvariable F_tol equal 1.0e-14\nvariable N_iter equal 1000000\nvariable N_eval equal 2500000\n\n");
    fprintf(fp,"units real\natom_style full\nboundary p p p\nread_data mol_%d.lammps\n",l+1);
    //fprintf(fp,"timestep 0.5\n");
    for(j=0;j<general_atom_types;++j)fprintf(fp,"# %d --> %s\n",j+1,species_intermed_array[j]);
    fprintf(fp,"group my_atoms type ");
    for(j=0;j<general_atom_types;++j)fprintf(fp,"%d ",j+1);fprintf(fp,"\n");
    fprintf(fp,"bond_style harmonic\n");
    for(j=0;j<intermed_bond_types;++j)
    {
        fprintf(fp,"# (%s)--%s--(%s)\n",global_bond_type_intermed_array[j][0],global_bond_type_intermed_array[j][1],global_bond_type_intermed_array[j][2]);
        fprintf(fp,"bond_coeff %d %lf %lf\n",j+1,0.5*bonds_FF[j][0],bonds_FF[j][1]);
    }
    fprintf(fp,"angle_style hybrid cosine/periodic fourier\n");
    for(j=0;j<intermed_angle_types;++j)
    {
        fprintf(fp,"# (%s)--%s--(%s)--%s--(%s)\n",global_angle_type_intermed_array[j][0],global_angle_type_intermed_array[j][1],global_angle_type_intermed_array[j][2],global_angle_type_intermed_array[j][3],global_angle_type_intermed_array[j][4]);
        if(angles_FF[j][1]==0)
        {
            fprintf(fp,"angle_coeff %d fourier %lf %lf %lf %lf\n",j+1,angles_FF[j][0],angles_FF[j][3],angles_FF[j][4],angles_FF[j][5]);
        }
        else
        {
            int_buffer=0;
            if(equal(angles_FF[j][1],1.0)==1)int_buffer=1;
            if(equal(angles_FF[j][1],2.0)==1)int_buffer=-1;
            if(equal(angles_FF[j][1],3.0)==1)int_buffer=-1;
            if(equal(angles_FF[j][1],4.0)==1)int_buffer=1;
            fprintf(fp,"angle_coeff %d cosine/periodic %lf %d %d\n",j+1,0.5*angles_FF[j][0],int_buffer,(int)angles_FF[j][1]);
        }
    }
    fprintf(fp,"dihedral_style harmonic\n");
    for(j=0;j<D_type_counter;++j)
    {
        fprintf(fp,"# %s\n",dihedrals_multi_text[j]);
        fprintf(fp,"dihedral_coeff %d %lf %d %d\n",j+1,0.5*dihedrals_FF_multi[j][0],(int)(-1.0*cos(dihedrals_FF_multi[j][1]*dihedrals_FF_multi[j][2]*pi/180.0)),(int)dihedrals_FF_multi[j][1]);
    }
    if(general_impropers>0)
    {
        fprintf(fp,"improper_style fourier\n");
        for(j=0;j<intermed_improper_types;++j)
        {
            fprintf(fp,"# %s -- %s -- %s -- %s\n",global_improper_type_intermed_array[j][0],global_improper_type_intermed_array[j][1],global_improper_type_intermed_array[j][2],global_improper_type_intermed_array[j][3]);
            //fprintf(fp,"improper_coeff %d %lf %lf %lf %lf %d\n",j+1,impropers_FF[j][0]/3.0,impropers_FF[j][2],impropers_FF[j][3],impropers_FF[j][4],0);
            fprintf(fp,"improper_coeff %d %lf %lf %lf %lf %d\n",j+1,impropers_FF[j][0]/3.0,impropers_FF[j][2],impropers_FF[j][3],impropers_FF[j][4],3);
        }
    }
    ljconvert=pow(2.0,-1.0/6.0);
    fprintf(fp,"pair_style lj/cut 12.0\n");
    for(j=0;j<general_atom_types;++j)
    {
        for(n=j;n<general_atom_types;++n)
        {
            
            fprintf(fp,"pair_coeff %d %d %lf %lf\n",j+1,n+1,sqrt(nb_FF[j][0]*nb_FF[n][0]),ljconvert*sqrt(nb_FF[j][1]*nb_FF[n][1]));
        }
    }
    fprintf(fp,"pair_modify mix geometric\nspecial_bonds lj 0.0 0.0 1.0\n\n");
    
    //
    sprintf(file_path,"%s/min_header.txt",current_folder);
    fp2=fopen(file_path,"r");
    if(fp2!=NULL)
    {
        while(fgets(buffer,cmax_length,fp2)!=NULL)
            fprintf(fp,"%s",buffer);
        fclose(fp2);
    }
    
    
    //
    
    fprintf(fp,
            "thermo 1\n"
            "thermo_style custom step pe ebond eangle edihed eimp evdwl\n"
            "min_style cg\n"
            "minimize ${E_tol} ${F_tol} ${N_iter} ${N_eval}\n"
            "reset_timestep 0\n\n");
    
    fprintf(fp,"# write snapshot\n"
            "dump coords my_atoms custom 1 snapshot_%d.dat id mol type q x y z ix iy iz\n"
            "dump_modify coords sort id\n\n",l+1);
    
    //dump xyz_coords all xyz ${xyz_step} xyz_cg_${system_name}/animation.xyz
    //dump_modify xyz_coords element C C C H C C
    
    
    
    fprintf(fp,"# write preview\n"
            "dump xyz_coords my_atoms xyz 1 animation.xyz\n"
            "dump_modify xyz_coords element ");
    for(j=0;j<general_atom_types;++j)fprintf(fp,"%c ",species_intermed_array[j][0]);
    fprintf(fp,"\n\n");
    
    fprintf(fp,
            "run 0\n\n"
            );
    
    fclose(fp);
    
    //sprintf(buffer,"cp %s/mol_%d.in %s/%d_mol_%d.in",current_folder,l+1,current_folder,counter,l+1);
    //system(buffer);
    
    //
    
    free(type_core);
    free(D_type);free(reg1);free(reg2);
    for(j=0;j<D_type_counter;++j)free(dihedrals_FF_multi[j]);free(dihedrals_FF_multi);
    for(j=0;j<D_type_counter;++j)free(dihedrals_multi_text[j]);free(dihedrals_multi_text);
    
    //free(*b_refine_array);free(*a_refine_array);free(*d_refine_array);free(*i_refine_array);

}
