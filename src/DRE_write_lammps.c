
#include"builder.h"

int equal(double A,double B);

void find_unique_d(int N, int M, double **array,int **out);

void DRE_write_lammps(char *current_folder,int counter,int l,
                      int general_atoms,int general_bonds,int general_angles,int general_dihedrals,int general_impropers,
                      int general_atom_types,int intermed_bond_types,int intermed_angle_types,int intermed_dihedral_types,int intermed_improper_types,
                      char **atom_type_core_FF,char **species_intermed_array,
                      double xlo,double xhi,double ylo,double yhi,double zlo,double zhi,
                      int *master_id_array,int *master_mol_array,double *master_q_array,double *master_x_array,double *master_y_array,double *master_z_array,
                      int **general_bonds_registry,int **general_angles_registry,int **general_dihedrals_registry,int **general_impropers_registry,
                      char ***global_bond_type_intermed_array,char ***global_angle_type_intermed_array,char ***global_dihedral_type_intermed_array,char ***global_improper_type_intermed_array,
                      double **nb_FF,double **bonds_FF,double **angles_FF,double **dihedrals_FF,double **impropers_FF,
                      int **b_refine_array,int **a_refine_array,int **d_refine_array,int **i_refine_array,
                      int *D_multi,
                      int equiv_B_N,char **equiv_B,int equiv_A_N,char **equiv_A,int equiv_D_N,char **equiv_D,
                      char ***general_bond_types_registry,char ***general_angle_types_registry,char ***general_dihedral_types_registry)
{
    FILE *fp;//,*fp2;
    char file_path[cmax_length];//,buffer[cmax_length];
    int j,k,n,D_type_counter,found,igl;//,int_buffer;
    int *type_core;
    int *D_type,*reg1,*reg2;
    char **dihedrals_multi_text;
    double **dihedrals_FF_multi;
    double mass,ljconvert;
    double lx,ly,lz;
    double nx,ny,nz;
    
    int m1,m2,m,o;
    
    int equiv_t,int_buffer;
    char el[7][cmax_length];
    char el2[7][cmax_length];
    char el3[7][cmax_length];
    double f1,f2,f3,f4,**dihedrals_FF_new;int **dihedrals_out;
    char **new_D_registry;int new_D_types;int new_D_type_array[general_dihedrals];
    
    int equiv_N;
    
    // use type defined dihedral multiplicities
    for(j=0;j<general_dihedrals;++j)
    {
        //printf("[%d]\t%d\t%d\t%s\t%s\n",j+1,general_dihedrals_registry[j][1],general_dihedrals_registry[j][2],
        //       atom_type_core_FF[general_dihedrals_registry[j][1]-1],atom_type_core_FF[general_dihedrals_registry[j][2]-1]);
        m1=1;
        if(strcmp(atom_type_core_FF[general_dihedrals_registry[j][2]-1],"H_")!=0)m1=atom_type_core_FF[general_dihedrals_registry[j][2]-1][2] - '0';
        m2=1;
        if(strcmp(atom_type_core_FF[general_dihedrals_registry[j][3]-1],"H_")!=0)m2=atom_type_core_FF[general_dihedrals_registry[j][3]-1][2] - '0';
        //printf("{%d %d}\n",m1,m2);
        D_multi[j]=m1*m2;
    }
    //getchar();
    
    new_D_types=0;
    
    for(j=0;j<general_dihedrals;++j)
    {
        // convert from SYBYL to DREIDING types; no overwrite!
        // stored in el
        for(k=0;k<7;++k)
        {
            if(strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][k],"C.3")==0)sprintf(el[k],"%s","C_3");
            else if(strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][k],"C.2")==0)sprintf(el[k],"%s","C_2");
            else if(strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][k],"O.3")==0)sprintf(el[k],"%s","O_3");
            else if(strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][k],"O.2")==0)sprintf(el[k],"%s","O_2");
            else if(strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][k],"H")==0)sprintf(el[k],"%s","H_");
            else if(strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][k],"1")==0)sprintf(el[k],"%s","1");
            else if(strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][k],"2")==0)sprintf(el[k],"%s","2");
            else if(strcmp(general_dihedral_types_registry[general_dihedrals_registry[j][0]-1][k],"ar")==0)sprintf(el[k],"%s","ar");
        }
        
        // find in equiv registry
        // locate the type
        for(k=0;k<equiv_D_N;++k)
        {
            // read and store to el2 and fs
            sscanf(equiv_D[k],"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf",
                   el2[0],el2[1],el2[2],el2[3],el2[4],el2[5],el2[6],&int_buffer,&int_buffer,&int_buffer,&f1,&f2,&f3,&f4);
            // check
            if((strcmp(el[0],el2[0])==0 &&
                strcmp(el[1],el2[1])==0 &&
                strcmp(el[2],el2[2])==0 &&
                strcmp(el[3],el2[3])==0 &&
                strcmp(el[4],el2[4])==0 &&
                strcmp(el[5],el2[5])==0 &&
                strcmp(el[6],el2[6])==0 ) ||
               
               (strcmp(el[0],el2[6])==0 &&
                strcmp(el[1],el2[5])==0 &&
                strcmp(el[2],el2[4])==0 &&
                strcmp(el[3],el2[3])==0 &&
                strcmp(el[4],el2[2])==0 &&
                strcmp(el[5],el2[1])==0 &&
                strcmp(el[6],el2[0])==0 )
               )
            {break;}
            
        }
        
        if(new_D_types==0)
        {
            // write the first entry in the NEW registry
            new_D_types=new_D_types+1;
            new_D_registry=(char**)malloc(new_D_types*sizeof(char*));
            new_D_registry[new_D_types-1]=(char*)malloc(cmax_length*sizeof(char));
            sprintf(new_D_registry[new_D_types-1],"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%lf\t%lf\t%lf\t%lf",
                    el2[0],el2[1],el2[2],el2[3],el2[4],el2[5],el2[6],f1/D_multi[j],f2,f3,f4);
            new_D_type_array[j]=new_D_types;
        }
        else
        {
            // try to locate the type in the NEW registry
            found=0;
            for(k=0;k<new_D_types;++k)
            {
                // read and store to el2 and fs
                sscanf(new_D_registry[k],"%s\t%s\t%s\t%s\t%s\t%s\t%s",
                       el3[0],el3[1],el3[2],el3[3],el3[4],el3[5],el3[6]);
                // check
                if((strcmp(el2[0],el3[0])==0 &&
                    strcmp(el2[1],el3[1])==0 &&
                    strcmp(el2[2],el3[2])==0 &&
                    strcmp(el2[3],el3[3])==0 &&
                    strcmp(el2[4],el3[4])==0 &&
                    strcmp(el2[5],el3[5])==0 &&
                    strcmp(el2[6],el3[6])==0 ) ||
                   
                   (strcmp(el2[0],el3[6])==0 &&
                    strcmp(el2[1],el3[5])==0 &&
                    strcmp(el2[2],el3[4])==0 &&
                    strcmp(el2[3],el3[3])==0 &&
                    strcmp(el2[4],el3[2])==0 &&
                    strcmp(el2[5],el3[1])==0 &&
                    strcmp(el2[6],el3[0])==0 )
                   )
                {
                    found=1;
                    
                    new_D_type_array[j]=k+1;
                    
                    break;
                }
            }
            // if not found:
            if(found==0)
            {
                // ADD!
                new_D_types=new_D_types+1;
                new_D_registry=(char**)realloc(new_D_registry,new_D_types*sizeof(char*));
                new_D_registry[new_D_types-1]=(char*)malloc(cmax_length*sizeof(char));
                sprintf(new_D_registry[new_D_types-1],"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%lf\t%lf\t%lf\t%lf",
                        el2[0],el2[1],el2[2],el2[3],el2[4],el2[5],el2[6],f1/D_multi[j],f2,f3,f4);
                
                new_D_type_array[j]=new_D_types;
                
            
            }
        
        }
    
    
    }
    
    dihedrals_FF_new=(double**)malloc(new_D_types*sizeof(double*));
    for(j=0;j<new_D_types;++j)dihedrals_FF_new[j]=(double*)malloc(4*sizeof(double));
    for(j=0;j<new_D_types;++j){sscanf(new_D_registry[j],"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%lf\t%lf\t%lf\t%lf",
                                      el2[0],el2[1],el2[2],el2[3],el2[4],el2[5],el2[6],
                                      &dihedrals_FF_new[j][0],&dihedrals_FF_new[j][1],&dihedrals_FF_new[j][2],&dihedrals_FF_new[j][3]);}
    
    // refine dihedrals
    dihedrals_out=(int**)malloc(new_D_types*sizeof(int*));
    for(j=0;j<new_D_types;++j)dihedrals_out[j]=(int*)malloc(3*sizeof(int));
    find_unique_d(new_D_types,4,dihedrals_FF_new,dihedrals_out);
    
    //for(j=0;j<equiv_D_N;++j)printf("[%d]\t%s\n",j+1,equiv_D[j]);
    
    //for(j=0;j<new_D_types;++j)printf("{%d}\t%s\t%d\t%d\t%d\n",j+1,new_D_registry[j],dihedrals_out[j][0],dihedrals_out[j][1],dihedrals_out[j][2]);
    
    
    
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
    /*
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
    for(j=0;j<D_type_counter;++j)dihedrals_FF_multi[j]=(double*)malloc(4*sizeof(double));
    
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
        dihedrals_FF_multi[j][3]=dihedrals_FF[reg1[j]-1][3];
        
        //printf("[%d]\t{%d}\t|%lf|\t%lf\t%lf\t%lf\t%lf\n",j+1,reg2[j],dihedrals_FF[reg1[j]-1][0],dihedrals_FF_multi[j][0],dihedrals_FF_multi[j][1],dihedrals_FF_multi[j][2],dihedrals_FF_multi[j][3]);
        
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
    */
    //
    
    //
    sprintf(file_path,"%s/mol_DREIDING_%d.lammps",current_folder,l+1);
    fp=fopen(file_path,"w+");
    fprintf(fp,"mol_DREIDING_%d_%d\n",counter+1,l+1);
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
//    fprintf(fp,"%d dihedral types\n",D_type_counter);
    k=0;for(j=0;j<new_D_types;++j)if(dihedrals_out[j][2]!=-1)k=k+1;
    fprintf(fp,"%d dihedral types\n",k);
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
        //fprintf(fp,"%d\t%d\t%d\t%d\t#\t(%s)--%s--(%s)\n",igl+1,(*b_refine_array)[general_bonds_registry[igl][0]-1],general_bonds_registry[igl][1],general_bonds_registry[igl][2],
        fprintf(fp,"%d\t%d\t%d\t%d\t#\t%s\t%s\t%s\n",igl+1,(*b_refine_array)[general_bonds_registry[igl][0]-1],general_bonds_registry[igl][1],general_bonds_registry[igl][2],
                general_bond_types_registry[general_bonds_registry[igl][0]-1][0],
                general_bond_types_registry[general_bonds_registry[igl][0]-1][1],
                general_bond_types_registry[general_bonds_registry[igl][0]-1][2]);
    fprintf(fp,"\n");
    fprintf(fp,"Angles\n");
    fprintf(fp,"\n");
    for(igl=0;igl<general_angles;++igl)
        //fprintf(fp,"%d\t%d\t%d\t%d\t%d\t#\t(%s)--%s--(%s)--%s--(%s)\n",igl+1,(*a_refine_array)[general_angles_registry[igl][0]-1],general_angles_registry[igl][1],general_angles_registry[igl][2],general_angles_registry[igl][3],
        fprintf(fp,"%d\t%d\t%d\t%d\t%d\t#\t%s\t%s\t%s\t%s\t%s\n",igl+1,(*a_refine_array)[general_angles_registry[igl][0]-1],general_angles_registry[igl][1],general_angles_registry[igl][2],general_angles_registry[igl][3],
                general_angle_types_registry[general_angles_registry[igl][0]-1][0],
                general_angle_types_registry[general_angles_registry[igl][0]-1][1],
                general_angle_types_registry[general_angles_registry[igl][0]-1][2],
                general_angle_types_registry[general_angles_registry[igl][0]-1][3],
                general_angle_types_registry[general_angles_registry[igl][0]-1][4]);
    fprintf(fp,"\n");
    fprintf(fp,"Dihedrals\n");
    fprintf(fp,"\n");
    for(igl=0;igl<general_dihedrals;++igl)
        //fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\n",igl+1,D_type[igl],general_dihedrals_registry[igl][1],general_dihedrals_registry[igl][2],general_dihedrals_registry[igl][3],general_dihedrals_registry[igl][4]);
        //fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\t#\t(%s)--%s--(%s)--%s--(%s)--%s--(%s)\n",igl+1,dihedrals_out[new_D_type_array[igl]-1][1],general_dihedrals_registry[igl][1],general_dihedrals_registry[igl][2],general_dihedrals_registry[igl][3],general_dihedrals_registry[igl][4],
        fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\t#\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",igl+1,dihedrals_out[new_D_type_array[igl]-1][1],general_dihedrals_registry[igl][1],general_dihedrals_registry[igl][2],general_dihedrals_registry[igl][3],general_dihedrals_registry[igl][4],
                general_dihedral_types_registry[general_dihedrals_registry[igl][0]-1][0],
                general_dihedral_types_registry[general_dihedrals_registry[igl][0]-1][1],
                general_dihedral_types_registry[general_dihedrals_registry[igl][0]-1][2],
                general_dihedral_types_registry[general_dihedrals_registry[igl][0]-1][3],
                general_dihedral_types_registry[general_dihedrals_registry[igl][0]-1][4],
                general_dihedral_types_registry[general_dihedrals_registry[igl][0]-1][5],
                general_dihedral_types_registry[general_dihedrals_registry[igl][0]-1][6]);
    fprintf(fp,"\n");
    if(general_impropers>0)
    {
        fprintf(fp,"Impropers\n");
        fprintf(fp,"\n");
        for(igl=0;igl<general_impropers;++igl)
            fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%d\n",igl+1,(*i_refine_array)[general_impropers_registry[igl][0]-1],general_impropers_registry[igl][1],general_impropers_registry[igl][2],general_impropers_registry[igl][3],general_impropers_registry[igl][4]);
        fprintf(fp,"\n");
    }
    fclose(fp);
    
    //
    
    sprintf(file_path,"%s/mol_DREIDING_%d.in",current_folder,l+1);
    fp=fopen(file_path,"w+");
    fprintf(fp,"\nvariable E_tol equal 1.0e-4\nvariable F_tol equal 1.0e-3\nvariable N_iter equal 10000\nvariable N_eval equal 25000\n\n");
    fprintf(fp,"units real\natom_style full\nboundary p p p\nread_data mol_DREIDING_%d.lammps\n",l+1);
    for(j=0;j<general_atom_types;++j)fprintf(fp,"# %d --> %s\n",j+1,species_intermed_array[j]);
    fprintf(fp,"group my_atoms type ");
    for(j=0;j<general_atom_types;++j)fprintf(fp,"%d ",j+1);fprintf(fp,"\n");
    fprintf(fp,"\n");
    fprintf(fp,"bond_style harmonic\n");
    fprintf(fp,"\n");
    for(j=0;j<intermed_bond_types;++j)
    {
        
        fprintf(fp,"bond_coeff %d %lf %lf\n",j+1,0.5*bonds_FF[j][0],bonds_FF[j][1]);
        
        equiv_N=0;
        for(k=0;k<equiv_B_N;++k)
        {
            sscanf(equiv_B[k],"%s\t%s\t%s\t%d\t%d",el[0],el[1],el[2],&int_buffer,&equiv_t);
            if(equiv_t==j+1)equiv_N=equiv_N+1;
            
        }
        fprintf(fp,"# %d equiv\n",equiv_N);
        
        for(k=0;k<equiv_B_N;++k)
        {
            sscanf(equiv_B[k],"%s\t%s\t%s\t%d\t%d",el[0],el[1],el[2],&int_buffer,&equiv_t);
            for(m=0;m<3;++m)
            {
                if(strcmp(el[m],"H_")==0)sprintf(el[m],"%s","H");
                for(o=0;o<cmax_length;++o){if(el[m][o]=='\0')break;if(el[m][o]=='_')el[m][o]='.';}
            }
            //if(equiv_t==j+1)fprintf(fp,"# (%s)--%s--(%s)\n",el[0],el[1],el[2]);
            if(equiv_t==j+1)fprintf(fp,"# %s\t%s\t%s\n",el[0],el[1],el[2]);
        }
        
        fprintf(fp,"\n");

    }
    fprintf(fp,"angle_style hybrid cosine/squared cosine\n");
    fprintf(fp,"\n");
    for(j=0;j<intermed_angle_types;++j)
    {
        
        if(equal(angles_FF[j][1],180.0)==1)
        {
            fprintf(fp,"angle_coeff %d cosine %lf\n",j+1,angles_FF[j][0]);
        }
        else
        {
            fprintf(fp,"angle_coeff %d cosine/squared %lf %lf\n",j+1,0.5*angles_FF[j][0],angles_FF[j][1]);
        }
        
        equiv_N=0;
        for(k=0;k<equiv_A_N;++k)
        {
            sscanf(equiv_A[k],"%s\t%s\t%s\t%s\t%s\t%d\t%d",el[0],el[1],el[2],el[3],el[4],&int_buffer,&equiv_t);
            if(equiv_t==j+1)equiv_N=equiv_N+1;
            
        }
        fprintf(fp,"# %d equiv\n",equiv_N);

        
        for(k=0;k<equiv_A_N;++k)
        {
            sscanf(equiv_A[k],"%s\t%s\t%s\t%s\t%s\t%d\t%d",el[0],el[1],el[2],el[3],el[4],&int_buffer,&equiv_t);
            for(m=0;m<5;++m)
            {
                if(strcmp(el[m],"H_")==0)sprintf(el[m],"%s","H");
                for(o=0;o<cmax_length;++o){if(el[m][o]=='\0')break;if(el[m][o]=='_')el[m][o]='.';}
            }
            //if(equiv_t==j+1)fprintf(fp,"# (%s)--%s--(%s)--%s--(%s)\n",el[0],el[1],el[2],el[3],el[4]);
            if(equiv_t==j+1)fprintf(fp,"# %s\t%s\t%s\t%s\t%s\n",el[0],el[1],el[2],el[3],el[4]);
            
        }
        
        fprintf(fp,"\n");
    }
    fprintf(fp,"dihedral_style harmonic\n");
    fprintf(fp,"\n");
    /*
    for(j=0;j<D_type_counter;++j)
    {
        fprintf(fp,"# %s\n",dihedrals_multi_text[j]);
        fprintf(fp,"dihedral_coeff %d %lf %d %d\n",j+1,dihedrals_FF_multi[j][0],(int)dihedrals_FF_multi[j][1],(int)dihedrals_FF_multi[j][2]);
    }
    */
    for(j=0;j<new_D_types;++j)
    {
        if(dihedrals_out[j][2]!=-1)
        {
            
            fprintf(fp,"dihedral_coeff %d %lf %d %d\n",dihedrals_out[j][1],dihedrals_FF_new[j][0],(int)dihedrals_FF_new[j][1],(int)dihedrals_FF_new[j][2]);
            
            equiv_N=0;
            for(k=0;k<new_D_types;++k)
            {
                if(dihedrals_out[k][1]==dihedrals_out[j][1])
                {
                    equiv_N=equiv_N+1;
                }
            }
            fprintf(fp,"# %d equiv\n",equiv_N);

            for(k=0;k<new_D_types;++k)
            {
                if(dihedrals_out[k][1]==dihedrals_out[j][1])
                {
                    sscanf(new_D_registry[k],"%s\t%s\t%s\t%s\t%s\t%s\t%s",el[0],el[1],el[2],el[3],el[4],el[5],el[6]);
                    for(m=0;m<7;++m)
                    {
                        if(strcmp(el[m],"H_")==0)sprintf(el[m],"%s","H");
                        for(o=0;o<cmax_length;++o){if(el[m][o]=='\0')break;if(el[m][o]=='_')el[m][o]='.';}
                    }
                    //fprintf(fp,"# (%s)--%s--(%s)--%s--(%s)--%s--(%s)\n",el[0],el[1],el[2],el[3],el[4],el[5],el[6]);
                    fprintf(fp,"# %s\t%s\t%s\t%s\t%s\t%s\t%s\n",el[0],el[1],el[2],el[3],el[4],el[5],el[6]);
                }
            }
            fprintf(fp,"\n");
        }
    }
    if(general_impropers>0)
    {
        fprintf(fp,"improper_style umbrella\n");
        for(j=0;j<intermed_improper_types;++j)
        {
            fprintf(fp,"# %s -- %s -- %s -- %s\n",global_improper_type_intermed_array[j][0],global_improper_type_intermed_array[j][1],global_improper_type_intermed_array[j][2],global_improper_type_intermed_array[j][3]);
            fprintf(fp,"improper_coeff %d %lf %lf\n",j+1,impropers_FF[j][0],impropers_FF[j][1]);
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
    
    fprintf(fp,
            "thermo 1\n"
            "thermo_style custom step pe ebond eangle edihed eimp evdwl\n"
            "min_style cg\n"
            "minimize ${E_tol} ${F_tol} ${N_iter} ${N_eval}\n"
            "reset_timestep 0\n\n");
    
    fprintf(fp,"# write snapshot\n"
            "dump coords my_atoms custom 1 snapshot_%d.dat id mol type q x y z ix iy iz\n"
            "dump_modify coords sort id\n\n",l+1);
    
    fprintf(fp,"# write preview\n"
            "dump xyz_coords my_atoms xyz 1 animation.xyz\n"
            "dump_modify xyz_coords element ");
    for(j=0;j<general_atom_types;++j)fprintf(fp,"%c ",species_intermed_array[j][0]);
    fprintf(fp,"\n\n");
    
    fprintf(fp,
            "run 0\n\n"
            );
    
    fclose(fp);
    
    
    //
    
    free(type_core);
    /*
    free(D_type);free(reg1);free(reg2);
    for(j=0;j<D_type_counter;++j)free(dihedrals_FF_multi[j]);free(dihedrals_FF_multi);
    for(j=0;j<D_type_counter;++j)free(dihedrals_multi_text[j]);free(dihedrals_multi_text);
    */
    //
    
    for(j=0;j<new_D_types;++j)free(dihedrals_out[j]);free(dihedrals_out);
    for(j=0;j<new_D_types;++j)free(dihedrals_FF_new[j]);free(dihedrals_FF_new);
    for(j=0;j<new_D_types;++j)free(new_D_registry[j]);free(new_D_registry);
    
    
}
