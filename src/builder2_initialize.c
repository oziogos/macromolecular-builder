
#include"builder.h"

void builder2_initialize(char *current_folder,char *argv1,int *molecules,int *mol_types,double *xlo,double *xhi,double *ylo,double *yhi,double *zlo,double *zhi,
                         char ***core_file,int **type_array,int **chains_array,int **mol_id_array,int **mol_type_array,int **com_flag_array,
                         double **X,double **Y,double **Z,double **alpha,double **beta,double **gamma,
                         int *total_chains,char ***quad_matrix,
                         int *master_growth_registry_rows,int *master_growth_registry_cols,int *master_species_rows,int *master_topo_rows,
                         int **MGR_start,int **MGR_stop,int **MA_start,int **MA_stop,int **MB_start,int **MB_stop,
                         int *entries_UFF,
                         char ***species_UFF_array,
                         double **D_UFF_array,double **x_UFF_array,double **r0_UFF_array,double **chi_UFF_array,double **Zstar_UFF_array,
                         double **theta0_UFF_array,double **Vtor_UFF_array,double **Utor_UFF_array,double **zeta_UFF_array,
                         int *entries_DRE,
                         char ***species_DRE_array,
                         double **R0_DRE_array,double **theta0_DRE_array,double **Rvdw0_DRE_array,double **D0_DRE_array,
                         int *backbones,char ***backbone_id,double **bb_dx,double ***bb_entries,int *bb_lines_max,
                         int *tacticity,char ***tacticity_id,char ***tacticity_type,
                         int *backbone_type_counter, char ***backbone_type_array,
                         char *mpi_from_config, char *lmp_from_config,
                         double **bb_pdf_dx,double ***bb_pdf_entries,int *bb_pdf_lines_max)
{
    FILE *fp,*fp2;
    char file_path[cmax_length],buffer[cmax_length],word[cmax_length],file_path2[cmax_length],buffer2[cmax_length];
    int i,j,k;
    int atoms_chain,bonds_chain;
    int found,index;
    double z,zmap;
    int *bb_lines;
    char D1[sub_length];
    char D2[sub_length];
    char D3[sub_length];
    char D4[sub_length];
    char **init_backbone_type_array;
    //int backbone_type_counter;
    //char **backbone_type_array;
    
    //==========================================================================
    //
    // A typical master file looks like this:
    //
    //--------------------------------------------------------------------------
    //
    // number_of_molecules  4                                                   stored in int molecules
    //
    // supercell            0.0 40.0 0.0 40.0 0.0 40.0                          stored in double xlo,xhi,ylo,yhi,zlo,zhi
    //
    // number_of_types      2                                                   stored in int mol_types
    //
    // type                 1                                                   (II)
    // core                 CH4.mol2                                            (I)
    // chains               1                                                   (III)
    // C6.mol2              1   1   1
    //
    // type                 2                                                   (II)
    // core                 benzene.mol2                                        (I)
    // chains               2                                                   (III)
    // phytyl.mol2          1   1   1
    // phytyl.mol2          1   4   1
    //
    // molecules
    // 1   1   10.0   10.0   10.0   0.0   0.0   0.0   0
    // 2   2   10.0   20.0   10.0   0.0   0.0   0.0   0
    // 3   2   20.0   10.0   10.0   0.0   0.0   0.0   0
    // 4   1   20.0   20.0   10.0   0.0   0.0   0.0   0
    //
    // ^   ^    ^      ^      ^      ^     ^     ^    ^
    // ^   ^    ^      ^      ^      ^     ^     ^    ^
    // A   B    C      D      E      F     G     H    I
    //
    //--------------------------------------------------------------------------
    //
    // int mol_types is used to preallocate:
    // - char **core_file  (I)
    // - int *type_array   (II)
    // - int *chains_array (III)
    // - int *MGR_start
    // - int *MGR_stop
    // - int *MA_start
    // - int *MA_stop
    // - int *MB_start
    // - int *MB_stop
    //
    // int molecules is used to preallocate:
    // - int *mol_id_array      (A)
    // - int *mol_type_array    (B)
    // - int *com_flag_array    (I)
    // - double *X              (C)
    // - double *Y              (D)
    // - double *Z              (E)
    // - double *alpha          (F)
    // - double *beta           (G)
    // - double *gamma          (H)
    //
    // int total_chains counts the number of chains from all molecular types and
    // is used to preallocate char **quad_matrix.
    //
    // char **quad_matrix stores the following information:
    // - molecular type
    // - side chain mol2 filename
    // - graft id
    // - host id
    // - initial bond order
    // - chain atoms
    // - chain bonds
    //
    // int variables master_growth_registry_rows, master_growth_registry_cols,
    // master_species_rows, master_topo_rows are initialized to zero.
    //
    // 'uff.prm' data file is accessed and read. The number of lines is stored
    // in int entries_UFF. Parameter arrays are preallocated and populated.
    //
    // uff.prm file contains the parameters of the Universal Force Field:
    // Sources:
    // - A. K. Rappe et al. J. Am. Chem. Soc. 114, 10024-10035 (1992)
    // - crucial corrections and clarifications from
    //   http://towhee.sourceforge.net/forcefields/uff.html
    // This is an example of an appropriate uff.prm file:
    //
    //  atomic
    //  6
    //  species	r0	theta0	x	D	zeta	Z*	chi	V_tor	U_tor
    //  H_	0.354	180.0	2.886	0.044	12.0	0.712	4.528	0.0	0.0
    //  C_3	0.757	109.47	3.851	0.105	12.73	1.912	5.343	2.119	2.0
    //  C_R	0.729	120.0	3.851	0.105	12.73	1.912	5.343	0.0	2.0
    //  C_2	0.732	120.0	3.851	0.105	12.73	1.912	5.343	0.0	2.0
    //  O_3	0.658	104.51	3.500	0.060	14.085	2.300	8.741	0.018	2.0
    //  O_2	0.634	120.0	3.500	0.060	14.085	2.300	8.741	0.0	2.0
    //
    // The function also checks is there is a 'backbones' section. If so, int
    // backbones stores the number of entries and is used to preallocate char
    // **backbone_id. This array stores the backbone identifiers. int backbones
    // preallocates also double *bb_dx and is used for the preallocation of
    // double **bb_entries. The former array holds the storage interval per
    // backbone type whereas the latter array stores the data read from the
    // appropriate files. int bb_lines_max holds the maximum number of lines
    // from all backbone types and dictates the number of rows for bb_entries.
    //
    //==========================================================================
    
    // read current directory
    getcwd(current_folder,cmax_length);
    // initialize the number of backbone types
    *backbones=0;
    // initialize GLOBAL counters
    *master_growth_registry_rows=0;      // rows of **master_growth_registry
    *master_growth_registry_cols=0;      // columns of **master_growth_registry
    *master_species_rows=0;              // length of **master_species_chain_array
    *master_topo_rows=0;                 // length of *master_B1_chain_array,*master_B2_chain_array,**master_B_type_chain_array
    // read the master configuration file; open file
    sprintf(file_path,"%s/%s",current_folder,argv1);
    fp=fopen(file_path,"r");
    if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
    // echo
    printf("\n$- input file ------------------------------------------\n\n");
    while(fgets(buffer,cmax_length,fp)!=NULL)printf("%s",buffer);rewind(fp);printf("\n");
    // search for mpi
    while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"mpi")==0)break;}
    sscanf(buffer,"%s\t%s",word,mpi_from_config);
    rewind(fp);
    // search for lmp
    while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"lmp")==0)break;}
    sscanf(buffer,"%s\t%s",word,lmp_from_config);
    rewind(fp);
    // search for the number of molecules
    while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"number_of_molecules")==0)break;}
    sscanf(buffer,"%s\t%d",word,&(*molecules));    // read molecules
    rewind(fp);
    // search for the number of molecule types
    while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"number_of_types")==0)break;}
    sscanf(buffer,"%s\t%d",word,&(*mol_types));    // read mol_types
    rewind(fp);
    // search for the supercell
    while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"supercell")==0)break;}
    sscanf(buffer,"%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",word,&(*xlo),&(*xhi),&(*ylo),&(*yhi),&(*zlo),&(*zhi));   // read orthogonal supercell parameters
    rewind(fp);
    // preallocations
    // mol_type related
    *core_file=(char**)malloc(*mol_types*sizeof(char*));      // core molecule mol2 filenames
    for(i=0;i<*mol_types;++i)(*core_file)[i]=(char*)malloc(sub_length*sizeof(char));
    *type_array=(int*)malloc(*mol_types*sizeof(int));         // molecular type index
    *chains_array=(int*)malloc(*mol_types*sizeof(int));       // number of chains per molecular type
    // molecules related
    *mol_id_array=(int*)malloc(*molecules*sizeof(int));       // molecular id index
    *mol_type_array=(int*)malloc(*molecules*sizeof(int));     // molecular type (directly linked to type_array)
    *com_flag_array=(int*)malloc(*molecules*sizeof(int));     // flag to calculate CoMs: if equal to 0, calculate CoM using all atoms per molecule; if larger than 0, use the defined atoms
    *X=(double*)malloc(*molecules*sizeof(double));            // CoM x position array
    *Y=(double*)malloc(*molecules*sizeof(double));            // CoM y position array
    *Z=(double*)malloc(*molecules*sizeof(double));            // CoM z position array
    *alpha=(double*)malloc(*molecules*sizeof(double));        // alpha extrinsic Euler angle
    *beta=(double*)malloc(*molecules*sizeof(double));         // beta extrinsic Euler angle
    *gamma=(double*)malloc(*molecules*sizeof(double));        // gamma extrinsic Euler angle
    // populate *type_array,**core_file,*chains_array
    for(i=0;i<*mol_types;++i)
    {
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"type")==0)break;}
        sscanf(buffer,"%s\t%d",word,&(*type_array)[i]);
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"core")==0)break;}
        sscanf(buffer,"%s\t%s",word,(*core_file)[i]);
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"chains")==0)break;}
        sscanf(buffer,"%s\t%d",word,&(*chains_array)[i]);
    }
    rewind(fp);
    // read chain info
    *total_chains=0;                                             // this is the total number of chains from every molecular type
    for(i=0;i<*mol_types;++i)*total_chains=*total_chains+(*chains_array)[i];
    *quad_matrix=(char**)malloc(*total_chains*sizeof(char*));     // stores all side chain associated info: mol2 filename, graft id, host id, init bo and number of atoms and bonds
    for(i=0;i<*total_chains;++i)(*quad_matrix)[i]=(char*)malloc(cmax_length*sizeof(char));
    // populate quad_matrix
    k=-1;
    for(i=0;i<*mol_types;++i)
    {
        // locate 'chains'
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"chains")==0)break;}
        // read data
        for(j=0;j<(*chains_array)[i];++j)
        {
            fgets(buffer,cmax_length,fp);
            // open associated side chain mol2 file and read number of atoms and bonds
            sprintf(word,"%s","");
            sscanf(buffer,"%s",word);
            sprintf(file_path2,"%s/%s",current_folder,word);
            fp2=fopen(file_path2,"r");
            if(fp2==NULL){printf("Could not locate %s\n",file_path2);exit(-1);}
            while(fgets(buffer2,cmax_length,fp2)!=NULL)if(strcmp(buffer2,"@<TRIPOS>MOLECULE\n")==0)break;
            fgets(buffer2,cmax_length,fp2);fgets(buffer2,cmax_length,fp2);sscanf(buffer2,"%d\t%d",&atoms_chain,&bonds_chain);
            fclose(fp2);
            buffer[strcspn(buffer,"\n")]='\0';
            // advance index and write data to quad_matrix
            k=k+1;
            sprintf((*quad_matrix)[k],"%d\t%s\t%d\t%d\n",(*type_array)[i],buffer,atoms_chain,bonds_chain);
        }
    }
    rewind(fp);
    // read molecules
    while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"molecules")==0)break;}
    for(i=0;i<*molecules;++i)
    {
        fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d",&(*mol_id_array)[i],&(*mol_type_array)[i],&(*X)[i],&(*Y)[i],&(*Z)[i],&(*alpha)[i],&(*beta)[i],&(*gamma)[i],&(*com_flag_array)[i]);
    }
    rewind(fp);
    // search for the number of backbones
    found=0;
    while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"backbones")==0){found=1;break;}}
    // found==1 means that there is an active 'backbones' block in the config file!
    if(found==1)
    {
        sscanf(buffer,"%s\t%d",word,&(*backbones));                 // read the number of backbone types
        *backbone_id=(char**)malloc((*backbones)*sizeof(char*));    // backbone_id stores each backbone dihedral type identified: D1-D2-D3-D4
        //for(i=0;i<*backbones;++i)(*backbone_id)[i]=(char*)malloc(sub_length*sizeof(char));
        for(i=0;i<*backbones;++i)(*backbone_id)[i]=(char*)malloc(cmax_length*sizeof(char));
        bb_lines=(int*)malloc((*backbones)*sizeof(int));            // bb_lines holds the number of invCDF entries per backbone D type
        *bb_dx=(double*)malloc((*backbones)*sizeof(double));        // bb_dx holds the record (storage) interval per invCDF type
        // here we open the D type invCDF file and resolve the number of entries via fgets()
        /*
        for(i=0;i<*backbones;++i)
        {
            fgets(buffer,cmax_length,fp);
            //sscanf(buffer,"%s\t%s",(*backbone_id)[i],word);
            sscanf(buffer,"%s\t%s\t%s\t%s\t%s",D1,D2,D3,D4,word);
            sprintf(file_path2,"%s/%s",current_folder,word);
            fp2=fopen(file_path2,"r");
            if(fp2==NULL){printf("Could not locate %s\n",file_path2);exit(-1);}
            bb_lines[i]=0;
            while(fgets(buffer2,cmax_length,fp2)!=NULL)bb_lines[i]=bb_lines[i]+1;
            fclose(fp2);
        }
        */
        // use '#invCDF'
        for(i=0;i<*backbones;++i)
        {
            fgets(buffer,cmax_length,fp);
            //sscanf(buffer,"%s\t%s",(*backbone_id)[i],word);
            sscanf(buffer,"%s\t%s\t%s\t%s\t%s",D1,D2,D3,D4,word);
            sprintf(file_path2,"%s/%s",current_folder,word);
            fp2=fopen(file_path2,"r");
            if(fp2==NULL){printf("Could not locate %s\n",file_path2);exit(-1);}
            //bb_lines[i]=0;
            //while(fgets(buffer2,cmax_length,fp2)!=NULL)bb_lines[i]=bb_lines[i]+1;
            // search for '#invCDF'
            while(fgets(buffer2,cmax_length,fp2)!=NULL){
                sprintf(word,"%s","");sscanf(buffer2,"%s",word);if(strcmp(word,"#invCDF")==0){sscanf(buffer2,"%s\t%d",word,&bb_lines[i]);break;}}
            fclose(fp2);
        }
        rewind(fp);
        *bb_lines_max=bb_lines[0];
        for(i=0;i<*backbones;++i)if(bb_lines[i]>*bb_lines_max)*bb_lines_max=bb_lines[i];
        *bb_entries=(double**)malloc((*bb_lines_max)*sizeof(double*));
        for(i=0;i<*bb_lines_max;++i)(*bb_entries)[i]=(double*)malloc((2*(*backbones))*sizeof(double));
        for(i=0;i<*bb_lines_max;++i)for(j=0;j<2*(*backbones);++j)(*bb_entries)[i][j]=0.0;
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"backbones")==0){break;}}
        for(i=0;i<*backbones;++i)
        {
            fgets(buffer,cmax_length,fp);
            //sscanf(buffer,"%s\t%s",(*backbone_id)[i],word);
            sscanf(buffer,"%s\t%s\t%s\t%s\t%s",D1,D2,D3,D4,word);
            sprintf((*backbone_id)[i],"%s\t%s\t%s\t%s",D1,D2,D3,D4);
            sprintf(file_path2,"%s/%s",current_folder,word);
            fp2=fopen(file_path2,"r");
            if(fp2==NULL){printf("Could not locate %s\n",file_path2);exit(-1);}
            /*
            j=-1;
            while(fgets(buffer2,cmax_length,fp2)!=NULL)
            {
                j=j+1;
                sscanf(buffer2,"%lf\t%lf",&(*bb_entries)[j][0+i*2],&(*bb_entries)[j][1+i*2]);
            }
            */
            while(fgets(buffer2,cmax_length,fp2)!=NULL){sprintf(word,"%s","");sscanf(buffer2,"%s",word);if(strcmp(word,"#invCDF")==0)break;}
            for(j=0;j<bb_lines[i];++j)
            {
                fgets(buffer2,cmax_length,fp2);
                sscanf(buffer2,"%lf\t%lf",&(*bb_entries)[j][0+i*2],&(*bb_entries)[j][1+i*2]);
            }
            fclose(fp2);
        }
        rewind(fp);
        
        // read #PDF; reuse bb_lines
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"backbones")==0){break;}}
        *bb_pdf_dx=(double*)malloc((*backbones)*sizeof(double));
        for(i=0;i<*backbones;++i)
        {
            fgets(buffer,cmax_length,fp);
            sscanf(buffer,"%s\t%s\t%s\t%s\t%s",D1,D2,D3,D4,word);
            sprintf(file_path2,"%s/%s",current_folder,word);
            fp2=fopen(file_path2,"r");
            if(fp2==NULL){printf("Could not locate %s\n",file_path2);exit(-1);}
            // search for '#PDF'
            while(fgets(buffer2,cmax_length,fp2)!=NULL){
                sprintf(word,"%s","");sscanf(buffer2,"%s",word);if(strcmp(word,"#PDF")==0){sscanf(buffer2,"%s\t%d",word,&bb_lines[i]);break;}}
            fclose(fp2);
        }
        rewind(fp);
        *bb_pdf_lines_max=bb_lines[0];
        for(i=0;i<*backbones;++i)if(bb_lines[i]>*bb_pdf_lines_max)*bb_pdf_lines_max=bb_lines[i];
        *bb_pdf_entries=(double**)malloc((*bb_pdf_lines_max)*sizeof(double*));
        for(i=0;i<*bb_pdf_lines_max;++i)(*bb_pdf_entries)[i]=(double*)malloc((2*(*backbones))*sizeof(double));
        for(i=0;i<*bb_pdf_lines_max;++i)for(j=0;j<2*(*backbones);++j)(*bb_pdf_entries)[i][j]=0.0;
        while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"backbones")==0){break;}}
        for(i=0;i<*backbones;++i)
        {
            fgets(buffer,cmax_length,fp);
            sscanf(buffer,"%s\t%s\t%s\t%s\t%s",D1,D2,D3,D4,word);
            sprintf(file_path2,"%s/%s",current_folder,word);
            fp2=fopen(file_path2,"r");
            if(fp2==NULL){printf("Could not locate %s\n",file_path2);exit(-1);}
            while(fgets(buffer2,cmax_length,fp2)!=NULL){sprintf(word,"%s","");sscanf(buffer2,"%s",word);if(strcmp(word,"#PDF")==0)break;}
            for(j=0;j<bb_lines[i];++j)
            {
                fgets(buffer2,cmax_length,fp2);
                sscanf(buffer2,"%lf\t%lf",&(*bb_pdf_entries)[j][0+i*2],&(*bb_pdf_entries)[j][1+i*2]);
            }
            fclose(fp2);
        }
    }
    /*
    for(i=0;i<*bb_pdf_lines_max;++i)
    {
        for(j=0;j<*backbones;++j)printf("%lf\t%lf\t",(*bb_pdf_entries)[i][0+j*2],(*bb_pdf_entries)[i][1+j*2]);
        printf("\n");
    }
    */
    // tacticity
    rewind(fp);
    *tacticity=0;
    found=0;
    while(fgets(buffer,cmax_length,fp)!=NULL){sprintf(word,"%s","");sscanf(buffer,"%s",word);if(strcmp(word,"tacticity")==0){found=1;break;}}
    if(found==1)
    {
        sscanf(buffer,"%s\t%d",word,&(*tacticity));    // read
        *tacticity_id=(char**)malloc((*tacticity)*sizeof(char*));
        for(i=0;i<*tacticity;++i)(*tacticity_id)[i]=(char*)malloc(sub_length*sizeof(char));
        *tacticity_type=(char**)malloc((*tacticity)*sizeof(char*));
        for(i=0;i<*tacticity;++i)(*tacticity_type)[i]=(char*)malloc(sub_length*sizeof(char));
        for(i=0;i<*tacticity;++i){fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%s",(*tacticity_id)[i],(*tacticity_type)[i]);}
        
    }
    
    fclose(fp);
    for(i=0;i<*backbones;++i)
    {
        (*bb_dx)[i]=(*bb_entries)[1][0+i*2]-(*bb_entries)[0][0+i*2];
        (*bb_pdf_dx)[i]=(*bb_pdf_entries)[1][0+i*2]-(*bb_pdf_entries)[0][0+i*2];
        sprintf(file_path,"%s/backbone_%d_stat.dat",current_folder,i+1);
        fp=fopen(file_path,"w+");
        for(k=0;k<1000;++k)
        {
            z=(double)rand()/RAND_MAX;
            index=(int)floor(z/(*bb_dx)[i]);
            zmap=(*bb_entries)[index][1+i*2]+((*bb_entries)[index+1][1+i*2]-(*bb_entries)[index][1+i*2])*(z-(*bb_entries)[index][0+i*2])/(*bb_dx)[i];
            //printf("%lf\n",zmap);
            fprintf(fp,"%lf\n",zmap);
        }
        fclose(fp);
    }
    // read uff.prm
    // open the file
    sprintf(file_path,"%s/uff.prm",current_folder);
    fp=fopen(file_path,"r");
    // file check
    if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
    // locate the atomic section
    while(fgets(buffer,cmax_length,fp)!=NULL)
        if(strcmp(buffer,"atomic\n")==0)
            break;
    // read number of entries
    fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&(*entries_UFF));
    // preallocate arrays to store UFF parameters
    *species_UFF_array=(char**)malloc(*entries_UFF*sizeof(char*));
    for(i=0;i<*entries_UFF;++i)(*species_UFF_array)[i]=(char*)malloc(sub_length*sizeof(char));
    *r0_UFF_array=(double*)malloc(*entries_UFF*sizeof(double));
    *theta0_UFF_array=(double*)malloc(*entries_UFF*sizeof(double));
    *x_UFF_array=(double*)malloc(*entries_UFF*sizeof(double));
    *D_UFF_array=(double*)malloc(*entries_UFF*sizeof(double));
    *zeta_UFF_array=(double*)malloc(*entries_UFF*sizeof(double));
    *Zstar_UFF_array=(double*)malloc(*entries_UFF*sizeof(double));
    *chi_UFF_array=(double*)malloc(*entries_UFF*sizeof(double));
    *Vtor_UFF_array=(double*)malloc(*entries_UFF*sizeof(double));
    *Utor_UFF_array=(double*)malloc(*entries_UFF*sizeof(double));
    // read data
    fgets(buffer,cmax_length,fp);
    for(i=0;i<*entries_UFF;++i)
    {
        fgets(buffer,cmax_length,fp);
        sscanf(buffer,"%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",(*species_UFF_array)[i],&(*r0_UFF_array)[i],&(*theta0_UFF_array)[i],&(*x_UFF_array)[i],&(*D_UFF_array)[i],&(*zeta_UFF_array)[i],&(*Zstar_UFF_array)[i],&(*chi_UFF_array)[i],&(*Vtor_UFF_array)[i],&(*Utor_UFF_array)[i]);
    }
    // close uff.prm
    fclose(fp);
    // preallocate arrays to store the upper and lower data boundaries
    *MGR_start=(int*)malloc(*mol_types*sizeof(int));  // GR: Growth Registry
    *MGR_stop=(int*)malloc(*mol_types*sizeof(int));
    *MA_start=(int*)malloc(*mol_types*sizeof(int));   // A: Atomic data (species)
    *MA_stop=(int*)malloc(*mol_types*sizeof(int));
    *MB_start=(int*)malloc(*mol_types*sizeof(int));   // B: Bond data (B1,B2,Btype)
    *MB_stop=(int*)malloc(*mol_types*sizeof(int));
    //
    if(*backbones>0){
        free(bb_lines);
        //free(bb_dx);
        //for(i=0;i<bb_lines_max;++i)free(bb_entries[i]);free(bb_entries);
    }
    
    if(*tacticity!=0)
    {
        sprintf(file_path,"%s/tacticity_history.dat",current_folder);
        fp=fopen(file_path,"w+");
        fclose(fp);
    }
    /*
    if(*backbones>0)
    {
        sprintf(file_path,"%s/Cn_history.dat",current_folder);
        fp=fopen(file_path,"w+");
        fclose(fp);
        
        sprintf(file_path,"%s/Cn_vs_build.dat",current_folder);
        fp=fopen(file_path,"w+");
        fclose(fp);
    }
     */
    //
    /*
    for(i=0;i<*backbones;++i)
    {
        printf("%s\n",(*backbone_id)[i]);
    }
    */
    if(*backbones>0)
    {
        init_backbone_type_array=(char**)malloc(4*(*backbones)*sizeof(char*));
        for(i=0;i<4*(*backbones);++i)init_backbone_type_array[i]=(char*)malloc(sub_length*sizeof(char));
        
        j=-1;
        for(i=0;i<*backbones;++i)
        {
            sscanf((*backbone_id)[i],"%s\t%s\t%s\t%s",D1,D2,D3,D4);
            j=j+1;sprintf(init_backbone_type_array[j],"%s",D1);
            j=j+1;sprintf(init_backbone_type_array[j],"%s",D2);
            j=j+1;sprintf(init_backbone_type_array[j],"%s",D3);
            j=j+1;sprintf(init_backbone_type_array[j],"%s",D4);
        }
        
        //for(i=0;i<4*(*backbones);++i)printf("%s\n",init_backbone_type_array[i]);
        
        for(i=0;i<4*(*backbones)-1;++i)
        {
            for(j=i+1;j<4*(*backbones);++j)
            {
                if(strcmp(init_backbone_type_array[i],init_backbone_type_array[j])==0)sprintf(init_backbone_type_array[j],"%s","___");
            }
        }
        
        //printf("-----\n");
        //for(i=0;i<4*(*backbones);++i)printf("%s\n",init_backbone_type_array[i]);
        
        k=0;
        for(i=0;i<4*(*backbones);++i)
            if(strcmp(init_backbone_type_array[i],"___")==0)k=k+1;
        *backbone_type_counter=4*(*backbones)-k;
        //printf("%d\n",*backbone_type_counter);
        
        *backbone_type_array=(char**)malloc(*backbone_type_counter*sizeof(char*));
        for(i=0;i<*backbone_type_counter;++i)(*backbone_type_array)[i]=(char*)malloc(sub_length*sizeof(char));
        
        j=-1;
        for(i=0;i<4*(*backbones);++i)
        {
            if(!(strcmp(init_backbone_type_array[i],"___")==0))
            {
                j=j+1;
                sprintf((*backbone_type_array)[j],"%s",init_backbone_type_array[i]);
            }
        }
        
        //for(i=0;i<*backbone_type_counter;++i)printf("%s\n",(*backbone_type_array)[i]);
        
        for(i=0;i<4*(*backbones);++i)free(init_backbone_type_array[i]);free(init_backbone_type_array);
    }
    
    
    // read dreiding.prm
    // open the file
    sprintf(file_path,"%s/dreiding.prm",current_folder);
    fp=fopen(file_path,"r");
    // file check
    if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
    // locate the atomic section
    while(fgets(buffer,cmax_length,fp)!=NULL)
        if(strcmp(buffer,"atomic\n")==0)
            break;
    // read number of entries
    fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&(*entries_DRE));
    // preallocate arrays to store DREIDING parameters
    *species_DRE_array=(char**)malloc(*entries_DRE*sizeof(char*));
    for(i=0;i<*entries_DRE;++i)(*species_DRE_array)[i]=(char*)malloc(sub_length*sizeof(char));
    *R0_DRE_array=(double*)malloc(*entries_DRE*sizeof(double));
    *theta0_DRE_array=(double*)malloc(*entries_DRE*sizeof(double));
    *Rvdw0_DRE_array=(double*)malloc(*entries_DRE*sizeof(double));
    *D0_DRE_array=(double*)malloc(*entries_DRE*sizeof(double));
    // read data
    fgets(buffer,cmax_length,fp);
    for(i=0;i<*entries_DRE;++i)
    {
        fgets(buffer,cmax_length,fp);
        sscanf(buffer,"%s\t%lf\t%lf\t%lf\t%lf",(*species_DRE_array)[i],&(*R0_DRE_array)[i],&(*theta0_DRE_array)[i],&(*Rvdw0_DRE_array)[i],&(*D0_DRE_array)[i]);
    }
    // close uff.prm
    fclose(fp);
    
    
    //getchar();
}
