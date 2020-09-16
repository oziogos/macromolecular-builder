
#include"builder.h"

void lb(FILE *fp,int start,int *branches_out,int *depth_out,int ***lb_output_matrix)
{
    char buffer[cmax_length];
    char **species;
    int i,j,k,l,atoms,bonds,bonds_f,atoms_f,max;
    int *id;
    int *b_id,*b1,*b2;
    int *b1_f,*b2_f;
    double *x,*y,*z;
    // neighbors
    int max_1_2;
    int *one_two_counter_array;
    int **one_two;
    int **nHneigh;
    //
    int ii,jj;
    int branches,start_neigh,eligible;
    int ecounter,master;
    int *aux;
    int **matrix1,**matrix2;
    
    char **label;
    
    fgets(buffer,cmax_length,fp);
    fgets(buffer,cmax_length,fp);
    fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d",&atoms,&bonds);
    
    id=(int*)malloc(atoms*sizeof(int));     // preallocations
    species=(char**)malloc(atoms*sizeof(char*));
    for(i=0;i<atoms;++i)species[i]=(char*)malloc(10*sizeof(char));
    x=(double*)malloc(atoms*sizeof(double));
    y=(double*)malloc(atoms*sizeof(double));
    z=(double*)malloc(atoms*sizeof(double));
    
    label=(char**)malloc(atoms*sizeof(char*));
    for(i=0;i<atoms;++i)label[i]=(char*)malloc(10*sizeof(char));
    
    while(fgets(buffer,cmax_length,fp)!=NULL)
    {
        if(strcmp(buffer,"@<TRIPOS>ATOM\n")==0)
            break;
    }
    for(i=0;i<atoms;++i)
    {
        fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%s\t%lf\t%lf\t%lf\t%s",&id[i],species[i],&x[i],&y[i],&z[i],label[i]);
    }
    
    b_id=(int*)malloc(bonds*sizeof(int));   // preallocations
    b1=(int*)malloc(bonds*sizeof(int));
    b2=(int*)malloc(bonds*sizeof(int));
    
    while(fgets(buffer,cmax_length,fp)!=NULL)
    {
        if(strcmp(buffer,"@<TRIPOS>BOND\n")==0)
            break;
    }
    for(i=0;i<bonds;++i)
    {
        fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d\t%d",&b_id[i],&b1[i],&b2[i]);
    }
    
    // eliminate X-H bonds; new bonds arrays: b1_f, b2_f
    
    bonds_f=bonds;
    for(i=0;i<bonds;++i)
    {
        if(strcmp(species[b1[i]-1],"H")==0 || strcmp(species[b2[i]-1],"H")==0)
        {
            bonds_f=bonds_f-1;
        }
    }
    b1_f=(int*)malloc(bonds_f*sizeof(int));
    b2_f=(int*)malloc(bonds_f*sizeof(int));
    
    j=-1;
    for(i=0;i<bonds;++i)
    {
        if((strcmp(species[b1[i]-1],"H")==0 || strcmp(species[b2[i]-1],"H")==0)==0)
        {
            j=j+1;
            b1_f[j]=b1[i];
            b2_f[j]=b2[i];
        }
    }
    
    //------------------------- - - - - - - --- - -- - - - - - -- -- --- ----- -
    // TOPOLOGY
    
    //--------------------------------------------------------------------------
    // 1-2 neighbors
    //--------------------------------------------------------------------------
    
    one_two_counter_array=(int*)malloc(atoms*sizeof(int));
    for(i=0;i<atoms;++i)one_two_counter_array[i]=0;
    for(i=0;i<bonds;++i)
    {
        one_two_counter_array[b1[i]-1]=one_two_counter_array[b1[i]-1]+1;
        one_two_counter_array[b2[i]-1]=one_two_counter_array[b2[i]-1]+1;
    }
    
    max_1_2=1;
    for(i=0;i<atoms;++i)
        if(one_two_counter_array[i]>max_1_2)max_1_2=one_two_counter_array[i];
    
    one_two=(int**)malloc(atoms*sizeof(int*));
    for(i=0;i<atoms;++i)one_two[i]=(int*)malloc(max_1_2*sizeof(int));
    for(i=0;i<atoms;++i)
        for(j=0;j<max_1_2;++j)
            one_two[i][j]=0;
    
    for(i=0;i<bonds;++i)
    {
        j=0;
        while(one_two[b1[i]-1][j]!=0)
            j=j+1;
        one_two[b1[i]-1][j]=b2[i];
        j=0;
        while(one_two[b2[i]-1][j]!=0)
            j=j+1;
        one_two[b2[i]-1][j]=b1[i];
    }
    
    // eliminate hydrogens
    for(i=0;i<atoms;++i)
    {
        for(j=0;j<max_1_2;++j)
            if(one_two[i][j]!=0)
                
                if(strcmp(species[one_two[i][j]-1],"H")==0)
                    one_two[i][j]=0;
    }
    
    // count non-hydrogen 1-2 neighbors for non-hydrogen species
    
    atoms_f=0;
    for(i=0;i<atoms;++i)
        if(strcmp(species[i],"H")!=0)
            atoms_f=atoms_f+1;
    nHneigh=(int**)malloc(atoms_f*sizeof(int*));
    for(i=0;i<atoms_f;++i)nHneigh[i]=(int*)malloc(2*sizeof(int));
    k=-1;
    for(i=0;i<atoms;++i)
    {
        if(strcmp(species[i],"H")!=0)
        {
            k=k+1;
            max=0;
            for(j=0;j<max_1_2;++j)
                if(one_two[i][j]!=0)
                    max=max+1;
            nHneigh[k][0]=i+1;
            nHneigh[k][1]=max;
        }
    }
    
    //
    
    // define the starting atom
    //start=1;
    
    // resolve number of branches
    branches=0;
    for(i=0;i<atoms_f;++i)
    {
        if(nHneigh[i][1]>1)
        {
            branches=branches+nHneigh[i][1]-2;
        }
    }
    for(i=0;i<atoms_f;++i)
    {
        if(nHneigh[i][0]==start)
        {
            start_neigh=nHneigh[i][1];
            break;
        }
    }
    if(start_neigh>1)
        branches=branches+2;
    else
        branches=branches+1;
    
    //printf("\nbranches = %d\n\n",branches);
    *branches_out=branches;
    
    //
    
    aux=(int*)malloc(atoms_f*sizeof(int));
    matrix1=(int**)malloc(branches*sizeof(int*));
    for(i=0;i<branches;++i)matrix1[i]=(int*)malloc(atoms_f*sizeof(int));
    for(i=0;i<branches;++i)
        for(j=0;j<atoms_f;++j)
            matrix1[i][j]=0;
    matrix2=(int**)malloc(branches*sizeof(int*));
    for(i=0;i<branches;++i)matrix2[i]=(int*)malloc(atoms_f*sizeof(int));
    for(i=0;i<branches;++i)
        for(j=0;j<atoms_f;++j)
            matrix2[i][j]=0;
    
    // initialize
    for(i=0;i<atoms_f;++i)
        aux[i]=0;
    // find the neighbors of the starting atom
    max=0;
    for(i=0;i<max_1_2;++i)
        if(one_two[start-1][i]!=0)
            max=max+1;
    j=-1;
    for(i=0;i<max_1_2;++i)
        if(one_two[start-1][i]!=0)
        {
            j=j+1;
            aux[j]=one_two[start-1][i];
        }
    for(i=0;i<max;++i)
    {
        matrix1[i][0]=start;
        matrix1[i][1]=aux[i];
    }
    
    //
    
    // scan lines that don't start with zero and find the outermost non-zero entry
    master=0;
    for(i=0;i<branches;++i)
    {
        if(matrix1[i][0]!=0)
        {
            for(j=0;j<atoms_f;++j)
            {
                if(matrix1[i][j]==0)
                {
                    k=j-1;
                    break;
                }
            }
            
            ecounter=0;
            for(j=0;j<max_1_2;++j)
            {
                if(one_two[matrix1[i][k]-1][j]!=0)
                {
                    eligible=1;
                    for(l=0;l<atoms_f;++l)
                    {
                        if(matrix1[i][l]==one_two[matrix1[i][k]-1][j])
                        {
                            eligible=0;
                            break;
                        }
                    }
                    if(eligible==1)
                    {
                        ecounter=ecounter+1;
                        
                        //
                        for(ii=0;ii<branches;++ii)
                        {
                            if(matrix2[ii][0]==0)
                            {
                                for(jj=0;jj<atoms_f;++jj)
                                {
                                    matrix2[ii][jj]=matrix1[i][jj];
                                }
                                for(jj=0;jj<atoms_f;++jj)
                                {
                                    if(matrix2[ii][jj]==0)
                                    {
                                        matrix2[ii][jj]=one_two[matrix1[i][k]-1][j];
                                        break;
                                    }
                                }
                                break;
                            }
                        }
                        //
                    }
                }
            }
            
            if(ecounter==0)
            {
                for(ii=0;ii<branches;++ii)
                {
                    if(matrix2[ii][0]==0)
                    {
                        for(jj=0;jj<atoms_f;++jj)
                        {
                            matrix2[ii][jj]=matrix1[i][jj];
                        }
                        
                        break;
                    }
                }
            }
            
            master=master+ecounter;
            
        }
        
    }
    
    while(master>0)
        //for(kk=0;kk<10;++kk)
    {
        
        // copy
        for(i=0;i<branches;++i)
            for(j=0;j<atoms_f;++j)
                matrix1[i][j]=matrix2[i][j];
        
        // erase
        for(i=0;i<branches;++i)
            for(j=0;j<atoms_f;++j)
                matrix2[i][j]=0;
        
        // repeat
        
        // scan lines that don't start with zero and find the outermost non-zero entry
        master=0;
        for(i=0;i<branches;++i)
        {
            if(matrix1[i][0]!=0)
            {
                for(j=0;j<atoms_f;++j)
                {
                    if(matrix1[i][j]==0)
                    {
                        k=j-1;
                        break;
                    }
                }
                
                ecounter=0;
                for(j=0;j<max_1_2;++j)
                {
                    if(one_two[matrix1[i][k]-1][j]!=0)
                    {
                        eligible=1;
                        for(l=0;l<atoms_f;++l)
                        {
                            if(matrix1[i][l]==one_two[matrix1[i][k]-1][j])
                            {
                                eligible=0;
                                break;
                            }
                        }
                        if(eligible==1)
                        {
                            ecounter=ecounter+1;
                            
                            //
                            for(ii=0;ii<branches;++ii)
                            {
                                if(matrix2[ii][0]==0)
                                {
                                    for(jj=0;jj<atoms_f;++jj)
                                    {
                                        matrix2[ii][jj]=matrix1[i][jj];
                                    }
                                    for(jj=0;jj<atoms_f;++jj)
                                    {
                                        if(matrix2[ii][jj]==0)
                                        {
                                            matrix2[ii][jj]=one_two[matrix1[i][k]-1][j];
                                            break;
                                        }
                                    }
                                    break;
                                }
                            }
                            //
                        }
                    }
                }
                
                if(ecounter==0)
                {
                    for(ii=0;ii<branches;++ii)
                    {
                        if(matrix2[ii][0]==0)
                        {
                            for(jj=0;jj<atoms_f;++jj)
                            {
                                matrix2[ii][jj]=matrix1[i][jj];
                            }
                            
                            break;
                        }
                    }
                }
                
                master=master+ecounter;
            }
            
        }
        
    }
    
    //
    
    *lb_output_matrix=(int**)malloc(branches*sizeof(int*));
    for(i=0;i<branches;++i)(*lb_output_matrix)[i]=(int*)malloc(atoms_f*sizeof(int));
    
    // out
    for(i=0;i<branches;++i)
        for(j=0;j<atoms_f;++j)
            (*lb_output_matrix)[i][j]=matrix2[i][j];
    
    max=1;
    for(i=0;i<branches;++i)
    {
        k=0;
        for(j=0;j<atoms_f;++j)
        {
            if(matrix2[i][j]!=0)
                k=k+1;
        }
        if(k>max)max=k;
    }
    //printf("\ndepth = %d\n\n",max);
    
    *depth_out=max;
    
    //
    
    free(id);for(i=0;i<atoms;++i)free(species[i]);free(species);
    free(b_id);free(b1);free(b2);
    free(b1_f);free(b2_f);
    
    free(one_two_counter_array);
    for(i=0;i<atoms;++i)free(one_two[i]);free(one_two);
    for(i=0;i<atoms_f;++i)free(nHneigh[i]);free(nHneigh);
    
    free(aux);
    for(i=0;i<branches;++i)free(matrix1[i]);free(matrix1);
    for(i=0;i<branches;++i)free(matrix2[i]);free(matrix2);
    
    free(x);free(y);free(z);
    
    for(i=0;i<atoms;++i)free(label[i]);free(label);
}
