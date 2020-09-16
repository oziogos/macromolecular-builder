
#include"builder.h"

void linked_list_cell_initialize(int cutoff,double xlo,double xhi,double ylo,double yhi,double zlo,double zhi,int general_atoms,
                                 int *Mx,int *My,int *Mz,int *M3,double *lx,double *ly,double *lz,double *lcx,double *lcy,double *lcz,int ***neighbors,int **head)
{
    int i,j,k,l,M2;
    int ghostsx,ghostsy,**cells,row,column,ii,jj;
    
    // calculate supercell edges
    *lx=xhi-xlo;*ly=yhi-ylo;*lz=zhi-zlo;
    
    // resolve Ms
    *Mx=floor(*lx/cutoff);
    *My=floor(*ly/cutoff);
    *Mz=floor(*lz/cutoff);
    if(*Mx<3){*Mx=3;}
    if(*My<3){*My=3;}
    if(*Mz<3){*Mz=3;}
    *M3=(*Mx)*(*My)*(*Mz);
    // cell edges
    *lcx=(*lx)/(*Mx);*lcy=(*ly)/(*My);*lcz=(*lz)/(*Mz);
    // preallocate and populate cell neighbors matrix
    *neighbors=(int**)malloc(26*sizeof(int*));
    for(i=0;i<26;++i){(*neighbors)[i]=(int*)malloc((*M3)*sizeof(int));}
    
    
    
    // cells
    M2=(*Mx)*(*My);
    
    // ghost cells
    ghostsx=*Mx+2;
    ghostsy=*My+2;
    cells=(int**)malloc(ghostsy*sizeof(int*));
    for (i=0;i<ghostsy;++i)
    {
        cells[i]=(int*)malloc(ghostsx*sizeof(int));
    }
    
    for (i=0;i<ghostsy;++i){for (j=0;j<ghostsx;++j){cells[i][j]=0;}}
    
    // index rule: subtract 1 from static values and keep running indices the same! reminder: place counters at the bottom!!
    l=0;
    for (i=*My-1;i>=0;--i)
    {
        for (j=0;j<=*Mx-1;++j)
        {
            l=l+1;
            cells[i+1][j+1]=l;
            
        }
    }
    
    for (j=1;j<=ghostsx-2;++j)
    {
        cells[0][j]=cells[ghostsy-2][j];
        cells[ghostsy-1][j]=cells[1][j];
    }
    
    for (i=0;i<=ghostsy-1;++i)
    {
        cells[i][ghostsx-1]=cells[i][1];
        cells[i][0]=cells[i][ghostsx-2];
    }
    
    for (i=0;i<26;++i){for (j=0;j<*M3;++j){(*neighbors)[i][j]=0;}}
    
    l=0;
    for (i=*My-1;i>=0;--i)
    {
        for (j=0;j<=*Mx-1;++j)
        {
            row=i+2;column=j+2;
            k=0;
            for (jj=0;jj<=2;++jj)
            {
                for (ii=0;ii<=2;++ii)
                {
                    if (ii!=1||jj!=1)
                    {
                        (*neighbors)[k][l]=cells[ii-2+row][jj-2+column];
                        k=k+1;
                    }
                }
            }
            l=l+1;
        }
    }
    
    for (j=M2;j<=*M3-1;++j)
    {
        for (i=0;i<=7;++i)
        {
            (*neighbors)[i][j]=(*neighbors)[i][j-M2]+M2;
        }
    }
    
    for (j=0;j<=*M3-1;++j)
    {
        for (i=0;i<=7;++i)
        {
            if ((*neighbors)[i][j]+M2<=*M3)
            {
                (*neighbors)[i+8][j]=(*neighbors)[i][j]+M2;
            }
            else
            {
                (*neighbors)[i+8][j]=(*neighbors)[i][j]-(*Mz-1)*M2;
            }
            if ((*neighbors)[i][j]-M2>0)
            {
                (*neighbors)[i+16][j]=(*neighbors)[i][j]-M2;
            }
            else
            {
                (*neighbors)[i+16][j]=(*neighbors)[i][j]+(*Mz-1)*M2;
            }
        }
        if (j+1+M2<=*M3)
        {
            (*neighbors)[24][j]=j+1+M2;
        }
        else
        {
            (*neighbors)[24][j]=j+1-(*Mz-1)*M2;
        }
        if (j+1-M2>0)
        {
            (*neighbors)[25][j]=j+1-M2;
        }
        else
        {
            (*neighbors)[25][j]=j+1+(*Mz-1)*M2;
        }
    }
    
    for (i=0;i<26;++i)
    {
        for (j=0;j<*M3;++j)
        {
            (*neighbors)[i][j]=(*neighbors)[i][j]-1;
        }
        
    }
    
    for (i=0;i<ghostsy;++i){free(cells[i]);}
    free(cells);
    
    // preallocate head and list
    *head=(int*)malloc((*M3)*sizeof(int));// head
}
