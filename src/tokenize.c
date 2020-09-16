
#include"builder.h"

void tokenize(char *buffer,int *ignore_out,int *elements_out)
{
    ////////////////////////////////////////////////////////////////////////////////////
    // This function reads a line, tokenizes and merges consecutive white spaces and  //
    // tabs, rewrites each line in a tab-delimited format and returns the number      //
    // of tabs.                                                                       //
    ////////////////////////////////////////////////////////////////////////////////////
    // final revision: 30/04/2015 OGZ
    //
    char second_buffer[cmax_length];
    int j,k,l;
    int end,ignore=0;
    int elements;
    // loop over buffer and tokenize spaces and tabs
    end=-1;
    for(j=0;j<cmax_length;++j)
    {
        end=end+1;
        if(buffer[j]=='\0')break;
        if(buffer[j]==' ' || buffer[j]=='\t' || buffer[j]=='\n')buffer[j]='$';
    }
    // erase multiple dollars
    for(j=0;j<end-1;++j)    // end-1 because we check pairs j,j+1
    {
        if(buffer[j]=='$' && buffer[j+1]=='$')  // if you find two consecutive $s
        {
            l=-1;                               // index for writing to second_buffer
            for(k=0;k<=j;++k)                   // copy the left segment (including the first $) to memory (in second_buffer)
            {
                l=l+1;
                second_buffer[l]=buffer[k];
            }
            for(k=j+2;k<=end;++k)               // copy the remaining segment, without the second $
            {
                l=l+1;
                second_buffer[l]=buffer[k];
            }
            sprintf(buffer,"%s",second_buffer); // overwrite buffer
            j=j-1;                              // retrack to account for multiple $s
            end=end-1;                          // shrink dimension by one position
        }
    }
    // strip leading token
    if(buffer[0]=='$')
    {
        sprintf(second_buffer,"%s",buffer);
        for(j=0;j<cmax_length;++j)
        {
            if(second_buffer[j]=='\0')break;
            buffer[j]=second_buffer[j+1];
        }
    }
    // count the number of final tokens
    elements=0;
    for(j=0;j<cmax_length;++j)
    {
        if(buffer[j]=='\0')break;
        if(buffer[j]=='$')elements=elements+1;
    }
    // replace tokens
    /*
     for(j=0;j<cmax_length;++j)
     {
     if(buffer[j]=='\0')break;
     if(buffer[j]=='$')buffer[j]='\t';
     }
     */
    *ignore_out=ignore;
    *elements_out=elements;
    //
}
