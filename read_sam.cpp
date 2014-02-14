#include "mex.h" 
#include "matrix.h"
#include <string.h>
#include <vector>
#include <fstream>


mxArray * getMexArrayInt(const std::vector<int>& v){
    mxArray * mx = mxCreateDoubleMatrix(1,v.size(), mxREAL);
    std::copy(v.begin(), v.end(), mxGetPr(mx));
    return mx;
}

void mexFunction(int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    char line[1000];
    int column;
    char *input_buf;
    input_buf = mxArrayToString(prhs[0]);
    std::vector<int> start; 
    std::vector<std::string> chrom; 
    FILE *file;
    mxArray *chrom_cell;
    
    file = fopen(input_buf, "r");
    if (file) {
        while (fgets(line, 1000, file) != NULL) {
            if (line[0] != '@'){
                char * pch;
                pch = strtok (line," \t");
                column = 0;
                while (pch != NULL)
                {            
                    pch = strtok (NULL, " \t");
                    if (column == 1 && pch == "*"){
                        break;
                    }
                    else if (column == 1) {
                        chrom.push_back (pch);
                    }
                    else if (column == 2) {
                        start.push_back(atoi(pch));
                    }
                    column = column ++;
                }
            }
        }
       
        fclose(file);
    }
    
    chrom_cell = mxCreateCellMatrix(start.size(), 1);

    for (int i=0; i<chrom.size(); i++){
        mxSetCell(chrom_cell, i, mxCreateString(chrom[i].c_str()));
    }
      
    plhs[0] = getMexArrayInt(start);
    plhs[1]= chrom_cell;   
    return;
}

