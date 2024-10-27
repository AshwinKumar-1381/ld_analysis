/*
 General utility functions to work with files
 
 Author         : Ashwin Kumar M
 Date created   : 26.10.24
 Last modified  : 26.10.24 
*/

#include "file_utils.h"

using namespace fileUtils;

// Make file path based on the parameters: 
// nr(file no.), type1(lmp/cpp), type2(eq/ss), extras(log/thermo/traj)

char* fileUtils::makeFilePath(int nr, char *type1, char *type2, string extras){

    char *fname = (char*) malloc(500*sizeof(char));
    const char *extra = extras.c_str();
    
    string str(filesystem::current_path());
    str.erase(str.begin() + str.find("/ld_analysis/traj"), str.end());
    sprintf(fname, "%s/LD/%s/Data%d", str.c_str(), type1, nr);
    
    if(strcmp(type1, "lmp") == 0){
        if(strcmp(extra, "traj") == 0) sprintf(fname, "%s/%s%s%d.xyz", fname, extra, type2, nr);
        else if(strcmp(extra, "traj_orig") == 0) sprintf(fname, "%s/traj%s%d_orig.xyz", fname, type2, nr);
        else sprintf(fname, "%s/%s.dat", fname, extra);
    }
    
    else{
        if(strcmp(extra, "traj") == 0) sprintf(fname, "%s/relaxa_%s/traj_%s_%d.xyz", fname, type2, type2, nr);
        else if(strcmp(extra, "traj_orig") == 0) sprintf(fname, "%s/relaxa_%s/traj_%s_%d_orig.xyz", fname, type2, type2, nr);
        else if(strcmp(extra, "log") == 0) sprintf(fname, "%s/%s.dat", fname, extra);
        else sprintf(fname, "%s/relaxa_%s/%s.dat", fname, type2, extra);
    }
    
    return(fname);
}

// Member function of struct log_param
void fileUtils::log_param::read_log(FILE *file){
    char *pipeString = (char*) malloc(500*sizeof(char));
    char *pipeChar = (char*) malloc(500*sizeof(char));
    int tmp;
          
    rewind(file);
    while(!feof(file))
    {
        fgets(pipeString, 500, file);
        sscanf(pipeString, "%s %d", pipeChar, &tmp);
        if(strcmp(pipeChar, name) == 0) val = tmp;
    }
    delete(pipeString, pipeChar);
}
