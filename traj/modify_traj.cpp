/*     
 -------------- Options --------------
 A) Merge the original trajectory file (nr) with a restart file (restart_nr)
 -------------------------------------

 -------------- Inputs ---------------
 option - A or B or C ...
 nr - file no. of original trajectory file
 type1 - lmp/cpp
 type2 - eq/ss
 
 Other arguments:
 
 option A - restart_nr
 option B - frameStart frameEnd newFrames
 -------------------------------------
    
 Program to work on trajectory files
 Author : Ashwin Kumar M
 Date created : 10.09.24
 Last modified : 25.10.24
*/

#include "file_utils.h"

using namespace fileUtils;

int countFrames(char *trajname, char *logname){
    FILE *logfile = fopen(logname, "r");
    fileUtils::log_param nAtoms = {0, 'N'};
    nAtoms.read_log(logfile);
    
    FILE *trajfile = fopen(trajname, "r");
    if(trajfile == NULL) {printf("Cannot open trajectory file. Exiting...\n"); exit(-1);}
    
    int nFrames = 0, temp = 0;
    char *pipeString = (char*) malloc(500*sizeof(char));
    while(!feof(trajfile)){
        fgets(pipeString, 500, trajfile);
        sscanf(pipeString, "%d", &temp);
        if(temp == nAtoms.val){
            nFrames++;
            for (int i = 0; i < nAtoms.val + 1; i++) fgets(pipeString, 500, trajfile);
        }
    }
    
    fclose(logfile);
    fclose(trajfile);
    delete(logfile, trajfile, pipeString);
    return(nFrames-1);
}

void mergeTraj(int nr, char *type1, char *type2, int restart_nr){

    FILE *file_ss, *file_res;
    char *cwd = (char*) malloc (500*sizeof(char));
    char *pipeString = (char*) malloc (500*sizeof(char));
    char *fname_ss = (char*) malloc (500*sizeof(char));
    char *fname_res = (char*) malloc (500*sizeof(char));
    
    string str(filesystem::current_path());
    str.erase(str.begin() + str.find("/code"), str.end());
    strcpy(cwd, str.c_str());
    
    sprintf(fname_res, "%s/LD/Data%d/restart_%d/traj_ss_%d.xyz", cwd, nr, restart_nr, nr);
    file_res = fopen(fname_res, "r");
    if(file_res == NULL)
    {
        printf("Restart files not found.\n");
        exit(-1);
    }
    else printf("Opening Restart file.\n");
    
    sprintf(fname_ss, "%s/LD/Data%d/relaxa_ss/traj_ss_%d.xyz", cwd, nr, nr);
    file_ss = fopen(fname_ss, "a+");
    if(file_ss == NULL)
    {
        printf("Original trajectory not found.\n");
        exit(-1);
    }
    else printf("Opening Trajectory file for writing.\n");
    
    while(!feof(file_res))
    {
        fgets(pipeString, 500, file_res);   
        fputs(pipeString, file_ss); 
    }
    
    fclose(file_res);
    fclose(file_ss);
}

void extractTraj(int nr, char *type1, char *type2, int frameStart, int frameEnd, int newFrames){
    
    string traj("traj"), log("log"), traj_orig("traj_orig");
    char *trajfileOld = fileUtils::makeFilePath(nr, type1, type2, traj_orig);
    char *trajfileNew = fileUtils::makeFilePath(nr, type1, type2, traj);
    char *logfile = fileUtils::makeFilePath(nr, type1, type2, log);
    char *pipeString = (char*) malloc(500*sizeof(char));
    
    sprintf(pipeString, "cp %s %s", )
    
    int nFrames = countFrames(trajfileOld, logfile);
    int nSkip = int(nFrames/newFrames);
    int frame_nr;
    printf("%d %d %d\n", nFrames, newFrames, nSkip);
    
    
    FILE *trajOld = fopen(trajfileOld, "r");
    printf("Opening original trajectory with %d frames...\n", nFrames);
    while(!feof(trajOld)){
        
    }
    
    
    fclose(trajOld);
}

int main(int argc, char* argv[]){

    char option = *argv[1];
    int nr = atoi(argv[2]);
    char *type1 = argv[3];
    char *type2 = argv[4];
    
    if(option == 'A'){
        int restart_nr = atoi(argv[5]);
        mergeTraj(nr, type1, type2, restart_nr);    
    }
    else if(option == 'B'){
        int frameStart = atoi(argv[5]);
        int frameEnd = atoi(argv[6]);
        int newFrames = atoi(argv[7]);
        extractTraj(nr, type1, type2, frameStart, frameEnd, newFrames);
    }
    
    return(0);
}


