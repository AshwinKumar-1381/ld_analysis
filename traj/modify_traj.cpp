/*       
 -------------- Options --------------
 A) Count the number of frames in trajectory file
 B) Merge the original trajectory file (nr) with a restart file (restart_nr)
 C) Extract certain frames from trajectory file

 -------------- Inputs ---------------
 option - A or B or C ...
 nr - file no. of original trajectory file
 type1 - lmp/cpp
 type2 - eq/ss
 
 Other arguments:
 
 option A - none
 option B - restart_nr
 option C - frameStart frameEnd newFrames
 -------------------------------------
  
 Program to work on trajectory files
 Author : Ashwin Kumar M
 Date created : 10.09.24
 Last modified : 27.10.24
*/

#include "library.h"
#include "file_utils.h"
#include "../../LD/LD-cpp/src/system.h"

using namespace fileUtils;

namespace fileUtils{
    int* countFrames(int nr, char *type1, char *type2, const char *extras); 
}

int* fileUtils::countFrames(int nr, char *type1, char *type2, const char *extras)
{
    string log("log");
    char *fname;
    fname = fileUtils::makeFilePath(nr, type1, type2, log.c_str());
    FILE *logFile = fopen(fname, "r");
    if(logFile == NULL){ printf("Log file not found. Error!\n"); exit(-1);}
    
    fileUtils::log_param nAtoms = {0, 'N'};
    nAtoms.read_log(logFile);
    fclose(logFile);
    
    fname = fileUtils::makeFilePath(nr, type1, type2, extras);
    FILE *trajFile = fopen(fname, "r");
    if(trajFile == NULL){ printf("Trajectory file not found. Error!\n"); exit(-1);}
    
    int nFrames = 0, temp = 0;
    char *pipeString = (char*) malloc(500*sizeof(char));
    int *vars = (int*) malloc(2*sizeof(int));
    
    rewind(trajFile);
    while(!feof(trajFile))
    {
        fgets(pipeString, 500, trajFile);
        sscanf(pipeString, "%d", &temp);
        if(temp == nAtoms.val)
        {
            nFrames++;
            for (int i = 0; i < nAtoms.val + 1; i++) fgets(pipeString, 500, trajFile);
        }
    }
    fclose(trajFile);
    
    vars[0] = nFrames-1;
    vars[1] = nAtoms.val;
    
    delete(pipeString);
    return(vars);
}

void mergeTraj(int nr, char *type1, char *type2, int restart_nr)
{
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
    if(file_res == NULL) { printf("Restart files not found.\n"); exit(-1);}
    else printf("Opening Restart file.\n");
    
    sprintf(fname_ss, "%s/LD/Data%d/relaxa_ss/traj_ss_%d.xyz", cwd, nr, nr);
    file_ss = fopen(fname_ss, "a+");
    if(file_ss == NULL){ printf("Original trajectory not found.\n"); exit(-1);}
    else printf("Opening Trajectory file for writing.\n");
    
    while(!feof(file_res))
    {
        fgets(pipeString, 500, file_res);   
        fputs(pipeString, file_ss); 
    }
    
    fclose(file_res);
    fclose(file_ss);
}

void extractTraj(int nr, char *type1, char *type2, int frameStart, int frameEnd, int newFrames)
{
    string traj("traj"), traj_orig("traj_orig");
    char *fname;
    
    fname = fileUtils::makeFilePath(nr, type1, type2, traj_orig.c_str());
    FILE *trajFileOld = fopen(fname, "r");
    if(trajFileOld == NULL){ printf("File %s can\'t be found. Error!", fname); exit(-1);}
    
    fname = fileUtils::makeFilePath(nr, type1, type2, traj.c_str());
    remove(fname);
    
    FILE *trajFileNew = fopen(fname, "a+");
    if(trajFileNew == NULL){ printf("File %s can\'t be created. Error!", fname); exit(-1);}
    
    char *pipeString = (char*) malloc(500*sizeof(char));
    int *vars = countFrames(nr, type1, type2, traj_orig.c_str());
    int nFrames = vars[0];
    int nSkip = int(nFrames/newFrames);
    int frame_nr = 0;
    
    printf("Number of frames in original trajectory = %d\n", vars[0]);
    
    rewind(trajFileOld);
    while(!feof(trajFileOld))
    {
        frame_nr++;
        if(frame_nr >= frameStart and frame_nr <= frameEnd and frame_nr%nSkip == 1)
        {
            for(int i = 0; i < vars[1] + 2; i++)
            {
                fgets(pipeString, 500, trajFileOld);
                fputs(pipeString, trajFileNew);    
            }
        }
        else
        {
            for(int i = 0; i < vars[1] + 2; i++) fgets(pipeString, 500, trajFileOld);
        }
    }
    fclose(trajFileNew);
    fclose(trajFileOld);
    
    vars = countFrames(nr, type1, type2, traj.c_str());
    printf("Number of frames in modified trajectory = %d\n", vars[0]);
}

int main(int argc, char* argv[])
{
    char option = *argv[1];
    int nr = atoi(argv[2]);
    char *type1 = argv[3];
    char *type2 = argv[4];
    
    if(option == 'A')
    {
        char *extras = argv[5];
        int *vars = fileUtils::countFrames(nr, type1, type2, extras);
        printf("The requested trajectory file has %d frames.\n", vars[0]);
    }
    else if(option == 'B')
    {
        int restart_nr = atoi(argv[5]);
        mergeTraj(nr, type1, type2, restart_nr);    
    }
    else if(option == 'C')
    {
        int frameStart = atoi(argv[5]);
        int frameEnd = atoi(argv[6]);
        int newFrames = atoi(argv[7]);
        extractTraj(nr, type1, type2, frameStart, frameEnd, newFrames);
    }
    
    return(0);
}
