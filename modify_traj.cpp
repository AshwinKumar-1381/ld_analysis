/*
    Program to work on trajectory files
    Author : Ashwin Kumar M
    Date created : 10.09.24
    Last modified : 10.09.24
    
    =======
    Options
    =======
    A) Merge the original trajectory file (nr) with a restart file (restart_nr)

    ======
    Inputs
    ======
    option - A or B or C ...
    nr - file no. of original trajectory file
    restart_nr - file no. of restart file
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<cstring>
#include<fstream>
#include<iostream>
#include<sys/stat.h>
#include<cmath>
#include<filesystem>

using namespace std;

void merge_traj(int nr, int restart_nr);

int main(int argc, char* argv[])
{
    char option = *argv[1];
    int nr = atoi(argv[2]);
    int restart_nr = atoi(argv[3]);
    
    if(option == 'A')
    {
        merge_traj(nr, restart_nr);
    }
    
    return(0);
}

void merge_traj(int nr, int restart_nr)
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
