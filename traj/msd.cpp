/*
 Code to read trajectory file and calculate
 Mean-Squared-Displacement (MSD) of particles

 ------- Inputs ---------
 nr         - file number
 t_start    - time to start with
 t_end      - time to end 
 frame_w    - time interval corresponding to each frame
 ------------------------
 
 Author        : Ashwin Kumar M
 Date created  : 11.10.24
 Last modified : 11.10.24
*/

#include<iostream>
#include<fstream>
#include<cstring>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<filesystem>

using namespace std;

struct atomCoords{
    float rx_t0 = 0.0, ry_t0 = 0.0;
    float rx_t1 = 0.0, ry_t1 = 0.0;
    float rx_t2 = 0.0, ry_t2 = 0.0;
    int jump_x = 0, jump_y = 0;
    void unwrap_pbc(int boxLength);  
};

void atomCoords::unwrap_pbc(int boxLength)
{
    if(rx_t2 - rx_t1 <= -0.5*boxLength) jump_x += 1;
    else if(rx_t2 - rx_t1 >= 0.5*boxLength) jump_x -= 1;
    if(ry_t2 - ry_t1 <= -0.5*boxLength) jump_y += 1;
    else if(ry_t2 - ry_t1 >= 0.5*boxLength) jump_y -= 1;
}

struct log_param{
    int val;
    char name[10];
    void read_log(FILE *file);
};

void log_param::read_log(FILE *file)
{
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

void computeMeanSquaredDisplacement(float *meanSquaredDisplacement, FILE *file, int frame_start, int nFrames, int nAtoms, int boxLength)
{
    for(int i = 0; i < nFrames; i++) meanSquaredDisplacement[i] = 0.0;
    
    atomCoords *atoms = (atomCoords*) malloc(nAtoms*sizeof(atomCoords));
    char* pipeString = (char*) malloc(500*sizeof(char));
    int frame_nr = 0;
    
    rewind(file);
    while(!feof(file))
    {
        fgets(pipeString, 500, file);
        fgets(pipeString, 500, file);
        frame_nr += 1;
        
        if(frame_nr < frame_start + 1)
        {
            for(int i = 0; i < nAtoms; i++) fgets(pipeString, 500, file);
        }
        
        if(frame_nr == frame_start + 1)
        {
            for(int i = 0; i < nAtoms; i++)
            {
                fgets(pipeString, 500, file);
                sscanf(pipeString, "%*c %f %f %*f", &atoms[i].rx_t0, &atoms[i].ry_t0);
                atoms[i].rx_t1 = atoms[i].rx_t0;
                atoms[i].ry_t1 = atoms[i].ry_t0;
            }
        }
        
        if(frame_nr > frame_start + 1 and frame_nr <= frame_start + nFrames)
        {
            for(int i = 0; i < nAtoms; i++)
            {
                fgets(pipeString, 500, file);
                sscanf(pipeString, "%*c %f %f %*f", &atoms[i].rx_t2, &atoms[i].ry_t2);
                atoms[i].unwrap_pbc(boxLength);
                
                float dx = (atoms[i].rx_t2 + boxLength*atoms[i].jump_x) - atoms[i].rx_t0;
                float dy = (atoms[i].ry_t2 + boxLength*atoms[i].jump_y) - atoms[i].ry_t0;
                
                meanSquaredDisplacement[frame_nr - (frame_start+1)] += dx*dx + dy*dy;
                
                atoms[i].rx_t1 = atoms[i].rx_t2;
                atoms[i].ry_t1 = atoms[i].ry_t2;
            }
            meanSquaredDisplacement[frame_nr - (frame_start+1)] /= nAtoms;
        }
    }
    delete(pipeString, atoms);
}

int main(int argc, char* argv[])
{   
    int nr = atoi(argv[1]);
    float t_start = atof(argv[2]);
    float t_end = atof(argv[3]);
    float frame_w = atof(argv[4]);
    char* type = argv[5];
    
    char* cwd = (char*) malloc(500*sizeof(char));
    string str(filesystem::current_path());
    str.erase(str.begin() + str.find("/ld_analysis"), str.end());
    sprintf(cwd, "%s/LD/Data%d", str.c_str(), nr);
    
    log_param nAtoms = {0, 'N'}, boxL = {0, 'L'};
    
    char* fname = (char*) malloc(500*sizeof(char));
    sprintf(fname, "%s/log.dat", cwd);
    FILE *file = fopen(fname, "r");
    if(file == NULL){ printf("Can not open log file. Exiting ...\n"); exit(-1);}
    else printf("Opening Log file to read params L and N...\n");
    
    nAtoms.read_log(file); 
    boxL.read_log(file);
    fclose(file);
    
    int N_frames = int((t_end - t_start)/frame_w) + 1;
    float* meanSquaredDisplacement = (float*) malloc(N_frames*sizeof(float));
    
    sprintf(fname, "%s/relaxa_%s/traj_%s_%d.xyz", cwd, type, type, nr);
    file = fopen(fname, "r");
    if(file == NULL) {printf("Can not open trajectory file. Exiting ...\n"); exit(-1);}
    else printf("Opening trajectory file to read atomic coordinates...\n");
    
    computeMeanSquaredDisplacement(meanSquaredDisplacement, file, int(t_start/frame_w+0.1), N_frames, nAtoms.val, boxL.val);
    fclose(file);
 
    sprintf(fname, "%s/relaxa_%s/msd%d.dat", cwd, type, nr);
    remove(fname);
    file = fopen(fname, "a+");
    if(file == NULL) {printf("Could not create msd dat file. Exiting...\n"); exit(-1);}
    else printf("Writing to msd%d.dat file...\n", nr); 
    for(int i = 0; i < N_frames; i++) fprintf(file, "%d %f\n", i, meanSquaredDisplacement[i]);
    fclose(file);
    
    return 0;
}
