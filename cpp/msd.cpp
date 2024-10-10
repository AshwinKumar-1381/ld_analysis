/*
 Code to read trajectory file and calculate
 Mean-Squared-Displacement (MSD) of particles

 ------- Inputs ---------
 nr         - file number
 t_start    - time to start with
 t_end      - time to end 
 frame_w    - time interval corresponding to each frame
 
 Author        : Ashwin Kumar M
 Date created  : 10.10.24
 Last modified : 10.10.24
*/

#include<iostream>
#include<fstream>
#include<cstring>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<filesystem>

using namespace std;

int main(int argc, char *argv[])
{
    int nr = atoi(argv[1]);
    int t_start = atoi(argv[2]);
    int t_end = atoi(argv[3]);
    float frame_w = atof(argv[4]);
    
    int L, N, tmp;
    int frame_start = int(t_start/frame_w + 0.1);
    int frame_end = int(t_end/frame_w + 0.1);
    
    FILE *file_log, *file_traj, *file_msd;
    char* pipeChar = (char*) malloc (500*sizeof(char));
    char* pipeString = (char*) malloc (500*sizeof(char));
    char* fname = (char*) malloc (500*sizeof(char));
    char* cwd = (char*) malloc (500*sizeof(char));
    
    float* rx_t0 = (float*) malloc(N*sizeof(float));
    float* ry_t0 = (float*) malloc(N*sizeof(float)); 
    float* rx_t1 = (float*) malloc(N*sizeof(float));
    float* ry_t1 = (float*) malloc(N*sizeof(float)); 
    float* rx_t2 = (float*) malloc(N*sizeof(float));
    float* ry_t2 = (float*) malloc(N*sizeof(float));
    int* jump_index_x = (int*) malloc(N*sizeof(int));
    int* jump_index_y = (int*) malloc(N*sizeof(int));
    float* msd = (float*) malloc((frame_end - frame_start + 1)*sizeof(float));
    
    for(int i = 0; i < N; i++)
    {
        jump_index_x[i] = 0;
        jump_index_y[i] = 0;
    }
    
    // Set the current working directory to Data folder
    string str(filesystem::current_path());
    str.erase(str.begin() + str.find("/code"), str.end());
    strcpy(cwd, str.c_str()); 
    sprintf(cwd, "%s/LD/Data%d", str.c_str(), nr);
    
    // Opening log file to read parameters
    sprintf(fname, "%s/log.dat", cwd);
    file_log = fopen(fname, "r");
    if(file_log == NULL)
    {
        printf("Could not open log file.\n");
        exit(-1);
    }
    
    while(!feof(file_log))
    {
        fgets(pipeString, 500, file_log);
        sscanf(pipeString, "%s %d", pipeChar, &tmp);
        if(strcmp(pipeChar,"L") == 0) L = tmp;
        if(strcmp(pipeChar,"N") == 0) N = tmp;
    }
    fclose(file_log);
    
    // read trajectory file
    sprintf(fname, "%s/relaxa_eq/traj_eq_%d.xyz", cwd, nr);
    file_traj = fopen(fname, "r");
    if(file_traj == NULL)
    {
        printf("Cannot open trajectory file.\n");
        exit(-1);
    }
    
    int frame_nr = -1;
    for(int i = 0; i <= frame_end - frame_start; i++) msd[i] = 0.0;
    printf("%d %d %d\n", frame_start, frame_end, frame_end - frame_start);
    
    
    while(!feof(file_traj))
    {
        fgets(pipeString, 500, file_traj);
        fgets(pipeString, 500, file_traj);
        frame_nr++;
        
        printf("%f\n", msd[int(frame_nr - frame_start)]);
        // storing particle coordinates of reference frame
        if(frame_nr == frame_start)
        {
            for(int i = 0; i < N; i++)
            {
                fgets(pipeString, 500, file_traj);
                sscanf(pipeString, "%*c %f %f %*f", &rx_t0[i], &ry_t0[i]);
                rx_t1[i] = rx_t0[i];
                ry_t1[i] = ry_t0[i];
            }            
        }
        
        else if(frame_nr > frame_start and frame_nr <= frame_end)
        {
            for(int i = 0; i < N; i++)
            {
                fgets(pipeString, 500, file_traj);
                sscanf(pipeString, "%*c %f %f %*f", &rx_t2[i], &ry_t2[i]);
            }
            
            // calculate MSD
            for (int i = 0; i < N; i++) 
            {
                float dx = rx_t2[i] - rx_t1[i];
                float dy = ry_t2[i] - ry_t1[i];
                
                // Un-wrapping PBCs
                if(dx <= -0.5*L) jump_index_x[i] += 1;
                else if(dx >= 0.5*L) jump_index_x[i] -= 1;
                dx = (rx_t2[i] + L*jump_index_x[i]) - rx_t0[i];
                
                if(dy <= -0.5*L) jump_index_y[i] += 1;
                else if(dy >= 0.5*L) jump_index_y[i] -= 1;
                dy = (ry_t2[i] + L*jump_index_y[i]) - ry_t0[i];
                
                msd[frame_nr - frame_start] += dx*dx + dy*dy;
                rx_t1[i] = rx_t2[i];
                ry_t1[i] = ry_t2[i];
            }
            
            msd[frame_nr - frame_start] /= N;
        }
        
        else
        {
            for(int i = 0; i < N; i++) fgets(pipeString, 500, file_traj);
        }
    }
    
    
    
    // Remove msd file if it already exists
    sprintf(fname, "%s/relaxa_eq/msd%d.dat", cwd, nr);
    remove(fname);
    
    file_msd = fopen(fname, "a+");
    if(file_msd == NULL)
    {
        printf("Cannot create MSD file.\n");
        exit(-1);
    }
    fprintf(file_msd, "frame_nr msd\n");
    for(int i = 0; i <= frame_end - frame_start; i++) fprintf(file_msd, "%d %f\n", i + 1, msd[i]);
    fclose(file_msd);
}
