/*
 Code to read trajectory file and calculate time-averaged  
 particle currents along y - direction
 
 ------------- Inputs ----------------
 nr         - file no. 
 type1      - lmp/cpp
 type2      - eq/ss
 t_ss_start - starting time for time-averaging 
 t_ss_end   - ending time for time-avergaing
 N_every    - perform averaging for every these many time units
 frame_w    - width at which frames are output in original traj file
 x_start    - coordinate position to start from
 x_end      - coordinate position to end
 bin_w      - width of particle bins
 --------------------------------------
 
 Author         : Ashwin Kumar M
 Date created   : 10.09.24
 Last modified  : 16.10.24
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

int main(int argc, char *argv[])
{
    int nr = atoi(argv[1]);
    char *type1 = argv[2];
    char *type2 = argv[3];
    float t_ss_start = atof(argv[4]);
    float t_ss_end = atof(argv[5]);
    int N_every = atoi(argv[6]);
    float frame_w = atof(argv[7]);
    float x_start = atof(argv[8]);
    float x_end = atof(argv[9]);
    float bin_w = atof(argv[10]);
    
    int i, L, N, frame_nr, N_avg;
    char id;
    float dely, tmp, PeA, PeB;
    
    FILE *file_log, *file_traj, *file_current;
    fpos_t pos;
    char *cwd = (char*) malloc (500*sizeof(char));
    char *fname = (char*) malloc (500*sizeof(char));
    char *pipeString = (char*) malloc (500*sizeof(char));
    char *pipeChar = (char*) malloc (500*sizeof(char));
    
    string str(filesystem::current_path());
    str.erase(str.begin() + str.find("/ld_analysis/traj"), str.end());
    strcpy(cwd, str.c_str());

    // Reading log file to obtain values of L and N 
    sprintf(fname, "%s/LD/%s/Data%d/log.dat", cwd, type1, nr);
    file_log = fopen(fname, "r");
    if (file_log == NULL)
    {
        printf("Error. Log File not found. \n");
        exit(-1);    
    }
    
    while (!feof(file_log))
    {
        fgets(pipeString, 500, file_log);
        sscanf(pipeString, "%s %f", pipeChar, &tmp);
        if(strcmp(pipeChar, "Pe_A") == 0) PeA = tmp;
        if(strcmp(pipeChar, "Pe_B") == 0) PeB = tmp;
        if(strcmp(pipeChar, "L") == 0) L = int(tmp);
        if(strcmp(pipeChar, "N") == 0) N = int(tmp);
    }
    fclose(file_log);

    // Bins to store particle positions
    int n_bin = int((x_end - x_start)/ bin_w);
    float *rx = (float*) malloc (N*sizeof(float));
    float *ry_t = (float*) malloc (N*sizeof(float));
    float *ry_t_dt = (float*) malloc (N*sizeof(float));
    int *bin_count = (int*) malloc (n_bin*sizeof(int));
    float *vel = (float*) malloc (n_bin*sizeof(float));
    float *vel_avg = (float*) malloc (n_bin*sizeof(float));
    for (i = 0; i < N; i++)
    {
        rx[i] = 0.0;
        ry_t[i] = 0.0;
        ry_t_dt[i] = 0.0;     
    }
    for (i = 0; i < n_bin; i++) 
    {
        bin_count[i] = 0;
        vel[i] = 0.0;
        vel_avg[i] = 0.0;
    }

    // Reading trajectory file for particle coordinates
    if(strcmp(type1, "lmp") == 0) sprintf(fname, "%s/LD/%s/Data%d/traj%d.xyz", cwd, type1, nr, nr);
    else sprintf(fname, "%s/LD/%s/Data%d/relaxa_%s/traj_%s_%d.xyz", cwd, type1, nr, type2, type2, nr);
    
    file_traj = fopen(fname, "r");
    if (file_traj == NULL)
    {
        printf("Error. Can't open traj file. \n");
        exit(-1);
    }

    frame_nr = 0;
    N_avg = 0;
    while(!feof(file_traj))
    {
        fgets(pipeString, 500, file_traj);
        fgets(pipeString, 500, file_traj);
        frame_nr++;
        
        if(frame_nr == int((t_ss_start + N_avg)/frame_w) and frame_nr < int(t_ss_end/frame_w))      
        {
            N_avg += N_every;
            // Reading frame at time t 
            for(i = 0; i < N; i++)
            {
                fgets(pipeString, 500, file_traj);
                sscanf(pipeString, "%*c %f %f %*f", &rx[i], &ry_t[i]);   
            }
            
            fgetpos (file_traj, &pos);
            
            int nSkip = int(N_every/frame_w) - 1;
            // Skipping nSkip frames 
            for(i = 0; i < nSkip*(N + 2); i++) fgets(pipeString, 500, file_traj);
            
            // Reading frame at time t + N_every
            fgets(pipeString, 500, file_traj);
            fgets(pipeString, 500, file_traj);
            for(i = 0; i < N; i++)
            {
                fgets(pipeString, 500, file_traj);
                sscanf(pipeString, "%*c %*f %f %*f", &ry_t_dt[i]);
            }
            
            fsetpos (file_traj, &pos);
            
            for(i = 0; i < N; i++)
            {
                // Using un-wrapped coordinates
                if(abs(ry_t_dt[i] - ry_t[i]) >= 0.5 * L)
                {
                    if(ry_t[i] > ry_t_dt[i]) ry_t_dt[i] += L;
                    else ry_t_dt[i] -= L; 
                }
                
                // Calculate particle displacements
                if(rx[i] >= x_start and rx[i] <= x_end)
                {
                    bin_count[int((rx[i] - x_start)/bin_w)] += 1;  
                    dely = ry_t_dt[i] - ry_t[i];
                    vel[int((rx[i] - x_start)/bin_w)] += dely;          
                }        
            }
            
            for(i = 0; i < n_bin; i++)
            {
                vel[i] /= bin_count[i];
                vel_avg[i] += vel[i];
                vel[i] = 0.0;
            }
        }
        else
        {
            for(i = 0; i < N; i++) fgets(pipeString, 500, file_traj);
        }
    }
    fclose(file_traj);    
    for (i = 0; i < n_bin; i++) vel_avg[i] = (vel_avg[i])/N_avg;
    
    // Output particle currents to .dat file
    if(strcmp(type1, "lmp") == 0) sprintf(fname, "%s/LD/%s/Data%d/current%d.dat", cwd, type1, nr, nr);
    else sprintf(fname, "%s/LD/%s/Data%d/relaxa_%s/current%d.dat", cwd, type1, nr, type2, nr);
    remove(fname);
    
    file_current = fopen(fname, "a+");
    if(file_current == NULL)
    {
        printf("Can't create file to output current data.\n");
        exit(-1);
    }
    
    fprintf(file_current, "%f %f %d %d %d %f %f\n", x_start, x_end, n_bin, L, N, PeA, PeB);
    fprintf(file_current, "bin currrent\n");
    for(i = 0; i < n_bin; i++) fprintf(file_current, "%d %f\n", i+1, vel_avg[i]);
    
    fclose(file_current);

    return(0);
}
