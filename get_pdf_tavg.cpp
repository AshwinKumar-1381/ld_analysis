/*
 Code to read trajectory file and calculate time-averaged  
 particle density function for particles A and B

 -------- Inputs ----------
 nr         - file no. 
 t_ss_start - starting time for time-averaging 
 t_ss_end   - ending time for time-avergaing
 N_every    - perform averaging for every these many time units
 frame_w    - width at which frames are output in original traj file
 bin_w      - width of particle bins

 Author        : Ashwin Kumar M
 Date created  : 07.09.24
 Last modified : 07.09.24
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
    float t_ss_start = atof(argv[2]);
    float t_ss_end = atof(argv[3]);
    int N_every = atoi(argv[4]);
    float frame_w = atof(argv[5]);
    float bin_w = atof(argv[6]);
    
    int i, L, N, frame_nr, N_avg;
    char id;
    float rx, ry, tmp, PeA, PeB;
    
    FILE *file_log, *file_traj, *file_pdf;
    char *cwd = (char*) malloc (500*sizeof(char));
    char *fname = (char*) malloc (500*sizeof(char));
    char *pipeString = (char*) malloc (500*sizeof(char));
    char *pipeChar = (char*) malloc (500*sizeof(char));
    
    string str(filesystem::current_path());
    str.erase(str.begin() + str.find("/code"), str.end());
    strcpy(cwd, str.c_str());

    // Reading log file to obtain values of L and N 
    sprintf(fname, "%s/LD/Data%d/log.dat", cwd, nr);
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
    
    // Bins to store particle densities
    int n_bin = int(L/ bin_w);
    float *bin_A = (float*) malloc (n_bin*sizeof(float));
    float *bin_B = (float*) malloc (n_bin*sizeof(float));
    float *bin_A_avg = (float*) malloc (n_bin*sizeof(float));
    float *bin_B_avg = (float*) malloc (n_bin*sizeof(float));
    for (i = 0; i < n_bin; i++) 
    {
        bin_A[i] = 0.0;
        bin_B[i] = 0.0;
        bin_A_avg[i] = 0.0;
        bin_B_avg[i] = 0.0;
    }
    
    // Remove pdf file if it already exists
    sprintf(fname, "%s/LD/Data%d/relaxa_ss/pdf_%d.dat", cwd, nr, nr);
    remove(fname);
    
    // Reading trajectory file for particle coordinates 
    sprintf(fname, "%s/LD/Data%d/relaxa_ss/traj_ss_%d.xyz", cwd, nr, nr);
    file_traj = fopen(fname, "r");
    if (file_traj == NULL)
    {
        printf("Error. Can't open traj file. \n");
        exit(-1);
    }
    else printf("Opening trajectory file for reading frames.\n\n");
    
    frame_nr = 0;
    N_avg = 0;
    while (!feof(file_traj))
    {
        fgets(pipeString, 500, file_traj);
        fgets(pipeString, 500, file_traj);

        if(frame_nr == int((t_ss_start + N_avg)/frame_w) and frame_nr <= int(t_ss_end/frame_w)) 
        {
            N_avg += N_every;
            // Binning particles 
            for(i = 0; i < N; i++)
            {
                fgets(pipeString, 500, file_traj);
                sscanf(pipeString, "%c %f %f %*f", &id, &rx, &ry);
                if(id == 'N') bin_A[int(rx/bin_w)] += 1;
                else if (id == 'O') bin_B[int(rx/bin_w)] += 1;   
            }
            
            // Calculating PDF
            for(i = 0; i < n_bin; i++)
            {
                bin_A[i] /= (bin_w*L);
                bin_B[i] /= (bin_w*L);
                bin_A_avg[i] += bin_A[i];
                bin_B_avg[i] += bin_B[i];
                bin_A[i] = 0.0;
                bin_B[i] = 0.0;
            }  
        } 
        else
        {
            for(i = 0; i < N; i++) fgets(pipeString, 500, file_traj);
        }
        frame_nr++; 
    }
    fclose(file_traj);
    // Averaging over N_avg frames
    for(i = 0; i < n_bin; i++)
    {
        bin_A_avg[i] /= N_avg;
        bin_B_avg[i] /= N_avg;
    }
    printf("Time averaging performed over %d frames spaced at %d time unit(s).\n", N_avg, N_every);
    
    // Output the avg. PDF to .dat file
    sprintf(fname, "%s/LD/Data%d/relaxa_ss/pdf_%d.dat", cwd, nr, nr);
    file_pdf = fopen(fname, "a+");
    
    if(file_pdf == NULL)
    {
        printf("Can't create new file.\n");
        exit(-1);
    }
    fprintf(file_pdf, "%d %d %d %f %f\n", n_bin, L, N, PeA, PeB);
    fprintf(file_pdf, "%s\n", "bin_no pdfA pdfB");
    
    for(i = 0; i < n_bin; i++) fprintf(file_pdf, "%d %f %f\n", i+1, bin_A_avg[i], bin_B_avg[i]);
    fclose(file_pdf);
     
    return(0);
} 
