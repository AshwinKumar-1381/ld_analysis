#include<iostream>
#include<stdlib.h>
#include<stdio.h>
#include<cstring>
#include<string.h>
#include<filesystem>
#include<fstream>

#define nTypes 2

struct particleBins{
    int nBin = 1;
    float *bins;  
    
    void initialize ();  
    void normalize (float factor);
    void copy (particleBins Bin);
};

void particleBins::initialize()
{
    bins = (float*) malloc(nBin*sizeof(float));
    for(int i = 0; i < nBin; i++) bins[i] = 0.0;
}

void particleBins::normalize(float factor)
{
    for (int i = 0; i < nBin; i++) bins[i] /= factor;
}

void particleBins::copy(particleBins Bin)
{
    for (int i = 0; i < nBin; i++) bins[i] += Bin.bins[i]; 
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

void computeParticleDistribution(FILE *file, int frame_start, int frame_end, int nFrames, int nAtoms, int boxL, float bin_w)
{   
    particleBins *localBins, *globalBins = (particleBins*) malloc(nTypes*sizeof(particleBins*));
    for(int i = 0; i < nTypes; i++)
    {
        localBins[i].nBin = int(boxL/bin_w);
        globalBins[i].nBin = localBins[i].nBin;
        localBins[i].initialize();
        globalBins[i].initialize();
    }
    
           
    char *pipeString = (char*) malloc(500*sizeof(char));
    int frame_nr = 0;
    float rx; char id;
    
    /*
    rewind(file);
    while (!feof(file))
    {
         fgets(pipeString, 500, file);
         fgets(pipeString, 500, file);   
         frame_nr ++;
         
         if (frame_nr > frame_start and frame_nr <= frame_end)
         {
            for (int i = 0; i < nAtoms; i++)
            {
                fgets(pipeString, 500, file);
                sscanf("%c %f %*f %*f", &id, &rx);
                if(id == 'N') localBins[0].bins[int(rx/bin_w)] += 1;
                else localBins[1].bins[int(rx/bin_w)] += 1;
            }
            
            for (int i = 0; i < nTypes; i++) 
            {
                localBins[i].normalize(boxL*bin_w);
                
                globalBins[i].copy(localBins[i]);
                localBins[i].initialize();
            }            
         }
         
         else
         {
            for (int i = 0; i < nAtoms; i++) fgets(pipeString, 500, file);
         }
    }
    
    for (int i = 0; i < globalBins[0].nBin; i++) printf("%d %f %f\n", i+1, globalBins[0].bins[i], globalBins[1].bins[i]);
    for (int i = 0; i < nTypes; i++) globalBins[i].normalize(nFrames);
    
    */
    
    delete(localBins);
}

using namespace std;

int main(int argc, char *argv[])
{
    int nr = atoi(argv[1]);
    float t_start = atof(argv[2]);
    float t_end = atof(argv[3]);
    float frame_w = atof(argv[4]);
    float bin_w = atof(argv[5]);
    char *type1 = argv[6];
    char *type2 = argv[7];
      
    char *fname = (char*) malloc(500*sizeof(char));
    string cwd(filesystem::current_path());
    cwd.erase(cwd.begin() + cwd.find("/ld_analysis"), cwd.end());
    
    sprintf(fname, "%s/LD/%s/Data%d/log.dat", cwd.c_str(), type1, nr);
    FILE *file = fopen(fname, "r");
    if(file == NULL) {printf("Cannot find log file. Exiting ...\n"); exit(-1);}
    
    log_param nAtoms = {0, 'N'}, boxL = {0, 'L'};
    nAtoms.read_log(file);
    boxL.read_log(file);
    fclose(file);    

    particleBins *globalBins = (particleBins*) malloc (nTypes* sizeof(particleBins*));   
    
    int frame_start = int(t_start/frame_w);
    int frame_end = int(t_end/frame_w);
    int nFrames = frame_end - frame_start + 1;
    
    if(strcmp(type1, "lmp") == 0) sprintf(fname, "%s/LD/%s/Data%d/traj%d.xyz", cwd.c_str(), type1, nr, nr);
    else sprintf(fname, "%s/LD/%s/Data%d/relaxa_%s/traj%d.xyz", cwd.c_str(), type1, nr, type2, nr);
    
    file = fopen(fname, "r");
    if(file == NULL) {printf("Cannot open trajectory file. Exiting...\n"); exit(-1);}
    else computeParticleDistribution(file, frame_start, frame_end, nFrames, nAtoms.val, boxL.val, bin_w);
    fclose(file);
    
    return(0);
}
