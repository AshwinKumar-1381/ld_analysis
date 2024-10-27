// Read an .xyz trajectory file, extract only certain time frames, and output to new .xyz file
// Date created : 10.07.24
// Last modified : 12.07.24

#include<stdio.h>
#include<iostream>
//#include<omp.h>

using namespace std;

int main()
{
    int N = 2000;                   // No. of particles
    int total_steps = 20001;      // Total no. of frames
    int intval = 10;                // extract frame for every
    int step = -1;                  
    
    char file_old[1000], file_new[1000], string[1000];
    sprintf(file_old, "/home/ethaya_lab_c1/Desktop/ashwin_md/LD/Data3/relaxa/reg/visual_relaxa.xyz");
    sprintf(file_new, "/home/ethaya_lab_c1/Desktop/ashwin_md/LD/Data3/relaxa/reg/visual_relaxa_stp10.xyz");
    
    // create input file object and open file to read
    FILE* dati;
    dati = fopen(file_old, "r");
    if (dati == NULL) printf("error\n");
    else printf("input file is open\n");
    
    // create output file object and create file to write
    FILE* dato;
    dato = fopen(file_new, "a+");
    if(dato == NULL) printf("error\n");
    else printf("output file is open\n");
    
    while(step <= total_steps and feof(dati) == 0)
    {
        step += 1;
        
        if(step % intval == 0)
        {
            for(int i = 1; i <= N + 2; i++)
            {
                fgets(string, 1000, dati);
                fputs(string, dato);
            }
        }
        
        else
        {
            for(int i = 1; i <= N + 2; i++) fgets(string, 1000, dati);
        } 
    }
       
    fclose(dati);
    fclose(dato);
}
