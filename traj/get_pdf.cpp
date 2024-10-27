// Code to read full trajectory file and output the density
// distribution for particles A and B at stipulated times
// Author : Ashwin Kumar M
// Date created : 16.08.24
// Last modified : 17.08.24

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string>
#include<fstream>
#include<iostream>
#include<sys/stat.h>
#include<cmath>

using namespace std;

char check_t(float t_out[], float t, int l);

int main(int argc, char *argv[])
{
    int i, N;
    char pdf_file[100], fname[100], str[1000];
    char type;
    float L, rx, tmp, t_relax_ss, t;
    double *bin_A, *bin_B;
    
    int nr = 6;                // file no.
    double w_bin = 2.0;         // bin width 
    float t_int = 0.1;          // time interval of frames in traj file
    
    float t_out[9] = {0.0, 0.1, 0.5, 1.0, 5.0, 10.0, 50.0, 100.0, 500.0};
    int l = int(sizeof(t_out)/sizeof(t_out[0]));                // size of t_out array
    // Converting time to frame no.
    for(i = 0; i < l; i++) t_out[i] = int(t_out[i]/t_int) + 1; 
    
    // Reading elog file to obtain values of L and N 
    sprintf(fname, "/home/ethaya_lab_c1/Desktop/ashwin_md/LD/Aug24/Data%d/elog.dat", nr);
    FILE *log;
    log = fopen(fname, "r");
    if (log == NULL)
    {
        printf("Error. Log File not found. \n");
        exit(-1);    
    }
    
    i = 0;
    while (!feof(log))
    {
        fscanf (log, "%f", &tmp);
        i += 1;
        if (i == 2) t_relax_ss = tmp;
        if (i == 9) L = tmp;
        if (i == 10) N = tmp;
    }
    fclose(log);
    
    // Bins to store particle densities
    int n_bin = int(L/ w_bin);
    bin_A = new double[n_bin];
    bin_B = new double[n_bin];
    
    // Remove pdf file if it already exists
    sprintf(pdf_file, "/home/ethaya_lab_c1/Desktop/ashwin_md/LD/Aug24/Data%d/relaxa_ss/pdf_%d.dat", nr, nr);
    remove(pdf_file);
    
    // Reading trajectory file
    sprintf(fname, "/home/ethaya_lab_c1/Desktop/ashwin_md/LD/Aug24/Data%d/relaxa_ss/traj_ss_%d.xyz", nr, nr);
    FILE *data;
    data = fopen(fname, "r");
    if(data == NULL)
    {
        printf("Error. Trajectory file not found.\n");
        exit(-1);
    }

    t = 1;
    while(!feof(data))
    {   
        if(check_t(t_out, t, l) == 'Y')
        {
            fgets(str, 1000, data);
            //cout << str;
            fgets(str, 1000, data);
        
            for (i = 0; i < n_bin; i++) 
            {
                bin_A[i] = 0.0;
                bin_B[i] = 0.0;
            }
            
            for(i = 1; i <= N; i++)
            {
                fgets(str, 1000, data);
                sscanf(str, "%s %12f %*12f %*3d", &type, &rx);
                if(type == 'N') bin_A[int(rx/w_bin)] += 1;
                else bin_B[int(rx/w_bin)] += 1;
            }
        
            for (i = 0; i < n_bin; i++) 
            {
                bin_A[i]/= (n_bin * w_bin);
                bin_B[i]/= (n_bin * w_bin);
            }
            
            FILE *data_out;
            data_out = fopen(pdf_file, "a+");
            if(data_out == NULL)
            {
                printf("Could not open pdf file\n");
                exit(-1);
            }
            else
            {
                //cout << t << endl;
                fprintf(data_out, "%f\n", (t-1)*t_int);
                for(i = 0; i < n_bin; i++) fprintf(data_out, "%d %6f %6f\n", i+1, bin_A[i], bin_B[i]);      
            }
            fclose(data_out);
            
        }
        else
        {  
            for (i = 1; i <= N + 2; i++) fgets(str, 1000, data);     
       
        }
        t += 1;
        
    }
    fclose(data);
    
    FILE *data_out;
    data_out = fopen(pdf_file, "a+");
    if(data_out == NULL)
    {
        printf("Could not open pdf file\n");
        exit(-1);
    }
    else fprintf(data_out, "%d", nr);
    fclose(data_out);
    
    return(0);
} 

char check_t(float t_out[], float t, int l)
{
    for (int j = 0; j < l; j++)
    {
        if(t == t_out[j]) return('Y');
    }       
    
    return('N');
}
