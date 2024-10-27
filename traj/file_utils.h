// file_utils.h

#include<iostream>
#include<fstream>
#include<cstring>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<filesystem>

using namespace std;

namespace fileUtils{
    char* makeFilePath(int nr, char *type1, char *type2, string extras);
    
    struct log_param{
        int val;
        char name[10];
        void read_log(FILE *file);
    };
}
