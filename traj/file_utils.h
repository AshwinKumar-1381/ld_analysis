// file_utils.h

namespace fileUtils{
    
    struct log_param{
        int val;
        char name[10];
        void read_log(FILE *file);
    };
    
    char* makeFilePath(int nr, char *type1, char *type2, const char *extras);
}
