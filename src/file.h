#ifndef FILE_H
#define FILE_H
#include "./headers.h"

FILE *openfile(string filepath, char * opt, bool quiet=true);
FILE *openfile64(string filepath, char * opt, bool quiet=true);

FILE *openfile(string filepath, char * opt, bool quiet)
{
	FILE *file = fopen(filepath.c_str(),opt);
	if (file == NULL){
		fprintf(stderr,"ERROR: fail to open %s.\n",filepath.c_str());
		exit(1);
	}
	if (quiet == false)
		printf("Open file %s.\n",filepath.c_str());
	return file;	
}

FILE *openfile64(string filepath, char * opt, bool quiet)
{
        FILE *file = fopen64(filepath.c_str(),opt);
        if (file == NULL){
                fprintf(stderr,"ERROR: fail to open %s.\n",filepath.c_str());
		exit(1);
        }
        if (quiet == false)
                printf("Open file %s.\n",filepath.c_str());
        return file;
}
#endif // FILE_H
