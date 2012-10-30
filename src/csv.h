#ifndef CSV_H
#define CSV_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXFLDS 50     /* maximum possible number of fields */
#define MAXFLDSIZE 1000   /* longest possible field + 1 = 31 byte field */
#define MAXFLANK 1000 
void parse( char *record, char *delim, char arr[][MAXFLDSIZE],int *fldcnt)
{
    	char*p=strtok(record,delim);
   	int fld=0;
	while(p!=NULL){
    		strcpy(arr[fld],p);
		fld++;
		p=strtok(NULL,delim);
	}		
	*fldcnt=fld;
}

#endif // CSV_H

