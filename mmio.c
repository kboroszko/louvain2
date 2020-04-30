#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include "mmio.h"

void validateFormat(MFormat form){
    int dat = form.format;
    int val = form.valueType;
    int sym = form.symmetry;
    if(sym == HERMITIAN && val != COMPLEX){
        THROW("validateFormat: hermitian must be complex!", 1);
    }
    if(val == PATTERN && sym != SYMMETRIC && sym != GENERAL && dat != COORDINATE){
        THROW("validateFormat: invalid data structure and/or symmetries for type pattern!", 2);
    }
}

MFormat readSignature(FILE* filePtr){
    MFormat ret;
    char formatStr[20];
    char valueTypeStr[20];
    char symmetryStr[20];
    int a = fscanf(filePtr,"%*s %*s %s",formatStr);
    int b = fscanf(filePtr,"%s", valueTypeStr)==1;
    int c = fscanf(filePtr,"%s\n",symmetryStr)==1;
    if(a==1 && b==1 && c==1){
        printf("readSignature: read format %s %s %s\n", formatStr, valueTypeStr, symmetryStr);
    } else {
        printf("readSignature: invalid file or format\n");
        exit(1);
    }

    if(strcmp("coordinate", formatStr) == 0){
         ret.format = COORDINATE;
    } else if(strcmp("array", formatStr) == 0){
        ret.format = ARRAY;
    } else {
        THROW("readSignature: invalid data structure!", 2);
    }

    if(strcmp("real", valueTypeStr) == 0){
        ret.valueType = REAL;
    } else if(strcmp("integer", valueTypeStr) == 0) {
        ret.valueType = INTEGER;
    } else if(strcmp("complex", valueTypeStr) == 0){
        ret.valueType = COMPLEX;
    }  else if(strcmp("pattern", valueTypeStr) == 0){
        ret.valueType = PATTERN;
    } else {
        THROW("readSignature: invalid value type!", 3);
    }

    if(strcmp("general",symmetryStr) == 0){
        ret.symmetry = GENERAL;
    } else if (strcmp("symmetric",symmetryStr)  == 0){
        ret.symmetry = SYMMETRIC;
    } else if(strcmp("skew-symmetric",symmetryStr)  == 0){
        ret.symmetry = SKEW;
    } else if(strcmp("hermitian",symmetryStr)  == 0){
        ret.symmetry = HERMITIAN;
    } else {
        THROW("readSignature: invalid symetries type", 4);
    }
    validateFormat(ret);
    return ret;
}

MData * readData(const char * filename){
    FILE* filePtr = fopen("mycielskian4.mtx","r");
    if (filePtr==NULL)
    {
        THROW("readData: no such file.", 10);
    }

    MFormat format = readSignature(filePtr);
    char line[1024];
    int num_cols, num_rows, num_lines;
//    fscanf(filePtr,"%s",line);
    while(fscanf(filePtr,"%1024[^\n]\n",line)!= EOF){
        if(line[0] != '%' && line[0] != '\000'){
            int success = -1;
            if(format.format == COORDINATE){
                success = sscanf(line, "%d %d %d", &num_rows, &num_cols, &num_lines);
            } else {
                success = sscanf(line, "%d %d", &num_rows, &num_cols);
                num_lines = num_rows * num_cols;
            }
            if(success != 3 && success != 2){
                THROW("readData: error reading matrix size!", 11);
            }
            break;
        }
    }

    printf("reading rows=%d  cols=%d  lines=%d", num_rows, num_cols, num_lines);
    return initMData(num_rows, num_cols, num_lines, format);
}


int main(){
    printf("hello world\n");

    MData * dat = readData("mycielskian4.mtx");

    destroyMData(dat);



    return 0;
}
