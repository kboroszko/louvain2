#include<stdlib.h>
#include<stdio.h>
#include<string.h>

#define COORDINATE 1
#define ARRAY 2

#define REAL 1
#define INTEGER 2
#define COMPLEX 3
#define PATTERN 4

#define GENERAL 1
#define SYMMETRIC 2
#define SKEW 3
#define HERMITIAN 4

#define THROW(msg, exit_code) printf("ERROR OCCURED\n%s\n", msg); exit(exit_code)

typedef struct {
    int format;
    int valueType;
    int symmetry;
} MFormat;

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
    char formatStr[200];
    char valueTypeStr[20];
    char symmetryStr[20];
    int a = fscanf(filePtr,"%*s %*s %s",formatStr);
    int b = fscanf(filePtr,"%s", valueTypeStr)==1;
    int c = fscanf(filePtr,"%s",symmetryStr)==1;
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




int main(){
    printf("hello world\n");

    FILE* ptr = fopen("mycielskian4.mtx","r");
    if (ptr==NULL)
    {
        printf("no such file.");
        return 1;
    }
    MFormat format = readSignature(ptr);
    printf("readSignature: read format %d %d %d\n", format.format, format.valueType, format.symmetry);



    return 0;
}
