//
// Created by kajetan on 30.04.2020.
//

#ifndef LOUVAIN2_MMIO_H
#define LOUVAIN2_MMIO_H

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

void validateFormat(MFormat form);

MFormat readSignature(FILE* filePtr);


#endif //LOUVAIN2_MMIO_H
