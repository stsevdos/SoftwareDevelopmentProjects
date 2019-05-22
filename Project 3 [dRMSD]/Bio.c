#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Bio.h"

float** getData(char *fileName, int *conformationCount, int *N)
{
    int i, j;
    float **data;
    char *buffer;
    FILE *inputFile;
    
    inputFile = fopen(fileName, "r");
    if (inputFile == NULL)
    {
        perror("Bio.c: getData: fopen failed");
        exit(-1);
    }
    
    buffer = malloc(sizeof(char)*200);
    if (buffer == NULL)
    {
        perror("Bio.c: getData: malloc failed");
        exit(-1);
    }
    /* Read first line containing number of conformations */
    fscanf(inputFile, "%s", buffer);
    *conformationCount = atoi(buffer);
    /* Read second line containing number of triads per conformation */
    fscanf(inputFile, "%s", buffer);
    *N = atoi(buffer);

    data = malloc(sizeof(float*)*(*conformationCount));
    if (data == NULL)
    {
        perror("Bio.c: getData: malloc failed");
        exit(-1);
    }
    
    for (i = 0; i < *conformationCount; i++)
    {
        data[i] = malloc(sizeof(float)*((*N)*3));
        if (data[i] == NULL)
        {
            perror("Bio.c: getData: malloc failed");
            exit(-1);
        }
        
        for (j = 0; j < (*N)*3; j++)
        {
            fscanf(inputFile, "%s", buffer);
            data[i][j] = atof(buffer);
        }
    }
    
    fclose(inputFile);  free(buffer);
    
    return data;
}

float getdRMSD(float *conf1, float *conf2, int confSize)
{
    int i, j;
    float dRMSD;
    
    dRMSD = 0.0;
    for (i = 0; i <confSize; i++)
    {
        for (j = 0; j < confSize; j++)
        {
            if (i == j)
                continue;
            dRMSD += (conf1[3*i] - conf2[3*i])*(conf1[3*i] - conf2[3*i]);
            dRMSD += (conf1[3*i+1] - conf2[3*i+1])*(conf1[3*i+1] - conf2[3*i+1]);
            dRMSD += (conf1[3*i+2] - conf2[3*i+2])*(conf1[3*i+2] - conf2[3*i+2]);
        }
        dRMSD = sqrt(dRMSD);
        dRMSD *= 1.0/(confSize*(confSize-1));
    }
        
    return dRMSD;
}