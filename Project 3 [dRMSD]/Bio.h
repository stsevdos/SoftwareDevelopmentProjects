#ifndef __bio__Bio__
#define __bio__Bio__

#include <stdio.h>

float** getData(char *fileName, int *conformationCount, int *N);
float getdRMSD(float *conf1, float *conf2, int confSize);

#endif /* defined(__bio__Bio__) */
