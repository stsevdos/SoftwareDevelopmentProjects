#ifndef __Recommendation__LSHFunctions__
#define __Recommendation__LSHFunctions__

#include <stdio.h>

int* selectKBits(int N, int k);

/* Distance functions */
double euclideanDistance(int user1, int user2, short **ratings, int dimension);
double cosineSimilarity(int user1, int user2, short **ratings, int dimension);

/* Hash functions */
int* euclideanHashFunctionH(vectorHash vh, short *user, int dimension, int k);
int euclideanHashFunctionF(vectorHash vh, int k, int *hFunctions);
short** cosineRandomVectors(int k, int dimension);
int cosineHashFunction(int k, int dimension, short **csVectors, short *ratings, float *userMean);

double initialRadius(short **ratings, int userCount, int musicCount, int userID, char metric);

#endif /* defined(__Recommendation__LSHFunctions__) */
