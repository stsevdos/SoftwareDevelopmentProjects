#ifndef __Recommendation__Recommendation__
#define __Recommendation__Recommendation__

#include <stdio.h>
#include "ADT.h"

int getUserCount(char *fileName);
int getMusicCount(char *fileName);
void getMusicRatings(char *fileName, short **ratings);
float* userRatingMean(short **ratings, int userCount, int musicCount);
double** recommendItems(int userCount, int musicCount, short **pNearest, short **ratings, float *r, char fileType, FILE *outputFile);
double crossValidation(short **pNearest, short **ratings, float *r, int userCount, int musicCount, double** pNNdist);


#endif /* defined(__Recommendation__Recommendation__) */
