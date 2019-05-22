
#ifndef __Clustering__ClusterFunctions__
#define __Clustering__ClusterFunctions__

#include <stdio.h>

int* selectKBits(int N, int k);
int sameMedoids(int *oldMedoids, int *newMedoids, int k);

/* INITIALIZATION FUNCTIONS */
int* initializationConcentrate(float **distanceMatrix, int itemCount, int k);

/* ASSIGNMENT FUNCTIONS */
int* assignmentPAM(float **distanceMatrix, int* medoids, int itemCount, int k, double *totalCost);

/* UPDATE FUNCTIONS */
int* updateClarans(float **distanceMatrix, int *clustersArray, int *medoids, int itemCount, int k, double *totalCost, int numLocal, int maxNeighbor);

double* silhouette(float **distanceMatrix, int *clustersArray, int *medoids, int itemCount, int k);

void SelectionSortTwoArrays(double *array, int *array2, int arrayLength);
void medoidDistance(float **distanceMatrix, double *minDistance, double *maxDistance, int *medoids, int k);

#endif /* defined(__Clustering__ClusterFunctions__) */
