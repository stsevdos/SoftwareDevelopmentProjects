
#ifndef __Clustering__ClusterFunctions__
#define __Clustering__ClusterFunctions__

#include <stdio.h>

double** calcDistanceMatrix(short **ratings, char fileType, int dimension, int itemCount);
int sameMedoids(int *oldMedoids, int *newMedoids, int k);

/* INITIALIZATION FUNCTIONS */
int* initializationConcentrate(double **distanceMatrix, int itemCount, int k);

/* ASSIGNMENT FUNCTIONS */
int* assignmentPAM(double **distanceMatrix, int* medoids, int itemCount, int k, double *totalCost);

/* UPDATE FUNCTIONS */
short** updateClarans(double **distanceMatrix, int *clustersArray, int *medoids, int itemCount, int k, double *totalCost, int numLocal, int maxNeighbor, int L, char fileType);
double* silhouette(double **distanceMatrix, int *clustersArray, int *medoids, int itemCount, int k);

void SelectionSortTwoArrays(double *array, int *array2, int arrayLength);
void medoidDistance(double **distanceMatrix, double *minDistance, double *maxDistance, int *medoids, int k);

#endif /* defined(__Clustering__ClusterFunctions__) */
