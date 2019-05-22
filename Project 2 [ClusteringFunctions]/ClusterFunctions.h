
#ifndef __Clustering__ClusterFunctions__
#define __Clustering__ClusterFunctions__

#include <stdio.h>

void readConfigFile(char *configFileName, int *k, int *hashK, int *L, int *maxNeighbor, int *numLocal);
int getItemCount(char *fileName, char *fileType);
int getDimension(char *fileName, char fileType);
double** readItems(char *fileName, char **itemNamesIndex, int dimension, int itemCount);
double** calcDistanceMatrix(double **itemsIndex, char fileType, int dimension, int itemCount);
int sameMedoids(int *oldMedoids, int *newMedoids, int k);

/* INITIALIZATION FUNCTIONS */
int* initializationConcentrate(double **distanceMatrix, int itemCount, int k);
int* initializationSpreadOut(double **distanceMatrix, int itemCount, int k);

/* ASSIGNMENT FUNCTIONS */
int* assignmentPAM(double **distanceMatrix, int* medoids, int itemCount, int k, double *totalCost);
int* assignmentLSH(double **distanceMatrix, int **bucketIndex, int *medoids, int itemCount, double *totalCost, int k, int hashK, int L, double radius, double threshold);

/* UPDATE FUNCTIONS */
int* updateALaLloyds(double **distanceMatrix, int *clustersArray, int *medoids, int itemCount, int k, double *totalCost, int assignment, int **bucketIndex, int hashK, int L, double radius, double threshold);
int* updateClarans(double **distanceMatrix, int *clustersArray, int *medoids, int itemCount, int k, double *totalCost, int numLocal, int maxNieghbor, int assignment, int **bucketIndex, int hashK, int L, double radius, double threshold);

double* silhouette(double **distanceMatrix, int *clustersArray, int *medoids, int itemCount, int k);

void SelectionSortTwoArrays(double *array, int *array2, int arrayLength);
void medoidDistance(double **distanceMatrix, double *minDistance, double *maxDistance, int *medoids, int k);
int* getClusterSize(int *clustersArray, int *medoids, int itemCount, int k);

#endif /* defined(__Clustering__ClusterFunctions__) */
