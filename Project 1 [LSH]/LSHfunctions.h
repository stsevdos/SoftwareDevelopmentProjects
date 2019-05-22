/* Locality Sensitive Hashing functions.h */

#ifndef LSH_LSHfunctions_h
#define LSH_LSHfunctions_h

/* Hamming Space */
int* selectKBits(int N, int k);
int hashFunctionG(item input, int* kBits, int k);
int readHamming(FILE* inputFile, int N, item* currentItem);
int getHammingDimension(char* fileName);
int hammingDistance(item input, item query, int dimension);

/* Vector Spaces ~ Euclidean */
int getVectorDimension(char* fileName, char* metric);
int* euclideanHashFunctionH(vectorHash vh, item vectorP, int dimension, int k);
int euclideanHashFunctionF(vectorHash vh, int k, int *hFunctions);
int readVector(FILE* inputFile, int N, item* currentItem);
double euclideanDistance(item input, item query, int dimension);

/* Distance Based Hashing */
int getDistanceMatrixDimension(char* fileName);
void readDistances(FILE* inputFile, int dimension, dMatrixItem *dmi);
void indexItems(item *temp, double **itemDistances, char **itemNames, int indexPosition);
double** dbhFucntionH(int *x1, int *x2, double **distancesMatrix, int dimension, int k);
void selectt1t2(dMatrixHash *dmh, double **hx1x2, int dimension, int k);
int* discritizeH(dMatrixHash dmh, double **hx1x2, int dimension, int k, int index);
int DBHhashFunctionG(int* kBits, int k);
int queryHashingDBH(double **distancesMatrix, int dimension, dMatrixHash dmh, int k, int hashTableNo);

/* Extras */
void SelectionSort(double* array, int arrayLength);

#endif