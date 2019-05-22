#include <stdio.h>
#include <stdlib.h>
#include "ClusterFunctions.h"
#include "Bio.h"

int main(int argc, const char * argv[])
{
    int i, j, k, bestK, conformationCount, N, *initMedoids, *clustersArray;
    float **data, **distanceMatrix;
    double totalCost, *silhouetteMatrix, bestSilhouette;
    FILE *outputFile;
    
    /* Get data from file */
    data = getData("bio_small_input.dat", &conformationCount, &N);
    distanceMatrix = malloc(sizeof(float*)*conformationCount);
    if (distanceMatrix == NULL)
    {
        perror("main.c: malloc failed");
        exit(-1);
    }
    /* Calculate distance matrix using dRMSD metric */
    for (i = 0; i < conformationCount; i++)
    {
        distanceMatrix[i] = malloc(sizeof(float)*conformationCount);
        if (distanceMatrix[i] == NULL)
        {
            perror("main.c: malloc failed");
            exit(-1);
        }
        for (j = 0; j < conformationCount; j++)
            distanceMatrix[i][j] = getdRMSD(data[i], data[j], conformationCount);
    }
    
    /* Find best k for clustering */
    bestK = 0;
    bestSilhouette = -1.0;
    for (k = 5; k < 20; k++)
    {
        initMedoids = initializationConcentrate(distanceMatrix, conformationCount, k);
        clustersArray = assignmentPAM(distanceMatrix, initMedoids, conformationCount, k, &totalCost);
        clustersArray = updateClarans(distanceMatrix, clustersArray, initMedoids, conformationCount, k, &totalCost, 10, 5);
        silhouetteMatrix = silhouette(distanceMatrix, clustersArray, initMedoids, conformationCount, k);
        if (silhouetteMatrix[k] > bestSilhouette)   /* if total silhouette (for clustering) is better than the one found before, swap */
        {
            bestK = k;
            bestSilhouette = silhouetteMatrix[k];
        }
        /*printf("k = %d\tsilhouette = %lf\n", k, silhouetteMatrix[k]);*/
        free(initMedoids);  free(clustersArray); free(silhouetteMatrix);
    }
    
    /* Perform clustering algorithms again, because arrays might have been deallocated above */
    initMedoids = initializationConcentrate(distanceMatrix, conformationCount, bestK);
    clustersArray = assignmentPAM(distanceMatrix, initMedoids, conformationCount, bestK, &totalCost);
    clustersArray = updateClarans(distanceMatrix, clustersArray, initMedoids, conformationCount, bestK, &totalCost, 10, 5);
    /* Write to output file */
    outputFile = fopen("conform.dat", "w");
    if (outputFile == NULL)
    {
        perror("main.c: malloc failed");
        exit(-1);
    }
    fprintf(outputFile, "k: %d\n", bestK);
    fprintf(outputFile, "s: %lf\n", bestSilhouette);
    for (i = 0; i < bestK; i++)     /* Print cluster items */
    {
        for (j = 0; j < conformationCount; j++)
        {
            /* if item in in cluster with medoid = initMedoids[i], print to output file */
            if (clustersArray[initMedoids[i]] == clustersArray[j])
                fprintf(outputFile, "%d ", j);
        }
        fprintf(outputFile, "\n");
    }
    
    fclose(outputFile);
    for (i = 0; i < conformationCount; i++)
    {
        free(data[i]);
        free(distanceMatrix[i]);
    }
    free(data); free(distanceMatrix);   free(initMedoids);  free(clustersArray);
    return 0;
}
