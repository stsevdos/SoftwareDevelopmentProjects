/* Clustering main */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include "LSHfunctions.h"
#include "ADT.h"
#include "LSHAlgorithm.h"
#include "ClusterFunctions.h"

int main(int argc, const char * argv[])
{
    int i, k, itemCount, dimension, *oldMedoids, *clustersArray, **bucketIndex, hashK, L, maxNeighbor, numLocal, *clusterSize;
    double **itemIndex, **distanceMatrix, totalCost, maxDistance, minDistance, radius, threshold, *evaluation, elapsedTime;
    char fileType, **itemNamesIndex, inputFileName[100], configFileName[100], outputFileName[100];
    FILE *outputFile;
    struct timeval time1, time2;

    /* Default Values */
    k = 5;
    hashK = 4;
    L = 5;
    maxNeighbor = 0;    /* Set 0 to calculate later if not in configuration file */
    numLocal = 2;
    
    srand((unsigned int)time(NULL));
    
    /* Read file names */
    if (argc <= 1)
    {
        printf("Too few arguments. Cannot execute program without all parameters passed.\n");
        exit(-1);
    }
    else
    {
        for (i = 1; i < argc; i += 2)
        {
            if (argv[i][1] == 'd')
                strcpy(inputFileName, argv[i+1]);
            else if (argv[i][1] == 'c')
                strcpy(configFileName, argv[i+1]);
            else if (argv[i][1] == 'o')
                strcpy(outputFileName, argv[i+1]);
            else
            {
                printf("Unknown flag %c\n", argv[i][1]);
                exit(-1);
            }
        }
    }

    /* Read parameters from configuration file */
    readConfigFile(configFileName, &k, &hashK, &L, &maxNeighbor, &numLocal);
    itemCount = getItemCount(inputFileName, &fileType);
    /* If clarans_set_fraction is not defined in config file, calculate now */
    if (maxNeighbor == 0)
    {
        if (250 > (0.12*k*(itemCount-k)))
            maxNeighbor = 250;
        else
            maxNeighbor = 0.12*k*(itemCount-k);
    }

    /* Create item index */
    itemNamesIndex = malloc(sizeof(char*)*itemCount);
    if (itemNamesIndex == NULL)
    {
        printf("main.c: main: Failed to allocate memory for itemNamesIndex\n");
        exit(-1);
    }
    
    for (i = 0; i < itemCount; i++)
    {
        itemNamesIndex[i] = malloc(sizeof(double)*100);
        if (itemNamesIndex[i] == NULL)
        {
            printf("main.c: main: Failed to allocate memory for itemNamesIndex\n");
            exit(-1);
        }

    }
    dimension = getDimension(inputFileName, fileType);
    itemIndex = readItems(inputFileName, itemNamesIndex, dimension, itemCount);
    distanceMatrix = calcDistanceMatrix(itemIndex, fileType, dimension, itemCount);
    bucketIndex = LSHAlgorithm(dimension, itemCount, fileType, hashK, L, itemIndex);
    outputFile = fopen(outputFileName, "w");
    if (outputFile == NULL)
    {
        printf("main.c: main: Failed to open (output) file %s\n", outputFileName);
        exit(-1);
    }

    /* Combination #1 */
    gettimeofday(&time1, NULL); /* START TIMER */
    oldMedoids = initializationConcentrate(distanceMatrix, itemCount, k);
    clustersArray = assignmentPAM(distanceMatrix, oldMedoids, itemCount, k, &totalCost);
    clustersArray = updateALaLloyds(distanceMatrix, clustersArray, oldMedoids, itemCount, k, &totalCost, 0, bucketIndex, hashK, L, radius, threshold);
    gettimeofday(&time2, NULL); /* STOP TIMER */
    elapsedTime = (time2.tv_sec - time1.tv_sec) * 1000.0;
    elapsedTime += (time2.tv_usec - time1.tv_usec) / 1000.0;
    clusterSize = getClusterSize(clustersArray, oldMedoids, itemCount, k);
    evaluation = silhouette(distanceMatrix, clustersArray, oldMedoids, itemCount, k);
    fprintf(outputFile, "Algorithm: I1A1U1\n");
    for (i = 0; i < k; i++)
        fprintf(outputFile, "CLUSTER-%d {size: %d, medoid: %s}\n", i+1, clusterSize[i], itemNamesIndex[oldMedoids[i]]);
    fprintf(outputFile, "clustering_time: %lf\n", elapsedTime);
    fprintf(outputFile, "Silhouette: [");
    for (i = 0; i < k; i++)
        fprintf(outputFile, "%lf, ", evaluation[i]);
    fprintf(outputFile, "%lf]\n", evaluation[k]);
    //free(evaluation);   free(clusterSize);  free(oldMedoids);   //free(clustersArray);

    /* Combination #2 */
    gettimeofday(&time1, NULL); /* START TIMER */
    oldMedoids = initializationConcentrate(distanceMatrix, itemCount, k);
    clustersArray = assignmentPAM(distanceMatrix, oldMedoids, itemCount, k, &totalCost);
    clustersArray = updateClarans(distanceMatrix, clustersArray, oldMedoids, itemCount, k, &totalCost, numLocal, maxNeighbor, 0, bucketIndex, hashK, L, radius, threshold);
    gettimeofday(&time2, NULL); /* STOP TIMER */
    elapsedTime = (time2.tv_sec - time1.tv_sec) * 1000.0;
    elapsedTime += (time2.tv_usec - time1.tv_usec) / 1000.0;
    clusterSize = getClusterSize(clustersArray, oldMedoids, itemCount, k);
    evaluation = silhouette(distanceMatrix, clustersArray, oldMedoids, itemCount, k);
    fprintf(outputFile, "Algorithm: I1A1U2\n");
    for (i = 0; i < k; i++)
        fprintf(outputFile, "CLUSTER-%d {size: %d, medoid: %s}\n", i+1, clusterSize[i], itemNamesIndex[oldMedoids[i]]);
    fprintf(outputFile, "clustering_time: %lf\n", elapsedTime);
    fprintf(outputFile, "Silhouette: [");
    for (i = 0; i < k; i++)
        fprintf(outputFile, "%lf, ", evaluation[i]);
    fprintf(outputFile, "%lf]\n", evaluation[k]);
    //free(evaluation);   free(clusterSize);  free(oldMedoids);   //free(clustersArray);

    /* Combination #3 */
    gettimeofday(&time1, NULL); /* START TIMER */
    oldMedoids = initializationConcentrate(distanceMatrix, itemCount, k);
    medoidDistance(distanceMatrix, &minDistance, &maxDistance, oldMedoids, k);
    radius = minDistance/2.0;
    threshold = 5.0 * maxDistance;
    clustersArray = assignmentLSH(distanceMatrix, bucketIndex, oldMedoids, itemCount, &totalCost, k, hashK, L, radius, threshold);
    clustersArray = updateALaLloyds(distanceMatrix, clustersArray, oldMedoids, itemCount, k, &totalCost, 1, bucketIndex, hashK, L, radius, threshold);
    gettimeofday(&time2, NULL); /* STOP TIMER */
    elapsedTime = (time2.tv_sec - time1.tv_sec) * 1000.0;
    elapsedTime += (time2.tv_usec - time1.tv_usec) / 1000.0;
    clusterSize = getClusterSize(clustersArray, oldMedoids, itemCount, k);
    evaluation = silhouette(distanceMatrix, clustersArray, oldMedoids, itemCount, k);
    fprintf(outputFile, "Algorithm: I1A2U1\n");
    for (i = 0; i < k; i++)
        fprintf(outputFile, "CLUSTER-%d {size: %d, medoid: %s}\n", i+1, clusterSize[i], itemNamesIndex[oldMedoids[i]]);
    fprintf(outputFile, "clustering_time: %lf\n", elapsedTime);
    fprintf(outputFile, "Silhouette: [");
    for (i = 0; i < k; i++)
        fprintf(outputFile, "%lf, ", evaluation[i]);
    fprintf(outputFile, "%lf]\n", evaluation[k]);
    //free(evaluation);   free(clusterSize);  free(oldMedoids);  // free(clustersArray);

    /* Combination #4 */
    gettimeofday(&time1, NULL); /* START TIMER */
    oldMedoids = initializationConcentrate(distanceMatrix, itemCount, k);
    medoidDistance(distanceMatrix, &minDistance, &maxDistance, oldMedoids, k);
    radius = minDistance/2.0;
    threshold = 5.0 * maxDistance;
    clustersArray = assignmentLSH(distanceMatrix, bucketIndex, oldMedoids, itemCount, &totalCost, k, hashK, L, radius, threshold);
    clustersArray = updateClarans(distanceMatrix, clustersArray, oldMedoids, itemCount, k, &totalCost, numLocal, maxNeighbor, 1, bucketIndex, hashK, L, radius, threshold);
    gettimeofday(&time2, NULL); /* STOP TIMER */
    elapsedTime = (time2.tv_sec - time1.tv_sec) * 1000.0;
    elapsedTime += (time2.tv_usec - time1.tv_usec) / 1000.0;
    clusterSize = getClusterSize(clustersArray, oldMedoids, itemCount, k);
    evaluation = silhouette(distanceMatrix, clustersArray, oldMedoids, itemCount, k);
    fprintf(outputFile, "Algorithm: I1A2U2\n");
    for (i = 0; i < k; i++)
        fprintf(outputFile, "CLUSTER-%d {size: %d, medoid: %s}\n", i+1, clusterSize[i], itemNamesIndex[oldMedoids[i]]);
    fprintf(outputFile, "clustering_time: %lf\n", elapsedTime);
    fprintf(outputFile, "Silhouette: [");
    for (i = 0; i < k; i++)
        fprintf(outputFile, "%lf, ", evaluation[i]);
    fprintf(outputFile, "%lf]\n", evaluation[k]);
    //free(evaluation);   free(clusterSize);  free(oldMedoids);   //free(clustersArray);

    /* Combination #5 */
    gettimeofday(&time1, NULL); /* START TIMER */
    oldMedoids = initializationSpreadOut(distanceMatrix, itemCount, k);
    clustersArray = assignmentPAM(distanceMatrix, oldMedoids, itemCount, k, &totalCost);
    clustersArray = updateALaLloyds(distanceMatrix, clustersArray, oldMedoids, itemCount, k, &totalCost, 0, bucketIndex, hashK, L, radius, threshold);
    gettimeofday(&time2, NULL); /* STOP TIMER */
    elapsedTime = (time2.tv_sec - time1.tv_sec) * 1000.0;
    elapsedTime += (time2.tv_usec - time1.tv_usec) / 1000.0;
    clusterSize = getClusterSize(clustersArray, oldMedoids, itemCount, k);
    evaluation = silhouette(distanceMatrix, clustersArray, oldMedoids, itemCount, k);
    fprintf(outputFile, "Algorithm: I2A1U1\n");
    for (i = 0; i < k; i++)
        fprintf(outputFile, "CLUSTER-%d {size: %d, medoid: %s}\n", i+1, clusterSize[i], itemNamesIndex[oldMedoids[i]]);
    fprintf(outputFile, "clustering_time: %lf\n", elapsedTime);
    fprintf(outputFile, "Silhouette: [");
    for (i = 0; i < k; i++)
        fprintf(outputFile, "%lf, ", evaluation[i]);
    fprintf(outputFile, "%lf]\n", evaluation[k]);
    //free(evaluation);   free(clusterSize);  free(oldMedoids);   //free(clustersArray);

    /* Combination #6 */
    gettimeofday(&time1, NULL); /* START TIMER */
    oldMedoids = initializationSpreadOut(distanceMatrix, itemCount, k);
    clustersArray = assignmentPAM(distanceMatrix, oldMedoids, itemCount, k, &totalCost);
    clustersArray = updateClarans(distanceMatrix, clustersArray, oldMedoids, itemCount, k, &totalCost, numLocal, maxNeighbor, 0, bucketIndex, hashK, L, radius, threshold);
    gettimeofday(&time2, NULL); /* STOP TIMER */
    elapsedTime = (time2.tv_sec - time1.tv_sec) * 1000.0;
    elapsedTime += (time2.tv_usec - time1.tv_usec) / 1000.0;
    clusterSize = getClusterSize(clustersArray, oldMedoids, itemCount, k);
    evaluation = silhouette(distanceMatrix, clustersArray, oldMedoids, itemCount, k);
    fprintf(outputFile, "Algorithm: I2A1U2\n");
    for (i = 0; i < k; i++)
        fprintf(outputFile, "CLUSTER-%d {size: %d, medoid: %s}\n", i+1, clusterSize[i], itemNamesIndex[oldMedoids[i]]);
    fprintf(outputFile, "clustering_time: %lf\n", elapsedTime);
    fprintf(outputFile, "Silhouette: [");
    for (i = 0; i < k; i++)
        fprintf(outputFile, "%lf, ", evaluation[i]);
    fprintf(outputFile, "%lf]\n", evaluation[k]);
    //free(evaluation);   free(clusterSize);  free(oldMedoids);   //free(clustersArray);

    /* Combination #7 */
    gettimeofday(&time1, NULL); /* START TIMER */
    oldMedoids = initializationSpreadOut(distanceMatrix, itemCount, k);
    medoidDistance(distanceMatrix, &minDistance, &maxDistance, oldMedoids, k);
    radius = minDistance/2.0;
    threshold = 5.0 * maxDistance;
    clustersArray = assignmentLSH(distanceMatrix, bucketIndex, oldMedoids, itemCount, &totalCost, k, hashK, L, radius, threshold);
    clustersArray = updateALaLloyds(distanceMatrix, clustersArray, oldMedoids, itemCount, k, &totalCost, 1, bucketIndex, hashK, L, radius, threshold);
    gettimeofday(&time2, NULL); /* STOP TIMER */
    elapsedTime = (time2.tv_sec - time1.tv_sec) * 1000.0;
    elapsedTime += (time2.tv_usec - time1.tv_usec) / 1000.0;
    clusterSize = getClusterSize(clustersArray, oldMedoids, itemCount, k);
    evaluation = silhouette(distanceMatrix, clustersArray, oldMedoids, itemCount, k);
    fprintf(outputFile, "Algorithm: I2A2U1\n");
    for (i = 0; i < k; i++)
        fprintf(outputFile, "CLUSTER-%d {size: %d, medoid: %s}\n", i+1, clusterSize[i], itemNamesIndex[oldMedoids[i]]);
    fprintf(outputFile, "clustering_time: %lf\n", elapsedTime);
    fprintf(outputFile, "Silhouette: [");
    for (i = 0; i < k; i++)
        fprintf(outputFile, "%lf, ", evaluation[i]);
    fprintf(outputFile, "%lf]\n", evaluation[k]);
    //free(evaluation);   free(clusterSize);  free(oldMedoids);  // free(clustersArray);

    /* Combination #8 */
    gettimeofday(&time1, NULL); /* START TIMER */
    oldMedoids = initializationSpreadOut(distanceMatrix, itemCount, k);
    medoidDistance(distanceMatrix, &minDistance, &maxDistance, oldMedoids, k);
    radius = minDistance/2.0;
    threshold = 5.0 * maxDistance;
    clustersArray = assignmentLSH(distanceMatrix, bucketIndex, oldMedoids, itemCount, &totalCost, k, hashK, L, radius, threshold);
    clustersArray = updateClarans(distanceMatrix, clustersArray, oldMedoids, itemCount, k, &totalCost, numLocal, maxNeighbor, 1, bucketIndex, hashK, L, radius, threshold);
    gettimeofday(&time2, NULL); /* STOP TIMER */
    elapsedTime = (time2.tv_sec - time1.tv_sec) * 1000.0;
    elapsedTime += (time2.tv_usec - time1.tv_usec) / 1000.0;
    clusterSize = getClusterSize(clustersArray, oldMedoids, itemCount, k);
    evaluation = silhouette(distanceMatrix, clustersArray, oldMedoids, itemCount, k);
    fprintf(outputFile, "Algorithm: I2A2U2\n");
    for (i = 0; i < k; i++)
        fprintf(outputFile, "CLUSTER-%d {size: %d, medoid: %s}\n", i+1, clusterSize[i], itemNamesIndex[oldMedoids[i]]);
    fprintf(outputFile, "clustering_time: %lf\n", elapsedTime);
    fprintf(outputFile, "Silhouette: [");
    for (i = 0; i < k; i++)
        fprintf(outputFile, "%lf, ", evaluation[i]);
    fprintf(outputFile, "%lf]\n", evaluation[k]);

    /* Deallocating memory */
    for (i = 0; i < itemCount; i++)
    {
        free(itemIndex[i]);
        free(itemNamesIndex[i]);
        if (fileType != 'm')    /* If input file is distance matrix, distanceMatrix = itemIndex */
            free(distanceMatrix[i]);
    }
    free(itemIndex);    free(itemNamesIndex);
    if (fileType != 'm')    /* If input file is distance matrix, distanceMatrix = itemIndex */
        free(distanceMatrix);
    for (i = 0; i < L; i++)
        free(bucketIndex[i]);
    free(evaluation);   free(clusterSize);  free(bucketIndex);  free(oldMedoids);
    fclose(outputFile);

    return 0;
}