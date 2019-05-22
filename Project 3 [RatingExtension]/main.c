#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include "ADT.h"
#include "LSHAlgorithm.h"
#include "LSHFunctions.h"
#include "ClusterFunctions.h"
#include "Recommendation.h"

#define k 3
#define L 2

int P;

int main(int argc, const char * argv[])
{
    int i, j, userCount, musicCount, *initMedoids, medoids, *clustersArray, validate;
    short *users, **ratings, **pNearest;
    float *r;
    double **distanceMatrix, totalCost, **pNNdist, MAE, elapsedTime;
    char *inputFileName, *outputFileName;
    FILE *outputFile;
    struct timeval time1, time2;
    
    if (argc > 6 || argc < 5)
    {
        printf("main.c: Need more arguments to run program.\n");
        exit(-1);
    }
    else
    {
        inputFileName = malloc(sizeof(char)*200);
        if (inputFileName == NULL)
        {
            perror("main.c: malloc failed");
            exit(-1);
        }
        outputFileName = malloc(sizeof(char)*200);
        if (outputFileName == NULL)
        {
            perror("main.c: malloc failed");
            exit(-1);
        }
        
        validate = 0;
        for (i = 1; i < argc; i += 2)
        {
            if (argv[i][1] == 'd')
                strcpy(inputFileName, argv[i+1]);
            else if (argv[i][1] == 'o')
                strcpy(outputFileName, argv[i+1]);
            else if (strcmp(argv[i], "-validate") == 0)
                validate = 1;
            else
            {
                printf("Unknown flag %s\n", argv[i]);
                exit(-1);
            }
            
        }
    }
    
    srand((unsigned int)time(NULL));
    
    userCount = getUserCount(inputFileName);
    users = malloc(sizeof(short)*userCount);
    if (users == NULL)
    {
        perror("main.c: malloc failed");
        exit(-1);
    }
    ratings = malloc(sizeof(short*)*userCount);
    if (ratings == NULL)
    {
        perror("main.c: malloc failed");
        exit(-1);
    }
    
    musicCount = getMusicCount(inputFileName);
    for (i = 0; i < userCount; i++)
    {
        ratings[i] = malloc(sizeof(short)*musicCount);
        if (ratings[i] == NULL)
        {
            perror("main.c: malloc failed");
            exit(-1);
        }
        for (j = 0; j < musicCount; j++)
            ratings[i][j] = 0;     /* Initialy, all musics are considered not rated */
    }
    getMusicRatings(inputFileName, ratings);    /* Get music ratings from input file */
    r = userRatingMean(ratings, userCount, musicCount);         /* Calculate mean rating for each user */
    distanceMatrix = calcDistanceMatrix(ratings, 'c', musicCount, userCount);
    
    outputFile = fopen(outputFileName, "w");
    if (outputFile == NULL)
    {
        perror("main.c: fopen failed");
        exit(-1);
    }
    
    /* NN LSH: Euclidean */
    fprintf(outputFile, "NN LSH (Euclidean)\n");
    gettimeofday(&time1, NULL); /* START TIMER */
    pNearest = LSHAlgorithm(musicCount, userCount, 'v', k, L, ratings, r);
    pNNdist = recommendItems(userCount, musicCount, pNearest, ratings, r, 'v', outputFile);
    gettimeofday(&time2, NULL); /* STOP TIMER */
    elapsedTime = (time2.tv_sec - time1.tv_sec) * 1000.0;
    elapsedTime += (time2.tv_usec - time1.tv_usec) / 1000.0;
    fprintf(outputFile, "Execution Time: %lf\n", elapsedTime);
    if (validate == 1)
    {
        MAE = crossValidation(pNearest, ratings, r, userCount, musicCount, pNNdist);
        fprintf(outputFile, "NN LSH (Euclidean) Recommendation MAE: %lf\n", MAE);
    }
    
    for (i = 0; i < userCount; i++)
    {
        free(pNearest[i]);
        free(pNNdist[i]);
    }
    free(pNearest); free(pNNdist);
    
    /* CLUSTERING: Euclidean */
    fprintf(outputFile, "Clustering (Euclidean)\n");
    medoids = 10;
    initMedoids = malloc(sizeof(int)*medoids);
    if (initMedoids == NULL)
    {
        perror("main.c: malloc failed");
        exit(-1);
    }
    initMedoids = initializationConcentrate(distanceMatrix, userCount, k);
    
    gettimeofday(&time1, NULL); /* START TIMER */
    assignmentPAM(distanceMatrix, initMedoids, userCount, k, &totalCost);
    pNearest = updateClarans(distanceMatrix, clustersArray, initMedoids, userCount, k, &totalCost, 5, 10, L, 'v');
    pNNdist = recommendItems(userCount, musicCount, pNearest, ratings, r, 'c', outputFile);
    gettimeofday(&time2, NULL); /* STOP TIMER */
    elapsedTime = (time2.tv_sec - time1.tv_sec) * 1000.0;
    elapsedTime += (time2.tv_usec - time1.tv_usec) / 1000.0;
    fprintf(outputFile, "Execution Time: %lf\n", elapsedTime);
    if (validate == 1)
    {
        MAE = crossValidation(pNearest, ratings, r, userCount, musicCount, pNNdist);
        fprintf(outputFile, "Clustering (Euclidean) Recommendation MAE: %lf\n", MAE);
    }
    
    for (i = 0; i < userCount; i++)
    {
        free(pNearest[i]);
        free(pNNdist[i]);
    }
    free(pNearest); free(pNNdist);
    fclose(outputFile);
    
    strcat(outputFileName, "2");
    outputFile = fopen(outputFileName, "w");
    if (outputFile == NULL)
    {
        perror("main.c: fopen failed");
        exit(-1);
    }
    /* NN LSH: Cosine */
    fprintf(outputFile, "NN LSH (Cosine)\n");
    gettimeofday(&time1, NULL); /* START TIMER */
    pNearest = LSHAlgorithm(musicCount, userCount, 'c', k, L, ratings, r);
    pNNdist = recommendItems(userCount, musicCount, pNearest, ratings, r, 'c', outputFile);
    gettimeofday(&time2, NULL); /* STOP TIMER */
    elapsedTime = (time2.tv_sec - time1.tv_sec) * 1000.0;
    elapsedTime += (time2.tv_usec - time1.tv_usec) / 1000.0;
    fprintf(outputFile, "Execution Time: %lf\n", elapsedTime);
    if (validate == 1)
    {
        MAE = crossValidation(pNearest, ratings, r, userCount, musicCount, pNNdist);
        fprintf(outputFile, "NN LSH (Cosine) Recommendation MAE: %lf\n", MAE);
    }

    for (i = 0; i < userCount; i++)
    {
        free(pNearest[i]);
        free(pNNdist[i]);
    }
    free(pNearest); free(pNNdist);
    
    /* CLUSTERING: Cosine */
    fprintf(outputFile, "Clustering (Cosine)\n");
    initMedoids = initializationConcentrate(distanceMatrix, userCount, k);
    gettimeofday(&time1, NULL); /* START TIMER */
    assignmentPAM(distanceMatrix, initMedoids, userCount, k, &totalCost);
    pNearest = updateClarans(distanceMatrix, clustersArray, initMedoids, userCount, k, &totalCost, 5, 10, L, 'v');
    pNNdist = recommendItems(userCount, musicCount, pNearest, ratings, r, 'c', outputFile);
    gettimeofday(&time2, NULL); /* STOP TIMER */
    elapsedTime = (time2.tv_sec - time1.tv_sec) * 1000.0;
    elapsedTime += (time2.tv_usec - time1.tv_usec) / 1000.0;
    fprintf(outputFile, "Execution Time: %lf\n", elapsedTime);
    if (validate == 1)
    {
        MAE = crossValidation(pNearest, ratings, r, userCount, musicCount, pNNdist);
        fprintf(outputFile, "Clustering (Cosine) Recommendation MAE: %lf\n", MAE);
    }
    
    for (i = 0; i < userCount; i++)
    {
        free(pNearest[i]);
        free(pNNdist[i]);
    }
    free(pNearest); free(pNNdist);
    fclose(outputFile);
    
    /* Deallocating memory */
    for (i = 0; i < userCount; i++)
        free(ratings[i]);
    
    free(r);    free(ratings);  free(users);    free(inputFileName);
    free(initMedoids);
    
    return 0;
}