#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "ClusterFunctions.h"
#include "ADT.h"
#include "LSHfunctions.h"

extern int P;

double** calcDistanceMatrix(short **ratings, char fileType, int dimension, int itemCount)
{
    int i, j;
    double **distanceMatrix;
    char *stringBuffer;
    
    stringBuffer = malloc(sizeof(char)*100);
    if (stringBuffer == NULL)
    {
        printf("ClusterFunctions.c: calcDistanceMatrix: Failed to allocate memory for stringBufffer\n");
        exit(-1);
    }
    distanceMatrix = malloc(sizeof(double*)*itemCount);
    if (distanceMatrix == NULL)
    {
        printf("ClusterFunctions.c: calcDistanceMatrix: Failed to allocate memory for distanceMatrix\n");
        exit(-1);
    }
    
    if (fileType == 'v')
    {
        /* If items are in Vector Space */
        for (i = 0; i < itemCount; i++)
        {
            distanceMatrix[i] = malloc(sizeof(double)*itemCount);
            if (distanceMatrix[i] == NULL)
            {
                printf("ClusterFunctions.c: calcDistanceMatrix: Failed to allocate memory for distanceMatrix[%d]\n", i);
                exit(-1);
            }
            /* Use euclideanDistnace to calculate distance for distance matrix */
            for (j = 0; j < itemCount; j++)
                distanceMatrix[i][j] = (double)euclideanDistance(i, j, ratings, dimension);
        }
    }
    else    /* Calculate distance matrix for cosine similarity */
    {
        for (i = 0; i < itemCount; i++)
        {
            distanceMatrix[i] = malloc(sizeof(double)*itemCount);
            if (distanceMatrix[i] == NULL)
            {
                printf("ClusterFunctions.c: calcDistanceMatrix: Failed to allocate memory for distanceMatrix[%d]\n", i);
                exit(-1);
            }
            /* Use cosineSimilarity to calculate distance for distance matrix */
            for (j = 0; j < itemCount; j++)
                distanceMatrix[i][j] = cosineSimilarity(i, j, ratings, dimension);
        }
    }
    
    free(stringBuffer);
    
    return distanceMatrix;
}

int sameMedoids(int *oldMedoids, int *newMedoids, int k)
{
    int i, j, notExists;
    
    notExists = 1;
    /* For each old medoid */
    for (i = 0; i < k; i++)
    {
        notExists = 1;
        /* Check if it has been selected before in any position of newMedoids */
        for (j = 0; j < k; j++)
        {
            if (oldMedoids[i] == newMedoids[j])
            {
                /* If found, continue checking for other old medoids */
                notExists = 0;
                break;
            }
        }
        /* If not found, arrays don't have the same elements, return 0 */
        if (notExists == 1)
            return 0;
    }
    return 1;
}

int* initializationConcentrate(double **distanceMatrix, int itemCount, int k)
{
    int i, j;
    int *initMedoids, *vIndexes;
    double *v, *distSums;
    
    v = malloc(sizeof(double)*itemCount);
    if (v == NULL)
    {
        printf("ClusterFunctions.c: initializationConcentrate: Failed to allocate memory for v array\n");
        exit(-1);
    }
    vIndexes = malloc(sizeof(int)*itemCount);
    if (vIndexes == NULL)
    {
        printf("ClusterFunctions.c: initializationConcentrate: Failed to allocate memory for vIndexes array\n");
        exit(-1);
    }
    distSums = malloc(sizeof(double)*itemCount);
    if (distSums == NULL)
    {
        printf("ClusterFunctions.c: initializationConcentrate: Failed to allocate memory for distSums array\n");
        exit(-1);
    }
    initMedoids = malloc(sizeof(int)*k);
    if (initMedoids == NULL)
    {
        printf("ClusterFunctions.c: initializationConcentrate: Failed to allocate memory for initMedoids\n");
        exit(-1);
    }
    for (i = 0; i < itemCount; i++)
    {
        /* Calculating denominator for vi values (for each item) */
        distSums[i] = 0;
        for (j = 0; j < itemCount; j++)
            distSums[i] += distanceMatrix[i][j];
    }
    for (i = 0; i < itemCount; i++)
    {
        v[i] = 0;
        vIndexes[i] = i;    /* Keep item index for sorting */
        /* Calculating vi values for each item */
        for (j = 0; j < itemCount; j++)
            v[i] += (double)(distanceMatrix[i][j]/distSums[j]);
    }
    /* Sort both arrays to insure valid vIndexes values */
    SelectionSortTwoArrays(v, vIndexes, itemCount);
    /* Select k items with smallest vi values */
    for (i = 0; i < k; i++)
        initMedoids[i] = vIndexes[i];
    
    free(v);
    free(vIndexes);
    free(distSums);
    
    return initMedoids;
}

void SelectionSortTwoArrays(double *array, int *array2, int arrayLength)
{
    int i, j;
    double temp;
    
    for (i = 0; i < arrayLength - 1; i++)
    {
        for (j = i+1; j < arrayLength; j++)
        {
            if (array[i] > array[j])
            {
                temp = array[i];
                array[i] = array[j];
                array[j] = temp;
                
                temp = array2[i];
                array2[i] = array2[j];
                array2[j] = temp;
                
            }
        }
    }
}

int* assignmentPAM(double **distanceMatrix, int *medoids, int itemCount, int k, double *totalCost)
{
    int i, j;
    int *clustersArray;
    double minDistance;
    
    clustersArray = malloc(sizeof(int)*itemCount);
    if (clustersArray == NULL)
    {
        printf("ClusterFunctions.c: assignmentPAM: Failed to allocate memory for clustersArray\n");
        exit(-1);
    }
    
    *totalCost = 0;
    for (i = 0; i < itemCount; i++) /* For each item */
    {
        /* Choose a large number as min distance */
        minDistance = (double)RAND_MAX;
        for (j = 0; j < k; j++) /* For each medoid */
        {
            /* Find medoid with min distance from item */
            if (minDistance > distanceMatrix[i][medoids[j]])
            {
                minDistance = distanceMatrix[i][medoids[j]];
                clustersArray[i] = j;   /* Item will go to medoids j cluster */
            }
        }
        *totalCost += minDistance; /* Calculate total cost */
    }
    return clustersArray;
}

short** updateClarans(double **distanceMatrix, int *clustersArray, int *medoids, int itemCount, int k, double *totalCost, int numLocal, int maxNeighbor, int L, char fileType)
{
    int i, j, m, *currentBestMedoids, *nonMedoids, temp, *assignm, whichChanges, whichGets, isMedoid, neighbourCounter, curUserMedoid, myCluster, *allNeighbours;
    short **pNearest;
    double optimalCost, localCost, currentCost, radius, initRadius, vicinity, *vicinities;
    
    if (fileType == 'v')
        initRadius = 3.74;
    else
        initRadius = 0.007;
    radius = initRadius;
    
    currentBestMedoids = malloc(sizeof(int)*k);
    if (currentBestMedoids == NULL)
    {
        perror("ClusterFunctions.c: updateClarans: malloc failed");
        exit(-1);
    }
    nonMedoids = malloc(sizeof(int)*(itemCount - k));
    if (nonMedoids == NULL)
    {
        perror("ClusterFunctions.c: updateClarans: malloc failed");
        exit(-1);
    }
    
    optimalCost = *totalCost;
    /* Creating an array that contains all non medoid items */
    for (i = 0, j = 0; i < itemCount; i++)
    {
        isMedoid = 0;
        for (m = 0; m < k; m++)
        {
            if (medoids[m] == i)
            {
                isMedoid = 1;
                break;
            }
        }
        
        if (isMedoid == 0)
        {
            nonMedoids[j] = i;
            j++;
        }
    }
    
    for (i = 0; i < numLocal; i++)
    {
        /* Randomly select k new medoids */
        currentBestMedoids =  selectKBits(itemCount, k);
        /* Make the assignment to update (local) cost */
        assignm = assignmentPAM(distanceMatrix, currentBestMedoids, itemCount, k, &localCost);
        //free(assignm);
        
        for (j = 0; j < maxNeighbor; j++)
        {
            /* Randomly select a medoid to change */
            whichChanges = rand()%k;
            /* Randomly select a non medoid item */
            whichGets = rand()%(itemCount - k);
            /* Interchange medoid with non medoid */
            temp = currentBestMedoids[whichChanges];
            currentBestMedoids[whichChanges] = nonMedoids[whichGets];
            nonMedoids[whichGets] = temp;
            /* Make the assignment to calculate new cost */
            
            assignm = assignmentPAM(distanceMatrix, currentBestMedoids, itemCount, k, &currentCost);
            //free(assignm);
            
            if (currentCost < localCost)
            {
                /* If newis lower than previous cost, restart */
                j = -1;
                localCost = currentCost;
            }
            else
            {
                /* Otherwise, revert back */
                temp = currentBestMedoids[whichChanges];
                currentBestMedoids[whichChanges] = nonMedoids[whichGets];
                nonMedoids[whichGets] = temp;
            }
        }
        
        /* If optimal cost is greater than local cost */
        if (optimalCost > localCost)
        {
            /* Update optimal cost */
            optimalCost = localCost;
            for (j = 0; j < k; j++)     /* Update medoids to new configuration */
            {
                medoids[j] = currentBestMedoids[j];
            }
            assignm = assignmentPAM(distanceMatrix, currentBestMedoids, itemCount, k, &localCost);
        }
    }
    clustersArray = assignmentPAM(distanceMatrix, medoids, itemCount, k, totalCost);
    
    pNearest = malloc(sizeof(short*)*itemCount);
    if (pNearest == NULL)
    {
        perror("ClusterFunctions.c: updateClarans: malloc failed");
        exit(-1);
    }
    for (i = 0; i < itemCount; i++)
    {
        pNearest[i] = malloc(sizeof(short)*P);
        if (pNearest[i] == NULL)
        {
            perror("ClusterFunctions.c: updateClarans: malloc failed");
            exit(-1);
        }
        for (j = 0; j < P; j++)
            pNearest[i][j] = -1.0;
        
        /* we find i-th user's P nearest neighbours from their cluster */
		curUserMedoid = clustersArray[i];
		myCluster=0;
		for (j = 0; j < itemCount; j++) /* Count user's cluster size */
		{
			if ((curUserMedoid == clustersArray[j]) && (i != j))
				myCluster++;
		}
		allNeighbours = malloc(sizeof(int)*myCluster);
        if (allNeighbours == NULL)
        {
            perror("ClusterFunctions.c: updateClarans: malloc failed");
            exit(-1);
        }
		vicinities = malloc(sizeof(double)*myCluster);
        if (vicinities == NULL)
        {
            perror("ClusterFunctions.c: updateClarans: malloc failed");
            exit(-1);
        }
				
		m = 0;
		for (j = 0; j < itemCount; j++)
		{
            /* if users are in same medoid */
			if ((curUserMedoid == clustersArray[j]) && (i != j))
			{
				allNeighbours[m] = j;
				vicinities[m] = distanceMatrix[i][j];
				m++;
			}
		}
		SelectionSortTwoArrays(vicinities, allNeighbours, myCluster); /* Sort both arrays to find user's P NN */
		for (j = 0; (j < myCluster) && (j < P); j++)
			pNearest[i][j] = allNeighbours[j];
		free(allNeighbours);			free(vicinities);
    }
    
    return pNearest;
}

void medoidDistance(double **distanceMatrix, double *minDistance, double *maxDistance, int *medoids, int k)
{
    int i, j;
    double tempDistance;
    
    *minDistance = (double)RAND_MAX;    /* Set min distance to a large number */
    *maxDistance = 0.0;
    
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < k; j++)
        {
            /* Do not calculate distance from itself */
            if (i == j)
                continue;
            
            tempDistance = distanceMatrix[medoids[i]][medoids[j]];
            if (tempDistance > *maxDistance)
                *maxDistance = tempDistance;
            else if (tempDistance < *minDistance)
                *minDistance = tempDistance;
            else
                continue;
        }
    }
}

double* silhouette(double **distanceMatrix, int *clustersArray, int *medoids, int itemCount, int k)
{
    int i, j, m, searchCluster, counter;
    double *a, *b, *s, sum, medoidDistance;
    
    a = malloc(sizeof(double)*itemCount);
    if (a == NULL)
    {
        perror("ClusterFunctions.c: silhouette: malloc failed");
        exit(-1);
    }
    b = malloc(sizeof(double)*itemCount);
    if (b == NULL)
    {
        perror("ClusterFunctions.c: silhouette: malloc failed");
        exit(-1);
    }
    s = malloc(sizeof(double)*(k+1));
    if (s == NULL)
    {
        perror("ClusterFunctions.c: silhouette: malloc failed");
        exit(-1);
    }
    for (i = 0; i < (k+1); i++)
        s[i] = 0.0;
    for (i = 0; i < itemCount; i++) /* For each item */
    {
        counter = 0;
        sum = 0.0;
        /* For each item, calculate sum of distances from other items in its cluster */
        for (j = 0; j <itemCount; j++)
        {
            if (clustersArray[i] == clustersArray[j])
            {
                counter++;
                sum += distanceMatrix[i][j];
            }
        }
        
        a[i] = sum/counter; /* Calculate mean distance from items in same cluster */
        counter = 0;
        sum = 0.0;
        medoidDistance = (double)RAND_MAX; /* Set min distance to neighbour medoids to large number */
        /* For each medoid, find medoid with min distance from item */
        for (m = 0; m < k; m++)
        {
            if (clustersArray[i] == m)  /* If medoid from same cluster, continue */
                continue;
            if ((medoidDistance > distanceMatrix[i][medoids[m]]) && (distanceMatrix[i][medoids[m]] != 0))
            {
                medoidDistance = distanceMatrix[i][medoids[m]];
                searchCluster = m;
            }
        }
        /* Calculate distance from items of searchCluster (cluster with min distance from item) */
        for (j = 0; j < itemCount; j++)
        {
            if (clustersArray[j] == searchCluster)
            {
                counter++;
                sum += distanceMatrix[i][j];
            }
        }
        b[i] = sum/counter; /* Calculate mean distance from items of neighbor cluster */
        
        if (b[i] > a[i])
            sum = (b[i] - a[i])/b[i];
        else
            sum = (b[i] - a[i])/a[i];
        s[i] = sum/counter;
        s[k] += sum;
    }
    
    s[k] /= itemCount;
    
    free(a);
    free(b);
    
    return s;
}