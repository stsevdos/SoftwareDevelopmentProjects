#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ADT.h"
#include "Recommendation.h"
#include "LSHFunctions.h"

extern int P;

int* selectKBits(int N, int k)
{
    int i, j, sameExists;
    
    int *indexes = malloc(k*sizeof(int));
    if (indexes == NULL)
    {
        perror("LSHfunctions.c: selectKBits: malloc failed");
        exit(-1);
    }
    
    indexes[0] = rand() % N;            /* Initializes first bit outside of loop so it can compare two of the for equality inside loop */
    for (i = 1; i < k;)
    {
        indexes[i] = rand() % N;
        
        sameExists = 0;
        for (j = 0; j < i; j++)         /* Checking if array indexes selected up to i have been selected before and if so randomly selects again */
        {
            if (indexes[i] == indexes[j])
                sameExists = 1;
        }
        if (sameExists == 0)           /* if same array index doesn't exist in already randomly selected k-bts, continue with next random bit */
            i++;
    }
    return indexes;                     /* returns the array of k-randomly selected bits (in array indexes, not actual hamming bits) to hash point in bucket */
}

double euclideanDistance(int user1, int user2, short **ratings, int dimension)
{
    int i;
    double distance = 0;
    
    for (i = 0; i < dimension; i++)
    {
        distance += (ratings[user1][i] - ratings[user2][i])*(ratings[user1][i] - ratings[user2][i]);
    }
    
    distance = sqrt(distance);
    
    return distance;
}

/* Returns array of all h(p) functions for a single point in space */
int* euclideanHashFunctionH(vectorHash vh, short *user, int dimension, int k)
{
    int i, j;
    
    int *hFunctions = malloc(sizeof(int)*k);
    if (hFunctions == NULL)
    {
        perror("LSHFunctions.c: euclideanHashFunctionH: malloc failed");
        exit(-1);
    }
    
    /* Calculating k-hi functions for one g hash function */
    for (i = 0; i < k; i++)
    {
        hFunctions[i] = 0;
        for (j = 0; j < dimension; j++)
        {
            if (vh.v[i][j] == 1)
                hFunctions[i] += user[j];
        }
        
        hFunctions[i] = (hFunctions[i] + vh.t)/W;
    }
    
    return hFunctions;
}

/* Returns index of hash table for a single point in space */
int euclideanHashFunctionF(vectorHash vh, int k, int *hFunctions)
{
    int i, f = 0, tableSize;
    
    long M;
    
    M = pow(2, 31) - 5;
    
    tableSize = 0.5 + pow(2, k);
    
    for (i = 0; i < k; i++)
    {
        f += (hFunctions[i]*vh.r[i]) % M;
    }
    
    f = f % tableSize;
    
    return f;
}

double cosineSimilarity(int user1, int user2, short **ratings, int dimension)
{
    int i;
    double similarity, sum1, sum2;
    
    similarity = 0; sum1 = 0;   sum2 = 0;

    for (i = 0; i < dimension; i++)
    {
        similarity += (ratings[user1][i]*ratings[user2][i]);
        sum1 += ratings[user1][i]*ratings[user1][i];
        sum2 += ratings[user2][i]*ratings[user2][i];
    }
    
    sum1 = sqrt(sum1);
    sum2 = sqrt(sum2);
    
    similarity /= sum1*sum2;
    if (similarity != 0.0)
        similarity = 1.0 - (similarity + 1.0)/2.0;
    else
        similarity = 1.0;
    
    return similarity;
}

short** cosineRandomVectors(int k, int dimension)
{
    int i, j;
    short **csVectors;
    
    csVectors = malloc(sizeof(short*)*k);
    if (csVectors == NULL)
    {
        perror("LSHFunctions.c: cosineRandomVectors: malloc failed");
        exit(-1);
    }
    
    for (i = 0; i < k; i++)
    {
        csVectors[i] = malloc(sizeof(short)*dimension);
        if (csVectors[i] == NULL)
        {
            perror("LSHFunctions.c: cosineRandomVectors: malloc failed");
            exit(-1);
        }
        
        for (j = 0; j < dimension; j++)
            csVectors[i][j] = rand()%2;
    }
    
    return csVectors;
}

int cosineHashFunction(int k, int dimension, short **csVectors, short *ratings, float *userMean)
{
    int *h, i, j, index;
    float dotProduct;
    
    h = malloc(sizeof(int)*k);
    if (h == NULL)
    {
        perror("LSHFunctions.c: cosineHashFunction: malloc failed");
        exit(-1);
    }
    for (i = 0; i < k; i++)
    {
        h[i] = 0;
        dotProduct = 0;
        for (j = 0; j < dimension; j++)
        {
            if ((ratings[j]>0)&&(csVectors[i][j] == 1))
                dotProduct += csVectors[i][j]*(ratings[j] - userMean[j]);
        }
        if (dotProduct >= 0)
            h[i] = 1;
        else
            h[i] = 0;
    }
    
    index = 0;
    for (i = 0; i < k; i++)
        index = (index * 2) + h[i];
    
    return index;
}

double initialRadius(short **ratings, int userCount, int musicCount, int userID, char metric)
{
	int i;
	double radius, temp;
    
	radius = 5000.0;
	for (i = 0; i < userCount; i++)
	{
		if(i==userID)
			continue;
			
		if (metric == 'v')
		{
			if (i == userID)
				continue;
					
			temp = euclideanDistance(userID, i, ratings, musicCount);
			if ((temp < radius) && (temp != 0.0))
				radius = temp;
		}
		else
		{
			if (i == userID)
				continue;
						
			temp = cosineSimilarity(userID, i, ratings, musicCount);
			if ((temp < radius) && (temp != 0.0))
				radius = temp;	
		}
	}

    return radius;
}