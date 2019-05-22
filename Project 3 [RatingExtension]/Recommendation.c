#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ADT.h"
#include "Recommendation.h"
#include "LSHFunctions.h"
#include "LSHAlgorithm.h"

extern int P;

int getUserCount(char *fileName)
{
    int userCount, lastUser;
    char *stringBuffer;
    FILE *inputFile;
    
    inputFile = fopen(fileName, "r");
    if (inputFile == NULL)
    {
        perror("Recommendation.c: getUserCount: fopen failed");
        exit(-1);
    }
    stringBuffer = malloc(sizeof(char)*50);
    if (stringBuffer == NULL)
    {
        perror("Recommendation.c: getUserCount: malloc failed");
        exit(-1);
    }
    
    /* Read first line from file (P: <int>) */
    fscanf(inputFile, "%s", stringBuffer);
    fscanf(inputFile, "%s", stringBuffer);
    P = atoi(stringBuffer);
    
    userCount = 0;
    while (!feof(inputFile))
    {
        /* Read user ID */
        fscanf(inputFile, "%s", stringBuffer);
        if (userCount > 0) /* If a user ID has been read before */
        {
            if (atoi(stringBuffer) != lastUser) /* And we read a different user ID */
            {
                userCount++;    /* increment users found */
                lastUser = userCount;
            }
        }
        else    /* if it's the first user ID we read from file */
        {
            userCount++;
            lastUser = userCount;
        }
        
        /* Read movie ID and movie rating from same line */
        fscanf(inputFile, "%s", stringBuffer);  /* Read movie ID */
        fscanf(inputFile, "%s", stringBuffer);  /* Read user rating */
    }
    
    userCount -= 1 ;
    
    fclose(inputFile);
    free(stringBuffer);
    
    return userCount;
}

int getMusicCount(char *fileName)
{
    int musicCount;
    char *stringBuffer;
    FILE *inputFile;
    
    inputFile = fopen(fileName, "r");
    if (inputFile == NULL)
    {
        perror("Recommendation.c: getMusicCount: fopen failed");
        exit(-1);
    }
    stringBuffer = malloc(sizeof(char)*50);
    if (stringBuffer == NULL)
    {
        perror("Recommendation.c: getMusicCount: malloc failed");
        exit(-1);
    }
    
    /* Read first line from file (P: <int>) */
    fscanf(inputFile, "%s", stringBuffer);
    fscanf(inputFile, "%s", stringBuffer);
    
    musicCount = -1;
    while (!feof(inputFile))
    {
        /* Read user ID */
        fscanf(inputFile, "%s", stringBuffer);
        /* Read music ID */
        fscanf(inputFile, "%s", stringBuffer);
        if (atoi(stringBuffer) > musicCount)
            musicCount = atoi(stringBuffer);    /* If new ID read is greater than last, set new ID as max */
        
        /* Read user rating */
        fscanf(inputFile, "%s", stringBuffer);
    }
    
    fclose(inputFile);
    free(stringBuffer);
    return musicCount;
}

void getMusicRatings(char *fileName, short **ratings)
{
    short userID, musicID, rating;
    char *stringBuffer;
    FILE *inputFile;
    
    inputFile = fopen(fileName, "r");
    if (inputFile == NULL)
    {
        perror("Recommendation.c: getMusicRatings: fopen failed");
        exit(-1);
    }
    stringBuffer = malloc(sizeof(char)*50);
    if (stringBuffer == NULL)
    {
        perror("Recommendation.c: getMusicRatings: malloc failed");
        exit(-1);
    }
    
    /* Read first line from file */
    fscanf(inputFile, "%s", stringBuffer);
     fscanf(inputFile, "%s", stringBuffer);
    
    while (!feof(inputFile))
    {
        /* Read user ID */
        fscanf(inputFile, "%s", stringBuffer);
        userID = atoi(stringBuffer) - 1;
        /* Read music ID */
        fscanf(inputFile, "%s", stringBuffer);
        musicID = atoi(stringBuffer) - 1;
        /* Read rating */
        fscanf(inputFile, "%s", stringBuffer);
        rating = atoi(stringBuffer);
        ratings[userID][musicID] = rating;
    }
    
    fclose(inputFile);
    free(stringBuffer);
}

float* userRatingMean(short **ratings, int userCount, int musicCount)
{
    int i, j, meanCount;
    float *userMean;
    
    userMean = malloc(sizeof(float)*userCount);
    if (userMean == NULL)
    {
        perror("Recommendation.c: userRatingMean: malloc failed");
        exit(-1);
    }
    
    for (i = 0; i < userCount; i++)
    {
        meanCount = 0;
        userMean[i] = 0;
        for (j = 0; j < musicCount; j++)
        {
            /* if music has been rated by user */
            if (ratings[i][j] != 0)
            {
                userMean[i] += ratings[i][j];
                meanCount++;
            }
        }
        userMean[i] /= meanCount;
    }
    
    return userMean;
}

double** recommendItems(int userCount, int musicCount, short **pNearest, short **ratings, float *r, char fileType, FILE *outputFile)
{
    int i, j, m, max, maxPos;
    float **recommendedItems;
    double **pNNVicinity, z;
    
    /*** pNearest neighbors similarity/distance ***/
    pNNVicinity = malloc(sizeof(double*)*userCount);
    if (pNNVicinity == NULL)
    {
        perror("Recommendation.c: recommendItems: malloc failed");
        exit(-1);
    }
    
    for (i = 0; i < userCount; i++)
    {
        pNNVicinity[i] = malloc(sizeof(double)*P);
        if (pNNVicinity[i] == NULL)
        {
            perror("Recommendation.c: recommendItems: malloc failed");
            exit(-1);
        }
        
        for (j = 0; j < P; j++)
        {
            if (pNearest[i][j] != -1)   /* if nearest neighbour exists (user has P nearest neighbours) */
            {
                if (fileType == 'v')
                    pNNVicinity[i][j] = euclideanDistance(i, pNearest[i][j], ratings, musicCount);
                else
                    pNNVicinity[i][j] = cosineSimilarity(i, pNearest[i][j], ratings, musicCount);
            }
        }
    }
    
    /*** Recommended items ***/
    recommendedItems = malloc(sizeof(float*)*userCount);
    if (recommendedItems == NULL)
    {
        perror("Recommendation.c: recommendItems: malloc failed");
        exit(-1);
    }
    
    for (i = 0; i < userCount; i++)
    {
        recommendedItems[i] = malloc(sizeof(float)*musicCount);
        if (recommendedItems[i] == NULL)
        {
            perror("Recommendation.c: recommendItems: malloc failed");
            exit(-1);
        }

        for (j = 0; j < musicCount; j++)
            recommendedItems[i][j] = 0; /* Initialize all items to 0 */
    }
    
    for (i = 0; i < userCount; i++) /* for each user */
    {
        for (j = 0; j < musicCount; j++)    /* for each item */
        {
            if (ratings[i][j] > 0)  /* if user has rated item, item can't be recommended to them, flag with -10 */
            {
                recommendedItems[i][j] = -10;
                continue;
            }
            
            z = 0.0;
            for (m = 0; m < P; m++) /* for each of the P nearest neighbours */
            {
                if (pNearest[i][m] >= 0 && ratings[pNearest[i][m]][j])  /* Calculate rating for item recommended by nearest neighbour */
                {
                    z += 1 - pNNVicinity[i][m];
                    recommendedItems[i][j] += (1-pNNVicinity[i][m]) * (ratings[pNearest[i][m]][j] - r[pNearest[i][m]]);
                }
            }
        }
        recommendedItems[i][j] /= (float)z;
    }
    
    for (i = 0; i < userCount; i++) /* for each user */
    {
        fprintf(outputFile, "<U%d> : ",i+1);
        for (j = 0; j < 5; j++) /* Find 5 top recommended items */
        {
            max = -10;
            for (m = 0; m < musicCount; m++)
            {
                if (recommendedItems[i][m] > max)
                {
                    max = recommendedItems[i][m];
                    maxPos = m;
                }
            }
            if (max == -10)     /* If item has been rated by user, skip */
                continue;
            /* and write to file */
            fprintf(outputFile, " item%d", maxPos+1);
            recommendedItems[i][maxPos] = -10;      /* Mark item as recommended */
        }
        fprintf(outputFile, "\n");
    }
    
    /* Deallocating memory */
    for (i = 0; i < userCount; i++)
        free(recommendedItems[i]);
    free(recommendedItems);
    
    return pNNVicinity;  
}

double crossValidation(short **pNearest, short **ratings, float *r, int userCount, int musicCount, double** pNNdist)
{
	int i, j, m, nearest;
	float sumWeight, predicted;
	int commonItems, checkCount = 0;
	double Pj, S, relativeRating, MAEi, MAEtotal = 0.0;
	
	for (i = 0; i < userCount; i++)
	{
		S = 0.0; commonItems = 0; Pj = 0.0;
		for (j = 0; j < musicCount; j++)
		{
			sumWeight = 0.0;
            predicted = 0.0;
			if (ratings[i][j] == 0)
				continue;
			
			for (m = 0; m < P; m++)
			{	
				nearest = pNearest[i][m];
				if (nearest == -1)
					break;
					
				if (ratings[nearest][j] > 0)
				{
					sumWeight += (1.0-pNNdist[i][m]);
					predicted += (1.0-pNNdist[i][m])*(ratings[nearest][j]-r[j]); 
				}
			}
		
			if (sumWeight > 0.0)
			{
				commonItems++;
				relativeRating = predicted/sumWeight;
				S += fabs(ratings[i][j] - ( relativeRating + r[i]));
			}	
		}
	
		if (commonItems > 0)
		{
			checkCount++;	
			MAEi = S/commonItems;
			MAEtotal += MAEi; 
		}
		else
			MAEi = -10;
        }
	
	return (MAEtotal/checkCount);
}