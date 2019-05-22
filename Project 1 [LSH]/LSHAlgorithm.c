/* LSH Algorithm.c */
#include "LSHAlgorithm.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include "ADT.h"
#include "LSHfunctions.h"

#define c 2

void LSHAlgorithm(int argc, const char * argv[])
{
    listNodePtr** hashTables, listStart;
    item temp;
    int **kDiffNums, i, itemIndex, N, hammingEOF, hashTableSize, j, m, vectorEOF, *hFunctions, *disc, pin, *duplicateCheck, itemCounter;
    int k, L;
    struct timeval time1, time2;
    double distance, radius, ***hx1x2, bestNNDistance, elapsedTime;
    vectorHash* vHash;
    dMatrixItem dmi;
    dMatrixHash* dmh;
    char tempChar[20], metric, typeOfFile, IfileName[30], QfileName[30], OfileName[30], *bestNNname;
    int newL, newK;
    FILE *queryFile;                /* Query File */
    FILE *inFile;                   /* Dataset file */
    FILE *outputFile;               /* Output file that contains results */
   
    /* Default values for k, L */
    k = 4;
    L = 5;
    
    /* If no arguments passed, ask user */
    if (argc <= 1)
    {
        printf("Please enter input file name: ");
        scanf("%s", IfileName);
        printf("Please enter query file name: ");
        scanf("%s", QfileName);
        printf("Please enter output file name: ");
        scanf("%s", OfileName);
        printf("Would you like to enter a new value for L? (1 = yes, 0 = no)");
        scanf("%d", &newL);
        if (newL == 1)
        {
            printf("Enter new L: ");
            scanf("%d", &L);
        }
        printf("Would you like to enter a new value for k? (1 = yes, 0 = no)");
        scanf("%d", &newK);
        if (newK == 1)
        {
            printf("Enter new k: ");
            scanf("%d", &k);
        }
    }
    else
    {
        for (i = 1; i < argc; i += 2)
        {
            if (argv[i][1] == 'L')
                L = atoi(argv[i+1]);
            else if (argv[i][1] == 'k')
                k = atoi(argv[i+1]);
            else if (argv[i][1] == 'd')
                strcpy(IfileName, argv[i+1]);
            else if (argv[i][1] == 'q')
                strcpy(QfileName, argv[i+1]);
            else if (argv[i][1] == 'o')
                strcpy(OfileName, argv[i+1]);
            else
            {
                printf("Warning! Wrong options passed. The program will now be terminated.\n");
                exit(-1);
            }
        }
    }
    
    inFile = fopen(IfileName, "r");
    if (inFile == NULL)
    {
        printf("LSHAlgorithm.c: failed to open file (inputFile) %s\n", IfileName);
        exit(-1);
    }
    srand((unsigned int)time(NULL));
    
    /* Read first line from input file */
    fscanf(inFile, "%s", tempChar);
    fscanf(inFile, "%s", tempChar);
    
    if (strcmp(tempChar, "hamming") == 0)
    {
        typeOfFile = 'h';
        /* Get hamming space's dimension */
        N = getHammingDimension(IfileName);
        /* Create L arrays of k-randomly selected bits */
        kDiffNums = malloc(L*sizeof(int*));
        if (kDiffNums == NULL)
        {
            printf("LSHAlgorithm.c: failed to allocate memory for kDiffNums\n");
            exit(-1);
        }
        /* For each of the L functions g, select the k-bits of hamming item to create each hash function g */
        for(i = 0; i < L; i++)
        {
            kDiffNums[i] = selectKBits(N, k);
        }
    }
    else if (strcmp(tempChar, "vector") == 0 || strcmp(tempChar, "euclidean") == 0)
    {
        typeOfFile = 'v';
        N = getVectorDimension(IfileName, &metric);
        fscanf(inFile, "%s", tempChar);
        fscanf(inFile, "%s", tempChar);
        
        /* Allocate memory for L vector hash structures, containing randomly selected v vectors in [0,1], t randomly selected in [0, w) and randomly selected ri integers for hash function g */
        vHash = malloc(L*sizeof(vectorHash));
        if (vHash == NULL)
        {
            printf("LSHAlgorithm.c: failed to allocate memory for vHash\n");
            exit(-1);
        }
        /* Creating L vector hash structures */
        for (i = 0; i < L; i++)
        {
            createVectorHash(&(vHash[i]), k, N);
        }
    }
    else /* distance matrix */
    {
        typeOfFile = 'm';
        N = getDistanceMatrixDimension(IfileName);      /* Get dimension of distance matrix */
        readDistances(inFile, N, &dmi);         /* Allocate 2D, NxN array that holds the distances between the items */
        
        /* Create L dMatrixHash structures, each holds data associated with only 1 hash table */
        dmh = malloc(L*sizeof(dMatrixHash));
        if (dmh == NULL)
        {
            printf("LSHAlgorithm.c: failed to allocate memory for dmh\n");
            exit(-1);
        }
        
        /* Create L hx1x2 2d arrays */
        hx1x2 = malloc(L*sizeof(double**));
        if (hx1x2 == NULL)
        {
            printf("LSHAlgorithm: failed to allocate memory for hx1x2\n");
            exit(-1);
        }
        
        /* For L hash functions g */
        for (i = 0; i < L; i++)
        {
            createDMatrixHash(&(dmh[i]), k, N);          /* Creates structure and selects k random x1, x2 items */
            hx1x2[i] = dbhFucntionH(dmh[i].x1, dmh[i].x2, dmi.distances, N, k);      /* Computes k-h[x1x2] for all items */
            selectt1t2(&dmh[i], hx1x2[i], N, k);          /* Randomly selects k pairs of t1, t2 */
        }
        
    }
    
    hashTableSize = (int)(0.5 + pow(2, k));		/* Creating L hash tables */
    hashTables = createHashTables(k, L);
    
    temp.indexNo = -1;
    
    i = 0;
    for (;;)
    {
        if (typeOfFile == 'h')
        {
            /* Read each item from file until end of file */
            hammingEOF = readHamming(inFile, N, &temp);
            if (hammingEOF == 0)
                break;
            
            temp.indexNo++;
            /* Hash it through L hashing functions */
            for(i = 0; i < L; i++)
            {
                itemIndex = hashFunctionG(temp, kDiffNums[i], k);           /* Get hash table bucket number/index */
                insertFirst(&(hashTables[i][itemIndex]), temp);             /* Insert into hash table i */
            }
        }
        else if (typeOfFile == 'v')
        {
            /* Read each item from file until end of file */
            vectorEOF = readVector(inFile, N, &temp);
            if (vectorEOF == 0)
                break;
            
            temp.indexNo++;
            /* Hash it through L hashing functions */
            for (i = 0; i < L; i++)
            {
                hFunctions = euclideanHashFunctionH(vHash[i], temp, N, k);          /* Pass item through k h functions */
                itemIndex = euclideanHashFunctionF(vHash[i], k, hFunctions);        /* Get hash table bucket number/index */
                insertFirst(&(hashTables[i][itemIndex]), temp);                     /* Insert into hash table i */
            }
        }
        else /* distance matrix */
        {
            
            if (i == N)     /* read all distances from input file */
                break;
            
            for (j = 0; j < L; j++)
            {
                indexItems(&temp, dmi.distances, dmi.dmItemName, i);        /* Index each item */
                disc = discritizeH(dmh[j], hx1x2[j], N, k, i);              /* Discritize each of the k h functions */
                itemIndex = DBHhashFunctionG(disc, k);                      /* Get hash table bucket number/index */
                insertFirst(&(hashTables[j][itemIndex]), temp);             /* Insert into hash table i */
            }
            i++;
        }
    }
    fclose(inFile);
    itemCounter = temp.indexNo + 1;
    
    /* QUERY ITEMS INPUT */
    
    queryFile = fopen(QfileName, "r");
    if (queryFile == NULL)
    {
        printf("LSHAlgorithm.c: failed to open file (queryFile) %s", QfileName);
        exit(-1);
    }
    
    outputFile = fopen(OfileName, "w");
    if (outputFile == NULL)
    {
        printf("LSHAlgorithm.c: failed to create and open file (outputFile) LSHoutput.txt\n");
        exit(-1);
    }
    /* Read first line from query file and get radius */
    fscanf(queryFile, "%s", tempChar);
    fscanf(queryFile, "%s", tempChar);
    radius = atof(tempChar);
    
    bestNNname = malloc(sizeof(char)*20);       /* Allocate memory for Best NN name */
    if (bestNNname == NULL)
    {
        printf("LSHAlgorithm.c: failed to allocate memory for bestNNname\n");
        exit(-1);
    }
    duplicateCheck = malloc(itemCounter*sizeof(int));
    if (duplicateCheck == NULL)
    {
        printf("LSHAlgorithm.c: failed to allocate memory for duplicateCheck\n");
        exit(-1);
    }
    
    for (;;)
    {
        gettimeofday(&time1, NULL);             /* Start timer */
        for (j = 0; j < itemCounter; j++)       /* Check for duplicates while writting to output file */
            duplicateCheck[j] = -1;
        
        if (typeOfFile == 'h')
        {
            hammingEOF = readHamming(inFile, N, &temp);
            if (hammingEOF == 0)
                break;
            
            fprintf(outputFile, "Query: %s\n", temp.itemName);
            fprintf(outputFile, "R-near neighbors:\n");
            bestNNDistance = RAND_MAX;      /* Set best distance as largest number possible, instead of infinity */
            strcpy(bestNNname, " ");        /* Set blank name */
            j = 0;
            
            for (i = 0; i < L; i++)
            {
                itemIndex = hashFunctionG(temp, kDiffNums[i], k);
                listStart = hashTables[i][itemIndex];
                
                while (listStart != NULL)
                {
                    distance = hammingDistance(listStart->hashedItem, temp, N);
                    if (distance < c*radius)       /* Check if distance suffices and also if item is duplicate (already printed) */
                    {
                        m = 0;
                        while (duplicateCheck[m] != (listStart->hashedItem).indexNo)
                        {
                            m++;
                            if (duplicateCheck[m] == -1)
                            {
                                duplicateCheck[m] = (listStart->hashedItem).indexNo;
                                fprintf(outputFile, "%s\n", (listStart->hashedItem).itemName);
                                break;
                            }
                        }
                        
                        if (distance < bestNNDistance)
                        {
                            strcpy(bestNNname, (listStart->hashedItem).itemName);
                            bestNNDistance = distance;
                        }
                    }
                    listStart = listStart->next;
                }
            }
            fprintf(outputFile, "Nearest neighbor: %s\n", bestNNname);
            fprintf(outputFile, "distanceLSH: %lf\n", bestNNDistance);
        }
        else if (typeOfFile == 'v')
        {
            /* Reading query items from query file */
            vectorEOF = readVector(inFile, N, &temp);
            if (vectorEOF == 0)
                break;
            
            fprintf(outputFile, "Query: %s\n", temp.itemName);
            fprintf(outputFile, "R-near neighbors:\n");
            bestNNDistance = RAND_MAX;
            strcpy(bestNNname, " ");
            
            for (i = 0; i < L; i++)
            {
                hFunctions = euclideanHashFunctionH(vHash[i], temp, N, k);      /* Pass item through k h functions */
                itemIndex = euclideanHashFunctionF(vHash[i], k, hFunctions);    /* Get hash table bucket number/index */
                listStart = hashTables[i][itemIndex];
                /* Start searching in bucket for nearest neighbors */
                while (listStart != NULL)
                {
                    distance = euclideanDistance(listStart->hashedItem, temp, N);
                    if (distance < c*radius)    /* Check if distance suffices and also if item is duplicate (already printed) */
                    {
                        m = 0;
                        while (duplicateCheck[m] != (listStart->hashedItem).indexNo)
                        {
                            m++;
                            if (duplicateCheck[m] == -1)
                            {
                                duplicateCheck[m] = (listStart->hashedItem).indexNo;
                                fprintf(outputFile, "%s\n", (listStart->hashedItem).itemName);
                                break;
                            }
                        }
                        
                        if (distance < bestNNDistance)  /* For approximate nearest neighbor */
                        {
                            strcpy(bestNNname, (listStart->hashedItem).itemName);
                            bestNNDistance = distance;
                        }
                    }
                    listStart = listStart->next;
                }
            }
            fprintf(outputFile, "Nearest neighbor: %s\n", bestNNname);
            fprintf(outputFile, "distanceLSH: %lf\n", bestNNDistance);
        }
        else /* matrix */
        {
            if (feof(queryFile))
                break;
            fscanf(queryFile, "%s", tempChar);
            
            /* Read query item from file and insert into dmi.distances[i][j]'s last index position */
            for (i = 0; i < N; i++)
                fscanf(queryFile, "%lf", &(dmi.distances[N][i]));
            
            fprintf(outputFile, "Query: %s\n", tempChar);
            fprintf(outputFile, "R-near neighbors:\n");
            bestNNDistance = RAND_MAX;
            strcpy(bestNNname, " ");
            
            for (i = 0; i < L; i++)
            {
                itemIndex = queryHashingDBH(dmi.distances, N, dmh[i], k, L);        /* Hash query item */
                listStart = hashTables[i][itemIndex];
                /* Searching for nearest neighbors */
                while (listStart != NULL)
                {
                    distance = dmi.distances[N][(listStart->hashedItem).indexNo];
                    if (distance < c*radius)    /* Check if distance suffices and also if item is duplicate (already printed) */
                    {
                        m = 0;
                        while (duplicateCheck[m] != (listStart->hashedItem).indexNo)
                        {
                            m++;
                            if (duplicateCheck[m] == -1)
                            {
                                duplicateCheck[m] = (listStart->hashedItem).indexNo;
                                fprintf(outputFile, "%s\n", (listStart->hashedItem).itemName);
                                break;
                            }
                        }
                        
                        if (distance < bestNNDistance)      /* For approximate nearest neighbor */
                        {
                            strcpy(bestNNname, (listStart->hashedItem).itemName);
                            bestNNDistance = distance;
                        }
                    }
                    listStart = listStart->next;
                }
            }
            fprintf(outputFile, "Nearest neighbor: %s\n", bestNNname);
            fprintf(outputFile, "distanceLSH: %lf\n", bestNNDistance);
        }
        
        gettimeofday(&time2, NULL);     /* Stop timer */
        elapsedTime = (time2.tv_sec - time1.tv_sec) * 1000.0;      /* sec to ms */
        elapsedTime += (time2.tv_usec - time1.tv_usec) / 1000.0;   /* us to ms */
        fprintf(outputFile, "tLSH: %g", elapsedTime);
        
        /* Exhaustive search for Nearest Neighbor */
        bestNNDistance = RAND_MAX;      /* set best NN distance to max possible */

        gettimeofday(&time1, NULL);
        for (i = 0; i < hashTableSize; i++)
        {
            listStart = hashTables[0][i];
            
            while (listStart != NULL)
            {
                if(typeOfFile == 'h')
                {
                    distance = hammingDistance(listStart->hashedItem, temp, N);
                    if (distance < c*radius)
                    {
                        if (distance < bestNNDistance)
                        {
                            strcpy(bestNNname, (listStart->hashedItem).itemName);
                            bestNNDistance = distance;
                        }
                    }
                }
                else if(typeOfFile == 'v')
                {
                    distance = euclideanDistance(listStart->hashedItem, temp, N);
                    if (distance < c*radius)
                    {
                        if (distance < bestNNDistance)
                        {
                            strcpy(bestNNname, (listStart->hashedItem).itemName);
                            bestNNDistance = distance;
                        }
                    }
                }
                else /* matrix */
                {
                    distance = dmi.distances[N][(listStart->hashedItem).indexNo];
                    if (distance < c*radius)
                    {
                        if (distance < bestNNDistance)
                        {
                            strcpy(bestNNname, (listStart->hashedItem).itemName);
                            bestNNDistance = distance;
                        }
                    }
                }
                
                listStart = listStart->next;
            }
        }
        gettimeofday(&time2, NULL);         /* Stop timer */
        elapsedTime = (time2.tv_sec - time1.tv_sec) * 1000.0;
        elapsedTime += (time2.tv_usec - time1.tv_usec) / 1000.0;
        fprintf(outputFile, "\ntTRUE: %g", elapsedTime);
        fprintf(outputFile, "\ndistanceTRUE: %f\n\n",bestNNDistance);
    }
    fclose(queryFile);
    fclose(outputFile);
    
    /* Deallocate memory used for any of the 3 files types */
    if (typeOfFile == 'h')      /* Deallocate memory used for Hamming space hashing parameters */
    {
        for (i = 0; i < L; i++)
        {
            free(kDiffNums[i]);
            kDiffNums[i] = NULL;
        }
        free(kDiffNums);
        kDiffNums = NULL;
    }
    else if (typeOfFile == 'v')     /* Deallocate memory used for Euclidean space hashing parameters */
    {
        free(hFunctions);
        hFunctions = NULL;
        
        destroyVectorHash(vHash, k);
    }
    else            /* Deallocate memory used for DBH parameters */
    {
        destroyDMatrixHash(dmh, k);
        free(disc);
        disc = NULL;
        
        for (i = 0; i < L; i++)
        {
            for (j = 0; j < N; j++)
            {
                free(hx1x2[i]);
                hx1x2[i] = NULL;
            }
        }
        free(hx1x2);
        hx1x2 = NULL;
        
        for (i = 0; i < N; i++)
        {
            free(dmi.distances[i]);
            dmi.distances[i] = NULL;
            free(dmi.dmItemName[i]);
            dmi.dmItemName[i] = NULL;
        }
        free(dmi.distances);
        dmi.distances = NULL;
        free(dmi.dmItemName);
        dmi.dmItemName = NULL;
    }
    
    /* Deallocate memory used for LSH algorithm, independent of file type */
    for (i = 0; i < L; i++)
    {
        for (j = 0; j < hashTableSize; j++)
            destroyList(&hashTables[i][j]);
        
        free(hashTables[i]);
        hashTables[i] = NULL;
    }
    free(hashTables);
    hashTables = NULL;
    
    free(duplicateCheck);
    duplicateCheck = NULL;
    
    free(bestNNname);
    bestNNname = NULL;
}