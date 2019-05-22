#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "LSHAlgorithm.h"

int main(int argc, const char * argv[])
{
    int userSelection;
    
    /* Run LSH algorithm for first time, with arguments passed to program */
    LSHAlgorithm(argc, argv);
    /* Ask user if they wish to perform another search */
    printf("Would you like to perform a new search with another set of files? (1 = yes, 0 = no)\n Selection: ");
    scanf("%d", &userSelection);
    /* While user says yes, run LSH algorithm, pass null arguments */
    while (userSelection == 1)
    {
        LSHAlgorithm(0, NULL);
        printf("Would you like to perform a new search with another set of files? (1 = yes, 0 = no)\n Selection: ");
        scanf("%d", &userSelection);
    }
    
    return 0;
}
