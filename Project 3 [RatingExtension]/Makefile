OBJS 	= main.o ADT.o LSHFunctions.o LSHAlgorithm.o ClusterFunctions.o	Recommendation.o
SOURCE	= main.c ADT.c LSHFunctions.c LSHAlgorithm.c ClusterFunctions.c	Recommendation.c
HEADER  = ADT.h LSHfunctions.h LSHAlgorithm.h ClusterFunctions.c Recommendation.h
OUT  	= recommendation
CC	= gcc
CFLAGS   = -g -c 

all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) -lm

main.o: main.c
	$(CC) $(CFLAGS) main.c

ADT.o: ADT.c
	$(CC) $(CFLAGS) ADT.c
	
LSHFunctions.o: LSHFunctions.c
	$(CC) $(CFLAGS) LSHFunctions.c

LSHAlgorithm.o: LSHAlgorithm.c
	$(CC) $(CFLAGS) LSHAlgorithm.c

ClusterFunctions.o: ClusterFunctions.c
	$(CC) $(CFLAGS) ClusterFunctions.c

Recommendation.o: Recommendation.c
	$(CC) $(CFLAGS)	Recommendation.c

clean:
	rm -f $(OBJS) $(OUT)

count:
	wc $(SOURCE) $(HEADER)