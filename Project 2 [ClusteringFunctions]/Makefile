OBJS 	= main.o ADT.o LSHfunctions.o LSHAlgorithm.o ClusterFunctions.o
SOURCE	= main.c ADT.c LSHfunctions.c LSHAlgorithm.c ClusterFunctions.c
HEADER  = ADT.h LSHfunctions.h LSHAlgorithm.h ClusterFunctions.c
OUT  	= medoids
CC	= gcc
CFLAGS   = -g -c 

all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) -lm

main.o: main.c
	$(CC) $(CFLAGS) main.c

ADT.o: ADT.c
	$(CC) $(CFLAGS) ADT.c
	
LSHfunctions.o: LSHfunctions.c
	$(CC) $(CFLAGS) LSHfunctions.c

LSHAlgorithm.o: LSHAlgorithm.c
	$(CC) $(CFLAGS) LSHAlgorithm.c

ClusterFunctions.o: ClusterFunctions.c
	$(CC) $(CFLAGS) ClusterFunctions.c

clean:
	rm -f $(OBJS) $(OUT)

count:
	wc $(SOURCE) $(HEADER)