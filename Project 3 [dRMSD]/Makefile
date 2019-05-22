OBJS 	= main.o ClusterFunctions.o Bio.o
SOURCE	= main.c ClusterFunctions.c Bio.c
HEADER  = ClusterFunctions.c Bio.h
OUT  	= bio
CC	= gcc
CFLAGS   = -g -c 

all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) -lm

main.o: main.c
	$(CC) $(CFLAGS) main.c

ClusterFunctions.o: ClusterFunctions.c
	$(CC) $(CFLAGS) ClusterFunctions.c

Bio.o: Bio.c
	$(CC) $(CFLAGS)	Bio.c

clean:
	rm -f $(OBJS) $(OUT)

count:
	wc $(SOURCE) $(HEADER)