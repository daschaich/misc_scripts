CC = gcc
CFLAGS = -O3 -Wall

EXECUTABLE = covarprop
OBJECT_FILES = covarprop.o

main: covarprop.o
	$(CC) $(CFLAGS) -o $(EXECUTABLE) $(OBJECT_FILES) $(LIBRARY_FLAGS)

clean:
	rm -f $(EXECUTABLE) *.o

covarprop.o: covarprop.c
	$(CC) $(CFLAGS) -c covarprop.c
