CC = gcc
CFLAGS = -O3 -Wall
LIBRARY_FLAGS = -lm

EXECUTABLES = meson_mass_co covarprop
MASS_FILES = meson_mass.o cofit.o minimize.o linalg.o
PROP_FILES = covarprop.o

main: meson_mass.o cofit.o minimize.o linalg.o
	$(CC) $(CFLAGS) -o meson_mass_co $(MASS_FILES) $(LIBRARY_FLAGS)

covarprop: covarprop.o
	$(CC) $(CFLAGS) -o covarprop $(PROP_FILES)

clean:
	rm -f $(EXECUTABLES) *.o

meson_mass.o: meson_mass.c
	$(CC) $(CFLAGS) -c meson_mass.c

cofit.o: cofit.c
	$(CC) $(CFLAGS) -c cofit.c

minimize.o: minimize.c
	$(CC) $(CFLAGS) -c minimize.c

linalg.o: linalg.c
	$(CC) $(CFLAGS) -c linalg.c

covarprop.o: covarprop.c
	$(CC) $(CFLAGS) -c covarprop.c
