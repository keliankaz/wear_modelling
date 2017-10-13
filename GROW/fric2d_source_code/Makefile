# Makefile for fric2d

CFLAGS =  -std=gnu99  -I. -DANSI
LIBRARIES = -lm

fric2d: fric2d.o getwords.o getoption.o ludcmp.o lubksb.o  nrutil.o
	cc -o fric2d fric2d.o getwords.o getoption.o ludcmp.o lubksb.o nrutil.o $(LIBRARIES)
