# Makefile for fric2d

CFLAGS =  -O -std=gnu99  -I. -DANSI
LIBRARIES = -lm

fric2d_3.2.8: fric2d_3.2.8.o getwords.o getoption.o ludcmp.o lubksb.o  nrutil.o
	gcc -o fric2d_3.2.8  fric2d_3.2.8.o getwords.o getoption.o ludcmp.o lubksb.o nrutil.o $(LIBRARIES)

clean:
	rm fric2d_3.2.8.o  getwords.o getoption.o ludcmp.o lubksb.o  nrutil.o

