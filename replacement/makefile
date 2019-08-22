CC = clang
CFLAGS = -Wall
OBJECTS = ../libs/matrix.o
LDLIBS = -lreadline

all: main
	make clean

main: $(OBJECTS)

clean:
	rm -rf *.o

purge: clean
	rm -f main

remake: purge
	make
	make clean
