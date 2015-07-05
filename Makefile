CC=g++
CFLAGS=-Wall -g -std=c++0x
LDFLAGS=-lopencv_core -lopencv_flann -L/usr/local/lib
EXEC=main

all: $(EXEC)

main: main.o parser.o arg.o method.o quilting.o
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

parser.o: parser.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

arg.o: arg.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

method.o: method.cpp arg.h
	$(CC) -o $@ -c $< $(CFLAGS)

quilting.o: quilting.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

main.o: main.cpp parser.h arg.h method.h quilting.h
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -rf *.o *.a *~ 

mrproper: clean
	rm -rf $(EXEC) 
