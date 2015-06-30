CC=g++
CFLAGS=-Wall -g -std=c++0x
LDFLAGS=-lopencv_core -lopencv_flann -L/usr/local/lib
EXEC=main

all: $(EXEC)

main: main.o parser.o arg.o
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

parser.o: parser.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

arg.o: arg.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

main.o: main.cpp parser.h arg.h
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -rf $(EXEC) *.o *.a *~
