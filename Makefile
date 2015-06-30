CC=g++
CFLAGS=-Wall -g
LDFLAGS=-lopencv_core -lopencv_flann -L/usr/local/lib
EXEC=main

all: $(EXEC)

main: main.o parser.o
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

parser.o: parser.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

main.o: main.cpp parser.h
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -rf $(EXEC) *~
