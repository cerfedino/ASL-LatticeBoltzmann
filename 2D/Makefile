# Compiles by default with the DEBUG flag
CFLAGS = -DDEBUG
.PHONY: all clean build run test

all: build run test
clean:
	find . -iname '*.o' -print -delete; rm -f output/timestamp.txt; rm -f output/reference/timestamp.txt


build: ./main.o
./main.o: main.cpp ./include/utils.o
	g++ $(CFLAGS) main.cpp ./include/utils.o -o ./main.o
./include/utils.o: include/utils.h include/utils.c
	g++ $(CFLAGS) -c include/utils.c -o ./include/utils.o

run: output/timestamp.txt
output/timestamp.txt: ./main.o
	mkdir -p output/; echo "ciao mamma" > output/timestamp.txt ; ./main.o

test: run output/reference/timestamp.txt
	python3 compare.py
output/reference/timestamp.txt: reference-python/latticeboltzmann_prints.py
	echo "ciao mamma" > output/reference/timestamp.txt; python3 reference-python/latticeboltzmann_prints.py