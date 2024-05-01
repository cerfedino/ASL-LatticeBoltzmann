# Compiles by default with the DEBUG flag
CFLAGS = -g -O0 -DDEBUG
.PHONY: all clean build run test leaks

all: build run test leaks
clean:
	find . -iname '*.o' -print -delete; rm -f output/timestamp.txt; rm -f output/reference/timestamp.txt

build: ./main_base.o ./main_optimized.o 
./main_optimized.o: main.cpp ./include/utils.o
	g++ -O3 $(CFLAGS) main.cpp ./include/utils.o -o ./main_optimized.o
./main_base.o: main_base.cpp ./include/utils.o
	g++ $(CFLAGS) main_base.cpp ./include/utils.o -o ./main_base.o
./include/utils.o: include/utils.h include/utils.c
	g++ -O3 $(CFLAGS) -c include/utils.c -o ./include/utils.o


run: output/timestamp.txt
output/timestamp.txt: ./main_optimized.o
	mkdir -p output/; ./main_optimized.o; echo "ciao mamma" > output/timestamp.txt

test: run output/reference/timestamp.txt
	python3 compare.py
output/reference/timestamp.txt: reference-python/latticeboltzmann_prints.py
	python3 reference-python/latticeboltzmann_prints.py; echo "ciao mamma" > output/reference/timestamp.txt; 

leaks: ./main_optimized.o
	valgrind --leak-check=full ./main_optimized.o