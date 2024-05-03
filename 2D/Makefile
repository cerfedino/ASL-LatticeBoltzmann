CFLAGS = -g
OPTIM_FLAGS = -Ofast
.PHONY: all clean build run_debug test leaks benchmark

all: build test leaks benchmark
clean:
	find . -iname '*.o' -print -delete; rm -f output/timestamp_debug.txt; rm -f output/reference/timestamp.txt

build: ./main_base.o ./main_base_debug.o ./main_optimized.o ./main_optimized_debug.o 
./main_optimized.o: main.cpp ./include/utils.c ./include/utils.h ./include/profile.c
	g++ $(CFLAGS) $(OPTIM_FLAGS) main.cpp ./include/utils.c ./include/profile.c -o ./main_optimized.o
./main_optimized_debug.o: main.cpp ./include/utils.c ./include/utils.h ./include/profile.c
	g++ $(CFLAGS) $(OPTIM_FLAGS) -DDEBUG main.cpp ./include/utils.c  ./include/profile.c -o ./main_optimized_debug.o
./main_optimized_profile.o: main.cpp ./include/utils.c ./include/utils.h ./include/profile.c
	g++ $(CFLAGS) $(OPTIM_FLAGS) -DPROFILE main.cpp ./include/utils.c  ./include/profile.c -o ./main_optimized_profile.o
./main_base.o: main_base.cpp ./include/utils.c ./include/utils.h ./include/profile.c
	g++ $(CFLAGS) main_base.cpp ./include/utils.c  ./include/profile.c -o ./main_base.o
./main_base_debug.o: main_base.cpp ./include/utils.c ./include/utils.h ./include/profile.c
	g++ $(CFLAGS) -DDEBUG main_base.cpp ./include/utils.c ./include/profile.c -o ./main_base_debug.o

run_debug: output/timestamp_debug.txt
output/timestamp_debug.txt: ./main_optimized_debug.o
	mkdir -p output/; ./main_optimized_debug.o; echo "ciao mamma" > output/timestamp_debug.txt

test: run_debug output/reference/timestamp.txt
	python3 compare.py
output/reference/timestamp.txt: reference-python/latticeboltzmann_prints.py
	python3 reference-python/latticeboltzmann_prints.py; echo "ciao mamma" > output/reference/timestamp.txt; 

leaks: ./main_optimized_debug.o
	valgrind --leak-check=full ./main_optimized_debug.o
benchmark: ./main_optimized.o ./main_base.o
	./main_optimized.o; ./main_base.o
profile: ./main_optimized_profile.o
	./main_optimized_profile.o