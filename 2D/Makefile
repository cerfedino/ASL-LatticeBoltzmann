CFLAGS = -g -march=native
OPTIM_FLAGS = -Ofast
LONG_RUN_PARAMS  = 400 200 5000
SHORT_RUN_PARAMS = 400 200 300

PAST_VERSIONS_OBJ = $(foreach file,$(wildcard previous_versions/*/main.cpp),$(file:.cpp=.o))
LONG_RUN_PARAMS_ESCAPED = $(shell echo $(LONG_RUN_PARAMS) | tr ' ' '_')
SHORT_RUN_PARAMS_ESCAPED = $(shell echo $(SHORT_RUN_PARAMS) | tr ' ' '_')

.PHONY: all clean build run_debug test leaks benchmark

all: build test leaks benchmark
clean:
	find . \( -iname '*.o' -o -iname 'timestamp_*.txt' \) -print -delete; rm -r artifacts/*

build: ./main_base.o ./main_base_debug.o ./main_optimized.o ./main_optimized_debug.o ./main_optimized_profile.o
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

# Recipe for building old versions
%.o: %.cpp ./include/utils.c ./include/utils.h ./include/profile.c
	g++ $(CFLAGS) $(OPTIM_FLAGS) -DPROFILE $< -o $@ ./include/utils.c ./include/utils.h ./include/profile.c
# Sometimes we use "timestamp" files as a way for Make to keep track of the last time a target was run
run_debug: output/00_latest_$(SHORT_RUN_PARAMS_ESCAPED)/timestamp_$(SHORT_RUN_PARAMS_ESCAPED).txt
output/00_latest_$(SHORT_RUN_PARAMS_ESCAPED)/timestamp_$(SHORT_RUN_PARAMS_ESCAPED).txt: ./main_optimized_debug.o
	mkdir -p output/; ./main_optimized_debug.o $(SHORT_RUN_PARAMS);

test: run_debug output/reference/timestamp_$(SHORT_RUN_PARAMS_ESCAPED).txt
	python3 compare.py --baseline=./output/reference/ --output=output/00_latest_$(SHORT_RUN_PARAMS_ESCAPED)/
output/reference/timestamp_$(SHORT_RUN_PARAMS_ESCAPED).txt: reference-python/latticeboltzmann_prints.py
	python3 reference-python/latticeboltzmann_prints.py $(SHORT_RUN_PARAMS)

leaks: ./main_optimized_debug.o
	mkdir -p output; valgrind --leak-check=full ./main_optimized_debug.o $(SHORT_RUN_PARAMS)



benchmark: ./main_optimized.o $(PAST_VERSIONS_OBJ)
	mkdir -p artifacts/; python3 benchmark.py
profile: ./main_optimized_profile.o
	./main_optimized_profile.o $(LONG_RUN_PARAMS)