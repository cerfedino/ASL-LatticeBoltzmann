CFLAGS = -g -march=native
OPTIM_FLAGS = -Ofast
LONG_RUN_PARAMS  = 400 200 300
SHORT_RUN_PARAMS = 400 200 300

PAST_VERSIONS_OBJ = $(foreach file,$(wildcard previous_versions/*/main.cpp),$(file:.cpp=_optimized_profile.o))
LONG_RUN_PARAMS_ESCAPED = $(shell echo $(LONG_RUN_PARAMS) | tr ' ' '_')
SHORT_RUN_PARAMS_ESCAPED = $(shell echo $(SHORT_RUN_PARAMS) | tr ' ' '_')

.PHONY: all clean build run_debug test leaks benchmark

all: build test leaks benchmark
clean:
	find . \( -iname '*.o' -o -iname 'timestamp_*.txt' \) -print -delete; rm -r artifacts/* output/*

build: ./src/main_optimized.o ./src/main_optimized_debug.o ./src/main_optimized_profile.o
%/main_optimized.o: %/main.cpp src/include/utils.c src/include/utils.h src/include/profile.c
	g++ $(CFLAGS) $(OPTIM_FLAGS) $< -o $@ src/include/utils.c src/include/utils.h src/include/profile.c
%/main_optimized_debug.o: %/main.cpp src/include/utils.c src/include/utils.h src/include/profile.c
	g++ $(CFLAGS) $(OPTIM_FLAGS) -DDEBUG $< -o $@ src/include/utils.c src/include/utils.h src/include/profile.c
%/main_optimized_profile.o: %/main.cpp src/include/utils.c src/include/utils.h src/include/profile.c
	g++ $(CFLAGS) $(OPTIM_FLAGS) -DPROFILE $< -o $@ src/include/utils.c src/include/utils.h src/include/profile.c


run_debug: output/00_latest_$(SHORT_RUN_PARAMS_ESCAPED)/timestamp_$(SHORT_RUN_PARAMS_ESCAPED).txt
output/00_latest_$(SHORT_RUN_PARAMS_ESCAPED)/timestamp_$(SHORT_RUN_PARAMS_ESCAPED).txt: ./src/main_optimized_debug.o
	mkdir -p output/; ./src/main_optimized_debug.o $(SHORT_RUN_PARAMS);

test: run_debug output/reference/timestamp_$(SHORT_RUN_PARAMS_ESCAPED).txt
	python3 compare.py --baseline=./output/reference/ --output=./output/00_latest_$(SHORT_RUN_PARAMS_ESCAPED)/
output/reference/timestamp_$(SHORT_RUN_PARAMS_ESCAPED).txt: ./reference-python/latticeboltzmann_prints.py
	python3 reference-python/latticeboltzmann_prints.py $(SHORT_RUN_PARAMS)

leaks: ./src/main_optimized_debug.o
	mkdir -p output; valgrind --leak-check=full ./src/main_optimized_debug.o $(SHORT_RUN_PARAMS)



benchmark: ./src/main_optimized.o $(PAST_VERSIONS_OBJ)
	mkdir -p artifacts/; python3 benchmark.py
profile: ./src/main_optimized_profile.o
	./src/main_optimized_profile.o $(LONG_RUN_PARAMS)