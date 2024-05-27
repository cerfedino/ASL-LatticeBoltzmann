CFLAGS = -g -march=native
OPTIM_FLAGS = -Ofast



ifeq ($(firstword $(MAKECMDGOALS)),benchmark)
    RUN_PARAMS = 400 200 5000
else ifeq ($(firstword $(MAKECMDGOALS)),profile)
    RUN_PARAMS = 400 200 5000
else
	# $(info 'USING DEFAULT RUN_PARAMS')
	RUN_PARAMS = 400 200 300
endif
RUN_PARAMS_ESCAPED =  $(shell echo $(RUN_PARAMS) | tr ' ' '_')
RUN_PARAMS_DEFINITIONS = -DMNx=$(word 1,$(RUN_PARAMS)) -DMNy=$(word 2,$(RUN_PARAMS)) -DMNt=$(word 3,$(RUN_PARAMS))




PAST_VERSIONS_OBJ = $(foreach file,$(wildcard previous_versions/*/main.cpp),$(file:.cpp=_optimized_profile_$(RUN_PARAMS_ESCAPED).o))

.PHONY: all clean test leaks benchmark profile


all: test leaks benchmark profile
clean:
	find . \( -iname '*.o' -o -iname 'timestamp_*.txt' \) -print -delete; rm -r artifacts/* output/*

./previous_versions/v0/main_optimized_profile_$(RUN_PARAMS_ESCAPED).o: ./previous_versions/v0/main.cpp src/include/utils.c src/include/profile.c
	g++ $(CFLAGS) -DPROFILE $< -o $@ src/include/profile.c
%/main_optimized_$(RUN_PARAMS_ESCAPED).o: %/main.cpp src/include/utils.c src/include/profile.c
	g++ $(CFLAGS) $(OPTIM_FLAGS) $(RUN_PARAMS_DEFINITIONS) $< -o $@ src/include/profile.c
%/main_optimized_debug_$(RUN_PARAMS_ESCAPED).o: %/main.cpp src/include/utils.c src/include/profile.c
	g++ $(CFLAGS) $(OPTIM_FLAGS) $(RUN_PARAMS_DEFINITIONS) -DDEBUG $< -o $@ src/include/profile.c
%/main_optimized_profile_$(RUN_PARAMS_ESCAPED).o: %/main.cpp src/include/utils.c src/include/profile.c
	g++ $(CFLAGS) $(OPTIM_FLAGS) $(RUN_PARAMS_DEFINITIONS) -DPROFILE $< -o $@ src/include/profile.c


output/00_latest_$(RUN_PARAMS_ESCAPED)/timestamp_$(RUN_PARAMS_ESCAPED).txt: ./src/main_optimized_debug_$(RUN_PARAMS_ESCAPED).o
	mkdir -p output/; ./src/main_optimized_debug_$(RUN_PARAMS_ESCAPED).o $(RUN_PARAMS);

test: output/00_latest_$(RUN_PARAMS_ESCAPED)/timestamp_$(RUN_PARAMS_ESCAPED).txt   output/reference/timestamp_$(RUN_PARAMS_ESCAPED).txt
	python3 compare.py --baseline=./output/reference/ --output=./output/00_latest_$(RUN_PARAMS_ESCAPED)/
output/reference/timestamp_$(RUN_PARAMS_ESCAPED).txt: ./reference-python/latticeboltzmann_prints.py
	python3 reference-python/latticeboltzmann_prints.py $(RUN_PARAMS)

leaks: ./src/main_optimized_debug_$(RUN_PARAMS_ESCAPED).o
	mkdir -p output; valgrind --leak-check=full ./src/main_optimized_debug_$(RUN_PARAMS_ESCAPED).o $(RUN_PARAMS)



benchmark: ./src/main_optimized_$(RUN_PARAMS_ESCAPED).o $(PAST_VERSIONS_OBJ)
	mkdir -p artifacts/; python3 benchmark.py $(RUN_PARAMS);
profile: ./src/main_optimized_profile._$(RUN_PARAMS_ESCAPED).o
	./src/main_optimized_profile_$(RUN_PARAMS_ESCAPED).o $(RUN_PARAMS)