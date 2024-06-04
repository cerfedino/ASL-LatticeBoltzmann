.PHONY: all clean build benchmark profile


RUN_PARAMS = 30 50 10 500
ifeq ($(firstword $(MAKECMDGOALS)),test)
    RUN_PARAMS = 30 50 10 100
endif

RUN_PARAMS_ESCAPED =  $(shell echo $(RUN_PARAMS) | tr ' ' '_')
RUN_PARAMS_DEFINITIONS = -DMNx=$(word 1,$(RUN_PARAMS)) -DMNy=$(word 2,$(RUN_PARAMS)) -DMNz=$(word 3,$(RUN_PARAMS)) -DMNt=$(word 4,$(RUN_PARAMS))
$(info 'Run params: $(RUN_PARAMS_DEFINITIONS)')

BASE_VERSION_OBJ = previous_versions/v0/src/cmdline/main_$(RUN_PARAMS_ESCAPED).o
CURRENT_VERSION_OBJ = src/cmdline/main_$(RUN_PARAMS_ESCAPED).o
PAST_VERSIONS_OBJ = $(foreach file,$(wildcard previous_versions/*/src/cmdline/main.cpp),$(file:.cpp=_$(RUN_PARAMS_ESCAPED).o))


all: build test benchmark profile

clean:
	@find . \( -iname '*.o' -o -iname 'timestamp.txt' \) -print -delete
	@find previous_versions -type d -print -exec bash -c 'rm -r {}/src/cmdline/{CMakeFiles,cmake_install.cmake,CMakeCache.txt,Makefile,output/*} 2> /dev/null' \;
	@find . -type d -exec bash -c 'rm -r src/cmdline/{CMakeFiles,cmake_install.cmake,CMakeCache.txt,Makefile,output/*} 2> /dev/null' \;

build: $(CURRENT_VERSION_OBJ) $(BASE_VERSION_OBJ) $(PAST_VERSIONS_OBJ) 

# Recipe for building any version
%/cmdline/main_$(RUN_PARAMS_ESCAPED).o: %/cmdline/main.cpp
	cd $(dir $<); pwd; cmake $(RUN_PARAMS_DEFINITIONS) ../ -B.; make
# Recipe for producing results for any version
%/cmdline/output/timestamp_$(RUN_PARAMS_ESCAPED).txt: %/cmdline/main_$(RUN_PARAMS_ESCAPED).o
	cd $(dir $<); ./main_$(RUN_PARAMS_ESCAPED).o $(RUN_PARAMS)


test: previous_versions/v0/src/cmdline/output/timestamp_$(RUN_PARAMS_ESCAPED).txt src/cmdline/output/timestamp_$(RUN_PARAMS_ESCAPED).txt
	python3 compare.py --baseline=previous_versions/v0/src/cmdline/output/ --output=src/cmdline/output/

benchmark: $(PAST_VERSIONS_OBJ)
	mkdir -p artifacts/; python3 benchmark.py $(RUN_PARAMS);
profile: $(CURRENT_VERSION_OBJ)
	cd $(dir $<); ./main_$(RUN_PARAMS_ESCAPED)_profiler.o $(RUN_PARAMS)