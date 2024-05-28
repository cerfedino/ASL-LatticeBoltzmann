$(info 'USING DEFAULT RUN_PARAMS')
RUN_PARAMS = 30 50 10 100
RUN_PARAMS_ESCAPED =  $(shell echo $(RUN_PARAMS) | tr ' ' '_')


BASE_VERSION_OBJ = previous_versions/v0/src/cmdline/main.o
CURRENT_VERSION_OBJ = src/cmdline/main.o
PAST_VERSIONS_OBJ = $(foreach file,$(wildcard previous_versions/*/src/cmdline/main.cpp),$(file:.cpp=.o))

.PHONY: all clean build 

all: build test
clean:
	@find . \( -iname '*.o' -o -iname 'timestamp.txt' \) -print -delete
	@find previous_versions/ -mindepth 1 -maxdepth 1 -type d -print -exec bash -c 'rm -r {}/src/cmdline/{CMakeFiles,cmake_install.cmake,CMakeCache.txt,Makefile}' \;
	@find . -mindepth 1 -maxdepth 1 -type d -print -exec bash -c 'rm -r src/cmdline/{CMakeFiles,cmake_install.cmake,CMakeCache.txt,Makefile}' \;

build: $(CURRENT_VERSION_OBJ) $(BASE_VERSION_OBJ) $(PAST_VERSIONS_OBJ) 


# Recipe for building any version
%/cmdline/main.o: %/cmdline/main.cpp
	cd $(dir $<); pwd; cmake ../; make
# Recipe for producing results for any version
%/cmdline/output/timestamp_$(RUN_PARAMS_ESCAPED).txt: %/cmdline/main_$(RUN_PARAMS_ESCAPED).o
	cd $(dir $<); ./main_$(RUN_PARAMS_ESCAPED).o $(RUN_PARAMS)


test: previous_versions/v0/src/cmdline/output/timestamp_$(RUN_PARAMS_ESCAPED).txt src/cmdline/output/timestamp_$(RUN_PARAMS_ESCAPED).txt
	python3 compare.py --baseline=previous_versions/v0/src/cmdline/output/ --output=src/cmdline/output/