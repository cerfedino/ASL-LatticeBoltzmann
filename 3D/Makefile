LONG_RUN_PARAMS=3 5 1 1000


BASE_VERSION_OBJ = previous_versions/v0/src/cmdline/main.o
CURRENT_VERSION_OBJ = src/cmdline/main.o
PAST_VERSIONS_OBJ = $(foreach file,$(wildcard previous_versions/*/src/cmdline/main.cpp),$(file:.cpp=.o))

.PHONY: all clean build 

all: build
clean:
	@find . \( -iname '*.o' -o -iname 'timestamp.txt' \) -print -delete
	@find previous_versions/ -mindepth 1 -maxdepth 1 -type d -print -exec bash -c 'rm -r {}/src/cmdline/{CMakeFiles,cmake_install.cmake,CMakeCache.txt,Makefile}' \;


build: $(CURRENT_VERSION_OBJ) $(BASE_VERSION_OBJ) $(PAST_VERSIONS_OBJ) 


# Recipe for building every version
%.o: %.cpp
	cd $(dir $<); pwd; cmake ../; make
