.SECONDEXPANSION:
.PHONY: all clean build benchmark profile


RUN_PARAMS = 20 20 20 100
ifeq ($(firstword $(MAKECMDGOALS)),test)
	RUN_PARAMS = 20 20 20 400
else ifeq ($(firstword $(MAKECMDGOALS)),sizevscyclesplot)
    RUN_PARAMS = 16 8 50
else ifeq ($(firstword $(MAKECMDGOALS)),sizevscyclesplot_report)
	RUN_PARAMS = 4 4 4 10
endif

RUN_PARAMS_ESCAPED =  $(shell echo $(RUN_PARAMS) | tr ' ' '_')
RUN_PARAMS_DEFINITIONS = -DMNx=$(word 1,$(RUN_PARAMS)) -DMNy=$(word 2,$(RUN_PARAMS)) -DMNz=$(word 3,$(RUN_PARAMS)) -DMNt=$(word 4,$(RUN_PARAMS))

BASE_VERSION_OBJ = previous_versions/v0/src/cmdline/main_$(RUN_PARAMS_ESCAPED).o
CURRENT_VERSION_OBJ = src/cmdline/main_$(RUN_PARAMS_ESCAPED).o

PAST_VERSIONS_CPP = $(wildcard previous_versions/*/src/cmdline/main.cpp)
PAST_VERSIONS_OBJ = $(foreach file,$(PAST_VERSIONS_CPP),$(file:.cpp=_$(RUN_PARAMS_ESCAPED).o))


all: build test profile

clean:
	@find . \( -iname '*.o' -o -iname 'timestamp.txt' \) -print -delete
	@find previous_versions -type d -print -exec bash -c 'rm -r {}/src/cmdline/{CMakeFiles,cmake_install.cmake,CMakeCache.txt,Makefile,output/*} 2> /dev/null' \;
	@find . -type d -exec bash -c 'rm -r src/cmdline/{CMakeFiles,cmake_install.cmake,CMakeCache.txt,Makefile,output/*} 2> /dev/null' \;

build: $(CURRENT_VERSION_OBJ) $(BASE_VERSION_OBJ) $(PAST_VERSIONS_OBJ) 


RANGE_I = $(shell seq 1 1 25)
multiply = $(shell echo $$(( $(1) * $(2) )))
NX_VALUES := $(foreach i,$(RANGE_I),$(eval NX_$(i)=$(call multiply, $(word 1,$(RUN_PARAMS)), $(i))))
NY_VALUES := $(foreach i,$(RANGE_I),$(eval NY_$(i)=$(call multiply, $(word 2,$(RUN_PARAMS)), $(i))))
NZ_VALUES := $(foreach i,$(RANGE_I),$(eval NZ_$(i)=$(call multiply, $(word 3,$(RUN_PARAMS)), $(i))))
PAST_VERSIONS_PROFILE_OBJ_RANGE = $(foreach i,$(RANGE_I),$(foreach file,$(PAST_VERSIONS_CPP),$(file:.cpp=_$(NX_$(i))_$(NY_$(i))_$(NZ_$(i))_$(word 4,$(RUN_PARAMS))_profiler.o)))


define PREPROCESS
$(eval MNx=$(shell echo "$(1)" | sed 's/_profiler//' | awk -F'[_.]' '{print $$(NF-4)}'))
$(eval MNy=$(shell echo "$(1)" | sed 's/_profiler//' | awk -F'[_.]' '{print $$(NF-3)}'))
$(eval MNz=$(shell echo "$(1)" | sed 's/_profiler//' | awk -F'[_.]' '{print $$(NF-2)}'))
$(eval MNt=$(shell echo "$(1)" | sed 's/_profiler//' | awk -F'[_.]' '{print $$(NF-1)}'))
echo "MNx: $(MNx) MNy: $(MNy) MNt: $(MNt)"
endef


%.o: $$(dir $$@)main.cpp
	$(call PREPROCESS,$@)
	cd $(dir $<); pwd; cmake -DMNx=$(MNx) -DMNy=$(MNy) -DMNz=$(MNz) -DMNt=$(MNt) ../; make


sizevscyclesplot: $(PAST_VERSIONS_PROFILE_OBJ_RANGE)
	mkdir -p artifacts/; python3 sizevscyclesplot.py $(RUN_PARAMS) $(RANGE_I)

sizevscyclesplot_report: $(PAST_VERSIONS_PROFILE_OBJ_RANGE)
	mkdir -p artifacts/; python3 sizevscyclesplot_report.py $(RUN_PARAMS) $(RANGE_I)


# Recipe for building any version
%/cmdline/main_$(RUN_PARAMS_ESCAPED).o: %/cmdline/main.cpp
	cd $(dir $<); pwd; cmake $(RUN_PARAMS_DEFINITIONS) ../ -B.; make

# %/cmdline/main_$(RUN_PARAMS_ESCAPED).o: %/cmdline/main.cpp
# 	cd $(dir $<); pwd; cmake $(RUN_PARAMS_DEFINITIONS) ../; make 

# Recipe for producing results for any version
%/cmdline/output/timestamp_$(RUN_PARAMS_ESCAPED).txt: %/cmdline/main_$(RUN_PARAMS_ESCAPED).o
	cd $(dir $<); ./main_$(RUN_PARAMS_ESCAPED).o $(RUN_PARAMS)


test: previous_versions/v0/src/cmdline/output/timestamp_$(RUN_PARAMS_ESCAPED).txt src/cmdline/output/timestamp_$(RUN_PARAMS_ESCAPED).txt
	python3 compare.py --baseline=previous_versions/v0/src/cmdline/output/ --output=src/cmdline/output/

benchmark: $(PAST_VERSIONS_OBJ)
	mkdir -p artifacts/; python3 benchmark.py $(RUN_PARAMS);
profile: $(CURRENT_VERSION_OBJ)
	cd $(dir $<); ./main_$(RUN_PARAMS_ESCAPED)_profiler.o $(RUN_PARAMS)
benchmark_cur: $(CURRENT_VERSION_OBJ)
	./src/cmdline/main_$(RUN_PARAMS_ESCAPED)_benchmark.o