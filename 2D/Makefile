CFLAGS = -g -march=native
OPTIM_FLAGS = -Ofast



ifeq ($(firstword $(MAKECMDGOALS)),benchmark)
    RUN_PARAMS = 400 200 300
else ifeq ($(firstword $(MAKECMDGOALS)),profile)
    RUN_PARAMS = 400 200 2000
else
	# $(info 'USING DEFAULT RUN_PARAMS')
	RUN_PARAMS = 400 200 300
endif
RUN_PARAMS_ESCAPED =  $(shell echo $(RUN_PARAMS) | tr ' ' '_')

PAST_VERSIONS_CPP = $(wildcard previous_versions/*/src/main.cpp)
PAST_VERSIONS_PROFILE_OBJ = $(foreach file,$(PAST_VERSIONS_CPP),$(file:.cpp=_optimized_profile_$(RUN_PARAMS_ESCAPED).o))


.PHONY: all clean test leaks benchmark profile


all: test leaks benchmark profile
clean:
	find . \( -iname '*.o' -o -iname 'timestamp_*.txt' \) -print -delete; rm -r artifacts/* output/*


define PREPROCESS
$(eval MNx=$(shell echo "$(1)" | awk -F'[_.]' '{print $$(NF-3)}'))
$(eval MNy=$(shell echo "$(1)" | awk -F'[_.]' '{print $$(NF-2)}'))
$(eval MNt=$(shell echo "$(1)" | awk -F'[_.]' '{print $$(NF-1)}'))
$(eval DDEBUG=$(if $(findstring _debug,$(1)),-DDEBUG,))
$(eval DPROFILE=$(if $(findstring _profile,$(1)),-DPROFILE,))
echo "MNx: $(MNx) MNy: $(MNy) MNt: $(MNt) DDEBUG: $(DDEBUG) DPROFILE: $(DPROFILE)"
endef

./previous_versions/v1/src/main_%.o: ./previous_versions/v1/src/main.cpp src/include/utils.c src/include/profile.c
	$(call PREPROCESS,$@)
	g++ $(CFLAGS) 							 -DMNx=$(MNx) -DMNy=$(MNy) -DMNt=$(MNt) $(DDEBUG) $(DPROFILE) $< -o $@ src/include/profile.c

%.o: ./src/main.cpp src/include/utils.c src/include/profile.c
	$(call PREPROCESS,$@)
	g++ $(CFLAGS) $(OPTIM_FLAGS) -DMNx=$(MNx) -DMNy=$(MNy) -DMNt=$(MNt) $(DDEBUG) $(DPROFILE) $< -o $@ src/include/profile.c



output/00_latest_$(RUN_PARAMS_ESCAPED)/timestamp_$(RUN_PARAMS_ESCAPED).txt: ./src/main_optimized_debug_$(RUN_PARAMS_ESCAPED).o
	mkdir -p output/; ./src/main_optimized_debug_$(RUN_PARAMS_ESCAPED).o $(RUN_PARAMS);

test: output/00_latest_$(RUN_PARAMS_ESCAPED)/timestamp_$(RUN_PARAMS_ESCAPED).txt output/reference/timestamp_$(RUN_PARAMS_ESCAPED).txt
	python3 compare.py --baseline=./output/reference/ --output=./output/00_latest_$(RUN_PARAMS_ESCAPED)/
output/reference/timestamp_%.txt:
	$(call PREPROCESS,$@)
	python3 reference-python/latticeboltzmann_prints.py $(MNx) $(MNy) $(MNt)

leaks: ./src/main_optimized_debug_$(RUN_PARAMS_ESCAPED).o
	mkdir -p output; valgrind --leak-check=full ./src/main_optimized_debug_$(RUN_PARAMS_ESCAPED).o $(RUN_PARAMS)



benchmark: $(PAST_VERSIONS_PROFILE_OBJ)
	mkdir -p artifacts/; python3 benchmark.py $(RUN_PARAMS);
profile: ./src/main_optimized_profile_$(RUN_PARAMS_ESCAPED).o
	./src/main_optimized_profile_$(RUN_PARAMS_ESCAPED).o $(RUN_PARAMS)
