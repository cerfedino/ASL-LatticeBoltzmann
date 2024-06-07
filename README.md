# Project topic: Lattice Boltzmann method
*Authors: Albert Cerfeda, Karlo Piskor, Dominic Wüst, Kamelia Ivanova*




### Using the Makefile
- `make clean`\
Removes every compiled binary, artifact and output folder

- `make test`\
Builds and executes the current version and the baseline, after which it invokes the `compare.py` to compare the reference and current version output folders.

- `make leaks`\
Builds the current version and checks for memory leaks using Valgrind

- `make benchmark`\
Builds all versions both with profiling and the Performance Application Programming Interface (PAPI) enabled. Benchmarks exact memory movements, cycles and flops and plots the main performance comparison plot.

- `make sizevscyclesplot`\
Builds all versions with varying input sizes, benchmarks them and plots the 'Size/Cycles vs Cycles' plot

- `make profile`\
Builds and profiles separate loops for the current version.

### Previous versions
Notice that every version is inside the `previous_versions` folder. Inside each folder there is a `.env` that contains important info about the version such as a title, description, commit hash and the equation for the work and allocated bytes.

Every folder under `previous_versions` has to be named like follows: `{_}v<int>` (e.g `v1`, `_v2`)


To specify compiler flags for a specific 3D version, open the `CMakeLists.txt` file  (e.g `previous_versions/v1/CMakeLists.txt`) and edit the `add_compile_options()` option.
