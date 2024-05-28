## Developer guide
Notice that every version is inside the `previous_versions` folder. To specify compiler flags for a specific version, open the `CMakeLists.txt` file  (e.g `previous_versions/v1/CMakeLists.txt`) and edit the `add_compile_options()` option.

### Using the Makefile

- `make build`\
Builds all current and previous versions

- `make test`\
Runs the current version and the baseline and invokes the `compare.py` to compare the reference and current version output folders.