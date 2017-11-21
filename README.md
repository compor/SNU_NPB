# NAS benchmarks

This repo contains the source code for the [NAS][1] Parallel Benchmarks [implemented in C][2].
Currently, this includes the serial and OMP parallel versions.

The main modifications performed in preparation for usage with [this][3] build harness were:

- move each program source to a `src` subdirectory
- remove existing build files (for `make` and `cmake`)


[1]: www.nas.nasa.gov/publications/npb.html
[2]: http://aces.snu.ac.kr/software/snu-npb/
[3]: https://github.com/compor/nauseous

