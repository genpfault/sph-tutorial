sph-tutorial
============

Experimenting with performance evaluation and improvements to Brandon Pelfrey's SPH fluid simulation tutorial.

Requirements:
* CMake >= 3.7.2
* C++11 compiler
* OpenGL support
* OpenMP support

Controls:
* Q/Escape: Exit
* Space: Add more particles
* Mouse button: Attract nearby particles

You can use the [`OMP_NUM_THREADS` environment variable](https://gcc.gnu.org/onlinedocs/libgomp/OMP_005fNUM_005fTHREADS.html#OMP_005fNUM_005fTHREADS) to limit the number of threads used by OpenMP.

Build quickstart:

    git clone https://github.com/genpfault/sph-tutorial.git
    cd sph-tutorial
    git submodule update --init --recursive
    mkdir build
    cd build
    cmake ../ -DCMAKE_BUILD_TYPE=RelWithDebInfo
    cmake --build .
