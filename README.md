This is a C++ library for computing shortest paths with higher-order properties like curvature and torsion taken into account. It implements the algorithms of our ICCV 2013 paper [1].

[![Build Status](https://travis-ci.org/PetterS/vessel.png)](https://travis-ci.org/PetterS/vessel)

Compilation
-----------
The following is required to compile the library:
* C++ compiler.
* CMake
* Spii installed; see https://github.com/PetterS/spii 

All tests pass with the following compilers:
* Visual Studio 2012
* GCC 4.8 (Cygwin)
* GCC 4.6 (Ubuntu)
Earlier compilers might not work.

You can check travis.yml for the commands used to build the library and run all tests on Ubuntu.

References
----------
1. ...
