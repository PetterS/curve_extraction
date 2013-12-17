This is a C++ library for computing shortest paths with higher-order properties like curvature and torsion taken into account. It implements the algorithms of our ICCV 2013 paper [1].

[![Build Status](https://travis-ci.org/PetterS/curve_extraction.png)](https://travis-ci.org/PetterS/curve_extraction)

Compilation
-----------
The following is required to compile the library:
* C++ compiler.
* CMake
* Spii installed; see https://github.com/PetterS/spii 

All tests pass with the following compilers:
* Visual Studio 2013
* GCC 4.8 (Cygwin)
* GCC 4.7 (Ubuntu)

Earlier compilers might not work.

You can check .travis.yml for the commands used to build the library and run all tests on Ubuntu.

Video
-----
 * http://www.youtube.com/watch?v=9qjQj3I5pBc

References
----------
1. Petter Strandmark, Johannes Ulén, Fredrik Kahl, Leo Grady. [Shortest Paths with Curvature and Torsion](http://www2.maths.lth.se/vision/publications/publications/view_paper.php?paper_id=582). International Conference on Computer Vision. 2013.
