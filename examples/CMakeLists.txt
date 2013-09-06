# Author: petter.strandmark@gmail.com (Petter Strandmark)

# Also depend on the header files so that they appear in IDEs.
FILE(GLOB CURVE_EXTRACTION_HEADERS *.h)

MACRO (CURVE_EXTRACTION_EXAMPLE NAME)
  ADD_EXECUTABLE(${NAME} ${NAME}.cpp)
  TARGET_LINK_LIBRARIES(${NAME} curve_extraction)
ENDMACRO (CURVE_EXTRACTION_EXAMPLE)

CURVE_EXTRACTION_EXAMPLE(curve_completion)
CURVE_EXTRACTION_EXAMPLE(demo_2d)
CURVE_EXTRACTION_EXAMPLE(demo_3d)
CURVE_EXTRACTION_EXAMPLE(tree_reconstruction)
