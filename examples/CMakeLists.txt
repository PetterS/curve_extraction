# Author: petter.strandmark@gmail.com (Petter Strandmark)

# Also depend on the header files so that they appear in IDEs.
FILE(GLOB VESSEL_HEADERS *.h)

MACRO (VESSEL_EXAMPLE NAME)
  ADD_EXECUTABLE(${NAME} ${NAME}.cpp)
  TARGET_LINK_LIBRARIES(${NAME} vessel)
ENDMACRO (VESSEL_EXAMPLE)

VESSEL_EXAMPLE(curve_completion)
VESSEL_EXAMPLE(demo_2d)
VESSEL_EXAMPLE(demo_3d)
VESSEL_EXAMPLE(tree_reconstruction)
