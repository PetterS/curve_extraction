# Author: petter.strandmark@gmail.com (Petter Strandmark)

SET(VESSEL_INTERNAL_SRC
	curvature.cpp
	data_term.cpp
	grid_mesh.cpp
	mesh.cpp
	shortest_path.cpp
	trilinear.cpp
)

# Also depend on the header files so that they appear in IDEs.
FILE(GLOB VESSEL_HEADERS *.h)

SET(VESSEL_LIBRARY_SOURCE
    ${VESSEL_INTERNAL_SRC}
    ${VESSEL_HEADERS})

ADD_LIBRARY(vessel STATIC ${VESSEL_LIBRARY_SOURCE})
TARGET_LINK_LIBRARIES(vessel ${VESSEL_LIBRARY_DEPENDENCIES})

INSTALL(TARGETS vessel
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)


