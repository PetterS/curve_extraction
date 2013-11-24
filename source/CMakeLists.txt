# Author: petter.strandmark@gmail.com (Petter Strandmark)

FILE(GLOB CURVE_EXTRACTION_INTERNAL_SRC *.cpp)

# Also depend on the header files so that they appear in IDEs.
SET(CURVE_EXTRACTION_LIBRARY_SOURCE
    ${CURVE_EXTRACTION_INTERNAL_SRC}
    ${CURVE_EXTRACTION_HEADERS})

ADD_LIBRARY(curve_extraction STATIC ${CURVE_EXTRACTION_LIBRARY_SOURCE})
TARGET_LINK_LIBRARIES(curve_extraction ${CURVE_EXTRACTION_LIBRARY_DEPENDENCIES})

INSTALL(TARGETS curve_extraction
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)


