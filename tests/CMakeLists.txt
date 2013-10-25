# Author: petter.strandmark@gmail.com (Petter Strandmark)

# Also depend on the header files so that they appear in IDEs.
FILE(GLOB CURVE_EXTRACTION_HEADERS *.h)



IF (${BUILD_TESTING})
	MACRO (CURVE_EXTRACTION_TEST NAME)
	  ADD_EXECUTABLE(${NAME} ${NAME}.cpp)
	  TARGET_LINK_LIBRARIES(${NAME} curve_extraction)

	  ADD_TEST(NAME ${NAME}
	           COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${NAME})
	ENDMACRO (CURVE_EXTRACTION_TEST)

	FILE(GLOB CURVE_EXTRACTION_TEST_FILES test_*.cpp)
	FOREACH (TEST_FILE ${CURVE_EXTRACTION_TEST_FILES})
		GET_FILENAME_COMPONENT(TEST_NAME ${TEST_FILE} NAME_WE)
		MESSAGE("-- Adding test: " ${TEST_NAME})
		CURVE_EXTRACTION_TEST(${TEST_NAME})
	ENDFOREACH()

	#
	# Valgrind (memcheck) testing.
	#
	find_program(VALGRIND NAMES valgrind PATH /usr/bin /usr/local/bin)
	if (VALGRIND)
		message("-- Valgrind found; configuring memcheck tests.")
		macro(VALGRIND_TEST NAME)
			add_test(NAME memcheck_${NAME}
			         COMMAND ${VALGRIND} --error-exitcode=1 --leak-check=full ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_${NAME})
		endmacro(VALGRIND_TEST)

		# Run shorter tests through memcheck.
		VALGRIND_TEST(curvature)
		VALGRIND_TEST(data_term)
		VALGRIND_TEST(shortest_path)
	endif (VALGRIND)

	#
	# Cppcheck testing.
	#
	find_program(CPPCHECK NAMES cppcheck PATH /usr/bin /usr/local/bin)
	if (CPPCHECK)
		message("-- Cppcheck found; configuring test.")
		add_test(NAME cppcheck
		         COMMAND ${CPPCHECK} ${CMAKE_SOURCE_DIR}/source/ -q --error-exitcode=1)
	endif (CPPCHECK)

ENDIF (${BUILD_TESTING})
