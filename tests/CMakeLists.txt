# Author: petter.strandmark@gmail.com (Petter Strandmark)

# Also depend on the header files so that they appear in IDEs.
FILE(GLOB VESSEL_HEADERS *.h)



IF (${BUILD_TESTING})
	MACRO (VESSEL_TEST NAME)
	  ADD_EXECUTABLE(test_${NAME} test_${NAME}.cpp)
	  TARGET_LINK_LIBRARIES(test_${NAME} vessel)

	  ADD_TEST(NAME test_${NAME}
	           COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_${NAME})
	ENDMACRO (VESSEL_TEST)

	VESSEL_TEST(curvature)
	VESSEL_TEST(data_term)
	VESSEL_TEST(grid_mesh)
	#VESSEL_TEST(indexed_heap)
	VESSEL_TEST(mesh)
	VESSEL_TEST(shortest_path)

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
