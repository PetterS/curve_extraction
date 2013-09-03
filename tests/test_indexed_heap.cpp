
// This will compile the internal gheap tests.
//#define GHEAP_CPP11 1
//#include <tests.cpp>

#include <algorithm>
#include <iostream>

#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()
#include <catch.hpp>

#include <vessel/indexed_heap.h>
using namespace vessel;

void check_consistency()
{
	REQUIRE(heap::index.size() == heap::priority.size());

	for (int e = 0; e < heap::index.size(); ++e) {
		if (heap::index[e] >= 0) {
			// Element is in heap.
			CHECK(heap::heap[heap::index[e]] == e);
		}
		else {
			// Pass
		}
	}
}

TEST_CASE("indexed_heap/insert_pop", "")
{
	const int n = 10;
	heap::initialize(n);

	heap::insert(0, 15.1f); check_consistency();
	heap::insert(1, 30.1f); check_consistency();
	heap::insert(2, 45.1f); check_consistency();
	heap::insert(3, 1.1f); check_consistency();
	heap::insert(4, 2.1f); check_consistency();
	heap::insert(5, 0.1f); check_consistency();
	heap::insert(6, 5.1f); check_consistency();
	heap::insert(7, 6.1f); check_consistency();
	heap::insert(8, 7.1f); check_consistency();
	heap::insert(9, 8.1f); check_consistency();

	REQUIRE(heap::smallest_element() == 5);
	REQUIRE(heap::smallest_priority() == 0.1f);
	heap::pop(); check_consistency();

	REQUIRE(heap::smallest_element() == 3);
	REQUIRE(heap::smallest_priority() == 1.1f);
	heap::pop(); check_consistency();

	REQUIRE(heap::smallest_element() == 4);
	REQUIRE(heap::smallest_priority() == 2.1f);
	heap::pop(); check_consistency();

	REQUIRE(heap::smallest_element() == 6);
	REQUIRE(heap::smallest_priority() == 5.1f);
	heap::pop(); check_consistency();

	REQUIRE(heap::smallest_element() == 7);
	REQUIRE(heap::smallest_priority() == 6.1f);
	heap::pop(); check_consistency();

	REQUIRE(heap::smallest_element() == 8);
	REQUIRE(heap::smallest_priority() == 7.1f);
	heap::pop(); check_consistency();

	REQUIRE(heap::smallest_element() == 9);
	REQUIRE(heap::smallest_priority() == 8.1f);
	heap::pop(); check_consistency();

	REQUIRE(heap::smallest_element() == 0);
	REQUIRE(heap::smallest_priority() == 15.1f);
	heap::pop(); check_consistency();

	REQUIRE(heap::smallest_priority() == 30.1f);
	REQUIRE(heap::smallest_element() == 1);
	heap::pop(); check_consistency();

	REQUIRE(heap::smallest_element() == 2);
	REQUIRE(heap::smallest_priority() == 45.1f);
	heap::pop(); check_consistency();
}

TEST_CASE("indexed_heap/decrease_priority", "")
{
	const int n = 10;
	heap::initialize(n);

	heap::insert(0, 15.1f); check_consistency();
	heap::insert(1, 30.1f); check_consistency();
	heap::insert(2, 45.1f); check_consistency();
	heap::insert(3, 1.1f); check_consistency();
	heap::insert(4, 2.1f); check_consistency();
	heap::insert(5, 0.1f); check_consistency();
	heap::insert(6, 5.1f); check_consistency();
	heap::insert(7, 6.1f); check_consistency();
	heap::insert(8, 7.1f); check_consistency();
	heap::insert(9, 8.1f); check_consistency();

	heap::decrease_priority(1, 0.0f); check_consistency();

	REQUIRE(heap::smallest_element() == 1);
	REQUIRE(heap::smallest_priority() == 0.0f);
	heap::pop(); check_consistency();

	REQUIRE(heap::smallest_element() == 5);
	REQUIRE(heap::smallest_priority() == 0.1f);
	heap::pop(); check_consistency();

	REQUIRE(heap::smallest_element() == 3);
	REQUIRE(heap::smallest_priority() == 1.1f);
	heap::pop(); check_consistency();

	REQUIRE(heap::smallest_element() == 4);
	REQUIRE(heap::smallest_priority() == 2.1f);
	heap::pop(); check_consistency();

	REQUIRE(heap::smallest_element() == 6);
	REQUIRE(heap::smallest_priority() == 5.1f);
	heap::pop(); check_consistency();
}
