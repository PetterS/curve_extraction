// Petter Strandmark 2013.

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <queue>
#include <set>
#include <stdexcept>

#include <curve_extraction/shortest_path.h>

namespace curve_extraction
{

double shortest_path(int n, const std::set<int>& start_set, const std::set<int>& end_set,
                     const std::function<void(int, std::vector<Neighbor>* neighbors)>& neighbors,
                     std::vector<int>* path, const std::function<double(int)>* get_lower_bound,
                     const ShortestPathOptions& options)
{
	// Datatype used for the internal storage. Using float saves memory
	// for really large problems.
	typedef float queue_cost;
	const queue_cost infinity = std::numeric_limits<queue_cost>::max();

	// The current distance from the start set to node i.
	std::vector<queue_cost> distance(n, infinity);

	// The estimated distance from node i to the end. This
	// storage in needed in order to erase entries from the
	// queue.
	std::vector<queue_cost> estimation;
	if (get_lower_bound) {
		estimation.resize(n, 0);
	}

	// The previous node in the shortest path from the start
	// to node i.
	std::vector<int> previous(n, -1);

	// The priority queue specifying the order in which to
	// process the nodes. Store costs as floats to save
	// memory.
	std::set<std::pair<queue_cost, int> > prio_queue;

	std::vector<Neighbor> neighbor_storage;
	// Allocate storage for 100 neighbors. Each call to clear() will
	// not deallocate the storage.
	neighbor_storage.reserve(100);

	if (options.store_visited) {
		options.visit_time.resize(0);
		options.visit_time.resize(n, -1);
	}

	// Check and put the start set into the queue.
	if (start_set.size() == 0) {
		throw std::runtime_error("shortest_path: empty start set");
	}
	for (auto itr = start_set.begin(); itr != start_set.end(); ++itr) {
		if (*itr < 0 || *itr >= n) {
			throw std::runtime_error("shortest_path: Invalid start set.");
		}
		distance[*itr] = 0;
		prio_queue.insert(std::make_pair(0, *itr));
	}
	// Check end_set.
	for (auto itr = end_set.begin(); itr != end_set.end(); ++itr) {
		if (*itr < 0 || *itr >= n) {
			throw std::runtime_error("shortest_path: Invalid end set.");
		}
	}

	int n_visited = 0;
	bool first_print = true;
	auto last_time = std::clock();

	// We have already stored the shortest path in the path vector.
	int end_node = -1;

	while (! prio_queue.empty()) {
		int i = prio_queue.begin()->second;
		prio_queue.erase(prio_queue.begin());

		if (options.store_visited || options.print_progress) {
			n_visited++;

			if (options.store_visited) {
				options.visit_time[i] = n_visited;
			}

			if (options.print_progress) {
				if (double(std::clock() - last_time) > 0.3 * double(CLOCKS_PER_SEC)) {
					last_time = std::clock();
					double fraction_done = double(n_visited) / double(n);
					if (!first_print) {
						std::fprintf(stderr, "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
					}
					first_print = false;
					std::fprintf(stderr, "%7.3f%% visited... ", 100.0 * fraction_done);
					std::fflush(stderr);
				}
			}
		}

		//std::cerr << "Vistiting distance[" << i << "] == " << distance[i] << "\n";

		if (distance[i] >= infinity) {
			throw std::runtime_error("shortest_path: Path too long.");
		}

		// Is node i a goal node? If so, we are done.
		if (end_node == -1 && end_set.find(i) != end_set.end()) {
			// Store the shortest path.
			end_node = i;
			int j = i;
			path->clear();
			while (start_set.find(j) == start_set.end()) {
				path->push_back(j);
				j = previous[j];
			}
			path->push_back(j);
			// Store the shortest path from the start to
			// the end.
			std::reverse(path->begin(), path->end());
			

			if (!options.compute_all_distances) {
				// We are satisfied with the shortest path only.
				return distance[i];
			}
		}

		// Get all neighbors of node i using the oracle.
		neighbor_storage.clear();
		neighbors(i, &neighbor_storage);

		for (auto itr = neighbor_storage.begin(); itr != neighbor_storage.end(); ++itr) {
			// Debug check.
			if (itr->distance < 0) {
				throw std::runtime_error("shortest_path: Negative const encountered.");
			}
			// Index of the neighbor.
			int j = itr->destination;
			// Distance from the start to j via node i.
			double new_dist = distance[i] + itr->distance;
			// Previously known best distance from the start
			// to node j.
			double old_dist = distance[itr->destination];

			// Did we find a better path to j?
			if (new_dist < old_dist) {
				// Remove j from the queue (if present).
				queue_cost est;
				if (get_lower_bound) {
					est = estimation[j];
				}
				else {
					// If lower bounds are not available, the estimated
					// total distance is just the distance from the
					// start to j
					est = distance[j];
				}
				prio_queue.erase(std::make_pair(est, j));
				// Update the best distance to j.
				distance[j] = new_dist;
				previous[j] = i;
				// Get an estimation of the best distance.
				est = new_dist;
				if (get_lower_bound) {
					// If a lower bound function is available, the
					// estimated distance is the distance from the
					// start to j plus the lower bound from j to
					// the end.
					est += (*get_lower_bound)(j);
					estimation[j] = est;
				}
				// Add j with the new priority.
				prio_queue.insert(std::make_pair(est, j));
				//prio_queue.emplace(est, j);

				if (options.maximum_queue_size > 0 &&
				    prio_queue.size() > options.maximum_queue_size) {
					throw std::runtime_error("shortest_path: Maximum queue size reached.");
				}
			}
		}
	}

	if (!options.compute_all_distances) {
		// We should have reached the end set by now.
		if (end_node == -1) {
			throw std::runtime_error("shortest_path: No path found.");
		}
		else {
			// The algorithm should have terminated by now.
			throw std::runtime_error("shortest_path: Internal error.");
		}
	}

	// Clear some temporary storage.
	neighbor_storage.reserve(0);
	previous.reserve(0);

	// Copy the distance to the output.
	options.distance = distance;

	// If we had an end set, return the distance to it.
	if (end_node >= 0)
		return distance[end_node];
	else
		// No end set was provided, but this is not an error when
		// computing all distances.
		return -1.0;
}

double shortest_path(int n, const std::set<int>& start_set, const std::set<int>& end_set,
                     const std::function<void(int, std::vector<Neighbor>* neighbors)>& neighbors,
                     std::vector<int>* path, const std::function<double(int)>& get_lower_bound,
                     const ShortestPathOptions& options)
{
	return shortest_path(n, start_set, end_set, neighbors, path, &get_lower_bound, options);
}


#if 0
double shortest_path_memory_efficient(
	int n,
	const std::set<int>& start_set,
	const std::set<int>& end_set,
	const std::function<void(int, std::vector<Neighbor>* neighbors)>& neighbors,
	std::vector<int>* path,
	const ShortestPathOptions& options)
{
	// Initialize the heap structure. 
	// This will allocate (sizeof(float) + sizeof(int)) * n bytes.
	// This is usually 8*n bytes.
	heap::initialize(n);

	// The previous node in the shortest path from the start
	// to node i.
	// This is usually 4*n bytes 
	std::vector<int> previous(n, -1);

	// Additionally, each element actually put in the queue will
	// occupy 4 additional bytes.

	//
	// TOTAL: This function will allocate no more than 16*n + O(1) bytes.
	//

	std::vector<Neighbor> neighbor_storage;
	// Allocate storage for 100 neighbors. Each call to clear() will
	// not deallocate the storage.
	neighbor_storage.reserve(100);

	// Check and put the start set into the queue.
	if (start_set.size() == 0) {
		throw std::runtime_error("shortest_path: empty start set");
	}
	for (auto itr = start_set.begin(); itr != start_set.end(); ++itr) {
		if (*itr < 0 || *itr >= n) {
			throw std::runtime_error("shortest_path: Invalid start set.");
		}
		heap::insert(*itr, 0.0f);
	}
	// Check end_set.
	for (auto itr = end_set.begin(); itr != end_set.end(); ++itr) {
		if (*itr < 0 || *itr >= n) {
			throw std::runtime_error("shortest_path: Invalid end set.");
		}
	}

	int n_visited = 0;
	bool first_print = true;
	auto last_time = std::clock();

	while (! heap::is_empty()) {
		int i = heap::smallest_element();
		double current_distance = heap::smallest_priority();
		heap::pop();

		if (options.print_progress) {
			n_visited++;
			if (double(std::clock() - last_time) > 0.3 * double(CLOCKS_PER_SEC)) {
				last_time = std::clock();
				double fraction_done = double(n_visited) / double(n);
				if (!first_print) {
					std::fprintf(stderr, "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
				}
				first_print = false;
				std::fprintf(stderr, "%7.3f%% visited... ", 100.0 * fraction_done);
				std::fflush(stderr);
			}
		}

		//std::cerr << "Vistiting distance[" << i << "] == " << current_distance << std::endl;

		if (current_distance >= heap::infinity) {
			throw std::runtime_error("shortest_path: Path too long.");
		}

		// Is node i a goal node? If so, we are done.
		if (end_set.find(i) != end_set.end()) {
			// Store the shortest path.
			int j = i;
			path->clear();
			while (start_set.find(j) == start_set.end()) {
				path->push_back(j);
				j = previous[j];
			}
			path->push_back(j);
			// Store the shortest path from the start to
			// the end.
			std::reverse(path->begin(), path->end());
			return heap::priority[i];
		}

		// Get all neighbors of node i using the oracle.
		neighbor_storage.clear();
		neighbors(i, &neighbor_storage);

		for (auto itr = neighbor_storage.begin(); itr != neighbor_storage.end(); ++itr) {
			// Debug check.
			if (itr->distance < 0) {
				throw std::runtime_error("shortest_path: Negative const encountered.");
			}
			// Index of the neighbor.
			int j = itr->destination;
			// Distance from the start to j via node i.
			double new_dist = current_distance + itr->distance;
			// Previously known best distance from the start
			// to node j.
			double old_dist = heap::priority[j];

			// Did we find a better path to j?
			if (new_dist < old_dist) {
				previous[j] = i;
				// Add j with the new priority.
				heap::decrease_priority(j, new_dist);
			}
		}
	}

	throw std::runtime_error("shortest_path: No path found.");
	return 0.0;
}
#endif

}  // namespace curve_extraction
