// Petter Strandmark 2013.
#ifndef CURVE_EXTRACTION_CURVATURE_H
#define CURVE_EXTRACTION_CURVATURE_H

namespace curve_extraction {

extern int curvature_cache_hits;
extern int curvature_cache_misses;
template<typename R>
R compute_curvature(R x1, R y1, R z1,
                    R x2, R y2, R z2,
                    R x3, R y3, R z3,
                    R power = 2.0,
                    bool writable_cache = true,
                    int n_approximation_points = 200
                    );

extern int torsion_cache_hits;
extern int torsion_cache_misses;

template<typename R>
R compute_torsion(R x1, R y1, R z1,
                  R x2, R y2, R z2,
                  R x3, R y3, R z3,
                  R x4, R y4, R z4,
                  R power = 2.0,
                  bool writable_cache = true,
                  int n_approximation_points = 200);

}  // namespace curve_extraction

#endif
