// Petter Strandmark 2013.
#ifndef CURVE_EXTRACTION_DATA_TERM_H
#define CURVE_EXTRACTION_DATA_TERM_H

#include <vector>

namespace curve_extraction {

class PieceWiseConstant
{
public:
	PieceWiseConstant(const double * unary,
	                  int M, int N, int O,
	                  const std::vector<double>& voxeldimensions);

	template<typename R>
	R evaluate(R x, R y, R z = 0.0) const;

	template<typename R>
	R evaluate_line_integral(R x1, R y1, R z1,
	                         R x2, R y2, R z2) const;
private:
	int xyz_to_ind(double x, double y, double z) const;
	template<typename R> 	bool inside_volume(R x, R y, R s) const;
	int M, N, O;
	const double* unary;
	const std::vector<double> voxeldimensions;
};

class TriLinear
{
public:
	TriLinear(const double * unary,
	          int M, int N, int O,
	          const std::vector<double>& voxeldimensions);

	template<typename R>
	R evaluate(R x, R y, R z = 0.0) const;

	template<typename R>
	R evaluate_line_integral(R x1, R y1, R z1,
	                         R x2, R y2, R z2) const;
private:
	int xyz_to_ind(double x, double y, double z) const;
	int M, N, O;
	const double* unary;
	const std::vector<double> voxeldimensions;
};

}  // namespace curve_extraction

#endif
