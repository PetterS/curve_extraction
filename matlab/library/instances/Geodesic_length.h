#pragma once
#include <assert.h>
#include "Boundary_points.h"

class Geodesic_length
{
  public:
    Geodesic_length (
      const matrix<double>& data,
      const InstanceSettings& settings)
      : penalty(settings.penalty[0]),
        data_dependent(true),
        boundary_points_calculator(data, settings.voxel_dimensions)
   {}

    template<typename R>
    inline R sqr(R x) const
    {
      return x*x;
    }

    template<typename R>
    R operator()(const R* const point1, const R* const point2) const
    {
      using std::abs;

      // For fixed connectivity relative boundary_points can be cached.
      Boundary_points<R> points = boundary_points_calculator( point1[0],point1[1],
                                                              point2[0],point2[1]);
      R cost = 0;

      for (int i = 0; i < points.num_pairs(); i++)
      {
        R x0,y0, dx,dy;
        double d11,d10,d01;
        tie(x0,y0,dx,dy, d11,d10,d01) = points.get_pair(i);

        R tolerance = 1e-8;
        if ( (abs(dx) < tolerance) || (abs(dy) < tolerance) || (abs(d11) < tolerance) )
        {
          cost += sqrt(sqr(dx)+sqr(dy)+sqr(dx*d10+dy*d01));
        } else
        {
          R a = 2*d11*dx*dy;
          R b = (y0*dx+x0*dy)/(2*dx*dy) + (d10*dx+d01*dy)/a;
          R c = (sqr(dx)+sqr(dy))/(sqr(a));
          R f = sqr(b)+c;
          R g = 1+2*b+f;

          R p1 = abs(a)/2;
          R p2 = c*(log ( abs( (b+1+sqrt(g))/(b+sqrt(f)) ) ) );
          R p3 = (b+1)*sqrt(g)-b*sqrt(f);

          cost += p1*(p2+p3);
        }
      }

      return penalty*cost;
    }

  double penalty;
  bool data_dependent;
  Boundary_points_calculator boundary_points_calculator;
};