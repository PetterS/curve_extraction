#pragma once
#include <assert.h>
#include "Boundary_points.h"
#include "Depth_corners.h"

class Geodesic_length
{
  public:
    Geodesic_length (
      const matrix<double>& data,
      const vector<double>& voxel_dimensions,
      double penalty)
      : data(data),
        vd(voxel_dimensions),
        penalty(penalty),
        data_depdent(true),
        boundary_points(voxel_dimensions, data.M, data.N),
        euclidean_length(data, voxel_dimensions, penalty),
        depth_corners(voxel_dimensions,data)
   {
   }

    inline double sqr(double x)
    {
      return x*x;
    }

    double operator () (double x0,double y0,double z0, double x1, double y1, double z1)
    {
      using std::abs;

      // For fixed connectivity relative boundary_points can be cached.
      std::vector<boundary_point> points = boundary_points(x0,y0,x1,y1);
      double cost = 0;

      auto prev = points.begin();
      auto next = prev;
      next++;

      double i00,i01,i10,i11;
      double dx,dy;
      double d11,d10,d01;
      double a,b,c,dsqr,f,g;
      double p1,p2,p3;
      
      for (; next != points.end(); prev++,next++)
      {
        std::tie(i00,i01,i10,i11) = depth_corners(*prev);

        x0 = prev->first;
        x1 = next->first;

        y0 = prev->second;
        y1 = next->second;

        dx = (x0-x1)*vd[0];
        dy = (y0-y1)*vd[1];

        // Fractional part
        x1 = (x1-std::floor(x0))*vd[0];
        y1 = (y1-std::floor(y0))*vd[1];

        x0 = (x0-std::floor(x0))*vd[0];
        y0 = (y0-std::floor(y0))*vd[1];

        d11 = (i00+i11-i10-i01)/(vd[0]*vd[1]);
        d10 = (i10-i00)/(vd[0]*vd[1]);
        d01 = (i01-i00)/(vd[0]*vd[1]);

        double tolerance = 1e-8;
        if ( (abs(dx) < tolerance) || (abs(dy) < tolerance) || (abs(d11) < tolerance) )
        {
          cost += sqrt(sqr(dx)+sqr(dy)+sqr(dx*d10+dy*d01));
        } else
        {
          a = 2*d11*dx*dy;
          b = (y0*dx+x0*dy)/(2*dx*dy) + (d01*dy+d10*dx)/a;
          c = (sqr(dx)+sqr(dy))/(sqr(a));
          f = sqr(b)+c;
          g = 1+2*b+f;

          p1 = abs(a)/2;
          p2 = c*(log ( abs( (b+1+sqrt(g))/(b+sqrt(f)) ) ) );
          p3 = (b+1)*sqrt(g)-b*sqrt(f);

          cost += p1*(p2+p3);
        }
      }

      return penalty*cost;
    }

  // Data
  bool data_depdent;
  const matrix<double> data;
  vector<double> vd;
  double penalty;

  // Functors
  Boundary_points boundary_points;
  Euclidean_length euclidean_length;
  Depth_corners depth_corners;
};