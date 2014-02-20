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
        depth_corners(data)
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

        i00 *= vd[2];
        i01 *= vd[2];
        i10 *= vd[2];
        i11 *= vd[2];

        x0 = prev->first;
        x1 = next->first;

        y0 = prev->second;
        y1 = next->second;

        dx = (x0-x1)*vd[0];
        dy = (y0-y1)*vd[1];

        // Fractional part
        x0 = (x0-std::floor(x0))*vd[0];
        y0 = (y0-std::floor(y0))*vd[1];

        // See paper for explanation.
        d11 = (i00+i11-i10-i01)/(vd[0]*vd[1]);
        d10 = (y0*(i01-i11) + y1*(i10-i00) )/(vd[0]*vd[1]);
        d01 = (x0*(i10-i11) + x1*(i01-i11) )/(vd[0]*vd[1]);

        a = 2*d11*dx*dy;
        b = 2*i00-i10-i01;
        c = b/(2*d11);
        dsqr = (sqr(dx)+sqr(dy))/(sqr(a));
        f = sqr(c)+dsqr;
        g = 1+2*c+f;

        //
        if ((dx == 0) || (dy == 0) || (d11 == 0))
        {
          cost += sqrt(sqr(dx)+sqr(dy)+sqr(dx*d10+dy*d01));
        } else
        {
          p1 = abs(a)/2;
          p2 = dsqr*(log ( abs( (c+1+sqrt(g))/(c+sqrt(f)) ) ) );
          p3 = (c+1)*sqrt(g)-c*sqrt(f);

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