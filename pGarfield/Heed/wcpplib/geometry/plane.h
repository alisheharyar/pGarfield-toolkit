#ifndef PLANE_H
#define PLANE_H

/* Copyright (c) 2000 Igor B. Smirnov

The file can be used, copied, modified, and distributed
according to the terms of GNU Lesser General Public License version 2.1
as published by the Free Software Foundation,
and provided that the above copyright notice, this permission notice,
and notices about any modifications of the original text
appear in all copies and in supporting documentation.
The file is provided "as is" without express or implied warranty.
*/

#include "wcpplib/geometry/vec.h"
#include "wcpplib/geometry/straight.h"

namespace Heed {

class polyline;

/// Plane, defined by defined by a point and a vector normal to the plane.

class plane : public absref {
 protected:
  /// Origin point, pivot.
  point piv;  
  /// Direction of normal, unit vector.
  vec dir;    

 public:
  point Gpiv() const { return piv; }
  vec Gdir() const { return dir; }

 protected:
  virtual absref_transmit get_components() override;
  static absref absref::* aref[2];

 public:
  plane() : piv(), dir() {}
  plane(const point& fpiv, const vec& fdir) : piv(fpiv), dir(unit_vec(fdir)) {}
  plane(const straight& sl, const point& pt);
  plane(const straight& sl1, const straight& sl2, vfloat prec);
  // Good if lines are crossed or if they are different not crossed parallel.
  // Otherwise vecerror != 0
  // Prec is used for crossing of lines.

  /// Copy constructor.
  plane(const plane& p) : piv(p.piv), dir(p.dir) {}
  /// Copy assignment operator.
  plane& operator=(const plane& fpl) {
    piv = fpl.piv;
    dir = fpl.dir;
    return *this;
  }

  friend int operator==(const plane& pl1, const plane& pl2);
  friend int operator!=(const plane& pl1, const plane& pl2) {
    return pl1 == pl2 ? 0 : 1;
  }
  friend bool apeq(const plane& pl1, const plane& pl2, vfloat prec);

  /// Return 1 if a point is in the plane (within precision prec).
  int check_point_in(const point& fp, vfloat prec) const;

  /// Figure out whether a straight line crosses the plane 
  /// and return the intersection point if it does.
  /// vecerror = 2: line is parallel to the plane.
  /// vecerror = 3: line is in the plane.
  point cross(const straight& sl) const;
  /// Determine the intersection with another plane.
  /// vecerror = 2: planes are parallel.
  /// vecerror = 3: planes are identical.
  straight cross(const plane& sl) const;

  int cross(const polyline& pll, point* crpt, int& qcrpt, polyline* crpll,
            int& qcrpll, vfloat prec) const;

  vfloat distance(const point& fpt) const;
  friend std::ostream& operator<<(std::ostream& file, const plane& s);
};

std::ostream& operator<<(std::ostream& file, const plane& s);
}

#endif
