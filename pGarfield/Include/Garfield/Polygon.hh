#ifndef G_POLYGON_H
#define G_POLYGON_H

#include <vector>

namespace Garfield {

namespace Polygon {

/// Determine whether the point (x, y) is located inside of the
/// polygon (xpl, ypl).
void Inside(const std::vector<double>& xpl, const std::vector<double>& ypl,
            const double x, const double y, bool& inside, bool& edge);

/// Determine the (signed) area of a polygon.
double Area(const std::vector<double>& xp, const std::vector<double>& yp);

/// Check whether a set of points builds a non-trivial polygon.
bool NonTrivial(const std::vector<double>& xp, const std::vector<double>& yp);

/// Try to eliminate "butterflies" (crossing of two adjacent segments  
/// of a polygon), by point exchanges.
void EliminateButterflies(std::vector<double>& xp, std::vector<double>& yp,
                          std::vector<double>& zp);
}

}

#endif
