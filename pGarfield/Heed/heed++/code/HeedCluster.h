#ifndef HEEDCLUSTER_H
#define HEEDCLUSTER_H

#include "wcpplib/geometry/vec.h"
#include "wcpplib/geometry/volume.h"

namespace Heed {

/// Cluster.
class HeedCluster {
 public:
  HeedCluster() = default;
  HeedCluster(double ftransferred_energy, const point& fpt,
              long fnatom, long fnshell)
      : transferred_energy(ftransferred_energy),
        pt(fpt),
        natom(fnatom),
        nshell(fnshell) {}
  /// Energy transfer in internal units.
  double transferred_energy = 0.;
  /// Coordinates in the global frame.
  point pt;
  long natom = 0;
  long nshell = 0;
  void print(std::ostream& file, int l) const;
};
}

#endif
