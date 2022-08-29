#include <iomanip>
#include "heed++/code/HeedCluster.h"

// 2003, I. Smirnov

namespace Heed {

void HeedCluster::print(std::ostream& file, int l) const {
  if (l <= 0) return;
  Ifile << "HeedCluster (l=" << l
        << "): transferred_energy=" << transferred_energy << " MeV\n";
  Ifile << "pt=" << pt << '\n';
  if (l > 1) {
    indn.n += 2;
    Ifile << "natom=" << natom << " nshell=" << nshell << '\n';
    indn.n -= 2;
  }
}
}
