#include <cmath>

#include "wcpplib/random/pois.h"
#include "wcpplib/random/rnorm.h"
#include "wcpplib/random/ranluxint.h"

namespace Heed {

long pois(const double amu) {
  // POISSON GENERATOR
  // CODED FROM LOS ALAMOS REPORT      LA-5061-MS
  // PROB(N)=EXP(-AMU)*AMU**N/FACT(N)
  // WHERE FACT(N) STANDS FOR FACTORIAL OF N

  if (amu <= 0.) return 0;
  if (amu > 100.) {
    return static_cast<long>(rnorm_improved() * sqrt(amu) + amu + 0.5);
  }
  double expma = exp(-amu);
  double pir = 1.;
  long n = -1;
  while (1) {
    ++n;
    pir *= SRANLUX();
    if (pir <= expma) break;
  }
  return n;
}
}
