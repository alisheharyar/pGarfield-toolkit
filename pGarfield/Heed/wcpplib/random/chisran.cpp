#include "wcpplib/random/chisran.h"
#include "wcpplib/util/FunNameStack.h"

// I. B. Smirnov, 2003.

namespace Heed {

double chispre(std::vector<double> &f, int s_allow_zero_f) {
  mfunnamep("double chispre(vector<double>& f, int s_allow_zero_f)");
  const size_t q = f.size();
  check_econd11(q, <= 0, mcerr);
  double r = 0;
  for (size_t i = 0; i < q; ++i) {
    if (s_allow_zero_f == 0) {
      check_econd11a(f[i], < 0.0, "i=" << i << '\n', mcerr);
    } else {
      if (f[i] < 0.0) {
        mcout << "Warning: f[i] < 0.0 in chispre\n";
        Iprint2n(mcout, i, f[i]);
        f[i] = 0.0;
      }
    }
    r += f[i];
    f[i] = r;
  }
  check_econd11(r, <= 0, mcerr);
  const double scale = 1. / r;
  for (size_t i = 0; i < q; ++i) f[i] *= scale;
  return r;
}

double chisran(double flat_random_number, const std::vector<double> &f) {
  mfunnamep("double chisran(double flat_random_number, vector<double>& f)");
  const long q = f.size();
  check_econd11(q, <= 0, mcerr);
  check_econd21(flat_random_number, < 0.0 ||, > 1.0, mcerr);
  if (flat_random_number == 0.0) {
    for (long n = 0; n < q; ++n) {
      if (f[n] > 0.0) return double(n);
    }
  } else {
    if (flat_random_number == 1.0) {
      for (long n = q - 1; n >= 0; n--) {
        if (f[n] < 1.0) return double(n + 1);
      }
    } else {
      if (flat_random_number <= f[0]) {
        return flat_random_number / f[0];
      } 
      long nl = 0;
      long nr = q - 1;
      long nc;
      while (nr - nl > 1) {
        nc = (nr + nl) / 2;
        if (flat_random_number < f[nc]) {
          nr = nc;
        } else {
          nl = nc;
        }
      }
      const double xl = double(nl + 1);
      const double xr = double(nr + 1);
      const double yl = f[nl];
      const double yr = f[nr];
      const double a = (xr - xl) / (yr - yl);
      const double b = xl;
      // Iprint3n(mcout, nl, nr, nc);
      // Iprint2n(mcout, xl, xr);
      // Iprint2n(mcout, yl, yr);
      // Iprint2n(mcout, a, b);
      return a * (flat_random_number - yl) + b;
    }
  }
  funnw.ehdr(mcerr);
  mcerr << "should never happen\n";
  spexit(mcerr);
  return 0.0;
}
}
