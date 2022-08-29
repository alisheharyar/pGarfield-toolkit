#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>

#include "Garfield/Numerics.hh"

namespace {

int deqnGen(const int n, std::vector<std::vector<double > >& a,
            std::vector<double>& b) {

  std::vector<int> ir(n, 0);
  double det = 0.;
  int ifail = 0, jfail = 0;
  Garfield::Numerics::CERNLIB::dfact(n, a, ir, ifail, det, jfail); 
  if (ifail != 0) return ifail;
  Garfield::Numerics::CERNLIB::dfeqn(n, a, ir, b);
  return 0;
}

/// Epsilon algorithm.
/// Determines the limit of a given sequence of approximations, 
/// by means of the epsilon algorithm of P. Wynn. 
/// An estimate of the absolute error is also given.
/// The condensed epsilon table is computed. Only those elements needed
/// for the computation of the next diagonal are preserved.
/// \param epstab elements of the two lower diagonals of the triangular
///               epsilon table. The elements are numbered starting at the
///               right-hand corner of the triangle.
/// \param n size of the epsilon table.
/// \param result resulting approximation to the integral.
/// \param abserr estimate of the absolute error computed from
///               result and the three previous results.
/// \param lastRes last three results.
/// \param nres number of calls to the function.
void qelg(unsigned int& n, std::array<double, 52>& epstab, 
          double& result, double& abserr, 
          std::array<double, 3>& lastRes, unsigned int& nres) {

  constexpr double eps = std::numeric_limits<double>::epsilon();

  ++nres;
  abserr = std::numeric_limits<double>::max(); 
  result = epstab[n - 1];
  if (n < 3) {
    abserr = std::max(abserr, 50. * eps * std::abs(result));
    return;
  }
  epstab[n + 1] = epstab[n - 1];
  epstab[n - 1] = std::numeric_limits<double>::max();
  // Number of elements to be computed in the new diagonal.
  const unsigned int nnew = (n - 1) / 2;
  const unsigned int nold = n;
  unsigned int k = n;
  for (unsigned int i = 1; i <= nnew; ++i) {
    double res = epstab[k + 1];
    // e0 - e3 are the four elements on which the computation of a new
    // element in the epsilon table is based.
    //                 e0
    //           e3    e1    new
    //                 e2
    const double e0 = epstab[k - 3];
    const double e1 = epstab[k - 2];
    const double e2 = res;
    const double delta2 = e2 - e1;
    const double err2 = std::abs(delta2);
    const double tol2 = std::max(std::abs(e2), std::abs(e1)) * eps;
    const double delta3 = e1 - e0;
    const double err3 = std::abs(delta3);
    const double tol3 = std::max(std::abs(e1), std::abs(e0)) * eps;
    if (err2 <= tol2 && err3 <= tol3) {
      // If e0, e1 and e2 are equal to within machine accuracy,
      // convergence is assumed.
      result = res;
      abserr = std::max(err2 + err3, 50. * eps * std::abs(result));
      return;
    }
    const double e3 = epstab[k - 1];
    epstab[k - 1] = e1;
    const double delta1 = e1 - e3;
    const double err1 = std::abs(delta1);
    const double tol1 = std::max(std::abs(e1), std::abs(e3)) * eps;
    // If two elements are very close to each other, omit
    // a part of the table by adjusting the value of n
    if (err1 <= tol1 || err2 <= tol2 || err3 <= tol3) {
      n = i + i - 1;
      break;
    }
    const double ss = 1. / delta1 + 1. / delta2 - 1./ delta3;
    // Test to detect irregular behaviour in the table, and
    // eventually omit a part of the table adjusting the value of n.
    if (std::abs(ss * e1) <= 1.e-4) {
      n = i + i - 1;
      break;
    }
    // Compute a new element and eventually adjust the value of result.
    res = e1 + 1. / ss;
    epstab[k - 1] = res;
    k -= 2;
    const double error = err2 + std::abs(res - e2) + err3;
    if (error <= abserr) {
      abserr = error;
      result = res;
    }
  }
  // Shift the table.
  constexpr unsigned int limexp = 50;
  if (n == limexp) n = 2 * (limexp / 2) - 1;
  unsigned int ib = (nold % 2 == 0) ? 1 : 0;
  for (unsigned int i = 0; i <= nnew; ++i) {
    epstab[ib] = epstab[ib + 2];
    ib += 2;
  } 
  if (nold != n) {
    for (unsigned int i = 0; i < n; ++i) {
      epstab[i] = epstab[nold - n + i];
    }
  }
  if (nres >= 4) {
    // Compute error estimate.
    abserr = std::abs(result - lastRes[2]) + 
             std::abs(result - lastRes[1]) +
             std::abs(result - lastRes[0]);
    lastRes[0] = lastRes[1];
    lastRes[1] = lastRes[2];
    lastRes[2] = result;
  } else {
    lastRes[nres - 1] = result;
    abserr = std::numeric_limits<double>::max();
  }
  abserr = std::max(abserr, 50. * eps * std::abs(result));
}

}

namespace Garfield {

namespace Numerics {

namespace QUADPACK {

void qagi(std::function<double(double)> f, double bound, const int inf, 
          const double epsabs, const double epsrel, 
          double& result, double& abserr, unsigned int& status) {

  status = 0;
  result = 0.;
  abserr = 0.;

  // Test on validity of parameters.
  constexpr double eps = std::numeric_limits<double>::epsilon();
  if (epsabs <= 0. && epsrel < std::max(50. * eps, 0.5e-28)) {
    status = 6;
    return;
  }
  // First approximation to the integral
  if (inf == 2) bound = 0.;
  double resabs0 = 0., resasc0 = 0.;
  qk15i(f, bound, inf, 0., 1., result, abserr, resabs0, resasc0);

  // Calculate the requested accuracy.
  double tol = std::max(epsabs, epsrel * std::abs(result));
  // Test on accuracy.
  if (abserr <= 100. * eps * resabs0 && abserr > tol) {
    // Roundoff error at the first attempt.
    status = 2;
    return;
  }
  // Test if the first approximation was good enough.
  if ((abserr <= tol && abserr != resasc0) || abserr == 0.) return;

  struct Interval {
    double a; ///< Left end point. 
    double b; ///< Right end point.
    double r; ///< Approximation to the integral over this interval.
    double e; ///< Error estimate.
  };
  std::vector<Interval> intervals(1);
  intervals[0].a = 0.;
  intervals[0].b = 1.;
  intervals[0].r = result;
  intervals[0].e = abserr;
  constexpr unsigned int nMaxIntervals = 500;
  unsigned int nIntervals = 1;
  // Interval to be bisected.
  auto it = intervals.begin();
  size_t nrmax = 0;

  // Initialize the epsilon table.
  std::array<double, 52> epstab;
  epstab[0] = result;
  // Count the number of elements currently in the epsilon table.
  unsigned int nEps = 2;
  // Keep track of the last three results.
  std::array<double, 3> lastRes = {0., 0., 0.};
  // Count the number of calls to the epsilon extrapolation function.
  unsigned int nRes = 0;
  // Flag denoting that we are attempting to perform extrapolation.
  bool extrap = false;
  // Flag indicating that extrapolation is no longer allowed.
  bool noext = false;
  unsigned int ktmin = 0;

  // Initialize the sum of the integrals over the subintervals.
  double area = result;
  // Initialize the sum of the errors over the subintervals.
  double errSum = abserr;
  // Length of the smallest interval considered up now, multiplied by 1.5. 
  double small = 0.375; 
  // Sum of the errors over the intervals larger than the smallest interval 
  // considered up to now.
  double errLarge = 0.;
  double errTest = 0.;
  // Error estimate of the interval with the largest error estimate.
  double errMax = abserr;
  abserr = std::numeric_limits<double>::max();

  // Count roundoff errors.
  std::array<unsigned int, 3> nRoundOff = {0, 0, 0};
  bool roundOffErrors = false;
  double correc = 0.;

  // Set flag whether the integrand is positive.
  const bool pos = (std::abs(result) >= (1. - 50. * eps) * resabs0);

  // Main loop.
  bool dosum = false;
  for (nIntervals = 2; nIntervals <= nMaxIntervals; ++nIntervals) {
    // Bisect the subinterval.
    const double a1 = (*it).a;
    const double b2 = (*it).b;
    const double b1 = 0.5 * (a1 + b2);
    const double a2 = b1;
    // Save the error on the interval before doing the subdivision.
    const double errLast = (*it).e;
    double area1 = 0., err1 = 0., resabs1 = 0., resasc1 = 0.;
    qk15i(f, bound, inf, a1, b1, area1, err1, resabs1, resasc1);
    double area2 = 0., err2 = 0., resabs2 = 0., resasc2 = 0.;
    qk15i(f, bound, inf, a2, b2, area2, err2, resabs2, resasc2);
    // Improve previous approximations to integral and error 
    // and test for accuracy.
    const double area12 = area1 + area2;
    const double err12 = err1 + err2;
    errSum += err12 - errMax;
    area += area12 - (*it).r;
    if (resasc1 != err1 && resasc2 != err2) {
      if (std::abs((*it).r - area12) <= 1.e-5 * std::abs(area12) && 
          err12 >= 0.99 * errMax) {
        if (extrap) {
          ++nRoundOff[1];
        } else {
          ++nRoundOff[0];
        }
      }   
      if (nIntervals > 10 && err12 > errMax) ++nRoundOff[2];
    }
    tol = std::max(epsabs, epsrel * std::abs(area));
    // Test for roundoff error and eventually set error flag.
    if (nRoundOff[0] + nRoundOff[1] >= 10 || nRoundOff[2] >= 20) status = 2;
    if (nRoundOff[1] >= 5) roundOffErrors = true;
    // Set error flag in the case that the number of subintervals equals limit.
    if (nIntervals == nMaxIntervals) status = 1;
    // Set error flag in the case of bad integrand behaviour
    // at some points of the integration range.
    constexpr double tol1 = 1. + 100. * eps;
    constexpr double tol2 = 1000. * std::numeric_limits<double>::min();
    if (std::max(std::abs(a1), std::abs(b2)) <= tol1 * (std::abs(a2) + tol2)) {
      status = 4;
    }
    // Append the newly-created intervals to the list.
    if (err2 > err1) {
      (*it).a = a2;
      (*it).r = area2;
      (*it).e = err2;
      Interval interval;
      interval.a = a1;
      interval.b = b1;
      interval.r = area1;
      interval.e = err1;
      intervals.push_back(std::move(interval)); 
    } else {
      (*it).b = b1;
      (*it).r = area1;
      (*it).e = err1;
      Interval interval;
      interval.a = a2;
      interval.b = b2;
      interval.r = area2;
      interval.e = err2;
      intervals.push_back(std::move(interval)); 
    }
    // Sort the intervals in descending order by error estimate. 
    std::sort(intervals.begin(), intervals.end(),
             [](const Interval& lhs, const Interval& rhs) {
               return (lhs.e > rhs.e);
             });
    // Select the subinterval to be bisected next.
    it = intervals.begin() + nrmax;
    errMax = (*it).e;
    if (errSum <= tol) {
      dosum = true;
      break;
    }
    if (status != 0) break;
    if (nIntervals == 2) {
      errLarge = errSum;
      errTest = tol;
      epstab[1] = area;
      continue;
    }
    if (noext) continue;
    errLarge -= errLast;
    if (std::abs(b1 - a1) > small) errLarge += err12;
    if (!extrap) {
      // Test whether the interval to be bisected next is the smallest one.
      if (std::abs((*it).b - (*it).a) > small) continue;
      extrap = true;
      nrmax = 1; 
    }
    // The smallest interval has the largest error.
    // Before bisecting decrease the sum of the errors over the
    // larger intervals (errLarge) and perform extrapolation.
    if (!roundOffErrors && errLarge > errTest) {
      const size_t k0 = nrmax;
      size_t k1 = nIntervals;
      if (nIntervals > (2 + nMaxIntervals / 2)) {
        k1 = nMaxIntervals + 3 - nIntervals;
      }
      bool found = false;
      for (unsigned int k = k0; k < k1; ++k) {
        it = intervals.begin() + nrmax;
        errMax = (*it).e;
        if (std::abs((*it).b - (*it).a) > small) {
          found = true;
          break;
        }
        ++nrmax;
      }
      if (found) continue;
    }
    // Perform extrapolation.
    epstab[nEps] = area;
    ++nEps;
    double resExtr = 0., errExtr = 0.;
    qelg(nEps, epstab, resExtr, errExtr, lastRes, nRes);
    ++ktmin;
    if (ktmin > 5 && abserr < 1.e-3 * errSum) status = 5;
    if (errExtr < abserr) {
      ktmin = 0;
      abserr = errExtr;
      result = resExtr;
      correc = errLarge;
      errTest = std::max(epsabs, epsrel * std::abs(resExtr));
      if (abserr <= errTest) break;
    }
    // Prepare bisection of the smallest interval.
    if (nEps == 1) noext = true;
    if (status == 5) break;
    it = intervals.begin();
    errMax = (*it).e;
    nrmax = 0;
    extrap = false;
    small *= 0.5;
    errLarge = errSum;
  }
  // Set final result and error estimate.
  if (abserr == std::numeric_limits<double>::max()) dosum = true;
  if (!dosum) {
    if ((status != 0 || roundOffErrors)) {
      if (roundOffErrors) {
        abserr += correc;
        if (status == 0) status = 3;
      }
      if (result != 0. && area != 0.) {
        if (abserr / std::abs(result) > errSum / std::abs(area)) dosum = true;
      } else {
        if (abserr > errSum) {
          dosum = true;
        } else if (area == 0.) {
          if (status > 2) --status;
          return;
        }
      }
    }
    // Test on divergence
    if (!dosum) {
      if (!pos && std::max(std::abs(result), std::abs(area)) <= resabs0 * 0.01) {
        if (status > 2) --status;
        return;
      }
      const double r = result / area;
      if (r < 0.01 || r > 100. || errSum > std::abs(area)) {
        status = 5;
        return;
      }
    }
  } else {
    // Compute global integral sum.
    result = 0.;
    for (const auto& interval : intervals) result += interval.r;
    abserr = errSum;
    if (status > 2) --status;
  }
}

void qk15i(std::function<double(double)> f, double bound, const int inf,
           const double a, const double b, double& result, double& abserr,
           double& resabs, double& resasc) {

  // The abscissae and weights are supplied for the interval (-1, 1).
  // Because of symmetry only the positive abscissae and
  // their corresponding weights are given.

  // Weights of the 7-point Gauss rule.
  constexpr double wg[8] = {0., 0.129484966168869693270611432679082,
                            0., 0.279705391489276667901467771423780,
                            0., 0.381830050505118944950369775488975,
                            0., 0.417959183673469387755102040816327};
  // Abscissae of the 15-point Kronrod rule.
  constexpr double xgk[8] = {
    0.991455371120812639206854697526329, 
    0.949107912342758524526189684047851,
    0.864864423359769072789712788640926, 
    0.741531185599394439863864773280788,
    0.586087235467691130294144838258730, 
    0.405845151377397166906606412076961,
    0.207784955007898467600689403773245, 
    0.};
  // Weights of the 15-point Kronrod rule.
  constexpr double wgk[8] = {
    0.022935322010529224963732008058970,
    0.063092092629978553290700663189204,
    0.104790010322250183839876322541518,
    0.140653259715525918745189590510238,
    0.169004726639267902826583426598550,
    0.190350578064785409913256402421014,
    0.204432940075298892414161999234649,
    0.209482141084727828012999174891714};

  const int dinf = std::min(1, inf);

  // Mid point of the interval.
  const double xc = 0.5 * (a + b);
  // Half-length of the interval.
  const double h = 0.5 * (b - a);
  // Transformed abscissa.
  const double tc = bound + dinf * (1. - xc) / xc;
  double fc = f(tc);
  if (inf == 2) fc += f(-tc);
  fc = (fc / xc) / xc;
  // Result of the 7-point Gauss formula.
  double resg = wg[7] * fc;
  // Result of the 15-point Kronrod formula.
  double resk = wgk[7] * fc;
  resabs = std::abs(resk);
  std::array<double, 7> fv1, fv2;
  for (unsigned int j = 0; j < 7; ++j) {
    const double x = h * xgk[j];
    const double x1 = xc - x;
    const double x2 = xc + x;
    const double t1 = bound + dinf * (1. - x1) / x1;
    const double t2 = bound + dinf * (1. - x2) / x2;
    double y1 = f(t1);
    double y2 = f(t2);
    if (inf == 2) {
      y1 += f(-t1);
      y2 += f(-t2);
    }
    y1 = (y1 / x1) / x1;
    y2 = (y2 / x2) / x2;
    fv1[j] = y1;
    fv2[j] = y2;
    const double fsum = y1 + y2;
    resg += wg[j] * fsum;
    resk += wgk[j] * fsum;
    resabs += wgk[j] * (std::abs(y1) + std::abs(y2));
  }
  // Approximation to the mean value of the transformed integrand over (a,b).
  const double reskh = resk * 0.5;
  resasc = wgk[7] * std::abs(fc - reskh);
  for (unsigned int j = 0; j < 7; ++j) {
    resasc += wgk[j] * (std::abs(fv1[j] - reskh) + std::abs(fv2[j] - reskh));
  }
  result = resk * h;
  resasc *= h;
  resabs *= h;
  abserr = std::abs((resk - resg) * h);
  if (resasc != 0. && abserr != 0.) {
    abserr = resasc * std::min(1., pow(200. * abserr / resasc, 1.5));
  }
  constexpr double eps = 50. * std::numeric_limits<double>::epsilon();
  if (resabs > std::numeric_limits<double>::min() / eps) {
    abserr = std::max(eps * resabs, abserr);
  }
}

void qk15(std::function<double(double)> f, const double a, const double b,
          double& result, double& abserr, double& resabs, double& resasc) {

  // Gauss quadrature weights and Kronron quadrature abscissae and weights
  // as evaluated with 80 decimal digit arithmetic by L. W. Fullerton,
  // Bell labs, Nov. 1981.

  // Weights of the 7-point Gauss rule.
  constexpr double wg[4] = {
    0.129484966168869693270611432679082,
    0.279705391489276667901467771423780,
    0.381830050505118944950369775488975,
    0.417959183673469387755102040816327};
  // Abscissae of the 15-point Kronrod rule.
  constexpr double xgk[8] = {
    0.991455371120812639206854697526329,
    0.949107912342758524526189684047851,
    0.864864423359769072789712788640926,
    0.741531185599394439863864773280788,
    0.586087235467691130294144838258730,
    0.405845151377397166906606412076961,
    0.207784955007898467600689403773245,
    0.};
  // Weights of the 15-point Kronrod rule.
  constexpr double wgk[8] = {
    0.022935322010529224963732008058970,
    0.063092092629978553290700663189204,
    0.104790010322250183839876322541518,
    0.140653259715525918745189590510238,
    0.169004726639267902826583426598550,
    0.190350578064785409913256402421014,
    0.204432940075298892414161999234649,
    0.209482141084727828012999174891714};

  // Mid point of the interval.
  const double xc = 0.5 * (a + b);
  // Half-length of the interval.
  const double h = 0.5 * (b - a);
  const double dh = std::abs(h);
  // Compute the 15-point Kronrod approximation to the integral, 
  // and estimate the absolute error.
  const double fc = f(xc);
  // Result of the 7-point Gauss formula.
  double resg = fc * wg[3];
  // Result of the 15-point Kronrod formula.
  double resk = fc * wgk[7];
  resabs = std::abs(resk);
  std::array<double, 7> fv1, fv2;
  for (unsigned int j = 0; j < 3; ++j) {
    const unsigned int k = j * 2 + 1;
    const double x = h * xgk[k];
    double y1 = f(xc - x);
    double y2 = f(xc + x);
    fv1[k] = y1;
    fv2[k] = y2;
    const double fsum = y1 + y2;
    resg += wg[j] * fsum;
    resk += wgk[k] * fsum;
    resabs += wgk[k] * (std::abs(y1) + std::abs(y2));
  }
  for (unsigned int j = 0; j < 4; ++j) {
    const unsigned int k = j * 2;
    const double x = h * xgk[k];
    const double y1 = f(xc - x);
    const double y2 = f(xc + x);
    fv1[k] = y1;
    fv2[k] = y2;
    const double fsum = y1 + y2;
    resk += wgk[k] * fsum;
    resabs += wgk[k] * (std::abs(y1) + std::abs(y2));
  }
  // Approximation to the mean value of f over (a,b), i.e. to i/(b-a).
  const double reskh = resk * 0.5;
  resasc = wgk[7] * std::abs(fc - reskh);
  for (unsigned int j = 0; j < 7; ++j) {
    resasc += wgk[j] * (std::abs(fv1[j] - reskh) + std::abs(fv2[j] - reskh));
  }
  result = resk * h;
  resabs *= dh;
  resasc *= dh;
  abserr = std::abs((resk - resg) * h);
  if (resasc != 0. && abserr != 0.) {
    abserr = resasc * std::min(1., pow(200. * abserr / resasc, 1.5));
  }
  constexpr double eps = 50. * std::numeric_limits<double>::epsilon();
  if (resabs > std::numeric_limits<double>::min() / eps) {
    abserr = std::max(eps * resabs, abserr);
  }
}

}

namespace CERNLIB {

int deqn(const int n, std::vector<std::vector<double> >& a,
         std::vector<double>& b) { 

  // REPLACES B BY THE SOLUTION X OF A*X=B, AFTER WHICH A IS UNDEFINED.

  if (n < 1) return 1;
  if (n == 1) {
    if (a[0][0] == 0.) return -1;
    const double s = 1. / a[0][0];
    b[0] *= s;
  } else if (n == 2) {
    // Cramer's rule.
    const double det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    if (det == 0.) return -1;
    const double s = 1. / det;
    const double b1 = b[0];
    b[0] = s * ( a[1][1] * b1 - a[0][1] * b[1]);
    b[1] = s * (-a[1][0] * b1 + a[0][0] * b[1]);
  } else if (n == 3) {
    // Factorize matrix A=L*U.
    // First pivot search.
    const double t1 = std::abs(a[0][0]);
    const double t2 = std::abs(a[1][0]);
    const double t3 = std::abs(a[2][0]);
    unsigned int m1 = 0, m2 = 0, m3 = 0;
    if (t1 < t2 && t3 < t2) {
      // Pivot is A21
      m1 = 1;
      m2 = 0;
      m3 = 2;
    } else if (t2 < t1 && t3 < t1) { 
      // Pivot is A11
      m1 = 0;
      m2 = 1;
      m3 = 2;
    } else {
      // Pivot is A31
      m1 = 2;
      m2 = 1;
      m3 = 0;
    }
    double temp = a[m1][0];
    if (temp == 0.) return deqnGen(n, a, b);
    const double l11 = 1. / temp;
    const double u12 = l11 * a[m1][1];
    const double u13 = l11 * a[m1][2];
    double l22 = a[m2][1] - a[m2][0] * u12;
    double l32 = a[m3][1] - a[m3][0] * u12;
    // Second pivot search.
    if (std::abs(l22) < std::abs(l32)) {
      std::swap(m2, m3);
      std::swap(l22, l32);
    }
    double l21 = a[m2][0];
    double l31 = a[m3][0];
    if (l22 == 0.) return deqnGen(n, a, b);
    l22 = 1. / l22;
    const double u23 = l22 * (a[m2][2] - l21 * u13);
    temp = a[m3][2] - l31 * u13 - l32 * u23;
    if (temp == 0.) return deqnGen(n, a, b);
    const double l33 = 1. / temp;
 
    // Solve L*Y=B and U*X=Y.
    const double y1 = l11 * b[m1];
    const double y2 = l22 * (b[m2] - l21 * y1);
    b[2] = l33 * (b[m3] - l31 * y1 - l32 * y2);
    b[1] = y2 - u23 * b[2];
    b[0] = y1 - u12 * b[1] - u13 * b[2];
  } else {
    // N > 3 cases. Factorize matrix and solve system.
    return deqnGen(n, a, b);
  }
  return 0;
}

void dfact(const int n, std::vector<std::vector<double> >& a,
           std::vector<int>& ir, int& ifail, double& det, int& jfail) {
  constexpr double g1 = 1.e-19;
  constexpr double g2 = 1.e-19;

  double t;
  int k;

  ifail = jfail = 0;

  int nxch = 0;
  det = 1.;

  for (int j = 1; j <= n; ++j) {
    k = j;
    double p = fabs(a[j - 1][j - 1]);
    if (j == n) {
      if (p <= 0.) {
        det = 0.;
        ifail = -1;
        jfail = 0;
        return;
      }
      det *= a[j - 1][j - 1];
      a[j - 1][j - 1] = 1. / a[j - 1][j - 1];
      t = fabs(det);
      if (t < g1) {
        det = 0.;
        if (jfail == 0) jfail = -1;
      } else if (t > g2) {
        det = 1.;
        if (jfail == 0) jfail = +1;
      }
      continue;
    }
    for (int i = j + 1; i <= n; ++i) {
      double q = std::abs(a[i - 1][j - 1]);
      if (q <= p) continue;
      k = i;
      p = q;
    }
    if (k != j) {
      for (int l = 1; l <= n; ++l) {
        std::swap(a[j - 1][l - 1], a[k - 1][l - 1]);
      }
      ++nxch;
      ir[nxch - 1] = j * 4096 + k;
    } else if (p <= 0.) {
      det = 0.;
      ifail = -1;
      jfail = 0;
      return;
    }
    det *= a[j - 1][j - 1];
    a[j - 1][j - 1] = 1. / a[j - 1][j - 1];
    t = fabs(det);
    if (t < g1) {
      det = 0.;
      if (jfail == 0) jfail = -1;
    } else if (t > g2) {
      det = 1.;
      if (jfail == 0) jfail = +1;
    }
    for (k = j + 1; k <= n; ++k) {
      double s11 = -a[j - 1][k - 1];
      double s12 = -a[k - 1][j];
      if (j == 1) {
        a[j - 1][k - 1] = -s11 * a[j - 1][j - 1];
        a[k - 1][j] = -(s12 + a[j - 1][j] * a[k - 1][j - 1]);
        continue;
      }
      for (int i = 1; i <= j - 1; ++i) {
        s11 += a[i - 1][k - 1] * a[j - 1][i - 1];
        s12 += a[i - 1][j] * a[k - 1][i - 1];
      }
      a[j - 1][k - 1] = -s11 * a[j - 1][j - 1];
      a[k - 1][j] = -a[j - 1][j] * a[k - 1][j - 1] - s12;
    }
  }

  if (nxch % 2 != 0) det = -det;
  if (jfail != 0) det = 0.;
  ir[n - 1] = nxch;
}

void dfeqn(const int n, std::vector<std::vector<double> >& a,
           std::vector<int>& ir, std::vector<double>& b) {
  if (n <= 0) return;

  int nxch = ir[n - 1];
  if (nxch != 0) {
    for (int m = 1; m <= nxch; ++m) {
      const int ij = ir[m - 1];
      const int i = ij / 4096;
      const int j = ij % 4096;
      std::swap(b[i - 1], b[j - 1]);
    }
  }

  b[0] *= a[0][0];
  if (n == 1) return;

  for (int i = 2; i <= n; ++i) {
    double s21 = -b[i - 1];
    for (int j = 1; j <= i - 1; ++j) {
      s21 += a[i - 1][j - 1] * b[j - 1];
    }
    b[i - 1] = -a[i - 1][i - 1] * s21;
  }

  for (int i = 1; i <= n - 1; ++i) {
    double s22 = -b[n - i - 1];
    for (int j = 1; j <= i; ++j) {
      s22 += a[n - i - 1][n - j] * b[n - j];
    }
    b[n - i - 1] = -s22;
  }
}

int dinv(const int n, std::vector<std::vector<double> >& a) {

  if (n < 1) return 1;
  if (n > 3) {
    // Factorize matrix and invert.
    double det = 0.;
    int ifail = 0;
    int jfail = 0;
    std::vector<int> ir(n, 0);
    dfact(n, a, ir, ifail, det, jfail);
    if (ifail != 0) return ifail;
    dfinv(n, a, ir);
  } else if (n == 3) {
    // Compute cofactors.
    const double c11 = a[1][1] * a[2][2] - a[1][2] * a[2][1];
    const double c12 = a[1][2] * a[2][0] - a[1][0] * a[2][2];
    const double c13 = a[1][0] * a[2][1] - a[1][1] * a[2][0];
    const double c21 = a[2][1] * a[0][2] - a[2][2] * a[0][1];
    const double c22 = a[2][2] * a[0][0] - a[2][0] * a[0][2];
    const double c23 = a[2][0] * a[0][1] - a[2][1] * a[0][0];
    const double c31 = a[0][1] * a[1][2] - a[0][2] * a[1][1];
    const double c32 = a[0][2] * a[1][0] - a[0][0] * a[1][2];
    const double c33 = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    const double t1 = fabs(a[0][0]);
    const double t2 = fabs(a[1][0]);
    const double t3 = fabs(a[2][0]);
    double det = 0.;
    double pivot = 0.;
    if (t2 < t1 && t3 < t1) {
      pivot = a[0][0];
      det = c22 * c33 - c23 * c32;
    } else if (t1 < t2 && t3 < t2) {
      pivot = a[1][0];
      det = c13 * c32 - c12 * c33;
    } else {
      pivot = a[2][0];
      det = c23 * c12 - c22 * c13;
    }
    // Set elements of inverse in A.
    if (det == 0.) return -1;
    double s = pivot / det;
    a[0][0] = s * c11;
    a[0][1] = s * c21;
    a[0][2] = s * c31;
    a[1][0] = s * c12;
    a[1][1] = s * c22;
    a[1][2] = s * c32;
    a[2][0] = s * c13;
    a[2][1] = s * c23;
    a[2][2] = s * c33;
  } else if (n == 2) {
    // Cramer's rule.
    const double det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    if (det == 0.) return -1;
    const double s = 1. / det;
    const double c11 = s * a[1][1];
    a[0][1] = -s * a[0][1];
    a[1][0] = -s * a[1][0];
    a[1][1] =  s * a[0][0];
    a[0][0] = c11;
  } else if (n == 1) {
    if (a[0][0] == 0.) return -1;
    a[0][0] = 1. / a[0][0];
  }
  return 0;
}

void dfinv(const int n, std::vector<std::vector<double> >& a,
           std::vector<int>& ir) {

  if (n <= 1) return;
  a[1][0] = -a[1][1] * a[0][0] * a[1][0];
  a[0][1] = -a[0][1];
  if (n > 2) {
    for (int i = 3; i <= n; ++i) {
      for (int j = 1; j <= i - 2; ++j) {
        double s31 = 0.;
        double s32 = a[j - 1][i - 1];
        for (int k = j; k <= i - 2; ++k) {
          s31 += a[k - 1][j - 1] * a[i - 1][k - 1];
          s32 += a[j - 1][k] * a[k][i - 1];
        }
        a[i - 1][j - 1] =
            -a[i - 1][i - 1] * (s31 + a[i - 2][j - 1] * a[i - 1][i - 2]);
        a[j - 1][i - 1] = -s32;
      }
      a[i - 1][i - 2] = -a[i - 1][i - 1] * a[i - 2][i - 2] * a[i - 1][i - 2];
      a[i - 2][i - 1] = -a[i - 2][i - 1];
    }
  }

  for (int i = 1; i <= n - 1; ++i) {
    for (int j = 1; j <= i; ++j) {
      double s33 = a[i - 1][j - 1];
      for (int k = 1; k <= n - i; ++k) {
        s33 += a[i + k - 1][j - 1] * a[i - 1][i + k - 1];
      }
      a[i - 1][j - 1] = s33;
    }
    for (int j = 1; j <= n - i; ++j) {
      double s34 = 0.;
      for (int k = j; k <= n - i; ++k) {
        s34 += a[i + k - 1][i + j - 1] * a[i - 1][i + k - 1];
      }
      a[i - 1][i + j - 1] = s34;
    }
  }

  int nxch = ir[n - 1];
  if (nxch == 0) return;

  for (int m = 1; m <= nxch; ++m) {
    int k = nxch - m + 1;
    int ij = ir[k - 1];
    int i = ij / 4096;
    int j = ij % 4096;
    for (k = 1; k <= n; ++k) {
      std::swap(a[k - 1][i - 1], a[k - 1][j - 1]);
    }
  }
}

int deqinv(const int n, std::vector<std::vector<double> >& a, 
           std::vector<double>& b) {

  // Test for parameter errors.
  if (n < 1) return 1;

  double det = 0.;
  if (n > 3) {
    // n > 3 cases. Factorize matrix, invert and solve system.
    std::vector<int> ir(n, 0);
    int ifail = 0, jfail = 0;
    dfact(n, a, ir, ifail, det, jfail);
    if (ifail != 0) return ifail;
    dfeqn(n, a, ir, b);
    dfinv(n, a, ir);
  } else if (n == 3) {
    // n = 3 case. Compute cofactors.
    const double c11 = a[1][1] * a[2][2] - a[1][2] * a[2][1];
    const double c12 = a[1][2] * a[2][0] - a[1][0] * a[2][2];
    const double c13 = a[1][0] * a[2][1] - a[1][1] * a[2][0];
    const double c21 = a[2][1] * a[0][2] - a[2][2] * a[0][1];
    const double c22 = a[2][2] * a[0][0] - a[2][0] * a[0][2];
    const double c23 = a[2][0] * a[0][1] - a[2][1] * a[0][0];
    const double c31 = a[0][1] * a[1][2] - a[0][2] * a[1][1];
    const double c32 = a[0][2] * a[1][0] - a[0][0] * a[1][2];
    const double c33 = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    const double t1 = fabs(a[0][0]);
    const double t2 = fabs(a[1][0]);
    const double t3 = fabs(a[2][0]);

    // Set temp = pivot and det = pivot * det.
    double temp = 0.;
    if (t2 < t1 && t3 < t1) {
      temp = a[0][0];
      det = c22 * c33 - c23 * c32;
    } else if (t1 < t2 && t3 < t2) {
      temp = a[1][0];
      det = c13 * c32 - c12 * c33;
    } else {
      temp = a[2][0];
      det = c23 * c12 - c22 * c13;
    }

    // Set elements of inverse in A.
    if (det == 0.) return -1;
    const double s = temp / det;
    a[0][0] = s * c11;
    a[0][1] = s * c21;
    a[0][2] = s * c31;
    a[1][0] = s * c12;
    a[1][1] = s * c22;
    a[1][2] = s * c32;
    a[2][0] = s * c13;
    a[2][1] = s * c23;
    a[2][2] = s * c33;

    // Replace b by Ainv * b.
    const double b1 = b[0];
    const double b2 = b[1];
    b[0] = a[0][0] * b1 + a[0][1] * b2 + a[0][2] * b[2];
    b[1] = a[1][0] * b1 + a[1][1] * b2 + a[1][2] * b[2];
    b[2] = a[2][0] * b1 + a[2][1] * b2 + a[2][2] * b[2];
  } else if (n == 2) {
    // n = 2 case by Cramer's rule.
    det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    if (det == 0.) return -1;
    const double s = 1. / det;
    const double c11 = s * a[1][1];
    a[0][1] = -s * a[0][1];
    a[1][0] = -s * a[1][0];
    a[1][1] = s * a[0][0];
    a[0][0] = c11;

    const double b1 = b[0];
    b[0] = c11 * b1 + a[0][1] * b[1];
    b[1] = a[1][0] * b1 + a[1][1] * b[1];
  } else {
    // n = 1 case.
    if (a[0][0] == 0.) return -1;
    a[0][0] = 1. / a[0][0];
    b[0] = a[0][0] * b[0];
  }
  return 0;
}

void cfact(const int n, std::vector<std::vector<std::complex<double> > >& a,
           std::vector<int>& ir, int& ifail, std::complex<double>& det,
           int& jfail) {
  constexpr double g1 = 1.e-19;
  constexpr double g2 = 1.e-19;

  ifail = jfail = 0;

  int nxch = 0;
  det = std::complex<double>(1., 0.);

  for (int j = 1; j <= n; ++j) {
    int k = j;
    double p = std::max(fabs(real(a[j - 1][j - 1])), fabs(imag(a[j - 1][j - 1])));
    if (j == n) {
      if (p <= 0.) {
        det = std::complex<double>(0., 0.);
        ifail = -1;
        jfail = 0;
        return;
      }
      det *= a[j - 1][j - 1];
      a[j - 1][j - 1] = std::complex<double>(1., 0.) / a[j - 1][j - 1];
      const double t = std::max(fabs(real(det)), fabs(imag(det)));
      if (t < g1) {
        det = std::complex<double>(0., 0.);
        if (jfail == 0) jfail = -1;
      } else if (t > g2) {
        det = std::complex<double>(1., 0.);
        if (jfail == 0) jfail = +1;
      }
      continue;
    }
    for (int i = j + 1; i <= n; ++i) {
      double q = std::max(fabs(real(a[i - 1][j - 1])), fabs(imag(a[i - 1][j - 1])));
      if (q <= p) continue;
      k = i;
      p = q;
    }
    if (k != j) {
      for (int l = 1; l <= n; ++l) {
        const auto tf = a[j - 1][l - 1];
        a[j - 1][l - 1] = a[k - 1][l - 1];
        a[k - 1][l - 1] = tf;
      }
      ++nxch;
      ir[nxch - 1] = j * 4096 + k;
    } else if (p <= 0.) {
      det = std::complex<double>(0., 0.);
      ifail = -1;
      jfail = 0;
      return;
    }
    det *= a[j - 1][j - 1];
    a[j - 1][j - 1] = 1. / a[j - 1][j - 1];
    const double t = std::max(fabs(real(det)), fabs(imag(det)));
    if (t < g1) {
      det = std::complex<double>(0., 0.);
      if (jfail == 0) jfail = -1;
    } else if (t > g2) {
      det = std::complex<double>(1., 0.);
      if (jfail == 0) jfail = +1;
    }
    for (k = j + 1; k <= n; ++k) {
      auto s11 = -a[j - 1][k - 1];
      auto s12 = -a[k - 1][j];
      if (j == 1) {
        a[j - 1][k - 1] = -s11 * a[j - 1][j - 1];
        a[k - 1][j] = -(s12 + a[j - 1][j] * a[k - 1][j - 1]);
        continue;
      }
      for (int i = 1; i <= j - 1; ++i) {
        s11 += a[i - 1][k - 1] * a[j - 1][i - 1];
        s12 += a[i - 1][j] * a[k - 1][i - 1];
      }
      a[j - 1][k - 1] = -s11 * a[j - 1][j - 1];
      a[k - 1][j] = -a[j - 1][j] * a[k - 1][j - 1] - s12;
    }
  }

  if (nxch % 2 != 0) det = -det;
  if (jfail != 0) det = std::complex<double>(0., 0.);
  ir[n - 1] = nxch;
}

void cfinv(const int n, std::vector<std::vector<std::complex<double> > >& a,
           std::vector<int>& ir) {

  if (n <= 1) return;
  a[1][0] = -a[1][1] * a[0][0] * a[1][0];
  a[0][1] = -a[0][1];
  if (n > 2) {
    for (int i = 3; i <= n; ++i) {
      for (int j = 1; j <= i - 2; ++j) {
        auto s31 = std::complex<double>(0., 0.);
        auto s32 = a[j - 1][i - 1];
        for (int k = j; k <= i - 2; ++k) {
          s31 += a[k - 1][j - 1] * a[i - 1][k - 1];
          s32 += a[j - 1][k] * a[k][i - 1];
        }
        a[i - 1][j - 1] =
            -a[i - 1][i - 1] * (s31 + a[i - 2][j - 1] * a[i - 1][i - 2]);
        a[j - 1][i - 1] = -s32;
      }
      a[i - 1][i - 2] = -a[i - 1][i - 1] * a[i - 2][i - 2] * a[i - 1][i - 2];
      a[i - 2][i - 1] = -a[i - 2][i - 1];
    }
  }

  for (int i = 1; i <= n - 1; ++i) {
    for (int j = 1; j <= i; ++j) {
      auto s33 = a[i - 1][j - 1];
      for (int k = 1; k <= n - i; ++k) {
        s33 += a[i + k - 1][j - 1] * a[i - 1][i + k - 1];
      }
      a[i - 1][j - 1] = s33;
    }
    for (int j = 1; j <= n - i; ++j) {
      std::complex<double> s34(0., 0.);
      for (int k = j; k <= n - i; ++k) {
        s34 += a[i + k - 1][i + j - 1] * a[i - 1][i + k - 1];
      }
      a[i - 1][i + j - 1] = s34;
    }
  }

  int nxch = ir[n - 1];
  if (nxch == 0) return;

  for (int m = 1; m <= nxch; ++m) {
    int k = nxch - m + 1;
    const int ij = ir[k - 1];
    const int i = ij / 4096;
    const int j = ij % 4096;
    for (k = 1; k <= n; ++k) {
      std::swap(a[k - 1][i - 1], a[k - 1][j - 1]);
    }
  }
}

int cinv(const int n, std::vector<std::vector<std::complex<double> > >& a) {

  // Test for parameter errors.
  if (n < 1) return 1;

  std::complex<double> det(0., 0.);
  if (n > 3) {
    // n > 3 cases. Factorize matrix and invert.
    std::vector<int> ir(n, 0);
    int ifail = 0, jfail = 0;
    cfact(n, a, ir, ifail, det, jfail);
    if (ifail != 0) return ifail;
    cfinv(n, a, ir);
  } else if (n == 3) {
    // n = 3 case. Compute cofactors.
    const auto c11 = a[1][1] * a[2][2] - a[1][2] * a[2][1];
    const auto c12 = a[1][2] * a[2][0] - a[1][0] * a[2][2];
    const auto c13 = a[1][0] * a[2][1] - a[1][1] * a[2][0];
    const auto c21 = a[2][1] * a[0][2] - a[2][2] * a[0][1];
    const auto c22 = a[2][2] * a[0][0] - a[2][0] * a[0][2];
    const auto c23 = a[2][0] * a[0][1] - a[2][1] * a[0][0];
    const auto c31 = a[0][1] * a[1][2] - a[0][2] * a[1][1];
    const auto c32 = a[0][2] * a[1][0] - a[0][0] * a[1][2];
    const auto c33 = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    const double t1 = fabs(real(a[0][0])) + fabs(imag(a[0][0]));
    const double t2 = fabs(real(a[1][0])) + fabs(imag(a[1][0]));
    const double t3 = fabs(real(a[2][0])) + fabs(imag(a[2][0]));

    // Set temp = pivot and det = pivot * det.
    std::complex<double> temp(0., 0.);
    if (t2 < t1 && t3 < t1) {
      temp = a[0][0];
      det = c22 * c33 - c23 * c32;
    } else if (t1 < t2 && t3 < t2) {
      temp = a[1][0];
      det = c13 * c32 - c12 * c33;
    } else {
      temp = a[2][0];
      det = c23 * c12 - c22 * c13;
    }
    // Set elements of inverse in A.
    if (real(det) == 0. && imag(det) == 0.) return -1;
    const auto s = temp / det;
    a[0][0] = s * c11;
    a[0][1] = s * c21;
    a[0][2] = s * c31;
    a[1][0] = s * c12;
    a[1][1] = s * c22;
    a[1][2] = s * c32;
    a[2][0] = s * c13;
    a[2][1] = s * c23;
    a[2][2] = s * c33;
  } else if (n == 2) {
    // n=2 case by Cramer's rule.
    det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    if (real(det) == 0. && imag(det) == 0.) return -1;
    const auto s = std::complex<double>(1., 0.) / det;
    const auto c11 = s * a[1][1];
    a[0][1] = -s * a[0][1];
    a[1][0] = -s * a[1][0];
    a[1][1] = s * a[0][0];
    a[0][0] = c11;
  } else {
    // n = 1 case.
    if (real(a[0][0]) == 0. && imag(a[0][0]) == 0.) return -1;
    a[0][0] = std::complex<double>(1., 0.) / a[0][0];
  }
  return 0;
}

void cfft(std::vector<std::complex<double> >& a, const int msign) {

  if (msign == 0) return;
  const int m = std::abs(msign);
  const int n = pow(2, m);
  // Bit reversal.
  int j = 1;
  for (int i = 1; i <= n - 1; ++i) {
    if (i < j) std::swap(a[j - 1], a[i - 1]);
    int k = n / 2;
    while (k < j) {
      j -= k;
      k /= 2;
    }
    j += k;
  }
  for (int i = 1; i <= n; i += 2) {
    const auto t = a[i];
    a[i] = a[i - 1] - t;
    a[i - 1] += t;
  }
  if (m == 1) return;
  double c = 0.;
  double s = msign >= 0 ? 1. : -1.;
  int le = 2;
  for (int l = 2; l <= m; ++l) {
    std::complex<double> w(c, s);
    std::complex<double> u = w;
    c = sqrt(0.5 * c + 0.5);
    s = imag(w) / (c + c);
    int le1 = le;
    le = le1 + le1;
    for (int i = 1; i <= n; i += le) {
      const auto t = a[i + le1 - 1];
      a[i + le1 - 1] = a[i - 1] - t;
      a[i - 1] += t;
    }
    for (j = 2; j <= le1; ++j) {
      for (int i = j; i <= n; i += le) {
        const auto t = a[i + le1 - 1] * u;
        a[i + le1 - 1] = a[i - 1] - t;
        a[i - 1] += t;
      }
      u *= w;
    }
  }
}

}

double Divdif(const std::vector<double>& f, const std::vector<double>& a,
              const int nn, const double x, const int mm) {

  double t[20], d[20];

  // Check the arguments.
  if (nn < 2) {
    std::cerr << "Divdif: Array length < 2.\n";
    return 0.;
  }
  if (mm < 1) {
    std::cerr << "Divdif: Interpolation order < 1.\n";
    return 0.;
  }

  // Deal with the case that X is located at first or last point.
  const double tol1 = 1.e-6 * fabs(fabs(a[1]) - fabs(a[0]));
  const double tol2 = 1.e-6 * fabs(fabs(a[nn-1]) - fabs(a[nn-2]));
  if (fabs(x - a[0]) < tol1) return f[0];
  if (fabs(x - a[nn - 1]) < tol2) return f[nn - 1];

  // Find subscript IX of X in array A.
  constexpr int mmax = 10;
  const int m = std::min({mm, mmax, nn - 1});
  const int mplus = m + 1;
  int ix = 0;
  int iy = nn + 1;
  if (a[0] > a[nn - 1]) {
    // Search decreasing arguments.
    do {
      const int mid = (ix + iy) / 2;
      if (x > a[mid - 1]) {
        iy = mid;
      } else {
        ix = mid;
      }
    } while (iy - ix > 1);
  } else {
    // Search increasing arguments.
    do {
      const int mid = (ix + iy) / 2;
      if (x < a[mid - 1]) {
        iy = mid;
      } else {
        ix = mid;
      }
    } while (iy - ix > 1);
  }
  //  Copy reordered interpolation points into (T[I],D[I]), setting
  //  EXTRA to True if M+2 points to be used.
  int npts = m + 2 - (m % 2);
  int ip = 0;
  int l = 0;
  do {
    const int isub = ix + l;
    if ((1 > isub) || (isub > nn)) {
      // Skip point.
      npts = mplus;
    } else {
      // Insert point.
      ip++;
      t[ip - 1] = a[isub - 1];
      d[ip - 1] = f[isub - 1];
    }
    if (ip < npts) {
      l = -l;
      if (l >= 0) ++l;
    }
  } while (ip < npts);

  const bool extra = npts != mplus;
  // Replace d by the leading diagonal of a divided-difference table,
  // supplemented by an extra line if EXTRA is True.
  for (l = 1; l <= m; l++) {
    if (extra) {
      const int isub = mplus - l;
      d[m + 1] = (d[m + 1] - d[m - 1]) / (t[m + 1] - t[isub - 1]);
    }
    int i = mplus;
    for (int j = l; j <= m; j++) {
      const int isub = i - l;
      d[i - 1] = (d[i - 1] - d[i - 2]) / (t[i - 1] - t[isub - 1]);
      i--;
    }
  }
  // Evaluate the Newton interpolation formula at X, averaging two values
  // of last difference if EXTRA is True.
  double sum = d[mplus - 1];
  if (extra) {
    sum = 0.5 * (sum + d[m + 1]);
  }
  int j = m;
  for (l = 1; l <= m; l++) {
    sum = d[j - 1] + (x - t[j - 1]) * sum;
    j--;
  }
  return sum;
}

bool Boxin2(const std::vector<std::vector<double> >& value,
            const std::vector<double>& xAxis, const std::vector<double>& yAxis,
            const int nx, const int ny, const double x, const double y,
            double& f, const int iOrder) {
  //-----------------------------------------------------------------------
  //   BOXIN2 - Interpolation of order 1 and 2 in an irregular rectangular
  //            2-dimensional grid.
  //-----------------------------------------------------------------------
  int iX0 = 0, iX1 = 0;
  int iY0 = 0, iY1 = 0;
  std::array<double, 3> fX;
  std::array<double, 3> fY;
  f = 0.;
  // Ensure we are in the grid.
  if ((xAxis[nx - 1] - x) * (x - xAxis[0]) < 0 ||
      (yAxis[ny - 1] - y) * (y - yAxis[0]) < 0) {
    std::cerr << "Boxin2: Point not in the grid; no interpolation.\n";
    return false;
  }
  // Make sure we have enough points.
  if (iOrder < 0 || iOrder > 2) {
    std::cerr << "Boxin2: Incorrect order; no interpolation.\n";
    return false;
  } else if (nx < 1 || ny < 1) {
    std::cerr << "Boxin2: Incorrect number of points; no interpolation.\n";
    return false;
  }
  if (iOrder == 0 || nx <= 1) {
    // Zeroth order interpolation in x.
    // Find the nearest node.
    double dist = fabs(x - xAxis[0]);
    int iNode = 0;
    for (int i = 1; i < nx; i++) {
      if (fabs(x - xAxis[i]) < dist) {
        dist = fabs(x - xAxis[i]);
        iNode = i;
      }
    }
    // Set the summing range.
    iX0 = iNode;
    iX1 = iNode;
    // Establish the shape functions.
    fX = {1., 0., 0.};
  } else if (iOrder == 1 || nx <= 2) {
    // First order interpolation in x.
    // Find the grid segment containing this point.
    int iGrid = 0;
    for (int i = 1; i < nx; i++) {
      if ((xAxis[i - 1] - x) * (x - xAxis[i]) >= 0.) {
        iGrid = i;
      }
    }
    const double x0 = xAxis[iGrid - 1];
    const double x1 = xAxis[iGrid];
    // Ensure there won't be divisions by zero.
    if (x1 == x0) {
      std::cerr << "Boxin2: Incorrect grid; no interpolation.\n";
      return false;
    }
    // Compute local coordinates.
    const double xL = (x - x0) / (x1 - x0);
    // Set the summing range.
    iX0 = iGrid - 1;
    iX1 = iGrid;
    // Set the shape functions.
    fX = {1. - xL, xL, 0.};
  } else if (iOrder == 2) {
    // Second order interpolation in x.
    // Find the nearest node and the grid segment.
    double dist = fabs(x - xAxis[0]);
    int iNode = 0;
    for (int i = 1; i < nx; ++i) {
      if (fabs(x - xAxis[i]) < dist) {
        dist = fabs(x - xAxis[i]);
        iNode = i;
      }
    }
    // Find the nearest fitting 2x2 matrix.
    int iGrid = std::max(1, std::min(nx - 2, iNode));
    const double x0 = xAxis[iGrid - 1];
    const double x1 = xAxis[iGrid];
    const double x2 = xAxis[iGrid + 1];
    // Ensure there won't be divisions by zero.
    if (x2 == x0) {
      std::cerr << "Boxin2: Incorrect grid; no interpolation.\n";
      return false;
    }
    // Compute the alpha and local coordinate for this grid segment.
    const double xAlpha = (x1 - x0) / (x2 - x0);
    const double xL = (x - x0) / (x2 - x0);
    // Ensure there won't be divisions by zero.
    if (xAlpha <= 0 || xAlpha >= 1) {
      std::cerr << "Boxin2: Incorrect grid; no interpolation.\n";
      return false;
    }
    // Set the summing range.
    iX0 = iGrid - 1;
    iX1 = iGrid + 1;
    // Set the shape functions.
    const double xL2 = xL * xL;
    fX[0] = xL2 / xAlpha - xL * (1. + xAlpha) / xAlpha + 1.;
    fX[1] = (xL2 - xL) / (xAlpha * xAlpha - xAlpha);
    fX[2] = (xL2 - xL * xAlpha) / (1. - xAlpha);
  }
  if (iOrder == 0 || ny <= 1) {
    // Zeroth order interpolation in y.
    // Find the nearest node.
    double dist = fabs(y - yAxis[0]);
    int iNode = 0;
    for (int i = 1; i < ny; i++) {
      if (fabs(y - yAxis[i]) < dist) {
        dist = fabs(y - yAxis[i]);
        iNode = i;
      }
    }
    // Set the summing range.
    iY0 = iNode;
    iY1 = iNode;
    // Establish the shape functions.
    fY = {1., 0., 0.};
  } else if (iOrder == 1 || ny <= 2) {
    // First order interpolation in y.
    // Find the grid segment containing this point.
    int iGrid = 0;
    for (int i = 1; i < ny; ++i) {
      if ((yAxis[i - 1] - y) * (y - yAxis[i]) >= 0) {
        iGrid = i;
      }
    }
    const double y0 = yAxis[iGrid - 1];
    const double y1 = yAxis[iGrid];
    // Ensure there won't be divisions by zero.
    if (y1 == y0) {
      std::cerr << "Boxin2: Incorrect grid; no interpolation.\n";
      return false;
    }
    // Compute local coordinates.
    const double yL = (y - y0) / (y1 - y0);
    // Set the summing range.
    iY0 = iGrid - 1;
    iY1 = iGrid;
    // Set the shape functions.
    fY = {1. - yL, yL, 0.};
  } else if (iOrder == 2) {
    // Second order interpolation in y.
    // Find the nearest node and the grid segment.
    double dist = fabs(y - yAxis[0]);
    int iNode = 0;
    for (int i = 1; i < ny; ++i) {
      if (fabs(y - yAxis[i]) < dist) {
        dist = fabs(y - yAxis[i]);
        iNode = i;
      }
    }
    // Find the nearest fitting 2x2 matrix.
    int iGrid = std::max(1, std::min(ny - 2, iNode));
    const double y0 = yAxis[iGrid - 1];
    const double y1 = yAxis[iGrid];
    const double y2 = yAxis[iGrid + 1];
    // Ensure there won't be divisions by zero.
    if (y2 == y0) {
      std::cerr << "Boxin2: Incorrect grid; no interpolation.\n";
      return false;
    }
    // Compute the alpha and local coordinate for this grid segment.
    const double yAlpha = (y1 - y0) / (y2 - y0);
    const double yL = (y - y0) / (y2 - y0);
    // Ensure there won't be divisions by zero.
    if (yAlpha <= 0 || yAlpha >= 1) {
      std::cerr << "Boxin2: Incorrect grid; no interpolation.\n";
      return false;
    }
    // Set the summing range.
    iY0 = iGrid - 1;
    iY1 = iGrid + 1;
    // Set the shape functions.
    const double yL2 = yL * yL;
    fY[0] = yL2 / yAlpha - yL * (1. + yAlpha) / yAlpha + 1.;
    fY[1] = (yL2 - yL) / (yAlpha * yAlpha - yAlpha);
    fY[2] = (yL2 - yL * yAlpha) / (1. - yAlpha);
  }

  // Sum the shape functions.
  for (int i = iX0; i <= iX1; ++i) {
    for (int j = iY0; j <= iY1; ++j) {
      f += value[i][j] * fX[i - iX0] * fY[j - iY0];
    }
  }
  return true;
}

bool Boxin3(const std::vector<std::vector<std::vector<double> > >& value,
            const std::vector<double>& xAxis, const std::vector<double>& yAxis,
            const std::vector<double>& zAxis, const int nx, const int ny,
            const int nz, const double xx, const double yy, const double zz,
            double& f, const int iOrder) {
  //-----------------------------------------------------------------------
  //   BOXIN3 - interpolation of order 1 and 2 in an irregular rectangular
  //            3-dimensional grid.
  //-----------------------------------------------------------------------
  int iX0 = 0, iX1 = 0;
  int iY0 = 0, iY1 = 0;
  int iZ0 = 0, iZ1 = 0;
  std::array<double, 4> fX;
  std::array<double, 4> fY;
  std::array<double, 4> fZ;

  f = 0.;
  // Ensure we are in the grid.
  const double x = std::min(std::max(xx, std::min(xAxis[0], xAxis[nx - 1])),
                            std::max(xAxis[0], xAxis[nx - 1]));
  const double y = std::min(std::max(yy, std::min(yAxis[0], yAxis[ny - 1])),
                            std::max(yAxis[0], yAxis[ny - 1]));
  const double z = std::min(std::max(zz, std::min(zAxis[0], zAxis[nz - 1])),
                            std::max(zAxis[0], zAxis[nz - 1]));

  // Make sure we have enough points.
  if (iOrder < 0 || iOrder > 2) {
    std::cerr << "Boxin3: Incorrect order; no interpolation.\n";
    return false;
  } else if (nx < 1 || ny < 1 || nz < 1) {
    std::cerr << "Boxin3: Incorrect number of points; no interpolation.\n";
    return false;
  }
  if (iOrder == 0 || nx == 1) {
    // Zeroth order interpolation in x.
    // Find the nearest node.
    double dist = fabs(x - xAxis[0]);
    int iNode = 0;
    for (int i = 1; i < nx; i++) {
      if (fabs(x - xAxis[i]) < dist) {
        dist = fabs(x - xAxis[i]);
        iNode = i;
      }
    }
    // Set the summing range.
    iX0 = iNode;
    iX1 = iNode;
    // Establish the shape functions.
    fX = {1., 0., 0., 0.};
  } else if (iOrder == 1 || nx == 2) {
    // First order interpolation in x.
    // Find the grid segment containing this point.
    int iGrid = 0;
    for (int i = 1; i < nx; i++) {
      if ((xAxis[i - 1] - x) * (x - xAxis[i]) >= 0.) {
        iGrid = i;
      }
    }
    const double x0 = xAxis[iGrid - 1];
    const double x1 = xAxis[iGrid];
    // Ensure there won't be divisions by zero.
    if (x1 == x0) {
      std::cerr << "Boxin3: Incorrect grid; no interpolation.\n";
      return false;
    }
    // Compute local coordinates.
    const double xL = (x - x0) / (x1 - x0);
    // Set the summing range.
    iX0 = iGrid - 1;
    iX1 = iGrid;
    // Set the shape functions.
    fX = {1. - xL, xL, 0., 0.};
  } else if (iOrder == 2) {
    // Second order interpolation in x.
    // Find the grid segment containing this point.
    int iGrid = 0;
    for (int i = 1; i < nx; i++) {
      if ((xAxis[i - 1] - x) * (x - xAxis[i]) >= 0.) {
        iGrid = i;
      }
    }
    if (iGrid == 1) {
      iX0 = iGrid - 1;
      iX1 = iGrid + 1;
      const double x0 = xAxis[iX0];
      const double x1 = xAxis[iX0 + 1];
      const double x2 = xAxis[iX0 + 2];
      if (x0 == x1 || x0 == x2 || x1 == x2) {
        std::cerr << "Boxin3: One or more grid points in x coincide.\n"
                  << "    No interpolation.\n";
        return false;
      }
      fX[0] = (x - x1) * (x - x2) / ((x0 - x1) * (x0 - x2));
      fX[1] = (x - x0) * (x - x2) / ((x1 - x0) * (x1 - x2));
      fX[2] = (x - x0) * (x - x1) / ((x2 - x0) * (x2 - x1));
    } else if (iGrid == nx - 1) {
      iX0 = iGrid - 2;
      iX1 = iGrid;
      const double x0 = xAxis[iX0];
      const double x1 = xAxis[iX0 + 1];
      const double x2 = xAxis[iX0 + 2];
      if (x0 == x1 || x0 == x2 || x1 == x2) {
        std::cerr << "Boxin3: One or more grid points in x coincide.\n"
                  << "    No interpolation.\n";
        return false;
      }
      fX[0] = (x - x1) * (x - x2) / ((x0 - x1) * (x0 - x2));
      fX[1] = (x - x0) * (x - x2) / ((x1 - x0) * (x1 - x2));
      fX[2] = (x - x0) * (x - x1) / ((x2 - x0) * (x2 - x1));
    } else {
      iX0 = iGrid - 2;
      iX1 = iGrid + 1;
      const double x0 = xAxis[iX0];
      const double x1 = xAxis[iX0 + 1];
      const double x2 = xAxis[iX0 + 2];
      const double x3 = xAxis[iX0 + 3];
      if (x0 == x1 || x0 == x2 || x0 == x3 || 
          x1 == x2 || x1 == x3 || x2 == x3) {
        std::cerr << "Boxin3: One or more grid points in x coincide.\n"
                  << "    No interpolation.\n";
        return false;
      }
      // Compute the local coordinate for this grid segment.
      const double xL = (x - x1) / (x2 - x1);
      fX[0] = ((x - x1) * (x - x2) / ((x0 - x1) * (x0 - x2))) * (1. - xL);
      fX[1] = ((x - x0) * (x - x2) / ((x1 - x0) * (x1 - x2))) * (1. - xL) +
              ((x - x2) * (x - x3) / ((x1 - x2) * (x1 - x3))) * xL;
      fX[2] = ((x - x0) * (x - x1) / ((x2 - x0) * (x2 - x1))) * (1. - xL) + 
              ((x - x1) * (x - x3) / ((x2 - x1) * (x2 - x3))) * xL;
      fX[3] = ((x - x1) * (x - x2) / ((x3 - x1) * (x3 - x2))) * xL;
    }
  }

  if (iOrder == 0 || ny == 1) {
    // Zeroth order interpolation in y.
    // Find the nearest node.
    double dist = fabs(y - yAxis[0]);
    int iNode = 0;
    for (int i = 1; i < ny; i++) {
      if (fabs(y - yAxis[i]) < dist) {
        dist = fabs(y - yAxis[i]);
        iNode = i;
      }
    }
    // Set the summing range.
    iY0 = iNode;
    iY1 = iNode;
    // Establish the shape functions.
    fY = {1., 0., 0., 0.};
  } else if (iOrder == 1 || ny == 2) {
    // First order interpolation in y.
    // Find the grid segment containing this point.
    int iGrid = 0;
    for (int i = 1; i < ny; i++) {
      if ((yAxis[i - 1] - y) * (y - yAxis[i]) >= 0.) {
        iGrid = i;
      }
    }
    // Ensure there won't be divisions by zero.
    const double y0 = yAxis[iGrid - 1];
    const double y1 = yAxis[iGrid];
    if (y1 == y0) {
      std::cerr << "Boxin3: Incorrect grid; no interpolation.\n";
      return false;
    }
    // Compute local coordinates.
    const double yL = (y - y0) / (y1 - y0);
    // Set the summing range.
    iY0 = iGrid - 1;
    iY1 = iGrid;
    // Set the shape functions.
    fY = {1. - yL, yL, 0., 0.};
  } else if (iOrder == 2) {
    // Second order interpolation in y.
    // Find the grid segment containing this point.
    int iGrid = 0;
    for (int i = 1; i < ny; i++) {
      if ((yAxis[i - 1] - y) * (y - yAxis[i]) >= 0.) {
        iGrid = i;
      }
    }
    if (iGrid == 1) {
      iY0 = iGrid - 1;
      iY1 = iGrid + 1;
      const double y0 = yAxis[iY0];
      const double y1 = yAxis[iY0 + 1];
      const double y2 = yAxis[iY0 + 2];
      if (y0 == y1 || y0 == y2 || y1 == y2) {
        std::cerr << "Boxin3: One or more grid points in y coincide.\n"
                  << "    No interpolation.\n";
        return false;
      }
      fY[0] = (y - y1) * (y - y2) / ((y0 - y1) * (y0 - y2));
      fY[1] = (y - y0) * (y - y2) / ((y1 - y0) * (y1 - y2));
      fY[2] = (y - y0) * (y - y1) / ((y2 - y0) * (y2 - y1));
    } else if (iGrid == ny - 1) {
      iY0 = iGrid - 2;
      iY1 = iGrid;
      const double y0 = yAxis[iY0];
      const double y1 = yAxis[iY0 + 1];
      const double y2 = yAxis[iY0 + 2];
      if (y0 == y1 || y0 == y2 || y1 == y2) {
        std::cerr << "Boxin3: One or more grid points in y coincide.\n"
                  << "    No interpolation.\n";
        return false;
      }
      fY[0] = (y - y1) * (y - y2) / ((y0 - y1) * (y0 - y2));
      fY[1] = (y - y0) * (y - y2) / ((y1 - y0) * (y1 - y2));
      fY[2] = (y - y0) * (y - y1) / ((y2 - y0) * (y2 - y1));
    } else {
      iY0 = iGrid - 2;
      iY1 = iGrid + 1;
      const double y0 = yAxis[iY0];
      const double y1 = yAxis[iY0 + 1];
      const double y2 = yAxis[iY0 + 2];
      const double y3 = yAxis[iY0 + 3];
      if (y0 == y1 || y0 == y2 || y0 == y3 || 
          y1 == y2 || y1 == y3 || y2 == y3) {
        std::cerr << "Boxin3: One or more grid points in y coincide.\n"
                  << "    No interpolation.\n";
        return false;
      }
      // Compute the local coordinate for this grid segment.
      const double yL = (y - y1) / (y2 - y1);
      fY[0] = ((y - y1) * (y - y2) / ((y0 - y1) * (y0 - y2))) * (1. - yL);
      fY[1] = ((y - y0) * (y - y2) / ((y1 - y0) * (y1 - y2))) * (1. - yL) +
              ((y - y2) * (y - y3) / ((y1 - y2) * (y1 - y3))) * yL;
      fY[2] = ((y - y0) * (y - y1) / ((y2 - y0) * (y2 - y1))) * (1. - yL) +
              ((y - y1) * (y - y3) / ((y2 - y1) * (y2 - y3))) * yL;
      fY[3] = ((y - y1) * (y - y2) / ((y3 - y1) * (y3 - y2))) * yL;
    }
  }

  if (iOrder == 0 || nz == 1) {
    // Zeroth order interpolation in z.
    // Find the nearest node.
    double dist = fabs(z - zAxis[0]);
    int iNode = 0;
    for (int i = 1; i < nz; i++) {
      if (fabs(z - zAxis[i]) < dist) {
        dist = fabs(z - zAxis[i]);
        iNode = i;
      }
    }
    // Set the summing range.
    iZ0 = iNode;
    iZ1 = iNode;
    // Establish the shape functions.
    fZ = {1., 0., 0., 0.};
  } else if (iOrder == 1 || nz == 2) {
    // First order interpolation in z.
    // Find the grid segment containing this point.
    int iGrid = 0;
    for (int i = 1; i < nz; i++) {
      if ((zAxis[i - 1] - z) * (z - zAxis[i]) >= 0.) {
        iGrid = i;
      }
    }
    const double z0 = zAxis[iGrid - 1];
    const double z1 = zAxis[iGrid];
    // Ensure there won't be divisions by zero.
    if (z1 == z0) {
      std::cerr << "Boxin3: Incorrect grid; no interpolation.\n";
      return false;
    }
    // Compute local coordinates.
    const double zL = (z - z0) / (z1 - z0);
    // Set the summing range.
    iZ0 = iGrid - 1;
    iZ1 = iGrid;
    // Set the shape functions.
    fZ = {1. - zL, zL, 0., 0.};
  } else if (iOrder == 2) {
    // Second order interpolation in z.
    // Find the grid segment containing this point.
    int iGrid = 0;
    for (int i = 1; i < nz; i++) {
      if ((zAxis[i - 1] - z) * (z - zAxis[i]) >= 0.) {
        iGrid = i;
      }
    }
    if (iGrid == 1) {
      iZ0 = iGrid - 1;
      iZ1 = iGrid + 1;
      const double z0 = zAxis[iZ0];
      const double z1 = zAxis[iZ0 + 1];
      const double z2 = zAxis[iZ0 + 2];
      if (z0 == z1 || z0 == z2 || z1 == z2) {
        std::cerr << "Boxin3: One or more grid points in z coincide.\n"
                  << "    No interpolation.\n";
        return false;
      }
      fZ[0] = (z - z1) * (z - z2) / ((z0 - z1) * (z0 - z2));
      fZ[1] = (z - z0) * (z - z2) / ((z1 - z0) * (z1 - z2));
      fZ[2] = (z - z0) * (z - z1) / ((z2 - z0) * (z2 - z1));
    } else if (iGrid == nz - 1) {
      iZ0 = iGrid - 2;
      iZ1 = iGrid;
      const double z0 = zAxis[iZ0];
      const double z1 = zAxis[iZ0 + 1];
      const double z2 = zAxis[iZ0 + 2];
      if (z0 == z1 || z0 == z2 || z1 == z2) {
        std::cerr << "Boxin3: One or more grid points in z coincide.\n"
                  << "    No interpolation.\n";
        return false;
      }
      fZ[0] = (z - z1) * (z - z2) / ((z0 - z1) * (z0 - z2));
      fZ[1] = (z - z0) * (z - z2) / ((z1 - z0) * (z1 - z2));
      fZ[2] = (z - z0) * (z - z1) / ((z2 - z0) * (z2 - z1));
    } else {
      iZ0 = iGrid - 2;
      iZ1 = iGrid + 1;
      const double z0 = zAxis[iZ0];
      const double z1 = zAxis[iZ0 + 1];
      const double z2 = zAxis[iZ0 + 2];
      const double z3 = zAxis[iZ0 + 3];
      if (z0 == z1 || z0 == z2 || z0 == z3 || 
          z1 == z2 || z1 == z3 || z2 == z3) {
        std::cerr << "Boxin3: One or more grid points in z coincide.\n"
                  << "    No interpolation.\n";
        return false;
      }
      // Compute the local coordinate for this grid segment.
      const double zL = (z - z1) / (z2 - z1);
      fZ[0] = ((z - z1) * (z - z2) / ((z0 - z1) * (z0 - z2))) * (1. - zL);
      fZ[1] = ((z - z0) * (z - z2) / ((z1 - z0) * (z1 - z2))) * (1. - zL) +
              ((z - z2) * (z - z3) / ((z1 - z2) * (z1 - z3))) * zL;
      fZ[2] = ((z - z0) * (z - z1) / ((z2 - z0) * (z2 - z1))) * (1. - zL) +
              ((z - z1) * (z - z3) / ((z2 - z1) * (z2 - z3))) * zL;
      fZ[3] = ((z - z1) * (z - z2) / ((z3 - z1) * (z3 - z2))) * zL;
    }
  }

  for (int i = iX0; i <= iX1; ++i) {
    for (int j = iY0; j <= iY1; ++j) {
      for (int k = iZ0; k <= iZ1; ++k) {
        f += value[i][j][k] * fX[i - iX0] * fY[j - iY0] * fZ[k - iZ0];
      }
    }
  }
  return true;
}

bool LeastSquaresFit(
    std::function<double(double, const std::vector<double>&)> f, 
    std::vector<double>& par, std::vector<double>& epar,
    const std::vector<double>& x, const std::vector<double>& y,
    const std::vector<double>& ey, const unsigned int nMaxIter,
    const double diff, double& chi2, const double eps, 
    const bool debug, const bool verbose) {

 //-----------------------------------------------------------------------
 //   LSQFIT - Subroutine fitting the parameters A in the routine F to
 //            the data points (X,Y) using a least squares method.
 //            Translated from an Algol routine written by Geert Jan van
 //            Oldenborgh and Rob Veenhof, based on Stoer + Bulirsch.
 //   VARIABLES : F( . ,A,VAL) : Subroutine to be fitted.
 //               (X,Y)        : Input data.
 //               D            : Derivative matrix.
 //               R            : Difference vector between Y and F(X,A).
 //               S            : Correction vector for A.
 //               EPSDIF       : Used for differentiating.
 //               EPS          : Numerical resolution.
 //   (Last updated on 23/ 5/11.)
 //-----------------------------------------------------------------------

  const unsigned int n = par.size();
  const unsigned int m = x.size();
  // Make sure that the # degrees of freedom < the number of data points.
  if (n > m) {
    std::cerr << "LeastSquaresFit: Number of parameters to be varied\n"
              << "                 exceeds the number of data points.\n";
    return false;
  }

  // Check the errors.
  if (*std::min_element(ey.cbegin(), ey.cend()) <= 0.) {
    std::cerr << "LeastSquaresFit: Not all errors are > 0; no fit done.\n";
    return false;
  }
  chi2 = 0.;
  // Largest difference.
  double diffc = -1.;
  // Initialise the difference vector R.
  std::vector<double> r(m, 0.);
  for (unsigned int i = 0; i < m; ++i) {
    // Compute initial residuals.
    r[i] = (y[i] - f(x[i], par)) / ey[i];
    // Compute initial maximum difference.
    diffc = std::max(std::abs(r[i]), diffc);
    // And compute initial chi2.
    chi2 += r[i] * r[i];
  }
  if (debug) {
    std::cout << "  Input data points:\n"
              << "                 X              Y          Y - F(X)\n";
    for (unsigned int i = 0; i < m; ++i) {
      std::printf(" %9u %15.8e %15.8e %15.8e\n", i, x[i], y[i], r[i]);
    }
    std::cout << "  Initial values of the fit parameters:\n"
              << "    Parameter            Value\n";
    for (unsigned int i = 0; i < n; ++i) {
      std::printf("    %9u  %15.8e\n", i, par[i]);
    }
  }
  if (verbose) {
    std::cout << "  MINIMISATION SUMMARY\n"
              << "  Initial situation:\n";
    std::printf("    Largest difference between fit and target: %15.8e\n",
                diffc);
    std::printf("    Sum of squares of these differences:       %15.8e\n",
                chi2);
  }
  // Start optimising loop.
  bool converged = false;
  double chi2L = 0.;
  for (unsigned int iter = 1; iter <= nMaxIter; ++iter) {
    // Check the stopping criteria: (1) max norm, (2) change in chi-squared.
    if ((diffc < diff) || (iter > 1 && std::abs(chi2L - chi2) < eps * chi2)) {
      if (debug || verbose) {
        if (diffc < diff) {
          std::cout << "  The maximum difference stopping criterion "
                    << "is satisfied.\n";
        } else {
          std::cout << "  The relative change in chi-squared has dropped "
                    << "below the threshold.\n";
        }
      }
      converged = true;
      break;
    } 
    // Calculate the derivative matrix.
    std::vector<std::vector<double> > d(n, std::vector<double>(m, 0.));
    for (unsigned int i = 0; i < n; ++i) {
      const double epsdif = eps * (1. + std::abs(par[i]));
      par[i] += 0.5 * epsdif;
      for (unsigned int j = 0; j < m; ++j) d[i][j] = f(x[j], par);
      par[i] -= epsdif;
      for (unsigned int j = 0; j < m; ++j) {
        d[i][j] = (d[i][j] - f(x[j], par)) / (epsdif * ey[j]);
      }
      par[i] += 0.5 * epsdif;
    }
    // Invert the matrix in Householder style.
    std::vector<double> colsum(n, 0.);
    std::vector<int> pivot(n, 0);
    for (unsigned int i = 0; i < n; ++i) {
      colsum[i] = std::inner_product(d[i].cbegin(), d[i].cend(), 
                                     d[i].cbegin(), 0.);
      pivot[i] = i;
    }
    // Decomposition.
    std::vector<double> alpha(n, 0.);
    bool singular = false;
    for (unsigned int k = 0; k < n; ++k) {
      double sigma = colsum[k];
      unsigned int jbar = k;
      for (unsigned int j = k + 1; j < n; ++j) {
        if (sigma < colsum[j]) {
          sigma = colsum[j];
          jbar = j;
        }
      }
      if (jbar != k) {
        // Interchange columns.
        std::swap(pivot[k], pivot[jbar]);
        std::swap(colsum[k], colsum[jbar]);
        std::swap(d[k], d[jbar]);
      }
      sigma = 0.;
      for (unsigned int i = k; i < m; ++i) sigma += d[k][i] * d[k][i]; 
      if (sigma == 0. || sqrt(sigma) < 1.e-8 * std::abs(d[k][k])) {
        singular = true;
        break;
      }
      alpha[k] = d[k][k] < 0. ? sqrt(sigma) : -sqrt(sigma);
      const double beta = 1. / (sigma - d[k][k] * alpha[k]);
      d[k][k] -= alpha[k];
      std::vector<double> b(n, 0.);
      for (unsigned int j = k + 1; j < n; ++j) {
        for (unsigned int i = k; i < n; ++i) b[j] += d[k][i] * d[j][i];
        b[j] *= beta;
      }
      for (unsigned int j = k + 1; j < n; ++j) {
        for (unsigned int i = k; i < m; ++i) {
          d[j][i] -= d[k][i] * b[j];
          colsum[j] -= d[j][k] * d[j][k];
        }
      }
    }
    if (singular) {
      std::cerr << "LeastSquaresFit: Householder matrix (nearly) singular;\n"
                << "    no further optimisation.\n"
                << "    Ensure the function depends on the parameters\n"
                << "    and try to supply reasonable starting values.\n";
      break;
    }
    // Solve.
    for (unsigned int j = 0; j < n; ++j) {
      double gamma = 0.;
      for (unsigned int i = j; i < m; ++i) gamma += d[j][i] * r[i];
      gamma *= 1. / (alpha[j] * d[j][j]);
      for (unsigned int i = j; i < m; ++i) r[i] += gamma * d[j][i];
    }
    std::vector<double> z(n, 0.);
    z[n - 1] = r[n - 1] / alpha[n - 1];
    for (int i = n - 1; i >= 1; --i) {
      double sum = 0.;
      for (unsigned int j = i + 1; j <= n; ++j) {
        sum += d[j - 1][i - 1] * z[j - 1];
      } 
      z[i - 1] = (r[i - 1] - sum) / alpha[i - 1]; 
    }
    // Correction vector.
    std::vector<double> s(n, 0.);
    for (unsigned int i = 0; i < n; ++i) s[pivot[i]] = z[i];
    // Generate some debugging output.
    if (debug) {
      std::cout << "  Correction vector in iteration " << iter << ":\n";
      for (unsigned int i = 0; i < n; ++i) {
        std::printf("    %5u  %15.8e\n", i, s[i]);
      }
    }
    // Add part of the correction vector to the estimate to improve chi2.
    chi2L = chi2;
    chi2 *= 2;
    for (unsigned int i = 0; i < n; ++i) par[i] += s[i] * 2;
    for (unsigned int i = 0; i <= 10; ++i) {
      if (chi2 <= chi2L) break;
      if (std::abs(chi2L - chi2) < eps * chi2) {
        if (debug) {
          std::cout << "    Too little improvement; reduction loop halted.\n";
        }
        break;
      }
      chi2 = 0.;
      const double scale = 1. / pow(2, i);
      for (unsigned int j = 0; j < n; ++j) par[j] -= s[j] * scale;
      for (unsigned int j = 0; j < m; ++j) {
        r[j] = (y[j] - f(x[j], par)) / ey[j];
        chi2 += r[j] * r[j];
      }
      if (debug) {
        std::printf("    Reduction loop %3i: chi2 = %15.8e\n", i, chi2);
      }
    }
    // Calculate the max. norm.
    diffc = std::abs(r[0]);
    for (unsigned int i = 1; i < m; ++i) {
      diffc = std::max(std::abs(r[i]), diffc);
    }
    // Print some debugging output.
    if (debug) {
      std::cout << "  Values of the fit parameters after iteration " 
                << iter << "\n    Parameter            Value\n";
      for (unsigned int i = 0; i < n; ++i) {
        std::printf("    %9u  %15.8e\n", i, par[i]);
      }
      std::printf("  for which chi2 = %15.8e and diff = %15.8e\n", 
                  chi2, diffc);
    } else if (verbose) {
      std::printf("  Step %3u: largest deviation = %15.8e, chi2 = %15.8e\n",
                  iter, diffc, chi2);
    }
  }
  // End of fit, perform error calculation.
  if (!converged) {
    std::cerr << "LeastSquaresFit: Maximum number of iterations reached.\n";
  }
  // Calculate the derivative matrix for the final settings.
  std::vector<std::vector<double> > d(n, std::vector<double>(m, 0.));
  for (unsigned int i = 0; i < n; ++i) {
    const double epsdif = eps * (1. + std::abs(par[i]));
    par[i] += 0.5 * epsdif;
    for (unsigned int j = 0; j < m; ++j) d[i][j] = f(x[j], par);
    par[i] -= epsdif;
    for (unsigned int j = 0; j < m; ++j) {
      d[i][j] = (d[i][j] - f(x[j], par)) / (epsdif * ey[j]);
    }
    par[i] += 0.5 * epsdif;
  }
  // Calculate the error matrix.
  std::vector<std::vector<double> > cov(n, std::vector<double>(n, 0.));
  for (unsigned int i = 0; i < n; ++i) {
    for (unsigned int j = 0; j < n; ++j) {
      cov[i][j] = std::inner_product(d[i].cbegin(), d[i].cend(), 
                                     d[j].cbegin(), 0.);
    }
  }
  // Compute the scaling factor for the errors.
  double scale = m > n ? chi2 / (m - n) : 1.;
  // Invert it to get the covariance matrix.
  epar.assign(n, 0.);
  if (Garfield::Numerics::CERNLIB::dinv(n, cov) != 0) {
    std::cerr << "LeastSquaresFit: Singular covariance matrix; "
              << "no error calculation.\n";
  } else {
    for (unsigned int i = 0; i < n; ++i) {
      for (unsigned int j = 0; j < n; ++j) cov[i][j] *= scale;
      epar[i] = sqrt(std::max(0., cov[i][i]));
    }
  }
  // Print results.
  if (debug) {
    std::cout << "  Comparison between input and fit:\n"
              << "            X            Y      F(X)\n";
    for (unsigned int i = 0; i < m; ++i) {
      std::printf(" %5u %15.8e %15.8e %15.8e\n", i, x[i], y[i], f(x[i], par));
    }
  }
  if (verbose) {
    std::cout << "  Final values of the fit parameters:\n"
              << "    Parameter            Value            Error\n";
    for (unsigned int i = 0; i < n; ++i) {
      std::printf("    %9u  %15.8e  %15.8e\n", i, par[i], epar[i]);
    }
    std::cout << "  The errors have been scaled by a factor of "
              << sqrt(scale) << ".\n";
    std::cout << "  Covariance matrix:\n";
    for (unsigned int i = 0; i < n; ++i) {
      for (unsigned int j = 0; j < n; ++j) {
        std::printf(" %15.8e", cov[i][j]);
      }
      std::cout << "\n";
    }
    std::cout << "  Correlation matrix:\n";
    for (unsigned int i = 0; i < n; ++i) {
      for (unsigned int j = 0; j < n; ++j) {
        double cor = 0.;
        if (cov[i][i] > 0. && cov[j][j] > 0.) {
          cor = cov[i][j] / sqrt(cov[i][i] * cov[j][j]);
        }
        std::printf(" %15.8e", cor);
      }
      std::cout << "\n";
    }
    std::cout << "  Minimisation finished.\n";
  }
  return converged;
}

}
}
