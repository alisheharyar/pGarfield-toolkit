/*
(c) 2005 Supratik Mukhopadhyay, Nayana Majumdar
*/

#define DEFINE_ISLESGLOBAL

#include "Isles.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_sf.h>

#define SHIFT 2.0
#define ARMAX 10000.0  // Maximum aspect ratio for an element
#define ARMIN 0.0001   // Minimum aspect ratio for an element

// #define XNSegApprox 10	// much less time but inaccurate results
// #define ZNSegApprox 10	// much less time but inaccurate results
#define XNSegApprox 100  // barely acceptable results
#define ZNSegApprox 100  // barely acceptable result
// #define XNSegApprox 1000	// much better results but lot more time
// #define ZNSegApprox 1000	// much better results but lot more time

#ifdef __cplusplus
#include <cmath>
#ifdef _GLIBCXX_HAVE_OBSOLETE_ISINF_ISNAN
using std::isinf;
using std::isnan;
#endif
namespace neBEM {
#endif

// Exact potential and flux for a rectangular element
// Potential needs all the calculations, whereas, the fluxes need part of the
// information generated during the computation for potential. This is the
// reason why we have included flux computation along with the potential
// computation.
// Expressions from:
int ExactRecSurf(double X, double Y, double Z, double xlo, double zlo,
                 double xhi, double zhi, double *Potential, Vector3D *Flux) {
  if (DebugISLES) printf("In ExactRecSurf ...\n");

  ++IslesCntr;
  ++ExactCntr;

  ApproxFlag = 0;

  if ((fabs(xhi - xlo) < 3.0 * MINDIST) || (fabs(zhi - zlo) < 3.0 * MINDIST)) {
    fprintf(stdout, "Element size too small! ... returning ...\n");
    return -1;
  }

  double ar = fabs((zhi - zlo) / (xhi - xlo));
  if (ar > ARMAX || ar < ARMIN) {
    fprintf(stdout, "Element too thin! ... returning ...\n");
    return -2;
  }

  if (fabs(X) < MINDIST) X = 0.0;
  if (fabs(Y) < MINDIST) Y = 0.0;
  if (fabs(Z) < MINDIST) Z = 0.0;

  double dxlo = X - xlo;  // zero at the X=xlo edge
  if (fabs(dxlo) < MINDIST) dxlo = 0.0;
  double dxhi = X - xhi;  // zero at the X=xhi edge
  if (fabs(dxhi) < MINDIST) dxhi = 0.0;
  double dzlo = Z - zlo;  // zero at the Z=zlo edge
  if (fabs(dzlo) < MINDIST) dzlo = 0.0;
  double dzhi = Z - zhi;  // zero at the Z=zhi edge
  if (fabs(dzhi) < MINDIST) dzhi = 0.0;

  // These four parameters can never be zero except at the four corners where
  // one of them can become zero. For example, at X=xlo, Y=0, Z=zlo, D11
  // is zero but the others are nonzero.
  double D11 = sqrt(dxlo * dxlo + Y * Y + dzlo * dzlo);
  if (fabs(D11) < MINDIST) D11 = 0.0;
  double D21 = sqrt(dxhi * dxhi + Y * Y + dzlo * dzlo);
  if (fabs(D21) < MINDIST) D21 = 0.0;
  double D12 = sqrt(dxlo * dxlo + Y * Y + dzhi * dzhi);
  if (fabs(D12) < MINDIST) D12 = 0.0;
  double D22 = sqrt(dxhi * dxhi + Y * Y + dzhi * dzhi);
  if (fabs(D22) < MINDIST) D22 = 0.0;

  // Parameters related to the Y terms
  int S1 = Sign(dzlo);
  int S2 = Sign(dzhi);
  double modY = fabs(Y);
  int SY = Sign(Y);
  double I1 = dxlo * modY;
  double I2 = dxhi * modY;
  double R1 = Y * Y + dzlo * dzlo;
  double R2 = Y * Y + dzhi * dzhi;
  if (fabs(I1) < MINDIST2) I1 = 0.0;
  if (fabs(I2) < MINDIST2) I2 = 0.0;
  if (fabs(R1) < MINDIST2) R1 = 0.0;
  if (fabs(R2) < MINDIST2) R2 = 0.0;

  if (DebugISLES) {
    fprintf(stdout, "X: %.16lg, Y: %.16lg, Z: %.16lg\n", X, Y, Z);
    fprintf(stdout, "xlo: %.16lg, zlo: %.16lg, xhi: %.16lg, zhi: %.16lg\n", xlo,
            zlo, xhi, zhi);
    fprintf(stdout, "dxlo: %.16lg, dzlo: %.16lg, dxhi: %.16lg, dzhi: %.16lg\n",
            dxlo, dzlo, dxhi, dzhi);
    fprintf(stdout, "D11: %.16lg, D12: %.16lg, D21: %.16lg, D22: %.16lg\n", D11,
            D12, D21, D22);
    fprintf(stdout, "S1: %d, S2: %d, modY: %.16lg\n", S1, S2, modY);
    fprintf(stdout, "I1: %.16lg, I2: %.16lg, R1: %.16lg, R2: %.16lg\n", I1, I2,
            R1, R2);
    fprintf(stdout, "MINDIST: %.16lg, MINDIST2: %.16lg, SHIFT: %.16lg\n",
            MINDIST, MINDIST2, SHIFT);
    fflush(stdout);
  }

  // Check for possible numerical difficuties and take care.
  // Presently the idea is to shift the field point slightly to a 'safe'
  // position. Note that a shift in Y does not work because the singularities
  // are associated with D11-s and dxlo-s and their combinations. A Y shift can
  // sometimes alleviate the problem, but it cannot gurantee a permanent1s
  // solution.
  // A better approach is to evaluate analytic expressions truly valid for
  // these critical regions, especially to take care of the >= and <=
  // comparisons. For the == case, both of the previous comparisons are true and
  // it is hard to justify one choice over the other. On the other hand, there
  // are closed form analytic solutions for the == cases, and the problem does
  // not stem from round-off errors as in the > or <, but not quite == cases.
  // Note that shifts that are exactly equal to MINDIST, can give rise to
  // un-ending recursions, leading to Segmentation fault.
  // Corners
  // Four corners
  // Averages over four point evaluations may be carried out (skipped at
  // present) This, however, may be difficult since we'll have to make sure that
  // the shift towards the element does not bring the point too close to the
  // same, or another difficult-to-evaluate situation. One of the ways to ensure
  // this is to make SHIFT large enough, but that is unreasnoable and will
  // introduce large amount of error.
  if ((fabs(D11) <= MINDIST)) {
    // close to xlo, 0, zlo
    if (DebugISLES) printf("fabs(D11) <= MINDIST ... ");
    double X1 = X;
    double Z1 = Z;
    if ((X >= xlo) && (Z >= zlo)) {
      // point on the element
      if (DebugISLES) printf("Case 1\n");
      X1 = X + SHIFT * MINDIST;
      Z1 = Z + SHIFT * MINDIST;
    } else if ((X <= xlo) && (Z >= zlo)) {
      // field point outside the element
      if (DebugISLES) printf("Case 2 ...\n");
      X1 = X - SHIFT * MINDIST;
      Z1 = Z + SHIFT * MINDIST;
    } else if ((X >= xlo) && (Z <= zlo)) {
      // field point outside the element
      if (DebugISLES) printf("Case 3 ...\n");
      X1 = X + SHIFT * MINDIST;
      Z1 = Z - SHIFT * MINDIST;
    } else if ((X <= xlo) && (Z <= zlo)) {
      // field point outside the element
      if (DebugISLES) printf("Case 4 ...\n");
      X1 = X - SHIFT * MINDIST;
      Z1 = Z - SHIFT * MINDIST;
    }
    double Pot1;
    Vector3D Flux1;
    ExactRecSurf(X1, Y, Z1, xlo, zlo, xhi, zhi, &Pot1, &Flux1);
    *Potential = Pot1;
    Flux->X = Flux1.X;
    Flux->Y = Flux1.Y;
    Flux->Z = Flux1.Z;
    return 0;
  }
  if ((fabs(D21) <= MINDIST)) {
    // close to xhi, 0, zlo
    if (DebugISLES) printf("fabs(D21) <= MINDIST ... ");
    double X1 = X;
    double Z1 = Z;
    if ((X >= xhi) && (Z >= zlo)) {
      // point outside the element
      if (DebugISLES) printf("Case 1 ...\n");
      X1 = X + SHIFT * MINDIST;
      Z1 = Z + SHIFT * MINDIST;
    } else if ((X <= xhi) && (Z >= zlo)) {
      // point on the element
      if (DebugISLES) printf("Case 2 ...\n");
      X1 = X - SHIFT * MINDIST;
      Z1 = Z + SHIFT * MINDIST;
    } else if ((X >= xhi) && (Z <= zlo)) {
      // field point outside the element
      if (DebugISLES) printf("Case 3 ...\n");
      X1 = X + SHIFT * MINDIST;
      Z1 = Z - SHIFT * MINDIST;
    } else if ((X <= xhi) && (Z <= zlo)) {
      // field point outside the element
      if (DebugISLES) printf("Case 4 ...\n");
      X1 = X - SHIFT * MINDIST;
      Z1 = Z - SHIFT * MINDIST;
    }
    double Pot1;
    Vector3D Flux1;
    ExactRecSurf(X1, Y, Z1, xlo, zlo, xhi, zhi, &Pot1, &Flux1);
    *Potential = Pot1;
    Flux->X = Flux1.X;
    Flux->Y = Flux1.Y;
    Flux->Z = Flux1.Z;
    return 0;
  }
  if ((fabs(D12) <= MINDIST)) {
    // close to xlo, 0, zhi
    if (DebugISLES) printf("fabs(D12) <= MINDIST ... ");
    double X1 = X;
    double Z1 = Z;
    if ((X >= xlo) && (Z >= zhi)) {
      // point outside the element
      if (DebugISLES) printf("Case 1 ...\n");
      X1 = X + SHIFT * MINDIST;
      Z1 = Z + SHIFT * MINDIST;
    } else if ((X <= xlo) && (Z >= zhi)) {
      // field point outside the element
      if (DebugISLES) printf("Case 2 ...\n");
      X1 = X - SHIFT * MINDIST;
      Z1 = Z + SHIFT * MINDIST;
    } else if ((X >= xlo) && (Z <= zhi)) {
      // field point on the element
      if (DebugISLES) printf("Case 3 ...\n");
      X1 = X + SHIFT * MINDIST;
      Z1 = Z - SHIFT * MINDIST;
    } else if ((X <= xlo) && (Z <= zhi)) {
      // field point outside the element
      if (DebugISLES) printf("Case 4 ...\n");
      X1 = X - SHIFT * MINDIST;
      Z1 = Z - SHIFT * MINDIST;
    }
    double Pot1;
    Vector3D Flux1;
    ExactRecSurf(X1, Y, Z1, xlo, zlo, xhi, zhi, &Pot1, &Flux1);
    *Potential = Pot1;
    Flux->X = Flux1.X;
    Flux->Y = Flux1.Y;
    Flux->Z = Flux1.Z;
    return 0;
  }
  if ((fabs(D22) <= MINDIST)) {
    // close to xhi, 0, zhi
    if (DebugISLES) printf("fabs(D22) <= MINDIST ... ");
    double X1 = X;
    double Z1 = Z;
    if ((X >= xhi) && (Z >= zhi)) {
      // point outside the element
      if (DebugISLES) printf("Case 1 ...\n");
      X1 = X + SHIFT * MINDIST;
      Z1 = Z + SHIFT * MINDIST;
    } else if ((X <= xhi) && (Z >= zhi)) {
      // field point outside the element
      if (DebugISLES) printf("Case 2 ...\n");
      X1 = X - SHIFT * MINDIST;
      Z1 = Z + SHIFT * MINDIST;
    } else if ((X >= xhi) && (Z <= zhi)) {
      // field point outside the element
      if (DebugISLES) printf("Case 3 ...\n");
      X1 = X + SHIFT * MINDIST;
      Z1 = Z - SHIFT * MINDIST;
    } else if ((X <= xhi) && (Z <= zhi)) {
      // field point on the element
      if (DebugISLES) printf("Case 4 ...\n");
      X1 = X - SHIFT * MINDIST;
      Z1 = Z - SHIFT * MINDIST;
    }
    double Pot1;
    Vector3D Flux1;
    ExactRecSurf(X1, Y, Z1, xlo, zlo, xhi, zhi, &Pot1, &Flux1);
    *Potential = Pot1;
    Flux->X = Flux1.X;
    Flux->Y = Flux1.Y;
    Flux->Z = Flux1.Z;
    return 0;
  }
  // Four edges	- average over two points on two sides of the edge may be ok.
  // Here also we'll have to make sure that the shift
  // towards the element does not bring the point too close to the same, or
  // another difficult-to-evaluate situation. One of the ways to ensure this
  // is to make SHIFT large enough, but that is unreasonable and will introduce
  // large amount of error.
  if (fabs(dxlo) < MINDIST) {
    // edge at x=xlo || to Z - axis
    if (DebugISLES) printf("fabs(dxlo) < MINDIST ... ");
    double X1 = X;
    double X2 = X;
    if (X >= xlo) {
      // field point on +ve side of YZ plane
      if (DebugISLES) printf("Case 1 ...\n");
      X1 += SHIFT * MINDIST;
      X2 -= SHIFT * MINDIST;
    } else {
      // field point on -ve side of YZ plane
      if (DebugISLES) printf("Case 2 ...\n");
      X1 -= SHIFT * MINDIST;
      X2 += SHIFT * MINDIST;
    }
    double Pot1, Pot2;
    Vector3D Flux1, Flux2;
    ExactRecSurf(X1, Y, Z, xlo, zlo, xhi, zhi, &Pot1, &Flux1);
    ExactRecSurf(X2, Y, Z, xlo, zlo, xhi, zhi, &Pot2, &Flux2);
    *Potential = 0.5 * (Pot1 + Pot2);
    Flux->X = 0.5 * (Flux1.X + Flux2.X);
    Flux->Y = 0.5 * (Flux1.Y + Flux2.Y);
    Flux->Z = 0.5 * (Flux1.Z + Flux2.Z);
    return 0;
  }
  if (fabs(dzlo) < MINDIST) {
    // edge at z=zlo, || to X axis
    if (DebugISLES) printf("fabs(dzlo) < MINDIST ... ");
    double Z1 = Z;
    double Z2 = Z;
    if (Z >= zlo) {
      // field point on +ve side of XY plane
      if (DebugISLES) printf("Case 1 ...\n");
      Z1 += SHIFT * MINDIST;
      Z2 -= SHIFT * MINDIST;
    } else {
      // field point on -ve side of XY plane
      if (DebugISLES) printf("Case 2 ...\n");
      Z1 -= SHIFT * MINDIST;
      Z2 += SHIFT * MINDIST;
    }
    double Pot1, Pot2;
    Vector3D Flux1, Flux2;
    ExactRecSurf(X, Y, Z1, xlo, zlo, xhi, zhi, &Pot1, &Flux1);
    ExactRecSurf(X, Y, Z2, xlo, zlo, xhi, zhi, &Pot2, &Flux2);
    *Potential = 0.5 * (Pot1 + Pot2);
    Flux->X = 0.5 * (Flux1.X + Flux2.X);
    Flux->Y = 0.5 * (Flux1.Y + Flux2.Y);
    Flux->Z = 0.5 * (Flux1.Z + Flux2.Z);
    return 0;
  }
  if (fabs(dxhi) < MINDIST) {
    // edge at x=xhi, || to Z axis
    if (DebugISLES) printf("fabs(dxhi) < MINDIST ... ");
    double X1 = X;
    double X2 = X;
    if (X >= xhi) {
      // field point on +ve side of YZ plane
      if (DebugISLES) printf("Case 1 ...\n");
      X1 += SHIFT * MINDIST;
      X2 -= SHIFT * MINDIST;
    } else {
      // field point on -ve side of YZ plane
      if (DebugISLES) printf("Case 2 ...\n");
      X1 -= SHIFT * MINDIST;
      X2 += SHIFT * MINDIST;
    }
    double Pot1, Pot2;
    Vector3D Flux1, Flux2;
    ExactRecSurf(X1, Y, Z, xlo, zlo, xhi, zhi, &Pot1, &Flux1);
    ExactRecSurf(X2, Y, Z, xlo, zlo, xhi, zhi, &Pot2, &Flux2);
    *Potential = 0.5 * (Pot1 + Pot2);
    Flux->X = 0.5 * (Flux1.X + Flux2.X);
    Flux->Y = 0.5 * (Flux1.Y + Flux2.Y);
    Flux->Z = 0.5 * (Flux1.Z + Flux2.Z);
    return 0;
  }
  if (fabs(dzhi) < MINDIST) {
    // edge at z=zhi || to X axis
    if (DebugISLES) printf("fabs(dzhi) < MINDIST ... ");
    double Z1 = Z;
    double Z2 = Z;
    if (Z >= zhi) {
      // field point on +ve side of XY plane
      if (DebugISLES) printf("Case 1 ...\n");
      Z1 += SHIFT * MINDIST;
      Z2 -= SHIFT * MINDIST;
    } else {
      // field point on -ve side of XY plane
      if (DebugISLES) printf("Case 2 ...\n");
      Z1 -= SHIFT * MINDIST;
      Z2 += SHIFT * MINDIST;
    }
    double Pot1, Pot2;
    Vector3D Flux1, Flux2;
    ExactRecSurf(X, Y, Z1, xlo, zlo, xhi, zhi, &Pot1, &Flux1);
    ExactRecSurf(X, Y, Z2, xlo, zlo, xhi, zhi, &Pot2, &Flux2);
    *Potential = 0.5 * (Pot1 + Pot2);
    Flux->X = 0.5 * (Flux1.X + Flux2.X);
    Flux->Y = 0.5 * (Flux1.Y + Flux2.Y);
    Flux->Z = 0.5 * (Flux1.Z + Flux2.Z);
    return 0;
  }

  // Logarithmic weak singularities are possible.
  // Checks to be performed for 0 or -ve denominators and also
  // 0 and +ve numerators.
  // Interestingly, 0/0 does not cause a problem.
  double DZTerm1 = log((D11 - dzlo) / (D12 - dzhi));
  double DZTerm2 = log((D21 - dzlo) / (D22 - dzhi));
  double DXTerm1 = log((D11 - dxlo) / (D21 - dxhi));
  double DXTerm2 = log((D12 - dxlo) / (D22 - dxhi));

  if (DebugISLES) {
    fprintf(
        stdout,
        "DZTerm1: %.16lg, DZTerm2: %.16lg, DXTerm1: %.16lg, DXTerm2: %.16lg\n",
        DZTerm1, DZTerm2, DXTerm1, DXTerm2);
  }
  // Four conditions based on the logarithmic terms
  if (isnan(DZTerm1) || isinf(DZTerm1)) {
    ++FailureCntr;
    --ExactCntr;
    ApproxFlag = 1;
    fprintf(fIsles, "DZTerm1 problem ... approximating: %d.\n", ApproxFlag);
    if (DebugISLES)
      fprintf(stdout, "DZTerm1 problem ... approximating: %d.\n", ApproxFlag);
    return (ApproxRecSurf(X, Y, Z, xlo, zlo, xhi, zhi, XNSegApprox, ZNSegApprox,
                          Potential, Flux));
  }
  if (isnan(DZTerm2) || isinf(DZTerm2)) {
    ++FailureCntr;
    --ExactCntr;
    ApproxFlag = 2;
    fprintf(fIsles, "DZTerm2 problem ... approximating: %d.\n", ApproxFlag);
    if (DebugISLES)
      fprintf(stdout, "DZTerm2 problem ... approximating: %d.\n", ApproxFlag);
    return (ApproxRecSurf(X, Y, Z, xlo, zlo, xhi, zhi, XNSegApprox, ZNSegApprox,
                          Potential, Flux));
  }
  if (isnan(DXTerm1) || isinf(DXTerm1)) {
    ++FailureCntr;
    --ExactCntr;
    ApproxFlag = 3;
    fprintf(fIsles, "DXTerm1 problem ... approximating: %d.\n", ApproxFlag);
    if (DebugISLES)
      fprintf(stdout, "DXTerm1 problem ... approximating: %d.\n", ApproxFlag);
    return (ApproxRecSurf(X, Y, Z, xlo, zlo, xhi, zhi, XNSegApprox, ZNSegApprox,
                          Potential, Flux));
  }
  if (isnan(DXTerm2) || isinf(DXTerm2)) {
    ++FailureCntr;
    --ExactCntr;
    ApproxFlag = 4;
    fprintf(fIsles, "DXTerm2 problem ... approximating: %d.\n", ApproxFlag);
    if (DebugISLES)
      fprintf(stdout, "DXTerm2 problem ... approximating: %d.\n", ApproxFlag);
    return (ApproxRecSurf(X, Y, Z, xlo, zlo, xhi, zhi, XNSegApprox, ZNSegApprox,
                          Potential, Flux));
  }
  // Four criticalities based on the tanhyperbolic terms
  // Since gsl_complex_arctanh_real can have any real number as its argument,
  // all the criticalities are related to gsl_complex_arctanh where the
  // imaginary component is present. So, fabs(I1) > mindist and
  // fabs(I2) > mindist are only being tested.
  if (S1 != 0) {
    if (fabs(I1) > MINDIST2) {
      if (D11 * fabs(dzlo) < MINDIST2) {
        ++FailureCntr;
        --ExactCntr;
        ApproxFlag = 5;
        fprintf(fIsles, "S1-I1 problem ... approximating: %d.\n", ApproxFlag);
        if (DebugISLES)
          fprintf(stdout, "S1-I1 problem ... approximating: %d.\n", ApproxFlag);
        return (ApproxRecSurf(X, Y, Z, xlo, zlo, xhi, zhi, XNSegApprox,
                              ZNSegApprox, Potential, Flux));
      }
    }
    if (fabs(I2) > MINDIST2) {
      if (D21 * fabs(dzlo) < MINDIST2) {
        ++FailureCntr;
        --ExactCntr;
        ApproxFlag = 6;
        fprintf(fIsles, "S1-I2 problem ... approximating: %d.\n", ApproxFlag);
        if (DebugISLES)
          fprintf(stdout, "S1-I2 problem ... approximating: %d.\n", ApproxFlag);
        return (ApproxRecSurf(X, Y, Z, xlo, zlo, xhi, zhi, XNSegApprox,
                              ZNSegApprox, Potential, Flux));
      }
    }
  }
  if (S2 != 0) {
    if (fabs(I1) > MINDIST2) {
      if (D12 * fabs(dzhi) < MINDIST2) {
        ++FailureCntr;
        --ExactCntr;
        ApproxFlag = 7;
        fprintf(fIsles, "S2-I1 problem ... approximating: %d.\n", ApproxFlag);
        if (DebugISLES)
          fprintf(stdout, "S2-I1 problem ... approximating: %d.\n", ApproxFlag);
        return (ApproxRecSurf(X, Y, Z, xlo, zlo, xhi, zhi, XNSegApprox,
                              ZNSegApprox, Potential, Flux));
      }
    }
    if (fabs(I2) > MINDIST2) {
      if (D22 * fabs(dzhi) < MINDIST2) {
        ++FailureCntr;
        --ExactCntr;
        ApproxFlag = 8;
        fprintf(fIsles, "S2-I2 problem ... approximating: %d.\n", ApproxFlag);
        if (DebugISLES)
          fprintf(stdout, "S2-I2 problem ... approximating: %d.\n", ApproxFlag);
        return (ApproxRecSurf(X, Y, Z, xlo, zlo, xhi, zhi, XNSegApprox,
                              ZNSegApprox, Potential, Flux));
      }
    }
  }

  double sumTanTerms = 0.;
  // The possibility of singularities for dzhi or dzlo (division by zero)
  // is overridden by the fact that S1 or S2 becomes zero in such cases
  // and the singularity is avoided.
  if (S1 != 0) {
    if (fabs(I1) > MINDIST2) {
      double a = D11 * D11 * dzlo * dzlo;
      double b = R1 * R1 + I1 * I1;
      double tmp = -S1 * atan(2 * I1 * D11 * fabs(dzlo) / (a - b));
      if (b > a) {
        if ((X > xlo && Z > zlo) || (X < xlo && Z < zlo)) {
          tmp -= ST_PI;
        } else if ((X < xlo && Z > zlo) || (X > xlo && Z < zlo)) {
          tmp += ST_PI;
        }
      }
      sumTanTerms += tmp;
    }

    if (fabs(I2) > MINDIST2) {
      double a = D21 * D21 * dzlo * dzlo;
      double b = R1 * R1 + I2 * I2;
      double tmp = -S1 * atan(2 * I2 * D21 * fabs(dzlo) / (a - b));
      if (b > a) {
        if ((X > xhi && Z > zlo) || (X < xhi && Z < zlo)) {
          tmp -= ST_PI;
        } else if ((X < xhi && Z > zlo) || (X > xhi && Z < zlo)) {
          tmp += ST_PI;
        }
      }
      sumTanTerms -= tmp;
    }
  }

  if (S2 != 0) {
    if (fabs(I1) > MINDIST2) {
      double a = D12 * D12 * dzhi * dzhi;
      double b = R2 * R2 + I1 * I1;
      double tmp = -S2 * atan(2 * I1 * D12 * fabs(dzhi) / (a - b));
      if (b > a) {
        if ((X > xlo && Z > zhi) || (X < xlo && Z < zhi)) {
          tmp -= ST_PI;
        } else if ((X < xlo && Z > zhi) || (X > xlo && Z < zhi)) {
          tmp += ST_PI;
        }
      }
      sumTanTerms -= tmp;
    }

    if (fabs(I2) > MINDIST2) {
      double a = D22 * D22 * dzhi * dzhi;
      double b = R2 * R2 + I2 * I2;
      double tmp = -S2 * atan(2 * I2 * D22 * fabs(dzhi) / (a - b));
      if (b > a) {
        if ((X > xhi && Z > zhi) || (X < xhi && Z < zhi)) {
          tmp -= ST_PI;
        } else if ((X < xhi && Z > zhi) || (X > xhi && Z < zhi)) {
          tmp += ST_PI;
        }
      }
      sumTanTerms += tmp;
    }
  }

  sumTanTerms *= -0.5;
  double Pot = -dxlo * DZTerm1 + dxhi * DZTerm2 + modY * sumTanTerms -
               dzlo * DXTerm1 + dzhi * DXTerm2;
  double Fx = DZTerm1 - DZTerm2;
  double Fy = -SY * sumTanTerms;
  double Fz = DXTerm1 - DXTerm2;
  if (DebugISLES) {
    printf("XTerms: %.16lg, YTerms: %.16lg, ZTerms: %.16lg\n",
           -dxlo * DZTerm1 + dxhi * DZTerm2, modY * sumTanTerms,
           -dzlo * DXTerm1 + dzhi * DXTerm2);
    printf("Pot: %lf, Fx: %lf, Fy: %lf, Fz: %lf\n", Pot, Fx, Fy, Fz);
    fflush(stdout);
  }

  // constants of integration
  // The only logic for the Fy constant seems to be the fact that the
  // potential has a negative of this constant
  if (((X > (xlo + MINDIST)) && (X < (xhi - MINDIST))) &&
      ((Z > (zlo + MINDIST)) && (Z < (zhi - MINDIST)))) {
    Pot -= 2.0 * modY * ST_PI;
    if (SY != 0)
      Fy += 2.0 * (double)SY * ST_PI;
    else
      Fy = 2.0 * ST_PI;
  }
  if (DebugISLES) {
    printf("Constants of integration added for potential and Fy.\n");
    printf("Pot: %lf, Fx: %lf, Fy: %lf, Fz: %lf\n", Pot, Fx, Fy, Fz);
    fflush(stdout);
  }

  // Error situations handled before returning the values
  if ((Pot < 0.0) || (isnan(Pot) || isinf(Pot))) {
    fprintf(fIsles, "\n--- Approximation in ExactRecSurf ---\n");
    fprintf(fIsles, "Negative, nan or inf potential ... approximating!\n");
    if (DebugISLES) {
      fprintf(stdout, "\n--- Approximation in ExactRecSurf ---\n");
      fprintf(stdout, "Negative, nan or inf potential ... approximating!\n");
    }
    fprintf(fIsles, "X: %.16lg, Y: %.16lg, Z: %.16lg\n", X, Y, Z);
    fprintf(fIsles, "xlo: %.16lg, zlo: %.16lg, xhi: %.16lg, zhi: %.16lg\n", xlo,
            zlo, xhi, zhi);
    fprintf(fIsles, "dxlo: %.16lg, dzlo: %.16lg, dxhi: %.16lg, dzhi: %.16lg\n",
            dxlo, dzlo, dxhi, dzhi);
    fprintf(fIsles, "D11: %.16lg, D12: %.16lg, D21: %.16lg, D22: %.16lg\n", D11,
            D12, D21, D22);
    fprintf(fIsles, "modY: %.16lg\n", modY);
    fprintf(fIsles, "I1: %.16lg, I2: %.16lg, R1: %.16lg, R2: %.16lg\n", I1, I2,
            R1, R2);
    fprintf(fIsles, "S1: %d, S2: %d\n", S1, S2);
    fprintf(fIsles, "Computed Pot: %.16lg\n", Pot);
    ApproxFlag = 13;
    ++FailureCntr;
    --ExactCntr;
    fprintf(fIsles, "Approximating: %d.\n", ApproxFlag);
    return (ApproxRecSurf(X, Y, Z, xlo, zlo, xhi, zhi, XNSegApprox, ZNSegApprox,
                          Potential, Flux));
  }
  if (isnan(Fx) || isinf(Fx)) {
    fprintf(fIsles, "\n--- Approximation in ExactRecSurf ---\n");
    fprintf(fIsles, "Nan or inf Fx ... approximating!\n");
    if (DebugISLES) {
      fprintf(stdout, "\n--- Approximation in ExactRecSurf ---\n");
      fprintf(stdout, "Nan or inf Fx ... approximating!\n");
    }
    fprintf(fIsles, "X: %.16lg, Y: %.16lg, Z: %.16lg\n", X, Y, Z);
    fprintf(fIsles, "xlo: %.16lg, zlo: %.16lg, xhi: %.16lg, zhi: %.16lg\n", xlo,
            zlo, xhi, zhi);
    fprintf(fIsles, "dxlo: %.16lg, dzlo: %.16lg, dxhi: %.16lg, dzhi: %.16lg\n",
            dxlo, dzlo, dxhi, dzhi);
    fprintf(fIsles, "D11: %.16lg, D12: %.16lg, D21: %.16lg, D22: %.16lg\n", D11,
            D12, D21, D22);
    fprintf(fIsles, "Computed Fx: %.16lg\n", Fx);
    ApproxFlag = 14;
    ++FailureCntr;
    --ExactCntr;
    fprintf(fIsles, "Approximating: %d.\n", ApproxFlag);
    return (ApproxRecSurf(X, Y, Z, xlo, zlo, xhi, zhi, XNSegApprox, ZNSegApprox,
                          Potential, Flux));
  }
  if (isnan(Fy) || isinf(Fy)) {
    fprintf(fIsles, "\n--- Approximation in ExactRecSurf ---\n");
    fprintf(fIsles, "Nan or inf Fy ... approximating!\n");
    if (DebugISLES) {
      fprintf(stdout, "\n--- Approximation in ExactRecSurf ---\n");
      fprintf(stdout, "Nan or inf Fy ... approximating!\n");
    }
    fprintf(fIsles, "X: %.16lg, Y: %.16lg, Z: %.16lg\n", X, Y, Z);
    fprintf(fIsles, "xlo: %.16lg, zlo: %.16lg, xhi: %.16lg, zhi: %.16lg\n", xlo,
            zlo, xhi, zhi);
    fprintf(fIsles, "dxlo: %.16lg, dzlo: %.16lg, dxhi: %.16lg, dzhi: %.16lg\n",
            dxlo, dzlo, dxhi, dzhi);
    fprintf(fIsles, "D11: %.16lg, D12: %.16lg, D21: %.16lg, D22: %.16lg\n", D11,
            D12, D21, D22);
    fprintf(fIsles, "S1: %d, S2: %d, SY: %d, modY: %.16lg\n", S1, S2, SY, modY);
    fprintf(fIsles, "I1: %.16lg, I2: %.16lg, R1: %.16lg, R2: %.16lg\n", I1, I2,
            R1, R2);
    fprintf(fIsles, "Computed Fy: %.16lg\n", Fy);
    ApproxFlag = 15;
    ++FailureCntr;
    --ExactCntr;
    fprintf(fIsles, "Approximating: %d.\n", ApproxFlag);
    return (ApproxRecSurf(X, Y, Z, xlo, zlo, xhi, zhi, XNSegApprox, ZNSegApprox,
                          Potential, Flux));
  }
  if (isnan(Fz) || isinf(Fz)) {
    fprintf(fIsles, "\n--- Approximation in ExactRecSurf ---\n");
    fprintf(fIsles, "Nan or inf Fz ... approximating!\n");
    if (DebugISLES) {
      fprintf(stdout, "\n--- Approximation in ExactRecSurf ---\n");
      fprintf(stdout, "Nan or inf Fz ... approximating!\n");
    }
    fprintf(fIsles, "X: %.16lg, Y: %.16lg, Z: %.16lg\n", X, Y, Z);
    fprintf(fIsles, "xlo: %.16lg, zlo: %.16lg, xhi: %.16lg, zhi: %.16lg\n", xlo,
            zlo, xhi, zhi);
    fprintf(fIsles, "dxlo: %.16lg, dzlo: %.16lg, dxhi: %.16lg, dzhi: %.16lg\n",
            dxlo, dzlo, dxhi, dzhi);
    fprintf(fIsles, "D11: %.16lg, D12: %.16lg, D21: %.16lg, D22: %.16lg\n", D11,
            D12, D21, D22);
    fprintf(fIsles, "Computed Fz: %.16lg\n", Fz);
    ApproxFlag = 16;
    ++FailureCntr;
    --ExactCntr;
    fprintf(fIsles, "Approximating: %d.\n", ApproxFlag);
    return (ApproxRecSurf(X, Y, Z, xlo, zlo, xhi, zhi, XNSegApprox, ZNSegApprox,
                          Potential, Flux));
  }

  *Potential = Pot;
  Flux->X = Fx;
  Flux->Y = Fy;
  Flux->Z = Fz;

  if (DebugISLES) printf("Going out of ExactRecSurf ...\n");

  return 0;
}  // ExactRecSurf ends

int ApproxRecSurf(double X, double Y, double Z, double xlo, double zlo,
                  double xhi, double zhi, int xseg, int zseg, double *Potential,
                  Vector3D *Flux) {
  if (DebugISLES) {
    printf("In ApproxRecSurf ...\n");
  }

  ++ApproxCntr;

  double dx = (xhi - xlo) / xseg;
  double dz = (zhi - zlo) / zseg;
  double xel = (xhi - xlo) / xseg;
  double zel = (zhi - zlo) / zseg;
  double diag = sqrt(dx * dx + dz * dz);
  double area = xel * zel;

  double Pot = 0., XFlux = 0., YFlux = 0., ZFlux = 0.;

  if (area > MINDIST2) {  // else not necessary
    for (int i = 1; i <= xseg; ++i) {
      double xi = xlo + (dx / 2.0) + (i - 1) * dx;
      for (int k = 1; k <= zseg; ++k) {
        double zk = zlo + (dz / 2.0) + (k - 1) * dz;

        double dist = sqrt((X - xi) * (X - xi) + Y * Y + (Z - zk) * (Z - zk));
        if (DebugISLES) printf("dist: %lg\n", dist);
        if (dist >= diag) {
          Pot += area / dist;
        } else if (dist <= MINDIST) {
          // Self influence
          Pot += 2.0 * (xel * log((zel + sqrt(xel * xel + zel * zel)) / xel) +
                        zel * log((xel + sqrt(xel * xel + zel * zel)) / zel));
        } else {
          // in the intermediate region where diag > dist > MINDIST
          Pot += area / diag;  // replace by expression of self-influence
          if (DebugISLES) printf("Special Pot: %lg\n", area / diag);
        }

        if (dist >= diag) {
          double f = area / (dist * dist * dist);
          XFlux += f * (X - xi);
          YFlux += f * Y;
          ZFlux += f * (Z - zk);
        } else {
          double f = area / (diag * diag * diag);
          XFlux += f * (X - xi);
          YFlux += f * Y;
          ZFlux += f * (Z - zk);
          if (DebugISLES) {
            printf("Special XFlux: %lg, YFlux: %lg, ZFlux: %lg\n", f * (X - xi),
                   f * Y, f * (Z - zk));
          }
        }  // else dist >= diag
      }    // zseg
    }      // xseg
  }        // if area > MINDIST2

  *Potential = Pot;
  Flux->X = XFlux;
  Flux->Y = YFlux;
  Flux->Z = ZFlux;

  return 0;
}  // ApproxRecSurf ends

// Exact potential for a triangular element (a right-angled
// triangle having (0,0), (1,0), (0,zMax) as its vertices).
// It is assumed that the affine transformation has been carried out separately
// and necessary adjustments to X,Y,Z have been carried out and supplied to
// the following function.
// Potential needs all the calculations, whereas, the fluxes need part of the
// information generated during the computation for potential. This is the
// reason why we have included flux computation along with the potential
// computation.
// Parameters like dxlo (X), dxhi (X-1.0), dzlo (Z), dzhi (Z-zMax) can be used
// as done in ExactRectSurf
// Expressions from:
int ExactTriSurf(double zMax, double X, double Y, double Z, double *Potential,
                 Vector3D *Flux) {
  if (DebugISLES) printf("In ExactTriSurf ...\n");

  ++IslesCntr;
  ++ExactCntr;
  // Reset the flag to indicate approximate evaluation of potential.
  ApproxFlag = 0;

  // We do not need to check for X, since element extent is always 0 - 1.
  if (zMax < 3.0 * SHIFT * MINDIST) {
    // should allow enough space for Z corrections
    // One SHIFT should not lead to another criticality
    fprintf(stdout, "Element size too small! ... returning ...\n");
    return -1;
  }

  if (zMax > ARMAX || zMax < ARMIN) {
    fprintf(stdout, "Element too thin! ... returning ...\n");
    return -1;
  }

  // These three parameters can never be zero except at the three corners where
  // one of them can become zero. For example, at X=0, Y=0, Z=0, D11
  // is zero but the others are nonzero.
  if (fabs(X) < MINDIST) X = 0.0;
  if (fabs(Y) < MINDIST) Y = 0.0;
  if (fabs(Z) < MINDIST) Z = 0.0;
  double modY = fabs(Y);
  if (modY < MINDIST) modY = 0.0;
  int S1 = Sign(Z);
  int SY = Sign(Y);

  // distances from corners (0,0,0), (1,0,0) and (0,0,zMax)
  double D11 = sqrt(X * X + Y * Y + Z * Z);
  if (D11 < MINDIST) D11 = 0.0;
  double D21 = sqrt((X - 1.0) * (X - 1.0) + Y * Y + Z * Z);
  if (D21 < MINDIST) D21 = 0.0;
  double D12 = sqrt(X * X + Y * Y + (Z - zMax) * (Z - zMax));
  if (D12 < MINDIST) D12 = 0.0;

  double G = zMax * (X - 1.0) + Z;
  if (fabs(G) < MINDIST) G = 0.0;
  double E1 = (X - zMax * (Z - zMax));
  if (fabs(E1) < MINDIST) E1 = 0.0;
  double E2 = (X - 1.0 - zMax * Z);
  if (fabs(E2) < MINDIST) E2 = 0.0;
  double H1 = Y * Y + (Z - zMax) * G;
  if (fabs(H1) < MINDIST2) H1 = 0.0;
  double H2 = Y * Y + Z * G;
  if (fabs(H2) < MINDIST2) H2 = 0.0;
  double R1 = Y * Y + Z * Z;
  if (fabs(R1) < MINDIST2) R1 = 0.0;
  double I1 = modY * X;
  if (fabs(I1) < MINDIST2) I1 = 0.0;
  double I2 = modY * (X - 1.0);
  if (fabs(I2) < MINDIST2) I2 = 0.0;
  double Hypot = sqrt(1.0 + zMax * zMax);
  if (Hypot < MINDIST) {  // superfluous
    fprintf(stdout, "Hypotenuse too small! ... returning ...\n");
    return -1;
  }
  double Zhyp = zMax * (1.0 - X);  // Z on hypotenuse (extended) for given X

  if (DebugISLES) {
    printf("\n\nzMax: %.16lg, X: %.16lg, Y: %.16lg, Z: %.16lg\n", zMax, X, Y,
           Z);
    printf("D11: %.16lg, D21: %.16lg, D12: %.16lg, Hypot: %.16lg\n", D11, D21,
           D12, Hypot);
    printf("modY: %.16lg, G: %.16lg\n", modY, G);
    printf("E1: %.16lg, E2: %.16lg, H1: %.16lg, H2: %.16lg\n", E1, E2, H1, H2);
    printf("S1: %d, SY: %d, R1: %.16lg, I1: %.16lg, I2: %.16lg\n", S1, SY, R1,
           I1, I2);
    printf("H1+G*D12: %.16lg, E1-zMax*D12: %.16lg\n", H1 + G * D12,
           E1 - zMax * D12);
    printf("H2+G*D21: %.16lg, E2-zMax*D21: %.16lg\n", H2 + G * D21,
           E2 - zMax * D21);
    printf("D11*fabs(Z): %.16lg, D21*fabs(Z): %.16lg\n", D11 * fabs(Z),
           D21 * fabs(Z));
    printf("Hypot*D12 - E1: %.16lg\n", Hypot * D12 - E1);
    printf("Hypot*D21 - E2: %.16lg\n\n", Hypot * D21 - E2);
    fprintf(stdout, "MINDIST: %.16lg, MINDIST2: %.16lg, SHIFT: %.16lg\n",
            MINDIST, MINDIST2, SHIFT);
    fflush(stdout);
  }

  // Check for possible numerical difficulties and take care.
  // Presently the idea is to shift the field point slightly to a 'safe'
  // position. Note that a shift in Y does not work because the singularities
  // are associated with D11-s and dxlo-s and their combinations. A Y shift can
  // sometimes alleviate the problem, but it cannot gurantee a permanent1s
  // solution.
  // A better approach is to evaluate analytic expressions truly valid for
  // these critical regions, especially to take care of the >= and <=
  // comparisons. For the == case, both of the previous comparisons are true and
  // it is hard to justify one choice over the other. On the other hand, there
  // are closed form analytic solutions for the == cases, and the problem does
  // not stem from round-off errors as in the > or <, but not quite == cases.
  // Possible singularity at D21 corner where X=1, Z=0
  // Possible singularity at D12 corner where X=0, Z=zMax
  // Check for possible numerical difficuties and take care
  // Note that modY = 0 cannot be a condition of criticality since considering
  // that would mean omitting even the barycenter properties.
  // Also note that shifts that are exactly equal to MINDIST, can give rise to
  // un-ending recursions, leading to Segmentation fault.
  // Special points and combinations

  // Three corners
  if ((fabs(D11) <= MINDIST)) {
    if (DebugISLES) printf("D11 <= MINDIST\n");
    double X1 = X;
    double Z1 = Z;
    if ((X >= 0.0) && (Z >= 0.0)) {
      // point on the element
      if (DebugISLES) printf("Case1=>\n");
      X1 = X + SHIFT * MINDIST;
      Z1 = Z + SHIFT * MINDIST;
    } else if ((X <= 0.0) && (Z >= 0.0)) {
      // field point outside the element
      if (DebugISLES) printf("Case2=>\n");
      X1 = X - SHIFT * MINDIST;
      Z1 = Z + SHIFT * MINDIST;
    } else if ((X >= 0.0) && (Z <= 0.0)) {
      // field point outside the element
      if (DebugISLES) printf("Case3=>\n");
      X1 = X + SHIFT * MINDIST;
      Z1 = Z - SHIFT * MINDIST;
    } else if ((X <= 0.0) && (Z <= 0.0)) {
      // field point outside the element
      if (DebugISLES) printf("Case4=>\n");
      X1 = X - SHIFT * MINDIST;
      Z1 = Z - SHIFT * MINDIST;
    }
    double Pot1;
    Vector3D Flux1;
    ExactTriSurf(zMax, X1, Y, Z1, &Pot1, &Flux1);
    *Potential = Pot1;
    Flux->X = Flux1.X;
    Flux->Y = Flux1.Y;
    Flux->Z = Flux1.Z;
    return 0;
  }
  if ((fabs(D21) <= MINDIST)) {
    if (DebugISLES) printf("D21 <= MINDIST\n");
    double X1 = X;
    double Z1 = Z;
    if ((X >= 1.0) && (Z >= 0.0)) {
      // point outside the element
      if (DebugISLES) printf("Case1=>\n");
      X1 = X + SHIFT * MINDIST;
      Z1 = Z + SHIFT * MINDIST;
    } else if ((X <= 1.0) && (Z >= 0.0)) {
      // field point on the element
      // difficult to decide (chk figure) - chk whether Z is beyond Zhyp
      if (DebugISLES) printf("Case2=>\n");
      X1 = X - SHIFT * MINDIST;
      Z1 = Z + SHIFT * MINDIST;
    } else if ((X >= 1.0) && (Z <= 0.0)) {
      // field point outside the element
      if (DebugISLES) printf("Case3=>\n");
      X1 = X + SHIFT * MINDIST;
      Z1 = Z - SHIFT * MINDIST;
    } else if ((X <= 1.0) && (Z <= 0.0)) {
      // field point outside the element
      if (DebugISLES) printf("Case4=>\n");
      X1 = X - SHIFT * MINDIST;
      Z1 = Z - SHIFT * MINDIST;
    }
    double Pot1;
    Vector3D Flux1;
    ExactTriSurf(zMax, X1, Y, Z1, &Pot1, &Flux1);
    *Potential = Pot1;
    Flux->X = Flux1.X;
    Flux->Y = Flux1.Y;
    Flux->Z = Flux1.Z;
    return 0;
  }
  if ((fabs(D12) <= MINDIST)) {
    if (DebugISLES) printf("D12 <= MINDIST\n");
    double X1 = X;
    double Z1 = Z;
    if ((X >= 0.0) && (Z >= zMax)) {
      // point outside the element
      if (DebugISLES) printf("Case1=>\n");
      X1 = X + SHIFT * MINDIST;
      Z1 = Z + SHIFT * MINDIST;
    } else if ((X <= 0.0) && (Z >= zMax)) {
      // field point outside the element
      if (DebugISLES) printf("Case2=>\n");
      X1 = X - SHIFT * MINDIST;
      Z1 = Z + SHIFT * MINDIST;
    } else if ((X >= 0.0) && (Z <= zMax)) {
      // field point on the element
      // can create problem for small element - chk whether Z is beyond Zhyp
      if (DebugISLES) printf("Case3=>\n");
      X1 = X + SHIFT * MINDIST;
      Z1 = Z - SHIFT * MINDIST;
    } else if ((X <= 0.0) && (Z <= zMax)) {
      // field point outside the element
      if (DebugISLES) printf("Case4=>\n");
      X1 = X - SHIFT * MINDIST;
      Z1 = Z - SHIFT * MINDIST;
    }
    double Pot1;
    Vector3D Flux1;
    ExactTriSurf(zMax, X1, Y, Z1, &Pot1, &Flux1);
    *Potential = Pot1;
    Flux->X = Flux1.X;
    Flux->Y = Flux1.Y;
    Flux->Z = Flux1.Z;
    return 0;
  }
  // Three Y lines at corners
  if ((fabs(X) < MINDIST) && (fabs(Z) < MINDIST)) {
    // Y line at (0,0,0) corner
    if (DebugISLES) printf("Y line at (0,0,0) corner\n");
    double X1 = X;
    double Z1 = Z;
    if ((X >= 0.0) && (Z >= 0.0)) {
      // point on the element
      if (DebugISLES) printf("Case1=>\n");
      X1 = X + SHIFT * MINDIST;
      Z1 = Z + SHIFT * MINDIST;
    } else if ((X <= 0.0) && (Z >= 0.0)) {
      // field point outside the element
      if (DebugISLES) printf("Case2=>\n");
      X1 = X - SHIFT * MINDIST;
      Z1 = Z + SHIFT * MINDIST;
    } else if ((X >= 0.0) && (Z <= 0.0)) {
      // field point outside the element
      if (DebugISLES) printf("Case3=>\n");
      X1 = X + SHIFT * MINDIST;
      Z1 = Z - SHIFT * MINDIST;
    } else if ((X <= 0.0) && (Z <= 0.0)) {
      // field point outside the element
      if (DebugISLES) printf("Case4=>\n");
      X1 = X - SHIFT * MINDIST;
      Z1 = Z - SHIFT * MINDIST;
    }
    double Pot1;
    Vector3D Flux1;
    ExactTriSurf(zMax, X1, Y, Z1, &Pot1, &Flux1);
    *Potential = Pot1;
    Flux->X = Flux1.X;
    Flux->Y = Flux1.Y;
    Flux->Z = Flux1.Z;
    return 0;
  }
  if ((fabs(X - 1.0) < MINDIST) && (fabs(Z) < MINDIST)) {
    // Y line at (1,0) corner
    if (DebugISLES) printf("Y line at (1,0,0) corner\n");
    double X1 = X;
    double Z1 = Z;
    if ((X >= 1.0) && (Z >= 0.0)) {
      // point on the element
      if (DebugISLES) printf("Case1=>\n");
      X1 = X + SHIFT * MINDIST;
      Z1 = Z + SHIFT * MINDIST;
    } else if ((X <= 1.0) && (Z >= 0.0)) {
      // field point outside the element
      if (DebugISLES) printf("Case2=>\n");
      X1 = X - SHIFT * MINDIST;
      Z1 = Z + SHIFT * MINDIST;
    } else if ((X >= 1.0) && (Z <= 0.0)) {
      // field point outside the element
      if (DebugISLES) printf("Case3=>\n");
      X1 = X + SHIFT * MINDIST;
      Z1 = Z - SHIFT * MINDIST;
    } else if ((X <= 1.0) && (Z <= 0.0)) {
      // field point outside the element
      if (DebugISLES) printf("Case4=>\n");
      X1 = X - SHIFT * MINDIST;
      Z1 = Z - SHIFT * MINDIST;
    }
    double Pot1;
    Vector3D Flux1;
    ExactTriSurf(zMax, X1, Y, Z1, &Pot1, &Flux1);
    *Potential = Pot1;
    Flux->X = Flux1.X;
    Flux->Y = Flux1.Y;
    Flux->Z = Flux1.Z;
    return 0;
  }
  if ((fabs(X) < MINDIST) && (fabs(Z - zMax) < MINDIST)) {
    // Y line at (0,zMax)
    if (DebugISLES) printf("Y line at (0,0,zMax) corner\n");
    double X1 = X;
    double Z1 = Z;
    if ((X >= 0.0) && (Z >= zMax)) {
      // point on the element
      if (DebugISLES) printf("Case1=>\n");
      X1 = X + SHIFT * MINDIST;
      Z1 = Z + SHIFT * MINDIST;
    } else if ((X <= 0.0) && (Z >= zMax)) {
      // field point outside the element
      if (DebugISLES) printf("Case2=>\n");
      X1 = X - SHIFT * MINDIST;
      Z1 = Z + SHIFT * MINDIST;
    } else if ((X >= 0.0) && (Z <= zMax)) {
      // field point outside the element
      if (DebugISLES) printf("Case3=>\n");
      X1 = X + SHIFT * MINDIST;
      Z1 = Z - SHIFT * MINDIST;
    } else if ((X <= 0.0) && (Z <= zMax)) {
      // field point outside the element
      if (DebugISLES) printf("Case4=>\n");
      X1 = X - SHIFT * MINDIST;
      Z1 = Z - SHIFT * MINDIST;
    }
    double Pot1;
    Vector3D Flux1;
    ExactTriSurf(zMax, X1, Y, Z1, &Pot1, &Flux1);
    *Potential = Pot1;
    Flux->X = Flux1.X;
    Flux->Y = Flux1.Y;
    Flux->Z = Flux1.Z;
    return 0;
  }
  // Three edges outside the extent of the real element.
  // Plus two edges delineating the virtual right-triangle that complements
  // the real one to create a rectangle.
  // Special X=1 condition which turns R1+-I2 terms in dire straits
  // Similar problem occurs for X=0, where R1+-I1 is NaN but this has been taken
  // care of earlier.
  // But this particular problem has a remedy - check note on complex tanh-1
  // below. There is the problem of dealing with gsl_complex_log_real(1.0)
  // though!
  if (fabs(X) < MINDIST) {
    // edge along X - axis
    if (DebugISLES) printf("edge along X-axis\n");
    double X1 = X, X2 = X;
    if (X >= 0.0) {
      // field point on +ve side of YZ plane
      if (DebugISLES) printf("Case1=>\n");
      X1 = X + SHIFT * MINDIST;
      X2 = X - SHIFT * MINDIST;
    } else if (X <= 0.0) {
      // field point on -ve side of YZ plane
      if (DebugISLES) printf("Case2=>\n");
      X1 = X - SHIFT * MINDIST;
      X2 = X + SHIFT * MINDIST;
    }
    double Pot1, Pot2;
    Vector3D Flux1, Flux2;
    ExactTriSurf(zMax, X1, Y, Z, &Pot1, &Flux1);
    ExactTriSurf(zMax, X2, Y, Z, &Pot2, &Flux2);
    *Potential = 0.5 * (Pot1 + Pot2);
    Flux->X = 0.5 * (Flux1.X + Flux2.X);
    Flux->Y = 0.5 * (Flux1.Y + Flux2.Y);
    Flux->Z = 0.5 * (Flux1.Z + Flux2.Z);
    return 0;
  }
  /* do not erase these blocks as yet
  if (fabs(Z) < MINDIST) {
    // edge along Z - axis
    if (DebugISLES) printf("edge along Z-axis\n");
    double Z1 = Z;
    if (Z >= 0.0) {
      // field point on +ve side of XY plane
      if (DebugISLES) printf("Case1=>\n");
      Z1 = Z + SHIFT*MINDIST;
    } else if (Z <= 0.0) {
      // field point on -ve side of XY plane
      if (DebugISLES) printf("Case2=>\n");
      Z1 = Z - SHIFT*MINDIST;
    }
    ExactTriSurf(zMax, X, Y, Z1, &Pot1, &Flux1);
    *Potential = Pot1;
    Flux->X = Flux1.X;
    Flux->Y = Flux1.Y;
    Flux->Z = Flux1.Z;
    return 0;
  }
  */
  if (fabs(X - 1.0) < MINDIST) {
    // edge along X=1.0
    if (DebugISLES) printf("edge along X = 1.\n");
    double X1 = X, X2 = X;
    if (X <= 1.0) {
      // field point on +ve side of YZ plane
      if (DebugISLES) printf("Case1=>\n");
      X1 = X - SHIFT * MINDIST;
      X2 = X + SHIFT * MINDIST;
    } else if (X >= 1.0) {
      // field point on -ve side of YZ plane
      if (DebugISLES) printf("Case2=>\n");
      X1 = X + SHIFT * MINDIST;
      X2 = X + SHIFT * MINDIST;
    }
    double Pot1, Pot2;
    Vector3D Flux1, Flux2;
    ExactTriSurf(zMax, X1, Y, Z, &Pot1, &Flux1);
    ExactTriSurf(zMax, X2, Y, Z, &Pot2, &Flux2);
    *Potential = 0.5 * (Pot1 + Pot2);
    Flux->X = 0.5 * (Flux1.X + Flux2.X);
    Flux->Y = 0.5 * (Flux1.Y + Flux2.Y);
    Flux->Z = 0.5 * (Flux1.Z + Flux2.Z);
    return 0;
  }
  if (fabs(Z - zMax) < MINDIST) {
    // edge along Z=zMax
    if (DebugISLES) printf("edge along Z = zMax\n");
    double Z1 = Z, Z2 = Z;
    if (Z >= zMax) {
      // field point on +ve side of XY plane
      if (DebugISLES) printf("Case1=>\n");
      Z1 = Z + SHIFT * MINDIST;
      Z2 = Z - SHIFT * MINDIST;
    } else if (Z <= zMax) {
      // field point on -ve side of XY plane
      if (DebugISLES) printf("Case2=>\n");
      Z1 = Z - SHIFT * MINDIST;
      Z2 = Z + SHIFT * MINDIST;
    }
    double Pot1, Pot2;
    Vector3D Flux1, Flux2;
    ExactTriSurf(zMax, X, Y, Z1, &Pot1, &Flux1);
    ExactTriSurf(zMax, X, Y, Z2, &Pot2, &Flux2);
    *Potential = 0.5 * (Pot1 + Pot2);
    Flux->X = 0.5 * (Flux1.X + Flux2.X);
    Flux->Y = 0.5 * (Flux1.Y + Flux2.Y);
    Flux->Z = 0.5 * (Flux1.Z + Flux2.Z);
    return 0;
  }
  if (fabs(Z - Zhyp) < MINDIST) {
    // edge along the hypotenuse
    if (DebugISLES) printf("edge along Hypotenuse\n");
    double X1 = X, Z1 = Z;
    if (Z <= Zhyp) {
      // towards element origin
      if (DebugISLES) printf("Case1=>\n");
      X1 = X - SHIFT * MINDIST;
      Z1 = Z - SHIFT * MINDIST;
    } else if (Z >= Zhyp) {
      // going further away from the element origin
      if (DebugISLES) printf("Case2=>\n");
      X1 = X + SHIFT * MINDIST;
      Z1 = Z + SHIFT * MINDIST;
    }
    double Pot1;
    Vector3D Flux1;
    ExactTriSurf(zMax, X1, Y, Z1, &Pot1, &Flux1);
    *Potential = Pot1;
    Flux->X = Flux1.X;
    Flux->Y = Flux1.Y;
    Flux->Z = Flux1.Z;
    return 0;
  }

  // Related to logarithmic terms (LTerm1 and LTerm2)
  if ((fabs(G) <= MINDIST) && (modY <= MINDIST)) {
    ApproxFlag = 10;
    ++FailureCntr;
    --ExactCntr;
    fprintf(
        fIsles,
        "(fabs(G) <= MINDIST) && (modY <= MINDIST) ... approximating: %d.\n",
        ApproxFlag);
    return (ApproxTriSurf(zMax, X, Y, Z, XNSegApprox, ZNSegApprox, Potential,
                          Flux));
  }
  /* do not erase these blocks as yet
  if ((fabs(X) <= MINDIST) && (modY <= MINDIST)) {
    // denominator zero for LP1 and LM1
    ApproxFlag = 11; ++FailureCntr; --ExactCntr;
    fprintf(fIsles, "denominator zero for LP1 and LM1 ... approximating: %d.\n",
            ApproxFlag);
    return (ApproxTriSurf(zMax, X, Y, Z, XNSegApprox, ZNSegApprox, Potential,
                          Flux));
  }
  if ((fabs(H1 + D12 * G) <= MINDIST2) && (modY <= MINDIST)) {
    // numerator zero for LP1 and LM1
    ApproxFlag = 12; ++FailureCntr; --ExactCntr;
    fprintf(fIsles, "numerator zero for LP1 and LM1 ... approximating: %d.\n",
            ApproxFlag);
    return(ApproxTriSurf(zMax, X, Y, Z, XNSegApprox, ZNSegApprox, Potential,
                         Flux));
  }
  if ((fabs(1.0 - X) <= MINDIST) && (modY <= MINDIST) )	{
    // denominator zero for LP2 and LM2
    ApproxFlag = 13; ++FailureCntr; --ExactCntr;
    fprintf(fIsles, "denominator zero for LP2 and LM2 ... approximating: %d.\n",
            ApproxFlag);
    return(ApproxTriSurf(zMax, X, Y, Z, XNSegApprox, ZNSegApprox, Potential,
                         Flux));
  }
  if ((fabs(H2 + D21 * G) <= MINDIST2) && (modY <= MINDIST)) {
    // numerator zero for LP2 and LM2
    ApproxFlag = 14; ++FailureCntr; --ExactCntr;
    fprintf(fIsles, "numerator zero for LP2 and LM2 ... approximating: %d.\n",
            ApproxFlag);
    return (ApproxTriSurf(zMax, X, Y, Z, XNSegApprox, ZNSegApprox, Potential,
                          Flux));
  }
  */

  // Related to complex inverse tan hyperbolic terms - TanhTerm1 and TanhTerm2
  /* do not erase these blocks as yet
  if (D11 * fabs(Z) <= MINDIST2) {
    ApproxFlag = 15; ++FailureCntr; --ExactCntr;
    fprintf(fIsles, "D11*fabs(Z) zero for TanhTerms ... approximating: %d.\n",
            ApproxFlag);
    return(ApproxTriSurf(zMax, X, Y, Z, XNSegApprox, ZNSegApprox, Potential,
  Flux));
  }
  if (D21 * fabs(Z) <= MINDIST2) {
    ApproxFlag = 16; ++FailureCntr; --ExactCntr;
    fprintf(fIsles, "D21*fabs(Z) zero for TanhTerms ... approximating: %d.\n",
            ApproxFlag);
    return(ApproxTriSurf(zMax, X, Y, Z, XNSegApprox, ZNSegApprox, Potential,
  Flux));
  }
  */

  // Exact computations begin here, at last!
  // DTerm1 and DTerm2
  double DblTmp1 = (Hypot * D12 - E1) / (Hypot * D21 - E2);
  if (DebugISLES) {
    printf("DblTmp1: %.16lg\n", DblTmp1);
    fflush(stdout);
  }
  if (DblTmp1 < MINDIST2) DblTmp1 = MINDIST2;
  double DTerm1 = log(DblTmp1);
  if (isnan(DTerm1) || isinf(DTerm1)) {
    ApproxFlag = 18;
    ++FailureCntr;
    --ExactCntr;
    fprintf(fIsles, "DTerm1 nan or inf ... approximating: %d.\n", ApproxFlag);
    return (ApproxTriSurf(zMax, X, Y, Z, XNSegApprox, ZNSegApprox, Potential,
                          Flux));
  }
  DTerm1 /= Hypot;

  double DblTmp2 = (D11 - X) / (D21 - X + 1.0);
  if (DebugISLES) {
    printf("DblTmp2: %.16lg\n", DblTmp2);
    fflush(stdout);
  }
  if (DblTmp2 < MINDIST2) DblTmp2 = MINDIST2;
  double DTerm2 = log(DblTmp2);
  if (isnan(DTerm2) || isinf(DTerm2)) {
    ApproxFlag = 18;
    ++FailureCntr;
    --ExactCntr;
    fprintf(fIsles, "DTerm2 nan or inf ... approximating: %d.\n", ApproxFlag);
    return (ApproxTriSurf(zMax, X, Y, Z, XNSegApprox, ZNSegApprox, Potential,
                          Flux));
  }

  if (DebugISLES) {
    printf("DTerm1: %.16lg, DTerm2: %.16lg\n", DTerm1, DTerm2);
    fflush(stdout);
  }

  // Logarithmic terms
  double deltaLog = 0.;
  double deltaArg = 0.;
  if (modY > MINDIST) {
    const double s1 = 1. / (X * X + modY * modY);
    double Re1 = ((-X) * (H1 + G * D12) + modY * modY * (E1 - zMax * D12)) * s1;
    double Im1 = ((-X) * modY * (E1 - zMax * D12) - (H1 + G * D12) * modY) * s1;
    const double s2 = 1. / ((1. - X) * (1. - X) + modY * modY);
    double Re2 =
        ((1. - X) * (H2 + G * D21) + modY * modY * (E2 - zMax * D21)) * s2;
    double Im2 =
        ((1. - X) * modY * (E2 - zMax * D21) - (H2 + G * D21) * modY) * s2;
    deltaLog = 0.5 * log((Re1 * Re1 + Im1 * Im1) / (Re2 * Re2 + Im2 * Im2));
    deltaArg = atan2(Im1, Re1) - atan2(Im2, Re2);
  } else {
    double f1 = fabs((H1 + G * D12) / (-X));
    if (f1 < MINDIST) f1 = MINDIST;
    double f2 = fabs((H2 + G * D21) / (1. - X));
    if (f2 < MINDIST) f2 = MINDIST;
    deltaLog = log(f1 / f2);
  }
  const double c1 = 2. / (G * G + zMax * zMax * modY * modY);
  double LTerm1 = c1 * (G * deltaLog - zMax * modY * deltaArg);
  double LTerm2 = c1 * (G * deltaArg + zMax * modY * deltaLog);
  if (DebugISLES) {
    printf("LTerm1: %.16lg, LTerm2: %.16lg\n", LTerm1, LTerm2);
    fflush(stdout);
  }

  // Computations for estimating TanhTerm1, TanhTerm2
  // Possible singularities - D11, D21 corners and Z=0 line / surface
  // Corners have been taken care of, as well as Z=0 (latter, through S1)
  // Incorrectly implemented - CHECK!!! CORRECTED!
  // If the argument of atanh is real and greater than 1.0, it is
  // perfectly computable. The value returned is a complex number,
  // however. This computation has to be carried out by invoking
  // gsl_complex_arctanh_real(double z).
  // The branch cut needs to be enforced only if the
  // argument is real and is equal to 1.0. Then the returned value
  // is undefined. But this happens only if Y=0 besides X=1.0!
  // Similar implementation issues are also there for X=0.
  // Coincides with X=1 since I2 turns out to be zero in such
  // cases. CmTmp3.dat[0] can be >= 1.0 in a variety of situations
  // since in the denominator D21, sqrt( (X-1)^2 + Y^2 + Z^2 ) is
  // sqrt(Y^2+Z^2) for X=1. This is multiplied by |Z| while on the
  // numerator we have Y^2+Z^2.

  double TanhTerm1 = 0.;
  double TanhTerm2 = 0.;
  if (abs(S1) > 0) {
    if (modY < MINDIST) {
      const double f1 = R1 / (D11 * fabs(Z));
      const double f2 = R1 / (D21 * fabs(Z));
      TanhTerm1 = log1p(f1) - log1p(-f1) - log1p(f2) + log1p(-f2);
    } else {
      TanhTerm1 = 1.;
      if (fabs(R1 - D11 * fabs(Z)) > MINDIST2 || true) {
        const double rp = D11 * fabs(Z) + R1;
        const double rm = D11 * fabs(Z) - R1;
        TanhTerm1 *= (I1 * I1 + rp * rp) / (I1 * I1 + rm * rm);
      }
      if (fabs(R1 - D21 * fabs(Z)) > MINDIST2 || true) {
        const double rp = D21 * fabs(Z) + R1;
        const double rm = D21 * fabs(Z) - R1;
        TanhTerm1 *= (I2 * I2 + rm * rm) / (I2 * I2 + rp * rp);
      }
      TanhTerm1 = 0.5 * log(TanhTerm1);
    }

    if (fabs(I1) > MINDIST2) {
      double a = D11 * D11 * Z * Z;
      double b = R1 * R1 + I1 * I1;
      double tmp = atan(2 * I1 * D11 * fabs(Z) / (a - b));
      if (b > a) {
        if (X > 0.) {
          tmp += ST_PI;
        } else {
          tmp -= ST_PI;
        }
      }
      TanhTerm2 += tmp;
    }
    if (fabs(I2) > MINDIST2) {
      double a = D21 * D21 * Z * Z;
      double b = R1 * R1 + I2 * I2;
      double tmp = atan(2 * I2 * D21 * fabs(Z) / (a - b));
      if (b > a) {
        if (X > 1.) {
          tmp += ST_PI;
        } else {
          tmp -= ST_PI;
        }
      }
      TanhTerm2 -= tmp;
    }
  }

  if (DebugISLES) {
    printf("TanhTerm1: %.16lg, TanhTerm2: %.16lg\n", TanhTerm1, TanhTerm2);
    fflush(stdout);
  }

  double Pot = (zMax * Y * Y - X * G) * LTerm1 -
               (zMax * X + G) * modY * LTerm2 + (S1 * X * TanhTerm1) +
               S1 * modY * TanhTerm2 + 2. * (G * DTerm1 - Z * DTerm2);
  Pot *= 0.5;
  double Fx =
      G * LTerm1 + zMax * modY * LTerm2 - S1 * TanhTerm1 - 2.0 * zMax * DTerm1;
  Fx *= 0.5;
  double Fy = -zMax * Y * LTerm1 + G * SY * LTerm2 - S1 * SY * TanhTerm2;
  Fy *= 0.5;

  double Fz = DTerm2 - DTerm1;
  if (DebugISLES) {
    printf("Pot: %.16lg, Fx: %.16lg, Fy: %.16lg, Fz: %.16lg\n", Pot, Fx, Fy,
           Fz);
    fflush(stdout);
  }

  // Final adjustments
  // Constants (?) of integration - Carry out only one of the options
  // Depends criticially on > or >=; similarly < or <=
  // Needs further investigation
  // As far as Fy is concerned, conditions 1 & 3, and 2 & 4 can be combined.
  // So, instead of the triangular area, it seems, that the rectangular bound
  // is more important. In such an event, Zhyp will be redundant.
  if ((X >= 0.0) && (X <= 1.0)) {
    // Possibility of point within element bounds
    int ConstAdd = 0;
    if ((Z >= 0.0) && (Z <= Zhyp)) {
      // within the element
      ConstAdd = 1;
      Pot -= modY * ST_PI;
      if (fabs(Y) < MINDIST)  // on the element surface
        Fy = ST_PI;
      else if (Y > 0.0)  // above or below the element surface
        Fy += ST_PI;
      else if (Y < 0.0)
        Fy -= ST_PI;
    } else if ((Z <= 0.0) && (Z >= -Zhyp)) {
      // within -ve element
      ConstAdd = 2;
      Pot += modY * ST_PI;
      if (fabs(Y) < MINDIST)
        Fy = 0.0;
      else if (Y > 0.0)
        Fy -= ST_PI;
      else if (Y < 0.0)
        Fy += ST_PI;
    } else if (((Z > Zhyp) && (Z <= zMax))) {
      // within +rect bounds
      // merge with +ve shadow?
      ConstAdd = 3;
      Pot -= modY * ST_PI;
      if (fabs(Y) < MINDIST)
        Fy = 0.0;
      else if (Y > 0.0)
        Fy += ST_PI;
      else if (Y < 0.0)
        Fy -= ST_PI;
    } else if ((Z < -Zhyp) && (Z >= -zMax)) {
      // within -rect bounds
      // merge with -ve shadow?
      ConstAdd = 4;
      Pot += modY * ST_PI;
      if (fabs(Y) < MINDIST)
        Fy = 0.0;
      else if (Y > 0.0)
        Fy -= ST_PI;
      else if (Y < 0.0)
        Fy += ST_PI;
    } else if (Z > zMax) {
      // +ve shadow of the triangle - WHY?
      ConstAdd = 5;
      Pot -= modY * ST_PI;
      if (fabs(Y) < MINDIST)
        Fy = 0.0;
      else if (Y > 0.0)
        Fy += ST_PI;
      else if (Y < 0.0)
        Fy -= ST_PI;
    } else if (Z < -zMax) {
      // -ve shadow of the triangle - WHY?
      ConstAdd = 6;
      Pot += modY * ST_PI;
      if (fabs(Y) < MINDIST)
        Fy = 0.0;
      else if (Y > 0.0)
        Fy -= ST_PI;
      else if (Y < 0.0)
        Fy += ST_PI;
    }

    /*
    if ((fabs(X - 1.0) < MINDIST) && (fabs(Z-zMax) < MINDIST)) {
      // (1,zMax) corner
      if (Y > 0.0) Fy += 0.5*ST_PI;
      if (Y < 0.0) Fy -= 0.5*ST_PI;
    } else if ((fabs(X - 1.0) < MINDIST) && (fabs(Z + zMax) < MINDIST)) {
      // (1,-zMax) corner
      if (Y > 0.0) Fy -= 0.5*ST_PI;
      if (Y < 0.0) Fy += 0.5*ST_PI;
    } else if ((fabs(X) < MINDIST) && (Z < 0.0)) {
      // X = 0, Z < 0 - off ele
      // along -ve Z axis
      if (Y > 0.0) Fy -= 0.5*ST_PI;
      if (Y < 0.0) Fy += 0.5*ST_PI;
    } else if ((fabs(X) < MINDIST) && (Z > 0.0)) {
      // X = 0, Z > 0 - on & off ele
      // along +ve Z axis
      if (Y > 0.0) Fy += 0.5*ST_PI;
      if (Y < 0.0) Fy -= 0.5*ST_PI;
    } else if (Z > 0.0) {
      // within rectangular bounds - changed for INPC
      if (Y > 0.0) Fy += ST_PI;
      if (Y < 0.0) Fy -= ST_PI;
    } else if (Z < 0.0) {
      // within -ve rectangular bounds - changed for INPC
      if (Y > 0.0) Fy -= ST_PI;
      if (Y < 0.0) Fy += ST_PI;
    }
  */
    // two other conditions, one for Z > zMax and another for Z < -zMax expected

    if (DebugISLES) {
      printf("After constant addition %d\n", ConstAdd);
      printf("Pot: %.16lg, Fx: %.16lg, Fy: %.16lg, Fz: %.16lg\n", Pot, Fx, Fy,
             Fz);
      fflush(stdout);
    }
  }  // point within element bounds

  // Error situations handled before returning the potential value
  if ((Pot < 0) || isnan(Pot) || isinf(Pot)) {
    fprintf(fIsles, "\n---Approximation in ExactTriSurf---\n");
    fprintf(fIsles, "Problem with potential ... approximating!\n");
    fprintf(fIsles, "zMax: %.16lg, X: %.16lg, Y: %.16lg, Z: %.16lg\n", zMax, X,
            Y, Z);
    fprintf(fIsles, "D11: %.16lg, D21: %.16lg, D12: %.16lg, Hypot: %.16lg\n",
            D11, D21, D12, Hypot);
    fprintf(fIsles, "modY: %.16lg, G: %.16lg\n", modY, G);
    fprintf(fIsles, "E1: %.16lg, E2: %.16lg, H1: %.16lg, H2: %.16lg\n", E1, E2,
            H1, H2);
    fprintf(fIsles, "S1: %d, R1: %.16lg, I1: %.16lg, I2: %.16lg\n", S1, R1, I1,
            I2);
    fprintf(fIsles, "Pot: %.16lg\n", Pot);
    ApproxFlag = 23;
    ++FailureCntr;
    --ExactCntr;
    return (ApproxTriSurf(zMax, X, Y, Z, XNSegApprox, ZNSegApprox, Potential,
                          Flux));
  }

  if (isnan(Fx) || isinf(Fx)) {
    fprintf(fIsles, "\n---Approximation in ExactTriSurf---\n");
    fprintf(fIsles, "Problem with Fx ... approximating!\n");
    fprintf(fIsles, "zMax: %.16lg, X: %.16lg, Y: %.16lg, Z: %.16lg\n", zMax, X,
            Y, Z);
    fprintf(fIsles, "D11: %.16lg, D21: %.16lg, D12: %.16lg, Hypot: %.16lg\n",
            D11, D21, D12, Hypot);
    fprintf(fIsles, "modY: %.16lg, G: %.16lg\n", modY, G);
    fprintf(fIsles, "E1: %.16lg, E2: %.16lg, H1: %.16lg, H2: %.16lg\n", E1, E2,
            H1, H2);
    fprintf(fIsles, "S1: %d, R1: %.16lg, I1: %.16lg, I2: %.16lg\n", S1, R1, I1,
            I2);
    fprintf(fIsles, "Fx: %.16lg\n", Fx);
    ApproxFlag = 24;
    ++FailureCntr;
    --ExactCntr;
    return (ApproxTriSurf(zMax, X, Y, Z, XNSegApprox, ZNSegApprox, Potential,
                          Flux));
  }

  if (isnan(Fy) || isinf(Fy)) {
    fprintf(fIsles, "\n---Approximation in ExactTriSurf---\n");
    fprintf(fIsles, "Problem with Fy ... approximating!\n");
    fprintf(fIsles, "zMax: %.16lg, X: %.16lg, Y: %.16lg, Z: %.16lg\n", zMax, X,
            Y, Z);
    fprintf(fIsles, "D11: %.16lg, D21: %.16lg, D12: %.16lg, Hypot: %.16lg\n",
            D11, D21, D12, Hypot);
    fprintf(fIsles, "modY: %.16lg, G: %.16lg\n", modY, G);
    fprintf(fIsles, "E1: %.16lg, E2: %.16lg, H1: %.16lg, H2: %.16lg\n", E1, E2,
            H1, H2);
    fprintf(fIsles, "S1: %d, R1: %.16lg, I1: %.16lg, I2: %.16lg\n", S1, R1, I1,
            I2);
    fprintf(fIsles, "Fy: %.16lg\n", Fy);
    ApproxFlag = 25;
    ++FailureCntr;
    --ExactCntr;
    return (ApproxTriSurf(zMax, X, Y, Z, XNSegApprox, ZNSegApprox, Potential,
                          Flux));
  }

  if (isnan(Fz) || isinf(Fz)) {
    fprintf(fIsles, "\n---Approximation in ExactTriSurf---\n");
    fprintf(fIsles, "Problem with Fz ... approximating!\n");
    fprintf(fIsles, "zMax: %.16lg, X: %.16lg, Y: %.16lg, Z: %.16lg\n", zMax, X,
            Y, Z);
    fprintf(fIsles, "D11: %.16lg, D21: %.16lg, D12: %.16lg, Hypot: %.16lg\n",
            D11, D21, D12, Hypot);
    fprintf(fIsles, "modY: %.16lg\n", modY);
    fprintf(fIsles, "E1: %.16lg, E2: %.16lg\n", E1, E2);
    fprintf(fIsles, "Fz: %.16lg\n", Fz);
    ApproxFlag = 26;
    ++FailureCntr;
    --ExactCntr;
    return (ApproxTriSurf(zMax, X, Y, Z, XNSegApprox, ZNSegApprox, Potential,
                          Flux));
  }

  // Return the computed value of potential and flux
  *Potential = Pot;
  Flux->X = Fx;
  Flux->Y = Fy;
  Flux->Z = Fz;

  if (DebugISLES) printf("Going out of ExactTriSurf ...\n\n");
  return 0;
}  // ExactTriSurf ends

// Following function assumes that the rt triangle is (0,0), (1,0) and (0,zMax)
// Barycenter criticalities:
// It should be noted here that reasonably correct results for potential
// is obtained using a low discretization (as low as 100 by 100). The same is
// not at all correct for forces, where correct results starts appearing for
// number of segments as high as 3000 by 3000!
// Note that the centroid and the distance of a field point from it has been
// taken care of properly for all the thre shapes (triangle, rectangle and
// trapezoid) have been properly taken care of in this implementation.
int ApproxTriSurf(double zMax, double X, double Y, double Z, int nbxseg,
                  int nbzseg, double *Potential, Vector3D *Flux) {
  if (DebugISLES) printf("In ApproxTriSurf ...\n");
  ++ApproxCntr;

  if (DebugISLES) {
    printf("zMax: %lg, X: %lg, Y: %lg, Z: %lg\n", zMax, X, Y, Z);
    printf("nbxseg: %d, nbzseg: %d\n", nbxseg, nbzseg);
  }

  double dx = (1.0 - 0.0) / nbxseg;
  double dz = (zMax - 0.0) / nbzseg;
  double diag = sqrt(dx * dx + dz * dz);
  if (DebugISLES) printf("dx: %lg, dz: %lg, diag: %lg\n", dx, dz, diag);
  if ((dx < MINDIST) || (dz < MINDIST)) {
    printf("sub-element size too small in ApproxTriSurf.\n");
    return -1;
  }

  double grad = (zMax - 0.0) / (1.0 - 0.0);
  if (DebugISLES) printf("grad: %lg\n", grad);

  double Pot = 0., XFlux = 0., YFlux = 0., ZFlux = 0.;
  for (int i = 1; i <= nbxseg; ++i) {
    double xbgn = (i - 1) * dx;
    double zlimit_xbgn = zMax - grad * xbgn;
    double xend = i * dx;
    double zlimit_xend = zMax - grad * xend;
    if (DebugISLES)
      printf("i: %d, xbgn: %lg, zlimit_xbgn: %lg, xend: %lg, zlimit_xend:%lg\n",
             i, xbgn, zlimit_xbgn, xend, zlimit_xend);

    for (int k = 1; k <= nbzseg; ++k) {
      double zbgn = (k - 1) * dz;
      double zend = k * dz;
      if (DebugISLES) printf("k: %d, zbgn: %lg, zend: %lg\n", k, zbgn, zend);
      int type_subele = 0;
      double area = 0.;
      double xi = 0., zk = 0.;
      if (zbgn >= zlimit_xbgn) {
        // completely outside the element - no effect of sub-element
        type_subele = 0;
        area = 0.0;
      } else if (zend <= zlimit_xend) {
        // completely within the element - rectangular sub-element
        type_subele = 1;
        xi = xbgn + 0.5 * dx;
        zk = zbgn + 0.5 * dz;
        area = dx * dz;
      } else if ((zbgn <= zlimit_xend) && (zend >= zlimit_xend)) {
        // partially inside the triangle
        type_subele = 2;
        // refer to the trapezoid figure
        double a = (zlimit_xend - zbgn);
        double b = (zlimit_xbgn - zbgn);
        double h = dx;

        if (fabs(a) <= MINDIST) {
          // centroid of a triangle
          type_subele = 3;
          xi = xbgn + (h / 3.0);
          zk = zbgn + (b / 3.0);
          area = 0.5 * b * h;
        } else {
          // centroid of the trapezoid
          xi = xbgn + (h * (2. * a + b)) / (3. * (a + b));
          zk = zbgn + (a * a + a * b + b * b) / (3. * (a + b));
          area = 0.5 * h * (a + b);
        }
      } else {
        // takes care of round-off issues
        type_subele = 4;
        area = 0.0;
      }
      if (DebugISLES) printf("type_subele: %d, area: %lg\n", type_subele, area);

      if (area <= MINDIST2) continue;
      double dist = sqrt((X - xi) * (X - xi) + Y * Y + (Z - zk) * (Z - zk));
      if (DebugISLES) printf("dist: %lg\n", dist);
      if (dist >= diag) {
        Pot += area / dist;
        double f = area / (dist * dist * dist);
        XFlux += f * (X - xi);
        YFlux += f * Y;
        ZFlux += f * (Z - zk);
      } else {
        Pot += area / diag;  // replace by expression of self-influence
        if (DebugISLES) printf("Special Pot: %lg\n", area / diag);
        double f = area / (diag * diag * diag);
        XFlux += f * (X - xi);
        YFlux += f * Y;
        ZFlux += f * (Z - zk);
        if (DebugISLES) {
          printf("Special XFlux: %lg, YFlux: %lg, ZFlux: %lg\n", f * (X - xi),
                 f * Y, f * (Z - zk));
        }
      }
    }  // nbzseg
  }    // nbxseg

  *Potential = Pot;
  Flux->X = XFlux;
  Flux->Y = YFlux;
  Flux->Z = ZFlux;

  return 0;
}  // ApproxTriSurf ends

// The wire segment is assumed to be along Z.
// The first three functions, ExactCentroidalP_W, ExactAxialP_W and
// ExactAxialFZ_W will have to remain independent in order to retain
// backward compatibility.
// Combined `PF' functions for ApproxWire, ImprovedWire and ExactThinWire
// have been developed to save some computation and time.

// Potential at the centroidal point of the wire
double ExactCentroidalP_W(double rW, double lW) {
  // Self-Influence: wire5.m
  if (DebugISLES) printf("In ExactCentroidalP_W ...\n");
  const double h = 0.5 * lW;
  double dtmp1 = hypot(rW, h);
  dtmp1 = log((dtmp1 + h) / (dtmp1 - h));
  return 2.0 * ST_PI * dtmp1 * rW;
}  // ExactCentroidalP ends

// Potential along the axis of the wire
double ExactAxialP_W(double rW, double lW, double Z) {
  if (DebugISLES) printf("In ExactAxialP_W ...\n");
  // Expressions from PF00Z.m
  const double h = 0.5 * lW;
  const double r2 = rW * rW;
  const double a = h + Z;
  const double b = h - Z;
  double Pot = log((a + sqrt(a * a + r2)) * (b + sqrt(b * b + r2)) / r2);
  return 2. * ST_PI * rW * Pot;
}  // ExactAxialP ends

// Axial field along the axis of the wire
// Repeated execution of similar expression - can be optimized
double ExactAxialFZ_W(double rW, double lW, double Z) {
  if (DebugISLES) printf("In ExactAxialFZ_W ...\n");
  double h = 0.5 * lW;
  double Fz = 2.0 * ST_PI *
              (sqrt(h * h + 2.0 * Z * h + Z * Z + rW * rW) -
               sqrt(h * h - 2.0 * Z * h + Z * Z + rW * rW)) /
              sqrt(h * h - 2.0 * Z * h + Z * Z + rW * rW) /
              sqrt(h * h + 2.0 * Z * h + Z * Z + rW * rW);
  return (rW * Fz);
}  // ExactAxialFZ ends

double ApproxP_W(double rW, double lW, double X, double Y, double Z, int zseg) {
  // Approximate potential due to a wire segment
  if (DebugISLES) printf("In ApproxP_W ...\n");
  ++ApproxCntr;

  double dz = lW / zseg;
  double area = 2.0 * ST_PI * rW * dz;

  double Pot = 0.0;
  double z0 = -0.5 * lW + 0.5 * dz;
  for (int k = 1; k <= zseg; ++k) {
    double zk = z0 + (k - 1) * dz;
    double dist = sqrt(X * X + Y * Y + (Z - zk) * (Z - zk));
    if (fabs(dist) >= MINDIST) {
      Pot += area / dist;
    }
  }  // zseg

  return (Pot);
}  // ApproxP_W ends

double ApproxFX_W(double rW, double lW, double X, double Y, double Z,
                  int zseg) {
  if (DebugISLES) printf("In ApproxFX_W ...\n");
  ++ApproxCntr;

  double dz = lW / zseg;
  double area = 2.0 * ST_PI * rW * dz;

  double Fx = 0.0;
  double z0 = -0.5 * lW + 0.5 * dz;
  for (int k = 1; k <= zseg; ++k) {
    double zk = z0 + (k - 1) * dz;
    double dist = sqrt(X * X + Y * Y + (Z - zk) * (Z - zk));
    double dist3 = pow(dist, 3.0);
    if (fabs(dist) >= MINDIST) {
      Fx += (area * X / dist3);
    }
  }  // zseg

  return (Fx);
}  // ApproxFX_W ends

double ApproxFY_W(double rW, double lW, double X, double Y, double Z,
                  int zseg) {
  if (DebugISLES) printf("In ApproxFY_W ...\n");
  ++ApproxCntr;

  double dz = lW / zseg;
  double area = 2.0 * ST_PI * rW * dz;

  double Fy = 0.0;
  double z0 = -0.5 * lW + 0.5 * dz;
  for (int k = 1; k <= zseg; ++k) {
    double zk = z0 + (k - 1) * dz;
    double dist = sqrt(X * X + Y * Y + (Z - zk) * (Z - zk));
    double dist3 = pow(dist, 3.0);
    if (fabs(dist) >= MINDIST) {
      Fy += (area * X / dist3);
    }
  }  // zseg

  return (Fy);
}  // ApproxFY_W ends

double ApproxFZ_W(double rW, double lW, double X, double Y, double Z,
                  int zseg) {
  if (DebugISLES) printf("In ApproxFZ_W ...\n");
  ++ApproxCntr;

  double dz = lW / zseg;
  double area = 2.0 * ST_PI * rW * dz;

  double Fz = 0.0;
  double z0 = -0.5 * lW + 0.5 * dz;
  for (int k = 1; k <= zseg; ++k) {
    double zk = z0 + (k - 1) * dz;
    double dist = sqrt(X * X + Y * Y + (Z - zk) * (Z - zk));
    double dist3 = pow(dist, 3.0);
    if (fabs(dist) >= MINDIST) {
      Fz += (area * X / dist3);
    }
  }  // zseg

  return (Fz);
}  // ApproxFZ_W ends

int ApproxWire(double rW, double lW, double X, double Y, double Z, int zseg,
               double *potential, Vector3D *Flux) {
  if (DebugISLES) printf("In ApproxWire ...\n");

  ++ApproxCntr;

  double dz = lW / zseg;
  double area = 2.0 * ST_PI * rW * dz;

  double Pot = 0., Fx = 0., Fy = 0., Fz = 0.;
  double z0 = -0.5 * lW + 0.5 * dz;
  for (int k = 1; k <= zseg; ++k) {
    double zk = z0 + (k - 1) * dz;
    double dist = sqrt(X * X + Y * Y + (Z - zk) * (Z - zk));
    double dist3 = pow(dist, 3.0);
    if (fabs(dist) >= MINDIST) {
      Pot += area / dist;
      Fx += (area * X / dist3);
      Fy += (area * Y / dist3);
      Fz += (area * Z / dist3);
    }
  }  // zseg

  *potential = Pot;
  Flux->X = Fx;
  Flux->Y = Fy;
  Flux->Z = Fz;

  return 0;
}

double ImprovedP_W(double rW, double lW, double X, double Y, double Z) {
  // Improved potential at far away points:
  // wire2.m applicable for thin wires
  if (DebugISLES) printf("In ImprovedP_W ...\n");

  double dz = 0.5 * lW;  // half length of the wire segment
  double dtmp1 = (X * X + Y * Y + (Z - dz) * (Z - dz));
  dtmp1 = sqrt(dtmp1);
  dtmp1 = dtmp1 - (Z - dz);
  double dtmp2 = (X * X + Y * Y + (Z + dz) * (Z + dz));
  dtmp2 = sqrt(dtmp2);
  dtmp2 = dtmp2 - (Z + dz);
  double dtmp3 = log(dtmp1 / dtmp2);
  return (2.0 * ST_PI * rW * dtmp3);
}  // ImprovedP_W ends

double ImprovedFX_W(double rW, double lW, double X, double Y, double Z) {
  if (DebugISLES) printf("In ImprovedFX_W ...\n");

  double dist = sqrt(X * X + Y * Y + Z * Z);
  if (dist < MINDIST) {
    // distance less than MINDIST
    return (0.0);
  } else if ((fabs(X) < MINDIST) && (fabs(Y) < MINDIST)) {
    // point on the axis of the wire element
    return (0.0);
  } else {
    // point far away from the wire element
    double C = (Z) - (lW / 2.0);
    double D = (Z) + (lW / 2.0);
    double A = sqrt(X * X + Y * Y + C * C);
    double B = sqrt(X * X + Y * Y + D * D);

    double tmp1 = 2.0 * ST_PI * rW;
    double tmp2 = (X / (A * (B - D))) - ((A - C) * X / (B * (B - D) * (B - D)));
    double tmp3 = (B - D) / (A - C);
    return (-1.0 * tmp1 * tmp2 * tmp3);
  }  // dist ... if ... else if ... else
}  // ImprovedFX_W ends

double ImprovedFY_W(double rW, double lW, double X, double Y, double Z) {
  if (DebugISLES) printf("In ImprovedFY_W ...\n");

  double dist = sqrt(X * X + Y * Y + Z * Z);
  if (dist < MINDIST) {
    // distance less than MINDIST
    return (0.0);
  } else if ((fabs(X) < MINDIST) && (fabs(Y) < MINDIST)) {
    // point on the axis of the wire element
    return (0.0);
  } else {
    // point far away from the wire element
    double C = (Z) - (lW / 2.0);
    double D = (Z) + (lW / 2.0);
    double A = sqrt(X * X + Y * Y + C * C);
    double B = sqrt(X * X + Y * Y + D * D);

    double tmp1 = 2.0 * ST_PI * rW;
    double tmp2 = (Y / (A * (B - D))) - ((A - C) * Y / (B * (B - D) * (B - D)));
    double tmp3 = (B - D) / (A - C);
    return (-1.0 * tmp1 * tmp2 * tmp3);
  }  // dist ... if ... else if ... else
}  // ImprovedFY_W ends

double ImprovedFZ_W(double rW, double lW, double X, double Y, double Z) {
  if (DebugISLES) printf("In ImprovedFZ_W ...\n");

  double dist = sqrt(X * X + Y * Y + Z * Z);
  if (dist < MINDIST) {
    // distance less than MINDIST
    return (0.0);
  } else if ((fabs(X) < MINDIST) && (fabs(Y) < MINDIST)) {
    // point on the axis of the wire element
    double C = Z - (lW / 2.0);
    double D = Z + (lW / 2.0);
    double A = sqrt(X * X + Y * Y + C * C);
    double B = sqrt(X * X + Y * Y + D * D);

    double tmp1 = 2.0 * ST_PI * rW;
    double tmp2 = (1.0 / B) - (1.0 / A);
    return (-1.0 * tmp1 * tmp2);
  } else {
    // point far away from the wire element
    double C = Z - (lW / 2.0);
    double D = Z + (lW / 2.0);
    double A = sqrt(X * X + Y * Y + C * C);
    double B = sqrt(X * X + Y * Y + D * D);

    double tmp1 = 2.0 * ST_PI * rW;
    double tmp2 = (1.0 / B) - (1.0 / A);
    return (-1.0 * tmp1 * tmp2);
  }  // dist ... if ... else if ... else
}  // ImprovedFZ_W ends

int ImprovedWire(double rW, double lW, double X, double Y, double Z,
                 double *potential, Vector3D *Flux) {
  if (DebugISLES) printf("In ImprovedWire ...\n");

  double dz = 0.5 * lW;  // half length of the wire segment

  double dtmp1 = (X * X + Y * Y + (Z - dz) * (Z - dz));
  dtmp1 = sqrt(dtmp1);
  dtmp1 = dtmp1 - (Z - dz);
  double dtmp2 = (X * X + Y * Y + (Z + dz) * (Z + dz));
  dtmp2 = sqrt(dtmp2);
  dtmp2 = dtmp2 - (Z + dz);
  double dtmp3 = log(dtmp1 / dtmp2);
  *potential = 2.0 * ST_PI * rW * dtmp3;

  double dist = sqrt(X * X + Y * Y + Z * Z);
  double Fx = 0.0, Fy = 0.0, Fz = 0.0;

  if (dist < MINDIST) {
    // distance less than MINDIST from centroid
    Fx = 0.0;
    Fy = 0.0;
    Fz = 0.0;
  } else if ((fabs(X) < MINDIST) && (fabs(Y) < MINDIST)) {
    // point on the axis of the wire element
    Fx = 0.0;
    Fy = 0.0;

    double C = Z - (lW / 2.0);
    double D = Z + (lW / 2.0);
    double A = sqrt(X * X + Y * Y + C * C);
    double B = sqrt(X * X + Y * Y + D * D);

    double tmp1 = 2.0 * ST_PI * rW;
    double tmp2 = (1.0 / B) - (1.0 / A);
    Fz = -1.0 * tmp1 * tmp2;
  } else {
    // point far away from the wire element
    double C = (Z) - (lW / 2.0);
    double D = (Z) + (lW / 2.0);
    double A = sqrt(X * X + Y * Y + C * C);
    double B = sqrt(X * X + Y * Y + D * D);

    double tmp1 = 2.0 * ST_PI * rW;
    double tmp2 = (X / (A * (B - D))) - ((A - C) * X / (B * (B - D) * (B - D)));
    double tmp3 = (B - D) / (A - C);
    Fx = -1.0 * tmp1 * tmp2 * tmp3;

    tmp2 = (Y / (A * (B - D))) - ((A - C) * Y / (B * (B - D) * (B - D)));
    Fy = -1.0 * tmp1 * tmp2 * tmp3;

    tmp2 = (1.0 / B) - (1.0 / A);
    Fz = -1.0 * tmp1 * tmp2;
  }  // dist ... if ... else if ... else

  Flux->X = Fx;
  Flux->Y = Fy;
  Flux->Z = Fz;

  return 0;
}

// Exact Potential due to thin wire at an arbitrary location
// Z axis is along the length of the wire
double ExactThinP_W(double rW, double lW, double X, double Y, double Z) {
  if (DebugISLES) {
    printf("In ExactThinP_W ...\n");
    printf("rW: %lg, lW: %lg, X: %lg, Y: %lg, Z: %lg\n", rW, lW, X, Y, Z);
  }
  // Expressions from PFXYZ_thin.m
  const double h = 0.5 * lW;
  const double r2 = X * X + Y * Y;
  const double a = h + Z;
  const double b = h - Z;
  double Pot = log((a + sqrt(r2 + a * a)) * (b + sqrt(r2 + b * b)) / r2);
  return 2. * ST_PI * rW * Pot;
}  // ExactThinP_W ends

// Exact FX due to thin wire at an arbitrary location
// Repeated execution of similar expression - can be optimized
double ExactThinFX_W(double rW, double lW, double X, double Y, double Z) {
  if (DebugISLES) {
    printf("In ExactThinFX_W ...\n");
    printf("rW: %lg, lW: %lg, X: %lg, Y: %lg, Z: %lg\n", rW, lW, X, Y, Z);
  }
  // Expressions from PFXYZ_thin.m
  double h = 0.5 * lW;
  double Fx = 2.0 * X *
              (sqrt(X * X + Y * Y + Z * Z + h * h + 2.0 * Z * h) * h -
               sqrt(X * X + Y * Y + Z * Z + h * h + 2.0 * Z * h) * Z +
               sqrt(X * X + Y * Y + Z * Z - 2.0 * Z * h + h * h) * h +
               sqrt(X * X + Y * Y + Z * Z - 2.0 * Z * h + h * h) * Z) /
              (X * X + Y * Y) /
              sqrt(X * X + Y * Y + Z * Z - 2.0 * Z * h + h * h) /
              sqrt(X * X + Y * Y + Z * Z + h * h + 2.0 * Z * h) * ST_PI;
  return (rW * Fx);
}  // ExactThinFX_W ends

// Exact FY due to thin wire at an arbitrary location
// Repeated execution of similar expression - can be optimized
double ExactThinFY_W(double rW, double lW, double X, double Y, double Z) {
  if (DebugISLES) {
    printf("In ExactThinFY_W ...\n");
    printf("rW: %lg, lW: %lg, X: %lg, Y: %lg, Z: %lg\n", rW, lW, X, Y, Z);
  }
  // Expressions from PFXYZ_thin.m
  double h = 0.5 * lW;
  double Fy = 2.0 * Y *
              (sqrt(X * X + Y * Y + Z * Z + h * h + 2.0 * Z * h) * h -
               sqrt(X * X + Y * Y + Z * Z + h * h + 2.0 * Z * h) * Z +
               sqrt(X * X + Y * Y + Z * Z - 2.0 * Z * h + h * h) * h +
               sqrt(X * X + Y * Y + Z * Z - 2.0 * Z * h + h * h) * Z) /
              (X * X + Y * Y) /
              sqrt(X * X + Y * Y + Z * Z - 2.0 * Z * h + h * h) /
              sqrt(X * X + Y * Y + Z * Z + h * h + 2.0 * Z * h) * ST_PI;
  return (rW * Fy);
}  // ExactThinFY_W ends

// Exact FZ (axial field) due to thin wire at an arbitrary location
// Repeated execution of similar expression - can be optimized
double ExactThinFZ_W(double rW, double lW, double X, double Y, double Z) {
  if (DebugISLES) {
    printf("In ExactThinFZ_W ...\n");
    printf("rW: %lg, lW: %lg, X: %lg, Y: %lg, Z: %lg\n", rW, lW, X, Y, Z);
  }

  // Expressions from PFXYZ_thin.m
  double h = 0.5 * lW;
  double Fz = 2.0 *
              (sqrt(X * X + Y * Y + Z * Z + h * h + 2.0 * Z * h) -
               sqrt(X * X + Y * Y + Z * Z - 2.0 * Z * h + h * h)) /
              sqrt(X * X + Y * Y + Z * Z - 2.0 * Z * h + h * h) /
              sqrt(X * X + Y * Y + Z * Z + h * h + 2.0 * Z * h) * ST_PI;
  return (rW * Fz);
}  // ExactThinFZ_W ends

int ExactThinWire(double rW, double lW, double X, double Y, double Z,
                  double *potential, Vector3D *Flux) {
  if (DebugISLES) {
    printf("In ExactThinWire_W ...\n");
    printf("rW: %lg, lW: %lg, X: %lg, Y: %lg, Z: %lg\n", rW, lW, X, Y, Z);
  }
  const double h = 0.5 * lW;
  const double r2 = X * X + Y * Y;
  const double a = h + Z;
  const double b = h - Z;
  const double c = sqrt(r2 + a * a);
  const double d = sqrt(r2 + b * b);
  const double s = 2. * ST_PI * rW;
  *potential = s * log((a + c) * (b + d) / r2);
  double f = s * (c * b + d * a) / (r2 * c * d);
  Flux->X = f * X;
  Flux->Y = f * Y;
  Flux->Z = s * (c - d) / (c * d);
  return 0;
}

// Radius of ring is 'a'
// The point is passed in the local coordinate system, in the ECS of the ring
// In the ECS, plane of the ring is r,theta (XY) and z is the ordinate
// up and down of this plane.
// The charge density is assumed to be one.
// Reference: Fonte 2013 - Survey of physical modelling in Resistive Plate
// Chambers.pdf, equation (4.52)

ISLESGLOBAL int ExactRingPF(double a, Point3D localPt, double *potential,
                            Vector3D *Flux) {
  int dbgFn = 0;

  // Coordinate transformation as performed for lines is necessary
  // Origin of the system is the center of the ring
  // The plane of the ring will need to be identified
  // An axis vertical to the plane will be the Y-axis
  // An axis orthogonal to the Y-axis will be X-axis (as done for lines)
  // Cross between the X and Y units will yield unit Z
  // Input parameter will need to be added - we now have only the radius
  // A way of identifying the plane of the ring is needed. It is quite likel
  // that the normal to the ring plane will be provided as an input. That
  // will immediately give us unit Y. The rest can be easily deduced.

  double roe = localPt.X * localPt.X + localPt.Y * localPt.Y;
  roe = sqrt(roe);
  double z = localPt.Z;
  double r1dot = ((a + roe) * (a + roe)) + z * z;
  r1dot = sqrt(r1dot);
  double u = 2.0 * sqrt(a * roe) / r1dot;

  // K1 and K2 are complete elliptic integrals
  // K implies first kind, according to GSL (and Wikipedia) convention
  gsl_mode_t mode = GSL_PREC_DOUBLE;
  double K1, K2;
  // Following are the routines with debugging information
  if (dbgFn) {
    gsl_sf_result rK1, rK2;
    gsl_sf_ellint_Kcomp_e(u, mode, &rK1);
    gsl_sf_ellint_Ecomp_e(u, mode, &rK2);
    K1 = rK1.val;
    K2 = rK2.val;
  } else {  // If debugging information is not necessary
    K1 = gsl_sf_ellint_Kcomp(u, mode);
    K2 = gsl_sf_ellint_Ecomp(u, mode);
  }

  double Vring = (a / ST_PI) * K1 / r1dot;
  // field in the radial direction
  double Eringroe = (a / ST_PI) * ((1.0 / (2.0 * r1dot * roe)) *
                (K1 - (((a * a - roe * roe + z * z) * K2) /
                       (r1dot * r1dot * (1.0 - u * u)))));
  // field in the z direction
  double Eringz = (a / ST_PI) * 
                  ((z * K2) / ((r1dot * r1dot * r1dot) * (1.0 - u * u)));

  *potential = Vring;
  Flux->X = Eringroe * localPt.X / roe;
  Flux->Y = Eringroe * localPt.Y / roe;
  Flux->Z = Eringz;

  return 0;
}  // ExactRingPF ends

// Exact potential and GCS flux due to a point source
double PointKnChPF(Point3D SourcePt, Point3D FieldPt, Vector3D *globalF) {
  double dist = GetDistancePoint3D(&SourcePt, &FieldPt);
  double dist3 = pow(dist, 3.0);

  if (dist3 < MINDIST3) {
    globalF->X = 0.0;
    globalF->Y = 0.0;
    globalF->Z = 0.0;
  } else {
    globalF->X = (FieldPt.X - SourcePt.X) / dist3;
    globalF->Y = (FieldPt.Y - SourcePt.Y) / dist3;
    globalF->Z = (FieldPt.Z - SourcePt.Z) / dist3;
  }

  if (dist < MINDIST)
    return (0.0);
  else
    return (1.0 / dist);
}  // PointKnChPF ends

// Pass start and stop points of the line and get back potential and GCS flux
// at field point assuming unit linear charge density.
double LineKnChPF(Point3D LineStart, Point3D LineStop, Point3D FieldPt,
                  Vector3D *globalF) {
  int debugFn = 0;

  if (debugFn) {
    printf("In LineKnChPF:\n");
    printf("LineStart: %lg %lg %lg\n", LineStart.X, LineStart.Y, LineStart.Z);
    printf("LineStop: %lg %lg %lg\n", LineStop.X, LineStop.Y, LineStop.Z);
    printf("FieldPt: %lg %lg %lg\n", FieldPt.X, FieldPt.Y, FieldPt.Z);
  }

  double xvert[2], yvert[2], zvert[2];
  xvert[0] = LineStart.X;
  xvert[1] = LineStop.X;
  yvert[0] = LineStart.Y;
  yvert[1] = LineStop.Y;
  zvert[0] = LineStart.Z;
  zvert[1] = LineStop.Z;

  double xfld = FieldPt.X;
  double yfld = FieldPt.Y;
  double zfld = FieldPt.Z;

  double xorigin = 0.5 * (xvert[1] + xvert[0]);
  double yorigin = 0.5 * (yvert[1] + yvert[0]);
  double zorigin = 0.5 * (zvert[1] + zvert[0]);
  double LZ = GetDistancePoint3D(&LineStop, &LineStart);
  if (debugFn) {
    printf("xorigin: %lg, yorigin: %lg, zorigin: %lg\n", xorigin, yorigin,
           zorigin);
    printf("LZ: %lg\n", LZ);
  }

  // Create a local coordinate system for which the z axis is along the length
  // of the line
  // FieldPt needs to be transformed to localPt

  // Find direction cosines of the line charge to initiate the transformation.
  DirnCosn3D DirCos;

  // Copied from the DiscretizeWire function of neBEM/src/PreProcess/ReTrim.c
  // Direction cosines along the wire - note difference from surface primitives!
  // The direction along the wire is considered to be the z axis of the LCS.
  // So, let us fix that axial vector first
  DirCos.ZUnit.X = (xvert[1] - xvert[0]) / LZ;  // useful
  DirCos.ZUnit.Y = (yvert[1] - yvert[0]) / LZ;
  DirCos.ZUnit.Z = (zvert[1] - zvert[0]) / LZ;  // useful
  // Now direction cosines for the X and Y axes.
  {
    Vector3D XUnit, YUnit, ZUnit;
    ZUnit.X = DirCos.ZUnit.X;
    ZUnit.Y = DirCos.ZUnit.Y;
    ZUnit.Z = DirCos.ZUnit.Z;

    // Find any vector that is orthogonal to the Z unit vector, and make that
    // the X-axis
    int caseComp = 1;
    double maxComp = fabs(ZUnit.X);
    if (fabs(ZUnit.Y) > maxComp) {
      caseComp = 2;
      maxComp = fabs(ZUnit.Y);
    }
    if (fabs(ZUnit.Z) > maxComp) {
      caseComp = 3;
      maxComp = fabs(ZUnit.Z);
    }

    switch (caseComp) {
      case 1:  // ZUnit.X is larger than the other direction cosines
        XUnit.Y = 1.0;
        XUnit.Z = 1.0;
        XUnit.X = (ZUnit.Y + ZUnit.Z) / ZUnit.X;
        break;
      case 2:  // ZUnit.Y is larger than the other direction cosines
        XUnit.X = 1.0;
        XUnit.Z = 1.0;
        XUnit.Y = (ZUnit.X + ZUnit.Z) / ZUnit.Y;
        break;
      case 3:  // ZUnit.Z is larger than the other direction cosines
        XUnit.X = 1.0;
        XUnit.Y = 1.0;
        XUnit.Z = (ZUnit.X + ZUnit.Y) / ZUnit.Z;
        break;
      default:;
    }                              // switch ends
    XUnit = UnitVector3D(&XUnit);  // direction cosines for X is created

    // y-Axis: vectorial product of axes 1 and 2.
    YUnit = Vector3DCrossProduct(&ZUnit, &XUnit);
    YUnit = UnitVector3D(&YUnit);
    // end of replacement

    DirCos.XUnit.X = XUnit.X;
    DirCos.XUnit.Y = XUnit.Y;
    DirCos.XUnit.Z = XUnit.Z;
    DirCos.YUnit.X = YUnit.X;
    DirCos.YUnit.Y = YUnit.Y;
    DirCos.YUnit.Z = YUnit.Z;
  }  // X and Y direction cosines computed

  Point3D localPt;
  // Through InitialVector[], field point gets translated to ECS origin.
  // Axes direction are, however, still global which when rotated to ECS
  // system, yields FinalVector[].
  {  // Rotate point3D from global to local system to get localPt.
    double InitialVector[3] = {xfld - xorigin, yfld - yorigin, zfld - zorigin};
    double TransformationMatrix[3][3] = {{0.0, 0.0, 0.0},
                                         {0.0, 0.0, 0.0},
                                         {0.0, 0.0, 0.0}};
    TransformationMatrix[0][0] = DirCos.XUnit.X;
    TransformationMatrix[0][1] = DirCos.XUnit.Y;
    TransformationMatrix[0][2] = DirCos.XUnit.Z;
    TransformationMatrix[1][0] = DirCos.YUnit.X;
    TransformationMatrix[1][1] = DirCos.YUnit.Y;
    TransformationMatrix[1][2] = DirCos.YUnit.Z;
    TransformationMatrix[2][0] = DirCos.ZUnit.X;
    TransformationMatrix[2][1] = DirCos.ZUnit.Y;
    TransformationMatrix[2][2] = DirCos.ZUnit.Z;
    double FinalVector[3];

    for (int i = 0; i < 3; ++i) {
      FinalVector[i] = 0.0;
      for (int j = 0; j < 3; ++j) {
        FinalVector[i] += TransformationMatrix[i][j] * InitialVector[j];
      }
    }

    localPt.X = FinalVector[0];
    localPt.Y = FinalVector[1];
    localPt.Z = FinalVector[2];
  }  // Point3D rotated

  double Pot;
  Vector3D localF;
  double xpt = localPt.X;
  double ypt = localPt.Y;
  double zpt = localPt.Z;

  if (debugFn) {
    // printf("LineStart: %lg %lg %lg\n", LineStart.X, LineStart.Y,
    // LineStart.Z); printf("LineStop: %lg %lg %lg\n", LineStop.X, LineStop.Y,
    // LineStop.Z);
    printf("localPt: %lg %lg %lg\n", localPt.X, localPt.Y, localPt.Z);
    printf("LZ: %lg\n", LZ);
  }

  // field point from element centroid
  double zptplus = zpt + 0.5 * LZ;
  double zptminus = zpt - 0.5 * LZ;

  // NOTE: MatLab15, Maple Tanay's expressions are the same, only rearranged in
  // different ways. In MatLab, there are several possibilities indicated, but
  // Maple produces only the reasonable ones.

  /*
  // Maple expressions: DefLineElementTry1.mw
  gsl_complex tmpcmp;
  tmpcmp.dat[0] = zptplus / sqrt(xpt*xpt + ypt*ypt); tmpcmp.dat[1] = 0.0;
  gsl_complex Pot1 = gsl_complex_arcsinh(tmpcmp);
  tmpcmp.dat[0] = zptminus / sqrt(xpt*xpt + ypt*ypt); tmpcmp.dat[1] = 0.0;
  gsl_complex Pot2 = gsl_complex_arcsinh(tmpcmp);
  Pot = Pot1.dat[0] - Pot2.dat[0];
  // The following expressions need to be rewritten after consulting the Maple
  // print-out.
  localF.X = xpt *((zptplus*sqrt(xpt*xpt + ypt*ypt + zptminus*zptminus) -
                    zptminus*sqrt(xpt*xpt + ypt*ypt + zptplus*zptplus)) / 
                   (sqrt(xpt*xpt + ypt*ypt + zptminus*zptminus) * 
                    sqrt(xpt*xpt + ypt*ypt + zptplus*zptplus)) * 
                   (1.0/sqrt(xpt*xpt + ypt*ypt));
  localF.Y = ypt *((zptplus*sqrt(xpt*xpt + ypt*ypt + zptminus*zptminus) -
                    zptminus*sqrt(xpt*xpt + ypt*ypt + zptplus*zptplus)) / 
                   (sqrt(xpt*xpt + ypt*ypt + zptminus*zptminus) * 
                    sqrt(xpt*xpt + ypt*ypt + zptplus*zptplus) * 
                    sqrt(xpt*xpt + ypt*ypt));
  localF.Z = (sqrt(xpt*xpt + ypt*ypt + zptplus*zptplus) - 
              sqrt(xpt*xpt + ypt*ypt + zptminus*zptminus)) / 
             (sqrt(xpt*xpt + ypt*ypt + zptplus*zptplus) * 
              sqrt(xpt*xpt + ypt*ypt + zptminus*zptminus));

  if (debugFn) {
    printf("Using Maple expressions =>\n");
    printf("Pot1: %lg + i %lg\n", Pot1.dat[0], Pot1.dat[1]);
    printf("Pot2: %lg + i %lg\n", Pot2.dat[0], Pot2.dat[1]);
    printf("Pot: %lg\n", Pot);
    printf("localF: %lg %lg %lg\n", localF.X, localF.Y, localF.Z);
  }
  */

  // MatLab expressions: PFXYZ_MatLab15_4.out
  double tmpxy = xpt * xpt + ypt * ypt;
  double tmpzplus = sqrt(zptplus * zptplus + xpt * xpt + ypt * ypt);
  double tmpzminus = sqrt(xpt * xpt + ypt * ypt + zptminus * zptminus);
  Pot =
      log(zptplus + tmpzplus) - log(zptminus + tmpzminus);
  localF.X = xpt * ((zptplus / tmpxy / tmpzplus) -
                    (zptminus / tmpxy / tmpzminus));
  localF.Y = ypt * ((zptplus / tmpxy / tmpzplus) -
                    (zptminus / tmpxy / tmpzminus));
  localF.Z = (1.0 / tmpzminus) - (1.0 / tmpzplus);
  if (debugFn) {
    printf("Using MatLab expressions =>\n");
    printf("Pot: %lg\n", Pot);
    printf("localF: %lg %lg %lg\n", localF.X, localF.Y, localF.Z);
  }

  /*
  // Tanay expressions - here the line is parallel to Z-axis, rather than Y.
  // According to Tanay's paper, the lines are parallel to Y.
  Pot = 0.0;	// Tanay does not have an expression for the potential.
  // borrow from Matlab
  // Pot =  log(zptplus + sqrt(zptplus*zptplus+xpt*xpt+ypt*ypt)) -
  //        log(zptminus + sqrt(xpt*xpt+ypt*ypt+zptminus*zptminus)); 
  localF.X = (xpt / (xpt*xpt + ypt*ypt)) * 
             ((zptplus / sqrt(zptplus*zptplus + xpt*xpt + ypt*ypt)) - 
              (zptminus / sqrt(zptminus*zptminus + xpt*xpt + ypt*ypt))); 
  localF.Y = (ypt / (xpt*xpt + ypt*ypt)) * 
             ((zptplus / sqrt(zptplus*zptplus + xpt*xpt + ypt*ypt)) -
              (zptminus / sqrt(zptminus*zptminus + xpt*xpt + ypt*ypt))); 
  localF.Z = - ((1.0 / sqrt(zptplus*zptplus + xpt*xpt + ypt*ypt))	-
                (1.0 / sqrt(zptminus*zptminus + xpt*xpt + ypt*ypt))); 
  if (debugFn) {
    printf("Using Tanay expressions =>\n");
    printf("Pot: %lg\n", Pot);
    printf("localF: %lg %lg %lg\n", localF.X, localF.Y, localF.Z);
  }
  */

  if (debugFn) {
    printf("Pot: %lg\n", Pot);
    printf("localF: %lg %lg %lg\n", localF.X, localF.Y, localF.Z);
  }

  (*globalF) = RotateVector3D(&localF, &DirCos, local2global);
  if (debugFn) {
    printf("globalF: %lg %lg %lg\n", globalF->X, globalF->Y, globalF->Z);
    exit(-1);
  }

  return (Pot);
}  // LineKnChPF ends

// Pass start and stop points of the wire and get back potential and GCS flux
// at field point assuming unit surface charge density.
double WireKnChPF(Point3D WireStart, Point3D WireStop, double radius,
                  Point3D FieldPt, Vector3D *globalF) {
  int debugFn = 1;

  if (debugFn) {
    printf("In WireKnChPF:\n");
    printf("WireStart: %lg %lg %lg\n", WireStart.X, WireStart.Y, WireStart.Z);
    printf("WireStop: %lg %lg %lg\n", WireStop.X, WireStop.Y, WireStop.Z);
    printf("radius: %lg\n", radius);
    printf("FieldPt: %lg %lg %lg\n", FieldPt.X, FieldPt.Y, FieldPt.Z);
  }

  double xvert[2], yvert[2], zvert[2];
  xvert[0] = WireStart.X;
  xvert[1] = WireStop.X;
  yvert[0] = WireStart.Y;
  yvert[1] = WireStop.Y;
  zvert[0] = WireStart.Z;
  zvert[1] = WireStop.Z;

  double xfld = FieldPt.X;
  double yfld = FieldPt.Y;
  double zfld = FieldPt.Z;

  double xorigin = 0.5 * (xvert[1] + xvert[0]);
  double yorigin = 0.5 * (yvert[1] + yvert[0]);
  double zorigin = 0.5 * (zvert[1] + zvert[0]);
  double LZ = GetDistancePoint3D(&WireStop, &WireStart);
  if (debugFn) {
    printf("xorigin: %lg, yorigin: %lg, zorigin: %lg\n", xorigin, yorigin,
           zorigin);
    printf("LZ: %lg\n", LZ);
  }

  // Create a local coordinate system for which the z axis is along the length
  // of the wire
  // FieldPt needs to be transformed to localPt

  // Find direction cosines of the line charge to initiate the transformation.
  DirnCosn3D DirCos;

  // Copied from the DiscretizeWire function of neBEM/src/PreProcess/ReTrim.c
  // Direction cosines along the wire - note difference from surface primitives!
  // The direction along the wire is considered to be the z axis of the LCS.
  // So, let us first fix that as the axial vector.
  // Check whether this method is appropriate. If found not ok, we may need
  // to change ReTrim.c as well.
  // Finding the unit X and unit Y could work as in line elements.
  DirCos.ZUnit.X = (xvert[1] - xvert[0]) / LZ;  // useful
  DirCos.ZUnit.Y = (yvert[1] - yvert[0]) / LZ;
  DirCos.ZUnit.Z = (zvert[1] - zvert[0]) / LZ;  // useful
  // Now direction cosines for the X and Y axes.
  {
    Vector3D XUnit, YUnit, ZUnit;
    ZUnit.X = DirCos.ZUnit.X;
    ZUnit.Y = DirCos.ZUnit.Y;
    ZUnit.Z = DirCos.ZUnit.Z;

    if (fabs(ZUnit.X) >= fabs(ZUnit.Z) && fabs(ZUnit.Y) >= fabs(ZUnit.Z)) {
      XUnit.X = -ZUnit.Y;
      XUnit.Y = ZUnit.X;
      XUnit.Z = 0.0;
    } else if (fabs(ZUnit.X) >= fabs(ZUnit.Y) &&
               fabs(ZUnit.Z) >= fabs(ZUnit.Y)) {
      XUnit.X = -ZUnit.Z;
      XUnit.Y = 0.0;
      XUnit.Z = ZUnit.X;
    } else {
      XUnit.X = 0.0;
      XUnit.Y = ZUnit.Z;
      XUnit.Z = -ZUnit.Y;
    }
    XUnit = UnitVector3D(&XUnit);

    // y-Axis: vectorial product of axes 1 and 2.
    YUnit = Vector3DCrossProduct(&ZUnit, &XUnit);
    YUnit = UnitVector3D(&YUnit);
    // end of replacement

    DirCos.XUnit.X = XUnit.X;
    DirCos.XUnit.Y = XUnit.Y;
    DirCos.XUnit.Z = XUnit.Z;
    DirCos.YUnit.X = YUnit.X;
    DirCos.YUnit.Y = YUnit.Y;
    DirCos.YUnit.Z = YUnit.Z;
  }  // X and Y direction cosines computed

  Point3D localPt;
  // Through InitialVector[], field point gets translated to ECS origin.
  // Axes direction are, however, still global which when rotated to ECS
  // system, yields FinalVector[].
  {  // Rotate point3D from global to local system to get localPt.
    double InitialVector[3] = {xfld - xorigin, yfld - yorigin, zfld - zorigin};
    double TransformationMatrix[3][3] = {{0.0, 0.0, 0.0},
                                         {0.0, 0.0, 0.0},
                                         {0.0, 0.0, 0.0}};
    TransformationMatrix[0][0] = DirCos.XUnit.X;
    TransformationMatrix[0][1] = DirCos.XUnit.Y;
    TransformationMatrix[0][2] = DirCos.XUnit.Z;
    TransformationMatrix[1][0] = DirCos.YUnit.X;
    TransformationMatrix[1][1] = DirCos.YUnit.Y;
    TransformationMatrix[1][2] = DirCos.YUnit.Z;
    TransformationMatrix[2][0] = DirCos.ZUnit.X;
    TransformationMatrix[2][1] = DirCos.ZUnit.Y;
    TransformationMatrix[2][2] = DirCos.ZUnit.Z;
    double FinalVector[3];

    for (int i = 0; i < 3; ++i) {
      FinalVector[i] = 0.0;
      for (int j = 0; j < 3; ++j) {
        FinalVector[i] += TransformationMatrix[i][j] * InitialVector[j];
      }
    }

    localPt.X = FinalVector[0];
    localPt.Y = FinalVector[1];
    localPt.Z = FinalVector[2];
  }  // Point3D rotated

  double Pot;
  Vector3D localF;
  double xpt = localPt.X;
  double ypt = localPt.Y;
  double zpt = localPt.Z;

  if (radius < MINDIST) radius = MINDIST;
  double rW = radius;
  double lW = LZ;

  if (debugFn) {
    // printf("WireStart: %lg %lg %lg\n", WireStart.X, WireStart.Y,
    // WireStart.Z); printf("WireStop: %lg %lg %lg\n", WireStop.X, WireStop.Y,
    // WireStop.Z);
    printf("localPt: %lg %lg %lg\n", localPt.X, localPt.Y, localPt.Z);
    printf("rW: %lg, lW: %lg\n", rW, lW);
  }

  // field point from element centroid
  double dist = sqrt(xpt * xpt + ypt * ypt + zpt * zpt);

  if (dist >= FarField * lW) {
    double dA = 2.0 * ST_PI * rW * lW;
    Pot = dA / dist;
    double dist3 = dist * dist * dist;
    localF.X = dA * xpt / dist3;
    localF.Y = dA * ypt / dist3;
    localF.Z = dA * zpt / dist3;
  } else if ((fabs(xpt) < MINDIST) && (fabs(ypt) < MINDIST) &&
             (fabs(zpt) < MINDIST)) {
    Pot = ExactCentroidalP_W(rW, lW);
    localF.X = 0.0;  // CHECK - these flux values need to be confirmed
    localF.Y = 0.0;
    localF.Z = 0.0;
  } else if ((fabs(xpt) < MINDIST) && (fabs(ypt) < MINDIST)) {
    Pot = ExactAxialP_W(rW, lW, zpt);
    localF.X = localF.Y = 0.0;
    localF.Z = ExactThinFZ_W(rW, lW, xpt, ypt, zpt);
  } else {
    Pot = ExactThinP_W(rW, lW, xpt, ypt, zpt);
    localF.X = ExactThinFX_W(rW, lW, xpt, ypt, zpt);
    localF.Y = ExactThinFY_W(rW, lW, xpt, ypt, zpt);
    localF.Z = ExactThinFZ_W(rW, lW, xpt, ypt, zpt);
  }

  if (debugFn) {
    printf("Pot: %lg\n", Pot);
    printf("localF: %lg %lg %lg\n", localF.X, localF.Y, localF.Z);
  }

  (*globalF) = RotateVector3D(&localF, &DirCos, local2global);
  if (debugFn) {
    printf("globalF: %lg %lg %lg\n", globalF->X, globalF->Y, globalF->Z);
    exit(-1);
  }

  return (Pot);
}  // WireKnChPF ends

// Pass four vertices of the area and get back potential and GCS flux
// at field point assuming unit charge density.
double AreaKnChPF(int NbVertices, Point3D *Vertex, Point3D FieldPt,
                  Vector3D *globalF) {
  if (NbVertices != 4) {
    printf(
        "Only rectangular areas with known charges allowed at present ...\n");
    exit(-1);
  }

  double xvert[4], yvert[4], zvert[4];
  xvert[0] = Vertex[0].X;
  xvert[1] = Vertex[1].X;
  xvert[2] = Vertex[2].X;
  xvert[3] = Vertex[3].X;
  yvert[0] = Vertex[0].Y;
  yvert[1] = Vertex[1].Y;
  yvert[2] = Vertex[2].Y;
  yvert[3] = Vertex[3].Y;
  zvert[0] = Vertex[0].Z;
  zvert[1] = Vertex[1].Z;
  zvert[2] = Vertex[2].Z;
  zvert[3] = Vertex[3].Z;

  double xfld = FieldPt.X;
  double yfld = FieldPt.Y;
  double zfld = FieldPt.Z;

  double xorigin = 0.25 * (xvert[3] + xvert[2] + xvert[1] + xvert[0]);
  double yorigin = 0.25 * (yvert[3] + yvert[2] + yvert[1] + yvert[0]);
  double zorigin = 0.25 * (zvert[3] + zvert[2] + zvert[1] + zvert[0]);

  // lengths of the sides
  double SurfLX = sqrt((xvert[1] - xvert[0]) * (xvert[1] - xvert[0]) +
                       (yvert[1] - yvert[0]) * (yvert[1] - yvert[0]) +
                       (zvert[1] - zvert[0]) * (zvert[1] - zvert[0]));
  double SurfLZ = sqrt((xvert[2] - xvert[1]) * (xvert[2] - xvert[1]) +
                       (yvert[2] - yvert[1]) * (yvert[2] - yvert[1]) +
                       (zvert[2] - zvert[1]) * (zvert[2] - zvert[1]));

  // Direction cosines - CHECK!
  //
  //   3      2
  //   ________
  //   |      |
  //   |      |
  //   |      |
  //   --------    -> +ve x-axis
  //   0      1
  //   |
  //   |
  //   V +ve z-axis
  //   The resulting +ve y-axis is out of the paper towards the reader.
  DirnCosn3D DirCos;
  DirCos.XUnit.X = (xvert[1] - xvert[0]) / SurfLX;
  DirCos.XUnit.Y = (yvert[1] - yvert[0]) / SurfLX;
  DirCos.XUnit.Z = (zvert[1] - zvert[0]) / SurfLX;
  DirCos.ZUnit.X = (xvert[0] - xvert[3]) / SurfLZ;
  DirCos.ZUnit.Y = (yvert[0] - yvert[3]) / SurfLZ;
  DirCos.ZUnit.Z = (zvert[0] - zvert[3]) / SurfLZ;
  DirCos.YUnit = Vector3DCrossProduct(&DirCos.ZUnit, &DirCos.XUnit);

  Point3D localPt;
  // Through InitialVector[], field point gets translated to ECS origin.
  // Axes direction are, however, still global which when rotated to ECS
  // system, yields FinalVector[].
  {  // Rotate point3D from global to local system to get localPt.
    double InitialVector[3] = {xfld - xorigin, yfld - yorigin, zfld - zorigin};
    double TransformationMatrix[3][3] = {{0.0, 0.0, 0.0},
                                         {0.0, 0.0, 0.0},
                                         {0.0, 0.0, 0.0}};
    TransformationMatrix[0][0] = DirCos.XUnit.X;
    TransformationMatrix[0][1] = DirCos.XUnit.Y;
    TransformationMatrix[0][2] = DirCos.XUnit.Z;
    TransformationMatrix[1][0] = DirCos.YUnit.X;
    TransformationMatrix[1][1] = DirCos.YUnit.Y;
    TransformationMatrix[1][2] = DirCos.YUnit.Z;
    TransformationMatrix[2][0] = DirCos.ZUnit.X;
    TransformationMatrix[2][1] = DirCos.ZUnit.Y;
    TransformationMatrix[2][2] = DirCos.ZUnit.Z;
    double FinalVector[3];

    for (int i = 0; i < 3; ++i) {
      FinalVector[i] = 0.0;
      for (int j = 0; j < 3; ++j) {
        FinalVector[i] += TransformationMatrix[i][j] * InitialVector[j];
      }
    }

    localPt.X = FinalVector[0];
    localPt.Y = FinalVector[1];
    localPt.Z = FinalVector[2];
  }  // Point3D rotated

  double Pot;
  Vector3D localF;
  double xpt = localPt.X;
  double ypt = localPt.Y;
  double zpt = localPt.Z;

  double a = SurfLX;
  double b = SurfLZ;
  double diag = sqrt(a * a + b * b);  // diagonal of the element

  // distance of field point from element centroid
  double dist = sqrt(xpt * xpt + ypt * ypt + zpt * zpt);
  double dist3 = dist * dist * dist;

  if (dist >= FarField * diag)  // all are distances and, hence, +ve
  {
    double dA = a * b;
    Pot = dA / dist;
    localF.X = xpt * dA / dist3;
    localF.Y = ypt * dA / dist3;
    localF.Z = zpt * dA / dist3;
  } else {
    // normalize distances by `a' while sending - likely to improve accuracy
    int fstatus =
        ExactRecSurf(xpt / a, ypt / a, zpt / a, -1.0 / 2.0, -(b / a) / 2.0,
                     1.0 / 2.0, (b / a) / 2.0, &Pot, &localF);
    if (fstatus) {  // non-zero
      printf("problem in computing Potential of rectangular element ... \n");
      printf("a: %lg, b: %lg, X: %lg, Y: %lg, Z: %lg\n", a, b, xpt, ypt, zpt);
      // printf("returning ...\n");
      // return -1; void function at present
    }
    Pot *= a;  // rescale Potential - cannot be done outside because of if(dist)
  }

  (*globalF) = RotateVector3D(&localF, &DirCos, local2global);
  return (Pot);
}  // AreaKnChPF ends

// Pass vertices of the volume and get back approx potential and GCS flux
// at field point assuming unit charge density.
double ApproxVolumeKnChPF(int /*NbPts*/, Point3D* /*SourcePt*/, 
                          Point3D /*FieldPt*/,
                          Vector3D *globalF) {
  printf("ApproxVolumeKnChPF not implemented yet ... returning zero flux\n");
  globalF->X = 0.0;
  globalF->Y = 0.0;
  globalF->Z = 0.0;

  printf("ApproxVolumeKnChPF not implemented yet ... returning 0.0\n");
  return 0.0;
}  // ApproxVolumeKnChPF ends

// Pass vertices of the volume and get back potential and GCS flux
// at field point assuming unit charge density.
double VolumeKnChPF(int /*NbPts*/, Point3D* /*SourcePt*/, Point3D /*FieldPt*/,
                    Vector3D *globalF) {
  printf("VolumeKnChPF not implemented yet ... returning zero flux\n");
  globalF->X = 0.0;
  globalF->Y = 0.0;
  globalF->Z = 0.0;

  printf("VolumeKnChPF not implemented yet ... returning 0.0\n");
  return 0.0;
}  // VolumeKnChPF ends

// Return the sign of the argument or zero if the argument is smaller than
// a given value
int Sign(double x) {
  if (fabs(x) < MINDIST)
    return (0);
  else
    return (x < 0 ? -1 : 1);
}  // Sign ends

#ifdef __cplusplus
}  // namespace
#endif
