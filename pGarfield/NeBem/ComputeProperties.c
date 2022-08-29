/*
(c) 2005, Supratik Mukhopadhayay, Nayana Majumdar
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Isles.h"
#include "NR.h"
#include "Vector.h"
#include "neBEM.h"
#include "neBEMInterface.h"

#ifdef __cplusplus
namespace {
static constexpr double InvFourPiEps0 = 1. / MyFACTOR;
}

namespace neBEM {
#endif

// Weighting field function (WtFldPFAtPoint) has been modified!
// Later it be merged with PFAtPoint function.
// Check the notes written ahead of the weighting field function.

// Potential per unit charge density on an element
double GetPotential(int ele, Point3D *localP) {
  double value;

  switch ((EleArr + ele - 1)->G.Type) {
    case 4:  // rectangular element
      value = RecPot(ele, localP);
      break;
    case 3:  // triangular element
      value = TriPot(ele, localP);
      break;
    case 2:  // linear (wire) element
      value = WirePot(ele, localP);
      break;
    default:
      printf("Geometrical type out of range! ... exiting ...\n");
      exit(-1);
      break;  // never comes here
  }           // switch over gtsrc ends

  return (value);
}  // end of GetPotential

// Potential due to unit charge density on this element
double RecPot(int ele, Point3D *localP) {
  if (DebugLevel == 301) {
    printf("In RecPot ...\n");
  }

  double Pot;
  Vector3D Field;
  double xpt = localP->X;
  double ypt = localP->Y;
  double zpt = localP->Z;

  double a = (EleArr + ele - 1)->G.LX;
  double b = (EleArr + ele - 1)->G.LZ;
  double diag = sqrt(a * a + b * b);  // diagonal of the element

  // distance of field point from element centroid
  double dist = sqrt(xpt * xpt + ypt * ypt + zpt * zpt);

  if (dist >= FarField * diag) {
    double dA = a * b;
    Pot = dA / dist;
  } else {
    // normalize distances by `a' while sending - likely to improve accuracy
    int fstatus = ExactRecSurf(xpt / a, ypt / a, zpt / a, -0.5, -(b / a) / 2.0,
                               0.5, (b / a) / 2.0, &Pot, &Field);
    if (fstatus) { // non-zero
      printf("problem in computing Potential of rectangular element ... \n");
      printf("a: %lg, b: %lg, X: %lg, Y: %lg, Z: %lg\n", a, b, xpt, ypt, zpt);
      // printf("returning ...\n");
      // return -1; void function at present
    }
    Pot *= a;  // rescale Potential - cannot be done outside because of the `if'
  }

#ifdef __cplusplus
  return Pot * InvFourPiEps0;
#else
  return (Pot / MyFACTOR);
#endif
}  // end of RecPot

// Potential due to unit charge density on a triangular element
double TriPot(int ele, Point3D *localP) {
  if (DebugLevel == 301) {
    printf("In TriPot ...\n");
  }

  double Pot;
  Vector3D Field;
  double xpt = localP->X;
  double ypt = localP->Y;
  double zpt = localP->Z;

  // distance of field point from element centroid
  double a = (EleArr + ele - 1)->G.LX;
  double b = (EleArr + ele - 1)->G.LZ;
  // largest side (hypotenuse) of the element
  double diag = sqrt(a * a + b * b);

  const double xm = xpt - a / 3.;
  const double zm = zpt - b / 3.;
  double dist = sqrt(xm * xm + ypt * ypt + zm * zm);

  if (dist >= FarField * diag) {
    double dA = 0.5 * a * b;  // area of the triangular element
    Pot = dA / dist;
  } else {
    int fstatus = ExactTriSurf(b / a, xpt / a, ypt / a, zpt / a, &Pot, &Field);
    if (fstatus) {  // non-zero
      printf("problem in computing Potential of triangular element ... \n");
      printf("a: %lg, b: %lg, X: %lg, Y: %lg, Z: %lg\n", a, b, xpt, ypt, zpt);
      // printf("returning ...\n");
      // return -1; void function at present
    }
    Pot *= a;  // rescale Potential
  }

#ifdef __cplusplus
  return Pot * InvFourPiEps0;
#else
  return (Pot / MyFACTOR);
#endif
}  // end of TriPot

// Potential due to unit charge density on this element
// arguments not normalized yet
double WirePot(int ele, Point3D *localP) {
  if (DebugLevel == 301) {
    printf("In WirePot ...\n");
  }

  double Pot;
  double xpt = localP->X;
  double ypt = localP->Y;
  double zpt = localP->Z;

  double rW = (EleArr + ele - 1)->G.LX;
  double lW = (EleArr + ele - 1)->G.LZ;

  // field point from element centroid
  double dist = sqrt(xpt * xpt + ypt * ypt + zpt * zpt);

  if (dist >= FarField * lW) {
    double dA = 2.0 * ST_PI * rW * lW;
    // Pot = ApproxP_W(rW, lW, X, Y, Z, 1);
    Pot = dA / dist;
  } else if ((fabs(xpt) < MINDIST) && (fabs(ypt) < MINDIST) &&
             (fabs(zpt) < MINDIST)) {
    Pot = ExactCentroidalP_W(rW, lW);
  } else if ((fabs(xpt) < MINDIST) && (fabs(ypt) < MINDIST)) {
    Pot = ExactAxialP_W(rW, lW, localP->Z);
  } else {
    Pot = ExactThinP_W(rW, lW, xpt, ypt, zpt);
  }

#ifdef __cplusplus
  return Pot * InvFourPiEps0;
#else
  return (Pot / MyFACTOR);
#endif
}  // end of WirePot

// Flux per unit charge density on an element returned as globalF
// in the global coordiante system
void GetFluxGCS(int ele, Point3D *localP, Vector3D *globalF) {
  Vector3D localF;

  switch ((EleArr + ele - 1)->G.Type) {
    case 4:  // rectangular element
      RecFlux(ele, localP, &localF);
      break;
    case 3:  // triangular element
      TriFlux(ele, localP, &localF);
      break;
    case 2:  // linear (wire) element
      WireFlux(ele, localP, &localF);
      break;
    default:
      printf("Geometrical type out of range! ... exiting ...\n");
      exit(-1);
      break;  // never comes here
  }           // switch over gtsrc ends

  (*globalF) = RotateVector3D(&localF, &(EleArr + ele - 1)->G.DC, local2global);
}  // end of GetFluxGCS

// Flux per unit charge density on an element returned as localF
// in the local coordiante system
void GetFlux(int ele, Point3D *localP, Vector3D *localF) {
  switch ((EleArr + ele - 1)->G.Type) {
    case 4:  // rectangular element
      RecFlux(ele, localP, localF);
      break;
    case 3:  // triangular element
      TriFlux(ele, localP, localF);
      break;
    case 2:  // linear (wire) element
      WireFlux(ele, localP, localF);
      break;
    default:
      printf("Geometrical type out of range! ... exiting ...\n");
      exit(-1);
      break;  // never comes here
  }           // switch over gtsrc ends
}  // end of GetFlux

// local coord flux localF per unit charge density on one rectangular element
void RecFlux(int ele, Point3D *localP, Vector3D *localF) {
  if (DebugLevel == 301) {
    printf("In RecFlux ...\n");
  }

  double xpt = localP->X;
  double ypt = localP->Y;
  double zpt = localP->Z;

  double a = (EleArr + ele - 1)->G.LX;
  double b = (EleArr + ele - 1)->G.LZ;
  double diag = sqrt(a * a + b * b);  // diagonal of the element

  // distance of field point from element centroid
  double dist = sqrt(xpt * xpt + ypt * ypt + zpt * zpt);

  // no re-scaling necessary for `E'
  if (dist >= FarField * diag) {
    const double f = a * b / (dist * dist * dist);
    localF->X = xpt * f;
    localF->Y = ypt * f;
    localF->Z = zpt * f;
  } else {
    double Pot;
    int fstatus = ExactRecSurf(xpt / a, ypt / a, zpt / a, -0.5, -(b / a) / 2.0,
                               0.5, (b / a) / 2.0, &Pot, localF);
    if (fstatus) {  // non-zero
      printf("problem in computing flux of rectangular element ... \n");
      // printf("returning ...\n");
      // return -1; void function at present
    }
  }

#ifdef __cplusplus
  localF->X *= InvFourPiEps0;
  localF->Y *= InvFourPiEps0;
  localF->Z *= InvFourPiEps0;
#else
  localF->X /= MyFACTOR;
  localF->Y /= MyFACTOR;
  localF->Z /= MyFACTOR;
#endif
}  // end of RecFlux

// local coord flux per unit charge density on a triangluar element
void TriFlux(int ele, Point3D *localP, Vector3D *localF) {
  if (DebugLevel == 301) {
    printf("In TriFlux ...\n");
  }

  double xpt = localP->X;
  double ypt = localP->Y;
  double zpt = localP->Z;

  double a = (EleArr + ele - 1)->G.LX;
  double b = (EleArr + ele - 1)->G.LZ;
  double diag = sqrt(a * a + b * b);  // diagonal of the element

  // printf("In TriFlux\n");
  // printf("a: %lg, b: %lg, X: %lg, Y: %lg, Z: %lg\n", a, b, X, Y, Z);

  // distance of field point from element centroid
  const double xm = xpt - a / 3.;
  const double zm = zpt - b / 3.;
  double dist = sqrt(xm * xm + ypt * ypt + zm * zm);

  if (dist >= FarField * diag) {
    const double f = 0.5 * a * b / (dist * dist * dist);
    localF->X = xpt * f;
    localF->Y = ypt * f;
    localF->Z = zpt * f;
  } else {
    double Pot;
    int fstatus = ExactTriSurf(b / a, xpt / a, ypt / a, zpt / a, &Pot, localF);
    // fstatus = ApproxTriSurf(b/a, X/a, Y/a, Z/a, 5000, 5000, &Pot, &Flux);
    if (fstatus) {  // non-zero
      printf("problem in computing flux of triangular element ... \n");
      // printf("returning ...\n");
      // return -1; void function at present
    }
  }

#ifdef __cplusplus
  localF->X *= InvFourPiEps0;
  localF->Y *= InvFourPiEps0;
  localF->Z *= InvFourPiEps0;
#else
  localF->X /= MyFACTOR;
  localF->Y /= MyFACTOR;
  localF->Z /= MyFACTOR;
#endif
  // printf("Pot: %lg, Ex: %lg, Ey: %lg, Ez: %lg\n",
  // Pot, localF.X, localF.Y, localF.Z);
  // printf("Out of TriFlux\n");
}  // end of TriFlux

// local coord flux per unit charge density on a wire element
void WireFlux(int ele, Point3D *localP, Vector3D *localF) {
  if (DebugLevel == 301) {
    printf("In WireFlux ...\n");
  }

  double xpt = localP->X;
  double ypt = localP->Y;
  double zpt = localP->Z;
  double rW = (EleArr + ele - 1)->G.LX;
  double lW = (EleArr + ele - 1)->G.LZ;

  // field point from element centroid
  double dist = sqrt(xpt * xpt + ypt * ypt + zpt * zpt);

  if (dist >= FarField * lW) {
    const double f = 2.0 * ST_PI * rW * lW / (dist * dist * dist);
    localF->X = xpt * f;
    localF->Y = ypt * f;
    localF->Z = zpt * f;
  } else {
    if ((fabs(xpt) < MINDIST) && (fabs(ypt) < MINDIST)) {
      localF->X = localF->Y = 0.0;
    } else {
      localF->X = ExactThinFX_W(rW, lW, xpt, ypt, zpt);

      localF->Y = ExactThinFY_W(rW, lW, xpt, ypt, zpt);
    }

    // Ez
    localF->Z = ExactThinFZ_W(rW, lW, xpt, ypt, zpt);
  }

#ifdef __cplusplus
  localF->X *= InvFourPiEps0;
  localF->Y *= InvFourPiEps0;
  localF->Z *= InvFourPiEps0;
#else
  localF->X /= MyFACTOR;
  localF->Y /= MyFACTOR;
  localF->Z /= MyFACTOR;
#endif
}  // end of WireFlux

// Gives three components of the total Potential and flux in the global
// coordinate system due to all the interface elements and all known charges.
// It may be interesting to inspect the relative influence of these two factors.
int PFAtPoint(Point3D *globalP, double *Potential, Vector3D *globalF) {
  double ElePot;
  Vector3D EleglobalF;
  int fstatus = ElePFAtPoint(globalP, &ElePot, &EleglobalF);
  if (fstatus) {
    printf(
        "Problem in ElePFAtPoint being called from PFAtPoint ... returning\n");
    return (-1);
  }
  *Potential = ElePot;
  globalF->X = EleglobalF.X;
  globalF->Y = EleglobalF.Y;
  globalF->Z = EleglobalF.Z;

  if (OptKnCh) {
    double KnChPot;
    Vector3D KnChglobalF;
    fstatus = KnChPFAtPoint(globalP, &KnChPot, &KnChglobalF);
    if (fstatus) {
      printf(
          "Problem in KnChPFAtPoint being called from PFAtPoint ... "
          "returning\n");
      return (-1);
    }
    *Potential += KnChPot;
    globalF->X += KnChglobalF.X;
    globalF->Y += KnChglobalF.Y;
    globalF->Z += KnChglobalF.Z;
  }

  return 0;
}  // PFAtPoint ends

// Gives three components of the total Potential and flux in the global
// coordinate system only due to all the interface elements.
// Multi-threading implemented in the following routine
int ElePFAtPoint(Point3D *globalP, double *Potential, Vector3D *globalF) {
  int dbgFn = 0;

  const double xfld = globalP->X;
  const double yfld = globalP->Y;
  const double zfld = globalP->Z;

  // Compute Potential and field at different locations
  *Potential = globalF->X = globalF->Y = globalF->Z = 0.0;

  // Effects due to base primitives and their repetitions are considered in the
  // local coordinate system of the primitive (or element), while effects due to
  // mirror elements and their repetitions are considered in the global
  // coordinate system (GCS). This works because the direction cosines of a
  // primitive (and its elements) and those of its repetitions are the same.
  // As a result, we can do just one transformation from local to global at the
  // end of calculations related to a primitive. This can save substantial
  // computation if a discretized version of the primitive is being used since
  // we avoid one unnecessary transformation for each element that comprises a
  // primitive.
  // Begin with primitive description of the device

  // Scope in OpenMP: Variables in the global data space are accessible to all
  // threads, while variables in a thread's private space is accessible to the
  // thread only (there are several variations - copying outside region etc)
  // Field point remains the same - kept outside private
  // source point changes with change in primitive - private
  // TransformationMatrix changes - kept within private (Note: matrices with
  // fixed dimensions can be maintained, but those with dynamic allocation
  // can not).
  double *pPot = dvector(1, NbPrimitives);
  // Field components in LCS for a primitive and its other incarnations.
  double *plFx = dvector(1, NbPrimitives);
  double *plFy = dvector(1, NbPrimitives);
  double *plFz = dvector(1, NbPrimitives);

  for (int prim = 1; prim <= NbPrimitives; ++prim) {
    pPot[prim] = plFx[prim] = plFy[prim] = plFz[prim] = 0.0;
  }

#ifdef _OPENMP
  int tid = 0, nthreads = 1;
#pragma omp parallel private(tid, nthreads)
#endif
  {
#ifdef _OPENMP
    if (dbgFn) {
      tid = omp_get_thread_num();
      if (tid == 0) {
        nthreads = omp_get_num_threads();
        printf("PFAtPoint computation with %d threads\n", nthreads);
      }
    }
#endif
// by default, nested parallelization is off in C
#ifdef _OPENMP
#pragma omp for
#endif
    for (int primsrc = 1; primsrc <= NbPrimitives; ++primsrc) {
      if (dbgFn) {
        printf("Evaluating effect of primsrc %d using on %lg, %lg, %lg\n",
               primsrc, xfld, yfld, zfld);
        fflush(stdout);
      }

      const double xpsrc = PrimOriginX[primsrc];
      const double ypsrc = PrimOriginY[primsrc];
      const double zpsrc = PrimOriginZ[primsrc];

      // Field in the local frame.
      double lFx = 0.;
      double lFy = 0.;
      double lFz = 0.;

      // Set up transform matrix for this primitive, which is also the same
      // for all the elements belonging to this primitive.
      double TransformationMatrix[3][3];
      TransformationMatrix[0][0] = PrimDC[primsrc].XUnit.X;
      TransformationMatrix[0][1] = PrimDC[primsrc].XUnit.Y;
      TransformationMatrix[0][2] = PrimDC[primsrc].XUnit.Z;
      TransformationMatrix[1][0] = PrimDC[primsrc].YUnit.X;
      TransformationMatrix[1][1] = PrimDC[primsrc].YUnit.Y;
      TransformationMatrix[1][2] = PrimDC[primsrc].YUnit.Z;
      TransformationMatrix[2][0] = PrimDC[primsrc].ZUnit.X;
      TransformationMatrix[2][1] = PrimDC[primsrc].ZUnit.Y;
      TransformationMatrix[2][2] = PrimDC[primsrc].ZUnit.Z;

      // The total influence is due to primitives on the basic device and due to
      // virtual primitives arising out of repetition, reflection etc and not
      // residing on the basic device

      // Influence due to basic primitive

      // Evaluate possibility whether primitive influence is accurate enough.
      // This could be based on localPP and the subtended solid angle.
      // If 1, then only primitive influence will be considered.
      int PrimOK = 0;
      // consider primitive representation accurate enough if it is
      // repeated and beyond PrimAfter repetitions.
      if (PrimAfter < 0) { // If PrimAfter < 0, PrimOK is zero
        PrimOK = 0;
      } else if (PrimAfter == 0) {  // only this is necessary
        PrimOK = 1;
      } else if (PrimAfter > 0) {
        PrimOK = 1;
      }
      if (PrimOK) {
        // Only primitive influence will be considered
        // Potential and flux (local system) due to base primitive
        double tmpPot;
        Vector3D tmpF;
        // Rotate point from global to local system
        double InitialVector[3] = {xfld - xpsrc, yfld - ypsrc, zfld - zpsrc};
        double FinalVector[3] = {0., 0., 0.};
        for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
            FinalVector[i] += TransformationMatrix[i][j] * InitialVector[j];
          }
        }
        Point3D localPP;
        localPP.X = FinalVector[0];
        localPP.Y = FinalVector[1];
        localPP.Z = FinalVector[2];
        GetPrimPF(primsrc, &localPP, &tmpPot, &tmpF);
        const double qpr = AvChDen[primsrc] + AvAsgndChDen[primsrc];
        pPot[primsrc] += qpr * tmpPot;
        lFx += qpr * tmpF.X;
        lFy += qpr * tmpF.Y;
        lFz += qpr * tmpF.Z;
        // if(DebugLevel == 301)
        if (dbgFn) {
          printf("PFAtPoint base primitive =>\n");
          printf("primsrc: %d, xlocal: %lg, ylocal: %lg, zlocal %lg\n", primsrc,
                 localPP.X, localPP.Y, localPP.Z);
          printf("primsrc: %d, Pot: %lg, Fx: %lg, Fx: %lg, Fz: %lg\n", primsrc,
                 tmpPot, tmpF.X, tmpF.Y, tmpF.Z);
          printf("primsrc: %d, pPot: %lg, lFx: %lg, lFy: %lg, lFz: %lg\n",
                 primsrc, pPot[primsrc], lFx, lFy, lFz);
          fflush(stdout);
          // exit(-1);
        }
      } else {
        // Need to consider element influence.
        double tPot;
        Vector3D tF;
        double ePot = 0.;
        Vector3D eF;
        eF.X = 0.0;
        eF.Y = 0.0;
        eF.Z = 0.0;
        const int eleMin = ElementBgn[primsrc];
        const int eleMax = ElementEnd[primsrc];
        for (int ele = eleMin; ele <= eleMax; ++ele) {
          const double xsrc = (EleArr + ele - 1)->G.Origin.X;
          const double ysrc = (EleArr + ele - 1)->G.Origin.Y;
          const double zsrc = (EleArr + ele - 1)->G.Origin.Z;
          // Rotate from global to local system; matrix as for primitive
          double vG[3] = {xfld - xsrc, yfld - ysrc, zfld - zsrc};
          double vL[3] = {0., 0., 0.};
          for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
              vL[i] += TransformationMatrix[i][j] * vG[j];
            }
          }
          // Potential and flux (local system) due to base primitive
          const int type = (EleArr + ele - 1)->G.Type;
          const double a = (EleArr + ele - 1)->G.LX;
          const double b = (EleArr + ele - 1)->G.LZ;
          GetPF(type, a, b, vL[0], vL[1], vL[2], &tPot, &tF);
          const double qel =
              (EleArr + ele - 1)->Solution + (EleArr + ele - 1)->Assigned;
          ePot += qel * tPot;
          eF.X += qel * tF.X;
          eF.Y += qel * tF.Y;
          eF.Z += qel * tF.Z;
          // if(DebugLevel == 301)
          if (dbgFn) {
            printf("PFAtPoint base primitive:%d\n", primsrc);
            printf("ele: %d, xlocal: %lg, ylocal: %lg, zlocal %lg\n", ele,
                   vL[0], vL[1], vL[2]);
            printf(
                "ele: %d, tPot: %lg, tFx: %lg, tFy: %lg, tFz: %lg, Solution: "
                "%g\n",
                ele, tPot, tF.X, tF.Y, tF.Z, qel);
            printf("ele: %d, ePot: %lg, eFx: %lg, eFy: %lg, eFz: %lg\n", ele,
                   ePot, eF.X, eF.Y, eF.Z);
            fflush(stdout);
          }
        }  // for all the elements on this primsrc primitive

        pPot[primsrc] += ePot;
        lFx += eF.X;
        lFy += eF.Y;
        lFz += eF.Z;
        if (dbgFn) {
          printf("prim%d, ePot: %lg, eFx: %lg, eFy: %lg, eFz: %lg\n", primsrc,
                 ePot, eF.X, eF.Y, eF.Z);
          printf("prim%d, pPot: %lg, lFx: %lg, lFy: %lg, lFz: %lg\n", primsrc,
                 pPot[primsrc], lFx, lFy, lFz);
          fflush(stdout);
        }
      }  // else elements influence

      // if(DebugLevel == 301)
      if (dbgFn) {
        printf("basic primitive\n");
        printf("primsrc: %d, pPot: %lg, lFx: %lg, lFy: %lg, lFz: %lg\n",
               primsrc, pPot[primsrc], lFx, lFy, lFz);
        fflush(stdout);
      }

      if (MirrorTypeX[primsrc] || MirrorTypeY[primsrc] ||
          MirrorTypeZ[primsrc]) {  // Mirror effect of base primitives
        printf("Mirror may not be correctly implemented ...\n");
        exit(0);
      }  // Mirror effect ends

      // Flux due to repeated primitives
      if ((PeriodicTypeX[primsrc] == 1) || (PeriodicTypeY[primsrc] == 1) ||
          (PeriodicTypeZ[primsrc] == 1)) {
        const int perx = PeriodicInX[primsrc];
        const int pery = PeriodicInY[primsrc];
        const int perz = PeriodicInZ[primsrc];
        if (perx || pery || perz) {
          for (int xrpt = -perx; xrpt <= perx; ++xrpt) {
            const double xShift = XPeriod[primsrc] * (double)xrpt;
            const double XPOfRpt = xpsrc + xShift;
            for (int yrpt = -pery; yrpt <= pery; ++yrpt) {
              const double yShift = YPeriod[primsrc] * (double)yrpt;
              const double YPOfRpt = ypsrc + yShift;
              for (int zrpt = -perz; zrpt <= perz; ++zrpt) {
                const double zShift = ZPeriod[primsrc] * (double)zrpt;
                const double ZPOfRpt = zpsrc + zShift;
                // Skip the basic primitive.
                if ((xrpt == 0) && (yrpt == 0) && (zrpt == 0)) continue;

                // Basic primitive repeated
                int repPrimOK = 0;
                // consider primitive representation accurate enough if it is
                // repeated and beyond PrimAfter repetitions.
                if (PrimAfter < 0) { // If PrimAfter <0, repPrimOK is zero
                  repPrimOK = 0;
                } else if ((abs(xrpt) >= PrimAfter)
                            && (abs(yrpt) >= PrimAfter)) {
                  repPrimOK = 1;
                }
                if (repPrimOK) {
                  // Use primitive representation
                  // Rotate point from global to local system
                  double InitialVector[3] = {xfld - XPOfRpt, yfld - YPOfRpt,
                                             zfld - ZPOfRpt};
                  double FinalVector[3] = {0., 0., 0.};
                  for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                      FinalVector[i] +=
                          TransformationMatrix[i][j] * InitialVector[j];
                    }
                  }
                  Point3D localPPR;
                  localPPR.X = FinalVector[0];
                  localPPR.Y = FinalVector[1];
                  localPPR.Z = FinalVector[2];
                  // Potential and flux (local system) due to repeated
                  // primitive
                  double tmpPot;
                  Vector3D tmpF;
                  GetPrimPF(primsrc, &localPPR, &tmpPot, &tmpF);
                  const double qpr = AvChDen[primsrc] + AvAsgndChDen[primsrc];
                  pPot[primsrc] += qpr * tmpPot;
                  lFx += qpr * tmpF.X;
                  lFy += qpr * tmpF.Y;
                  lFz += qpr * tmpF.Z;
                  // if(DebugLevel == 301)
                  if (dbgFn) {
                    printf(
                        "primsrc: %d, xlocal: %lg, ylocal: %lg, zlocal: "
                        "%lg\n",
                        primsrc, localPPR.X, localPPR.Y, localPPR.Z);
                    printf("primsrc: %d, Pot: %lg, Fx: %lg, Fy: %lg, Fz: %lg\n",
                           primsrc, tmpPot * qpr, tmpF.X * qpr, tmpF.Y * qpr,
                           tmpF.Z * qpr);
                    printf(
                        "primsrc: %d, pPot: %lg, lFx: %lg, lFy: %lg, lFz: "
                        "%lg\n",
                        primsrc, pPot[primsrc], lFx, lFy, lFz);
                    fflush(stdout);
                  }
                } else {
                  // Use discretized representation of a repeated primitive
                  double tPot;
                  Vector3D tF;
                  double erPot = 0.0;
                  Vector3D erF;
                  erF.X = 0.0;
                  erF.Y = 0.0;
                  erF.Z = 0.0;
                  const int eleMin = ElementBgn[primsrc];
                  const int eleMax = ElementEnd[primsrc];
                  for (int ele = eleMin; ele <= eleMax; ++ele) {
                    const double xrsrc = (EleArr + ele - 1)->G.Origin.X;
                    const double yrsrc = (EleArr + ele - 1)->G.Origin.Y;
                    const double zrsrc = (EleArr + ele - 1)->G.Origin.Z;

                    const double XEOfRpt = xrsrc + xShift;
                    const double YEOfRpt = yrsrc + yShift;
                    const double ZEOfRpt = zrsrc + zShift;
                    // Rotate from global to local system
                    double vG[3] = {xfld - XEOfRpt, yfld - YEOfRpt,
                                    zfld - ZEOfRpt};
                    double vL[3] = {0., 0., 0.};
                    for (int i = 0; i < 3; ++i) {
                      for (int j = 0; j < 3; ++j) {
                        vL[i] += TransformationMatrix[i][j] * vG[j];
                      }
                    }
                    // Allowed, because all the local coordinates have the
                    // same orientations. Only the origins are mutually
                    // displaced along a line.
                    const int type = (EleArr + ele - 1)->G.Type;
                    const double a = (EleArr + ele - 1)->G.LX;
                    const double b = (EleArr + ele - 1)->G.LZ;
                    GetPF(type, a, b, vL[0], vL[1], vL[2], &tPot, &tF);
                    const double qel = (EleArr + ele - 1)->Solution +
                                       (EleArr + ele - 1)->Assigned;
                    erPot += qel * tPot;
                    erF.X += qel * tF.X;
                    erF.Y += qel * tF.Y;
                    erF.Z += qel * tF.Z;
                    // if(DebugLevel == 301)
                    if (dbgFn) {
                      printf("PFAtPoint base primitive:%d\n", primsrc);
                      printf("ele: %d, xlocal: %lg, ylocal: %lg, zlocal %lg\n",
                             ele, vL[0], vL[1], vL[2]);
                      printf(
                          "ele: %d, tPot: %lg, tFx: %lg, tFy: %lg, tFz: %lg, "
                          "Solution: %g\n",
                          ele, tPot, tF.X, tF.Y, tF.Z, qel);
                      printf(
                          "ele: %d, ePot: %lg, eFx: %lg, eFy: %lg, eFz: %lg\n",
                          ele, erPot, erF.X, erF.Y, erF.Z);
                      fflush(stdout);
                    }
                  }  // for all the elements on this primsrc repeated
                     // primitive

                  pPot[primsrc] += erPot;
                  lFx += erF.X;
                  lFy += erF.Y;
                  lFz += erF.Z;
                }  // else discretized representation of this primitive

                // if(DebugLevel == 301)
                if (dbgFn) {
                  printf("basic repeated xrpt: %d. yrpt: %d, zrpt: %d\n", xrpt,
                         yrpt, zrpt);
                  printf(
                      "primsrc: %d, pPot: %lg, lFx: %lg, lFy: %lg, lFz: %lg\n",
                      primsrc, pPot[primsrc], lFx, lFy, lFz);
                  fflush(stdout);
                }

                if (MirrorTypeX[primsrc] || MirrorTypeY[primsrc] ||
                    MirrorTypeZ[primsrc]) {  // Mirror effect of repeated
                                             // primitives - not parallelized
                  printf(
                      "Mirror not correctly implemented in this version of "
                      "neBEM ...\n");
                  exit(0);

                  double tmpPot;
                  Vector3D tmpF;
                  Point3D srcptp;
                  Point3D localPPRM;  // point primitive repeated mirrored
                  DirnCosn3D DirCos;

                  srcptp.X = XPOfRpt;
                  srcptp.Y = YPOfRpt;
                  srcptp.Z = ZPOfRpt;

                  Point3D fldpt;
                  fldpt.X = xfld;
                  fldpt.Y = yfld;
                  fldpt.Z = zfld;
                  if (MirrorTypeX[primsrc]) {
                    MirrorTypeY[primsrc] = 0;
                    MirrorTypeZ[primsrc] = 0;
                  }
                  if (MirrorTypeY[primsrc]) MirrorTypeZ[primsrc] = 0;

                  if (MirrorTypeX[primsrc]) {
                    localPPRM = ReflectPrimitiveOnMirror(
                        'X', primsrc, srcptp, fldpt,
                        MirrorDistXFromOrigin[primsrc], &DirCos);

                    // check whether primitive description is good enough
                    int mirrPrimOK = 0;
                    if (mirrPrimOK) {
                      GetPrimPFGCS(primsrc, &localPPRM, &tmpPot, &tmpF,
                                   &DirCos);
                      const double qpr =
                          AvChDen[primsrc] + AvAsgndChDen[primsrc];
                      if (MirrorTypeX[primsrc] == 1) {
                        // opposite charge density
                        pPot[primsrc] -= qpr * tmpPot;
                        lFx -= qpr * tmpF.X;
                        lFy -= qpr * tmpF.Y;
                        lFz -= qpr * tmpF.Z;
                      } else if (MirrorTypeX[primsrc] == 2) {
                        // same charge density
                        pPot[primsrc] += qpr * tmpPot;
                        lFx += qpr * tmpF.X;
                        lFy += qpr * tmpF.Y;
                        lFz += qpr * tmpF.Z;
                      }
                    } else {              // consider element representation
                      Point3D localPERM;  // point element repeated mirrored
                      Point3D srcpte;

                      const int eleMin = ElementBgn[primsrc];
                      const int eleMax = ElementEnd[primsrc];
                      for (int ele = eleMin; ele <= eleMax; ++ele) {
                        const double xsrc = (EleArr + ele - 1)->G.Origin.X;
                        const double ysrc = (EleArr + ele - 1)->G.Origin.Y;
                        const double zsrc = (EleArr + ele - 1)->G.Origin.Z;

                        const double XEOfRpt = xsrc + xShift;
                        const double YEOfRpt = ysrc + yShift;
                        const double ZEOfRpt = zsrc + zShift;

                        srcpte.X = XEOfRpt;
                        srcpte.Y = YEOfRpt;
                        srcpte.Z = ZEOfRpt;

                        localPERM = ReflectOnMirror(
                            'X', ele, srcpte, fldpt,
                            MirrorDistXFromOrigin[primsrc], &DirCos);
                        const int type = (EleArr + ele - 1)->G.Type;
                        const double a = (EleArr + ele - 1)->G.LX;
                        const double b = (EleArr + ele - 1)->G.LZ;
                        GetPFGCS(type, a, b, &localPERM, &tmpPot, &tmpF,
                                 &DirCos);  // force?
                        const double qel = (EleArr + ele - 1)->Solution +
                                           (EleArr + ele - 1)->Assigned;
                        if (MirrorTypeX[primsrc] == 1) {
                          // opposite charge density
                          pPot[primsrc] -= qel * tmpPot;
                          lFx -= qel * tmpF.X;
                          lFy -= qel * tmpF.Y;
                          lFz -= qel * tmpF.Z;
                        } else if (MirrorTypeX[primsrc] == 2) {
                          // same charge density
                          pPot[primsrc] += qel * tmpPot;
                          lFx += qel * tmpF.X;
                          lFy += qel * tmpF.Y;
                          lFz += qel * tmpF.Z;
                        }
                      }  // loop for all elements on the primsrc primitive
                    }    // else element representation
                  }      // MirrorTypeX

                  if (MirrorTypeY[primsrc]) {
                    localPPRM = ReflectOnMirror('Y', primsrc, srcptp, fldpt,
                                                MirrorDistYFromOrigin[primsrc],
                                                &DirCos);

                    // check whether primitive description is good enough
                    int mirrPrimOK = 0;
                    if (mirrPrimOK) {
                      GetPrimPFGCS(primsrc, &localPPRM, &tmpPot, &tmpF,
                                   &DirCos);
                      const double qpr =
                          AvChDen[primsrc] + AvAsgndChDen[primsrc];
                      if (MirrorTypeY[primsrc] == 1) {
                        // opposite charge density
                        pPot[primsrc] -= qpr * tmpPot;
                        lFx -= qpr * tmpF.X;
                        lFy -= qpr * tmpF.Y;
                        lFz -= qpr * tmpF.Z;
                      } else if (MirrorTypeY[primsrc] == 2) {
                        // same charge density
                        pPot[primsrc] += qpr * tmpPot;
                        lFx += qpr * tmpF.X;
                        lFy += qpr * tmpF.Y;
                        lFz += qpr * tmpF.Z;
                      }
                    } else {  // consider element representation
                      Point3D localPERM;
                      Point3D srcpte;

                      const int eleMin = ElementBgn[primsrc];
                      const int eleMax = ElementEnd[primsrc];
                      for (int ele = eleMin; ele <= eleMax; ++ele) {
                        const double xsrc = (EleArr + ele - 1)->G.Origin.X;
                        const double ysrc = (EleArr + ele - 1)->G.Origin.Y;
                        const double zsrc = (EleArr + ele - 1)->G.Origin.Z;

                        const double XEOfRpt = xsrc + xShift;
                        const double YEOfRpt = ysrc + yShift;
                        const double ZEOfRpt = zsrc + zShift;

                        srcpte.X = XEOfRpt;
                        srcpte.Y = YEOfRpt;
                        srcpte.Z = ZEOfRpt;

                        localPERM = ReflectOnMirror(
                            'Y', ele, srcpte, fldpt,
                            MirrorDistYFromOrigin[primsrc], &DirCos);
                        const int type = (EleArr + ele - 1)->G.Type;
                        const double a = (EleArr + ele - 1)->G.LX;
                        const double b = (EleArr + ele - 1)->G.LZ;
                        GetPFGCS(type, a, b, &localPERM, &tmpPot, &tmpF,
                                 &DirCos);
                        const double qel = (EleArr + ele - 1)->Solution +
                                           (EleArr + ele - 1)->Assigned;
                        if (MirrorTypeY[primsrc] == 1) {
                          // opposite charge density
                          pPot[primsrc] -= qel * tmpPot;
                          lFx -= qel * tmpF.X;
                          lFy -= qel * tmpF.Y;
                          lFz -= qel * tmpF.Z;
                        } else if (MirrorTypeY[primsrc] == 2) {
                          // same charge density
                          pPot[primsrc] += qel * tmpPot;
                          lFx += qel * tmpF.X;
                          lFy += qel * tmpF.Y;
                          lFz += qel * tmpF.Z;
                        }
                      }  // loop for all elements on the primsrc primitive
                    }    // else element representations
                  }      // MirrorTypeY

                  if (MirrorTypeZ[primsrc]) {
                    localPPRM = ReflectOnMirror('Z', primsrc, srcptp, fldpt,
                                                MirrorDistZFromOrigin[primsrc],
                                                &DirCos);

                    // check whether primitive description is good enough
                    int mirrPrimOK = 0;
                    if (mirrPrimOK) {
                      GetPrimPFGCS(primsrc, &localPPRM, &tmpPot, &tmpF,
                                   &DirCos);
                      const double qpr =
                          AvChDen[primsrc] + AvAsgndChDen[primsrc];
                      if (MirrorTypeZ[primsrc] == 1) {
                        // opposite charge density
                        pPot[primsrc] -= qpr * tmpPot;
                        lFx -= qpr * tmpF.X;
                        lFy -= qpr * tmpF.Y;
                        lFz -= qpr * tmpF.Z;
                      } else if (MirrorTypeZ[primsrc] == 2) {
                        // same charge density
                        pPot[primsrc] += qpr * tmpPot;
                        lFx += qpr * tmpF.X;
                        lFy += qpr * tmpF.Y;
                        lFz += qpr * tmpF.Z;
                      }
                    } else {
                      // elements to be considered
                      Point3D localPERM;
                      Point3D srcpte;

                      const int eleMin = ElementBgn[primsrc];
                      const int eleMax = ElementEnd[primsrc];
                      for (int ele = eleMin; ele <= eleMax; ++ele) {
                        const double xsrc = (EleArr + ele - 1)->G.Origin.X;
                        const double ysrc = (EleArr + ele - 1)->G.Origin.Y;
                        const double zsrc = (EleArr + ele - 1)->G.Origin.Z;

                        const double XEOfRpt = xsrc + xShift;
                        const double YEOfRpt = ysrc + yShift;
                        const double ZEOfRpt = zsrc + zShift;

                        srcpte.X = XEOfRpt;
                        srcpte.Y = YEOfRpt;
                        srcpte.Z = ZEOfRpt;

                        localPERM = ReflectOnMirror(
                            'Z', ele, srcpte, fldpt,
                            MirrorDistZFromOrigin[primsrc], &DirCos);
                        const int type = (EleArr + ele - 1)->G.Type;
                        const double a = (EleArr + ele - 1)->G.LX;
                        const double b = (EleArr + ele - 1)->G.LZ;
                        GetPFGCS(type, a, b, &localPERM, &tmpPot, &tmpF,
                                 &DirCos);
                        const double qel = (EleArr + ele - 1)->Solution +
                                           (EleArr + ele - 1)->Assigned;
                        if (MirrorTypeZ[primsrc] == 1) {
                          // opposite charge density
                          pPot[primsrc] -= qel * tmpPot;
                          lFx -= qel * tmpF.X;
                          lFy -= qel * tmpF.Y;
                          lFz -= qel * tmpF.Z;
                        } else if (MirrorTypeZ[primsrc] == 2) {
                          // same charge density
                          pPot[primsrc] += qel * tmpPot;
                          lFx += qel * tmpF.X;
                          lFy += qel * tmpF.Y;
                          lFz += qel * tmpF.Z;
                        }
                      }  // loop for all elements on the primsrc primitive
                    }    // else consider element representation
                  }      // MirrorTypeZ
                }        // Mirror effect for repeated primitives ends

              }  // for zrpt
            }    // for yrpt
          }      // for xrpt
        }        // PeriodicInX || PeriodicInY || PeriodicInZ
      }          // PeriodicType == 1
      Vector3D localF;
      localF.X = lFx;
      localF.Y = lFy;
      localF.Z = lFz;
      Vector3D tmpF = RotateVector3D(&localF, &PrimDC[primsrc], local2global);
      plFx[primsrc] = tmpF.X;
      plFy[primsrc] = tmpF.Y;
      plFz[primsrc] = tmpF.Z;
    }  // for all primitives: basic device, mirror reflections and repetitions
  }    // pragma omp parallel

  double totPot = 0.0;
  Vector3D totF;
  totF.X = totF.Y = totF.Z = 0.0;
  for (int prim = 1; prim <= NbPrimitives; ++prim) {
    totPot += pPot[prim];
    totF.X += plFx[prim];
    totF.Y += plFy[prim];
    totF.Z += plFz[prim];
  }

  // This is done at the end of the function - before freeing memory
#ifdef __cplusplus
  *Potential = totPot * InvFourPiEps0;
  globalF->X = totF.X * InvFourPiEps0;
  globalF->Y = totF.Y * InvFourPiEps0;
  globalF->Z = totF.Z * InvFourPiEps0;
#else
  *Potential = totPot / MyFACTOR;
  globalF->X = totF.X / MyFACTOR;
  globalF->Y = totF.Y / MyFACTOR;
  globalF->Z = totF.Z / MyFACTOR;
#endif
  (*Potential) += VSystemChargeZero;  // respect total system charge constraint

  if (dbgFn) {
    printf("Final values due to all primitives: ");
    // printf("xfld\tyfld\tzfld\tPot\tFx\tFy\tFz\n");	// refer, do not
    // uncomment
    printf("%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n\n", xfld, yfld, zfld,
           (*Potential), globalF->X, globalF->Y, globalF->Z);
    fflush(stdout);
  }

  free_dvector(pPot, 1, NbPrimitives);
  free_dvector(plFx, 1, NbPrimitives);
  free_dvector(plFy, 1, NbPrimitives);
  free_dvector(plFz, 1, NbPrimitives);

  return (0);
}  // end of ElePFAtPoint

// At present implemented without OpenMP
// Evaluate effects due to known charge distributions
// Since there is no intermediate function that interfaces PointKnChPF etc,
// division by MyFACTOR is necessary - carried out just before end of function.
// CHECK OpenMP / GPU possibilities:
// Do parallelize before using these known charges - points or distributions
int KnChPFAtPoint(Point3D *globalP, double *Potential, Vector3D *globalF) {
  int dbgFn = 0;
  double tmpPot;
  Point3D fldPt;
  fldPt.X = globalP->X;
  fldPt.Y = globalP->Y;
  fldPt.Z = globalP->Z;
  Vector3D tmpF;

  // initialize the values for each point being evaluated.
  *Potential = 0.0;
  globalF->X = globalF->Y = globalF->Z = 0.0;

  for (int point = 1; point <= NbPointsKnCh; ++point) {
    Point3D srcPt;
    srcPt.X = PointKnChArr[point].P.X;
    srcPt.Y = PointKnChArr[point].P.Y;
    srcPt.Z = PointKnChArr[point].P.Z;

    tmpPot = PointKnChPF(srcPt, fldPt, &tmpF);
    (*Potential) += PointKnChArr[point].Assigned * tmpPot;
    globalF->X += PointKnChArr[point].Assigned * tmpF.X;
    globalF->Y += PointKnChArr[point].Assigned * tmpF.Y;
    globalF->Z += PointKnChArr[point].Assigned * tmpF.Z;
  }

  if (dbgFn) {
    printf("Final values due to all known point charges (*MyFACTOR): ");
    // printf("xfld\tyfld\tzfld\tPot\tFx\tFy\tFz\n"); // refer, do not uncomment
    printf("%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n\n", fldPt.X, fldPt.Y, fldPt.Z,
           (*Potential), globalF->X, globalF->Y, globalF->Z);
    fflush(stdout);
  }

  for (int line = 1; line <= NbLinesKnCh; ++line) {
    Point3D startPt, stopPt;
    startPt.X = LineKnChArr[line].Start.X;
    startPt.Y = LineKnChArr[line].Start.Y;
    startPt.Z = LineKnChArr[line].Start.Z;
    stopPt.X = LineKnChArr[line].Stop.X;
    stopPt.Y = LineKnChArr[line].Stop.Y;
    stopPt.Z = LineKnChArr[line].Stop.Z;
    tmpPot = LineKnChPF(startPt, stopPt, fldPt, &tmpF);
    (*Potential) += LineKnChArr[line].Assigned * tmpPot;
    globalF->X += LineKnChArr[line].Assigned * tmpF.X;
    globalF->Y += LineKnChArr[line].Assigned * tmpF.Y;
    globalF->Z += LineKnChArr[line].Assigned * tmpF.Z;
  }

  if (dbgFn) {
    printf("Final values due to all known line charges (*MyFACTOR): ");
    // printf("xfld\tyfld\tzfld\tPot\tFx\tFy\tFz\n"); // refer, do not uncomment
    printf("%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n\n", fldPt.X, fldPt.Y, fldPt.Z,
           (*Potential), globalF->X, globalF->Y, globalF->Z);
    fflush(stdout);
  }

  // Conversion to simpler structure may be necessary, as done for points and
  // lines, above.
  for (int area = 1; area <= NbAreasKnCh; ++area) {
    tmpPot = AreaKnChPF((AreaKnChArr + area - 1)->NbVertices,
                        ((AreaKnChArr + area - 1)->Vertex), fldPt, &tmpF);
    (*Potential) += (AreaKnChArr + area - 1)->Assigned * tmpPot;
    globalF->X += (AreaKnChArr + area - 1)->Assigned * tmpF.X;
    globalF->Y += (AreaKnChArr + area - 1)->Assigned * tmpF.Y;
    globalF->Z += (AreaKnChArr + area - 1)->Assigned * tmpF.Z;
  }

  if (dbgFn) {
    printf("Final values due to all known area charges (*MyFACTOR): ");
    // printf("xfld\tyfld\tzfld\tPot\tFx\tFy\tFz\n"); // refer, do not uncomment
    printf("%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n\n", fldPt.X, fldPt.Y, fldPt.Z,
           (*Potential), globalF->X, globalF->Y, globalF->Z);
    fflush(stdout);
  }

  for (int vol = 1; vol <= NbVolumesKnCh; ++vol) {
    tmpPot = VolumeKnChPF((VolumeKnChArr + vol - 1)->NbVertices,
                          ((VolumeKnChArr + vol - 1)->Vertex), fldPt, &tmpF);
    (*Potential) += (VolumeKnChArr + vol - 1)->Assigned * tmpPot;
    globalF->X += (VolumeKnChArr + vol - 1)->Assigned * tmpF.X;
    globalF->Y += (VolumeKnChArr + vol - 1)->Assigned * tmpF.Y;
    globalF->Z += (VolumeKnChArr + vol - 1)->Assigned * tmpF.Z;
  }

  if (dbgFn) {
    printf("Final values due to all known volume charges (*MyFACTOR): ");
    // printf("xfld\tyfld\tzfld\tPot\tFx\tFy\tFz\n"); // for reference, do not
    // uncomment
    printf("%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n\n", fldPt.X, fldPt.Y, fldPt.Z,
           (*Potential), globalF->X, globalF->Y, globalF->Z);
    fflush(stdout);
  }

  (*Potential) /= MyFACTOR;
  globalF->X /= MyFACTOR;
  globalF->Y /= MyFACTOR;
  globalF->Z /= MyFACTOR;

  return 0;
}  // KnChPFAtPoint ends

// Compute voxelized data for export to Garfield++
int VoxelFPR(void) {
  int dbgFn = 0;

  printf(
      "\nPhysical potential and field computation for voxelized data export\n");

  char VoxelFile[256];
  strcpy(VoxelFile, BCOutDir);
  strcat(VoxelFile, "/VoxelFPR.out");
  FILE *fVoxel = fopen(VoxelFile, "w");
  if (fVoxel == NULL) {
    neBEMMessage("VoxelFPR - VoxelFile");
    return -1;
  }
  fprintf(
      fVoxel,
      "# X(cm)\tY(cm)\tZ(cm)\tFX(V/cm)\tFY(V/cm)\tFZ(V/cm)\tPot(V)\tRegion\n");

  if (dbgFn) {
    printf("VoxelFPR.out created ...\n");
    fflush(stdout);
  }

  int nbXCells = Voxel.NbXCells;
  int nbYCells = Voxel.NbYCells;
  int nbZCells = Voxel.NbZCells;
  double startX = Voxel.Xmin;
  double startY = Voxel.Ymin;
  double startZ = Voxel.Zmin;
  double delX = (Voxel.Xmax - Voxel.Xmin) / nbXCells;
  double delY = (Voxel.Ymax - Voxel.Ymin) / nbYCells;
  double delZ = (Voxel.Zmax - Voxel.Zmin) / nbZCells;

  int ivol;  // relates XYZ position to volume number
  double *VoxelFX = dvector(0, nbZCells + 1);
  double *VoxelFY = dvector(0, nbZCells + 1);
  double *VoxelFZ = dvector(0, nbZCells + 1);
  double *VoxelP = dvector(0, nbZCells + 1);

  if (dbgFn) {
    printf("nbXCells, nbYCells, nbZCells: %d, %d, %d\n", nbXCells, nbYCells,
           nbZCells);
    printf("startX, startY, startZ: %le, %le, %le\n", startX, startY, startZ);
    printf("delX, delY, delZ: %le, %le, %le\n", delX, delY, delZ);
    fflush(stdout);
  }

  for (int i = 1; i <= nbXCells + 1; ++i) {
    for (int j = 1; j <= nbYCells + 1; ++j) {
      if (dbgFn) {
        printf("VoxelFPR => i: %4d, j: %4d", i, j);
        fflush(stdout);
      }

      Point3D point;
#ifdef _OPENMP
      int nthreads = 1, tid = 0;
#pragma omp parallel private(nthreads, tid)
#endif
      {
#ifdef _OPENMP
        if (dbgFn) {
          tid = omp_get_thread_num();
          if (tid == 0) {
            nthreads = omp_get_num_threads();
            printf("Starting voxel computation with %d threads\n", nthreads);
          }
        }
#endif
        int k;
        Vector3D field;
        double potential;
#ifdef _OPENMP
#pragma omp for private(k, point, potential, field)
#endif
        for (k = 1; k <= nbZCells + 1; ++k) {
          potential = 0.0;
          field.X = field.Y = field.Z = 0.0;

          point.X = startX + (i - 1) * delX;  // all 3 components need to be
          point.Y = startY + (j - 1) * delY;  // evaluated after pragma omp
          point.Z = startZ + (k - 1) * delZ;

          if (dbgFn) {
            printf("i, j, k: %d, %d, %d\n", i, j, k);
            printf("point X, Y, Z: %.8lg\t%.8lg\t%.8lg\n",
                   point.X / LengthScale, point.Y / LengthScale,
                   point.Z / LengthScale);
            fflush(stdout);
          }

          int fstatus = PFAtPoint(&point, &potential, &field);
          if (fstatus != 0) {
            neBEMMessage("wrong PFAtPoint return value in VoxelFPR\n");
            // return -1;
          }

          if (dbgFn) {
            printf("%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\n",
                   point.X / LengthScale, point.Y / LengthScale,
                   point.Z / LengthScale, field.X, field.Y, field.Z,
                   potential / LengthScale);
            fflush(stdout);
          }

          VoxelFX[k] = field.X;
          VoxelFY[k] = field.Y;
          VoxelFZ[k] = field.Z;
          VoxelP[k] = potential;
        }  // loop k
      }    // pragma omp parallel

      for (int k = 1; k <= nbZCells + 1; ++k) {  // file output
        point.X = startX + (i - 1) * delX;
        point.Y = startY + (j - 1) * delY;
        point.Z = startZ + (k - 1) * delZ;

        ivol = neBEMVolumePoint(point.X, point.Y, point.Z);
        /*
        volMaterial[ivol];	// region linked to material
        neBEMVolumeDescription(ivol, &vshp, &vmat, &veps, &vpot, &vq, &vtype);
        if (dbgFn) {
          printf("volref: %d\n", ivol);
          printf("shape: %d,  material: %d\n", volShape[ivol], volMaterial[ivol]);
          printf("eps: %d,  pot: %d\n", volEpsilon[ivol], volPotential[ivol]);
          printf("q: %d,  type: %d\n", volCharge[ivol], volBoundaryType[ivol]);
          printf("shape: %d,  material: %d\n", vshp, vmat);
          printf("eps: %d,  pot: %d\n", veps, vpot);
          printf("q: %d,  type: %d\n", vq, vtype);
        }
        */

        fprintf(fVoxel,
                "%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%4d\n",
                100.0 * point.X / LengthScale, 100.0 * point.Y / LengthScale,
                100.0 * point.Z / LengthScale, VoxelFX[k] / 100.0,
                VoxelFY[k] / 100.0, VoxelFZ[k] / 100.0, VoxelP[k] / LengthScale,
                ivol + 1);
      }
      fflush(fVoxel);  // file output over

      // printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
    }  // loop j
  }    // loop i

  fclose(fVoxel);

  free_dvector(VoxelFX, 0, nbZCells + 1);
  free_dvector(VoxelFY, 0, nbZCells + 1);
  free_dvector(VoxelFZ, 0, nbZCells + 1);
  free_dvector(VoxelP, 0, nbZCells + 1);

  return 0;
}  // end of VoxelFPR

// Compute 3dMap data for export to Garfield++.
// The 3dMap data for weighting field can be generated using the same function.
// After creation, the exported file has to be properly
// named and then imported by ComponentNeBem3dMap.
int MapFPR(void) {
  int dbgFn = 0;

  printf("\nPhysical potential and field computation for 3dMap data export\n");

  char MapInfoFile[256];
  strcpy(MapInfoFile, BCOutDir);
  strcat(MapInfoFile, "/MapInfo.out");
  FILE *fMapInfo = fopen(MapInfoFile, "w");
  if (fMapInfo == NULL) {
    neBEMMessage("MapFPR - MapInfoFile");
    return -1;
  }
  if (dbgFn) {
    printf("MapInfoFile.out created ...\n");
    fflush(stdout);
  }

  // In certain versions, we may have only the version number in the header and
  // nothing more. In that case, it is unsafe to assume that OptMap or
  // OptStaggerMap will be at all present in the output file. This decision
  // may need to be taken immediately after reading the MapVersion value.
  fprintf(fMapInfo, "%s\n", MapVersion);
  fprintf(fMapInfo, "%d\n", OptMap);
  fprintf(fMapInfo, "%d\n", OptStaggerMap);
  fprintf(fMapInfo, "%d\n", Map.NbXCells + 1);
  fprintf(fMapInfo, "%d\n", Map.NbYCells + 1);
  fprintf(fMapInfo, "%d\n", Map.NbZCells + 1);
  fprintf(fMapInfo, "%le %le\n", Map.Xmin * 100.0, Map.Xmax * 100.0);
  fprintf(fMapInfo, "%le %le\n", Map.Ymin * 100.0, Map.Ymax * 100.0);
  fprintf(fMapInfo, "%le %le\n", Map.Zmin * 100.0, Map.Zmax * 100.0);
  fprintf(fMapInfo, "%le\n", Map.XStagger * 100.0);
  fprintf(fMapInfo, "%le\n", Map.YStagger * 100.0);
  fprintf(fMapInfo, "%le\n", Map.ZStagger * 100.0);
  fprintf(fMapInfo, "MapFPR.out\n");
  // if (OptStaggerMap) fprintf(fMapInfo, "StgrMapFPR.out\n"); /// not being
  // read
  fclose(fMapInfo);

  char MapFile[256];
  strcpy(MapFile, BCOutDir);
  strcat(MapFile, "/MapFPR.out");
  FILE *fMap = fopen(MapFile, "w");
  if (fMap == NULL) {
    neBEMMessage("MapFPR - MapFile");
    return -1;
  }
  if (dbgFn) {
    printf("MapFPR.out created ...\n");
    fflush(stdout);
  }

  fprintf(
      fMap,
      "# X(cm)\tY(cm)\tZ(cm)\tFX(V/cm)\tFY(V/cm)\tFZ(V/cm)\tPot(V)\tRegion\n");

  int nbXCells = Map.NbXCells;
  int nbYCells = Map.NbYCells;
  int nbZCells = Map.NbZCells;
  double startX = Map.Xmin;
  double startY = Map.Ymin;
  double startZ = Map.Zmin;
  double delX = (Map.Xmax - Map.Xmin) / nbXCells;
  double delY = (Map.Ymax - Map.Ymin) / nbYCells;
  double delZ = (Map.Zmax - Map.Zmin) / nbZCells;

  double *MapFX = dvector(0, nbZCells + 1);
  double *MapFY = dvector(0, nbZCells + 1);
  double *MapFZ = dvector(0, nbZCells + 1);
  double *MapP = dvector(0, nbZCells + 1);

  if (dbgFn) {
    printf("nbXCells, nbYCells, nbZCells: %d, %d, %d\n", nbXCells, nbYCells,
           nbZCells);
    printf("startX, startY, startZ: %le, %le, %le\n", startX, startY, startZ);
    printf("delX, delY, delZ: %le, %le, %le\n", delX, delY, delZ);
    fflush(stdout);
  }

  for (int i = 1; i <= nbXCells + 1; ++i) {
    for (int j = 1; j <= nbYCells + 1; ++j) {
      if (dbgFn) {
        printf("MapFPR => i: %4d, j: %4d", i, j);
        fflush(stdout);
      }

      Point3D point;
#ifdef _OPENMP
      int nthreads = 1, tid = 0;
#pragma omp parallel private(nthreads, tid)
#endif
      {
#ifdef _OPENMP
        if (dbgFn) {
          tid = omp_get_thread_num();
          if (tid == 0) {
            nthreads = omp_get_num_threads();
            printf("Starting voxel computation with %d threads\n", nthreads);
          }
        }
#endif
        int k;
        Vector3D field;
        double potential;
#ifdef _OPENMP
#pragma omp for private(k, point, potential, field)
#endif
        for (k = 1; k <= nbZCells + 1; ++k) {
          point.X = startX + (i - 1) * delX;  // all 3 components need to be
          point.Y = startY + (j - 1) * delY;  // evaluated after pragma omp
          point.Z = startZ + (k - 1) * delZ;

          if (dbgFn) {
            printf("i, j, k: %d, %d, %d\n", i, j, k);
            printf("point X, Y, Z: %.8lg\t%.8lg\t%.8lg\n",
                   point.X / LengthScale, point.Y / LengthScale,
                   point.Z / LengthScale);
            fflush(stdout);
          }

          if (OptReadFastPF) {
            int fstatus = FastPFAtPoint(&point, &potential, &field);
            if (fstatus != 0) {
              neBEMMessage("wrong FastPFAtPoint return value in MapFPR\n");
              // return -1;
            }
          } else {
            neBEMMessage(
                "Suggestion: Use of FastVol can expedite generation of Map.\n");
            int fstatus = PFAtPoint(&point, &potential, &field);
            if (fstatus != 0) {
              neBEMMessage("wrong PFAtPoint return value in MapFPR\n");
              // return -1;
            }
          }

          if (dbgFn) {
            printf("%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\n",
                   point.X / LengthScale, point.Y / LengthScale,
                   point.Z / LengthScale, field.X, field.Y, field.Z,
                   potential / LengthScale);
            fflush(stdout);
          }

          MapFX[k] = field.X;
          MapFY[k] = field.Y;
          MapFZ[k] = field.Z;
          MapP[k] = potential;
        }  // loop k
      }    // pragma omp parallel

      for (int k = 1; k <= nbZCells + 1; ++k) {
        // file output
        point.X = startX + (i - 1) * delX;
        point.Y = startY + (j - 1) * delY;
        point.Z = startZ + (k - 1) * delZ;

        int ivol = neBEMVolumePoint(point.X, point.Y, point.Z);
        /*
        volMaterial[ivol];	// region linked to material
        neBEMVolumeDescription(ivol, &vshp, &vmat, &veps, &vpot, &vq, &vtype);
        if (dbgFn) {
          printf("volref: %d\n", ivol);
          printf("shape: %d,  material: %d\n", volShape[ivol],
        volMaterial[ivol]); printf("eps: %d,  pot: %d\n", volEpsilon[ivol],
        volPotential[ivol]); printf("q: %d,  type: %d\n", volCharge[ivol],
        volBoundaryType[ivol]); printf("shape: %d,  material: %d\n", vshp,
        vmat); printf("eps: %d,  pot: %d\n", veps, vpot); printf("q: %d,  type:
        %d\n", vq, vtype);
        }
        */

        fprintf(fMap, "%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%4d\n",
                100.0 * point.X / LengthScale, 100.0 * point.Y / LengthScale,
                100.0 * point.Z / LengthScale, MapFX[k] / 100.0,
                MapFY[k] / 100.0, MapFZ[k] / 100.0, MapP[k] / LengthScale,
                ivol + 1);
      }
      fflush(fMap);  // file output over
    }                // loop j
  }                  // loop i

  fclose(fMap);

  free_dvector(MapFX, 0, nbZCells + 1);
  free_dvector(MapFY, 0, nbZCells + 1);
  free_dvector(MapFZ, 0, nbZCells + 1);
  free_dvector(MapP, 0, nbZCells + 1);

  // If staggered map
  if (OptStaggerMap) {
    char StgrMapFile[256];
    strcpy(StgrMapFile, BCOutDir);
    strcat(StgrMapFile, "/StgrMapFPR.out");
    fMap = fopen(StgrMapFile, "w");
    if (fMap == NULL) {
      neBEMMessage("StgrMapFPR - Staggered MapFile");
      return -1;
    }
    if (dbgFn) {
      printf("StgrMapFPR.out created ...\n");
      fflush(stdout);
    }

    fprintf(
        fMap,
        "# "
        "X(cm)\tY(cm)\tZ(cm)\tFX(V/cm)\tFY(V/cm)\tFZ(V/cm)\tPot(V)\tRegion\n");
    // Very static stagger where X-shift is one map long
    double LX = (Map.Xmax - Map.Xmin);
    Map.Xmin = Map.Xmax;
    Map.Xmax = Map.Xmin + LX;
    double LY = (Map.Ymax - Map.Ymin);
    Map.Ymin = Map.Ymin + Map.YStagger;
    Map.Ymax = Map.Ymin + LY;
    nbXCells = Map.NbXCells;
    nbYCells = Map.NbYCells;
    nbZCells = Map.NbZCells;
    startX = Map.Xmin;
    // and y-shift is of the presecribed amount
    startY = Map.Ymin + Map.YStagger;
    startZ = Map.Zmin;
    delX = (Map.Xmax - Map.Xmin) / nbXCells;
    delY = (Map.Ymax - Map.Ymin) / nbYCells;
    delZ = (Map.Zmax - Map.Zmin) / nbZCells;

    MapFX = dvector(0, nbZCells + 1);
    MapFY = dvector(0, nbZCells + 1);
    MapFZ = dvector(0, nbZCells + 1);
    MapP = dvector(0, nbZCells + 1);

    if (dbgFn) {
      printf("nbXCells, nbYCells, nbZCells: %d, %d, %d\n", nbXCells, nbYCells,
             nbZCells);
      printf("startX, startY, startZ: %le, %le, %le\n", startX, startY, startZ);
      printf("delX, delY, delZ: %le, %le, %le\n", delX, delY, delZ);
      fflush(stdout);
    }

    for (int i = 1; i <= nbXCells + 1; ++i) {
      for (int j = 1; j <= nbYCells + 1; ++j) {
        if (dbgFn) {
          printf("StgrMapFPR => i: %4d, j: %4d", i, j);
          fflush(stdout);
        }

        Point3D point;
#ifdef _OPENMP
        int nthreads = 1, tid = 0;
#pragma omp parallel private(nthreads, tid)
#endif
        {
#ifdef _OPENMP
          if (dbgFn) {
            tid = omp_get_thread_num();
            if (tid == 0) {
              nthreads = omp_get_num_threads();
              printf("Starting voxel computation with %d threads\n", nthreads);
            }
          }  // if dbgFn
#endif
          int k;
          Vector3D field;
          double potential;
#ifdef _OPENMP
#pragma omp for private(k, point, potential, field)
#endif
          for (k = 1; k <= nbZCells + 1; ++k) {
            point.X = startX + (i - 1) * delX;  // all 3 components need to be
            point.Y = startY + (j - 1) * delY;  // evaluated after pragma omp
            point.Z = startZ + (k - 1) * delZ;

            if (dbgFn) {
              printf("i, j, k: %d, %d, %d\n", i, j, k);
              printf("point X, Y, Z: %.8lg\t%.8lg\t%.8lg\n",
                     point.X / LengthScale, point.Y / LengthScale,
                     point.Z / LengthScale);
              fflush(stdout);
            }

            if (OptReadFastPF) {
              int fstatus = FastPFAtPoint(&point, &potential, &field);
              if (fstatus != 0) {
                neBEMMessage("wrong FastPFAtPoint return value in MapFPR\n");
                // return -1;
              }
            } else {
              neBEMMessage(
                  "Suggestion: Use of FastVol can expedite generation of "
                  "Map.\n");
              int fstatus = PFAtPoint(&point, &potential, &field);
              if (fstatus != 0) {
                neBEMMessage("wrong PFAtPoint return value in MapFPR\n");
                // return -1;
              }
            }

            if (dbgFn) {
              printf("%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\n",
                     point.X / LengthScale, point.Y / LengthScale,
                     point.Z / LengthScale, field.X, field.Y, field.Z,
                     potential / LengthScale);
              fflush(stdout);
            }

            MapFX[k] = field.X;
            MapFY[k] = field.Y;
            MapFZ[k] = field.Z;
            MapP[k] = potential;
          }  // loop k
        }    // pragma omp parallel

        for (int k = 1; k <= nbZCells + 1; ++k) {
          // file output
          point.X = startX + (i - 1) * delX;
          point.Y = startY + (j - 1) * delY;
          point.Z = startZ + (k - 1) * delZ;

          int ivol = neBEMVolumePoint(point.X, point.Y, point.Z);
          /*
          volMaterial[ivol];	// region linked to material
          neBEMVolumeDescription(ivol, &vshp, &vmat, &veps, &vpot, &vq, &vtype);
          if (dbgFn) {
            printf("volref: %d\n", ivol);
            printf("shape: %d,  material: %d\n", volShape[ivol], volMaterial[ivol]); 
            printf("eps: %d,  pot: %d\n", volEpsilon[ivol], volPotential[ivol]); 
            printf("q: %d,  type: %d\n", volCharge[ivol], volBoundaryType[ivol]); 
            printf("shape: %d,  material: %d\n", vshp, vmat); 
            printf("eps: %d,  pot: %d\n", veps, vpot); 
            printf("q: %d, type: %d\n", vq, vtype);
          }
          */

          fprintf(
              fMap, "%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%4d\n",
              100.0 * point.X / LengthScale, 100.0 * point.Y / LengthScale,
              100.0 * point.Z / LengthScale, MapFX[k] / 100.0, MapFY[k] / 100.0,
              MapFZ[k] / 100.0, MapP[k] / LengthScale, ivol + 1);
        }              // for k <= nbZCells
        fflush(fMap);  // file output over

      }  // loop j
    }    // loop i

    fclose(fMap);

    free_dvector(MapFX, 0, nbZCells + 1);
    free_dvector(MapFY, 0, nbZCells + 1);
    free_dvector(MapFZ, 0, nbZCells + 1);
    free_dvector(MapP, 0, nbZCells + 1);
  }  // If staggered map

  return 0;
}  // end of MapFPR

// Compute potential and field in a mesh within the Fast Volume
// Note: PFAtPoint itself is evaluated in two steps: one for charges on the
//       elements, and the other for known charges.
int CreateFastVolPF(void) {
  // The following may be necessary only during the first time step / iteration
  // At present, FastVolElePF() considers both element and KnCh effects
  // and create one combined fast volume.
  int fstatus = CreateFastVolElePF();
  if (fstatus) {
    printf(
        "Problem in FastVolElePF being called from FastVolPF ... returning\n");
    return -1;
  }

  /*
  // The following is likely to change throughout the computation and necessary
  // at all time steps. However, in order to achieve computational economy, it
  // may be prudent to carry out the following only after several time steps.
  if (OptKnCh) {
    fstatus = CreateFastVolKnChPF();
    if (fstatus) {
      printf("Problem in FastVolKnChPF being called from FastVolPF...
  returning\n"); return -1;
    }
  }
  */

  return 0;
}  // CreateFastVolPF ends

// Compute potential and field in a mesh within the Fast Volume
// Possible pitfall: evaluation of n-skips
// As the name implies, this function uses the ElePFAtPoint function only.
// As a result, KnChPFAtPoint is not included in the resulting values.
// The effects due to known charge distributions can be included separately
// as if they are perturbations on the background system.
int CreateFastVolElePF(void) {
  int dbgFn = 0;
  int fstatus;

  int nbXCells;
  int nbYCells;
  int nbZCells;
  double startX;
  double startY;
  double startZ;
  double delX;
  double delY;
  double delZ;

  printf(
      "\nPhysical potential and field computation within basic fast volume\n");
  int bskip = 0, iskip = 0, jskip = 0, kskip = 0;

  // calculate n-skips based on NbPtSkip
  if (NbPtSkip) {
    int volptcnt = 0, endskip = 0;

    for (int block = 1; block <= FastVol.NbBlocks; ++block) {
      nbXCells = BlkNbXCells[block];
      nbYCells = BlkNbYCells[block];
      nbZCells = BlkNbZCells[block];
      for (int i = 1; i <= nbXCells + 1; ++i) {
        for (int j = 1; j <= nbYCells + 1; ++j) {
          for (int k = 1; k <= nbZCells + 1; ++k) {
            ++volptcnt;

            if (volptcnt == NbPtSkip) {
              bskip = block - 1;
              iskip = i - 1;
              jskip = j - 1;
              kskip = k;
              endskip = 1;
            }

            if (endskip) break;
          }
          if (endskip) break;
        }
        if (endskip) break;
      }
      if (endskip) break;
    }
    if (dbgFn) {
      printf(
          "Basic fast volume => bskip, iskip, jskip, kskip: %d, %d, %d, %d\n",
          bskip, iskip, jskip, kskip);
    }
  }  // NbPtSkip

  char FastVolPFFile[256];
  strcpy(FastVolPFFile, BCOutDir);
  strcat(FastVolPFFile, "/FastVolPF.out");
  FILE *fFastVolPF = fopen(FastVolPFFile, "w");
  if (fFastVolPF == NULL) {
    neBEMMessage("FastVolPF - FastVolPFFile");
    return -1;
  }
  fprintf(fFastVolPF, "#block\tX\tY\tZ\tPot\tFX\tFY\tFZ\n");

  if (dbgFn) {
    printf("FastVolPF.out created ...\n");
    fflush(stdout);
  }

  for (int block = 1 + bskip; block <= FastVol.NbBlocks; ++block) {
    nbXCells = BlkNbXCells[block];
    nbYCells = BlkNbYCells[block];
    nbZCells = BlkNbZCells[block];
    startX = FastVol.CrnrX;
    startY = FastVol.CrnrY;
    startZ = BlkCrnrZ[block];
    delX = FastVol.LX / nbXCells;
    delY = FastVol.LY / nbYCells;
    delZ = BlkLZ[block] / nbZCells;
    printf(
        "NbBlocks: %d, block: %d, nbXCells: %d, nbYCells: %d, nbZCells: %d\n",
        FastVol.NbBlocks, block, nbXCells, nbYCells, nbZCells);

    if (dbgFn) {
      printf("block: %d\n", block);
      printf("nbXCells, nbYCells, nbZCells: %d, %d, %d\n", nbXCells, nbYCells,
             nbZCells);
      printf("startX, startY, startZ: %le, %le, %le\n", startX, startY, startZ);
      printf("delX, delY, delZ: %le, %le, %le\n", delX, delY, delZ);
      printf("bskip, iskip, jskip, kskip: %d, %d, %d, %d\n", bskip, iskip,
             jskip, kskip);
      fflush(stdout);
    }
    // total number of points in a given block
    // int blktotpt = (nbXCells + 1) * (nbYCells + 1) * (nbZCells + 1);
    for (int i = 1 + iskip; i <= nbXCells + 1; ++i) {
      for (int j = 1 + jskip; j <= nbYCells + 1; ++j) {
        printf("Fast volume => block: %3d, i: %4d, j: %4d", block, i, j);
        fflush(stdout);

        Point3D point;
#ifdef _OPENMP
        int nthreads = 1, tid = 0;
#pragma omp parallel private(nthreads, tid)
#endif
        {
#ifdef _OPENMP
          if (dbgFn) {
            tid = omp_get_thread_num();
            if (tid == 0) {
              nthreads = omp_get_num_threads();
              printf("Starting fast volume computation with %d threads\n",
                     nthreads);
            }
          }
#endif
          int k;
          int omitFlag;
          double potential;
          Vector3D field;
#ifdef _OPENMP
#pragma omp for private(k, point, omitFlag, potential, field)
#endif
          for (k = 1 + kskip; k <= nbZCells + 1; ++k) {
            potential = 0.0;
            field.X = field.Y = field.Z = 0.0;

            point.X = startX + (i - 1) * delX;
            point.Y = startY + (j - 1) * delY;
            point.Z = startZ + (k - 1) * delZ;

            // Check whether the point falls within a volume that should be
            // ignored
            omitFlag = 0;
            for (int omit = 1; omit <= FastVol.NbOmitVols; ++omit) {
              if ((point.X > OmitVolCrnrX[omit]) &&
                  (point.X < OmitVolCrnrX[omit] + OmitVolLX[omit]) &&
                  (point.Y > OmitVolCrnrY[omit]) &&
                  (point.Y < OmitVolCrnrY[omit] + OmitVolLY[omit]) &&
                  (point.Z > OmitVolCrnrZ[omit]) &&
                  (point.Z < OmitVolCrnrZ[omit] + OmitVolLZ[omit])) {
                omitFlag = 1;
                break;
              }
            }  // loop over omitted volumes

            if (dbgFn) {
              printf("block, i, j, k: %d, %d, %d, %d\n", block, i, j, k);
              printf("point X, Y, Z: %.8lg\t%.8lg\t%.8lg\n",
                     point.X / LengthScale, point.Y / LengthScale,
                     point.Z / LengthScale);
              printf("omitFlag: %d\n", omitFlag);
              fflush(stdout);
            }

            if (omitFlag) {
              potential = field.X = field.Y = field.Z = 0.0;
            } else {
              // fstatus = ElePFAtPoint(&point, &potential, &field);
              fstatus =
                  PFAtPoint(&point, &potential, &field);  // both ele and KnCh
              if (fstatus != 0) {
                neBEMMessage(
                    "wrong ElePFAtPoint return value in FastVolElePF.\n");
                // return -1;
              }
            }  // else omitFlag
            if (dbgFn) {
              printf("%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\n",
                     point.X / LengthScale, point.Y / LengthScale,
                     point.Z / LengthScale, potential / LengthScale, field.X,
                     field.Y, field.Z);
              fflush(stdout);
            }

            FastPot[block][i][j][k] = potential;
            FastFX[block][i][j][k] = field.X;
            FastFY[block][i][j][k] = field.Y;
            FastFZ[block][i][j][k] = field.Z;
          }  // loop k
        }    // pragma omp parallel

        for (int k = 1 + kskip; k <= nbZCells + 1; ++k)  // file output
        {
          point.X = startX + (i - 1) * delX;
          point.Y = startY + (j - 1) * delY;
          point.Z = startZ + (k - 1) * delZ;

          fprintf(fFastVolPF,
                  "%4d\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\n",
                  block, point.X / LengthScale, point.Y / LengthScale,
                  point.Z / LengthScale, FastPot[block][i][j][k],
                  FastFX[block][i][j][k], FastFY[block][i][j][k],
                  FastFZ[block][i][j][k]);
        }
        fflush(fFastVolPF);  // file output over

        printf(
            "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
            "\b\b\b\b\b\b\b\b\b\b");
      }  // loop j
    }    // loop i
  }      // loop block

  fclose(fFastVolPF);

  if (OptStaggerFastVol) {
    printf("Potential and field computation within staggered fast volume\n");

    bskip = iskip = jskip = kskip = 0;

    // calculate n-skips based on NbStgPtSkip
    if (NbStgPtSkip) {
      int volptcnt = 0, endskip = 0;

      for (int block = 1; block <= FastVol.NbBlocks; ++block) {
        nbXCells = BlkNbXCells[block];
        nbYCells = BlkNbYCells[block];
        nbZCells = BlkNbZCells[block];
        for (int i = 1; i <= nbXCells + 1; ++i) {
          for (int j = 1; j <= nbYCells + 1; ++j) {
            for (int k = 1; k <= nbZCells + 1; ++k) {
              ++volptcnt;

              if (volptcnt == NbStgPtSkip) {
                bskip = block - 1;
                iskip = i - 1;
                jskip = j - 1;
                kskip = k;
                endskip = 1;
              }

              if (endskip) break;
            }
            if (endskip) break;
          }
          if (endskip) break;
        }
        if (endskip) break;
      }
      if (dbgFn) {
        printf(
            "Staggered volume => bskip, iskip, jskip, kskip: %d, %d, %d, %d\n",
            bskip, iskip, jskip, kskip);
      }
    }  // NbStgPtSkip

    char StgFastVolPFFile[256];
    FILE *fStgFastVolPF;
    strcpy(StgFastVolPFFile, BCOutDir);
    strcat(StgFastVolPFFile, "/StgFastVolPF.out");
    fStgFastVolPF = fopen(StgFastVolPFFile, "w");
    if (fStgFastVolPF == NULL) {
      neBEMMessage("FastVolPF - StgFastVolPFFile");
      return -1;
    }
    fprintf(fStgFastVolPF, "#block\tX\tY\tZ\tPot\tFX\tFY\tFZ\n");

    if (dbgFn) {
      printf("StgFastVolPF.out created ...\n");
      fflush(stdout);
    }

    for (int block = 1 + bskip; block <= FastVol.NbBlocks; ++block) {
      nbXCells = BlkNbXCells[block];
      nbYCells = BlkNbYCells[block];
      nbZCells = BlkNbZCells[block];
      startX = FastVol.CrnrX + FastVol.LX;
      startY = FastVol.CrnrY + FastVol.YStagger;
      startZ = BlkCrnrZ[block];
      delX = FastVol.LX / nbXCells;
      delY = FastVol.LY / nbYCells;
      delZ = BlkLZ[block] / nbZCells;
      printf(
          "NbBlocks: %d, block: %d, nbXCells: %d, nbYCells: %d, nbZCells: %d\n",
          FastVol.NbBlocks, block, nbXCells, nbYCells, nbZCells);

      if (dbgFn) {
        printf("block: %d\n", block);
        printf("nbXCells, nbYCells, nbZCells: %d, %d, %d\n", nbXCells, nbYCells,
               nbZCells);
        printf("startX, startY, startZ: %le, %le, %le\n", startX, startY,
               startZ);
        printf("delX, delY, delZ: %le, %le, %le\n", delX, delY, delZ);
        printf("bskip, iskip, jskip, kskip: %d, %d, %d, %d\n", bskip, iskip,
               jskip, kskip);
        fflush(stdout);
      }

      // int blktotpt = (nbXCells + 1) * (nbYCells + 1) * (nbZCells + 1);
      for (int i = 1 + iskip; i <= nbXCells + 1; ++i) {
        for (int j = 1 + jskip; j <= nbYCells + 1; ++j) {
          printf("Fast volume => block: %3d, i: %4d, j: %4d", block, i, j);
          fflush(stdout);

          Point3D point;
#ifdef _OPENMP
          int nthreads = 1, tid = 0;
#pragma omp parallel private(nthreads, tid)
#endif
          {
#ifdef _OPENMP
            if (dbgFn) {
              tid = omp_get_thread_num();
              if (tid == 0) {
                nthreads = omp_get_num_threads();
                printf(
                    "Starting staggered fast volume computation with %d "
                    "threads\n",
                    nthreads);
              }
            }
#endif
            int k;
            int omitFlag;
            double potential;
            Vector3D field;
#ifdef _OPENMP
#pragma omp for private(k, point, omitFlag, potential, field)
#endif
            for (k = 1 + kskip; k <= nbZCells + 1; ++k) {
              potential = 0.0;
              field.X = field.Y = field.Z = 0.0;

              point.X = startX + (i - 1) * delX;
              point.Y = startY + (j - 1) * delY;
              point.Z = startZ + (k - 1) * delZ;

              // Check whether point falls within a volume that should be
              // ignored
              omitFlag = 0;
              for (int omit = 1; omit <= FastVol.NbOmitVols; ++omit) {
                if ((point.X > OmitVolCrnrX[omit] + FastVol.LX) &&
                    (point.X <
                     OmitVolCrnrX[omit] + OmitVolLX[omit] + FastVol.LX) &&
                    (point.Y > OmitVolCrnrY[omit] + FastVol.YStagger) &&
                    (point.Y <
                     OmitVolCrnrY[omit] + OmitVolLY[omit] + FastVol.YStagger) &&
                    (point.Z > OmitVolCrnrZ[omit]) &&
                    (point.Z < OmitVolCrnrZ[omit] + OmitVolLZ[omit])) {
                  omitFlag = 1;
                  break;
                }
              }  // loop over omitted volumes

              if (dbgFn) {
                printf("point X, Y, Z: %.8lg\t%.8lg\t%.8lg\n",
                       point.X / LengthScale, point.Y / LengthScale,
                       point.Z / LengthScale);
                printf("omitFlag: %d\n", omitFlag);
                fflush(stdout);
              }

              if (omitFlag) {
                potential = field.X = field.Y = field.Z = 0.0;
              } else {
                // fstatus = ElePFAtPoint(&point, &potential, &field);
                fstatus =
                    PFAtPoint(&point, &potential, &field);  // both ele & KnCh
                if (fstatus != 0) {
                  neBEMMessage(
                      "wrong PFAtPoint return value in FastVolElePF.\n");
                  // return -1;
                }
              } // else omitFlag
              if (dbgFn) {
                printf("%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\n",
                       point.X / LengthScale, point.Y / LengthScale,
                       point.Z / LengthScale, potential / LengthScale, field.X,
                       field.Y, field.Z);
                fflush(stdout);
              }

              StgFastPot[block][i][j][k] = potential;
              StgFastFX[block][i][j][k] = field.X;
              StgFastFY[block][i][j][k] = field.Y;
              StgFastFZ[block][i][j][k] = field.Z;
            }  // loop k
          }    // pragma omp

          for (int k = 1 + kskip; k <= nbZCells + 1; ++k)  // file output
          {
            point.X = startX + (i - 1) * delX;
            point.Y = startY + (j - 1) * delY;
            point.Z = startZ + (k - 1) * delZ;

            fprintf(fStgFastVolPF,
                    "%4d\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\n",
                    block, point.X / LengthScale, point.Y / LengthScale,
                    point.Z / LengthScale, StgFastPot[block][i][j][k],
                    StgFastFX[block][i][j][k], StgFastFY[block][i][j][k],
                    StgFastFZ[block][i][j][k]);
          }
          fflush(fStgFastVolPF);  // file output over

          printf(
              "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
              "\b\b\b\b\b\b\b\b\b\b\b");
        }  // loop j
      }    // loop i
    }      // loop block

    fclose(fStgFastVolPF);
  }  // if OptStaggerFastVol

  return 0;
}  // CreateFastVolElePF ends

// Gives three components of the total Potential and flux in the global
// coordinate system due to all the elements using the results stored in
// the FAST volume mesh.
int FastPFAtPoint(Point3D *globalP, double *Potential, Vector3D *globalF) {
  int dbgFn = 0;
  double Xpt = globalP->X;
  double Ypt = globalP->Y;
  double Zpt = globalP->Z;
  double RptVolLX = FastVol.LX;
  double RptVolLY = FastVol.LY;
  double RptVolLZ = FastVol.LZ;
  double CornerX = FastVol.CrnrX;
  double CornerY = FastVol.CrnrY;
  double CornerZ = FastVol.CrnrZ;
  double TriLin(double xd, double yd, double zd, double c000, double c100,
                double c010, double c001, double c110, double c101, double c011,
                double c111);

  // First of all, check how the point in question should be treated ...

  // Check whether the point falls within a volume that is not regarded as
  // FastVol
  for (int ignore = 1; ignore <= FastVol.NbIgnoreVols; ++ignore) {
    if ((Xpt >= (IgnoreVolCrnrX[ignore])) &&
        (Xpt <= (IgnoreVolCrnrX[ignore] + IgnoreVolLX[ignore])) &&
        (Ypt >= (IgnoreVolCrnrY[ignore])) &&
        (Ypt <= (IgnoreVolCrnrY[ignore] + IgnoreVolLY[ignore])) &&
        (Zpt >= (IgnoreVolCrnrZ[ignore])) &&
        (Zpt <= (IgnoreVolCrnrZ[ignore] + IgnoreVolLZ[ignore]))) {
      if (dbgFn)
        neBEMMessage("In FastPFAtPoint: point in an ignored volume!\n");

      int fstatus = PFAtPoint(globalP, Potential, globalF);
      if (fstatus != 0) {
        neBEMMessage("wrong PFAtPoint return value in FastVolPF.\n");
        return -1;
      } else
        return 0;
    }
  }  // loop over ignored volumes

  // If not ignored, the point qualifies for FastVol evaluation ...

  // for a staggered fast volume, the volume repeated in X is larger
  if (OptStaggerFastVol) {
    RptVolLX += FastVol.LX;
  }
  if (dbgFn) {
    printf("\nin FastPFAtPoint\n");
    printf("x, y, z: %g, %g, %g\n", Xpt, Ypt, Zpt);
    printf("RptVolLX, RptVolLY, RptVolLZ: %g, %g, %g\n", RptVolLX, RptVolLY,
           RptVolLZ);
    printf("CornerX, CornerY, CornerZ: %g, %g, %g\n", CornerX, CornerY,
           CornerZ);
    printf("Nb of blocks: %d\n", FastVol.NbBlocks);
    for (int block = 1; block <= FastVol.NbBlocks; ++block) {
      printf("NbOfXCells: %d\n", BlkNbXCells[block]);
      printf("NbOfYCells: %d\n", BlkNbYCells[block]);
      printf("NbOfZCells: %d\n", BlkNbZCells[block]);
      printf("LZ: %le\n", BlkLZ[block]);
      printf("CornerZ: %le\n", BlkCrnrZ[block]);
    }
  }

  // Find equivalent position inside the basic / staggered volume.
  // real distance from volume corner
  double dx = Xpt - CornerX;
  double dy = Ypt - CornerY;
  double dz = Zpt - CornerZ;
  if (dbgFn)
    printf("real dx, dy, dz from volume corner: %g, %g, %g\n", dx, dy, dz);

  int NbFastVolX = (int)(dx / RptVolLX);
  if (dx < 0.0) --NbFastVolX;
  int NbFastVolY = (int)(dy / RptVolLY);
  if (dy < 0.0) --NbFastVolY;
  int NbFastVolZ = (int)(dz / RptVolLZ);
  if (dz < 0.0) --NbFastVolZ;
  if (dbgFn)
    printf("Volumes in x, y, z: %d, %d, %d\n", NbFastVolX, NbFastVolY,
           NbFastVolZ);

  // equivalent distances from fast volume corner
  dx -= NbFastVolX * RptVolLX;
  dy -= NbFastVolY * RptVolLY;
  dz -= NbFastVolZ * RptVolLZ;
  // The following conditions should never happen - generate an error message
  if (dx < 0.0) {
    dx = 0.0;
    neBEMMessage("equiv dx < 0.0 - not correct!\n");
  }
  if (dy < 0.0) {
    dy = 0.0;
    neBEMMessage("equiv dy < 0.0 - not correct!\n");
  }
  if (dz < 0.0) {
    dz = 0.0;
    neBEMMessage("equiv dz < 0.0 - not correct!\n");
  }
  if (dx > RptVolLX) {
    dx = RptVolLX;
    neBEMMessage("equiv dx > RptVolLX - not correct!\n");
  }
  if (dy > RptVolLY) {
    dy = RptVolLY;
    neBEMMessage("equiv dy > RptVolLY - not correct!\n");
  }
  if (dz > RptVolLZ) {
    dz = RptVolLZ;
    neBEMMessage("equiv dz > RptVolLZ - not correct!\n");
  }
  if (dbgFn)
    printf("equivalent dist from corner - dx, dy, dz: %g, %g, %g\n", dx, dy,
           dz);

  // Take care of possible trouble-makers
  if (dx < MINDIST) dx = MINDIST;  // -ve dx has already been made equal to 0
  if (dy < MINDIST) dy = MINDIST;
  if (dz < MINDIST) dz = MINDIST;
  if ((RptVolLX - dx) < MINDIST)
    dx = RptVolLX - MINDIST;  // dx > RptVolLX taken care of
  if ((RptVolLY - dy) < MINDIST) dy = RptVolLY - MINDIST;
  if ((RptVolLZ - dz) < MINDIST) dz = RptVolLZ - MINDIST;
  // For staggered volumes, there is another plane where difficulties may occur
  if ((dx <= FastVol.LX) && (FastVol.LX - dx) < MINDIST)
    dx = FastVol.LX - MINDIST;
  // else if((dx > FastVol.LX) && (fabs(FastVol.LX-dx) < MINDIST))
  else if ((dx > FastVol.LX) && (dx - FastVol.LX) < MINDIST)  // more inttuitive
    dx = FastVol.LX + MINDIST;
  if (dbgFn)
    printf("equivalent dist adjusted - dx, dy, dz: %g, %g, %g\n", dx, dy, dz);

  // If volume is staggered, we have a few more things to do before finalizing
  // the values of equivalent distance
  // sector identification
  // _................__________________
  // |    .     .     |    Sector 3    |
  // |     .   .      |                |
  // |       .        |       .        |
  // |                |     .   .      |
  // |    Sector 2    |    .     .     |
  // |----------------|     .   .      |
  // |                |       .        |
  // |       .        |                |
  // |     .   .      |                |
  // |    .     .     |----------------|
  // |     .   .      |    Sector 4    |
  // |       .        |       .        |
  // |                |     .   .      |
  // |    Sector 1    |    .     .     |
  // |----------------|................|

  int sector = 1;  // kept outside `if' since this is necessary further below
  if (OptStaggerFastVol) {
    if ((dx >= 0.0) && (dx <= FastVol.LX) && (dy >= 0.0) &&
        (dy <= FastVol.LY)) {
      // point lies in sector 1, everything remains unchanged
      sector = 1;
    } else if ((dx >= 0.0) && (dx <= FastVol.LX) && (dy > FastVol.LY) &&
               (dy <= FastVol.LY + FastVol.YStagger)) {
      // point lies in sector 2, move basic volume one step up
      sector = 2;
      ++NbFastVolY;
      CornerY += FastVol.LY;  // repeat length in Y is LY
      dy -= FastVol.LY;
    } else if ((dx > FastVol.LX) && (dx <= 2.0 * FastVol.LX) &&
               (dy >= FastVol.YStagger) &&
               (dy <= FastVol.LY + FastVol.YStagger)) {
      // point lies in sector 3, pt in staggered vol, change corner coords
      sector = 3;
      CornerX += FastVol.LX;
      CornerY += FastVol.YStagger;
      dx -= FastVol.LX;
      dy -= FastVol.YStagger;
    } else if ((dx > FastVol.LX) && (dx <= 2.0 * FastVol.LX) && (dy >= 0.0) &&
               (dy < FastVol.YStagger)) {
      // point lies in sector 4, move basic volume one step down and consider
      // staggered fast volume
      sector = 4;
      --NbFastVolY;
      CornerX += FastVol.LX;  // in the staggered part of the repeated volume
      CornerY -= (FastVol.LY - FastVol.YStagger);
      dx -= FastVol.LX;
      dy += (FastVol.LY - FastVol.YStagger);
    } else {
      neBEMMessage("FastPFAtPoint: point in none of the sectors!\n");
    }
    if (dbgFn) printf("stagger modified dx, dy, dz: %g, %g, %g\n", dx, dy, dz);
  }

  // Take care of possible trouble-makers - once more
  if (dx < MINDIST) dx = MINDIST;  // -ve dx has already been made equal to 0
  if (dy < MINDIST) dy = MINDIST;
  if (dz < MINDIST) dz = MINDIST;
  if ((RptVolLX - dx) < MINDIST)
    dx = RptVolLX - MINDIST;  // dx > RptVolLX taken care of
  if ((RptVolLY - dy) < MINDIST) dy = RptVolLY - MINDIST;
  if ((RptVolLZ - dz) < MINDIST) dz = RptVolLZ - MINDIST;
  // For staggered volumes, there is another plane where difficulties may occur
  if ((dx <= FastVol.LX) && (FastVol.LX - dx) < MINDIST)
    dx = FastVol.LX - MINDIST;
  // else if((dx > FastVol.LX) && (fabs(FastVol.LX-dx) < MINDIST))
  else if ((dx > FastVol.LX) && (dx - FastVol.LX) < MINDIST)  // more intuitive
    dx = FastVol.LX + MINDIST;
  if (dbgFn)
    printf("equivalent dist adjusted for staggered: %g, %g, %g\n", dx, dy, dz);

  // Check whether the point falls within a volume that is omitted
  for(int omit = 1; omit <= FastVol.NbOmitVols; ++omit) {
    if ((dx >= (OmitVolCrnrX[omit]-FastVol.CrnrX)) && 
        (dx <= (OmitVolCrnrX[omit]+OmitVolLX[omit]-FastVol.CrnrX)) && 
        (dy >= (OmitVolCrnrY[omit]-FastVol.CrnrY)) && 
        (dy <= (OmitVolCrnrY[omit]+OmitVolLY[omit]-FastVol.CrnrY)) && 
        (dz >= (OmitVolCrnrZ[omit]-FastVol.CrnrZ)) && 
        (dz <= (OmitVolCrnrZ[omit]+OmitVolLZ[omit]-FastVol.CrnrZ))) {
      neBEMMessage("In FastPFAtPoint: point in an omitted volume!\n"); 
      *Potential = 0.0; 
      globalF->X = 0.0; globalF->Y = 0.0; globalF->Z = 0.0;
    }
  }	// loop over omitted volumes

  // Find the block in which the point lies
  int thisBlock = 0;
  for (int block = 1; block <= FastVol.NbBlocks; ++block) {
    double blkBtmZ = BlkCrnrZ[block] - CornerZ;  // since CornerZ has been
    double blkTopZ = blkBtmZ + BlkLZ[block];     // subtracted from dz already
    if (dbgFn) {
      printf("block, dz, blkBtmZ, blkTopZ: %d, %lg, %lg, %lg\n", block, dz,
             blkBtmZ, blkTopZ);
    }

    // take care of difficult situations
    if ((dz <= blkBtmZ) && ((blkBtmZ - dz) < MINDIST)) dz = blkBtmZ - MINDIST;
    if ((dz >= blkBtmZ) && ((dz - blkBtmZ) < MINDIST)) dz = blkBtmZ + MINDIST;
    if ((dz <= blkTopZ) && ((blkTopZ - dz) < MINDIST)) dz = blkTopZ - MINDIST;
    if ((dz >= blkTopZ) && ((dz - blkTopZ) < MINDIST)) dz = blkTopZ + MINDIST;

    if ((dz >= blkBtmZ) && (dz <= blkTopZ)) {
      thisBlock = block;
      break;
    }
  }
  if (!thisBlock) {
    neBEMMessage("FastPFAtPoint: point in none of the blocks!\n");
  }

  int nbXCells = BlkNbXCells[thisBlock];
  int nbYCells = BlkNbYCells[thisBlock];
  int nbZCells = BlkNbZCells[thisBlock];
  double delX = FastVol.LX / nbXCells;
  double delY = FastVol.LY / nbYCells;
  double delZ = BlkLZ[thisBlock] / nbZCells;
  dz -= (BlkCrnrZ[thisBlock] - CornerZ);  // distance from the block corner

  if (dbgFn) {
    printf("thisBlock: %d\n", thisBlock);
    printf("nbXCells, nbYCells, nbZCells: %d, %d, %d\n", nbXCells, nbYCells,
           nbZCells);
    printf("BlkCrnrZ: %lg\n", BlkCrnrZ[thisBlock]);
    printf("delX, delY, delZ: %le, %le, %le\n", delX, delY, delZ);
    printf("dz: %lg\n", dz);
    fflush(stdout);
  }

  // Find cell in block of basic / staggered volume within which the point lies
  int celli = (int)(dx / delX) + 1;  // Find cell in which the point lies
  if (celli < 1) {
    celli = 1;
    dx = 0.5 * delX;
    neBEMMessage("FastPFAtPoint - celli < 1\n");
  }
  if (celli > nbXCells) {
    celli = nbXCells;
    dx = FastVol.LX - 0.5 * delX;
    neBEMMessage("FastPFAtPoint - celli > nbXCells\n");
  }
  int cellj = (int)(dy / delY) + 1;
  if (cellj < 1) {
    cellj = 1;
    dy = 0.5 * delY;
    neBEMMessage("FastPFAtPoint - cellj < 1\n");
  }
  if (cellj > nbYCells) {
    cellj = nbYCells;
    dy = FastVol.LY - 0.5 * delY;
    neBEMMessage("FastPFAtPoint - cellj > nbYCells\n");
  }
  int cellk = (int)(dz / delZ) + 1;
  if (cellk < 1) {
    cellk = 1;
    dz = 0.5 * delX;
    neBEMMessage("FastPFAtPoint - cellk < 1\n");
  }
  if (cellk > nbZCells) {
    cellk = nbZCells;
    dz = FastVol.LZ - 0.5 * delZ;
    neBEMMessage("FastPFAtPoint - cellk > nbZCells\n");
  }
  if (dbgFn) printf("Cells in x, y, z: %d, %d, %d\n", celli, cellj, cellk);

  // Interpolate potential and field at the point using the corner values of
  // of the cell and, if necessary, of the neighbouring cells
  // These gradients can also be calculated while computing the potential and
  // field at the cells and stored in memory, provided enough memory is
  // available

  // distances from cell corner
  double dxcellcrnr = dx - (double)(celli - 1) * delX;
  double dycellcrnr = dy - (double)(cellj - 1) * delY;
  double dzcellcrnr = dz - (double)(cellk - 1) * delZ;
  if (dbgFn)
    printf("cell crnr dx, dy, dz: %g, %g, %g\n", dxcellcrnr, dycellcrnr,
           dzcellcrnr);

  // normalized distances
  double xd = dxcellcrnr / delX;  // xd = (x-x0)/(x1-x0)
  double yd = dycellcrnr / delY;  // etc
  double zd = dzcellcrnr / delZ;
  if (xd <= 0.0) xd = 0.0;
  if (yd <= 0.0) yd = 0.0;
  if (zd <= 0.0) zd = 0.0;
  if (xd >= 1.0) xd = 1.0;
  if (yd >= 1.0) yd = 1.0;
  if (zd >= 1.0) zd = 1.0;

  // corner values of potential and field
  double P000 = FastPot[thisBlock][celli][cellj][cellk];  // lowest corner
  double FX000 = FastFX[thisBlock][celli][cellj][cellk];
  double FY000 = FastFY[thisBlock][celli][cellj][cellk];
  double FZ000 = FastFZ[thisBlock][celli][cellj][cellk];
  double P100 = FastPot[thisBlock][celli + 1][cellj][cellk];
  double FX100 = FastFX[thisBlock][celli + 1][cellj][cellk];
  double FY100 = FastFY[thisBlock][celli + 1][cellj][cellk];
  double FZ100 = FastFZ[thisBlock][celli + 1][cellj][cellk];
  double P010 = FastPot[thisBlock][celli][cellj + 1][cellk];
  double FX010 = FastFX[thisBlock][celli][cellj + 1][cellk];
  double FY010 = FastFY[thisBlock][celli][cellj + 1][cellk];
  double FZ010 = FastFZ[thisBlock][celli][cellj + 1][cellk];
  double P001 = FastPot[thisBlock][celli][cellj][cellk + 1];
  double FX001 = FastFX[thisBlock][celli][cellj][cellk + 1];
  double FY001 = FastFY[thisBlock][celli][cellj][cellk + 1];
  double FZ001 = FastFZ[thisBlock][celli][cellj][cellk + 1];
  double P110 = FastPot[thisBlock][celli + 1][cellj + 1][cellk];
  double FX110 = FastFX[thisBlock][celli + 1][cellj + 1][cellk];
  double FY110 = FastFY[thisBlock][celli + 1][cellj + 1][cellk];
  double FZ110 = FastFZ[thisBlock][celli + 1][cellj + 1][cellk];
  double P101 = FastPot[thisBlock][celli + 1][cellj][cellk + 1];
  double FX101 = FastFX[thisBlock][celli + 1][cellj][cellk + 1];
  double FY101 = FastFY[thisBlock][celli + 1][cellj][cellk + 1];
  double FZ101 = FastFZ[thisBlock][celli + 1][cellj][cellk + 1];
  double P011 = FastPot[thisBlock][celli][cellj + 1][cellk + 1];
  double FX011 = FastFX[thisBlock][celli][cellj + 1][cellk + 1];
  double FY011 = FastFY[thisBlock][celli][cellj + 1][cellk + 1];
  double FZ011 = FastFZ[thisBlock][celli][cellj + 1][cellk + 1];
  double P111 = FastPot[thisBlock][celli + 1][cellj + 1][cellk + 1];
  double FX111 = FastFX[thisBlock][celli + 1][cellj + 1][cellk + 1];
  double FY111 = FastFY[thisBlock][celli + 1][cellj + 1][cellk + 1];
  double FZ111 = FastFZ[thisBlock][celli + 1][cellj + 1][cellk + 1];
  if (OptStaggerFastVol) {
    if (sector == 1) {  // nothing to be done
    }
    if (sector == 2) {  // volume shifted up but point not in the staggered part
    }
    if (sector == 3) {  // staggered volume
      P000 = StgFastPot[thisBlock][celli][cellj][cellk];
      FX000 = StgFastFX[thisBlock][celli][cellj][cellk];
      FY000 = StgFastFY[thisBlock][celli][cellj][cellk];
      FZ000 = StgFastFZ[thisBlock][celli][cellj][cellk];
      P100 = StgFastPot[thisBlock][celli + 1][cellj][cellk];
      FX100 = StgFastFX[thisBlock][celli + 1][cellj][cellk];
      FY100 = StgFastFY[thisBlock][celli + 1][cellj][cellk];
      FZ100 = StgFastFZ[thisBlock][celli + 1][cellj][cellk];
      P010 = StgFastPot[thisBlock][celli][cellj + 1][cellk];
      FX010 = StgFastFX[thisBlock][celli][cellj + 1][cellk];
      FY010 = StgFastFY[thisBlock][celli][cellj + 1][cellk];
      FZ010 = StgFastFZ[thisBlock][celli][cellj + 1][cellk];
      P001 = StgFastPot[thisBlock][celli][cellj][cellk + 1];
      FX001 = StgFastFX[thisBlock][celli][cellj][cellk + 1];
      FY001 = StgFastFY[thisBlock][celli][cellj][cellk + 1];
      FZ001 = StgFastFZ[thisBlock][celli][cellj][cellk + 1];
      P110 = StgFastPot[thisBlock][celli + 1][cellj + 1][cellk];
      FX110 = StgFastFX[thisBlock][celli + 1][cellj + 1][cellk];
      FY110 = StgFastFY[thisBlock][celli + 1][cellj + 1][cellk];
      FZ110 = StgFastFZ[thisBlock][celli + 1][cellj + 1][cellk];
      P101 = StgFastPot[thisBlock][celli + 1][cellj][cellk + 1];
      FX101 = StgFastFX[thisBlock][celli + 1][cellj][cellk + 1];
      FY101 = StgFastFY[thisBlock][celli + 1][cellj][cellk + 1];
      FZ101 = StgFastFZ[thisBlock][celli + 1][cellj][cellk + 1];
      P011 = StgFastPot[thisBlock][celli][cellj + 1][cellk + 1];
      FX011 = StgFastFX[thisBlock][celli][cellj + 1][cellk + 1];
      FY011 = StgFastFY[thisBlock][celli][cellj + 1][cellk + 1];
      FZ011 = StgFastFZ[thisBlock][celli][cellj + 1][cellk + 1];
      P111 = StgFastPot[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FX111 = StgFastFX[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FY111 = StgFastFY[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FZ111 = StgFastFZ[thisBlock][celli + 1][cellj + 1][cellk + 1];
    }
    if (sector == 4) {  // volume shifted down and point in the staggered part
      P000 = StgFastPot[thisBlock][celli][cellj][cellk];
      FX000 = StgFastFX[thisBlock][celli][cellj][cellk];
      FY000 = StgFastFY[thisBlock][celli][cellj][cellk];
      FZ000 = StgFastFZ[thisBlock][celli][cellj][cellk];
      P100 = StgFastPot[thisBlock][celli + 1][cellj][cellk];
      FX100 = StgFastFX[thisBlock][celli + 1][cellj][cellk];
      FY100 = StgFastFY[thisBlock][celli + 1][cellj][cellk];
      FZ100 = StgFastFZ[thisBlock][celli + 1][cellj][cellk];
      P010 = StgFastPot[thisBlock][celli][cellj + 1][cellk];
      FX010 = StgFastFX[thisBlock][celli][cellj + 1][cellk];
      FY010 = StgFastFY[thisBlock][celli][cellj + 1][cellk];
      FZ010 = StgFastFZ[thisBlock][celli][cellj + 1][cellk];
      P001 = StgFastPot[thisBlock][celli][cellj][cellk + 1];
      FX001 = StgFastFX[thisBlock][celli][cellj][cellk + 1];
      FY001 = StgFastFY[thisBlock][celli][cellj][cellk + 1];
      FZ001 = StgFastFZ[thisBlock][celli][cellj][cellk + 1];
      P110 = StgFastPot[thisBlock][celli + 1][cellj + 1][cellk];
      FX110 = StgFastFX[thisBlock][celli + 1][cellj + 1][cellk];
      FY110 = StgFastFY[thisBlock][celli + 1][cellj + 1][cellk];
      FZ110 = StgFastFZ[thisBlock][celli + 1][cellj + 1][cellk];
      P101 = StgFastPot[thisBlock][celli + 1][cellj][cellk + 1];
      FX101 = StgFastFX[thisBlock][celli + 1][cellj][cellk + 1];
      FY101 = StgFastFY[thisBlock][celli + 1][cellj][cellk + 1];
      FZ101 = StgFastFZ[thisBlock][celli + 1][cellj][cellk + 1];
      P011 = StgFastPot[thisBlock][celli][cellj + 1][cellk + 1];
      FX011 = StgFastFX[thisBlock][celli][cellj + 1][cellk + 1];
      FY011 = StgFastFY[thisBlock][celli][cellj + 1][cellk + 1];
      FZ011 = StgFastFZ[thisBlock][celli][cellj + 1][cellk + 1];
      P111 = StgFastPot[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FX111 = StgFastFX[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FY111 = StgFastFY[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FZ111 = StgFastFZ[thisBlock][celli + 1][cellj + 1][cellk + 1];
    }
  }  // if OptStaggerFastVol

  double intP =
      TriLin(xd, yd, zd, P000, P100, P010, P001, P110, P101, P011, P111);
  double intFX = TriLin(xd, yd, zd, FX000, FX100, FX010, FX001, FX110, FX101,
                        FX011, FX111);
  double intFY = TriLin(xd, yd, zd, FY000, FY100, FY010, FY001, FY110, FY101,
                        FY011, FY111);
  double intFZ = TriLin(xd, yd, zd, FZ000, FZ100, FZ010, FZ001, FZ110, FZ101,
                        FZ011, FZ111);

  *Potential = intP;
  globalF->X = intFX;
  globalF->Y = intFY;
  globalF->Z = intFZ;

  if (dbgFn) {
    printf("Cell corner values:\n");
    printf("Potential: %g, %g, %g, %g\n", P000, P100, P010, P001);
    printf("Potential: %g, %g, %g, %g\n", P110, P101, P011, P111);
    printf("FastFX: %g, %g, %g, %g\n", FX000, FX100, FX010, FX001);
    printf("FastFX: %g, %g, %g, %g\n", FX110, FX101, FX011, FX111);
    printf("FastFY: %g, %g, %g, %g\n", FY000, FY100, FY010, FY001);
    printf("FastFY: %g, %g, %g, %g\n", FY110, FY101, FY011, FY111);
    printf("FastFZ: %g, %g, %g, %g\n", FZ000, FZ100, FZ010, FZ001);
    printf("FastFZ: %g, %g, %g, %g\n", FZ110, FZ101, FZ011, FZ111);
    printf("Pot, FX, FY, FZ: %g, %g, %g, %g\n", *Potential, globalF->X,
           globalF->Y, globalF->Z);
  }

  if (dbgFn) {
    printf("out FastPFAtPoint\n");
    fflush(stdout);
  }

  return 0;
}  // FastPFAtPoint ends

// There could be a function
int FastElePFAtPoint(Point3D * /*globalP*/, double * /*Potential*/,
                     Vector3D * /*globalF*/) {
  return 0;
}  // FastElePFAtPoint ends

// Gives three components of the total Potential and flux in the global
// coordinate system due to all the known charges using the results stored in
// the FAST volume KnCh mesh.
int FastKnChPFAtPoint(Point3D *globalP, double *Potential, Vector3D *globalF) {
  int dbgFn = 0;
  double Xpt = globalP->X;
  double Ypt = globalP->Y;
  double Zpt = globalP->Z;
  double RptVolLX = FastVol.LX;
  double RptVolLY = FastVol.LY;
  double RptVolLZ = FastVol.LZ;
  double CornerX = FastVol.CrnrX;
  double CornerY = FastVol.CrnrY;
  double CornerZ = FastVol.CrnrZ;
  double TriLin(double xd, double yd, double zd, double c000, double c100,
                double c010, double c001, double c110, double c101, double c011,
                double c111);

  // First of all, check how the point in question should be treated ...

  // Check whether the point falls within a volume that is not regarded as
  // FastVol
  for (int ignore = 1; ignore <= FastVol.NbIgnoreVols; ++ignore) {
    if ((Xpt >= (IgnoreVolCrnrX[ignore])) &&
        (Xpt <= (IgnoreVolCrnrX[ignore] + IgnoreVolLX[ignore])) &&
        (Ypt >= (IgnoreVolCrnrY[ignore])) &&
        (Ypt <= (IgnoreVolCrnrY[ignore] + IgnoreVolLY[ignore])) &&
        (Zpt >= (IgnoreVolCrnrZ[ignore])) &&
        (Zpt <= (IgnoreVolCrnrZ[ignore] + IgnoreVolLZ[ignore]))) {
      if (dbgFn)
        neBEMMessage("In FastKnChPFAtPoint: point in an ignored volume!\n");

      int fstatus = KnChPFAtPoint(globalP, Potential, globalF);
      if (fstatus != 0) {
        neBEMMessage("wrong KnChPFAtPoint return value in FastVolKnChPF.\n");
        return -1;
      } else
        return 0;
    }
  }  // loop over ignored volumes

  // If not ignored, the point qualifies for FastVol evaluation ...

  // for a staggered fast volume, the volume repeated in X is larger
  if (OptStaggerFastVol) {
    RptVolLX += FastVol.LX;
  }
  if (dbgFn) {
    printf("\nin FastKnChPFAtPoint\n");
    printf("x, y, z: %g, %g, %g\n", Xpt, Ypt, Zpt);
    printf("RptVolLX, RptVolLY, RptVolLZ: %g, %g, %g\n", RptVolLX, RptVolLY,
           RptVolLZ);
    printf("CornerX, CornerY, CornerZ: %g, %g, %g\n", CornerX, CornerY,
           CornerZ);
    printf("Nb of blocks: %d\n", FastVol.NbBlocks);
    for (int block = 1; block <= FastVol.NbBlocks; ++block) {
      printf("NbOfXCells: %d\n", BlkNbXCells[block]);
      printf("NbOfYCells: %d\n", BlkNbYCells[block]);
      printf("NbOfZCells: %d\n", BlkNbZCells[block]);
      printf("LZ: %le\n", BlkLZ[block]);
      printf("CornerZ: %le\n", BlkCrnrZ[block]);
    }
  }

  // Find equivalent position inside the basic / staggered volume.
  // real distance from volume corner
  double dx = Xpt - CornerX;
  double dy = Ypt - CornerY;
  double dz = Zpt - CornerZ;
  if (dbgFn)
    printf("real dx, dy, dz from volume corner: %g, %g, %g\n", dx, dy, dz);

  int NbFastVolX = (int)(dx / RptVolLX);
  if (dx < 0.0) --NbFastVolX;
  int NbFastVolY = (int)(dy / RptVolLY);
  if (dy < 0.0) --NbFastVolY;
  int NbFastVolZ = (int)(dz / RptVolLZ);
  if (dz < 0.0) --NbFastVolZ;
  if (dbgFn)
    printf("Volumes in x, y, z: %d, %d, %d\n", NbFastVolX, NbFastVolY,
           NbFastVolZ);

  // equivalent distances from fast volume corner
  dx -= NbFastVolX * RptVolLX;
  dy -= NbFastVolY * RptVolLY;
  dz -= NbFastVolZ * RptVolLZ;
  // The following conditions should never happen - generate an error message
  if (dx < 0.0) {
    dx = 0.0;
    neBEMMessage("equiv dx < 0.0 - not correct!\n");
  }
  if (dy < 0.0) {
    dy = 0.0;
    neBEMMessage("equiv dy < 0.0 - not correct!\n");
  }
  if (dz < 0.0) {
    dz = 0.0;
    neBEMMessage("equiv dz < 0.0 - not correct!\n");
  }
  if (dx > RptVolLX) {
    dx = RptVolLX;
    neBEMMessage("equiv dx > RptVolLX - not correct!\n");
  }
  if (dy > RptVolLY) {
    dy = RptVolLY;
    neBEMMessage("equiv dy > RptVolLY - not correct!\n");
  }
  if (dz > RptVolLZ) {
    dz = RptVolLZ;
    neBEMMessage("equiv dz > RptVolLZ - not correct!\n");
  }
  if (dbgFn)
    printf("equivalent dist from corner - dx, dy, dz: %g, %g, %g\n", dx, dy,
           dz);

  // Take care of possible trouble-makers
  if (dx < MINDIST) dx = MINDIST;  // -ve dx has already been made equal to 0
  if (dy < MINDIST) dy = MINDIST;
  if (dz < MINDIST) dz = MINDIST;
  if ((RptVolLX - dx) < MINDIST)
    dx = RptVolLX - MINDIST;  // dx > RptVolLX taken care of
  if ((RptVolLY - dy) < MINDIST) dy = RptVolLY - MINDIST;
  if ((RptVolLZ - dz) < MINDIST) dz = RptVolLZ - MINDIST;
  // For staggered volumes, there is another plane where difficulties may occur
  if ((dx <= FastVol.LX) && (FastVol.LX - dx) < MINDIST)
    dx = FastVol.LX - MINDIST;
  // else if((dx > FastVol.LX) && (fabs(FastVol.LX-dx) < MINDIST))
  else if ((dx > FastVol.LX) && (dx - FastVol.LX) < MINDIST)  // more inttuitive
    dx = FastVol.LX + MINDIST;
  if (dbgFn)
    printf("equivalent dist adjusted - dx, dy, dz: %g, %g, %g\n", dx, dy, dz);

  // If volume is staggered, we have a few more things to do before finalizing
  // the values of equivalent distance
  // sector identification
  // _................__________________
  // |    .     .     |    Sector 3    |
  // |     .   .      |                |
  // |       .        |       .        |
  // |                |     .   .      |
  // |    Sector 2    |    .     .     |
  // |----------------|     .   .      |
  // |                |       .        |
  // |       .        |                |
  // |     .   .      |                |
  // |    .     .     |----------------|
  // |     .   .      |    Sector 4    |
  // |       .        |       .        |
  // |                |     .   .      |
  // |    Sector 1    |    .     .     |
  // |----------------|................|

  int sector = 1;  // kept outside `if' since this is necessary further below
  if (OptStaggerFastVol) {
    if ((dx >= 0.0) && (dx <= FastVol.LX) && (dy >= 0.0) &&
        (dy <= FastVol.LY)) {
      // point lies in sector 1, everything remains unchanged
      sector = 1;
    } else if ((dx >= 0.0) && (dx <= FastVol.LX) && (dy > FastVol.LY) &&
               (dy <= FastVol.LY + FastVol.YStagger)) {
      // point lies in sector 2, move basic volume one step up
      sector = 2;
      ++NbFastVolY;
      CornerY += FastVol.LY;  // repeat length in Y is LY
      dy -= FastVol.LY;
    } else if ((dx > FastVol.LX) && (dx <= 2.0 * FastVol.LX) &&
               (dy >= FastVol.YStagger) &&
               (dy <= FastVol.LY + FastVol.YStagger)) {
      // point lies in sector 3, pt in staggered vol, change corner coords
      sector = 3;
      CornerX += FastVol.LX;
      CornerY += FastVol.YStagger;
      dx -= FastVol.LX;
      dy -= FastVol.YStagger;
    } else if ((dx > FastVol.LX) && (dx <= 2.0 * FastVol.LX) && (dy >= 0.0) &&
               (dy < FastVol.YStagger)) {
      // point lies in sector 4, move basic volume one step down and consider
      // staggered fast volume
      sector = 4;
      --NbFastVolY;
      CornerX += FastVol.LX;  // in the staggered part of the repeated volume
      CornerY -= (FastVol.LY - FastVol.YStagger);
      dx -= FastVol.LX;
      dy += (FastVol.LY - FastVol.YStagger);
    } else {
      neBEMMessage("FastKnChPFAtPoint: point in none of the sectors!\n");
    }
    if (dbgFn) printf("stagger modified dx, dy, dz: %g, %g, %g\n", dx, dy, dz);
  }

  // Take care of possible trouble-makers - once more
  if (dx < MINDIST) dx = MINDIST;  // -ve dx has already been made equal to 0
  if (dy < MINDIST) dy = MINDIST;
  if (dz < MINDIST) dz = MINDIST;
  if ((RptVolLX - dx) < MINDIST)
    dx = RptVolLX - MINDIST;  // dx > RptVolLX taken care of
  if ((RptVolLY - dy) < MINDIST) dy = RptVolLY - MINDIST;
  if ((RptVolLZ - dz) < MINDIST) dz = RptVolLZ - MINDIST;
  // For staggered volumes, there is another plane where difficulties may occur
  if ((dx <= FastVol.LX) && (FastVol.LX - dx) < MINDIST)
    dx = FastVol.LX - MINDIST;
  // else if((dx > FastVol.LX) && (fabs(FastVol.LX-dx) < MINDIST))
  else if ((dx > FastVol.LX) && (dx - FastVol.LX) < MINDIST)  // more intuitive
    dx = FastVol.LX + MINDIST;
  if (dbgFn)
    printf("equivalent dist adjusted for staggered: %g, %g, %g\n", dx, dy, dz);

  // Check whether the point falls within a volume that is omitted
  for(int omit = 1; omit <= FastVol.NbOmitVols; ++omit) {
    if((dx >= (OmitVolCrnrX[omit]-FastVol.CrnrX)) && 
       (dx <= (OmitVolCrnrX[omit]+OmitVolLX[omit]-FastVol.CrnrX)) && 
       (dy >= (OmitVolCrnrY[omit]-FastVol.CrnrY)) && 
       (dy <= (OmitVolCrnrY[omit]+OmitVolLY[omit]-FastVol.CrnrY)) && 
       (dz >= (OmitVolCrnrZ[omit]-FastVol.CrnrZ)) && 
       (dz <= (OmitVolCrnrZ[omit]+OmitVolLZ[omit]-FastVol.CrnrZ))) {
      neBEMMessage("In FastKnChPFAtPoint: point in an omitted volume!\n");
      *Potential = 0.0; 
      globalF->X = 0.0; globalF->Y = 0.0; globalF->Z = 0.0;
    }
  }	// loop over omitted volumes

  int thisBlock = 0;
  for (int block = 1; block <= FastVol.NbBlocks; ++block) {
    double blkBtmZ = BlkCrnrZ[block] - CornerZ;  // since CornerZ has been
    double blkTopZ = blkBtmZ + BlkLZ[block];     // subtracted from dz already
    if (dbgFn) {
      printf("block, dz, blkBtmZ, blkTopZ: %d, %lg, %lg, %lg\n", block, dz,
             blkBtmZ, blkTopZ);
    }

    // take care of difficult situations
    if ((dz <= blkBtmZ) && ((blkBtmZ - dz) < MINDIST)) dz = blkBtmZ - MINDIST;
    if ((dz >= blkBtmZ) && ((dz - blkBtmZ) < MINDIST)) dz = blkBtmZ + MINDIST;
    if ((dz <= blkTopZ) && ((blkTopZ - dz) < MINDIST)) dz = blkTopZ - MINDIST;
    if ((dz >= blkTopZ) && ((dz - blkTopZ) < MINDIST)) dz = blkTopZ + MINDIST;

    if ((dz >= blkBtmZ) && (dz <= blkTopZ)) {
      thisBlock = block;
      break;
    }
  }
  if (!thisBlock) {
    neBEMMessage("FastKnChPFAtPoint: point in none of the blocks!\n");
  }

  int nbXCells = BlkNbXCells[thisBlock];
  int nbYCells = BlkNbYCells[thisBlock];
  int nbZCells = BlkNbZCells[thisBlock];
  double delX = FastVol.LX / nbXCells;
  double delY = FastVol.LY / nbYCells;
  double delZ = BlkLZ[thisBlock] / nbZCells;
  dz -= (BlkCrnrZ[thisBlock] - CornerZ);  // distance from the block corner

  if (dbgFn) {
    printf("thisBlock: %d\n", thisBlock);
    printf("nbXCells, nbYCells, nbZCells: %d, %d, %d\n", nbXCells, nbYCells,
           nbZCells);
    printf("BlkCrnrZ: %lg\n", BlkCrnrZ[thisBlock]);
    printf("delX, delY, delZ: %le, %le, %le\n", delX, delY, delZ);
    printf("dz: %lg\n", dz);
    fflush(stdout);
  }

  // Find cell in block of basic / staggered volume within which the point lies
  int celli = (int)(dx / delX) + 1;  // Find cell in which the point lies
  if (celli < 1) {
    celli = 1;
    dx = 0.5 * delX;
    neBEMMessage("FastKnChPFAtPoint - celli < 1\n");
  }
  if (celli > nbXCells) {
    celli = nbXCells;
    dx = FastVol.LX - 0.5 * delX;
    neBEMMessage("FastKnChPFAtPoint - celli > nbXCells\n");
  }
  int cellj = (int)(dy / delY) + 1;
  if (cellj < 1) {
    cellj = 1;
    dy = 0.5 * delY;
    neBEMMessage("FastKnChPFAtPoint - cellj < 1\n");
  }
  if (cellj > nbYCells) {
    cellj = nbYCells;
    dy = FastVol.LY - 0.5 * delY;
    neBEMMessage("FastKnChPFAtPoint - cellj > nbYCells\n");
  }
  int cellk = (int)(dz / delZ) + 1;
  if (cellk < 1) {
    cellk = 1;
    dz = 0.5 * delX;
    neBEMMessage("FastKnChPFAtPoint - cellk < 1\n");
  }
  if (cellk > nbZCells) {
    cellk = nbZCells;
    dz = FastVol.LZ - 0.5 * delZ;
    neBEMMessage("FastKnChPFAtPoint - cellk > nbZCells\n");
  }
  if (dbgFn) printf("Cells in x, y, z: %d, %d, %d\n", celli, cellj, cellk);

  // Interpolate potential and field at the point using the corner values of
  // of the cell and, if necessary, of the neighbouring cells
  // These gradients can also be calculated while computing the potential and
  // field at the cells and stored in memory, provided enough memory is
  // available

  // distances from cell corner
  double dxcellcrnr = dx - (double)(celli - 1) * delX;
  double dycellcrnr = dy - (double)(cellj - 1) * delY;
  double dzcellcrnr = dz - (double)(cellk - 1) * delZ;
  if (dbgFn)
    printf("cell crnr dx, dy, dz: %g, %g, %g\n", dxcellcrnr, dycellcrnr,
           dzcellcrnr);

  // normalized distances
  double xd = dxcellcrnr / delX;  // xd = (x-x0)/(x1-x0)
  double yd = dycellcrnr / delY;  // etc
  double zd = dzcellcrnr / delZ;
  if (xd <= 0.0) xd = 0.0;
  if (yd <= 0.0) yd = 0.0;
  if (zd <= 0.0) zd = 0.0;
  if (xd >= 1.0) xd = 1.0;
  if (yd >= 1.0) yd = 1.0;
  if (zd >= 1.0) zd = 1.0;

  // corner values of potential and field
  double P000 = FastPotKnCh[thisBlock][celli][cellj][cellk];  // lowest corner
  double FX000 = FastFXKnCh[thisBlock][celli][cellj][cellk];
  double FY000 = FastFYKnCh[thisBlock][celli][cellj][cellk];
  double FZ000 = FastFZKnCh[thisBlock][celli][cellj][cellk];
  double P100 = FastPotKnCh[thisBlock][celli + 1][cellj][cellk];
  double FX100 = FastFXKnCh[thisBlock][celli + 1][cellj][cellk];
  double FY100 = FastFYKnCh[thisBlock][celli + 1][cellj][cellk];
  double FZ100 = FastFZKnCh[thisBlock][celli + 1][cellj][cellk];
  double P010 = FastPotKnCh[thisBlock][celli][cellj + 1][cellk];
  double FX010 = FastFXKnCh[thisBlock][celli][cellj + 1][cellk];
  double FY010 = FastFYKnCh[thisBlock][celli][cellj + 1][cellk];
  double FZ010 = FastFZKnCh[thisBlock][celli][cellj + 1][cellk];
  double P001 = FastPotKnCh[thisBlock][celli][cellj][cellk + 1];
  double FX001 = FastFXKnCh[thisBlock][celli][cellj][cellk + 1];
  double FY001 = FastFYKnCh[thisBlock][celli][cellj][cellk + 1];
  double FZ001 = FastFZKnCh[thisBlock][celli][cellj][cellk + 1];
  double P110 = FastPotKnCh[thisBlock][celli + 1][cellj + 1][cellk];
  double FX110 = FastFXKnCh[thisBlock][celli + 1][cellj + 1][cellk];
  double FY110 = FastFYKnCh[thisBlock][celli + 1][cellj + 1][cellk];
  double FZ110 = FastFZKnCh[thisBlock][celli + 1][cellj + 1][cellk];
  double P101 = FastPotKnCh[thisBlock][celli + 1][cellj][cellk + 1];
  double FX101 = FastFXKnCh[thisBlock][celli + 1][cellj][cellk + 1];
  double FY101 = FastFYKnCh[thisBlock][celli + 1][cellj][cellk + 1];
  double FZ101 = FastFZKnCh[thisBlock][celli + 1][cellj][cellk + 1];
  double P011 = FastPotKnCh[thisBlock][celli][cellj + 1][cellk + 1];
  double FX011 = FastFXKnCh[thisBlock][celli][cellj + 1][cellk + 1];
  double FY011 = FastFYKnCh[thisBlock][celli][cellj + 1][cellk + 1];
  double FZ011 = FastFZKnCh[thisBlock][celli][cellj + 1][cellk + 1];
  double P111 = FastPotKnCh[thisBlock][celli + 1][cellj + 1][cellk + 1];
  double FX111 = FastFXKnCh[thisBlock][celli + 1][cellj + 1][cellk + 1];
  double FY111 = FastFYKnCh[thisBlock][celli + 1][cellj + 1][cellk + 1];
  double FZ111 = FastFZKnCh[thisBlock][celli + 1][cellj + 1][cellk + 1];
  if (OptStaggerFastVol) {
    if (sector == 1) {  // nothing to be done
    }
    if (sector == 2) {  // volume shifted up but point not in the staggered part
    }
    if (sector == 3) {  // staggered volume
      P000 = StgFastPotKnCh[thisBlock][celli][cellj][cellk];
      FX000 = StgFastFXKnCh[thisBlock][celli][cellj][cellk];
      FY000 = StgFastFYKnCh[thisBlock][celli][cellj][cellk];
      FZ000 = StgFastFZKnCh[thisBlock][celli][cellj][cellk];
      P100 = StgFastPotKnCh[thisBlock][celli + 1][cellj][cellk];
      FX100 = StgFastFXKnCh[thisBlock][celli + 1][cellj][cellk];
      FY100 = StgFastFYKnCh[thisBlock][celli + 1][cellj][cellk];
      FZ100 = StgFastFZKnCh[thisBlock][celli + 1][cellj][cellk];
      P010 = StgFastPotKnCh[thisBlock][celli][cellj + 1][cellk];
      FX010 = StgFastFXKnCh[thisBlock][celli][cellj + 1][cellk];
      FY010 = StgFastFYKnCh[thisBlock][celli][cellj + 1][cellk];
      FZ010 = StgFastFZKnCh[thisBlock][celli][cellj + 1][cellk];
      P001 = StgFastPotKnCh[thisBlock][celli][cellj][cellk + 1];
      FX001 = StgFastFXKnCh[thisBlock][celli][cellj][cellk + 1];
      FY001 = StgFastFYKnCh[thisBlock][celli][cellj][cellk + 1];
      FZ001 = StgFastFZKnCh[thisBlock][celli][cellj][cellk + 1];
      P110 = StgFastPotKnCh[thisBlock][celli + 1][cellj + 1][cellk];
      FX110 = StgFastFXKnCh[thisBlock][celli + 1][cellj + 1][cellk];
      FY110 = StgFastFYKnCh[thisBlock][celli + 1][cellj + 1][cellk];
      FZ110 = StgFastFZKnCh[thisBlock][celli + 1][cellj + 1][cellk];
      P101 = StgFastPotKnCh[thisBlock][celli + 1][cellj][cellk + 1];
      FX101 = StgFastFXKnCh[thisBlock][celli + 1][cellj][cellk + 1];
      FY101 = StgFastFYKnCh[thisBlock][celli + 1][cellj][cellk + 1];
      FZ101 = StgFastFZKnCh[thisBlock][celli + 1][cellj][cellk + 1];
      P011 = StgFastPotKnCh[thisBlock][celli][cellj + 1][cellk + 1];
      FX011 = StgFastFXKnCh[thisBlock][celli][cellj + 1][cellk + 1];
      FY011 = StgFastFYKnCh[thisBlock][celli][cellj + 1][cellk + 1];
      FZ011 = StgFastFZKnCh[thisBlock][celli][cellj + 1][cellk + 1];
      P111 = StgFastPotKnCh[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FX111 = StgFastFXKnCh[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FY111 = StgFastFYKnCh[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FZ111 = StgFastFZKnCh[thisBlock][celli + 1][cellj + 1][cellk + 1];
    }
    if (sector == 4) {  // volume shifted down and point in the staggered part
      P000 = StgFastPotKnCh[thisBlock][celli][cellj][cellk];
      FX000 = StgFastFXKnCh[thisBlock][celli][cellj][cellk];
      FY000 = StgFastFYKnCh[thisBlock][celli][cellj][cellk];
      FZ000 = StgFastFZKnCh[thisBlock][celli][cellj][cellk];
      P100 = StgFastPotKnCh[thisBlock][celli + 1][cellj][cellk];
      FX100 = StgFastFXKnCh[thisBlock][celli + 1][cellj][cellk];
      FY100 = StgFastFYKnCh[thisBlock][celli + 1][cellj][cellk];
      FZ100 = StgFastFZKnCh[thisBlock][celli + 1][cellj][cellk];
      P010 = StgFastPotKnCh[thisBlock][celli][cellj + 1][cellk];
      FX010 = StgFastFXKnCh[thisBlock][celli][cellj + 1][cellk];
      FY010 = StgFastFYKnCh[thisBlock][celli][cellj + 1][cellk];
      FZ010 = StgFastFZKnCh[thisBlock][celli][cellj + 1][cellk];
      P001 = StgFastPotKnCh[thisBlock][celli][cellj][cellk + 1];
      FX001 = StgFastFXKnCh[thisBlock][celli][cellj][cellk + 1];
      FY001 = StgFastFYKnCh[thisBlock][celli][cellj][cellk + 1];
      FZ001 = StgFastFZKnCh[thisBlock][celli][cellj][cellk + 1];
      P110 = StgFastPotKnCh[thisBlock][celli + 1][cellj + 1][cellk];
      FX110 = StgFastFXKnCh[thisBlock][celli + 1][cellj + 1][cellk];
      FY110 = StgFastFYKnCh[thisBlock][celli + 1][cellj + 1][cellk];
      FZ110 = StgFastFZKnCh[thisBlock][celli + 1][cellj + 1][cellk];
      P101 = StgFastPotKnCh[thisBlock][celli + 1][cellj][cellk + 1];
      FX101 = StgFastFXKnCh[thisBlock][celli + 1][cellj][cellk + 1];
      FY101 = StgFastFYKnCh[thisBlock][celli + 1][cellj][cellk + 1];
      FZ101 = StgFastFZKnCh[thisBlock][celli + 1][cellj][cellk + 1];
      P011 = StgFastPotKnCh[thisBlock][celli][cellj + 1][cellk + 1];
      FX011 = StgFastFXKnCh[thisBlock][celli][cellj + 1][cellk + 1];
      FY011 = StgFastFYKnCh[thisBlock][celli][cellj + 1][cellk + 1];
      FZ011 = StgFastFZKnCh[thisBlock][celli][cellj + 1][cellk + 1];
      P111 = StgFastPotKnCh[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FX111 = StgFastFXKnCh[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FY111 = StgFastFYKnCh[thisBlock][celli + 1][cellj + 1][cellk + 1];
      FZ111 = StgFastFZKnCh[thisBlock][celli + 1][cellj + 1][cellk + 1];
    }
  }  // if OptStaggerFastVol

  double intP =
      TriLin(xd, yd, zd, P000, P100, P010, P001, P110, P101, P011, P111);
  double intFX = TriLin(xd, yd, zd, FX000, FX100, FX010, FX001, FX110, FX101,
                        FX011, FX111);
  double intFY = TriLin(xd, yd, zd, FY000, FY100, FY010, FY001, FY110, FY101,
                        FY011, FY111);
  double intFZ = TriLin(xd, yd, zd, FZ000, FZ100, FZ010, FZ001, FZ110, FZ101,
                        FZ011, FZ111);

  *Potential = intP;
  globalF->X = intFX;
  globalF->Y = intFY;
  globalF->Z = intFZ;

  if (dbgFn) {
    printf("Cell corner values:\n");
    printf("Potential: %g, %g, %g, %g\n", P000, P100, P010, P001);
    printf("Potential: %g, %g, %g, %g\n", P110, P101, P011, P111);
    printf("FastFX: %g, %g, %g, %g\n", FX000, FX100, FX010, FX001);
    printf("FastFX: %g, %g, %g, %g\n", FX110, FX101, FX011, FX111);
    printf("FastFY: %g, %g, %g, %g\n", FY000, FY100, FY010, FY001);
    printf("FastFY: %g, %g, %g, %g\n", FY110, FY101, FY011, FY111);
    printf("FastFZ: %g, %g, %g, %g\n", FZ000, FZ100, FZ010, FZ001);
    printf("FastFZ: %g, %g, %g, %g\n", FZ110, FZ101, FZ011, FZ111);
    printf("Pot, FX, FY, FZ: %g, %g, %g, %g\n", *Potential, globalF->X,
           globalF->Y, globalF->Z);
  }

  if (dbgFn) {
    printf("out FastKnChPFAtPoint\n");
    fflush(stdout);
  }

  return 0;
}  // FastKnChPFAtPoint ends

// Gives three components of weighting field in the global coordinate system
// due to all the elements
// Note that local evaluation of influence and additional influences have not
// been incorporated here. Iff local evaluation show a substantial advantage
// over the cleaner function call, we'll implement the former in this function.
// This function may be merged with PFAtPoint since the only change is in
// the use of weighting field charge density instead of the physical charge
// denstiy. However, care should be taken to check the last to points mentioned
// in this function - VSystemChargeZero and effects of known charge densities.
// Multi-threading implemented in the following routine.
// Gives three components of the total Potential and flux in the global
// coordinate system due to all the elements
int WtFldPFAtPoint(Point3D *globalP, double *Potential, Vector3D *globalF,
                   int IdWtField) {
  int dbgFn = 0;

  const double xfld = globalP->X;
  const double yfld = globalP->Y;
  const double zfld = globalP->Z;

  // Compute Potential and field at different locations
  *Potential = globalF->X = globalF->Y = globalF->Z = 0.0;

  // Effects due to base primitives and their repetitions are considered in the
  // local coordinate system of the primitive (or element), while effects due to
  // mirror elements and their repetitions are considered in the global
  // coordinate system (GCS). This works because the direction cosines of a
  // primitive (and its elements) and those of its repetitions are the same.
  // As a result, we can do just one transformation from local to global at the
  // end of calculations related to a primitive. This can save substantial
  // computation if a discretized version of the primitive is being used since
  // we avoid one unnecessary transformation for each element that comprises a
  // primitive.
  // Begin with primitive description of the device

  // Scope in OpenMP: Variables in the global data space are accessible to all
  // threads, while variables in a thread's private space is accessible to the
  // thread only (there are several variations - copying outside region etc)
  // Field point remains the same - kept outside private
  // source point changes with change in primitive - private
  // TransformationMatrix changes - kept within private (Note: matrices with
  // fixed dimensions can be maintained, but those with dynamic allocation
  // can not).
  double *pPot = dvector(1, NbPrimitives);
  double *plFx = dvector(1, NbPrimitives);  // field components in LCS
  double *plFy = dvector(1, NbPrimitives);  // for a primitive
  double *plFz = dvector(1, NbPrimitives);  // and its other incarnations

  for (int prim = 1; prim <= NbPrimitives; ++prim) {
    pPot[prim] = plFx[prim] = plFy[prim] = plFz[prim] = 0.0;
  }

#ifdef _OPENMP
  int tid = 0, nthreads = 1;
#pragma omp parallel private(tid, nthreads)
#endif
  {
#ifdef _OPENMP
    if (dbgFn) {
      tid = omp_get_thread_num();
      if (tid == 0) {
        nthreads = omp_get_num_threads();
        printf("PFAtPoint computation with %d threads\n", nthreads);
      }
    }
#endif
// by default, nested parallelization is off in C
#ifdef _OPENMP
#pragma omp for
#endif
    for (int primsrc = 1; primsrc <= NbPrimitives; ++primsrc) {
      if (dbgFn) {
        printf("Evaluating effect of primsrc %d using on %lg, %lg, %lg\n",
               primsrc, xfld, yfld, zfld);
        fflush(stdout);
      }

      const double xpsrc = PrimOriginX[primsrc];
      const double ypsrc = PrimOriginY[primsrc];
      const double zpsrc = PrimOriginZ[primsrc];

      // Field in the local frame.
      double lFx = 0.;
      double lFy = 0.;
      double lFz = 0.;

      // Set up transform matrix for this primitive, which is also the same
      // for all the elements belonging to this primitive.
      double TransformationMatrix[3][3];
      TransformationMatrix[0][0] = PrimDC[primsrc].XUnit.X;
      TransformationMatrix[0][1] = PrimDC[primsrc].XUnit.Y;
      TransformationMatrix[0][2] = PrimDC[primsrc].XUnit.Z;
      TransformationMatrix[1][0] = PrimDC[primsrc].YUnit.X;
      TransformationMatrix[1][1] = PrimDC[primsrc].YUnit.Y;
      TransformationMatrix[1][2] = PrimDC[primsrc].YUnit.Z;
      TransformationMatrix[2][0] = PrimDC[primsrc].ZUnit.X;
      TransformationMatrix[2][1] = PrimDC[primsrc].ZUnit.Y;
      TransformationMatrix[2][2] = PrimDC[primsrc].ZUnit.Z;

      // The total influence is due to primitives on the basic device and due to
      // virtual primitives arising out of repetition, reflection etc and not
      // residing on the basic device

      {  // basic primitive
        // point translated to the ECS origin, but axes direction global
        Point3D localPP;
        {  // Rotate point from global to local system
          double InitialVector[3] = {xfld - xpsrc, yfld - ypsrc, zfld - zpsrc};
          double FinalVector[3] = {0., 0., 0.};
          for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
              FinalVector[i] += TransformationMatrix[i][j] * InitialVector[j];
            }
          }
          localPP.X = FinalVector[0];
          localPP.Y = FinalVector[1];
          localPP.Z = FinalVector[2];
        }  // Point3D rotated

        // evaluate possibility whether primitive influence is accurate enough
        // This could be based on localPP and the subtended solid angle
        // If 1, then only primitive influence will be considered
        int PrimOK = 0;
        // consider primitive representation accurate enough if it is
        // repeated and beyond WtFldPrimAfter repetitions.
        if (WtFldPrimAfter < 0) { // If WtFldPrimAfter <0, PrimOK is zero
          PrimOK = 0;
        } else if (WtFldPrimAfter == 0) {  // only this is necessary
          PrimOK = 1;
        } else if (WtFldPrimAfter > 0) {
          PrimOK = 1;
        }
        if (PrimOK) {
          // Potential and flux (local system) due to base primitive
          double tmpPot;
          Vector3D tmpF;
          GetPrimPF(primsrc, &localPP, &tmpPot, &tmpF);
          const double qpr = AvWtChDen[IdWtField][primsrc];
          pPot[primsrc] += qpr * tmpPot;
          lFx += qpr * tmpF.X;
          lFy += qpr * tmpF.Y;
          lFz += qpr * tmpF.Z;
          // if(DebugLevel == 301)
          if (dbgFn) {
            printf("PFAtPoint base primitive =>\n");
            printf("primsrc: %d, xlocal: %lg, ylocal: %lg, zlocal %lg\n",
                   primsrc, localPP.X, localPP.Y, localPP.Z);
            printf("primsrc: %d, Pot: %lg, Fx: %lg, Fx: %lg, Fz: %lg\n",
                   primsrc, tmpPot, tmpF.X, tmpF.Y, tmpF.Z);
            printf("primsrc: %d, pPot: %lg, lFx: %lg, lFy: %lg, lFz: %lg\n",
                   primsrc, pPot[primsrc], lFx, lFy, lFz);
            fflush(stdout);
            // exit(-1);
          }
        } else {
          // element influence
          double tPot;
          Vector3D tF;
          double ePot = 0.;
          Vector3D eF;
          eF.X = 0.0;
          eF.Y = 0.0;
          eF.Z = 0.0;

          const int eleMin = ElementBgn[primsrc];
          const int eleMax = ElementEnd[primsrc];
          for (int ele = eleMin; ele <= eleMax; ++ele) {
            const double xsrc = (EleArr + ele - 1)->G.Origin.X;
            const double ysrc = (EleArr + ele - 1)->G.Origin.Y;
            const double zsrc = (EleArr + ele - 1)->G.Origin.Z;
            // Rotate vector from global to local system; matrix as for
            // primitive
            double vG[3] = {xfld - xsrc, yfld - ysrc, zfld - zsrc};
            double vL[3] = {0., 0., 0.};
            for (int i = 0; i < 3; ++i) {
              for (int j = 0; j < 3; ++j) {
                vL[i] += TransformationMatrix[i][j] * vG[j];
              }
            }
            // Potential and flux (local system) due to base primitive
            const int type = (EleArr + ele - 1)->G.Type;
            const double a = (EleArr + ele - 1)->G.LX;
            const double b = (EleArr + ele - 1)->G.LZ;
            GetPF(type, a, b, vL[0], vL[1], vL[2], &tPot, &tF);
            const double qel = WtFieldChDen[IdWtField][ele];
            ePot += qel * tPot;
            eF.X += qel * tF.X;
            eF.Y += qel * tF.Y;
            eF.Z += qel * tF.Z;
            // if(DebugLevel == 301)
            if (dbgFn) {
              printf("PFAtPoint base primitive:%d\n", primsrc);
              printf("ele: %d, xlocal: %lg, ylocal: %lg, zlocal %lg\n", ele,
                     vL[0], vL[1], vL[2]);
              printf(
                  "ele: %d, tPot: %lg, tFx: %lg, tFy: %lg, tFz: %lg, Solution: "
                  "%g\n",
                  ele, tPot, tF.X, tF.Y, tF.Z, qel);
              printf("ele: %d, ePot: %lg, eFx: %lg, eFy: %lg, eFz: %lg\n", ele,
                     ePot, eF.X, eF.Y, eF.Z);
              fflush(stdout);
            }
          }  // for all the elements on this primsrc primitive

          pPot[primsrc] += ePot;
          lFx += eF.X;
          lFy += eF.Y;
          lFz += eF.Z;
          if (dbgFn) {
            printf("prim%d, ePot: %lg, eFx: %lg, eFy: %lg, eFz: %lg\n", primsrc,
                   ePot, eF.X, eF.Y, eF.Z);
            printf("prim%d, pPot: %lg, lFx: %lg, lFy: %lg, lFz: %lg\n", primsrc,
                   pPot[primsrc], lFx, lFy, lFz);
            fflush(stdout);
          }
        }  // else elements influence

        // if(DebugLevel == 301)
        if (dbgFn) {
          printf("basic primtive\n");
          printf("primsrc: %d, pPot: %lg, lFx: %lg, lFy: %lg, lFz: %lg\n",
                 primsrc, pPot[primsrc], lFx, lFy, lFz);
          fflush(stdout);
        }
      }  // basic primitive ends

      if (MirrorTypeX[primsrc] || MirrorTypeY[primsrc] ||
          MirrorTypeZ[primsrc]) {  // Mirror effect of base primitives
        printf("Mirror may not be correctly implemented ...\n");
        exit(0);
      }  // Mirror effect ends

      // Flux due to repeated primitives
      if ((PeriodicTypeX[primsrc] == 1) || (PeriodicTypeY[primsrc] == 1) ||
          (PeriodicTypeZ[primsrc] == 1)) {
        const int perx = PeriodicInX[primsrc];
        const int pery = PeriodicInY[primsrc];
        const int perz = PeriodicInZ[primsrc];
        if (perx || pery || perz) {
          for (int xrpt = -perx; xrpt <= perx; ++xrpt) {
            const double xShift = XPeriod[primsrc] * (double)xrpt;
            double XPOfRpt = xpsrc + xShift;
            for (int yrpt = -pery; yrpt <= pery; ++yrpt) {
              const double yShift = YPeriod[primsrc] * (double)yrpt;
              double YPOfRpt = ypsrc + yShift;
              for (int zrpt = -perz; zrpt <= perz; ++zrpt) {
                const double zShift = ZPeriod[primsrc] * (double)zrpt;
                double ZPOfRpt = zpsrc + zShift;
                // Skip the base device.
                if ((xrpt == 0) && (yrpt == 0) && (zrpt == 0)) continue;
                {  // basic primitive repeated
                  Point3D localPPR;
                  {  // Rotate point from global to local system
                    double InitialVector[3] = {xfld - XPOfRpt, yfld - YPOfRpt,
                                               zfld - ZPOfRpt};
                    double FinalVector[3] = {0., 0., 0.};
                    for (int i = 0; i < 3; ++i) {
                      for (int j = 0; j < 3; ++j) {
                        FinalVector[i] +=
                            TransformationMatrix[i][j] * InitialVector[j];
                      }
                    }
                    localPPR.X = FinalVector[0];
                    localPPR.Y = FinalVector[1];
                    localPPR.Z = FinalVector[2];
                  }  // Point3D rotated

                  int repPrimOK = 0;

                  // consider primitive representation accurate enough if it is
                  // repeated and beyond WtFldPrimAfter repetitions.
                  if (WtFldPrimAfter < 0) {//WtFldPrimAfter <0 => repPrimOK = 0
                    repPrimOK = 0;
                  } else if ((abs(xrpt) >= WtFldPrimAfter) &&
                             (abs(yrpt) >= WtFldPrimAfter)) {
                    repPrimOK = 1;
                  }
                  // enforce primitive representation since it is unlikely
                  // that the weighting field will be modified much due to
                  // such an approximation for the repeated primitives.
                  if (repPrimOK) {  // use primitive representation
                    // Potential and flux (local system) due to repeated
                    // primitive
                    double tmpPot;
                    Vector3D tmpF;
                    GetPrimPF(primsrc, &localPPR, &tmpPot, &tmpF);
                    const double qpr = AvWtChDen[IdWtField][primsrc];
                    pPot[primsrc] += qpr * tmpPot;
                    lFx += qpr * tmpF.X;
                    lFy += qpr * tmpF.Y;
                    lFz += qpr * tmpF.Z;
                    // if(DebugLevel == 301)
                    if (dbgFn) {
                      printf(
                          "primsrc: %d, xlocal: %lg, ylocal: %lg, zlocal: "
                          "%lg\n",
                          primsrc, localPPR.X, localPPR.Y, localPPR.Z);
                      printf(
                          "primsrc: %d, Pot: %lg, Fx: %lg, Fy: %lg, Fz: %lg\n",
                          primsrc, tmpPot * qpr, tmpF.X * qpr, tmpF.Y * qpr,
                          tmpF.Z * qpr);
                      printf(
                          "primsrc: %d, pPot: %lg, lFx: %lg, lFy: %lg, lFz: "
                          "%lg\n",
                          primsrc, pPot[primsrc], lFx, lFy, lFz);
                      fflush(stdout);
                    }
                  } else {
                    // use discretized representation of a repeated primitive
                    double tPot;
                    Vector3D tF;
                    double erPot = 0.;
                    Vector3D erF;
                    erF.X = 0.0;
                    erF.Y = 0.0;
                    erF.Z = 0.0;

                    const int eleMin = ElementBgn[primsrc];
                    const int eleMax = ElementEnd[primsrc];
                    for (int ele = eleMin; ele <= eleMax; ++ele) {
                      const double xrsrc = (EleArr + ele - 1)->G.Origin.X;
                      const double yrsrc = (EleArr + ele - 1)->G.Origin.Y;
                      const double zrsrc = (EleArr + ele - 1)->G.Origin.Z;

                      const double XEOfRpt = xrsrc + xShift;
                      const double YEOfRpt = yrsrc + yShift;
                      const double ZEOfRpt = zrsrc + zShift;

                      // Rotate point from global to local system.
                      double vG[3] = {xfld - XEOfRpt, yfld - YEOfRpt,
                                      zfld - ZEOfRpt};
                      double vL[3] = {0., 0., 0.};
                      for (int i = 0; i < 3; ++i) {
                        for (int j = 0; j < 3; ++j) {
                          vL[i] += TransformationMatrix[i][j] * vG[j];
                        }
                      }
                      // Allowed, because all the local coordinates have the
                      // same orientations. Only the origins are mutually
                      // displaced along a line.
                      const int type = (EleArr + ele - 1)->G.Type;
                      const double a = (EleArr + ele - 1)->G.LX;
                      const double b = (EleArr + ele - 1)->G.LZ;
                      GetPF(type, a, b, vL[0], vL[1], vL[2], &tPot, &tF);
                      const double qel = WtFieldChDen[IdWtField][ele];
                      erPot += qel * tPot;
                      erF.X += qel * tF.X;
                      erF.Y += qel * tF.Y;
                      erF.Z += qel * tF.Z;
                      // if(DebugLevel == 301)
                      if (dbgFn) {
                        printf("PFAtPoint base primitive:%d\n", primsrc);
                        printf(
                            "ele: %d, xlocal: %lg, ylocal: %lg, zlocal %lg\n",
                            ele, vL[0], vL[1], vL[2]);
                        printf(
                            "ele: %d, tPot: %lg, tFx: %lg, tFy: %lg, tFz: %lg, "
                            "Solution: %g\n",
                            ele, tPot, tF.X, tF.Y, tF.Z, qel);
                        printf(
                            "ele: %d, ePot: %lg, eFx: %lg, eFy: %lg, eFz: "
                            "%lg\n",
                            ele, erPot, erF.X, erF.Y, erF.Z);
                        fflush(stdout);
                      }
                    }  // for all the elements on this primsrc repeated
                       // primitive

                    pPot[primsrc] += erPot;
                    lFx += erF.X;
                    lFy += erF.Y;
                    lFz += erF.Z;
                  }  // else discretized representation of this primitive

                  // if(DebugLevel == 301)
                  if (dbgFn) {
                    printf("basic repeated xrpt: %d. yrpt: %d, zrpt: %d\n",
                           xrpt, yrpt, zrpt);
                    printf(
                        "primsrc: %d, pPot: %lg, lFx: %lg, lFy: %lg, lFz: "
                        "%lg\n",
                        primsrc, pPot[primsrc], lFx, lFy, lFz);
                    fflush(stdout);
                  }
                }  // repetition of basic primitive

                if (MirrorTypeX[primsrc] || MirrorTypeY[primsrc] ||
                    MirrorTypeZ[primsrc]) {
                  // Mirror effect of repeated primitives - not parallelized
                  printf(
                      "Mirror not correctly implemented in this version of "
                      "neBEM ...\n");
                  exit(0);
                }  // Mirror effect for repeated primitives ends

              }  // for zrpt
            }    // for yrpt
          }      // for xrpt
        }        // PeriodicInX || PeriodicInY || PeriodicInZ
      }          // PeriodicType == 1
      Vector3D localF;
      localF.X = lFx;
      localF.Y = lFy;
      localF.Z = lFz;
      Vector3D tmpF = RotateVector3D(&localF, &PrimDC[primsrc], local2global);
      plFx[primsrc] = tmpF.X;  // local fluxes lFx, lFy, lFz in GCS
      plFy[primsrc] = tmpF.Y;
      plFz[primsrc] = tmpF.Z;
    }  // for all primitives: basic device, mirror reflections and repetitions
  }    // pragma omp parallel

  double totPot = 0.0;
  Vector3D totF;
  totF.X = totF.Y = totF.Z = 0.0;
  for (int prim = 1; prim <= NbPrimitives; ++prim) {
    totPot += pPot[prim];
    totF.X += plFx[prim];
    totF.Y += plFy[prim];
    totF.Z += plFz[prim];
  }

  // This should be done at the end of the function - before freeing memory
#ifdef __cplusplus
  *Potential = totPot * InvFourPiEps0;
  globalF->X = totF.X * InvFourPiEps0;
  globalF->Y = totF.Y * InvFourPiEps0;
  globalF->Z = totF.Z * InvFourPiEps0;
#else
  *Potential = totPot / MyFACTOR;
  globalF->X = totF.X / MyFACTOR;
  globalF->Y = totF.Y / MyFACTOR;
  globalF->Z = totF.Z / MyFACTOR;
#endif

  /*
  For weighting field, effect of KnCh is possibly zero.
  Similarly, there is no reason to respect constraint on total system charge.
  */

  if (dbgFn) {
    printf("Final values due to all primitives and other influences: ");
    // printf("xfld\tyfld\tzfld\tPot\tFx\tFy\tFz\n");	// refer, do not
    // uncomment
    printf("%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n\n", xfld, yfld, zfld,
           (*Potential), globalF->X, globalF->Y, globalF->Z);
    fflush(stdout);
  }

  free_dvector(pPot, 1, NbPrimitives);
  free_dvector(plFx, 1, NbPrimitives);
  free_dvector(plFy, 1, NbPrimitives);
  free_dvector(plFz, 1, NbPrimitives);

  return (0);
}  // end of WtFldPFAtPoint

// Compute weighting potential and field in a given Fast Volume
// Possible pitfall: evaluation of n-skips
int CreateWtFldFastVolPF(int IdWtField) {
  if (IdWtField >= MAXWtFld) {
    printf(
        "neBEMPrepareWeightingField: reached MaxWtField (%d) weighting "
        "fields.\n",
        MAXWtFld);
    return -1;
  }

  int dbgFn = 0;
  int fstatus;

  int nbXCells;
  int nbYCells;
  int nbZCells;
  double startX;
  double startY;
  double startZ;
  double delX;
  double delY;
  double delZ;

  printf("\nComputing weighting potential & field within basic fast volume\n");
  int bskip = 0, iskip = 0, jskip = 0, kskip = 0;

  // calculate n-skips based on NbPtSkip
  if (NbPtSkip) {
    int volptcnt = 0, endskip = 0;

    for (int block = 1; block <= WtFldFastVol[IdWtField].NbBlocks; ++block) {
      nbXCells = WtFldBlkNbXCells[IdWtField][block];
      nbYCells = WtFldBlkNbYCells[IdWtField][block];
      nbZCells = WtFldBlkNbZCells[IdWtField][block];
      for (int i = 1; i <= nbXCells + 1; ++i) {
        for (int j = 1; j <= nbYCells + 1; ++j) {
          for (int k = 1; k <= nbZCells + 1; ++k) {
            ++volptcnt;

            if (volptcnt == WtFldNbPtSkip[IdWtField]) {
              bskip = block - 1;
              iskip = i - 1;
              jskip = j - 1;
              kskip = k;
              endskip = 1;
            }

            if (endskip) break;
          }
          if (endskip) break;
        }
        if (endskip) break;
      }
      if (endskip) break;
    }
    if (dbgFn) {
      printf(
          "Basic fast volume => bskip, iskip, jskip, kskip: %d, %d, %d, %d\n",
          bskip, iskip, jskip, kskip);
    }
  }  // WtFldNbPtSkip

  // stringify the integer
  char stringIdWtField[16];
  sprintf(stringIdWtField, "%d", IdWtField);

  char WtFldFastVolPFFile[256];
  strcpy(WtFldFastVolPFFile, BCOutDir);
  strcat(WtFldFastVolPFFile, "/WtFldFastVolPF_");
  strcat(WtFldFastVolPFFile, stringIdWtField);
  strcat(WtFldFastVolPFFile, ".out");
  FILE *fWtFldFastVolPF = fopen(WtFldFastVolPFFile, "w");
  if (fWtFldFastVolPF == NULL) {
    neBEMMessage("CreateWtFldFastVolPF - WtFldFastVolPFFile");
    return -1;
  }
  fprintf(fWtFldFastVolPF, "#block\tX\tY\tZ\tPot\tFX\tFY\tFZ\n");

  if (dbgFn) {
    printf("WtFldFastVolPF.out created ...\n");
    fflush(stdout);
  }

  for (int block = 1 + bskip; block <= WtFldFastVol[IdWtField].NbBlocks;
       ++block) {
    nbXCells = WtFldBlkNbXCells[IdWtField][block];
    nbYCells = WtFldBlkNbYCells[IdWtField][block];
    nbZCells = WtFldBlkNbZCells[IdWtField][block];
    startX = WtFldFastVol[IdWtField].CrnrX;
    startY = WtFldFastVol[IdWtField].CrnrY;
    startZ = WtFldBlkCrnrZ[IdWtField][block];
    delX = WtFldFastVol[IdWtField].LX / nbXCells;
    delY = WtFldFastVol[IdWtField].LY / nbYCells;
    delZ = WtFldBlkLZ[IdWtField][block] / nbZCells;
    printf(
        "WtFldNbBlocks: %d, block: %d, nbXCells: %d, nbYCells: %d, nbZCells: "
        "%d\n",
        WtFldFastVol[IdWtField].NbBlocks, block, nbXCells, nbYCells, nbZCells);

    if (dbgFn) {
      printf("block: %d\n", block);
      printf("nbXCells, nbYCells, nbZCells: %d, %d, %d\n", nbXCells, nbYCells,
             nbZCells);
      printf("startX, startY, startZ: %le, %le, %le\n", startX, startY, startZ);
      printf("delX, delY, delZ: %le, %le, %le\n", delX, delY, delZ);
      printf("bskip, iskip, jskip, kskip: %d, %d, %d, %d\n", bskip, iskip,
             jskip, kskip);
      fflush(stdout);
    }
    // total number of points in a given block
    // int blktotpt = (nbXCells + 1) * (nbYCells + 1) * (nbZCells + 1);
    for (int i = 1 + iskip; i <= nbXCells + 1; ++i) {
      for (int j = 1 + jskip; j <= nbYCells + 1; ++j) {
        printf("Fast volume => block: %3d, i: %4d, j: %4d", block, i, j);
        fflush(stdout);

        Point3D point;
#ifdef _OPENMP
        int nthreads = 1, tid = 0;
#pragma omp parallel private(nthreads, tid)
#endif
        {
#ifdef _OPENMP
          if (dbgFn) {
            tid = omp_get_thread_num();
            if (tid == 0) {
              nthreads = omp_get_num_threads();
              printf("Starting fast volume computation with %d threads\n",
                     nthreads);
            }
          }
#endif
          int k;
          int omitFlag;
          double potential;
          Vector3D field;
#ifdef _OPENMP
#pragma omp for private(k, point, omitFlag, potential, field)
#endif
          for (k = 1 + kskip; k <= nbZCells + 1; ++k) {
            potential = 0.0;
            field.X = field.Y = field.Z = 0.0;

            point.X = startX + (i - 1) * delX;
            point.Y = startY + (j - 1) * delY;
            point.Z = startZ + (k - 1) * delZ;

            // Check whether the point falls within a volume that should be
            // ignored
            omitFlag = 0;
            for (int omit = 1; omit <= WtFldFastVol[IdWtField].NbOmitVols;
                 ++omit) {
              if ((point.X > WtFldOmitVolCrnrX[IdWtField][omit]) &&
                  (point.X < WtFldOmitVolCrnrX[IdWtField][omit] +
                                 WtFldOmitVolLX[IdWtField][omit]) &&
                  (point.Y > WtFldOmitVolCrnrY[IdWtField][omit]) &&
                  (point.Y < WtFldOmitVolCrnrY[IdWtField][omit] +
                                 WtFldOmitVolLY[IdWtField][omit]) &&
                  (point.Z > WtFldOmitVolCrnrZ[IdWtField][omit]) &&
                  (point.Z < WtFldOmitVolCrnrZ[IdWtField][omit] +
                                 WtFldOmitVolLZ[IdWtField][omit])) {
                omitFlag = 1;
                break;
              }
            }  // loop over omitted volumes

            if (dbgFn) {
              printf("block, i, j, k: %d, %d, %d, %d\n", block, i, j, k);
              printf("point X, Y, Z: %.8lg\t%.8lg\t%.8lg\n",
                     point.X / LengthScale, point.Y / LengthScale,
                     point.Z / LengthScale);
              printf("omitFlag: %d\n", omitFlag);
              fflush(stdout);
            }

            if (omitFlag) {
              potential = field.X = field.Y = field.Z = 0.0;
            } else {
              fstatus = WtFldPFAtPoint(&point, &potential, &field, IdWtField);
              if (fstatus != 0) {
                neBEMMessage("wrong return from WtFldPFAtPoint.\n");
                // return -1;
              }
            } // else omitFlag
            if (dbgFn) {
              printf("%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\n",
                     point.X / LengthScale, point.Y / LengthScale,
                     point.Z / LengthScale, potential / LengthScale, field.X,
                     field.Y, field.Z);
              fflush(stdout);
            }

            WtFldFastPot[IdWtField][block][i][j][k] = potential;
            WtFldFastFX[IdWtField][block][i][j][k] = field.X;
            WtFldFastFY[IdWtField][block][i][j][k] = field.Y;
            WtFldFastFZ[IdWtField][block][i][j][k] = field.Z;
          }  // loop k
        }    // pragma omp parallel

        for (int k = 1 + kskip; k <= nbZCells + 1; ++k)  // file output
        {
          point.X = startX + (i - 1) * delX;
          point.Y = startY + (j - 1) * delY;
          point.Z = startZ + (k - 1) * delZ;

          fprintf(fWtFldFastVolPF,
                  "%4d\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\n",
                  block, point.X / LengthScale, point.Y / LengthScale,
                  point.Z / LengthScale,
                  WtFldFastPot[IdWtField][block][i][j][k],
                  WtFldFastFX[IdWtField][block][i][j][k],
                  WtFldFastFY[IdWtField][block][i][j][k],
                  WtFldFastFZ[IdWtField][block][i][j][k]);
        }
        fflush(fWtFldFastVolPF);  // file output over

        printf(
            "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
            "\b\b\b\b\b\b\b\b\b\b");
      }  // loop j
    }    // loop i
  }      // loop block

  fclose(fWtFldFastVolPF);

  if (OptStaggerWtFldFastVol[IdWtField]) {
    printf(
        "\nComputing weighting potential & field within staggered fast "
        "volume\n");

    bskip = iskip = jskip = kskip = 0;

    // calculate n-skips based on StgWtFldNbPtSkip
    if (StgWtFldNbPtSkip[IdWtField]) {
      int volptcnt = 0, endskip = 0;

      for (int block = 1; block <= WtFldFastVol[IdWtField].NbBlocks; ++block) {
        nbXCells = WtFldBlkNbXCells[IdWtField][block];
        nbYCells = WtFldBlkNbYCells[IdWtField][block];
        nbZCells = WtFldBlkNbZCells[IdWtField][block];
        for (int i = 1; i <= nbXCells + 1; ++i) {
          for (int j = 1; j <= nbYCells + 1; ++j) {
            for (int k = 1; k <= nbZCells + 1; ++k) {
              ++volptcnt;

              if (volptcnt == StgWtFldNbPtSkip[IdWtField]) {
                bskip = block - 1;
                iskip = i - 1;
                jskip = j - 1;
                kskip = k;
                endskip = 1;
              }

              if (endskip) break;
            }
            if (endskip) break;
          }
          if (endskip) break;
        }
        if (endskip) break;
      }
      if (dbgFn) {
        printf(
            "Staggered volume => bskip, iskip, jskip, kskip: %d, %d, %d, %d\n",
            bskip, iskip, jskip, kskip);
      }
    }  // StgWtFldNbPtSkip

    char StgWtFldFastVolPFFile[256];
    strcpy(StgWtFldFastVolPFFile, BCOutDir);
    strcat(StgWtFldFastVolPFFile, "/StgWtFldFastVolPF_");
    strcat(StgWtFldFastVolPFFile, stringIdWtField);
    strcat(StgWtFldFastVolPFFile, ".out");
    FILE *fStgWtFldFastVolPF = fopen(StgWtFldFastVolPFFile, "w");
    if (fStgWtFldFastVolPF == NULL) {
      neBEMMessage("WtFldFastVolPF - StgWtFldFastVolPFFile");
      return -1;
    }
    fprintf(fStgWtFldFastVolPF, "#block\tX\tY\tZ\tPot\tFX\tFY\tFZ\n");

    if (dbgFn) {
      printf("StgWtFldFastVolPF.out created ...\n");
      fflush(stdout);
    }

    for (int block = 1 + bskip; block <= WtFldFastVol[IdWtField].NbBlocks;
         ++block) {
      nbXCells = WtFldBlkNbXCells[IdWtField][block];
      nbYCells = WtFldBlkNbYCells[IdWtField][block];
      nbZCells = WtFldBlkNbZCells[IdWtField][block];
      startX = WtFldFastVol[IdWtField].CrnrX + WtFldFastVol[IdWtField].LX;
      startY = WtFldFastVol[IdWtField].CrnrY + WtFldFastVol[IdWtField].YStagger;
      startZ = WtFldBlkCrnrZ[IdWtField][block];
      delX = WtFldFastVol[IdWtField].LX / nbXCells;
      delY = WtFldFastVol[IdWtField].LY / nbYCells;
      delZ = WtFldBlkLZ[IdWtField][block] / nbZCells;
      printf(
          "NbBlocks: %d, block: %d, nbXCells: %d, nbYCells: %d, nbZCells: %d\n",
          WtFldFastVol[IdWtField].NbBlocks, block, nbXCells, nbYCells,
          nbZCells);

      if (dbgFn) {
        printf("block: %d\n", block);
        printf("nbXCells, nbYCells, nbZCells: %d, %d, %d\n", nbXCells, nbYCells,
               nbZCells);
        printf("startX, startY, startZ: %le, %le, %le\n", startX, startY,
               startZ);
        printf("delX, delY, delZ: %le, %le, %le\n", delX, delY, delZ);
        printf("bskip, iskip, jskip, kskip: %d, %d, %d, %d\n", bskip, iskip,
               jskip, kskip);
        fflush(stdout);
      }

      // int blktotpt = (nbXCells + 1) * (nbYCells + 1) * (nbZCells + 1);
      for (int i = 1 + iskip; i <= nbXCells + 1; ++i) {
        for (int j = 1 + jskip; j <= nbYCells + 1; ++j) {
          printf("Fast volume => block: %3d, i: %4d, j: %4d", block, i, j);
          fflush(stdout);

          Point3D point;
#ifdef _OPENMP
          int nthreads = 1, tid = 0;
#pragma omp parallel private(nthreads, tid)
#endif
          {
#ifdef _OPENMP
            if (dbgFn) {
              tid = omp_get_thread_num();
              if (tid == 0) {
                nthreads = omp_get_num_threads();
                printf(
                    "Starting staggered fast volume computation with %d "
                    "threads\n",
                    nthreads);
              }
            }
#endif
            int k;
            int omitFlag;
            double potential;
            Vector3D field;
#ifdef _OPENMP
#pragma omp for private(k, point, omitFlag, potential, field)
#endif
            for (k = 1 + kskip; k <= nbZCells + 1; ++k) {
              potential = 0.0;
              field.X = field.Y = field.Z = 0.0;

              point.X = startX + (i - 1) * delX;
              point.Y = startY + (j - 1) * delY;
              point.Z = startZ + (k - 1) * delZ;

              // Check whether point falls within a volume that should be
              // ignored
              omitFlag = 0;
              for (int omit = 1; omit <= WtFldFastVol[IdWtField].NbOmitVols;
                   ++omit) {
                if ((point.X > WtFldOmitVolCrnrX[IdWtField][omit] +
                                   WtFldFastVol[IdWtField].LX) &&
                    (point.X < WtFldOmitVolCrnrX[IdWtField][omit] +
                                   WtFldOmitVolLX[IdWtField][omit] +
                                   WtFldFastVol[IdWtField].LX) &&
                    (point.Y > WtFldOmitVolCrnrY[IdWtField][omit] +
                                   WtFldFastVol[IdWtField].YStagger) &&
                    (point.Y < WtFldOmitVolCrnrY[IdWtField][omit] +
                                   WtFldOmitVolLY[IdWtField][omit] +
                                   WtFldFastVol[IdWtField].YStagger) &&
                    (point.Z > WtFldOmitVolCrnrZ[IdWtField][omit]) &&
                    (point.Z < WtFldOmitVolCrnrZ[IdWtField][omit] +
                                   WtFldOmitVolLZ[IdWtField][omit])) {
                  omitFlag = 1;
                  break;
                }
              }  // loop over omitted volumes

              if (dbgFn) {
                printf("point X, Y, Z: %.8lg\t%.8lg\t%.8lg\n",
                       point.X / LengthScale, point.Y / LengthScale,
                       point.Z / LengthScale);
                printf("omitFlag: %d\n", omitFlag);
                fflush(stdout);
              }

              if (omitFlag) {
                potential = field.X = field.Y = field.Z = 0.0;
              } else {
                fstatus = WtFldPFAtPoint(&point, &potential, &field, IdWtField);
                if (fstatus != 0) {
                  neBEMMessage(
                      "wrong WtFldPFAtPoint return value in staggered part of "
                      "CreateWtFldFastVolElePF.\n");
                  // return -1;
                }
              } // else omitFlag
              if (dbgFn) {
                printf("%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\n",
                       point.X / LengthScale, point.Y / LengthScale,
                       point.Z / LengthScale, potential / LengthScale, field.X,
                       field.Y, field.Z);
                fflush(stdout);
              }

              StgWtFldFastPot[IdWtField][block][i][j][k] = potential;
              StgWtFldFastFX[IdWtField][block][i][j][k] = field.X;
              StgWtFldFastFY[IdWtField][block][i][j][k] = field.Y;
              StgWtFldFastFZ[IdWtField][block][i][j][k] = field.Z;
            }  // loop k
          }    // pragma omp

          for (int k = 1 + kskip; k <= nbZCells + 1; ++k) {
            // file output
            point.X = startX + (i - 1) * delX;
            point.Y = startY + (j - 1) * delY;
            point.Z = startZ + (k - 1) * delZ;

            fprintf(fStgWtFldFastVolPF,
                    "%4d\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\t%.8lg\n",
                    block, point.X / LengthScale, point.Y / LengthScale,
                    point.Z / LengthScale,
                    StgWtFldFastPot[IdWtField][block][i][j][k],
                    StgWtFldFastFX[IdWtField][block][i][j][k],
                    StgWtFldFastFY[IdWtField][block][i][j][k],
                    StgWtFldFastFZ[IdWtField][block][i][j][k]);
          }
          fflush(fStgWtFldFastVolPF);  // file output over

          printf(
              "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
              "\b\b\b\b\b\b\b\b\b\b\b");
        }  // loop j
      }    // loop i
    }      // loop block

    fclose(fStgWtFldFastVolPF);
  }  // if OptStaggerWtFldFastVol

  return 0;
}  // CreateWtFldFastVol ends

// Gives three components of the total Potential and flux in the global
// coordinate system due to all the elements using the results stored in
// the FAST volume mesh. The Fast volume is generated in the normal manner
// but by making necessary changes in the boundary conditions. This Fast
// volume is then renamed. The same is true for the data in staggered volume.
// These names are provided to the code by the neBEMWtFldFastVol.inp
int WtFldFastPFAtPoint(Point3D *globalP, double *Potential, Vector3D *globalF,
                       int IdWtField) {
  if (IdWtField >= MAXWtFld) {
    printf(
        "neBEMPrepareWeightingField: reached MaxWtField (%d) weighting "
        "fields.\n",
        MAXWtFld);
    return -1;
  }

  int dbgFn = 0;
  double Xpt = globalP->X;
  double Ypt = globalP->Y;
  double Zpt = globalP->Z;
  double RptVolLX = WtFldFastVol[IdWtField].LX;
  double RptVolLY = WtFldFastVol[IdWtField].LY;
  double RptVolLZ = WtFldFastVol[IdWtField].LZ;
  double CornerX = WtFldFastVol[IdWtField].CrnrX;
  double CornerY = WtFldFastVol[IdWtField].CrnrY;
  double CornerZ = WtFldFastVol[IdWtField].CrnrZ;
  double TriLin(double xd, double yd, double zd, double c000, double c100,
                double c010, double c001, double c110, double c101, double c011,
                double c111);

  // First of all, check how the point in question should be treated ...

  // Check whether the point falls within a volume that is not regarded as
  // FastVol
  for (int ignore = 1; ignore <= WtFldFastVol[IdWtField].NbIgnoreVols;
       ++ignore) {
    if ((Xpt >= (WtFldIgnoreVolCrnrX[IdWtField][ignore])) &&
        (Xpt <= (WtFldIgnoreVolCrnrX[IdWtField][ignore] +
                 WtFldIgnoreVolLX[IdWtField][ignore])) &&
        (Ypt >= (WtFldIgnoreVolCrnrY[IdWtField][ignore])) &&
        (Ypt <= (WtFldIgnoreVolCrnrY[IdWtField][ignore] +
                 WtFldIgnoreVolLY[IdWtField][ignore])) &&
        (Zpt >= (WtFldIgnoreVolCrnrZ[IdWtField][ignore])) &&
        (Zpt <= (WtFldIgnoreVolCrnrZ[IdWtField][ignore] +
                 WtFldIgnoreVolLZ[IdWtField][ignore]))) {
      if (dbgFn)
        neBEMMessage("In WtFldFastPFAtPoint: point in an ignored volume!\n");

      // KnCh does not have any effect
      int fstatus = ElePFAtPoint(globalP, Potential, globalF);
      if (fstatus != 0) {
        neBEMMessage("wrong WtFldPFAtPoint return value in FastVolPF.\n");
        return -1;
      } else
        return 0;
    }
  }  // loop over ignored volumes

  // If not ignored, the point qualifies for FastVol evaluation ...

  // for a staggered fast volume, the volume repeated in X is larger
  if (OptStaggerWtFldFastVol[IdWtField]) {
    RptVolLX += WtFldFastVol[IdWtField].LX;
  }
  if (dbgFn) {
    printf("\nin WtFldFastPFAtPoint\n");
    printf("x, y, z: %g, %g, %g\n", Xpt, Ypt, Zpt);
    printf("RptVolLX, RptVolLY, RptVolLZ: %g, %g, %g\n", RptVolLX, RptVolLY,
           RptVolLZ);
    printf("CornerX, CornerY, CornerZ: %g, %g, %g\n", CornerX, CornerY,
           CornerZ);
    printf("Nb of blocks: %d\n", WtFldFastVol[IdWtField].NbBlocks);
    for (int block = 1; block <= WtFldFastVol[IdWtField].NbBlocks; ++block) {
      printf("NbOfXCells: %d\n", WtFldBlkNbXCells[IdWtField][block]);
      printf("NbOfYCells: %d\n", WtFldBlkNbYCells[IdWtField][block]);
      printf("NbOfZCells: %d\n", WtFldBlkNbZCells[IdWtField][block]);
      printf("LZ: %le\n", WtFldBlkLZ[IdWtField][block]);
      printf("CornerZ: %le\n", WtFldBlkCrnrZ[IdWtField][block]);
    }
  }

  // Find equivalent position inside the basic / staggered volume.
  // real distance from volume corner
  double dx = Xpt - CornerX;
  double dy = Ypt - CornerY;
  double dz = Zpt - CornerZ;
  if (dbgFn)
    printf("real dx, dy, dz from volume corner: %g, %g, %g\n", dx, dy, dz);

  int NbFastVolX = (int)(dx / RptVolLX);
  if (dx < 0.0) --NbFastVolX;
  int NbFastVolY = (int)(dy / RptVolLY);
  if (dy < 0.0) --NbFastVolY;
  int NbFastVolZ = (int)(dz / RptVolLZ);
  if (dz < 0.0) --NbFastVolZ;
  if (dbgFn)
    printf("Volumes in x, y, z: %d, %d, %d\n", NbFastVolX, NbFastVolY,
           NbFastVolZ);

  // equivalent distances from fast volume corner
  dx -= NbFastVolX * RptVolLX;
  dy -= NbFastVolY * RptVolLY;
  dz -= NbFastVolZ * RptVolLZ;
  // The following conditions should never happen - generate an error message
  if (dx < 0.0) {
    dx = 0.0;
    neBEMMessage("equiv dx < 0.0 - not correct!\n");
  }
  if (dy < 0.0) {
    dy = 0.0;
    neBEMMessage("equiv dy < 0.0 - not correct!\n");
  }
  if (dz < 0.0) {
    dz = 0.0;
    neBEMMessage("equiv dz < 0.0 - not correct!\n");
  }
  if (dx > RptVolLX) {
    dx = RptVolLX;
    neBEMMessage("equiv dx > RptVolLX - not correct!\n");
  }
  if (dy > RptVolLY) {
    dy = RptVolLY;
    neBEMMessage("equiv dy > RptVolLY - not correct!\n");
  }
  if (dz > RptVolLZ) {
    dz = RptVolLZ;
    neBEMMessage("equiv dz > RptVolLZ - not correct!\n");
  }
  if (dbgFn)
    printf("equivalent dist from corner - dx, dy, dz: %g, %g, %g\n", dx, dy,
           dz);

  // Take care of possible trouble-makers
  if (dx < MINDIST) dx = MINDIST;  // -ve dx has already been made equal to 0
  if (dy < MINDIST) dy = MINDIST;
  if (dz < MINDIST) dz = MINDIST;
  if ((RptVolLX - dx) < MINDIST)
    dx = RptVolLX - MINDIST;  // dx > RptVolLX taken care of
  if ((RptVolLY - dy) < MINDIST) dy = RptVolLY - MINDIST;
  if ((RptVolLZ - dz) < MINDIST) dz = RptVolLZ - MINDIST;
  // For staggered volumes, there is another plane where difficulties may occur
  if ((dx <= WtFldFastVol[IdWtField].LX) &&
      (WtFldFastVol[IdWtField].LX - dx) < MINDIST)
    dx = WtFldFastVol[IdWtField].LX - MINDIST;
  else if ((dx > WtFldFastVol[IdWtField].LX) &&
           (fabs(WtFldFastVol[IdWtField].LX - dx) < MINDIST))
    dx = WtFldFastVol[IdWtField].LX + MINDIST;
  if (dbgFn)
    printf("equivalent dist adjusted - dx, dy, dz: %g, %g, %g\n", dx, dy, dz);

  // If volume is staggered, we have a few more things to do before finalizing
  // the values of equivalent distance
  // sector identification
  // _................__________________
  // |    .     .     |    Sector 3    |
  // |     .   .      |                |
  // |       .        |       .        |
  // |                |     .   .      |
  // |    Sector 2    |    .     .     |
  // |----------------|     .   .      |
  // |                |       .        |
  // |       .        |                |
  // |     .   .      |                |
  // |    .     .     |----------------|
  // |     .   .      |    Sector 4    |
  // |       .        |       .        |
  // |                |     .   .      |
  // |    Sector 1    |    .     .     |
  // |----------------|................|

  int sector = 1;  // kept outside `if' since this is necessary further below
  if (OptStaggerWtFldFastVol[IdWtField]) {
    if ((dx >= 0.0) && (dx <= WtFldFastVol[IdWtField].LX) && (dy >= 0.0) &&
        (dy <= WtFldFastVol[IdWtField].LY)) {
      // point lies in sector 1, everything remains unchanged
      sector = 1;
    } else if ((dx >= 0.0) && (dx <= WtFldFastVol[IdWtField].LX) &&
               (dy > WtFldFastVol[IdWtField].LY) &&
               (dy <= WtFldFastVol[IdWtField].LY +
                          WtFldFastVol[IdWtField].YStagger)) {
      // point lies in sector 2, move basic volume one step up
      sector = 2;
      ++NbFastVolY;
      CornerY += WtFldFastVol[IdWtField].LY;  // repeat length in Y is LY
      dy -= WtFldFastVol[IdWtField].LY;
    } else if ((dx > WtFldFastVol[IdWtField].LX) &&
               (dx <= 2.0 * WtFldFastVol[IdWtField].LX) &&
               (dy >= WtFldFastVol[IdWtField].YStagger) &&
               (dy <= WtFldFastVol[IdWtField].LY +
                          WtFldFastVol[IdWtField].YStagger)) {
      // point lies in sector 3, pt in staggered vol, change corner coords
      sector = 3;
      CornerX += WtFldFastVol[IdWtField].LX;
      CornerY += WtFldFastVol[IdWtField].YStagger;
      dx -= WtFldFastVol[IdWtField].LX;
      dy -= WtFldFastVol[IdWtField].YStagger;
    } else if ((dx > WtFldFastVol[IdWtField].LX) &&
               (dx <= 2.0 * WtFldFastVol[IdWtField].LX) && (dy >= 0.0) &&
               (dy < WtFldFastVol[IdWtField].YStagger)) {
      // point lies in sector 4, move basic volume one step down and consider
      // staggered fast volume
      sector = 4;
      --NbFastVolY;
      CornerX += WtFldFastVol[IdWtField]
                     .LX;  // in the staggered part of the repeated volume
      CornerY -=
          (WtFldFastVol[IdWtField].LY - WtFldFastVol[IdWtField].YStagger);
      dx -= WtFldFastVol[IdWtField].LX;
      dy += (WtFldFastVol[IdWtField].LY - WtFldFastVol[IdWtField].YStagger);
    } else {
      neBEMMessage("WtFldFastPFAtPoint: point in none of the sectors!\n");
    }
    if (dbgFn) printf("stagger modified dx, dy, dz: %g, %g, %g\n", dx, dy, dz);
  }

  // Take care of possible trouble-makers - once more
  if (dx < MINDIST) dx = MINDIST;  // -ve dx has already been made equal to 0
  if (dy < MINDIST) dy = MINDIST;
  if (dz < MINDIST) dz = MINDIST;
  if ((RptVolLX - dx) < MINDIST)
    dx = RptVolLX - MINDIST;  // dx > RptVolLX taken care of
  if ((RptVolLY - dy) < MINDIST) dy = RptVolLY - MINDIST;
  if ((RptVolLZ - dz) < MINDIST) dz = RptVolLZ - MINDIST;
  // For staggered volumes, there is another plane where difficulties may occur
  if ((dx <= WtFldFastVol[IdWtField].LX) &&
      (WtFldFastVol[IdWtField].LX - dx) < MINDIST)
    dx = WtFldFastVol[IdWtField].LX - MINDIST;
  else if ((dx > WtFldFastVol[IdWtField].LX) &&
           (fabs(WtFldFastVol[IdWtField].LX - dx) < MINDIST))
    dx = WtFldFastVol[IdWtField].LX + MINDIST;
  if (dbgFn)
    printf("equivalent dist adjusted for staggered: %g, %g, %g\n", dx, dy, dz);

  // Find the block in which the point lies
  int thisBlock = 0;
  for (int block = 1; block <= WtFldFastVol[IdWtField].NbBlocks; ++block) {
    double blkBtmZ =
        WtFldBlkCrnrZ[IdWtField][block] - CornerZ;  // since CornerZ has been
    double blkTopZ =
        blkBtmZ + WtFldBlkLZ[IdWtField][block];  // subtracted from dz already
    if (dbgFn) {
      printf("block, dz, blkBtmZ, blkTopZ: %d, %lg, %lg, %lg\n", block, dz,
             blkBtmZ, blkTopZ);
    }

    // take care of difficult situations
    if ((dz <= blkBtmZ) && ((blkBtmZ - dz) < MINDIST)) dz = blkBtmZ - MINDIST;
    if ((dz >= blkBtmZ) && ((dz - blkBtmZ) < MINDIST)) dz = blkBtmZ + MINDIST;
    if ((dz <= blkTopZ) && ((blkTopZ - dz) < MINDIST)) dz = blkTopZ - MINDIST;
    if ((dz >= blkTopZ) && ((dz - blkTopZ) < MINDIST)) dz = blkTopZ + MINDIST;

    if ((dz >= blkBtmZ) && (dz <= blkTopZ)) {
      thisBlock = block;
      break;
    }
  }
  if (!thisBlock) {
    neBEMMessage("WtFldFastPFAtPoint: point in none of the blocks!\n");
  }

  int nbXCells = WtFldBlkNbXCells[IdWtField][thisBlock];
  int nbYCells = WtFldBlkNbYCells[IdWtField][thisBlock];
  int nbZCells = WtFldBlkNbZCells[IdWtField][thisBlock];
  double delX = WtFldFastVol[IdWtField].LX / nbXCells;
  double delY = WtFldFastVol[IdWtField].LY / nbYCells;
  double delZ = WtFldBlkLZ[IdWtField][thisBlock] / nbZCells;
  dz -= (WtFldBlkCrnrZ[IdWtField][thisBlock] -
         CornerZ);  // distance from the block corner

  if (dbgFn) {
    printf("thisBlock: %d\n", thisBlock);
    printf("nbXCells, nbYCells, nbZCells: %d, %d, %d\n", nbXCells, nbYCells,
           nbZCells);
    printf("WtFldBlkCrnrZ: %lg\n", WtFldBlkCrnrZ[IdWtField][thisBlock]);
    printf("delX, delY, delZ: %le, %le, %le\n", delX, delY, delZ);
    printf("dz: %lg\n", dz);
    fflush(stdout);
  }

  // Find cell in block of basic / staggered volume within which the point lies
  int celli = (int)(dx / delX) + 1;  // Find cell in which the point lies
  if (celli < 1) {
    celli = 1;
    dx = 0.5 * delX;
    neBEMMessage("WtFldFastPFAtPoint - celli < 1\n");
  }
  if (celli > nbXCells) {
    celli = nbXCells;
    dx = WtFldFastVol[IdWtField].LX - 0.5 * delX;
    neBEMMessage("WtFldFastPFAtPoint - celli > nbXCells\n");
  }
  int cellj = (int)(dy / delY) + 1;
  if (cellj < 1) {
    cellj = 1;
    dy = 0.5 * delY;
    neBEMMessage("WtFldFastPFAtPoint - cellj < 1\n");
  }
  if (cellj > nbYCells) {
    cellj = nbYCells;
    dy = WtFldFastVol[IdWtField].LY - 0.5 * delY;
    neBEMMessage("WtFldFastPFAtPoint - cellj > nbYCells\n");
  }
  int cellk = (int)(dz / delZ) + 1;
  if (cellk < 1) {
    cellk = 1;
    dz = 0.5 * delX;
    neBEMMessage("WtFldFastPFAtPoint - cellk < 1\n");
  }
  if (cellk > nbZCells) {
    cellk = nbZCells;
    dz = WtFldFastVol[IdWtField].LZ - 0.5 * delZ;
    neBEMMessage("WtFldFastPFAtPoint - cellk > nbZCells\n");
  }
  if (dbgFn) printf("Cells in x, y, z: %d, %d, %d\n", celli, cellj, cellk);

  // Interpolate potential and field at the point using the corner values of
  // of the cell and, if necessary, of the neighbouring cells
  // These gradients can also be calculated while computing the potential and
  // field at the cells and stored in memory, provided enough memory is
  // available

  // distances from cell corner
  double dxcellcrnr = dx - (double)(celli - 1) * delX;
  double dycellcrnr = dy - (double)(cellj - 1) * delY;
  double dzcellcrnr = dz - (double)(cellk - 1) * delZ;
  if (dbgFn)
    printf("cell crnr dx, dy, dz: %g, %g, %g\n", dxcellcrnr, dycellcrnr,
           dzcellcrnr);

  // normalized distances
  double xd = dxcellcrnr / delX;  // xd = (x-x0)/(x1-x0)
  double yd = dycellcrnr / delY;  // etc
  double zd = dzcellcrnr / delZ;
  if (xd <= 0.0) xd = 0.0;
  if (yd <= 0.0) yd = 0.0;
  if (zd <= 0.0) zd = 0.0;
  if (xd >= 1.0) xd = 1.0;
  if (yd >= 1.0) yd = 1.0;
  if (zd >= 1.0) zd = 1.0;

  // corner values of potential and field
  double P000 =
      WtFldFastPot[IdWtField][thisBlock][celli][cellj][cellk];  // lowest corner
  double FX000 = WtFldFastFX[IdWtField][thisBlock][celli][cellj][cellk];
  double FY000 = WtFldFastFY[IdWtField][thisBlock][celli][cellj][cellk];
  double FZ000 = WtFldFastFZ[IdWtField][thisBlock][celli][cellj][cellk];
  double P100 = WtFldFastPot[IdWtField][thisBlock][celli + 1][cellj][cellk];
  double FX100 = WtFldFastFX[IdWtField][thisBlock][celli + 1][cellj][cellk];
  double FY100 = WtFldFastFY[IdWtField][thisBlock][celli + 1][cellj][cellk];
  double FZ100 = WtFldFastFZ[IdWtField][thisBlock][celli + 1][cellj][cellk];
  double P010 = WtFldFastPot[IdWtField][thisBlock][celli][cellj + 1][cellk];
  double FX010 = WtFldFastFX[IdWtField][thisBlock][celli][cellj + 1][cellk];
  double FY010 = WtFldFastFY[IdWtField][thisBlock][celli][cellj + 1][cellk];
  double FZ010 = WtFldFastFZ[IdWtField][thisBlock][celli][cellj + 1][cellk];
  double P001 = WtFldFastPot[IdWtField][thisBlock][celli][cellj][cellk + 1];
  double FX001 = WtFldFastFX[IdWtField][thisBlock][celli][cellj][cellk + 1];
  double FY001 = WtFldFastFY[IdWtField][thisBlock][celli][cellj][cellk + 1];
  double FZ001 = WtFldFastFZ[IdWtField][thisBlock][celli][cellj][cellk + 1];
  double P110 = WtFldFastPot[IdWtField][thisBlock][celli + 1][cellj + 1][cellk];
  double FX110 = WtFldFastFX[IdWtField][thisBlock][celli + 1][cellj + 1][cellk];
  double FY110 = WtFldFastFY[IdWtField][thisBlock][celli + 1][cellj + 1][cellk];
  double FZ110 = WtFldFastFZ[IdWtField][thisBlock][celli + 1][cellj + 1][cellk];
  double P101 = WtFldFastPot[IdWtField][thisBlock][celli + 1][cellj][cellk + 1];
  double FX101 = WtFldFastFX[IdWtField][thisBlock][celli + 1][cellj][cellk + 1];
  double FY101 = WtFldFastFY[IdWtField][thisBlock][celli + 1][cellj][cellk + 1];
  double FZ101 = WtFldFastFZ[IdWtField][thisBlock][celli + 1][cellj][cellk + 1];
  double P011 = WtFldFastPot[IdWtField][thisBlock][celli][cellj + 1][cellk + 1];
  double FX011 = WtFldFastFX[IdWtField][thisBlock][celli][cellj + 1][cellk + 1];
  double FY011 = WtFldFastFY[IdWtField][thisBlock][celli][cellj + 1][cellk + 1];
  double FZ011 = WtFldFastFZ[IdWtField][thisBlock][celli][cellj + 1][cellk + 1];
  double P111 =
      WtFldFastPot[IdWtField][thisBlock][celli + 1][cellj + 1][cellk + 1];
  double FX111 =
      WtFldFastFX[IdWtField][thisBlock][celli + 1][cellj + 1][cellk + 1];
  double FY111 =
      WtFldFastFY[IdWtField][thisBlock][celli + 1][cellj + 1][cellk + 1];
  double FZ111 =
      WtFldFastFZ[IdWtField][thisBlock][celli + 1][cellj + 1][cellk + 1];
  if (OptStaggerWtFldFastVol[IdWtField]) {
    if (sector == 1) {  // nothing to be done
    }
    if (sector == 2) {  // volume shifted up but point not in the staggered part
    }
    if (sector == 3) {  // staggered volume
      P000 = StgWtFldFastPot[IdWtField][thisBlock][celli][cellj][cellk];
      FX000 = StgWtFldFastFX[IdWtField][thisBlock][celli][cellj][cellk];
      FY000 = StgWtFldFastFY[IdWtField][thisBlock][celli][cellj][cellk];
      FZ000 = StgWtFldFastFZ[IdWtField][thisBlock][celli][cellj][cellk];
      P100 = StgWtFldFastPot[IdWtField][thisBlock][celli + 1][cellj][cellk];
      FX100 = StgWtFldFastFX[IdWtField][thisBlock][celli + 1][cellj][cellk];
      FY100 = StgWtFldFastFY[IdWtField][thisBlock][celli + 1][cellj][cellk];
      FZ100 = StgWtFldFastFZ[IdWtField][thisBlock][celli + 1][cellj][cellk];
      P010 = StgWtFldFastPot[IdWtField][thisBlock][celli][cellj + 1][cellk];
      FX010 = StgWtFldFastFX[IdWtField][thisBlock][celli][cellj + 1][cellk];
      FY010 = StgWtFldFastFY[IdWtField][thisBlock][celli][cellj + 1][cellk];
      FZ010 = StgWtFldFastFZ[IdWtField][thisBlock][celli][cellj + 1][cellk];
      P001 = StgWtFldFastPot[IdWtField][thisBlock][celli][cellj][cellk + 1];
      FX001 = StgWtFldFastFX[IdWtField][thisBlock][celli][cellj][cellk + 1];
      FY001 = StgWtFldFastFY[IdWtField][thisBlock][celli][cellj][cellk + 1];
      FZ001 = StgWtFldFastFZ[IdWtField][thisBlock][celli][cellj][cellk + 1];
      P110 = StgWtFldFastPot[IdWtField][thisBlock][celli + 1][cellj + 1][cellk];
      FX110 = StgWtFldFastFX[IdWtField][thisBlock][celli + 1][cellj + 1][cellk];
      FY110 = StgWtFldFastFY[IdWtField][thisBlock][celli + 1][cellj + 1][cellk];
      FZ110 = StgWtFldFastFZ[IdWtField][thisBlock][celli + 1][cellj + 1][cellk];
      P101 = StgWtFldFastPot[IdWtField][thisBlock][celli + 1][cellj][cellk + 1];
      FX101 = StgWtFldFastFX[IdWtField][thisBlock][celli + 1][cellj][cellk + 1];
      FY101 = StgWtFldFastFY[IdWtField][thisBlock][celli + 1][cellj][cellk + 1];
      FZ101 = StgWtFldFastFZ[IdWtField][thisBlock][celli + 1][cellj][cellk + 1];
      P011 = StgWtFldFastPot[IdWtField][thisBlock][celli][cellj + 1][cellk + 1];
      FX011 = StgWtFldFastFX[IdWtField][thisBlock][celli][cellj + 1][cellk + 1];
      FY011 = StgWtFldFastFY[IdWtField][thisBlock][celli][cellj + 1][cellk + 1];
      FZ011 = StgWtFldFastFZ[IdWtField][thisBlock][celli][cellj + 1][cellk + 1];
      P111 = StgWtFldFastPot[IdWtField][thisBlock][celli + 1][cellj + 1]
                            [cellk + 1];
      FX111 =
          StgWtFldFastFX[IdWtField][thisBlock][celli + 1][cellj + 1][cellk + 1];
      FY111 =
          StgWtFldFastFY[IdWtField][thisBlock][celli + 1][cellj + 1][cellk + 1];
      FZ111 =
          StgWtFldFastFZ[IdWtField][thisBlock][celli + 1][cellj + 1][cellk + 1];
    }
    if (sector == 4) {  // volume shifted down and point in the staggered part
      P000 = StgWtFldFastPot[IdWtField][thisBlock][celli][cellj][cellk];
      FX000 = StgWtFldFastFX[IdWtField][thisBlock][celli][cellj][cellk];
      FY000 = StgWtFldFastFY[IdWtField][thisBlock][celli][cellj][cellk];
      FZ000 = StgWtFldFastFZ[IdWtField][thisBlock][celli][cellj][cellk];
      P100 = StgWtFldFastPot[IdWtField][thisBlock][celli + 1][cellj][cellk];
      FX100 = StgWtFldFastFX[IdWtField][thisBlock][celli + 1][cellj][cellk];
      FY100 = StgWtFldFastFY[IdWtField][thisBlock][celli + 1][cellj][cellk];
      FZ100 = StgWtFldFastFZ[IdWtField][thisBlock][celli + 1][cellj][cellk];
      P010 = StgWtFldFastPot[IdWtField][thisBlock][celli][cellj + 1][cellk];
      FX010 = StgWtFldFastFX[IdWtField][thisBlock][celli][cellj + 1][cellk];
      FY010 = StgWtFldFastFY[IdWtField][thisBlock][celli][cellj + 1][cellk];
      FZ010 = StgWtFldFastFZ[IdWtField][thisBlock][celli][cellj + 1][cellk];
      P001 = StgWtFldFastPot[IdWtField][thisBlock][celli][cellj][cellk + 1];
      FX001 = StgWtFldFastFX[IdWtField][thisBlock][celli][cellj][cellk + 1];
      FY001 = StgWtFldFastFY[IdWtField][thisBlock][celli][cellj][cellk + 1];
      FZ001 = StgWtFldFastFZ[IdWtField][thisBlock][celli][cellj][cellk + 1];
      P110 = StgWtFldFastPot[IdWtField][thisBlock][celli + 1][cellj + 1][cellk];
      FX110 = StgWtFldFastFX[IdWtField][thisBlock][celli + 1][cellj + 1][cellk];
      FY110 = StgWtFldFastFY[IdWtField][thisBlock][celli + 1][cellj + 1][cellk];
      FZ110 = StgWtFldFastFZ[IdWtField][thisBlock][celli + 1][cellj + 1][cellk];
      P101 = StgWtFldFastPot[IdWtField][thisBlock][celli + 1][cellj][cellk + 1];
      FX101 = StgWtFldFastFX[IdWtField][thisBlock][celli + 1][cellj][cellk + 1];
      FY101 = StgWtFldFastFY[IdWtField][thisBlock][celli + 1][cellj][cellk + 1];
      FZ101 = StgWtFldFastFZ[IdWtField][thisBlock][celli + 1][cellj][cellk + 1];
      P011 = StgWtFldFastPot[IdWtField][thisBlock][celli][cellj + 1][cellk + 1];
      FX011 = StgWtFldFastFX[IdWtField][thisBlock][celli][cellj + 1][cellk + 1];
      FY011 = StgWtFldFastFY[IdWtField][thisBlock][celli][cellj + 1][cellk + 1];
      FZ011 = StgWtFldFastFZ[IdWtField][thisBlock][celli][cellj + 1][cellk + 1];
      P111 = StgWtFldFastPot[IdWtField][thisBlock][celli + 1][cellj + 1]
                            [cellk + 1];
      FX111 =
          StgWtFldFastFX[IdWtField][thisBlock][celli + 1][cellj + 1][cellk + 1];
      FY111 =
          StgWtFldFastFY[IdWtField][thisBlock][celli + 1][cellj + 1][cellk + 1];
      FZ111 =
          StgWtFldFastFZ[IdWtField][thisBlock][celli + 1][cellj + 1][cellk + 1];
    }
  }

  double intP =
      TriLin(xd, yd, zd, P000, P100, P010, P001, P110, P101, P011, P111);
  double intFX = TriLin(xd, yd, zd, FX000, FX100, FX010, FX001, FX110, FX101,
                        FX011, FX111);
  double intFY = TriLin(xd, yd, zd, FY000, FY100, FY010, FY001, FY110, FY101,
                        FY011, FY111);
  double intFZ = TriLin(xd, yd, zd, FZ000, FZ100, FZ010, FZ001, FZ110, FZ101,
                        FZ011, FZ111);

  *Potential = intP;
  globalF->X = intFX;
  globalF->Y = intFY;
  globalF->Z = intFZ;

  if (dbgFn) {
    printf("WtFldCell corner values:\n");
    printf("Potential: %g, %g, %g, %g\n", P000, P100, P010, P001);
    printf("Potential: %g, %g, %g, %g\n", P110, P101, P011, P111);
    printf("FastFX: %g, %g, %g, %g\n", FX000, FX100, FX010, FX001);
    printf("FastFX: %g, %g, %g, %g\n", FX110, FX101, FX011, FX111);
    printf("FastFY: %g, %g, %g, %g\n", FY000, FY100, FY010, FY001);
    printf("FastFY: %g, %g, %g, %g\n", FY110, FY101, FY011, FY111);
    printf("FastFZ: %g, %g, %g, %g\n", FZ000, FZ100, FZ010, FZ001);
    printf("FastFZ: %g, %g, %g, %g\n", FZ110, FZ101, FZ011, FZ111);
    printf("Pot, FX, FY, FZ: %g, %g, %g, %g\n", *Potential, globalF->X,
           globalF->Y, globalF->Z);
  }

  if (dbgFn) {
    printf("out WtFldFastPFAtPoint\n");
    fflush(stdout);
  }

  return 0;
}  // WtFldFastPFAtPoint ends

// Potential and flux per unit charge density on an element returned as
// Pot, Fx, Fy, Fz components
// in the global coordinate system
void GetPFGCS(int type, double a, double b, Point3D *localP, double *Potential,
              Vector3D *globalF, DirnCosn3D *DirCos) {
  Vector3D localF;
  const double x = localP->X;
  const double y = localP->Y;
  const double z = localP->Z;
  switch (type) {
    case 4:  // rectangular element
      RecPF(a, b, x, y, z, Potential, &localF);
      break;
    case 3:  // triangular element
      TriPF(a, b, x, y, z, Potential, &localF);
      break;
    case 2:  // linear (wire) element
      WirePF(a, b, x, y, z, Potential, &localF);
      break;
    default:
      printf("Geometrical type out of range! ... exiting ...\n");
      exit(-1);
      break;  // never comes here
  }           // switch over gtsrc ends

  (*globalF) = RotateVector3D(&localF, DirCos, local2global);
}  // end of GetPFGCS

// Potential and flux per unit charge density on an element returned as
// Pot, Fx, Fy, Fz components in the local coordinate system
void GetPF(int type, double a, double b, double x, double y, double z,
           double *Potential, Vector3D *localF) {
  switch (type) {
    case 4:  // rectangular element
      RecPF(a, b, x, y, z, Potential, localF);
      break;
    case 3:  // triangular element
      TriPF(a, b, x, y, z, Potential, localF);
      break;
    case 2:  // linear (wire) element
      WirePF(a, b, x, y, z, Potential, localF);
      break;
    default:
      printf("Geometrical type out of range! ... exiting ...\n");
      exit(-1);
      break;  // never comes here
  }           // switch over gtsrc ends
}  // end of GetPF

// Flux per unit charge density on a rectangular element
// Are X and Z directions the same as obtained using the direction cosines?
void RecPF(double a, double b, double x, double y, double z, double *Potential,
           Vector3D *localF) {
  const double d2 = x * x + y * y + z * z;
  if (d2 >= FarField2 * (a * a + b * b)) {
    (*Potential) = a * b / sqrt(d2);
    const double f = (*Potential) / d2;
    localF->X = x * f;
    localF->Y = y * f;
    localF->Z = z * f;
  } else {
    int fstatus = ExactRecSurf(x / a, y / a, z / a, -0.5, -(b / a) / 2.0, 0.5,
                               (b / a) / 2.0, Potential, localF);
    if (fstatus) {  // non-zero
      printf("problem in RecPF ... \n");
      // printf("returning ...\n");
      // return -1; void function at present
    }
    (*Potential) *= a;  // rescale - cannot be done outside because of the `if'
  }

}  // end of RecPF

// Flux per unit charge density on a triangular element
void TriPF(double a, double b, double x, double y, double z, double *Potential,
           Vector3D *localF) {
  // printf("In TriPF\n");
  const double xm = x - a / 3.;
  const double zm = z - b / 3.;
  const double d2 = xm * xm + y * y + zm * zm;

  if (d2 >= FarField2 * (a * a + b * b)) {
    (*Potential) = 0.5 * a * b / sqrt(d2);
    const double f = (*Potential) / d2;
    localF->X = x * f;
    localF->Y = y * f;
    localF->Z = z * f;
  } else {
    int fstatus = ExactTriSurf(b / a, x / a, y / a, z / a, Potential, localF);
    if (fstatus) {  // non-zero
      printf("problem in TriPF ... \n");
      // printf("returning ...\n");
      // return -1; void function at present
    }
    (*Potential) *= a;  // rescale - cannot be done outside because of the `if'
  }
  // printf("Out of TriPF\n");
}  // end of TriPF

// Flux per unit charge density on a wire element
void WirePF(double rW, double lW, double x, double y, double z,
            double *Potential, Vector3D *localF) {
  const double d2 = x * x + y * y + z * z;
  if (d2 >= FarField2 * lW * lW) {
    double dA = 2.0 * ST_PI * rW * lW;
    const double dist = sqrt(d2);
    (*Potential) = dA / dist;
    double f = dA / (dist * d2);
    localF->X = x * f;
    localF->Y = y * f;
    localF->Z = z * f;
  } else {
    if ((fabs(x) < MINDIST) && (fabs(y) < MINDIST)) {
      if (fabs(z) < MINDIST)
        (*Potential) = ExactCentroidalP_W(rW, lW);
      else
        (*Potential) = ExactAxialP_W(rW, lW, z);

      localF->X = localF->Y = 0.0;
      localF->Z = ExactThinFZ_W(rW, lW, x, y, z);
    } else {
      ExactThinWire(rW, lW, x, y, z, Potential, localF);
    }
  }
}  // end of WirePF

// Potential and flux per unit charge density on an element returned as
// Pot, Fx, Fy, Fz components
// in the global coordiante system
void GetPrimPFGCS(int prim, Point3D *localP, double *Potential,
                  Vector3D *globalF, DirnCosn3D *DirCos) {
  Vector3D localF;

  switch (PrimType[prim]) {
    case 4:  // rectangular primitive
      RecPrimPF(prim, localP, Potential, &localF);
      break;
    case 3:  // triangular primitive
      TriPrimPF(prim, localP, Potential, &localF);
      break;
    case 2:  // linear (wire) primitive
      WirePrimPF(prim, localP, Potential, &localF);
      break;
    default:
      printf("Geometrical type out of range! ... exiting ...\n");
      exit(-1);
      break;  // never comes here
  }           // switch over gtsrc ends

  (*globalF) = RotateVector3D(&localF, DirCos, local2global);
}  // end of GetPrimPFGCS

// Potential and flux per unit charge density on an element returned as
// Pot, Fx, Fy, Fz components in the local coordinate system
void GetPrimPF(int prim, Point3D *localP, double *Potential, Vector3D *localF) {
  switch (PrimType[prim]) {
    case 4:  // rectangular primitive
      RecPrimPF(prim, localP, Potential, localF);
      break;
    case 3:  // triangular primitive
      TriPrimPF(prim, localP, Potential, localF);
      break;
    case 2:  // linear (wire) primitive
      WirePrimPF(prim, localP, Potential, localF);
      break;
    default:
      printf("Geometrical type out of range! ... exiting ...\n");
      exit(-1);
      break;  // never comes here
  }           // switch over gtsrc ends
}  // end of GetPrimPF

// Flux per unit charge density on a rectangular primitive
void RecPrimPF(int prim, Point3D *localP, double *Potential, Vector3D *localF) {
  double xpt = localP->X;
  double ypt = localP->Y;
  double zpt = localP->Z;
  double dist = sqrt(xpt * xpt + ypt * ypt + zpt * zpt);

  double a = PrimLX[prim];
  double b = PrimLZ[prim];
  double diag = sqrt(a * a + b * b);  // diagonal

  if (dist >= FarField * diag) {
    double dA = a * b;  // area
    (*Potential) = dA / dist;
    const double f = dA / (dist * dist * dist);
    localF->X = xpt * f;
    localF->Y = ypt * f;
    localF->Z = zpt * f;
  } else {
    int fstatus = ExactRecSurf(xpt / a, ypt / a, zpt / a, -0.5, -(b / a) / 2.0,
                               0.5, (b / a) / 2.0, Potential, localF);
    if (fstatus) {  // non-zero
      printf("problem in RecPrimPF ... \n");
      // printf("returning ...\n");
      // return -1; void function at present
    }
    (*Potential) *= a;  // rescale - cannot be done outside because of the `if'
  }

}  // end of RecPrimPF

// Flux per unit charge density on a triangluar primitive
// Note that vertex[1] is the right angle corner
void TriPrimPF(int prim, Point3D *localP, double *Potential, Vector3D *localF) {
  double xpt = localP->X;
  double ypt = localP->Y;
  double zpt = localP->Z;
  double a = PrimLX[prim];
  double b = PrimLZ[prim];
  // longest side (hypotenuse)
  double diag = sqrt(a * a + b * b);
  const double xm = xpt - a / 3.;
  const double zm = zpt - b / 3.;
  double dist = sqrt(xm * xm + ypt * ypt + zm * zm);

  if (dist >= FarField * diag) {
    double dA = 0.5 * a * b;  // area
    (*Potential) = dA / dist;
    double f = dA / (dist * dist * dist);
    localF->X = xpt * f;
    localF->Y = ypt * f;
    localF->Z = zpt * f;
  } else {
    int fstatus =
        ExactTriSurf(b / a, xpt / a, ypt / a, zpt / a, Potential, localF);
    if (fstatus) { // non-zero
      printf("problem in TriPrimPF ... \n");
    }
    (*Potential) *= a;  // rescale - cannot be done outside because of the `if'
  }

  // printf("Out of TriPrimPF\n");
}  // end of TriPrimPF

// Flux per unit charge density on a wire primitive
void WirePrimPF(int prim, Point3D *localP, double *Potential,
                Vector3D *localF) {
  double xpt = localP->X;
  double ypt = localP->Y;
  double zpt = localP->Z;
  double rW = Radius[prim];
  double lW = PrimLZ[prim];
  double dist = sqrt(xpt * xpt + ypt * ypt + zpt * zpt);  

  if (dist >= FarField * lW) {
    double dA = 2.0 * ST_PI * rW * lW;
    (*Potential) = dA / dist;
    double f = dA / (dist * dist * dist);
    localF->X = xpt * f;
    localF->Y = ypt * f;
    localF->Z = zpt * f;
  } else {
    if ((fabs(xpt) < MINDIST) && (fabs(ypt) < MINDIST)) {
      if (fabs(zpt) < MINDIST)
        (*Potential) = ExactCentroidalP_W(rW, lW);
      else
        (*Potential) = ExactAxialP_W(rW, lW, zpt);

      localF->X = localF->Y = 0.0;
      localF->Z = ExactThinFZ_W(rW, lW, xpt, ypt, zpt);
    } else {
      ExactThinWire(rW, lW, xpt, ypt, zpt, Potential, localF);
    }
  }
}  // end of WirePrimPF

double TriLin(double xd, double yd, double zd, double c000, double c100,
              double c010, double c001, double c110, double c101, double c011,
              double c111) {
  double c00 = c000 * (1.0 - xd) + c100 * xd;
  double c10 = c010 * (1.0 - xd) + c110 * xd;
  double c01 = c001 * (1.0 - xd) + c101 * xd;
  double c11 = c011 * (1.0 - xd) + c111 * xd;
  double c0 = c00 * (1.0 - yd) + c10 * yd;
  double c1 = c01 * (1.0 - yd) + c11 * yd;
  return (c0 * (1.0 - zd) + c1 * zd);
}

#ifdef __cplusplus
}  // namespace
#endif
