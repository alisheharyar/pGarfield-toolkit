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

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

#include "Isles.h"
#include "NR.h"
#include "Vector.h"
#include "neBEM.h"
#include "neBEMInterface.h"

#ifdef __cplusplus
namespace neBEM {
#endif

int InfluenceMatrixFlag;

// At the end of this function, one should have the charge density distribution
// This function is not called when only post-processing is being carried out,
// i.e., when NewModel == NewMesh == NewBC == 0.
// In such an event, the solution is directly read using ReadSolution
int ComputeSolution(void) {
  int LHMatrix(void), RHVector(void), Solve(void);
  int InvertMatrix(void);
  int ReadInvertedMatrix(void);
#ifdef _OPENMP
  double time_begin = 0., time_end = 0.;
#endif
  printf(
      "----------------------------------------------------------------\n\n");
  printf("ComputeSolution: neBEM solution begins ...\n");
  printf("                 TimeStep: %d\n", TimeStep);
  printf("                 NewModel: %d, NewMesh: %d, NewBC: %d, NewPP: %d\n",
         NewModel, NewMesh, NewBC, NewPP);
  printf(
      "                 ModelCntr: %d, MeshCntr: %d, BCCntr: %d, PPCntr: %d\n",
      ModelCntr, MeshCntr, BCCntr, PPCntr);
  fflush(stdout);
  DebugISLES = 0;  // an integer declared in Isles header file

  NbEqns = NbUnknowns = NbElements;

  switch (OptInvMatProc) {
    case 0:
      OptLU = 1;
      OptSVD = 0;
      OptGSL = 0;
      break;
    case 1:
      OptLU = 0;
      OptSVD = 1;
      OptGSL = 0;
      break;
    case 2:
      OptLU = 0;
      OptSVD = 0;
      OptGSL = 1;
      break;
    default:
      OptLU = 1;
      OptSVD = 0;
      OptGSL = 0;
  }
  if ((OptSVD == 0) && (OptLU == 0) && (OptGSL == 0)) {
    printf("Cannot proceed with OptSVD, OptLU and OptGSL zero.\n");
    printf("Assuming the safer option OptSVD = 1.\n");
    OptLU = 0;
    OptSVD = 1;
    OptGSL = 0;
  }
  if ((OptSVD == 1) && (OptLU == 1) && (OptGSL == 1)) {
    printf("Cannot proceed with all OptSVD, OptLU and OptGSL one.\n");
    printf("Assuming the safer option OptSVD = 1.\n");
    OptLU = 0;
    OptSVD = 1;
    OptGSL = 0;
  }

  NbConstraints = 0;
  if (OptSystemChargeZero && NbFloatingConductors) {
    printf(
        "ComputeSolution: Simultaneous presence of OptSystemChargeZero && "
        "NbFloatingConductors!\n");
    printf("                 Returning ...\n");
    return (-1);
  }
  if (OptSystemChargeZero)  // Constraint making total charge on the sysyem,
                            // zero
  {
    ++NbConstraints;
    NbUnknowns = NbElements + NbConstraints;
    NbEqns = NbElements + NbConstraints;
    NbSystemChargeZero =
        NbUnknowns;          // which equation & unknown relates to this
  }                          // constraint
  if (NbFloatingConductors)  // Number of floating conductors now restricted to
                             // one
  {
    if (NbFloatingConductors > 1) {
      printf("Number of floating conductors > 1! ... not yet implemented.\n");
      printf("Returning\n");
      return -1;
    }
    ++NbConstraints;
    NbUnknowns = NbElements + NbConstraints;
    NbEqns = NbElements + NbConstraints;
    NbFloatCon = NbUnknowns;  // which equation and unknown relates to this
  }                           // floating conductor

  if (NewModel || NewMesh) {
    OptValidateSolution = 1;
    InfluenceMatrixFlag = 1;
  } else {
    InfluenceMatrixFlag = 0;
  }
  // Computation of the influence coefficient matrix and inversion of the same
  // is necessary only when we are considering a NewModel, and / or a NewMesh,
  // i.e., InfluenceMatrixFlag is true.
  // When OptChargingUp is true, then influence matrix, if not already
  // available, is computed in the EffectChargingUp function. For an existing
  // model and unchanged mesh, we can simply read in the relevant inverted
  // matrix.
  if (TimeStep == 1) {
    if (InfluenceMatrixFlag) {
      startClock = clock();
      printf("ComputeSolution: LHMatrix ... ");
      fflush(stdout);
#ifdef _OPENMP
      time_begin = omp_get_wtime();
#endif
      int fstatus = LHMatrix();
#ifdef _OPENMP
      time_end = omp_get_wtime();
      printf("Elapsed time: %lg\n", time_end - time_begin);
#endif
      if (fstatus != 0) {
        neBEMMessage("ComputeSolution - LHMatrix");
        return -1;
      }
      printf("ComputeSolution: LHMatrix done!\n");
      fflush(stdout);
      stopClock = clock();
      neBEMTimeElapsed(startClock, stopClock);
      printf("to setup influence matrix.\n");

      startClock = clock();
      printf("ComputeSolution: Inverting influence matrix ...\n");
      fflush(stdout);
      fstatus = InvertMatrix();
      if (fstatus != 0) {
        neBEMMessage("ComputeSolution - InvertMatrix");
        return -1;
      }
      printf("ComputeSolution: Matrix inversion over.\n");
      stopClock = clock();
      neBEMTimeElapsed(startClock, stopClock);
      printf("to invert influence matrix.\n");
    }  // if InfluenceMatrixFlag

    if ((!InfluenceMatrixFlag) && NewBC) {
      if (OptReadInvMatrix) {
        startClock = clock();
        printf(
            "ComputeSolution: Reading inverted matrix ... will take time ...");
        int fstatus = ReadInvertedMatrix();
        if (fstatus != 0) {
          neBEMMessage("ComputeSolution - ReadInvertedMatrix");
          return -1;
        }
        printf("                 done!\n");
        stopClock = clock();
        neBEMTimeElapsed(startClock, stopClock);
        printf("to read inverted influence matrix.\n");
      } else {
        neBEMMessage("ComputeSolution - NewBC but no InvMat ... ");
        neBEMMessage("don't know how to proceed!\n");
        return -1;
      }
    }  // if (!InfluenceMatrixFlag) && NewBC
  }    // if TimeStep == 1

  // If TimeStep > 1, use the Infl (LHS) and InvMat existing in memory

  int fstatus = 0;

  // Update known charges
  if (OptKnCh) {
    startClock = clock();
    printf("ComputeSolution: UpdateKnownCharges ... ");
    fflush(stdout);
    fstatus = UpdateKnownCharges();
    if (fstatus != 0) {
      neBEMMessage("ComputeSolution - UpdateKnownCharges");
      return -1;
    }
    printf("ComputeSolution: UpdateKnownCharges done!\n");
    fflush(stdout);
    stopClock = clock();
    neBEMTimeElapsed(startClock, stopClock);
    printf("to set up UpdateKnownCharges.\n");
  }  // if OptKnCh

  // Update charging up
  if (OptChargingUp) {
    startClock = clock();
    printf("ComputeSolution: UpdateChargingUp ... ");
    fflush(stdout);
    fstatus = UpdateChargingUp();
    if (fstatus != 0) {
      neBEMMessage("ComputeSolution - UpdateChargingUp");
      return -1;
    }
    printf("ComputeSolution: UpdateChargingUp done!\n");
    fflush(stdout);
    stopClock = clock();
    neBEMTimeElapsed(startClock, stopClock);
    printf("to set up UpdateChargingUp.\n");
  }  // if OptChargingUp

  // RHS
  startClock = clock();
  printf("ComputeSolution: RHVector ... ");
  fflush(stdout);
  fstatus = RHVector();
  if (fstatus != 0) {
    neBEMMessage("ComputeSolution - RHVector");
    return -1;
  }
  printf("ComputeSolution: RHVector done!\n");
  fflush(stdout);
  stopClock = clock();
  neBEMTimeElapsed(startClock, stopClock);
  printf("to set up RH vector.\n");

  // Solve
  // The Solve routine simply involves the matrix
  // multiplication of inverted matrix and the current RHVector.
  startClock = clock();
  printf("ComputeSolution: Solve ... ");
  fflush(stdout);
#ifdef _OPENMP
  time_begin = omp_get_wtime();
#endif
  fstatus = Solve();
#ifdef _OPENMP
  time_end = omp_get_wtime();
  printf("Elapsed time: %lg\n", time_end - time_begin);
#endif
  if (fstatus != 0) {
    neBEMMessage("ComputeSolution - Solve");
    return -1;
  }
  printf("ComputeSolution: Solve done!\n");
  fflush(stdout);
  stopClock = clock();
  neBEMTimeElapsed(startClock, stopClock);
  printf("to compute solution.\n");

  // The other locations where this matrix is freed is in
  // InvertMatrix (if OptValidateSolution is false!),
  // or in
  // Solve (due to a combination of OptValidateSolution true and
  // InfluenceMatrixFlag false)
  // due to RepeatLHMatrix
  // or OptStoreInflMatrix.
  if (InfluenceMatrixFlag && OptValidateSolution && EndOfTime)
    free_dmatrix(Inf, 1, NbEqns, 1, NbUnknowns);

  printf("ComputeSolution: neBEM solution ends ...\a\n");
  fflush(stdout);
  return (0);
}  // end of neBEM

// The coordinates of and distances from the barycentre of the source element
// have been included in the computations and are being passed as parameters
// of the ComputeInfluence (and to several successivey called functions).
// This is because, while the influence of triangular elements are computed
// based on a local coordinate system that has its origin at the right-angle
// corner of the triangle, the approximate influence of these elements needs
// to be computed based on a bary-centric co-ordinate system of the source
// element. This additional computation and extra burden of passing parameters
// can be reduced if the exact influence due to triangular elements can be
// computed based on a bary-centric co-ordinate system.
// Note that the co-ordinate transformations, that require the direction
// cosines of the local co-ordinate system is always the one already defined
// for the element. This implies that the directions cosines are assumed to
// remain unchanged when shifted from the right-corner of the triangle to its
// bary-centre.
// The problem does not arise for a rectangular element since the exact
// influence due to a rectangular element is always computed in terms of its
// centroidal co-ordinate system.
// The same is true for wire elements - origin of the LCS and the centroid are
// collocated.
int LHMatrix(void) {
#ifdef _OPENMP
  int dbgFn = 0;
#endif
  printf(
      "\nLHMatrix: The size of the Influence coefficient matrix is %d X %d\n",
      NbEqns, NbUnknowns);
  fflush(stdout);

  // The influence coefficient matrix is created only when the
  // InfluenceMatrixFlag
  // is true. The other eventualities when this can get created is when
  // this flag is false but the OptValidateSolution flag is true (through
  // RepeatLHMatrix (here!)
  // or1s
  // OptStoreInflMatrix (in function Solve)).
  Inf = dmatrix(1, NbEqns, 1, NbUnknowns);

  // Influence coefficient matrix
  // For each field point where the boundary condition is known (collocation
  // point), influence from all the source elements needs to be summed up
  // The field points are followed using elefld (field counter) and the
  // source elements are followed using elesrc (source counter)
  // printf("field point: ");	// do not remove
  printf("Computing influence coefficient matrix ... will take time ...\n");

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
        printf("Starting influence matrix computation with %d threads\n",
               nthreads);
      }
    }
#endif
    // printf("Field point:");
    // fflush(stdout);

    for (int elefld = 1; elefld <= NbElements; ++elefld) {
      // printf("%6d", elefld);
      // fflush(stdout);	// do not remove

      // Retrieve element properties at the field point
      // boundary condn applied at collocation point
      const double xfld = (EleArr + elefld - 1)->BC.CollPt.X;
      const double yfld = (EleArr + elefld - 1)->BC.CollPt.Y;
      const double zfld = (EleArr + elefld - 1)->BC.CollPt.Z;

#ifdef _OPENMP
#pragma omp for
#endif
      for (int elesrc = 1; elesrc <= NbElements; ++elesrc) {
        if (DebugLevel == 301) {
          printf("\n\nelefld: %d, elesrc: %d\n", elefld, elesrc);
        }

        // Retrieve element properties at the field point
        const int primsrc = (EleArr + elesrc - 1)->PrimitiveNb;
        const double xsrc = (EleArr + elesrc - 1)->G.Origin.X;
        const double ysrc = (EleArr + elesrc - 1)->G.Origin.Y;
        const double zsrc = (EleArr + elesrc - 1)->G.Origin.Z;

        if ((EleArr + elesrc - 1)->E.Type == 0) {
          printf(
              "LHMatrix: Wrong EType for elesrc %d element %d on %dth "
              "primitive!\n",
              elesrc, (EleArr + elesrc - 1)->Id, primsrc);
          exit(-1);
        }

        // The total influence is due to elements on the basic device and due to
        // virtual elements arising out of repetition, reflection etc and not
        // residing on the basic device

        // Influence due to elements belonging to the basic device
        {
          Point3D localP;
          // Through InitialVector[], field point gets translated to ECS origin
          // Axes direction are, however, still global which when rotated to ECS
          // system, yields FinalVector[].
          {  // Rotate point3D from global to local system
            double InitialVector[3] = {xfld - xsrc, yfld - ysrc, zfld - zsrc};
            double TransformationMatrix[3][3] = {
                {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
            DirnCosn3D *DirCos = &(EleArr + elesrc - 1)->G.DC;
            TransformationMatrix[0][0] = DirCos->XUnit.X;
            TransformationMatrix[0][1] = DirCos->XUnit.Y;
            TransformationMatrix[0][2] = DirCos->XUnit.Z;
            TransformationMatrix[1][0] = DirCos->YUnit.X;
            TransformationMatrix[1][1] = DirCos->YUnit.Y;
            TransformationMatrix[1][2] = DirCos->YUnit.Z;
            TransformationMatrix[2][0] = DirCos->ZUnit.X;
            TransformationMatrix[2][1] = DirCos->ZUnit.Y;
            TransformationMatrix[2][2] = DirCos->ZUnit.Z;
            double FinalVector[3] = {0., 0., 0.};

            for (int i = 0; i < 3; ++i) {
              for (int j = 0; j < 3; ++j) {
                FinalVector[i] += TransformationMatrix[i][j] * InitialVector[j];
              }
            }

            localP.X = FinalVector[0];
            localP.Y = FinalVector[1];
            localP.Z = FinalVector[2];
          }

          // Initiate debugging, if necessary
          if ((elefld == 0) && (elesrc == 0))
            DebugISLES = 1;
          else
            DebugISLES = 0;

          Inf[elefld][elesrc] = ComputeInfluence(elefld, elesrc, &localP,
                                                 &(EleArr + elesrc - 1)->G.DC);
          if (DebugLevel == 301) {
            printf("elefld: %d, elesrc: %d, Influence: %.16lg\n", elefld,
                   elesrc, Inf[elefld][elesrc]);
          }

          // Take care of reflections	of basic device
          // At present, reflection on a single mirror is allowed

          if (MirrorTypeX[primsrc] || MirrorTypeY[primsrc] ||
              MirrorTypeZ[primsrc]) {
            printf(
                "Mirror not correctly implemented in this version of neBEM "
                "...\n");
            exit(0);

            Point3D fldpt, srcpt;
            DirnCosn3D DirCos;

            fldpt.X = xfld;
            fldpt.Y = yfld;
            fldpt.Z = zfld;
            srcpt.X = xsrc;
            srcpt.Y = ysrc;
            srcpt.Z = zsrc;

            if (MirrorTypeX[primsrc]) {
              MirrorTypeY[primsrc] = 0;
              MirrorTypeZ[primsrc] = 0;
            }
            if (MirrorTypeY[primsrc])  // no point checking MirrorTypeX
              MirrorTypeZ[primsrc] = 0;

            // If the reflection is other than that of an element (a known space
            // charge, for example), elesrc should be 0 and no question of
            // reflection of DC would arise. However, if the space charge is
            // itself distributed on a surface or in a volume, reflection of DC
            // etc will be important. What happens when reflection of an wire is
            // considered? At a later stage, reflection and periodicity should
            // become a property of an element and should be computed through a
            // single call of ComputeInfluence
            if (MirrorTypeX[primsrc]) {
              localP = ReflectOnMirror('X', elesrc, srcpt, fldpt,
                                       MirrorDistXFromOrigin[primsrc], &DirCos);
              double AddnalInfl =
                  ComputeInfluence(elefld, elesrc, &localP, &DirCos);

              if (MirrorTypeX[primsrc] ==
                  1)  // element having opposite charge density
                Inf[elefld][elesrc] -= AddnalInfl;  // classical image charge
              if (MirrorTypeX[primsrc] ==
                  2)  // element having same charge density
                Inf[elefld][elesrc] += AddnalInfl;
            }

            if (MirrorTypeY[primsrc]) {
              localP = ReflectOnMirror('Y', elesrc, srcpt, fldpt,
                                       MirrorDistYFromOrigin[primsrc], &DirCos);
              double AddnalInfl =
                  ComputeInfluence(elefld, elesrc, &localP, &DirCos);

              if (MirrorTypeY[primsrc] ==
                  1)  // element having opposite charge density
                Inf[elefld][elesrc] -= AddnalInfl;  // classical image charge
              if (MirrorTypeY[primsrc] ==
                  2)  // element having same charge density
                Inf[elefld][elesrc] += AddnalInfl;
            }

            if (MirrorTypeZ[primsrc]) {
              localP = ReflectOnMirror('Z', elesrc, srcpt, fldpt,
                                       MirrorDistZFromOrigin[primsrc], &DirCos);
              double AddnalInfl =
                  ComputeInfluence(elefld, elesrc, &localP, &DirCos);

              if (MirrorTypeZ[primsrc] ==
                  1)  // element having opposite charge density
                Inf[elefld][elesrc] -= AddnalInfl;  // classical image charge
              if (MirrorTypeZ[primsrc] ==
                  2)  // element having same charge density
                Inf[elefld][elesrc] += AddnalInfl;
            }

            if (DebugLevel == 301) {
              printf("After reflection of basic device =>\n");
              printf("elefld: %d, elesrc: %d, Influence: %.16lg\n", elefld,
                     elesrc, Inf[elefld][elesrc]);
            }
          }  // reflections of basic device, taken care of

          DebugISLES =
              0;  // Stop beyond basic device - may need changes. CHECK!
        }  // end of influence due to elements belonging to the basic device

        {  // Influence due to virtual elements

          // If the source element is repeated due to periodicity (the field
          // points remain unchanged), the influence at the field point will
          // change due to additional influence of the extra elements. The
          // repeated elements have properties identical to the real element
          // (including the solved value of charge density - otherwise we would
          // not be able to simply add up the influence at the field point!),
          // except that its location is determined by direction of periodicty
          // (at present assumed to be along one of the coordinate axes, but
          // this constraint can be easily relaxed later - we can have a
          // direction cosine along which the elements may be repeated) and the
          // distance of repeatition. Thus, for each repeated element, we need
          // to compute its position and the rest remains identical as above.
          // The influence, evaluated as a temporary double, can be added to the
          // base value computed above. PeriodicInX etc are either zero or +ve

          if ((PeriodicTypeX[primsrc] == 1) || (PeriodicTypeY[primsrc] == 1) ||
              (PeriodicTypeZ[primsrc] == 1)) {
            if (PeriodicInX[primsrc] || PeriodicInY[primsrc] ||
                PeriodicInZ[primsrc]) {
              for (int xrpt = -PeriodicInX[primsrc];
                   xrpt <= PeriodicInX[primsrc]; ++xrpt) {
                double XOfRpt = xsrc + XPeriod[primsrc] * (double)xrpt;

                for (int yrpt = -PeriodicInY[primsrc];
                     yrpt <= PeriodicInY[primsrc]; ++yrpt) {
                  double YOfRpt = ysrc + YPeriod[primsrc] * (double)yrpt;

                  for (int zrpt = -PeriodicInZ[primsrc];
                       zrpt <= PeriodicInZ[primsrc]; ++zrpt) {
                    double ZOfRpt = zsrc + ZPeriod[primsrc] * (double)zrpt;

                    if ((xrpt == 0) && (yrpt == 0) && (zrpt == 0))
                      continue;  // this is the base device

                    Point3D localP;
                    // axis direction in the global system
                    {  // Rotate point3D from global to local system
                      // Vector in the GCS
                      double InitialVector[3] = {xfld - XOfRpt, yfld - YOfRpt,
                                                 zfld - ZOfRpt};
                      double TransformationMatrix[3][3] = {
                          {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
                      DirnCosn3D *DirCos = &(EleArr + elesrc - 1)->G.DC;
                      TransformationMatrix[0][0] = DirCos->XUnit.X;
                      TransformationMatrix[0][1] = DirCos->XUnit.Y;
                      TransformationMatrix[0][2] = DirCos->XUnit.Z;
                      TransformationMatrix[1][0] = DirCos->YUnit.X;
                      TransformationMatrix[1][1] = DirCos->YUnit.Y;
                      TransformationMatrix[1][2] = DirCos->YUnit.Z;
                      TransformationMatrix[2][0] = DirCos->ZUnit.X;
                      TransformationMatrix[2][1] = DirCos->ZUnit.Y;
                      TransformationMatrix[2][2] = DirCos->ZUnit.Z;
                      double FinalVector[3] = {0., 0., 0.};

                      for (int i = 0; i < 3; ++i) {
                        for (int j = 0; j < 3; ++j) {
                          FinalVector[i] +=
                              TransformationMatrix[i][j] * InitialVector[j];
                        }
                      }

                      localP.X = FinalVector[0];  // Vector in the ECS
                      localP.Y = FinalVector[1];
                      localP.Z = FinalVector[2];
                    }  // Point3D rotated

                    // Direction cosines remain unchanged for a regular
                    // repetition
                    double AddnalInfl = ComputeInfluence(
                        elefld, elesrc, &localP, &(EleArr + elesrc - 1)->G.DC);
                    Inf[elefld][elesrc] += AddnalInfl;

                    if (DebugLevel == 301) {
                      printf("REPEATED\n");
                      printf("elefld: %d, elesrc: %d\n", elefld, elesrc);
                      printf("xsrc: %lg, ysrc: %lg, zsrc: %lg\n", xsrc, ysrc,
                             zsrc);
                      printf("xrpt: %d, yrpt: %d, zrpt: %d\n", xrpt, yrpt,
                             zrpt);
                      printf("XOfRpt: %lg, YOfRpt: %lg, ZOfRpt: %lg\n", XOfRpt,
                             YOfRpt, ZOfRpt);
                      printf("AddnalInfl: %lg\n", AddnalInfl);
                      printf("Inf: %lg\n", Inf[elefld][elesrc]);
                    }

                    // Take care of reflections of repetitions
                    // At present, reflection on a single mirror is allowed
                    if (MirrorTypeX[primsrc] || MirrorTypeY[primsrc] ||
                        MirrorTypeZ[primsrc]) {
                      printf(
                          "Mirror not correctly implemented in this version of "
                          "neBEM ...\n");
                      exit(0);

                      Point3D fldpt, srcpt;
                      DirnCosn3D DirCos;

                      fldpt.X = xfld;
                      fldpt.Y = yfld;
                      fldpt.Z = zfld;
                      srcpt.X = XOfRpt;
                      srcpt.Y = YOfRpt;
                      srcpt.Z = ZOfRpt;

                      if (MirrorTypeX[primsrc]) {
                        MirrorTypeY[primsrc] = 0;
                        MirrorTypeZ[primsrc] = 0;
                      }
                      if (MirrorTypeY[primsrc]) MirrorTypeZ[primsrc] = 0;

                      if (MirrorTypeX[primsrc]) {
                        localP = ReflectOnMirror('X', elesrc, srcpt, fldpt,
                                                 MirrorDistXFromOrigin[primsrc],
                                                 &DirCos);
                        AddnalInfl =
                            ComputeInfluence(elefld, elesrc, &localP, &DirCos);

                        if (MirrorTypeX[primsrc] ==
                            1)  // opposite charge density
                          Inf[elefld][elesrc] -=
                              AddnalInfl;  // classical image charge
                        if (MirrorTypeX[primsrc] == 2)  // same charge density
                          Inf[elefld][elesrc] += AddnalInfl;
                      }

                      if (MirrorTypeY[primsrc]) {
                        localP = ReflectOnMirror('Y', elesrc, srcpt, fldpt,
                                                 MirrorDistYFromOrigin[primsrc],
                                                 &DirCos);
                        AddnalInfl =
                            ComputeInfluence(elefld, elesrc, &localP, &DirCos);

                        if (MirrorTypeY[primsrc] ==
                            1)  // opposite charge density
                          Inf[elefld][elesrc] -=
                              AddnalInfl;  // classical image charge
                        if (MirrorTypeY[primsrc] == 2)  // same charge density
                          Inf[elefld][elesrc] += AddnalInfl;
                      }

                      if (MirrorTypeZ[primsrc]) {
                        localP = ReflectOnMirror('Z', elesrc, srcpt, fldpt,
                                                 MirrorDistZFromOrigin[primsrc],
                                                 &DirCos);
                        AddnalInfl =
                            ComputeInfluence(elefld, elesrc, &localP, &DirCos);

                        if (MirrorTypeZ[primsrc] ==
                            1)  // opposite charge density
                          Inf[elefld][elesrc] -=
                              AddnalInfl;  // classical image charge
                        if (MirrorTypeZ[primsrc] == 2)  // same charge density
                          Inf[elefld][elesrc] += AddnalInfl;
                      }

                      if (DebugLevel == 301) {
                        printf("REPEATED and reflected\n");
                        printf("elefld: %d, elesrc: %d\n", elefld, elesrc);
                        printf("xsrc: %lg, ysrc: %lg, zsrc: %lg\n", xsrc, ysrc,
                               zsrc);
                        printf("xrpt: %d, yrpt: %d, zrpt: %d\n", xrpt, yrpt,
                               zrpt);
                        printf("XOfRpt: %lg, YOfRpt: %lg, ZOfRpt: %lg\n",
                               XOfRpt, YOfRpt, ZOfRpt);
                        printf("AddnalInfl: %lg\n", AddnalInfl);
                        printf("Inf: %lg\n", Inf[elefld][elesrc]);
                      }
                    }  // reflections of repetitions taken care of

                  }  // for zrpt
                }    // for yrpt
              }      // for xrpt
            }        // PeriodicInX || PeriodicInY || PeriodicInZ
          }          // PeriodicType == 1
        }            // end of influence due to virtual elements

      }  // loop for elesrc, source element (influencing)

      // printf("\b\b\b\b\b\b");
    }  // loop for elefld, field element (influenced)
  }    // pragma omp parallel

  // Enforce total charge on the system to be zero.
  // All the voltages in the system need to be shifted by an unknown amount
  // V_shift
  // V_shift is obtained by introducing one additional row and column and
  // impposing the constraint that the sum of all the charges in the system is
  // zero.
  // Please note that charge = Element Charge Density * Element Area
  if (OptSystemChargeZero) {
    // an additional column
    for (int row = 1; row <= NbEqns; ++row) {
      if (((EleArr + row - 1)->E.Type == 1) ||
          ((EleArr + row - 1)->E.Type == 3))  // The
        Inf[row][NbUnknowns] = 1.0;  // VSystemChargeZero is subtracted only
      else                           // from potentials
        Inf[row][NbUnknowns] = 0.0;
    }

    // an additional row
    for (int col = 1; col <= NbUnknowns; ++col)
      Inf[NbEqns][col] =
          (EleArr + col - 1)->G.dA;  // if charge density is computed
    // Inf[NbEqns][col] = 1.0;	// if charge is computed

    // the last element
    Inf[NbEqns][NbUnknowns] = 0.0;
  }  // if(OptSystemChargeZero)
  else {
    VSystemChargeZero = 0.0;
  }

  if (NbFloatingConductors)  // assume only one floating conductor
  {
    // an additional column
    for (int row = 1; row <= NbEqns; ++row) {
      int etfld = (EleArr + row - 1)->E.Type;
      if (etfld == 3)  // element of a floating conductor
        Inf[row][NbUnknowns] = -1.0;
      else
        Inf[row][NbUnknowns] = 0.0;
    }  // additional column

    // an additional row
    for (int col = 1; col <= NbUnknowns; ++col) {
      int etfld = (EleArr + col - 1)->E.Type;
      if (etfld == 3)  // element of a floating conductor
      {
        Inf[NbEqns][col] =
            (EleArr + col - 1)->G.dA;  // if charge density is computed
        // Inf[NbEqns][col] = 1.0;	// if charge is computed
      } else {
        Inf[NbEqns][col] = 0.0;  // if charge density is computed
        // Inf[NbEqns][col] = 0.0;	// if charge is computed
      }
    }  // additional row

    // the last element
    Inf[NbEqns][NbUnknowns] = 0.0;
  }  // if NbFloatingConductors

  if (OptStoreInflMatrix && OptFormattedFile) {
    printf("storing the influence matrix in a formatted file ...\n");
    fflush(stdout);

    char InflFile[256];
    strcpy(InflFile, MeshOutDir);
    strcat(InflFile, "/Infl.out");
    FILE *fInf = fopen(InflFile, "w+");
    if (fInf == NULL) {
      neBEMMessage("LHMatrix - InflFile");
      return -1;
    }
    fprintf(fInf, "%d %d\n", NbEqns, NbUnknowns);

    for (int elefld = 1; elefld <= NbEqns; ++elefld) {
      for (int elesrc = 1; elesrc <= NbUnknowns; ++elesrc)
        fprintf(fInf, "%.16lg\n", Inf[elefld][elesrc]);
      fprintf(fInf, "\n");
    }

    fclose(fInf);
  }  // if OptStoreInflMatrix && OptFormattedFile

  if (OptStoreInflMatrix &&
      OptUnformattedFile)  // Raw cannot be implemented now.
  {  // It may be because of the memory allocation using the NR routines -
    neBEMMessage(
        "LHMatrix - Binary write of Infl matrix not implemented yet.\n");
    return -1;

    char InflFile[256];
    strcpy(InflFile, MeshOutDir);
    strcat(InflFile, "/RawInfl.out");
    FILE *fInf = fopen(InflFile, "wb");
    if (fInf == NULL) {
      neBEMMessage("LHMatrix - RawInflFile");
      return -1;
    }
    printf("\nfInf: %p\n", (void *)fInf);
    int rfw = fwrite(Inf, sizeof(double), NbEqns * NbUnknowns, fInf);
    fclose(fInf);
    printf("\nNb of items successfully written in raw mode in %s is %d\n",
           InflFile, rfw);

    /* following block used to check raw reading from a file written in raw
    double **RawInf;
    RawInf = dmatrix(1,NbEqns,1,NbUnknowns);
    strcpy(InflFile, MeshOutDir); strcat(InflFile, "/RawInfl.out");
    fInf = fopen(InflFile, "rb");
    // assert(fInf != NULL);
          if(fInf == NULL)
                  {
                  neBEMMessage("LHMatrix - RawInflFile");
                  return -1;
                  }
    rfw = fread(RawInf, sizeof(double), NbEqns*NbUnknowns, fInf);
    fclose(fInf);
    printf("Nb of items successfully read in raw mode from %s is %d\n",
                                  InflFile, rfw);
    for(int unknown = 1; unknown <= NbUnknowns; ++unknown)
          for(int eqn = 1; eqn <= NbEqns; ++eqn)
                  printf("Unknown:%d, Eqn:%d => diff Inf: %lg, RawInf: %lg is
    %lg\n", unknown, eqn, Inf[unknown][eqn], RawInf[unknown][eqn],
                                                  fabs(Inf[unknown][eqn] -
    RawInf[unknown][eqn])); Please do not delete */
  }  // if OptStoreInflMatrix && OptUnformattedFile

  neBEMState = 6;

  return (0);
}  // end of LHMatrix

/* To be tried later
// Generate a row vector of influence (in terms of potential or continutity) due
to all the elements at a given field point
// This row vector, when multiplied by the charge density on all the elements
(how about non-elemental charges?) gives rise to the
// total potential, or field (normal component for a given coordinate system) at
that point
// If no specific direction is provided, the GCS is used as default
Vector1D* InflVec(Point3D fldpt, DirnCosn3D *fldDC, int Pot0Cont1) {
  int dbgFn = 0;
  int primsrc;
  double xsrc, ysrc, zsrc;
  Point3D localP;

  printf("\nLHMatrix: The size of the Influence coefficient vector is %d\n",
         NbUnknowns); fflush(stdout);

  // Influence coefficient vector
  // For each field point influences from all the source elements
  // need to be summed up.
  // The source elements are followed using elesrc (source counter)
  printf("Computing influence coefficient vector ... will take some time
...\n");

  int nthreads = 1, tid = 0;
#pragma omp parallel private(nthreads, tid)
  {
    if(dbgFn) {
      tid = omp_get_thread_num();
      if (tid == 0) {
        nthreads = omp_get_num_threads();
        printf("Starting influence matrix computation with %d threads\n",
               nthreads);
      }
    }

    double xfld = fldpt.X;
    double yfld = fldpt.Y;
    double zfld = fldpt.Z;

#pragma omp for private(primsrc, xsrc, ysrc, zsrc, localP)
    for (int elesrc = 1; elesrc <= NbElements; ++elesrc) {
      if (DebugLevel == 301) {
        printf("\n\nelesrc: %d\n", elesrc);
      }

    // Retrieve element properties at the field point
    int primsrc = (EleArr+elesrc-1)->PrimitiveNb;
    double xsrc = (EleArr+elesrc-1)->G.Origin.X;
    double ysrc = (EleArr+elesrc-1)->G.Origin.Y;
    double zsrc = (EleArr+elesrc-1)->G.Origin.Z;

                // The total influence is due to elements on the basic device
and due to
                // virtual elements arising out of repetition, reflection etc
and not
                // residing on the basic device

                // Influence due to elements belonging to the basic device
                {
                // point translated to the ECS origin, but axes direction global
        { // Rotate point3D from global to local system
        double InitialVector[4];
        double TransformationMatrix[4][4] = {{0.0, 0.0, 0.0, 0.0},
                                                {0.0, 0.0, 0.0, 0.0},
                                                {0.0, 0.0, 0.0, 0.0},
                                                {0.0, 0.0, 0.0, 1.0}};
        DirnCosn3D *DirCos = &(EleArr+elesrc-1)->G.DC;
        double FinalVector[4];

        InitialVector[0] = xfld - xsrc; InitialVector[1] = yfld - ysrc;
        InitialVector[2] = zfld - zsrc; InitialVector[3] = 1.0;

        TransformationMatrix[0][0] = DirCos->XUnit.X;
        TransformationMatrix[0][1] = DirCos->XUnit.Y;
        TransformationMatrix[0][2] = DirCos->XUnit.Z;
        TransformationMatrix[1][0] = DirCos->YUnit.X;
        TransformationMatrix[1][1] = DirCos->YUnit.Y;
        TransformationMatrix[1][2] = DirCos->YUnit.Z;
        TransformationMatrix[2][0] = DirCos->ZUnit.X;
        TransformationMatrix[2][1] = DirCos->ZUnit.Y;
        TransformationMatrix[2][2] = DirCos->ZUnit.Z;

        for(int i = 0; i < 4; ++i)
        {
        FinalVector[i] = 0.0;
        for(int j = 0 ; j < 4; ++j)
        {
        FinalVector[i] += TransformationMatrix[i][j] * InitialVector[j];
        }
        }

        localP.X = FinalVector[0];
        localP.Y = FinalVector[1];
        localP.Z = FinalVector[2];
        } // Point3D rotated

                // Initiate debugging, if necessary
                if(elesrc == 0)
      DebugISLES = 1;
                else
                        DebugISLES = 0;

                InfVec[elesrc] = ComputeEleInf(Pot0Cont1, elesrc, &localP,
&(EleArr+elesrc-1)->G.DC); if(DebugLevel == 301)
                        {
                        printf("elesrc: %d, Influence: %.16lg\n", elesrc,
InfVec[elesrc]);
                        }

                // Take care of reflections	of basic device
                // At present, reflection on a single mirror is allowed

          if(MirrorTypeX[primsrc] || MirrorTypeY[primsrc] ||
MirrorTypeZ[primsrc])
                {
    printf("Mirror not correctly implemented in this version of neBEM ...\n");
    exit(0);

                Point3D srcpt;
                DirnCosn3D DirCos;

                srcpt.X = xsrc; srcpt.Y = ysrc; srcpt.Z = zsrc;

                if(MirrorTypeX[primsrc])
                        { MirrorTypeY[primsrc] = 0; MirrorTypeZ[primsrc] = 0; }
                if(MirrorTypeY[primsrc])	// no point checking MirrorTypeX
                        MirrorTypeZ[primsrc] = 0;

                // If the reflection is other than that of an element (a known
space charge,
                // for example), elesrc should be 0 and no question of
reflection of DC
                // would arise. However, if the space charge is itself
distributed on a
                // surface or in a volume, reflection of DC etc will be
important. What
                // happens when reflection of an wire is considered?
                // At a later stage, reflection and periodicity should become a
property
                // of an element and should be computed through a single call of
                // ComputeInfluence
                if(MirrorTypeX[primsrc])
                        {
                        localP = ReflectOnMirror('X', elesrc, srcpt, fldpt,
                                                                                                                                MirrorDistXFromOrigin[primsrc], &DirCos);
                        double AddnalInfl = ComputeEleInf(Pot0Cont1, elesrc,
&localP, &DirCos);

                        if(MirrorTypeX[primsrc] == 1)	// element having
opposite charge density Inf[elefld][elesrc] -= AddnalInfl;	// classical
image charge if(MirrorTypeX[primsrc] == 2)	// element having same charge
density Inf[elefld][elesrc] += AddnalInfl;
                        }

                if(MirrorTypeY[primsrc])
                        {
                        localP = ReflectOnMirror('Y', elesrc, srcpt, fldpt,
                                                                                                                                MirrorDistYFromOrigin[primsrc], &DirCos);
                        double AddnalInfl = ComputeEleInf(Pot0Cont1, elesrc,
&localP, &DirCos);

                        if(MirrorTypeY[primsrc] == 1)	// element having
opposite charge density Inf[elefld][elesrc] -= AddnalInfl;	// classical
image charge if(MirrorTypeY[primsrc] == 2)	// element having same charge
density Inf[elefld][elesrc] += AddnalInfl;
                        }

                if(MirrorTypeZ[primsrc])
                        {
                        localP = ReflectOnMirror('Z', elesrc, srcpt, fldpt,
                                                                                                                                MirrorDistZFromOrigin[primsrc], &DirCos);
                        double AddnalInfl = ComputeEleInf(Pot0Cont1, elesrc,
&localP, &DirCos);

                        if(MirrorTypeZ[primsrc] == 1)	// element having
opposite charge density InfVec[elesrc] -= AddnalInfl;	// classical image
charge if(MirrorTypeZ[primsrc] == 2)	// element having same charge density
                                InfVec[elesrc] += AddnalInfl;
                        }

                if(DebugLevel == 301)
                        {
                        printf("After reflection of basic device =>\n");
                        printf("elesrc: %d, Influence: %.16lg\n", elesrc,
InfVec[elesrc]);
                        }
                }	// reflections of basic device, taken care of

                DebugISLES = 0;	// Stop beyond basic device - may need changes.
CHECK! }	// end of influence due to elements belonging to the basic
device

                {	// Influence due to virtual elements

                // If the source element is repeated due to periodicity (the
field points
                // remain unchanged), the influence at the field point will
change due to
                // additional influence of the extra elements. The repeated
elements
                // have properties identical to the real element (including the
solved
                // value of charge density - otherwise we would not be able to
simply add
                // up the influence at the field point!), except that its
location is
                // determined by direction of periodicty (at present assumed to
be along
                // one of the coordinate axes, but this constraint can be easily
relaxed
                // later - we can have a direction cosine along which the
elements may be
                // repeated) and the distance of repeatition. Thus, for each
repeated
                // element, we need to compute its position and the rest remains
identical
                // as above. The influence, evaluated as a temporary double, can
be added
                // to the base value computed above.
                // PeriodicInX etc are either zero or +ve

                if((PeriodicTypeX[primsrc] == 1) || (PeriodicTypeY[primsrc] ==
1)
                                || (PeriodicTypeZ[primsrc] == 1))
                        {
                        if(PeriodicInX[primsrc] || PeriodicInY[primsrc] ||
PeriodicInZ[primsrc])
                                {
                                double AddnalInfl=0.0,
                                                                XOfRpt, YOfRpt,
ZOfRpt;

                                for(int xrpt = -PeriodicInX[primsrc];
                                                        xrpt <=
PeriodicInX[primsrc]; ++xrpt)
                                        {
                                        XOfRpt = xsrc + XPeriod[primsrc] *
(double)xrpt;

                                        for(int yrpt = -PeriodicInY[primsrc];
                                                                yrpt <=
PeriodicInY[primsrc]; ++yrpt)
                                                {
                                                YOfRpt = ysrc + YPeriod[primsrc]
* (double)yrpt;

                                                for(int zrpt =
-PeriodicInZ[primsrc]; zrpt <= PeriodicInZ[primsrc]; ++zrpt)
                                                        {
                                                        ZOfRpt = zsrc +
ZPeriod[primsrc] * (double)zrpt;

                                                        if( (xrpt == 0) && (yrpt
== 0) && (zrpt == 0) ) continue;	// this is the base device

                                                        // axis direction in the
global system { // Rotate point3D from global to local system double
InitialVector[4]; double TransformationMatrix[4][4] = {{0.0, 0.0, 0.0, 0.0},
                                                                                {0.0, 0.0, 0.0, 0.0},
                                                                                {0.0, 0.0, 0.0, 0.0},
                                                                                {0.0, 0.0, 0.0, 1.0}};
                                        DirnCosn3D *DirCos =
&(EleArr+elesrc-1)->G.DC; double FinalVector[4];

                                        InitialVector[0] = xfld - XOfRpt;
// Vector in the GCS InitialVector[1] = yfld - YOfRpt; InitialVector[2] = zfld -
ZOfRpt; InitialVector[3] = 1.0;

                                        TransformationMatrix[0][0] =
DirCos->XUnit.X; TransformationMatrix[0][1] = DirCos->XUnit.Y;
                                        TransformationMatrix[0][2] =
DirCos->XUnit.Z; TransformationMatrix[1][0] = DirCos->YUnit.X;
                                        TransformationMatrix[1][1] =
DirCos->YUnit.Y; TransformationMatrix[1][2] = DirCos->YUnit.Z;
                                        TransformationMatrix[2][0] =
DirCos->ZUnit.X; TransformationMatrix[2][1] = DirCos->ZUnit.Y;
                                        TransformationMatrix[2][2] =
DirCos->ZUnit.Z;

                                        for(int i = 0; i < 4; ++i)
                                        {
                                        FinalVector[i] = 0.0;
                                        for(int j = 0 ; j < 4; ++j)
                                                {
                                                FinalVector[i] +=
TransformationMatrix[i][j]
                                                                                                                                                * InitialVector[j];
                                                }
                                        }

                                        localP.X = FinalVector[0];
// Vector in the ECS localP.Y = FinalVector[1]; localP.Z = FinalVector[2]; } //
Point3D rotated

                                                        // Direction cosines
remain unchanged for a regular repetition AddnalInfl = ComputeEleInf(Pot0Cont1,
elesrc, &localP, &(EleArr+elesrc-1)->G.DC); InfVec[elesrc] += AddnalInfl;

                                                        if(DebugLevel == 301)
                                                                {
                                                                printf("REPEATED\n");
                                                                printf("elesrc:
%d\n", elefld, elesrc); printf("xsrc: %lg, ysrc: %lg, zsrc: %lg\n", xsrc, ysrc,
zsrc); printf("xrpt: %d, yrpt: %d, zrpt: %d\n", xrpt, yrpt, zrpt);
                                                                printf("XOfRpt:
%lg, YOfRpt: %lg, ZOfRpt: %lg\n", XOfRpt, YOfRpt, ZOfRpt); printf("AddnalInfl:
%lg\n", AddnalInfl); printf("InfVec: %lg\n", InfVec[elesrc]);
                                                                }

                                                        // Take care of
reflections of repetitions
                                                        // At present,
reflection on a single mirror is allowed if(MirrorTypeX[primsrc] ||
MirrorTypeY[primsrc]
                                                                        ||
MirrorTypeZ[primsrc])
                                                        {
                                        printf("Mirror not correctly implemented
in this version of neBEM ...\n"); exit(0);

                                                        Point3D srcpt;
                                                        DirnCosn3D DirCos;

                                                        srcpt.X = XOfRpt;
srcpt.Y = YOfRpt; srcpt.Z = ZOfRpt;

                                                        if(MirrorTypeX[primsrc])
                                                                {
MirrorTypeY[primsrc] = 0; MirrorTypeZ[primsrc] = 0; } if(MirrorTypeY[primsrc])
MirrorTypeZ [primsrc]= 0;

                                                        if(MirrorTypeX[primsrc])
                                                                {
                                                                localP =
ReflectOnMirror('X', elesrc, srcpt, fldpt, MirrorDistXFromOrigin[primsrc],
                                                                                                                                                                        &DirCos);
                                                                AddnalInfl =
ComputeEleInf(Pot0Cont1, elesrc, &localP, &DirCos);

                                                                if(MirrorTypeX[primsrc]
== 1)	// opposite charge density InfVec[elesrc] -= AddnalInfl;	//
classical image charge if(MirrorTypeX[primsrc] == 2)	// same charge density
                                                                        InfVec[elesrc]
+= AddnalInfl;
                                                                }

                                                        if(MirrorTypeY[primsrc])
                                                                {
                                                                localP =
ReflectOnMirror('Y', elesrc, srcpt, fldpt, MirrorDistYFromOrigin[primsrc],
                                                                                                                                                                        &DirCos);
                                                                AddnalInfl =
ComputeEleInf(Pot0Cont1, elesrc, &localP, &DirCos);

                                                                if(MirrorTypeY[primsrc]
== 1)	// opposite charge density InfVec[elesrc] -= AddnalInfl;	//
classical image charge if(MirrorTypeY[primsrc] == 2)	// same charge density
                                                                        InfVec[elesrc]
+= AddnalInfl;
                                                                }

                                                        if(MirrorTypeZ[primsrc])
                                                                {
                                                                localP =
ReflectOnMirror('Z', elesrc, srcpt, fldpt, MirrorDistZFromOrigin[primsrc],
                                                                                                                                                                        &DirCos);
                                                                AddnalInfl =
ComputeEleInf(Pot0Cont1, elesrc, &localP, &DirCos);

                                                                if(MirrorTypeZ[primsrc]
== 1)	// opposite charge density InfVec[elesrc] -= AddnalInfl;	//
classical image charge if(MirrorTypeZ[primsrc] == 2)	// same charge density
                                                                        InfVec[elesrc]
+= AddnalInfl;
                                                                }

                                                        if(DebugLevel == 301)
                                                                {
                                                                printf("REPEATED
and reflected\n"); printf("elesrc: %d\n", elesrc); printf("xsrc: %lg, ysrc: %lg,
zsrc: %lg\n", xsrc, ysrc, zsrc); printf("xrpt: %d, yrpt: %d, zrpt: %d\n", xrpt,
yrpt, zrpt); printf("XOfRpt: %lg, YOfRpt: %lg, ZOfRpt: %lg\n", XOfRpt, YOfRpt,
ZOfRpt); printf("AddnalInfl: %lg\n", AddnalInfl); printf("InfVec: %lg\n",
InfVec[elesrc]);
                                                                }
                                                        }	// reflections
of repetitions taken care of

                                                        }	// for zrpt
                                                }	// for yrpt
                                        }	// for xrpt
                                }	// PeriodicInX || PeriodicInY ||
PeriodicInZ }	// PeriodicType == 1 }	// end of influence due to virtual
elements

                }	// loop for elesrc, source element (influencing)

        // printf("\b\b\b\b\b\b");
}	// pragma omp parallel

if(OptStoreInfVec && OptFormattedFile)
        {
        printf("storing the influence matrix in a formatted file ...\n");
        fflush(stdout);

        char InfVecFile[256];
        strcpy(InfVecFile, MeshOutDir); strcat(InfVecFile, "/InfVec.out");
        FILE *fInfVec = fopen(InfVecFile, "w+");
        if(fInfVec == NULL)
                {
                neBEMMessage("InfluenceVector - InfVecFile");
                return -1;
                }
        fprintf(fInfVec, "%d\n", NbUnknowns);

        for(int elesrc = 1; elesrc <= NbUnknowns; ++elesrc)
                fprintf(fInfVec, "%.16lg\n", InfVec[elesrc]);
        fprintf(fInfVec, "\n");

        fclose(fInfVec);
        }	// if OptStoreInflMatrix && OptFormattedFile

if(OptStoreInfVec && OptUnformattedFile)	// Raw cannot be implemented
now.
  {		// It may be because of the memory allocation using the NR
routines - neBEMMessage("LHMatrix - Binary write of Infl matrix not implemented
yet.\n"); return -1;

  char InfVecFile[256];
  strcpy(InfVecFile, MeshOutDir); strcat(InfVecFile, "/RawInfVec.out");
  FILE *fInfVec = fopen(InfVecFile, "wb");
        if(fInfVec == NULL)
                {
                neBEMMessage("InfluenceVector - RawInfVecFile");
                return -1;
                }
  printf("\nfInfVec: %p\n", fInfVec);
  int rfw = fwrite(InfVec, sizeof(double), NbUnknowns, fInfVec);
  fclose(fInfVec);
  printf("\nNb of items successfully written in raw mode in %s is %d\n",
          InfVecFile, rfw);

  SLASH-STAR following block used to check raw reading from a file written in
raw double **RawInf; RawInf = dmatrix(1,NbEqns,1,NbUnknowns); strcpy(InflFile,
MeshOutDir); strcat(InflFile, "/RawInfl.out"); fInf = fopen(InflFile, "rb");
  // assert(fInf != NULL);
        if(fInf == NULL)
                {
                neBEMMessage("LHMatrix - RawInflFile");
                return -1;
                }
  rfw = fread(RawInf, sizeof(double), NbEqns*NbUnknowns, fInf);
  fclose(fInf);
  printf("Nb of items successfully read in raw mode from %s is %d\n",
                                InflFile, rfw);
  for(int unknown = 1; unknown <= NbUnknowns; ++unknown)
        for(int eqn = 1; eqn <= NbEqns; ++eqn)
                printf("Unknown:%d, Eqn:%d => diff Inf: %lg, RawInf: %lg is
%lg\n", unknown, eqn, Inf[unknown][eqn], RawInf[unknown][eqn],
                                                fabs(Inf[unknown][eqn] -
RawInf[unknown][eqn])); Please do not delete STAR-SLASH }	// if
OptStoreInflMatrix && OptUnformattedFile

neBEMState = 6;

return(0);
}	// end of Influence
To be tried later */

// Invert a matrix of size NEquations X NbUnknowns
// Note that slightly modified influence coefficient matrix can be solved
// without going for a fresh re-inversion
int InvertMatrix(void) {
  int DecomposeMatrixSVD(double **SVDInf, double *SVDw, double **SVDv);

  InvMat = dmatrix(1, NbUnknowns, 1, NbEqns);

  if (OptGSL) {
    printf("InvertMatrix: matrix decomposition using GSL ... ");
    printf("no OpenMP implementation ...");
    fflush(stdout);

    int s;  // signum for LU decomposition

    gsl_matrix *m = gsl_matrix_alloc(NbUnknowns, NbEqns);
    gsl_matrix *inverse = gsl_matrix_alloc(NbUnknowns, NbEqns);
    gsl_permutation *perm = gsl_permutation_alloc(NbUnknowns);

    for (int i = 0; i < NbUnknowns; ++i)
      for (int j = 0; j < NbEqns; ++j)
        gsl_matrix_set(m, i, j, Inf[i + 1][j + 1]);

    // Make LU decomposition of matrix m
    gsl_linalg_LU_decomp(m, perm, &s);

    // Invert the matrix m
    gsl_linalg_LU_invert(m, perm, inverse);

    for (int i = 0; i < NbUnknowns; ++i)
      for (int j = 0; j < NbEqns; ++j)
        InvMat[i + 1][j + 1] = gsl_matrix_get(inverse, i, j);

    gsl_matrix_free(m);
    gsl_matrix_free(inverse);
    printf("InvertMatrix: ... completed using GSL\n");
  }  // if OptGSL

  if (OptSVD) {
    printf("InvertMatrix: matrix decomposition using SVD ... ");
    printf("no OpenMP implementation ...");
    fflush(stdout);

    clock_t SVDstartClock = clock();
    printf("ComputeSolution: Decomposing influence matrix ...\n");
    fflush(stdout);

    // These may as well be declared in neBEM.h because a signicant amount
    // of information can be extracted from them in a later version.
    double **SVDInf, *SVDw, **SVDv;
    SVDInf = dmatrix(1, NbEqns, 1, NbUnknowns);
    SVDw = dvector(1, NbUnknowns);
    SVDv = dmatrix(1, NbUnknowns, 1, NbUnknowns);  // a square matrix!

    // Keep the original Inf[] safe
    for (int i = 1; i <= NbEqns; i++) {
      int j;
#ifdef _OPENMP
#pragma omp parallel for private(j)
#endif
      for (j = 1; j <= NbUnknowns; j++)
        SVDInf[i][j] = Inf[i][j];  // end of omp parallel for
    }

    int fstatus = DecomposeMatrixSVD(SVDInf, SVDw, SVDv);
    if (fstatus != 0) {
      neBEMMessage("ComputeSolution - DecomposeMatrixSVD");
      return -1;
    }
    printf("ComputeSolution: Matrix decomposition over.\n");
    clock_t SVDstopClock = clock();
    neBEMTimeElapsed(SVDstartClock, SVDstopClock);
    printf("to singular value decompose the influence matrix.\n");

    // Find the pseudo-inverse of the influence coefficient matrix
    double **tmpmat;
    tmpmat = dmatrix(1, NbEqns, 1, NbUnknowns);

    // calculate w+ (transpose)u
    for (int j = 1; j <= NbUnknowns;
         j++) {  // w+ is obtained by replacing every non-zero diagonal entry of
                 // [W]
      int i;
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
      for (i = 1; i <= NbEqns; i++)  // (note W, not w) by its reciprocal and
      {                              // transposing the resulting matrix
        if (SVDw[j])                 // nonzero result only if w[i] is nonzero
        {
          tmpmat[j][i] = SVDInf[i][j] / SVDw[j];
        } else {
          tmpmat[j][i] = 0.0;
        }
      }  // end of omp parallel for
    }

    // multiply by [V] to get answer
    for (int i = 1; i <= NbUnknowns; i++)
      for (int j = 1; j <= NbEqns; j++) {
        InvMat[i][j] = 0.0;

        int k;
        double sum = 0.0;
#ifdef _OPENMP
#pragma omp parallel for private(k) reduction(+ : sum)
#endif
        for (k = 1; k <= NbUnknowns; k++)
          sum += SVDv[i][k] * tmpmat[k][j];  // end of omp parallel for

        InvMat[i][j] = sum;
      }

    free_dmatrix(tmpmat, 1, NbEqns, 1, NbUnknowns);
    // free SVD matrices? may need to be manitained for later use
    free_dmatrix(SVDInf, 1, NbEqns, 1, NbUnknowns);
    free_dvector(SVDw, 1, NbUnknowns);
    free_dmatrix(SVDv, 1, NbUnknowns, 1, NbUnknowns);
    printf("InvertMatrix: completed using SVD ...\n");
    fflush(stdout);
  }  // if OptSVD

  if (OptLU) {
    int *index;
    double d, *col, **y;
    y = dmatrix(1, NbUnknowns, 1, NbUnknowns);
    col = dvector(1, NbUnknowns);
    index = ivector(1, NbUnknowns);

    double **tmpInf;  // temporary influence matrix necessary for ludcmp
    tmpInf = dmatrix(1, NbEqns, 1, NbUnknowns);
    for (int i = 1; i <= NbEqns; i++)  // Keep original Inf[] safe
    {
      int j;
#ifdef _OPENMP
#pragma omp parallel for private(j)
#endif
      for (j = 1; j <= NbUnknowns; j++)
        tmpInf[i][j] = Inf[i][j];  // end of omp parallel for
    }

    printf("InvertMatrix: matrix decomposition using LU ... ");
    fflush(stdout);
    ludcmp(tmpInf, NbUnknowns, index, &d);  // The tmpInf matrix over-written

    for (int j = 1; j <= NbUnknowns; j++)  // Find inverse by columns.
    {
      int i;
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
      for (i = 1; i <= NbUnknowns; i++)
        col[i] = 0.0;  // end of omp parallel for

      col[j] = 1.0;

      lubksb(tmpInf, NbUnknowns, index, col);  // changed avatar of tmpInf used
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
      for (i = 1; i <= NbEqns; i++) {
        y[i][j] = col[i];
        InvMat[i][j] = y[i][j];
      }  // end of omp parallel for
      // printf("\b\b\b\b\b\b");
    }

    free_ivector(index, 1, NbUnknowns);
    free_dvector(col, 1, NbUnknowns);
    free_dmatrix(y, 1, NbUnknowns, 1, NbUnknowns);
    free_dmatrix(tmpInf, 1, NbEqns, 1, NbUnknowns);

    printf("InvertMatrix: completed using LU ...\n");
    fflush(stdout);
  }  // if OptLU

  // Since this is occurring within the InvertMatrix function, it is obvious
  // that the InfluenceMatrixFlag is true and the LHMatrix has already been
  // computed giving rise to a valid Inf[].
  // Free Inf[] if validation is not opted for, despite going through a fresh
  // solution attempt - highly deprecated!
  if (!OptValidateSolution) free_dmatrix(Inf, 1, NbEqns, 1, NbUnknowns);

  // It is necessary to write this file always if we want to avoid
  // matrix inversion for analyzing the same device with a different BC
  printf("OptStoreInvMatrix: %d, OptFormattedFile: %d\n", OptStoreInvMatrix,
         OptFormattedFile);

  if (OptStoreInvMatrix && OptFormattedFile) {
    printf("storing the inverted matrix in a formatted file ...\n");
    fflush(stdout);

    FILE *fInv;  // can be a very large file - change to raw and zipped format
    char InvMFile[256];
    strcpy(InvMFile, MeshOutDir);
    strcat(InvMFile, "/InvMat.out");
    fInv = fopen(InvMFile, "w");

    // following line may be removed after the dimension issue is resolved.
    fprintf(fInv, "%d %d\n", NbEqns, NbUnknowns);
    for (int i = 1; i <= NbEqns; i++)
      for (int j = 1; j <= NbUnknowns; j++)
        fprintf(fInv, "%.16le\n", InvMat[i][j]);

    fclose(fInv);
  }  // if OptStoreInvMatrix && OptFormattedFile

  if (OptStoreInvMatrix && OptUnformattedFile) {
    // not implemented yet
    // not implemented
    neBEMMessage("InvertMatrix - Binary write not yet implemented.");
    return (-1);
  }

  neBEMState = 7;

  return (0);
}  // end of InvertMatrix

// Decompose a matrix of size NbEqns X NbUnknowns using SVD
// NbEqns has to be equal to or greater than NbUnknowns for SVD to work.
int DecomposeMatrixSVD(double **SVDInf, double *SVDw, double **SVDv) {
  // The following needs optimization - wmin, wmax etc.
  double wmin, wmax;

  printf("DecomposeMatrix: matrix decomposition using SVD ... ");
  fflush(stdout);
  svdcmp(SVDInf, NbEqns, NbUnknowns, SVDw, SVDv);  // SVDInf matrix over-written
  printf("DecomposeMatrix: decomposition completed ...\n");
  fflush(stdout);

  wmax = SVDw[1];  // Will be the maximum singular value obtained - changed 0.0
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    int j;
#ifdef _OPENMP
#pragma omp for private(j)
#endif
    for (j = 1; j <= NbUnknowns; j++) {
      if (SVDw[j] > wmax) wmax = SVDw[j];
    }
  }  // omp parallel

  // This is where we set the threshold for singular values allowed to be
  // nonzero. The constant is typical, but not universal. You have to experiment
  // with your own application.
  wmin = wmax * 1.0e-12;
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    int j;
#ifdef _OPENMP
#pragma omp for private(j)
#endif
    for (j = 1; j <= NbUnknowns; j++) {
      if (SVDw[j] < wmin) SVDw[j] = 0.0;
    }
  }  // omp parallel
  // printf("wmin: %le, wmax: %le\n", wmin, wmax);

  neBEMState = 7;

  return (0);
}  // end of DecomposeMatrixSVD

// Read an inverted matrix of size NbEquations X NbUnknowns
// Called only when OptReadInvMatrix is true. Hence, no need to check the
// same all over again.
int ReadInvertedMatrix(void) {
  // Computation of solution using LU only in this version
  if (OptFormattedFile || OptUnformattedFile) {
    // InvMat = dmatrix(1,NbEqns,1,NbUnknowns); // may be reinstated after
    // the dimension issue is resolvedint ReadInvertedMatrix(void)
  } else {
    printf(
        "ReadInvertedMatrix: OptFormattedFile and OptUnformattedFile, both are "
        "false ... ");
    printf(
        "                    Can not read inverted matrix ... returning ...\n");
    return (-1);
  }

  // We need to provide two inputs to implement this
  // 	1. The fact that we are only changing the BC; MotherInputFile and object
  //     files should be the same as the original - in such an event, it is
  //     obvious where the following file is to be found.
  //  2. The new set of BCs and the new output file - can be a subdirectory
  //     of the original (may be default-ed to BC0).
  //  3. BC0, if existing, should never be overwritten / removed by the code.
  //     If at all, it should be removed manually.

  if (OptFormattedFile) {
    FILE *fInv;  // can be a very large file - change to raw and zipped format

    char InvMFile[256];
    strcpy(InvMFile, MeshOutDir);
    strcat(InvMFile, "/InvMat.out");
    fInv = fopen(InvMFile, "r");
    // assert(fInv != NULL);
    if (fInv == NULL) {
      neBEMMessage("ReadInvertedMatrix - inverted matrix not found.");
      return (-1);
    }

    // Retrieve the dimensions and create matrix
    // These lines may be removed once the dimension issue is resolved
    // Better still, these numbers can be used to cross-check that there is
    // a possibility we are reading the correct inverted matrix - at least the
    // dimensions match!
    int chkNbEqns, chkNbUnknowns;
    fscanf(fInv, "%d %d\n", &chkNbEqns, &chkNbUnknowns);
    if ((chkNbEqns != NbEqns) || (chkNbUnknowns != NbUnknowns)) {
      neBEMMessage(
          "ReadInvertedMatrix - inverted matrix imension do not match!");
      return (-1);
    }

    printf("ReadInvertedMatrix: Matrix dimensions: %d equations, %d unknowns\n",
           NbEqns, NbUnknowns);
    InvMat = dmatrix(1, NbEqns, 1, NbUnknowns);

    for (int i = 1; i <= NbEqns; i++) {
      printf("%6d", i);  // fflush(stdout);
      for (int j = 1; j <= NbUnknowns; j++) {
        fscanf(fInv, "%le\n", &InvMat[i][j]);
      }
      printf("\b\b\b\b\b\b");
    }

    fclose(fInv);
  }                             // if OptFormattedFile
  else if (OptUnformattedFile)  // both can not be true!
  {
    // not implemented
    neBEMMessage("ReadInvertedMatrix - Binary read not yet implemented.");
    return (-1);
  }

  neBEMState = 7;

  return (0);
}  // end of ReadInvertedMatrix

double ComputeInfluence(int elefld, int elesrc, Point3D *localP,
                        DirnCosn3D *DirCos) {
  if (DebugLevel == 301) {
    printf("In ComputeInfluence ...\n");
  }

  double value;  // influence coefficient

  if (0) {
    printf("\nContinuity satisfaction using following parameters ...\n");
    printf("gtsrc: %d, lxsrc: %lg, lzsrc% lg, dA: %lg\n",
           (EleArr + elesrc - 1)->G.Type, (EleArr + elesrc - 1)->G.LX,
           (EleArr + elesrc - 1)->G.LZ, (EleArr + elesrc - 1)->G.dA);
    printf("xlocal: %lg, ylocal: %lg, zlocal: %lg\n", localP->X, localP->Y,
           localP->Z);
  }

  switch ((EleArr + elefld - 1)
              ->E.Type)  // depending on the etype at the field point
  {                      // different boundary conditions need to be applied
    case 1:              // conductor with known potential
      value = SatisfyValue(elesrc, localP);
      return (value);
      break;

    case 2:  // Conductor with a specified charge - has to have a BC as well!
      printf("Conductors with specific charge not implemented yet.\n");
      return -1;  // The potential has to be the same for the complete
                  // component.
      break;  // If the value of pot is known, the situation is same as above.
              // if not, it is similar to a floating conductor.

    case 3:  // floating conductor
      value = SatisfyValue(elesrc, localP);
      return (value);
      // printf("Floating conductors not implemented yet.\n");
      // return -1;
      break;

      // normal component of the displacement vector is continuous across each
      // dielectric-to-dielectric interface
    case 4:  // DD interface
      value = SatisfyContinuity(elefld, elesrc, localP, DirCos);
      return (value);
      break;

    case 5:  // Dielectric with surface charge density known; same as above BC?
      value = SatisfyContinuity(elefld, elesrc, localP, DirCos);
      return (value);
      // The BC has to be the same and some part of the charge
      break;  // density necessary to satisfy the BC needs to be solved for.

    case 6:  // Symmetry boundary, E parallel
      printf("Symmetry boundary, E parallel not implemented yet. \n");
      return -1;
      break;

    case 7:  // Symmetry boundary, E perpendicular (not likely to be used)
      printf("Symmetry boundary, E perpendicular not implemented yet. \n");
      return -1;
      break;

    default:
      printf("Electric type %d out of range! ... exiting.\n",
             (EleArr + elefld - 1)->E.Type);
      return (-1);
      break;  // unreachable
  }           // switch on etfld ends
}  // end of ComputeInfluence

/* To be tried later
// Compute influence at a point due to a source element / known charged entities
// Pot0Cont1 should better be replaced by (EleArr+elefld-1)->E.Type
double ComputeEleInf(int srcEle, DirnCosn3D *srcDirCos, Point3D *fldPt,
DirnCosn3D *fldDirCos, int Pot0Cont1)
{
if(DebugLevel == 301) { printf("In ComputeEleInf ...\n"); }

double value;	// influence coefficient

if(0)
        {
        printf("\nContinuity satisfaction using following parameters ...\n");
        printf("gtsrc: %d, lxsrc: %lg, lzsrc% lg, dA: %lg\n",
                                        (EleArr+srcEle-1)->G.Type,
(EleArr+srcEle-1)->G.LX, (EleArr+srcEle-1)->G.LZ, (EleArr+srcEle-1)->G.dA);
        printf("xlocal: %lg, ylocal: %lg, zlocal: %lg\n",
                                        localP->X, localP->Y, localP->Z);
        }

switch(Pot0Cont1)	// depending on the etype at the field point
        {			// different boundary conditions need to be
applied case 0:		// conductor with known potential value =
EleSatisfyValue(srcEle, fldPt); return(value); break;

// normal component of the displacement vector is continuous across each
// dielectric-to-dielectric interface
        case 1:		// DD interface
                value = EleSatisfyContinuity(srcEle, srcDirCos, fldPt,
fldDirCos); return(value); break;

        default:
        }	// switch on Pot0Cont1 ends
}	// end of ComputeEleInf
To be tried later */

// Satisfy Dirichlet condition
double SatisfyValue(int elesrc, Point3D *localP) {
  if (DebugLevel == 301) {
    printf("In SatisfyValue ...\n");
  }

  double value;

  switch ((EleArr + elesrc - 1)->G.Type) {
    case 4:  // rectangular element
      value = RecPot(elesrc, localP);
      return (value);
      // return(value/dA);
      break;
    case 3:  // triangular element
      value = TriPot(elesrc, localP);
      return (value);
      // return(value/dA);
      break;
    case 2:  // linear (wire) element
      value = WirePot(elesrc, localP);
      return (value);
      // return(value/dA);
      break;
    default:
      printf("Geometrical type out of range! ... exiting ...\n");
      return (-1);
      break;  // never comes here
  }           // switch over gtsrc ends
}  // end of SatisfyValue

// Satisfy Neumann condition.
// Normal displacement field is continuous across the interface.
// The potential is also continuous, implying continuity of the tangential
// field, but that is not explicitly imposed.
// Slightly modified in V1.8 in comparison to previous versions
// Further modification in V1.8.3:
// Check equation (18) in Bardhan's paper. Note that the D* operator is
// single-layer negative electric field operator (equation 12)
double SatisfyContinuity(int elefld, int elesrc, Point3D *localP,
                         DirnCosn3D *DirCos) {
  if (DebugLevel == 301) {
    printf("In SatisfyContinuity ...\n");
  }

  double value = 0.0;
  Vector3D localF, globalF;

  // For self-influence, lmsrc is equal to lmfld which allows the following.
  // Additional conditions on distances ensure that periodic copies are not
  // treated as `self'.
  // Here, we have the additional problem of correctly treating the
  // self-influence of traingular elements for which the co-ordinate origin
  // does not coincide with the barycentre and thus elesrc and elefld
  // are separated from each other by lx/3 and lz/3 in X and Z,
  // respectively. The Y coordinates match, however, since the element
  // local coordinate system has the Y=0 plane coinciding with the element
  // on which the element is lying.
  // Since elefld and elesrc are the same for self-influence, etsrc can either
  // be 4 or 5. So, the check on the value of etsrc is superfluous, and two
  // separate if blocks are merged into one.
  // Check for other "special" cases!
  if ((elefld == elesrc) &&
      (fabs(localP->X) < (EleArr + elesrc - 1)->G.LX / 2.0) &&
      (fabs(localP->Y) < MINDIST) &&
      (fabs(localP->Z) <
       (EleArr + elesrc - 1)->G.LZ / 2.0))  // self-inf for DD intrfc
  {  // consistent with eqn 18 of Bardhan's paper where lmsrc is inverse
    value = 1.0 / (2.0 * EPS0 * (EleArr + elesrc - 1)->E.Lambda);
  }  // of the multiplying factor of roe(r). EPS0 arises due to electrostatics.
  else {
    value = 0.0;
    // Following fluxes in the influencing ECS
    switch ((EleArr + elesrc - 1)->G.Type) {
      case 4:  // rectangular element
        RecFlux(elesrc, localP, &localF);
        break;
      case 3:  // triangular element
        TriFlux(elesrc, localP, &localF);
        break;
      case 2:  // linear (wire) element
        WireFlux(elesrc, localP, &localF);
        break;
      default:
        printf("Geometrical type out of range! ... exiting ...\n");
        exit(-1);
        break;  // never comes here
    }           // switch over gtsrc ends

    // flux in the global co-ordinate system
    // globalF = RotateVector3D(&localF, &(EleArr+elesrc-1)->G.DC,
    // local2global);
    globalF = RotateVector3D(&localF, DirCos, local2global);

    // Fluxes in the influenced ECS co-ordinate system
    localF =
        RotateVector3D(&globalF, &(EleArr + elefld - 1)->G.DC, global2local);

    value = -localF.Y;
    // normal derivative of Green's function is - normal force
    // (changed from -Fy to +Fy on 02/11/11 - CHECK!!!);
    // From = to +=, further change on 04/11/11; Reverted + to - on 05/11/11
  }  // else self-influence

  return (value);
}  // end of SatisfyContinuity

/* Two functions to be tried later
// Satisfy Dirichlet condition due to a source element at an arbitraty field
point double EleSatisfyValue(int srcEle, Point3D *fldPt)
{
if(DebugLevel == 301) { printf("In EleSatisfyValue ...\n"); }

double value;

switch((EleArr+srcEle-1)->G.Type)
        {
        case 4:		// rectangular element
                value = RecPot(srcEle, fldPt);
                return(value);
                // return(value/dA);
                break;
        case 3:		// triangular element
                value = TriPot(srcEle, fldPt);
                return(value);
                // return(value/dA);
                break;
        case 2:		// linear (wire) element
                value = WirePot(srcEle, fldPt);
                return(value);
                // return(value/dA);
                break;
        default:
                printf("Geometrical type out of range! ... exiting ...\n");
                return(-1);
                break;	// never comes here
        }	// switch over gtsrc ends
}	// end of EleSatisfyValue


// Satisfy Neumann condition due to a source element at an arbitrary field point
// Normal displacement field is continuous across the interface.
// The potential is also continuous, implying continuity of the tangential
// field, but that is not explicitly imposed.
// Slightly modified in V1.8 in comparison to previous versions
// Further modification in V1.8.3:
// Check equation (18) in Bardhan's paper. Note that the D* operator is
// single-layer negative electric field operator (equation 12)
double EleSatisfyContinuity(int srcEle, DirnCosn3D *eleDirCos, Point3D *fldPt,
DirnCosn3D *fldDirCos)
{
if(DebugLevel == 301) { printf("In SatisfyContinuity ...\n"); }

double value = 0.0;
Vector3D localF, globalF;

// For self-influence, lmsrc is equal to lmfld which allows the following.
// Additional conditions on distances ensure that periodic copies are not
// treated as `self'.
// Here, we have the additional problem of correctly treating the
// self-influence of traingular elements for which the co-ordinate origin
// does not coincide with the barycentre and thus srcEle and elefld
// are separated from each other by lx/3 and lz/3 in X and Z,
// respectively. The Y coordinates match, however, since the element
// local coordinate system has the Y=0 plane coinciding with the element
// on which the element is lying.
// Since elefld and srcEle are the same for self-influence, etsrc can either be
// 4 or 5. So, the check on the value of etsrc is superfluous, and two
// separate if blocks are merged into one.
// Check for other "special" cases!
if((fabs(fldPt->X) < (EleArr+srcEle-1)->G.LX/2.0)
                && (fabs(fldPt->Y) < MINDIST)
                && (fabs(fldPt->Z) < (EleArr+srcEle-1)->G.LZ/2.0))// self-inf
for DD intrfc {	// consistent with eqn 18 of Bardhan's paper where lmsrc is
inverse value = 1.0 / (2.0*EPS0*(EleArr+srcEle-1)->E.Lambda); }	// of the
multiplying factor of roe(r). EPS0 arises due to electrostatics. else
        {
        value = 0.0;
        // Following fluxes in the influencing ECS
        switch((EleArr+srcEle-1)->G.Type)
                {
                case 4:		// rectangular element
                        RecFlux(srcEle, fldPt, &localF);
                break;
                case 3:		// triangular element
                        TriFlux(srcEle, fldPt, &localF);
                break;
                case 2:		// linear (wire) element
                        WireFlux(srcEle, fldPt, &localF);
                break;
                default:
                        printf("Geometrical type out of range! ... exiting
...\n"); exit(-1); break;	// never comes here }	// switch over gtsrc
ends

        // flux in the global co-ordinate system
        // globalF = RotateVector3D(&localF, &(EleArr+srcEle-1)->G.DC,
local2global); globalF = RotateVector3D(&localF, srcDirCos, local2global);

        // Fluxes in the influenced ECS co-ordinate system
        localF = RotateVector3D(&globalF, fldDirCos, global2local);

        value = -localF.Y;
        // normal derivative of Green's function is - normal force
        // (changed from -Fy to +Fy on 02/11/11 - CHECK!!!);
        // From = to +=, further change on 04/11/11; Reverted + to - on 05/11/11
        }	// else self-influence

return(value);
}	// end of EleSatisfyContinuity
Two functions to be tried later */

// The RHVector is specified by the boundary conditions at the field points
// Ideally, there should be a flag related to the presence of known charges.
// If no known charges are there, it is not necessary to try out functions
// ValueKnCh and ContinuityKnCh. However, now these functions are always used.
// Similar functions for ChargingUp are also being used.
// Check the note in neBEMInterface.c regarding these two functions.
int RHVector(void) {
  double value, valueKnCh, valueChUp;
  char outfile[256];
  FILE *fout;

  if (TimeStep == 1) RHS = dvector(1, NbEqns);

  strcpy(outfile, BCOutDir);
  strcat(outfile, "/BCondns.out");
  fout = fopen(outfile, "w");
  fprintf(fout, "#BCondn Vector\n");
  fprintf(fout, "#elefld\tAssigned\tBC\tKnCh\tChUp\tRHValue\n");

  printf("created BCondns.out file ...\n");
  fflush(stdout);

  for (int elefld = 1; elefld <= NbElements; ++elefld) {
    if (0) printf("\nIn RHVector, elefld: %d\n", elefld);
    value = valueKnCh = valueChUp = 0.0;
    value = (EleArr + elefld - 1)
                ->BC.Value;  // previouly this line was within case 1

    switch ((EleArr + elefld - 1)->E.Type) {
      case 1:  // Conducting surfaces
        // value = (EleArr+elefld-1)->BC.Value;
        if (OptKnCh) {
          valueKnCh = ValueKnCh(elefld);  // effect of all known charges
          if (isnan(valueKnCh)) exit(-1);
          if (isinf(valueKnCh)) exit(-1);
        }
        if (OptChargingUp) {
          valueChUp = ValueChUp(elefld);  // effect of device charging up
          if (isnan(valueChUp)) exit(-1);
          if (isinf(valueChUp)) exit(-1);
        }
        RHS[elefld] = value - valueKnCh - valueChUp;
        break;
      case 2:  // Conducting surfaces with known charge
        printf("Conducting surface with charge not implemented.\n");
        return -1;
        break;  // NOTE: no RHVector
      case 3:   // Floating conducting surfaces
        if (OptKnCh) {
          valueKnCh = ValueKnCh(elefld);
          if (isnan(valueKnCh)) exit(-1);
          if (isinf(valueKnCh)) exit(-1);
        }
        if (OptChargingUp) {
          valueChUp = ValueChUp(elefld);
          if (isnan(valueChUp)) exit(-1);
          if (isinf(valueChUp)) exit(-1);
        }
        RHS[elefld] = value - valueKnCh - valueChUp;
        break;
      case 4:  // Dielectric interfaces
        if (OptKnCh) {
          valueKnCh = ContinuityKnCh(elefld);
          if (isnan(valueKnCh)) exit(-1);
          if (isinf(valueKnCh)) exit(-1);
        }
        if (OptChargingUp) {
          valueChUp = ContinuityChUp(elefld);
          if (isnan(valueChUp)) exit(-1);
          if (isinf(valueChUp)) exit(-1);
        }
        RHS[elefld] = value - valueKnCh - valueChUp;
        RHS[elefld] +=
            (EleArr + elefld - 1)->Assigned;  // effect due to charging up
        break;
      case 5:  // Dielectric interfaces with known charge
        if (OptKnCh) {
          valueKnCh = ContinuityKnCh(elefld);
          if (isnan(valueKnCh)) exit(-1);
          if (isinf(valueKnCh)) exit(-1);
        }
        if (OptChargingUp) {
          valueChUp = ContinuityChUp(elefld);
          if (isnan(valueChUp)) exit(-1);
          if (isinf(valueChUp)) exit(-1);
        }
        RHS[elefld] = value - valueKnCh - valueChUp;  // Check Bardhan's eqn 16
        RHS[elefld] +=
            (EleArr + elefld - 1)->Assigned;  // effect due to assigned
        // charge; what happens when assigned elements are charged up in
        // addition?
        break;
      case 6:  // E parallel symmetry boundary
        printf("Symmetry boundary, E parallel not implemented yet.\n");
        return -1;
        break;
      case 7:  // E perpendicular symmetry boundary
        printf("Symmetry boundary, E perpendicular not implemented yet.\n");
        return -1;
        break;
      default:
        printf("elefld in RHVector out of range ... returning\n");
        return -1;
    }
    fprintf(fout, "%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n", elefld,
            (EleArr + elefld - 1)->G.Origin.X,
            (EleArr + elefld - 1)->G.Origin.Y,
            (EleArr + elefld - 1)->G.Origin.Z, (EleArr + elefld - 1)->Assigned,
            value, valueKnCh, valueChUp, RHS[elefld]);
  }  // for elefld ends

  // Floating conductor are long neglected. Needs close inspection.
  if (NbConstraints) {
    for (int eqn = NbElements + 1; eqn <= NbEqns; ++eqn) {
      if (eqn ==
          NbSystemChargeZero)  // if row corresponds to zero-charge condition
      {
        double assigned, dA, SumAssigned = 0.0;
        // can be parallelized
        for (int ele = 1; ele <= NbElements; ++ele) {
          assigned = (EleArr + ele - 1)->Assigned;
          dA = (EleArr + ele - 1)->G.dA;
          SumAssigned += assigned * dA;  // charge density * element area
        }
        RHS[eqn] =
            -SumAssigned;  // applied charge considered - too small change
        RHS[eqn] =
            0.0;  // applied charge not considered while meeting constraint
      }           // if eqn == NbSystemChargeZero
      else {
        RHS[eqn] = 0.0;
      }
    }
  }  // if NbConstraints

  printf("computations for RHVector completed ...\n");
  fflush(stdout);

  fclose(fout);

  neBEMState = 5;

  return (0);
}  // RHVector

// Dirichlet contribution due to known charge distributions.
// For the following two to work, all the geometrical attributes (repetitions,
// reflections etc.) need to be considered. This is true for the point, line
// and surface charges, as well. But can these known charges be assumed to
// repeat at all?
double ValueKnCh(int elefld) {
  int dbgFn = 0;

  if (dbgFn) printf("\nelefld: %d\n", elefld);

  double value = 0.0;
  double assigned = 0.0;
  double xfld = (EleArr + elefld - 1)->BC.CollPt.X;
  double yfld = (EleArr + elefld - 1)->BC.CollPt.Y;
  double zfld = (EleArr + elefld - 1)->BC.CollPt.Z;

  // Retrieve element properties at the field point
  // Location needed for Dirichlet (potential)

  // OMPCheck - parallelize

  // Effect of known charges on the interface elements
  for (int elesrc = 1; elesrc <= NbElements; ++elesrc) {
    Point3D localP;

    assigned = (EleArr + elesrc - 1)->Assigned;
    if (fabs(assigned) <= 1.0e-16) continue;

    // Retrieve element properties from the structure
    if ((EleArr + elesrc - 1)->E.Type == 0) {
      printf("Wrong EType for elesrc %d element %d on %dth primitive!\n",
             elesrc, (EleArr + elesrc - 1)->Id,
             (EleArr + elesrc - 1)->PrimitiveNb);
      exit(-1);
    }

    Point3D *pOrigin = &(EleArr + elesrc - 1)->G.Origin;

    {  // Rotate point3D from global to local system
      double InitialVector[3] =
             {xfld - pOrigin->X, yfld - pOrigin->Y, zfld - pOrigin->Z};
      double TransformationMatrix[3][3] = {{0.0, 0.0, 0.0},
                                           {0.0, 0.0, 0.0},
                                           {0.0, 0.0, 0.0}};
      DirnCosn3D *DirCos = &(EleArr + elesrc - 1)->G.DC;
      TransformationMatrix[0][0] = DirCos->XUnit.X;
      TransformationMatrix[0][1] = DirCos->XUnit.Y;
      TransformationMatrix[0][2] = DirCos->XUnit.Z;
      TransformationMatrix[1][0] = DirCos->YUnit.X;
      TransformationMatrix[1][1] = DirCos->YUnit.Y;
      TransformationMatrix[1][2] = DirCos->YUnit.Z;
      TransformationMatrix[2][0] = DirCos->ZUnit.X;
      TransformationMatrix[2][1] = DirCos->ZUnit.Y;
      TransformationMatrix[2][2] = DirCos->ZUnit.Z;
      double FinalVector[3];

      for (int i = 0; i < 3; ++i) {
        FinalVector[i] = 0.0;
        for (int j = 0; j < 3; ++j) {
          FinalVector[i] += TransformationMatrix[i][j] * InitialVector[j];
        }
      }

      localP.X = FinalVector[0];
      localP.Y = FinalVector[1];
      localP.Z = FinalVector[2];
    }  // Point3D rotated

    switch ((EleArr + elesrc - 1)->G.Type) {
      case 4:  // rectangular element
        value += assigned * RecPot(elesrc, &localP);
        // return(value/dA);
        break;
      case 3:  // triangular element
        value += assigned * TriPot(elesrc, &localP);
        // return(value/dA);
        break;
      case 2:  // linear (wire) element
        value += assigned * WirePot(elesrc, &localP);
        // return(value/dA);
        break;
      default:
        printf("Geometrical type out of range! ... exiting ...\n");
        exit(-1);
        break;  // never comes here
    }           // switch over gtsrc ends
  }             // for all source elements - elesrc

  value += MyFACTOR;  // in order reduce divisions by MyFACTOR later
  if (dbgFn) {
    printf("value after known charges on elements (* MyFACTOR): %g\n", value);
  }

  // globalP is now required with a different definition
  Point3D fieldPt;
  fieldPt.X = xfld;
  fieldPt.Y = yfld;
  fieldPt.Z = zfld;

  // Contribution due to known point charges
  Vector3D tmpglobalF;  // flux not being used to evaluate Dirichlet value.
  for (int point = 1; point <= NbPointsKnCh; ++point) {
    Point3D sourcePt;
    sourcePt.X = PointKnChArr[point].P.X;
    sourcePt.Y = PointKnChArr[point].P.Y;
    sourcePt.Z = PointKnChArr[point].P.Z;
    value += PointKnChArr[point].Assigned *
             PointKnChPF(sourcePt, fieldPt, &tmpglobalF);
  }  // for all points
  if (dbgFn) {
    printf("value after known charges at points: %g\n", value);
  }

  // Contribution due to known line charges
  for (int line = 1; line <= NbLinesKnCh; ++line) {
    Point3D startPt, stopPt;
    startPt.X = LineKnChArr[line].Start.X;
    startPt.Y = LineKnChArr[line].Start.Y;
    startPt.Z = LineKnChArr[line].Start.Z;
    stopPt.X = LineKnChArr[line].Stop.X;
    stopPt.Y = LineKnChArr[line].Stop.Y;
    stopPt.Z = LineKnChArr[line].Stop.Z;
    double radius = LineKnChArr[line].Radius;
    if (dbgFn) {
      printf(
          "In ValueKnCh, Nb: %d, X1Y1Z1: %lg %lg %lg X2Y2Z2 %lg %lg %lg R: %lg "
          "Lmda:%lg\n",
          LineKnChArr[line].Nb, startPt.X, startPt.Y, startPt.Z, stopPt.X,
          stopPt.Y, stopPt.Z, radius, LineKnChArr[line].Assigned);
    }
    value += LineKnChArr[line].Assigned *
             LineKnChPF(startPt, stopPt, fieldPt, &tmpglobalF);
  }  // for lines
  if (dbgFn) {
    printf("value after known charges on lines: %g\n", value);
  }

  // Contribution due to known area charges
  // We may need to convert the information to simpler structures, as done
  // for points and lines, above.
  for (int area = 1; area <= NbAreasKnCh; ++area) {
    value +=
        (AreaKnChArr + area - 1)->Assigned *
        AreaKnChPF((AreaKnChArr + area - 1)->NbVertices,
                   ((AreaKnChArr + area - 1)->Vertex), fieldPt, &tmpglobalF);
  }  // for areas
  if (dbgFn) {
    printf("value after known charges on areas: %g\n", value);
  }

  // Contribution due to known volume charges
  for (int vol = 1; vol <= NbVolumesKnCh; ++vol) {
    value += (VolumeKnChArr + vol - 1)->Assigned *
             VolumeKnChPF((VolumeKnChArr + vol - 1)->NbVertices,
                          ((VolumeKnChArr + vol - 1)->Vertex), fieldPt,
                          &tmpglobalF);
  }  // for vols

  value /= MyFACTOR;
  if (dbgFn) {
    printf("value after known charges in volumes (/ MyFACTOR): %g\n", value);
  }

  return (value);
}  // end of ValueKnCh

// Neumann contribution due to known charge distribution.
// Note that the return value is negated from the BC which was set up without
// considering the presence of existing charges.
// Check dielectric-dielectric formulation in
// (NumSolnOfBIEforMolES_JPBardhan.pdf):
// Numerical solution of boundary-integral equations for molecular
// electrostatics,
// by Jaydeep P. Bardhan,
// THE JOURNAL OF CHEMICAL PHYSICS 130, 094102 (2009)
double ContinuityKnCh(int elefld) {
  int dbgFn = 0;

  if (dbgFn) printf("\nelefld: %d\n", elefld);

  double value = 0.0;
  double assigned = 0.0;
  double xfld = (EleArr + elefld - 1)->BC.CollPt.X;
  double yfld = (EleArr + elefld - 1)->BC.CollPt.Y;
  double zfld = (EleArr + elefld - 1)->BC.CollPt.Z;

  Vector3D localF, globalF;

  // OMPCheck - parallelize

  // Effect of known charges on interface elements
  for (int elesrc = 1; elesrc <= NbElements; ++elesrc) {
    Point3D localP;

    assigned = (EleArr + elesrc - 1)->Assigned;
    if (fabs(assigned) <= 1.0e-16) continue;

    Point3D *pOrigin = &(EleArr + elesrc - 1)->G.Origin;

    // Retrieve element properties from the structure
    if ((EleArr + elesrc - 1)->E.Type == 0) {
      printf("Wrong EType for elesrc %d element %d on %dth primitive!\n",
             elesrc, (EleArr + elesrc - 1)->Id,
             (EleArr + elesrc - 1)->PrimitiveNb);
      exit(-1);
    }

    {  // Rotate point3D from global to local system
      double InitialVector[3]
                = {xfld - pOrigin->X, yfld - pOrigin->Y, zfld - pOrigin->Z};
      double TransformationMatrix[3][3] = {{0.0, 0.0, 0.0},
                                           {0.0, 0.0, 0.0},
                                           {0.0, 0.0, 0.0}};
      DirnCosn3D *DirCos = &(EleArr + elesrc - 1)->G.DC;
      TransformationMatrix[0][0] = DirCos->XUnit.X;
      TransformationMatrix[0][1] = DirCos->XUnit.Y;
      TransformationMatrix[0][2] = DirCos->XUnit.Z;
      TransformationMatrix[1][0] = DirCos->YUnit.X;
      TransformationMatrix[1][1] = DirCos->YUnit.Y;
      TransformationMatrix[1][2] = DirCos->YUnit.Z;
      TransformationMatrix[2][0] = DirCos->ZUnit.X;
      TransformationMatrix[2][1] = DirCos->ZUnit.Y;
      TransformationMatrix[2][2] = DirCos->ZUnit.Z;
      double FinalVector[3];

      for (int i = 0; i < 3; ++i) {
        FinalVector[i] = 0.0;
        for (int j = 0; j < 3; ++j) {
          FinalVector[i] += TransformationMatrix[i][j] * InitialVector[j];
        }
      }

      localP.X = FinalVector[0];
      localP.Y = FinalVector[1];
      localP.Z = FinalVector[2];
    }  // Point3D rotated

    // if((((EleArr+elesrc-1)->E.Type == 4) && (elefld == elesrc))  //
    // self-influence
    // || (((EleArr+elesrc-1)->E.Type == 5) && (elefld == elesrc))) // DD intfce
    // For self-influence, lmsrc is equal to lmfld which allows the following.
    if ((elefld == elesrc) &&
        (fabs(localP.X) < (EleArr + elesrc - 1)->G.LX / 2.0) &&
        (fabs(localP.Y) < MINDIST) &&
        (fabs(localP.Z) < (EleArr + elesrc - 1)->G.LZ / 2.0)) {
      value += assigned * 1.0 / (2.0 * EPS0 * (EleArr + elesrc - 1)->E.Lambda);
    } else {
      // Retrieve element properties from the structure
      if ((EleArr + elesrc - 1)->E.Type == 0) {
        printf("Wrong EType for elesrc %d element %d on %dth primitive!\n",
               elesrc, (EleArr + elesrc - 1)->Id,
               (EleArr + elesrc - 1)->PrimitiveNb);
        exit(-1);
      }

      // Following fluxes in the influencing ECS
      switch ((EleArr + elesrc - 1)->G.Type) {
        case 4:  // rectangular element
          RecFlux(elesrc, &localP, &localF);
          break;
        case 3:  // triangular element
          TriFlux(elesrc, &localP, &localF);
          break;
        case 2:  // linear (wire) element
          WireFlux(elesrc, &localP, &localF);
          break;
        default:
          printf("Geometrical type out of range! ... exiting ...\n");
          exit(-1);
          break;  // never comes here
      }           // switch over gtsrc ends

      // in GCS - mirror points?!
      globalF =
          RotateVector3D(&localF, &(EleArr + elesrc - 1)->G.DC, local2global);
      // in ECS (local to the field element)
      localF =
          RotateVector3D(&globalF, &(EleArr + elefld - 1)->G.DC, global2local);

      value -= assigned * localF.Y;  // +ve gradient of Green's is -ve normal
                                     // force (changed from -= to += on 02/11/11
                                     // - CHECK!!! - and back to -= on 05/11/11)
    }                                // else self-influence
  }                                  // for elesrc

  value *= MyFACTOR; // so that later MyFACTOR divisions are reduced.
  if (dbgFn) {
    printf("value (* MyFACTOR): %g\n", value);
  }

  // Contributions from points, lines areas and volumes carrying known charge
  // (density).
  Point3D fieldPt;
  fieldPt.X = xfld;
  fieldPt.Y = yfld;  // the global position of the field point necessary now
  fieldPt.Z = zfld;
  for (int point = 1; point <= NbPointsKnCh; ++point) {
    Point3D sourcePt;
    sourcePt.X = PointKnChArr[point].P.X;
    sourcePt.Y = PointKnChArr[point].P.Y;
    sourcePt.Z = PointKnChArr[point].P.Z;
    // potential not being used to evaluate Neumann continuity
    (void)PointKnChPF(sourcePt, fieldPt, &globalF);
    localF =
        RotateVector3D(&globalF, &(EleArr + elefld - 1)->G.DC, global2local);
    value -= PointKnChArr[point].Assigned * localF.Y;
  }  // loop over points NbPointsKnCh

  for (int line = 1; line <= NbLinesKnCh; ++line) {
    Point3D startPt, stopPt;
    startPt.X = LineKnChArr[line].Start.X;
    startPt.Y = LineKnChArr[line].Start.Y;
    startPt.Z = LineKnChArr[line].Start.Z;
    stopPt.X = LineKnChArr[line].Stop.X;
    stopPt.Y = LineKnChArr[line].Stop.Y;
    stopPt.Z = LineKnChArr[line].Stop.Z;
    double radius = LineKnChArr[line].Radius;
    if (0) {
      printf(
          "In ContinuityKnCh, Nb: %d, X1Y1Z1: %lg %lg %lg X2Y2Z2 %lg %lg %lg "
          "R: %lg Lmda:%lg\n",
          LineKnChArr[line].Nb, startPt.X, startPt.Y, startPt.Z, stopPt.X,
          stopPt.Y, stopPt.Z, radius, LineKnChArr[line].Assigned);
    }
    // potential not being used to evaluate Neumann continuity
    (void)LineKnChPF(startPt, stopPt, fieldPt, &globalF);
    localF =
        RotateVector3D(&globalF, &(EleArr + elefld - 1)->G.DC, global2local);
    value -= LineKnChArr[line].Assigned * localF.Y;
  }  // loop over lines
  // Simplifying conversions, similar to points and lines may be necessary.
  for (int area = 1; area <= NbAreasKnCh; ++area) {
    // potential not being used to evaluate Neumann continuity
    (void)AreaKnChPF((AreaKnChArr + area - 1)->NbVertices,
                   ((AreaKnChArr + area - 1)->Vertex), fieldPt, &globalF);
    localF =
        RotateVector3D(&globalF, &(EleArr + elefld - 1)->G.DC, global2local);
    value -= (AreaKnChArr + area - 1)->Assigned * localF.Y;
  }  // loop over areas

  // Simplifying conversions, similar to points and lines may be necessary.
  for (int vol = 1; vol <= NbVolumesKnCh; ++vol) {
    // potential not being used to evaluate Neumann continuity
    (void)VolumeKnChPF((VolumeKnChArr + vol - 1)->NbVertices,
                     ((VolumeKnChArr + vol - 1)->Vertex), fieldPt, &globalF);
    localF =
        RotateVector3D(&globalF, &(EleArr + elefld - 1)->G.DC, global2local);
    value -= (VolumeKnChArr + vol - 1)->Assigned * localF.Y;
  }  // loop over volumes

  value /= MyFACTOR; // factored in
  return (value);
}  // end of ContinuityKnCh

// Effect of charging up on the Dirichlet boundary conditions
double ValueChUp(int elefld) {
  int dbgFn = 0;

  if (dbgFn) printf("\nelefld: %d\n", elefld);
  if (dbgFn) {
    printf("In ValueChUp ...\n");
  }

  double value = 0.0;
  double assigned = 0.0;
  double xfld = (EleArr + elefld - 1)->BC.CollPt.X;
  double yfld = (EleArr + elefld - 1)->BC.CollPt.Y;
  double zfld = (EleArr + elefld - 1)->BC.CollPt.Z;

  // Retrieve element properties at the field point
  // Location needed for Dirichlet (potential)

  // OMPCheck - parallelize

  // Effect of known charges on the interface elements
  for (int elesrc = 1; elesrc <= NbElements; ++elesrc) {
    Point3D localP;

    assigned = (EleArr + elesrc - 1)->Assigned;
    if (fabs(assigned) <= 1.0e-16) continue;

    // Retrieve element properties from the structure
    if ((EleArr + elesrc - 1)->E.Type == 0) {
      printf("Wrong EType for elesrc %d element %d on %dth primitive!\n",
             elesrc, (EleArr + elesrc - 1)->Id,
             (EleArr + elesrc - 1)->PrimitiveNb);
      exit(-1);
    }

    Point3D *pOrigin = &(EleArr + elesrc - 1)->G.Origin;

    {  // Rotate point3D from global to local system
      double InitialVector[3]
                  = {xfld - pOrigin->X, yfld - pOrigin->Y, zfld - pOrigin->Z};
      double TransformationMatrix[3][3] = {{0.0, 0.0, 0.0},
                                           {0.0, 0.0, 0.0},
                                           {0.0, 0.0, 0.0}};
      DirnCosn3D *DirCos = &(EleArr + elesrc - 1)->G.DC;
      TransformationMatrix[0][0] = DirCos->XUnit.X;
      TransformationMatrix[0][1] = DirCos->XUnit.Y;
      TransformationMatrix[0][2] = DirCos->XUnit.Z;
      TransformationMatrix[1][0] = DirCos->YUnit.X;
      TransformationMatrix[1][1] = DirCos->YUnit.Y;
      TransformationMatrix[1][2] = DirCos->YUnit.Z;
      TransformationMatrix[2][0] = DirCos->ZUnit.X;
      TransformationMatrix[2][1] = DirCos->ZUnit.Y;
      TransformationMatrix[2][2] = DirCos->ZUnit.Z;
      double FinalVector[3];

      for (int i = 0; i < 3; ++i) {
        FinalVector[i] = 0.0;
        for (int j = 0; j < 3; ++j) {
          FinalVector[i] += TransformationMatrix[i][j] * InitialVector[j];
        }
      }

      localP.X = FinalVector[0];
      localP.Y = FinalVector[1];
      localP.Z = FinalVector[2];
    }  // Point3D rotated

    switch ((EleArr + elesrc - 1)->G.Type) {
      case 4:  // rectangular element
        value += assigned * RecPot(elesrc, &localP);
        // return(value/dA);
        break;
      case 3:  // triangular element
        value += assigned * TriPot(elesrc, &localP);
        // return(value/dA);
        break;
      case 2:  // linear (wire) element
        value += assigned * WirePot(elesrc, &localP);
        // return(value/dA);
        break;
      default:
        printf("Geometrical type out of range! ... exiting ...\n");
        exit(-1);
        break;  // never comes here
    }           // switch over gtsrc ends
  }             // for all source elements - elesrc

  value *= MyFACTOR;
  if (dbgFn) {
    printf("value after known charges on elements (*MyFACTOR): %g\n", value);
  }

  // globalP is now required with a different definition
  Point3D globalP;
  globalP.X = xfld;
  globalP.Y = yfld;
  globalP.Z = zfld;

  // Contribution due to known point charges
  Vector3D tmpglobalF;  // flux not being used to evaluate Dirichlet value.
  for (int point = 1; point <= NbPointsKnCh; ++point) {
    value += (PointKnChArr + point - 1)->Assigned *
             PointKnChPF((PointKnChArr + point - 1)->P, globalP, &tmpglobalF);
  }  // for all points
  if (dbgFn) {
    printf("value after known charges at points: %g\n", value);
  }

  // Contribution due to known line charges
  for (int line = 1; line <= NbLinesKnCh; ++line) {
    value += (LineKnChArr + line - 1)->Assigned *
             LineKnChPF((LineKnChArr + line - 1)->Start,
                        (LineKnChArr + line - 1)->Stop, globalP, &tmpglobalF);
  }  // for lines
  if (dbgFn) {
    printf("value after known charges on lines: %g\n", value);
  }

  // Contribution due to known area charges
  for (int area = 1; area <= NbAreasKnCh; ++area) {
    value +=
        (AreaKnChArr + area - 1)->Assigned *
        AreaKnChPF((AreaKnChArr + area - 1)->NbVertices,
                   ((AreaKnChArr + area - 1)->Vertex), globalP, &tmpglobalF);
  }  // for areas
  if (dbgFn) {
    printf("value after known charges on areas: %g\n", value);
  }

  // Contribution due to known volume charges
  for (int vol = 1; vol <= NbVolumesKnCh; ++vol) {
    value += (VolumeKnChArr + vol - 1)->Assigned *
             VolumeKnChPF((VolumeKnChArr + vol - 1)->NbVertices,
                          ((VolumeKnChArr + vol - 1)->Vertex), globalP,
                          &tmpglobalF);
  }  // for vols

  value /= MyFACTOR;
  if (dbgFn) {
    printf("value after known charges in volumes (/ MyFACTOR): %g\n", value);
    printf("Exiting ValueChUp ...\n");
  }

  return (value);
}  // ValueChUp ends

// Effect of charging up on the Neumann boundary condition
double ContinuityChUp(int elefld) {
  int dbgFn = 0;

  if (dbgFn) printf("\nelefld: %d\n", elefld);
  if (dbgFn) {
    printf("In ContinuityChUp ...\n");
  }

  double value = 0.0;
  double assigned = 0.0;
  double xfld = (EleArr + elefld - 1)->BC.CollPt.X;
  double yfld = (EleArr + elefld - 1)->BC.CollPt.Y;
  double zfld = (EleArr + elefld - 1)->BC.CollPt.Z;

  Vector3D localF, globalF;

  // OMPCheck - parallelize

  // Effect of known charges on interface elements
  for (int elesrc = 1; elesrc <= NbElements; ++elesrc) {
    Point3D localP;

    assigned = (EleArr + elesrc - 1)->Assigned;
    if (fabs(assigned) <= 1.0e-16) continue;

    Point3D *pOrigin = &(EleArr + elesrc - 1)->G.Origin;

    // Retrieve element properties from the structure
    if ((EleArr + elesrc - 1)->E.Type == 0) {
      printf("Wrong EType for elesrc %d element %d on %dth primitive!\n",
             elesrc, (EleArr + elesrc - 1)->Id,
             (EleArr + elesrc - 1)->PrimitiveNb);
      exit(-1);
    }

    {  // Rotate point3D from global to local system
      double InitialVector[3]
                = {xfld - pOrigin->X, yfld - pOrigin->Y, zfld - pOrigin->Z};
      double TransformationMatrix[3][3] = {{0.0, 0.0, 0.0},
                                           {0.0, 0.0, 0.0},
                                           {0.0, 0.0, 0.0}};
      DirnCosn3D *DirCos = &(EleArr + elesrc - 1)->G.DC;
      TransformationMatrix[0][0] = DirCos->XUnit.X;
      TransformationMatrix[0][1] = DirCos->XUnit.Y;
      TransformationMatrix[0][2] = DirCos->XUnit.Z;
      TransformationMatrix[1][0] = DirCos->YUnit.X;
      TransformationMatrix[1][1] = DirCos->YUnit.Y;
      TransformationMatrix[1][2] = DirCos->YUnit.Z;
      TransformationMatrix[2][0] = DirCos->ZUnit.X;
      TransformationMatrix[2][1] = DirCos->ZUnit.Y;
      TransformationMatrix[2][2] = DirCos->ZUnit.Z;
      double FinalVector[3];

      for (int i = 0; i < 3; ++i) {
        FinalVector[i] = 0.0;
        for (int j = 0; j < 3; ++j) {
          FinalVector[i] += TransformationMatrix[i][j] * InitialVector[j];
        }
      }

      localP.X = FinalVector[0];
      localP.Y = FinalVector[1];
      localP.Z = FinalVector[2];
    }  // Point3D rotated

    // if((((EleArr+elesrc-1)->E.Type == 4) && (elefld == elesrc))  //
    // self-influence
    // || (((EleArr+elesrc-1)->E.Type == 5) && (elefld == elesrc))) // DD intfce
    // For self-influence, lmsrc is equal to lmfld which allows the following.
    if ((elefld == elesrc) &&
        (fabs(localP.X) < (EleArr + elesrc - 1)->G.LX / 2.0) &&
        (fabs(localP.Y) < MINDIST) &&
        (fabs(localP.Z) < (EleArr + elesrc - 1)->G.LZ / 2.0)) {
      value += assigned * 1.0 / (2.0 * EPS0 * (EleArr + elesrc - 1)->E.Lambda);
    } else {
      // Retrieve element properties from the structure
      if ((EleArr + elesrc - 1)->E.Type == 0) {
        printf("Wrong EType for elesrc %d element %d on %dth primitive!\n",
               elesrc, (EleArr + elesrc - 1)->Id,
               (EleArr + elesrc - 1)->PrimitiveNb);
        exit(-1);
      }

      switch ((EleArr + elesrc - 1)->G.Type) {
        case 4:  // rectangular element
          RecFlux(elesrc, &localP, &localF);
          break;
        case 3:  // triangular element
          TriFlux(elesrc, &localP, &localF);
          break;
        case 2:  // linear (wire) element
          WireFlux(elesrc, &localP, &localF);
          break;
        default:
          printf("Geometrical type out of range! ... exiting ...\n");
          exit(-1);
          break;  // never comes here
      }           // switch over gtsrc ends

      // in GCS - mirror points?!
      globalF =
          RotateVector3D(&localF, &(EleArr + elesrc - 1)->G.DC, local2global);
      // in ECS (local to the field element)
      localF =
          RotateVector3D(&globalF, &(EleArr + elefld - 1)->G.DC, global2local);

      value -= assigned * localF.Y;  // +ve gradient of Green's is -ve normal
                                     // force (changed from -= to += on 02/11/11
                                     // - CHECK!!! - and back to -= on 05/11/11)
    }                                // else self-influence
  }

  value *= MyFACTOR;
  if (dbgFn) {
    printf("value (* MyFACTOR): %g\n", value);
  }

  // Contributions from points, lines areas and volumes carrying known charge
  // (density).
  Point3D globalP;
  globalP.X = xfld;
  globalP.Y = yfld;  // the global position of the field point necessary now
  globalP.Z = zfld;
  for (int point = 1; point <= NbPointsKnCh; ++point) {
    // potential not being used to evaluate Neumann continuity
    (void)PointKnChPF((PointKnChArr + point - 1)->P, globalP, &globalF);
    localF =
        RotateVector3D(&globalF, &(EleArr + elefld - 1)->G.DC, global2local);
    value -= (PointKnChArr + point - 1)->Assigned * localF.Y;
  }  // loop over points NbPointsKnCh

  for (int line = 1; line <= NbLinesKnCh; ++line) {
    // potential not being used to evaluate Neumann continuity
    (void)LineKnChPF((LineKnChArr + line - 1)->Start,
                   (LineKnChArr + line - 1)->Stop, globalP, &globalF);
    localF =
        RotateVector3D(&globalF, &(EleArr + elefld - 1)->G.DC, global2local);
    value -= (LineKnChArr + line - 1)->Assigned * localF.Y;
  }  // loop over lines

  for (int area = 1; area <= NbAreasKnCh; ++area) {
    // potential not being used to evaluate Neumann continuity
    (void)AreaKnChPF((AreaKnChArr + area - 1)->NbVertices,
                   ((AreaKnChArr + area - 1)->Vertex), globalP, &globalF);
    localF =
        RotateVector3D(&globalF, &(EleArr + elefld - 1)->G.DC, global2local);
    value -= (AreaKnChArr + area - 1)->Assigned * localF.Y;
  }  // loop over areas

  for (int vol = 1; vol <= NbVolumesKnCh; ++vol) {
    // potential not being used to evaluate Neumann continuity
    (void)VolumeKnChPF((VolumeKnChArr + vol - 1)->NbVertices,
                     ((VolumeKnChArr + vol - 1)->Vertex), globalP, &globalF);
    localF =
        RotateVector3D(&globalF, &(EleArr + elefld - 1)->G.DC, global2local);
    value -= (VolumeKnChArr + vol - 1)->Assigned * localF.Y;
  }  // loop over volumes

  value /= MyFACTOR;
  if (dbgFn) {
    printf("value after known charges in volumes (/ MyFACTOR): %g\n", value);
    printf("Exiting ContinuityChUp ...\n");
  }

  return (value);
}  // ContinuityChUp ends

// An error estimate should be carried out irrespective of the DebugLevel
int Solve(void) {
  double **RawInf = NULL;
  int LHMatrix(void);

  Solution = dvector(1, NbUnknowns);
  // printf("Solution array allocated\n"); fflush(stdout);

  printf("\ncomputing solution for each element: ");

  for (int i = 1; i <= NbUnknowns; i++) {
    printf("%6d ...", i);
    Solution[i] = 0.0;

    double sum = 0.0;
    int j;
#ifdef _OPENMP
#pragma omp parallel for private(j) reduction(+ : sum)
#endif
    for (j = 1; j <= NbUnknowns; j++) {
      sum += InvMat[i][j] * RHS[j];  // new
    }

    Solution[i] = sum;
    printf("\b\b\b\b\b\b\b\b\b\b");
  }
  fflush(stdout);

  printf("\nsolution over for all elements ...\n");
  fflush(stdout);

  if (NbConstraints) {
    if (OptSystemChargeZero) {
      VSystemChargeZero = Solution[NbSystemChargeZero];
      printf("\nsolution over for system charge zero constraint ...\n");
      fflush(stdout);
    }

    if (NbFloatingConductors) {
      VFloatCon = Solution[NbFloatCon];
      printf(
          "\nsolution over for floating conductor charge zero constraint "
          "...\n");
      fflush(stdout);
    }
  }  // if NbConstraints

  // Update the element structure array and write the solution in a file
  char SolnFile[256];
  strcpy(SolnFile, BCOutDir);
  strcat(SolnFile, "/Soln.out");
  FILE *fSoln = fopen(SolnFile, "w");
  if (fSoln == NULL) {
    neBEMMessage("Solve - SolnFile");
    return -1;
  }
  fprintf(fSoln, "#EleNb\tX\tY\tZ\tSolution\tAssigned\tTotal\n");
  for (int ele = 1; ele <= NbElements; ++ele) {
    (EleArr + ele - 1)->Solution = Solution[ele];
    fprintf(fSoln, "%d\t%lg\t%lg\t%lg\t%.16lg\t%.16lg\t%.16lg\n", ele,
            (EleArr + ele - 1)->G.Origin.X, (EleArr + ele - 1)->G.Origin.Y,
            (EleArr + ele - 1)->G.Origin.Z, (EleArr + ele - 1)->Solution,
            (EleArr + ele - 1)->Assigned,
            ((EleArr + ele - 1)->Solution + (EleArr + ele - 1)->Assigned));
  }

  if (NbConstraints) {
    if (OptSystemChargeZero) {
      fprintf(fSoln, "#NbSystemChargeZero\tVSystemChargeZero\n");
      fprintf(fSoln, "# %d\t%lg\n", NbSystemChargeZero, VSystemChargeZero);
    }
    if (NbFloatingConductors) {
      fprintf(fSoln, "#NbFloatCon\tVFloatCon\n");
      fprintf(fSoln, "# %d\t%lg\n", NbFloatCon, VFloatCon);
    }
  }  // if NbConstraints

  fclose(fSoln);

  // Find primitive related charge densities
  char PrimSolnFile[256];
  strcpy(PrimSolnFile, BCOutDir);
  strcat(PrimSolnFile, "/PrimSoln.out");
  FILE *fPrimSoln = fopen(PrimSolnFile, "w");
  if (fPrimSoln == NULL) {
    neBEMMessage("Solve - PrimSolnFile");
    return -1;
  }
  fprintf(fPrimSoln,
          "#PrimNb\tEleBgn\tEleEnd\tX\tY\tZ\tAvChDen\tAvAsgndChDen\n");
  // OMPCheck - may be parallelized
  for (int prim = 1; prim <= NbPrimitives; ++prim) {
    double area = 0.0;  // need area of the primitive as well!
    AvChDen[prim] = 0.0;
    AvAsgndChDen[prim] = 0.0;

    for (int ele = ElementBgn[prim]; ele <= ElementEnd[prim]; ++ele) {
      area += (EleArr + ele - 1)->G.dA;
      AvChDen[prim] += (EleArr + ele - 1)->Solution * (EleArr + ele - 1)->G.dA;
      AvAsgndChDen[prim] +=
          (EleArr + ele - 1)->Assigned * (EleArr + ele - 1)->G.dA;
    }

    AvChDen[prim] /= area;
    AvAsgndChDen[prim] /= area;

    fprintf(fPrimSoln, "%d\t%d\t%d\t%lg\t%lg\t%lg\t%.16lg\t%.16g\n", prim,
            ElementBgn[prim], ElementEnd[prim], PrimOriginX[prim],
            PrimOriginY[prim], PrimOriginZ[prim], AvChDen[prim],
            AvAsgndChDen[prim]);
  }

  fclose(fPrimSoln);

  neBEMState = 9;

  // Validate the obtained solution, i.e., estimate the residue at the
  // collocation points. The validation is expected to be carried out by the
  // influence matrix already available in the memory. If, however, for some
  // reason (such as a repeat calculation), the influence coefficient matrix is
  // not available, it can be retrieved by carrying out the computation once
  // again, reading it from a formatted or an unformatted file.
  printf("OptValidateSolution: %d\n", OptValidateSolution);
  if (OptValidateSolution) {
    printf("Computing solution at the collocation points for comparison.\n");
    fflush(stdout);

    if (InfluenceMatrixFlag) {
      printf("Influence matrix in memory ...\n");
    }

    // Only when InfluenceMatrixFlag is false, the influence coefficient
    // reamins uncomputed. Then the question of re-computation or reading it
    // from a file (formatted / unformatted) arises.
    if (!InfluenceMatrixFlag) {
      if (TimeStep != 1) {
        // influence matrix to be computed only in the first time step
        printf("Influence matrix in memory ...\n");
      } else {
        printf("Influence matrix NOT in memory ...\n");

        if (OptRepeatLHMatrix) {
          printf("repeating influence coefficient matrix computation ...\n");

          int fstatus = LHMatrix();
          // assert(fstatus == 0);
          if (fstatus != 0) {
            neBEMMessage("Solve - LHMatrix in OptRepeatLHMatrix");
            return -1;
          }
        }

        if (OptStoreInflMatrix && OptFormattedFile) {
          printf(
              "reading influence coefficient matrix from formatted file...\n");

          char InflFile[256];
          strcpy(InflFile, MeshOutDir);
          strcat(InflFile, "/Infl.out");
          FILE *fInf = fopen(InflFile, "r");
          // assert(fInf != NULL);
          if (fInf == NULL) {
            neBEMMessage("Solve - InflFile in OptValidate.");
            return 1;
          }

          int chkNbEqns, chkNbUnknowns;
          fscanf(fInf, "%d %d\n", &chkNbEqns, &chkNbUnknowns);
          if ((chkNbEqns != NbEqns) || (chkNbUnknowns != NbUnknowns)) {
            neBEMMessage("Solve - matrix dimension do not match!");
            return (-1);
          }

          printf("Solve: Matrix dimensions: %d equations, %d unknowns\n",
                 NbEqns, NbUnknowns);

          Inf = dmatrix(1, NbEqns, 1, NbUnknowns);

          for (int elefld = 1; elefld <= NbEqns; ++elefld)
            for (int elesrc = 1; elesrc <= NbUnknowns; ++elesrc)
              fscanf(fInf, "%le\n", &Inf[elefld][elesrc]);
          fclose(fInf);
        }  // stored influence matrix and formatted file

        if (OptStoreInflMatrix && OptUnformattedFile) {
          neBEMMessage(
              "Solve - Binary read of Infl matrix not implemented yet.\n");
          return -1;

          RawInf = dmatrix(1, NbEqns, 1, NbUnknowns);

          printf(
              "reading influence coefficient matrix from unformatted file "
              "...\n");

          char InflFile[256];
          strcpy(InflFile, MeshOutDir);
          strcat(InflFile, "/RawInfl.out");
          printf("\nread from file %s\n", InflFile);
          fflush(stdout);
          FILE *fInf = fopen(InflFile, "rb");
          // assert(fInf != NULL);
          if (fInf == NULL) {
            neBEMMessage("Solve - RawInflFile in OptValidate");
            return -1;
          }
          printf("\nfInf: %p\n", (void *)fInf);
          int rfw = fread(RawInf, sizeof(double), NbEqns * NbUnknowns, fInf);
          fclose(fInf);
          printf("\nNb of items successfully read from raw mode in %s is %d\n",
                 InflFile, rfw);
          fflush(stdout);

          for (int unknown = 1; unknown <= NbUnknowns; ++unknown)
            for (int eqn = 1; eqn <= NbEqns; ++eqn)
              printf(
                  "Unknown:%d, Eqn:%d => diff Inf: %lg, RawInf: %lg is %lg\n",
                  unknown, eqn, Inf[unknown][eqn], RawInf[unknown][eqn],
                  fabs(Inf[unknown][eqn] - RawInf[unknown][eqn]));
        }  // if OptStoreInflMatrix and Unformatted file
      }    // else TimeStep != 1
    }      // if(!InfluenceMatrixFlag)

    // Used for all validations except where re-computation is forced which
    // is taken care of by else of this if
    if (Inf || RawInf) {
      double XChk;
      char Chkfile[256];
      strcpy(Chkfile, BCOutDir);
      strcat(Chkfile, "/XChk.out");
      FILE *fChk = fopen(Chkfile, "w");  // assert(fChk != NULL);
      if (fChk == NULL) {
        neBEMMessage("Solve - ChkFile");
        return -1;
      }
      fprintf(fChk, "#Row\tGivenPot\tError\n");

      int ElementOfMaxError = 1;
      double *Error, MaxError = 0.0;
      Error = dvector(1, NbEqns);
      int elesrc;
#ifdef _OPENMP
#pragma omp parallel for private(elesrc)
#endif
      for (int elefld = 1; elefld <= NbEqns; elefld++) {
        XChk = 0.0;
        for (elesrc = 1; elesrc <= NbUnknowns; elesrc++) {
          if ((!InfluenceMatrixFlag) && OptStoreInflMatrix &&
              OptUnformattedFile)
            XChk += RawInf[elefld][elesrc] * Solution[elesrc];
          else
            XChk += Inf[elefld][elesrc] * Solution[elesrc];
        }  // for elesrc
        Error[elefld] = fabs(RHS[elefld] - XChk);

        if (Error[elefld] > MaxError) {
          MaxError = Error[elefld];  // update maximum error related info
          ElementOfMaxError = elefld;
        }
      }  // for elefld
      for (int elefld = 1; elefld <= NbEqns; elefld++)
        fprintf(fChk, "%d\t%lg\t%lg\n", elefld, RHS[elefld], Error[elefld]);
      free_dvector(Error, 1, NbEqns);

      printf("Computed values at the collocation points for comparison.\n");
      printf("Error maximum on element %d and its magnitude is %lg.\n",
             ElementOfMaxError, MaxError);
      fflush(stdout);

      fprintf(fChk, "#Error maximum on element %d and its magnitude is %lg.\n",
              ElementOfMaxError, MaxError);
      fclose(fChk);

      if ((!InfluenceMatrixFlag) && OptRepeatLHMatrix)
        free_dmatrix(Inf, 1, NbEqns, 1, NbUnknowns);
      if ((!InfluenceMatrixFlag) && OptStoreInflMatrix && OptFormattedFile)
        free_dmatrix(Inf, 1, NbEqns, 1, NbUnknowns);
      if ((!InfluenceMatrixFlag) && OptStoreInflMatrix && OptUnformattedFile)
        free_dmatrix(RawInf, 1, NbEqns, 1, NbUnknowns);
    }  // if(Inf || RawInf)
    else {
      if (OptForceValidation) {
        neBEMMessage(
            "Solve - Infl matrix not available, re-computation forced.\n");

        int fstatus = LHMatrix();
        if (fstatus != 0) {
          neBEMMessage("Solve - LHMatrix in forced OptRepeatLHMatrix");
          return -1;
        }

        double XChk;
        char Chkfile[256];
        FILE *fChk;
        strcpy(Chkfile, BCOutDir);
        strcat(Chkfile, "/XChk.out");
        fChk = fopen(Chkfile, "w");
        if (fChk == NULL) {
          neBEMMessage("Solve - ChkFile in forced recomputation");
          return -1;
        }
        fprintf(fChk, "#Row\tGivenPot\tComputedPot\tDiff\n");

        int ElementOfMaxError = 1;
        double Error, MaxError = 0.0;
        for (int elefld = 1; elefld <= NbEqns; elefld++) {
          XChk = 0.0;
          for (int elesrc = 1; elesrc <= NbUnknowns; elesrc++) {
            XChk += Inf[elefld][elesrc] * Solution[elesrc];
          }
          Error = fabs(RHS[elefld] - XChk);
          fprintf(fChk, "%d\t%lg\t%lg\t%lg\n", elefld, RHS[elefld], XChk,
                  Error);

          if (Error > MaxError) {
            MaxError = Error;  // update maximum error related info
            ElementOfMaxError = elefld;
          }
        }

        printf("Computed values at the collocation points for comparison.\n");
        printf("Error maximum on element %d and its magnitude is %lg.\n",
               ElementOfMaxError, MaxError);
        fflush(stdout);

        fprintf(fChk,
                "#Error maximum on element %d and its magnitude is %lg.\n",
                ElementOfMaxError, MaxError);
        fclose(fChk);

        free_dmatrix(Inf, 1, NbEqns, 1, NbUnknowns);
      } else {  // this is not an error, though
        neBEMMessage("Solve - Infl matrix not available, no validation.\n");
      }
    }  // else (Inf || RawInf)
  }    // if(OptValidateSolution)

  /*
  Find out the error at relevant points. The points are chosen such that the
  mesh can be refined based on the errors obtained. For example, if the error
  at a particular point is found to be more than acceptable, an additional
  element (including necessary adjustment of the elements on the primitive
  in question) can easily be inserted to improve the accuracy of the solution.
  The validation is NOT expected to be carried out by the influence matrix
  already available in the memory, since the evaluation points no longer
  coincides with the collocation (qualocation?) points used to compute the
  solution. There is no question of repeating the computation, as while
  estimating the RESIDUE (Validation block immediately above this paragraph).
  The relevant points of rectangular elements are the centroid of four
  possible quadrants, whereas, for the triangular elements, these are the
  centroid of the rectangular and two barycentre of the two possible right
  triangles.
  In the following figures, X-s represent the collocation points, while x-s
  represent the points introduced to estimate error. For the triangle,
  there is a possibility of overlap of X and x of the rectangle.
           ----------------------------
           |             |            |
           |      x      |      x     |
           |             |            |
           --------------X-------------
           |             |            |
           |      x      |      x     |
           |             |            |
           ----------------------------
           \
           |\
           | \
           |x \
           ----\
           |   |\
           | x | \
           |   |x \
           --------\
  */

  OptEstimateError = 0;  // temporary measure
  printf("OptEstimateError: %d\n", OptEstimateError);
  if (OptEstimateError) {
    printf(
        "Properties at non-collocation points on element for estimating "
        "error.\n");
    fflush(stdout);

    double Err;
    int PrimitiveOfMaxError = 1;
    int ElementOfMaxError = 1;
    double MaxError = 0.0;
    double xerrMax = 0.0, yerrMax = 0.0, zerrMax = 0.0;

    char Errfile[256];
    strcpy(Errfile, BCOutDir);
    strcat(Errfile, "/Errors.out");
    FILE *fErr = fopen(Errfile, "w");
    if (fErr == NULL) {
      neBEMMessage("Solve - ErrFile");
      return -1;
    }
    fprintf(fErr,
            "#Prim\tEle\tGType\tIType\txerr\tyerr\tzerr\tGivenBC\tComputVal\tDi"
            "ff\n");

    // loop over all the elements - at present wire elements are NOT
    // considered!! We now also have ElementBgn and ElementEnd for each
    // primitive. So, we can loop over primtives and, within it, the necessary
    // elements can be considered. It turns out that going from primitive to
    // primitive can result into lesser computation. For example, E.Type for an
    // entire primitive is the same, and does not need to be evaluated from
    // element to element. conductor-dielectric: 1, conductor with known charge:
    // 2 conductor with floating potential: 3, dielectric-dielectric: 4
    // dielectric with known charge: 5
    // At present, error-accounting is maintained only for two types of
    // interfaces, C-D interfaces (specified voltage), D-D interfaces
    // (continuity of displacement current implied)
    // OMPCheck - may be parallelized
    for (int prim = 1; prim <= NbPrimitives; ++prim)
      for (int ele = ElementBgn[prim]; ele <= ElementEnd[prim]; ++ele) {
        double normdisp = 1.0e-8;  // displacement in the normal direction

        if ((prim == 91) && (ele == 337))  // debug switches are turned on
        {
          // DebugLevel = 1001;
          DebugLevel = 0;
        } else  // turn off switches
        {
          DebugLevel = 0;
        }

        // find relevant points for this element; compute PF and cross-check
        // with BC
        if ((EleArr + ele - 1)->G.Type == 2) continue;

        if ((EleArr + ele - 1)->G.Type == 3)  // triangle
        {  // 3 points need to be considered for estimating error (above fig)
          Point3D globalP;
          double Potential;
          Vector3D globalF, localF;
          double x0 = (EleArr + ele - 1)->G.Vertex[0].X;
          double y0 = (EleArr + ele - 1)->G.Vertex[0].Y;
          double z0 = (EleArr + ele - 1)->G.Vertex[0].Z;
          double x1 = (EleArr + ele - 1)->G.Vertex[1].X;
          double y1 = (EleArr + ele - 1)->G.Vertex[1].Y;
          double z1 = (EleArr + ele - 1)->G.Vertex[1].Z;
          double x2 = (EleArr + ele - 1)->G.Vertex[2].X;
          double y2 = (EleArr + ele - 1)->G.Vertex[2].Y;
          double z2 = (EleArr + ele - 1)->G.Vertex[2].Z;
          double xb = (x0 + x1 + x2) / 3.0;  // b for barycentre
          double yb = (y0 + y1 + y2) / 3.0;  // b for barycentre
          double zb = (z0 + z1 + z2) / 3.0;  // b for barycentre

          // check values at the barycentre
          globalP.X = xb;
          globalP.Y = yb;
          globalP.Z = zb;
          if (InterfaceType[prim] == 1) {
            PFAtPoint(&globalP, &Potential, &globalF);
            Err = ApplPot[prim] - Potential;
            fprintf(fErr, "%d\t%d\t%d\t%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                    prim, ele, 3, 1, xb, yb, zb, ApplPot[prim], Potential, Err);
            if (fabs(Err) > fabs(MaxError)) {
              MaxError = Err;
              PrimitiveOfMaxError = prim;
              ElementOfMaxError = ele;
              xerrMax = xb;
              yerrMax = yb;
              zerrMax = zb;
            }
          }
          if (InterfaceType[prim] ==
              4) {  // compute displacement currents in the two dielectrics
            double xplus = xb + (EleArr + ele - 1)->G.DC.XUnit.X * normdisp;
            xplus += (EleArr + ele - 1)->G.DC.YUnit.X * normdisp;
            xplus += (EleArr + ele - 1)->G.DC.ZUnit.X * normdisp;
            double yplus = yb + (EleArr + ele - 1)->G.DC.XUnit.Y * normdisp;
            yplus += (EleArr + ele - 1)->G.DC.YUnit.Y * normdisp;
            yplus += (EleArr + ele - 1)->G.DC.ZUnit.Y * normdisp;
            double zplus = zb + (EleArr + ele - 1)->G.DC.XUnit.Z * normdisp;
            zplus += (EleArr + ele - 1)->G.DC.YUnit.Z * normdisp;
            zplus += (EleArr + ele - 1)->G.DC.ZUnit.Z * normdisp;
            globalP.X = xplus;
            globalP.Y = yplus;
            globalP.Z = zplus;
            PFAtPoint(&globalP, &Potential, &globalF);
            localF  // Flux in the ECS
                = RotateVector3D(&globalF, &(EleArr + ele - 1)->G.DC,
                                 global2local);
            double value1 = -localF.Y;
            double xminus = xb - (EleArr + ele - 1)->G.DC.XUnit.X * normdisp;
            xminus -= (EleArr + ele - 1)->G.DC.YUnit.X * normdisp;
            xminus -= (EleArr + ele - 1)->G.DC.ZUnit.X * normdisp;
            double yminus = yb - (EleArr + ele - 1)->G.DC.XUnit.Y * normdisp;
            yminus -= (EleArr + ele - 1)->G.DC.YUnit.Y * normdisp;
            yminus -= (EleArr + ele - 1)->G.DC.ZUnit.Y * normdisp;
            double zminus = zb - (EleArr + ele - 1)->G.DC.XUnit.Z * normdisp;
            zminus -= (EleArr + ele - 1)->G.DC.YUnit.Z * normdisp;
            zminus -= (EleArr + ele - 1)->G.DC.ZUnit.Z * normdisp;
            globalP.X = xminus;
            globalP.Y = yminus;
            globalP.Z = zminus;
            PFAtPoint(&globalP, &Potential, &globalF);
            localF  // Flux in the ECS
                = RotateVector3D(&globalF, &(EleArr + ele - 1)->G.DC,
                                 global2local);
            double value2 = -localF.Y;
            double epsratio = (Epsilon2[prim] / Epsilon1[prim]);
            Err = epsratio - (value1 / value2);
            fprintf(fErr, "%d\t%d\t%d\t%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                    prim, ele, 3, 4, xb, yb, zb, epsratio, (value1 / value2),
                    Err);
            if (fabs(Err) > fabs(MaxError)) {
              MaxError = Err;
              PrimitiveOfMaxError = prim;
              ElementOfMaxError = ele;
              xerrMax = xb;
              yerrMax = yb;
              zerrMax = zb;
            }
          }

          // values at the centroid of the rectangular part
          double xerr =
              0.25 * (x0 + 0.5 * (x0 + x1) + 0.5 * (x1 + x2) + 0.5 * (x0 + x2));
          double yerr =
              0.25 * (y0 + 0.5 * (y0 + y1) + 0.5 * (y1 + y2) + 0.5 * (y0 + y2));
          double zerr =
              0.25 * (z0 + 0.5 * (z0 + z1) + 0.5 * (z1 + z2) + 0.5 * (z0 + z2));
          globalP.X = xerr;
          globalP.Y = yerr;
          globalP.Z = zerr;
          // At present, error-accounting is maintained only for two types of
          // interfaces, C-D interfaces (specified voltage), D-D interfaces
          // (continuity of displacement current implied)
          if (InterfaceType[prim] == 1) {
            PFAtPoint(&globalP, &Potential, &globalF);
            Err = ApplPot[prim] - Potential;
            fprintf(fErr, "%d\t%d\t%d\t%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                    prim, ele, 3, 1, xerr, yerr, zerr, ApplPot[prim], Potential,
                    Err);
            if (fabs(Err) > fabs(MaxError)) {
              MaxError = Err;
              PrimitiveOfMaxError = prim;
              ElementOfMaxError = ele;
              xerrMax = xerr;
              yerrMax = yerr;
              zerrMax = zerr;
            }
          }
          if (InterfaceType[prim] ==
              4) {  // compute displacement currents in the two dielectrics
            double xplus = xerr + (EleArr + ele - 1)->G.DC.XUnit.X * normdisp;
            xplus += (EleArr + ele - 1)->G.DC.YUnit.X * normdisp;
            xplus += (EleArr + ele - 1)->G.DC.ZUnit.X * normdisp;
            double yplus = yerr + (EleArr + ele - 1)->G.DC.XUnit.Y * normdisp;
            yplus += (EleArr + ele - 1)->G.DC.YUnit.Y * normdisp;
            yplus += (EleArr + ele - 1)->G.DC.ZUnit.Y * normdisp;
            double zplus = zerr + (EleArr + ele - 1)->G.DC.XUnit.Z * normdisp;
            zplus += (EleArr + ele - 1)->G.DC.YUnit.Z * normdisp;
            zplus += (EleArr + ele - 1)->G.DC.ZUnit.Z * normdisp;
            globalP.X = xplus;
            globalP.Y = yplus;
            globalP.Z = zplus;
            PFAtPoint(&globalP, &Potential, &globalF);
            localF  // Flux in the ECS
                = RotateVector3D(&globalF, &(EleArr + ele - 1)->G.DC,
                                 global2local);
            double value1 = -localF.Y;
            double xminus = xerr - (EleArr + ele - 1)->G.DC.XUnit.X * normdisp;
            xminus -= (EleArr + ele - 1)->G.DC.YUnit.X * normdisp;
            xminus -= (EleArr + ele - 1)->G.DC.ZUnit.X * normdisp;
            double yminus = yerr - (EleArr + ele - 1)->G.DC.XUnit.Y * normdisp;
            yminus -= (EleArr + ele - 1)->G.DC.YUnit.Y * normdisp;
            yminus -= (EleArr + ele - 1)->G.DC.ZUnit.Y * normdisp;
            double zminus = zerr - (EleArr + ele - 1)->G.DC.XUnit.Z * normdisp;
            zminus -= (EleArr + ele - 1)->G.DC.YUnit.Z * normdisp;
            zminus -= (EleArr + ele - 1)->G.DC.ZUnit.Z * normdisp;
            globalP.X = xminus;
            globalP.Y = yminus;
            globalP.Z = zminus;
            PFAtPoint(&globalP, &Potential, &globalF);
            localF  // Flux in the ECS
                = RotateVector3D(&globalF, &(EleArr + ele - 1)->G.DC,
                                 global2local);
            double value2 = -localF.Y;
            double epsratio = (Epsilon2[prim] / Epsilon1[prim]);
            Err = epsratio - (value1 / value2);
            fprintf(fErr, "%d\t%d\t%d\t%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                    prim, ele, 3, 4, xerr, yerr, zerr, epsratio,
                    (value1 / value2), Err);
            if (fabs(Err) > fabs(MaxError)) {
              MaxError = Err;
              PrimitiveOfMaxError = prim;
              ElementOfMaxError = ele;
              xerrMax = xerr;
              yerrMax = yerr;
              zerrMax = zerr;
            }
          }

          // barycentre of the triangular part towards +ve X-axis of ECS
          xerr = (0.5 * (x0 + x1) + 0.5 * (x1 + x2) + x1) / 3.0;
          yerr = (0.5 * (y0 + y1) + 0.5 * (y1 + y2) + y1) / 3.0;
          zerr = (0.5 * (z0 + z1) + 0.5 * (z1 + z2) + z1) / 3.0;
          globalP.X = xerr;
          globalP.Y = yerr;
          globalP.Z = zerr;
          // At present, error-accounting is maintained only for two types of
          // interfaces, C-D interfaces (specified voltage), D-D interfaces
          // (continuity of displacement current implied)
          if (InterfaceType[prim] == 1) {
            PFAtPoint(&globalP, &Potential, &globalF);
            Err = ApplPot[prim] - Potential;
            fprintf(fErr, "%d\t%d\t%d\t%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                    prim, ele, 3, 1, xerr, yerr, zerr, ApplPot[prim], Potential,
                    Err);
            if (fabs(Err) > fabs(MaxError)) {
              MaxError = Err;
              PrimitiveOfMaxError = prim;
              ElementOfMaxError = ele;
              xerrMax = xerr;
              yerrMax = yerr;
              zerrMax = zerr;
            }
          }
          if (InterfaceType[prim] ==
              4) {  // compute displacement currents in the two dielectrics
            double xplus = xerr + (EleArr + ele - 1)->G.DC.XUnit.X * normdisp;
            xplus += (EleArr + ele - 1)->G.DC.YUnit.X * normdisp;
            xplus += (EleArr + ele - 1)->G.DC.ZUnit.X * normdisp;
            double yplus = yerr + (EleArr + ele - 1)->G.DC.XUnit.Y * normdisp;
            yplus += (EleArr + ele - 1)->G.DC.YUnit.Y * normdisp;
            yplus += (EleArr + ele - 1)->G.DC.ZUnit.Y * normdisp;
            double zplus = zerr + (EleArr + ele - 1)->G.DC.XUnit.Z * normdisp;
            zplus += (EleArr + ele - 1)->G.DC.YUnit.Z * normdisp;
            zplus += (EleArr + ele - 1)->G.DC.ZUnit.Z * normdisp;
            globalP.X = xplus;
            globalP.Y = yplus;
            globalP.Z = zplus;
            PFAtPoint(&globalP, &Potential, &globalF);
            localF  // Flux in the ECS
                = RotateVector3D(&globalF, &(EleArr + ele - 1)->G.DC,
                                 global2local);
            double value1 = -localF.Y;
            double xminus = xerr - (EleArr + ele - 1)->G.DC.XUnit.X * normdisp;
            xminus -= (EleArr + ele - 1)->G.DC.YUnit.X * normdisp;
            xminus -= (EleArr + ele - 1)->G.DC.ZUnit.X * normdisp;
            double yminus = yerr - (EleArr + ele - 1)->G.DC.XUnit.Y * normdisp;
            yminus -= (EleArr + ele - 1)->G.DC.YUnit.Y * normdisp;
            yminus -= (EleArr + ele - 1)->G.DC.ZUnit.Y * normdisp;
            double zminus = zerr - (EleArr + ele - 1)->G.DC.XUnit.Z * normdisp;
            zminus -= (EleArr + ele - 1)->G.DC.YUnit.Z * normdisp;
            zminus -= (EleArr + ele - 1)->G.DC.ZUnit.Z * normdisp;
            globalP.X = xminus;
            globalP.Y = yminus;
            globalP.Z = zminus;
            PFAtPoint(&globalP, &Potential, &globalF);
            localF  // Flux in the ECS
                = RotateVector3D(&globalF, &(EleArr + ele - 1)->G.DC,
                                 global2local);
            double value2 = -localF.Y;
            double epsratio = (Epsilon2[prim] / Epsilon1[prim]);
            Err = epsratio - (value1 / value2);
            fprintf(fErr, "%d\t%d\t%d\t%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                    prim, ele, 3, 4, xerr, yerr, zerr, epsratio,
                    (value1 / value2), Err);
            if (fabs(Err) > fabs(MaxError)) {
              MaxError = Err;
              PrimitiveOfMaxError = prim;
              ElementOfMaxError = ele;
              xerrMax = xerr;
              yerrMax = yerr;
              zerrMax = zerr;
            }
          }

          // barycentre of the triangular part towards +ve Z-axis of ECS
          xerr = (0.5 * (x0 + x2) + 0.5 * (x1 + x2) + x2) / 3.0;
          yerr = (0.5 * (y0 + y2) + 0.5 * (y1 + y2) + y2) / 3.0;
          zerr = (0.5 * (z0 + z2) + 0.5 * (z1 + z2) + z2) / 3.0;
          globalP.X = xerr;
          globalP.Y = yerr;
          globalP.Z = zerr;
          // At present, error-accounting is maintained only for two types of
          // interfaces, C-D interfaces (specified voltage), D-D interfaces
          // (continuity of displacement current implied)
          if (InterfaceType[prim] == 1) {
            PFAtPoint(&globalP, &Potential, &globalF);
            Err = ApplPot[prim] - Potential;
            fprintf(fErr, "%d\t%d\t%d\t%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                    prim, ele, 3, 1, xerr, yerr, zerr, ApplPot[prim], Potential,
                    Err);
            if (fabs(Err) > fabs(MaxError)) {
              MaxError = Err;
              PrimitiveOfMaxError = prim;
              ElementOfMaxError = ele;
              xerrMax = xerr;
              yerrMax = yerr;
              zerrMax = zerr;
            }
          }
          if (InterfaceType[prim] ==
              4) {  // compute displacement currents in the two dielectrics
            double xplus = xerr + (EleArr + ele - 1)->G.DC.XUnit.X * normdisp;
            xplus += (EleArr + ele - 1)->G.DC.YUnit.X * normdisp;
            xplus += (EleArr + ele - 1)->G.DC.ZUnit.X * normdisp;
            double yplus = yerr + (EleArr + ele - 1)->G.DC.XUnit.Y * normdisp;
            yplus += (EleArr + ele - 1)->G.DC.YUnit.Y * normdisp;
            yplus += (EleArr + ele - 1)->G.DC.ZUnit.Y * normdisp;
            double zplus = zerr + (EleArr + ele - 1)->G.DC.XUnit.Z * normdisp;
            zplus += (EleArr + ele - 1)->G.DC.YUnit.Z * normdisp;
            zplus += (EleArr + ele - 1)->G.DC.ZUnit.Z * normdisp;
            globalP.X = xplus;
            globalP.Y = yplus;
            globalP.Z = zplus;
            PFAtPoint(&globalP, &Potential, &globalF);
            localF  // Flux in the ECS
                = RotateVector3D(&globalF, &(EleArr + ele - 1)->G.DC,
                                 global2local);
            double value1 = -localF.Y;
            double xminus = xerr - (EleArr + ele - 1)->G.DC.XUnit.X * normdisp;
            xminus -= (EleArr + ele - 1)->G.DC.YUnit.X * normdisp;
            xminus -= (EleArr + ele - 1)->G.DC.ZUnit.X * normdisp;
            double yminus = yerr - (EleArr + ele - 1)->G.DC.XUnit.Y * normdisp;
            yminus -= (EleArr + ele - 1)->G.DC.YUnit.Y * normdisp;
            yminus -= (EleArr + ele - 1)->G.DC.ZUnit.Y * normdisp;
            double zminus = zerr - (EleArr + ele - 1)->G.DC.XUnit.Z * normdisp;
            zminus -= (EleArr + ele - 1)->G.DC.YUnit.Z * normdisp;
            zminus -= (EleArr + ele - 1)->G.DC.ZUnit.Z * normdisp;
            globalP.X = xminus;
            globalP.Y = yminus;
            globalP.Z = zminus;
            PFAtPoint(&globalP, &Potential, &globalF);
            localF  // Flux in the ECS
                = RotateVector3D(&globalF, &(EleArr + ele - 1)->G.DC,
                                 global2local);
            double value2 = -localF.Y;
            double epsratio = (Epsilon2[prim] / Epsilon1[prim]);
            Err = epsratio - (value1 / value2);
            fprintf(fErr, "%d\t%d\t%d\t%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                    prim, ele, 4, 4, xerr, yerr, zerr, epsratio,
                    (value1 / value2), Err);
            if (fabs(Err) > fabs(MaxError)) {
              MaxError = Err;
              PrimitiveOfMaxError = prim;
              ElementOfMaxError = ele;
              xerrMax = xerr;
              yerrMax = yerr;
              zerrMax = zerr;
            }
          }
        }  // if triangle

        if ((EleArr + ele - 1)->G.Type == 4)  // rectangle
        {  // 4 points need to be considered for estimating error (above fig)
          Point3D globalP;
          double Potential;
          Vector3D globalF, localF;

          double x0 = (EleArr + ele - 1)->G.Vertex[0].X;
          double y0 = (EleArr + ele - 1)->G.Vertex[0].Y;
          double z0 = (EleArr + ele - 1)->G.Vertex[0].Z;
          double x1 = (EleArr + ele - 1)->G.Vertex[1].X;
          double y1 = (EleArr + ele - 1)->G.Vertex[1].Y;
          double z1 = (EleArr + ele - 1)->G.Vertex[1].Z;
          double x2 = (EleArr + ele - 1)->G.Vertex[2].X;
          double y2 = (EleArr + ele - 1)->G.Vertex[2].Y;
          double z2 = (EleArr + ele - 1)->G.Vertex[2].Z;
          double x3 = (EleArr + ele - 1)->G.Vertex[3].X;
          double y3 = (EleArr + ele - 1)->G.Vertex[3].Y;
          double z3 = (EleArr + ele - 1)->G.Vertex[3].Z;
          double xo = 0.25 * (x0 + x1 + x2 + x3);  // o for origin
          double yo = 0.25 * (y0 + y1 + y2 + y3);  // o for origin
          double zo = 0.25 * (z0 + z1 + z2 + z3);  // o for origin

          // check values at the centroid of the element
          globalP.X = xo;
          globalP.Y = yo;
          globalP.Z = zo;
          if (InterfaceType[prim] == 1) {
            PFAtPoint(&globalP, &Potential, &globalF);
            Err = ApplPot[prim] - Potential;
            fprintf(fErr, "%d\t%d\t%d\t%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                    prim, ele, 4, 1, xo, yo, zo, ApplPot[prim], Potential, Err);
            if (fabs(Err) > fabs(MaxError)) {
              MaxError = Err;
              PrimitiveOfMaxError = prim;
              ElementOfMaxError = ele;
              xerrMax = xo;
              yerrMax = yo;
              zerrMax = zo;
            }
          }
          if (InterfaceType[prim] == 4) {
            // compute displacement currents in the two dielectrics
            double xplus = xo + (EleArr + ele - 1)->G.DC.XUnit.X * normdisp;
            xplus += (EleArr + ele - 1)->G.DC.YUnit.X * normdisp;
            xplus += (EleArr + ele - 1)->G.DC.ZUnit.X * normdisp;
            double yplus = yo + (EleArr + ele - 1)->G.DC.XUnit.Y * normdisp;
            yplus += (EleArr + ele - 1)->G.DC.YUnit.Y * normdisp;
            yplus += (EleArr + ele - 1)->G.DC.ZUnit.Y * normdisp;
            double zplus = zo + (EleArr + ele - 1)->G.DC.XUnit.Z * normdisp;
            zplus += (EleArr + ele - 1)->G.DC.YUnit.Z * normdisp;
            zplus += (EleArr + ele - 1)->G.DC.ZUnit.Z * normdisp;
            globalP.X = xplus;
            globalP.Y = yplus;
            globalP.Z = zplus;
            PFAtPoint(&globalP, &Potential, &globalF);
            localF  // Flux in the ECS
                = RotateVector3D(&globalF, &(EleArr + ele - 1)->G.DC,
                                 global2local);
            double dispfld1 = Epsilon1[prim] * localF.Y;
            double xminus = xo - (EleArr + ele - 1)->G.DC.XUnit.X * normdisp;
            xminus -= (EleArr + ele - 1)->G.DC.YUnit.X * normdisp;
            xminus -= (EleArr + ele - 1)->G.DC.ZUnit.X * normdisp;
            double yminus = yo - (EleArr + ele - 1)->G.DC.XUnit.Y * normdisp;
            yminus -= (EleArr + ele - 1)->G.DC.YUnit.Y * normdisp;
            yminus -= (EleArr + ele - 1)->G.DC.ZUnit.Y * normdisp;
            double zminus = zo - (EleArr + ele - 1)->G.DC.XUnit.Z * normdisp;
            zminus -= (EleArr + ele - 1)->G.DC.YUnit.Z * normdisp;
            zminus -= (EleArr + ele - 1)->G.DC.ZUnit.Z * normdisp;
            globalP.X = xminus;
            globalP.Y = yminus;
            globalP.Z = zminus;
            PFAtPoint(&globalP, &Potential, &globalF);
            localF  // Flux in the ECS
                = RotateVector3D(&globalF, &(EleArr + ele - 1)->G.DC,
                                 global2local);
            double dispfld2 = Epsilon2[prim] * localF.Y;
            globalP.X = xo;
            globalP.Y = yo;
            globalP.Z = zo;
            PFAtPoint(&globalP, &Potential, &globalF);
            localF  // Flux in the ECS
                = RotateVector3D(&globalF, &(EleArr + ele - 1)->G.DC,
                                 global2local);
            double dispfldo = Epsilon1[prim] * localF.Y;
            Err = (dispfld2 - dispfld1) /
                  dispfldo;  // - (&(EleArr+ele-1)->Assigned);
            fprintf(fErr, "%d\t%d\t%d\t%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                    prim, ele, 4, 4, xo, yo, zo, dispfld1, dispfld2, Err);
            if ((prim == 91) && (ele == 337)) {
              // debug switches are turned on
              double TotalInfl = 0.0;
              for (int elesrc = 1; elesrc <= NbElements; ++elesrc) {
                TotalInfl += Inf[ele][elesrc] * (EleArr + elesrc - 1)->Solution;
              }
              printf("TotalInfl: %lg\n", TotalInfl);
              printf(
                  "%d\t%d\t%d\t%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%"
                  "lg\n",
                  prim, ele, 4, 4, xo, yo, zo, Epsilon1[prim], Epsilon2[prim],
                  dispfld1, dispfld2, dispfldo, Err);
            }
            if (fabs(Err) > fabs(MaxError)) {
              MaxError = Err;
              PrimitiveOfMaxError = prim;
              ElementOfMaxError = ele;
              xerrMax = xo;
              yerrMax = yo;
              zerrMax = zo;
            }
          }

          // centroid of bottom-left rectangle
          double xerr = 0.25 * (x0 + 0.5 * (x1 + x0) + xo + 0.5 * (x0 + x3));
          double yerr = 0.25 * (y0 + 0.5 * (y1 + y0) + yo + 0.5 * (y0 + y3));
          double zerr = 0.25 * (z0 + 0.5 * (z1 + z0) + zo + 0.5 * (z0 + z3));
          globalP.X = xerr;
          globalP.Y = yerr;
          globalP.Z = zerr;
          if (InterfaceType[prim] == 1) {
            PFAtPoint(&globalP, &Potential, &globalF);
            Err = ApplPot[prim] - Potential;
            fprintf(fErr, "%d\t%d\t%d\t%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                    prim, ele, 4, 1, xerr, yerr, zerr, ApplPot[prim], Potential,
                    Err);
            if (fabs(Err) > fabs(MaxError)) {
              MaxError = Err;
              PrimitiveOfMaxError = prim;
              ElementOfMaxError = ele;
              xerrMax = xerr;
              yerrMax = yerr;
              zerrMax = zerr;
            }
          }
          if (InterfaceType[prim] ==
              4) {  // compute displacement currents in the two dielectrics
            double xplus = xerr + (EleArr + ele - 1)->G.DC.XUnit.X * normdisp;
            xplus += (EleArr + ele - 1)->G.DC.YUnit.X * normdisp;
            xplus += (EleArr + ele - 1)->G.DC.ZUnit.X * normdisp;
            double yplus = yerr + (EleArr + ele - 1)->G.DC.XUnit.Y * normdisp;
            yplus += (EleArr + ele - 1)->G.DC.YUnit.Y * normdisp;
            yplus += (EleArr + ele - 1)->G.DC.ZUnit.Y * normdisp;
            double zplus = zerr + (EleArr + ele - 1)->G.DC.XUnit.Z * normdisp;
            zplus += (EleArr + ele - 1)->G.DC.YUnit.Z * normdisp;
            zplus += (EleArr + ele - 1)->G.DC.ZUnit.Z * normdisp;
            globalP.X = xplus;
            globalP.Y = yplus;
            globalP.Z = zplus;
            PFAtPoint(&globalP, &Potential, &globalF);
            localF  // Flux in the ECS
                = RotateVector3D(&globalF, &(EleArr + ele - 1)->G.DC,
                                 global2local);
            double value1 = -localF.Y;
            double xminus = xerr - (EleArr + ele - 1)->G.DC.XUnit.X * normdisp;
            xminus -= (EleArr + ele - 1)->G.DC.YUnit.X * normdisp;
            xminus -= (EleArr + ele - 1)->G.DC.ZUnit.X * normdisp;
            double yminus = yerr - (EleArr + ele - 1)->G.DC.XUnit.Y * normdisp;
            yminus -= (EleArr + ele - 1)->G.DC.YUnit.Y * normdisp;
            yminus -= (EleArr + ele - 1)->G.DC.ZUnit.Y * normdisp;
            double zminus = zerr - (EleArr + ele - 1)->G.DC.XUnit.Z * normdisp;
            zminus -= (EleArr + ele - 1)->G.DC.YUnit.Z * normdisp;
            zminus -= (EleArr + ele - 1)->G.DC.ZUnit.Z * normdisp;
            globalP.X = xminus;
            globalP.Y = yminus;
            globalP.Z = zminus;
            PFAtPoint(&globalP, &Potential, &globalF);
            localF  // Flux in the ECS
                = RotateVector3D(&globalF, &(EleArr + ele - 1)->G.DC,
                                 global2local);
            double value2 = -localF.Y;
            double epsratio = (Epsilon2[prim] / Epsilon1[prim]);
            Err = epsratio - (value1 / value2);
            fprintf(fErr, "%d\t%d\t%d\t%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                    prim, ele, 4, 4, xerr, yerr, zerr, epsratio,
                    (value1 / value2), Err);
            if (fabs(Err) > fabs(MaxError)) {
              MaxError = Err;
              PrimitiveOfMaxError = prim;
              ElementOfMaxError = ele;
              xerrMax = xerr;
              yerrMax = yerr;
              zerrMax = zerr;
            }
          }

          // centroid of bottom-right rectangle
          xerr = 0.25 * (0.5 * (x1 + x0) + x1 + 0.5 * (x1 + x2) + xo);
          yerr = 0.25 * (0.5 * (y1 + y0) + y1 + 0.5 * (y1 + y2) + yo);
          zerr = 0.25 * (0.5 * (z1 + z0) + z1 + 0.5 * (z1 + z2) + zo);
          globalP.X = xerr;
          globalP.Y = yerr;
          globalP.Z = zerr;
          // At present, error-accounting is maintained only for two types of
          // interfaces, C-D interfaces (specified voltage), D-D interfaces
          // (continuity of displacement current implied)
          if (InterfaceType[prim] == 1) {
            PFAtPoint(&globalP, &Potential, &globalF);
            Err = ApplPot[prim] - Potential;
            fprintf(fErr, "%d\t%d\t%d\t%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                    prim, ele, 4, 1, xerr, yerr, zerr, ApplPot[prim], Potential,
                    Err);
            if (fabs(Err) > fabs(MaxError)) {
              MaxError = Err;
              PrimitiveOfMaxError = prim;
              ElementOfMaxError = ele;
              xerrMax = xerr;
              yerrMax = yerr;
              zerrMax = zerr;
            }
          }
          if (InterfaceType[prim] ==
              4) {  // compute displacement currents in the two dielectrics
            double xplus = xerr + (EleArr + ele - 1)->G.DC.XUnit.X * normdisp;
            xplus += (EleArr + ele - 1)->G.DC.YUnit.X * normdisp;
            xplus += (EleArr + ele - 1)->G.DC.ZUnit.X * normdisp;
            double yplus = yerr + (EleArr + ele - 1)->G.DC.XUnit.Y * normdisp;
            yplus += (EleArr + ele - 1)->G.DC.YUnit.Y * normdisp;
            yplus += (EleArr + ele - 1)->G.DC.ZUnit.Y * normdisp;
            double zplus = zerr + (EleArr + ele - 1)->G.DC.XUnit.Z * normdisp;
            zplus += (EleArr + ele - 1)->G.DC.YUnit.Z * normdisp;
            zplus += (EleArr + ele - 1)->G.DC.ZUnit.Z * normdisp;
            globalP.X = xplus;
            globalP.Y = yplus;
            globalP.Z = zplus;
            PFAtPoint(&globalP, &Potential, &globalF);
            localF  // Flux in the ECS
                = RotateVector3D(&globalF, &(EleArr + ele - 1)->G.DC,
                                 global2local);
            double value1 = -localF.Y;
            double xminus = xerr - (EleArr + ele - 1)->G.DC.XUnit.X * normdisp;
            xminus -= (EleArr + ele - 1)->G.DC.YUnit.X * normdisp;
            xminus -= (EleArr + ele - 1)->G.DC.ZUnit.X * normdisp;
            double yminus = yerr - (EleArr + ele - 1)->G.DC.XUnit.Y * normdisp;
            yminus -= (EleArr + ele - 1)->G.DC.YUnit.Y * normdisp;
            yminus -= (EleArr + ele - 1)->G.DC.ZUnit.Y * normdisp;
            double zminus = zerr - (EleArr + ele - 1)->G.DC.XUnit.Z * normdisp;
            zminus -= (EleArr + ele - 1)->G.DC.YUnit.Z * normdisp;
            zminus -= (EleArr + ele - 1)->G.DC.ZUnit.Z * normdisp;
            globalP.X = xminus;
            globalP.Y = yminus;
            globalP.Z = zminus;
            PFAtPoint(&globalP, &Potential, &globalF);
            localF  // Flux in the ECS
                = RotateVector3D(&globalF, &(EleArr + ele - 1)->G.DC,
                                 global2local);
            double value2 = -localF.Y;
            double epsratio = (Epsilon2[prim] / Epsilon1[prim]);
            Err = epsratio - (value1 / value2);
            fprintf(fErr, "%d\t%d\t%d\t%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                    prim, ele, 4, 4, xerr, yerr, zerr, epsratio,
                    (value1 / value2), Err);
            if (fabs(Err) > fabs(MaxError)) {
              MaxError = Err;
              PrimitiveOfMaxError = prim;
              ElementOfMaxError = ele;
              xerrMax = xerr;
              yerrMax = yerr;
              zerrMax = zerr;
            }
          }

          // centroid of top-right rectangle
          xerr = 0.25 * (xo + 0.5 * (x1 + x2) + x2 + 0.5 * (x3 + x2));
          yerr = 0.25 * (yo + 0.5 * (y1 + y2) + y2 + 0.5 * (y3 + y2));
          zerr = 0.25 * (zo + 0.5 * (z1 + z2) + z2 + 0.5 * (z3 + z2));
          globalP.X = xerr;
          globalP.Y = yerr;
          globalP.Z = zerr;
          // At present, error-accounting is maintained only for two types of
          // interfaces, C-D interfaces (specified voltage), D-D interfaces
          // (continuity of displacement current implied)
          if (InterfaceType[prim] == 1) {
            PFAtPoint(&globalP, &Potential, &globalF);
            Err = ApplPot[prim] - Potential;
            fprintf(fErr, "%d\t%d\t%d\t%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                    prim, ele, 4, 1, xerr, yerr, zerr, ApplPot[prim], Potential,
                    Err);
            if (fabs(Err) > fabs(MaxError)) {
              MaxError = Err;
              PrimitiveOfMaxError = prim;
              ElementOfMaxError = ele;
              xerrMax = xerr;
              yerrMax = yerr;
              zerrMax = zerr;
            }
          }
          if (InterfaceType[prim] ==
              4) {  // compute displacement currents in the two dielectrics
            double xplus = xerr + (EleArr + ele - 1)->G.DC.XUnit.X * normdisp;
            xplus += (EleArr + ele - 1)->G.DC.YUnit.X * normdisp;
            xplus += (EleArr + ele - 1)->G.DC.ZUnit.X * normdisp;
            double yplus = yerr + (EleArr + ele - 1)->G.DC.XUnit.Y * normdisp;
            yplus += (EleArr + ele - 1)->G.DC.YUnit.Y * normdisp;
            yplus += (EleArr + ele - 1)->G.DC.ZUnit.Y * normdisp;
            double zplus = zerr + (EleArr + ele - 1)->G.DC.XUnit.Z * normdisp;
            zplus += (EleArr + ele - 1)->G.DC.YUnit.Z * normdisp;
            zplus += (EleArr + ele - 1)->G.DC.ZUnit.Z * normdisp;
            globalP.X = xplus;
            globalP.Y = yplus;
            globalP.Z = zplus;
            PFAtPoint(&globalP, &Potential, &globalF);
            localF  // Flux in the ECS
                = RotateVector3D(&globalF, &(EleArr + ele - 1)->G.DC,
                                 global2local);
            double value1 = -localF.Y;
            double xminus = xerr - (EleArr + ele - 1)->G.DC.XUnit.X * normdisp;
            xminus -= (EleArr + ele - 1)->G.DC.YUnit.X * normdisp;
            xminus -= (EleArr + ele - 1)->G.DC.ZUnit.X * normdisp;
            double yminus = yerr - (EleArr + ele - 1)->G.DC.XUnit.Y * normdisp;
            yminus -= (EleArr + ele - 1)->G.DC.YUnit.Y * normdisp;
            yminus -= (EleArr + ele - 1)->G.DC.ZUnit.Y * normdisp;
            double zminus = zerr - (EleArr + ele - 1)->G.DC.XUnit.Z * normdisp;
            zminus -= (EleArr + ele - 1)->G.DC.YUnit.Z * normdisp;
            zminus -= (EleArr + ele - 1)->G.DC.ZUnit.Z * normdisp;
            globalP.X = xminus;
            globalP.Y = yminus;
            globalP.Z = zminus;
            PFAtPoint(&globalP, &Potential, &globalF);
            localF  // Flux in the ECS
                = RotateVector3D(&globalF, &(EleArr + ele - 1)->G.DC,
                                 global2local);
            double value2 = -localF.Y;
            double epsratio = (Epsilon2[prim] / Epsilon1[prim]);
            Err = epsratio - (value1 / value2);
            fprintf(fErr, "%d\t%d\t%d\t%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                    prim, ele, 4, 4, xerr, yerr, zerr, epsratio,
                    (value1 / value2), Err);
            if (fabs(Err) > fabs(MaxError)) {
              MaxError = Err;
              PrimitiveOfMaxError = prim;
              ElementOfMaxError = ele;
              xerrMax = xerr;
              yerrMax = yerr;
              zerrMax = zerr;
            }
          }

          // centroid of top-left rectangle
          xerr = 0.25 * (0.5 * (x0 + x3) + xo + 0.5 * (x3 + x2) + x3);
          yerr = 0.25 * (0.5 * (y0 + y3) + yo + 0.5 * (y3 + y2) + y3);
          zerr = 0.25 * (0.5 * (z0 + z3) + zo + 0.5 * (z3 + z2) + z3);
          globalP.X = xerr;
          globalP.Y = yerr;
          globalP.Z = zerr;
          // At present, error-accounting is maintained only for two types of
          // interfaces, C-D interfaces (specified voltage), D-D interfaces
          // (continuity of displacement current implied)
          if (InterfaceType[prim] == 1) {
            PFAtPoint(&globalP, &Potential, &globalF);
            Err = ApplPot[prim] - Potential;
            fprintf(fErr, "%d\t%d\t%d\t%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                    prim, ele, 4, 1, xerr, yerr, zerr, ApplPot[prim], Potential,
                    Err);
            if (fabs(Err) > fabs(MaxError)) {
              MaxError = Err;
              PrimitiveOfMaxError = prim;
              ElementOfMaxError = ele;
              xerrMax = xerr;
              yerrMax = yerr;
              zerrMax = zerr;
            }
          }
          if (InterfaceType[prim] ==
              4) {  // compute displacement currents in the two dielectrics
            double xplus = xerr + (EleArr + ele - 1)->G.DC.XUnit.X * normdisp;
            xplus += (EleArr + ele - 1)->G.DC.YUnit.X * normdisp;
            xplus += (EleArr + ele - 1)->G.DC.ZUnit.X * normdisp;
            double yplus = yerr + (EleArr + ele - 1)->G.DC.XUnit.Y * normdisp;
            yplus += (EleArr + ele - 1)->G.DC.YUnit.Y * normdisp;
            yplus += (EleArr + ele - 1)->G.DC.ZUnit.Y * normdisp;
            double zplus = zerr + (EleArr + ele - 1)->G.DC.XUnit.Z * normdisp;
            zplus += (EleArr + ele - 1)->G.DC.YUnit.Z * normdisp;
            zplus += (EleArr + ele - 1)->G.DC.ZUnit.Z * normdisp;
            globalP.X = xplus;
            globalP.Y = yplus;
            globalP.Z = zplus;
            PFAtPoint(&globalP, &Potential, &globalF);
            localF  // Flux in the ECS
                = RotateVector3D(&globalF, &(EleArr + ele - 1)->G.DC,
                                 global2local);
            double value1 = -localF.Y;
            double xminus = xerr - (EleArr + ele - 1)->G.DC.XUnit.X * normdisp;
            xminus -= (EleArr + ele - 1)->G.DC.YUnit.X * normdisp;
            xminus -= (EleArr + ele - 1)->G.DC.ZUnit.X * normdisp;
            double yminus = yerr - (EleArr + ele - 1)->G.DC.XUnit.Y * normdisp;
            yminus -= (EleArr + ele - 1)->G.DC.YUnit.Y * normdisp;
            yminus -= (EleArr + ele - 1)->G.DC.ZUnit.Y * normdisp;
            double zminus = zerr - (EleArr + ele - 1)->G.DC.XUnit.Z * normdisp;
            zminus -= (EleArr + ele - 1)->G.DC.YUnit.Z * normdisp;
            zminus -= (EleArr + ele - 1)->G.DC.ZUnit.Z * normdisp;
            globalP.X = xminus;
            globalP.Y = yminus;
            globalP.Z = zminus;
            PFAtPoint(&globalP, &Potential, &globalF);
            localF  // Flux in the ECS
                = RotateVector3D(&globalF, &(EleArr + ele - 1)->G.DC,
                                 global2local);
            double value2 = -localF.Y;
            double epsratio = (Epsilon2[prim] / Epsilon1[prim]);
            Err = epsratio - (value1 / value2);
            fprintf(fErr, "%d\t%d\t%d\t%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                    prim, ele, 4, 4, xerr, yerr, zerr, epsratio,
                    (value1 / value2), Err);
            if (fabs(Err) > fabs(MaxError)) {
              MaxError = Err;
              PrimitiveOfMaxError = prim;
              ElementOfMaxError = ele;
              xerrMax = xerr;
              yerrMax = yerr;
              zerrMax = zerr;
            }
          }
        }  // if rectangle

        if (DebugLevel == 1001) exit(1);
      }  // loop over primitives and elements

    fprintf(fErr,
            "#Error maximum on element %d  on primitive %d with x,y,z %lg, "
            "%lg, %lg and its magnitude is %lg.\n",
            ElementOfMaxError, PrimitiveOfMaxError, xerrMax, yerrMax, zerrMax,
            MaxError);
    fclose(fErr);
  }  // if(OptEstimateError)

  return (0);
}  // Solve ends here

// In case we are PostProcessing only the solution needs to be read in
// from an existing Solution file present in the BCOutDir
int ReadSolution(void) {
  char SolnFile[256];
  strcpy(SolnFile, BCOutDir);
  strcat(SolnFile, "/Soln.out");
  FILE *fSoln = fopen(SolnFile, "r");
  // assert(fSoln != NULL);
  if (fSoln == NULL) {
    neBEMMessage("ReadSoln - unable to open solution file.");
    return -1;
  }

  int itmp;
  double dtmp, sol;
  char instr[256];
  fgets(instr, 256, fSoln);
  for (int ele = 1; ele <= NbElements; ++ele) {
    fscanf(fSoln, "%d %lg %lg %lg %lg\n", &itmp, &dtmp, &dtmp, &dtmp, &sol);
    // assert(ele == itmp);
    if (ele != itmp) {
      neBEMMessage("ReadSolution - ele_itmp in ReadSolution");
      return -1;
    }
    (EleArr + ele - 1)->Solution = sol;
  }
  printf("\nReadSolution: Solution read in for all elements ...\n");
  fflush(stdout);

  if (NbConstraints) {
    if (OptSystemChargeZero) {
      fgets(instr, 256, fSoln);
      fscanf(fSoln, "%d %lg\n", &NbSystemChargeZero, &VSystemChargeZero);
      printf(
          "ReadSolution: Read in voltage shift to ensure system charge "
          "zero.\n");
    }
    if (NbFloatingConductors) {
      fgets(instr, 256, fSoln);
      fscanf(fSoln, "%d %lg\n", &NbFloatCon, &VFloatCon);
      printf("ReadSolution: Read in voltage on floating conductor.\n");
    }
    fflush(stdout);
  }  // if NbConstraints

  fclose(fSoln);

  // Find primitive related charge densities
  // OMPCheck - may be parallelized
  for (int prim = 1; prim <= NbPrimitives; ++prim) {
    double area = 0.0;  // need area of the primitive as well!
    AvChDen[prim] = 0.0;
    AvAsgndChDen[prim] = 0.0;

    for (int ele = ElementBgn[prim]; ele <= ElementEnd[prim]; ++ele) {
      const double dA = (EleArr + ele - 1)->G.dA;
      area += dA;
      AvChDen[prim] += (EleArr + ele - 1)->Solution * dA;
      AvAsgndChDen[prim] += (EleArr + ele - 1)->Assigned * dA;
    }

    AvChDen[prim] /= area;
    AvAsgndChDen[prim] /= area;
  }

  neBEMState = 9;

  return (0);
}  // end of neBEMReadSolution

// Error estimate is easy here - check that the integration yields unity!
// TryWtField is a script that elucidates the idea.
int WeightingFieldSolution(int NbPrimsWtField, int PrimListWtField[],
                           double solnarray[]) {
  // Check for the inverted matrix
  if (!InvMat) {
    printf(
        "WeightingFieldSolution: Capacitance matrix not in memory, can not "
        "calculate weighting charges.\n");
    return (-1);
  }

  for (int i = 1; i <= NbUnknowns; i++) solnarray[i] = 0.0;

  for (int ele = 1, InList; ele <= NbElements; ++ele) {
    int prim = (EleArr + ele - 1)->PrimitiveNb;

    InList = 0;  // assume that this prim is not in the list
    for (int primwtfl = 0; primwtfl < NbPrimsWtField; ++primwtfl) {
      if (prim == PrimListWtField[primwtfl]) {
        InList = 1;
        break;  // get out of the for loop
      }
    }  // for primwtfl

    if (InList) {
      for (int i = 1; i <= NbUnknowns; ++i) {
        solnarray[i] += InvMat[i][ele];
      }
    }
  }  // for ele

  return (0);
}  // WtFieldSolution ends

// Create a function for Reflect to get the effect of given mirrors - this
// function should provide the localP. Rest can be done within the calling
// routines.
// Inputs for the function: Ele array, count ele, MirrorTypes, MirrorDists,
// xfld, yfld, zfld, xsrc, ysrc, zsrc, primsrc (if necessary)
// Rest can be computed within the function: p1, n, Image, dist, MirroredDC,
// globalP
// Output from the function: localP (distance from the reflected point to
// the field point where the influence is being evaluated)
//
// Shift the reflected point to take care of mirrors not at origin
// The equation to such a plane is (n.X)*x + (n.Y)*y + (n.Z)*z = d
// where d is the perpendicular distance from the origin
// So, we refine mirror definition using this additional information
// which, according to Garfield, has to be deduced from the device geom
// In neBEM, this distance has been deduced in neBEMInterface.

// the `handed-ness' of the element can get altered whenever elments
// are not perpendicular to the mirror
// Signs for the `X' componets need to be changed.

Point3D ReflectPrimitiveOnMirror(char Axis, int primsrc, Point3D srcpt,
                                 Point3D fldpt, double distance,
                                 DirnCosn3D *MirroredDC) {
  Vector3D n;     // define mirror by a bi-vector perpendicular to it
  Point3D Image;  // reflected point by mirror assumed at origin

  switch (Axis) {
    case 'x':
    case 'X':
      n.X = 1.0;
      n.Y = 0.0;
      n.Z = 0.0;

      Image = ReflectPoint3DByMirrorAtOrigin(&srcpt, &n);
      Image.X += 2.0 * distance;

      MirroredDC->XUnit.X = -PrimDC[primsrc].XUnit.X;
      MirroredDC->XUnit.Y = PrimDC[primsrc].XUnit.Y;
      MirroredDC->XUnit.Z = PrimDC[primsrc].XUnit.Z;
      MirroredDC->YUnit.X = -PrimDC[primsrc].YUnit.X;
      MirroredDC->YUnit.Y = PrimDC[primsrc].YUnit.Y;
      MirroredDC->YUnit.Z = PrimDC[primsrc].YUnit.Z;
      MirroredDC->ZUnit.X = -PrimDC[primsrc].ZUnit.X;
      MirroredDC->ZUnit.Y = PrimDC[primsrc].ZUnit.Y;
      MirroredDC->ZUnit.Z = PrimDC[primsrc].ZUnit.Z;
      break;

    case 'y':
    case 'Y':
      n.X = 0.0;
      n.Y = 1.0;
      n.Z = 0.0;

      Image = ReflectPoint3DByMirrorAtOrigin(&srcpt, &n);
      Image.Y += 2.0 * distance;

      MirroredDC->XUnit.X = PrimDC[primsrc].XUnit.X;
      MirroredDC->XUnit.Y = -PrimDC[primsrc].XUnit.Y;
      MirroredDC->XUnit.Z = PrimDC[primsrc].XUnit.Z;
      MirroredDC->YUnit.X = PrimDC[primsrc].YUnit.X;
      MirroredDC->YUnit.Y = -PrimDC[primsrc].YUnit.Y;
      MirroredDC->YUnit.Z = PrimDC[primsrc].YUnit.Z;
      MirroredDC->ZUnit.X = PrimDC[primsrc].ZUnit.X;
      MirroredDC->ZUnit.Y = -PrimDC[primsrc].ZUnit.Y;
      MirroredDC->ZUnit.Z = PrimDC[primsrc].ZUnit.Z;
      break;

    case 'z':
    case 'Z':
      n.X = 0.0;
      n.Y = 0.0;
      n.Z = 1.0;

      Image = ReflectPoint3DByMirrorAtOrigin(&srcpt, &n);
      Image.Z += 2.0 * distance;

      MirroredDC->XUnit.X = PrimDC[primsrc].XUnit.X;
      MirroredDC->XUnit.Y = PrimDC[primsrc].XUnit.Y;
      MirroredDC->XUnit.Z = -PrimDC[primsrc].XUnit.Z;
      MirroredDC->YUnit.X = PrimDC[primsrc].YUnit.X;
      MirroredDC->YUnit.Y = PrimDC[primsrc].YUnit.Y;
      MirroredDC->YUnit.Z = -PrimDC[primsrc].YUnit.Z;
      MirroredDC->ZUnit.X = PrimDC[primsrc].ZUnit.X;
      MirroredDC->ZUnit.Y = PrimDC[primsrc].ZUnit.Y;
      MirroredDC->ZUnit.Z = -PrimDC[primsrc].ZUnit.Z;
      break;

    default:
      printf("Axis not chosen properly!!! No reflection occurred!\n");
      Image.X = srcpt.X;
      Image.Y = srcpt.Y;
      Image.Z = srcpt.Z;
  }  // switch on Axis ends

  // printf("Axis: %c, distance: %lg\n", Axis, distance);
  // printf("srcpt: %lg, %lg, %lg\n", srcpt.X, srcpt.Y, srcpt.Z);
  // printf("Image: %lg, %lg, %lg\n", Image.X, Image.Y, Image.Z);
  // getchar();

  // traslated to ECS origin, but axes direction as in global system
  Point3D globalP, localP;
  globalP.X = fldpt.X - Image.X;
  globalP.Y = fldpt.Y - Image.Y;
  globalP.Z = fldpt.Z - Image.Z;

  // entirely in ECS
  return (localP = RotatePoint3D(&globalP, MirroredDC, global2local));
}  // ReflectPrimtiveOnMirror ends here

Point3D ReflectOnMirror(char Axis, int elesrc, Point3D srcpt, Point3D fldpt,
                        double distance, DirnCosn3D *MirroredDC) {
  Vector3D n;     // define mirror by a bi-vector perpendicular to it
  Point3D Image;  // reflected point by mirror assumed at origin

  switch (Axis) {
    case 'x':
    case 'X':
      n.X = 1.0;
      n.Y = 0.0;
      n.Z = 0.0;

      Image = ReflectPoint3DByMirrorAtOrigin(&srcpt, &n);
      Image.X += 2.0 * distance;

      MirroredDC->XUnit.X = -(EleArr + elesrc - 1)->G.DC.XUnit.X;
      MirroredDC->XUnit.Y = (EleArr + elesrc - 1)->G.DC.XUnit.Y;
      MirroredDC->XUnit.Z = (EleArr + elesrc - 1)->G.DC.XUnit.Z;
      MirroredDC->YUnit.X = -(EleArr + elesrc - 1)->G.DC.YUnit.X;
      MirroredDC->YUnit.Y = (EleArr + elesrc - 1)->G.DC.YUnit.Y;
      MirroredDC->YUnit.Z = (EleArr + elesrc - 1)->G.DC.YUnit.Z;
      MirroredDC->ZUnit.X = -(EleArr + elesrc - 1)->G.DC.ZUnit.X;
      MirroredDC->ZUnit.Y = (EleArr + elesrc - 1)->G.DC.ZUnit.Y;
      MirroredDC->ZUnit.Z = (EleArr + elesrc - 1)->G.DC.ZUnit.Z;
      break;

    case 'y':
    case 'Y':
      n.X = 0.0;
      n.Y = 1.0;
      n.Z = 0.0;

      Image = ReflectPoint3DByMirrorAtOrigin(&srcpt, &n);
      Image.Y += 2.0 * distance;

      MirroredDC->XUnit.X = (EleArr + elesrc - 1)->G.DC.XUnit.X;
      MirroredDC->XUnit.Y = -(EleArr + elesrc - 1)->G.DC.XUnit.Y;
      MirroredDC->XUnit.Z = (EleArr + elesrc - 1)->G.DC.XUnit.Z;
      MirroredDC->YUnit.X = (EleArr + elesrc - 1)->G.DC.YUnit.X;
      MirroredDC->YUnit.Y = -(EleArr + elesrc - 1)->G.DC.YUnit.Y;
      MirroredDC->YUnit.Z = (EleArr + elesrc - 1)->G.DC.YUnit.Z;
      MirroredDC->ZUnit.X = (EleArr + elesrc - 1)->G.DC.ZUnit.X;
      MirroredDC->ZUnit.Y = -(EleArr + elesrc - 1)->G.DC.ZUnit.Y;
      MirroredDC->ZUnit.Z = (EleArr + elesrc - 1)->G.DC.ZUnit.Z;
      break;

    case 'z':
    case 'Z':
      n.X = 0.0;
      n.Y = 0.0;
      n.Z = 1.0;

      Image = ReflectPoint3DByMirrorAtOrigin(&srcpt, &n);
      Image.Z += 2.0 * distance;

      MirroredDC->XUnit.X = (EleArr + elesrc - 1)->G.DC.XUnit.X;
      MirroredDC->XUnit.Y = (EleArr + elesrc - 1)->G.DC.XUnit.Y;
      MirroredDC->XUnit.Z = -(EleArr + elesrc - 1)->G.DC.XUnit.Z;
      MirroredDC->YUnit.X = (EleArr + elesrc - 1)->G.DC.YUnit.X;
      MirroredDC->YUnit.Y = (EleArr + elesrc - 1)->G.DC.YUnit.Y;
      MirroredDC->YUnit.Z = -(EleArr + elesrc - 1)->G.DC.YUnit.Z;
      MirroredDC->ZUnit.X = (EleArr + elesrc - 1)->G.DC.ZUnit.X;
      MirroredDC->ZUnit.Y = (EleArr + elesrc - 1)->G.DC.ZUnit.Y;
      MirroredDC->ZUnit.Z = -(EleArr + elesrc - 1)->G.DC.ZUnit.Z;
      break;

    default:
      printf("Axis not chosen properly!!! No reflection occurred!\n");
      Image.X = srcpt.X;
      Image.Y = srcpt.Y;
      Image.Z = srcpt.Z;
  }  // switch on Axis ends

  // printf("Axis: %c, distance: %lg\n", Axis, distance);
  // printf("srcpt: %lg, %lg, %lg\n", srcpt.X, srcpt.Y, srcpt.Z);
  // printf("Image: %lg, %lg, %lg\n", Image.X, Image.Y, Image.Z);
  // getchar();

  // traslated to ECS origin, but axes direction as in global system
  Point3D globalP, localP;
  globalP.X = fldpt.X - Image.X;
  globalP.Y = fldpt.Y - Image.Y;
  globalP.Z = fldpt.Z - Image.Z;

  // entirely in ECS
  return (localP = RotatePoint3D(&globalP, MirroredDC, global2local));
}  // ReflectOnMirror ends

int UpdateKnownCharges(void) { return 0; }  // UpdateKnownCharges ends

int UpdateChargingUp(void) { return 0; }  // UpdateChargingUp ends

#ifdef __cplusplus
}  // namespace
#endif