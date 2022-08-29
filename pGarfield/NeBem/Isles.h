/*
(c) 2005, Supratik Mukhopadhayay, Nayana Majumdar
*/
#ifndef _Isles_H_
#define _Isles_H_

#ifdef DEFINE_ISLESGLOBAL
#define ISLESGLOBAL
#else
#define ISLESGLOBAL extern
#endif

#include <math.h>
#include <stdio.h>

#include "Vector.h"

#define MINDIST 1.0e-8
#define MINDIST2 1.0e-16
#define MINDIST3 1.0e-20

#define FarField 10.0  // Beyond FarField*ElementSize, use nodal expressions
#define FarField2 100.0

#define ST_PI 3.14159265358979323846

#ifdef __cplusplus
namespace neBEM {
#endif

ISLESGLOBAL char ISLESVersion[10];

// Nb of times these get called: FailureCntr keeps a tab of how many times
// the exact expressions failed to evaluate.
ISLESGLOBAL int IslesCntr, ExactCntr, ApproxCntr, FailureCntr;
ISLESGLOBAL int ApproxFlag, DebugISLES;
ISLESGLOBAL FILE *fIsles;

// ECS: Element coordinate system
// GCS: Global coordinate system

// X, Y, Z defines the field point in the ECS whose origin is at the
// centre of the rectangular element, the element being in the XZ plane.
// lo and hi define the extent of the element.
ISLESGLOBAL int ExactRecSurf(double X, double Y, double Z, double xlo,
                             double zlo, double xhi, double zhi,
                             double *Potential, Vector3D *Flux);

ISLESGLOBAL int ApproxRecSurf(double X, double Y, double Z, double xlo,
                              double zlo, double xhi, double zhi, int xseg,
                              int zseg, double *Potential, Vector3D *Flux);

// X, Y, Z defines the field point in the ECS whose origin is at the right
// angle corner of a triangular element, the element being on the XZ plane.
// x extent is always from 0 to 1, while z is from 0 to zmax
ISLESGLOBAL int ExactTriSurf(double zMax, double X, double Y, double Z,
                             double *Potential, Vector3D *Flux);

ISLESGLOBAL int ApproxTriSurf(double zMax, double X, double Y, double Z,
                              int xseg, int zseg, double *Potential,
                              Vector3D *Flux);

// potential at the surfce of the thin wire
ISLESGLOBAL double ExactCentroidalP_W(double rW, double lW);

// potential along the axis of the thin wire
ISLESGLOBAL double ExactAxialP_W(double rW, double lW, double Z);

// axial field along the axis of the wire
ISLESGLOBAL double ExactAxialFZ_W(double rW, double lW, double Z);

ISLESGLOBAL double ApproxP_W(double rW, double lW, double X, double Y, double Z,
                             int zseg);

ISLESGLOBAL double ApproxFX_W(double rW, double lW, double X, double Y,
                              double Z, int zseg);

ISLESGLOBAL double ApproxFY_W(double rW, double lW, double X, double Y,
                              double Z, int zseg);

ISLESGLOBAL double ApproxFZ_W(double rW, double lW, double X, double Y,
                              double Z, int zseg);

ISLESGLOBAL int ApproxWire(double rW, double lW, double X, double Y, double Z,
                           int zseg, double *potential, Vector3D *Flux);

ISLESGLOBAL double ImprovedP_W(double rW, double lW, double X, double Y,
                               double Z);

ISLESGLOBAL double ImprovedFX_W(double rW, double lW, double X, double Y,
                                double Z);

ISLESGLOBAL double ImprovedFY_W(double rW, double lW, double X, double Y,
                                double Z);

ISLESGLOBAL double ImprovedFZ_W(double rW, double lW, double X, double Y,
                                double Z);

ISLESGLOBAL int ImprovedWire(double rW, double lW, double X, double Y, double Z,
                             double *potential, Vector3D *Flux);

// Potential and fluxes at a field point in the ECS due to a thin wire.
// The four functions yield the mentioned four parameters separately.
ISLESGLOBAL double ExactThinP_W(double rW, double lW, double X, double Y,
                                double Z);

ISLESGLOBAL double ExactThinFX_W(double rW, double lW, double X, double Y,
                                 double Z);

ISLESGLOBAL double ExactThinFY_W(double rW, double lW, double X, double Y,
                                 double Z);

ISLESGLOBAL double ExactThinFZ_W(double rW, double lW, double X, double Y,
                                 double Z);

ISLESGLOBAL int ExactThinWire(double rW, double lW, double X, double Y,
                              double Z, double *potential, Vector3D *Flux);

// The point is passed in the local coordinate system, in the ECS of the ring
// In the ECS, plane of the ring is r,theta (XY) and z is the ordinate
// up and down of this plane.
ISLESGLOBAL int ExactRingPF(double rW, Point3D localPt, double *potential,
                            Vector3D *Flux);

// The point is passed in the local coordinate system, in the ECS of the disc
// In the ECS, plane of the ring is r,theta (XY) and z is the ordinate
// up and down of this plane.
ISLESGLOBAL int ExactDiscPF(double rW, Point3D localPt, double *potential,
                            Vector3D *Flux);

// All coordinates and vectors are in the global co-ordinate system.
// In case of area (both triangular and rectangular), it is assumed that
// the +ve normal direction is obtained by
// \vec{P1-P0} \cross \vec{P2-P1}.
// Note that, according to our convention, the second point of a right
// triangle \vec{P1} is the right corner of the triangle.
ISLESGLOBAL double PointKnChPF(Point3D SourcePt, Point3D FieldPt,
                               Vector3D *globalF);

ISLESGLOBAL double LineKnChPF(Point3D LineStart, Point3D LineStop,
                              Point3D FieldPt, Vector3D *globalF);

ISLESGLOBAL double WireKnChPF(Point3D WireStart, Point3D WireStop,
                              double radius, Point3D FieldPt,
                              Vector3D *globalF);

ISLESGLOBAL double AreaKnChPF(int NbVertices, Point3D *Vertices,
                              Point3D FieldPt, Vector3D *globalF);

ISLESGLOBAL double ApproxVolumeKnChPF(int NbPts, Point3D *SourcePt,
                                      Point3D FieldPt, Vector3D *globalF);

ISLESGLOBAL double VolumeKnChPF(int NbPts, Point3D *Vertices, Point3D FieldPt,
                                Vector3D *globalF);

ISLESGLOBAL int Sign(double x);

#ifdef __cplusplus
}  // namespace
#endif

#endif
