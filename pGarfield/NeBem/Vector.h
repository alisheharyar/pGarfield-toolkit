/*
(c) 2005, Supratik Mukhopadhayay, Nayana Majumdar
*/
#ifndef _Vector_H_
#define _Vector_H_

#ifdef DEFINE_VGLOBAL
#define VGLOBAL
#else
#define VGLOBAL extern
#endif

#define global2local 1
#define local2global -1

#ifdef __cplusplus
namespace neBEM {
#endif

// three coordinates
typedef struct {
  double X;
  double Y;
  double Z;
} Point3D;

// three vector components
typedef struct {
  double X;
  double Y;
  double Z;
} Vector3D;

// direction cosine
typedef struct {
  Vector3D XUnit;
  Vector3D YUnit;
  Vector3D ZUnit;
} DirnCosn3D;

// Creates a 3D point {x,y,z} - note sequence
VGLOBAL Point3D CreatePoint3D(double x, double y, double z);

// Get the distance between two 3D points
VGLOBAL double GetDistancePoint3D(Point3D *a, Point3D *b);

// Create a 3D distance vector from point a to point b - note sequence
VGLOBAL Vector3D CreateDistanceVector3D(Point3D *a, Point3D *b);

// Compute magnitude of a 3D vector
VGLOBAL double MagVector3D(Vector3D *);

// Create the unit vector of vector a
VGLOBAL Vector3D UnitVector3D(Vector3D *a);

// prints out co-ordinates of a 3D point
VGLOBAL int PrintPoint3D(Point3D);

// prints out components of a 3D vector
VGLOBAL int PrintVector3D(Vector3D A);

// prints out components of 3D direction cosines
VGLOBAL int PrintDirnCosn3D(DirnCosn3D A);

// Rotates a vector
VGLOBAL void VectorRotate_Rect3D(double Xin, double Yin, double Zin,
                                 double RotX, double RotY, double RotZ,
                                 int Opt,  // 1 forward, -1 backward
                                 double *Xout, double *Yout, double *Zout);

// Transforms a point
VGLOBAL void CoordRotate_Rect3D(double Xin, double Yin, double Zin, double RotX,
                                double RotY, double RotZ,
                                int Opt,  // 1 forward, -1 backward
                                double *Xout, double *Yout, double *Zout);

// Carries out a vector dot-product
VGLOBAL double Vector3DDotProduct(Vector3D *, Vector3D *);

// Carries out a vector cross-product - be sure that the order is correct!
VGLOBAL Vector3D Vector3DCrossProduct(Vector3D *, Vector3D *);

// Translates a point to a new origin specified in terms of the existing system
VGLOBAL Point3D TranslatePoint3D(Point3D *A, Point3D *Origin, int Sense);

// Position of a point in a coord system with direction cosines specified
// in terms of the existing system
VGLOBAL Point3D RotatePoint3D(Point3D *A, DirnCosn3D *Origin, int Sense);

// Vector in a coord system with direction cosines specified
// in terms of the existing system
VGLOBAL Vector3D RotateVector3D(Vector3D *A, DirnCosn3D *Origin, int Sense);

// Transform point: Get the new coordinates of a point described in the
// original//                  coordinate system. The new system of coordinates
// is also
//                  described in terms of the original coordiate system, the
//                  origin being NewOrigin and new direction cosines being
//                  termed as NewDirns.
//                  Translate-Rotate sequence
// Inputs: point in the original coordinate system, origin of the new system
//         w.r.t the original coord system, direction cosines of the new
//         system w.r.t the original coord system
VGLOBAL Point3D TransformPoint3D(Point3D *initial, Point3D *NewOrigin,
                                 DirnCosn3D *NewDirns);

// Reflect a point using a mirror passing through the origin
VGLOBAL Point3D ReflectPoint3DByMirrorAtOrigin(Point3D *p1, Vector3D *n);

// Product of two matrices
VGLOBAL double **MatrixProduct(int NbRows1, int NbCols1, double **Matrix1,
                               int NbRows2, int NbCols2, double **Matrix2);

#ifdef __cplusplus
}  // namespace
#endif

#endif
