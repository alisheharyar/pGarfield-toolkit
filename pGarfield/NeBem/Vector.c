/*
(c) 2005, Supratik Mukhopadhayay, Nayana Majumdar
*/
// Note that Opt = 1 yields a forward rotation while -1 yields a -ve rotation.
#define DEFINE_VGLOBAL

#include "Vector.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
namespace neBEM {
#endif

void VectorRotate_Rect3D(double Xin, double Yin, double Zin, double RotX,
                         double RotY, double RotZ,
                         int Opt,  // 1 forward, -1 backward
                         double *Xout, double *Yout, double *Zout) {
  double R11, R12, R13,  // Follow Numerical Methods Using Matlab
      R21, R22, R23,     // J.H.Mathews and K.D.Fink
      R31, R32, R33;     // 4th Edition, Prentice-Hall of India Pvt. Ltd.
  double X_X, Y_X, Z_X,  // New Delhi, 2004
      X_Y, Y_Y, Z_Y;     // p.115

  if ((fabs(RotX) < 1.0e-12)  // A most happy unrotated situation
      && (fabs(RotY) < 1.0e-12) && (fabs(RotZ) < 1.0e-12)) {
    *Xout = Xin;
    *Yout = Yin;
    *Zout = Zin;
    return;
  }

  // Rotation in X
  if (fabs(RotX) < 1.0e-12) {
    X_X = Xin;
    Y_X = Yin;
    Z_X = Zin;
  } else {
    R11 = 1.0;
    R12 = 0.0;
    R13 = 0.0;
    R21 = 0.0;
    R22 = cos(Opt * RotX);
    R23 = -sin(Opt * RotX);
    R31 = 0.0;
    R32 = sin(Opt * RotX);
    R33 = cos(Opt * RotX);
    X_X = R11 * Xin + R12 * Yin + R13 * Zin;  // separate matrix multiplication
    Y_X = R21 * Xin + R22 * Yin + R23 * Zin;  // function avoided intentionally.
    Z_X = R31 * Xin + R32 * Yin + R33 * Zin;
  }

  // Rotation in Y
  if (fabs(RotY) < 1.0e-12) {
    X_Y = X_X;
    Y_Y = Y_X;
    Z_Y = Z_X;
  } else {
    R11 = cos(Opt * RotY);
    R12 = 0.0;
    R13 = sin(Opt * RotY);
    R21 = 0.0;
    R22 = 1.0;
    R23 = 0.0;
    R31 = -sin(Opt * RotY);
    R32 = 0.0;
    R33 = cos(Opt * RotY);
    X_Y = R11 * X_X + R12 * Y_X + R13 * Z_X;  // separate matrix multiplication
    Y_Y = R21 * X_X + R22 * Y_X + R23 * Z_X;  // function avoided intentionally.
    Z_Y = R31 * X_X + R32 * Y_X + R33 * Z_X;
  }

  // Rotation in Z - final rotation
  if (fabs(RotZ) < 1.0e-12) {
    *Xout = X_Y;
    *Yout = Y_Y;
    *Zout = Z_Y;
  } else {
    R11 = cos(Opt * RotZ);
    R12 = -sin(Opt * RotZ);
    R13 = 0.0;
    R21 = sin(Opt * RotZ);
    R22 = cos(Opt * RotZ);
    R23 = 0.0;
    R31 = 0.0;
    R32 = 0.0;
    R33 = 1.0;
    *Xout =
        R11 * X_Y + R12 * Y_Y + R13 * Z_Y;  // separate matrix multiplication
    *Yout =
        R21 * X_Y + R22 * Y_Y + R23 * Z_Y;  // function avoided intentionally.
    *Zout = R31 * X_Y + R32 * Y_Y + R33 * Z_Y;
  }
}

// Rotation of coordinate system is effectively oppsite to rotating the vector
// keeping the coordinate system the same, as in the above routine.
void CoordRotate_Rect3D(double Xin, double Yin, double Zin, double RotX,
                        double RotY, double RotZ,
                        int Opt,  // 1 forward, -1 backward
                        double *Xout, double *Yout, double *Zout) {
  double R11, R12, R13,  // Follow Numerical Methods Using Matlab
      R21, R22, R23,     // J.H.Mathews and K.D.Fink
      R31, R32, R33;     // 4th Edition, Prentice-Hall of India Pvt. Ltd.
  double X_X, Y_X, Z_X,  // New Delhi, 2004
      X_Y, Y_Y, Z_Y;     // p.115

  if ((fabs(RotX) < 1.0e-12)  // A most happy unrotated situation
      && (fabs(RotY) < 1.0e-12) && (fabs(RotZ) < 1.0e-12)) {
    *Xout = Xin;
    *Yout = Yin;
    *Zout = Zin;
    return;
  }

  // Rotation in X
  if (fabs(RotX) < 1.0e-12) {
    X_X = Xin;
    Y_X = Yin;
    Z_X = Zin;
  } else {
    R11 = 1.0;
    R12 = 0.0;
    R13 = 0.0;
    R21 = 0.0;
    R22 = cos(Opt * RotX);
    R23 = sin(Opt * RotX);
    R31 = 0.0;
    R32 = -sin(Opt * RotX);
    R33 = cos(Opt * RotX);
    X_X = R11 * Xin + R12 * Yin + R13 * Zin;  // separate matrix multiplication
    Y_X = R21 * Xin + R22 * Yin + R23 * Zin;  // function avoided intentionally.
    Z_X = R31 * Xin + R32 * Yin + R33 * Zin;
  }

  // Rotation in Y
  if (fabs(RotY) < 1.0e-12) {
    X_Y = X_X;
    Y_Y = Y_X;
    Z_Y = Z_X;
  } else {
    R11 = cos(Opt * RotY);
    R12 = 0.0;
    R13 = -sin(Opt * RotY);
    R21 = 0.0;
    R22 = 1.0;
    R23 = 0.0;
    R31 = sin(Opt * RotY);
    R32 = 0.0;
    R33 = cos(Opt * RotY);
    X_Y = R11 * X_X + R12 * Y_X + R13 * Z_X;  // separate matrix multiplication
    Y_Y = R21 * X_X + R22 * Y_X + R23 * Z_X;  // function avoided intentionally.
    Z_Y = R31 * X_X + R32 * Y_X + R33 * Z_X;
  }

  // Rotation in Z - final rotation
  if (fabs(RotZ) < 1.0e-12) {
    *Xout = X_Y;
    *Yout = Y_Y;
    *Zout = Z_Y;
  } else {
    R11 = cos(Opt * RotZ);
    R12 = sin(Opt * RotZ);
    R13 = 0.0;
    R21 = -sin(Opt * RotZ);
    R22 = cos(Opt * RotZ);
    R23 = 0.0;
    R31 = 0.0;
    R32 = 0.0;
    R33 = 1.0;
    *Xout =
        R11 * X_Y + R12 * Y_Y + R13 * Z_Y;  // separate matrix multiplication
    *Yout =
        R21 * X_Y + R22 * Y_Y + R23 * Z_Y;  // function avoided intentionally.
    *Zout = R31 * X_Y + R32 * Y_Y + R33 * Z_Y;
  }
}

// Create a 3D point
Point3D CreatePoint3D(double x, double y, double z) {
  Point3D p;
  p.X = x;
  p.Y = y;
  p.Z = z;

  return (p);
}

// Get distance between two 3D poionts
double GetDistancePoint3D(Point3D *a, Point3D *b) {
  return (sqrt((b->X - a->X) * (b->X - a->X) + (b->Y - a->Y) * (b->Y - a->Y) +
               (b->Z - a->Z) * (b->Z - a->Z)));
}

// Create a 3D distance vector from point a to point b
Vector3D CreateDistanceVector3D(Point3D *a, Point3D *b) {
  Vector3D v;

  v.X = b->X - a->X;
  v.Y = b->Y - a->Y;
  v.Z = b->Z - a->Z;

  return (v);
}

// Computes magnitude of a 3D vector
double MagVector3D(Vector3D *A) {
  double mag;

  mag = sqrt(A->X * A->X + A->Y * A->Y + A->Z * A->Z);

  return (mag);
}

// Computes vector dot product
double Vector3DDotProduct(Vector3D *A, Vector3D *B) {
  double product;

  product = A->X * B->X + A->Y * B->Y + A->Z * B->Z;

  return (product);
}

// Computes unit vector for a given vector
Vector3D UnitVector3D(Vector3D *v) {
  Vector3D u;

  double mag = MagVector3D(v);
  if (fabs(mag) <= 1.0e-12) {
    printf("UnitVector3D: magnitude smaller than 1.0e-12; no normalization.\n");
    u.X = v->X;
    u.Y = v->Y;
    u.Z = v->Z;
  } else {
    u.X = v->X / mag;
    u.Y = v->Y / mag;
    u.Z = v->Z / mag;
  }

  return (u);
}

// Computes vector cross product
Vector3D Vector3DCrossProduct(Vector3D *A, Vector3D *B) {
  Vector3D product;

  product.X = A->Y * B->Z - A->Z * B->Y;
  product.Y = A->Z * B->X - A->X * B->Z;
  product.Z = A->X * B->Y - A->Y * B->X;

  return (product);
}

// prints out co-ordinates of a 3D point
int PrintPoint3D(Point3D A) {
  printf("%lg %lg %lg", A.X, A.Y, A.Z);
  return (0);
}

// prints out components of a 3D vector
int PrintVector3D(Vector3D A) {
  printf("%lg %lg %lg", A.X, A.Y, A.Z);
  return (0);
}

// prints out components of 3D direction cosines
int PrintDirnCosn3D(DirnCosn3D A) {
  printf("XUnit: ");
  PrintVector3D(A.XUnit);
  printf("\n");
  printf("YUnit: ");
  PrintVector3D(A.YUnit);
  printf("\n");
  printf("ZUnit: ");
  PrintVector3D(A.ZUnit);
  printf("\n");
  return (0);
}

// Translates a point to a new origin specified in terms of the existing system
// Translating the axes by (xT, yT, zT) is equivalent to translating the point
// by (-xT, -yT, -zT)
Point3D TranslatePoint3D(Point3D *A, Point3D *Origin, int Sense) {
  double InitialVector[4];
  double TranslationMatrix[4][4] = {{1.0, 0.0, 0.0, 0.0},
                                    {0.0, 1.0, 0.0, 0.0},
                                    {0.0, 0.0, 1.0, 0.0},
                                    {0.0, 0.0, 0.0, 1.0}};
  double FinalVector[4];
  Point3D TranslatedPt;

  InitialVector[0] = A->X;
  InitialVector[1] = A->Y;
  InitialVector[2] = A->Z;
  InitialVector[3] = 1.0;

  switch (Sense) {
    case 1:
      TranslationMatrix[0][3] = -Origin->X;
      TranslationMatrix[1][3] = -Origin->Y;
      TranslationMatrix[2][3] = -Origin->Z;
      break;

    case -1:
      TranslationMatrix[0][3] = Origin->X;
      TranslationMatrix[1][3] = Origin->Y;
      TranslationMatrix[2][3] = Origin->Z;
      break;

    default:
      printf("Only forward and inverse senses are allowed ...\n");
      exit(-1);
  }

  for (int i = 0; i < 4; ++i) {
    FinalVector[i] = 0.0;
    for (int j = 0; j < 4; ++j) {
      FinalVector[i] += TranslationMatrix[i][j] * InitialVector[j];
    }
  }

  TranslatedPt.X = FinalVector[0];
  TranslatedPt.Y = FinalVector[1];
  TranslatedPt.Z = FinalVector[2];
  return (TranslatedPt);
}

// Position of a point in a coord system with direction cosines specified
// in terms of the existing system
// Theory from Rotation representation (Wikipedia)
// Sense +1 implies we are moving towards the new coordinate system whose
// direction cosines we know in terms of the original coordinate system.
// This option can be invoked by choosing 'global2local' defined in Vector.h
// Sense -1 implies we are moving towards the original coordinate system in
// which we know the direction cosines of the new coordinate system.
// This option can be invoked by choosing 'local2global' defined in Vector.h
Point3D RotatePoint3D(Point3D *A, DirnCosn3D *DC, int Sense) {
  double TransformationMatrix[3][3] = {
      {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  switch (Sense) {
    case 1:
      TransformationMatrix[0][0] = DC->XUnit.X;
      TransformationMatrix[0][1] = DC->XUnit.Y;
      TransformationMatrix[0][2] = DC->XUnit.Z;
      TransformationMatrix[1][0] = DC->YUnit.X;
      TransformationMatrix[1][1] = DC->YUnit.Y;
      TransformationMatrix[1][2] = DC->YUnit.Z;
      TransformationMatrix[2][0] = DC->ZUnit.X;
      TransformationMatrix[2][1] = DC->ZUnit.Y;
      TransformationMatrix[2][2] = DC->ZUnit.Z;
      break;

    case -1:
      TransformationMatrix[0][0] = DC->XUnit.X;
      TransformationMatrix[0][1] = DC->YUnit.X;
      TransformationMatrix[0][2] = DC->ZUnit.X;
      TransformationMatrix[1][0] = DC->XUnit.Y;
      TransformationMatrix[1][1] = DC->YUnit.Y;
      TransformationMatrix[1][2] = DC->ZUnit.Y;
      TransformationMatrix[2][0] = DC->XUnit.Z;
      TransformationMatrix[2][1] = DC->YUnit.Z;
      TransformationMatrix[2][2] = DC->ZUnit.Z;
      break;

    default:
      printf("Only forward and inverse senses are allowed ...\n");
      exit(-1);
  }

  double InitialVector[3] = {A->X, A->Y, A->Z};
  double FinalVector[3] = {0., 0., 0.};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      FinalVector[i] += TransformationMatrix[i][j] * InitialVector[j];
    }
  }
  Point3D RotatedPt;
  RotatedPt.X = FinalVector[0];
  RotatedPt.Y = FinalVector[1];
  RotatedPt.Z = FinalVector[2];
  return (RotatedPt);
}

// Vector in a coord system with direction cosines specified
// in terms of the existing system
// Theory from Rotation representation (Wikipedia)
// Sense +1 implies we are moving towards the new coordinate system whose
// direction cosines we know in terms of the original coordinate system.
// This option can be invoked by choosing 'global2local' defined in Vector.h
// Sense -1 implies we are moving towards the original coordinate system in
// which we know the direction cosines of the new coordinate system.
// This option can be invoked by choosing 'local2global' defined in Vector.h
Vector3D RotateVector3D(Vector3D *A, DirnCosn3D *DC, int Sense) {
  double TransformationMatrix[4][4] = {{0.0, 0.0, 0.0, 0.0},
                                       {0.0, 0.0, 0.0, 0.0},
                                       {0.0, 0.0, 0.0, 0.0},
                                       {0.0, 0.0, 0.0, 1.0}};
  switch (Sense) {
    case 1:
      TransformationMatrix[0][0] = DC->XUnit.X;
      TransformationMatrix[0][1] = DC->XUnit.Y;
      TransformationMatrix[0][2] = DC->XUnit.Z;
      TransformationMatrix[1][0] = DC->YUnit.X;
      TransformationMatrix[1][1] = DC->YUnit.Y;
      TransformationMatrix[1][2] = DC->YUnit.Z;
      TransformationMatrix[2][0] = DC->ZUnit.X;
      TransformationMatrix[2][1] = DC->ZUnit.Y;
      TransformationMatrix[2][2] = DC->ZUnit.Z;
      break;

    case -1:
      TransformationMatrix[0][0] = DC->XUnit.X;
      TransformationMatrix[0][1] = DC->YUnit.X;
      TransformationMatrix[0][2] = DC->ZUnit.X;
      TransformationMatrix[1][0] = DC->XUnit.Y;
      TransformationMatrix[1][1] = DC->YUnit.Y;
      TransformationMatrix[1][2] = DC->ZUnit.Y;
      TransformationMatrix[2][0] = DC->XUnit.Z;
      TransformationMatrix[2][1] = DC->YUnit.Z;
      TransformationMatrix[2][2] = DC->ZUnit.Z;
      break;

    default:
      printf("Only forward and inverse senses are allowed ...\n");
      exit(-1);
  }

  double InitialVector[3] = {A->X, A->Y, A->Z};
  double FinalVector[3] = {0., 0., 0.};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      FinalVector[i] += TransformationMatrix[i][j] * InitialVector[j];
    }
  }
  Vector3D RotatedVector;
  RotatedVector.X = FinalVector[0];
  RotatedVector.Y = FinalVector[1];
  RotatedVector.Z = FinalVector[2];
  return (RotatedVector);
}

// Transform point: Get the new coordinates of a point described in the original
//                  coordinate system. The new system of coordinates is also
//                  described in terms of the original coordiate system, the
//                  origin being NewOrigin and new direction cosines being
//                  termed as NewDirns.
//                  Translate-Rotate sequence
// Inputs: point in the original coordinate system, origin of the new system
//         w.r.t the original coord system, direction cosines of the new
//         system w.r.t the original coord system
Point3D TransformPoint3D(Point3D *initial, Point3D *NewOrigin,
                         DirnCosn3D *NewDirns) {
  Point3D TmpPoint, final;

  TmpPoint = TranslatePoint3D(initial, NewOrigin, 1);
  final = RotatePoint3D(&TmpPoint, NewDirns, 1);
  return (final);
}

// Reflect a 3D point (p1) on a mirror which is defined by a bi-vector (n)
// perpendicular to the mirror
// It is assumed that the mirror passes through the origin
Point3D ReflectPoint3DByMirrorAtOrigin(Point3D *p1, Vector3D *n) {
  double matrix[3][3];
  matrix[0][0] = -n->X * n->X + n->Y * n->Y + n->Z * n->Z;
  matrix[0][1] = -2.0 * n->X * n->Y;
  matrix[0][2] = -2.0 * n->X * n->Z;
  matrix[1][0] = -2.0 * n->X * n->Y;
  matrix[1][1] = n->X * n->X - n->Y * n->Y + n->Z * n->Z;
  matrix[1][2] = -2.0 * n->Y * n->Z;
  matrix[2][0] = -2.0 * n->X * n->Z;
  matrix[2][1] = -2.0 * n->Y * n->Z;
  matrix[2][2] = n->X * n->X + n->Y * n->Y - n->Z * n->Z;

  Point3D p2;  // reflected point
  p2.X = matrix[0][0] * p1->X + matrix[0][1] * p1->Y + matrix[0][2] * p1->Z;
  p2.Y = matrix[1][0] * p1->X + matrix[1][1] * p1->Y + matrix[1][2] * p1->Z;
  p2.Z = matrix[2][0] * p1->X + matrix[2][1] * p1->Y + matrix[2][2] * p1->Z;

  return (p2);
}  // ReflectPoint3DByMirrorAtOrigin ends

#ifdef __cplusplus
}  // namespace
#endif
