/*
(c) 2005, Supratik Mukhopadhayay, Nayana Majumdar
*/

#include <stdio.h>
#include <time.h>
#include <unistd.h>

#include "Isles.h"
#include "NR.h"
#include "Vector.h"
#include "neBEM.h"
#include "neBEMInterface.h"

#define MyPI 3.14159265358979323846

#ifdef __cplusplus
namespace neBEM {
#endif

// GCS: global coordinate system
// PCS: primitive coordinate system
// ECS: element coordinate system

// How do we decide on the number of elements, in each direction, appropriate
// for a given surface?
// Since no linear element is being considered, no assumption is being made
// regarding a linear element for the time being.
// shall we return the pointer to the element array here? That seems to be
// more intuitive!
// Note that for the right triangle, the second vertex (vertex[1], since the
// vector begins from 0) is the 90 degree corner.
int SurfaceElements(int prim, int nvertex, double xvert[], double yvert[],
                    double zvert[], double xnorm, double ynorm, double znorm,
                    int volref1, int volref2, int inttype, double potential,
                    double charge, double lambda, int NbSegX, int NbSegZ) {
  // Decide the geometry of this primitive - if it is a rectangle, our job is
  // greatly simplified. To begin with, we check the number of vertices to
  // take the decision automatically.
  // Note that a triangle is the next best bet. All other primitive will have to
  // be reduced to rectangles (as many as possible) and triangles.
  // Incidentally, a PrimType (determined in neBEMReadGeom (neBEMInterface.c))
  // determines whether the primitive is a surface (2D) or a wire (1D),
  // for a given surface, SurfShape determines whether it is a triangle,
  // rectangle, or any other polygon besides these two (square is, of course, a
  // special case of rectangle)
  int fstatus;
  switch (nvertex) {
    case 3:  // triangle
      fstatus = DiscretizeTriangle(prim, nvertex, xvert, yvert, zvert, xnorm,
                                   ynorm, znorm, volref1, volref2, inttype,
                                   potential, charge, lambda, NbSegX, NbSegZ);
      // assert(fstatus == 0);
      if (fstatus != 0) {
        neBEMMessage("SurfaceElements - DiscretizeTriangle");
        return -1;
      }
      break;

    case 4:  // rectangle
      fstatus = DiscretizeRectangle(prim, nvertex, xvert, yvert, zvert, xnorm,
                                    ynorm, znorm, volref1, volref2, inttype,
                                    potential, charge, lambda, NbSegX, NbSegZ);
      // assert(fstatus == 0);
      if (fstatus != 0) {
        neBEMMessage("SurfaceElements - DiscretizeRectangle");
        return -1;
      }
      break;

    default:
      printf("nvertex out of bounds in SurfaceElements ... exiting ...\n");
      exit(-1);
  }

  return (0);
}  // end of SurfaceElements

// Analyze wire and break up into smaller wire elements
// How do we decide on the number of elements? Currently, it is based on user
// inputs that need to be made automatic and adaptive
int WireElements(int prim, int nvertex, double xvert[], double yvert[],
                 double zvert[], double radius, int volref1, int volref2,
                 int inttype, double potential, double charge, double lambda,
                 int NbSegs) {
  int fstatus;

  switch (nvertex) {
    case 2:  // wire
      fstatus =
          DiscretizeWire(prim, nvertex, xvert, yvert, zvert, radius, volref1,
                         volref2, inttype, potential, charge, lambda, NbSegs);
      // assert(fstatus == 0);
      if (fstatus != 0) {
        neBEMMessage("WireElements - DiscretizeWire");
        return -1;
      }
      break;

    default:
      printf("nvertex out of bounds in WireElements ... exiting ...\n");
      exit(-1);
  }

  return (0);
}  // end of WireElement

// Try to set up elements on this primitive such that the average element area
// is close to what has been requested.
// If the number of elements is too small (<5), over-ride and have 5*5 elements
// If the size of the element is too small (<MINDIST2), or a side of the element
// is too small (<MINDIST), over-ride 5*5 and create elements (3*3) and (1*1)
// such that they are of acceptable size and shape (the latter, if possible)
// If the primitive is smaller than MINDIST * MINDIST, report and quit
// If the number of elements is more than a limit, make the element size
// larger than that requested. Report the incidences.
// Since, for a surface primitive, the larger number between Coord1Seg and
// Coord2Seg is used for the longer arm (look for OverSmart in
// DiscretizeTriangle and DiscretizeRectagnle), the aspect ratio problem is
// expected to be taken care of to a large extent.
int AnalyzePrimitive(int prim, int *NbSegCoord1, int *NbSegCoord2) {
  int fstatus;

  switch (NbVertices[prim]) {
    case 2:  // wire
      fstatus = AnalyzeWire(prim, NbSegCoord1);
      *NbSegCoord2 = 0;
      // assert(fstatus == 0);
      if (fstatus != 0) {
        neBEMMessage("AnalyzePrimitive - AnalyzeWire");
        return -1;
      }
      return (2);
      break;
    case 3:  // triangle
      fstatus = AnalyzeSurface(prim, NbSegCoord1, NbSegCoord2);
      // assert(fstatus == 0);
      if (fstatus != 0) {
        neBEMMessage("AnalyzePrimitive - AnalyzeSurface");
        return -1;
      }
      return (3);
      break;
    case 4:  // rectangle
      fstatus = AnalyzeSurface(prim, NbSegCoord1, NbSegCoord2);
      // assert(fstatus == 0);
      if (fstatus != 0) {
        neBEMMessage("AnalyzePrimitive - AnalyzeSurface");
        return -1;
      }
      return (4);
      break;
    default:  // no shape!
      return (0);
  }
}  // end of AnalyzePrimitive

int AnalyzeWire(int prim, int *NbSeg) {
  int nb = *NbSeg;

  if (nb < 1)  // absurd! use the trio: target, min, max
  {
    double lWire = (XVertex[prim][1] - XVertex[prim][0]) *
                       (XVertex[prim][1] - XVertex[prim][0]) +
                   (YVertex[prim][1] - YVertex[prim][0]) *
                       (YVertex[prim][1] - YVertex[prim][0]) +
                   (ZVertex[prim][1] - ZVertex[prim][0]) *
                       (ZVertex[prim][1] - ZVertex[prim][0]);
    lWire = sqrt(lWire);

    nb = (int)(lWire / ElementLengthRqstd);

    if ((nb > MinNbElementsOnLength) &&
        (nb < MaxNbElementsOnLength)) {  // nothing to be done
    }
    // Check whether the length of the wire primitive is long enough
    else if (lWire < MINDIST) {
      nb = 1;
      fprintf(fMeshLog, "Wire element too small on primitive %d!\n", prim);
    }  // if lWire < MINDIST
    // Need to have at least MinNbElementsOnLength elements per wire primitive
    else if (nb < MinNbElementsOnLength) {
      nb = MinNbElementsOnLength;
      double ellength = lWire / (double)nb;
      if (ellength <
          MINDIST)  // which may not be possible if the length is small
      {
        nb = (int)(lWire / MINDIST);
        if (nb < 1)  // However, it is necessary to have at least one element!
        {
          nb = 1;
          fprintf(fMeshLog, "Wire element very small on primitive %d!\n", prim);
        }
      }
    }  // if nb < MinNbElementsOnLength
    else if (nb > MaxNbElementsOnLength) {
      nb = MaxNbElementsOnLength;
      fprintf(fMeshLog, "Too many elements on wire primitive %d!\n", prim);
      fprintf(fMeshLog, "Number of elements reduced to maximum allowed %d\n",
              MaxNbElementsOnLength);
    }  // if nb > MaxNbElementsOnLength

    *NbSeg = nb;

    fprintf(fMeshLog, "Number of elements on wire primitive %d is %d.\n\n",
            prim, *NbSeg);
  } else {  // number of dicretization specified by user
    double lWire = (XVertex[prim][1] - XVertex[prim][0]) *
                       (XVertex[prim][1] - XVertex[prim][0]) +
                   (YVertex[prim][1] - YVertex[prim][0]) *
                       (YVertex[prim][1] - YVertex[prim][0]) +
                   (ZVertex[prim][1] - ZVertex[prim][0]) *
                       (ZVertex[prim][1] - ZVertex[prim][0]);
    lWire = sqrt(lWire);

    double ellength = lWire / (double)nb;

    if (lWire < MINDIST)  // at least one element is a necessity
    {
      nb = 1;
      fprintf(fMeshLog, "Fatal: Wire element too small on primitive %d!\n",
              prim);
    }                             // if lWire < MINDIST
    else if (ellength < MINDIST)  // element length more than twice MINDIST
    {
      nb = (int)(lWire / (2.0 * MINDIST));
      if (nb < 1) {
        nb = 1;
        fprintf(fMeshLog, "Fatal: Wire element too small on primitive %d!\n",
                prim);
      }
    }  // if ellength < MINDIST

    *NbSeg = nb;

    fprintf(fMeshLog, "Number of elements on wire primitive %d is %d.\n\n",
            prim, *NbSeg);
  }

  if (nb)
    return 0;
  else
    return -1;
}  // AnalyzeWire ends

int AnalyzeSurface(int prim, int *NbSegCoord1, int *NbSegCoord2) {
  int nb1 = *NbSegCoord1, nb2 = *NbSegCoord2;

  if ((nb1 < 1) || (nb2 < 1))  // absurd! use the trio: target, min, max
  {
    // Triangle primitives have their right angle on vertex 1
    double l1 = (XVertex[prim][0] - XVertex[prim][1]) *
                    (XVertex[prim][0] - XVertex[prim][1]) +
                (YVertex[prim][0] - YVertex[prim][1]) *
                    (YVertex[prim][0] - YVertex[prim][1]) +
                (ZVertex[prim][0] - ZVertex[prim][1]) *
                    (ZVertex[prim][0] - ZVertex[prim][1]);
    l1 = sqrt(l1);
    double l2 = (XVertex[prim][2] - XVertex[prim][1]) *
                    (XVertex[prim][2] - XVertex[prim][1]) +
                (YVertex[prim][2] - YVertex[prim][1]) *
                    (YVertex[prim][2] - YVertex[prim][1]) +
                (ZVertex[prim][2] - ZVertex[prim][1]) *
                    (ZVertex[prim][2] - ZVertex[prim][1]);
    l2 = sqrt(l2);

    // We can use the lengths independently and forget about area
    // for the time being

    // double area = l1 * l2;

    nb1 = (int)(l1 / ElementLengthRqstd);
    if ((nb1 > MinNbElementsOnLength) && (nb1 < MaxNbElementsOnLength)) {
      // nothing to be done
    } else if (l1 < MINDIST) {
      fprintf(fMeshLog, "Length1 too small on primitive %d!\n", prim);
      nb1 = 1;
    } else if (nb1 < MinNbElementsOnLength) {
      // Need to have at least MinNbElementsOnLength elements per wire primitive
      nb1 = MinNbElementsOnLength;
      double ellength = l1 / (double)nb1;
      // which may not be possible if the length is small
      if (ellength < MINDIST) {
        nb1 = (int)(l1 / MINDIST);
        // However, it is necessary to have at least one element!
        if (nb1 < 1) {
          fprintf(fMeshLog, "Length1 very small on primitive %d!\n", prim);
          nb1 = 1;
        }
      }
    }  // else if nb1 < MinNbElementsOnLength

    if (nb1 > MaxNbElementsOnLength) {
      fprintf(fMeshLog, "Too many elements on Length1 for primitive %d!\n",
              prim);
      fprintf(fMeshLog, "Number of elements reduced to maximum allowed %d\n",
              MaxNbElementsOnLength);
      nb1 = MaxNbElementsOnLength;
    }

    nb2 = (int)(l2 / ElementLengthRqstd);
    if ((nb2 > MinNbElementsOnLength) && (nb2 < MaxNbElementsOnLength)) {
      // nothing to be done
    } else if (l2 < MINDIST) {
      fprintf(fMeshLog, "Length2 element too small on primitive %d!\n", prim);
      nb2 = 1;
    } else if (nb2 < MinNbElementsOnLength) {
      // Need to have at least MinNbElementsOnLength elements per wire primitive
      nb2 = MinNbElementsOnLength;
      double ellength = l2 / (double)nb2;
      // which may not be possible if the length is small
      if (ellength < MINDIST) {
        // However, it is necessary to have at least one element!
        nb2 = (int)(l2 / MINDIST);
        if (nb2 < 1) {
          fprintf(fMeshLog, "Length2 element very small on primitive %d!\n",
                  prim);
          nb2 = 1;
        }
      }
    }  // else if nb2 < MinNbElementsOnLength

    if (nb2 > MaxNbElementsOnLength) {
      fprintf(fMeshLog, "Too many elements on Length2 of primitive %d!\n",
              prim);
      fprintf(fMeshLog, "Number of elements reduced to maximum allowed %d\n",
              MaxNbElementsOnLength);
      nb2 = MaxNbElementsOnLength;
    }

    *NbSegCoord1 = nb1;
    *NbSegCoord2 = nb2;

    fprintf(fMeshLog,
            "Number of elements on surface primitive %d is %d X %d.\n\n", prim,
            *NbSegCoord1, *NbSegCoord2);
  } else {  // number of discretization specified by the user
    // Triangle primitives have their right angle on the vertex 1
    double l1 = (XVertex[prim][0] - XVertex[prim][1]) *
                    (XVertex[prim][0] - XVertex[prim][1]) +
                (YVertex[prim][0] - YVertex[prim][1]) *
                    (YVertex[prim][0] - YVertex[prim][1]) +
                (ZVertex[prim][0] - ZVertex[prim][1]) *
                    (ZVertex[prim][0] - ZVertex[prim][1]);
    l1 = sqrt(l1);
    double l2 = (XVertex[prim][2] - XVertex[prim][1]) *
                    (XVertex[prim][2] - XVertex[prim][1]) +
                (YVertex[prim][2] - YVertex[prim][1]) *
                    (YVertex[prim][2] - YVertex[prim][1]) +
                (ZVertex[prim][2] - ZVertex[prim][1]) *
                    (ZVertex[prim][2] - ZVertex[prim][1]);
    l2 = sqrt(l2);

    if (l1 > l2) {
      if (nb2 > nb1) {
        // swap numbers to allow larger number to larger side
        int tmpnb = nb1;
        nb1 = nb2;
        nb2 = tmpnb;
      }
    }  // if l1 > l2

    double ellength1 = l1 / (double)nb1;
    double ellength2 = l2 / (double)nb2;

    if (l1 < MINDIST) {
      nb1 = 1;
      fprintf(fMeshLog, "Fatal: Side length l1 too small! prim: %d\n", prim);
    } else if (ellength1 < MINDIST)  // element length more than twice MINDIST
    {
      nb1 = (int)(l1 / (2.0 * MINDIST));
      if (nb1 < 1) {
        nb1 = 1;
        fprintf(fMeshLog, "Fatal: Side length l1 too small on primitive %d!\n",
                prim);
      }
    }  // if ellength1 < MINDIST

    if (l2 < MINDIST) {
      nb2 = 1;
      fprintf(fMeshLog, "Fatal: Side length l2 too small! prim: %d\n", prim);
    } else if (ellength2 < MINDIST)  // element length more than twice MINDIST
    {
      nb2 = (int)(l2 / (2.0 * MINDIST));
      if (nb2 < 1) {
        nb2 = 1;
        fprintf(fMeshLog, "Fatal: Side length l2 too small on primitive %d!\n",
                prim);
      }
    }  // if ellength2 < MINDIST

    *NbSegCoord1 = nb1;
    *NbSegCoord2 = nb2;

    fprintf(fMeshLog,
            "Number of elements on surface primitive %d is %d X %d.\n\n", prim,
            *NbSegCoord1, *NbSegCoord2);
  }

  if ((nb1 > 0) && (nb2 > 0))
    return 0;
  else
    return -1;
}  // AnalyzeSurface ends

// Discretize wire into linear wire elements
int DiscretizeWire(int prim, int nvertex, double xvert[], double yvert[],
                   double zvert[], double radius, int volref1, int volref2,
                   int inttype, double potential, double charge, double lambda,
                   int NbSegs) {
  int WireParentObj, WireEType;
  double WireR, WireL;
  double WireLambda, WireV;
  double WireElX, WireElY, WireElZ, WireElL;
  DirnCosn3D PrimDirnCosn;  // direction cosine of the current primitive
  char primstr[10];
  char gpElem[256];
  FILE *fPrim, *fElem, *fgpPrim, *fgpElem;

  // Check inputs
  if (PrimType[prim] != 2) {
    neBEMMessage("DiscretizeWire - PrimType in DiscretizeWire");
    return -1;
  }
  if (nvertex != 2) {
    neBEMMessage("DiscretizeWire - nvertex in DiscretizeWire");
    return -1;
  }
  if (radius < MINDIST) {
    neBEMMessage("DiscretizeWire - radius in DiscretizeWire");
    return -1;
  }
  if (NbSegs <= 0) {
    neBEMMessage("DiscretizeWire - NbSegs in DiscretizeWire");
    return -1;
  }

  if (OptPrintVertexAndNormal) {
    printf("nvertex: %d\n", nvertex);
    for (int vert = 0; vert < nvertex; ++vert) {
      printf("vert: %d, x: %lg, y: %lg, z: %lg\n", vert, xvert[vert],
             yvert[vert], zvert[vert]);
    }
    printf("radius: %lg\n", radius);
  }  // if OptPrintVertexAndNormal

  // necessary for separating filenames
  sprintf(primstr, "%d", prim);

  // in order to avoid warning messages
  fPrim = NULL;
  fElem = NULL;
  fgpElem = NULL;

  WireParentObj = 1;  // ParentObj not being used now

  WireL = sqrt((xvert[1] - xvert[0]) * (xvert[1] - xvert[0])  // length of wire
               + (yvert[1] - yvert[0]) * (yvert[1] - yvert[0]) +
               (zvert[1] - zvert[0]) * (zvert[1] - zvert[0]));
  WireR = radius;
  WireElL = WireL / NbSegs;  // length of each wire element

  // Direction cosines along the wire - note difference from surface primitives!
  // The direction along the wire is considered to be the z axis of the LCS
  // So, let us fix that axial vector first
  PrimDirnCosn.ZUnit.X = (xvert[1] - xvert[0]) / WireL;  // useful
  PrimDirnCosn.ZUnit.Y = (yvert[1] - yvert[0]) / WireL;
  PrimDirnCosn.ZUnit.Z = (zvert[1] - zvert[0]) / WireL;  // useful
  // Next, let us find out the coefficients of a plane that passes through the
  // wire centroid and is normal to the axis of the wire. This is basically the
  // mid-plane of the cylindrical wire
  // Any vector OR on the plane normal to the axial direction satisfies
  // \vec{OR} . \vec{OA} = 0
  // where O is the wire centroid, A is a point on the axis and R is a point on
  // the cylindrical surface of the wire. \vec{OA} can be easily replaced by the
  // axial vector that is equivalent to the vector PrimDirnCosn.ZUnit
  // The equation of the plane can be shown to be:
  // XCoef * X + YCoef * Y + ZCoef * Z = Const
  // double XCoef, YCoef, ZCoef, Const; - not needed any more
  // double rnorm;
  // Point3D O, R;
  // Vector3D OR;
  // O = CreatePoint3D(WireX, WireY, WireZ);
  {
    Vector3D XUnit, YUnit, ZUnit;
    ZUnit.X = PrimDirnCosn.ZUnit.X;
    ZUnit.Y = PrimDirnCosn.ZUnit.Y;
    ZUnit.Z = PrimDirnCosn.ZUnit.Z;

    /* old code	- abs instead of fabs??!!
    XCoef = ZUnit.X;
    YCoef = ZUnit.Y;
    ZCoef = ZUnit.Z;
    double WireX = 0.5 * (xvert[1] + xvert[0]);
    double WireY = 0.5 * (yvert[1] + yvert[0]);
    double WireZ = 0.5 * (zvert[1] + zvert[0]);
    Const = WireX * ZUnit.X + WireY * ZUnit.Y + WireZ * ZUnit.Z;
    if(abs(XCoef) < 1.0e-12)	// X can be anything!
            {
            XUnit.X = 1.0;
            XUnit.Y = 0.0;
            XUnit.Z = 0.0;
            YUnit = Vector3DCrossProduct(ZUnit, XUnit);
            }
    else
            {
            // For a point on the above surface where both Y and Z are zero
            O = CreatePoint3D(WireX, WireY, WireZ);
            R = CreatePoint3D(Const, 0, 0);
            // Create the vector joining O and R; find X and Y unit vectors
            OR = CreateDistanceVector3D(O,R);
            XUnit = UnitVector3D(OR);
            YUnit = Vector3DCrossProduct(ZUnit, XUnit);
            }
    old code */

    // replaced following Rob's suggestions (used functions instead of direct
    // evaluation, although the latter is probably faster)
    // x-Axis: orthogonal in the 2 largest components.
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

    PrimDirnCosn.XUnit.X = XUnit.X;
    PrimDirnCosn.XUnit.Y = XUnit.Y;
    PrimDirnCosn.XUnit.Z = XUnit.Z;
    PrimDirnCosn.YUnit.X = YUnit.X;
    PrimDirnCosn.YUnit.Y = YUnit.Y;
    PrimDirnCosn.YUnit.Z = YUnit.Z;
  }  // X and Y direction cosines computed

  // primitive direction cosine assignments
  PrimDC[prim].XUnit.X = PrimDirnCosn.XUnit.X;
  PrimDC[prim].XUnit.Y = PrimDirnCosn.XUnit.Y;
  PrimDC[prim].XUnit.Z = PrimDirnCosn.XUnit.Z;
  PrimDC[prim].YUnit.X = PrimDirnCosn.YUnit.X;
  PrimDC[prim].YUnit.Y = PrimDirnCosn.YUnit.Y;
  PrimDC[prim].YUnit.Z = PrimDirnCosn.YUnit.Z;
  PrimDC[prim].ZUnit.X = PrimDirnCosn.ZUnit.X;
  PrimDC[prim].ZUnit.Y = PrimDirnCosn.ZUnit.Y;
  PrimDC[prim].ZUnit.Z = PrimDirnCosn.ZUnit.Z;

  // primitive origin: also the barycenter for a wire element
  PrimOriginX[prim] = 0.5 * (xvert[0] + xvert[1]);
  PrimOriginY[prim] = 0.5 * (yvert[0] + yvert[1]);
  PrimOriginZ[prim] = 0.5 * (zvert[0] + zvert[1]);
  PrimLX[prim] = WireR;  // radius for wire
  PrimLZ[prim] = WireL;  // length of wire

  WireEType = inttype;
  WireV = potential;
  WireLambda = lambda;

  // file output for a primitive
  if (OptPrimitiveFiles) {
    char OutPrim[256];
    strcpy(OutPrim, ModelOutDir);
    strcat(OutPrim, "/Primitives/Primitive");
    strcat(OutPrim, primstr);
    strcat(OutPrim, ".out");
    fPrim = fopen(OutPrim, "w");
    if (fPrim == NULL) {
      neBEMMessage("DiscretizeWire - OutPrim");
      return -1;
    }
    fprintf(fPrim, "#prim: %d, nvertex: %d\n", prim, nvertex);
    fprintf(fPrim, "Node1: %lg\t%lg\t%lg\n", xvert[0], yvert[0], zvert[0]);
    fprintf(fPrim, "Node2: %lg\t%lg\t%lg\n", xvert[1], yvert[1], zvert[1]);
    fprintf(fPrim, "PrimOrigin: %lg\t%lg\t%lg\n", PrimOriginX[prim],
            PrimOriginY[prim], PrimOriginZ[prim]);
    fprintf(fPrim, "Primitive lengths: %lg\t%lg\n", PrimLX[prim], PrimLZ[prim]);
    fprintf(fPrim, "#DirnCosn: \n");
    fprintf(fPrim, "%lg, %lg, %lg\n", PrimDirnCosn.XUnit.X,
            PrimDirnCosn.XUnit.Y, PrimDirnCosn.XUnit.Z);
    fprintf(fPrim, "%lg, %lg, %lg\n", PrimDirnCosn.YUnit.X,
            PrimDirnCosn.YUnit.Y, PrimDirnCosn.YUnit.Z);
    fprintf(fPrim, "%lg, %lg, %lg\n", PrimDirnCosn.ZUnit.X,
            PrimDirnCosn.ZUnit.Y, PrimDirnCosn.ZUnit.Z);
    fprintf(fPrim, "#volref1: %d, volref2: %d\n", volref1, volref2);
    fprintf(fPrim, "#NbSegs: %d\n", NbSegs);
    fprintf(fPrim, "#ParentObj: %d\tEType: %d\n", WireParentObj, WireEType);
    fprintf(fPrim, "#WireR: %lg\tWireL: %lg\n", WireR, WireL);
    fprintf(fPrim, "#SurfLambda: %lg\tSurfV: %lg\n", WireLambda, WireV);
  }

  // necessary for gnuplot
  if (OptGnuplot && OptGnuplotPrimitives) {
    char gpPrim[256];
    strcpy(gpPrim, MeshOutDir);
    strcat(gpPrim, "/GViewDir/gpPrim");
    strcat(gpPrim, primstr);
    strcat(gpPrim, ".out");
    fgpPrim = fopen(gpPrim, "w");
    if (fgpPrim == NULL) {
      neBEMMessage("DiscretizeWire - OutgpPrim");
      return -1;
    }
    fprintf(fgpPrim, "%g\t%g\t%g\n\n", xvert[0], yvert[0], zvert[0]);
    fprintf(fgpPrim, "%g\t%g\t%g\n\n", xvert[1], yvert[1], zvert[1]);
    fclose(fgpPrim);

    if (prim == 1)
      fprintf(fgnuPrim, " '%s\' w l", gpPrim);
    else
      fprintf(fgnuPrim, ", \\\n \'%s\' w l", gpPrim);
  }

  // file outputs for elements on primitive
  if (OptElementFiles) {
    char OutElem[256];
    strcpy(OutElem, MeshOutDir);
    strcat(OutElem, "/Elements/ElemOnPrim");
    strcat(OutElem, primstr);
    strcat(OutElem, ".out");
    fElem = fopen(OutElem, "w");
    if (fElem == NULL) {
      neBEMMessage("DiscretizeWire - OutElem");
      return -1;
    }
  }
  // gnuplot friendly file outputs for elements on primitive
  if (OptGnuplot && OptGnuplotElements) {
    strcpy(gpElem, MeshOutDir);
    strcat(gpElem, "/GViewDir/gpElemOnPrim");
    strcat(gpElem, primstr);
    strcat(gpElem, ".out");
    fgpElem = fopen(gpElem, "w");
    if (fgpElem == NULL) {
      neBEMMessage("DiscretizeWire - OutgpElem");
      if (fElem) fclose(fElem);
      return -1;
    }
  }

  double xincr = (xvert[1] - xvert[0]) / (double)NbSegs;
  double yincr = (yvert[1] - yvert[0]) / (double)NbSegs;
  double zincr = (zvert[1] - zvert[0]) / (double)NbSegs;

  ElementBgn[prim] = EleCntr + 1;
  double xv0, yv0, zv0, xv1, yv1, zv1;
  for (int seg = 1; seg <= NbSegs; ++seg) {
    xv0 = xvert[0] + ((double)seg - 1.0) * xincr;
    yv0 = yvert[0] + ((double)seg - 1.0) * yincr;
    zv0 = zvert[0] + ((double)seg - 1.0) * zincr;
    xv1 = xvert[0] + ((double)seg) * xincr;
    yv1 = yvert[0] + ((double)seg) * yincr;
    zv1 = zvert[0] + ((double)seg) * zincr;
    WireElX = xvert[0] + ((double)seg - 1.0) * xincr + 0.5 * xincr;
    WireElY = yvert[0] + ((double)seg - 1.0) * yincr + 0.5 * yincr;
    WireElZ = zvert[0] + ((double)seg - 1.0) * zincr + 0.5 * zincr;

    // Assign element values and write in the file
    // If element counter exceeds the maximum allowed number of elements, warn!
    ++EleCntr;
    if (EleCntr > NbElements) {
      neBEMMessage("DiscretizeWire - EleCntr more than NbElements!");
      if (fgpElem) fclose(fgpElem);
      if (fElem) fclose(fElem);
      return -1;
    }

    (EleArr + EleCntr - 1)->DeviceNb =
        1;  // At present, there is only one device
    (EleArr + EleCntr - 1)->ComponentNb = WireParentObj;
    (EleArr + EleCntr - 1)->PrimitiveNb = prim;
    (EleArr + EleCntr - 1)->Id = EleCntr;
    (EleArr + EleCntr - 1)->G.Type = 2;  // linear (wire) here
    (EleArr + EleCntr - 1)->G.Origin.X = WireElX;
    (EleArr + EleCntr - 1)->G.Origin.Y = WireElY;
    (EleArr + EleCntr - 1)->G.Origin.Z = WireElZ;
    (EleArr + EleCntr - 1)->G.Vertex[0].X = xv0;
    (EleArr + EleCntr - 1)->G.Vertex[0].Y = yv0;
    (EleArr + EleCntr - 1)->G.Vertex[0].Z = zv0;
    (EleArr + EleCntr - 1)->G.Vertex[1].X = xv1;
    (EleArr + EleCntr - 1)->G.Vertex[1].Y = yv1;
    (EleArr + EleCntr - 1)->G.Vertex[1].Z = zv1;
    (EleArr + EleCntr - 1)->G.LX = WireR;    // radius of the wire element
    (EleArr + EleCntr - 1)->G.LZ = WireElL;  // wire element length
    (EleArr + EleCntr - 1)->G.dA = 2.0 * MyPI * (EleArr + EleCntr - 1)->G.LX *
                                   (EleArr + EleCntr - 1)->G.LZ;
    (EleArr + EleCntr - 1)->G.DC.XUnit.X = PrimDirnCosn.XUnit.X;
    (EleArr + EleCntr - 1)->G.DC.XUnit.Y = PrimDirnCosn.XUnit.Y;
    (EleArr + EleCntr - 1)->G.DC.XUnit.Z = PrimDirnCosn.XUnit.Z;
    (EleArr + EleCntr - 1)->G.DC.YUnit.X = PrimDirnCosn.YUnit.X;
    (EleArr + EleCntr - 1)->G.DC.YUnit.Y = PrimDirnCosn.YUnit.Y;
    (EleArr + EleCntr - 1)->G.DC.YUnit.Z = PrimDirnCosn.YUnit.Z;
    (EleArr + EleCntr - 1)->G.DC.ZUnit.X = PrimDirnCosn.ZUnit.X;
    (EleArr + EleCntr - 1)->G.DC.ZUnit.Y = PrimDirnCosn.ZUnit.Y;
    (EleArr + EleCntr - 1)->G.DC.ZUnit.Z = PrimDirnCosn.ZUnit.Z;
    (EleArr + EleCntr - 1)->E.Type = WireEType;
    (EleArr + EleCntr - 1)->E.Lambda = WireLambda;
    (EleArr + EleCntr - 1)->Solution = 0.0;
    (EleArr + EleCntr - 1)->Assigned = charge;
    (EleArr + EleCntr - 1)->BC.NbOfBCs = 1;  // assume one BC per element
    (EleArr + EleCntr - 1)->BC.CollPt.X =
        (EleArr + EleCntr - 1)->G.Origin.X;  // modify
    (EleArr + EleCntr - 1)->BC.CollPt.Y =
        (EleArr + EleCntr - 1)->G.Origin.Y;  // to be on
    (EleArr + EleCntr - 1)->BC.CollPt.Z =
        (EleArr + EleCntr - 1)->G.Origin.Z;  // surface?

    // File operations begin
    // rfw = fwrite(&Ele, sizeof(Element), 1, fpEle);
    // printf("Return of fwrite is %d\n", rfw);

    if (OptElementFiles) {
      fprintf(fElem, "##Element Counter: %d\n", EleCntr);
      fprintf(fElem, "#DevNb\tCompNb\tPrimNb\tId\n");
      fprintf(fElem, "%d\t%d\t%d\t%d\n", (EleArr + EleCntr - 1)->DeviceNb,
              (EleArr + EleCntr - 1)->ComponentNb,
              (EleArr + EleCntr - 1)->PrimitiveNb, (EleArr + EleCntr - 1)->Id);
      fprintf(fElem, "#GType\tX\tY\tZ\tLX\tLZ\tdA\n");
      fprintf(fElem, "%d\t%.16lg\t%.16lg\t%.16lg\t%.16lg\t%.16lg\t%.16lg\n",
              (EleArr + EleCntr - 1)->G.Type,
              (EleArr + EleCntr - 1)->G.Origin.X,
              (EleArr + EleCntr - 1)->G.Origin.Y,
              (EleArr + EleCntr - 1)->G.Origin.Z, (EleArr + EleCntr - 1)->G.LX,
              (EleArr + EleCntr - 1)->G.LZ, (EleArr + EleCntr - 1)->G.dA);
      fprintf(fElem, "#DirnCosn: \n");
      fprintf(fElem, "%lg, %lg, %lg\n", (EleArr + EleCntr - 1)->G.DC.XUnit.X,
              (EleArr + EleCntr - 1)->G.DC.XUnit.Y,
              (EleArr + EleCntr - 1)->G.DC.XUnit.Z);
      fprintf(fElem, "%lg, %lg, %lg\n", (EleArr + EleCntr - 1)->G.DC.YUnit.X,
              (EleArr + EleCntr - 1)->G.DC.YUnit.Y,
              (EleArr + EleCntr - 1)->G.DC.YUnit.Z);
      fprintf(fElem, "%lg, %lg, %lg\n", (EleArr + EleCntr - 1)->G.DC.ZUnit.X,
              (EleArr + EleCntr - 1)->G.DC.ZUnit.Y,
              (EleArr + EleCntr - 1)->G.DC.ZUnit.Z);
      fprintf(fElem, "#EType\tLambda\n");
      fprintf(fElem, "%d\t%lg\n", (EleArr + EleCntr - 1)->E.Type,
              (EleArr + EleCntr - 1)->E.Lambda);
      fprintf(fElem, "#NbBCs\tCPX\tCPY\tCPZ\tValue\n");
      fprintf(fElem, "%d\t%.16lg\t%.16lg\t%.16lg\t%lg\n",
              (EleArr + EleCntr - 1)->BC.NbOfBCs,
              (EleArr + EleCntr - 1)->BC.CollPt.X,
              (EleArr + EleCntr - 1)->BC.CollPt.Y,
              (EleArr + EleCntr - 1)->BC.CollPt.Z,
              (EleArr + EleCntr - 1)->BC.Value);
    }  // if OptElementFiles
       //
    // mark centroid
    if (OptGnuplot && OptGnuplotElements) {
      fprintf(fgpElem, "%g\t%g\t%g\n", (EleArr + EleCntr - 1)->BC.CollPt.X,
              (EleArr + EleCntr - 1)->BC.CollPt.Y,
              (EleArr + EleCntr - 1)->BC.CollPt.Z);
    }  // if OptElementFiles
       // File operations end
  }    // seg loop for wire elements
  ElementEnd[prim] = EleCntr;
  NbElmntsOnPrim[prim] = ElementEnd[prim] - ElementBgn[prim] + 1;

  if (OptPrimitiveFiles) {
    fprintf(fPrim, "Element begin: %d, Element end: %d\n", ElementBgn[prim],
            ElementEnd[prim]);
    fprintf(fPrim, "Number of elements on primitive: %d\n",
            NbElmntsOnPrim[prim]);
    fclose(fPrim);
  }

  if (OptElementFiles) fclose(fElem);
  if (OptGnuplot && OptGnuplotElements) fclose(fgpElem);

  return (0);
}  // end of DiscretizeWire

// NbSegX is considered to be the number of rectangular + one triangular
// elements into which the base of the triangular primitive is divided.
// NbSegZ is the total number of elements into which the height is divided.
int DiscretizeTriangle(int prim, int nvertex, double xvert[], double yvert[],
                       double zvert[], double xnorm, double ynorm, double znorm,
                       int volref1, int volref2, int inttype, double potential,
                       double charge, double lambda, int NbSegX, int NbSegZ) {
  int SurfParentObj, SurfEType;
  double SurfX, SurfY, SurfZ, SurfLX, SurfLZ;
  double SurfElX, SurfElY, SurfElZ, SurfElLX, SurfElLZ;
  DirnCosn3D PrimDirnCosn;  // direction cosine of the current primitive
  char primstr[10];
  char gpElem[256], gpMesh[256];
  FILE *fPrim, *fElem, *fgpPrim, *fgpElem, *fgpMesh;

  // Check inputs
  if ((NbSegX <= 0) || (NbSegZ <= 0)) {
    printf("segmentation input wrong in DiscretizeTriangle ...\n");
    return -1;
  }

  if (OptPrintVertexAndNormal) {
    printf("nvertex: %d\n", nvertex);
    for (int vert = 0; vert < nvertex; ++vert) {
      printf("vert: %d, x: %lg, y: %lg, z: %lg\n", vert, xvert[vert],
             yvert[vert], zvert[vert]);
    }
    printf("Normal: %lg, %lg, %lg\n", xnorm, ynorm, znorm);
  }  // if OptPrintVertexAndNormal

  // necessary for separating filenames
  sprintf(primstr, "%d", prim);

  // in order to avoid warning messages
  fPrim = NULL;
  fElem = NULL;
  fgpMesh = NULL;
  fgpElem = NULL;

  // Compute all the properties of this surface
  // Boundary types from 1 to 7 have been defined
  SurfParentObj = 1;
  SurfEType = inttype;
  if ((SurfEType <= 0) || (SurfEType >= 8)) {
    printf("Wrong SurfEType for prim %d\n", prim);
    exit(-1);
  }
  // Origin of the local coordinate center is at the right angle corner
  SurfX = xvert[1];
  SurfY = yvert[1];
  SurfZ = zvert[1];

  // Find the proper direction cosines first - little more tricky that in the
  // rectangular case
  int flagDC = 0;  // DC has not been ascertained as yet
  // We begin with trial 1: one of the possible orientations
  // Intially, the lengths are necessary
  // lengths of the sides - note that right angle corner is [1]-th element
  SurfLX = sqrt((xvert[0] - xvert[1]) * (xvert[0] - xvert[1]) +
                (yvert[0] - yvert[1]) * (yvert[0] - yvert[1]) +
                (zvert[0] - zvert[1]) * (zvert[0] - zvert[1]));
  SurfLZ = sqrt((xvert[2] - xvert[1]) * (xvert[2] - xvert[1]) +
                (yvert[2] - yvert[1]) * (yvert[2] - yvert[1]) +
                (zvert[2] - zvert[1]) * (zvert[2] - zvert[1]));
  // Direction cosines - note that right angle corner is [1]-th element
  PrimDirnCosn.XUnit.X = (xvert[0] - xvert[1]) / SurfLX;
  PrimDirnCosn.XUnit.Y = (yvert[0] - yvert[1]) / SurfLX;
  PrimDirnCosn.XUnit.Z = (zvert[0] - zvert[1]) / SurfLX;
  PrimDirnCosn.ZUnit.X = (xvert[2] - xvert[1]) / SurfLZ;
  PrimDirnCosn.ZUnit.Y = (yvert[2] - yvert[1]) / SurfLZ;
  PrimDirnCosn.ZUnit.Z = (zvert[2] - zvert[1]) / SurfLZ;
  PrimDirnCosn.YUnit =
      Vector3DCrossProduct(&PrimDirnCosn.ZUnit, &PrimDirnCosn.XUnit);
  if ((fabs(PrimDirnCosn.YUnit.X - xnorm) <= 1.0e-3) &&
      (fabs(PrimDirnCosn.YUnit.Y - ynorm) <= 1.0e-3) &&
      (fabs(PrimDirnCosn.YUnit.Z - znorm) <= 1.0e-3))
    flagDC = 1;
  if (DebugLevel == 202) {
    printf("First attempt: \n");
    PrintDirnCosn3D(PrimDirnCosn);
    printf("\n");
  }

  if (!flagDC)  // if DC search is unsuccessful, try the other orientation
  {
    SurfLX = sqrt((xvert[2] - xvert[1]) * (xvert[2] - xvert[1]) +
                  (yvert[2] - yvert[1]) * (yvert[2] - yvert[1]) +
                  (zvert[2] - zvert[1]) * (zvert[2] - zvert[1]));
    SurfLZ = sqrt((xvert[0] - xvert[1]) * (xvert[0] - xvert[1]) +
                  (yvert[0] - yvert[1]) * (yvert[0] - yvert[1]) +
                  (zvert[0] - zvert[1]) * (zvert[0] - zvert[1]));
    // Direction cosines - note that right angle corner is [1]-th element
    PrimDirnCosn.XUnit.X = (xvert[2] - xvert[1]) / SurfLX;
    PrimDirnCosn.XUnit.Y = (yvert[2] - yvert[1]) / SurfLX;
    PrimDirnCosn.XUnit.Z = (zvert[2] - zvert[1]) / SurfLX;
    PrimDirnCosn.ZUnit.X = (xvert[0] - xvert[1]) / SurfLZ;
    PrimDirnCosn.ZUnit.Y = (yvert[0] - yvert[1]) / SurfLZ;
    PrimDirnCosn.ZUnit.Z = (zvert[0] - zvert[1]) / SurfLZ;
    PrimDirnCosn.YUnit =
        Vector3DCrossProduct(&PrimDirnCosn.ZUnit, &PrimDirnCosn.XUnit);
    if ((fabs(PrimDirnCosn.YUnit.X - xnorm) <= 1.0e-3) &&
        (fabs(PrimDirnCosn.YUnit.Y - ynorm) <= 1.0e-3) &&
        (fabs(PrimDirnCosn.YUnit.Z - znorm) <= 1.0e-3))
      flagDC = 2;
    if (DebugLevel == 202) {
      printf("Second attempt: \n");
      PrintDirnCosn3D(PrimDirnCosn);
      printf("\n");
    }
  }

  if (!flagDC)  // No other possibility, DC search failed!!!
  {
    printf("Triangle DC problem ... returning ...\n");
    // getchar();
    return -1;
  }

  // primitive direction cosine assignments
  PrimDC[prim].XUnit.X = PrimDirnCosn.XUnit.X;
  PrimDC[prim].XUnit.Y = PrimDirnCosn.XUnit.Y;
  PrimDC[prim].XUnit.Z = PrimDirnCosn.XUnit.Z;
  PrimDC[prim].YUnit.X = PrimDirnCosn.YUnit.X;
  PrimDC[prim].YUnit.Y = PrimDirnCosn.YUnit.Y;
  PrimDC[prim].YUnit.Z = PrimDirnCosn.YUnit.Z;
  PrimDC[prim].ZUnit.X = PrimDirnCosn.ZUnit.X;
  PrimDC[prim].ZUnit.Y = PrimDirnCosn.ZUnit.Y;
  PrimDC[prim].ZUnit.Z = PrimDirnCosn.ZUnit.Z;

  // primitive origin - for a triangle, origin is at the right corner
  PrimOriginX[prim] = SurfX;
  PrimOriginY[prim] = SurfY;
  PrimOriginZ[prim] = SurfZ;
  PrimLX[prim] = SurfLX;
  PrimLZ[prim] = SurfLZ;
  if (flagDC == 1) {
    int tmpVolRef1 = VolRef1[prim];
    VolRef1[prim] = VolRef2[prim];
    VolRef2[prim] = tmpVolRef1;
    int tmpEpsilon1 = Epsilon1[prim];
    Epsilon1[prim] = Epsilon2[prim];
    Epsilon2[prim] = tmpEpsilon1;
  }

  double SurfLambda = lambda;
  double SurfV = potential;

  // file output for a primitive
  if (OptPrimitiveFiles) {
    char OutPrim[256];
    strcpy(OutPrim, ModelOutDir);
    strcat(OutPrim, "/Primitives/Primitive");
    strcat(OutPrim, primstr);
    strcat(OutPrim, ".out");
    fPrim = fopen(OutPrim, "w");
    if (fPrim == NULL) {
      neBEMMessage("DiscretizeTriangle - OutPrim");
      return -1;
    }
    fprintf(fPrim, "#prim: %d, nvertex: %d\n", prim, nvertex);
    fprintf(fPrim, "Node1: %lg\t%lg\t%lg\n", xvert[0], yvert[0], zvert[0]);
    fprintf(fPrim, "Node2: %lg\t%lg\t%lg\n", xvert[1], yvert[1], zvert[1]);
    fprintf(fPrim, "Node3: %lg\t%lg\t%lg\n", xvert[2], yvert[2], zvert[2]);
    fprintf(fPrim, "PrimOrigin: %lg\t%lg\t%lg\n", PrimOriginX[prim],
            PrimOriginY[prim], PrimOriginZ[prim]);
    fprintf(fPrim, "Primitive lengths: %lg\t%lg\n", PrimLX[prim], PrimLZ[prim]);
    fprintf(fPrim, "Norm: %lg\t%lg\t%lg\n", xnorm, ynorm, znorm);
    fprintf(fPrim, "#volref1: %d, volref2: %d\n", volref1, volref2);
    fprintf(fPrim, "#NbSegX: %d, NbSegZ: %d (check note!)\n", NbSegX, NbSegZ);
    fprintf(fPrim, "#ParentObj: %d\tEType: %d\n", SurfParentObj, SurfEType);
    fprintf(fPrim, "#SurfX\tSurfY\tSurfZ\tSurfLX\tSurfLZ (Rt. Corner)\n");
    fprintf(fPrim, "%lg\t%lg\t%lg\t%lg\t%lg\n", SurfX, SurfY, SurfZ, SurfLX,
            SurfLZ);
    fprintf(fPrim, "#DirnCosn: \n");
    fprintf(fPrim, "%lg, %lg, %lg\n", PrimDirnCosn.XUnit.X,
            PrimDirnCosn.XUnit.Y, PrimDirnCosn.XUnit.Z);
    fprintf(fPrim, "%lg, %lg, %lg\n", PrimDirnCosn.YUnit.X,
            PrimDirnCosn.YUnit.Y, PrimDirnCosn.YUnit.Z);
    fprintf(fPrim, "%lg, %lg, %lg\n", PrimDirnCosn.ZUnit.X,
            PrimDirnCosn.ZUnit.Y, PrimDirnCosn.ZUnit.Z);
    fprintf(fPrim, "#SurfLambda: %lg\tSurfV: %lg\n", SurfLambda, SurfV);
  }

  // necessary for gnuplot
  if (OptGnuplot && OptGnuplotPrimitives) {
    char gpPrim[256];
    strcpy(gpPrim, MeshOutDir);
    strcat(gpPrim, "/GViewDir/gpPrim");
    strcat(gpPrim, primstr);
    strcat(gpPrim, ".out");
    fgpPrim = fopen(gpPrim, "w");
    if (fgpPrim == NULL) {
      neBEMMessage("DiscretizeTriangle - OutgpPrim");
      return -1;
    }
    fprintf(fgpPrim, "%g\t%g\t%g\n\n", xvert[0], yvert[0], zvert[0]);
    fprintf(fgpPrim, "%g\t%g\t%g\n\n", xvert[1], yvert[1], zvert[1]);
    fprintf(fgpPrim, "%g\t%g\t%g\n\n", xvert[2], yvert[2], zvert[2]);
    fprintf(fgpPrim, "%g\t%g\t%g\n", xvert[0], yvert[0], zvert[0]);
    fclose(fgpPrim);

    if (prim == 1)
      fprintf(fgnuPrim, " '%s\' w l", gpPrim);
    else
      fprintf(fgnuPrim, ", \\\n \'%s\' w l", gpPrim);
  }

  // file outputs for elements on primitive
  if (OptElementFiles) {
    char OutElem[256];
    strcpy(OutElem, MeshOutDir);
    strcat(OutElem, "/Elements/ElemOnPrim");
    strcat(OutElem, primstr);
    strcat(OutElem, ".out");
    fElem = fopen(OutElem, "w");
    if (fElem == NULL) {
      neBEMMessage("DiscretizeTriangle - OutElem");
      return -1;
    }
  }
  // gnuplot friendly file outputs for elements on primitive
  if (OptGnuplot && OptGnuplotElements) {
    strcpy(gpElem, MeshOutDir);
    strcat(gpElem, "/GViewDir/gpElemOnPrim");
    strcat(gpElem, primstr);
    strcat(gpElem, ".out");
    fgpElem = fopen(gpElem, "w");
    // assert(fgpElem != NULL);
    if (fgpElem == NULL) {
      neBEMMessage("DiscretizeTriangle - OutgpElem");
      if (fElem) fclose(fElem);
      return -1;
    }
    // gnuplot friendly file outputs for elements on primitive
    strcpy(gpMesh, MeshOutDir);
    strcat(gpMesh, "/GViewDir/gpMeshOnPrim");
    strcat(gpMesh, primstr);
    strcat(gpMesh, ".out");
    fgpMesh = fopen(gpMesh, "w");
    if (fgpMesh == NULL) {
      neBEMMessage("DiscretizeTriangle - OutgpMesh");
      fclose(fgpElem);
      if (fElem) fclose(fElem);
      return -1;
    }
  }

  // Compute element positions (CGs) in primitive local coordinate system (PCS).
  // Then map these CGs to the global coordinate system.
  // (xav, 0, zav) is the CG of an element wrt the primitive coordinate system
  // From this, we find the offset of the element from the primitive CG in the
  // global coordinate system by just rotating the displacement vector
  // (0,0,0) -> (xav,0,zav) using the known DCM of the surface.
  // This rotated vector (now in the global system) is then added to the
  // position vector of the primitive CG (always in the global system) to get
  // the CG of the element in the global system.

  // Consult note for clarifications on the discretization algorithm.
  // SurfElLX is true for the lowest row  of elements, but indicative of the
  // other rows
  // Make sure that the larger number is being used to discretize the longer
  // side - can be OverSmart in certain cases - there may be cases where
  // smaller number of elements suffice for a longer side
  if (NbSegX == NbSegZ) {
    SurfElLX = SurfLX / NbSegX;  // element sizes
    SurfElLZ = SurfLZ / NbSegZ;
  } else if (NbSegX > NbSegZ) {
    if (SurfLX > SurfLZ) {
      SurfElLX = SurfLX / NbSegX;  // element sizes
      SurfElLZ = SurfLZ / NbSegZ;
    } else  // interchange NbSegX and NbSegZ
    {
      int tmp = NbSegZ;
      NbSegZ = NbSegX;
      NbSegX = tmp;
      SurfElLX = SurfLX / NbSegX;  // element sizes
      SurfElLZ = SurfLZ / NbSegZ;
    }
  }     // NbSegX > NbSegZ
  else  // NbSegX < NbSegZ
  {
    if (SurfLX < SurfLZ) {
      SurfElLX = SurfLX / NbSegX;  // element sizes
      SurfElLZ = SurfLZ / NbSegZ;
    } else  // interchange NbSegX and NbSegZ
    {
      int tmp = NbSegZ;
      NbSegZ = NbSegX;
      NbSegX = tmp;
      SurfElLX = SurfLX / NbSegX;  // element sizes
      SurfElLZ = SurfLZ / NbSegZ;
    }
  }

  // Analysis of element aspect ratio.
  // Note that we can afford only to reduce the total number of elements.
  // Otherwise, we'll have to realloc `EleArr' array.
  // Later, we'll make the corrections such that the total number of elements
  // remain close to the originally intended value.
  double AR = SurfElLX / SurfElLZ;  // indicative element aspect ratio
  if (OptPrimitiveFiles) {
    fprintf(fPrim,
            "Using the input, the aspect ratio of the elements on prim: %d\n",
            prim);
    fprintf(fPrim,
            "NbSegX: %d, SurfElLX: %lg, NbSegZ: %d, SurfElLZ: %lg, AR: %lg\n",
            NbSegX, SurfElLX, NbSegZ, SurfElLZ, AR);
  }
  if (AR > 10.0)  // NbSegZ needs to be reduced
  {
    double tmpElLZ = SurfElLX / 10.0;
    NbSegZ = (int)(SurfLZ / tmpElLZ);
    if (NbSegZ <= 0) NbSegZ = 1;
    SurfElLZ = SurfLZ / NbSegZ;
  }
  if (AR < 0.1)  // NbSegX need to be reduced
  {
    double tmpElLX = SurfElLZ * 0.1;
    NbSegX = (int)(SurfLX / tmpElLX);
    if (NbSegX <= 0) NbSegX = 1;
    SurfElLX = SurfLX / NbSegX;
  }
  AR = SurfElLX / SurfElLZ;  // indicative element aspect ratio
  if (OptPrimitiveFiles) {
    fprintf(fPrim, "After analyzing the likely aspect ratio of the elements\n");
    fprintf(fPrim,
            "NbSegX: %d, SurfElLX: %lg, NbSegZ: %d, SurfElLZ: %lg, AR: %lg\n",
            NbSegX, SurfElLX, NbSegZ, SurfElLZ, AR);
  }

  // The top-most right angle triangle is a lone element on the highest row, its
  // properties being determined irrespective of the others. Despite that, this
  // element is being treated as one among the others, especially to facilitate
  // element indexing, file output etc.
  // Each row, thus, has at least one triangle, and possibly several rectangular
  // elements. All the elements in any given row has the same SurfElLZ as has
  // been determined above.
  // First we create the triangular element and then move on to create the
  // rectangular ones.
  double xv0, yv0, zv0, xv1, yv1, zv1, xv2, yv2, zv2;
  ElementBgn[prim] = EleCntr + 1;
  for (int k = 1; k <= NbSegZ; ++k)  // consider the k-th row
  {
    double grad = (SurfLZ / SurfLX);
    double zlopt = (k - 1) * SurfElLZ;
    double zhipt = (k)*SurfElLZ;
    double xlopt = (SurfLZ - zlopt) / grad;  // used again on 21 Feb 2014
    double xhipt = (SurfLZ - zhipt) / grad;

    // the triangular element on the k-th row can now be specified in PCS
    double xtorigin = xhipt;
    double ytorigin = 0.0;
    double ztorigin = zlopt;
    {  // Separate block for position rotation - local2global
      Point3D localDisp, globalDisp;

      localDisp.X = xtorigin;
      localDisp.Y = ytorigin;
      localDisp.Z = ztorigin;
      globalDisp = RotatePoint3D(&localDisp, &PrimDirnCosn, local2global);
      SurfElX =
          SurfX + globalDisp.X;  // these are the coords in GCS of origin of
      SurfElY =
          SurfY + globalDisp.Y;  // the triangluar element under consideration
      SurfElZ = SurfZ + globalDisp.Z;
    }  // vector rotation over

    // Assign element values and write in the file
    ++EleCntr;
    if (EleCntr > NbElements) {
      neBEMMessage("DiscretizeTriangle - EleCntr more than NbElements 1!");
      if (fgpElem) fclose(fgpElem);
      if (fgpMesh) fclose(fgpMesh);
      return -1;
    }

    (EleArr + EleCntr - 1)->DeviceNb =
        1;  // At present, there is only one device
    (EleArr + EleCntr - 1)->ComponentNb = SurfParentObj;
    (EleArr + EleCntr - 1)->PrimitiveNb = prim;
    (EleArr + EleCntr - 1)->Id = EleCntr;
    (EleArr + EleCntr - 1)->G.Type = 3;  // triangular here
    (EleArr + EleCntr - 1)->G.Origin.X = SurfElX;
    (EleArr + EleCntr - 1)->G.Origin.Y = SurfElY;
    (EleArr + EleCntr - 1)->G.Origin.Z = SurfElZ;
    // (EleArr+EleCntr-1)->G.LX = SurfElLX;	// previously written as xlopt -
    // xhipt;
    (EleArr + EleCntr - 1)->G.LX =
        xlopt - xhipt;  // back to old ways on 21 Feb 2014
    (EleArr + EleCntr - 1)->G.LZ = SurfElLZ;
    (EleArr + EleCntr - 1)->G.LZ =
        zhipt - zlopt;  // to be on the safe side, 21/2/14
    (EleArr + EleCntr - 1)->G.dA =
        0.5 * (EleArr + EleCntr - 1)->G.LX * (EleArr + EleCntr - 1)->G.LZ;
    // Safe to use the direction cosines obtained for the triangular primitive
    // since they are bound to remain unchanged for the rectangular sub-elements
    (EleArr + EleCntr - 1)->G.DC.XUnit.X = PrimDirnCosn.XUnit.X;
    (EleArr + EleCntr - 1)->G.DC.XUnit.Y = PrimDirnCosn.XUnit.Y;
    (EleArr + EleCntr - 1)->G.DC.XUnit.Z = PrimDirnCosn.XUnit.Z;
    (EleArr + EleCntr - 1)->G.DC.YUnit.X = PrimDirnCosn.YUnit.X;
    (EleArr + EleCntr - 1)->G.DC.YUnit.Y = PrimDirnCosn.YUnit.Y;
    (EleArr + EleCntr - 1)->G.DC.YUnit.Z = PrimDirnCosn.YUnit.Z;
    (EleArr + EleCntr - 1)->G.DC.ZUnit.X = PrimDirnCosn.ZUnit.X;
    (EleArr + EleCntr - 1)->G.DC.ZUnit.Y = PrimDirnCosn.ZUnit.Y;
    (EleArr + EleCntr - 1)->G.DC.ZUnit.Z = PrimDirnCosn.ZUnit.Z;
    (EleArr + EleCntr - 1)->E.Type = SurfEType;
    (EleArr + EleCntr - 1)->E.Lambda = SurfLambda;
    (EleArr + EleCntr - 1)->Solution = 0.0;
    (EleArr + EleCntr - 1)->Assigned = charge;
    // Boundary condition is applied at the barycenter, not at the origin
    // of the element coordinate system (ECS) which is at the right corner
    (EleArr + EleCntr - 1)->BC.NbOfBCs = 1;  // assume one BC per element

    xv0 = (EleArr + EleCntr - 1)->G.Origin.X;
    yv0 = (EleArr + EleCntr - 1)->G.Origin.Y;
    zv0 = (EleArr + EleCntr - 1)->G.Origin.Z;
    xv1 = (EleArr + EleCntr - 1)->G.Origin.X +
          (EleArr + EleCntr - 1)->G.LX * (EleArr + EleCntr - 1)->G.DC.XUnit.X;
    yv1 = (EleArr + EleCntr - 1)->G.Origin.Y +
          (EleArr + EleCntr - 1)->G.LX * (EleArr + EleCntr - 1)->G.DC.XUnit.Y;
    zv1 = (EleArr + EleCntr - 1)->G.Origin.Z +
          (EleArr + EleCntr - 1)->G.LX * (EleArr + EleCntr - 1)->G.DC.XUnit.Z;
    xv2 = (EleArr + EleCntr - 1)->G.Origin.X +
          (EleArr + EleCntr - 1)->G.LZ * (EleArr + EleCntr - 1)->G.DC.ZUnit.X;
    yv2 = (EleArr + EleCntr - 1)->G.Origin.Y +
          (EleArr + EleCntr - 1)->G.LZ * (EleArr + EleCntr - 1)->G.DC.ZUnit.Y;
    zv2 = (EleArr + EleCntr - 1)->G.Origin.Z +
          (EleArr + EleCntr - 1)->G.LZ * (EleArr + EleCntr - 1)->G.DC.ZUnit.Z;
    // assign vertices of the element
    (EleArr + EleCntr - 1)->G.Vertex[0].X = xv0;
    (EleArr + EleCntr - 1)->G.Vertex[0].Y = yv0;
    (EleArr + EleCntr - 1)->G.Vertex[0].Z = zv0;
    (EleArr + EleCntr - 1)->G.Vertex[1].X = xv1;
    (EleArr + EleCntr - 1)->G.Vertex[1].Y = yv1;
    (EleArr + EleCntr - 1)->G.Vertex[1].Z = zv1;
    (EleArr + EleCntr - 1)->G.Vertex[2].X = xv2;
    (EleArr + EleCntr - 1)->G.Vertex[2].Y = yv2;
    (EleArr + EleCntr - 1)->G.Vertex[2].Z = zv2;
    (EleArr + EleCntr - 1)->G.Vertex[3].X = 0.0;
    (EleArr + EleCntr - 1)->G.Vertex[3].Y = 0.0;
    (EleArr + EleCntr - 1)->G.Vertex[3].Z = 0.0;

    if (DebugLevel == 201) {
      printf("Primitive nb: %d\n", (EleArr + EleCntr - 1)->PrimitiveNb);
      printf("Element id: %d\n", (EleArr + EleCntr - 1)->Id);
      printf("Element X, Y, Z: %lg %lg %lg\n",
             (EleArr + EleCntr - 1)->G.Origin.X,
             (EleArr + EleCntr - 1)->G.Origin.Y,
             (EleArr + EleCntr - 1)->G.Origin.Z);
      printf("Element LX, LZ: %lg %lg\n", (EleArr + EleCntr - 1)->G.LX,
             (EleArr + EleCntr - 1)->G.LZ);
      printf("Element (primitive) X axis dirn cosines: %lg, %lg, %lg\n",
             PrimDirnCosn.XUnit.X, PrimDirnCosn.XUnit.Y, PrimDirnCosn.XUnit.Z);
      printf("Element (primitive) Y axis dirn cosines: %lg, %lg, %lg\n",
             PrimDirnCosn.YUnit.X, PrimDirnCosn.YUnit.Y, PrimDirnCosn.YUnit.Z);
      printf("Element (primitive) Z axis dirn cosines: %lg, %lg, %lg\n",
             PrimDirnCosn.ZUnit.X, PrimDirnCosn.ZUnit.Y, PrimDirnCosn.ZUnit.Z);
    }
    // Following are the location in the ECS
    double dxl = (EleArr + EleCntr - 1)->G.LX / 3.0;
    double dyl = 0.0;
    double dzl = (EleArr + EleCntr - 1)->G.LZ / 3.0;
    {  // Separate block for position rotation - local2global
      Point3D localDisp, globalDisp;

      localDisp.X = dxl;
      localDisp.Y = dyl;
      localDisp.Z = dzl;
      if (DebugLevel == 201) {
        printf("Element dxl, dxy, dxz: %lg %lg %lg\n", localDisp.X, localDisp.Y,
               localDisp.Z);
      }

      globalDisp = RotatePoint3D(&localDisp, &PrimDirnCosn, local2global);
      (EleArr + EleCntr - 1)->BC.CollPt.X =
          (EleArr + EleCntr - 1)->G.Origin.X + globalDisp.X;
      (EleArr + EleCntr - 1)->BC.CollPt.Y =
          (EleArr + EleCntr - 1)->G.Origin.Y + globalDisp.Y;
      (EleArr + EleCntr - 1)->BC.CollPt.Z =
          (EleArr + EleCntr - 1)->G.Origin.Z + globalDisp.Z;
      if (DebugLevel == 201) {
        printf("Element global dxl, dxy, dxz: %lg %lg %lg\n", globalDisp.X,
               globalDisp.Y, globalDisp.Z);
        printf("Element BCX, BCY, BCZ: %lg %lg %lg\n",
               (EleArr + EleCntr - 1)->BC.CollPt.X,
               (EleArr + EleCntr - 1)->BC.CollPt.Y,
               (EleArr + EleCntr - 1)->BC.CollPt.Z);
      }
    }  // vector rotation over
    // (EleArr+EleCntr-1)->BC.Value = SurfV; // assigned in BoundaryConditions

    if (OptElementFiles) {
      fprintf(fElem, "##Element Counter: %d\n", EleCntr);
      fprintf(fElem, "#DevNb\tCompNb\tPrimNb\tId\n");
      fprintf(fElem, "%d\t%d\t%d\t%d\n", (EleArr + EleCntr - 1)->DeviceNb,
              (EleArr + EleCntr - 1)->ComponentNb,
              (EleArr + EleCntr - 1)->PrimitiveNb, (EleArr + EleCntr - 1)->Id);
      fprintf(fElem, "#GType\tX\tY\tZ\tLX\tLZ\tdA\n");
      fprintf(fElem, "%d\t%.16lg\t%.16lg\t%.16lg\t%.16lg\t%.16lg\t%.16lg\n",
              (EleArr + EleCntr - 1)->G.Type,
              (EleArr + EleCntr - 1)->G.Origin.X,
              (EleArr + EleCntr - 1)->G.Origin.Y,
              (EleArr + EleCntr - 1)->G.Origin.Z, (EleArr + EleCntr - 1)->G.LX,
              (EleArr + EleCntr - 1)->G.LZ, (EleArr + EleCntr - 1)->G.dA);
      fprintf(fElem, "#DirnCosn: \n");
      fprintf(fElem, "%lg, %lg, %lg\n", (EleArr + EleCntr - 1)->G.DC.XUnit.X,
              (EleArr + EleCntr - 1)->G.DC.XUnit.Y,
              (EleArr + EleCntr - 1)->G.DC.XUnit.Z);
      fprintf(fElem, "%lg, %lg, %lg\n", (EleArr + EleCntr - 1)->G.DC.YUnit.X,
              (EleArr + EleCntr - 1)->G.DC.YUnit.Y,
              (EleArr + EleCntr - 1)->G.DC.YUnit.Z);
      fprintf(fElem, "%lg, %lg, %lg\n", (EleArr + EleCntr - 1)->G.DC.ZUnit.X,
              (EleArr + EleCntr - 1)->G.DC.ZUnit.Y,
              (EleArr + EleCntr - 1)->G.DC.ZUnit.Z);
      fprintf(fElem, "#EType\tLambda\n");
      fprintf(fElem, "%d\t%lg\n", (EleArr + EleCntr - 1)->E.Type,
              (EleArr + EleCntr - 1)->E.Lambda);
      fprintf(fElem, "#NbBCs\tCPX\tCPY\tCPZ\tValue\n");
      fprintf(fElem, "%d\t%.16lg\t%.16lg\t%.16lg\t%lg\n",
              (EleArr + EleCntr - 1)->BC.NbOfBCs,
              (EleArr + EleCntr - 1)->BC.CollPt.X,
              (EleArr + EleCntr - 1)->BC.CollPt.Y,
              (EleArr + EleCntr - 1)->BC.CollPt.Z,
              (EleArr + EleCntr - 1)->BC.Value);
    }  // if OptElementFiles

    // mark bary-center and draw mesh
    if (OptGnuplot && OptGnuplotElements) {
      fprintf(fgpElem, "%g\t%g\t%g\n", (EleArr + EleCntr - 1)->BC.CollPt.X,
              (EleArr + EleCntr - 1)->BC.CollPt.Y,
              (EleArr + EleCntr - 1)->BC.CollPt.Z);

      // draw mesh
      // assign vertices of the element
      xv0 = (EleArr + EleCntr - 1)->G.Vertex[0].X;
      yv0 = (EleArr + EleCntr - 1)->G.Vertex[0].Y;
      zv0 = (EleArr + EleCntr - 1)->G.Vertex[0].Z;
      xv1 = (EleArr + EleCntr - 1)->G.Vertex[1].X;
      yv1 = (EleArr + EleCntr - 1)->G.Vertex[1].Y;
      zv1 = (EleArr + EleCntr - 1)->G.Vertex[1].Z;
      xv2 = (EleArr + EleCntr - 1)->G.Vertex[2].X;
      yv2 = (EleArr + EleCntr - 1)->G.Vertex[2].Y;
      zv2 = (EleArr + EleCntr - 1)->G.Vertex[2].Z;

      fprintf(fgpMesh, "%g\t%g\t%g\n", xv0, yv0, zv0);
      fprintf(fgpMesh, "%g\t%g\t%g\n", xv1, yv1, zv1);
      fprintf(fgpMesh, "%g\t%g\t%g\n", xv2, yv2, zv2);
      fprintf(fgpMesh, "%g\t%g\t%g\n\n", xv0, yv0, zv0);
    }  // if OptGnuplotElements

    if (k == NbSegZ)  // no rectangular element on this row
      continue;

    // determine NbSegXOnThisRow and ElLXOnThisRow for the rectagnular elements
    // and then loop.
    double RowLX = xhipt;  // the triangular portion is left outside the slice
    int NbSegXOnThisRow;
    double ElLXOnThisRow;
    if (RowLX <= SurfElLX) {
      NbSegXOnThisRow = 1;
      ElLXOnThisRow = RowLX;
    } else {
      NbSegXOnThisRow = (int)(RowLX / SurfElLX);
      ElLXOnThisRow = RowLX / NbSegXOnThisRow;
    }
    double x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3;
    for (int i = 1; i <= NbSegXOnThisRow; ++i) {
      double xorigin = (i - 1) * ElLXOnThisRow + 0.5 * ElLXOnThisRow;  // PCS
      double yorigin = 0.0;  // centroid of the rectagnular element
      double zorigin = 0.5 * (zlopt + zhipt);
      // printf("k: %d, i: %d, xo: %lg, yo: %lg, zo: %lg\n", k, i,
      // xorigin, yorigin, zorigin);
      // printf("xlopt: %lg, zlopt: %lg, xhipt: %lg, zhipt: %lg\n",
      // xlopt, zlopt, xhipt, zhipt);

      {  // Separate block for vector rotation
        Point3D localDisp, globalDisp;

        localDisp.X = xorigin;
        localDisp.Y = yorigin;
        localDisp.Z = zorigin;
        globalDisp = RotatePoint3D(&localDisp, &PrimDirnCosn, local2global);
        SurfElX = SurfX + globalDisp.X;  // GCS
        SurfElY = SurfY + globalDisp.Y;
        SurfElZ = SurfZ + globalDisp.Z;
      }  // vector rotation over
      // printf("SurfX: %lg, SurfY: %lg, SurfZ: %lg\n",
      // SurfX, SurfY, SurfZ);
      // printf("SurfElX: %lg, SurfElY: %lg, SurfElZ: %lg\n",
      // SurfElX, SurfElY, SurfElZ);

      // Assign element values and write in the file
      ++EleCntr;
      if (EleCntr > NbElements) {
        neBEMMessage("DiscretizeTriangle - EleCntr more than NbElements 2!");
        return -1;
      }

      (EleArr + EleCntr - 1)->DeviceNb =
          1;  // At present, there is only one device
      (EleArr + EleCntr - 1)->ComponentNb = SurfParentObj;
      (EleArr + EleCntr - 1)->PrimitiveNb = prim;
      (EleArr + EleCntr - 1)->Id = EleCntr;
      (EleArr + EleCntr - 1)->G.Type = 4;  // rectagnular here
      (EleArr + EleCntr - 1)->G.Origin.X = SurfElX;
      (EleArr + EleCntr - 1)->G.Origin.Y = SurfElY;
      (EleArr + EleCntr - 1)->G.Origin.Z = SurfElZ;
      (EleArr + EleCntr - 1)->G.LX = ElLXOnThisRow;
      (EleArr + EleCntr - 1)->G.LZ = SurfElLZ;
      (EleArr + EleCntr - 1)->G.LZ =
          zhipt - zlopt;  // to be on the safe side! 21/2/14
      (EleArr + EleCntr - 1)->G.dA =
          (EleArr + EleCntr - 1)->G.LX * (EleArr + EleCntr - 1)->G.LZ;
      // Safe to use the direction cosines obtained for the triangular primitive
      // since they are bound to remain unchanged for the triangular
      // sub-elements
      (EleArr + EleCntr - 1)->G.DC.XUnit.X = PrimDirnCosn.XUnit.X;
      (EleArr + EleCntr - 1)->G.DC.XUnit.Y = PrimDirnCosn.XUnit.Y;
      (EleArr + EleCntr - 1)->G.DC.XUnit.Z = PrimDirnCosn.XUnit.Z;
      (EleArr + EleCntr - 1)->G.DC.YUnit.X = PrimDirnCosn.YUnit.X;
      (EleArr + EleCntr - 1)->G.DC.YUnit.Y = PrimDirnCosn.YUnit.Y;
      (EleArr + EleCntr - 1)->G.DC.YUnit.Z = PrimDirnCosn.YUnit.Z;
      (EleArr + EleCntr - 1)->G.DC.ZUnit.X = PrimDirnCosn.ZUnit.X;
      (EleArr + EleCntr - 1)->G.DC.ZUnit.Y = PrimDirnCosn.ZUnit.Y;
      (EleArr + EleCntr - 1)->G.DC.ZUnit.Z = PrimDirnCosn.ZUnit.Z;
      (EleArr + EleCntr - 1)->E.Type = SurfEType;
      (EleArr + EleCntr - 1)->E.Lambda = SurfLambda;
      (EleArr + EleCntr - 1)->Solution = 0.0;
      (EleArr + EleCntr - 1)->Assigned = charge;
      // Boundary condition is applied at the origin for this rectangular
      // element coordinate system (ECS)
      (EleArr + EleCntr - 1)->BC.NbOfBCs = 1;  // assume one BC per element
      // Following are the location in the ECS
      (EleArr + EleCntr - 1)->BC.CollPt.X = (EleArr + EleCntr - 1)->G.Origin.X;
      (EleArr + EleCntr - 1)->BC.CollPt.Y = (EleArr + EleCntr - 1)->G.Origin.Y;
      (EleArr + EleCntr - 1)->BC.CollPt.Z = (EleArr + EleCntr - 1)->G.Origin.Z;
      // find element vertices
      // 1) displacement vector in the ECS is first identified
      // 2) this vector is transformed to the GCS
      // 3) the global displacement vector, when added to the centroid in GCS,
      // gives the node positions in GCS.
      x0 = -0.5 *
           (EleArr + EleCntr - 1)->G.LX;  // xyz displacement of node wrt ECS
      y0 = 0.0;
      z0 = -0.5 * (EleArr + EleCntr - 1)->G.LZ;
      {  // Separate block for position rotation - local2global
        Point3D localDisp, globalDisp;

        localDisp.X = x0;
        localDisp.Y = y0;
        localDisp.Z = z0;  // displacement in GCS
        globalDisp = RotatePoint3D(&localDisp, &PrimDirnCosn, local2global);
        x0 = (EleArr + EleCntr - 1)->G.Origin.X +
             globalDisp.X;  // xyz position in GCS
        y0 = (EleArr + EleCntr - 1)->G.Origin.Y + globalDisp.Y;
        z0 = (EleArr + EleCntr - 1)->G.Origin.Z + globalDisp.Z;
      }  // position rotation over

      x1 = 0.5 * (EleArr + EleCntr - 1)->G.LX;
      y1 = 0.0;
      z1 = -0.5 * (EleArr + EleCntr - 1)->G.LZ;
      {  // Separate block for position rotation - local2global
        Point3D localDisp, globalDisp;

        localDisp.X = x1;
        localDisp.Y = y1;
        localDisp.Z = z1;
        globalDisp = RotatePoint3D(&localDisp, &PrimDirnCosn, local2global);
        x1 = (EleArr + EleCntr - 1)->G.Origin.X + globalDisp.X;
        y1 = (EleArr + EleCntr - 1)->G.Origin.Y + globalDisp.Y;
        z1 = (EleArr + EleCntr - 1)->G.Origin.Z + globalDisp.Z;
      }  // position rotation over

      x2 = 0.5 * (EleArr + EleCntr - 1)->G.LX;
      y2 = 0.0;
      z2 = 0.5 * (EleArr + EleCntr - 1)->G.LZ;
      {  // Separate block for position rotation - local2global
        Point3D localDisp, globalDisp;

        localDisp.X = x2;
        localDisp.Y = y2;
        localDisp.Z = z2;
        globalDisp = RotatePoint3D(&localDisp, &PrimDirnCosn, local2global);
        x2 = (EleArr + EleCntr - 1)->G.Origin.X + globalDisp.X;
        y2 = (EleArr + EleCntr - 1)->G.Origin.Y + globalDisp.Y;
        z2 = (EleArr + EleCntr - 1)->G.Origin.Z + globalDisp.Z;
      }  // position rotation over

      x3 = -0.5 * (EleArr + EleCntr - 1)->G.LX;
      y3 = 0.0;
      z3 = 0.5 * (EleArr + EleCntr - 1)->G.LZ;
      {  // Separate block for position rotation - local2global
        Point3D localDisp, globalDisp;

        localDisp.X = x3;
        localDisp.Y = y3;
        localDisp.Z = z3;
        globalDisp = RotatePoint3D(&localDisp, &PrimDirnCosn, local2global);
        x3 = (EleArr + EleCntr - 1)->G.Origin.X + globalDisp.X;
        y3 = (EleArr + EleCntr - 1)->G.Origin.Y + globalDisp.Y;
        z3 = (EleArr + EleCntr - 1)->G.Origin.Z + globalDisp.Z;
      }  // position rotation over

      // assign vertices of the element
      (EleArr + EleCntr - 1)->G.Vertex[0].X = x0;
      (EleArr + EleCntr - 1)->G.Vertex[0].Y = y0;
      (EleArr + EleCntr - 1)->G.Vertex[0].Z = z0;
      (EleArr + EleCntr - 1)->G.Vertex[1].X = x1;
      (EleArr + EleCntr - 1)->G.Vertex[1].Y = y1;
      (EleArr + EleCntr - 1)->G.Vertex[1].Z = z1;
      (EleArr + EleCntr - 1)->G.Vertex[2].X = x2;
      (EleArr + EleCntr - 1)->G.Vertex[2].Y = y2;
      (EleArr + EleCntr - 1)->G.Vertex[2].Z = z2;
      (EleArr + EleCntr - 1)->G.Vertex[3].X = x3;
      (EleArr + EleCntr - 1)->G.Vertex[3].Y = y3;
      (EleArr + EleCntr - 1)->G.Vertex[3].Z = z3;

      if (OptElementFiles) {
        fprintf(fElem, "##Element Counter: %d\n", EleCntr);
        fprintf(fElem, "#DevNb\tCompNb\tPrimNb\tId\n");
        fprintf(fElem, "%d\t%d\t%d\t%d\n", (EleArr + EleCntr - 1)->DeviceNb,
                (EleArr + EleCntr - 1)->ComponentNb,
                (EleArr + EleCntr - 1)->PrimitiveNb,
                (EleArr + EleCntr - 1)->Id);
        fprintf(fElem, "#GType\tX\tY\tZ\tLX\tLZ\tdA\n");
        fprintf(
            fElem, "%d\t%.16lg\t%.16lg\t%.16lg\t%.16lg\t%.16lg\t%.16lg\n",
            (EleArr + EleCntr - 1)->G.Type, (EleArr + EleCntr - 1)->G.Origin.X,
            (EleArr + EleCntr - 1)->G.Origin.Y,
            (EleArr + EleCntr - 1)->G.Origin.Z, (EleArr + EleCntr - 1)->G.LX,
            (EleArr + EleCntr - 1)->G.LZ, (EleArr + EleCntr - 1)->G.dA);
        fprintf(fElem, "#DirnCosn: \n");
        fprintf(fElem, "%lg, %lg, %lg\n", (EleArr + EleCntr - 1)->G.DC.XUnit.X,
                (EleArr + EleCntr - 1)->G.DC.XUnit.Y,
                (EleArr + EleCntr - 1)->G.DC.XUnit.Z);
        fprintf(fElem, "%lg, %lg, %lg\n", (EleArr + EleCntr - 1)->G.DC.YUnit.X,
                (EleArr + EleCntr - 1)->G.DC.YUnit.Y,
                (EleArr + EleCntr - 1)->G.DC.YUnit.Z);
        fprintf(fElem, "%lg, %lg, %lg\n", (EleArr + EleCntr - 1)->G.DC.ZUnit.X,
                (EleArr + EleCntr - 1)->G.DC.ZUnit.Y,
                (EleArr + EleCntr - 1)->G.DC.ZUnit.Z);
        fprintf(fElem, "#EType\tLambda\n");
        fprintf(fElem, "%d\t%lg\n", (EleArr + EleCntr - 1)->E.Type,
                (EleArr + EleCntr - 1)->E.Lambda);
        fprintf(fElem, "#NbBCs\tCPX\tCPY\tCPZ\tValue\n");
        fprintf(fElem, "%d\t%.16lg\t%.16lg\t%.16lg\t%lg\n",
                (EleArr + EleCntr - 1)->BC.NbOfBCs,
                (EleArr + EleCntr - 1)->BC.CollPt.X,
                (EleArr + EleCntr - 1)->BC.CollPt.Y,
                (EleArr + EleCntr - 1)->BC.CollPt.Z,
                (EleArr + EleCntr - 1)->BC.Value);
      }  // if OptElementFiles

      // draw centroid and mesh
      if (OptGnuplot && OptGnuplotElements) {
        fprintf(fgpElem, "%g\t%g\t%g\n", (EleArr + EleCntr - 1)->BC.CollPt.X,
                (EleArr + EleCntr - 1)->BC.CollPt.Y,
                (EleArr + EleCntr - 1)->BC.CollPt.Z);
      }  // if OptGnuplot && OptGnuplotElements

      if (OptGnuplot && OptGnuplotElements) {
        fprintf(fgpMesh, "%g\t%g\t%g\n\n", x0, y0, z0);
        fprintf(fgpMesh, "%g\t%g\t%g\n\n", x1, y1, z1);
        fprintf(fgpMesh, "%g\t%g\t%g\n\n", x2, y2, z2);
        fprintf(fgpMesh, "%g\t%g\t%g\n\n", x3, y3, z3);
        fprintf(fgpMesh, "%g\t%g\t%g\n\n", x0, y0, z0);
      }  // if OptGnuplot && OptGnuplotElements
    }    // for i
  }      // for k
  ElementEnd[prim] = EleCntr;
  NbElmntsOnPrim[prim] = ElementEnd[prim] - ElementBgn[prim] + 1;

  if (OptPrimitiveFiles) {
    fprintf(fPrim, "Element begin: %d, Element end: %d\n", ElementBgn[prim],
            ElementEnd[prim]);
    fprintf(fPrim, "Number of elements on primitive: %d\n",
            NbElmntsOnPrim[prim]);
    fclose(fPrim);
  }

  if (OptGnuplot && OptGnuplotElements) {
    if (prim == 1)
      fprintf(fgnuElem, " '%s\' w p", gpElem);
    else
      fprintf(fgnuElem, ", \\\n \'%s\' w p", gpElem);
    if (prim == 1) {
      fprintf(fgnuMesh, " '%s\' w l", gpMesh);
      fprintf(fgnuMesh, ", \\\n \'%s\' w p ps 1", gpElem);
    } else {
      fprintf(fgnuMesh, ", \\\n \'%s\' w l", gpMesh);
      fprintf(fgnuMesh, ", \\\n \'%s\' w p ps 1", gpElem);
    }

    fclose(fgpElem);
    fclose(fgpMesh);
  }  // if OptGnuplot && OptGnuplotElements

  if (OptElementFiles) fclose(fElem);

  return (0);
}  // end of DiscretizeTriangles

// It may be noted here that the direction cosines of a given primitive are
// the same as those of the elements residing on it. There is only a chage
// of origin associated here.
int DiscretizeRectangle(int prim, int nvertex, double xvert[], double yvert[],
                        double zvert[], double xnorm, double ynorm,
                        double znorm, int volref1, int volref2, int inttype,
                        double potential, double charge, double lambda,
                        int NbSegX, int NbSegZ) {
  int SurfParentObj, SurfEType;
  double SurfX, SurfY, SurfZ, SurfLX, SurfLZ;
  double SurfElX, SurfElY, SurfElZ, SurfElLX, SurfElLZ;
  DirnCosn3D PrimDirnCosn;  // direction cosine of the current primitive
  char primstr[10];
  char gpElem[256], gpMesh[256];
  FILE *fPrim, *fElem, *fgpPrim, *fgpElem, *fgpMesh;

  // Check inputs
  if ((NbSegX <= 0) || (NbSegZ <= 0)) {
    printf("segmentation input wrong in DiscretizeRectangle ...\n");
    return -1;
  }

  if (OptPrintVertexAndNormal) {
    printf("nvertex: %d\n", nvertex);
    for (int vert = 0; vert < nvertex; ++vert) {
      printf("vert: %d, x: %lg, y: %lg, z: %lg\n", vert, xvert[vert],
             yvert[vert], zvert[vert]);
    }
    printf("Normal: %lg, %lg, %lg\n", xnorm, ynorm, znorm);
  }  // if OptPrintVertexAndNormal

  // necessary for separating filenames
  sprintf(primstr, "%d", prim);

  // in order to avoid warning messages
  fPrim = NULL;
  fElem = NULL;
  fgpMesh = NULL;
  fgpElem = NULL;

  // Get volume information, to begin with
  // int shape, material, boundarytype;
  // double eps, potential, charge;
  // neBEMVolumeDescription(volref1, &shape, &material,
  // &eps, &potential, &charge, &boundarytype);

  // compute all the properties of this surface
  SurfParentObj = 1;
  SurfEType = inttype;
  if (SurfEType == 0) {
    printf("Wrong SurfEType for prim %d\n", prim);
    exit(-1);
  }
  // centroid of the local coordinate system
  SurfX = 0.25 * (xvert[0] + xvert[1] + xvert[2] + xvert[3]);
  SurfY = 0.25 * (yvert[0] + yvert[1] + yvert[2] + yvert[3]);
  SurfZ = 0.25 * (zvert[0] + zvert[1] + zvert[2] + zvert[3]);
  // lengths of the sides
  SurfLX = sqrt((xvert[1] - xvert[0]) * (xvert[1] - xvert[0]) +
                (yvert[1] - yvert[0]) * (yvert[1] - yvert[0]) +
                (zvert[1] - zvert[0]) * (zvert[1] - zvert[0]));
  SurfLZ = sqrt((xvert[2] - xvert[1]) * (xvert[2] - xvert[1]) +
                (yvert[2] - yvert[1]) * (yvert[2] - yvert[1]) +
                (zvert[2] - zvert[1]) * (zvert[2] - zvert[1]));
  // Direction cosines
  PrimDirnCosn.XUnit.X = (xvert[1] - xvert[0]) / SurfLX;
  PrimDirnCosn.XUnit.Y = (yvert[1] - yvert[0]) / SurfLX;
  PrimDirnCosn.XUnit.Z = (zvert[1] - zvert[0]) / SurfLX;
  PrimDirnCosn.YUnit.X = xnorm;
  PrimDirnCosn.YUnit.Y = ynorm;
  PrimDirnCosn.YUnit.Z = znorm;
  PrimDirnCosn.ZUnit =
      Vector3DCrossProduct(&PrimDirnCosn.XUnit, &PrimDirnCosn.YUnit);

  // primitive direction cosine assignments
  PrimDC[prim].XUnit.X = PrimDirnCosn.XUnit.X;
  PrimDC[prim].XUnit.Y = PrimDirnCosn.XUnit.Y;
  PrimDC[prim].XUnit.Z = PrimDirnCosn.XUnit.Z;
  PrimDC[prim].YUnit.X = PrimDirnCosn.YUnit.X;
  PrimDC[prim].YUnit.Y = PrimDirnCosn.YUnit.Y;
  PrimDC[prim].YUnit.Z = PrimDirnCosn.YUnit.Z;
  PrimDC[prim].ZUnit.X = PrimDirnCosn.ZUnit.X;
  PrimDC[prim].ZUnit.Y = PrimDirnCosn.ZUnit.Y;
  PrimDC[prim].ZUnit.Z = PrimDirnCosn.ZUnit.Z;

  // primitive origin: also the barcenter for a rectangular element
  PrimOriginX[prim] = SurfX;
  PrimOriginY[prim] = SurfY;
  PrimOriginZ[prim] = SurfZ;
  PrimLX[prim] = SurfLX;
  PrimLZ[prim] = SurfLZ;

  double SurfLambda = lambda;
  double SurfV = potential;

  // file output for a primitive
  if (OptPrimitiveFiles) {
    char OutPrim[256];
    strcpy(OutPrim, ModelOutDir);
    strcat(OutPrim, "/Primitives/Primitive");
    strcat(OutPrim, primstr);
    strcat(OutPrim, ".out");
    fPrim = fopen(OutPrim, "w");
    if (fPrim == NULL) {
      neBEMMessage("DiscretizeRectangle - OutPrim");
      return -1;
    }
    fprintf(fPrim, "#prim: %d, nvertex: %d\n", prim, nvertex);
    fprintf(fPrim, "Node1: %lg\t%lg\t%lg\n", xvert[0], yvert[0], zvert[0]);
    fprintf(fPrim, "Node2: %lg\t%lg\t%lg\n", xvert[1], yvert[1], zvert[1]);
    fprintf(fPrim, "Node3: %lg\t%lg\t%lg\n", xvert[2], yvert[2], zvert[2]);
    fprintf(fPrim, "Node4: %lg\t%lg\t%lg\n", xvert[3], yvert[3], zvert[3]);
    fprintf(fPrim, "PrimOrigin: %lg\t%lg\t%lg\n", PrimOriginX[prim],
            PrimOriginY[prim], PrimOriginZ[prim]);
    fprintf(fPrim, "Primitive lengths: %lg\t%lg\n", PrimLX[prim], PrimLZ[prim]);
    fprintf(fPrim, "Norm: %lg\t%lg\t%lg\n", xnorm, ynorm, znorm);
    fprintf(fPrim, "#volref1: %d, volref2: %d\n", volref1, volref2);
    fprintf(fPrim, "#NbSegX: %d, NbSegZ: %d\n", NbSegX, NbSegZ);
    fprintf(fPrim, "#ParentObj: %d\tEType: %d\n", SurfParentObj, SurfEType);
    fprintf(fPrim, "#SurfX\tSurfY\tSurfZ\tSurfLZ\tSurfLZ\n");
    fprintf(fPrim, "%lg\t%lg\t%lg\t%lg\t%lg\n", SurfX, SurfY, SurfZ, SurfLX,
            SurfLZ);
    // fprintf(fPrim, "#SurfRX: %lg\tSurfRY: %lg\tSurfRZ: %lg\n",
    // SurfRX, SurfRY, SurfRZ);
    fprintf(fPrim, "#DirnCosn: \n");
    fprintf(fPrim, "%lg, %lg, %lg\n", PrimDirnCosn.XUnit.X,
            PrimDirnCosn.XUnit.Y, PrimDirnCosn.XUnit.Z);
    fprintf(fPrim, "%lg, %lg, %lg\n", PrimDirnCosn.YUnit.X,
            PrimDirnCosn.YUnit.Y, PrimDirnCosn.YUnit.Z);
    fprintf(fPrim, "%lg, %lg, %lg\n", PrimDirnCosn.ZUnit.X,
            PrimDirnCosn.ZUnit.Y, PrimDirnCosn.ZUnit.Z);
    fprintf(fPrim, "#SurfLambda: %lg\tSurfV: %lg\n", SurfLambda, SurfV);
  }  // if OptPrimitiveFiles

  // necessary for gnuplot
  if (OptGnuplot && OptGnuplotPrimitives) {
    char gpPrim[256];
    strcpy(gpPrim, MeshOutDir);
    strcat(gpPrim, "/GViewDir/gpPrim");
    strcat(gpPrim, primstr);
    strcat(gpPrim, ".out");
    fgpPrim = fopen(gpPrim, "w");
    if (fgpPrim == NULL) {
      neBEMMessage("DiscretizeRectangle - OutgpPrim");
      return -1;
    }
    fprintf(fgpPrim, "%g\t%g\t%g\n\n", xvert[0], yvert[0], zvert[0]);
    fprintf(fgpPrim, "%g\t%g\t%g\n\n", xvert[1], yvert[1], zvert[1]);
    fprintf(fgpPrim, "%g\t%g\t%g\n\n", xvert[2], yvert[2], zvert[2]);
    fprintf(fgpPrim, "%g\t%g\t%g\n\n", xvert[3], yvert[3], zvert[3]);
    fprintf(fgpPrim, "%g\t%g\t%g\n", xvert[0], yvert[0], zvert[0]);
    fclose(fgpPrim);

    if (prim == 1)
      fprintf(fgnuPrim, " '%s\' w l", gpPrim);
    else
      fprintf(fgnuPrim, ", \\\n \'%s\' w l", gpPrim);
  }  // if OptGnuplot && OptGnuplotPrimitives

  // file outputs for elements on primitive
  if (OptElementFiles) {
    char OutElem[256];
    strcpy(OutElem, MeshOutDir);
    strcat(OutElem, "/Elements/ElemOnPrim");
    strcat(OutElem, primstr);
    strcat(OutElem, ".out");
    fElem = fopen(OutElem, "w");
    // assert(fElem != NULL);
    if (fElem == NULL) {
      neBEMMessage("DiscretizeRectangle - OutElem");
      return -1;
    }
  }  // if OptElementFiles

  // gnuplot friendly file outputs for elements on primitive
  if (OptGnuplot && OptGnuplotElements) {
    strcpy(gpElem, MeshOutDir);
    strcat(gpElem, "/GViewDir/gpElemOnPrim");
    strcat(gpElem, primstr);
    strcat(gpElem, ".out");
    fgpElem = fopen(gpElem, "w");
    // assert(fgpElem != NULL);
    if (fgpElem == NULL) {
      neBEMMessage("DiscretizeRectangle - OutgpElem");
      if (fElem) fclose(fElem);
      return -1;
    }
    // gnuplot friendly file outputs for elements on primitive
    strcpy(gpMesh, MeshOutDir);
    strcat(gpMesh, "/GViewDir/gpMeshOnPrim");
    strcat(gpMesh, primstr);
    strcat(gpMesh, ".out");
    fgpMesh = fopen(gpMesh, "w");
    if (fgpMesh == NULL) {
      neBEMMessage("DiscretizeRectangle - OutgpMesh");
      fclose(fgpElem);
      return -1;
    }
  }  // if OptGnuplot && OptElements

  // Compute element positions (CGs) in the primitive local coordinate system.
  // Then map these CGs to the global coordinate system.
  // (xav, 0, zav) is the CG of an element wrt the primitive coordinate system
  // From this, we find the offset of the element from the primitive CG in the
  // global coordinate system by just rotating the displacement vector
  // (0,0,0) -> (xav,0,zav) using the known DCM of the surface.
  // This rotated vector (now in the global system) is then added to the
  // position vector of the primitive CG (always in the global system) to get
  // the CG of the element in the global system. Make sure that the larger
  // number is being used to discretize the longer side - can be OverSmart in
  // certain cases - there may be cases where smaller number of elements suffice
  // for a longer side
  if (NbSegX == NbSegZ) {
    SurfElLX = SurfLX / NbSegX;  // element sizes
    SurfElLZ = SurfLZ / NbSegZ;
  } else if (NbSegX > NbSegZ) {
    if (SurfLX > SurfLZ) {
      SurfElLX = SurfLX / NbSegX;  // element sizes
      SurfElLZ = SurfLZ / NbSegZ;
    } else  // interchange NbSegX and NbSegZ
    {
      int tmp = NbSegZ;
      NbSegZ = NbSegX;
      NbSegX = tmp;
      SurfElLX = SurfLX / NbSegX;  // element sizes
      SurfElLZ = SurfLZ / NbSegZ;
    }
  }     // NbSegX > NbSegZ
  else  // NbSegX < NbSegZ
  {
    if (SurfLX < SurfLZ) {
      SurfElLX = SurfLX / NbSegX;  // element sizes
      SurfElLZ = SurfLZ / NbSegZ;
    } else  // interchange NbSegX and NbSegZ
    {
      int tmp = NbSegZ;
      NbSegZ = NbSegX;
      NbSegX = tmp;
      SurfElLX = SurfLX / NbSegX;  // element sizes
      SurfElLZ = SurfLZ / NbSegZ;
    }
  }

  // Analysis of element aspect ratio.
  // Note that we can afford only to reduce the total number of elements.
  // Otherwise, we'll have to realloc `EleArr' array.
  // Later, we'll make the corrections such that the total number of elements
  // remain close to the originally intended value.
  double AR = SurfElLX / SurfElLZ;  // indicative element aspect ratio
  if (OptPrimitiveFiles) {
    fprintf(fPrim,
            "Using the input, the aspect ratio of the elements on prim: %d\n",
            prim);
    fprintf(fPrim,
            "NbSegX: %d, SurfElLX: %lg, NbSegZ: %d, SurfElLZ: %lg, AR: %lg\n",
            NbSegX, SurfElLX, NbSegZ, SurfElLZ, AR);
  }
  if (AR > 10.0)  // NbSegZ needs to be reduced
  {
    double tmpElLZ = SurfElLX / 10.0;
    NbSegZ = (int)(SurfLZ / tmpElLZ);
    if (NbSegZ <= 0) NbSegZ = 1;
    SurfElLZ = SurfLZ / NbSegZ;
  }
  if (AR < 0.1)  // NbSegX need to be reduced
  {
    double tmpElLX = SurfElLZ * 0.1;
    NbSegX = (int)(SurfLX / tmpElLX);
    if (NbSegX <= 0) NbSegX = 1;
    SurfElLX = SurfLX / NbSegX;
  }
  AR = SurfElLX / SurfElLZ;  // indicative element aspect ratio
  if (OptPrimitiveFiles) {
    fprintf(fPrim, "After analyzing the aspect ratio of the elements\n");
    fprintf(fPrim,
            "NbSegX: %d, SurfElLX: %lg, NbSegZ: %d, SurfElLZ: %lg, AR: %lg\n",
            NbSegX, SurfElLX, NbSegZ, SurfElLZ, AR);
  }

  double x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3, xav, zav;
  ElementBgn[prim] = EleCntr + 1;
  for (int i = 1; i <= NbSegX; ++i) {
    x1 = -SurfLX / 2.0 +
         (double)(i - 1) * SurfElLX;  // assuming centroid at 0,0,0
    x2 = -SurfLX / 2.0 + (double)(i)*SurfElLX;
    xav = 0.5 * (x1 + x2);

    for (int k = 1; k <= NbSegZ; ++k) {
      z1 = -SurfLZ / 2.0 + (double)(k - 1) * SurfElLZ;
      z2 = -SurfLZ / 2.0 + (double)(k)*SurfElLZ;
      zav = 0.5 * (z1 + z2);

      {  // Separate block for position rotation - local2global
        Point3D localDisp, globalDisp;

        localDisp.X = xav;
        localDisp.Y = 0.0;
        localDisp.Z = zav;
        globalDisp = RotatePoint3D(&localDisp, &PrimDirnCosn, local2global);
        SurfElX = SurfX + globalDisp.X;
        SurfElY = SurfY + globalDisp.Y;
        SurfElZ = SurfZ + globalDisp.Z;
      }  // vector rotation over

      // Assign element values and write in the file
      ++EleCntr;
      if (EleCntr > NbElements) {
        neBEMMessage("DiscretizeRectangle - EleCntr more than NbElements!");
        if (fgpMesh) fclose(fgpMesh);
        return -1;
      }

      (EleArr + EleCntr - 1)->DeviceNb =
          1;  // At present, there is only one device
      (EleArr + EleCntr - 1)->ComponentNb = SurfParentObj;
      (EleArr + EleCntr - 1)->PrimitiveNb = prim;
      (EleArr + EleCntr - 1)->Id = EleCntr;
      (EleArr + EleCntr - 1)->G.Type = 4;  // rectangular here
      (EleArr + EleCntr - 1)->G.Origin.X = SurfElX;
      (EleArr + EleCntr - 1)->G.Origin.Y = SurfElY;
      (EleArr + EleCntr - 1)->G.Origin.Z = SurfElZ;
      (EleArr + EleCntr - 1)->G.LX = SurfElLX;
      (EleArr + EleCntr - 1)->G.LZ = SurfElLZ;
      (EleArr + EleCntr - 1)->G.dA =
          (EleArr + EleCntr - 1)->G.LX * (EleArr + EleCntr - 1)->G.LZ;
      (EleArr + EleCntr - 1)->G.DC.XUnit.X = PrimDirnCosn.XUnit.X;
      (EleArr + EleCntr - 1)->G.DC.XUnit.Y = PrimDirnCosn.XUnit.Y;
      (EleArr + EleCntr - 1)->G.DC.XUnit.Z = PrimDirnCosn.XUnit.Z;
      (EleArr + EleCntr - 1)->G.DC.YUnit.X = PrimDirnCosn.YUnit.X;
      (EleArr + EleCntr - 1)->G.DC.YUnit.Y = PrimDirnCosn.YUnit.Y;
      (EleArr + EleCntr - 1)->G.DC.YUnit.Z = PrimDirnCosn.YUnit.Z;
      (EleArr + EleCntr - 1)->G.DC.ZUnit.X = PrimDirnCosn.ZUnit.X;
      (EleArr + EleCntr - 1)->G.DC.ZUnit.Y = PrimDirnCosn.ZUnit.Y;
      (EleArr + EleCntr - 1)->G.DC.ZUnit.Z = PrimDirnCosn.ZUnit.Z;
      (EleArr + EleCntr - 1)->E.Type = SurfEType;
      (EleArr + EleCntr - 1)->E.Lambda = SurfLambda;
      (EleArr + EleCntr - 1)->Solution = 0.0;
      (EleArr + EleCntr - 1)->Assigned = charge;
      (EleArr + EleCntr - 1)->BC.NbOfBCs = 1;  // assume one BC per element
      (EleArr + EleCntr - 1)->BC.CollPt.X = (EleArr + EleCntr - 1)->G.Origin.X;
      (EleArr + EleCntr - 1)->BC.CollPt.Y = (EleArr + EleCntr - 1)->G.Origin.Y;
      (EleArr + EleCntr - 1)->BC.CollPt.Z = (EleArr + EleCntr - 1)->G.Origin.Z;

      // find element vertices
      // 1) displacement vector in the ECS is first identified
      // 2) this vector is transformed to the GCS
      // 3) the global displacement vector, when added to the centroid in GCS,
      // gives the node positions in GCS.
      x0 = -0.5 *
           (EleArr + EleCntr - 1)->G.LX;  // xyz displacement of node wrt ECS
      y0 = 0.0;
      z0 = -0.5 * (EleArr + EleCntr - 1)->G.LZ;
      {  // Separate block for position rotation - local2global
        Point3D localDisp, globalDisp;

        localDisp.X = x0;
        localDisp.Y = y0;
        localDisp.Z = z0;  // displacement in GCS
        globalDisp = RotatePoint3D(&localDisp, &PrimDirnCosn, local2global);
        x0 = (EleArr + EleCntr - 1)->G.Origin.X +
             globalDisp.X;  // xyz position in GCS
        y0 = (EleArr + EleCntr - 1)->G.Origin.Y + globalDisp.Y;
        z0 = (EleArr + EleCntr - 1)->G.Origin.Z + globalDisp.Z;
      }  // position rotation over

      x1 = 0.5 * (EleArr + EleCntr - 1)->G.LX;
      y1 = 0.0;
      z1 = -0.5 * (EleArr + EleCntr - 1)->G.LZ;
      {  // Separate block for position rotation - local2global
        Point3D localDisp, globalDisp;

        localDisp.X = x1;
        localDisp.Y = y1;
        localDisp.Z = z1;
        globalDisp = RotatePoint3D(&localDisp, &PrimDirnCosn, local2global);
        x1 = (EleArr + EleCntr - 1)->G.Origin.X + globalDisp.X;
        y1 = (EleArr + EleCntr - 1)->G.Origin.Y + globalDisp.Y;
        z1 = (EleArr + EleCntr - 1)->G.Origin.Z + globalDisp.Z;
      }  // position rotation over

      x2 = 0.5 * (EleArr + EleCntr - 1)->G.LX;
      y2 = 0.0;
      z2 = 0.5 * (EleArr + EleCntr - 1)->G.LZ;
      {  // Separate block for position rotation - local2global
        Point3D localDisp, globalDisp;

        localDisp.X = x2;
        localDisp.Y = y2;
        localDisp.Z = z2;
        globalDisp = RotatePoint3D(&localDisp, &PrimDirnCosn, local2global);
        x2 = (EleArr + EleCntr - 1)->G.Origin.X + globalDisp.X;
        y2 = (EleArr + EleCntr - 1)->G.Origin.Y + globalDisp.Y;
        z2 = (EleArr + EleCntr - 1)->G.Origin.Z + globalDisp.Z;
      }  // position rotation over

      x3 = -0.5 * (EleArr + EleCntr - 1)->G.LX;
      y3 = 0.0;
      z3 = 0.5 * (EleArr + EleCntr - 1)->G.LZ;
      {  // Separate block for position rotation - local2global
        Point3D localDisp, globalDisp;

        localDisp.X = x3;
        localDisp.Y = y3;
        localDisp.Z = z3;
        globalDisp = RotatePoint3D(&localDisp, &PrimDirnCosn, local2global);
        x3 = (EleArr + EleCntr - 1)->G.Origin.X + globalDisp.X;
        y3 = (EleArr + EleCntr - 1)->G.Origin.Y + globalDisp.Y;
        z3 = (EleArr + EleCntr - 1)->G.Origin.Z + globalDisp.Z;
      }  // position rotation over

      // assign vertices of the element
      (EleArr + EleCntr - 1)->G.Vertex[0].X = x0;
      (EleArr + EleCntr - 1)->G.Vertex[0].Y = y0;
      (EleArr + EleCntr - 1)->G.Vertex[0].Z = z0;
      (EleArr + EleCntr - 1)->G.Vertex[1].X = x1;
      (EleArr + EleCntr - 1)->G.Vertex[1].Y = y1;
      (EleArr + EleCntr - 1)->G.Vertex[1].Z = z1;
      (EleArr + EleCntr - 1)->G.Vertex[2].X = x2;
      (EleArr + EleCntr - 1)->G.Vertex[2].Y = y2;
      (EleArr + EleCntr - 1)->G.Vertex[2].Z = z2;
      (EleArr + EleCntr - 1)->G.Vertex[3].X = x3;
      (EleArr + EleCntr - 1)->G.Vertex[3].Y = y3;
      (EleArr + EleCntr - 1)->G.Vertex[3].Z = z3;

      if (OptElementFiles) {
        fprintf(fElem, "##Element Counter: %d\n", EleCntr);
        fprintf(fElem, "#DevNb\tCompNb\tPrimNb\tId\n");
        fprintf(fElem, "%d\t%d\t%d\t%d\n", (EleArr + EleCntr - 1)->DeviceNb,
                (EleArr + EleCntr - 1)->ComponentNb,
                (EleArr + EleCntr - 1)->PrimitiveNb,
                (EleArr + EleCntr - 1)->Id);
        fprintf(fElem, "#GType\tX\tY\tZ\tLX\tLZ\tdA\n");
        fprintf(
            fElem, "%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
            (EleArr + EleCntr - 1)->G.Type, (EleArr + EleCntr - 1)->G.Origin.X,
            (EleArr + EleCntr - 1)->G.Origin.Y,
            (EleArr + EleCntr - 1)->G.Origin.Z, (EleArr + EleCntr - 1)->G.LX,
            (EleArr + EleCntr - 1)->G.LZ, (EleArr + EleCntr - 1)->G.dA);
        fprintf(fElem, "#DirnCosn: \n");
        fprintf(fElem, "%lg, %lg, %lg\n", (EleArr + EleCntr - 1)->G.DC.XUnit.X,
                (EleArr + EleCntr - 1)->G.DC.XUnit.Y,
                (EleArr + EleCntr - 1)->G.DC.XUnit.Z);
        fprintf(fElem, "%lg, %lg, %lg\n", (EleArr + EleCntr - 1)->G.DC.YUnit.X,
                (EleArr + EleCntr - 1)->G.DC.YUnit.Y,
                (EleArr + EleCntr - 1)->G.DC.YUnit.Z);
        fprintf(fElem, "%lg, %lg, %lg\n", (EleArr + EleCntr - 1)->G.DC.ZUnit.X,
                (EleArr + EleCntr - 1)->G.DC.ZUnit.Y,
                (EleArr + EleCntr - 1)->G.DC.ZUnit.Z);
        fprintf(fElem, "#EType\tLambda\n");
        fprintf(fElem, "%d\t%lg\n", (EleArr + EleCntr - 1)->E.Type,
                (EleArr + EleCntr - 1)->E.Lambda);
        fprintf(fElem, "#NbBCs\tCPX\tCPY\tCPZ\tValue\n");
        fprintf(fElem, "%d\t%lg\t%lg\t%lg\t%lg\n",
                (EleArr + EleCntr - 1)->BC.NbOfBCs,
                (EleArr + EleCntr - 1)->BC.CollPt.X,
                (EleArr + EleCntr - 1)->BC.CollPt.Y,
                (EleArr + EleCntr - 1)->BC.CollPt.Z,
                (EleArr + EleCntr - 1)->BC.Value);
      }  // if OptElementFiles

      // mark centroid
      if (OptGnuplot && OptGnuplotElements) {
        fprintf(fgpElem, "%g\t%g\t%g\n", (EleArr + EleCntr - 1)->BC.CollPt.X,
                (EleArr + EleCntr - 1)->BC.CollPt.Y,
                (EleArr + EleCntr - 1)->BC.CollPt.Z);
      }  // if OptGnuplot && OptGnuplotElements

      if (OptGnuplot && OptGnuplotElements) {
        fprintf(fgpMesh, "%g\t%g\t%g\n", x0, y0, z0);
        fprintf(fgpMesh, "%g\t%g\t%g\n", x1, y1, z1);
        fprintf(fgpMesh, "%g\t%g\t%g\n", x2, y2, z2);
        fprintf(fgpMesh, "%g\t%g\t%g\n", x3, y3, z3);
        fprintf(fgpMesh, "%g\t%g\t%g\n\n", x0, y0, z0);
        /*
        fprintf(fgpMesh, "%g\t%g\t%g\t%d\n", x0, y0, z0, 1);
        fprintf(fgpMesh, "%g\t%g\t%g\t%d\n\n", x1, y1, z1, 2);
        fprintf(fgpMesh, "%g\t%g\t%g\t%d\n", x2, y2, z2, 1);
        fprintf(fgpMesh, "%g\t%g\t%g\t%d\n\n\n", x2, y2, z2, 2);
        fprintf(fgpMesh, "%g\t%g\t%g\t%d\n", x0, y0, z0, 1);
        fprintf(fgpMesh, "%g\t%g\t%g\t%d\n\n", x3, y3, z3, 2);
        fprintf(fgpMesh, "%g\t%g\t%g\t%d\n", x2, y2, z2, 1);
        fprintf(fgpMesh, "%g\t%g\t%g\t%d\n\n\n", x2, y2, z2, 2);
        */
      }  // if(OptGnuplot && OptGnuplotElements)
    }    // for k
  }      // for i
  ElementEnd[prim] = EleCntr;
  NbElmntsOnPrim[prim] = ElementEnd[prim] - ElementBgn[prim] + 1;

  if (OptPrimitiveFiles) {
    fprintf(fPrim, "Element begin: %d, Element end: %d\n", ElementBgn[prim],
            ElementEnd[prim]);
    fprintf(fPrim, "Number of elements on primitive: %d\n",
            NbElmntsOnPrim[prim]);
    fclose(fPrim);
  }

  if (OptElementFiles) fclose(fElem);

  if (OptGnuplot && OptGnuplotElements) {
    if (prim == 1)
      fprintf(fgnuElem, " '%s\' w p", gpElem);
    else
      fprintf(fgnuElem, ", \\\n \'%s\' w p", gpElem);
    if (prim == 1) {
      fprintf(fgnuMesh, " '%s\' w l", gpMesh);
      fprintf(fgnuMesh, ", \\\n \'%s\' w p ps 1", gpElem);
    } else {
      fprintf(fgnuMesh, ", \\\n \'%s\' w l", gpMesh);
      fprintf(fgnuMesh, ", \\\n \'%s\' w p ps 1", gpElem);
    }

    fclose(fgpElem);
    fclose(fgpMesh);
  }

  return (0);
}  // end of DiscretizeRectangle

int DiscretizePolygon(int /*prim*/, int /*nvertex*/, double /*xvert*/[],
                      double /*yvert*/[], double /*zvert*/[], double /*xnorm*/,
                      double /*ynorm*/, double /*znorm*/, int /*volref1*/,
                      int /*volref2*/, int /*inttype*/, double /*potential*/,
                      double /*charge*/, double /*lambda*/, int NbSegX,
                      int NbSegZ) {
  // Check inputs
  if ((NbSegX <= 0) || (NbSegZ <= 0)) {
    printf("segmentation input wrong in DiscretizePolygon ...\n");
    return -1;
  }

  return (0);
}  // end of DiscretizePolygon

int BoundaryConditions(void) {
  // Loop over all the elements of all the primitives.
  // Some of the primitives do not have a boundary condition to be satisfied
  // We'll omit these primitives but before doing that we need to know to
  // which primitive a particular element belongs to.
  // The RHS will also need to be modified to accommodate the presence of
  // floating conductors and charged (known) substances.
  // Looping on primitives (rather than elements) can save precious time!
  for (int ele = 1; ele <= NbElements; ++ele) {
    int prim = (EleArr + ele - 1)->PrimitiveNb;
    // Note this condition is pure geometry, not electric!
    switch (PrimType[prim]) {
      case 2:
      case 3:  // for floating conductor and for dielectric-dielectric intfc
      case 4:  // the BC is zero (may be modified due to known charges, later)
        (EleArr + ele - 1)->BC.Value = ApplPot[prim];
        break;

      default:
        printf("Primitive out of range in BoundaryConditions ... returning\n");
        return (-1);
    }  // switch on PrimType over
  }    // ele loop over

  return (0);
}  // end of BoundaryConditions

// Set up initial conditions (C styling to be fixed)
int InitialConditions(void) {
  int fstatus = 0;

  // Known charges
  if (OptKnCh) {
    startClock = clock();
    printf("InitialConditions: InitKnownCharges ... ");
    fflush(stdout);
    fstatus = InitKnownCharges();
    if (fstatus != 0) {
      neBEMMessage("InitialConditions - InitKnownCharges");
      return -1;
    }
    printf("InitialConditions: InitKnownCharges done!\n");
    fflush(stdout);
    stopClock = clock();
    neBEMTimeElapsed(startClock, stopClock);
    printf("to InitKnownCharges.\n");
  }  // if OptKnCh

  // Charging up
  if (OptChargingUp) {
    startClock = clock();
    printf("InitialConditions: InitChargingUp ... ");
    fflush(stdout);
    fstatus = InitChargingUp();
    if (fstatus != 0) {
      neBEMMessage("InitialConditions - InitChargingUp");
      return -1;
    }
    printf("InitialConditions: InitChargingUp done!\n");
    fflush(stdout);
    stopClock = clock();
    neBEMTimeElapsed(startClock, stopClock);
    printf("to InitChargingUp.\n");
  }  // if OptChargingUp

  return 0;
}  // InitialConditions ends

// Very important note (Vin):
// There are two approaches by which existing charge distributions can be
// considered. The first one, neBEMKnownCharges, is the simpler approach
// in which the external charges are left as they are. In the second one,
// neBEMChargingUp, the charges that have crossed an element of a dielectric
// interface, is considered to be sitting on the element. All such accumulated
// charges are assigned to that element by suitably modifying
// (EleArr+elesrc-1)->Assigned.
// However, possible pitfalls such as considering the effects twice, exist.
// Hence the code needs to be carefully exampined.
// int neBEMKnownCharges(int (*Pt2UserFunction)(void)) - apt for UpdateKnCh
int InitKnownCharges(void) {
  int debugFn = 0;

  /*
  // How to check that the pointer points to a valid function?
  if(Pt2UserFunction == NULL)
          {
          printf("Not a valid function ... returning ...\n");
          return(-1);
          }
  else
          {
          // printf("Pt2UserFunction points to %p\n", Pt2UserFunction);
          }

  // Status of the known charges conditions is meaningful only after the
  // element discretization has been completed, i.e., beyond the 5th state.
  if(neBEMState >= 5)
          {	// the following function is declared in the Interface.h.
          int fstatus = (*Pt2UserFunction)();	// user to supply function
          if(fstatus != 0)
                  {
                  neBEMMessage("neBEMKnownCharges - Pt2UserFunction");
                  return -1;
                  }
          if(neBEMState > 5)	// assume LHS and inversion to be over?
                  neBEMState = 8;
          }
  else
          {
          printf("Known charges are meaningful only beyond state 4 ...\n");
          printf("returning ...\n");
          return(-1);
          }
  */

  // Set up parameters related to known charge calculations
  // Electrons and ions can be distributed as points, lines, areas or volumes.
  {
    FILE *KnChInpFile = fopen("neBEMInp/neBEMKnCh.inp", "r");
    if (KnChInpFile == NULL) {
      OptKnCh = 0;
      printf(
          "neBEMKnCh.inp absent ... assuming absence of known charges ...\n");
    } else {
      fscanf(KnChInpFile, "OptKnCh: %d\n", &OptKnCh);
      if (1) printf("OptKnCh: %d\n", OptKnCh);

      if (!OptKnCh) printf("OptKnCh = 0 ... assuming no known charges ...\n");

      if (OptKnCh) {
        char PointKnChFile[256];
        char LineKnChFile[256];
        char AreaKnChFile[256];
        char VolumeKnChFile[256];
        double LScaleFactor;
        double KnChFactor;

        fscanf(KnChInpFile, "NbPointsKnCh: %d\n", &NbPointsKnCh);
        fscanf(KnChInpFile, "PointKnChFile: %s\n", PointKnChFile);
        fscanf(KnChInpFile, "NbLinesKnCh: %d\n", &NbLinesKnCh);
        fscanf(KnChInpFile, "LineKnChFile: %s\n", LineKnChFile);
        fscanf(KnChInpFile, "NbAreasKnCh: %d\n", &NbAreasKnCh);
        fscanf(KnChInpFile, "AreaKnChFile: %s\n", AreaKnChFile);
        fscanf(KnChInpFile, "NbVolumesKnCh: %d\n", &NbVolumesKnCh);
        fscanf(KnChInpFile, "VolumeKnChFile: %s\n", VolumeKnChFile);
        fscanf(KnChInpFile, "LScaleFactor: %lg\n", &LScaleFactor);
        fscanf(KnChInpFile, "KnChFactor: %lg\n", &KnChFactor);
        if (1) {
          printf("NbPointsKnCh: %d\n", NbPointsKnCh);
          printf("PointKnChFile: %s\n", PointKnChFile);
          printf("NbLinesKnCh: %d\n", NbLinesKnCh);
          printf("LineKnChFile: %s\n", LineKnChFile);
          printf("NbAreasKnCh: %d\n", NbAreasKnCh);
          printf("AreaKnChFile: %s\n", AreaKnChFile);
          printf("NbVolumesKnCh: %d\n", NbVolumesKnCh);
          printf("VolumeKnChFile: %s\n", VolumeKnChFile);
          printf("LScaleFactor: %lg\n", LScaleFactor);
          printf("KnChFactor: %lg\n", KnChFactor);
        }

        if (NbPointsKnCh) {  // Points
          if (debugFn)
            printf("No. of points with known charges: %d\n", NbPointsKnCh);

          FILE *fptrPointKnChFile = fopen(PointKnChFile, "r");
          if (fptrPointKnChFile == NULL) {
            neBEMMessage("PointKnCh file absent ... returning\n");
            return -10;
          }

          PointKnChArr =
              (PointKnCh *)malloc((NbPointsKnCh + 1) * sizeof(PointKnCh));

          for (int point = 0; point <= NbPointsKnCh; ++point) {  
            // CHECK!!! ele limits start from 0, but all else from 1 to ...
            PointKnChArr[point].Nb = 0;
            PointKnChArr[point].P.X = 0.0;
            PointKnChArr[point].P.Y = 0.0;
            PointKnChArr[point].P.Z = 0.0;
            PointKnChArr[point].Assigned = 0.0;
          }

          char header[256];
          fgets(header, 256, fptrPointKnChFile);  // header

          for (int point = 1; point <= NbPointsKnCh; ++point) {
            fscanf(fptrPointKnChFile, "%d %lg %lg %lg %lg\n",
                   &PointKnChArr[point].Nb, &PointKnChArr[point].P.X,
                   &PointKnChArr[point].P.Y, &PointKnChArr[point].P.Z,
                   &PointKnChArr[point].Assigned);

            PointKnChArr[point].P.X /= LScaleFactor;     // convert Garfield cm
            PointKnChArr[point].P.Y /= LScaleFactor;     // to SI meter
            PointKnChArr[point].P.Z /= LScaleFactor;     // and similar other
            PointKnChArr[point].Assigned *= KnChFactor;  // requirements

            if (debugFn) {
              printf("Nb: %d , X: %lg, Y:%lg, Z:%lg, Assigned:%lg\n",
                     PointKnChArr[point].Nb, PointKnChArr[point].P.X,
                     PointKnChArr[point].P.Y, PointKnChArr[point].P.Z,
                     PointKnChArr[point].Assigned);
            }
          }

          fclose(fptrPointKnChFile);
          if (debugFn) printf("done for all points\n");
        }  // if NbPointsKnCh: Inputs and calculations for points ends

        if (NbLinesKnCh) {  // Lines
          if (debugFn)
            printf("No. of lines with known charges: %d\n", NbLinesKnCh);

          FILE *fptrLineKnChFile = fopen(LineKnChFile, "r");
          if (fptrLineKnChFile == NULL) {
            neBEMMessage("LineKnCh file absent ... returning\n");
            return -10;
          }

          LineKnChArr =
              (LineKnCh *)malloc((NbLinesKnCh + 1) * sizeof(LineKnCh));

          for (int line = 0; line <= NbLinesKnCh; ++line) {  
            // CHECK!!! ele limits start from 0, but all else from 1 to ...
            LineKnChArr[line].Nb = 0;
            LineKnChArr[line].Start.X = 0.0;
            LineKnChArr[line].Start.Y = 0.0;
            LineKnChArr[line].Start.Z = 0.0;
            LineKnChArr[line].Stop.X = 0.0;
            LineKnChArr[line].Stop.Y = 0.0;
            LineKnChArr[line].Stop.Z = 0.0;
            LineKnChArr[line].Radius = 0.0;
            LineKnChArr[line].Assigned = 0.0;
          }

          char header[256];
          fgets(header, 256, fptrLineKnChFile);  // header

          for (int line = 1; line <= NbLinesKnCh; ++line) {
            fscanf(fptrLineKnChFile, "%d %lg %lg %lg %lg %lg %lg %lg %lg\n",
                   &LineKnChArr[line].Nb, &LineKnChArr[line].Start.X,
                   &LineKnChArr[line].Start.Y, &LineKnChArr[line].Start.Z,
                   &LineKnChArr[line].Stop.X, &LineKnChArr[line].Stop.Y,
                   &LineKnChArr[line].Stop.Z, &LineKnChArr[line].Radius,
                   &LineKnChArr[line].Assigned);

            LineKnChArr[line].Start.X /= LScaleFactor;  // cm to m
            LineKnChArr[line].Start.Y /= LScaleFactor;  // and similar other
            LineKnChArr[line].Start.Z /= LScaleFactor;
            LineKnChArr[line].Stop.X /= LScaleFactor;
            LineKnChArr[line].Stop.Y /= LScaleFactor;
            LineKnChArr[line].Stop.Z /= LScaleFactor;
            LineKnChArr[line].Radius /= LScaleFactor;
            LineKnChArr[line].Assigned *= LScaleFactor;  // Q/cm to Q/m
            LineKnChArr[line].Assigned *= KnChFactor;

            if (debugFn) {
              printf(
                  "Nb: %d, X1Y1Z1: %lg %lg %lg X2Y2Z2 %lg %lg %lg R: %lg "
                  "Lmda:%lg\n",
                  LineKnChArr[line].Nb, LineKnChArr[line].Start.X,
                  LineKnChArr[line].Start.Y, LineKnChArr[line].Start.Z,
                  LineKnChArr[line].Stop.X, LineKnChArr[line].Stop.Y,
                  LineKnChArr[line].Stop.Z, LineKnChArr[line].Radius,
                  LineKnChArr[line].Assigned);
            }
          }

          fclose(fptrLineKnChFile);
          if (debugFn) printf("done for all lines\n");
        }  // if NbLinesKnCh: Inputs and calculations for lines ends

        if (NbAreasKnCh) {  // Areas
          if (debugFn)
            printf("No. of areas with known charges: %d\n", NbAreasKnCh);

          FILE *fptrAreaKnChFile = fopen(AreaKnChFile, "r");
          if (fptrAreaKnChFile == NULL) {
            neBEMMessage("AreaKnCh file absent ... returning\n");
            return -10;
          }

          AreaKnChArr =
              (AreaKnCh *)malloc((NbAreasKnCh + 1) * sizeof(AreaKnCh));

          for (int area = 0; area <= NbAreasKnCh; ++area) {  
            // CHECK!!! ele limits start from 0, but all else from 1 to ...
            AreaKnChArr[area].Nb = 0;
            AreaKnChArr[area].NbVertices = 0;
            for (int vert = 1; vert <= (AreaKnChArr + area - 1)->NbVertices;
                 ++vert) {
              (AreaKnChArr + area - 1)->Vertex[vert].X = 0.0;
              (AreaKnChArr + area - 1)->Vertex[vert].Y = 0.0;
              (AreaKnChArr + area - 1)->Vertex[vert].Z = 0.0;
            }
            AreaKnChArr[area].Assigned = 0.0;
          }

          char header[256];
          fgets(header, 256, fptrAreaKnChFile);  // header

          for (int area = 1; area <= NbAreasKnCh; ++area) {
            fscanf(fptrAreaKnChFile, "%d %d %le\n",
                   &(AreaKnChArr + area - 1)->Nb,
                   &(AreaKnChArr + area - 1)->NbVertices,
                   &(AreaKnChArr + area - 1)->Assigned);
            for (int vert = 1; vert <= (AreaKnChArr + area - 1)->NbVertices;
                 ++vert) {
              fscanf(fptrAreaKnChFile, "%le %le %le\n",
                     &(AreaKnChArr + area - 1)->Vertex[vert].X,
                     &(AreaKnChArr + area - 1)->Vertex[vert].Y,
                     &(AreaKnChArr + area - 1)->Vertex[vert].Z);
            }

            for (int vert = 1; vert <= (AreaKnChArr + area - 1)->NbVertices;
                 ++vert) {  // cm to m
              (AreaKnChArr + area - 1)->Vertex[vert].X /= LScaleFactor;
              (AreaKnChArr + area - 1)->Vertex[vert].Y /= LScaleFactor;
              (AreaKnChArr + area - 1)->Vertex[vert].Z /= LScaleFactor;
            }
            AreaKnChArr[area].Assigned *=
                (LScaleFactor * LScaleFactor);  // cm^2 to m^2
            AreaKnChArr[area].Assigned *= KnChFactor;
          }

          fclose(fptrAreaKnChFile);
          if (debugFn) printf("done for all areas\n");
        }  // if AreasKnCh: Inputs and calculations for areas ends

        if (NbVolumesKnCh) {  // Volumes
          if (debugFn)
            printf("No. of volumes with known charges: %d\n", NbVolumesKnCh);

          FILE *fptrVolumeKnChFile = fopen(VolumeKnChFile, "r");
          if (fptrVolumeKnChFile == NULL) {
            neBEMMessage("VolumeKnCh file absent ... returning\n");
            return -10;
          }

          VolumeKnChArr =
              (VolumeKnCh *)malloc((NbVolumesKnCh + 1) * sizeof(VolumeKnCh));

          for (int volume = 0; volume <= NbVolumesKnCh; ++volume) {  
            // CHECK!!! ele limits start from 0, but all else from 1 to ...
            VolumeKnChArr[volume].Nb = 0;
            VolumeKnChArr[volume].NbVertices = 0;
            for (int vert = 1; vert <= (VolumeKnChArr + volume - 1)->NbVertices;
                 ++vert) {
              (VolumeKnChArr + volume - 1)->Vertex[vert].X = 0.0;
              (VolumeKnChArr + volume - 1)->Vertex[vert].Y = 0.0;
              (VolumeKnChArr + volume - 1)->Vertex[vert].Z = 0.0;
            }
            VolumeKnChArr[volume].Assigned *=
                (LScaleFactor * LScaleFactor * LScaleFactor);  // cm^3 to m^3
            VolumeKnChArr[volume].Assigned = 0.0;
          }

          char header[256];
          fgets(header, 256, fptrVolumeKnChFile);  // header

          for (int volume = 1; volume <= NbVolumesKnCh; ++volume) {
            fscanf(fptrVolumeKnChFile, "%d %d %le\n",
                   &(VolumeKnChArr + volume - 1)->Nb,
                   &(VolumeKnChArr + volume - 1)->NbVertices,
                   &(VolumeKnChArr + volume - 1)->Assigned);
            for (int vert = 1; vert <= (VolumeKnChArr + volume - 1)->NbVertices;
                 ++vert) {
              fscanf(fptrVolumeKnChFile, "%le %le %le\n",
                     &(VolumeKnChArr + volume - 1)->Vertex[vert].X,
                     &(VolumeKnChArr + volume - 1)->Vertex[vert].Y,
                     &(VolumeKnChArr + volume - 1)->Vertex[vert].Z);
            }

            for (int vert = 1; vert <= (VolumeKnChArr + volume - 1)->NbVertices;
                 ++vert) {  // cm to m
              (VolumeKnChArr + volume - 1)->Vertex[vert].X /= LScaleFactor;
              (VolumeKnChArr + volume - 1)->Vertex[vert].Y /= LScaleFactor;
              (VolumeKnChArr + volume - 1)->Vertex[vert].Z /= LScaleFactor;
            }
            VolumeKnChArr[volume].Assigned *= KnChFactor;
          }

          fclose(fptrVolumeKnChFile);
          if (debugFn) printf("done for all volumes\n");
        }  // if NbVolumesKnCh: Inputs and calculations for volumes ends

      }  // if OptKnCh

      fclose(KnChInpFile);
    }  // else KnChInpFile
  }    // parameters related to known charge calculations

  return (0);
}  // InitKnownCharges ends

// It is quite possible that the elements had a charge assgined to them
// even before they are being considered for getting charged up.
// Moreover, following the sequence in which the algorithm has been developed,
// the same element accumulates the available electrons and then the ions.
// The resultant of all three (prior charge, electrons and ions) turns out
// to be the assigned charge on the elements, after the execution of this
// function.
int InitChargingUp(void) {
  int debugFn = 0;

  // status of elements before being charged up
  if (debugFn) {
    for (int ele = 1; ele <= NbElements; ++ele) {
      printf("ele, Assigned charge: %d, %lg\n", ele,
             (EleArr + ele - 1)->Assigned);
    }
  }

  // Set up parameters related to charging-up calculations
  // The plan is to distribute the electrons and ions ending in dielectric
  // volumes to the elements on the volumes
  {
    FILE *ChargingUpInpFile = fopen("neBEMInp/neBEMChargingUp.inp", "r");
    if (ChargingUpInpFile == NULL) {
      printf(
          "neBEMChargingUp.inp absent ... assuming no charging up effect "
          "...\n");
      // assign NbChUpEEle and NbChUpIEle and Prims to be zeros?
    } else {
      fscanf(ChargingUpInpFile, "OptChargingUp: %d\n", &OptChargingUp);
      if (!OptChargingUp)
        printf("OptChargingUp = 0 ... assuming no charging up effect ...\n");
      if (1) printf("OptChargingUp: %d\n", OptChargingUp);

      if (OptChargingUp) {
        char ChargingUpEFile[256];
        char ChargingUpIFile[256];
        double LScaleFactor;
        double ChUpFactor;
        int *NbChUpEonEle, *NbChUpIonEle;

        fscanf(ChargingUpInpFile, "ChargingUpEFile: %s\n", ChargingUpEFile);
        fscanf(ChargingUpInpFile, "ChargingUpIFile: %s\n", ChargingUpIFile);
        fscanf(ChargingUpInpFile, "LScaleFactor: %lg\n", &LScaleFactor);
        fscanf(ChargingUpInpFile, "ChUpFactor: %lg\n", &ChUpFactor);
        if (1) {
          printf("ChargingUpEFile: %s\n", ChargingUpEFile);
          printf("ChargingUpIFile: %s\n", ChargingUpIFile);
          printf("LScaleFactor: %lg\n", LScaleFactor);
          printf("ChUpFactor: %lg\n", ChUpFactor);
        }

        {  // Calculation for electrons
          FILE *fptrChargingUpEFile = fopen(ChargingUpEFile, "r");
          if (fptrChargingUpEFile == NULL) {
            neBEMMessage("ChargingUpE file absent ... returning\n");
            return -10;
          }
          int NbOfE = neBEMGetNbOfLines(ChargingUpEFile);
          if (NbOfE <= 1) {
            neBEMMessage("Too few lines in ChargingUpE ... returning\n");
            return -11;
          } else {  // initialize
            NbChUpEonEle = (int *)malloc((NbElements + 1) * sizeof(int));
            for (int ele = 0; ele <= NbElements; ++ele) {  
              // CHECK!!! ele limits start from 0, but all else from 1 to ...
              NbChUpEonEle[ele] = 0;
            }
          }

          // read the header line
          char header[256];
          fgets(header, 256, fptrChargingUpEFile);

          --NbOfE;  // one line was for the header
          if (debugFn) printf("No. of electrons: %d\n", NbOfE);
          char tmpEFile[256];
          strcpy(tmpEFile, "/tmp/ElectronTempFile.out");
          FILE *ftmpEF = fopen(tmpEFile, "w");
          if (ftmpEF == NULL) {
            printf("cannot open temporary output file ... returning ...\n");
            return -100;
          }
          FILE *fPtEChUpMap = fopen("PtEChUpMap.out", "w");
          if (fPtEChUpMap == NULL) {
            printf("cannot open PtEChUpMap.out file for writing ...\n");
            return 110;
          }

          char label;
          int vol, enb;  // label, volume and electron number
          double xlbend, ylbend, zlbend, xend, yend,
              zend;  // lbend == Last But END
          Point3D
              ptintsct;  // each electron likely to have an intersection point
          for (int electron = 1; electron <= NbOfE; ++electron) {
            fscanf(fptrChargingUpEFile, "%c %d %d %lg %lg %lg %lg %lg %lg\n",
                   &label, &vol, &enb, &xlbend, &xend, &ylbend, &yend, &zlbend,
                   &zend);
            xlbend /= LScaleFactor;  // cm to metre
            xend /= LScaleFactor;    // and similar other
            ylbend /= LScaleFactor;
            yend /= LScaleFactor;
            zlbend /= LScaleFactor;
            zend /= LScaleFactor;
            ptintsct.X = 0.0;
            ptintsct.Y = 0.0;
            ptintsct.Z = 0.0;  // initialize

            // find the parametric equation of this last segment
            // if xend < xlbend, swap the directions
            // This has not been mentioned as mandatory in Vince's book
            // "Geometry for Computer Graphics", but implied in the book "A
            // Programmer's Geometry"
            if (xend < xlbend)  // not implemented now
            {
            }
            double lseg = (xend - xlbend) * (xend - xlbend) +
                          (yend - ylbend) * (yend - ylbend) +
                          (zend - zlbend) * (zend - zlbend);
            lseg = sqrt(lseg);
            double xgrd =
                (xend - xlbend) / lseg;  // normalized direction vector
            double ygrd = (yend - ylbend) / lseg;
            double zgrd = (zend - zlbend) / lseg;
            if (debugFn) {
              printf("\nelectron: %d\n", electron);
              printf("xlbend: %lg, ylbend: %lg, zlbend: %lg\n", xlbend, ylbend,
                     zlbend);
              printf("xend: %lg, yend: %lg, zend: %lg\n", xend, yend, zend);
              printf("xgrd: %lg, ygrd: %lg, zgrd: %lg\n", xgrd, ygrd, zgrd);
              fprintf(ftmpEF, "#e: %d, label: %c, vol: %d\n", electron, label,
                      vol);
              fprintf(ftmpEF, "%lg %lg %lg\n", xlbend, ylbend, zlbend);
              fprintf(ftmpEF, "%lg %lg %lg\n", xend, yend, zend);
              fprintf(ftmpEF, "#xgrd: %lg, ygrd: %lg, zgrd: %lg\n", xgrd, ygrd,
                      zgrd);
              fprintf(ftmpEF, "\n");
            }

            // determine which element gets this electron
            // At present, the logic is as follow:
            // Using the information on the last segment, find out which
            // primitive is pierced by it and at which point From the
            // intersection point, find out the element number on the primitive
            // in a local sense Using the information (start and end elements of
            // a given primitive) identify the element in a global sense. This
            // approach should be lot more efficient than checking intersection
            // element by element.
            // The intersection point is computed following algorithm
            // implemented in a matlab code (plane_imp_line_par_int_3d.m) Also
            // check which primitive in the list is the closet to the end point

            double SumOfAngles;
            int PrimIntsctd = -1,
                EleIntsctd = -1;   // intersected primitive & element
            int nearestprim = -1;  // absurd value
            double dist = 1.0e6, mindist = 1.0e6;  // absurdly high numbers
            // check all primitives
            for (int prim = 1; prim <= NbPrimitives; ++prim) { 
              if (InterfaceType[prim] != 4)
                continue;  // primitive not a dielectric

              int intersect = 0, extrasect = 1;  // worst of conditions
              int InPrim = 0, InEle = 0;
              if (debugFn)
                printf("prim: %d, mindist: %lg, nearestprim: %d\n", prim,
                       mindist, nearestprim);

              // Use two nodes at a time to get two independent vectors on
              // primitive Get cross-product of these two vector - normal to the
              // plane Note that the normal is already associated with the
              // primitive of type 3 and 4; 2 is wire and does not have any
              // associated normal
              if ((PrimType[prim] == 3) || (PrimType[prim] == 4)) {
                if (debugFn) {
                  printf("prim: %d\n", prim);
                  printf("vertex0: %lg, %lg, %lg\n", XVertex[prim][0],
                         YVertex[prim][0], ZVertex[prim][0]);
                  printf("vertex1: %lg, %lg, %lg\n", XVertex[prim][1],
                         YVertex[prim][1], ZVertex[prim][1]);
                  printf("vertex2: %lg, %lg, %lg\n", XVertex[prim][2],
                         YVertex[prim][2], ZVertex[prim][2]);
                  if (PrimType[prim] == 4) {
                    printf("vertex3: %lg, %lg, %lg\n", XVertex[prim][3],
                           YVertex[prim][3], ZVertex[prim][3]);
                  }
                  fprintf(ftmpEF, "#prim: %d\n", prim);
                  fprintf(ftmpEF, "%lg %lg %lg\n", XVertex[prim][0],
                          YVertex[prim][0], ZVertex[prim][0]);
                  fprintf(ftmpEF, "%lg %lg %lg\n", XVertex[prim][1],
                          YVertex[prim][1], ZVertex[prim][1]);
                  fprintf(ftmpEF, "%lg %lg %lg\n", XVertex[prim][2],
                          YVertex[prim][2], ZVertex[prim][2]);
                  if (PrimType[prim] == 4) {
                    fprintf(ftmpEF, "%lg %lg %lg\n", XVertex[prim][3],
                            YVertex[prim][3], ZVertex[prim][3]);
                  }
                  fprintf(ftmpEF, "%lg %lg %lg\n", XVertex[prim][0],
                          YVertex[prim][0], ZVertex[prim][0]);
                  fprintf(ftmpEF, "\n");
                  fflush(stdout);
                }  // debugFn

                // use a, b, c (normal is ai + bj + ck) at one of the nodes to
                // get d ax + by + cz + d = 0 is the equation of the plane
                double a = XNorm[prim];
                double b = YNorm[prim];
                double c = ZNorm[prim];
                double d = -a * XVertex[prim][0] - b * YVertex[prim][0] -
                           c * ZVertex[prim][0];

                // distance of the end point to this primitve is
                dist = (xend * a + yend * b + zend * c + d) /
                       sqrt(a * a + b * b + c * c);
                dist = fabs(dist);  // if only magnitude is required
                if (prim == 1) {
                  mindist = dist;
                  nearestprim = prim;
                }
                if ((prim == 1) && debugFn)
                  printf(
                      "after prim == 1 check mindist: %lg, nearestprim: %d\n",
                      mindist, nearestprim);

                // Point of intersection
                // Algo as on p62 (pdf 81) of Vince - Geometry for Computer
                // Graphics 1.5.13 Intersection of a line and a plane Algorithm
                // as implemented in plne_imp_line_par_int_3d.m a (nx), b (ny),
                // c (nz), d are a, b, c, d vx, vy, vz are xgrd, ygrd and zgrd
                // tx, ty, tz are xlbend, ylbend, zlbend
                // In the present case, n and v are unit vectors
                double norm1 = sqrt(a * a + b * b + c * c);
                double norm2 = sqrt(xgrd * xgrd + ygrd * ygrd + zgrd * zgrd);
                double denom =
                    a * xgrd + b * ygrd + c * zgrd;  // (vec)n.(vec)v; if 0, ||
                double tol =
                    1.0e-16;  // CHECK: -8 in original code; sizes small here
                intersect = extrasect = 0;

                if (debugFn) {
                  printf("a, b, c, d, dist: %lg, %lg, %lg, %lg, %lg\n", a, b, c,
                         d, dist);
                  printf("vector n: ai + bj + ck\n");
                  printf("vector v: xgrd, ygrd, zgrd: %lg, %lg, %lg\n", xgrd,
                         ygrd, zgrd);
                  printf("norm1, norm2, (vec n . vec v) denom: %lg, %lg, %lg\n",
                         norm1, norm2, denom);
                  printf("if vec n . vec v == 0, line and plane parallel\n");
                  fflush(stdout);
                }

                if (fabs(denom) < tol * norm1 * norm2) { 
                  // line parallel to the plane
                  if (fabs(a * xlbend + b * ylbend + c * zlbend + d) <=
                      1.0e-16) {  // CHECK: was == 0.0 in original code
                    intersect = 1;
                    extrasect = 0;  // line ends on the plane
                    ptintsct.X = xlbend;
                    ptintsct.Y = ylbend;
                    ptintsct.Z = zlbend;
                  } else {
                    intersect = 0;
                    extrasect = 1;     // both wrong
                    ptintsct.X = 0.0;  // Wrong to assign 0 values
                    ptintsct.Y =
                        0.0;  // However, they are never going to be used
                    ptintsct.Z = 0.0;  // since intersect is 0
                  }
                  if (debugFn) {
                    printf("line and plane parallel ...\n");
                    printf("intersect: %d, extrasect: %d\n", intersect,
                           extrasect);
                    printf("intersection point: %lg, %lg, %lg\n", ptintsct.X,
                           ptintsct.Y, ptintsct.Z);
                  }       // if line and plane are parallel
                } else {  // if they are not parallel, they must intersect
                  intersect = 1;
                  double t =
                      -(a * xlbend + b * ylbend + c * zlbend + d) / denom;

                  // check whether t is less than the length of the segment
                  // and in the correct direction
                  // If not, then an extrapolated intersection is not of
                  // interest
                  if ((t < 0.0) ||
                      (fabs(t) > fabs(lseg)))  // wrong dirn or beyond end
                  {
                    extrasect = 1;
                    ptintsct.X = xlbend + t * xgrd;
                    ptintsct.Y = ylbend + t * ygrd;
                    ptintsct.Z = zlbend + t * zgrd;
                  } else {
                    extrasect = 0;
                    ptintsct.X = xlbend + t * xgrd;
                    ptintsct.Y = ylbend + t * ygrd;
                    ptintsct.Z = zlbend + t * zgrd;
                  }
                  if (debugFn) {
                    printf("line and plane NOT parallel ...\n");
                    printf("intersect: %d, extrasect: %d\n", intersect,
                           extrasect);
                    printf("intersection point: %lg, %lg, %lg\n", ptintsct.X,
                           ptintsct.Y, ptintsct.Z);
                    printf("t, lseg: %lg, %lg\n", t, lseg);
                    printf(
                        "for an interesting intersection, lseg > t > 0.0 "
                        "...\n\n");
                    fflush(stdout);
                  }   // must intersect
                }     // if not parallel
              }       // if PrimType is 3 or 4
              else {  // this is a wire primitive - assume no charging up issues
                dist = -1.0;  // an absurd negative distance
                intersect = 0;
                extrasect = 0;
                continue;
              }  // else PrimType 3 or 4

              if (dist < mindist) {
                mindist = dist;
                nearestprim = prim;
              }
              if (debugFn)
                printf("nearestprim: %d, mindist: %lg\n\n", nearestprim,
                       mindist);

              // implicit assumption: the first primitive that gets pierced by
              // the ray is the one that we want. There can be other primitives
              // that are pierced by the same ray. So, this logic should be
              // refined further
              if ((intersect == 1) && (extrasect == 0)) {
                // check whether the intersection point is within primitive
                // polygon
                int nvert = PrimType[prim];
                Point3D polynode[4];
                polynode[0].X = XVertex[prim][0];
                polynode[0].Y = YVertex[prim][0];
                polynode[0].Z = ZVertex[prim][0];
                polynode[1].X = XVertex[prim][1];
                polynode[1].Y = YVertex[prim][1];
                polynode[1].Z = ZVertex[prim][1];
                polynode[2].X = XVertex[prim][2];
                polynode[2].Y = YVertex[prim][2];
                polynode[2].Z = ZVertex[prim][2];
                if (PrimType[prim] == 4) {
                  polynode[3].X = XVertex[prim][3];
                  polynode[3].Y = YVertex[prim][3];
                  polynode[3].Z = ZVertex[prim][3];
                }
                // printf("neBEMChkInPoly for primitive %d\n", prim);
                SumOfAngles = neBEMChkInPoly(nvert, polynode, ptintsct);
                if (fabs(fabs(SumOfAngles) - neBEMtwopi) <= 1.0e-8) {
                  InPrim = 1;
                  PrimIntsctd = prim;
                }
                if (debugFn) {
                  // print polynode and InPrim
                  printf("Prim: %d\n", prim);
                  printf("ptintsct: %lg, %lg, %lg\n", ptintsct.X, ptintsct.Y,
                         ptintsct.Z);
                  printf("nvert: %d\n", nvert);
                  printf("polynode0: %lg, %lg, %lg\n", polynode[0].X,
                         polynode[0].Y, polynode[0].Z);
                  printf("polynode1: %lg, %lg, %lg\n", polynode[1].X,
                         polynode[1].Y, polynode[1].Z);
                  printf("polynode2: %lg, %lg, %lg\n", polynode[2].X,
                         polynode[2].Y, polynode[2].Z);
                  if (nvert == 4) {
                    printf("polynode3: %lg, %lg, %lg\n", polynode[3].X,
                           polynode[3].Y, polynode[3].Z);
                  }
                  printf("SumOfAngles: %lg, InPrim: %d\n", SumOfAngles, InPrim);
                  fflush(stdout);
                }

                if (!InPrim && (prim != NbPrimitives)) {
                  continue;  // check next primitive
                }

                // Once identified, check in which element belonging to this
                // primitive contains the point of intersection
                if (InPrim) {
                  InEle = 0;
                  for (int ele = ElementBgn[prim]; ele <= ElementEnd[prim];
                       ++ele) {
                    nvert = 0;
                    if ((EleArr + ele - 1)->G.Type == 3) nvert = 3;
                    if ((EleArr + ele - 1)->G.Type == 4) nvert = 4;
                    if (!nvert) {
                      neBEMMessage(
                          "no vertex in element! ... neBEMKnownCharges ...\n");
                      return -20;
                    }

                    polynode[0].X = (EleArr + ele - 1)->G.Vertex[0].X;
                    polynode[0].Y = (EleArr + ele - 1)->G.Vertex[0].Y;
                    polynode[0].Z = (EleArr + ele - 1)->G.Vertex[0].Z;
                    polynode[1].X = (EleArr + ele - 1)->G.Vertex[1].X;
                    polynode[1].Y = (EleArr + ele - 1)->G.Vertex[1].Y;
                    polynode[1].Z = (EleArr + ele - 1)->G.Vertex[1].Z;
                    polynode[2].X = (EleArr + ele - 1)->G.Vertex[2].X;
                    polynode[2].Y = (EleArr + ele - 1)->G.Vertex[2].Y;
                    polynode[2].Z = (EleArr + ele - 1)->G.Vertex[2].Z;
                    if (nvert == 4) {
                      polynode[3].X = (EleArr + ele - 1)->G.Vertex[3].X;
                      polynode[3].Y = (EleArr + ele - 1)->G.Vertex[3].Y;
                      polynode[3].Z = (EleArr + ele - 1)->G.Vertex[3].Z;
                    }

                    // printf("neBEMChkInPoly for element %d\n", ele);
                    SumOfAngles = neBEMChkInPoly(nvert, polynode, ptintsct);
                    if (fabs(fabs(SumOfAngles) - neBEMtwopi) <= 1.0e-8)
                      InEle = 1;
                    if (debugFn) {
                      // print polynode and InEle
                      printf("Ele: %d\n", ele);
                      printf("ptintsct: %lg, %lg, %lg\n", ptintsct.X,
                             ptintsct.Y, ptintsct.Z);
                      printf("nvert: %d\n", nvert);
                      printf("polynode0: %lg, %lg, %lg\n", polynode[0].X,
                             polynode[0].Y, polynode[0].Z);
                      printf("polynode1: %lg, %lg, %lg\n", polynode[1].X,
                             polynode[1].Y, polynode[1].Z);
                      printf("polynode2: %lg, %lg, %lg\n", polynode[2].X,
                             polynode[2].Y, polynode[2].Z);
                      if (nvert == 4) {
                        printf("polynode3: %lg, %lg, %lg\n", polynode[3].X,
                               polynode[3].Y, polynode[3].Z);
                      }
                      printf("SumOfAngles: %lg, InEle: %d\n", SumOfAngles,
                             InEle);
                      fflush(stdout);
                    }
                    if (InEle) {
                      ptintsct.X = (EleArr + ele - 1)->G.Origin.X;
                      ptintsct.Y = (EleArr + ele - 1)->G.Origin.Y;
                      ptintsct.Z = (EleArr + ele - 1)->G.Origin.Z;
                      // Associate this electron to the identified element
                      EleIntsctd = ele;
                      NbChUpEonEle[ele]++;
                      fprintf(fPtEChUpMap, "%d %lg %lg %lg %d %d %d %d\n",
                              electron, ptintsct.X, ptintsct.Y, ptintsct.Z,
                              prim, InPrim, ele, InEle);

                      if (debugFn) {
                        printf("# electron: %d\n", electron);
                        printf("%lg %lg %lg\n", xlbend, ylbend, zlbend);
                        printf("%lg %lg %lg\n", xend, yend, zend);
                        printf("%lg, %lg, %lg\n", ptintsct.X, ptintsct.Y,
                               ptintsct.Z);
                        printf("# Associated primitive: %d\n", prim);
                        printf(
                            "# Associated element and origin: %d, %lg, %lg, "
                            "%lg\n",
                            ele, (EleArr + ele - 1)->G.Origin.X,
                            (EleArr + ele - 1)->G.Origin.Y,
                            (EleArr + ele - 1)->G.Origin.Z);
                        printf("#NbChUpEonEle on element: %d\n",
                               NbChUpEonEle[ele]);
                        fprintf(ftmpEF, "#Element: %d\n", ele);
                        fprintf(ftmpEF, "%lg %lg %lg\n", polynode[0].X,
                                polynode[0].Y, polynode[0].Z);
                        fprintf(ftmpEF, "%lg %lg %lg\n", polynode[1].X,
                                polynode[1].Y, polynode[1].Z);
                        fprintf(ftmpEF, "%lg %lg %lg\n", polynode[2].X,
                                polynode[2].Y, polynode[2].Z);
                        if (nvert == 4) {
                          fprintf(ftmpEF, "%lg %lg %lg\n", polynode[3].X,
                                  polynode[3].Y, polynode[3].Z);
                        }
                        fprintf(ftmpEF, "%lg %lg %lg\n", polynode[0].X,
                                polynode[0].Y, polynode[0].Z);
                        fprintf(ftmpEF, "\n");
                        fflush(stdout);
                      }
                      break;  // desired element has been found!
                    }
                  }  // for all elements on this primitive

                  if (InEle)
                    break;
                  else {
                    neBEMMessage(
                        "Element not identified ... neBEMKnownCharges\n");
                    return -2;
                  }
                }  // if InPrim
              }    // if intersection and no extrasection

              if ((InPrim) && (intersect) && (!extrasect) &&
                  (InEle))  // all satisfied
                break;  // do not check any further primtive for this electron

              // If, after checking all the primitives, no interstion is found
              // valid
              if (prim ==
                  (NbPrimitives))  // end of the list and no intersection
              {
                int nvert;
                Point3D polynode[4];
                int nearestele = ElementBgn[nearestprim];
                double distele = 1.0e6,
                       mindistele = 1.0e6;  // absurdly high value

                if (debugFn) {
                  printf("prim == (NbPrimitives) ... checking nearest ...\n");
                  printf("nearestprim: %d, mindist: %lg\n", nearestprim,
                         mindist);
                }

                if (mindist <= 10.0e-6) {
                  PrimIntsctd = nearestprim;
                  InPrim = 1;
                } else {
                  InPrim = 0;
                  InEle = 0;
                  break;
                }

                for (int ele = ElementBgn[nearestprim];  // check all elements
                     ele <= ElementEnd[nearestprim]; ++ele) {
                  nvert = 0;
                  if ((EleArr + ele - 1)->G.Type == 3) nvert = 3;
                  if ((EleArr + ele - 1)->G.Type == 4) nvert = 4;
                  if (!nvert) {
                    neBEMMessage(
                        "no vertex element! ... neBEMKnownCharges ...\n");
                    return -20;
                  }

                  /*
                  polynode[0].X = (EleArr+ele-1)->G.Vertex[0].X;
                  polynode[0].Y = (EleArr+ele-1)->G.Vertex[0].Y;
                  polynode[0].Z = (EleArr+ele-1)->G.Vertex[0].Z;
                  polynode[1].X = (EleArr+ele-1)->G.Vertex[1].X;
                  polynode[1].Y = (EleArr+ele-1)->G.Vertex[1].Y;
                  polynode[1].Z = (EleArr+ele-1)->G.Vertex[1].Z;
                  polynode[2].X = (EleArr+ele-1)->G.Vertex[2].X;
                  polynode[2].Y = (EleArr+ele-1)->G.Vertex[2].Y;
                  polynode[2].Z = (EleArr+ele-1)->G.Vertex[2].Z;
                  if(nvert == 4)
                          {
                          polynode[3].X = (EleArr+ele-1)->G.Vertex[3].X;
                          polynode[3].Y = (EleArr+ele-1)->G.Vertex[3].Y;
                          polynode[3].Z = (EleArr+ele-1)->G.Vertex[3].Z;
                          }

                  Vector3D v01, v12, elenorm, unitelenorm;
                  v01.X = polynode[1].X - polynode[0].X;
                  v01.Y = polynode[1].Y - polynode[0].Y;
                  v01.Z = polynode[1].Z - polynode[0].Z;
                  v12.X = polynode[2].X - polynode[1].X;
                  v12.Y = polynode[2].Y - polynode[1].Y;
                  v12.Z = polynode[2].Z - polynode[1].Z;
                  elenorm = Vector3DCrossProduct(&v01, &v12);
                  unitelenorm = UnitVector3D(&elenorm);

                  if((nvert == 3) || (nvert == 4))
                          {
                          if(debugFn)
                                  {
                                  printf("nearestprim: %d, element: %d\n",
  nearestprim, ele); printf("vertex0: %lg, %lg, %lg\n", polynode[0].X,
  polynode[0].Y, polynode[0].Z); printf("vertex1: %lg, %lg, %lg\n",
                                                  polynode[1].X, polynode[1].Y,
  polynode[1].Z); printf("vertex2: %lg, %lg, %lg\n", polynode[2].X,
  polynode[2].Y, polynode[2].Z); if(PrimType[prim] == 4)
                                          {
                                          printf("vertex3: %lg, %lg, %lg\n",
                                                  polynode[3].X, polynode[3].Y,
  polynode[3].Z);
                                          }
                                  fprintf(ftmpEF, "#nearestprim: %d, element:
  %d\n", nearestprim, ele); fprintf(ftmpEF, "%lg %lg %lg\n", polynode[0].X,
  polynode[0].Y, polynode[0].Z); fprintf(ftmpEF, "%lg %lg %lg\n", polynode[1].X,
  polynode[1].Y, polynode[1].Z); fprintf(ftmpEF, "%lg %lg %lg\n", polynode[2].X,
  polynode[2].Y, polynode[2].Z); if(PrimType[prim] == 4)
                                          {
                                          fprintf(ftmpEF, "%lg %lg %lg\n",
                                                  polynode[3].X, polynode[3].Y,
  polynode[3].Z);
                                          }
                                  fprintf(ftmpEF, "%lg %lg %lg\n",
                                                  polynode[0].X, polynode[0].Y,
  polynode[0].Z); fprintf(ftmpEF, "\n"); fflush(stdout); }	// debugFn

  // use a, b, c (normal is ai + bj + ck) at one of the nodes to get d
  // ax + by + cz + d = 0 is the equation of the plane
                          double a = unitelenorm.X;
                          double b = unitelenorm.Y;
                          double c = unitelenorm.Z;
                          double d = - a*polynode[0].X - b*polynode[0].Y
                                                                          -
  c*polynode[0].Z;

                          // distance of the end point to this primitve is
                          distele = (xend*a + yend*b + zend*c + d)
                                                                                  / sqrt(a*a + b*b + c*c);
                          distele = fabs(distele);	// if only magnitude is
  required
                          */

                  Vector3D eleOrigin;
                  eleOrigin.X = (EleArr + ele - 1)->G.Origin.X;
                  eleOrigin.Y = (EleArr + ele - 1)->G.Origin.Y;
                  eleOrigin.Z = (EleArr + ele - 1)->G.Origin.Z;
                  distele = (eleOrigin.X - xend) * (eleOrigin.X - xend) +
                            (eleOrigin.Y - yend) * (eleOrigin.Y - yend) +
                            (eleOrigin.Z - zend) * (eleOrigin.Z - zend);
                  distele = sqrt(distele);

                  if (ele == ElementBgn[nearestprim]) {
                    mindistele = distele;
                    nearestele = ele;
                  }
                  if (distele < mindistele) {
                    mindistele = distele;
                    nearestele = ele;
                  }

                  if (debugFn) {
                    // printf("a, b, c, d, dist: %lg, %lg, %lg, %lg, %lg\n",
                    // a, b, c, d,  dist);
                    // printf("vector n: ai + bj + ck\n");
                    // printf("vector v: xgrd, ygrd, zgrd: %lg, %lg, %lg\n",
                    // xgrd, ygrd, zgrd);
                    printf(
                        "distele: %lg, mindistele: %lg,from nearest ele "
                        "origin: %d\n",
                        distele, mindistele, nearestele);
                    fflush(stdout);
                  }

                  // }	// if PrimType is 3 or 4
                }  // for elements in nearestprim

                // if(mindistele <= 10.0e-6)
                // {
                EleIntsctd = nearestele;
                InEle = 1;
                ptintsct.X = (EleArr + EleIntsctd - 1)->G.Origin.X;
                ptintsct.Y = (EleArr + EleIntsctd - 1)->G.Origin.Y;
                ptintsct.Z = (EleArr + EleIntsctd - 1)->G.Origin.Z;
                NbChUpEonEle[EleIntsctd]++;

                fprintf(fPtEChUpMap, "%d %lg %lg %lg %d %d %d %d\n", electron,
                        ptintsct.X, ptintsct.Y, ptintsct.Z, PrimIntsctd, InPrim,
                        EleIntsctd, InEle);
                // }	// if mindistele

                if (debugFn) {
                  printf("# electron: %d\n", electron);
                  printf("%lg %lg %lg\n", xlbend, ylbend, zlbend);
                  printf("%lg %lg %lg\n", xend, yend, zend);
                  printf("%lg, %lg, %lg\n", ptintsct.X, ptintsct.Y, ptintsct.Z);
                  printf("# Associated primitive: %d\n", PrimIntsctd);
                  printf("# Associated element and origin: %d, %lg, %lg, %lg\n",
                         EleIntsctd, (EleArr + EleIntsctd - 1)->G.Origin.X,
                         (EleArr + EleIntsctd - 1)->G.Origin.Y,
                         (EleArr + EleIntsctd - 1)->G.Origin.Z);
                  printf("#NbChUpEonEle on element: %d\n",
                         NbChUpEonEle[EleIntsctd]);
                  fflush(stdout);

                  fprintf(ftmpEF, "#Element: %d\n", EleIntsctd);
                  polynode[0].X = (EleArr + EleIntsctd - 1)->G.Vertex[0].X;
                  polynode[0].Y = (EleArr + EleIntsctd - 1)->G.Vertex[0].Y;
                  polynode[0].Z = (EleArr + EleIntsctd - 1)->G.Vertex[0].Z;
                  polynode[1].X = (EleArr + EleIntsctd - 1)->G.Vertex[1].X;
                  polynode[1].Y = (EleArr + EleIntsctd - 1)->G.Vertex[1].Y;
                  polynode[1].Z = (EleArr + EleIntsctd - 1)->G.Vertex[1].Z;
                  polynode[2].X = (EleArr + EleIntsctd - 1)->G.Vertex[2].X;
                  polynode[2].Y = (EleArr + EleIntsctd - 1)->G.Vertex[2].Y;
                  polynode[2].Z = (EleArr + EleIntsctd - 1)->G.Vertex[2].Z;
                  if (nvert == 4) {
                    polynode[3].X = (EleArr + EleIntsctd - 1)->G.Vertex[3].X;
                    polynode[3].Y = (EleArr + EleIntsctd - 1)->G.Vertex[3].Y;
                    polynode[3].Z = (EleArr + EleIntsctd - 1)->G.Vertex[3].Z;
                  }
                  fprintf(ftmpEF, "%lg %lg %lg\n", polynode[0].X, polynode[0].Y,
                          polynode[0].Z);
                  fprintf(ftmpEF, "%lg %lg %lg\n", polynode[1].X, polynode[1].Y,
                          polynode[1].Z);
                  fprintf(ftmpEF, "%lg %lg %lg\n", polynode[2].X, polynode[2].Y,
                          polynode[2].Z);
                  if (nvert == 4) {
                    fprintf(ftmpEF, "%lg %lg %lg\n", polynode[3].X,
                            polynode[3].Y, polynode[3].Z);
                  }
                  fprintf(ftmpEF, "%lg %lg %lg\n", polynode[0].X, polynode[0].Y,
                          polynode[0].Z);
                  fprintf(ftmpEF, "\n");
                }  // debug
              }    // if prim == NbPrimitives

            }  // for all primitives // just not those on the volume

            if (debugFn)
              printf("writing file for checking electron positions\n");

            if (debugFn)  // check electron positions, volume primitives and
                          // elements
            {
              char elecposdbg[256], enbstr[10];
              sprintf(enbstr, "%d", electron);
              strcpy(elecposdbg, "/tmp/Electron");
              strcat(elecposdbg, enbstr);
              strcat(elecposdbg, ".out");
              FILE *fepd = fopen(elecposdbg, "w");
              if (fepd == NULL) {
                printf(
                    "cannot open writable file to debug electron positions "
                    "...\n");
                printf("returning ...\n");
                return -111;
              }
              // write electron number, end, lbend, volume, primitive, elements,
              // intxn
              fprintf(fepd, "#electron: %d %d\n", enb,
                      electron);  // should print same
              fprintf(fepd, "#last but end position:\n");
              fprintf(fepd, "%lg %lg %lg\n", xlbend, ylbend, zlbend);
              fprintf(fepd, "#end position:\n");
              fprintf(fepd, "%lg %lg %lg\n\n", xend, yend, zend);
              fprintf(fepd, "#intersected primitive number: %d\n", PrimIntsctd);
              if (PrimIntsctd >= 1) {
                fprintf(fepd, "#PrimType: %d\n", PrimType[PrimIntsctd]);
                fprintf(fepd, "#prim vertices:\n");
                fprintf(fepd, "%lg %lg %lg\n", XVertex[PrimIntsctd][0],
                        YVertex[PrimIntsctd][0], ZVertex[PrimIntsctd][0]);
                fprintf(fepd, "%lg %lg %lg\n", XVertex[PrimIntsctd][1],
                        YVertex[PrimIntsctd][1], ZVertex[PrimIntsctd][1]);
                fprintf(fepd, "%lg %lg %lg\n", XVertex[PrimIntsctd][2],
                        YVertex[PrimIntsctd][2], ZVertex[PrimIntsctd][2]);
                if (PrimType[PrimIntsctd] == 4) {
                  fprintf(fepd, "%lg %lg %lg\n", XVertex[PrimIntsctd][3],
                          YVertex[PrimIntsctd][3], ZVertex[PrimIntsctd][3]);
                }
                fprintf(fepd, "%lg %lg %lg\n", XVertex[PrimIntsctd][0],
                        YVertex[PrimIntsctd][0], ZVertex[PrimIntsctd][0]);
                fprintf(fepd, "\n");

                fprintf(fepd, "#ptintsct:\n");
                fprintf(fepd, "%lg %lg %lg\n", ptintsct.X, ptintsct.Y,
                        ptintsct.Z);
                fprintf(fepd, "\n");
              }
              fprintf(fepd, "#intersected element number: %d\n", EleIntsctd);
              if (EleIntsctd >= 1) {
                int gtype = (EleArr + EleIntsctd - 1)->G.Type;
                fprintf(fepd, "#EleType: %d\n", gtype);
                fprintf(fepd, "#element vertices:\n");
                double x0 = (EleArr + EleIntsctd - 1)->G.Vertex[0].X;
                double y0 = (EleArr + EleIntsctd - 1)->G.Vertex[0].Y;
                double z0 = (EleArr + EleIntsctd - 1)->G.Vertex[0].Z;
                double x1 = (EleArr + EleIntsctd - 1)->G.Vertex[1].X;
                double y1 = (EleArr + EleIntsctd - 1)->G.Vertex[1].Y;
                double z1 = (EleArr + EleIntsctd - 1)->G.Vertex[1].Z;
                double x2 = (EleArr + EleIntsctd - 1)->G.Vertex[2].X;
                double y2 = (EleArr + EleIntsctd - 1)->G.Vertex[2].Y;
                double z2 = (EleArr + EleIntsctd - 1)->G.Vertex[2].Z;
                fprintf(fepd, "%lg %lg %lg\n", x0, y0, z0);
                fprintf(fepd, "%lg %lg %lg\n", x1, y1, z1);
                fprintf(fepd, "%lg %lg %lg\n", x2, y2, z2);
                if (gtype == 4) {
                  double x3 = (EleArr + EleIntsctd - 1)->G.Vertex[3].X;
                  double y3 = (EleArr + EleIntsctd - 1)->G.Vertex[3].Y;
                  double z3 = (EleArr + EleIntsctd - 1)->G.Vertex[3].Z;
                  fprintf(fepd, "%lg %lg %lg\n", x3, y3, z3);
                }
                fprintf(fepd, "%lg %lg %lg\n", x0, y0, z0);
                fprintf(fepd, "\n");

                fprintf(fepd, "#ptintsct:\n");
                fprintf(fepd, "%lg %lg %lg\n", ptintsct.X, ptintsct.Y,
                        ptintsct.Z);
                fprintf(fepd, "\n");
              }

              fclose(fepd);
            }  // if 1
            if (debugFn)
              printf("done writing file for checking electron positions\n");
          }  // for all the electrons
          fclose(fPtEChUpMap);
          if (debugFn) printf("done for all electrons\n");

          FILE *fEleEChUpMap = fopen("EleEChUpMap.out", "w");
          if (fEleEChUpMap == NULL) {
            printf("cannot open EleEChUpMap.out file for writing ...\n");
            return 111;
          }
          for (int ele = 1; ele <= NbElements; ++ele) {
            (EleArr + ele - 1)->Assigned +=
                ChUpFactor * Q_E * NbChUpEonEle[ele] / (EleArr + ele - 1)->G.dA;
            fprintf(fEleEChUpMap, "%d %lg %lg %lg %d %lg\n", ele,
                    (EleArr + ele - 1)->G.Origin.X,
                    (EleArr + ele - 1)->G.Origin.Y,
                    (EleArr + ele - 1)->G.Origin.Z, NbChUpEonEle[ele],
                    (EleArr + ele - 1)->Assigned);
          }
          fclose(fEleEChUpMap);

          fclose(ftmpEF);
          free(NbChUpEonEle);
        }  // Calculation for electrons ends

        {  // Calculation for ions
          FILE *fptrChargingUpIFile = fopen(ChargingUpIFile, "r");
          if (fptrChargingUpIFile == NULL) {
            neBEMMessage("ChargingUpI file absent ... returning\n");
            return -10;
          }
          int NbOfI = neBEMGetNbOfLines(ChargingUpIFile);
          if (NbOfI <= 1) {
            neBEMMessage("Too few lines in ChargingUpI ... returning\n");
            return -11;
          } else {  // initialize
            NbChUpIonEle = (int *)malloc((NbElements + 1) * sizeof(int));
            for (int ele = 0; ele <= NbElements; ++ele) {  
              // CHECK!!! ele limit starts from 0 but all other from 1 to ...
              NbChUpIonEle[ele] = 0;
            }
          }

          // read the header line
          char header[256];
          fgets(header, 256, fptrChargingUpIFile);

          --NbOfI;  // one line was for the header
          if (debugFn) printf("No. of ions: %d\n", NbOfI);
          char tmpIFile[256];
          strcpy(tmpIFile, "/tmp/IonTempFile.out");
          FILE *ftmpIF = fopen(tmpIFile, "w");
          if (ftmpIF == NULL) {
            printf("cannot open temporary ion output file ... returning ...\n");
            return -100;
          }
          FILE *fPtIChUpMap = fopen("PtIChUpMap.out", "w");
          if (fPtIChUpMap == NULL) {
            printf("cannot open PtIChUpMap.out file for writing ...\n");
            return 110;
          }

          char label;
          int inb, vol;  // label, volume and ion number
          double xlbend, ylbend, zlbend, xend, yend,
              zend;          // lbend == Last But END
          Point3D ptintsct;  // each ion likely to have an intersection point
          for (int ion = 1; ion <= NbOfI; ++ion) {
            fscanf(fptrChargingUpIFile, "%c %d %d %lg %lg %lg %lg %lg %lg\n",
                   &label, &vol, &inb, &xlbend, &xend, &ylbend, &yend, &zlbend,
                   &zend);
            xlbend /= LScaleFactor;  // cm to metre
            xend /= LScaleFactor;    // and similar other
            ylbend /= LScaleFactor;
            yend /= LScaleFactor;
            zlbend /= LScaleFactor;
            zend /= LScaleFactor;
            ptintsct.X = 0.0;
            ptintsct.Y = 0.0;
            ptintsct.Z = 0.0;  // initialize

            // find the parametric equation of this last segment
            // if xend < xlbend, swap the directions
            // This has not been mentioned as mandatory in Vince's book
            // "Geometry for Computer Graphics", but implied in the book "A
            // Programmer's Geometry"
            if (xend < xlbend)  // not implemented now
            {
            }
            double lseg = (xend - xlbend) * (xend - xlbend) +
                          (yend - ylbend) * (yend - ylbend) +
                          (zend - zlbend) * (zend - zlbend);
            lseg = sqrt(lseg);
            double xgrd =
                (xend - xlbend) / lseg;  // normalized direction vector
            double ygrd = (yend - ylbend) / lseg;
            double zgrd = (zend - zlbend) / lseg;
            if (debugFn) {
              printf("\nion: %d\n", ion);
              printf("xlbend: %lg, ylbend: %lg, zlbend: %lg\n", xlbend, ylbend,
                     zlbend);
              printf("xend: %lg, yend: %lg, zend: %lg\n", xend, yend, zend);
              printf("xgrd: %lg, ygrd: %lg, zgrd: %lg\n", xgrd, ygrd, zgrd);
              fprintf(ftmpIF, "#e: %d, label: %c, vol: %d\n", ion, label, vol);
              fprintf(ftmpIF, "%lg %lg %lg\n", xlbend, ylbend, zlbend);
              fprintf(ftmpIF, "%lg %lg %lg\n", xend, yend, zend);
              fprintf(ftmpIF, "#xgrd: %lg, ygrd: %lg, zgrd: %lg\n", xgrd, ygrd,
                      zgrd);
              fprintf(ftmpIF, "\n");
            }

            // determine which element gets this electron
            // At present, the logic is as follow:
            // Using the information on the last segment, find out which
            // primitive is pierced by it and at which point From the
            // intersection point, find out the element number on the primitive
            // in a local sense Using the information (start and end elements of
            // a given primitive) identify the element in a global sense. This
            // approach should be lot more efficient than checking intersection
            // element by element.
            // The intersection point is computed following algorithm
            // implemented in a matlab code (plane_imp_line_par_int_3d.m) Also
            // check which primitive in the list is the closet to the end point

            int PrimIntsctd = -1,
                EleIntsctd = -1;   // intersected primitive & element
            int nearestprim = -1;  // absurd value
            double dist = 1.0e6, mindist = 1.0e6;  // absurdly high numbers
            double SumOfAngles;
            // check all primitives
            for (int prim = 1; prim <= NbPrimitives; ++prim) { 
              if (InterfaceType[prim] != 4)
                continue;  // primitive not a dielectric

              int intersect = 0, extrasect = 1;  // worst of conditions
              int InPrim = 0, InEle = 0;
              if (debugFn)
                printf("prim: %d, mindist: %lg, nearestprim: %d\n", prim,
                       mindist, nearestprim);

              // get the primitive nodes

              // Use two nodes at a time to get two independent vectors on
              // primitive Get cross-product of these two vector - normal to the
              // plane Note that the normal is already associated with the
              // primitive of type 3 and 4; 2 is wire and does not have any
              // associated normal
              if ((PrimType[prim] == 3) || (PrimType[prim] == 4)) {
                if (debugFn) {
                  printf("prim: %d\n", prim);
                  printf("vertex0: %lg, %lg, %lg\n", XVertex[prim][0],
                         YVertex[prim][0], ZVertex[prim][0]);
                  printf("vertex1: %lg, %lg, %lg\n", XVertex[prim][1],
                         YVertex[prim][1], ZVertex[prim][1]);
                  printf("vertex2: %lg, %lg, %lg\n", XVertex[prim][2],
                         YVertex[prim][2], ZVertex[prim][2]);
                  if (PrimType[prim] == 4) {
                    printf("vertex3: %lg, %lg, %lg\n", XVertex[prim][3],
                           YVertex[prim][3], ZVertex[prim][3]);
                  }
                  fprintf(ftmpIF, "#prim: %d\n", prim);
                  fprintf(ftmpIF, "%lg %lg %lg\n", XVertex[prim][0],
                          YVertex[prim][0], ZVertex[prim][0]);
                  fprintf(ftmpIF, "%lg %lg %lg\n", XVertex[prim][1],
                          YVertex[prim][1], ZVertex[prim][1]);
                  fprintf(ftmpIF, "%lg %lg %lg\n", XVertex[prim][2],
                          YVertex[prim][2], ZVertex[prim][2]);
                  if (PrimType[prim] == 4) {
                    fprintf(ftmpIF, "%lg %lg %lg\n", XVertex[prim][3],
                            YVertex[prim][3], ZVertex[prim][3]);
                  }
                  fprintf(ftmpIF, "%lg %lg %lg\n", XVertex[prim][0],
                          YVertex[prim][0], ZVertex[prim][0]);
                  fprintf(ftmpIF, "\n");
                  fflush(stdout);
                }  // debugFn

                // use a, b, c (normal is ai + bj + ck) at one of the nodes to
                // get d ax + by + cz + d = 0 is the equation of the plane
                double d = -XNorm[prim] * XVertex[prim][0] -
                           YNorm[prim] * YVertex[prim][0] -
                           ZNorm[prim] * ZVertex[prim][0];

                // distance of the end point to this primitve is
                dist =
                    (xend * XNorm[prim] + yend * YNorm[prim] +
                     zend * ZNorm[prim] + d) /
                    sqrt(XNorm[prim] * XNorm[prim] + YNorm[prim] * YNorm[prim] +
                         ZNorm[prim] * ZNorm[prim]);
                dist = fabs(dist);  // if only magnitude is required
                if (prim == 1) {
                  mindist = dist;
                  nearestprim = prim;
                }
                if ((prim == 1) && debugFn)
                  printf(
                      "after prim == 1 check mindist: %lg, nearestprim: %d\n",
                      mindist, nearestprim);

                // Point of intersection
                // Algo as on p62 (pdf 81) of Vince - Geometry for Computer
                // Graphics 1.5.13 Intersection of a line and a plane Algorithm
                // as implemented in plne_imp_line_par_int_3d.m a (nx), b (ny),
                // c (nz), d are a, b, c, d vx, vy, vz are xgrd, ygrd and zgrd
                // tx, ty, tz are xlbend, ylbend, zlbend
                // In the present case, n and v are unit vectors
                double a = XNorm[prim];
                double b = YNorm[prim];
                double c = ZNorm[prim];
                double norm1 = sqrt(a * a + b * b + c * c);
                double norm2 = sqrt(xgrd * xgrd + ygrd * ygrd + zgrd * zgrd);
                double denom =
                    a * xgrd + b * ygrd + c * zgrd;  // (vec)n.(vec)v; if 0, ||
                double tol =
                    1.0e-12;  // CHECK: -8 in original code; sizes small here
                intersect = extrasect = 0;

                if (debugFn) {
                  printf("a, b, c, d, dist: %lg, %lg, %lg, %lg, %lg\n", a, b, c,
                         d, dist);
                  printf("vector n: ai + bj + ck\n");
                  printf("vector v: xgrd, ygrd, zgrd: %lg, %lg, %lg\n", xgrd,
                         ygrd, zgrd);
                  printf("norm1, norm2, (vec n . vec v) denom: %lg, %lg, %lg\n",
                         norm1, norm2, denom);
                  printf("if vec n . vec v == 0, line and plane parallel\n");
                  fflush(stdout);
                }

                if (fabs(denom) <
                    tol * norm1 * norm2)  // line parallel to the plane
                {
                  if (a * xlbend + b * ylbend + c * zlbend + d ==
                      0.0)  // CHECK == for float
                  {
                    intersect = 1;
                    extrasect = 0;
                    ptintsct.X = xlbend;
                    ptintsct.Y = ylbend;
                    ptintsct.Z = zlbend;
                  } else {
                    intersect = 0;
                    extrasect = 0;
                    ptintsct.X = 0.0;  // Wrong to assign 0 values
                    ptintsct.Y =
                        0.0;  // However, they are never going to be used
                    ptintsct.Z = 0.0;  // since intersect is 0
                  }
                  if (debugFn) {
                    printf("line and plane parallel ...\n");
                    printf("intersect: %d, extrasect: %d\n", intersect,
                           extrasect);
                    printf("intersection point: %lg, %lg, %lg\n", ptintsct.X,
                           ptintsct.Y, ptintsct.Z);
                  }       // if line and plane are parallel
                } else {  // if they are not parallel, they must intersect
                  intersect = 1;
                  double t =
                      -(a * xlbend + b * ylbend + c * zlbend + d) / denom;

                  // check whether t is less than the length of the segment
                  // and in the correct direction
                  // If not, then an extrapolated intersection is not of
                  // interest
                  if ((t < 0.0) || (fabs(t) > fabs(lseg))) {
                    extrasect = 1;
                    ptintsct.X = xlbend + t * xgrd;
                    ptintsct.Y = ylbend + t * ygrd;
                    ptintsct.Z = zlbend + t * zgrd;
                  } else {
                    extrasect = 0;
                    ptintsct.X = xlbend + t * xgrd;
                    ptintsct.Y = ylbend + t * ygrd;
                    ptintsct.Z = zlbend + t * zgrd;
                  }
                  if (debugFn) {
                    printf("line and plane NOT parallel ...\n");
                    printf("intersect: %d, extrasect: %d\n", intersect,
                           extrasect);
                    printf("intersection point: %lg, %lg, %lg\n", ptintsct.X,
                           ptintsct.Y, ptintsct.Z);
                    printf("t, lseg: %lg, %lg\n", t, lseg);
                    printf(
                        "for an interesting intersection, lseg > t > 0.0 "
                        "...\n\n");
                    fflush(stdout);
                  }   // must intersect
                }     // if not parallel
              }       // if PrimType is 3 or 4
              else {  // this is a wire primitive - assume no charging up issues
                dist = -1.0;  // an absurd negative distance
                intersect = 0;
                extrasect = 0;
              }  // else PrimType 3 or 4

              if (dist < mindist) {
                mindist = dist;
                nearestprim = prim;
              }

              // implicit assumption: the first primitive that gets pierced by
              // the ray is the one that we want. There can be other primitives
              // that are pierced by the same ray. So, this logic should be
              // refined further
              if ((intersect == 1) && (extrasect == 0)) {
                // check whether the intersection point is within primitive
                // polygon
                int nvert = PrimType[prim];
                Point3D polynode[4];
                polynode[0].X = XVertex[prim][0];
                polynode[0].Y = YVertex[prim][0];
                polynode[0].Z = ZVertex[prim][0];
                polynode[1].X = XVertex[prim][1];
                polynode[1].Y = YVertex[prim][1];
                polynode[1].Z = ZVertex[prim][1];
                polynode[2].X = XVertex[prim][2];
                polynode[2].Y = YVertex[prim][2];
                polynode[2].Z = ZVertex[prim][2];
                if (PrimType[prim] == 4) {
                  polynode[3].X = XVertex[prim][3];
                  polynode[3].Y = YVertex[prim][3];
                  polynode[3].Z = ZVertex[prim][3];
                }
                // printf("neBEMChkInPoly for primitive %d\n", prim);
                SumOfAngles = neBEMChkInPoly(nvert, polynode, ptintsct);
                if (fabs(fabs(SumOfAngles) - neBEMtwopi) <= 1.0e-8) {
                  InPrim = 1;
                  PrimIntsctd = prim;
                }
                if (debugFn) {
                  // print polynode and InPrim
                  printf("Prim: %d\n", prim);
                  printf("ptintsct: %lg, %lg, %lg\n", ptintsct.X, ptintsct.Y,
                         ptintsct.Z);
                  printf("nvert: %d\n", nvert);
                  printf("polynode0: %lg, %lg, %lg\n", polynode[0].X,
                         polynode[0].Y, polynode[0].Z);
                  printf("polynode1: %lg, %lg, %lg\n", polynode[1].X,
                         polynode[1].Y, polynode[1].Z);
                  printf("polynode2: %lg, %lg, %lg\n", polynode[2].X,
                         polynode[2].Y, polynode[2].Z);
                  if (nvert == 4) {
                    printf("polynode3: %lg, %lg, %lg\n", polynode[3].X,
                           polynode[3].Y, polynode[3].Z);
                  }
                  printf("SumOfAngles: %lg, InPrim: %d\n", SumOfAngles, InPrim);
                  fflush(stdout);
                }
                if (!InPrim) continue;  // check next primitive

                // Once identified, check in which element belonging to this
                // primitive contains the point of intersection
                InEle = 0;
                for (int ele = ElementBgn[prim]; ele <= ElementEnd[prim];
                     ++ele) {
                  nvert = 0;
                  if ((EleArr + ele - 1)->G.Type == 3) nvert = 3;
                  if ((EleArr + ele - 1)->G.Type == 4) nvert = 4;
                  if (!nvert) {
                    neBEMMessage(
                        "no vertex in element! ... neBEMKnownCharges ...\n");
                    return -20;
                  }

                  polynode[0].X = (EleArr + ele - 1)->G.Vertex[0].X;
                  polynode[0].Y = (EleArr + ele - 1)->G.Vertex[0].Y;
                  polynode[0].Z = (EleArr + ele - 1)->G.Vertex[0].Z;
                  polynode[1].X = (EleArr + ele - 1)->G.Vertex[1].X;
                  polynode[1].Y = (EleArr + ele - 1)->G.Vertex[1].Y;
                  polynode[1].Z = (EleArr + ele - 1)->G.Vertex[1].Z;
                  polynode[2].X = (EleArr + ele - 1)->G.Vertex[2].X;
                  polynode[2].Y = (EleArr + ele - 1)->G.Vertex[2].Y;
                  polynode[2].Z = (EleArr + ele - 1)->G.Vertex[2].Z;
                  if (nvert == 4) {
                    polynode[3].X = (EleArr + ele - 1)->G.Vertex[3].X;
                    polynode[3].Y = (EleArr + ele - 1)->G.Vertex[3].Y;
                    polynode[3].Z = (EleArr + ele - 1)->G.Vertex[3].Z;
                  }

                  // printf("neBEMChkInPoly for element %d\n", ele);
                  SumOfAngles = neBEMChkInPoly(nvert, polynode, ptintsct);
                  if (fabs(fabs(SumOfAngles) - neBEMtwopi) <= 1.0e-8) InEle = 1;
                  if (debugFn) {
                    // print polynode and InEle
                    printf("Ele: %d\n", ele);
                    printf("ptintsct: %lg, %lg, %lg\n", ptintsct.X, ptintsct.Y,
                           ptintsct.Z);
                    printf("nvert: %d\n", nvert);
                    printf("polynode0: %lg, %lg, %lg\n", polynode[0].X,
                           polynode[0].Y, polynode[0].Z);
                    printf("polynode1: %lg, %lg, %lg\n", polynode[1].X,
                           polynode[1].Y, polynode[1].Z);
                    printf("polynode2: %lg, %lg, %lg\n", polynode[2].X,
                           polynode[2].Y, polynode[2].Z);
                    if (nvert == 4) {
                      printf("polynode3: %lg, %lg, %lg\n", polynode[3].X,
                             polynode[3].Y, polynode[3].Z);
                    }
                    printf("SumOfAngles: %lg, InEle: %d\n", SumOfAngles, InEle);
                    fflush(stdout);
                  }
                  if (InEle) {
                    ptintsct.X = (EleArr + ele - 1)->G.Origin.X;
                    ptintsct.Y = (EleArr + ele - 1)->G.Origin.Y;
                    ptintsct.Z = (EleArr + ele - 1)->G.Origin.Z;
                    EleIntsctd = ele;
                    // Associate this electron to the identified element
                    NbChUpIonEle[ele]++;
                    fprintf(fPtIChUpMap, "%d %lg %lg %lg %d %d %d %d\n", ion,
                            ptintsct.X, ptintsct.Y, ptintsct.Z, prim, InPrim,
                            ele, InEle);

                    if (debugFn) {
                      printf("# ion: %d\n", ion);
                      printf("%lg %lg %lg\n", xlbend, ylbend, zlbend);
                      printf("%lg %lg %lg\n", xend, yend, zend);
                      printf("%lg, %lg, %lg\n", ptintsct.X, ptintsct.Y,
                             ptintsct.Z);
                      printf("# Associated primitive: %d\n", prim);
                      printf(
                          "# Associated element and origin: %d, %lg, %lg, "
                          "%lg\n",
                          ele, (EleArr + ele - 1)->G.Origin.X,
                          (EleArr + ele - 1)->G.Origin.Y,
                          (EleArr + ele - 1)->G.Origin.Z);
                      printf("#NbChUpIonEle on element: %d\n",
                             NbChUpIonEle[ele]);
                      fprintf(ftmpIF, "#Element: %d\n", ele);
                      fprintf(ftmpIF, "%lg %lg %lg\n", polynode[0].X,
                              polynode[0].Y, polynode[0].Z);
                      fprintf(ftmpIF, "%lg %lg %lg\n", polynode[1].X,
                              polynode[1].Y, polynode[1].Z);
                      fprintf(ftmpIF, "%lg %lg %lg\n", polynode[2].X,
                              polynode[2].Y, polynode[2].Z);
                      if (nvert == 4) {
                        fprintf(ftmpIF, "%lg %lg %lg\n", polynode[3].X,
                                polynode[3].Y, polynode[3].Z);
                      }
                      fprintf(ftmpIF, "%lg %lg %lg\n", polynode[0].X,
                              polynode[0].Y, polynode[0].Z);
                      fprintf(ftmpIF, "\n");
                      fflush(stdout);
                    }
                    break;  // desired element has been found!
                  }         // if InEle
                }           // for all elements on this primitive

                if (InEle)
                  break;
                else {
                  neBEMMessage(
                      "Element cannot be identified ... neBEMKnownCharges\n");
                  return -2;
                }
              }  // if proper intersection and no extrasection

              if ((InPrim) && (intersect) && (!extrasect) &&
                  (InEle))  // all satisfied
                break;  // do not check any further primtive for this electron

              // If, after checking all the primitives, no interstion is found
              // valid
              if (prim ==
                  (NbPrimitives))  // end of the list and no intersection
              {
                int nvert;
                Point3D polynode[4];
                int nearestele = ElementBgn[nearestprim];
                double distele = 1.0e6,
                       mindistele = 1.0e6;  // absurdly high value

                if (debugFn) {
                  printf("prim == (NbPrimitives) ... checking nearest ...\n");
                  printf("nearestprim: %d, mindist: %lg\n", nearestprim,
                         mindist);
                }

                if (mindist <= 10.0e-6) {
                  PrimIntsctd = nearestprim;
                  InPrim = 1;
                } else {
                  InPrim = 0;
                  InEle = 0;
                  break;
                }

                for (int ele = ElementBgn[nearestprim];  // check all elements
                     ele <= ElementEnd[nearestprim]; ++ele) {
                  nvert = 0;
                  if ((EleArr + ele - 1)->G.Type == 3) nvert = 3;
                  if ((EleArr + ele - 1)->G.Type == 4) nvert = 4;
                  if (!nvert) {
                    neBEMMessage(
                        "no vertex element! ... neBEMKnownCharges ...\n");
                    return -20;
                  }

                  /*
                  polynode[0].X = (EleArr+ele-1)->G.Vertex[0].X;
                  polynode[0].Y = (EleArr+ele-1)->G.Vertex[0].Y;
                  polynode[0].Z = (EleArr+ele-1)->G.Vertex[0].Z;
                  polynode[1].X = (EleArr+ele-1)->G.Vertex[1].X;
                  polynode[1].Y = (EleArr+ele-1)->G.Vertex[1].Y;
                  polynode[1].Z = (EleArr+ele-1)->G.Vertex[1].Z;
                  polynode[2].X = (EleArr+ele-1)->G.Vertex[2].X;
                  polynode[2].Y = (EleArr+ele-1)->G.Vertex[2].Y;
                  polynode[2].Z = (EleArr+ele-1)->G.Vertex[2].Z;
                  if(nvert == 4)
                          {
                          polynode[3].X = (EleArr+ele-1)->G.Vertex[3].X;
                          polynode[3].Y = (EleArr+ele-1)->G.Vertex[3].Y;
                          polynode[3].Z = (EleArr+ele-1)->G.Vertex[3].Z;
                          }

                  Vector3D v01, v12, elenorm, unitelenorm;
                  v01.X = polynode[1].X - polynode[0].X;
                  v01.Y = polynode[1].Y - polynode[0].Y;
                  v01.Z = polynode[1].Z - polynode[0].Z;
                  v12.X = polynode[2].X - polynode[1].X;
                  v12.Y = polynode[2].Y - polynode[1].Y;
                  v12.Z = polynode[2].Z - polynode[1].Z;
                  elenorm = Vector3DCrossProduct(&v01, &v12);
                  unitelenorm = UnitVector3D(&elenorm);

                  if((nvert == 3) || (nvert == 4))
                          {
                          if(debugFn)
                                  {
                                  printf("nearestprim: %d, element: %d\n",
          nearestprim, ele); printf("vertex0: %lg, %lg, %lg\n", polynode[0].X,
          polynode[0].Y, polynode[0].Z); printf("vertex1: %lg, %lg, %lg\n",
                                                          polynode[1].X,
          polynode[1].Y, polynode[1].Z); printf("vertex2: %lg, %lg, %lg\n",
                                                  polynode[2].X, polynode[2].Y,
          polynode[2].Z); if(PrimType[prim] == 4)
                                          {
                                          printf("vertex3: %lg, %lg, %lg\n",
                                                          polynode[3].X,
          polynode[3].Y, polynode[3].Z);
                                          }
                                  fprintf(ftmpIF, "#nearestprim: %d, element:
          %d\n", nearestprim, ele); fprintf(ftmpIF, "%lg %lg %lg\n",
                                                          polynode[0].X,
          polynode[0].Y, polynode[0].Z); fprintf(ftmpIF, "%lg %lg %lg\n",
                                                          polynode[1].X,
          polynode[1].Y, polynode[1].Z); fprintf(ftmpIF, "%lg %lg %lg\n",
                                                          polynode[2].X,
          polynode[2].Y, polynode[2].Z); if(PrimType[prim] == 4)
                                          {
                                          fprintf(ftmpIF, "%lg %lg %lg\n",
                                                  polynode[3].X, polynode[3].Y,
          polynode[3].Z);
                                          }
                                  fprintf(ftmpIF, "%lg %lg %lg\n",
                                                          polynode[0].X,
          polynode[0].Y, polynode[0].Z); fprintf(ftmpIF, "\n"); fflush(stdout);
                                  }	// debugFn

          // use a, b, c (normal is ai + bj + ck) at one of the nodes to get d
          // ax + by + cz + d = 0 is the equation of the plane
                          double a = unitelenorm.X;
                          double b = unitelenorm.Y;
                          double c = unitelenorm.Z;
                          double d = - unitelenorm.X*polynode[0].X
                                                                                  - unitelenorm.Y*polynode[0].Y
                                                                                  - unitelenorm.Z*polynode[0].Z;

                          // distance of the end point to this primitve is
                          distele = (xend * unitelenorm.X + yend * unitelenorm.Y
                                                                                  + zend * unitelenorm.Z + d)
                                                                                          /
                                                                                  sqrt(unitelenorm.X*unitelenorm.X
                                                                                                          + unitelenorm.Y*unitelenorm.Y
                                                                                                          + unitelenorm.Z*unitelenorm.Z);
                          distele = fabs(distele);	// if only magnitude is
          required
                          */

                  Vector3D eleOrigin;
                  eleOrigin.X = (EleArr + ele - 1)->G.Origin.X;
                  eleOrigin.Y = (EleArr + ele - 1)->G.Origin.Y;
                  eleOrigin.Z = (EleArr + ele - 1)->G.Origin.Z;
                  distele = (eleOrigin.X - xend) * (eleOrigin.X - xend) +
                            (eleOrigin.Y - yend) * (eleOrigin.Y - yend) +
                            (eleOrigin.Z - zend) * (eleOrigin.Z - zend);
                  distele = sqrt(distele);

                  if (ele == ElementBgn[nearestprim]) {
                    mindistele = distele;
                    nearestele = ele;
                  }
                  if (distele < mindistele) {
                    mindistele = distele;
                    nearestele = ele;
                  }

                  if (debugFn) {
                    // printf("a, b, c, d, dist: %lg, %lg, %lg, %lg, %lg\n",
                    // a, b, c, d,  dist);
                    // printf("vector n: ai + bj + ck\n");
                    // printf("vector v: xgrd, ygrd, zgrd: %lg, %lg, %lg\n",
                    // xgrd, ygrd, zgrd);
                    printf(
                        "distele: %lg, mindist: %lg,  from nearest ele: %d\n",
                        distele, mindistele, nearestele);
                    fflush(stdout);
                  }

                  // }	// if PrimType is 3 or 4
                }  // for elements in nearestprim

                // if(mindistele <= 10.0e-6)
                // {
                EleIntsctd = nearestele;
                InEle = 1;
                ptintsct.X = (EleArr + EleIntsctd - 1)->G.Origin.X;
                ptintsct.Y = (EleArr + EleIntsctd - 1)->G.Origin.Y;
                ptintsct.Z = (EleArr + EleIntsctd - 1)->G.Origin.Z;
                NbChUpIonEle[EleIntsctd]++;

                fprintf(fPtIChUpMap, "%d %lg %lg %lg %d %d %d %d\n", ion,
                        ptintsct.X, ptintsct.Y, ptintsct.Z, PrimIntsctd, InPrim,
                        EleIntsctd, InEle);
                // }

                if (debugFn) {
                  printf("# ion: %d\n", ion);
                  printf("%lg %lg %lg\n", xlbend, ylbend, zlbend);
                  printf("%lg %lg %lg\n", xend, yend, zend);
                  printf("%lg, %lg, %lg\n", ptintsct.X, ptintsct.Y, ptintsct.Z);
                  printf("# Associated primitive: %d\n", PrimIntsctd);
                  printf("# Associated element and origin: %d, %lg, %lg, %lg\n",
                         EleIntsctd, (EleArr + EleIntsctd - 1)->G.Origin.X,
                         (EleArr + EleIntsctd - 1)->G.Origin.Y,
                         (EleArr + EleIntsctd - 1)->G.Origin.Z);
                  printf("#NbChUpIonEle on element: %d\n",
                         NbChUpIonEle[EleIntsctd]);
                  fprintf(ftmpIF, "#Element: %d\n", EleIntsctd);
                  polynode[0].X = (EleArr + EleIntsctd - 1)->G.Vertex[0].X;
                  polynode[0].Y = (EleArr + EleIntsctd - 1)->G.Vertex[0].Y;
                  polynode[0].Z = (EleArr + EleIntsctd - 1)->G.Vertex[0].Z;
                  polynode[1].X = (EleArr + EleIntsctd - 1)->G.Vertex[1].X;
                  polynode[1].Y = (EleArr + EleIntsctd - 1)->G.Vertex[1].Y;
                  polynode[1].Z = (EleArr + EleIntsctd - 1)->G.Vertex[1].Z;
                  polynode[2].X = (EleArr + EleIntsctd - 1)->G.Vertex[2].X;
                  polynode[2].Y = (EleArr + EleIntsctd - 1)->G.Vertex[2].Y;
                  polynode[2].Z = (EleArr + EleIntsctd - 1)->G.Vertex[2].Z;
                  if (nvert == 4) {
                    polynode[3].X = (EleArr + EleIntsctd - 1)->G.Vertex[3].X;
                    polynode[3].Y = (EleArr + EleIntsctd - 1)->G.Vertex[3].Y;
                    polynode[3].Z = (EleArr + EleIntsctd - 1)->G.Vertex[3].Z;
                  }
                  fprintf(ftmpIF, "%lg %lg %lg\n", polynode[0].X, polynode[0].Y,
                          polynode[0].Z);
                  fprintf(ftmpIF, "%lg %lg %lg\n", polynode[1].X, polynode[1].Y,
                          polynode[1].Z);
                  fprintf(ftmpIF, "%lg %lg %lg\n", polynode[2].X, polynode[2].Y,
                          polynode[2].Z);
                  if (nvert == 4) {
                    fprintf(ftmpIF, "%lg %lg %lg\n", polynode[3].X,
                            polynode[3].Y, polynode[3].Z);
                  }
                  fprintf(ftmpIF, "%lg %lg %lg\n", polynode[0].X, polynode[0].Y,
                          polynode[0].Z);
                  fprintf(ftmpIF, "\n");
                  fflush(stdout);
                }  // debug
              }    // if prim == NbPrimitives

            }  // for all primitives // just not those on the volume

            if (debugFn)  // check ion positions, volume primitives and elements
            {
              char ionposdbg[256], inbstr[10];
              sprintf(inbstr, "%d", ion);
              strcpy(ionposdbg, "/tmp/Ion");
              strcat(ionposdbg, inbstr);
              strcat(ionposdbg, ".out");
              FILE *fipd = fopen(ionposdbg, "w");
              if (fipd == NULL) {
                printf(
                    "cannot open writable file to debug ion positions ...\n");
                printf("returning ...\n");
                return -111;
              }
              // write electron number, end, lbend, volume, primitive, elements,
              // intxn
              fprintf(fipd, "#ion: %d %d\n", inb, ion);  // should print same
              fprintf(fipd, "#last but end position:\n");
              fprintf(fipd, "%lg %lg %lg\n", xlbend, ylbend, zlbend);
              fprintf(fipd, "#end position:\n");
              fprintf(fipd, "%lg %lg %lg\n\n", xend, yend, zend);

              fprintf(fipd, "#intersected primitive number: %d\n", PrimIntsctd);
              if (PrimIntsctd >= 1) {
                fprintf(fipd, "#PrimType: %d\n", PrimType[PrimIntsctd]);
                fprintf(fipd, "#prim vertices:\n");
                fprintf(fipd, "%lg %lg %lg\n", XVertex[PrimIntsctd][0],
                        YVertex[PrimIntsctd][0], ZVertex[PrimIntsctd][0]);
                fprintf(fipd, "%lg %lg %lg\n", XVertex[PrimIntsctd][1],
                        YVertex[PrimIntsctd][1], ZVertex[PrimIntsctd][1]);
                fprintf(fipd, "%lg %lg %lg\n", XVertex[PrimIntsctd][2],
                        YVertex[PrimIntsctd][2], ZVertex[PrimIntsctd][2]);
                if (PrimType[PrimIntsctd] == 4) {
                  fprintf(fipd, "%lg %lg %lg\n", XVertex[PrimIntsctd][3],
                          YVertex[PrimIntsctd][3], ZVertex[PrimIntsctd][3]);
                }
                fprintf(fipd, "%lg %lg %lg\n", XVertex[PrimIntsctd][0],
                        YVertex[PrimIntsctd][0], ZVertex[PrimIntsctd][0]);
                fprintf(fipd, "\n");

                fprintf(fipd, "#ptintsct:\n");
                fprintf(fipd, "%lg %lg %lg\n", ptintsct.X, ptintsct.Y,
                        ptintsct.Z);
                fprintf(fipd, "\n");
              }

              fprintf(fipd, "#intersected element number: %d\n", EleIntsctd);
              if (EleIntsctd >= 1) {
                int gtype = (EleArr + EleIntsctd - 1)->G.Type;
                fprintf(fipd, "#EleType: %d\n", gtype);
                fprintf(fipd, "#element vertices:\n");
                double x0 = (EleArr + EleIntsctd - 1)->G.Vertex[0].X;
                double y0 = (EleArr + EleIntsctd - 1)->G.Vertex[0].Y;
                double z0 = (EleArr + EleIntsctd - 1)->G.Vertex[0].Z;
                double x1 = (EleArr + EleIntsctd - 1)->G.Vertex[1].X;
                double y1 = (EleArr + EleIntsctd - 1)->G.Vertex[1].Y;
                double z1 = (EleArr + EleIntsctd - 1)->G.Vertex[1].Z;
                double x2 = (EleArr + EleIntsctd - 1)->G.Vertex[2].X;
                double y2 = (EleArr + EleIntsctd - 1)->G.Vertex[2].Y;
                double z2 = (EleArr + EleIntsctd - 1)->G.Vertex[2].Z;
                fprintf(fipd, "%lg %lg %lg\n", x0, y0, z0);
                fprintf(fipd, "%lg %lg %lg\n", x1, y1, z1);
                fprintf(fipd, "%lg %lg %lg\n", x2, y2, z2);
                if (gtype == 4) {
                  double x3 = (EleArr + EleIntsctd - 1)->G.Vertex[3].X;
                  double y3 = (EleArr + EleIntsctd - 1)->G.Vertex[3].Y;
                  double z3 = (EleArr + EleIntsctd - 1)->G.Vertex[3].Z;
                  fprintf(fipd, "%lg %lg %lg\n", x3, y3, z3);
                }
                fprintf(fipd, "%lg %lg %lg\n", x0, y0, z0);
                fprintf(fipd, "\n");

                fprintf(fipd, "#ptintsct:\n");
                fprintf(fipd, "%lg %lg %lg\n", ptintsct.X, ptintsct.Y,
                        ptintsct.Z);
                fprintf(fipd, "\n");
              }
              fclose(fipd);
            }  // if 1
          }    // for all the ions
          fclose(fPtIChUpMap);

          // This file contains information about number of ions (I)
          // and total (E+I) charge deposition on each element
          FILE *fEleEIChUpMap = fopen("EleE+IChUpMap.out", "w");
          if (fEleEIChUpMap == NULL) {
            printf("cannot open EleE+IChUpMap.out file for writing ...\n");
            return 111;
          }
          for (int ele = 1; ele <= NbElements; ++ele) {
            (EleArr + ele - 1)->Assigned +=
                ChUpFactor * Q_I * NbChUpIonEle[ele] / (EleArr + ele - 1)->G.dA;
            fprintf(fEleEIChUpMap, "%d %lg %lg %lg %d %lg\n", ele,
                    (EleArr + ele - 1)->G.Origin.X,
                    (EleArr + ele - 1)->G.Origin.Y,
                    (EleArr + ele - 1)->G.Origin.Z, NbChUpIonEle[ele],
                    (EleArr + ele - 1)->Assigned);
          }
          fclose(fEleEIChUpMap);

          fclose(ftmpIF);
          free(NbChUpIonEle);
        }  // Calculation for ions ends

      }  // OptChargingUp

      fclose(ChargingUpInpFile);
    }  // else ChargingUpInpFile

    if (debugFn) {
      // print all the elements and their number of charging up e-s and i-s
    }
  }  // charging up parameters set up

  return (0);
}  // InitChargingUp ends

#ifdef __cplusplus
}  // namespace
#endif
