#ifndef _neBEMInterface_H_
#define _neBEMInterface_H_

#ifdef DEFINE_INTFACEGLOBAL
#define INTFACEGLOBAL
#else
#define INTFACEGLOBAL extern
#endif

#include <stdio.h>
#include <time.h>

#include "Vector.h"

#ifdef __cplusplus
namespace neBEM {
#endif

INTFACEGLOBAL clock_t startClock, stopClock;  // time related variables

INTFACEGLOBAL double neBEMTimeElapsed(clock_t, clock_t);  // elapsed time

INTFACEGLOBAL int neBEMState;  // the state in which neBEM is

// Print relevant messages to stdout device
INTFACEGLOBAL int neBEMMessage(const char *message);

// return number of lines in a file
INTFACEGLOBAL int neBEMGetNbOfLines(const char *fname);

// check whether a 3D point is within a 3D plane
// #define neBEMepsilon 1.0e-12
#define neBEMtwopi 6.283185307179586476925287
#define neBEMrtod 57.2957795
#define MODULUS(p) (sqrt(p.X * p.X + p.Y * p.Y + p.Z * p.Z))
INTFACEGLOBAL double neBEMChkInPoly(int nvert, Point3D node[], Point3D pt);

// load default values of important global variables
INTFACEGLOBAL int neBEMSetDefaults(void);
INTFACEGLOBAL int ReadInitFile(char filename[]);

// Geometry details from files
// Maintains compatilibity with the older approach
INTFACEGLOBAL int OptDeviceFile;
INTFACEGLOBAL char DeviceInputFile[256];
INTFACEGLOBAL int neBEMGetInputsFromFiles(void);

// Number of threads requested
INTFACEGLOBAL int NbThreads;

INTFACEGLOBAL int neBEMInitialize(void);
INTFACEGLOBAL int OptPrintPrimaryDetails, OptPrintVolumeDetails;

INTFACEGLOBAL int OptGnuplot, OptGnuplotPrimitives, OptGnuplotElements;
INTFACEGLOBAL FILE *fgnuPrim, *fgnuElem, *fgnuMesh;  // gnu plot files

// Geometry from Garfield or some other Pre-processor
INTFACEGLOBAL int neBEMReadGeometry(void);

INTFACEGLOBAL int neBEMGetNbPrimitives(void);

INTFACEGLOBAL int OptPrintVertexAndNormal, OptPrimitiveFiles;
INTFACEGLOBAL int neBEMGetPrimitive(int prim, int *nvertex, double xvert[],
                                    double yvert[], double zvert[],
                                    double *xnorm, double *ynorm, double *znorm,
                                    int *volref1, int *volref2);

INTFACEGLOBAL int neBEMVolumeDescription(int ivol, int *shape, int *material,
                                         double *epsilon, double *potential,
                                         double *charge, int *boundarytype);

INTFACEGLOBAL int neBEMVolumePoint(double x, double y, double z);

INTFACEGLOBAL void neBEMVolumePrimitives(int volume, int *nprim,
                                         int primlist[]);

// Periodicity information from Garfield
// int prim introduced since V1.03 to allow different periodicities for
// different primitives
INTFACEGLOBAL int neBEMGetPeriodicities(int prim, int *ix, int *jx, double *sx,
                                        int *iy, int *jy, double *sy, int *iz,
                                        int *jz, double *sz);

// Mirror information from Garfield
// different mirror properties for different primtives are allowed
INTFACEGLOBAL int neBEMGetMirror(int prim, int *ix, int *jx, double *sx,
                                 int *iy, int *jy, double *sy, int *iz, int *jz,
                                 double *sz);

// Bounding plane information from Garfield
INTFACEGLOBAL int neBEMGetBoundingPlanes(
    int *ixmin, double *cxmin, double *vxmin, int *ixmax, double *cxmax,
    double *vxmax, int *iymin, double *cymin, double *vymin, int *iymax,
    double *cymax, double *vymax, int *izmin, double *czmin, double *vzmin,
    int *izmax, double *czmax, double *vzmax);

// argument: array containing segment nbs in each dirn
INTFACEGLOBAL int OptElementFiles;
INTFACEGLOBAL int neBEMDiscretize(int **elementNbs);

// Set boundary condition
INTFACEGLOBAL int neBEMBoundaryInitialConditions(void);

// argument: status of matrix inversion (0: not inverted; 1: inverted)
INTFACEGLOBAL int neBEMSolve(void);

// arguments: a 3D point, potential and flux vector, the last two to be
// returned by the function
INTFACEGLOBAL int neBEMPFCallCntr;
INTFACEGLOBAL int neBEMPF(Point3D *, double *, Vector3D *);

// returns the total charge in a volume
// the volume identifier is the argument
INTFACEGLOBAL double neBEMVolumeCharge(int volume);

// Weighting field calculation preparation
// arguments: number of primitives considered for this weighting field
// computation and the related list; returns the identification tag for the
// solution for this weighting field set-up
INTFACEGLOBAL int neBEMPrepareWeightingField(int NbPrimsWtField,
                                             int PrimListWtField[]);

// Deallocates memory reserved for a weighting field
INTFACEGLOBAL void neBEMDeleteWeightingField(int IdWtField);

// Deallocates all memory reserved for all weighting fields
INTFACEGLOBAL void neBEMDeleteAllWeightingFields(void);

// Get weighting field and potential at a specific point
// arguments: evaluation position, field vector and the identification tag
// to indicate the necessary weighting field configuration; returns potential
// values evaluation failure returns DBL_MAX for the value of potential
INTFACEGLOBAL double neBEMWeightingField(Point3D *point, Vector3D *field,
                                         int IdWtField);

// no argument
INTFACEGLOBAL int neBEMEnd(void);

// no argument
// const int OptReuseDir = 1;
// #define OptReuseDir 1
INTFACEGLOBAL int OptReuseDir;
INTFACEGLOBAL int CreateDirStr(void);

// Read-Write Primitives and Elements
INTFACEGLOBAL int WritePrimitives(void);
INTFACEGLOBAL int ReadPrimitives(void);
INTFACEGLOBAL int WriteElements(void);
INTFACEGLOBAL int ReadElements(void);

#ifdef __cplusplus
}  // namespace
#endif

#endif
