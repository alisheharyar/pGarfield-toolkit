/*
(c) 2005, Supratik Mukhopadhayay, Nayana Majumdar
*/
#ifndef _neBEM_H_
#define _neBEM_H_

#ifdef DEFINE_neBEMGLOBAL
#define neBEMGLOBAL
#else
#define neBEMGLOBAL extern
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Vector.h"

#define EPS0 8.854187817e-12       // in C2/Nm2 (SI) equivalent to pF/m
#define MyFACTOR 111.26500547e-12  // 4 pi eps_0 is ofn reqrd for normalization
// #define EPS0 1.0						// relevant for other
// physics problems
// #define MyFACTOR 1.0				// relevant for other physics
// problems

#define Q_E -1.60217646e-19  // charge of electron in SI units (Coulomb)
#define Q_I 1.60217646e-19   // charge of ion in SI units (Coulomb)

#define MAXWtFld 11

#ifdef __cplusplus
namespace neBEM {
#endif

neBEMGLOBAL char neBEMVersion[10];

// options to decide whether or not storage of influence and inverted matrices
// are desired
// The following, especially OptValidateSolution should, prefereably, be 1
// Storing geometry and other info consumes disk space
// RepeatLHMatrix consumes time
neBEMGLOBAL int OptInvMatProc;
neBEMGLOBAL int OptValidateSolution;
neBEMGLOBAL int OptForceValidation;
neBEMGLOBAL int OptEstimateError;
neBEMGLOBAL int OptStoreInflMatrix;
neBEMGLOBAL int OptReadInflMatrix;
neBEMGLOBAL int OptStoreInvMatrix;
neBEMGLOBAL int OptReadInvMatrix;
neBEMGLOBAL int OptStorePrimitives;
neBEMGLOBAL int OptReadPrimitives;
neBEMGLOBAL int OptRmPrim;
neBEMGLOBAL int OptStoreElements;
neBEMGLOBAL int OptReadElements;
neBEMGLOBAL int OptFormattedFile;
neBEMGLOBAL int OptUnformattedFile;
neBEMGLOBAL int OptRepeatLHMatrix;
neBEMGLOBAL int OptKnCh;
neBEMGLOBAL int OptChargingUp;

// Geometry variables
neBEMGLOBAL int NbVolumes;
neBEMGLOBAL int NbPrimitives;
neBEMGLOBAL int OrgnlNbPrimitives;  // original nmbr may be less than effective
neBEMGLOBAL int MaxNbVertices;  // maximum allowed number (only 4, at present)
                                // allocation from 0 to MaxNbVertices - 1

// Related to volumes
neBEMGLOBAL int *volRef, *volShape, *volMaterial, *volBoundaryType;
neBEMGLOBAL double *volEpsilon, *volPotential, *volCharge;

// Related to primitives
neBEMGLOBAL int *PrimType,
    *InterfaceType;  // removed redundant *PrimCnt in V1.7.2
neBEMGLOBAL int **OrgnlToEffPrim;
neBEMGLOBAL int *NbVertices;
neBEMGLOBAL double **XVertex, **YVertex, **ZVertex, *XNorm, *YNorm, *ZNorm,
    *PrimLX, *PrimLZ, *Radius;
neBEMGLOBAL double *PrimOriginX, *PrimOriginY, *PrimOriginZ;
neBEMGLOBAL DirnCosn3D *PrimDC;
neBEMGLOBAL double *Epsilon1, *Epsilon2, *Lambda, *ApplPot, *ApplCh;
neBEMGLOBAL int *VolRef1, *VolRef2;
neBEMGLOBAL int VolMax;
neBEMGLOBAL int *PeriodicTypeX, *PeriodicTypeY, *PeriodicTypeZ;
neBEMGLOBAL int *PeriodicInX, *PeriodicInY, *PeriodicInZ;
neBEMGLOBAL double *XPeriod, *YPeriod, *ZPeriod;
neBEMGLOBAL int *MirrorTypeX, *MirrorTypeY, *MirrorTypeZ;
neBEMGLOBAL double *MirrorDistXFromOrigin, *MirrorDistYFromOrigin,
    *MirrorDistZFromOrigin;
neBEMGLOBAL int *BndPlaneInXMin, *BndPlaneInYMin, *BndPlaneInZMin;
neBEMGLOBAL int *BndPlaneInXMax, *BndPlaneInYMax, *BndPlaneInZMax;
neBEMGLOBAL double *XBndPlaneInXMin, *YBndPlaneInYMin, *ZBndPlaneInZMin;
neBEMGLOBAL double *XBndPlaneInXMax, *YBndPlaneInYMax, *ZBndPlaneInZMax;
neBEMGLOBAL double *VBndPlaneInXMin, *VBndPlaneInYMin, *VBndPlaneInZMin;
neBEMGLOBAL double *VBndPlaneInXMax, *VBndPlaneInYMax, *VBndPlaneInZMax;
neBEMGLOBAL double *AvChDen, *AvAsgndChDen;

// Related to both primitives and elements
// Beginning and ending element numbers on a griven primitive
neBEMGLOBAL int *NbElmntsOnPrim, *ElementBgn, *ElementEnd;

// Related to surfaces and wires
neBEMGLOBAL int NbSurfs, NbWires;
neBEMGLOBAL int *NbSurfSegX, *NbSurfSegZ;
neBEMGLOBAL int *NbWireSeg;

// Related to elements
// minimum number of elements allowed along Length
neBEMGLOBAL int MinNbElementsOnLength;
// maximum number of elements allowed along Length
neBEMGLOBAL int MaxNbElementsOnLength;
// user requested length of along Length
neBEMGLOBAL double ElementLengthRqstd;

// int MinNbElementsOnSurface;	// minimum number of elements allowed on a
// surface int MaxNbElementsOnSurface;	// maximum number of elements allowed on
// a surface double ElementAreaRqstd;	// user requested area of each element
neBEMGLOBAL int EleCntr;     // Element counter
neBEMGLOBAL int NbElements;  // total number of elements
neBEMGLOBAL FILE *fMeshLog;

// Related to solution constraints
// Whether total charge in the system is zero
neBEMGLOBAL int OptSystemChargeZero;
// which eqn / unknown is related to this constraint
neBEMGLOBAL int NbSystemChargeZero;
// voltage shift necessary to enforce zero charge
neBEMGLOBAL double VSystemChargeZero;

neBEMGLOBAL int NbFloatingConductors;  // Number of floating conductors
// which eqn / unknown is related to this constraint
neBEMGLOBAL int NbFloatCon;
// value of floating potential on the conductor
neBEMGLOBAL double VFloatCon;

typedef struct {
  short int Type;     // 4: rectangular, 3: triangular, 2: linear (wire)
  Point3D Origin;     // centroid / barycenter / axis-center (local origin)
  Point3D Vertex[4];  // element vertex begins with index 0 and goes to max 3
  double LX, LZ;      // length, breadth / base, height / radius, length
  double dA;          // area
  DirnCosn3D DC;      // Direction cosines
} GeomProp;

// 1: conductor at known potential, 2: charged conductor, 3: floating conductor
// 4: DD interface satisfying continuity, 5: charged DD interface
// 6: E parallel symmetry, 7: E perpendicular symmetry.
typedef struct {
  short int Type;
  double Lambda;  // ratio of dielectric permiitivites
} ElecProp;

typedef struct {
  short int NbOfBCs;  // nb of boundary conditions on this element
  Point3D CollPt;     // Collocation (only one, for the time being)
  double Value;       // potential / charge density
} BCProp;

// we need a reference to the volumes (volref1 and volref2) that an element
// belongs to
typedef struct {
  short int DeviceNb;  // each setup can be made of several devices
  int ComponentNb;     // each device made of several components
  int PrimitiveNb;     // each component can be made of several primitives
  int InterfaceId;
  int Id;  // element id number - each made of several elements
  // Point3D Vertex[4];	// since we consider only upto rectangles: within G
  GeomProp G;  // geomtype, origin, vertex, lengths, area, direction cosines
  ElecProp E;  // electype, BC value
  BCProp BC;   // boundary condn properties (should this BC thing be freed?)
  double Solution;  // accumulated charge, or similar solution
  double Assigned;  // assigned charge, or similar property
} Element;

neBEMGLOBAL Element *EleArr;  // represents an array of elements

// Various flags that need to be passed on to the solver to decide from what
// level the solution should proceed
neBEMGLOBAL int NewModel, NewMesh, NewBC, NewPP;
// Various counters to keep track of the various possible solutions for a
// given device.
// The first of these variables is being maintained for possible future use.
neBEMGLOBAL int ModelCntr, MeshCntr, BCCntr, PPCntr;

// Computation parameter
// Variables used throughout the solution procedure
// Unless we have an over- or under-determined system, the number of equations
// and unknowns will be the same.
// The number of equations will be as follows:
// One for each element that carries an unknown charge density
// One for each unknown related to each constraint equation
neBEMGLOBAL int NbConstraints;  // Arising of different physical considerations
neBEMGLOBAL int NbEqns, NbUnknowns,  // rows and columns in [Inf]
    DebugLevel;
neBEMGLOBAL int OptSVD, OptLU, OptGSL;  // option SVD, LU and GSL decompositions
neBEMGLOBAL double **Inf, **InvMat, *RHS, *Solution;
neBEMGLOBAL double LengthScale;

// Variables to facilitate time stepping
neBEMGLOBAL int TimeStep, EndOfTime;
neBEMGLOBAL char TimeStr[256];

// Variables related to geometry viewing
neBEMGLOBAL char GnuplotTmpDir[256], GnuplotScriptFile[256];  // Dirs and files

// Outputs are written in various subdirectories.
// DeviceOutDir:
// The top level directory for a device under study.
//  ModelOutDir:
// Files related to a unique model used to define a given device are kept here.
// MeshOutDir :
// Files down to elements used to discretize primitves device are kept here.
// This directory is also the place to store the inverted matrix file.
// BCOutDir :
// Files down to RHS, solution vectors, XChk are kept here.
// PPDir:
// Rest of the files, i.e., those related to post-processing are kept here.
// In addition, we have the GViewDir directory (as a subdirectory of each
// MeshOutDir) wherein the gnuplot files needed to view the element geometry,
// primitives, elements and mesh are stored.
neBEMGLOBAL char DeviceOutDir[256], ModelOutDir[256], NativeOutDir[256],
    NativePrimDir[256], MeshOutDir[256], BCOutDir[256], PPOutDir[256];

// Source code: PreProcess/ReTrim.c
neBEMGLOBAL int SurfaceElements(int prim, int nvertex, double xvert[],
                                double yvert[], double zvert[], double xnorm,
                                double ynorm, double znorm, int volref1,
                                int volref2, int inttype, double potential,
                                double charge, double lambda, int NbSegX,
                                int NbSegZ);
neBEMGLOBAL int WireElements(int prim, int nvertex, double xvert[],
                             double yvert[], double zvert[], double radius,
                             int volref1, int volref2, int inttype,
                             double potential, double charge, double lambda,
                             int NbWireSeg);
neBEMGLOBAL int BoundaryConditions(void);
neBEMGLOBAL int InitialConditions(void);
// Initiate known charge(s) / charge density (ies) within the device.
neBEMGLOBAL int InitKnownCharges(void);
neBEMGLOBAL int InitChargingUp(void);

// Number of points, lines, areas, volumes with known property (density) values
neBEMGLOBAL int NbPointsKnCh, NbLinesKnCh, NbAreasKnCh, NbVolumesKnCh;
// Details of each of the above entities with known property
// Point with known charge
typedef struct {
  int Nb;
  Point3D P;  // available from Vector.h
  double Assigned;
} PointKnCh;

neBEMGLOBAL PointKnCh *PointKnChArr;

// Line with known linear charge density
typedef struct {
  int Nb;
  Point3D Start;  // the line extends from Start to Stop
  Point3D Stop;
  double Radius;
  double Assigned;
} LineKnCh;

neBEMGLOBAL LineKnCh *LineKnChArr;

// Surface with known charge density
// Surfaces can be right-triangular or rectangular
// Needs to be extended to arbitrary polygons
typedef struct {
  int Nb;
  int NbVertices;  // number of vertices defining the surface (restricted to 4)
  Point3D Vertex[5];  // the surface is defined by Vertex[1-4]
  double Assigned;
} AreaKnCh;

neBEMGLOBAL AreaKnCh *AreaKnChArr;

// Volumes with known charge density
// Volumes can be tetrahedral or rectangular
// Needs to be extended to arbitrarly-shaped volumes
typedef struct {
  int Nb;
  int NbVertices;  // number of vertices defining the volume (restricted to 8)
  Point3D Vertex[9];
  double Assigned;
} VolumeKnCh;

neBEMGLOBAL VolumeKnCh *VolumeKnChArr;

neBEMGLOBAL int AnalyzePrimitive(int, int *, int *);
neBEMGLOBAL int AnalyzeWire(int, int *);
neBEMGLOBAL int AnalyzeSurface(int, int *, int *);
neBEMGLOBAL int DiscretizeWire(int prim, int nvertex, double xvert[],
                               double yvert[], double zvert[], double radius,
                               int volref1, int volref2, int inttype,
                               double potential, double charge, double lambda,
                               int NbSegs);
neBEMGLOBAL int DiscretizeTriangle(int prim, int nvertex, double xvert[],
                                   double yvert[], double zvert[], double xnorm,
                                   double ynorm, double znorm, int volref1,
                                   int volref2, int inttype, double potential,
                                   double charge, double lambda, int NbSegX,
                                   int NbSegZ);
neBEMGLOBAL int DiscretizeRectangle(int prim, int nvertex, double xvert[],
                                    double yvert[], double zvert[],
                                    double xnorm, double ynorm, double znorm,
                                    int volref1, int volref2, int inttype,
                                    double potential, double charge,
                                    double lambda, int NbSegX, int NbSegZ);
neBEMGLOBAL int DiscretizePolygon(int prim, int nvertex, double xvert[],
                                  double yvert[], double zvert[], double xnorm,
                                  double ynorm, double znorm, int volref1,
                                  int volref2, int inttype, double potential,
                                  double charge, double lambda, int NbSegX,
                                  int NbSegZ);

// Source code: Solve/neBEM.c
neBEMGLOBAL int ComputeSolution(void);
neBEMGLOBAL int ReadSolution(void);
neBEMGLOBAL int UpdateKnownCharges(void);
neBEMGLOBAL int UpdateChargingUp(void);

// Apparently, localPt need not be passed to ComputeInfluence. This is true if
// there are no repetitions / reflections. With these and similar possibilities,
// localPt may not strictly be computed using only the details of the field
// and source elements. It would further require information on repetition and
// reflection which can lead to a quick increase in the amount of arguments
// being passed to the function and, as a result, detrimental to the efficiency.
neBEMGLOBAL double ComputeInfluence(int fld, int src, Point3D *localPt,
                                    DirnCosn3D *DirCos);
neBEMGLOBAL double SatisfyValue(int src, Point3D *localPt);
neBEMGLOBAL double SatisfyContinuity(int fld, int src, Point3D *localPt,
                                     DirnCosn3D *DirCos);
// neBEMGLOBAL double EffectChUp(int caseid, int fld);
neBEMGLOBAL double EffectChUp(int fld);
neBEMGLOBAL double ValueChUp(int fld);
neBEMGLOBAL double ContinuityChUp(int fld);
// neBEMGLOBAL double EffectKnCh(int caseid, int fld);
neBEMGLOBAL double EffectKnCh(int fld);
neBEMGLOBAL double ValueKnCh(int fld);
neBEMGLOBAL double ContinuityKnCh(int fld);

// Weighting field charge density solution
// arguments: boundary condition array, and the solution (charge density, for
// electrostatic problems) array; returns success (0) or failure (non-zero)
neBEMGLOBAL double **WtFieldChDen, **AvWtChDen;
neBEMGLOBAL int WeightingFieldSolution(int NbPrimsWtField,
                                       int PrimListWtField[],
                                       double WtFieldChDen[]);
neBEMGLOBAL Point3D ReflectPrimitiveOnMirror(char Axis, int prim, Point3D srcpt,
                                             Point3D fieldpt, double distance,
                                             DirnCosn3D *DirCos);
neBEMGLOBAL Point3D ReflectOnMirror(char Axis, int elesrc, Point3D srcpt,
                                    Point3D fieldpt, double distance,
                                    DirnCosn3D *DirCos);

// Source code: Solve/ComputeProperties.c
// Compute potential and flux components at globalPt due to all elements
neBEMGLOBAL int PFAtPoint(Point3D *globalPt, double *Pot, Vector3D *Flux);
neBEMGLOBAL int ElePFAtPoint(Point3D *globalPt, double *Pot, Vector3D *Flux);
neBEMGLOBAL int KnChPFAtPoint(Point3D *globalPt, double *Pot, Vector3D *Flux);

// Choose between element and primitive representations
neBEMGLOBAL int PrimAfter;       // for physical field computations.
neBEMGLOBAL int WtFldPrimAfter;  // for weighting field computations.

// Variables and functions related to voxelized output for Garfield++
neBEMGLOBAL int OptVoxel;
neBEMGLOBAL int OptStaggerVoxel;
// inputs for the voxelized volume
typedef struct {
  double Xmin;
  double Xmax;
  double Ymin;
  double Ymax;
  double Zmin;
  double Zmax;
  double XStagger;
  double YStagger;
  double ZStagger;
  int NbXCells;
  int NbYCells;
  int NbZCells;
} VoxelVol;
neBEMGLOBAL VoxelVol Voxel;

// Compute flux components, potential and regions within the voxelized
// volume and store them in an external file
neBEMGLOBAL int VoxelFPR(void);

// Variables and functions related to output 3dMap for Garfield++
neBEMGLOBAL int OptMap;
neBEMGLOBAL int OptStaggerMap;
neBEMGLOBAL char MapVersion[10];
// inputs for the voxelized volume
typedef struct {
  double Xmin;
  double Xmax;
  double Ymin;
  double Ymax;
  double Zmin;
  double Zmax;
  double XStagger;
  double YStagger;
  double ZStagger;
  int NbXCells;
  int NbYCells;
  int NbZCells;
} MapVol;
neBEMGLOBAL MapVol Map;

// Compute flux components, potential and regions within the 3dMap volume
// and store them in an external file
neBEMGLOBAL int MapFPR(void);

// Variables and functions related to the FAST algorithm (physical potential
// and fields)
neBEMGLOBAL int OptFastVol;
neBEMGLOBAL int OptStaggerFastVol;
neBEMGLOBAL int OptCreateFastPF;
neBEMGLOBAL int OptReadFastPF;
neBEMGLOBAL int VersionFV;   // although not being used, these two are
neBEMGLOBAL int NbBlocksFV;  // assigned values through ComponentNeBem3d
neBEMGLOBAL int NbPtSkip;
neBEMGLOBAL int NbStgPtSkip;
// inputs for the fast volume
typedef struct {
  double LX;
  double LY;
  double LZ;
  double CrnrX;
  double CrnrY;
  double CrnrZ;
  double YStagger;
  int NbBlocks;
  int NbOmitVols;
  int NbIgnoreVols;
} FastAlgoVol;
neBEMGLOBAL FastAlgoVol FastVol;
// inputs for each block within the fast volume
neBEMGLOBAL int *BlkNbXCells;
neBEMGLOBAL int *BlkNbYCells;
neBEMGLOBAL int *BlkNbZCells;
neBEMGLOBAL double *BlkLZ;
neBEMGLOBAL double *BlkCrnrZ;
neBEMGLOBAL double *OmitVolLX;
neBEMGLOBAL double *OmitVolLY;
neBEMGLOBAL double *OmitVolLZ;
neBEMGLOBAL double *OmitVolCrnrX;
neBEMGLOBAL double *OmitVolCrnrY;
neBEMGLOBAL double *OmitVolCrnrZ;
neBEMGLOBAL double *IgnoreVolLX;
neBEMGLOBAL double *IgnoreVolLY;
neBEMGLOBAL double *IgnoreVolLZ;
neBEMGLOBAL double *IgnoreVolCrnrX;
neBEMGLOBAL double *IgnoreVolCrnrY;
neBEMGLOBAL double *IgnoreVolCrnrZ;
// The following could have been members of a structure
// (see above for examples of creating arrays of structures)
neBEMGLOBAL double ****FastPot;
neBEMGLOBAL double ****FastFX, ****FastFY, ****FastFZ;
neBEMGLOBAL double ****StgFastPot;
neBEMGLOBAL double ****StgFastFX, ****StgFastFY, ****StgFastFZ;
neBEMGLOBAL double ****FastPotKnCh;
neBEMGLOBAL double ****FastFXKnCh, ****FastFYKnCh, ****FastFZKnCh;
neBEMGLOBAL double ****StgFastPotKnCh;
neBEMGLOBAL double ****StgFastFXKnCh, ****StgFastFYKnCh, ****StgFastFZKnCh;

// Create Fast volumes with potential and flux components
neBEMGLOBAL int CreateFastVolPF(void);
neBEMGLOBAL int CreateFastVolElePF(void);
neBEMGLOBAL int CreateFastVolKnChPF(void);

// Evaluate potential and flux components at globalPt using FAST algorithm
neBEMGLOBAL int FastPFAtPoint(Point3D *globalPt, double *Pot, Vector3D *Flux);
neBEMGLOBAL int FastElePFAtPoint(Point3D *globalPt, double *Pot,
                                 Vector3D *Flux);
neBEMGLOBAL int FastKnChPFAtPoint(Point3D *globalPt, double *Pot,
                                  Vector3D *Flux);

// identify the pickup electrode
neBEMGLOBAL int IdPkupElektrd;

// variables related to FIXED weighting field
neBEMGLOBAL int OptFixedWtField[MAXWtFld];
neBEMGLOBAL double FixedWtPotential[MAXWtFld];
neBEMGLOBAL double FixedWtFieldX[MAXWtFld];
neBEMGLOBAL double FixedWtFieldY[MAXWtFld];
neBEMGLOBAL double FixedWtFieldZ[MAXWtFld];

// Variable related to the weighting field FAST algorithm
neBEMGLOBAL int OptWtFldFastVol[MAXWtFld];
neBEMGLOBAL int OptStaggerWtFldFastVol[MAXWtFld];
neBEMGLOBAL int OptCreateWtFldFastPF[MAXWtFld];
neBEMGLOBAL int OptReadWtFldFastPF[MAXWtFld];
neBEMGLOBAL int VersionWtFldFV[MAXWtFld];   // although not being used, these
neBEMGLOBAL int NbBlocksWtFldFV[MAXWtFld];  // are assigned values using Set fns
neBEMGLOBAL int WtFldNbPtSkip[MAXWtFld];
neBEMGLOBAL int StgWtFldNbPtSkip[MAXWtFld];
// inputs for the fast volume
typedef struct {
  double LX;
  double LY;
  double LZ;
  double CrnrX;
  double CrnrY;
  double CrnrZ;
  double YStagger;
  int NbBlocks;
  int NbOmitVols;
  int NbIgnoreVols;
} WtFldFastAlgoVol;
neBEMGLOBAL WtFldFastAlgoVol WtFldFastVol[MAXWtFld];
// inputs for each block within the fast volume
neBEMGLOBAL int *WtFldBlkNbXCells[MAXWtFld];
neBEMGLOBAL int *WtFldBlkNbYCells[MAXWtFld];
neBEMGLOBAL int *WtFldBlkNbZCells[MAXWtFld];
neBEMGLOBAL double *WtFldBlkLZ[MAXWtFld];
neBEMGLOBAL double *WtFldBlkCrnrZ[MAXWtFld];
neBEMGLOBAL double *WtFldOmitVolLX[MAXWtFld];
neBEMGLOBAL double *WtFldOmitVolLY[MAXWtFld];
neBEMGLOBAL double *WtFldOmitVolLZ[MAXWtFld];
neBEMGLOBAL double *WtFldOmitVolCrnrX[MAXWtFld];
neBEMGLOBAL double *WtFldOmitVolCrnrY[MAXWtFld];
neBEMGLOBAL double *WtFldOmitVolCrnrZ[MAXWtFld];
neBEMGLOBAL double *WtFldIgnoreVolLX[MAXWtFld];
neBEMGLOBAL double *WtFldIgnoreVolLY[MAXWtFld];
neBEMGLOBAL double *WtFldIgnoreVolLZ[MAXWtFld];
neBEMGLOBAL double *WtFldIgnoreVolCrnrX[MAXWtFld];
neBEMGLOBAL double *WtFldIgnoreVolCrnrY[MAXWtFld];
neBEMGLOBAL double *WtFldIgnoreVolCrnrZ[MAXWtFld];
// The following could have been members of a structure
// (see above for examples of creating arrays of structures)
neBEMGLOBAL double ****WtFldFastPot[MAXWtFld];
neBEMGLOBAL double ****WtFldFastFX[MAXWtFld], ****WtFldFastFY[MAXWtFld],
    ****WtFldFastFZ[MAXWtFld];
neBEMGLOBAL double ****StgWtFldFastPot[MAXWtFld];
neBEMGLOBAL double ****StgWtFldFastFX[MAXWtFld], ****StgWtFldFastFY[MAXWtFld],
    ****StgWtFldFastFZ[MAXWtFld];

// Create Fast volumes with weighting potential and flux components
neBEMGLOBAL int CreateWtFldFastVolPF(int Id);

// Evaluate weighting potential and flux components at globalPt using FAST
// algorithm
neBEMGLOBAL int WtFldFastPFAtPoint(Point3D *globalPt, double *Pot,
                                   Vector3D *Flux, int Id);

// Weighting field values
neBEMGLOBAL int WtFldPFAtPoint(Point3D *globalPt, double *Pot, Vector3D *Flux,
                               int Id);

// Compute potential at xlocal, ylocal, zlocal due to an element defined by
// gtsrc (type), lxsrc (dimension), lzsrc (dimension), (length), dA (area)
// xlocal, ylocal, zlocal are measured in the element local coordinate system
// and the charge on the element is assumed to be unity
neBEMGLOBAL double GetPotential(int src, Point3D *localPt);

// Flux components at xlocal, ylocal, zlocal due to an element defined by
// gtsrc (type), lxsrc (dimension), lzsrc (dimension), (length), dA (area)
// xlocal, ylocal, zlocal and flux components are in the element local
// coordinate system and the charge on the element is assumed to be unity
neBEMGLOBAL void GetFluxGCS(int src, Point3D *localPt, Vector3D *Flux);
neBEMGLOBAL void GetFlux(int src, Point3D *localPt, Vector3D *Flux);
neBEMGLOBAL double RecPot(int src, Point3D *localPt);
neBEMGLOBAL double TriPot(int src, Point3D *localPt);
neBEMGLOBAL double WirePot(int src, Point3D *localPt);
neBEMGLOBAL void RecFlux(int src, Point3D *localPt, Vector3D *Flux);
neBEMGLOBAL void TriFlux(int src, Point3D *localPt, Vector3D *Flux);
neBEMGLOBAL void WireFlux(int src, Point3D *localPt, Vector3D *Flux);

// get both potential and flux
neBEMGLOBAL void GetPFGCS(int type, double a, double b, Point3D *localPt,
                          double *Pot, Vector3D *Flux, DirnCosn3D *DirCos);
neBEMGLOBAL void GetPF(int type, double a, double b, double x, double y,
                       double z, double *Pot, Vector3D *Flux);
neBEMGLOBAL void RecPF(double a, double b, double x, double y, double z,
                       double *Pot, Vector3D *Flux);
neBEMGLOBAL void TriPF(double a, double b, double x, double y, double z,
                       double *Pot, Vector3D *Flux);
neBEMGLOBAL void WirePF(double rW, double lW, double x, double y, double z,
                        double *Pot, Vector3D *Flux);
neBEMGLOBAL void GetPrimPFGCS(int src, Point3D *localPt, double *Pot,
                              Vector3D *Flux, DirnCosn3D *DirCos);
neBEMGLOBAL void GetPrimPF(int src, Point3D *localPt, double *Pot,
                           Vector3D *Flux);
neBEMGLOBAL void RecPrimPF(int src, Point3D *localPt, double *Pot,
                           Vector3D *Flux);
neBEMGLOBAL void TriPrimPF(int src, Point3D *localPt, double *Pot,
                           Vector3D *Flux);
neBEMGLOBAL void WirePrimPF(int src, Point3D *localPt, double *Pot,
                            Vector3D *Flux);

#ifdef __cplusplus
}  // namespace
#endif

#endif
