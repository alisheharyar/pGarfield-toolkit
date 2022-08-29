/*
(c) 2009 Supratik Mukhopadhyay, Nayana Majumdar
*/
// Gets number of primitives (polygonal surfaces (2D) and wires (1D)) and
// their details
// Assumptions:
// MaxNbVertices: 4 (maximum nb of vertices a polygon can have)
// Vertices have been counted starting from zero - it may be better to count
// them starting with 1
// Elements have a count from 0 to NbElement - 1: can create confusion

#define DEFINE_INTFACEGLOBAL
#define DEFINE_neBEMGLOBAL
#define DEFINE_NRGLOBAL

#include <float.h>
#include <stdio.h>
#include <sys/stat.h>  // use of stat function
#include <unistd.h>

#ifdef __cplusplus
#include <vector>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Isles.h"
#include "NR.h"
#include "neBEM.h"
#include "neBEMInterface.h"

#ifdef __cplusplus
namespace neBEM {
#endif

// Called from a code requesting neBEM services
int neBEMInitialize(void) {
  // Version information
  strcpy(neBEMVersion, "1.9.09");
  strcpy(ISLESVersion, "1.4.9");
  printf("Using neBEM version %s and ISLES version %s\n", neBEMVersion,
         ISLESVersion);

  // The following is not an absolute necessity, but can ease the burden of the
  // user by setting default values of some of the important variable. Please
  // take a look at the src/Interface/ExampleDev2neBEM.c to find out what values
  // are being set, and what they mean.
  int fstatus = neBEMSetDefaults();
  if (fstatus != 0) {
    neBEMMessage("neBEMInitialize - neBEMSetDefaults");
    return -1;
  }

  // Change some of the global variables according to requirements.
  // Note that the solution flags and counters are set before the neBEMSolve
  // is invoked.
  LengthScale = 1.0;
  DebugLevel = 0;

  // The following function allows input files to be read and the geometry
  // set according to these files, as was customary in the previous versions
  if (OptDeviceFile) {
    printf("Reading geometry details from %s\n", DeviceInputFile);
    fstatus = neBEMGetInputsFromFiles();
    if (fstatus != 0) {
      neBEMMessage("neBEMInitialize - neBEMGetInputFromFiles");
      return -1;
    }
  }

  // creation of the mother output directory should be necessary only once
  if (neBEMState == 0) {
    fstatus = CreateDirStr();
    if (fstatus != 0) {
      neBEMMessage("neBEMInitialize - CreateDirStr");
      return -1;
    }
  }

  // Create Isles log file for keeping track of approximations (numerical
  // quadrature) in case evaluation of algebraic expressions fails.
  char IslesFile[256];
  strcpy(IslesFile, PPOutDir);
  strcat(IslesFile, "/Isles.log");
  fIsles = fopen(IslesFile, "w");
  if (fIsles == NULL) {
    neBEMMessage("neBEMInitialize - IslesFile");
    return -1;
  }

  // following integers keep track of the success / failure of the exact
  // expressions in estimating realistic physical properties.
  IslesCntr = ExactCntr = FailureCntr = ApproxCntr = 0;

  // Set up parameters related to neBEM computations
  // The file will be removed soon
  int RqstdThreads = 1;
  if (NbThreads > 0) RqstdThreads = NbThreads;
  /*
  FILE *processFile = fopen("neBEMInp/neBEMProcess.inp", "r");
  if (processFile == NULL) {
    printf("neBEMProcess.inp absent ... assuming defaults ...\n");
    PrimAfter = 0;
    if (NbThreads > 0) RqstdThreads = NbThreads;
  } else {
    fscanf(processFile, "PrimAfter: %d\n", &PrimAfter);
    fscanf(processFile, "RqstdThreads: %d\n", &RqstdThreads);
    fclose(processFile);
  }
  */
  // PrimAfter is assigned a value through "SetPrimAfter" member function of
  // ComponentNeBem3d class.
  // Default value of PrimAfter is a negative integer.

#ifdef _OPENMP
  int MaxProcessors = omp_get_num_procs();
  if (RqstdThreads > 1) {
    if (RqstdThreads < MaxProcessors) {
      // one processor left alone
      omp_set_num_threads(RqstdThreads);
    } else {
      printf("RqstdThreads: %d\n", RqstdThreads);
      RqstdThreads = MaxProcessors - 1;
      omp_set_num_threads(RqstdThreads);
      printf("Adjusted RqstdThreads: %d\n", RqstdThreads);
    }
  } else {
    // Work with one thread
    RqstdThreads = 1;  // cannot be zero or negative!
    omp_set_num_threads(RqstdThreads);
    printf("RqstdThreads: %d => No Multi-threading ...\n", RqstdThreads);
  }

  // OpenMP related information
  printf("PrimAfter: %d\n", PrimAfter);
  printf("RqstdThreads: %d, MaxProcessors: %d\n", RqstdThreads, MaxProcessors);
  printf("Maximum number of threads to be used for parallelization: %d\n",
         omp_get_max_threads());
  printf("Number of threads used for neBEMInitialize: %d\n",
         omp_get_num_threads());
#endif

  // Set up parameters related to voxelized data export for Garfield++
  // To be removed soon
  FILE *voxelInpFile = fopen("neBEMInp/neBEMVoxel.inp", "r");
  if (voxelInpFile == NULL) {
    printf("neBEMVoxel.inp absent ... assuming OptVoxel = 0 ...\n");
    OptVoxel = 0;
    OptStaggerVoxel = 0;
  } else {
    fscanf(voxelInpFile, "OptVoxel: %d\n", &OptVoxel);
    fscanf(voxelInpFile, "OptStaggerVoxel: %d\n", &OptStaggerVoxel);
    fscanf(voxelInpFile, "Xmin: %le\n", &Voxel.Xmin);
    fscanf(voxelInpFile, "Xmax: %le\n", &Voxel.Xmax);
    fscanf(voxelInpFile, "Ymin: %le\n", &Voxel.Ymin);
    fscanf(voxelInpFile, "Ymax: %le\n", &Voxel.Ymax);
    fscanf(voxelInpFile, "Zmin: %le\n", &Voxel.Zmin);
    fscanf(voxelInpFile, "Zmax: %le\n", &Voxel.Zmax);
    fscanf(voxelInpFile, "XStagger: %le\n", &Voxel.XStagger);
    fscanf(voxelInpFile, "YStagger: %le\n", &Voxel.YStagger);
    fscanf(voxelInpFile, "ZStagger: %le\n", &Voxel.ZStagger);
    fscanf(voxelInpFile, "NbOfXCells: %d\n", &Voxel.NbXCells);
    fscanf(voxelInpFile, "NbOfYCells: %d\n", &Voxel.NbYCells);
    fscanf(voxelInpFile, "NbOfZCells: %d\n", &Voxel.NbZCells);
    fclose(voxelInpFile);
  }  // inputs for Voxel

  // Set up parameters related to 3dMap data export for Garfield++
  // To be removed soon
  FILE *mapInpFile = fopen("neBEMInp/neBEMMap.inp", "r");
  if (mapInpFile == NULL) {
    printf("neBEMMap.inp absent ... assuming OptMap = 0 ...\n");
    OptMap = 0;
    OptStaggerMap = 0;
  } else {
    // While reading the input, OptMap and OptStaggerMap have to be read
    // first since that will decide whether there is a map and its version.
    fscanf(mapInpFile, "OptMap: %d\n", &OptMap);
    fscanf(mapInpFile, "OptStaggerMap: %d\n", &OptStaggerMap);
    fscanf(mapInpFile, "MapVersion: %9s\n", MapVersion);
    fscanf(mapInpFile, "Xmin: %le\n", &Map.Xmin);
    fscanf(mapInpFile, "Xmax: %le\n", &Map.Xmax);
    fscanf(mapInpFile, "Ymin: %le\n", &Map.Ymin);
    fscanf(mapInpFile, "Ymax: %le\n", &Map.Ymax);
    fscanf(mapInpFile, "Zmin: %le\n", &Map.Zmin);
    fscanf(mapInpFile, "Zmax: %le\n", &Map.Zmax);
    fscanf(mapInpFile, "XStagger: %le\n", &Map.XStagger);
    fscanf(mapInpFile, "YStagger: %le\n", &Map.YStagger);
    fscanf(mapInpFile, "ZStagger: %le\n", &Map.ZStagger);
    fscanf(mapInpFile, "NbOfXCells: %d\n", &Map.NbXCells);
    fscanf(mapInpFile, "NbOfYCells: %d\n", &Map.NbYCells);
    fscanf(mapInpFile, "NbOfZCells: %d\n", &Map.NbZCells);
    fclose(mapInpFile);
  }  // inputs for 3dMap

  // Set up parameters related to fast volume
  if (OptFastVol) {
    FILE *fastInpFile = fopen("neBEMInp/neBEMFastVol.inp", "r");
    if (fastInpFile == NULL) {
      printf("neBEMFastVol.inp absent ... assuming OptFastVol = 0 ...\n");
      OptFastVol = 0;
      OptStaggerFastVol = 0;
      OptCreateFastPF = 0;
      OptReadFastPF = 0;
      FastVol.NbBlocks = 0;
      FastVol.NbOmitVols = 0;
      FastVol.NbIgnoreVols = 0;
    } else {
      fscanf(fastInpFile, "OptFastVol: %d\n", &OptFastVol);
      fscanf(fastInpFile, "OptStaggerFastVol: %d\n", &OptStaggerFastVol);
      fscanf(fastInpFile, "OptCreateFastPF: %d\n", &OptCreateFastPF);
      fscanf(fastInpFile, "OptReadFastPF: %d\n", &OptReadFastPF);
      fscanf(fastInpFile, "NbPtSkip: %d\n", &NbPtSkip);
      fscanf(fastInpFile, "NbStgPtSkip: %d\n", &NbStgPtSkip);
      fscanf(fastInpFile, "LX: %le\n", &FastVol.LX);
      fscanf(fastInpFile, "LY: %le\n", &FastVol.LY);
      fscanf(fastInpFile, "LZ: %le\n", &FastVol.LZ);
      fscanf(fastInpFile, "CornerX: %le\n", &FastVol.CrnrX);
      fscanf(fastInpFile, "CornerY: %le\n", &FastVol.CrnrY);
      fscanf(fastInpFile, "CornerZ: %le\n", &FastVol.CrnrZ);
      fscanf(fastInpFile, "YStagger: %le\n", &FastVol.YStagger);
      if (!OptStaggerFastVol)
        FastVol.YStagger = 0.0;  // ignore any non-zero value
      fscanf(fastInpFile, "NbOfBlocks: %d\n", &FastVol.NbBlocks);
      BlkNbXCells = ivector(1, FastVol.NbBlocks);
      BlkNbYCells = ivector(1, FastVol.NbBlocks);
      BlkNbZCells = ivector(1, FastVol.NbBlocks);
      BlkLZ = dvector(1, FastVol.NbBlocks);
      BlkCrnrZ = dvector(1, FastVol.NbBlocks);
      for (int block = 1; block <= FastVol.NbBlocks; ++block) {
        fscanf(fastInpFile, "NbOfXCells: %d\n", &BlkNbXCells[block]);
        fscanf(fastInpFile, "NbOfYCells: %d\n", &BlkNbYCells[block]);
        fscanf(fastInpFile, "NbOfZCells: %d\n", &BlkNbZCells[block]);
        fscanf(fastInpFile, "LZ: %le\n", &BlkLZ[block]);
        fscanf(fastInpFile, "CornerZ: %le\n", &BlkCrnrZ[block]);
      }  // inputs for blocks
      fscanf(fastInpFile, "NbOfOmitVols: %d\n", &FastVol.NbOmitVols);
      if (FastVol.NbOmitVols) {
        OmitVolLX = dvector(1, FastVol.NbOmitVols);
        OmitVolLY = dvector(1, FastVol.NbOmitVols);
        OmitVolLZ = dvector(1, FastVol.NbOmitVols);
        OmitVolCrnrX = dvector(1, FastVol.NbOmitVols);
        OmitVolCrnrY = dvector(1, FastVol.NbOmitVols);
        OmitVolCrnrZ = dvector(1, FastVol.NbOmitVols);
        for (int omit = 1; omit <= FastVol.NbOmitVols; ++omit) {
          fscanf(fastInpFile, "OmitVolLX: %le\n", &OmitVolLX[omit]);
          fscanf(fastInpFile, "OmitVolLY: %le\n", &OmitVolLY[omit]);
          fscanf(fastInpFile, "OmitVolLZ: %le\n", &OmitVolLZ[omit]);
          fscanf(fastInpFile, "OmitVolCornerX: %le\n", &OmitVolCrnrX[omit]);
          fscanf(fastInpFile, "OmitVolCornerY: %le\n", &OmitVolCrnrY[omit]);
          fscanf(fastInpFile, "OmitVolCornerZ: %le\n", &OmitVolCrnrZ[omit]);
        }  // inputs for OmitVols
      }    // inputs for OmitVols
      fscanf(fastInpFile, "NbOfIgnoreVols: %d\n", &FastVol.NbIgnoreVols);
      if (FastVol.NbIgnoreVols) {
        IgnoreVolLX = dvector(1, FastVol.NbIgnoreVols);
        IgnoreVolLY = dvector(1, FastVol.NbIgnoreVols);
        IgnoreVolLZ = dvector(1, FastVol.NbIgnoreVols);
        IgnoreVolCrnrX = dvector(1, FastVol.NbIgnoreVols);
        IgnoreVolCrnrY = dvector(1, FastVol.NbIgnoreVols);
        IgnoreVolCrnrZ = dvector(1, FastVol.NbIgnoreVols);
        for (int ignore = 1; ignore <= FastVol.NbIgnoreVols; ++ignore) {
          fscanf(fastInpFile, "IgnoreVolLX: %le\n", &IgnoreVolLX[ignore]);
          fscanf(fastInpFile, "IgnoreVolLY: %le\n", &IgnoreVolLY[ignore]);
          fscanf(fastInpFile, "IgnoreVolLZ: %le\n", &IgnoreVolLZ[ignore]);
          fscanf(fastInpFile, "IgnoreVolCornerX: %le\n",
                 &IgnoreVolCrnrX[ignore]);
          fscanf(fastInpFile, "IgnoreVolCornerY: %le\n",
                 &IgnoreVolCrnrY[ignore]);
          fscanf(fastInpFile, "IgnoreVolCornerZ: %le\n",
                 &IgnoreVolCrnrZ[ignore]);
        }  // inputs for IgnoreVols
      }    // inputs for IgnoreVols
      for (int ignore = 1; ignore <= FastVol.NbIgnoreVols; ++ignore) {
        printf("IgnoreVolLX: %le\n", IgnoreVolLX[ignore]);
        printf("IgnoreVolLY: %le\n", IgnoreVolLY[ignore]);
        printf("IgnoreVolLZ: %le\n", IgnoreVolLZ[ignore]);
        printf("IgnoreVolCornerX: %le\n", IgnoreVolCrnrX[ignore]);
        printf("IgnoreVolCornerY: %le\n", IgnoreVolCrnrY[ignore]);
        printf("IgnoreVolCornerZ: %le\n", IgnoreVolCrnrZ[ignore]);
      }  // inputs for IgnoreVols
      fclose(fastInpFile);
    }  // else fastInpFile
  }    // if OptFastVol

  printf("neBEM initialized ...\n");
  fflush(stdout);
  sleep(3);  // wait for three seconds so that the user gets time to react

  neBEMState = 1;  // state 1 implied initialization of neBEM completed

  // announce success - later, add the name of the calling code
  return (0);
}  // neBEMInitialize ends

// Reads geometry details
// Note that reflection and periodicity (reflection) information is gathered
// using the neBEMGetPeriodicities() function, called from here.
// Repetition variables were introduced to facilitate GarfieldInterface.
// Now Garfield can pass parameters directly to neBEM and, hence, these
// variables have become redundant.
// They are likely to be removed in near future.
int neBEMReadGeometry(void) {
  int dbgFn = 0;
  int fstatus;

  startClock = clock();

  // For a model that was defined before and for which data was stored in a file
  if ((!NewModel) && (!NewBC) && (OptStorePrimitives)) {
    fstatus = ReadPrimitives();
    if (fstatus) {
      neBEMMessage("neBEMReadGeometry - problem reading stored Primitives.\n");
      return -1;
    }
    neBEMState = 3;  // primitives read in after initialization and Nbs
    return 0;
  }

  printf("geometry inputs ...\n");
  if (neBEMState != 1) {
    printf("reading geometry possible only after initialization ...\n");
    return -1;
  }

  NbPrimitives = neBEMGetNbPrimitives();
  OrgnlNbPrimitives = NbPrimitives;
  if (NbPrimitives == 0) {
    // nothing to do - return control to calling routine
    neBEMMessage("neBEMReadGeometry - no primitive.\n");
    return (-1);  // for the time being
  }

  NbSurfs = 0;
  NbWires = 0;
  neBEMState = 2;

  // Allocate memory for storing the geometry primitives till the elements are
  // created.
  // Explicit storage of these variables related to primitives may not be
  // necessary if the elements are created immediately after a primitive is read
  // in.
  // MaxNbVertices = 4; // value specified through SetDefaults or init files

  // neBEM has been initialized, NbPrimitives set
  PrimType = ivector(1, NbPrimitives);
  NbVertices = ivector(1, NbPrimitives);
  OrgnlToEffPrim = imatrix(1, NbPrimitives, 0, 2);  // 0 init, 1 intrfc, 2 rmv
  XVertex = dmatrix(1, NbPrimitives, 0, MaxNbVertices - 1);
  YVertex = dmatrix(1, NbPrimitives, 0, MaxNbVertices - 1);
  ZVertex = dmatrix(1, NbPrimitives, 0, MaxNbVertices - 1);
  XNorm = dvector(1, NbPrimitives);
  YNorm = dvector(1, NbPrimitives);
  ZNorm = dvector(1, NbPrimitives);
  PrimLX = dvector(1, NbPrimitives);
  PrimLZ = dvector(1, NbPrimitives);
  Radius = dvector(1, NbPrimitives);  // can lead to a little memory misuse
  PrimOriginX = dvector(1, NbPrimitives);
  PrimOriginY = dvector(1, NbPrimitives);
  PrimOriginZ = dvector(1, NbPrimitives);
  PrimDC = (DirnCosn3D *)malloc(NbPrimitives * sizeof(DirnCosn3D));
  VolRef1 = ivector(1, NbPrimitives);
  VolRef2 = ivector(1, NbPrimitives);
  NbSurfSegX = ivector(1, NbPrimitives);
  NbSurfSegZ = ivector(1, NbPrimitives);
  NbWireSeg = ivector(1, NbPrimitives);  // little memory misuse
  InterfaceType = ivector(1, NbPrimitives);
  Epsilon1 = dvector(1, NbPrimitives);
  Epsilon2 = dvector(1, NbPrimitives);
  Lambda = dvector(1, NbPrimitives);
  ApplPot = dvector(1, NbPrimitives);
  ApplCh = dvector(1, NbPrimitives);
  PeriodicTypeX = ivector(1, NbPrimitives);
  PeriodicTypeY = ivector(1, NbPrimitives);
  PeriodicTypeZ = ivector(1, NbPrimitives);
  PeriodicInX = ivector(1, NbPrimitives);
  PeriodicInY = ivector(1, NbPrimitives);
  PeriodicInZ = ivector(1, NbPrimitives);
  XPeriod = dvector(1, NbPrimitives);
  YPeriod = dvector(1, NbPrimitives);
  ZPeriod = dvector(1, NbPrimitives);
  MirrorTypeX = ivector(1, NbPrimitives);
  MirrorTypeY = ivector(1, NbPrimitives);
  MirrorTypeZ = ivector(1, NbPrimitives);
  MirrorDistXFromOrigin = dvector(1, NbPrimitives);
  MirrorDistYFromOrigin = dvector(1, NbPrimitives);
  MirrorDistZFromOrigin = dvector(1, NbPrimitives);
  BndPlaneInXMin = ivector(1, NbPrimitives);
  BndPlaneInXMax = ivector(1, NbPrimitives);
  BndPlaneInYMin = ivector(1, NbPrimitives);
  BndPlaneInYMax = ivector(1, NbPrimitives);
  BndPlaneInZMin = ivector(1, NbPrimitives);
  BndPlaneInZMax = ivector(1, NbPrimitives);
  XBndPlaneInXMin = dvector(1, NbPrimitives);
  XBndPlaneInXMax = dvector(1, NbPrimitives);
  YBndPlaneInYMin = dvector(1, NbPrimitives);
  YBndPlaneInYMax = dvector(1, NbPrimitives);
  ZBndPlaneInZMin = dvector(1, NbPrimitives);
  ZBndPlaneInZMax = dvector(1, NbPrimitives);
  VBndPlaneInXMin = dvector(1, NbPrimitives);
  VBndPlaneInXMax = dvector(1, NbPrimitives);
  VBndPlaneInYMin = dvector(1, NbPrimitives);
  VBndPlaneInYMax = dvector(1, NbPrimitives);
  VBndPlaneInZMin = dvector(1, NbPrimitives);
  VBndPlaneInZMax = dvector(1, NbPrimitives);
  NbElmntsOnPrim = ivector(1, NbPrimitives);
  ElementBgn = ivector(1, NbPrimitives);
  ElementEnd = ivector(1, NbPrimitives);
  AvChDen = dvector(1, NbPrimitives);
  AvAsgndChDen = dvector(1, NbPrimitives);

  // Loop over the primitives - major loop
  int nvertex, volref1, volref2, volmax = 0;
#ifdef __cplusplus
  std::vector<double> xvert(MaxNbVertices, 0.);
  std::vector<double> yvert(MaxNbVertices, 0.);
  std::vector<double> zvert(MaxNbVertices, 0.);
#else
  double xvert[MaxNbVertices], yvert[MaxNbVertices], zvert[MaxNbVertices];
#endif
  double xnorm, ynorm, znorm;  // in case of wire , radius is read as xnorm
  for (int prim = 1; prim <= NbPrimitives; ++prim) {
#ifdef __cplusplus
    fstatus = neBEMGetPrimitive(prim, &nvertex, xvert.data(), yvert.data(),
                                zvert.data(), &xnorm, &ynorm, &znorm, &volref1,
                                &volref2);
#else
    fstatus = neBEMGetPrimitive(prim, &nvertex, xvert, yvert, zvert,  // arrays
                                &xnorm, &ynorm, &znorm, &volref1, &volref2);
#endif
    if (fstatus != 0) {
      neBEMMessage("neBEMReadGeometry - neBEMGetPrimitve");
      return -1;
    }
    if (volmax < volref1) {
      volmax = volref1;
    }  // maxm nb of volumes
    if (volmax < volref2) {
      volmax = volref2;
    }  // maxm nb of volumes

    if (nvertex > MaxNbVertices) {
      printf("Number of vertices for primitive %d exceeds %d!\n", prim,
             MaxNbVertices);
      printf("Returning to garfield ...\n");
      return (-1);
    }

    PrimType[prim] = nvertex;  // wire:2, triangle:3, rectangle:4
    NbVertices[prim] = nvertex;
    for (int vert = 0; vert < NbVertices[prim]; ++vert) {
      XVertex[prim][vert] = xvert[vert];
      YVertex[prim][vert] = yvert[vert];
      ZVertex[prim][vert] = zvert[vert];
    }
    if (PrimType[prim] == 2) {
      // wire
      XNorm[prim] = 0.0;  // modulus not 1 - an absurd trio!
      YNorm[prim] = 0.0;
      ZNorm[prim] = 0.0;
      Radius[prim] = xnorm;
    }
    if ((PrimType[prim] == 3) || (PrimType[prim] == 4)) {
      XNorm[prim] = xnorm;
      YNorm[prim] = ynorm;
      ZNorm[prim] = znorm;
      Radius[prim] = 0.0;  // absurd radius!
    }
    VolRef1[prim] = volref1;
    VolRef2[prim] = volref2;

    // feedback for user begins - suppress later
    if (OptPrintPrimaryDetails) {
      printf("neBEM:\tprimitive %d between volumes %d, %d has %d vertices\n",
             prim, volref1, volref2, nvertex);
      for (int ivertex = 0; ivertex < nvertex; ivertex++) {
        printf("\tnode %d (%g,%g,%g)\n", ivertex, xvert[ivertex],
               yvert[ivertex], zvert[ivertex]);
      }
      printf("\tnormal vector: (%g,%g,%g)\n", xnorm, ynorm, znorm);
    }
    // feedback for user ends - suppress later

    // Now look for the volume related information for this primitve
    // This is obtained from the specified volume references
    // volref1 refers to the volume itself
    // volref2 describes the volume in the direction of the +ve normal
    // Note that materials from 1 to 10 are conductors and
    // 										 from 11 to 20
    // are dielectrics
    int shape1, material1, boundarytype1;
    double epsilon1, potential1, charge1;
    if (volref1 == -1) {
      // Must be an error, since no device is made of vacuum
      neBEMMessage("neBEMReadGeometry - volref1 = -1!");
      return -1;
    } else {
      neBEMVolumeDescription(volref1, &shape1, &material1, &epsilon1,
                             &potential1, &charge1, &boundarytype1);
    }
    if (OptPrintVolumeDetails) {
      printf("\tvolref1: %d\n", volref1);
      printf("\t\tboundarytype1: %d, shape1: %d, material1: %d\n",
             boundarytype1, shape1, material1);
      printf("\t\tepsilon1: %lg, potential1: %lg, charge1: %lg\n", epsilon1,
             potential1, charge1);
    }
    // in the -ve normal direction - properties of the external volume
    int shape2, material2, boundarytype2;
    double epsilon2, potential2, charge2;
    if (volref2 == -1) {
      shape2 = 0;
      material2 = 11;
      epsilon2 = 1.0;
      potential2 = 0.0;
      charge2 = 0.0;
      boundarytype2 = 0;
    } else {
      neBEMVolumeDescription(volref2, &shape2, &material2, &epsilon2,
                             &potential2, &charge2, &boundarytype2);
    }
    if (OptPrintVolumeDetails) {
      printf("\tvolref2: %d\n", volref2);
      printf("\t\tboundarytype2: %d, shape2: %d, material2: %d\n",
             boundarytype2, shape2, material2);
      printf("\t\tepsilon2: %lg, potential2: %lg, charge2: %lg\n", epsilon2,
             potential2, charge2);
    }

    // Put default values to variables that depend on the interface type
    // Is there any risk involved in putting these defaults?
    // At present, they even seem necessary. For example, for floating
    // conductors or dielectric-dielectric interface, the formulation requires
    // that the RHS is zero (may be modified by the effect of known charges).
    Epsilon1[prim] = epsilon1;
    Epsilon2[prim] = epsilon2;  // 1: self, 2: external
    ApplPot[prim] = 0.0;
    Lambda[prim] = 0.0;
    ApplCh[prim] = 0.0;

    // BoundaryTypes:
    // --------------
    // Vacuum: 0
    // Conductor at specified potential: 1
    // Conductor with a specified charge: 2
    // Floating conductor (zero charge, perpendicular E): 3
    // Dielectric intrface (plastic-plastic) without "manual" charge: 4
    // Dielectric with surface charge (plastic-gas, typically): 5
    // Symmetry boundary, E parallel: 6 (may not be necessary)
    // Symmetry boundary, E perpendicular: 7 (may not be necessary)

    // InterfaceType:
    // --------------
    // To be skipped: 0
    // Conductor-dielectric: 1
    // Conductor with known charge: 2
    // Conductor at floating potential: 3
    // Dielectric-dielectric: 4
    // Dielectric with given surface charge: 5
    // Check dielectric-dielectric formulation in
    // (NumSolnOfBIEforMolES_JPBardhan.pdf):
    // Numerical solution of boundary-integral equations for molecular
    // electrostatics,
    // by Jaydeep P. Bardhan,
    // THE JOURNAL OF CHEMICAL PHYSICS 130, 094102 (2009)

    switch (boundarytype1) {  // the volume itself is volref1
      case 1:                 // conductor at specified potential
        if (boundarytype2 == 0 || boundarytype2 == 4) {
          // dielectric-conductor
          InterfaceType[prim] = 1;
          ApplPot[prim] = potential1;
        } else if (boundarytype2 == 1) {
          // conductor-conductor
          if (fabs(potential1 - potential2)  // same potential
              < 1e-6 * (1 + fabs(potential1) + fabs(potential2))) {
            printf("neBEMReadGeometry: identical potentials; skipped.\n");
            printf("Primitive skipped: #%d\n", prim);
            InterfaceType[prim] = 0;
          } else {
            // different potentials
            printf("neBEMReadGeometry: different potentials; rejected.\n");
            return -1;
          }
        } else {
          // conductor-unknown
          printf(
              "neBEMReadGeometry: conductor at given potential; rejected.\n");
          return -1;
        }
        break;

      case 2:  // conductor with a specified charge
        if ((boundarytype2 == 0) || (boundarytype2 == 4)) {
          // conductor-dielectric
          InterfaceType[prim] = 2;
          ApplCh[prim] = charge1;
        } else {
          printf("neBEMReadGeometry: charged conductor; rejected.\n");
          return -1;
        }
        break;

      case 3:  // floating conductor (zero charge, perpendicular E)
        if ((boundarytype2 == 0) || (boundarytype2 == 4)) {
          // conductor-dielectric
          InterfaceType[prim] = 3;
          if (!NbFloatingConductors)   // assuming only one floating conductor
            NbFloatingConductors = 1;  // in the system
        } else {
          printf("neBEMReadGeometry: floating conductor; rejected.\n");
          return -1;
        }
        break;

      case 4:  // dielectric interface (plastic-plastic) without "manual" charge
        if (boundarytype2 == 0) {
          // dielectric-vacuum
          // epsilon1 is self dielectric-constant
          // epsilon2 is towards positive normal
          InterfaceType[prim] = 4;
          Lambda[prim] = (epsilon1 - epsilon2) / (epsilon1 + epsilon2);
          // consistent with Bardhan's eqn 16 where (1 / (2*Lambda)) is used
        } else if (boundarytype2 == 1) {
          // dielectric-conductor
          InterfaceType[prim] = 1;  // conductor at known potential
          ApplPot[prim] = potential2;
        } else if (boundarytype2 == 2) {
          // dielectric-conductor
          InterfaceType[prim] = 2;  // conductor with known charge
          ApplCh[prim] = charge2;
        } else if (boundarytype2 == 3) {
          // dielectric-conductor
          InterfaceType[prim] = 3;  // conductor at floating potential
        } else if (boundarytype2 == 4) {
          // dielectric-dielectric
          if (fabs(epsilon1 - epsilon2) <
              1e-6 * (1 + fabs(epsilon1) + fabs(epsilon2))) {
            // identical dielectrica
            printf(
                "neBEMReadGeometry: between identical dielectrica; skipd.\n");
            printf("Primitive skipped: #%d\n", prim);
            InterfaceType[prim] = 0;
          } else {
            // distinctly different dielectrica
            // epsilon1 is self dielectric-constant
            // epsilon2 towards positive normal
            InterfaceType[prim] = 4;
            Lambda[prim] = (epsilon1 - epsilon2) / (epsilon1 + epsilon2);
            // consistent with Bardhan's paper (1 / Lambda)
          }
        } else if (boundarytype2 == 5) {
          // dielectric-dielectric with charge
          if (fabs(epsilon1 - epsilon2)  // identical dielectrica
              < 1e-6 * (1 + fabs(epsilon1) + fabs(epsilon2))) {
            printf(
                "neBEMReadGeometry: between identical dielectrica; skipped.\n");
            printf("Primitive skipped: #%d\n", prim);
            InterfaceType[prim] = 0;
          } else {
            // distinctly different dielectrica
            InterfaceType[prim] = 5;  // epsilon2 towards positive normal
            ApplCh[prim] = charge2;
            Lambda[prim] = (epsilon1 - epsilon2) / (epsilon1 + epsilon2);
          }
        }       // if-else if boundarytypes 0 and 4
        else {  // dielectric-unknown
          printf("neBEMReadGeometry: unknown dielectric; rejected.\n");
          return -1;
        }
        break;

      case 5:  // dielectric with surface charge (plastic-gas, typically)
        if (boundarytype2 == 0) {   // dielectric-vacuum
          InterfaceType[prim] = 5;  // epsilon2 is towards +ve normal
          ApplCh[prim] = charge1;
          Lambda[prim] = (epsilon1 - epsilon2) / (epsilon1 + epsilon2);
          // consistent with Bardhan's paper (1 / Lambda)
        } else if (boundarytype2 == 4) {  // dielectric-dielectric
          if (fabs(epsilon1 - epsilon2) <
              1e-6 * (1 + fabs(epsilon1) + fabs(epsilon2))) {
            // identical dielectrica
            printf(
                "neBEMReadGeometry: between identical dielectrica; skipd.\n");
            printf("Primitive skipped: #%d\n", prim);
            InterfaceType[prim] = 0;
          } else {
            // distinctly different dielectrica
            InterfaceType[prim] = 5;  // epsilon2 towards positive normal
            ApplCh[prim] = charge1;
            Lambda[prim] = (epsilon1 - epsilon2) / (epsilon1 + epsilon2);
            // consistent with Bardhan's paper (1 / Lambda)
          }
        }  // if-else if boundarytypes 0 and 4
        else {
          printf(
              "neBEMReadGeometry: charged dielectric adjacent to a conductor; "
              "rejected.\n");
          return -1;
        }
        break;

      case 6:  // symmetry boundary, E parallel
        if (boundarytype2 == 0) {
          InterfaceType[prim] = 6;
        } else {
          printf("neBEMReadGeometry: E-parallel symmetry; rejected.\n");
          return -1;
        }
        break;

      case 7:  // symmetry boundary, E perpendicular
        if (boundarytype2 == 0) {
          InterfaceType[prim] = 7;
        } else {
          printf("neBEMReadGeometry: E-perpendicular symmetry; rejected.\n");
          return -1;
        }
        break;

      default:
        printf("neBEMReadGeometry: Boundary type 1: %d\n", boundarytype1);
        printf("neBEMReadGeometry: Boundary type 2: %d\n", boundarytype2);
        printf("neBEMReadGeometry:          out of range ... exiting.\n");
        return -1;
    }  // switch boundarytype1 ends

    if (OptPrintVolumeDetails) {
      printf(
          "\tType: %d, ApplPot: %lg, Epsilon1: %lg, Epsilon2: %lg, Lambda: "
          "%lg, ApplCh: %lg\n",
          InterfaceType[prim], ApplPot[prim], Epsilon1[prim], Epsilon2[prim],
          Lambda[prim], ApplCh[prim]);
    }

    // Read the periodicities
    // Note that mirror has been taken care of separately below
    // ix: PeriodicTypeX (1 - simple, 2 - mirror, 3 - axial, 4 - rotation)
    // jx: PeriodicInX (Number of copies internal to neBEM)
    // sx: XPeriod
    // NOTE: A change in this part of the algorithm is likely to affect
    // src/Solve/neBEM.c (LHMatrix) and
    // src/Solve/ComputeProperties.c (PFAtPoint and WtFldPFAtPoint)
    {
      int ix, iy, iz;
      int jx, jy, jz;
      double sx, sy, sz;
      fstatus = neBEMGetPeriodicities(prim, &ix, &jx, &sx, &iy, &jy, &sy, &iz,
                                      &jz, &sz);
      if (fstatus != 0) {
        neBEMMessage("neBEMReadGeometry - neBEMGetPeriodicities");
        return -1;
      }
      if (jx < 0) jx = 0;
      if (jy < 0) jy = 0;
      if (jz < 0) jz = 0;

      PeriodicTypeX[prim] = ix;
      PeriodicTypeY[prim] = iy;
      PeriodicTypeZ[prim] = iz;
      if (0) {
        printf("For primitive: %d\n", prim);
        printf("\tPeriodicTypeX: %d, PeriodicTypeY: %d, PeriodicTypeZ: %d\n",
               ix, iy, iz);
        printf("\tPeriodicInX: %d, PeriodicInY: %d, PeriodicInZ: %d\n", jx, jy,
               jz);
        printf("\tXPeriod: %lg, YPeriod: %lg, ZPeriod: %lg\n", sx, sy, sz);
      }
      if (ix > 0) {
        // These checks need to be done separately. Otherwise, there is
        // a possibility that non-zero values of PeriodicIn* and *Period
        // are used throughout the code despite PeriodicType* is 0
        PeriodicInX[prim] = jx;
        XPeriod[prim] = sx;
      } else {
        PeriodicInX[prim] = 0;
        XPeriod[prim] = 0.0;
      }
      if (iy > 0) {
        PeriodicInY[prim] = jy;
        YPeriod[prim] = sy;
      } else {
        PeriodicInY[prim] = 0;
        YPeriod[prim] = 0.0;
      }
      if (iz > 0) {
        PeriodicInZ[prim] = jz;
        ZPeriod[prim] = sz;
      } else {
        PeriodicInZ[prim] = 0;
        ZPeriod[prim] = 0.0;
      }
    }  // read periodicity information

    // Read mirror information
    // Mirror can be in perpendicular to any of the three cartesian axes
    // There can be more than one mirror at the same time
    // MirrorType (1 - charge density -ve of original, equivalent to method of
    // images, 2 - charge density same as original)
    {
      int ix, iy, iz;
      int jx, jy, jz;  // not used at present
      double sx, sy, sz;
      fstatus =
          neBEMGetMirror(prim, &ix, &jx, &sx, &iy, &jy, &sy, &iz, &jz, &sz);
      if (fstatus != 0) {
        neBEMMessage("neBEMReadGeometry - neBEMGetMirror");
        return -1;
      }
      if (jx < 0) jx = 0;
      if (jy < 0) jy = 0;
      if (jz < 0) jz = 0;

      MirrorTypeX[prim] = ix;
      MirrorTypeY[prim] = iy;
      MirrorTypeZ[prim] = iz;
      if (0) {
        printf("For primitive: %d\n", prim);
        printf("\tMirrorTypeX: %d, MirrorTypeY: %d, MirrorTypeZ: %d\n", ix, iy,
               iz);
        printf("\tNOT USED ==> MirrorInX: %d, MirrorInY: %d, MirrorInZ: %d\n",
               jx, jy, jz);
        printf("\tMirrorDistX: %lg, MirrorDistY: %lg, MirrorDistZ: %lg\n", sx,
               sy, sz);
        getchar();
      }
      if (ix > 0) {
        // printf("neBEMReadGeometry: Mirror have been requested.\n");
        MirrorDistXFromOrigin[prim] = sx;  // assumed to pass through the origin
      } else {
        MirrorDistXFromOrigin[prim] = 0.0;  // pass through the origin
      }
      if (iy > 0) {
        // printf("neBEMReadGeometry: Mirror have been requested.\n");
        MirrorDistYFromOrigin[prim] = sy;
      } else {
        MirrorDistYFromOrigin[prim] = 0.0;
      }
      if (iz > 0) {
        // printf("neBEMReadGeometry: Mirror have been requested.\n");
        MirrorDistZFromOrigin[prim] = sz;
      } else {
        MirrorDistZFromOrigin[prim] = 0.0;
      }
    }  // read mirror information

    // Information on bounding planes
    // ixmin=0: lower x-plane does not exist
    // ixmin=1: lower x-plane does exist
    // cxmin: coordinate of lower x-plane
    // vxmin: potential of lower x-plane
    // Similar for ixmax, iymin, iymax, izmin, izmax
    int ixmin, ixmax, iymin, iymax, izmin, izmax;
    double cxmin, cxmax, cymin, cymax, czmin, czmax;
    double vxmin, vxmax, vymin, vymax, vzmin, vzmax;
    fstatus = neBEMGetBoundingPlanes(
        &ixmin, &cxmin, &vxmin, &ixmax, &cxmax, &vxmax, &iymin, &cymin, &vymin,
        &iymax, &cymax, &vymax, &izmin, &czmin, &vzmin, &izmax, &czmax, &vzmax);
    if (fstatus != 0) {
      neBEMMessage("neBEMReadGeometry - neBEMGetBoundingPlanes");
      return -1;
    }
    BndPlaneInXMin[prim] = ixmin;
    BndPlaneInXMax[prim] = ixmax;
    BndPlaneInYMin[prim] = iymin;
    BndPlaneInYMax[prim] = iymax;
    BndPlaneInZMin[prim] = izmin;
    BndPlaneInZMax[prim] = izmax;
    if (ixmin) {
      XBndPlaneInXMin[prim] = cxmin;
      VBndPlaneInXMin[prim] = vxmin;
    } else {
      XBndPlaneInXMin[prim] = 0.0;
      VBndPlaneInXMin[prim] = 0.0;
    }
    if (ixmax > 0) {
      XBndPlaneInXMax[prim] = cxmax;
      VBndPlaneInXMax[prim] = vxmax;
    } else {
      XBndPlaneInXMax[prim] = 0.0;
      VBndPlaneInXMax[prim] = 0.0;
    }
    if (iymin > 0) {
      YBndPlaneInYMin[prim] = cymin;
      VBndPlaneInYMin[prim] = vymin;
    } else {
      YBndPlaneInYMin[prim] = 0.0;
      VBndPlaneInYMin[prim] = 0.0;
    }
    if (iymax > 0) {
      YBndPlaneInYMax[prim] = cymax;
      VBndPlaneInYMax[prim] = vymax;
    } else {
      YBndPlaneInYMax[prim] = 0.0;
      VBndPlaneInYMax[prim] = 0.0;
    }
    if (izmin > 0) {
      ZBndPlaneInZMin[prim] = czmin;
      VBndPlaneInZMin[prim] = vzmin;
    } else {
      ZBndPlaneInZMin[prim] = 0.0;
      VBndPlaneInZMin[prim] = 0.0;
    }
    if (izmax > 0) {
      ZBndPlaneInZMax[prim] = czmax;
      VBndPlaneInZMax[prim] = vzmax;
    } else {
      ZBndPlaneInZMax[prim] = 0.0;
      VBndPlaneInZMax[prim] = 0.0;
    }
  }  // loop over the primitives - major loop

  VolMax = volmax;              // Maximum nb of volumes in the problem
  volRef = ivector(0, VolMax);  // variables to store volume related infomration
  volShape = ivector(0, VolMax);
  volMaterial = ivector(0, VolMax);
  volEpsilon = dvector(0, VolMax);
  volPotential = dvector(0, VolMax);
  volCharge = dvector(0, VolMax);
  volBoundaryType = ivector(0, VolMax);
  for (int volref = 0; volref <= VolMax; ++volref) {
    neBEMVolumeDescription(volref, &volShape[volref], &volMaterial[volref],
                           &volEpsilon[volref], &volPotential[volref],
                           &volCharge[volref], &volBoundaryType[volref]);
    if (dbgFn) {
      printf("volref: %d\n", volref);
      printf("shape: %d,  material: %d\n", volShape[volref],
             volMaterial[volref]);
      printf("eps: %lg,  pot: %lg\n", volEpsilon[volref], volPotential[volref]);
      printf("q: %lg,  type: %d\n", volCharge[volref], volBoundaryType[volref]);
    }
  }

  // Ignore unnecessary primitives from the final count
  // Ideally, all the removal conditions for a primitive should be checked in
  // one loop and the list should be updated in one single go.
  {
    for (int prim = 1; prim <= NbPrimitives; ++prim) {
      OrgnlToEffPrim[prim][0] = prim;
      OrgnlToEffPrim[prim][1] = prim;
      OrgnlToEffPrim[prim][2] = prim;
    }

    {  // Skip primitive
      // Remove skipped primitives having InterfaceType == 0.
      // Also remove primitives having too small dimensions.
      int NbSkipped = 0, effprim;
      double DVertex[4], minDVertex = 0.0;  // maximum number of vertices is 4
      for (int prim = 1; prim <= NbPrimitives; ++prim) {
        effprim = prim - NbSkipped;

        // Check dimensions of the primitive
        for (int vert = 0; vert < NbVertices[prim] - 1; ++vert) {
          DVertex[vert] =
              sqrt(((XVertex[prim][vert + 1] - XVertex[prim][vert]) *
                    (XVertex[prim][vert + 1] - XVertex[prim][vert])) +
                   ((YVertex[prim][vert + 1] - YVertex[prim][vert]) *
                    (YVertex[prim][vert + 1] - YVertex[prim][vert])) +
                   ((ZVertex[prim][vert + 1] - ZVertex[prim][vert]) *
                    (ZVertex[prim][vert + 1] - ZVertex[prim][vert])));
          if (vert == 0)
            minDVertex = DVertex[vert];
          else {
            if (DVertex[vert] < minDVertex) minDVertex = DVertex[vert];
          }
        }

        if ((InterfaceType[prim]) && (minDVertex > MINDIST)) {
          OrgnlToEffPrim[prim][1] = effprim;
          OrgnlToEffPrim[prim][2] = effprim;
          PrimType[effprim] = PrimType[prim];
          NbVertices[effprim] = NbVertices[prim];
          for (int vert = 0; vert < NbVertices[effprim]; ++vert) {
            XVertex[effprim][vert] = XVertex[prim][vert];
            YVertex[effprim][vert] = YVertex[prim][vert];
            ZVertex[effprim][vert] = ZVertex[prim][vert];
          }
          if (PrimType[effprim] == 2)  // wire
          {
            XNorm[effprim] = 0.0;  // modulus not 1 - an absurd trio!
            YNorm[effprim] = 0.0;
            ZNorm[effprim] = 0.0;
            Radius[effprim] = Radius[prim];
          }
          if ((PrimType[effprim] == 3) || (PrimType[effprim] == 4)) {
            XNorm[effprim] = XNorm[prim];
            YNorm[effprim] = YNorm[prim];
            ZNorm[effprim] = ZNorm[prim];
            Radius[effprim] = 0.0;  // absurd radius!
          }
          VolRef1[effprim] = VolRef1[prim];
          VolRef2[effprim] = VolRef2[prim];

          InterfaceType[effprim] = InterfaceType[prim];
          Epsilon1[effprim] = Epsilon1[prim];
          Epsilon2[effprim] = Epsilon2[prim];
          Lambda[effprim] = Lambda[prim];
          ApplPot[effprim] = ApplPot[prim];
          ApplCh[effprim] = ApplCh[prim];
          PeriodicTypeX[effprim] = PeriodicTypeX[prim];
          PeriodicTypeY[effprim] = PeriodicTypeY[prim];
          PeriodicTypeZ[effprim] = PeriodicTypeZ[prim];
          PeriodicInX[effprim] = PeriodicInX[prim];
          PeriodicInY[effprim] = PeriodicInY[prim];
          PeriodicInZ[effprim] = PeriodicInZ[prim];
          XPeriod[effprim] = XPeriod[prim];
          YPeriod[effprim] = YPeriod[prim];
          ZPeriod[effprim] = ZPeriod[prim];
          MirrorTypeX[effprim] = MirrorTypeX[prim];
          MirrorTypeY[effprim] = MirrorTypeY[prim];
          MirrorTypeZ[effprim] = MirrorTypeZ[prim];
          MirrorDistXFromOrigin[effprim] = MirrorDistXFromOrigin[prim];
          MirrorDistYFromOrigin[effprim] = MirrorDistYFromOrigin[prim];
          MirrorDistZFromOrigin[effprim] = MirrorDistZFromOrigin[prim];
          BndPlaneInXMin[effprim] = BndPlaneInXMin[prim];
          BndPlaneInXMax[effprim] = BndPlaneInXMax[prim];
          BndPlaneInYMin[effprim] = BndPlaneInYMin[prim];
          BndPlaneInYMax[effprim] = BndPlaneInYMax[prim];
          BndPlaneInZMin[effprim] = BndPlaneInZMin[prim];
          BndPlaneInZMax[effprim] = BndPlaneInZMax[prim];
          XBndPlaneInXMin[effprim] = XBndPlaneInXMin[prim];
          XBndPlaneInXMax[effprim] = XBndPlaneInXMax[prim];
          YBndPlaneInYMin[effprim] = YBndPlaneInYMin[prim];
          YBndPlaneInYMax[effprim] = YBndPlaneInYMax[prim];
          ZBndPlaneInZMin[effprim] = ZBndPlaneInZMin[prim];
          ZBndPlaneInZMax[effprim] = ZBndPlaneInZMax[prim];
          VBndPlaneInXMin[effprim] = VBndPlaneInXMin[prim];
          VBndPlaneInXMax[effprim] = VBndPlaneInXMax[prim];
          VBndPlaneInYMin[effprim] = VBndPlaneInYMin[prim];
          VBndPlaneInYMax[effprim] = VBndPlaneInYMax[prim];
          VBndPlaneInZMin[effprim] = VBndPlaneInZMin[prim];
          VBndPlaneInZMax[effprim] = VBndPlaneInZMax[prim];
        }  // InterfaceType
        else {
          OrgnlToEffPrim[prim][1] = 0;  // removed from the list
          OrgnlToEffPrim[prim][2] = 0;
          ++NbSkipped;
          if (DebugLevel == 101) {
            printf("Skipped primitive %d, InterfaceType: %d, minDVertex: %lg\n",
                   prim, InterfaceType[prim], minDVertex);
          }
        }
      }  // loop over primitives to remove the skipped primitives
      NbPrimitives -= NbSkipped;
      printf("Number of primitives skipped: %d, Effective NbPrimitives: %d\n",
             NbSkipped, NbPrimitives);
    }  // Skip primitives

    if (OptRmPrim) {
      int NbRmPrims;
      FILE *rmprimFile = fopen("neBEMInp/neBEMRmPrim.inp", "r");
      if (rmprimFile == NULL) {
        printf("neBEMRmPrim.inp absent ... assuming defaults ...\n");
        NbRmPrims = 0;
      } else {
        fscanf(rmprimFile, "NbRmPrims: %d\n", &NbRmPrims);
        if (NbRmPrims) {
          int tint;
#ifdef __cplusplus
          std::vector<double> rmXNorm(NbRmPrims + 1, 0.);
          std::vector<double> rmYNorm(NbRmPrims + 1, 0.);
          std::vector<double> rmZNorm(NbRmPrims + 1, 0.);
          std::vector<double> rmXVert(NbRmPrims + 1, 0.);
          std::vector<double> rmYVert(NbRmPrims + 1, 0.);
          std::vector<double> rmZVert(NbRmPrims + 1, 0.);
#else
          double rmXNorm[NbRmPrims + 1], rmYNorm[NbRmPrims + 1];
          double rmZNorm[NbRmPrims + 1];
          double rmXVert[NbRmPrims + 1], rmYVert[NbRmPrims + 1];
          double rmZVert[NbRmPrims + 1];
#endif
          for (int rmprim = 1; rmprim <= NbRmPrims; ++rmprim) {
            fscanf(rmprimFile, "Prim: %d\n", &tint);
            fscanf(rmprimFile, "rmXNorm: %le\n", &rmXNorm[rmprim]);
            fscanf(rmprimFile, "rmYNorm: %le\n", &rmYNorm[rmprim]);
            fscanf(rmprimFile, "rmZNorm: %le\n", &rmZNorm[rmprim]);
            fscanf(rmprimFile, "rmXVert: %le\n", &rmXVert[rmprim]);
            fscanf(rmprimFile, "rmYVert: %le\n", &rmYVert[rmprim]);
            fscanf(rmprimFile, "rmZVert: %le\n", &rmZVert[rmprim]);
            printf(
                "rmprim: %d, rmXNorm: %lg, rmYNorm: %lg, rmZNorm: %lg, "
                "rmXVert: %lg, rmYVert: %lg, rmZVert: %lg\n",
                rmprim, rmXNorm[rmprim], rmYNorm[rmprim], rmZNorm[rmprim],
                rmXVert[rmprim], rmYVert[rmprim], rmZVert[rmprim]);
          }
#ifdef __cplusplus
          std::vector<int> remove(NbPrimitives + 1, 0);
#else
          int remove[NbPrimitives + 1];
#endif
          // Check updated prim list
          for (int prim = 1; prim <= NbPrimitives; ++prim) {
            remove[prim] = 0;
            if (dbgFn) {
              printf("\n\nprim: %d, XVertex: %lg, YVertex: %lg, ZVertex: %lg\n",
                     prim, XVertex[prim][0], YVertex[prim][0],
                     ZVertex[prim][0]);
              printf("XNorm: %lg, YNorm: %lg, ZNorm: %lg\n", XNorm[prim],
                     YNorm[prim], ZNorm[prim]);
            }

            for (int rmprim = 1; rmprim <= NbRmPrims; ++rmprim) {
              if (dbgFn) {
                printf(
                    "rmprim: %d, rmXVertex: %lg, rmYVertex: %lg, rmZVertex: "
                    "%lg\n",
                    rmprim, rmXVert[rmprim], rmYVert[rmprim], rmZVert[rmprim]);
                printf("rmXNorm: %lg, rmYNorm: %lg, rmZNorm: %lg\n",
                       rmXNorm[rmprim], rmYNorm[rmprim], rmZNorm[rmprim]);
              }

              // check the normal
              if ((fabs(fabs(XNorm[prim]) - fabs(rmXNorm[rmprim])) <=
                   MINDIST) &&
                  (fabs(fabs(YNorm[prim]) - fabs(rmYNorm[rmprim])) <=
                   MINDIST) &&
                  (fabs(fabs(ZNorm[prim]) - fabs(rmZNorm[rmprim])) <=
                   MINDIST)) {  // prim and rmprim are parallel
                // coplanarity check to be implemented later.
                // For the time-being, we will assume that the planes to be
                // removed have their normals parallel to a given axis. So, we
                // only check that and remove the primitive if the distace along
                // that axis match. Possible pitfall => the primitives may be
                // coplanar but non-overlapping!
                if (fabs(fabs(XNorm[prim]) - 1.0) <= 1.0e-12) {
                  // primitive || to YZ
                  if (fabs(XVertex[prim][0] - rmXVert[rmprim]) <= MINDIST) {
                    remove[prim] = 1;
                  }
                }
                if (fabs(fabs(YNorm[prim]) - 1.0) <= 1.0e-12) {
                  // primitive || to XZ
                  if (fabs(YVertex[prim][0] - rmYVert[rmprim]) <= MINDIST) {
                    remove[prim] = 1;
                  }
                }
                if (fabs(fabs(ZNorm[prim]) - 1.0) <= 1.0e-12) {
                  // primitive || to XY
                  if (fabs(ZVertex[prim][0] - rmZVert[rmprim]) <= MINDIST) {
                    remove[prim] = 1;
                  }
                }
              }  // case where prim and rmprim are parallel
              if (dbgFn) {
                printf("prim: %d, rmprim: %d, remove: %d\n", prim, rmprim,
                       remove[prim]);
              }
              if (remove[prim] == 1)
                break;  // once removed, no point checking others
            }           // for rmprim - loop over all removal specification

          }  // for prim loop over all primitives

          int NbRemoved = 0;
          char RmPrimFile[256];
          strcpy(RmPrimFile, NativePrimDir);
          strcat(RmPrimFile, "/RmPrims.info");
          FILE *fprrm = fopen(RmPrimFile, "w");
          if (fprrm == NULL) {
            printf(
                "error opening RmPrims.info file in write mode ... "
                "returning\n");
            return (-1);
          }
          // Note that some of the original primitives have already been removed
          // based on interface and dimension considerations
          int orgnlNb = 0;
          for (int prim = 1; prim <= NbPrimitives; ++prim) {
            // identify primitive number in the original list
            for (int orgnlprim = 1; orgnlprim <= OrgnlNbPrimitives;
                 ++orgnlprim) {
              if (OrgnlToEffPrim[orgnlprim][1] ==
                  prim)  // number updated for intrfc
              {
                orgnlNb = orgnlprim;
                break;
              }
            }  // loop for finding out its position in the original list

            if (remove[prim] == 1) {
              ++NbRemoved;
              OrgnlToEffPrim[orgnlNb][2] = 0;
              fprintf(fprrm, "NbRemoved: %d, Removed primitive: %d\n",
                      NbRemoved, prim);
              fprintf(fprrm, "PrimType: %d\n", PrimType[prim]);
              fprintf(fprrm, "NbVertices: %d\n", NbVertices[prim]);
              for (int vert = 0; vert < NbVertices[prim]; ++vert) {
                fprintf(fprrm, "Vertx %d: %lg, %lg, %lg\n", vert,
                        XVertex[prim][vert], YVertex[prim][vert],
                        ZVertex[prim][vert]);
              }
              fprintf(fprrm, "Normals: %lg, %lg, %lg\n", XNorm[prim],
                      YNorm[prim], ZNorm[prim]);
              continue;
            }       // if remove
            else {  // keep this one in the updated list of primitives
              int effprim = prim - NbRemoved;

              OrgnlToEffPrim[orgnlNb][2] = effprim;
              PrimType[effprim] = PrimType[prim];
              NbVertices[effprim] = NbVertices[prim];
              for (int vert = 0; vert < NbVertices[effprim]; ++vert) {
                XVertex[effprim][vert] = XVertex[prim][vert];
                YVertex[effprim][vert] = YVertex[prim][vert];
                ZVertex[effprim][vert] = ZVertex[prim][vert];
              }
              if (PrimType[effprim] == 2) {
                // wire
                XNorm[effprim] = 0.0;  // modulus not 1 - an absurd trio!
                YNorm[effprim] = 0.0;
                ZNorm[effprim] = 0.0;
                Radius[effprim] = Radius[prim];
              }
              if ((PrimType[effprim] == 3) || (PrimType[effprim] == 4)) {
                XNorm[effprim] = XNorm[prim];
                YNorm[effprim] = YNorm[prim];
                ZNorm[effprim] = ZNorm[prim];
                Radius[effprim] = 0.0;  // absurd radius!
              }
              VolRef1[effprim] = VolRef1[prim];
              VolRef2[effprim] = VolRef2[prim];

              InterfaceType[effprim] = InterfaceType[prim];
              Epsilon1[effprim] = Epsilon1[prim];
              Epsilon2[effprim] = Epsilon2[prim];
              Lambda[effprim] = Lambda[prim];
              ApplPot[effprim] = ApplPot[prim];
              ApplCh[effprim] = ApplCh[prim];
              PeriodicTypeX[effprim] = PeriodicTypeX[prim];
              PeriodicTypeY[effprim] = PeriodicTypeY[prim];
              PeriodicTypeZ[effprim] = PeriodicTypeZ[prim];
              PeriodicInX[effprim] = PeriodicInX[prim];
              PeriodicInY[effprim] = PeriodicInY[prim];
              PeriodicInZ[effprim] = PeriodicInZ[prim];
              XPeriod[effprim] = XPeriod[prim];
              YPeriod[effprim] = YPeriod[prim];
              ZPeriod[effprim] = ZPeriod[prim];
              MirrorTypeX[effprim] = MirrorTypeX[prim];
              MirrorTypeY[effprim] = MirrorTypeY[prim];
              MirrorTypeZ[effprim] = MirrorTypeZ[prim];
              MirrorDistXFromOrigin[effprim] = MirrorDistXFromOrigin[prim];
              MirrorDistYFromOrigin[effprim] = MirrorDistYFromOrigin[prim];
              MirrorDistZFromOrigin[effprim] = MirrorDistZFromOrigin[prim];
              BndPlaneInXMin[effprim] = BndPlaneInXMin[prim];
              BndPlaneInXMax[effprim] = BndPlaneInXMax[prim];
              BndPlaneInYMin[effprim] = BndPlaneInYMin[prim];
              BndPlaneInYMax[effprim] = BndPlaneInYMax[prim];
              BndPlaneInZMin[effprim] = BndPlaneInZMin[prim];
              BndPlaneInZMax[effprim] = BndPlaneInZMax[prim];
              XBndPlaneInXMin[effprim] = XBndPlaneInXMin[prim];
              XBndPlaneInXMax[effprim] = XBndPlaneInXMax[prim];
              YBndPlaneInYMin[effprim] = YBndPlaneInYMin[prim];
              YBndPlaneInYMax[effprim] = YBndPlaneInYMax[prim];
              ZBndPlaneInZMin[effprim] = ZBndPlaneInZMin[prim];
              ZBndPlaneInZMax[effprim] = ZBndPlaneInZMax[prim];
              VBndPlaneInXMin[effprim] = VBndPlaneInXMin[prim];
              VBndPlaneInXMax[effprim] = VBndPlaneInXMax[prim];
              VBndPlaneInYMin[effprim] = VBndPlaneInYMin[prim];
              VBndPlaneInYMax[effprim] = VBndPlaneInYMax[prim];
              VBndPlaneInZMin[effprim] = VBndPlaneInZMin[prim];
              VBndPlaneInZMax[effprim] = VBndPlaneInZMax[prim];
            }  // else remove == 0
          }    // loop over primitives to remove the primitives tagged to be
               // removed
          fclose(fprrm);

          NbPrimitives -= NbRemoved;
          printf(
              "Number of primitives removed: %d, Effective NbPrimitives: %d\n",
              NbRemoved, NbPrimitives);
          fflush(stdout);
        }  // if NbRmPrims true, implying primitives need to be removed
        fclose(rmprimFile);
      }  // if the rmprimFile is not NULL, prepare to remove primitives
    }    // if OptRmPrim: remove primitives as desired by the user

    // Information about primitives which are being ignored
    char IgnorePrimFile[256];
    strcpy(IgnorePrimFile, NativePrimDir);
    strcat(IgnorePrimFile, "/IgnorePrims.info");
    FILE *fignore = fopen(IgnorePrimFile, "w");
    if (fignore == NULL) {
      printf(
          "error opening IgnorePrims.info file in write mode ... returning\n");
      return (-1);
    }

    for (int prim = 1; prim <= OrgnlNbPrimitives; ++prim) {
      fprintf(fignore, "%d %d %d\n", OrgnlToEffPrim[prim][0],
              OrgnlToEffPrim[prim][1], OrgnlToEffPrim[prim][2]);
    }

    fclose(fignore);
  }  // Ignore unnecessary primitives from the final count

  // Reduced-Order Modelling information
  printf("ROM: switch to primitive representation after %d repetitions.\n",
         PrimAfter);

  // Store model data in native neBEM format
  char NativeInFile[256];

  strcpy(NativeInFile, NativeOutDir);
  strcat(NativeInFile, "/neBEMNative.inp");
  FILE *fNativeInFile = fopen(NativeInFile, "w");
  fprintf(fNativeInFile, "#====>Input directory\n");
  fprintf(fNativeInFile, "%s\n", NativePrimDir);
  fprintf(fNativeInFile, "#====>No. of primitives:\n");
  fprintf(fNativeInFile, "%d\n", NbPrimitives);
  fprintf(fNativeInFile, "#====>No. of volumes:\n");
  fprintf(fNativeInFile, "%d\n", VolMax);
  for (int prim = 1; prim <= NbPrimitives; ++prim) {
    char NativePrimFile[256];
    char strPrimNb[11];
    sprintf(strPrimNb, "%d", prim);
    strcpy(NativePrimFile, "Primitive");
    strcat(NativePrimFile, strPrimNb);
    strcat(NativePrimFile, ".inp");

    fprintf(fNativeInFile, "#Input filename for primitive:\n");
    fprintf(fNativeInFile, "%s\n", NativePrimFile);

    strcpy(NativePrimFile, NativePrimDir);
    strcat(NativePrimFile, "/Primitive");
    strcat(NativePrimFile, strPrimNb);
    strcat(NativePrimFile, ".inp");

    FILE *fNativePrim = fopen(NativePrimFile, "w");

    fprintf(fNativePrim, "#Number of vertices:\n");
    fprintf(fNativePrim, "%d\n", NbVertices[prim]);
    for (int vertex = 0; vertex < NbVertices[prim]; ++vertex) {
      fprintf(fNativePrim, "#Vertex %d:\n", vertex);
      fprintf(fNativePrim, "%lg %lg %lg\n", XVertex[prim][vertex],
              YVertex[prim][vertex], ZVertex[prim][vertex]);
    }  // for vertex
    fprintf(fNativePrim, "#Outward normal:\n");
    fprintf(fNativePrim, "%lg %lg %lg\n", XNorm[prim], YNorm[prim],
            ZNorm[prim]);
    fprintf(fNativePrim, "#Volume references:\n");
    fprintf(fNativePrim, "%d %d\n", VolRef1[prim], VolRef2[prim]);
    fprintf(fNativePrim, "#Maximum number of segments:\n");
    fprintf(fNativePrim, "%d %d\n", 0, 0);  // use the trio target, min, max
    fprintf(fNativePrim, "#Information on X periodicity:\n");
    fprintf(fNativePrim, "%d %d %lg\n", PeriodicTypeX[prim], PeriodicInX[prim],
            XPeriod[prim]);
    fprintf(fNativePrim, "#Information on Y periodicity:\n");
    fprintf(fNativePrim, "%d %d %lg\n", PeriodicTypeY[prim], PeriodicInY[prim],
            YPeriod[prim]);
    fprintf(fNativePrim, "#Information on Z periodicity:\n");
    fprintf(fNativePrim, "%d %d %lg\n", PeriodicTypeZ[prim], PeriodicInZ[prim],
            ZPeriod[prim]);

    fclose(fNativePrim);
  }  // for prim

  fprintf(fNativeInFile, "#====>No. of boundary conditions per element:\n");
  fprintf(fNativeInFile, "%d\n", 1);  // CHECK!!!
  fprintf(fNativeInFile, "#====>Device output directory name:\n");
  fprintf(fNativeInFile, "NativeResults\n");
  fprintf(fNativeInFile, "#====>Map input file:\n");
  fprintf(fNativeInFile, "MapFile.inp\n");
  fclose(fNativeInFile);

  for (int volume = 0; volume <= VolMax; ++volume) {
    // Note that materials from 1 to 10 are conductors and
    //                     from 11 to 20 are dielectrics
    int shape, material, boundarytype;
    double epsilon, potential, charge;
    neBEMVolumeDescription(volume, &shape, &material, &epsilon, &potential,
                           &charge, &boundarytype);

    char NativeVolFile[256];
    char strVolNb[11];
    sprintf(strVolNb, "%d", volume);

    strcpy(NativeVolFile, NativePrimDir);
    strcat(NativeVolFile, "/Volume");
    strcat(NativeVolFile, strVolNb);
    strcat(NativeVolFile, ".inp");

    FILE *fNativeVol = fopen(NativeVolFile, "w");

    fprintf(fNativeVol, "#Shape of the volume:\n");
    fprintf(fNativeVol, "%d\n", shape);
    fprintf(fNativeVol, "#Material of the volume:\n");
    fprintf(fNativeVol, "%d\n", material);
    fprintf(fNativeVol, "#Relative permittivity:\n");
    fprintf(fNativeVol, "%lg\n", epsilon);
    fprintf(fNativeVol, "#Potential:\n");
    fprintf(fNativeVol, "%lg\n", potential);
    fprintf(fNativeVol, "#Charge:\n");
    fprintf(fNativeVol, "%lg\n", charge);
    fprintf(fNativeVol, "#Boundary type:\n");
    fprintf(fNativeVol, "%d\n", boundarytype);

    fclose(fNativeVol);
  }

  // create a dummy map file
  {
    char NativeMapFile[256];

    strcpy(NativeMapFile, NativePrimDir);
    strcat(NativeMapFile, "/MapFile.inp");

    FILE *fNativeMap = fopen(NativeMapFile, "w");

    fprintf(fNativeMap, "#Number of maps\n");
    fprintf(fNativeMap, "%d\n", 0);
    fprintf(fNativeMap, "#Number of lines\n");
    fprintf(fNativeMap, "%d\n", 0);
    fprintf(fNativeMap, "#Number of points\n");
    fprintf(fNativeMap, "%d\n", 0);

    fclose(fNativeMap);
  }

  // Store primitive related data in a file for a new model, if opted for
  if (NewModel && OptStorePrimitives) {
    if (OptFormattedFile) {
      fstatus = WritePrimitives();
      if (fstatus) {
        neBEMMessage("neBEMReadGeometry - problem writing Primtives.\n");
        return -1;
      }
    }  // formatted file

    if (OptUnformattedFile) {
      neBEMMessage(
          "neBEMReadGeometry - unformatted write not inplemented yet.\n");
      return -1;
    }  // unformatted file
  }    // store primitives

  printf("neBEMReadGeometry: Geometry read!\n");

  neBEMState = 3;  // info about primitives read in after initialization and Nbs

  stopClock = clock();
  neBEMTimeElapsed(startClock, stopClock);
  printf("to read geometry\n");

  return (0);  // Success!
}  // neBEMReadGeometry ends

// Discretize the primitives using discretization information supplied by the
// user.
int neBEMDiscretize(int **NbElemsOnPrimitives) {
  // int fstatus;

  // For a model and a mesh that were defined before and for which data were
  // stored in two files - one for primitives and one for elements
  // The following operation does not assume that the primitives have been read
  // in from a stored file. In essence, these two read-in-s are maintained
  // independent of each other. However, for a stored element file to be useful,
  // both the model and mesh have to be old (not new).
  if ((!NewModel) && (!NewMesh) && (!NewBC) && (OptStoreElements)) {
    int fstatus = ReadElements();
    if (fstatus) {
      neBEMMessage("neBEMDiscretize - problem reading stored Elements.\n");
      return -1;
    }
    neBEMState = 4;
    return 0;
  }

  // Otherwise, continue with fresh discretization
  if (neBEMState != 3) {
    printf("discretization can continue only in State 3 ...\n");
    return -1;
  }

  // Only the number of primitives has been ascertained.
  // All the rest globally important numbers will be determined in this function
  NbSurfs = 0;
  NbWires = 0;
  NbElements = 0;

  // Here, the primitive can be analyzed and the elements necessary to
  // discretize it, may be determined. A user hint may help that can be supplied
  // during the setting up of the device. The hint may be as naive as gross,
  // coarse, medium, regular, fine depending on which the element size (in %
  // of the device / primitive) may be decided upon.
  // Earlier method of specifying the number of primitives on a priori has been
  // over-ridden.
  // Now we prescribe the length / area of each element and discretize each
  // primitive such that each element has an length / area close to the
  // suggested value.
  // Since the number of elements are being decided here, all the shapes will
  // have to be one among wire, right triangle and restangle. Any decomposition
  // of arbitrary polygons into rectangles (as many as possible) and right
  // triangles will have to be done earlier. All the counts have been
  // incremented by one to be on the safe side! The additional memory can be
  // minimized through a more careful computation of the required allocation for
  // each type of primitive.
  char MeshLogFile[256];

  strcpy(MeshLogFile, MeshOutDir);
  strcat(MeshLogFile, "/MeshLog.out");
  fMeshLog = fopen(MeshLogFile, "w");
  fprintf(fMeshLog, "Details of primitive discretization\n");

  for (int prim = 1; prim <= NbPrimitives; ++prim) {
    if (NbVertices[prim] == 4) {
      NbSurfSegX[prim] = NbElemsOnPrimitives[prim][1];
      NbSurfSegZ[prim] = NbElemsOnPrimitives[prim][2];
      int fstatus =
          AnalyzePrimitive(prim, &NbSurfSegX[prim], &NbSurfSegZ[prim]);
      if (fstatus == 0) {
        neBEMMessage("neBEMDiscretize - AnalyzePrimitve");
        return -1;
      }
      NbElements += (NbSurfSegX[prim] + 1) * (NbSurfSegZ[prim] + 1);
    }
    if (NbVertices[prim] == 3) {
      NbSurfSegX[prim] = NbElemsOnPrimitives[prim][1];
      NbSurfSegZ[prim] = NbElemsOnPrimitives[prim][2];
      int fstatus =
          AnalyzePrimitive(prim, &NbSurfSegX[prim], &NbSurfSegZ[prim]);
      if (fstatus == 0) {
        neBEMMessage("neBEMDiscretize - AnalyzePrimitive");
        return -1;
      }
      NbElements += (NbSurfSegX[prim] + 1) * (NbSurfSegZ[prim] + 1);
    }
    if (NbVertices[prim] == 2) {
      int itmp;
      NbWireSeg[prim] = NbElemsOnPrimitives[prim][1];
      int fstatus = AnalyzePrimitive(prim, &NbWireSeg[prim], &itmp);
      if (fstatus == 0) {
        neBEMMessage("neBEMDiscretize - AnalyzePrimitive");
        return -1;
      }
      NbElements += (NbWireSeg[prim] + 1);
    }
    if (DebugLevel == 101) {
      if (NbVertices[prim] == 2) {
        printf("Primitive %d to be discretized into %d elements.\n", prim,
               NbWireSeg[prim]);
      } else {
        printf("Primitive %d to be discretized into %d X %d elements.\n", prim,
               NbSurfSegX[prim], NbSurfSegZ[prim]);
      }
    }
  }
  printf("Memory allocated for maximum %d elements.\n", NbElements);
  fclose(fMeshLog);

  // Allocate enough space to store all the elements
  if (neBEMState == 3) {
    printf("neBEMDiscretize: NbElements = %d, sizeof(Element) = %zu\n",
           NbElements, sizeof(Element));
    if (EleArr) {
      Element *tmp = (Element *)realloc(EleArr, NbElements * sizeof(Element));
      if (tmp != NULL) {
        EleArr = tmp;
        EleCntr = 0;
      } else {
        free(EleArr);
        printf("neBEMDiscretize: Re-allocating EleArr failed.\n");
        return (1);
      }
      printf("neBEMDiscretize: Re-allocated EleArr.\n");
    }  // if EleArr => re-allocation
    else {
      EleArr = (Element *)malloc(NbElements * sizeof(Element));
      if (EleArr == NULL) {
        neBEMMessage("neBEMDiscretize - EleArr malloc");
        return -1;
      }
    }  // else EleArr => fresh allocation
  }    // neBEMState == 3

  // Prepare a data file that will contain the plotting information of the
  // the primitives and the elements
  if (OptGnuplot) {
    char GnuFile[256];
    strcpy(GnuFile, MeshOutDir);
    strcat(GnuFile, "/GViewDir/gPrimView.gp");
    fgnuPrim = fopen(GnuFile, "w");
    fprintf(fgnuPrim, "set title \"neBEM primitives in gnuplot VIEWER\"\n");
    // fprintf(fgnu, "#set label 1 \'LengthScale = %d\', LengthScale, right\n");
    fprintf(fgnuPrim, "#set pm3d\n");
    fprintf(fgnuPrim, "#set style data pm3d\n");
    fprintf(fgnuPrim, "#set palette model CMY\n");
    fprintf(fgnuPrim, "set hidden3d\n");
    fprintf(fgnuPrim, "set nokey\n");
    fprintf(fgnuPrim, "set xlabel \"X\"\n");
    fprintf(fgnuPrim, "set ylabel \"Y\"\n");
    fprintf(fgnuPrim, "set zlabel \"Z\"\n");
    fprintf(fgnuPrim, "set view 70, 335, 1, 1\n");
    fprintf(fgnuPrim, "\nsplot \\\n");

    strcpy(GnuFile, MeshOutDir);
    strcat(GnuFile, "/GViewDir/gElemView.gp");
    fgnuElem = fopen(GnuFile, "w");
    fprintf(fgnuElem, "set title \"neBEM elements in gnuplot VIEWER\"\n");
    // fprintf(fgnu, "#set label 1 \'LengthScale = %d\', LengthScale, right\n");
    fprintf(fgnuElem, "#set pm3d\n");
    fprintf(fgnuElem, "#set style data pm3d\n");
    fprintf(fgnuElem, "#set palette model CMY\n");
    fprintf(fgnuElem, "set hidden3d\n");
    fprintf(fgnuElem, "set nokey\n");
    fprintf(fgnuElem, "set xlabel \"X\"\n");
    fprintf(fgnuElem, "set ylabel \"Y\"\n");
    fprintf(fgnuElem, "set zlabel \"Z\"\n");
    fprintf(fgnuElem, "set view 70, 335, 1, 1\n");
    fprintf(fgnuElem, "\nsplot \\\n");

    strcpy(GnuFile, MeshOutDir);
    strcat(GnuFile, "/GViewDir/gMeshView.gp");
    fgnuMesh = fopen(GnuFile, "w");
    fprintf(fgnuMesh, "set title \"neBEM mesh in gnuplot VIEWER\"\n");
    // fprintf(fgnu, "#set label 1 \'LengthScale = %d\', LengthScale, right\n");
    fprintf(fgnuMesh, "#set pm3d\n");
    fprintf(fgnuMesh, "#set style data pm3d\n");
    fprintf(fgnuMesh, "#set palette model CMY\n");
    fprintf(fgnuMesh, "set hidden3d\n");
    fprintf(fgnuMesh, "set nokey\n");
    fprintf(fgnuMesh, "set xlabel \"X\"\n");
    fprintf(fgnuMesh, "set ylabel \"Y\"\n");
    fprintf(fgnuMesh, "set zlabel \"Z\"\n");
    fprintf(fgnuMesh, "set view 70, 335, 1, 1\n");
    fprintf(fgnuMesh, "\nsplot \\\n");
  }

  for (int prim = 1; prim <= NbPrimitives; ++prim) {
    switch (PrimType[prim]) {
      int fstatus;
      case 3:  // triangular surface
      case 4:  // rectangular surface
        ++NbSurfs;
        fstatus = SurfaceElements(
            prim, NbVertices[prim], XVertex[prim], YVertex[prim], ZVertex[prim],
            XNorm[prim], YNorm[prim], ZNorm[prim], VolRef1[prim], VolRef2[prim],
            InterfaceType[prim], ApplPot[prim], ApplCh[prim], Lambda[prim],
            NbSurfSegX[prim], NbSurfSegZ[prim]);
        if (fstatus != 0) {
          neBEMMessage("neBEMDiscretize - SurfaceElements");
          return -1;
        }
        break;
      case 2:       // wire - a wire presumably has only 2 vertices
        ++NbWires;  // it has one radius and one segmentation information
        fstatus = WireElements(
            prim, NbVertices[prim], XVertex[prim], YVertex[prim], ZVertex[prim],
            Radius[prim], VolRef1[prim], VolRef2[prim], InterfaceType[prim],
            ApplPot[prim], ApplCh[prim], Lambda[prim], NbWireSeg[prim]);
        if (fstatus != 0) {
          neBEMMessage("neBEMDiscretize - WireElements");
          return -1;
        }
        break;

      default:
        printf("PrimType out of range in CreateElements ... exiting ...\n");
        exit(-1);
    }  // switch PrimType ends
  }    // loop on prim number ends

  if (OptGnuplot) {
    fprintf(fgnuPrim, "\n\npause-1");
    fclose(fgnuPrim);
    fprintf(fgnuElem, "\n\npause-1");
    fclose(fgnuElem);
    fprintf(fgnuMesh, "\n\npause-1");
    fclose(fgnuMesh);
  }

  // If the required memory exceeds the maximum allowed number of elements
  if (EleCntr > NbElements) {
    neBEMMessage("neBEMDiscretize - EleCntr more than NbElements!");
    return -1;
  }

  // Check whether collocation points overlap
  {
    int startcntr = 1, cntr1, cntr2, NbCollPtOverlaps = 0;
    Point3D pt1, pt2;
    double dist;
    for (cntr1 = startcntr; cntr1 <= EleCntr; ++cntr1) {
      pt1.X = (EleArr + cntr1 - 1)->BC.CollPt.X;
      pt1.Y = (EleArr + cntr1 - 1)->BC.CollPt.Y;
      pt1.Z = (EleArr + cntr1 - 1)->BC.CollPt.Z;
      for (cntr2 = cntr1 + 1; cntr2 <= EleCntr; ++cntr2) {
        pt2.X = (EleArr + cntr2 - 1)->BC.CollPt.X;
        pt2.Y = (EleArr + cntr2 - 1)->BC.CollPt.Y;
        pt2.Z = (EleArr + cntr2 - 1)->BC.CollPt.Z;

        dist = GetDistancePoint3D(&pt1, &pt2);
        if (dist <=
            MINDIST)  // we need a linked-list here so that the overlapped
        {  // element is easily deleted and the rest upgraded immediately
          ++NbCollPtOverlaps;

          // Upgrade the element array manually, starting from cntr2 and restart
          // the overlap check. At present it is only a warning to the user with
          // some relevant information.
          // Find the primitives and volumes for the overlapping elements
          // The element structure should also maintain information on the
          // volumes that an element belongs to.
          int prim1 = (EleArr + cntr1 - 1)->PrimitiveNb;
          int volele1 = VolRef1[prim1];
          int prim2 = (EleArr + cntr2 - 1)->PrimitiveNb;
          int volele2 = VolRef1[prim2];

          neBEMMessage("neBEMDiscretize - Overlapping collocation points!");
          printf("Element %d, primitive %d, volume %d overlaps with\n", cntr1,
                 prim1, volele1);
          printf("\telement %d, primitive %d, volume %d.\n", cntr2, prim2,
                 volele2);
          printf("\tposition 1: (%g , %g , %g) micron,\n", 1e6 * pt1.X,
                 1e6 * pt1.Y, 1e6 * pt1.Z);
          printf("\tposition 2: (%g , %g , %g) micron.\n", 1e6 * pt2.X,
                 1e6 * pt2.Y, 1e6 * pt2.Z);
          printf("Please redo the geometry.\n");
          return -1;
        }  // if dist <= MINDIST
      }    // for cntr2
    }      // for cntr1
  }        // check collocation point overlap

  NbElements = EleCntr;  // the final number of elements
  printf("Total final number of elements: %d\n", NbElements);

  // Store element related data in a file for a new mesh created, if opted for
  if (NewMesh && OptStoreElements) {
    if (OptFormattedFile) {
      int fstatus = WriteElements();
      if (fstatus) {
        neBEMMessage("neBEMDiscretize - problem writing Elements.\n");
        return -1;
      }
    }  // formatted file

    if (OptUnformattedFile) {
      neBEMMessage(
          "neBEMDiscretize - unformatted write not inplemented yet.\n");
      return -1;
    }  // unformatted file
  }    // store elements

  neBEMState = 4;
  stopClock = clock();
  neBEMTimeElapsed(startClock, stopClock);
  printf("to complete discretization\n");

  return (0);
}  // neBEMDiscretize ends

int neBEMBoundaryInitialConditions(void) {
  startClock = clock();

  // The boundary conditions can be set only with neBEMState == 4 or 7
  // This state is assigned either after element discretization has been
  // completed or in a condition when we are looking for modifying only the
  // boundary condition for a device having same geometry (hence, the same
  // inverted influence coefficient matrix)
  if ((neBEMState == 4) || (neBEMState == 7)) {
    int fstatus = BoundaryConditions();
    if (fstatus != 0) {
      neBEMMessage("neBEMBondaryInitialConditions - BoundaryConditions");
      return -1;
    }
    fstatus = InitialConditions();
    if (fstatus != 0) {
      neBEMMessage("neBEMBondaryInitialConditions - InitialConditions");
      return -1;
    }
    if (neBEMState == 4) neBEMState = 5;  // create LHMatrix, invert etc
    if (neBEMState == 7) neBEMState = 8;  // assume LHS and inversion to be over
  } else {
    printf("Boundary & initial conditions can be set in states 4 / 7 ...\n");
    printf("returning ...\n");
    return (-1);
  }

  stopClock = clock();
  neBEMTimeElapsed(startClock, stopClock);
  printf("to setup boundary and initial conditions.\n");

  return (0);
}  // neBEMBoundaryInitialConditions ends

// There are several flags associated with this crucial step.
// These flags have a hieracrchy, the first mentioned having the highest
// priority.
// NewModel: 1 implies a fresh calculation.
// NewMesh: 1 implies a new mesh for a new / old device.
// NewBC: 1 implies new RHS for the same LHS; skips matrix inversion.
// NewPP: 1 implies the use of the same solution;
//             skips matrix inversion, as well as, computing the solution.
// It should be noted that a given device can be modeled in different manners.
// The same model for a device can have different discretization,
// same mesh different boundary conditions leading to different solutions
// and the same solution used to carry out different post-processes, we maintain
// four counters as well.
// ModelCntr: keeps track of the model for a given device,
// MeshCntr: keeps track of the mesh for a given model
// BCCntr: keeps track of the association of the boundary condition and
//               its solution. This has to maintained by the user manually and
//               supplied, for example, while carrying out a post-processing
//               for a solution that was computed before.
// PPCntr: numbers different post-processes for a given solution resulting from
//         a given set of preceding conditions.
int neBEMSolve(void) {
  int dbgFn = 0;

  clock_t startSolveClock = clock();

  if (TimeStep < 1) {
    neBEMMessage("neBEMSolve - TimeStep cannot be less than one!;\n");
    neBEMMessage("             Please put TimeStep = 1 for static problems.\n");
  }

  if ((neBEMState == 5) || (neBEMState == 8)) {
    if (neBEMState == 8)  // neBEMState 8 must have inverted flag on
    {  // so it must be a case looking for solution with a NewBC
      if (NewBC == 0) {
        neBEMMessage("neBEMSolve - NewBC zero for neBEMState = 8!");
        neBEMMessage("           - Nothing to be done ... returning.");
        return -1;
      }
    }

    if (NewModel) {  // effectively, NewMesh = NewBC = NewPP = 1;
      int fstatus = ComputeSolution();
      if (fstatus != 0) {
        neBEMMessage("neBEMSolve - NewModel");
        return -1;
      }
    } else {  // NewModel == 0
      if (NewMesh) {
        // effectively, NewBC = NewPP = 1;
        int fstatus = ComputeSolution();
        if (fstatus != 0) {
          neBEMMessage("neBEMSolve - NewMesh");
          return -1;
        }
      } else {        // NewModel == NewMesh == 0
        if (NewBC) {  // effectively, NewPP = 1;
          int fstatus = ComputeSolution();
          if (fstatus != 0) {
            neBEMMessage("neBEMSolve - Failure computing new solution");
            return -1;
          }
        } else {  // NewBC == 0
          if (NewPP) {
            int fstatus = ReadSolution();
            if (fstatus != 0) {
              neBEMMessage("neBEMSolve - Failure reading solution");
              return (-1);
            }
          } else {  // NewPP == 0
            printf("neBEMSolve: Nothing to do ... returning ...\n");
            return (-1);
          }  // NewPP == 0
        }    // NewBC == 0
      }      // NewModel == NewDiscretization == 0
    }        // NewModel == 0

    neBEMState = 9;
  } else {
    printf("neBEMSolve: neBEMSolve can be called only in state 5 / 8 ...\n");
    printf("returning ...\n");
    return (-1);
  }

  if (FailureCntr) {
    printf(
        "neBEMSolve: Approximations were made while computing the influence "
        "coefficients.\n");
    printf("            Please check the \"%s/Isles.log\" file.\n", PPOutDir);
  }

  clock_t stopSolveClock = clock();
  neBEMTimeElapsed(startSolveClock, stopSolveClock);
  printf("to complete solve\n");

  // Prepare voxelized data that will be exported to Garfield++
  if (OptVoxel) {
    clock_t startVoxelClock = clock();

    int fstatus = VoxelFPR();
    if (fstatus != 0) {
      neBEMMessage("neBEMSolve - Failure computing VoxelFPR");
      return -1;
    }

    clock_t stopVoxelClock = clock();
    neBEMTimeElapsed(startVoxelClock, stopVoxelClock);
    printf("to compute VoxelFPR\n");
  }

  // Prepare 3dMap data that will be exported to Garfield++
  if (OptMap) {
    clock_t startMapClock = clock();

    int fstatus = MapFPR();
    if (fstatus != 0) {
      neBEMMessage("neBEMSolve - Failure computing MapFPR");
      return -1;
    }

    clock_t stopMapClock = clock();
    neBEMTimeElapsed(startMapClock, stopMapClock);
    printf("to compute MapFPR\n");
  }

  // allocate memory for potential and field components within FAST volume mesh
  // and compute / read FastVol data
  // Similar allocation, computation and reading may be necessary for the KnCh
  // effects.
  // The other approach could be to create fast volumes that always have both
  // the influences (elements + known charges) added together. This approach
  // seems more managable now and is being followed.
  // Thus, if we want to inspect the effects of elements and known charges
  // separately, we will have to generate one fast volume with OptKnCh = 0,
  // and another with OptKnCh = 1. Subtraction of these two fast volumes will
  // provide us with the effect of KnCh.
  if (OptFastVol) {
    int MaxXCells = BlkNbXCells[1];
    int MaxYCells = BlkNbYCells[1];
    int MaxZCells = BlkNbZCells[1];
    clock_t startFastClock = clock();

    // find maximum number of Xcells etc in all the blocks
    // simplifies memory allocation using nrutils but hogs memory!
    for (int block = 1; block <= FastVol.NbBlocks; ++block) {
      if (block == 1) {
        MaxXCells = BlkNbXCells[1];
        MaxYCells = BlkNbYCells[1];
        MaxZCells = BlkNbZCells[1];
      } else {
        if (MaxXCells < BlkNbXCells[block]) MaxXCells = BlkNbXCells[block];
        if (MaxYCells < BlkNbYCells[block]) MaxYCells = BlkNbYCells[block];
        if (MaxZCells < BlkNbZCells[block]) MaxZCells = BlkNbZCells[block];
      }
    }  // loop block for finding maxm cells among all the blocks

    if (dbgFn) {
      printf("OptFastVol: %d\n", OptFastVol);
      printf("NbPtSkip: %d\n", NbPtSkip);
      printf("OptStaggerFastVol: %d\n", OptStaggerFastVol);
      printf("NbStgPtSkip: %d\n", NbStgPtSkip);
      printf("OptReadFastPF: %d\n", OptReadFastPF);
      printf("LX: %le\n", FastVol.LX);
      printf("LY: %le\n", FastVol.LY);
      printf("LZ: %le\n", FastVol.LZ);
      printf("CornerX: %le\n", FastVol.CrnrX);
      printf("CornerY: %le\n", FastVol.CrnrY);
      printf("CornerZ: %le\n", FastVol.CrnrZ);
      printf("YStagger: %le\n", FastVol.YStagger);
      printf("NbOfBlocks: %d\n", FastVol.NbBlocks);
      for (int block = 1; block <= FastVol.NbBlocks; ++block) {
        printf("NbOfXCells[%d]: %d\n", block, BlkNbXCells[block]);
        printf("NbOfYCells[%d]: %d\n", block, BlkNbYCells[block]);
        printf("NbOfZCells[%d]: %d\n", block, BlkNbZCells[block]);
        printf("LZ[%d]: %le\n", block, BlkLZ[block]);
        printf("CornerZ[%d]: %le\n", block, BlkCrnrZ[block]);
      }
      printf("NbOfOmitVols: %d\n", FastVol.NbOmitVols);
      if (FastVol.NbOmitVols) {
        for (int omit = 1; omit <= FastVol.NbOmitVols; ++omit) {
          printf("OmitVolLX[%d]: %le\n", omit, OmitVolLX[omit]);
          printf("OmitVolLY[%d]: %le\n", omit, OmitVolLY[omit]);
          printf("OmitVolLZ[%d]: %le\n", omit, OmitVolLZ[omit]);
          printf("OmitCrnrX[%d]: %le\n", omit, OmitVolCrnrX[omit]);
          printf("OmitCrnrY[%d]: %le\n", omit, OmitVolCrnrY[omit]);
          printf("OmitCrnrZ[%d]: %le\n", omit, OmitVolCrnrZ[omit]);
        }
      }
      printf("MaxXCells, MaxYCells, MaxZCells: %d, %d, %d\n", MaxXCells,
             MaxYCells, MaxZCells);
    }  // dbgFn

    // Following memory allocations are necessary even for creating
    // the fast volumes.
    // In fact, since these tensors are already in place, the fast volumes
    // can be used to evaluate potential and fields immediately after they are
    // created.
    /* Memory wastage!!! Improve as soon as possible. */
    FastPot = d4tensor(1, FastVol.NbBlocks, 1, MaxXCells + 1, 1, MaxYCells + 1,
                       1, MaxZCells + 1);
    FastFX = d4tensor(1, FastVol.NbBlocks, 1, MaxXCells + 1, 1, MaxYCells + 1,
                      1, MaxZCells + 1);
    FastFY = d4tensor(1, FastVol.NbBlocks, 1, MaxXCells + 1, 1, MaxYCells + 1,
                      1, MaxZCells + 1);
    FastFZ = d4tensor(1, FastVol.NbBlocks, 1, MaxXCells + 1, 1, MaxYCells + 1,
                      1, MaxZCells + 1);

    if (OptStaggerFastVol) {
      /* Memory wastage!!! Improve as soon as possible. */
      StgFastPot = d4tensor(1, FastVol.NbBlocks, 1, MaxXCells + 1, 1,
                            MaxYCells + 1, 1, MaxZCells + 1);
      StgFastFX = d4tensor(1, FastVol.NbBlocks, 1, MaxXCells + 1, 1,
                           MaxYCells + 1, 1, MaxZCells + 1);
      StgFastFY = d4tensor(1, FastVol.NbBlocks, 1, MaxXCells + 1, 1,
                           MaxYCells + 1, 1, MaxZCells + 1);
      StgFastFZ = d4tensor(1, FastVol.NbBlocks, 1, MaxXCells + 1, 1,
                           MaxYCells + 1, 1, MaxZCells + 1);
    }  // if OptStaggerFastVol

    if (OptCreateFastPF) {
      int fstatus = CreateFastVolPF();
      if (fstatus != 0) {
        neBEMMessage("neBEMSolve - Failure computing FastVolPF");
        return -1;
      }

      clock_t stopFastClock = clock();
      neBEMTimeElapsed(startFastClock, stopFastClock);
      printf("to compute FastVolPF\n");
    }  // if OptCreateFastPF

    if (OptReadFastPF) {  // reading from file
      int nbXCells, nbYCells, nbZCells;
      int tmpblk;
      double xpt, ypt, zpt;

      char FastVolPFFile[256];
      strcpy(FastVolPFFile, BCOutDir);
      strcat(FastVolPFFile, "/FastVolPF.out");
      FILE *fFastVolPF = fopen(FastVolPFFile, "r");
      if (fFastVolPF == NULL) {
        neBEMMessage("in neBEMSolve - FastVolPFFile");
        return -1;
      }

      fscanf(fFastVolPF, "#block\tX\tY\tZ\tPot\tFX\tFY\tFZ\n");

      for (int block = 1; block <= FastVol.NbBlocks; ++block) {
        nbXCells = BlkNbXCells[block];
        nbYCells = BlkNbYCells[block];
        nbZCells = BlkNbZCells[block];

        for (int i = 1; i <= nbXCells + 1; ++i) {
          for (int j = 1; j <= nbYCells + 1; ++j) {
            for (int k = 1; k <= nbZCells + 1; ++k) {
              fscanf(fFastVolPF, "%d\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n",
                     &tmpblk, &xpt, &ypt, &zpt, &FastPot[block][i][j][k],
                     &FastFX[block][i][j][k], &FastFY[block][i][j][k],
                     &FastFZ[block][i][j][k]);
            }  // loop k
          }    // loop j
        }      // loop i
      }        // loop block
      fclose(fFastVolPF);

      if (OptStaggerFastVol) {
        char StgFastVolPFFile[256];
        strcpy(StgFastVolPFFile, BCOutDir);
        strcat(StgFastVolPFFile, "/StgFastVolPF.out");
        FILE *fStgFastVolPF = fopen(StgFastVolPFFile, "r");
        if (fStgFastVolPF == NULL) {
          neBEMMessage("in neBEMSolve - StgFastVolPFFile");
          return -1;
        }

        fscanf(fStgFastVolPF, "#block\tX\tY\tZ\tPot\tFX\tFY\tFZ\n");

        for (int block = 1; block <= FastVol.NbBlocks; ++block) {
          nbXCells = BlkNbXCells[block];
          nbYCells = BlkNbYCells[block];
          nbZCells = BlkNbZCells[block];

          for (int i = 1; i <= nbXCells + 1; ++i) {
            for (int j = 1; j <= nbYCells + 1; ++j) {
              for (int k = 1; k <= nbZCells + 1; ++k) {
                fscanf(fStgFastVolPF, "%d\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n",
                       &tmpblk, &xpt, &ypt, &zpt, &StgFastPot[block][i][j][k],
                       &StgFastFX[block][i][j][k], &StgFastFY[block][i][j][k],
                       &StgFastFZ[block][i][j][k]);
              }  // loop k
            }    // loop j
          }      // loop i
        }        // loop block
        fclose(fStgFastVolPF);
      }  // if OptStaggerFastVol

      clock_t stopFastClock = clock();
      neBEMTimeElapsed(startFastClock, stopFastClock);
      printf("to read FastVolPF\n");
    }  // if OptReadFastPF
  }    // if OptFastVol

  return (0);
}  // neBEMSolve ends

// Get potential and field at a given point
int neBEMPF(Point3D *point, double *potential, Vector3D *field) {
  if (neBEMState < 9) {
    printf("neBEMPF cannot be called before reaching state 9.\n");
    return (-1);
  }

  // printf("neBEMPF called %8d times", ++neBEMPFCallCntr);

  double Pot;
  int fstatus;
  if (OptFastVol)  // Note: this is not the Create or Read option
  {
    fstatus = FastPFAtPoint(point, &Pot, field);
    if (fstatus != 0) {
      neBEMMessage("neBEMPF - FastPFAtPoint");
      return -1;
    }
  } else {
    fstatus = PFAtPoint(point, &Pot, field);
    if (fstatus != 0) {
      neBEMMessage("neBEMPF - PFAtPoint");
      return -1;
    }
  }

  *potential = Pot;

  // printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");

  return (0);
}  // neBEMPF ends

// Actual preparation of the weighting field, including those related to
// corresponding weighting field fast volumes.
// The return value identifies the weighting field. Error: id < 0.
// The list contains all primitives that are part of this particular
// read-out group. These can come from several volumes, but it is
// not anticipated that only some of the primitives of one volume
// are listed.
// Weighitng field boundary conditions are useful only when the inverted
// matrix is available.
// This state is assigned either after element discretization has been
// completed or in a condition when we are looking for modifying only the
// boundary condition for a device having same geometry (hence, the same
// inverted influence coefficient matrix)
int neBEMPrepareWeightingField(int nprim, int primlist[]) {
  static int IdWtField = 0;

  int dbgFn = 0;
  int fstatus = 0;

  if (neBEMState < 7) {
    printf(
        "neBEMPrepareWeightingField: Weighting computations only meaningful "
        "beyond neBEMState 7 ...\n");
    return -1;
  }

  // Find first free slot
  const int MaxWtField =
      MAXWtFld;  // used also while deallocating these memories
  if (WtFieldChDen == NULL)
    WtFieldChDen = (double **)malloc(MaxWtField * sizeof(double *));
  if (AvWtChDen == NULL)
    AvWtChDen = (double **)malloc(MaxWtField * sizeof(double *));

  ++IdWtField;
  if (IdWtField >= MaxWtField) {
    printf(
        "neBEMPrepareWeightingField: reached MaxWtField (%d) weighting "
        "fields.\n",
        MAXWtFld);
    return -1;
  } else {
    printf("\nPreparing weighting field for %d-th set.\n", IdWtField);
  }  // else within MaxWtField

  // Allocate a new column to store this solution set
  WtFieldChDen[IdWtField] = (double *)malloc((NbElements + 2) * sizeof(double));
  AvWtChDen[IdWtField] = (double *)malloc((NbPrimitives + 2) * sizeof(double));

  fstatus = WeightingFieldSolution(nprim, primlist, WtFieldChDen[IdWtField]);
  if (fstatus) {
    neBEMMessage("neBEMPrepareWeightingField - WeightingFieldSolution");
    return -1;
  } else {
    printf("Computed weighting field solution\n");
  }

  if (!fstatus) {  // estimate primitive related avrg wt field charge densities
    // OMPCheck - may be parallelized
    for (int prim = 1; prim <= NbPrimitives; ++prim) {
      double area = 0.0;  // need area of the primitive as well!
      AvWtChDen[IdWtField][prim] = 0.0;

      for (int ele = ElementBgn[prim]; ele <= ElementEnd[prim]; ++ele) {
        area += (EleArr + ele - 1)->G.dA;
        AvWtChDen[IdWtField][prim] +=
            WtFieldChDen[IdWtField][ele] * (EleArr + ele - 1)->G.dA;
      }

      AvWtChDen[IdWtField][prim] /= area;
    }
    printf("Computed primitive-averaged weighting field solutions\n");
  }  // if status flag not raised

  // stringify the integer
  char strIdWtField[5];
  sprintf(strIdWtField, "%d", IdWtField);
  // printf("strIdWtField: %s\n", strIdWtField);

  // Set up parameters related to fixed specification of weighting field
  // printf("OptFixedWtField: %d\n", OptFixedWtField[IdWtField]);fflush(stdout);
  if (OptFixedWtField[IdWtField]) {
    char fileWtField[256];
    strcpy(fileWtField, "neBEMInp/neBEMFixedWtField_");
    strcat(fileWtField, strIdWtField);
    strcat(fileWtField, ".inp");
    FILE *fixedWtInpFile = fopen(fileWtField, "r");
    if (fixedWtInpFile == NULL) {
      printf(
          "neBEMFixedWtField.inp absent ... assuming OptFixedWtField = 0 "
          "...\n");
      OptFixedWtField[IdWtField] = 0;
      FixedWtPotential[IdWtField] = 0.0;
      FixedWtFieldX[IdWtField] = 0.0;
      FixedWtFieldY[IdWtField] = 0.0;
      FixedWtFieldZ[IdWtField] = 0.0;
    } else {
      fscanf(fixedWtInpFile, "OptFixedWtField: %d\n",
             &OptFixedWtField[IdWtField]);
      fscanf(fixedWtInpFile, "FixedWtPotential: %lg\n",
             &FixedWtPotential[IdWtField]);
      fscanf(fixedWtInpFile, "FixedWtFieldX: %lg\n", &FixedWtFieldX[IdWtField]);
      fscanf(fixedWtInpFile, "FixedWtFieldY: %lg\n", &FixedWtFieldY[IdWtField]);
      fscanf(fixedWtInpFile, "FixedWtFieldZ: %lg\n", &FixedWtFieldZ[IdWtField]);
      fclose(fixedWtInpFile);
    }  // else fixedWtInpFile
  }    // if OptFixedWtField

  // Weighting field fast volume related computations
  // Set up parameters related to weighting field fast volumes.
  // There can be multiple weighting fields depending on readout electrodes.
  // As a result, there can be multiple versions of this input file, each
  // reflecting the requirement of one set of electrodes.
  // printf("OptWtFldFastVol: %d\n", OptWtFldFastVol[IdWtField]);fflush(stdout);
  if (OptWtFldFastVol[IdWtField]) {
    char fileWtField[256];
    strcpy(fileWtField, "neBEMInp/neBEMWtFldFastVol_");
    strcat(fileWtField, strIdWtField);
    strcat(fileWtField, ".inp");
    FILE *fastWtFldInpFile = fopen(fileWtField, "r");
    // FILE *fastWtFldInpFile = fopen("neBEMInp/neBEMWtFldFastVol.inp", "r");
    if (fastWtFldInpFile == NULL) {
      printf(
          "neBEMWtFldFastVol.inp absent ... assuming OptWtFldFastVol = 0 "
          "...\n");
      OptWtFldFastVol[IdWtField] = 0;
      OptStaggerWtFldFastVol[IdWtField] = 0;
      OptCreateWtFldFastPF[IdWtField] = 0;
      OptReadWtFldFastPF[IdWtField] = 0;
      WtFldFastVol[IdWtField].NbBlocks = 0;
      WtFldFastVol[IdWtField].NbOmitVols = 0;
      WtFldFastVol[IdWtField].NbIgnoreVols = 0;
    } else {
      fscanf(fastWtFldInpFile, "OptFastVol: %d\n", &OptWtFldFastVol[IdWtField]);
      fscanf(fastWtFldInpFile, "OptStaggerFastVol: %d\n",
             &OptStaggerWtFldFastVol[IdWtField]);
      fscanf(fastWtFldInpFile, "OptCreateFastPF: %d\n",
             &OptCreateWtFldFastPF[IdWtField]);
      fscanf(fastWtFldInpFile, "OptReadFastPF: %d\n",
             &OptReadWtFldFastPF[IdWtField]);
      fscanf(fastWtFldInpFile, "NbPtSkip: %d\n", &WtFldNbPtSkip[IdWtField]);
      fscanf(fastWtFldInpFile, "NbStgPtSkip: %d\n",
             &StgWtFldNbPtSkip[IdWtField]);
      fscanf(fastWtFldInpFile, "LX: %le\n", &WtFldFastVol[IdWtField].LX);
      fscanf(fastWtFldInpFile, "LY: %le\n", &WtFldFastVol[IdWtField].LY);
      fscanf(fastWtFldInpFile, "LZ: %le\n", &WtFldFastVol[IdWtField].LZ);
      fscanf(fastWtFldInpFile, "CornerX: %le\n",
             &WtFldFastVol[IdWtField].CrnrX);
      fscanf(fastWtFldInpFile, "CornerY: %le\n",
             &WtFldFastVol[IdWtField].CrnrY);
      fscanf(fastWtFldInpFile, "CornerZ: %le\n",
             &WtFldFastVol[IdWtField].CrnrZ);
      fscanf(fastWtFldInpFile, "YStagger: %le\n",
             &WtFldFastVol[IdWtField].YStagger);
      if (!OptStaggerWtFldFastVol[IdWtField])
        WtFldFastVol[IdWtField].YStagger = 0.0;  // ignore any non-zero value
      fscanf(fastWtFldInpFile, "NbOfBlocks: %d\n",
             &WtFldFastVol[IdWtField].NbBlocks);
      WtFldBlkNbXCells[IdWtField] =
          ivector(1, WtFldFastVol[IdWtField].NbBlocks);
      WtFldBlkNbYCells[IdWtField] =
          ivector(1, WtFldFastVol[IdWtField].NbBlocks);
      WtFldBlkNbZCells[IdWtField] =
          ivector(1, WtFldFastVol[IdWtField].NbBlocks);
      WtFldBlkLZ[IdWtField] = dvector(1, WtFldFastVol[IdWtField].NbBlocks);
      WtFldBlkCrnrZ[IdWtField] = dvector(1, WtFldFastVol[IdWtField].NbBlocks);
      for (int block = 1; block <= WtFldFastVol[IdWtField].NbBlocks; ++block) {
        fscanf(fastWtFldInpFile, "NbOfXCells: %d\n",
               &WtFldBlkNbXCells[IdWtField][block]);
        fscanf(fastWtFldInpFile, "NbOfYCells: %d\n",
               &WtFldBlkNbYCells[IdWtField][block]);
        fscanf(fastWtFldInpFile, "NbOfZCells: %d\n",
               &WtFldBlkNbZCells[IdWtField][block]);
        fscanf(fastWtFldInpFile, "LZ: %le\n", &WtFldBlkLZ[IdWtField][block]);
        fscanf(fastWtFldInpFile, "CornerZ: %le\n",
               &WtFldBlkCrnrZ[IdWtField][block]);
      }  // inputs for blocks
      fscanf(fastWtFldInpFile, "NbOfOmitVols: %d\n",
             &WtFldFastVol[IdWtField].NbOmitVols);
      if (WtFldFastVol[IdWtField].NbOmitVols) {
        WtFldOmitVolLX[IdWtField] =
            dvector(1, WtFldFastVol[IdWtField].NbOmitVols);
        WtFldOmitVolLY[IdWtField] =
            dvector(1, WtFldFastVol[IdWtField].NbOmitVols);
        WtFldOmitVolLZ[IdWtField] =
            dvector(1, WtFldFastVol[IdWtField].NbOmitVols);
        WtFldOmitVolCrnrX[IdWtField] =
            dvector(1, WtFldFastVol[IdWtField].NbOmitVols);
        WtFldOmitVolCrnrY[IdWtField] =
            dvector(1, WtFldFastVol[IdWtField].NbOmitVols);
        WtFldOmitVolCrnrZ[IdWtField] =
            dvector(1, WtFldFastVol[IdWtField].NbOmitVols);
        for (int omit = 1; omit <= WtFldFastVol[IdWtField].NbOmitVols; ++omit) {
          fscanf(fastWtFldInpFile, "OmitVolLX: %le\n",
                 &WtFldOmitVolLX[IdWtField][omit]);
          fscanf(fastWtFldInpFile, "OmitVolLY: %le\n",
                 &WtFldOmitVolLY[IdWtField][omit]);
          fscanf(fastWtFldInpFile, "OmitVolLZ: %le\n",
                 &WtFldOmitVolLZ[IdWtField][omit]);
          fscanf(fastWtFldInpFile, "OmitVolCornerX: %le\n",
                 &WtFldOmitVolCrnrX[IdWtField][omit]);
          fscanf(fastWtFldInpFile, "OmitVolCornerY: %le\n",
                 &WtFldOmitVolCrnrY[IdWtField][omit]);
          fscanf(fastWtFldInpFile, "OmitVolCornerZ: %le\n",
                 &WtFldOmitVolCrnrZ[IdWtField][omit]);
        }  // for loop inputs for OmitVols
      }    // inputs for OmitVols
      fscanf(fastWtFldInpFile, "NbOfIgnoreVols: %d\n",
             &WtFldFastVol[IdWtField].NbIgnoreVols);
      if (WtFldFastVol[IdWtField].NbIgnoreVols) {
        WtFldIgnoreVolLX[IdWtField] =
            dvector(1, WtFldFastVol[IdWtField].NbIgnoreVols);
        WtFldIgnoreVolLY[IdWtField] =
            dvector(1, WtFldFastVol[IdWtField].NbIgnoreVols);
        WtFldIgnoreVolLZ[IdWtField] =
            dvector(1, WtFldFastVol[IdWtField].NbIgnoreVols);
        WtFldIgnoreVolCrnrX[IdWtField] =
            dvector(1, WtFldFastVol[IdWtField].NbIgnoreVols);
        WtFldIgnoreVolCrnrY[IdWtField] =
            dvector(1, WtFldFastVol[IdWtField].NbIgnoreVols);
        WtFldIgnoreVolCrnrZ[IdWtField] =
            dvector(1, WtFldFastVol[IdWtField].NbIgnoreVols);
        for (int ignore = 1; ignore <= WtFldFastVol[IdWtField].NbIgnoreVols;
             ++ignore) {
          fscanf(fastWtFldInpFile, "IgnoreVolLX: %le\n",
                 &WtFldIgnoreVolLX[IdWtField][ignore]);
          fscanf(fastWtFldInpFile, "IgnoreVolLY: %le\n",
                 &WtFldIgnoreVolLY[IdWtField][ignore]);
          fscanf(fastWtFldInpFile, "IgnoreVolLZ: %le\n",
                 &WtFldIgnoreVolLZ[IdWtField][ignore]);
          fscanf(fastWtFldInpFile, "IgnoreVolCornerX: %le\n",
                 &WtFldIgnoreVolCrnrX[IdWtField][ignore]);
          fscanf(fastWtFldInpFile, "IgnoreVolCornerY: %le\n",
                 &WtFldIgnoreVolCrnrY[IdWtField][ignore]);
          fscanf(fastWtFldInpFile, "IgnoreVolCornerZ: %le\n",
                 &WtFldIgnoreVolCrnrZ[IdWtField][ignore]);
        }  // for loop inputs for IgnoreVols
      }    // inputs for IgnoreVols
      if (dbgFn) {
        for (int ignore = 1; ignore <= WtFldFastVol[IdWtField].NbIgnoreVols;
             ++ignore) {
          printf("WtFldIgnoreVolLX: %le\n",
                 WtFldIgnoreVolLX[IdWtField][ignore]);
          printf("WtFldIgnoreVolLY: %le\n",
                 WtFldIgnoreVolLY[IdWtField][ignore]);
          printf("WtFldIgnoreVolLZ: %le\n",
                 WtFldIgnoreVolLZ[IdWtField][ignore]);
          printf("WtFldIgnoreVolCornerX: %le\n",
                 WtFldIgnoreVolCrnrX[IdWtField][ignore]);
          printf("WtFldIgnoreVolCornerY: %le\n",
                 WtFldIgnoreVolCrnrY[IdWtField][ignore]);
          printf("WtFldIgnoreVolCornerZ: %le\n",
                 WtFldIgnoreVolCrnrZ[IdWtField][ignore]);
        }  // inputs for IgnoreVols
      }
      fclose(fastWtFldInpFile);
    }  // else fastWtFldInpFile

    int MaxXCells = WtFldBlkNbXCells[IdWtField][1];
    int MaxYCells = WtFldBlkNbYCells[IdWtField][1];
    int MaxZCells = WtFldBlkNbZCells[IdWtField][1];
    clock_t startFastClock = clock();

    // find maximum number of Xcells etc in all the blocks
    // simplifies memory allocation using nrutils but hogs memory!
    for (int block = 1; block <= WtFldFastVol[IdWtField].NbBlocks; ++block) {
      if (block == 1) {
        MaxXCells = WtFldBlkNbXCells[IdWtField][1];
        MaxYCells = WtFldBlkNbYCells[IdWtField][1];
        MaxZCells = WtFldBlkNbZCells[IdWtField][1];
      } else {
        if (MaxXCells < WtFldBlkNbXCells[IdWtField][block])
          MaxXCells = WtFldBlkNbXCells[IdWtField][block];
        if (MaxYCells < WtFldBlkNbYCells[IdWtField][block])
          MaxYCells = WtFldBlkNbYCells[IdWtField][block];
        if (MaxZCells < WtFldBlkNbZCells[IdWtField][block])
          MaxZCells = WtFldBlkNbZCells[IdWtField][block];
      }
    }  // loop block for finding maxm cells among all the blocks

    if (dbgFn) {
      printf("OptWtFldFastVol: %d\n", OptWtFldFastVol[IdWtField]);
      printf("OptStaggerWtFldFastVol: %d\n", OptStaggerWtFldFastVol[IdWtField]);
      printf("OptCreateWtFldFastPF: %d\n", OptCreateWtFldFastPF[IdWtField]);
      printf("OptReadWtFldFastPF: %d\n", OptReadWtFldFastPF[IdWtField]);
      printf("WtFldNbPtSkip: %d\n", WtFldNbPtSkip[IdWtField]);
      printf("StgWtFldNbPtSkip: %d\n", StgWtFldNbPtSkip[IdWtField]);
      printf("LX: %le\n", WtFldFastVol[IdWtField].LX);
      printf("LY: %le\n", WtFldFastVol[IdWtField].LY);
      printf("LZ: %le\n", WtFldFastVol[IdWtField].LZ);
      printf("CornerX: %le\n", WtFldFastVol[IdWtField].CrnrX);
      printf("CornerY: %le\n", WtFldFastVol[IdWtField].CrnrY);
      printf("CornerZ: %le\n", WtFldFastVol[IdWtField].CrnrZ);
      printf("YStagger: %le\n", WtFldFastVol[IdWtField].YStagger);
      printf("NbOfBlocks: %d\n", WtFldFastVol[IdWtField].NbBlocks);
      for (int block = 1; block <= WtFldFastVol[IdWtField].NbBlocks; ++block) {
        printf("NbOfXCells[%d]: %d\n", block,
               WtFldBlkNbXCells[IdWtField][block]);
        printf("NbOfYCells[%d]: %d\n", block,
               WtFldBlkNbYCells[IdWtField][block]);
        printf("NbOfZCells[%d]: %d\n", block,
               WtFldBlkNbZCells[IdWtField][block]);
        printf("LZ[%d]: %le\n", block, WtFldBlkLZ[IdWtField][block]);
        printf("CornerZ[%d]: %le\n", block, WtFldBlkCrnrZ[IdWtField][block]);
      }
      printf("NbOfOmitVols: %d\n", WtFldFastVol[IdWtField].NbOmitVols);
      if (WtFldFastVol[IdWtField].NbOmitVols) {
        for (int omit = 1; omit <= WtFldFastVol[IdWtField].NbOmitVols; ++omit) {
          printf("OmitVolLX[%d]: %le\n", omit, WtFldOmitVolLX[IdWtField][omit]);
          printf("OmitVolLY[%d]: %le\n", omit, WtFldOmitVolLY[IdWtField][omit]);
          printf("OmitVolLZ[%d]: %le\n", omit, WtFldOmitVolLZ[IdWtField][omit]);
          printf("OmitCrnrX[%d]: %le\n", omit,
                 WtFldOmitVolCrnrX[IdWtField][omit]);
          printf("OmitCrnrY[%d]: %le\n", omit,
                 WtFldOmitVolCrnrY[IdWtField][omit]);
          printf("OmitCrnrZ[%d]: %le\n", omit,
                 WtFldOmitVolCrnrZ[IdWtField][omit]);
        }
      }
      printf("MaxXCells, MaxYCells, MaxZCells: %d, %d, %d\n", MaxXCells,
             MaxYCells, MaxZCells);
      printf("NbOfIgnoreVols: %d\n", WtFldFastVol[IdWtField].NbIgnoreVols);
    }  // dbgFn

    // Following memory allocations are necessary even for creating
    // the fast volumes.
    // In fact, since these tensors are already in place, the fast volumes
    // can be used to evaluate potential and fields immediately after they are
    // created.
    /* Memory wastage!!! Improve as soon as possible. */
    WtFldFastPot[IdWtField] =
        d4tensor(1, WtFldFastVol[IdWtField].NbBlocks, 1, MaxXCells + 1, 1,
                 MaxYCells + 1, 1, MaxZCells + 1);
    WtFldFastFX[IdWtField] =
        d4tensor(1, WtFldFastVol[IdWtField].NbBlocks, 1, MaxXCells + 1, 1,
                 MaxYCells + 1, 1, MaxZCells + 1);
    WtFldFastFY[IdWtField] =
        d4tensor(1, WtFldFastVol[IdWtField].NbBlocks, 1, MaxXCells + 1, 1,
                 MaxYCells + 1, 1, MaxZCells + 1);
    WtFldFastFZ[IdWtField] =
        d4tensor(1, WtFldFastVol[IdWtField].NbBlocks, 1, MaxXCells + 1, 1,
                 MaxYCells + 1, 1, MaxZCells + 1);

    if (OptStaggerWtFldFastVol[IdWtField]) {
      /* Memory wastage!!! Improve as soon as possible. */
      StgWtFldFastPot[IdWtField] =
          d4tensor(1, WtFldFastVol[IdWtField].NbBlocks, 1, MaxXCells + 1, 1,
                   MaxYCells + 1, 1, MaxZCells + 1);
      StgWtFldFastFX[IdWtField] =
          d4tensor(1, WtFldFastVol[IdWtField].NbBlocks, 1, MaxXCells + 1, 1,
                   MaxYCells + 1, 1, MaxZCells + 1);
      StgWtFldFastFY[IdWtField] =
          d4tensor(1, WtFldFastVol[IdWtField].NbBlocks, 1, MaxXCells + 1, 1,
                   MaxYCells + 1, 1, MaxZCells + 1);
      StgWtFldFastFZ[IdWtField] =
          d4tensor(1, WtFldFastVol[IdWtField].NbBlocks, 1, MaxXCells + 1, 1,
                   MaxYCells + 1, 1, MaxZCells + 1);
    }  // if OptStaggerWtFldFastVol

    if (OptCreateWtFldFastPF[IdWtField]) {
      fstatus = CreateWtFldFastVolPF(IdWtField);

      clock_t stopFastClock = clock();
      neBEMTimeElapsed(startFastClock, stopFastClock);
      printf("to compute WtFldFastVolPF\n");
      // neBEMMessage(
      // "neBEMSolve - Failure computing WtFldFastVolPF: not implemented");
      // return -1;
    }  // if OptCreateWtFldFastPF

    if (OptReadWtFldFastPF[IdWtField]) {  // reading from file
      int nbXCells, nbYCells, nbZCells;
      int tmpblk;
      double xpt, ypt, zpt;

      // stringify the integer
      char stringIdWtField[5];
      sprintf(stringIdWtField, "%d", IdWtField);
      char FastVolPFFile[256];
      strcpy(FastVolPFFile, BCOutDir);
      strcat(FastVolPFFile, "/WtFldFastVolPF_");
      strcat(FastVolPFFile, stringIdWtField);
      strcat(FastVolPFFile, ".out");
      FILE *fFastVolPF = fopen(FastVolPFFile, "r");
      if (fFastVolPF == NULL) {
        neBEMMessage("in neBEMSolve - WtFldFastVolPFFile");
        return -1;
      }

      fscanf(fFastVolPF, "#block\tX\tY\tZ\tPot\tFX\tFY\tFZ\n");

      for (int block = 1; block <= WtFldFastVol[IdWtField].NbBlocks; ++block) {
        nbXCells = WtFldBlkNbXCells[IdWtField][block];
        nbYCells = WtFldBlkNbYCells[IdWtField][block];
        nbZCells = WtFldBlkNbZCells[IdWtField][block];

        for (int i = 1; i <= nbXCells + 1; ++i) {
          for (int j = 1; j <= nbYCells + 1; ++j) {
            for (int k = 1; k <= nbZCells + 1; ++k) {
              fscanf(fFastVolPF, "%d\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n",
                     &tmpblk, &xpt, &ypt, &zpt,
                     &WtFldFastPot[IdWtField][block][i][j][k],
                     &WtFldFastFX[IdWtField][block][i][j][k],
                     &WtFldFastFY[IdWtField][block][i][j][k],
                     &WtFldFastFZ[IdWtField][block][i][j][k]);
            }  // loop k
          }    // loop j
        }      // loop i
      }        // loop block
      fclose(fFastVolPF);

      if (OptStaggerWtFldFastVol[IdWtField]) {
        // stringify the integer
        sprintf(stringIdWtField, "%d", IdWtField);
        char StgFastVolPFFile[256];
        FILE *fStgFastVolPF;
        strcpy(StgFastVolPFFile, BCOutDir);
        // strcat(StgFastVolPFFile, "/WtFldStgFastVolPF.out");
        strcat(StgFastVolPFFile, "/StgWtFldFastVolPF_");
        strcat(StgFastVolPFFile, stringIdWtField);
        strcat(StgFastVolPFFile, ".out");
        fStgFastVolPF = fopen(StgFastVolPFFile, "r");

        if (fStgFastVolPF == NULL) {
          neBEMMessage("in neBEMSolve - StgWtFldFastVolPFFile");
          return -1;
        }

        fscanf(fStgFastVolPF, "#block\tX\tY\tZ\tPot\tFX\tFY\tFZ\n");

        for (int block = 1; block <= WtFldFastVol[IdWtField].NbBlocks;
             ++block) {
          nbXCells = WtFldBlkNbXCells[IdWtField][block];
          nbYCells = WtFldBlkNbYCells[IdWtField][block];
          nbZCells = WtFldBlkNbZCells[IdWtField][block];

          for (int i = 1; i <= nbXCells + 1; ++i) {
            for (int j = 1; j <= nbYCells + 1; ++j) {
              for (int k = 1; k <= nbZCells + 1; ++k) {
                fscanf(fStgFastVolPF, "%d\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n",
                       &tmpblk, &xpt, &ypt, &zpt,
                       &StgWtFldFastPot[IdWtField][block][i][j][k],
                       &StgWtFldFastFX[IdWtField][block][i][j][k],
                       &StgWtFldFastFY[IdWtField][block][i][j][k],
                       &StgWtFldFastFZ[IdWtField][block][i][j][k]);
              }  // loop k
            }    // loop j
          }      // loop i
        }        // loop block
        fclose(fStgFastVolPF);
      }  // if OptStaggerWtFldFastVol

      clock_t stopFastClock = clock();
      neBEMTimeElapsed(startFastClock, stopFastClock);
      printf("to read WtFldFastVolPF\n");
    }  // if OptReadWtFldFastPF
  }    // if OptWtFldFastVol

  return IdWtField;
}  // neBEMPrepareWeightingField ends

// Deallocates memory reserved for a weighting field
void neBEMDeleteWeightingField(int IdWtField) {
  free(WtFieldChDen[IdWtField]);
  free(AvWtChDen[IdWtField]);
}

// Deallocates all memory reserved for all weighting fields
void neBEMDeleteAllWeightingFields(void) {
  const int MaxWtField = MAXWtFld;  // being used while allocating memory
  for (int id = 1; id < MaxWtField; ++id) {  // count from 1
    free(WtFieldChDen[id]);
    free(AvWtChDen[id]);
  }
  free(WtFieldChDen);
  free(AvWtChDen);
}

// Get weighting field (potential also) at a specific point
// returns DBL_MAX as the value of potential when something goes wrong.
double neBEMWeightingField(Point3D *point, Vector3D *field, int IdWtField) {
  double potential;

  if (neBEMState < 9) {
    printf("neBEMWeightingField cannot be called before reaching state 9.\n");
    return (-1);
  }

  if (OptFixedWtField[IdWtField]) {
    // minimum computation, too restricted! Does not even consider variation
    // through IdWtField
    potential = FixedWtPotential[IdWtField];
    field->X = FixedWtFieldX[IdWtField];
    field->Y = FixedWtFieldY[IdWtField];
    field->Z = FixedWtFieldZ[IdWtField];
  } else if (OptWtFldFastVol[IdWtField]) {
    // bit more computation, lot more flexibility
    // Note: this is not the Create or Read option
    int fstatus = WtFldFastPFAtPoint(point, &potential, field, IdWtField);
    if (fstatus != 0) {
      neBEMMessage("neBEMWeightingField - WtFldFastPFAtPoint");
      return DBL_MAX;
    }
  } else {
    int fstatus = WtFldPFAtPoint(point, &potential, field, IdWtField);
    if (fstatus != 0) {
      neBEMMessage("neBEMWeightingField - WtFldPFAtPoint");
      return DBL_MAX;
    }
  }

  return potential;
}  // neBEMWeightingField ends

double neBEMVolumeCharge(int volume) {
  // Initialise the sum
  double sumcharge = 0.0;

  // Loop over all elements
  for (int elem = 1; elem <= NbElements; ++elem) {
    // Find the primitive number for the element
    int prim = (EleArr + elem - 1)->PrimitiveNb;

    // Find out to which volume this belongs to
    int volref = VolRef1[prim];

    // Skip the rest if the volume is not right
    if (volref + 1 != volume) {
      continue;
    }

    // Calculate the periodic weight of the primitive
    int rptCnt = 0;
    if (PeriodicInX[prim] || PeriodicInY[prim] || PeriodicInZ[prim]) {
      for (int xrpt = -PeriodicInX[prim]; xrpt <= PeriodicInX[prim]; ++xrpt)
        for (int yrpt = -PeriodicInY[prim]; yrpt <= PeriodicInY[prim]; ++yrpt)
          for (int zrpt = -PeriodicInZ[prim]; zrpt <= PeriodicInZ[prim];
               ++zrpt) {
            if ((xrpt == 0) && (yrpt == 0) && (zrpt == 0))
              continue;
            else
              ++rptCnt;
          }
    } else {
      rptCnt = 1;
    }

    // Add the charge
    // printf("Element: %d, volume: %d, charge: %g\n", elem, volref,
    //    (EleArr+elem-1)->Solution * (EleArr+elem-1)->G.dA);
    sumcharge +=
        rptCnt * (EleArr + elem - 1)->Solution * (EleArr + elem - 1)->G.dA;
  }  // loop over elements

  // Return the result
  // printf("Charge on volume %d: %g C\n", volume, sumcharge);
  return sumcharge;
}  // end of neBEMVolumeCharge

int neBEMEnd(void) {
  fprintf(fIsles,
          "IslesCntr: %d, ExactCntr: %d, FailureCntr: %d, ApproxCntr: %d\n",
          IslesCntr, ExactCntr, FailureCntr, ApproxCntr);
  fclose(fIsles);
  fIsles = NULL;
  printf("neBEM ends ... bye!\n");

  return 0;
}  // neBEMEnd ends

// In a given problem, the DeviceOutDir is the one that is unique for a given
// device. Within it, there are several sub-directories, one related to each
// device counter; within each device counter directory, there can be
// several sub-directories related to each Mesh specification;
// within each mesh sub-dir, there can be several sub-directories, one for each
// set of boundary conditions and finally, for each boundary condition,
// several related to each post-processing counter.
int CreateDirStr(void) {
  char strModelCntr[10], strMeshCntr[10], strBCCntr[10], strPPCntr[10];
  int CreateOrUseDir(char[]);
  int CreateDirOrQuit(char[]);

  sprintf(strModelCntr, "/Model%d", ModelCntr);
  sprintf(strMeshCntr, "/M%d", MeshCntr);
  sprintf(strBCCntr, "/BC%d", BCCntr);
  sprintf(strPPCntr, "/PP%d", PPCntr);

  strcpy(ModelOutDir, DeviceOutDir);
  strcat(ModelOutDir, strModelCntr);
  strcpy(NativeOutDir, ModelOutDir);
  strcat(NativeOutDir, "/neBEMNatives/");
  strcpy(NativePrimDir, NativeOutDir);
  strcat(NativePrimDir, "Primitives/");
  strcpy(MeshOutDir, ModelOutDir);
  strcat(MeshOutDir, strMeshCntr);
  strcpy(BCOutDir, MeshOutDir);
  strcat(BCOutDir, strBCCntr);
  strcpy(PPOutDir, BCOutDir);
  strcat(PPOutDir, strPPCntr);

  // Create DeviceOutDir, if necessary
  int fstatus = CreateOrUseDir(DeviceOutDir);
  if (fstatus != 0) {
    neBEMMessage("CreateDirStr - CreateOrUseDir");
    return -1;
  }

  if (NewModel) {
    // create ModelOutDir
    if (OptReuseDir) {
      fstatus = CreateOrUseDir(ModelOutDir);
      fstatus = CreateOrUseDir(NativeOutDir);
      fstatus = CreateOrUseDir(NativePrimDir);
    } else {
      fstatus = CreateDirOrQuit(ModelOutDir);
      fstatus = CreateDirOrQuit(NativeOutDir);
      fstatus = CreateDirOrQuit(NativePrimDir);
    }
    if (fstatus != 0) {
      neBEMMessage("CreateDirStr - ModelOutDir");
      return -1;
    }
  }

  if (NewMesh) {
    // create MeshOutDir
    if (OptReuseDir)
      fstatus = CreateOrUseDir(MeshOutDir);
    else
      fstatus = CreateDirOrQuit(MeshOutDir);
    if (fstatus != 0) {
      neBEMMessage("CreateDirStr - MeshOutDir");
      return -1;
    }
  }

  if (NewBC) {
    // create BCOutDir
    if (OptReuseDir)
      fstatus = CreateOrUseDir(BCOutDir);
    else
      fstatus = CreateDirOrQuit(BCOutDir);
    if (fstatus != 0) {
      neBEMMessage("CreateDirStr - BCOutDir");
      return -1;
    }
  }

  if (NewPP) {
    // create PPOutDir
    if (OptReuseDir)
      fstatus = CreateOrUseDir(PPOutDir);
    else
      fstatus = CreateDirOrQuit(PPOutDir);
    if (fstatus != 0) {
      neBEMMessage("CreateDirStr - PPOutDir");
      return -1;
    }
  }

  // Create other relevant sub-directories
  char subdir[256];

  strcpy(subdir, ModelOutDir);
  strcat(subdir, "/Primitives/");
  if (OptReuseDir)
    fstatus = CreateOrUseDir(subdir);
  else
    fstatus = CreateDirOrQuit(subdir);

  strcpy(subdir, MeshOutDir);
  strcat(subdir, "/Elements/");
  if (OptReuseDir)
    fstatus = CreateOrUseDir(subdir);
  else
    fstatus = CreateDirOrQuit(subdir);

  strcpy(subdir, MeshOutDir);
  strcat(subdir, "/GViewDir/");
  if (OptReuseDir)
    fstatus = CreateOrUseDir(subdir);
  else
    fstatus = CreateDirOrQuit(subdir);

  return (0);
}  // CreateDirStr ends

// Create or use directory dirname
int CreateOrUseDir(char dirname[]) {
  char dirstr[256];
  struct stat st;

  if (stat(dirname, &st) == 0)  // feel safe to use an existing directory
  {
    printf("Previous %s exists ... using the existing directory ... \n",
           dirname);
  } else {
    sprintf(dirstr, "mkdir -p %s", dirname);
    if (system(dirstr))  // returns 0 if successful
    {
      printf("Cannot create dirname %s ... returning ...\n", dirname);
      return (-1);
    }
  }

  return (0);
}  // CreateOrUseDir ends

// Create directory dirname or quit reporting failure
int CreateDirOrQuit(char dirname[]) {
  char dirstr[256];
  struct stat st;

  if (stat(dirname, &st) == 0)  // not safe to use an existing directory
  {
    printf("Previous %s exists ... please check inputs and counters ... \n",
           dirname);
    return (-1);
  } else {
    sprintf(dirstr, "mkdir -p %s", dirname);
    if (system(dirstr))  // returns 0 if successful
    {
      printf("Cannot create dirname %s ... returning ...\n", dirname);
      return (-1);
    }
  }

  return (0);
}  // CreateOrQuitDir ends

int neBEMMessage(const char *message) {
  fprintf(stdout, "neBEMMessage: %s\n", message);

  return 0;
}  // neBEMMessage ends

int WritePrimitives(void) {
  char PrimitiveFile[256];

  strcpy(PrimitiveFile, ModelOutDir);
  strcat(PrimitiveFile, "/Primitives/StorePrims.out");

  FILE *fStrPrm = fopen(PrimitiveFile, "w");
  if (fStrPrm == NULL) {
    neBEMMessage("WritePrimitives - Could not create file to store primitives");
    return -1;
  }

  fprintf(fStrPrm, "%d %d\n", NbVolumes, VolMax);
  fprintf(fStrPrm, "%d\n", NbPrimitives);
  fprintf(fStrPrm, "%d\n", MaxNbVertices);

  for (int prim = 1; prim <= NbPrimitives; ++prim) {
    fprintf(fStrPrm, "%d\n", PrimType[prim]);
    fprintf(fStrPrm, "%d\n", InterfaceType[prim]);
    fprintf(fStrPrm, "%d\n", NbVertices[prim]);

    for (int vert = 0; vert < NbVertices[prim]; ++vert) {
      fprintf(fStrPrm, "%le %le %le\n", XVertex[prim][vert],
              YVertex[prim][vert], ZVertex[prim][vert]);
    }  // vert loop

    fprintf(fStrPrm, "%le %le %le\n", XNorm[prim], YNorm[prim], ZNorm[prim]);
    fprintf(fStrPrm, "%le\n", Radius[prim]);

    fprintf(fStrPrm, "%le %le %le %le %le\n", Epsilon1[prim], Epsilon2[prim],
            Lambda[prim], ApplPot[prim], ApplCh[prim]);

    fprintf(fStrPrm, "%d %d\n", VolRef1[prim], VolRef2[prim]);

    fprintf(fStrPrm, "%d %d %d\n", PeriodicTypeX[prim], PeriodicTypeY[prim],
            PeriodicTypeZ[prim]);
    fprintf(fStrPrm, "%d %d %d\n", PeriodicInX[prim], PeriodicInY[prim],
            PeriodicInZ[prim]);
    fprintf(fStrPrm, "%le %le %le\n", XPeriod[prim], YPeriod[prim],
            ZPeriod[prim]);
    fprintf(fStrPrm, "%le %le %le\n", MirrorDistXFromOrigin[prim],
            MirrorDistYFromOrigin[prim], MirrorDistZFromOrigin[prim]);
  }  // prim loop

  fclose(fStrPrm);

  return 0;
}  // WritePrimitives ends

int WriteElements(void) {
  char ElementFile[256];

  strcpy(ElementFile, MeshOutDir);
  strcat(ElementFile, "/Elements/StoreElems.out");

  FILE *fStrEle = fopen(ElementFile, "w");
  if (fStrEle == NULL) {
    neBEMMessage("WriteElements - Could not create file to store elements");
    return -1;
  }

  fprintf(fStrEle, "%d %d\n", NbSurfs, NbWires);

  for (int prim = 1; prim <= NbPrimitives; ++prim) {
    fprintf(fStrEle, "%d %d\n", NbSurfSegX[prim], NbSurfSegZ[prim]);
    fprintf(fStrEle, "%d\n", NbWireSeg[prim]);
  }

  fprintf(fStrEle, "%d\n", NbElements);

  for (int ele = 1; ele <= NbElements; ++ele) {
    fprintf(fStrEle, "%d %d %d %d %d\n", (EleArr + ele - 1)->DeviceNb,
            (EleArr + ele - 1)->ComponentNb, (EleArr + ele - 1)->PrimitiveNb,
            (EleArr + ele - 1)->InterfaceId, (EleArr + ele - 1)->Id);
    fprintf(fStrEle, "%d %le %le %le %le %le %le\n", (EleArr + ele - 1)->G.Type,
            (EleArr + ele - 1)->G.Origin.X, (EleArr + ele - 1)->G.Origin.Y,
            (EleArr + ele - 1)->G.Origin.Z, (EleArr + ele - 1)->G.LX,
            (EleArr + ele - 1)->G.LZ, (EleArr + ele - 1)->G.dA);
    fprintf(fStrEle, "%le %le %le\n", (EleArr + ele - 1)->G.DC.XUnit.X,
            (EleArr + ele - 1)->G.DC.XUnit.Y, (EleArr + ele - 1)->G.DC.XUnit.Z);
    fprintf(fStrEle, "%le %le %le\n", (EleArr + ele - 1)->G.DC.YUnit.X,
            (EleArr + ele - 1)->G.DC.YUnit.Y, (EleArr + ele - 1)->G.DC.YUnit.Z);
    fprintf(fStrEle, "%le %le %le\n", (EleArr + ele - 1)->G.DC.ZUnit.X,
            (EleArr + ele - 1)->G.DC.ZUnit.Y, (EleArr + ele - 1)->G.DC.ZUnit.Z);
    fprintf(fStrEle, "%d %le\n", (EleArr + ele - 1)->E.Type,
            (EleArr + ele - 1)->E.Lambda);
    fprintf(fStrEle, "%d %le %le %le %le\n", (EleArr + ele - 1)->BC.NbOfBCs,
            (EleArr + ele - 1)->BC.CollPt.X, (EleArr + ele - 1)->BC.CollPt.Y,
            (EleArr + ele - 1)->BC.CollPt.Z, (EleArr + ele - 1)->BC.Value);
    fprintf(fStrEle, "%le %le\n", (EleArr + ele - 1)->Solution,
            (EleArr + ele - 1)->Assigned);
  }

  fprintf(fStrEle, "%d %d %d %d\n", NbPointsKnCh, NbLinesKnCh, NbAreasKnCh,
          NbVolumesKnCh);

  for (int pt = 1; pt <= NbPointsKnCh; ++pt) {
    fprintf(fStrEle, "%d %le\n", (PointKnChArr + pt - 1)->Nb,
            (PointKnChArr + pt - 1)->Assigned);
    fprintf(fStrEle, "%le %le %le\n", (PointKnChArr + pt - 1)->P.X,
            (PointKnChArr + pt - 1)->P.Y, (PointKnChArr + pt - 1)->P.Z);
  }

  for (int line = 1; line <= NbLinesKnCh; ++line) {
    fprintf(fStrEle, "%d %le %le\n", (LineKnChArr + line - 1)->Nb,
            (LineKnChArr + line - 1)->Radius,
            (LineKnChArr + line - 1)->Assigned);
    fprintf(fStrEle, "%le %le %le\n", (LineKnChArr + line - 1)->Start.X,
            (LineKnChArr + line - 1)->Start.Y,
            (LineKnChArr + line - 1)->Start.Z);
    fprintf(fStrEle, "%le %le %le\n", (LineKnChArr + line - 1)->Stop.X,
            (LineKnChArr + line - 1)->Stop.Y, (LineKnChArr + line - 1)->Stop.Z);
  }

  for (int area = 1; area <= NbAreasKnCh; ++area) {
    fprintf(fStrEle, "%d %d %le\n", (AreaKnChArr + area - 1)->Nb,
            (AreaKnChArr + area - 1)->NbVertices,
            (AreaKnChArr + area - 1)->Assigned);
    for (int vert = 1; vert <= (AreaKnChArr + area - 1)->NbVertices; ++vert) {
      fprintf(fStrEle, "%le %le %le\n",
              (AreaKnChArr + area - 1)->Vertex[vert].X,
              (AreaKnChArr + area - 1)->Vertex[vert].Y,
              (AreaKnChArr + area - 1)->Vertex[vert].Z);
    }
  }

  for (int vol = 1; vol <= NbVolumesKnCh; ++vol) {
    fprintf(fStrEle, "%d %d %le\n", (VolumeKnChArr + vol - 1)->Nb,
            (VolumeKnChArr + vol - 1)->NbVertices,
            (VolumeKnChArr + vol - 1)->Assigned);
    for (int vert = 1; vert <= (VolumeKnChArr + vol - 1)->NbVertices; ++vert) {
      fprintf(fStrEle, "%le %le %le\n",
              (VolumeKnChArr + vol - 1)->Vertex[vert].X,
              (VolumeKnChArr + vol - 1)->Vertex[vert].Y,
              (VolumeKnChArr + vol - 1)->Vertex[vert].Z);
    }
  }

  fclose(fStrEle);

  return 0;
}  // WriteElements ends

int ReadPrimitives(void) {
  int dbgFn = 0;

  char PrimitiveFile[256];

  strcpy(PrimitiveFile, ModelOutDir);
  strcat(PrimitiveFile, "/Primitives/StorePrims.out");

  FILE *fStrPrm = fopen(PrimitiveFile, "r");
  if (fStrPrm == NULL) {
    neBEMMessage("ReadPrimitives - Could not open file to read primitives");
    return -1;
  }

  fscanf(fStrPrm, "%d %d\n", &NbVolumes, &VolMax);
  fscanf(fStrPrm, "%d\n", &NbPrimitives);
  fscanf(fStrPrm, "%d\n", &MaxNbVertices);

  // assign neBEMState and allocate memory
  neBEMState = 2;
  PrimType = ivector(1, NbPrimitives);
  NbVertices = ivector(1, NbPrimitives);
  XVertex = dmatrix(1, NbPrimitives, 0, MaxNbVertices - 1);
  YVertex = dmatrix(1, NbPrimitives, 0, MaxNbVertices - 1);
  ZVertex = dmatrix(1, NbPrimitives, 0, MaxNbVertices - 1);
  XNorm = dvector(1, NbPrimitives);
  YNorm = dvector(1, NbPrimitives);
  ZNorm = dvector(1, NbPrimitives);
  Radius = dvector(1, NbPrimitives);  // can lead to a little memory misuse
  VolRef1 = ivector(1, NbPrimitives);
  VolRef2 = ivector(1, NbPrimitives);
  NbSurfSegX = ivector(1, NbPrimitives);
  NbSurfSegZ = ivector(1, NbPrimitives);
  NbWireSeg = ivector(1, NbPrimitives);  // little memory misuse
  InterfaceType = ivector(1, NbPrimitives);
  Lambda = dvector(1, NbPrimitives);
  ApplPot = dvector(1, NbPrimitives);
  ApplCh = dvector(1, NbPrimitives);
  PeriodicTypeX = ivector(1, NbPrimitives);
  PeriodicTypeY = ivector(1, NbPrimitives);
  PeriodicTypeZ = ivector(1, NbPrimitives);
  PeriodicInX = ivector(1, NbPrimitives);
  PeriodicInY = ivector(1, NbPrimitives);
  PeriodicInZ = ivector(1, NbPrimitives);
  XPeriod = dvector(1, NbPrimitives);
  YPeriod = dvector(1, NbPrimitives);
  ZPeriod = dvector(1, NbPrimitives);
  MirrorDistXFromOrigin = dvector(1, NbPrimitives);
  MirrorDistYFromOrigin = dvector(1, NbPrimitives);
  MirrorDistZFromOrigin = dvector(1, NbPrimitives);

  for (int prim = 1; prim <= NbPrimitives; ++prim) {
    fscanf(fStrPrm, "%d\n", &PrimType[prim]);
    fscanf(fStrPrm, "%d\n", &InterfaceType[prim]);
    fscanf(fStrPrm, "%d\n", &NbVertices[prim]);

    for (int vert = 0; vert < NbVertices[prim]; ++vert) {
      fscanf(fStrPrm, "%le %le %le\n", &XVertex[prim][vert],
             &YVertex[prim][vert], &ZVertex[prim][vert]);
    }  // vert loop

    fscanf(fStrPrm, "%le %le %le\n", &XNorm[prim], &YNorm[prim], &ZNorm[prim]);
    fscanf(fStrPrm, "%le\n", &Radius[prim]);

    fscanf(fStrPrm, "%le %le %le %le %le\n", &Epsilon1[prim], &Epsilon2[prim],
           &Lambda[prim], &ApplPot[prim], &ApplCh[prim]);

    fscanf(fStrPrm, "%d %d\n", &VolRef1[prim], &VolRef2[prim]);

    fscanf(fStrPrm, "%d %d %d\n", &PeriodicTypeX[prim], &PeriodicTypeY[prim],
           &PeriodicTypeZ[prim]);
    fscanf(fStrPrm, "%d %d %d\n", &PeriodicInX[prim], &PeriodicInY[prim],
           &PeriodicInZ[prim]);
    fscanf(fStrPrm, "%le %le %le\n", &XPeriod[prim], &YPeriod[prim],
           &ZPeriod[prim]);
    fscanf(fStrPrm, "%le %le %le\n", &MirrorDistXFromOrigin[prim],
           &MirrorDistYFromOrigin[prim], &MirrorDistZFromOrigin[prim]);
  }  // prim loop

  volRef = ivector(0, VolMax);
  volShape = ivector(0, VolMax);
  volMaterial = ivector(0, VolMax);
  volEpsilon = dvector(0, VolMax);
  volPotential = dvector(0, VolMax);
  volCharge = dvector(0, VolMax);
  volBoundaryType = ivector(0, VolMax);
  for (int volref = 0; volref <= VolMax; ++volref) {
    neBEMVolumeDescription(volref, &volShape[volref], &volMaterial[volref],
                           &volEpsilon[volref], &volPotential[volref],
                           &volCharge[volref], &volBoundaryType[volref]);
    if (dbgFn) {
      printf("volref: %d\n", volref);
      printf("shape: %d,  material: %d\n", volShape[volref],
             volMaterial[volref]);
      printf("eps: %lg,  pot: %lg\n", volEpsilon[volref], volPotential[volref]);
      printf("q: %lg,  type: %d\n", volCharge[volref], volBoundaryType[volref]);
    }
  }  // volume loop

  fclose(fStrPrm);

  return 0;
}  // ReadPrimitives ends

int ReadElements(void) {
  char ElementFile[256];

  strcpy(ElementFile, MeshOutDir);
  strcat(ElementFile, "/Elements/StoreElems.out");

  FILE *fStrEle = fopen(ElementFile, "r");
  if (fStrEle == NULL) {
    neBEMMessage("ReadElements - Could not open file to read elements");
    return -1;
  }

  fscanf(fStrEle, "%d %d\n", &NbSurfs, &NbWires);

  for (int prim = 1; prim <= NbPrimitives; ++prim) {
    fscanf(fStrEle, "%d %d\n", &NbSurfSegX[prim], &NbSurfSegZ[prim]);
    fscanf(fStrEle, "%d\n", &NbWireSeg[prim]);
  }

  fscanf(fStrEle, "%d\n", &NbElements);

  if (neBEMState == 3) {
    if (EleArr) {
      Element *tmp = (Element *)realloc(EleArr, NbElements * sizeof(Element));
      if (tmp != NULL) {
        EleArr = tmp;
        EleCntr = 0;
      } else {
        free(EleArr);
        printf("neBEMDiscretize: Re-allocating EleArr failed.\n");
        return (1);
      }
      printf("neBEMDiscretize: Re-allocated EleArr.\n");
    }  // if EleArr => re-allocation
    else {
      EleArr = (Element *)malloc(NbElements * sizeof(Element));
      if (EleArr == NULL) {
        neBEMMessage("neBEMDiscretize - EleArr malloc");
        return -1;
      }
    }  // else EleArr => fresh allocation
  } else {
    neBEMMessage("neBEMDiscretize - EleArr malloc; neBEMState mismatch!");
    return -1;
  }  // else neBEMState == 3

  for (int ele = 1; ele <= NbElements; ++ele) {
    fscanf(fStrEle, "%hd %d %d %d %d\n", &(EleArr + ele - 1)->DeviceNb,
           &(EleArr + ele - 1)->ComponentNb, &(EleArr + ele - 1)->PrimitiveNb,
           &(EleArr + ele - 1)->InterfaceId, &(EleArr + ele - 1)->Id);
    fscanf(fStrEle, "%hd %le %le %le %le %le %le\n",
           &(EleArr + ele - 1)->G.Type, &(EleArr + ele - 1)->G.Origin.X,
           &(EleArr + ele - 1)->G.Origin.Y, &(EleArr + ele - 1)->G.Origin.Z,
           &(EleArr + ele - 1)->G.LX, &(EleArr + ele - 1)->G.LZ,
           &(EleArr + ele - 1)->G.dA);
    fscanf(fStrEle, "%le %le %le\n", &(EleArr + ele - 1)->G.DC.XUnit.X,
           &(EleArr + ele - 1)->G.DC.XUnit.Y,
           &(EleArr + ele - 1)->G.DC.XUnit.Z);
    fscanf(fStrEle, "%le %le %le\n", &(EleArr + ele - 1)->G.DC.YUnit.X,
           &(EleArr + ele - 1)->G.DC.YUnit.Y,
           &(EleArr + ele - 1)->G.DC.YUnit.Z);
    fscanf(fStrEle, "%le %le %le\n", &(EleArr + ele - 1)->G.DC.ZUnit.X,
           &(EleArr + ele - 1)->G.DC.ZUnit.Y,
           &(EleArr + ele - 1)->G.DC.ZUnit.Z);
    fscanf(fStrEle, "%hd %le\n", &(EleArr + ele - 1)->E.Type,
           &(EleArr + ele - 1)->E.Lambda);
    fscanf(fStrEle, "%hd %le %le %le %le\n", &(EleArr + ele - 1)->BC.NbOfBCs,
           &(EleArr + ele - 1)->BC.CollPt.X, &(EleArr + ele - 1)->BC.CollPt.Y,
           &(EleArr + ele - 1)->BC.CollPt.Z, &(EleArr + ele - 1)->BC.Value);
    fscanf(fStrEle, "%le %le\n", &(EleArr + ele - 1)->Solution,
           &(EleArr + ele - 1)->Assigned);
  }

  fscanf(fStrEle, "%d %d %d %d\n", &NbPointsKnCh, &NbLinesKnCh, &NbAreasKnCh,
         &NbVolumesKnCh);

  for (int pt = 1; pt <= NbPointsKnCh; ++pt) {
    fscanf(fStrEle, "%d %le\n", &(PointKnChArr + pt - 1)->Nb,
           &(PointKnChArr + pt - 1)->Assigned);
    fscanf(fStrEle, "%le %le %le\n", &(PointKnChArr + pt - 1)->P.X,
           &(PointKnChArr + pt - 1)->P.Y, &(PointKnChArr + pt - 1)->P.Z);
  }

  for (int line = 1; line <= NbLinesKnCh; ++line) {
    fscanf(fStrEle, "%d %le %le\n", &(LineKnChArr + line - 1)->Nb,
           &(LineKnChArr + line - 1)->Radius,
           &(LineKnChArr + line - 1)->Assigned);
    fscanf(fStrEle, "%le %le %le\n", &(LineKnChArr + line - 1)->Start.X,
           &(LineKnChArr + line - 1)->Start.Y,
           &(LineKnChArr + line - 1)->Start.Z);
    fscanf(fStrEle, "%le %le %le\n", &(LineKnChArr + line - 1)->Stop.X,
           &(LineKnChArr + line - 1)->Stop.Y,
           &(LineKnChArr + line - 1)->Stop.Z);
  }

  for (int area = 1; area <= NbAreasKnCh; ++area) {
    fscanf(fStrEle, "%d %d %le\n", &(AreaKnChArr + area - 1)->Nb,
           &(AreaKnChArr + area - 1)->NbVertices,
           &(AreaKnChArr + area - 1)->Assigned);
    for (int vert = 1; vert <= (AreaKnChArr + area - 1)->NbVertices; ++vert) {
      fscanf(fStrEle, "%le %le %le\n",
             &(AreaKnChArr + area - 1)->Vertex[vert].X,
             &(AreaKnChArr + area - 1)->Vertex[vert].Y,
             &(AreaKnChArr + area - 1)->Vertex[vert].Z);
    }
  }

  for (int vol = 1; vol <= NbVolumesKnCh; ++vol) {
    fscanf(fStrEle, "%d %d %le\n", &(VolumeKnChArr + vol - 1)->Nb,
           &(VolumeKnChArr + vol - 1)->NbVertices,
           &(VolumeKnChArr + vol - 1)->Assigned);
    for (int vert = 1; vert <= (VolumeKnChArr + vol - 1)->NbVertices; ++vert) {
      fscanf(fStrEle, "%le %le %le\n",
             &(VolumeKnChArr + vol - 1)->Vertex[vert].X,
             &(VolumeKnChArr + vol - 1)->Vertex[vert].Y,
             &(VolumeKnChArr + vol - 1)->Vertex[vert].Z);
    }
  }

  fclose(fStrEle);

  return 0;
}  // ReadElements ends

// returns the time elapsed between start and stop
double neBEMTimeElapsed(clock_t t0, clock_t t1) {
  double elapsedTime = ((double)(t1 - t0)) / CLOCKS_PER_SEC;
  printf("neBEMV%s TimeElapsed ===> %lg seconds ", neBEMVersion, elapsedTime);

  return (elapsedTime);
}

// returns number of lines in file fname
// Interface.h
int neBEMGetNbOfLines(const char fname[]) {
  unsigned int number_of_lines = 0;

  FILE *infile = fopen(fname, "r");

  int ch;
  while (EOF != (ch = getc(infile)))
    if ('\n' == ch) ++number_of_lines;

  return number_of_lines;
}

// check whether a 3D point is within a 3D polygon
// Downloads/ComNum/Geometry/Determining if a point lies on the interior of a
// polygon.htm
double neBEMChkInPoly(int n, Point3D *p, Point3D q) {
  double anglesum = 0.0;
  // double theta = 0.0;
  Point3D p1, p2;

  // printf("In neBEMChkInPoly ... \n");
  // printf("n: %d\n", n);

  for (int i = 0; i < n; i++) {
    // printf("i: %d\n", i);
    p1.X = p[i].X - q.X;
    p1.Y = p[i].Y - q.Y;
    p1.Z = p[i].Z - q.Z;
    // printf("p[%d]: %lg %lg %lg\n", i, p[i].X, p[i].Y, p[i].Z);
    // printf("q: %lg %lg %lg\n", q.X, q.Y, q.Z);

    if (i < n - 1) {
      p2.X = p[i + 1].X - q.X;
      p2.Y = p[i + 1].Y - q.Y;
      p2.Z = p[i + 1].Z - q.Z;
    } else {
      p2.X = p[0].X - q.X;
      p2.Y = p[0].Y - q.Y;
      p2.Z = p[0].Z - q.Z;
    }

    double m1 = MODULUS(p1);
    double m2 = MODULUS(p2);
    // printf("m1: %lg, m2: %lg, m1*m2: %lg", m1, m2, m1*m2);

    if (m1 * m2 <= 1.0e-12) {
      // vetors of 1 micron - we may need to reduce further
      return (neBEMtwopi); /* We are on a node, consider this inside */
    }
    double costheta = (p1.X * p2.X + p1.Y * p2.Y + p1.Z * p2.Z) / (m1 * m2);

    /*
    double oldtheta = theta;
    theta = acos(costheta);
    // printf("n: %d, i: %d, theta: %lg\n", n, i, neBEMrtod*theta);
    if (Sign(theta) != Sign(oldtheta)) {
      // polygon either non-covex, or the point is outside the polygon
      return(0.0);  // absurd value implying outside polygon
    }
    */
    anglesum += acos(costheta);
    // printf("n: %d, i: %d, anglesum: %lg %lg\n", n, i, anglesum,
    // neBEMrtod*anglesum);
  }

  return (anglesum);
}  // neBEMChkInPoly

#ifdef __cplusplus
}  // namespace
#endif