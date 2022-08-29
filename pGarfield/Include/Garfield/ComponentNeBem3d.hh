#ifndef G_COMPONENT_NEBEM_3D_H
#define G_COMPONENT_NEBEM_3D_H

#include <map>

#include "Component.hh"

namespace Garfield {

/// Interface to neBEM.

class ComponentNeBem3d : public Component {
 public:
  /// Constructor
  ComponentNeBem3d();
  /// Destructor
  ~ComponentNeBem3d() {}

  Medium* GetMedium(const double x, const double y, const double z) override;

  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, Medium*& m, int& status) override;
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, double& v, Medium*& m,
                     int& status) override;
  bool GetVoltageRange(double& vmin, double& vmax) override;

  void WeightingField(const double x, const double y, const double z,
                      double& wx, double& wy, double& wz,
                      const std::string& label) override;
  double WeightingPotential(const double x, const double y, const double z,
                            const std::string& label) override;

  /// Add a plane at constant x.
  void AddPlaneX(const double x, const double voltage);
  /// Add a plane at constant y.
  void AddPlaneY(const double y, const double voltage);
  /// Add a plane at constant z.
  void AddPlaneZ(const double z, const double voltage);
  /// Get the number of equipotential planes at constant x.
  unsigned int GetNumberOfPlanesX() const;
  /// Get the number of equipotential planes at constant y.
  unsigned int GetNumberOfPlanesY() const;
  /// Get the number of equipotential planes at constant z.
  unsigned int GetNumberOfPlanesZ() const;
  /// Retrieve the parameters of a plane at constant x.
  bool GetPlaneX(const unsigned int i, double& x, double& v) const;
  /// Retrieve the parameters of a plane at constant y.
  bool GetPlaneY(const unsigned int i, double& y, double& v) const;
  /// Retrieve the parameters of a plane at constant z.
  bool GetPlaneZ(const unsigned int i, double& z, double& v) const;

  unsigned int GetNumberOfPrimitives() const { return m_primitives.size(); }
  bool GetPrimitive(const unsigned int i, double& a, double& b, double& c,
                    std::vector<double>& xv, std::vector<double>& yv,
                    std::vector<double>& zv, int& interface, double& v,
                    double& q, double& lambda) const;
  bool GetPrimitive(const unsigned int i, double& a, double& b, double& c,
                    std::vector<double>& xv, std::vector<double>& yv,
                    std::vector<double>& zv, int& vol1, int& vol2) const;
  bool GetVolume(const unsigned int vol, int& shape, int& material, double& eps,
                 double& potential, double& charge, int& bc);
  int GetVolume(const double x, const double y, const double z);

  unsigned int GetNumberOfElements() const { return m_elements.size(); }
  bool GetElement(const unsigned int i, std::vector<double>& xv,
                  std::vector<double>& yv, std::vector<double>& zv,
                  int& interface, double& bc, double& lambda) const;

  /// Retrieve surface panels, remove contacts and cut polygons to rectangles
  /// and right-angle triangles.
  bool Initialise();

  /// Set the default value of the target linear size of the elements
  /// produced by neBEM's discretisation process.
  void SetTargetElementSize(const double length);
  /// Set the smallest and largest allowed number of elements along
  /// the lenght of a primitive.
  void SetMinMaxNumberOfElements(const unsigned int nmin,
                                 const unsigned int nmax);

  void SetNewModel(const unsigned int NewModel);
  void SetNewMesh(const unsigned int NewMesh);
  void SetNewBC(const unsigned int NewBC);
  void SetNewPP(const unsigned int NewPP);
  void SetModelOptions(const unsigned int NewModel, const unsigned int NewMesh,
                       const unsigned int NewBC, const unsigned int NewPP);

  /// Set storing options (OptStoreInflMatrix, OptStoreInvMatrix,
  /// OptStoreInvMatrix, OptStoreInvMatrix)
  /// OptStorePrimitives, OptStorePrimitives)
  /// OptStoreElements, OptStoreElements)
  /// OptFormattedFile, OptUnformattedFile)
  void SetStoreInflMatrix(const unsigned int OptStoreInflMatrix);
  void SetReadInflMatrix(const unsigned int OptReadInflMatrix);
  void SetStoreInvMatrix(const unsigned int OptStoreInvMatrix);
  void SetReadInvMatrix(const unsigned int OptReadInvMatrix);
  void SetStorePrimitives(const unsigned int OptStorePrimitives);
  void SetReadPrimitives(const unsigned int OptReadPrimitives);
  void SetStoreElements(const unsigned int OptStoreElements);
  void SetReadElements(const unsigned int OptReadElements);
  void SetFormattedFile(const unsigned int OptFormattedFile);
  void SetUnformattedFile(const unsigned int OptUnformattedFile);
  void SetStoreReadOptions(const unsigned int OptStoreInflMatrix,
                           const unsigned int OptReadInflMatrix,
                           const unsigned int OptStoreInvMatrix,
                           const unsigned int OptReadInvMatrix,
                           const unsigned int OptStorePrimitives,
                           const unsigned int OptReadPrimitives,
                           const unsigned int OptStoreElements,
                           const unsigned int OptReadElements,
                           const unsigned int OptFormattedFile,
                           const unsigned int OptUnformattedFile);

  // Set whether an older model is to be re-used. Expects an inverted matrix
  // stored during an earlier computation that had identical model and mesh.
  // This execution considers only a change in the boundary conditions.
  void SetReuseModel(void);

  /// Other functions to be, are
  /// void SetPlotOptions(OptGnuplot=0, OptGnuplotPrimitives=0,
  /// OptGnuplotElements=0,
  /// OptPrimitiveFiles=0, OptElementFiles=0)

  // Functions that set computation details and constraints
  void SetSystemChargeZero(const unsigned int OptSystemChargeZero);
  void SetValidateSolution(const unsigned int OptValidateSolution);
  void SetForceValidation(const unsigned int OptForceValidation);
  void SetRepeatLHMatrix(const unsigned int OptRepeatLHMatrix);
  void SetComputeOptions(const unsigned int OptSystemChargeZero,
                         const unsigned int OptValidateSolution,
                         const unsigned int OptForceValidation,
                         const unsigned int OptRepeatLHMatrix);

  // Fast volume related information (physical potential and fields)
  void SetFastVolOptions(const unsigned int OptFastVol,
                         const unsigned int OptCreateFastPF,
                         const unsigned int OptReadFastPF);
  void SetFastVolVersion(const unsigned int VersionFV);
  void SetFastVolBlocks(const unsigned int NbBlocksFV);

  // Needs to include IdWtField information for each of these WtFld functions
  // Weighting potential and field related Fast volume information
  void SetWtFldFastVolOptions(const unsigned int IdWtField,
                              const unsigned int OptWtFldFastVol,
                              const unsigned int OptCreateWtFldFastPF,
                              const unsigned int OptReadWtFldFastPF);
  void SetWtFldFastVolVersion(const unsigned int IdWtField,
                              const unsigned int VersionWtFldFV);
  void SetWtFldFastVolBlocks(const unsigned int IdWtField,
                             const unsigned int NbBlocksWtFldFV);

  // Known charge options
  void SetKnownChargeOptions(const unsigned int OptKnownCharge);

  // Charging up options
  void SetChargingUpOptions(const unsigned int OptChargingUp);

  /// Invert the influence matrix using lower-upper (LU) decomposition.
  void UseLUInversion() { m_inversion = Inversion::LU; }
  /// Invert the influence matrix using singular value decomposition.
  void UseSVDInversion() { m_inversion = Inversion::SVD; }

  /// Set the parameters \f$n_x, n_y, n_z\f$ defining the number of periodic
  /// copies that neBEM will use when dealing with periodic configurations.
  /// neBEM will use \f$2 \times n + 1\f$ copies (default: \f$n = 5\f$).
  void SetPeriodicCopies(const unsigned int nx, const unsigned int ny,
                         const unsigned int nz);
  /// Retrieve the number of periodic copies used by neBEM.
  void GetPeriodicCopies(unsigned int& nx, unsigned int& ny,
                         unsigned int& nz) const {
    nx = m_nCopiesX;
    ny = m_nCopiesY;
    nz = m_nCopiesZ;
  }
  /// Set the periodic length [cm] in the x-direction.
  void SetPeriodicityX(const double s);
  /// Set the periodic length [cm] in the y-direction.
  void SetPeriodicityY(const double s);
  /// Set the periodic length [cm] in the z-direction.
  void SetPeriodicityZ(const double s);
  /// Set the periodic length [cm] in the x-direction.
  void SetMirrorPeriodicityX(const double s);
  /// Set the periodic length [cm] in the y-direction.
  void SetMirrorPeriodicityY(const double s);
  /// Set the periodic length [cm] in the z-direction.
  void SetMirrorPeriodicityZ(const double s);
  /// Get the periodic length in the x-direction.
  bool GetPeriodicityX(double& s) const;
  /// Get the periodic length in the y-direction.
  bool GetPeriodicityY(double& s) const;
  /// Get the periodic length in the z-direction.
  bool GetPeriodicityZ(double& s) const;

  /// Set the number of threads to be used by neBEM.
  void SetNumberOfThreads(const unsigned int n) { m_nThreads = n > 0 ? n : 1; }

  /// Set the number of repetitions after which primitive properties are used
  /// for the physical field. 
  /// A negative value (default) implies all the elements are always evaluated.
  void SetPrimAfter(const int n) { m_primAfter = n; }

  /// Set the number of repetitions after which primitive properties are used
  /// for the weighting field.
  /// A negative value (default) implies all the elements are always evaluated.
  void SetWtFldPrimAfter(const int n) { m_wtFldPrimAfter = n; }

  /// Set option related to removal of primitives.
  void SetOptRmPrim(const unsigned int n) { m_optRmPrim = n; }

 protected:
  void Reset() override;
  void UpdatePeriodicity() override;

 private:
  struct Primitive {
    /// Perpendicular vector
    double a, b, c;
    /// X-coordinates of vertices
    std::vector<double> xv;
    /// Y-coordinates of vertices
    std::vector<double> yv;
    /// Z-coordinates of vertices
    std::vector<double> zv;
    /// Interface type.
    int interface;
    /// Potential
    double v;
    /// Charge
    double q;
    /// Ratio of dielectric constants
    double lambda;
    /// Target element size.
    double elementSize;
    /// Volumes.
    int vol1, vol2;
  };
  /// List of primitives.
  std::vector<Primitive> m_primitives;

  struct Element {
    /// Local origin.
    std::array<double, 3> origin;
    double lx;
    double lz;
    /// Area.
    double dA;
    /// Direction cosines.
    std::array<std::array<double, 3>, 3> dcos;
    /// X-coordinates of vertices
    std::vector<double> xv;
    /// Y-coordinates of vertices
    std::vector<double> yv;
    /// Z-coordinates of vertices
    std::vector<double> zv;
    /// Interface type.
    int interface;
    /// Ratio of dielectric permittivities.
    double lambda;
    /// Collocation point.
    std::array<double, 3> collocationPoint;
    /// Boundary condition.
    double bc;
    /// Fixed charge density.
    double assigned;
    /// Solution (accumulated charge).
    double solution;
  };
  /// List of elements.
  std::vector<Element> m_elements;

  /// Plane existence.
  std::array<bool, 6> m_ynplan{{false, false, false, false, false, false}};
  /// Plane coordinates.
  std::array<double, 6> m_coplan{{0., 0., 0., 0., 0., 0.}};
  /// Plane potentials.
  std::array<double, 6> m_vtplan{{0., 0., 0., 0., 0., 0.}};

  // Model specifications
  unsigned int m_newModel = 1;
  unsigned int m_newMesh = 1;
  unsigned int m_newBC = 1;
  unsigned int m_newPP = 1;

  // Store and read options
  unsigned int m_optStoreInflMatrix = 0;
  unsigned int m_optReadInflMatrix = 0;
  unsigned int m_optStoreInvMatrix = 1;
  unsigned int m_optReadInvMatrix = 0;
  unsigned int m_optStorePrimitives = 0;
  unsigned int m_optReadPrimitives = 0;
  unsigned int m_optStoreElements = 0;
  unsigned int m_optReadElements = 0;
  unsigned int m_optStoreFormatted = 1;
  unsigned int m_optStoreUnformatted = 0;

  // Plot options
  // unsigned int m_optGnuplotPrimitives = 0;
  // unsigned int m_optGnuplotElements = 0;
  // unsigned int m_optPrimitiveFiles = 0;
  // unsigned int m_optElementFiles = 0;

  // Compute options
  unsigned int m_optSystemChargeZero = 1;
  unsigned int m_optValidateSolution = 1;
  unsigned int m_optForceValidation = 0;
  unsigned int m_optRepeatLHMatrix = 0;

  // Fast volume information (physical potential and fields)
  unsigned int m_optFastVol = 0;
  unsigned int m_optCreateFastPF = 0;
  unsigned int m_optReadFastPF = 0;
  unsigned int m_versionFV = 0;
  unsigned int m_nbBlocksFV = 0;

  // Weighting potential and field related Fast volume information
  unsigned int m_idWtField = 0;
  unsigned int m_optWtFldFastVol[11];
  unsigned int m_optCreateWtFldFastPF[11];
  unsigned int m_optReadWtFldFastPF[11];
  unsigned int m_versionWtFldFV[11];
  unsigned int m_nbBlocksWtFldFV[11];

  // Known charge options
  unsigned int m_optKnownCharge = 0;

  // Charging up options
  unsigned int m_optChargingUp = 0;

  // Number of threads to be used by neBEM.
  unsigned int m_nThreads = 1;

  // Number of repetitions, after which only primitive properties are used.
  // a negative value implies elements are used always.
  int m_primAfter = -1;  

  // Number of repetitions, after which only primitive properties are used
  // for weighting field calculations.
  // A negative value implies only elements are used.
  int m_wtFldPrimAfter = -1;  

  // Option for removing primitives from a device geometry.
  // Zero implies none to be removed.
  unsigned int m_optRmPrim = 0;  

  static constexpr double MinDist = 1.e-6;
  /// Target size of elements [cm].
  double m_targetElementSize = 50.0e-4;
  /// Smallest number of elements produced along the axis of a primitive.
  unsigned int m_minNbElementsOnLength = 1;
  /// Largest number of elements produced along the axis of a primitive.
  unsigned int m_maxNbElementsOnLength = 100;
  /// Periodic lengths.
  std::array<double, 3> m_periodicLength{{0., 0., 0.}};
  /// Number of periodic copies along x.
  unsigned int m_nCopiesX = 5;
  /// Number of periodic copies along y.
  unsigned int m_nCopiesY = 5;
  /// Number of periodic copies along z.
  unsigned int m_nCopiesZ = 5;

  enum class Inversion { LU = 0, SVD };
  Inversion m_inversion = Inversion::LU;

  /// Electrode labels and corresponding neBEM weighting field indices.
  std::map<std::string, int> m_wfields;

  void InitValues();
  /// Reduce panels to the basic period.
  void ShiftPanels(std::vector<Panel>& panels) const;
  /// Isolate the parts of polygon 1 that are not hidden by 2 and vice versa.
  bool EliminateOverlaps(const Panel& panel1, const Panel& panel2,
                         std::vector<Panel>& panelsOut,
                         std::vector<int>& itypo);

  bool TraceEnclosed(const std::vector<double>& xl1,
                     const std::vector<double>& yl1,
                     const std::vector<double>& xl2,
                     const std::vector<double>& yl2, const Panel& originalPanel,
                     std::vector<Panel>& newPanels) const;

  void TraceNonOverlap(
      const std::vector<double>& xp1, const std::vector<double>& yp1,
      const std::vector<double>& xl1, const std::vector<double>& yl1,
      const std::vector<double>& xl2, const std::vector<double>& yl2,
      const std::vector<int>& flags1, const std::vector<int>& flags2,
      const std::vector<int>& links1, const std::vector<int>& links2,
      std::vector<bool>& mark1, int ip1, const Panel& originalPanel,
      std::vector<Panel>& newPanels) const;

  void TraceOverlap(
      const std::vector<double>& xp1, const std::vector<double>& yp1,
      const std::vector<double>& xp2, const std::vector<double>& yp2,
      const std::vector<double>& xl1, const std::vector<double>& yl1,
      const std::vector<double>& xl2, const std::vector<double>& yl2,
      const std::vector<int>& flags1, const std::vector<int>& links1,
      const std::vector<int>& links2, std::vector<bool>& mark1, int ip1,
      int ip2, const Panel& originalPanel, std::vector<Panel>& newPanels) const;

  /// Split a polygon into rectangles and right-angled triangles.
  bool MakePrimitives(const Panel& panelIn,
                      std::vector<Panel>& panelsOut) const;

  /// Check whether a polygon contains parallel lines.
  /// If it does, split it in rectangular and non-rectangular parts.
  bool SplitTrapezium(const Panel panelIn, std::vector<Panel>& stack,
                      std::vector<Panel>& panelsOut, const double epsang) const;

  unsigned int NbOfSegments(const double length, const double target) const;
  bool DiscretizeWire(const Primitive& primitive, const double targetSize,
                      std::vector<Element>& elements) const;
  bool DiscretizeTriangle(const Primitive& primitive, const double targetSize,
                          std::vector<Element>& elements) const;
  bool DiscretizeRectangle(const Primitive& prim, const double targetSize,
                           std::vector<Element>& elements) const;
  int InterfaceType(const Solid::BoundaryCondition bc) const;
};

extern ComponentNeBem3d* gComponentNeBem3d;

}  // namespace Garfield

#endif
