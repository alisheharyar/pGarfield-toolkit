#ifndef G_TRACK_SRIM_H
#define G_TRACK_SRIM_H
#include <vector>
#include "Track.hh"

namespace Garfield {

/// Generate tracks based on SRIM energy loss, range and straggling tables.
///  - http://www.srim.org

class TrackSrim : public Track {
 public:
  /// Constructor
  TrackSrim();
  /// Destructor
  virtual ~TrackSrim() {}

  /// Load data from a SRIM file.
  bool ReadFile(const std::string& file);
  /// Print the energy loss table. 
  void Print();
  /// Plot the electromagnetic, hadronic, and total energy loss
  /// as function of the projectile energy.
  void PlotEnergyLoss();
  /// Plot the projected range as function of the projectile energy.
  void PlotRange();
  /// Plot the transverse and longitudinal straggling as function
  /// of the projectile energy. 
  void PlotStraggling();

  /// Set the fluctuation model
  /// (0 = none, 1 = Landau, 2 = Vavilov, 3 = Gaussian, 4 = Combined).
  /// By default, the combined model (4) is used.
  void SetModel(const int m) { m_model = m; }
  /// Get the fluctuation model.
  int GetModel() const { return m_model; }

  /// Specify how many electrons should be grouped to a cluster.
  void SetTargetClusterSize(const int n) { m_nsize = n; }
  /// Retrieve the target cluster size.
  int GetTargetClusterSize() const { return m_nsize; }

  /// Set the max. number of clusters on a track.
  void SetClustersMaximum(const int n) { m_maxclusters = n; }
  /// Retrieve the max. number of clusters on a track.
  int GetClustersMaximum() const { return m_maxclusters; }

  /// Set the W value [eV].
  void SetWorkFunction(const double w) { m_work = w; }
  /// Get the W value [eV].
  double GetWorkFunction() const { return m_work; }
  /// Set the Fano factor.
  void SetFanoFactor(const double f) { 
    m_fano = f; 
    m_fset = true;
  }
  /// Use the default Fano factor.
  void UnsetFanoFactor() { m_fset = false; } 
  /// Get the Fano factor.
  double GetFanoFactor() const { return m_fano; }
  /// Set the density [g/cm3] of the target medium.
  void SetDensity(const double density) { m_rho = density; }
  /// Get the density [g/cm3] of the target medium.
  double GetDensity() const { return m_rho; }
  /// Set A and Z of the target medium.
  void SetAtomicMassNumbers(const double a, const double z) {
    m_a = a;
    m_z = z;
  }
  /// Get A and Z of the target medium.
  void GetAtomicMassMumbers(double& a, double& z) const {
    a = m_a;
    z = m_z;
  }

  /// Simulate transverse straggling (default: on).
  void EnableTransverseStraggling(const bool on = true) { 
    m_useTransStraggle = on; 
  }
  /// Simulate longitudinal straggling (default: off).
  void EnableLongitudinalStraggling(const bool on = true) { 
    m_useLongStraggle = on; 
  }

  bool NewTrack(const double x0, const double y0, const double z0,
                const double t0, const double dx0, const double dy0,
                const double dz0) override;
  bool GetCluster(double& xcls, double& ycls, double& zcls,
                  double& tcls, int& n, double& e, double& extra) override;

 protected:
  /// Include transverse straggling
  bool m_useTransStraggle = true;
  /// Include longitudinal straggling
  bool m_useLongStraggle = false;

  /// Has the charge been defined?
  bool m_chargeset = false;
  /// Charge of the projectile
  double m_qion = 1.;
  /// Mass [MeV] of the projectile
  double m_mion = -1.;
  /// Mass density [g/cm3] of the target
  double m_rho = -1.;
  /// Work function [eV] of the target
  double m_work = -1.;
  /// Has the Fano factor been set?
  bool m_fset = false;
  /// Fano factor [-] of the target
  double m_fano = -1.;
  /// Effective A of the target
  double m_a = -1.;
  /// Effective Z of the target 
  double m_z = -1.;

  /// Maximum number of clusters allowed (infinite if 0)
  int m_maxclusters = -1;
  /// Energy in energy loss table [MeV]
  std::vector<double> m_ekin;
  /// EM energy loss [MeV cm2/g]
  std::vector<double> m_emloss;
  /// Hadronic energy loss [MeV cm2/g]
  std::vector<double> m_hdloss;
  /// Projected range [cm]
  std::vector<double> m_range;
  /// Transverse straggling [cm]
  std::vector<double> m_transstraggle;
  /// Longitudinal straggling [cm]
  std::vector<double> m_longstraggle;

  /// Index of the next cluster to be returned
  size_t m_currcluster = 0;
  /// Fluctuation model (0 = none, 1 = Landau, 2 = Vavilov,
  ///                    3 = Gaussian, 4 = Combined)
  unsigned int m_model = 4;
  /// Targeted cluster size
  int m_nsize = -1;
  struct Cluster {
    double x, y, z, t;  ///< Cluster location and time
    double ec;          ///< Energy spent to make the cluster
    double kinetic;     ///< Ion energy when cluster was created
    int electrons;      ///< Number of electrons in this cluster
  };
  std::vector<Cluster> m_clusters;

  double Xi(const double x, const double beta2, const double edens) const;
  double DedxEM(const double e) const;
  double DedxHD(const double e) const;
  bool PreciseLoss(const double step, const double estart, double& deem,
                   double& dehd) const;
  bool EstimateRange(const double ekin, const double step, 
                     double& stpmax) const;
  bool SmallestStep(const double ekin, const double edens,
                    double de, double step, double& stpmin);
  Medium* GetMedium(const std::array<double, 3>& x) const;
  double Terminate(const std::array<double, 3>& x0,
                   const std::array<double, 3>& v0, const double step0) const;
  double TerminateBfield(const std::array<double, 3>& x0,
                         const std::array<double, 3>& v0,
                         const double dt0, const double vmag) const;
  double RndmEnergyLoss(const double ekin, const double de,
                        const double step, const double edens) const;
};
}

#endif
