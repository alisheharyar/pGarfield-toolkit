#ifndef G_MEDIUM_GAS_H
#define G_MEDIUM_GAS_H

#include <array>
#include <cmath>
#include <vector>
#include <bitset>

#include "Medium.hh"

namespace Garfield {

/// Base class for gas media.

class MediumGas : public Medium {
 public:
  /// Constructor
  MediumGas();
  /// Destructor
  virtual ~MediumGas() {}

  /// Set the gas mixture.
  bool SetComposition(const std::string& gas1, const double f1 = 1.,
                      const std::string& gas2 = "", const double f2 = 0.,
                      const std::string& gas3 = "", const double f3 = 0.,
                      const std::string& gas4 = "", const double f4 = 0.,
                      const std::string& gas5 = "", const double f5 = 0.,
                      const std::string& gas6 = "", const double f6 = 0.);
  /// Retrieve the gas mixture.
  void GetComposition(std::string& gas1, double& f1, 
                      std::string& gas2, double& f2, 
                      std::string& gas3, double& f3,
                      std::string& gas4, double& f4, 
                      std::string& gas5, double& f5, 
                      std::string& gas6, double& f6) const;

  /// Read table of gas properties (transport parameters) from file.
  bool LoadGasFile(const std::string& filename, const bool quiet = false);
  /// Save the present table of gas properties (transport parameters) to a file.
  bool WriteGasFile(const std::string& filename);
  /// Read table of gas properties from and merge with the existing dataset.
  bool MergeGasFile(const std::string& filename, const bool replaceOld);

  /// Switch on simulation of Penning transfers, using pre-implemented 
  /// parameterisations of the transfer probability (if available).
  virtual bool EnablePenningTransfer();
  /** Switch on simulation of Penning transfers by means of
    * transfer probabilities, for all excitation levels in the mixture.
    * \param r transfer probability [0, 1]
    * \param lambda parameter for sampling the distance of the Penning electron
             with respect to the excitation.
    */
  virtual bool EnablePenningTransfer(const double r, const double lambda);
  /// Switch on simulation of Penning transfers by means of
  /// transfer probabilities, for all excitations of a given component.
  virtual bool EnablePenningTransfer(const double r, const double lambda,
                                     std::string gasname);
  /// Switch the simulation of Penning transfers off globally.
  virtual void DisablePenningTransfer();
  /// Switch the simulation of Penning transfers off for a given component.
  virtual bool DisablePenningTransfer(std::string gasname);

  /// Print information about the present gas mixture and available data.
  virtual void PrintGas();
  /// Print a list of all available gases.
  static void PrintGases();

  /// Read a table of (positive) ion mobilities vs. electric field from file.
  bool LoadIonMobility(const std::string& filename, const bool quiet = false);
  /// Read a table of negative ion mobilities vs. electric field from file.
  bool LoadNegativeIonMobility(const std::string& filename, 
                               const bool quiet = false);

  /// Adjust the Townsend coefficient using the excitation and ionisation 
  /// rates stored in the gas table and the Penning transfer probabilities.
  bool AdjustTownsendCoefficient();

  /// Return the number of ionisation levels in the table.
  size_t GetNumberOfIonisationLevels() const { return m_ionLevels.size(); }
  /// Return the number of excitation levels in the table.
  size_t GetNumberOfExcitationLevels() const { return m_excLevels.size(); }
  /// Return the identifier and threshold of an ionisation level.
  void GetIonisationLevel(const size_t level, std::string& label, 
                          double& energy) const; 
  /// Return the identifier and energy of an excitation level.
  void GetExcitationLevel(const size_t level, std::string& label, 
                          double& energy) const; 
  /// Get an entry in the table of ionisation rates.
  bool GetElectronIonisationRate(const size_t level,
                                 const size_t ie, const size_t ib,
                                 const size_t ia, double& f) const;
  /// Get an entry in the table of excitation rates.
  bool GetElectronExcitationRate(const size_t level,
                                 const size_t ie, const size_t ib,
                                 const size_t ia, double& f) const;

  bool IsGas() const override { return true; }

  void GetComponent(const unsigned int i, std::string& label,
                    double& f) override;

  void SetAtomicNumber(const double z) override;
  double GetAtomicNumber() const override;
  void SetAtomicWeight(const double a) override;
  double GetAtomicWeight() const override;
  void SetNumberDensity(const double n) override;
  double GetNumberDensity() const override;
  void SetMassDensity(const double rho) override;
  double GetMassDensity() const override;

  void ResetTables() override;

  void SetExtrapolationMethodExcitationRates(const std::string& low,
                                             const std::string& high) {
    SetExtrapolationMethod(low, high, m_extrExc, "ExcitationRates");
  }
  void SetExtrapolationMethodIonisationRates(const std::string& low,
                                             const std::string& high) {
    SetExtrapolationMethod(low, high, m_extrIon, "IonisationRates");
  }
  void SetInterpolationMethodExcitationRates(const unsigned int intrp) {
    if (intrp > 0) m_intpExc = intrp;
  }
  void SetInterpolationMethodIonisationRates(const unsigned int intrp) {
    if (intrp > 0) m_intpIon = intrp;
  }

  // Scaling laws.
  // TODO: cache scaling factors.
  double ScaleElectricField(const double e) const override {
    return e * m_pressureTable / m_pressure;
  }
  double UnScaleElectricField(const double e) const override {
    return e * m_pressure / m_pressureTable;
  }
  double ScaleDiffusion(const double d) const override {
    return d * sqrt(m_pressureTable / m_pressure);
  }
  double ScaleDiffusionTensor(const double d) const override {
    return d * m_pressureTable / m_pressure;
  }
  double ScaleTownsend(const double alpha) const override {
    return alpha * m_pressure / m_pressureTable;
  }
  double ScaleAttachment(const double eta) const override {
    return eta * m_pressure / m_pressureTable;
  }
  double ScaleLorentzAngle(const double lor) const override {
    return lor * m_pressure / m_pressureTable;
  }

  bool GetPhotoAbsorptionCrossSection(const double e, double& sigma,
                                      const unsigned int i) override;

 protected:
  static constexpr unsigned int m_nMaxGases = 6;

  // Gas mixture
  std::array<std::string, m_nMaxGases> m_gas;
  std::array<double, m_nMaxGases> m_fraction;
  std::array<double, m_nMaxGases> m_atWeight;
  std::array<double, m_nMaxGases> m_atNum;

  // Penning transfer
  // Flag enabling/disabling Penning transfer
  bool m_usePenning = false;
  // Penning transfer probability
  double m_rPenningGlobal = 0.;
  // Mean distance of Penning ionisation
  double m_lambdaPenningGlobal = 0.;
  // Penning transfer probability per component
  std::array<double, m_nMaxGases> m_rPenningGas;
  // Penning transfer distance per component
  std::array<double, m_nMaxGases> m_lambdaPenningGas;

  // Pressure at which the transport parameter table was calculated
  double m_pressureTable;
  // Temperature at which the transport parameter table was calculated
  double m_temperatureTable;

  // Table of Townsend coefficients without Penning transfer
  std::vector<std::vector<std::vector<double> > > m_eAlp0;

  // Tables for excitation and ionisation rates
  std::vector<std::vector<std::vector<std::vector<double> > > > m_excRates;
  std::vector<std::vector<std::vector<std::vector<double> > > > m_ionRates;

  // Store excitation and ionization information
  struct ExcLevel {
    std::string label;
    double energy;
    double prob;
    double rms;
    double dt;
  };
  std::vector<ExcLevel> m_excLevels;

  struct IonLevel {
    std::string label;
    double energy;
  };
  std::vector<IonLevel> m_ionLevels;

  // Extrapolation/interpolation for excitation and ionisation rates.
  std::pair<unsigned int, unsigned int> m_extrExc = {0, 1};
  std::pair<unsigned int, unsigned int> m_extrIon = {0, 1};
  unsigned int m_intpExc = 2;
  unsigned int m_intpIon = 2;

  bool LoadMobility(const std::string& filename, const bool quiet,
                    const bool negative);
  bool ReadHeader(std::ifstream& gasfile, int& version,
                  std::bitset<20>& gasok, bool& is3d, 
                  std::vector<double>& mixture,
                  std::vector<double>& efields, std::vector<double>& bfields,
                  std::vector<double>& angles, std::vector<ExcLevel>& excLevels,
                  std::vector<IonLevel>& ionLevels);
  void ReadFooter(std::ifstream& gasfile,
                  std::array<unsigned int, 13>& extrapH,
                  std::array<unsigned int, 13>& extrapL,
                  std::array<unsigned int, 13>& interp, 
                  unsigned int& thrAlp, unsigned int& thrAtt, 
                  unsigned int& thrDis, 
                  double& ionDiffL, double& ionDiffT,
                  double& pgas, double& tgas);
  void ReadRecord3D(std::ifstream& gasfile, double& ve, double& vb, double& vx,
                    double& dl, double& dt, double& alpha, double& alpha0, 
                    double& eta, double& mu, double& lor,
                    double& dis, std::array<double, 6>& dif, 
                    std::vector<double>& rexc, std::vector<double>& rion);
  void ReadRecord1D(std::ifstream& gasfile, double& ve, double& vb, double& vx,
                    double& dl, double& dt, double& alpha, double& alpha0, 
                    double& eta, double& mu, double& lor,
                    double& dis, std::array<double, 6>& dif, 
                    std::vector<double>& rexc, std::vector<double>& rion);
  void InsertE(const int ie, const int ne, const int nb, const int na);
  void InsertB(const int ib, const int ne, const int nb, const int na);
  void InsertA(const int ia, const int ne, const int nb, const int na);
  void ZeroRowE(const int ie, const int nb, const int na);
  void ZeroRowB(const int ib, const int ne, const int na);
  void ZeroRowA(const int ia, const int ne, const int nb);
  bool GetMixture(const std::vector<double>& mixture, const int version,
                  std::vector<std::string>& gasnames,
                  std::vector<double>& percentages) const;
  void GetGasBits(std::bitset<20>& gasok) const;
 
  static bool GetGasInfo(const std::string& gasname, 
                         double& a, double& z, double& w, double& f);
  static std::string GetGasName(const int gasnumber, const int version);
  static std::string GetGasName(std::string input);
  static int GetGasNumberGasFile(const std::string& input);
  static const std::vector<std::string> GetAliases(const std::string& gas);

};
}

#endif
