#ifndef G_SHAPER_H
#define G_SHAPER_H

#include <string>
#include <cmath>

namespace Garfield {

/// Class for signal processing

class Shaper {
 public:
  /// Default constructor.
  Shaper() = delete;
  /// Constructor.
  Shaper(const unsigned int n, const double tau, const double g, 
         std::string shaperType);
  /// Destructor.
  ~Shaper() {}
  
  /// Evaluate the transfer function.
  double Shape(const double t) const;
  /// Transfer function for a unipolar shaper.
  double UnipolarShaper(const double t) const;
  /// Transfer function for a bipolar shaper.
  double BipolarShaper(const double t) const;
  /// Time for the transfer function to rise from zero to peak height.
  double PeakingTime() const { return m_tp; }

  /// Return the integral of the transfer function squared.
  double TransferFuncSq() const { return m_transfer_func_sq; }

  /// Is it a unipolar shaper?
  bool IsUnipolar() const { return (m_type == ShaperType::Unipolar); } 
  /// Is it a bipolar shaper?
  bool IsBipolar() const { return (m_type == ShaperType::Bipolar); }
  /// Retrieve the parameters.
  void GetParameters(unsigned int& n, double& tp) {
    n = m_n;
    tp = m_tp;
  }
 
 private:
  std::string m_className = "Shaper";

  // Shaper type.
  enum class ShaperType { Unipolar = 0, Bipolar };
  ShaperType m_type = ShaperType::Unipolar;
  // Order of the shaper.
  unsigned int m_n = 1;
  // Time constant.
  double m_tau = 1.;
  // Peaking time.
  double m_tp = 1.;
  // Normalization factor.
  double m_prefactor = 1.;
  // Gain.
  double m_g = 1.;
  // Integral of the transfer function squared.
  double m_transfer_func_sq = -1.;
};
}

#endif
