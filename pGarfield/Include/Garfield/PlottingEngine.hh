#ifndef G_PLOTTING_ENGINE_H
#define G_PLOTTING_ENGINE_H

#include <TStyle.h>

namespace Garfield {

/// Plotting style.

class PlottingEngine {
 public:
  /// Default constructor.
  PlottingEngine();
  /// Destructor
  ~PlottingEngine() = default;

  /// Use serif font.
  void SetSerif() { m_serif = true; }
  /// Use sans-serif font.
  void SetSansSerif() { m_serif = false; }

  /// Set the colour palette.
  void SetPalette(int ncol) { m_palette = ncol; }

  /// Apply the default Garfield ROOT style.
  void SetDefaultStyle();

 private:
  std::string m_className = "PlottingEngine";

  bool m_serif = false;
  int m_palette = 0;

  TStyle m_garfieldStyle;
};
}

#endif
