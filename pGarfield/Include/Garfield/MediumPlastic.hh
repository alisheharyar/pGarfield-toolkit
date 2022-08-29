#ifndef G_MEDIUM_PLASTIC_H
#define G_MEDIUM_PLASTIC_H

#include "Medium.hh"

namespace Garfield {

/// Plastic medium.

class MediumPlastic : public Medium {
 public:
  // Constructor
  MediumPlastic() : Medium() {
    m_className = "MediumPlastic";
    m_name = "Plastic";
  }
  // Destructor
  virtual ~MediumPlastic() {}

  void EnableDrift(const bool /*on*/) override {}
  void EnablePrimaryIonisation(const bool /*on*/) override {}
};
}

#endif
