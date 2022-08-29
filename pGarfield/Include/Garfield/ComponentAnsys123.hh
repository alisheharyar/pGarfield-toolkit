#ifndef G_COMPONENT_ANSYS123_H
#define G_COMPONENT_ANSYS123_H

#include "ComponentFieldMap.hh"

namespace Garfield {

/// Component for importing and interpolating three-dimensional ANSYS field maps.

class ComponentAnsys123 : public ComponentFieldMap {
 public:
  /// Constructor
  ComponentAnsys123();
  /// Destructor
  ~ComponentAnsys123() {}

  bool Initialise(std::string elist = "ELIST.lis",
                  std::string nlist = "NLIST.lis",
                  std::string mplist = "MPLIST.lis",
                  std::string prnsol = "PRNSOL.lis", std::string unit = "cm");

  bool SetWeightingField(std::string prnsol, std::string label);

};
}
#endif
