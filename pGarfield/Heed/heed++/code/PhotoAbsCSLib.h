#ifndef PHOTOABSCSLIB_H
#define PHOTOABSCSLIB_H

#include <map>

#include "heed++/code/PhotoAbsCS.h"

namespace Heed {

/// Library of photoabsorption cross sections for some frequently used atoms
/// and molecules.
/// 2004, I. Smirnov

class PhotoAbsCSLib {

 public:
  static AtomPhotoAbsCS* getAPACS(const std::string& name); 
 private:
  static std::map<std::string, ExAtomPhotoAbsCS> apacs;
  static std::map<std::string, SimpleAtomPhotoAbsCS> hpacs;
  static void initialise();
};

}

#endif
