#include <iomanip>
#include "wcpplib/matter/MatterDef.h"
#include "wcpplib/util/FunNameStack.h"
#include "wcpplib/clhep_units/WPhysicalConstants.h"

// 1998-2004 I. Smirnov

namespace Heed {

void MatterDef::calc_I_eff() { I_effh = Z_mean() * 12.0 * CLHEP::eV; }

MatterDef::MatterDef(const std::string& fname, const std::string& fnotation,
                     long fqatom, const std::vector<std::string>& fatom_not,
                     const std::vector<double>& fweight_quan, double fdensity,
                     double ftemperature)
    : AtomMixDef(fqatom, fatom_not, fweight_quan),
      nameh(fname),
      notationh(fnotation),
      temperatureh(ftemperature),
      densityh(fdensity) {
  mfunname("MatterDef::MatterDef(...)");
  calc_I_eff();
}

MatterDef::MatterDef(const std::string& fname, const std::string& fnotation,
                     const std::string& fatom_not, double fdensity,
                     double ftemperature) :
    MatterDef(fname, fnotation, 1, {fatom_not}, {1.}, fdensity, ftemperature) {

}

MatterDef::MatterDef(const std::string& fname, const std::string& fnotation,
                     const std::string& fatom_not1, double fweight_quan1,
                     const std::string& fatom_not2, double fweight_quan2,
                     double fdensity, double ftemperature) :
    MatterDef(fname, fnotation, 2, {fatom_not1, fatom_not2},
              {fweight_quan1, fweight_quan2}, fdensity, ftemperature) {

}

MatterDef::MatterDef(const std::string& fname, const std::string& fnotation,
                     const std::string& fatom_not1, double fweight_quan1,
                     const std::string& fatom_not2, double fweight_quan2,
                     const std::string& fatom_not3, double fweight_quan3,
                     double fdensity, double ftemperature) :
    MatterDef(fname, fnotation, 3, {fatom_not1, fatom_not2, fatom_not3},
              {fweight_quan1, fweight_quan2, fweight_quan3}, 
              fdensity, ftemperature) {

}

void MatterDef::print(std::ostream& file, int l) const {
  if (l > 0) file << (*this);
}

std::ostream& operator<<(std::ostream& file, const MatterDef& f) {
  mfunname("ostream& operator << (ostream& file, const MatterDef& f)");
  Ifile << "MatterDef: name=" << std::setw(10) << f.name()
        << " notation=" << std::setw(3) << f.notation() << '\n';
  indn.n += 2;
  Ifile << "density/(gram/cm3)=" << f.density() / (CLHEP::gram / CLHEP::cm3)
        << " temperature/kelvin=" << f.temperature() / CLHEP::kelvin
        << " I_eff/eV=" << f.I_eff() / CLHEP::eV << '\n';
  f.AtomMixDef::print(file, 1);
  indn.n -= 2;
  return file;
}

}
