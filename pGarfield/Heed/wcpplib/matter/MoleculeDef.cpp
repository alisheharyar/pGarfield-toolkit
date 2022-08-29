#include <iomanip>
#include "wcpplib/matter/MoleculeDef.h"
#include "wcpplib/util/FunNameStack.h"
#include "wcpplib/clhep_units/WPhysicalConstants.h"
#include "wcpplib/math/cubic.h"

// 1998-2004 I. Smirnov

namespace Heed {

using CLHEP::k_Boltzmann;
using CLHEP::Avogadro;
using CLHEP::cm3;
using CLHEP::gram;
using CLHEP::mole;
using CLHEP::bar;
using CLHEP::hep_pascal;
using CLHEP::kelvin;

VanDerWaals::VanDerWaals(double fPk, double fTk) : Pkh(fPk), Tkh(fTk) {
  // Rydberg constant
  const double R = k_Boltzmann * Avogadro;

  Vkh = R * 3.0 / 8.0 * Tkh / Pkh;
  ah = 3 * Pkh * Vkh * Vkh;
  bh = 1.0 / 3.0 * Vkh;
}

double VanDerWaals::volume_of_mole(double T, double p, int& s_not_single) {
  mfunname("VanDerWaals::volume_of_mole(...)");

  double Tr = T / Tkh;
  double Pr = p / Pkh;
  Iprint2n(mcout, Tr, Pr);
  Cubic cb(Pr, -1.0 / 3.0 * (Pr + 8 * Tr), 3, -1);
  double r[3];
  int q = cb.find_real_zero(r);
  check_econd11(q, <= 0, mcerr);
  double x = r[q - 1];   // this is the relative volume taken by one mole
  double res = x * Vkh;  // this is the absolute volume taken by one mole
  Iprint2n(mcout, x, res);
  s_not_single = q == 2 ? 1 : 0;
  return res;
}

VanDerWaals* VanDerWaals::copy() const { return new VanDerWaals(*this); }

std::ostream& operator<<(std::ostream& file, const VanDerWaals& f) {
  mfunname(
      "std::ostream& operator << (std::ostream& file, const VanDerWaals& f)");
  Ifile << "VanDerWaals:\n";
  indn.n += 2;
  Iprintn(file, f.Pk() / (CLHEP::atmosphere));
  Iprintn(file, f.Tk() / (CLHEP::kelvin));
  Iprintn(file, f.Vk() / (cm3));
  Ifile << "For comparison, the volume of a mole of ideal gas\n";
  Ifile << "at the same conditions takes\n";
  Iprintn(file, (k_Boltzmann * Avogadro * f.Tk() / f.Pk()) / (cm3 * mole));
  Iprintn(file, f.a() / (CLHEP::atmosphere * cm3 * cm3));
  Iprintn(file, f.b() / (cm3));
  indn.n -= 2;
  return file;
}

MoleculeDef::MoleculeDef(const std::string& fname, const std::string& fnotation,
                         long fqatom, const std::vector<std::string>& fatom_not,
                         const std::vector<long>& fqatom_ps,
                         std::shared_ptr<VanDerWaals> fvdw)
    : AtomMixDef(fqatom, fatom_not, fqatom_ps),
      nameh(fname),
      notationh(fnotation),
      qatom_psh(fqatom_ps) {
  mfunname("MoleculeDef::MoleculeDef(...)");
  m_vdw = std::move(fvdw);
  for (long n = 0; n < qatom(); n++) {
    Z_totalh += qatom_psh[n] * atom(n)->Z();
    A_totalh += qatom_psh[n] * atom(n)->A();
    tqatomh += qatom_psh[n];
    check_econd11(qatom_psh[n], <= 0, mcerr);
  }
}

// one atom in molecule
MoleculeDef::MoleculeDef(const std::string& fname, const std::string& fnotation,
                         const std::string& fatom_not, long fqatom_ps,
                         std::shared_ptr<VanDerWaals> fvdw) :
    MoleculeDef(fname, fnotation, 1, {fatom_not}, {fqatom_ps}, fvdw) { 
}

// two atoms
MoleculeDef::MoleculeDef(const std::string& fname, const std::string& fnotation,
                         const std::string& fatom_not1, long fqatom_ps1,
                         const std::string& fatom_not2, long fqatom_ps2,
                         std::shared_ptr<VanDerWaals> fvdw) : 
    MoleculeDef(fname, fnotation, 2, {fatom_not1, fatom_not2}, 
                {fqatom_ps1, fqatom_ps2}, fvdw) {

}

// three atoms
MoleculeDef::MoleculeDef(const std::string& fname, const std::string& fnotation,
                         const std::string& fatom_not1, long fqatom_ps1,
                         const std::string& fatom_not2, long fqatom_ps2,
                         const std::string& fatom_not3, long fqatom_ps3,
                         std::shared_ptr<VanDerWaals> fvdw) :
    MoleculeDef(fname, fnotation, 3, {fatom_not1, fatom_not2, fatom_not3}, 
                {fqatom_ps1, fqatom_ps2, fqatom_ps3}, fvdw) {

}

void MoleculeDef::print(std::ostream& file, int l) const {
  if (l > 0) file << (*this);
}

std::ostream& operator<<(std::ostream& file, const MoleculeDef& f) {
  mfunnamep("std::ostream& operator << (std::ostream&, const MoleculeDef&)");
  constexpr double gpm = gram / mole;
  Ifile << "MoleculeDef: name=" << std::setw(10) << f.name()
        << " notation=" << std::setw(3) << f.notation() << '\n';
  indn.n += 2;
  Ifile << "Z_total()=" << std::setw(3) << f.Z_total()
        << " A_total()/(gram/mole)=" << f.A_total() / gpm 
        << " tqatom()=" << f.tqatom() << '\n';
  Iprintn(file, f.qatom());
  indn.n += 2;
  for (long n = 0; n < f.qatom(); n++) {
    Ifile << "n=" << n << " atom(n)->notation=" << f.atom(n)->notation()
          << " qatom_ps(n)=" << f.qatom_ps(n) << '\n';
  }
  indn.n -= 2;
  f.AtomMixDef::print(file, 1);
  VanDerWaals* at = f.vdw().get();
  if (at) {
    Ifile << "Density at the crucial conditions for ideal gas (for debug):\n";
    double rydberg = k_Boltzmann * Avogadro;  // more precise
    // mcout<<"rydberg/(joule/(kelvin*mole)) ="
    //     << rydberg/(joule/(kelvin*mole))<<'\n';
    // double sa = f.A_total();
    Iprintn(mcout,
            f.A_total() * at->Pk() / (rydberg * at->Tk()) / (gram / cm3));
    Ifile << "For the Waals:\n";
    Iprintn(mcout, f.A_total() / at->Vk() / (gram / cm3));
  }
  indn.n -= 2;
  return file;
}

std::list<MoleculeDef> MoleculeDefs::molecules;

void MoleculeDefs::printMolecules(std::ostream& file) {
  Ifile << "MoleculeDefs::printMolecules:\n";
  for (const auto& molecule : getMolecules()) {
    file << molecule;
  }
}

const MoleculeDef* MoleculeDefs::getMolecule(const std::string& fnotation) {
  for (const auto& molecule : getMolecules()) {
    if (molecule.notation() == fnotation) return &molecule;
  }
  return nullptr;
}

const std::list<MoleculeDef>& MoleculeDefs::getMolecules() {
  if (!molecules.empty()) return molecules;
  molecules.emplace_back(MoleculeDef("Hydrogen", "H2", "H", 2));
  molecules.emplace_back(MoleculeDef("Helium", "He", "He", 1));
  molecules.emplace_back(MoleculeDef("Nitrogen", "N2", "N", 2));
  molecules.emplace_back(MoleculeDef("Oxygen", "O2", "O", 2));
  molecules.emplace_back(MoleculeDef("Neon", "Ne", "Ne", 1));
  // molecules.emplace_back(MoleculeDef("Argon_without_K", "Ar_without_K",
  //                                    "Ar_without_K", 1));
  molecules.emplace_back(MoleculeDef("Argon", "Ar", "Ar", 1,
      std::make_shared<VanDerWaals>(48.6 * bar, 150.7 * kelvin)));
  molecules.emplace_back(MoleculeDef("Krypton", "Kr", "Kr", 1,
      std::make_shared<VanDerWaals>(55.0 * bar, 209.4 * kelvin)));
  molecules.emplace_back(MoleculeDef("Xenon", "Xe", "Xe", 1,
      std::make_shared<VanDerWaals>(55.0 * bar, 209.4 * kelvin)));

  molecules.emplace_back(MoleculeDef("NH3", "NH3", "N", 1, "H", 3));
  molecules.emplace_back(MoleculeDef("N2O", "N2O", "N", 2, "O", 1));
  molecules.emplace_back(MoleculeDef("CO2", "CO2", "C", 1, "O", 2));
  molecules.emplace_back(MoleculeDef("CH4", "CH4", "C", 1, "H", 4,
      std::make_shared<VanDerWaals>(4.64e6 * hep_pascal, 
                                    (273.15 - 82.5) * kelvin)));
  molecules.emplace_back(MoleculeDef("CF4", "CF4", "C", 1, "F", 4,
      std::make_shared<VanDerWaals>(42.5 * bar, 369.8 * kelvin)));
  molecules.emplace_back(MoleculeDef("SF4", "SF4", "S", 1, "F", 4));
  molecules.emplace_back(MoleculeDef("SF6", "SF6", "S", 1, "F", 6));
  molecules.emplace_back(MoleculeDef("C2H2", "C2H2", "C", 2, "H", 2));
  molecules.emplace_back(MoleculeDef("C2H4", "C2H4", "C", 2, "H", 4));
  molecules.emplace_back(MoleculeDef("C2H6", "C2H6", "C", 2, "H", 6));
  molecules.emplace_back(MoleculeDef("C3H8", "C3H8", "C", 3, "H", 8,
      std::make_shared<VanDerWaals>(42.5 * bar, 369.8 * kelvin)));
  molecules.emplace_back(MoleculeDef("C4H10", "C4H10", "C", 4, "H", 10,
      std::make_shared<VanDerWaals>(40.0 * bar, 418.3 * kelvin)));
  molecules.emplace_back(MoleculeDef("C2H2F4", "C2H2F4", "C", 2, "F", 4, "H", 2));
  molecules.emplace_back(MoleculeDef("H2O", "H2O", "H", 2, "O", 1,
      std::make_shared<VanDerWaals>(22.9e6 * hep_pascal, (273.15 + 374.15) * kelvin)));
  molecules.emplace_back(MoleculeDef("Methylal", "Methylal", "O", 2, "C", 3, "H", 8,
      std::make_shared<VanDerWaals>(39.5 * bar, 480.6 * kelvin)));

  // Additional molecule definitions for compatibility with Magboltz
  molecules.emplace_back(MoleculeDef("C5H12", "C5H12", "C", 5, "H", 12));
  molecules.emplace_back(MoleculeDef("NO", "NO", "N", 1, "O", 1));
  molecules.emplace_back(MoleculeDef("CO", "CO", "C", 1, "O", 1));
  molecules.emplace_back(MoleculeDef("DME", "DME", "C", 2, "H", 6, "O", 1));
  molecules.emplace_back(MoleculeDef("C2F6", "C2F6", "C", 2, "F", 6));
  molecules.emplace_back(MoleculeDef("C3H6", "C3H6", "C", 3, "H", 6));
  molecules.emplace_back(MoleculeDef("CH3OH", "CH3OH", "C", 1, "H", 4, "O", 1));
  molecules.emplace_back(MoleculeDef("C2H5OH", "C2H5OH", "C", 2, "H", 6, "O", 1));
  molecules.emplace_back(MoleculeDef("C3H7OH", "C3H7OH", "C", 3, "H", 8, "O", 1));
  molecules.emplace_back(MoleculeDef("Cs", "Cs", "Cs", 1));
  molecules.emplace_back(MoleculeDef("F2", "F2", "F", 2));
  molecules.emplace_back(MoleculeDef("CS2", "CS2", "C", 1, "S", 2));
  molecules.emplace_back(MoleculeDef("COS", "COS", "C", 1, "O", 1, "S", 1));
  molecules.emplace_back(MoleculeDef("BF3", "BF3", "B", 1, "F", 3));
  molecules.emplace_back(MoleculeDef("C2HF5", "C2HF5", "C", 2, "H", 1, "F", 5));
  molecules.emplace_back(MoleculeDef("CHF3", "CHF3", "C", 1, "H", 1, "F", 3));
  molecules.emplace_back(MoleculeDef("CF3Br", "CF3Br", "C", 1, "F", 3, "Br", 1));
  molecules.emplace_back(MoleculeDef("C3F8", "C3F8", "C", 3, "F", 8));
  molecules.emplace_back(MoleculeDef("O3", "O3", "O", 3));
  molecules.emplace_back(MoleculeDef("Hg", "Hg", "Hg", 1));
  molecules.emplace_back(MoleculeDef("H2S", "H2S", "H", 2, "S", 1));
  molecules.emplace_back(MoleculeDef("GeH4", "GeH4", "Ge", 1, "H", 4));
  molecules.emplace_back(MoleculeDef("SiH4", "SiH4", "Si", 1, "H", 4));
  molecules.emplace_back(MoleculeDef("CCl4", "CCl4", "C", 1, "Cl", 4));
  return molecules;
}

}
