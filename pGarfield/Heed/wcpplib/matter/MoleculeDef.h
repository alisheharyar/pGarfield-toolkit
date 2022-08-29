#ifndef MOLECULE_DEF_H
#define MOLECULE_DEF_H

#include <memory>
#include "wcpplib/matter/AtomDef.h"

namespace Heed {

/// Helper class for Van-der-Waals equation.
class VanDerWaals {
 private:
  double ah;
  double bh;
  double Vkh;
  double Pkh;
  double Tkh;

 public:
  VanDerWaals(double fPk, double fTk);
  double a() const { return ah; }
  double b() const { return bh; }
  double Vk() const { return Vkh; }
  double Pk() const { return Pkh; }
  double Tk() const { return Tkh; }
  /*
  double pressure(double M, // the number of moles
                  double volume,
                  double T);
  double volume(double T,  // relative to T_k
                double p,  // relative to p_k
                int &s_not_single);
  */
  // Return number of moles in the unit volume
  double volume_of_mole(double T, double p, int& s_not_single);

  VanDerWaals* copy() const;
};
std::ostream& operator<<(std::ostream& file, const VanDerWaals& f);

/// Definition of molecule as a mixture of atoms.
/// Only the basic information: the name, the notation,
/// the mean charge and atomic weight and the parameters of mixture class.
///
/// 1998-2004 I. Smirnov

class MoleculeDef : public AtomMixDef {
  std::string nameh = "none";
  std::string notationh = "none";
  /// Number of atoms of particular sort in the molecule.
  /// Obviously it is not normalized to one, but instead
  /// the sum is equal to tqatomh
  std::vector<long> qatom_psh;
  long Z_totalh = 0;
  double A_totalh = 0.;
  /// Total number of atoms in molecule
  /// Attention: this is not the number of different sorts of atoms
  /// The latter is qatom() from AtomMixDef
  long tqatomh = 0;
  std::shared_ptr<VanDerWaals> m_vdw;

 public:
  const std::string& name() const { return nameh; }
  const std::string& notation() const { return notationh; }
  const std::vector<long>& qatom_ps() const { return qatom_psh; }
  long qatom_ps(long n) const { return qatom_psh[n]; }
  long Z_total() const { return Z_totalh; }
  double A_total() const { return A_totalh; }
  long tqatom() const { return tqatomh; }
  const std::shared_ptr<VanDerWaals>& vdw() const { return m_vdw; }
  MoleculeDef() = default;
  MoleculeDef(const std::string& fname, const std::string& fnotation,
              long fqatom, const std::vector<std::string>& fatom_not,
              const std::vector<long>& fqatom_ps,
              std::shared_ptr<VanDerWaals> fvdw = {});
  MoleculeDef(const std::string& fname, const std::string& fnotation,
              const std::string& fatom_not, long fqatom_ps,
              std::shared_ptr<VanDerWaals> fvdw = {});
  MoleculeDef(const std::string& fname, const std::string& fnotation,
              const std::string& fatom_not1, long fqatom_ps1,
              const std::string& fatom_not2, long fqatom_ps2,
              std::shared_ptr<VanDerWaals> fvdw = {});
  MoleculeDef(const std::string& fname, const std::string& fnotation,
              const std::string& fatom_not1, long fqatom_ps1,
              const std::string& fatom_not2, long fqatom_ps2,
              const std::string& fatom_not3, long fqatom_ps3,
              std::shared_ptr<VanDerWaals> fvdw = {});
  ~MoleculeDef() = default;

  void print(std::ostream& file, int l) const;

  MoleculeDef* copy() const { return new MoleculeDef(*this); }
};
std::ostream& operator<<(std::ostream& file, const MoleculeDef& f);

/// Library of molecules.

class MoleculeDefs {

public:
  static const std::list<MoleculeDef>& getMolecules();
  /// Return the address of the molecule with this name.
  /// If there is no molecule with this notation, the function returns NULL
  /// but does not terminate the program as that for AtomDef. Be careful.
  static const MoleculeDef* getMolecule(const std::string& fnotation);
  /// Print all registered molecules.
  static void printMolecules(std::ostream& file);
private:
  static std::list<MoleculeDef> molecules;
};

}

#endif
