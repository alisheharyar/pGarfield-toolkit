#ifndef ATOM_DEF_H
#define ATOM_DEF_H

#include <iostream>
#include <vector>
#include <list>

namespace Heed {

/// Definition of atoms.
/// Only the basic information: name, notation, atomic weight and charge.
///
/// In principle I am going to initiate all atoms from Mendeleev's table,
/// but I haven't finished yet. Only its first half is filled at the moment.
///
/// 1998-2004, I. Smirnov.

class AtomDef {
  std::string nameh = "none";
  std::string notationh = "none";
  /// Atomic number.
  int Zh = 0;
  /// Atomic mass in internal units. Transfer to gram/mole if need.
  double Ah = 0.;

 public:
  /// Default constructor
  AtomDef() = default;
  /// Constructor
  AtomDef(const std::string& fnameh, const std::string& fnotationh, int fZh,
          double fAh);
  /// Destructor
  ~AtomDef() = default;

  const std::string& name() const { return nameh; }
  const std::string& notation() const { return notationh; }
  int Z() const { return Zh; }
  double A() const { return Ah; }

  void print(std::ostream& file, int l = 0) const;
  AtomDef* copy() const { return new AtomDef(*this); } 
};
std::ostream& operator<<(std::ostream& file, const AtomDef& f);

/// Library of atoms.

class AtomDefs {

 public:
  static void addAtom(const std::string& name, const std::string& notation,
                      const int z, const double a);
  static const std::list<AtomDef>& getAtoms();

  /// Return the address of atom with this name if it is registered in system,
  /// or NULL otherwise
  static const AtomDef* getAtom(const std::string& fnotation);
  /// Return the atomic number corresponding to a given Z.
  /// If the atom is not registered, the program is terminated. Be careful!
  static double getA(int fZ);
  /// Return the address of atom corresponding to a given Z.
  /// If the atom is not registered, the program is terminated. Be careful!
  static const AtomDef* getAtom(int fZ);

  /// Print all registered atoms.
  static void printAtoms(std::ostream& file);
 private:
  static std::list<AtomDef> atoms;
};

/// Definition of atomic mixtures. Pointers to atoms, weights and
/// various mean parameters.

class AtomMixDef {
  /// Number of different atoms.
  long qatomh = 0;
  /// Constituent atoms.
  std::vector<const AtomDef*> atomh;
  std::vector<double> weight_quanh;  // sum is 1
  std::vector<double> weight_massh;  // sum is 1

  /// Weighted mean Z
  double Z_meanh = 0.;
  /// Weighted mean A (in internal units). Transfer to gram/mole if needed.
  double A_meanh = 0.;
  /// Weighted mean 1 / A (in internal units).
  double inv_A_meanh = 0.;
  /// Weighted mean ratio Z / A.
  double mean_ratio_Z_to_Ah = 0.;
  double NumberOfElectronsInGramh = 0.;

 public:
  /// Default constructor
  AtomMixDef() = default;
  /// Constructor from list of atoms and weights.
  AtomMixDef(unsigned long fqatom, const std::vector<std::string>& fatom_not,
             const std::vector<double>& fweight_quan);
  /// Constructor from list of atoms and number of atoms per molecule.
  AtomMixDef(unsigned long fqatom, const std::vector<std::string>& fatom_not,
             const std::vector<long>& fweight_quan);
  void print(std::ostream& file, int l) const;
  long qatom() const { return qatomh; }
  const std::vector<const AtomDef*>& atom() const { return atomh; }
  const AtomDef* atom(long n) const { return atomh[n]; }
  const std::vector<double>& weight_quan() const { return weight_quanh; }
  const std::vector<double>& weight_mass() const { return weight_massh; }
  double weight_quan(long n) const { return weight_quanh[n]; }
  double weight_mass(long n) const { return weight_massh[n]; }
  double Z_mean() const { return Z_meanh; }
  double A_mean() const { return A_meanh; }
  double inv_A_mean() const { return inv_A_meanh; }
  double mean_ratio_Z_to_A() const { return mean_ratio_Z_to_Ah; }
  double NumberOfElectronsInGram() const {
    return NumberOfElectronsInGramh;
  }
};
std::ostream& operator<<(std::ostream& file, const AtomMixDef& f);
}

#endif
