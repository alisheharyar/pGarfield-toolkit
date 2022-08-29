#include "heed++/code/PhotoAbsCSLib.h"
#include "heed++/code/PhysicalConstants.h"

#include <iostream>

// 2004, I. Smirnov

namespace {

std::string getDataBasePath() { 

  std::string path;
  // First try if the environment variable HEED_DATABASE is defined.
  char* heed_database = std::getenv("HEED_DATABASE");
  if (heed_database) {
    path = std::string(heed_database);
  } else {
    // If HEED_DATABASE is not defined, try GARFIELD_INSTALL.
    char* garfield_install = std::getenv("GARFIELD_INSTALL");
    if (garfield_install) {
      path = std::string(garfield_install) + "/share/Heed/database/";
    } else {
      // Try GARFIELD_HOME.
      char* garfield_home = std::getenv("GARFIELD_HOME");
      if (garfield_home) {
        path = std::string(garfield_home) + "/Heed/heed++/database/";
      }
    }
  }
  return path;
}

Heed::ExAtomPhotoAbsCS generate_Ar_PACS(const std::string& shelllist_dir,
                                        const std::string& pacs_table_dir) {

  Heed::ExAtomPhotoAbsCS Argon_PACS_mod_esc(18, shelllist_dir + "shelllist.dat",
                                            pacs_table_dir + "Ar.dat");

  // ExAtomPhotoAbsCS Argon_PACS_mod_esc(18,
  //                                     shelllist_dir + "shelllist.dat",
  //                                     shelllist_dir + "mw3.dat");

  // ExAtomPhotoAbsCS Argon_PACS_mod_esc(18, "argon",
  //                                     shelllist_dir + "ftbf18.dat", 2);

  Heed::AtomicSecondaryProducts* asp = Argon_PACS_mod_esc.get_asp(1);
  std::vector<double> electron_energy;
  std::vector<double> photon_energy;
  electron_energy.emplace_back(0.000200);
  // electron_energy.emplace_back(0.002670);
  asp->add_channel(0.65, electron_energy, photon_energy);
  electron_energy.resize(2);
  electron_energy[0] = 0.000050;
  electron_energy[1] = 0.000200;
  asp->add_channel(0.35, electron_energy, photon_energy, 1);
  // mcout<<"L1:\n";
  // asp->print(mcout, 2);

  asp = Argon_PACS_mod_esc.get_asp(2);
  electron_energy.resize(1);
  electron_energy[0] = 0.000200;
  asp->add_channel(1.0, electron_energy, photon_energy, 1);
  // mcout<<"L2:\n";
  // asp->print(mcout, 2);

  asp = Argon_PACS_mod_esc.get_asp(3);
  electron_energy.resize(1);
  electron_energy[0] = 0.000200;
  asp->add_channel(1.0, electron_energy, photon_energy, 1);
  // mcout<<"L3:\n";
  // asp->print(mcout, 2);

  return Argon_PACS_mod_esc;
}
}

namespace Heed {

using CLHEP::gram;
using CLHEP::mole;


std::map<std::string, ExAtomPhotoAbsCS> PhotoAbsCSLib::apacs;
std::map<std::string, SimpleAtomPhotoAbsCS> PhotoAbsCSLib::hpacs;

AtomPhotoAbsCS* PhotoAbsCSLib::getAPACS(const std::string& name) {
  if (apacs.empty()) initialise();
  if (name == "H" || name.find("H for") == 0) {
    return &hpacs[name];
  }
  if (apacs.count(name) > 0) return &apacs[name];
  return nullptr;
}

void PhotoAbsCSLib::initialise() {
  if (!hpacs.empty() || !apacs.empty()) return;
  // Hydrogen.
  hpacs.emplace("H", 
      SimpleAtomPhotoAbsCS(1, std::make_shared<HydrogenPhotoAbsCS>()));
  hpacs.emplace("H for H2", 
      SimpleAtomPhotoAbsCS(1, std::make_shared<PhenoPhotoAbsCS>("Hydrogen_for_H2", 1, 15.43e-6, 3.228)));
  hpacs.emplace("H for CH4", 
      SimpleAtomPhotoAbsCS(1, std::make_shared<PhenoPhotoAbsCS>("Hydrogen_for_CH4", 1, 12.65e-06, 3.228)));
  hpacs.emplace("H for NH4",
      SimpleAtomPhotoAbsCS(1, std::make_shared<PhenoPhotoAbsCS>("Hydrogen_for_NH4", 1, 10.0e-06, 3.228))); 
  // hpacs.emplace("H for CH4", 
  //     SimpleTablePhotoAbsCS("Hydrogen_for_CH4", 1, 12.65e-6,
  //                           dbpath + "H_for_CH4.dat");

  std::string dbpath = getDataBasePath();
  if (dbpath.empty()) {
    std::cerr << "Heed::PhotoAbsCsLib::initialise:\n"
              << "    Could not retrieve database path.\n";
  }
  if (dbpath[dbpath.size() - 1] != '/') dbpath.append("/");

  const std::string pacs_table_dir = dbpath + "henke/";

  // Argon.
  apacs.emplace("Ar", generate_Ar_PACS(dbpath, pacs_table_dir));

  // Other atoms.
  std::string shells = dbpath + "shelllist.dat"; 
  const std::map<std::string, int> atoms = {
    {"He",  2}, {"Li",  3}, {"Be",  4}, {"B",   5}, {"C",   6}, 
    {"O",   8}, {"F",   9}, {"Ne", 10}, {"Na", 11}, {"Mg", 12},
    {"Al", 13}, {"Si", 14}, {"P",  15}, {"S",  16}, {"Cl", 17}, 
    {"Ga", 31}, {"Ge", 32}, {"As", 33}, {"Br", 35}, {"Kr", 36}, 
    {"Cd", 48}, {"Te", 52}, {"Xe", 54}, {"Cs", 55}, {"Hg", 80}, 
    {"U",  92}};
  for (const auto& atom : atoms) {
    std::string pacstable = pacs_table_dir + atom.first + ".dat";
    apacs.emplace(atom.first, ExAtomPhotoAbsCS(atom.second, shells, pacstable));
  }
  // "Standard" Argon:
  // apacs.emplace("Ar", ExAtomPhotoAbsCS(18, shells,
  //                                      pacs_table_dir + "Ar.dat"));
  // Optional variants:
  // apacs.emplace("Ar", ExAtomPhotoAbsCS(18, shells, dbpath + "mw3.dat"));
  // Variant for debug, pointwise cross section
  // apacs.emplace("Ar", ExAtomPhotoAbsCS(18, "argon",
  //                                      dbpath + "ftbf18.dat", 2));
  // Variant for debug, fitted cross section
  // apacs.emplace("Ar", ExAtomPhotoAbsCS(18, "argon", dbpath + "shelltscf.dat",
  //                                      2, 0, 0.0);
  // Variant for debug, fitted cross section with replacement from Henke
  // apacs.emplace("Ar", ExAtomPhotoAbsCS(18, "argon", dbpath + "shelltscf.dat",
  //                                      pacs_table_dir + "Ar.dat",
  //                                      40.0e-6, 2, 0.0);
  // Another variant for debug, fitted cross section with replacement from
  // Marr and West, should be similar to old Fortran verion
  // apacs.emplace("Ar", ExAtomPhotoAbsCS(18, "argon", dbpath + "shelltscf.dat",
  //                                      dbpath + "mw3.dat", 40.0e-6, 2, 0.0);

  // For debug, FitBT
  // apacs.emplace("C", ExAtomPhotoAbsCS(6, "carbon", dbpath + "shelltscf.dat",
  //                                     2, 0, 0.0);
  // apacs.emplace("C for CH4", ExAtomPhotoAbsCS(6, shells,
  //                                             pacs_table_dir + "C.dat",
  //                                             "C_for_CH4", 12.65e-06);

  apacs.emplace("C for CH4", ExAtomPhotoAbsCS(6, shells,
                                              dbpath + "C_for_CH4.dat",
                                              "C_for_CH4", 12.65e-6));
  apacs.emplace("C for C2H4", ExAtomPhotoAbsCS(6, shells,
                                               pacs_table_dir + "C.dat",
                                               "C_for_C2H4", 10.51e-06));
  apacs.emplace("C for C2H6", ExAtomPhotoAbsCS(6, shells,
                                               pacs_table_dir + "C.dat",
                                               "C_for_C2H6", 11.52e-06));
  apacs.emplace("C for C4H10", ExAtomPhotoAbsCS(6, shells,
                                                pacs_table_dir + "C.dat",
                                                "C_for_C4H10", 10.55e-06));
  apacs.emplace("C for Methylal", ExAtomPhotoAbsCS(6, shells,
                                                   pacs_table_dir + "C.dat",
                                                   "C_for_Methylal", 10.0e-06));
  apacs.emplace("C for CF4", ExAtomPhotoAbsCS(6, shells,
                                              pacs_table_dir + "C.dat", 
                                              "C_for_CF4", 16.23e-06));
  apacs.emplace("C for CO2", ExAtomPhotoAbsCS(6, shells,
                                              pacs_table_dir + "C.dat", 
                                              "C_for_CO2", 13.79e-06));
  apacs.emplace("N", ExAtomPhotoAbsCS(7, shells, pacs_table_dir + "N.dat", 
                                      "N_for_N2", 15.581e-6));
  apacs.emplace("O for CO2", ExAtomPhotoAbsCS(8, shells,
                                              pacs_table_dir + "O.dat", 
                                              "O_for_CO2", 13.79e-6));

  std::string sshells = dbpath + "shelllist_solid.dat";

  apacs.emplace("Diamond", ExAtomPhotoAbsCS(6, sshells,
                                            pacs_table_dir + "C.dat", 
                                            "Diamond"));
  apacs.emplace("Si crystal", ExAtomPhotoAbsCS(14, sshells,
                                               pacs_table_dir + "Si.dat",
                                               "Si_crystal"));
  apacs.emplace("Ge crystal", ExAtomPhotoAbsCS(32, shells,
                                               pacs_table_dir + "Ge.dat",
                                               "Ge_crystal", 0.67e-06));
  apacs.emplace("Si G4", ExAtomPhotoAbsCS(14, sshells,                                                                     dbpath + "Si_G4.dat", "Si_G4"));
  apacs.emplace("Ga for GaAs", ExAtomPhotoAbsCS(31, sshells,
                                                pacs_table_dir + "Ga.dat",
                                                "Ga_for_GaAs"));
  apacs.emplace("As for GaAs", ExAtomPhotoAbsCS(33, sshells,
                                                pacs_table_dir + "As.dat",
                                                "As_for_GaAs"));
  apacs.emplace("Cd for CdTe", ExAtomPhotoAbsCS(48, sshells,
                                                pacs_table_dir + "Cd.dat",
                                                "Cd_for_CdTe"));
  apacs.emplace("Te for CdTe", ExAtomPhotoAbsCS(52, sshells,
                                                pacs_table_dir + "Te.dat",
                                                "Te_for_CdTe"));
}

}
