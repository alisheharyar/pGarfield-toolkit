#include "wcpplib/particle/particle_def.h"
#include "wcpplib/clhep_units/WPhysicalConstants.h"
#include "wcpplib/util/FunNameStack.h"

// 1998 - 2004,   I. Smirnov

namespace Heed {

using CLHEP::electron_mass_c2;
using CLHEP::proton_mass_c2;
using CLHEP::neutron_mass_c2;
using CLHEP::c_squared;
using CLHEP::electron_charge;
using CLHEP::eplus;
using CLHEP::MeV;
using CLHEP::GeV;

particle_def electron_def("electron", "e-", electron_mass_c2 / c_squared,
                          electron_charge, 0.5);
particle_def positron_def("positron", "e+", electron_def);

particle_def muon_minus_def("muon_minus", "mu-", 105.658367 * MeV / c_squared,
                            electron_charge, 0.5);
particle_def muon_plus_def("muon_plus", "mu+", muon_minus_def);

particle_def proton_def("proton", "p+", proton_mass_c2 / c_squared, eplus, 0.5);
particle_def anti_proton_def("", "p-", proton_def);

// light unflavored mesons
particle_def pi_plus_meson_def("pi_plus", "pi+",
                               139.56755 * MeV / c_squared, eplus, 0.0);
particle_def pi_minus_meson_def("pi_minus", "pi-",
                                139.56755 * MeV / c_squared, -eplus, 0.0);
particle_def K_plus_meson_def("K_plus", "K+",
                              493.677 * MeV / c_squared, 1, 0.0);
particle_def K_minus_meson_def("K_minus", "K-", K_plus_meson_def);

particle_def deuteron_def("deuteron", "d", 1875.613 * MeV / c_squared, eplus,
                          0.0);
particle_def alpha_particle_def("alpha_particle", "alpha",
                                3727.417 * MeV / c_squared, 2 * eplus, 0.);

particle_def::particle_def(const std::string& fname,
                           const std::string& fnotation, double fmass,
                           double fcharge, float fspin) : 
    name(fname), notation(fnotation), mass(fmass), 
    charge(fcharge), spin(fspin) {
}

particle_def::particle_def(const std::string& fname,
                           const std::string& fnotation, particle_def& p) {
  // creates anti-particle through the call of anti_particle(p)
  *this = anti_particle(p);
  if (!(fname == "" || fname == " ")) name = fname;
  if (!(fnotation == "" || fnotation == " ")) notation = fnotation;
}

particle_def particle_def::anti_particle(const particle_def& p) {
  std::string aname = "anti-" + p.name;
  std::string anot = "anti-" + p.notation;
  return particle_def(aname, anot, p.mass, -p.charge, -p.spin);
}

void particle_def::set_mass(const double m) { mass = m * MeV / c_squared; }

void particle_def::set_charge(const double z) { charge = z * eplus; }

void particle_def::print(std::ostream& file, int l) const {
  if (l > 0) file << (*this);
}

std::ostream& operator<<(std::ostream& file, const particle_def& f) {
  Ifile << "particle_def: name=" << f.name << " notation=" << f.notation
        << '\n';
  Ifile << "mass=" << f.mass
        << " mass/(GeV/c_squared)=" << f.mass / (GeV / c_squared)
        << " charge=" << f.charge << " charge/eplus=" << f.charge / eplus
        << '\n';
  Ifile << "spin=" << f.spin << '\n';
  return file;
}

}
