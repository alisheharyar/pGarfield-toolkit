#include "wcpplib/particle/eparticle.h"

// 1998 - 2004,   I. Smirnov

namespace Heed {

eparticle::eparticle(manip_absvol* primvol, const point& pt, const vec& vel,
                     vfloat ftime, particle_def* fpardef, HeedFieldMap* fieldmap)
    : mparticle(primvol, pt, vel, ftime, fpardef->mass), 
      m_pardef(fpardef), m_fieldMap(fieldmap) {
}

int eparticle::force(const point& pt, vec& f, vec& f_perp, vfloat& mrange) {
  vec efield(0., 0., 0.);
  vec hfield(0., 0., 0.);
  if (!m_fieldMap) {
    std::cerr << "Field map not defined.\n";
    return 1;
  }
  m_fieldMap->field_map(pt, efield, hfield, mrange);
  f = m_pardef->charge * efield;
  f_perp = m_pardef->charge * hfield;
  return 1;
}

void eparticle::print(std::ostream& file, int l) const {
  if (l < 0) return;
  Ifile << "eparticle: particle is ";
  if (!m_pardef) {
    file << "none";
  } else {
    file << m_pardef->notation;
  }
  file << '\n';
  mparticle::print(file, l);
}
}
