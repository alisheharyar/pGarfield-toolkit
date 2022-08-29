#ifndef HEEDPARTICLE_BGM_H
#define HEEDPARTICLE_BGM_H

#include <vector>
#include "HeedCluster.h"
#include "wcpplib/particle/eparticle.h"

namespace Heed {

/// Definition of the particle which can be traced through the geometry.
/// 2003, I. Smirnov

class HeedParticle_BGM : public eparticle {
 public:
  /// Default constructor.
  HeedParticle_BGM() : eparticle() {}
  /// Constructor.
  /// if fs_loss_only == true - only transfer energy and
  /// no other physics: no deposition of clusters,
  /// no generation of virtual photons.
  /// Thus it is just a PAI without even clusters
  HeedParticle_BGM(manip_absvol* primvol, const point& pt, const vec& vel,
                   vfloat time, particle_def* fpardef, HeedFieldMap* fieldmap,
                   bool fs_loss_only = false, bool fs_print_listing = false);
  /// Destructor
  virtual ~HeedParticle_BGM() {}

  void print(std::ostream& file, int l) const override;
  HeedParticle_BGM* copy() const override { 
    return new HeedParticle_BGM(*this); 
  }

 protected:
  void physics(std::vector<gparticle*>& secondaries) override;

 private:
  bool m_print_listing = false;
  bool m_loss_only = false;

  long m_particle_number = 0;

  double m_edep = 0.;

  std::vector<HeedCluster> m_clusterBank;
};
}

#endif
