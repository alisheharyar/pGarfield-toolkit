#ifndef G_GEOMETRY_ROOT_H
#define G_GEOMETRY_ROOT_H

#include <map>

#include <TGeoManager.h>
#include <TGeoMaterial.h>

#include "Geometry.hh"

namespace Garfield {

/// Use a geometry defined using the ROOT TGeo package.

class GeometryRoot : public Geometry {
 public:
  /// Constructor
  GeometryRoot();
  /// Destructor
  ~GeometryRoot() {}

  /// Set the geometry (pointer to ROOT TGeoManager).
  void SetGeometry(TGeoManager* geoman);

  Medium* GetMedium(const double x, const double y, const double z,
                    const bool tesselated = false) const override;

  /// Get the number of materials defined in the ROOT geometry.
  unsigned int GetNumberOfMaterials();
  /// Get a pointer to the ROOT material with a given index.
  TGeoMaterial* GetMaterial(const unsigned int i);
  /// Get a pointer to the ROOT material with a given name.
  TGeoMaterial* GetMaterial(const char* name);
  /// Associate a ROOT material with a Garfield medium.
  void SetMedium(const unsigned int imat, Medium* med);
  /// Associate a ROOT material with a Garfield medium.
  void SetMedium(const char* mat, Medium* med);

  bool IsInside(const double x, const double y, const double z,
                const bool tesselated = false) const override;
  bool GetBoundingBox(double& xmin, double& ymin, double& zmin, double& xmax,
                      double& ymax, double& zmax) override;

  /// Switch debugging and warning messages on/off.
  void EnableDebugging(const bool on = true) { m_debug = on; }

 protected:
  // ROOT geometry manager
  TGeoManager* m_geoManager = nullptr;

  // List of ROOT materials associated to Garfield media
  std::map<std::string, Medium*> m_materials;

  // Switch on/off debugging messages.
  bool m_debug = false;
  void PrintGeoNotDefined(const std::string& fcn) const;
};
}

#endif
