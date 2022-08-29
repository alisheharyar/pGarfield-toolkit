#ifndef G_TGEOTET
#define G_TGEOTET

#include <array>

#include "TGeoBBox.h"

class TGeoTet : public TGeoBBox {

public:
  TGeoTet() {}
  /// Constructor
  TGeoTet(const char* name,
          const std::array<std::array<double, 3>, 4>& vertices);
  /// Destructor
  virtual ~TGeoTet() {}

  void ComputeBBox() override;
  int DistancetoPrimitive(int, int) override { return 99999; }
  const TBuffer3D &GetBuffer3D(int reqSections, bool localFrame) const override;
  void GetMeshNumbers(int& nvert, int& nsegs, int& npols) const override {
    nvert =  4;
    nsegs = 12;
    npols =  4;
  }
  int GetNmeshVertices() const override { return 4; }
  void InspectShape() const override {}
  TBuffer3D* MakeBuffer3D() const override;
  void Print(Option_t *option = "") const override;
  void SavePrimitive(std::ostream &, Option_t *) override {}
  void SetPoints(double* points) const override;
  void SetPoints(float* points) const override;
  void SetSegsAndPols(TBuffer3D &buff) const override;
  void Sizeof3D() const override {}

  // ClassDef(TGeoTet, 1)

 private:
  std::array<std::array<double, 3>, 4> fVertices;

  TGeoTet(const TGeoTet&) = delete;
  TGeoTet& operator=(const TGeoTet&) = delete;

};

#endif
