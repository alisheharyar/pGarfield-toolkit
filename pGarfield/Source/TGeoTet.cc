#include <iostream>

#include "Garfield/TGeoTet.hh"

#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoVolume.h"
#include "TVirtualGeoPainter.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"

// ClassImp(TGeoTet)

using Vec = std::array<double, 3>;

namespace {

Vec Dir(const Vec& a, const Vec& b) {
  Vec c;
  for (size_t i = 0; i < 3; ++i) c[i] = b[i] - a[i];
  return c;
}

Vec Cross(const Vec& a, const Vec& b) {

  Vec c = {a[1] * b[2] - a[2] * b[1],
           a[2] * b[0] - a[0] * b[2],
           a[0] * b[1] - a[1] - b[0]};
  return c;
}

double Dot(const Vec& a, const Vec& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

}

TGeoTet::TGeoTet(const char *name,
                 const std::array<Vec, 4>& vertices) 
    : TGeoBBox(name, 0, 0, 0) {
  fVertices = vertices;

  Vec u;
  Vec v;
  Vec w;
  for (size_t i = 0; i < 3; ++i) {
    u[i] = fVertices[1][i] - fVertices[0][i];
    v[i] = fVertices[2][i] - fVertices[0][i];
    w[i] = fVertices[3][i] - fVertices[0][i];
  }
  Vec x = Cross(v, w);
  double det = Dot(x, u);
  if (det < 0.) std::swap(fVertices[0], fVertices[1]);
}

void TGeoTet::ComputeBBox() {
  const double kBig = TGeoShape::Big();
  double vmin[3] = {kBig, kBig, kBig};
  double vmax[3] = {-kBig, -kBig, -kBig};
  for (size_t i = 0; i < 4; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      vmin[j] = std::min(vmin[j], fVertices[i][j]);
      vmax[j] = std::max(vmax[j], fVertices[i][j]);
    }
  }
  fDX = 0.5 * (vmax[0] - vmin[0]);
  fDY = 0.5 * (vmax[1] - vmin[1]);
  fDZ = 0.5 * (vmax[2] - vmin[2]);
  for (size_t i = 0; i < 3; ++i) { 
    fOrigin[i] = 0.5 * (vmax[i] + vmin[i]);
  }
}

TBuffer3D *TGeoTet::MakeBuffer3D() const {
  // Number of vertices.
  constexpr int nv = 4;
  // Number of segments.
  // constexpr int ns = 6;
  constexpr int ns = 12;
  // Number of polygons.
  constexpr int np = 4;
  auto buff = new TBuffer3D(TBuffer3DTypes::kGeneric, nv, 3 * nv, ns, 3 * ns, np, 5 * np);
  if (buff) {
    SetPoints(buff->fPnts);
    SetSegsAndPols(*buff);
  }
  return buff;
}

void TGeoTet::Print(Option_t *) const {
  std::cout << "=== Tetrahedron " << GetName() << "\n";
}

void TGeoTet::SetSegsAndPols(TBuffer3D &buff) const {
  const int c = GetBasicColor();

  auto v01 = Dir(fVertices[0], fVertices[1]);
  auto v02 = Dir(fVertices[0], fVertices[2]);
  auto v03 = Dir(fVertices[0], fVertices[3]);

  auto v12 = Dir(fVertices[1], fVertices[2]);
  auto v13 = Dir(fVertices[1], fVertices[3]);
  auto v10 = Dir(fVertices[1], fVertices[0]); 

  auto v23 = Dir(fVertices[2], fVertices[3]);
  auto v20 = Dir(fVertices[2], fVertices[0]);
  auto v21 = Dir(fVertices[2], fVertices[1]);

  auto v30 = Dir(fVertices[3], fVertices[0]);
  auto v31 = Dir(fVertices[3], fVertices[1]);
  auto v32 = Dir(fVertices[3], fVertices[2]);

  std::vector<std::array<int, 3> > faces;
  if (Dot(Cross(v01, v02), v03) > 0) {
    faces.push_back({1, 0, 2});
  } else {
    faces.push_back({0, 1, 2});
  }

  if (Dot(Cross(v12, v13), v10) > 0) {
    faces.push_back({2, 1, 3});
  } else {
    faces.push_back({1, 2, 3});
  }

  if (Dot(Cross(v23, v20), v21) > 0) {
    faces.push_back({3, 2, 0});
  } else {
    faces.push_back({2, 3, 0});
  }

  if (Dot(Cross(v30, v31), v32) > 0) {
    faces.push_back({0, 3, 1});
  } else {
    faces.push_back({3, 0, 1});
  }

  size_t ind = 0;
  for (const auto& face : faces) {
    buff.fSegs[ind++] = c; 
    buff.fSegs[ind++] = face[0];
    buff.fSegs[ind++] = face[1];
    buff.fSegs[ind++] = c; 
    buff.fSegs[ind++] = face[1];
    buff.fSegs[ind++] = face[2];
    buff.fSegs[ind++] = c; 
    buff.fSegs[ind++] = face[2];
    buff.fSegs[ind++] = face[0];
  }

  ind = 0;
  for (size_t i = 0; i < 4; ++i) {
    buff.fPols[ind++] = c;
    buff.fPols[ind++] = 3;
    buff.fPols[ind++] = 3 * i;
    buff.fPols[ind++] = 3 * i + 1;
    buff.fPols[ind++] = 3 * i + 2;
  }
}

void TGeoTet::SetPoints(double *points) const {
  size_t ind = 0;
  for (const auto& vertex : fVertices) {
    points[ind++] = vertex[0];
    points[ind++] = vertex[1];
    points[ind++] = vertex[2];
  }
}

void TGeoTet::SetPoints(float *points) const {
  size_t ind = 0;
  for (const auto& vertex : fVertices) {
    points[ind++] = vertex[0];
    points[ind++] = vertex[1];
    points[ind++] = vertex[2];
  }
}

/// Fills a static 3D buffer and returns a reference.
const TBuffer3D& TGeoTet::GetBuffer3D(int reqSections, bool localFrame) const {
  static TBuffer3D buffer(TBuffer3DTypes::kGeneric);

  FillBuffer3D(buffer, reqSections, localFrame);

  constexpr int nv = 4;
  // constexpr int ns = 6;
  constexpr int ns = 12;
  constexpr int np = 4;
  if (reqSections & TBuffer3D::kRawSizes) {
    if (buffer.SetRawSizes(nv, 3 * nv, ns, 3 * ns, np, 5 * np)) {
      buffer.SetSectionsValid(TBuffer3D::kRawSizes);
    }
  }
  if ((reqSections & TBuffer3D::kRaw) && buffer.SectionsValid(TBuffer3D::kRawSizes)) {
    SetPoints(buffer.fPnts);
    if (!buffer.fLocalFrame) {
      TransformPoints(buffer.fPnts, buffer.NbPnts());
    }
    SetSegsAndPols(buffer);
    buffer.SetSectionsValid(TBuffer3D::kRaw);
  }
  return buffer;
}

