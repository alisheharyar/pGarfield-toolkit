#include <cmath>
#include <iostream>

#include <TGeoBBox.h>
#include <TGeoCone.h>
#include <TGeoArb8.h>
#include <TGeoXtru.h>
#include <TGeoBoolNode.h>
#include <TGeoCompositeShape.h>
#include <TPolyLine.h>

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/Solid.hh"
#include "Garfield/ViewGeometry.hh"

namespace Garfield {

ViewGeometry::ViewGeometry() : ViewBase("ViewGeometry") {} 

ViewGeometry::~ViewGeometry() {
  Reset();
}

void ViewGeometry::SetGeometry(GeometrySimple* geo) {
  if (!geo) {
    std::cerr << m_className << "::SetGeometry: Null pointer.\n";
    return;
  }

  m_geometry = geo;
}

void ViewGeometry::Plot(const bool twod) {

  if (twod) {
    Plot2d();
  } else {
    Plot3d();
  }
}

void ViewGeometry::Plot3d() {
  if (!m_geometry) {
    std::cerr << m_className << "::Plot3d: Geometry is not defined.\n";
    return;
  }

  const auto nSolids = m_geometry->GetNumberOfSolids();
  if (nSolids == 0) {
    std::cerr << m_className << "::Plot3d: Geometry is empty.\n";
    return;
  }

  // Get the bounding box.
  double xMin = 0., yMin = 0., zMin = 0.;
  double xMax = 0., yMax = 0., zMax = 0.;
  if (!m_geometry->GetBoundingBox(xMin, yMin, zMin, xMax, yMax, zMax)) {
    std::cerr << m_className << "::Plot3d: Cannot retrieve bounding box.\n";
    return;
  }
  auto pad = GetCanvas();
  pad->cd();
  gGeoManager = nullptr;
  m_geoManager.reset(new TGeoManager("ViewGeometryGeoManager", ""));
  TGeoMaterial* matVacuum = new TGeoMaterial("Vacuum", 0., 0., 0.);
  TGeoMedium* medVacuum = new TGeoMedium("Vacuum", 1, matVacuum);
  m_media.push_back(medVacuum);
  // Use silicon as "default" material.
  TGeoMaterial* matDefault = new TGeoMaterial("Default", 28.085, 14., 2.329);
  TGeoMedium* medDefault = new TGeoMedium("Default", 1, matDefault);
  TGeoVolume* world = m_geoManager->MakeBox(
      "World", medVacuum, std::max(fabs(xMin), fabs(xMax)),
      std::max(fabs(yMin), fabs(yMax)), std::max(fabs(zMin), fabs(zMax)));
  m_geoManager->SetTopVolume(world);
  m_volumes.push_back(world);

  for (size_t i = 0; i < nSolids; ++i) {
    auto solid = m_geometry->GetSolid(i);
    if (!solid) {
      std::cerr << m_className << "::Plot3d:\n"
                << "    Could not get solid " << i << " from geometry.\n";
      continue;
    }
    // Get the center coordinates.
    double x0 = 0., y0 = 0., z0 = 0.;
    if (!solid->GetCentre(x0, y0, z0)) {
      std::cerr << m_className << "::Plot3d: Could not determine solid centre.\n";
      continue;
    }
    // Get the rotation.
    double ctheta = 1., stheta = 0.;
    double cphi = 1., sphi = 0.;
    if (!solid->GetOrientation(ctheta, stheta, cphi, sphi)) {
      std::cerr << m_className << "::Plot3d:\n"
                << "    Could not determine solid orientation.\n";
      continue;
    }
    double matrix[9] = {cphi * ctheta, -sphi, cphi * stheta,
                        sphi * ctheta, cphi,  sphi * stheta,
                        -stheta,       0,     ctheta};
    TGeoVolume* volume = nullptr;
    if (solid->IsTube()) {
      const double rt = solid->GetRadius();
      const double lz = solid->GetHalfLengthZ();
      volume = m_geoManager->MakeTube("Tube", medDefault, 0., rt, lz);
    } else if (solid->IsWire()) {
      const double rw = solid->GetRadius();
      const double lz = solid->GetHalfLengthZ();
      volume = m_geoManager->MakeTube("Wire", medDefault, 0., rw, lz);
    } else if (solid->IsBox()) {
      const double dx = solid->GetHalfLengthX();
      const double dy = solid->GetHalfLengthY();
      const double dz = solid->GetHalfLengthZ();
      volume = m_geoManager->MakeBox("Box", medDefault, dx, dy, dz);
    } else if (solid->IsSphere()) {
      const double rmin = std::max(solid->GetInnerRadius(), 0.);
      const double rmax = solid->GetOuterRadius();
      volume = m_geoManager->MakeSphere("Sphere", medDefault, rmin, rmax);
    } else if (solid->IsHole()) {
      const double r1 = solid->GetLowerRadius();
      const double r2 = solid->GetUpperRadius();
      const double dx = solid->GetHalfLengthX();
      const double dy = solid->GetHalfLengthY();
      const double dz = solid->GetHalfLengthZ();
      const double lz = 10 * std::max({dx, dy, dz});
      const double rm = 0.5 * (r1 + r2);
      const double dr = 0.5 * (r2 - r1) * lz / dz;
      TGeoBBox* box = new TGeoBBox("HoleBox", dx, dy, dz);
      TGeoCone* cone = new TGeoCone("HoleCone", lz, 0, rm - dr, 0, rm + dr);
      TGeoCompositeShape* hole = new TGeoCompositeShape("Hole", 
        new TGeoSubtraction(box, cone));
      hole->RegisterYourself();
      volume = new TGeoVolume("Hole", hole, medDefault); 
    } else if (solid->IsRidge()) {
      const double dx = solid->GetHalfLengthX();
      const double dy = solid->GetHalfLengthY();
      const double dz = 0.5 * solid->GetRidgeHeight();
      const double xr = solid->GetRidgeOffset();
      volume = m_geoManager->MakeArb8("Ridge", medDefault, dz);
      auto arb = (TGeoArb8*)volume->GetShape();
      arb->SetVertex(0, -dx, -dy);
      arb->SetVertex(1, -dx, +dy);
      arb->SetVertex(2, +dx, +dy);
      arb->SetVertex(3, +dx, -dy);
      arb->SetVertex(4, xr, -dy);
      arb->SetVertex(5, xr, +dy);
      arb->SetVertex(6, xr, +dy);
      arb->SetVertex(7, xr, -dy);
      z0 += dz;
    } else if (solid->IsExtrusion()) {
      const double dz = solid->GetHalfLengthZ();
      std::vector<double> xp;
      std::vector<double> yp;
      if (!solid->GetProfile(xp, yp)) continue;
      volume = m_geoManager->MakeXtru("Extrusion", medDefault, 2);
      auto xtru = (TGeoXtru*)volume->GetShape();
      xtru->DefinePolygon(xp.size(), xp.data(), yp.data());
      xtru->DefineSection(0, -dz);
      xtru->DefineSection(1, +dz);
    } else {
      std::cerr << m_className << "::Plot3d: Unknown type of solid.\n";
      continue;
    }
    Medium* medium = m_geometry->GetMedium(x0, y0, z0);
    if (solid->GetColour() >= 0) {
      volume->SetLineColor(solid->GetColour());
    } else if (!medium) {
      volume->SetLineColor(kGreen + 2);
      volume->SetTransparency(50);
    } else if (medium->IsGas()) {
      volume->SetLineColor(kBlue + medium->GetId());
      volume->SetTransparency(50);
    } else if (medium->IsSemiconductor()) {
      volume->SetLineColor(kRed + medium->GetId());
      volume->SetTransparency(50);
    } else {
      volume->SetLineColor(kViolet + medium->GetId());
      volume->SetTransparency(0);
    }
    TGeoRotation r;
    r.SetMatrix(matrix);
    TGeoTranslation t(x0, y0, z0);
    TGeoCombiTrans* transform = new TGeoCombiTrans(t, r);
    m_volumes.push_back(volume);
    m_geoManager->GetTopVolume()->AddNode(volume, 1, transform);
  }
  m_geoManager->CloseGeometry();
  m_geoManager->GetTopNode()->Draw("gl");
}

void ViewGeometry::Plot2d() {

  if (!m_geometry) {
    std::cerr << m_className << "::Plot2d: Geometry is not defined.\n";
    return;
  }

  const auto nSolids = m_geometry->GetNumberOfSolids();
  if (nSolids == 0) {
    std::cerr << m_className << "::Plot2d: Geometry is empty.\n";
    return;
  }

  // Determine the plot limits.
  double x0 = 0., y0 = 0.;
  double x1 = 0., y1 = 0.;
  if (m_userPlotLimits) {
    x0 = m_xMinPlot;
    y0 = m_yMinPlot;
    x1 = m_xMaxPlot;
    y1 = m_yMaxPlot;
  } else if (m_userBox) {
    PlotLimitsFromUserBox(x0, y0, x1, y1);
  } else {
    std::array<double, 3> bbmin;
    std::array<double, 3> bbmax;
    if (!m_geometry->GetBoundingBox(bbmin[0], bbmin[1], bbmin[2], 
                                    bbmax[0], bbmax[1], bbmax[2])) {
      std::cerr << m_className << "::Plot2d: Cannot retrieve bounding box.\n";
      return;
    }
    PlotLimits(bbmin, bbmax, x0, y0, x1, y1);
  }

  auto canvas = GetCanvas();
  canvas->cd();
  canvas->SetTitle("Geometry");

  bool empty = !RangeSet(gPad);
  const double bm = canvas->GetBottomMargin();
  const double lm = canvas->GetLeftMargin();
  const double rm = canvas->GetRightMargin();
  const double tm = canvas->GetTopMargin();
  if (!empty) {
    TPad* pad = new TPad("geo", "", 0, 0, 1, 1);
    pad->SetFillStyle(0);
    pad->SetFrameFillStyle(0);
    pad->Draw();
    pad->cd();
  }
  gPad->Range(x0 - (x1 - x0) * (lm / (1. - rm - lm)),
              y0 - (y1 - y0) * (bm / (1. - tm - lm)),
              x1 + (x1 - x0) * (rm / (1. - rm - lm)),
              y1 + (y1 - y0) * (tm / (1. - tm - lm)));
  TPolyLine pl;
  pl.SetLineWidth(2);
  for (size_t i = 0; i < nSolids; ++i) {
    auto solid = m_geometry->GetSolid(i);
    if (!solid) continue;
    std::vector<Panel> panels;
    solid->Cut(m_proj[2][0], m_proj[2][1], m_proj[2][2],
               m_plane[0], m_plane[1], m_plane[2], panels);
    for (const auto& panel : panels) {
      const auto nv = panel.xv.size();
      if (nv < 3) continue;
      std::vector<double> xpl;
      std::vector<double> ypl;
      for (size_t j = 0; j < nv; ++j) {
        double u, v;
        ToPlane(panel.xv[j], panel.yv[j], panel.zv[j], u, v);
        xpl.push_back(u);
        ypl.push_back(v);
      }
      xpl.push_back(xpl[0]);
      ypl.push_back(ypl[0]);
      if (panel.colour < 0) {
        pl.SetLineColor(kBlack);
      } else {
        pl.SetLineColor(panel.colour);
      }
      pl.DrawPolyLine(nv + 1, xpl.data(), ypl.data(), "same");
    }
  }
  gPad->Update();
}

void ViewGeometry::Reset() {
  for (auto it = m_volumes.begin(), end = m_volumes.end(); it != end; ++it) {
    if (*it) {
      TGeoShape* shape = (*it)->GetShape();
      if (shape) delete shape;
      delete *it;
    }
  }
  m_volumes.clear();
  for (auto it = m_media.begin(), end = m_media.end(); it != end; ++it) {
    if (*it) {
      TGeoMaterial* material = (*it)->GetMaterial();
      if (material) delete material;
      delete *it;
    }
  }
  m_media.clear();

  m_geoManager.reset(nullptr);
}
 
}
