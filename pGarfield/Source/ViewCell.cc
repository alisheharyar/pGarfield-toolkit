#include <cmath>
#include <iostream>

#include <TEllipse.h>
#include <TGeoBBox.h>
#include <TGeoPgon.h>
#include <TLine.h>
#include <TMarker.h>
#include <TPolyLine.h>

#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/ComponentNeBem2d.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/ViewCell.hh"

namespace Garfield {

ViewCell::ViewCell() : ViewBase("ViewCell") {}

void ViewCell::SetComponent(ComponentAnalyticField* cmp) {
  if (!cmp) {
    std::cerr << m_className << "::SetComponent: Null pointer.\n";
    return;
  }
  m_component = cmp;
}

void ViewCell::SetComponent(ComponentNeBem2d* cmp) {
  if (!cmp) {
    std::cerr << m_className << "::SetComponent: Null pointer.\n";
    return;
  }
  m_nebem = cmp;
}

void ViewCell::Plot2d() {
  if (!Plot(true)) {
    std::cerr << m_className << "::Plot2d: Error creating plot.\n";
  }
}

void ViewCell::Plot3d() {
  if (!Plot(false)) {
    std::cerr << m_className << "::Plot3d: Error creating plot.\n";
  }
}

bool ViewCell::Plot(const bool twod) {
  if (!m_component && !m_nebem) {
    std::cerr << m_className << "::Plot: Component is not defined.\n";
    return false;
  }

  double pmin = 0., pmax = 0.;
  if (m_component) {
    if (!m_component->GetVoltageRange(pmin, pmax)) {
      std::cerr << m_className << "::Plot: Component is not ready.\n";
      return false;
    }
  } else {
    if (!m_nebem->GetVoltageRange(pmin, pmax)) {
      std::cerr << m_className << "::Plot: Component is not ready.\n";
      return false;
    }
  }

  // Get the bounding box.
  double x0 = m_xMinBox, y0 = m_yMinBox, z0 = m_zMinBox;
  double x1 = m_xMaxBox, y1 = m_yMaxBox, z1 = m_zMaxBox;
  if (twod && m_userPlotLimits) {
    x0 = m_xMinPlot;
    y0 = m_yMinPlot;
    x1 = m_xMaxPlot;
    y1 = m_yMaxPlot;
  } else if (!m_userBox) {
    if (m_component) {
      if (!m_component->GetBoundingBox(x0, y0, z0, x1, y1, z1)) {
        std::cerr << m_className << "::Plot:\n"
                  << "    Bounding box cannot be determined.\n"
                  << "    Call SetArea first.\n";
        return false;
      }
    } else {
      if (!m_nebem->GetBoundingBox(x0, y0, z0, x1, y1, z1)) {
        std::cerr << m_className << "::Plot:\n"
                  << "    Bounding box cannot be determined.\n"
                  << "    Call SetArea first.\n";
        return false;
      }
    }
  }
  const double dx = std::max(fabs(x0), fabs(x1));
  const double dy = std::max(fabs(y0), fabs(y1));
  const double dz = std::max(fabs(z0), fabs(z1));

  auto canvas = GetCanvas();
  canvas->cd();
  canvas->SetTitle("Cell layout");

  if (twod) {
    const bool empty = !RangeSet(gPad);
    const double bm = canvas->GetBottomMargin();
    const double lm = canvas->GetLeftMargin();
    const double rm = canvas->GetRightMargin();
    const double tm = canvas->GetTopMargin();
    if (!empty) {
      TPad* pad = new TPad("cell", "", 0, 0, 1, 1);
      pad->SetFillStyle(0);
      pad->SetFrameFillStyle(0);
      pad->Draw();
      pad->cd();
    }
    gPad->Range(x0 - (x1 - x0) * (lm / (1. - rm - lm)),
                y0 - (y1 - y0) * (bm / (1. - tm - lm)),
                x1 + (x1 - x0) * (rm / (1. - rm - lm)),
                y1 + (y1 - y0) * (tm / (1. - tm - lm)));
  }

  if (m_nebem) return PlotNeBem(twod);

  // Get the cell type.
  const std::string cellType = m_component->GetCellType();

  // Get the x/y periodicities.
  double sx = 0., sy = 0.;
  const bool perX = m_component->GetPeriodicityX(sx);
  const bool perY = m_component->GetPeriodicityY(sy);
  // Determine the number of periods present in the cell.
  const int nMinX = perX ? int(x0 / sx) - 1 : 0;
  const int nMaxX = perX ? int(x1 / sx) + 1 : 0;
  const int nMinY = perY ? int(y0 / sy) - 1 : 0;
  const int nMaxY = perY ? int(y1 / sy) + 1 : 0;
  // Get the phi periodicity.
  double sphi = 0.;
  const bool perPhi = m_component->GetPeriodicityPhi(sphi);
  const int nPhi = perPhi ? int(360. / sphi) : 0;
  sphi *= DegreeToRad; 
  if (!twod) SetupGeo(dx, dy, dz);
  const bool polar = m_component->IsPolar();

  // Get the number of wires.
  const unsigned int nWires = m_component->GetNumberOfWires();
  std::vector<std::string> wireTypes;
  // Loop over the wires.
  for (unsigned int i = 0; i < nWires; ++i) {
    double xw = 0., yw = 0., dw = 0., vw = 0., lw = 0., qw = 0.;
    std::string lbl;
    int nTrap;
    m_component->GetWire(i, xw, yw, dw, vw, lbl, lw, qw, nTrap);
    auto it = std::find(wireTypes.begin(), wireTypes.end(), lbl);
    const int type = std::distance(wireTypes.begin(), it);
    if (it == wireTypes.end()) wireTypes.push_back(lbl); 
    if (polar) {
      const double r = xw;
      for (int j = 0; j <= nPhi; ++j) {
        const double phi = yw * Garfield::DegreeToRad + j * sphi;
        const double x = r * cos(phi);
        const double y = r * sin(phi);
        if (twod) {
          PlotWire(x, y, dw, type);
        } else {
          PlotWire(x, y, dw, type, std::min(0.5 * lw, dz));
        }
      }
      continue;
    }
    for (int nx = nMinX; nx <= nMaxX; ++nx) {
      const double x = xw + nx * sx;
      if (x + 0.5 * dw <= x0 || x - 0.5 * dw >= x1) continue;
      for (int ny = nMinY; ny <= nMaxY; ++ny) {
        const double y = yw + ny * sy;
        if (y + 0.5 * dw <= y0 || y - 0.5 * dw >= y1) continue;
        if (twod) {
          PlotWire(x, y, dw, type);
        } else {
          PlotWire(x, y, dw, type, std::min(0.5 * lw, dz));
        }
      }
    }
  }

  // Draw the x planes.
  const unsigned int nPlanesX = m_component->GetNumberOfPlanesX();
  for (unsigned int i = 0; i < nPlanesX; ++i) {
    double xp = 0., vp = 0.;
    std::string lbl;
    m_component->GetPlaneX(i, xp, vp, lbl);
    for (int nx = nMinX; nx <= nMaxX; ++nx) {
      const double x = xp + nx * sx;
      if (x < x0 || x > x1) continue;
      if (twod) {
        PlotPlane(x, y0, x, y1);
      } else {
        const double width = std::min(0.01 * dx, 0.01 * dy);
        PlotPlane(width, dy, dz, x, 0);
      }
    }
  }

  // Draw the y planes.
  const unsigned int nPlanesY = m_component->GetNumberOfPlanesY();
  for (unsigned int i = 0; i < nPlanesY; ++i) {
    double yp = 0., vp = 0.;
    std::string lbl;
    m_component->GetPlaneY(i, yp, vp, lbl);
    for (int ny = nMinY; ny <= nMaxY; ++ny) {
      const double y = yp + ny * sy;
      if (y < y0 || y > y1) continue;
      if (twod) {
        PlotPlane(x0, y, x1, y);
      } else {
        const double width = std::min(0.01 * dx, 0.01 * dy);
        PlotPlane(dx, width, dz, 0, y);
      }
    }
  }

  // Draw the r and phi planes.
  const unsigned int nPlanesR = m_component->GetNumberOfPlanesR();
  const unsigned int nPlanesPhi = m_component->GetNumberOfPlanesPhi();
  if (nPlanesR > 0 || nPlanesPhi > 0) {
    double vp = 0.;
    std::string lbl;
    std::array<double, 2> phi = {0., 360.};
    double dphi = 360.;
    if (nPlanesPhi == 2 && !m_component->GetPeriodicityPhi(dphi)) {
      m_component->GetPlanePhi(0, phi[0], vp, lbl);
      m_component->GetPlanePhi(1, phi[1], vp, lbl);
    }
    double w = 0.01 * std::min(dx, dy);
    std::array<double, 2> r = {0., 0.};
    if (nPlanesR > 0) {
      m_component->GetPlaneR(0, r[0], vp, lbl);
      w = 0.02 * r[0];
    }
    if (nPlanesR == 2) {
      m_component->GetPlaneR(1, r[1], vp, lbl);
      w = 0.01 * (r[1] - r[0]);
    }
    for (unsigned int i = 0; i < nPlanesR; ++i) {
      if (twod) {
        TEllipse circle;
        circle.SetDrawOption("same");
        circle.SetFillStyle(0);
        circle.SetNoEdges();
        circle.DrawEllipse(0, 0, r[i], r[i], phi[0], phi[1], 0);
        continue;
      }
      const auto med = m_geo->GetMedium("Metal");
      const auto tubs = m_geo->MakeTubs("Plane", med, r[i] - w, r[i] + w,
                                        dz, phi[0], phi[1]);
      tubs->SetLineColor(kGreen - 5);
      tubs->SetTransparency(75);
      m_geo->GetTopVolume()->AddNode(tubs, 1);
    }
    if (nPlanesR == 1) {
      std::swap(r[0], r[1]);
    } else if (nPlanesR == 0) {
      r[1] = std::max(dx, dy); 
    }
    for (unsigned int i = 0; i < nPlanesPhi; ++i) {
      const double cp = cos(phi[i] * DegreeToRad);
      const double sp = sin(phi[i] * DegreeToRad);
      if (twod) {
        PlotPlane(r[0] * cp, r[0] * sp, r[1] * cp, r[1] * sp);
        continue;
      } 
      const auto med = m_geo->GetMedium("Metal");
      const double dr = 0.5 * (r[1] - r[0]);
      TGeoVolume* plane = m_geo->MakeBox("Plane", med, dr, w, dz);
      plane->SetLineColor(kGreen - 5);
      plane->SetTransparency(75);
      TGeoRotation* rot = new TGeoRotation("", phi[i], 0, 0);
      const double rm = 0.5 * (r[0] + r[1]); 
      TGeoCombiTrans* tr = new TGeoCombiTrans(cp * rm, sp * rm, 0, rot);
      m_geo->GetTopVolume()->AddNode(plane, 1, tr);
    } 
  }
  double rt = 0., vt = 0.;
  int nt = 0;
  std::string lbl;
  if (m_component->GetTube(rt, vt, nt, lbl)) {
    if (twod) {
      PlotTube(0., 0., rt, nt);
    } else {
      PlotTube(0., 0., 0.98 * rt, 1.02 * rt, nt, dz);
    }
  }

  if (!twod) {
    m_geo->CloseGeometry();
    m_geo->GetTopNode()->Draw("ogl");
  } else {
    gPad->Update();
  }

  return true;
}

bool ViewCell::PlotNeBem(const bool twod) {

  if (!twod) {
    std::cerr << m_className << "::PlotNeBem: 3D plot not implemented yet.\n";
    return false;
  }

  // Draw the regions.
  const unsigned int nRegions = m_nebem->GetNumberOfRegions();
  for (unsigned int i = nRegions; i-- > 0;) {
    std::vector<double> xv;
    std::vector<double> yv;
    Medium* medium = nullptr;
    unsigned int bctype = 1;
    double v = 0.;
    if (!m_nebem->GetRegion(i, xv, yv, medium, bctype, v)) continue;
    const unsigned int n = xv.size();
    if (n < 3) continue;
    TLine line;
    line.SetDrawOption("same");
    if (bctype == 4) {
      line.SetLineStyle(2);
    } else {
      line.SetLineStyle(1);
    }
    for (unsigned int j = 0; j < n; ++j) {
      const unsigned int k = j < n - 1 ? j + 1 : 0;
      line.DrawLine(xv[j], yv[j], xv[k], yv[k]);
    }
  }

  // Draw the wires.
  const unsigned int nWires = m_nebem->GetNumberOfWires();
  for (unsigned int i = 0; i < nWires; ++i) {
    double x = 0., y = 0., d = 0., v = 0., q = 0.;
    if (!m_nebem->GetWire(i, x, y, d, v, q)) continue;
    PlotWire(x, y, d, 0);
  }

  // Draw the straight-line segments.
  const unsigned int nSegments = m_nebem->GetNumberOfSegments();
  for (unsigned int i = 0; i < nSegments; ++i) {
    double x0 = 0., y0 = 0., x1 = 0., y1 = 0., v = 0.;
    if (!m_nebem->GetSegment(i, x0, y0, x1, y1, v)) continue;
    PlotPlane(x0, y0, x1, y1);
  }

  gPad->Update();
  return true;
}

void ViewCell::SetupGeo(const double dx, const double dy, const double dz) {

  if (!m_geo) {
    gGeoManager = nullptr;
    m_geo.reset(new TGeoManager("ViewCellGeoManager", "Cell layout"));
    TGeoMaterial* matVacuum = new TGeoMaterial("Vacuum", 0., 0., 0.);
    TGeoMaterial* matMetal = new TGeoMaterial("Metal", 63.546, 29., 8.92);
    TGeoMedium* medVacuum = new TGeoMedium("Vacuum", 0, matVacuum);
    TGeoMedium* medMetal = new TGeoMedium("Metal", 1, matMetal);
    m_geo->AddMaterial(matVacuum);
    m_geo->AddMaterial(medMetal->GetMaterial());
    TGeoVolume* world =
        m_geo->MakeBox("World", medVacuum, 1.05 * dx, 1.05 * dy, 1.05 * dz);
    m_geo->SetTopVolume(world);
  } else {
    TGeoVolume* top = m_geo->GetTopVolume();
    TGeoBBox* box = dynamic_cast<TGeoBBox*>(top);
    double halfLenghts[3] = {1.05 * dx, 1.05 * dy, 1.05 * dz};
    if (box) box->SetDimensions(halfLenghts);
  }
}

void ViewCell::PlotWire(const double x, const double y, const double d,
                        const int type) {
  if (m_useWireMarker) {
    int markerStyle = 1;
    if (type == 0) {
      markerStyle = kFullCircle;
    } else if (type == 1) {
      markerStyle = kOpenCircle;
    } else if (type == 2) {
      markerStyle = kFullSquare;
    } else if (type == 3) {
      markerStyle = kOpenSquare;
    } else {
      markerStyle = 26 + type;
    }
    TMarker marker;
    marker.SetMarkerStyle(markerStyle);
    marker.SetDrawOption("Psame");
    marker.DrawMarker(x, y);
    return;
  }

  TEllipse circle;
  circle.SetDrawOption("same");
  circle.SetFillColor(kWhite);
  circle.DrawEllipse(x, y, 0.5 * d, 0.5 * d, 0, 360, 0);
}

void ViewCell::PlotWire(const double x, const double y, const double d,
                        const int type, const double lz) {

  const auto medium = m_geo->GetMedium("Metal");
  TGeoVolume* wire = m_geo->MakeTube("Wire", medium, 0., 0.5 * d, lz);
  switch (type) {
    case 0:
      wire->SetLineColor(kGray + 2);
      break;
    case 1:
      wire->SetLineColor(kRed + 2);
      break;
    case 2:
      wire->SetLineColor(kPink + 3);
      break;
    case 3:
      wire->SetLineColor(kCyan + 3);
      break;
    default:
      wire->SetLineColor(kBlue + type);
      break;
  }
  m_geo->GetTopVolume()->AddNode(wire, 1, new TGeoTranslation(x, y, 0.));
}

void ViewCell::PlotTube(const double x0, const double y0, const double r,
                        const int n) {
  if (n <= 0) {
    TEllipse circle;
    circle.SetDrawOption("same");
    circle.SetFillStyle(0);
    circle.DrawEllipse(x0, y0, r, r, 0, 360, 0);
    return;
  }

  std::vector<double> x(n + 1);
  std::vector<double> y(n + 1);
  for (int i = 0; i <= n; ++i) {
    const double phi = i * TwoPi / double(n);
    x[i] = x0 + r * cos(phi);
    y[i] = y0 + r * sin(phi);
  }
  TPolyLine pline;
  pline.SetDrawOption("same");
  pline.DrawPolyLine(n + 1, x.data(), y.data());
}

void ViewCell::PlotTube(const double x0, const double y0, 
                        const double r1, const double r2, const int n,
                        const double lz) {

  TGeoVolume* tube = nullptr;
  if (n <= 0) {
    // Round tube.
    tube = m_geo->MakeTube("Tube", m_geo->GetMedium("Metal"), r1, r2, lz);
  } else {
    tube = m_geo->MakePgon("Tube", m_geo->GetMedium("Metal"), 0., 360., n, 2);
    TGeoPgon* pgon = dynamic_cast<TGeoPgon*>(tube->GetShape());
    if (pgon) {
      pgon->DefineSection(0, -lz, r1, r2);
      pgon->DefineSection(1, +lz, r1, r2);
    }
  }
  tube->SetLineColor(kGreen + 2);
  tube->SetTransparency(75);
  m_geo->GetTopVolume()->AddNode(tube, 1, new TGeoTranslation(x0, y0, 0));
}

void ViewCell::PlotPlane(const double x0, const double y0, 
                         const double x1, const double y1) {

  TLine line;
  line.SetDrawOption("same");
  line.DrawLine(x0, y0, x1, y1);
}

void ViewCell::PlotPlane(const double dx, const double dy, const double dz,
                         const double x0, const double y0) {

  const auto medium = m_geo->GetMedium("Metal");
  TGeoVolume* plane = m_geo->MakeBox("Plane", medium, dx, dy, dz);
  plane->SetLineColor(kGreen - 5);
  plane->SetTransparency(75);
  m_geo->GetTopVolume()->AddNode(plane, 1, new TGeoTranslation(x0, y0, 0));
}

}
