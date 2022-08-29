#include "Garfield/AvalancheGrid.hh"

#include <TF1.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>

#include "Garfield/Medium.hh"
#include "Garfield/Random.hh"

namespace Garfield {
void AvalancheGrid::SetZGrid(Grid &av, const double ztop, const double zbottom,
                             const int zsteps) {
  // Creating the z-coordinate grid.
  av.zsteps = zsteps;
  av.zStepSize = (ztop - zbottom) / zsteps;
  // Loop bassed grid creation.
  for (int i = 0; i < zsteps; i++) {
    av.zgrid.push_back(zbottom + i * av.zStepSize);
  }
}

void AvalancheGrid::SetYGrid(Grid &av, const double ytop, const double ybottom,
                             const int ysteps) {
  // Idem to SetZGrid for the x-coordinate grid.
  av.ysteps = ysteps;
  av.yStepSize = (ytop - ybottom) / ysteps;

  if (av.yStepSize == 0)
    av.yStepSize = 1;

  for (int i = 0; i < ysteps; i++) {
    av.ygrid.push_back(ybottom + i * av.yStepSize);
  }
}

void AvalancheGrid::SetXGrid(Grid &av, const double xtop, const double xbottom,
                             const int xsteps) {
  // Idem to SetZGrid for the x-coordinate grid.
  av.xsteps = xsteps;
  av.xStepSize = (xtop - xbottom) / xsteps;

  if (av.xStepSize == 0)
    av.xStepSize = 1;

  for (int i = 0; i < xsteps; i++) {
    av.xgrid.push_back(xbottom + i * av.xStepSize);
  }
}

void AvalancheGrid::SetGrid(const double xmin, const double xmax,
                            const int xsteps, const double ymin,
                            const double ymax, const int ysteps,
                            const double zmin, const double zmax,
                            const int zsteps) {
  m_avgrid.gridset = true;

  if (zmin >= zmax || zsteps <= 0 || xmin > xmax || xsteps <= 0 ||
      ymin > ymax || ysteps <= 0) {
    std::cerr << m_className
              << "::SetGrid:Error: Grid is not properly defined.\n";
    return;
  }

  // Setting grid

  SetZGrid(m_avgrid, zmax, zmin, zsteps);
  SetYGrid(m_avgrid, ymax, ymin, ysteps);
  SetXGrid(m_avgrid, xmax, xmin, xsteps);

  if (m_debug) {
    std::cerr << m_className << "::SetGrid: Grid created:\n";
    std::cerr << "       x range = (" << xmin << "," << xmax << ").\n";
    std::cerr << "       y range = (" << ymin << "," << ymax << ").\n";
    std::cerr << "       z range = (" << zmin << "," << zmax << ").\n";
  }
}

int AvalancheGrid::GetAvalancheSize(double dx, const int nsize,
                                    const double alpha, const double eta) {
  // Algorithm to get the size of the avalanche after it has propagated over a
  // distance dx.

  int newnsize = 0; // Holder for final size.

  const double k = eta / alpha;
  const double ndx =
      exp((alpha - eta) * dx); // Scaling Townsend and Attachment coef. to 1/mm.
  // If the size is higher than 1e3 the central limit theorem will be used to
  // describe the growth of the Townsend avalanche.
  if (nsize < 1e3) {
    // Running over all electrons in the avalanche.
    for (int i = 0; i < nsize; i++) {
      // Draw a random number from the uniform distribution (0,1).
      double s = RndmUniformPos();
      // Condition to which the random number will be compared. If the number is
      // smaller than the condition, nothing happens. Otherwise, the single
      // electron will be attached or retrieve additional electrons from the
      // gas.
      double condition = k * (ndx - 1) / (ndx - k);

      if (s >= condition)
        newnsize += (int)(1 + log((ndx - k) * (1 - s) / (ndx * (1 - k))) /
                                  log(1 - (1 - k) / (ndx - k)));
    }

  } else {
    // Central limit theorem.
    const double sigma = sqrt((1 + k) * ndx * (ndx - 1) / (1 - k));
    newnsize = RndmGaussian(nsize * ndx, sqrt(nsize) * sigma);
  }

  return newnsize;
}

bool AvalancheGrid::SnapToGrid(Grid &av, const double x, const double y,
                               const double z, const double /*v*/,
                               const int n) {
  // Snap electron from AvalancheMicroscopic to the predefined grid.
  if (!av.gridset) {
    std::cerr << m_className << "::SnapToGrid:Error: grid is not defined.\n";
    return false;
  }
  // Finding the z position on the grid.

  int indexX, indexY, indexZ = 0;

  // TODO: Snap must be dependent on the direction of drift.
  indexX = round((x - av.xgrid.front()) / av.xStepSize);
  indexY = floor((y - av.ygrid.front()) / av.yStepSize);
  indexZ = round((z - av.zgrid.front()) / av.zStepSize);

  if (m_debug)
    std::cerr << m_className << "::SnapToGrid: ix = " << indexX
              << ", iy = " << indexY << ", iz = " << indexZ << ".\n\n";
  if (indexX < 0 || indexX >= av.xsteps || indexY < 0 || indexY >= av.ysteps ||
      indexZ < 0 || indexZ >= av.zsteps) {
    if (m_debug)
      std::cerr << m_className << "::SnapToGrid: Point is outside the grid.\n";
    return false;
  }

  AvalancheNode newNode;
  newNode.ix = indexX;
  newNode.iy = indexY;
  newNode.iz = indexZ;
  if (!GetParameters(newNode)) {
    // TODO:Memory leak?
    if (m_debug)
      std::cerr << m_className
                << "::SnapToGrid: Could not retrieve parameters from sensor.\n";
    return false;
  }

  // When snapping the electron to the grid the distance traveled can yield
  // additional electrons or attachment.

  double step = z - av.zgrid[indexZ];

  if (newNode.velNormal[0] != 0) {
    step = x - av.xgrid[indexX];
  } else if (newNode.velNormal[1] != 0) {
    step = y - av.ygrid[indexY];
  }

  const int nholder =
      GetAvalancheSize(step, n, newNode.townsend, newNode.attachment);
  if (nholder == 0) {
    if (m_debug)
      std::cerr << m_className << "::SnapToGrid: n from 1 to 0 -> cancel.\n";
    return false;
  }

  newNode.n = nholder < m_MaxSize ? nholder : m_MaxSize;
  av.N += newNode.n;

  bool alreadyExists = false;

  for (AvalancheNode &existingNode : m_activeNodes) {
    if (existingNode.ix == newNode.ix && existingNode.iy == newNode.iy &&
        existingNode.iz == newNode.iz) {
      alreadyExists = true;
      existingNode.n += newNode.n;
    }
  }

  // TODO: What if time changes as you are importing avalanches?
  newNode.time = av.time;
  if (!alreadyExists)
    m_activeNodes.push_back(newNode);

  if (m_debug)
    std::cerr << m_className << "::SnapToGrid: n from 1 to " << nholder
              << ".\n";

  if (m_debug)
    std::cerr << m_className << "::SnapToGrid: Snapped to (x,y,z) = (" << x
              << " -> " << av.xgrid[indexX] << ", " << y << " -> "
              << av.ygrid[indexY] << ", " << z << " -> " << av.zgrid[indexZ]
              << ").\n";
  return true;
}

void AvalancheGrid::NextAvalancheGridPoint(Grid &av) {
  // This main function propagates the electrons and applies the avalanche
  // statistics.
  int Nholder = 0; // Holds the avalanche size before propagating it to the
  // next point in the grid.
  av.run = false;
  for (AvalancheNode &node : m_activeNodes) { // For every avalanche node

    if (!node.active) {
      continue;
    } else {
      av.run = true;
    }

    if (m_debug)
      std::cerr << m_className << "::NextAvalancheGridPoint:(ix,iy,iz) = ("
                << node.ix << "," << node.iy << "," << node.iz << ").\n";

    // Get avalanche size.
    Nholder = node.n;

    if (Nholder == 0)
      continue; // If empty go to next point.
    // If the total avalanche size is smaller than the set saturation
    // limit the GetAvalancheSize function is utilized to obtain the size
    // after its propagation to the next z-coordinate grid point. Else,
    // the size will be kept constant under the propagation.

    if (!m_layerIndix && av.N < m_MaxSize) {
      int holdnsize = GetAvalancheSize(node.stepSize, node.n, node.townsend,
                                       node.attachment);

      if (m_MaxSize - av.N < holdnsize - node.n)
        holdnsize = m_MaxSize - av.N + node.n;

      node.n = holdnsize;

    } else if (m_layerIndix && m_NLayer[node.layer - 1] < m_MaxSize) {
      int holdnsize = GetAvalancheSize(node.stepSize, node.n, node.townsend,
                                       node.attachment);

      if (m_MaxSize - m_NLayer[node.layer - 1] < holdnsize - node.n)
        holdnsize = m_MaxSize - m_NLayer[node.layer - 1] + node.n;

      node.n = holdnsize;

    } else {
      m_Saturated = true;

      if (m_SaturationTime == -1)
        m_SaturationTime = node.time + node.dt;
    }
    // Produce induced signal on readout electrodes.

    m_sensor->AddSignal(-(Nholder + node.n) / 2, node.time, node.time + node.dt,
                        av.xgrid[node.ix], av.ygrid[node.iy], av.zgrid[node.iz],
                        av.xgrid[node.ix + node.velNormal[0]],
                        av.ygrid[node.iy + node.velNormal[1]],
                        av.zgrid[node.iz + node.velNormal[2]], false, true);

    // Update total number of electrons.

    if (m_layerIndix)
      m_NLayer[node.layer - 1] += node.n - Nholder;

    av.N += node.n - Nholder;

    if (m_diffusion) {
      // TODO: to impliment
    }

    if (m_debug)
      std::cerr << "n = " << Nholder << " -> " << node.n << ".\n";

    // Update position index.

    node.ix += node.velNormal[0];
    node.iy += node.velNormal[1];
    node.iz += node.velNormal[2];

    // After all active grid points have propagated, update the time.
    if (m_debug)
      std::cerr << "t = " << node.time << " -> ";
    node.time += node.dt;
    if (m_debug)
      std::cerr << node.time << ".\n";

    DeactivateNode(node);
  }
  if (m_debug)
    std::cerr << "N = " << av.N << ".\n\n";
}

void AvalancheGrid::DeactivateNode(AvalancheNode &node) {

  if (node.n == 0)
    node.active = false;

  if (node.velNormal[2] != 0) {
    if ((node.velNormal[2] < 0 && node.iz == 0) ||
        (node.velNormal[2] > 0 && node.iz == m_avgrid.zsteps - 1))
      node.active = false;
  } else if (node.velNormal[1] != 0) {
    if ((node.velNormal[1] < 0 && node.iy == 0) ||
        (node.velNormal[1] > 0 && node.iy == m_avgrid.ysteps - 1))
      node.active = false;
  } else {
    if ((node.velNormal[0] < 0 && node.ix == 0) ||
        (node.velNormal[0] > 0 && node.ix == m_avgrid.xsteps - 1))
      node.active = false;
  }

  double e[3], v;
  int status;
  Medium *m = nullptr;
  m_sensor->ElectricField(m_avgrid.xgrid[node.ix], m_avgrid.ygrid[node.iy],
                          m_avgrid.zgrid[node.iz], e[0], e[1], e[2], v, m,
                          status);

  if (status == -5 || status == -6) {
    node.active = false; // If not inside a gas gap return false to terminate
  }

  if (m_debug && !node.active)
    std::cerr << m_className << "::DeactivateNode: Node deactivated.\n";
}

void AvalancheGrid::StartGridAvalanche() {
  // Start the AvalancheGrid algorithm.
  if ((!m_avmc && !m_driftAvalanche) || !m_sensor)
    return;

  if (!m_importAvalanche && m_avmc)
    GetElectronsFromAvalancheMicroscopic();

  std::cerr << m_className
            << "::StartGridAvalanche::Starting grid based simulation with "
            << m_avgrid.N << " initial electrons.\n";
  if (m_avgrid.N <= 0) {
    std::cerr << m_className << "::StartGridAvalanche::Cancelled.\n";
    return;
  }

  m_nestart = m_avgrid.N;

  // Main loop.
  while (m_avgrid.run == true) {
    if (m_debug)
      std::cerr
          << "============ \n"
          << m_className
          << "::StartGridAvalanche: Looping over nodes.\n ============ \n";
    NextAvalancheGridPoint(m_avgrid);
  }

  std::vector<double> tlist = {};
  for (AvalancheNode &node : m_activeNodes) {
    tlist.push_back(node.time);
  }
  double maxTime = *max_element(std::begin(tlist), std::end(tlist));

  if (m_Saturated)
    std::cerr << m_className
              << "::StartGridAvalanche::Avalanche maximum size of " << m_MaxSize
              << " electrons reached at " << m_SaturationTime << " ns.\n";

  std::cerr << m_className
            << "::StartGridAvalanche::Final avalanche size = " << m_avgrid.N
            << " ended at t = " << maxTime << " ns.\n";

  return;
}

void AvalancheGrid::CreateAvalanche(const double x, const double y,
                                    const double z, const double t,
                                    const int n) {
  m_driftAvalanche = true;

  if (m_avgrid.time == 0 && m_avgrid.time != t && m_debug)
    std::cerr << m_className
              << "::CreateAvalanche::Overwriting start time of avalanche for t "
                 "= 0 to "
              << t << ".\n";
  m_avgrid.time = t;

  if (SnapToGrid(m_avgrid, x, y, z, 0, n) && m_debug)
    std::cerr << m_className
              << "::CreateAvalanche::Electron added at (t,x,y,z) =  (" << t
              << "," << x << "," << y << "," << z << ").\n";
}

void AvalancheGrid::GetElectronsFromAvalancheMicroscopic() {
  // Get the information of the electrons from the AvalancheMicroscopic class.
  if (!m_avmc)
    return;

  if (!m_importAvalanche)
    m_importAvalanche = true;

  int np = m_avmc->GetNumberOfElectronEndpoints();

  if (np == 0)
    return;

  // Get initial positions of electrons
  double x1, y1, z1, t1;
  double x2, y2, z2, t2;
  double e1, e2;
  int status;

  double vel = 0.;

  for (int i = 0; i < np; ++i) {
    m_avmc->GetElectronEndpoint(i, x1, y1, z1, t1, e1, x2, y2, z2, t2, e2,
                                status);

    vel = (z2 - z1) / (t2 - t1);

    m_avgrid.time = t2;

    if (SnapToGrid(m_avgrid, x2, y2, z2, vel) && m_debug)
      std::cerr << m_className
                << "::GetElectronsFromAvalancheMicroscopic::Electron added at "
                   "(x,y,z) =  ("
                << x2 << "," << y2 << "," << z2 << ").\n";
  }
}

bool AvalancheGrid::GetParameters(AvalancheNode &node) {
  if (!m_sensor)
    return false;

  double x = m_avgrid.xgrid[node.ix];
  double y = m_avgrid.ygrid[node.iy];
  double z = m_avgrid.zgrid[node.iz];

  if (m_debug)
    std::cerr << m_className
              << "::GetParametersFromSensor::Getting parameters from "
                 "(x,y,z) =  ("
              << x << "," << y << "," << z << ").\n";

  double e[3], v;
  int status;
  Medium *m = nullptr;
  m_sensor->ElectricField(x, y, z, e[0], e[1], e[2], v, m, status);

  if (m_debug)
    std::cerr << m_className << "::GetParametersFromSensor::status = " << status
              << ".\n";

  if (status == -5 || status == -6)
    return false; // If not inside a gas gap return false to terminate

  if (m_Townsend >=
      0) { // If Townsend coef. is not set by user, take it from the sensor.
    node.townsend = m_Townsend;
  } else {
    m->ElectronTownsend(e[0], e[1], e[2], 0., 0., 0., node.townsend);
  }

  if (m_Attachment >=
      0) { // If attachment coef. is not set by user, take it from the sensor.
    node.attachment = m_Attachment;
  } else {
    m->ElectronAttachment(e[0], e[1], e[2], 0., 0., 0., node.attachment);
  }

  if (m_Velocity >
      0) { // If velocity is not set by user, take it from the sensor.
    node.velocity = m_Velocity;
    node.velNormal = m_velNormal;
  } else {
    double vx, vy, vz;
    m->ElectronVelocity(e[0], e[1], e[2], 0., 0., 0., vx, vy, vz);

    double vel = sqrt(vx * vx + vy * vy + vz * vz);
    if (vel == 0.)
      return false;
    if (vel != std::abs(vx) && vel != std::abs(vy) && vel != std::abs(vz))
      return false;
    int nx = (int)round(vx / vel);
    int ny = (int)round(vy / vel);
    int nz = (int)round(vz / vel);

    node.velNormal = {nx, ny, nz};
    node.velocity = -std::abs(vel);
  }

  if (node.velNormal[0] != 0) {
    node.stepSize = m_avgrid.xStepSize;
  } else if (node.velNormal[1] != 0) {
    node.stepSize = m_avgrid.yStepSize;
  } else {
    node.stepSize = m_avgrid.zStepSize;
  }

  if (m_debug)
    std::cerr << m_className
              << "::GetParametersFromSensor::stepSize = " << node.stepSize
              << "[cm].\n";

  if (m_debug)
    std::cerr << m_className << "::GetParametersFromSensor::velNormal = ("
              << node.velNormal[0] << ", " << node.velNormal[1] << ", "
              << node.velNormal[2] << ") [1].\n";

  node.dt = std::abs(node.stepSize / node.velocity);

  // print
  if (m_debug || !m_printPar)
    std::cerr << m_className << "::GetParametersFromSensor::Electric field = ("
              << e[0] / 1000 << ", " << e[1] / 1000 << ", " << e[2] / 1000
              << ") [kV/cm].\n";

  if (m_debug || !m_printPar)
    std::cerr << m_className
              << "::GetParametersFromSensor::Townsend = " << node.townsend
              << " [1/cm], Attachment = " << node.attachment
              << " [1/cm], Velocity = " << node.velocity << " [cm/ns].\n";

  if (m_debug)
    std::cerr << m_className << "::StartGridAvalanche::Time steps per loop "
              << node.dt << " ns.\n";
  m_printPar = true;
  return true;
}

void AvalancheGrid::Reset() {

  std::cerr << m_className << "::Reset::Resetting AvalancheGrid.\n";
  m_avgrid.time = 0;
  m_avgrid.N = 0;
  m_avgrid.run = true;

  m_Saturated = false;
  m_SaturationTime = -1;

  m_driftAvalanche = false;

  m_activeNodes.clear();
  m_layerIndix = false;
  m_NLayer.clear();
}

void AvalancheGrid::AsignLayerIndex(ComponentParallelPlate *RPC) {
  int im = 0;
  double epsM = 0;
  std::vector<double> nLayer(RPC->NumberOfLayers(), 0);
  for (AvalancheNode &node : m_activeNodes) {
    double y = m_avgrid.ygrid[node.iy];
    RPC->getLayer(y, im, epsM);
    node.layer = im;
    nLayer[im - 1] += node.n;
    // std::cerr << m_className << "::AsignLayerIndex::im = "<<im<<".\n";
  }

  m_NLayer = nLayer;
  m_layerIndix = true;
}

} // namespace Garfield
