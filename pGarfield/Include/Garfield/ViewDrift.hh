#ifndef G_VIEW_DRIFT
#define G_VIEW_DRIFT

#include <string>
#include <vector>
#include <array>
#include <utility>
#include <mutex>

#include <Rtypes.h>

#include "GarfieldConstants.hh"
#include "ViewBase.hh"

namespace Garfield {

/// Visualize drift lines and tracks.

class ViewDrift : public ViewBase {
 public:
  /// Constructor.
  ViewDrift();
  /// Destructor.
  ~ViewDrift() = default;

  /// Delete existing drift lines, tracks and markers.
  void Clear();

  /// Draw the drift lines.
  void Plot(const bool twod = false, const bool axis = true,
            const bool snapshot = false);
  /// Make a 2D plot of the drift lines in the current viewing plane.
  void Plot2d(const bool axis = true, const bool snapshot = false);
  /// Make a 3D plot of the drift lines.
  void Plot3d(const bool axis = true, const bool ogl = true,
              const bool snapshot = false);

  /// Draw markers (or not) at every collision along a track.
  void EnableClusterMarkers(const bool on = true) { m_drawClusters = on; }
  /// Set the size of the cluster markers (see TAttMarker).
  void SetClusterMarkerSize(const double size);
  /// Set the size of the collision markers (see TAttMarker).
  void SetCollisionMarkerSize(const double size);

  /// Set the colour with which to draw electron drift lines.
  void SetColourElectrons(const short col) { m_colElectron = col; } 
  /// Set the colour with which to draw hole drift lines.
  void SetColourHoles(const short col) { m_colHole = col; } 
  /// Set the colour with which to draw ion drift lines.
  void SetColourIons(const short col) { m_colIon = col; } 
  /// Set the colour with which to draw charged particle tracks.
  void SetColourTracks(const short col) { m_colTrack = col; } 
  /// Set the colour with which to draw photons.
  void SetColourPhotons(const short col) { m_colPhoton = col; } 
  /// Set the colour with which to draw excitation markers.
  void SetColourExcitations(const short col) { m_colExcitation = col; } 
  /// Set the colour with which to draw ionisation markers.
  void SetColourIonisations(const short col) { m_colIonisation = col; } 
  /// Set the colour with which to draw attachment markers.
  void SetColourAttachments(const short col) { m_colAttachment = col; } 

  /// Get the number of drift lines stored. 
  size_t GetNumberOfDriftLines() const { return m_driftLines.size(); }
  /// Retrieve the coordinates of a given drift line.
  void GetDriftLine(const size_t i, 
                    std::vector<std::array<float, 3> >& driftLine, 
                    bool& electron) const;

  // Functions used by the transport classes.
  void NewDriftLine(const Particle particle, const size_t np, size_t& id, 
                    const float x0, const float y0, const float z0);
  void NewChargedParticleTrack(const size_t np, size_t& id, const float x0,
                               const float y0, const float z0);

  void SetDriftLinePoint(const size_t iL, const size_t iP,
                         const float x, const float y, const float z);
  void AddDriftLinePoint(const size_t iL, const float x, const float y,
                         const float z);
  void SetTrackPoint(const size_t iL, const size_t iP,
                     const float x, const float y, const float z);
  void AddTrackPoint(const size_t iL, const float x, const float y,
                     const float z);
  void AddExcitation(const float x, const float y, const float z);
  void AddIonisation(const float x, const float y, const float z);
  void AddAttachment(const float x, const float y, const float z);

  void AddPhoton(const float x0, const float y0, const float z0,
                 const float x1, const float y1, const float z1);

  friend class ViewFEMesh;

 private:
  std::mutex m_mutex;

  std::vector<std::pair<std::vector<std::array<float, 3> >,
                        Particle> > m_driftLines;

  std::vector<std::vector<std::array<float, 3> > > m_tracks;
  std::vector<std::array<std::array<float, 3>, 2> > m_photons;

  std::vector<std::array<float, 3> > m_exc;
  std::vector<std::array<float, 3> > m_ion;
  std::vector<std::array<float, 3> > m_att;

  double m_markerSizeCluster = 0.01;
  double m_markerSizeCollision = 0.5;
  
  short m_colTrack = kGreen + 3;
  short m_colPhoton = kBlue + 1;
  short m_colElectron = kOrange - 3;
  short m_colHole = kRed + 1;
  short m_colIon = kRed + 1;
  short m_colExcitation = kGreen + 3;
  short m_colIonisation = kOrange - 3;
  short m_colAttachment = kCyan + 3;

  bool m_drawClusters = false;

  bool SetPlotLimits2d();
  bool SetPlotLimits3d();
  void DrawMarkers2d(const std::vector<std::array<float, 3> >& points,
                     const short col, const double size);
  void DrawMarkers3d(const std::vector<std::array<float, 3> >& points,
                     const short col, const double size);
};
}
#endif
