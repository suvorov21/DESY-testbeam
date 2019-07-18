#ifndef SRC_SPATIALRESOL_SPATIALRESOLANA_HXX_
#define SRC_SPATIALRESOL_SPATIALRESOLANA_HXX_

#include "AnalysisBase.hxx"
#include "CrossingReconstruction.hxx"

/// Spatial resolution analysis
class SpatialResolAna: public AnalysisBase {
 public:
  SpatialResolAna(int argc, char** argv);
  virtual ~SpatialResolAna() {;}

  /// Initialise histoes, input files, selections
  bool Initialize();
  /// Process the selection output called Event
  bool ProcessEvent(const Event event);
  /// Write output files (histos, trees)
  /** Specify only for the values that are not included in the vector */
  bool WriteOutput();
};

#endif  // SRC_SPATIALRESOL_SPATIALRESOLANA_HXX_
