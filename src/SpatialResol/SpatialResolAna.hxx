#ifndef SRC_SPATIALRESOL_SPATIALRESOLANA_HXX_
#define SRC_SPATIALRESOL_SPATIALRESOLANA_HXX_

#include "AnalysisBase.hxx"
#include "SelectionBase.hxx"

class SpatialResolAna: public AnalysisBase {
 public:
  SpatialResolAna(int argc, char** argv);
  virtual ~SpatialResolAna() {;}

  // Initialise histoes, input files, selections
  bool Initialize();
  // Process the selection output called Event
  bool ProcessEvent(const Event event);
  // write output files (histos, trees)
  bool WriteOutput();
};

#endif  // SRC_SPATIALRESOL_SPATIALRESOLANA_HXX_
