#ifndef SRC_RECONSTRUCTION_CROSSINGRECONSTRUCTION_HXX_
#define SRC_RECONSTRUCTION_CROSSINGRECONSTRUCTION_HXX_

#include "ReconstructionBase.hxx"

//! Reconstruction for passing through tracks

//! Was developed for CERN beam test. Was known as 3D Reconstruction
//! This Reconstruction is optimised for going through tracks.
//! The clusters at the beginning and at the end MicroMegas are selected first.
//! Then all possible cluster matchs in 3D space are studied.
//! If clusters are connected with hits this supposed to be a found track.
class CrossingReconstruction: public ReconstructionBase {
 public:
  CrossingReconstruction();

  virtual bool Initialize();
  virtual bool SelectEvent(const Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples],
                           TRawEvent* event
                           );

 private:
  const int cluster_offset = 3;

  const int cluster_threshold   = 50;
  const int cluster_time_gate   = 30;
  const int cluster_space_gate  = 3;

  const int time_gate = 30;

  const int box_width_t2    = time_gate;
  const int box_width_t1    = time_gate;
  const int scan_time_up    = time_gate;
  const int scan_time_down  = time_gate;
  const int box_width_j  = 3;

  const int track_collection_dist = 4;

  const int breaking_thr  = 0;  // max number of gaps
  const int column_thr    = 34; // min number of columns
  const int OOT_thr       = 20; // suppress cosmic pile up

  const uint cluster_lim  = 4;

  const bool pile_up_cut  = true;
};

#endif  // SRC_RECONSTRUCTION_CROSSINGRECONSTRUCTION_HXX_
