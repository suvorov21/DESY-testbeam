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
  virtual ~CrossingReconstruction() {;}

  virtual bool Initialize();
  virtual bool SelectEvent(const Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples], Event &event);

 private:

};

#endif  // SRC_RECONSTRUCTION_CROSSINGRECONSTRUCTION_HXX_
