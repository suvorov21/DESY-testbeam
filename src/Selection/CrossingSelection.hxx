#ifndef SRC_SELECTION_CROSSINGSELECTION_HXX_
#define SRC_SELECTION_CROSSINGSELECTION_HXX_

#include "SelectionBase.hxx"

//! Selection for passing through tracks

//! Was developed for CERN beam test. Was known as 3D selection
//! This selection is optimised for going through tracks.
//! The clusters at the beginning and at the end MicroMegas are selected first.
//! Then all possible cluster matchs in 3D space are studied.
//! If clusters are connected with hits this supposed to be a found track.
class CrossingSelection: public SelectionBase {
 public:
  CrossingSelection();
  virtual ~CrossingSelection() {;}

  virtual bool Initialize();
  virtual bool SelectEvent(const Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples], Event &event);

 private:

};

#endif  // SRC_SELECTION_CROSSINGSELECTION_HXX_
