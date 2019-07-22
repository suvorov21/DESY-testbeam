#ifndef SRC_BASE_RECONSTRUCTIONBASE_HXX_
#define SRC_BASE_RECONSTRUCTIONBASE_HXX_

/** @cond */
#include <vector>
#include <iostream>
#include <string>
/** @endcond */

#include "DataStorage.hxx"

// It is important to include Selection.hxx after the stucts have been defined.
// #include "Selection.hxx"

/// Template for the Reconstruction class
class ReconstructionBase {
 public:
  ReconstructionBase() {;}
  virtual ~ReconstructionBase() {;}

  virtual bool Initialize();
  virtual bool SelectEvent(const Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples],
                  TEvent* event);

  virtual std::vector<std::vector<Int_t>> GetEmptyEvent();
};

#endif  // SRC_BASE_RECONSTRUCTIONBASE_HXX_
