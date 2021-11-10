#ifndef SRC_BASE_RECONSTRUCTIONBASE_HXX_
#define SRC_BASE_RECONSTRUCTIONBASE_HXX_

/** @cond */
#include <vector>
#include <iostream>
#include <string>
#include <unistd.h>  // getopt on Mac
#include <getopt.h>
#include <fstream>   // read file lists
/** @endcond */

#include "TEvent.hxx"
#include "TCluster.hxx"
#include "Geom.hxx"

/// Template for the Reconstruction class
class ReconstructionBase {
 public:
  ReconstructionBase(): _verbose(0) {}
  virtual ~ReconstructionBase() = default;

  virtual bool Initialize(int verbose);
  virtual bool SelectEvent(const std::shared_ptr<TEvent>& event);

protected:
  int _verbose;
};

#endif  // SRC_BASE_RECONSTRUCTIONBASE_HXX_
