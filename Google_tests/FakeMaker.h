//
// Created by SERGEY SUVOROV on 28/11/2021.
//

#ifndef DESY_TESTBEAM_FAKEMAKER_H
#define DESY_TESTBEAM_FAKEMAKER_H

#include "TEvent.hxx"
#include "TCluster.hxx"

namespace fake {
  std::shared_ptr<TEvent> GetEvent();
  THitPtrVec GetTrack();
  TClusterPtrVec GetClusterdTrack();
  THitPtrVec GetCenteredCluster();
  THitPtrVec GetSideCLuster();
}

#endif // DESY_TESTBEAM_FAKEMAKER_H
