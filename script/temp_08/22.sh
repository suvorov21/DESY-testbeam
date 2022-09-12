#!/bin/bash
cd /afs/cern.ch/work/u/uyevarou/private/Desy21/dev_v2/desy_testbeam/build
export PATH=/cvmfs/sft.cern.ch/lcg/contrib/CMake/3.18.3/Linux-x86_64/bin/:${PATH}
. /cvmfs/sft.cern.ch/lcg/views/LCG_100/x86_64-centos7-clang11-opt/setup.sh
/afs/cern.ch/work/u/uyevarou/private/Desy21/dev_v2/desy_testbeam/build/SpatialResol.exe -b --param ../params/horiz_CERN22.ini -t 0 -i /afs/cern.ch/work/u/uyevarou/private/Desy21/dev_v2/desy_testbeam/script//temp_08/52.list -o /eos/experiment/neutplatform/np07/HAT/DESY_2021/vlada_reco/new//CERN22_output_v1/All_ERAMS_350V_412ns_p_0p5GeV_1_iter0.root; /afs/cern.ch/work/u/uyevarou/private/Desy21/dev_v2/desy_testbeam/build/SpatialResol.exe -b --param ../params/horiz_CERN22.ini -t 1 -i /afs/cern.ch/work/u/uyevarou/private/Desy21/dev_v2/desy_testbeam/script//temp_08/52.list -o /eos/experiment/neutplatform/np07/HAT/DESY_2021/vlada_reco/new//CERN22_output_v1/All_ERAMS_350V_412ns_p_0p5GeV_1_iter1.root; /afs/cern.ch/work/u/uyevarou/private/Desy21/dev_v2/desy_testbeam/build/SpatialResol.exe -b --param ../params/horiz_CERN22.ini -t 2 -i /afs/cern.ch/work/u/uyevarou/private/Desy21/dev_v2/desy_testbeam/script//temp_08/52.list -o /eos/experiment/neutplatform/np07/HAT/DESY_2021/vlada_reco/new//CERN22_output_v1/All_ERAMS_350V_412ns_p_0p5GeV_1_iter2.root; 
