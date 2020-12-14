#!bin/sh

source /cvmfs/sft.cern.ch/lcg/contrib/gcc/7.3.0binutils/x86_64-centos7-gcc7-opt/setup.sh
source /afs/cern.ch/work/s/ssuvorov/public/ROOT/root-6.18.00-build_gcc73/bin/thisroot.sh

export SOFTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )/"
echo "Project path:  $SOFTDIR"
