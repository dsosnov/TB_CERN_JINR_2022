#!/bin/env sh
source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc11-opt/setup.sh
fullpath=$2
folder=$1
pushd "$fullpath"
hadd -j -f "${folder}-merged.root" "${folder}/out_tiger_"*".root"
popd