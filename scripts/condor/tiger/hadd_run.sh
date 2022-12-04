#!/bin/env sh
source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc11-opt/setup.sh
fullpath=$2
folder=$1
selection=$3
pushd "$fullpath"
hadd -f "${folder}-merged${selection}.root" "${folder}/out_tiger_RUN_"*"-SubRUN_${selection}"*"_GEMROC_"*".root"
popd