#!/bin/env sh
source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc11-opt/setup.sh
fullpath=$3
filename=$1
folder=$2
pushd "$fullpath"
root -b -q tiger_tree_converter_bin.cxx'("'${filename}'", "'$folder'")'
popd