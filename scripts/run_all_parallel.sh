#!/usr/bin/env bash
N=12
(
i=0
for f in ../data/run_*.root; do
   ((i=i%N)); ((i++==0)) && wait
    filenameFull=${f##*/}
    filename="`echo $filenameFull | sed 's/.root//'`"
    root -b -q -e 'gROOT->ProcessLine(".L vmm.C"); gROOT->ProcessLine("(new vmm(\"'$filename'\"))->Loop()")' &
done
)
