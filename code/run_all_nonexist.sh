#!/usr/bin/env bash
for f in ../data/run_*.root; do
    filenameFull=${f##*/}
    filename="`echo $filenameFull | sed 's/.root//'`" #"${filenameFull//r/R}"
    outfile="../out/out_${filenameFull}"

    [[ ! -f "$outfile" ]] && root -b -q -e 'gROOT->ProcessLine(".L vmm.C"); gROOT->ProcessLine("(new vmm(\"'$filename'\"))->Loop()")'
done

