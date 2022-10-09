#!/bin/env bash

if [[ ! -f "../code/link.C" ]]; then
  rootcint -f ../code/link.C -c -p ../code/link.h ../code/LinkDef.h
fi
# root -b -q -e 'gROOT->ProcessLine(".L ../code/link.C"); gROOT->ProcessLine(".L apv.C"); gROOT->ProcessLine("(new apv(\"fe55_run17-data\"))->Loop()")'
# root -b -q -e 'gROOT->ProcessLine(".L ../code/link.C"); gROOT->ProcessLine(".L apv.C"); gROOT->ProcessLine("(new apv(\"run56\"))->Loop()")'
# root -b -q -e 'gROOT->ProcessLine(".L ../code/link.C"); gROOT->ProcessLine(".L apv.C"); gROOT->ProcessLine("(new apv(\"run67\"))->Loop()")'
# root -b -q -e 'gROOT->ProcessLine(".L ../code/link.C"); gROOT->ProcessLine(".L apv.C"); gROOT->ProcessLine("(new apv(\"run103\"))->Loop()")'
# root -b -q -e 'gROOT->ProcessLine(".L ../code/link.C"); gROOT->ProcessLine(".L apv.C"); gROOT->ProcessLine("(new apv(\"run140\"))->Loop()")'


# pedestal only
# root -b -q -e 'gROOT->ProcessLine(".L ../code/link.C"); gROOT->ProcessLine(".L apv.C"); gROOT->ProcessLine("(new apv(\"run137\"))->Loop()")'
# root -b -q -e 'gROOT->ProcessLine(".L ../code/link.C"); gROOT->ProcessLine(".L apv.C"); gROOT->ProcessLine("(new apv(\"run140\"))->Loop()")'

# root -b -q -e 'gROOT->ProcessLine(".L ../code/link.C"); gROOT->ProcessLine(".L apv.C"); gROOT->ProcessLine("(new apv(\"run166\"))->Loop()")'
# root -b -q -e 'gROOT->ProcessLine(".L ../code/link.C"); gROOT->ProcessLine(".L apv.C"); gROOT->ProcessLine("(new apv(\"run168\"))->Loop()")'

j=""
for i in `echo 1{6{1,2,3,4,5,6,7,8,9},7{0,1,2},8{0,1,2,3}}`; do j=$j' run'$i','; done
echo $j
# root -b -q -e 'gROOT->ProcessLine(".L ../code/link.C"); gROOT->ProcessLine(".L apv.C"); gROOT->ProcessLine("(new apv({\"run166\"}))->Loop()")'
# echo root -b -q -e 'gROOT->ProcessLine(".L ../code/link.C"); gROOT->ProcessLine(".L apv.C"); gROOT->ProcessLine("(new apv({'$j'}))->Loop()")'
echo root -b -q -e 'gROOT->ProcessLine(".L ../code/link.C"); gROOT->ProcessLine(".L apv.C"); gROOT->ProcessLine("(new apv({'$j'}))->Loop()")'
