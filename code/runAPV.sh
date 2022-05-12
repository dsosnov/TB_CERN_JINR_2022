if [[ ! -f "link.C" ]]; then
  rootcint -f link.C -c -p link.h LinkDef.h
fi
# root -b -q -e 'gROOT->ProcessLine(".L link.C"); gROOT->ProcessLine(".L apv.C"); gROOT->ProcessLine("(new apv(\"fe55_run17-data\"))->Loop()")'
root -b -q -e 'gROOT->ProcessLine(".L link.C"); gROOT->ProcessLine(".L apv.C"); gROOT->ProcessLine("(new apv(\"run56\"))->Loop()")'
