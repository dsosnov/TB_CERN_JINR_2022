#!/usr/bin/env bash
if [[ ! -f "link.C" ]]; then
  rootcint -f link.C -c -p link.h LinkDef.h
fi
root -b -q -e 'gROOT->ProcessLine(".L link.C"); gROOT->ProcessLine(".L vmm.C"); gROOT->ProcessLine(".L evBuilder.C"); gROOT->ProcessLine("(new evBuilder(\"run_0224\"))->Loop()")'
