#!/usr/bin/env bash
root -b -q -e 'gROOT->ProcessLine(".L ../code/vmm.C"); gROOT->ProcessLine("(new vmm(\"run_0057\"))->Loop()")'
