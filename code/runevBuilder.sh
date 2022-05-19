#!/usr/bin/env bash
root -b -q -e 'gROOT->ProcessLine(".L evBuilder.C"); gROOT->ProcessLine("(new evBuilder(\"run_0227\", \"g3_p25_s100\", \"map-20220515\"))->Loop()")'
