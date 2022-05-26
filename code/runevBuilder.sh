#!/usr/bin/env bash
root -b -q -e 'gROOT->ProcessLine(".L evBuilder.C"); gROOT->ProcessLine("(new evBuilder(\"run_0260\", \"g3_p25_s100-0&60\", \"map-20220515\"))->Loop()")'
