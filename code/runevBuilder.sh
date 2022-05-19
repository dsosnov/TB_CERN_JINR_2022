#!/usr/bin/env bash
root -b -q -e 'gROOT->ProcessLine(".L evBuilder.C"); gROOT->ProcessLine("(new evBuilder(\"run_0227\"))->Loop()")'
