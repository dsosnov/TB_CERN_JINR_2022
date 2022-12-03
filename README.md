# CERN (JINR + PNPI + RD51) TESTBEAM SW TOOLS

## Run commands:
### Working with Mu2E (VMM)
1. Copy or link files to data/vmm folder:
`
ln -s /mnt/data-1/rd51-straw/TestBeam-2022-October/data-verso/run_*.root data/vmm/
`
2. Move to "code" directory:
`
cd code
`
3. Run "Loop" script on file:
```c++
root -l
.L evBuilder.C
b = new evBuilder("run_0057", "g1_p25_s100-0&60", "map-20220515")
b->Loop()
.q
```
or, as oneliner:
`
root -b -q -e 'gROOT->ProcessLine(".L evBuilder.C"); gROOT->ProcessLine("(new evBuilder(\"run_0227\", \"g3_p25_s100-0&60\", \"map-20220515\"))->Loop()")'
`
### Working with APV
1. Copy or link files to data/apv folder:
`
ln -s /mnt/data-1/rd51-straw/TestBeam-2022-October/data-mmdaq3/run*.root data/apv/
`
2. Move to "code" directory:
`
cd code
`
3. Create dictionary file for used classes:
`
rootcint -f link.C -c link.h LinkDef.h
`
4. Run "Loop" script on file:
```c++
root -l
.L link.C
.L apv.C
a = new apv("run418")
a->Loop()
.q
```
or, as oneliner:
`
root -b -q -e 'gROOT->ProcessLine(".L link.C"); gROOT->ProcessLine(".L apv.C"); gROOT->ProcessLine("(new apv(\"run418\"))->Loop()")'
`
### Working with TIGER
1. Convert .dat files to root in RUN folder:
`
root -b -q tiger_tree_converter_bin.cxx'("SubRUN_10_GEMROC_0_TL.dat", "RUN_1")'
`
or
`
root -b -q tiger_tree_converter_bin.cxx'("RUN_1")'
`
2. Copy or link RUN directories to data/tiger folder:
`
ln -s /mnt/data-1/rd51-straw/TestBeam-2022-October/data-tiger/RUN_* data/tiger/
`
3. Move to "code" directory:
`
cd code
`
4. Run "Loop" script on file:
```c++
root l
.L tiger.C
t = new tiger("SubRUN_1_GEMROC_0_TL", "RUN_1", "20221020", "20221019", {})
t->Loop()
.q
```
or, as oneliner:
`
root -b -q -e 'gROOT->ProcessLine(".L tiger.C"); gROOT->ProcessLine("(new tiger(\"SubRUN_1_GEMROC_0_TL\", \"RUN_1\", \"20221020\", \"20221020\", {}))->Loop()")'
`
Options:
1. File name without ".root" ending
2. Directory name. If empty, files searched in data/tiger directory itself
3. Mapping file (txt). Can be set as "map-tiger-20221019.txt", "map-tiger-20221019" or even "20221019"
4. Calibration files (txt). Can be set as "20221019" for the files:
   - "../configs/tiger_tfine_calibration-20221019.txt" (calibration for tFine value)
   - "../configs/tiger_efine_calibration_SH-20221019.txt" (calibration for eFine values for hits saved in Sample-And-Hold mode)
   - "../configs/tiger_efine_calibration_ToT-20221019.txt" (calibration for eFine values for hits saved in Time-over-Threshold mode)
   - "../configs/tiger_efine_noise_limits-20221019.txt" (file with eFine for noise for hits in SH mode)
5. Vector set the energy modes for each detector (0 is for Sample-And-Hold, 1 is for Time-over-Threshold)
   First value in vector is "default", for detector not listed in vector. Value at position 1 - for detector 0, etc.
   Default value is 0 (Sample-And-Hold) for all detectors
