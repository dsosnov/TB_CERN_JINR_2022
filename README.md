# CERN (JINR + PNPI) TESTBEAM SPRING 2022

## fast analysis tools for the VMM3 testbeam setup

### СODE folder includes

- vmm.C & .h -- source and header files for the fast analysis
- calibration.C & .h -- source and header files for straw [TDO] to [ns] calibration 

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
1. Copy or link RUN directories to data/tiger folder:
`
ln -s /mnt/data-1/rd51-straw/TestBeam-2022-October/data-tiger/RUN_* data/tiger/
`
2. Move to "code" directory:
`
cd code
`
3. Run "Loop" script on file:
```c++
root l
.L tiger.C
t = new tiger("SubRUN_25_GEMROC_0_TL", "RUN_21", "map-tiger-20221013-2")
t->Loop()
.q
```
or, as oneliner:
`
root -b -q -e 'gROOT->ProcessLine(".L tiger.C"); gROOT->ProcessLine("(new tiger(\"SubRUN_25_GEMROC_0_TL\", \"RUN_21\", \"map-tiger-20221013-2\"))->Loop()")'
`
