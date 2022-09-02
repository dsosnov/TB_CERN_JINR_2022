#include "../code/tigerTree.h"

void tiger_fullTime_command(){
  auto mainTree = static_cast<TTree*>(gDirectory->Get("tigerTL"));
  // mainTree->Print();
  Int_t    tCoarse; mainTree->SetBranchAddress("tCoarse", &tCoarse);
  Short_t  tFine; mainTree->SetBranchAddress("tFine", &tFine);
  Int_t    frameCount; mainTree->SetBranchAddress("frameCount", &frameCount);
  Long64_t frameCountLoops; mainTree->SetBranchAddress("frameCountLoops", &frameCountLoops);
  mainTree->GetEntry(0);
  printf("tiger_fullTime(tCoarse,tFine,frameCount,frameCountLoops, %d,%d,%d,%lld)\n", tCoarse,tFine,frameCount,frameCountLoops);
  printf("tiger_fullTime_auto(tCoarse,tFine,frameCount,frameCountLoops)\n");
}

double tiger_fullTime(
  Int_t    tCoarse,
  Short_t  tFine,
  Int_t    frameCount,
  Long64_t frameCountLoops,
  Int_t    tCoarseFirst,
  Short_t  tFineFirst,
  Int_t    frameCountFirst,
  Long64_t frameCountLoopsFirst
  ){
  tigerHitTL currHit, firstHit;
  firstHit.tCoarse = tCoarseFirst;
  firstHit.tFine = tFineFirst;
  firstHit.frameCount = frameCountFirst;
  firstHit.frameCountLoops = frameCountLoopsFirst;
  currHit.tCoarse = tCoarse;
  currHit.tFine = tFine;
  currHit.frameCount = frameCount;
  currHit.frameCountLoops = frameCountLoops;
  double fullTime = timeDifferenceFineNS(currHit, firstHit);
  return fullTime;
}

double tiger_fullTime_auto(
  Int_t    tCoarse,
  Short_t  tFine,
  Int_t    frameCount,
  Long64_t frameCountLoops){
  static Int_t    tCoarseFirst = 0;
  static Short_t  tFineFirst = 0;
  static Int_t    frameCountFirst = 0;
  static Long64_t frameCountLoopsFirst = 0;
  if(!tCoarseFirst && !tFineFirst && !frameCountFirst && !frameCountLoopsFirst){
    auto mainTree = static_cast<TTree*>(gDirectory->Get("tigerTL"));
    mainTree->SetBranchAddress("tCoarse", &tCoarseFirst);
    mainTree->SetBranchAddress("tFine", &tFineFirst);
    mainTree->SetBranchAddress("frameCount", &frameCountFirst);
    mainTree->SetBranchAddress("frameCountLoops", &frameCountLoopsFirst);
    mainTree->GetEntry(0);
    mainTree->ResetBranchAddresses();
  }
  return tiger_fullTime(
    tCoarse,
    tFine,
    frameCount,
    frameCountLoops,
    tCoarseFirst,
    tFineFirst,
    frameCountFirst,
    frameCountLoopsFirst
    );
}
