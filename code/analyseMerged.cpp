
#include "apv_cluster.h"
#include "evBuilder.C"
#include "apv.C"

#include "hitsMapper.h"
#include "analyseMerged.h"

bool drLayerFromVMM = true;
void analyseMerged(string runVMM = "0832", string runAPV = "423", bool tight = false, bool timefix = true, string partialfileEnd = ""){
  string tightText = tight ? "_tight" : "";
  string fixTimeText = timefix ? "_timefix" : "";
  auto mergedFileName = TString("../out/runMerged_run_"+runVMM+"_run"+runAPV+tightText+fixTimeText+partialfileEnd+".root");

  auto tVMM = new TChain("vmm");
  tVMM->AddFile(mergedFileName);
  auto vmman = new evBuilder(tVMM, "g1_p25_s100-0&60", "map-20220605");
  vmman->useSyncSignal();
  auto tAPV = new TChain("apv_raw");
  tAPV->AddFile(mergedFileName);
  auto apvan = new apv(tAPV, nullptr);
  apvan->useSyncSignal();

  auto fOut = TFile::Open(Form("../out/_analysed_%s_%s%s.root", runVMM.c_str(), runAPV.c_str(), (tightText+fixTimeText+partialfileEnd).c_str()), "recreate");

  auto hL0Straw = new TH2F("hL0Straw", "hL0Straw; straw; L0 strip", 6, 24, 30, 206, 54, 260);
  auto hL1Straw = new TH2F("hL1Straw", "hL2Straw; straw; L1 strip", 6, 24, 30, 206, 54, 260);
  auto hL2Straw = new TH2F("hL2Straw", "hL2Straw; straw; L2 strip", 6, 24, 30, 206, 54, 260);
  auto estimatedStraw = new TH2F("estimatedStraw", "Estimated track position; straw; estimated position, mm", 6, 24, 30, 500, 0, 50);
  // // 1 in MM strips
  // auto estimatedStraw1 = new TH2F("estimatedStraw1", "Estimated track position; straw; estimated position, mm", 6, 24, 30, 206, 54, 260);

  map<pair<int,int>, TH2F*> dtHists;
  for(auto &straw: {25, 26, 27, 28, 29}){
    for(auto &l: {0, 1, 2}){
      dtHists.emplace(make_pair(l, straw), new TH2F(Form("hL%dStraw%dDT",l, straw), Form("hL%dStraw%dDT; MM strip; dT, ns",l, straw), 206, 54, 260, 200, -100, 100));
    }
  }
  map<int, TH2F*> dtHistsProj, dtHistsProj1;
  for(auto &straw: {25, 26, 27, 28, 29}){
    dtHistsProj.emplace(straw, new TH2F(Form("hProjecteddStraw%dDT",straw), Form("hProjecteddStraw%dDT; position, mm; dT, ns",straw), 500, 0, 50, 200, -100, 100));
    dtHistsProj1.emplace(straw, new TH2F(Form("hProjecteddStraw%dDT_1",straw), Form("hProjecteddStraw%dDT_1; position, strip; dT, ns",straw), 206, 54, 260, 200, -100, 100));
  }

  int nGood = 0;
  auto nEvents = apvan->GetEntries();
  for(auto event = 0; event < nEvents; event++){
    auto dataAPV = apvan->GetCentralHits2ROnlyData(event);
    auto dataVMM = vmman->getHits(event);
    map<int, pair<double,double>> positions;
    if(!dataAPV.hitsPerLayer.at(2).size())
      continue;
    for(auto i = 0; i < 3; i++){
      auto layerData = dataAPV.hitsPerLayer.at(i);
      auto meanPos = getMeanPosition(layerData, i);
      getMeanClusterPositions(layerData, i);
      if(!meanPos)
        continue;
      if(i < 2 || !drLayerFromVMM)
        positions.emplace(i, meanPos.value());
      // printf("Event %d, position for layer %d: %g\n", event, i, meanPos.value());
      TH2F* hist = nullptr;
      switch(i){
        case 0:
          hist = hL0Straw;
          break;
        case 1:
          hist = hL1Straw;
          break;
        case 2:
          hist = hL2Straw;
          break;
      }
      for(auto h: dataVMM){
        if(h.detector != 1)
          continue;
        if(hist != nullptr){
          for(auto &hAPV: layerData)
            hist->Fill(h.strip, hAPV.first);
        }
        for(auto &hAPV: layerData){
          // printf("Filling %d, %d\n",i, h.strip);
          dtHists.at(make_pair(i, h.strip))->Fill(hAPV.first, h.timeToScint);
        }
      }
    }
    if(positions.size() < 2)
      continue;
    // printf("Event %d\t", event);

    if(drLayerFromVMM){
      // Adding MM layer from VMM data
      map<int, int> hitsMMVMM;
      for(auto h: dataVMM){
        if(h.detector != 4) continue;
        // printf("MM VMM strip %d\n", h.strip);
        hitsMMVMM.emplace(h.strip, h.pdo);
      }
      // printf("Number of MM VMM strips %lu\n", hitsMMVMM.size());
      auto meanPosMMVMM = getMeanPosition(hitsMMVMM, 2, true);
      getMeanClusterPositions(hitsMMVMM, 2, true);
      if(!meanPosMMVMM)
        continue;
      positions.emplace(2, meanPosMMVMM.value());
    }

    auto trackParam = getEstimatedTrack(positions);
    auto estimatedCoord = estimatePositionInLayer(trackParam, 3);
    printf("Event %d -- estimated: %g\n", event, estimatedCoord);
    for(auto h: dataVMM){
      if(h.detector != 1)
        continue;
      estimatedStraw->Fill(h.strip, stripToCoord(estimatedCoord));
      dtHistsProj.at(h.strip)->Fill(stripToCoord(estimatedCoord), h.timeToScint);
      dtHistsProj1.at(h.strip)->Fill(estimatedCoord, h.timeToScint);
    }
    nGood++;
  }

  for(auto &h: {hL0Straw, hL1Straw, hL2Straw}){
    renormToUnityByY(h);
    renormToUnityByX(h);
  }

  fOut->Write();
  fOut->Close();
  
  printf("nEvents: %u - %u | %d\n", apvan->GetEntries(), vmman->GetEntries(), nGood);
}
