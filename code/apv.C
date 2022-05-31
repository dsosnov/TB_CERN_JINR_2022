#define apv_cxx
#include "apv.h"
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <limits>

map<unsigned long, analysisGeneral::mm2CenterHitParameters> apv::GetCentralHits(unsigned long long fromSec, unsigned long long toSec){
  printf("apv::GetCentralHits(%llu, %llu)\n", fromSec, toSec);

  map<unsigned long, analysisGeneral::mm2CenterHitParameters> outputData = {};

  if(!isChain())
    return outputData;
  Long64_t nentries = fChain->GetEntries();

  unsigned long long previousTimestamp = 0;
  unsigned long long hitsToPrev = 0;
  
  for (auto event = 0; event < nentries; event++){
    Long64_t ientry = LoadTree(event);
    if (ientry < 0) break;
    GetEntry(event);
    // printf("Entry: %d: error - %d\n", event, error);
    if(error)
      continue;
      
    if(fromSec > 0 && daqTimeSec < fromSec)
      continue;
    if(toSec > 0 && toSec >= fromSec && daqTimeSec > toSec)
      break;

    /* Checking if there is any mapped channels */
    bool notErr = false;
    for (auto &r: *mmReadout)
      notErr = notErr || (r != 'E');
    if(!notErr) continue;
    
    unsigned long long currentTimestamp = daqTimeSec * 1E6 + daqTimeMicroSec;
    previousTimestamp = currentTimestamp;

    hits.clear();
    for (int j = 0; j < max_q->size(); j++){
      // printf("Record inside entry: %d\n", j);
      auto readout = mmReadout->at(j);
      if(readout == 'E') //non-mapped channel
        continue;
      auto layer = mmLayer->at(j);
      auto strip = mmStrip->at(j);
      auto maxQ = max_q->at(j);
      
      auto maxTime = t_max_q->at(j);
      
      hits.push_back({layer, strip, maxQ, maxTime, raw_q->at(j)});      
    }

    /* Constructing clasters */
    constructClasters();
    
    auto tracks = constructTracks(clasters);
    bool trackIn2Center = false;
    apvTrack *highestTrack = nullptr;
    for(auto &t: tracks){
      auto x2 = get<2>(getHitsForTrack(t));
      if(x2 > 153 && x2 < 210){
        trackIn2Center = true;
        if(highestTrack == nullptr || highestTrack->maxQ() < t.maxQ())
          highestTrack = &t;
      }
    }

    hitsToPrev++;
    if(!trackIn2Center)
      continue;    

    analysisGeneral::mm2CenterHitParameters hit;
    hit.timeSec = daqTimeSec;
    hit.timeMSec = daqTimeMicroSec;
    hit.nHitsToPrev = hitsToPrev;
    
    hit.approximated = !highestTrack->isX2();
    if(highestTrack->isX2()){
      auto c = highestTrack->getX2Claster();
      hit.stripX = c->center();
      hit.pdo = c->maxQ();
      hit.time = c->maxQTime() * 25;
    } else {
      hit.stripX = get<2>(getHitsForTrack(*highestTrack));
      hit.pdo = highestTrack->maxQ();
      hit.time = (highestTrack->getClasters().begin())->maxQTime() * 25; // highestTrack->maxQTime() * 25;
    }

    // hit.approximated = true;
    // hit.stripX = get<2>(getHitsForTrack(*highestTrack));
    // hit.pdo = highestTrack->maxQ();
    // hit.time = (highestTrack->getClasters().begin())->maxQTime() * 25; // highestTrack->maxQTime() * 25;

    hit.pdoRelative = static_cast<double>(hit.pdo) / 2048;
    hitsToPrev = 0;
    

    outputData.emplace(event, hit);
  }
  return outputData;
}

void apv::LoopSecond(unsigned long long sec){
  printf("apv::LoopSecond(%llu)\n", sec);

  if(!isChain())
    return;
  Long64_t nentries = fChain->GetEntries();

  unsigned long long previousTimestamp = 0;
  
  for (auto event = 0; event < nentries; event++){
    Long64_t ientry = LoadTree(event);
    if (ientry < 0) break;
    GetEntry(event);
    // printf("Entry: %d: error - %d\n", event, error);
    if(error)
      continue;
      
    if(daqTimeSec < sec)
      continue;
    if(daqTimeSec > sec)
      break;

    /* Checking if there is any mapped channels */
    bool notErr = false;
    for (auto &r: *mmReadout)
      notErr = notErr || (r != 'E');
    if(!notErr) continue;
    
    unsigned long long currentTimestamp = daqTimeSec * 1E6 + daqTimeMicroSec;
    previousTimestamp = currentTimestamp;

    hits.clear();
    for (int j = 0; j < max_q->size(); j++){
      // printf("Record inside entry: %d\n", j);
      auto readout = mmReadout->at(j);
      if(readout == 'E') //non-mapped channel
        continue;
      auto layer = mmLayer->at(j);
      auto strip = mmStrip->at(j);
      auto maxQ = max_q->at(j);
      
      auto maxTime = t_max_q->at(j);
      
      hits.push_back({layer, strip, maxQ, maxTime, raw_q->at(j)});      
    }

    /* Constructing clasters */
    constructClasters();
    
    bool clasterInRange0 = false, clasterInRange1 = false, clasterInRange2 = false;
    for(auto &c: clasters){
      if(c.getLayer() == 0 && c.center() >= 124 && c.center() <= 168-16)
        clasterInRange0 = true;
      else if(c.getLayer() == 1 && c.center() >= 134 && c.center() <= 184-16)
        clasterInRange1 = true;
      else if(c.getLayer() == 2 && c.center() >= 150 && c.center() <= 209-16 && (c.center() < 170 || c.center() >173))
        clasterInRange2 = true;
    }

    auto tracks = constructTracks(clasters);
    bool trackIn2Center = false;
    for(auto &t: tracks){
      auto x2 = get<2>(getHitsForTrack(t));
      if(x2 > 153 && x2 < 210){
        trackIn2Center = true;
        break;
      }
    }
    if(trackIn2Center){
      printf("Evevt parameters: evt %lld, time: %d & %d, timestamp: %d, trigger: %d;", evt, daqTimeSec, daqTimeMicroSec, srsTimeStamp, srsTrigger);
      printf("\n");
      // for(auto &c: clasters){
      //   printf("     ");
      //   c.print();
      // }
      for(auto &t: tracks){
        auto x2 = get<2>(getHitsForTrack(t));
        if(x2 > 153 && x2 < 210){
          printf("  Track hits to %g on layer 2: {%g, %.2g} (%lu hits)\n", x2, t.intersect(), t.slope(), t.nClasters());
        }
      }
    }
  }

}

void apv::Loop(unsigned long n)
{
  clasterTree->Reset();

  printf("apv::Loop\n");
  if (isChain())
    fChain->GetEntry(0);

  TFile *out = new TFile("../out/out_apv_" + file + ending, "RECREATE"); // PATH where to save out_*.root file

  vector<TDirectory*> dirs;
  for(auto &i: {0, 1, 2, 3}){
    dirs.push_back(out->mkdir(Form("Layer %d", i)));
  }


  // fast check plots
  auto hevts = make_shared<TH1F>("evt", Form("Run %s: evt", file.Data()), 2500, 0, 2500);
  auto hdaqTimeSec = make_shared<TH1F>("daqTimeSec", Form("Run %s: daqTimeSec", file.Data()), 2500, daqTimeSec, daqTimeSec + 2500);
  auto hdaqTimeMSec = make_shared<TH1F>("daqTimeMSec", Form("Run %s: daqTimeMSec", file.Data()), 1000, 0, 1000E3);
  auto hsrsTrigger = make_shared<TH1F>("srsTrigger", Form("Run %s: srsTrigger", file.Data()), 500, srsTrigger, srsTrigger + 5000);

  auto hdaqTimeDifference = make_shared<TH1F>("hdaqTimeDifference", Form("Run %s: hdaqTimeDifference", file.Data()), 10000, 0, 10000);
  
  auto hClasterShiftBetweenLayers01 = make_shared<TH1F>("hClasterShiftBetweenLayers01", Form("Run %s: hClasterShiftBetweenLayers01", file.Data()), 200, -100, 100);
  auto hClasterShiftBetweenLayers02 = make_shared<TH1F>("hClasterShiftBetweenLayers02", Form("Run %s: hClasterShiftBetweenLayers02", file.Data()), 200, -100, 100);
  auto hClasterShiftBetweenLayers12 = make_shared<TH1F>("hClasterShiftBetweenLayers12", Form("Run %s: hClasterShiftBetweenLayers12", file.Data()), 200, -100, 100);
  
  vector<shared_ptr<TH1F>> hMaxQ, hMaxQTime, hProfile /*, hTriggerShiftByMaxQ*/;
  vector<shared_ptr<TH2F>> hPositionVSMaxQ, hPositionVSMaxQTime /*, hTriggerShiftByMaxQ*/;
  for(auto &i: {0, 1, 2, 3}){
    dirs.at(i)->cd();
    hMaxQ.push_back(make_shared<TH1F>(Form("l%d_maxQ", i), Form("Run %s: l%d_maxQ", file.Data(), i), 2500, 0, 2500));
    hMaxQTime.push_back(make_shared<TH1F>(Form("l%d_maxQTime", i), Form("Run %s: l%d_maxQTime", file.Data(), i), 30, 0, 30*25));
    hProfile.push_back(make_shared<TH1F>(Form("l%d_profile", i), Form("Run %s: l%d_profile", file.Data(), i), 361, 0, 361));
    hPositionVSMaxQ.push_back(make_shared<TH2F>(Form("l%d_hPositionVSMaxQ", i), Form("Run %s: l%d_hPositionVSMaxQ", file.Data(), i), 361, 0, 361, 2500, 0, 2500));
    hPositionVSMaxQTime.push_back(make_shared<TH2F>(Form("l%d_hPositionVSMaxQTime", i), Form("Run %s: l%d_hPositionVSMaxQTime", file.Data(), i), 361, 0, 361, 30, 0, 30*25));
    out->cd();
  }
  // auto hMaxQAll = make_shared<TH1F>(Form("maxQ"), Form("Run %s: maxQ", file.Data()), 2500, 0, 2500);
  // auto hMaxQTimeAll = make_shared<TH1F>(Form("maxQTime"), Form("Run %s: maxQTime", file.Data()), 30, 0, 30*25);

  vector<shared_ptr<TH1F>> hClasterMaxQ, hClasterQ, hClasterPosition, hClasterSize;
  vector<shared_ptr<TH2F>> hClasterPositionVSSize, hClasterPositionVSMaxQ, hClasterPositionVSQ;
  for(auto &i: {0, 1, 2, 3}){
    dirs.at(i)->cd();
    hClasterPosition.push_back(make_shared<TH1F>(Form("l%d_clasterPosition", i), Form("Run %s: l%d_clasterPosition", file.Data(), i), 361, 0, 361));
    hClasterMaxQ.push_back(make_shared<TH1F>(Form("l%d_hClasterMaxQ", i), Form("Run %s: l%d_hClasterMaxQ", file.Data(), i), 4096, 0, 4096));
    hClasterQ.push_back(make_shared<TH1F>(Form("l%d_hClasterQ", i), Form("Run %s: l%d_hClasterQ", file.Data(), i), 10240, 0, 10240));
    hClasterSize.push_back(make_shared<TH1F>(Form("l%d_clasterSize", i), Form("Run %s: l%d_clasterSize", file.Data(), i), 50, 0, 50));
    hClasterPositionVSMaxQ.push_back(make_shared<TH2F>(Form("l%d_clasterPositionVSMaxQ", i), Form("Run %s: l%d_clasterPositionVSMaxQ", file.Data(), i), 361, 0, 361, 4096, 0, 4096));
    hClasterPositionVSQ.push_back(make_shared<TH2F>(Form("l%d_clasterPositionVSQ", i), Form("Run %s: l%d_clasterPositionVSQ", file.Data(), i), 361, 0, 361, 10240, 0, 10240));
    hClasterPositionVSSize.push_back(make_shared<TH2F>(Form("l%d_clasterPositionVSSize", i), Form("Run %s: l%d_clasterPositionVSSize", file.Data(), i), 361, 0, 361, 40, 0, 40));
    out->cd();
  }
  // auto hClasterPositionAll = make_shared<TH1F>(Form("hClasterPosition"), Form("Run %s: hClasterPosition", file.Data()), 361, 0, 361);
  // auto hClasterMaxQAll = make_shared<TH1F>(Form("hClasterQ"), Form("Run %s: hClasterQ", file.Data()), 4096, 0, 4096);
  // auto hClasterQAll = make_shared<TH1F>(Form("hClasterQ"), Form("Run %s: hClasterQ", file.Data()), 10240, 0, 10240);
  // auto hClasterSizeAll = make_shared<TH1F>(Form("hClasterSize"), Form("Run %s: hClasterSize", file.Data()), 50, 0, 50);
  // auto hClasterPositionVSMaxQAll = make_shared<TH2F>(Form("hClasterPositionVSMaxQ"), Form("Run %s: hClasterPositionVSMaxQ", file.Data()), 361, 0, 361, 4096, 0, 4096);
  // auto hClasterPositionVSQAll = make_shared<TH2F>(Form("hClasterPositionVSQ"), Form("Run %s: hClasterPositionVSQ", file.Data()), 361, 0, 361, 10240, 0, 10240);
  // auto hClasterPositionVSSizeAll = make_shared<TH2F>(Form("hClasterPositionVSSize"), Form("Run %s: hClasterPositionVSSize", file.Data()), 361, 0, 361, 40, 0, 40);

  vector<shared_ptr<TH1F>> hPedMeanVal, hPedStdevVal, hPedSigmaVal, hPed;
  for(auto &i: {0, 1, 2, 3}){
    dirs.at(i)->cd();
    hPedMeanVal.push_back(make_shared<TH1F>(Form("l%d_hPedMeanVal", i), Form("Run %s: l%d_hPedMeanVal", file.Data(), i), 361, 0, 361));
    hPedStdevVal.push_back(make_shared<TH1F>(Form("l%d_hPedStdevVal", i), Form("Run %s: l%d_hPedStdevVal", file.Data(), i), 361, 0, 361));
    hPedSigmaVal.push_back(make_shared<TH1F>(Form("l%d_hPedSigmaVal", i), Form("Run %s: l%d_profilehPedSigmaVal", file.Data(), i), 361, 0, 361));
    hPed.push_back(make_shared<TH1F>(Form("l%d_pedestals", i), Form("Run %s: l%d_pedestals", file.Data(), i), 361, 0, 361));
    out->cd();
  }


  ulong nEventsWHitsTwoLayers = 0, nEventsWHitsThreeLayers = 0;
  /* Claster histograms */
  // =============================== TDO & distributions ===============================

  // printf("Chain tree: %p\n", fChain);
  if(isChain()){
    printf("Chain n events: %lld\n", fChain->GetEntries());
    Long64_t nentries = fChain->GetEntries();

    if(n > 0 && nentries > n)
      nentries = n;

    unsigned long long previousTimestamp = 0;
  
    // nentries = 2000;
    for (auto event = 0; event < nentries; event++){
      // for(auto event = 80330; event <  nentries; event++){
      Long64_t ientry = LoadTree(event);
      if (ientry < 0) break;
      GetEntry(event);
      // printf("Entry: %d: error - %d\n", event, error);
      if(error)
        continue;

      /* Checking if there is any mapped channels */
      bool notErr = false;
      for (auto &r: *mmReadout)
        notErr = notErr || (r != 'E');
      if(!notErr) continue;
    
      /* Filling entry histogram */
      hevts->Fill(evt);
      hdaqTimeSec->Fill(daqTimeSec);
      hdaqTimeMSec->Fill(daqTimeMicroSec);
      hsrsTrigger->Fill(srsTrigger);

      unsigned long long currentTimestamp = daqTimeSec * 1E6 + daqTimeMicroSec;
      if(previousTimestamp > 0)
        hdaqTimeDifference->Fill(currentTimestamp - previousTimestamp);
      previousTimestamp = currentTimestamp;

      printf("Evevt parameters: evt %lld, time: %d & %d, timestamp: %d, trigger: %d;", evt, daqTimeSec, daqTimeMicroSec, srsTimeStamp, srsTrigger);
      // printf(" Unique timestamp: %llu;", unique_srs_time_stamp(daqTimeSec, daqTimeMicroSec, srsTimeStamp));
      printf("\n");

      // printf("  Channels (%lu):", max_q->size());
      /* Per-channel */
      hits.clear();
      for (int j = 0; j < max_q->size(); j++){
        // printf("Record inside entry: %d\n", j);
        auto readout = mmReadout->at(j);
        if(readout == 'E') //non-mapped channel
          continue;
        auto layer = mmLayer->at(j);
        auto strip = mmStrip->at(j);
        hProfile.at(layer)->Fill(strip);
        auto maxQ = max_q->at(j);
        hMaxQ.at(layer)->Fill(maxQ);
      
        // printf(" %d-%d (%d)", layer, strip, maxQ);

        auto maxTime = t_max_q->at(j);
        hMaxQTime.at(layer)->Fill(maxTime*25);

        hPositionVSMaxQ.at(layer)->Fill(strip, maxQ);
        hPositionVSMaxQTime.at(layer)->Fill(strip, maxTime*25);
      
        hits.push_back({layer, strip, maxQ, maxTime, raw_q->at(j)});
      
        // // printf("Raw q pointer: %p\n", raw_q);
        // int maxADCBin = -1;
        // // printf("raw_q->at(%d).size(): %d\n", j, raw_q->at(j).size());
        // for(uint qBin = 0; qBin < raw_q->at(j).size(); qBin++){
        //   // if(raw_q->at(j).at(qBin) == maxQ){
        //   //   hTriggerShiftByMaxQ.at(layer)->Fill(qBin*25);
        //   // }
        //   if(maxADCBin < 0 || raw_q->at(j).at(qBin) > raw_q->at(j).at(maxADCBin))
        //     maxADCBin = qBin;
        // }
        // hChePerTrigger.at(layer)->Fill(srsTrigger, strip);
        // auto maxADC = raw_q->at(j).at(maxADCBin);
        // // hmaxQFullHistory.at(layer)->Fill(maxADC);
        // // hmaxQTimeFullHistory.at(layer)->Fill(maxADCBin*25);
      }
      // printf("\n");

      /* Constructing clasters */
      constructClasters();
    
      bool clasterInRange0 = false, clasterInRange1 = false, clasterInRange2 = false;
      for(auto &c: clasters){
        if(c.getLayer() == 0 && c.center() >= 124 && c.center() <= 168-16)
          clasterInRange0 = true;
        else if(c.getLayer() == 1 && c.center() >= 134 && c.center() <= 184-16)
          clasterInRange1 = true;
        else if(c.getLayer() == 2 && c.center() >= 150 && c.center() <= 209-16 && (c.center() < 170 || c.center() >173))
          clasterInRange2 = true;
      }
      if(clasterInRange0 && clasterInRange1){
        nEventsWHitsTwoLayers++;
        if(clasterInRange2)
          nEventsWHitsThreeLayers++;
      }

      for(auto &c: clasters){
        printf(" ");
        c.print();
      
        hClasterPosition.at(c.getLayer())->Fill(c.center());
        hClasterMaxQ.at(c.getLayer())->Fill(c.maxQ());
        hClasterQ.at(c.getLayer())->Fill(c.q());
        hClasterSize.at(c.getLayer())->Fill(c.nHits());

        hClasterPositionVSMaxQ.at(c.getLayer())->Fill(c.center(), c.maxQ());;
        hClasterPositionVSQ.at(c.getLayer())->Fill(c.center(), c.q());;
        hClasterPositionVSSize.at(c.getLayer())->Fill(c.center(), c.nHits());
      }

      if(clasters.size() == 3){
        if(clasters.at(0).getLayer() == 0 &&
           clasters.at(1).getLayer() == 1 &&
           clasters.at(2).getLayer() == 2 &&
           (clasters.at(2).center() <= 122 || clasters.at(2).center() > 238)){
          hClasterShiftBetweenLayers01->Fill(clasters.at(0).center() - clasters.at(1).center());
          hClasterShiftBetweenLayers02->Fill(clasters.at(0).center() - clasters.at(2).center());
          hClasterShiftBetweenLayers12->Fill(clasters.at(1).center() - clasters.at(2).center());
        }
      }
      // printf("::N events with hits in two   layers: %lu (of %lld events -> %.2f)\n", nEventsWHitsTwoLayers, event+1, static_cast<double>(nEventsWHitsTwoLayers) / static_cast<double>(event+1));  
      // printf("::N events with hits in three layers: %lu (of %lld events -> %.2f)\n", nEventsWHitsThreeLayers, event+1, static_cast<double>(nEventsWHitsThreeLayers) / static_cast<double>(event+1));  
    
      // clasterTree->Fill();
    }
  
    printf("N events with hits in two   layers to region double-readed wyth mu2E: %lu (of %lld events -> %.2f)\n", nEventsWHitsTwoLayers, nentries, static_cast<double>(nEventsWHitsTwoLayers) / static_cast<double>(nentries));
    printf("N events with hits in three layers to region double-readed wyth mu2E: %lu (of %lld events -> %.2f)\n", nEventsWHitsThreeLayers, nentries, static_cast<double>(nEventsWHitsThreeLayers) / static_cast<double>(nentries));
  }
  
  if(isChainPed()){
    /* Checking pedestals */
    for (int j = 0; j < mmReadoutPed->size(); j++){
      // printf("Record inside entry: %d\n", j);
      auto readout = mmReadoutPed->at(j);
      if(readout == 'E') //non-mapped channel
        continue;
      auto layer = mmLayerPed->at(j);
      auto strip = mmStripPed->at(j);
      auto meanPed = ped_meanPed->at(j);
      auto stdevPed = ped_stdevPed->at(j);
      auto sigmaPed = ped_sigmaPed->at(j);
      auto bin = hPed.at(layer)->GetXaxis()->FindBin(strip);
      hPedMeanVal.at(layer)->Fill(strip, meanPed);
      hPedMeanVal.at(layer)->SetBinError(bin, 0);
      hPedStdevVal.at(layer)->Fill(strip, stdevPed);
      hPedStdevVal.at(layer)->SetBinError(bin, 0);
      hPedSigmaVal.at(layer)->Fill(strip, sigmaPed);
      hPedSigmaVal.at(layer)->SetBinError(bin, 0);
      hPed.at(layer)->SetBinContent(bin, meanPed);
      hPed.at(layer)->SetBinError(bin, stdevPed);
    }
  }

  // straw31_vs_straw30_banana_bcid->Write();
  out->Write();
  // out->Print();
  out->Close();
}
