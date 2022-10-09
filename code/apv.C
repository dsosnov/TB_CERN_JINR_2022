#ifndef apv_cxx
#define apv_cxx
#include "apv.h"
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <limits>

apv::doubleReadoutHits apv::GetCentralHits2ROnlyData(unsigned long long event){
    Long64_t ientry = LoadTree(event);

    set<unsigned int> channelsAPVPulser = {};

    doubleReadoutHits hit;
    hit.sync = false;
    hit.timeSec = 0;
    hit.timeMSec = 0;
    hit.timeSrs = 0;
    hit.hits = {};
    hit.hitsSync = {};
    hit.hitsPerLayer = {};
    
    if (ientry < 0)
      return hit;
    GetEntry(event);
    if(error)
      return hit;
    
  /* Checking if there is any mapped channels */
    bool notErr = false;
    for (auto &r: *mmReadout)
      notErr = notErr || (r != 'E');
    if(!syncSignal && !notErr)
      return {};

    hit.timeSec = static_cast<unsigned int>(daqTimeSec);
    hit.timeMSec = static_cast<unsigned int>(daqTimeMicroSec);
    hit.timeSrs = srsTimeStamp;
    
    unsigned long long currentTimestamp = daqTimeSec * 1E6 + daqTimeMicroSec;

    hit.hits.clear();
    hit.hits.clear();
    hit.hitsSync.clear();
    channelsAPVPulser.clear();
    hit.hitsPerLayer.clear();
    for(auto i = 0; i < nAPVLayers; i++){
      hit.hitsPerLayer.emplace(i, map<int, int>());
    }
    for (int j = 0; j < max_q->size(); j++){
      // printf("Record inside entry: %d\n", j);
      auto maxQ = max_q->at(j);
      auto maxTime = t_max_q->at(j);
      if (syncSignal){
        auto chip = srsChip->at(j);
        auto chan = srsChan->at(j);
        if(chip == pulserAPV){
          channelsAPVPulser.emplace(chan);
          hit.hitsSync.push_back({5, static_cast<int>(chan), maxQ, maxTime, raw_q->at(j)});
        }
      }
      auto readout = mmReadout->at(j);
      if(readout == 'E') //non-mapped channel
        continue;
      auto layer = mmLayer->at(j);
      auto strip = mmStrip->at(j);

      hits.push_back({layer, strip, maxQ, maxTime, raw_q->at(j)});      
      if(layer == layerDoubleReadout && strip > 153 && strip < 210)
        hit.hits.push_back({layer, strip, maxQ, maxTime, raw_q->at(j)});
      hit.hitsPerLayer.at(layer).emplace(strip, maxQ);
    }

    hit.sync = channelsAPVPulser.size() >= 127;
    if(!hit.sync && (hit.hits.size() == 0))
      return hit;
    return hit;
}
map<unsigned long, apv::doubleReadoutHits> apv::GetCentralHits2ROnly(unsigned long long fromSec, unsigned long long toSec, bool saveOnly){
  printf("apv::GetCentralHits2ROnly(%llu, %llu)\n", fromSec, toSec);

  auto out = new TFile("../out/out_apv_" + file + "_centralHits" + ending, "RECREATE");
  auto outTree = new TTree("apv_event", "apv_event");
  outTree->AutoSave("10000");
  pair<unsigned long, apv::doubleReadoutHits> eventData;
  outTree->Branch("event", &eventData);

  map<unsigned long, apv::doubleReadoutHits> outputData = {};
  vector<apvHit> hitsDR;
  vector<apvHit> hitsSync;

  if(!isChain())
    return outputData;
  Long64_t nentries = fChain->GetEntries();

  set<unsigned int> channelsAPVPulser = {};
  map<int, map<int, int>> hitsPerLayer;
  
  for (auto event = 0; event < nentries; event++){
    auto hit = GetCentralHits2ROnlyData(event);
      
    if(fromSec > 0 && hit.timeSec < fromSec)
      continue;
    if(toSec > 0 && toSec >= fromSec && hit.timeSec > toSec)
      break;

    if(!hit.sync && (hit.hits.size() == 0))
      continue;

    if(!saveOnly)
      outputData.emplace(event, hit);
    eventData = make_pair(event, hit);
    outTree->Fill();
  }
  outTree->Write();
  out->Close();
  return outputData;
}

map<unsigned long, analysisGeneral::mm2CenterHitParameters> apv::GetCentralHits(unsigned long long fromSec, unsigned long long toSec, bool saveOnly){
  printf("apv::GetCentralHits(%llu, %llu)\n", fromSec, toSec);

  map<unsigned long, analysisGeneral::mm2CenterHitParameters> outputData = {};

  if(!isChain())
    return outputData;
  Long64_t nentries = fChain->GetEntries();

  long long previousSyncTimestamp = -1;
  unsigned long long hitsToPrev = 0;
  set<unsigned int> channelsAPVPulser = {};
  unsigned long long previousSync = 0;
  map<int, map<int, int>> hitsPerLayer;
  
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
    if(!syncSignal && !notErr) continue;
    
    unsigned long long currentTimestamp = daqTimeSec * 1E6 + daqTimeMicroSec;

    hits.clear();
    channelsAPVPulser.clear();
    hitsPerLayer.clear();
    for(auto i = 0; i < nAPVLayers; i++)
      hitsPerLayer.emplace(i, map<int, int>());
    for (int j = 0; j < max_q->size(); j++){
      // printf("Record inside entry: %d\n", j);
      if (syncSignal){
        auto chip = srsChip->at(j);
        auto chan = srsChan->at(j);
        if(chip == pulserAPV) channelsAPVPulser.emplace(chan);
      }
      auto readout = mmReadout->at(j);
      if(readout == 'E') //non-mapped channel
        continue;
      auto layer = mmLayer->at(j);
      auto strip = mmStrip->at(j);
      auto maxQ = max_q->at(j);
      
      auto maxTime = t_max_q->at(j);
      
      hits.push_back({layer, strip, maxQ, maxTime, raw_q->at(j)});
      hitsPerLayer.at(layer).emplace(strip, maxQ);
    }

    /* Constructing clusters */
    constructClusters();
    
    auto tracks = constructTracks(clusters);
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
    bool isSyncSignal = channelsAPVPulser.size() >= 127;
    if(!isSyncSignal && !trackIn2Center)
      continue;

    analysisGeneral::mm2CenterHitParameters hit;
    hit.timeSec = daqTimeSec;
    hit.timeMSec = daqTimeMicroSec;
    hit.nHitsToPrev = hitsToPrev;
    hit.srsT = srsTimeStamp;
    
    hit.timeSinceSync = (previousSyncTimestamp < 0) ? -1.0 : static_cast<double>(currentTimestamp - previousSyncTimestamp);

    hit.sync = false;
    hit.approx = false;
    hit.signal = highestTrack != nullptr;

    if(isSyncSignal)
      hit.sync = true;

    if(highestTrack != nullptr){
      hit.approx = !highestTrack->isX2();
      if(highestTrack->isX2()){
        auto c = highestTrack->getX2Cluster();
        hit.stripX = c->center();
        hit.pdo = c->maxQ();
        hit.time = c->maxQTime() * 25;
        for(auto h: c->getHits()){
          hit.hitsX.emplace(h.strip, h.max_q);
        }
      } else {
        hit.stripX = get<2>(getHitsForTrack(*highestTrack));
        hit.pdo = highestTrack->maxQ();
        hit.time = (highestTrack->getClusters().begin())->maxQTime() * 25; // highestTrack->maxQTime() * 25;
      }
    } else {
      hit.stripX = 0;
      hit.pdo = 0;
      hit.time = 0;
    }

    // hit.approximated = true;
    // hit.stripX = get<2>(getHitsForTrack(*highestTrack));
    // hit.pdo = highestTrack->maxQ();
    // hit.time = (highestTrack->getClusters().begin())->maxQTime() * 25; // highestTrack->maxQTime() * 25;

    hit.pdoRelative = static_cast<double>(hit.pdo) / 4096;
    hitsToPrev = 0;

    hit.previousSync = previousSync;

    hit.hitsPerLayer = hitsPerLayer;

    outputData.emplace(event, hit);
    if(isSyncSignal){
      previousSyncTimestamp = currentTimestamp;
      previousSync = ientry;
    }

  }
  return outputData;
}

void apv::Loop(unsigned long n)
{
  clusterTree->Reset();

  printf("apv::Loop\n");
  if (isChain())
    fChain->GetEntry(0);

  TFile *out = new TFile("../out/out_apv_" + file + ending, "RECREATE"); // PATH where to save out_*.root file

  vector<TDirectory*> dirs;
  for(auto i = 0; i < nAPVLayers; i++){
    dirs.push_back(out->mkdir(Form("Layer %d", i)));
  }

  // fast check plots
  auto hevts = make_shared<TH1F>("evt", Form("Run %s: evt", file.Data()), 2500, 0, 2500);
  auto hdaqTimeSec = make_shared<TH1F>("daqTimeSec", Form("Run %s: daqTimeSec", file.Data()), 2500, daqTimeSec, daqTimeSec + 2500);
  auto hdaqTimeMSec = make_shared<TH1F>("daqTimeMSec", Form("Run %s: daqTimeMSec", file.Data()), 1000, 0, 1000E3);
  auto hsrsTrigger = make_shared<TH1F>("srsTrigger", Form("Run %s: srsTrigger", file.Data()), 500, srsTrigger, srsTrigger + 5000);

  auto hdaqTimeDifference = make_shared<TH1F>("hdaqTimeDifference", Form("Run %s: hdaqTimeDifference; #Delta T_{trig}, #mus", file.Data()), 10000, 0, 100000);
  auto hdaqTimeDifferenceSync = make_shared<TH1F>("hdaqTimeDifferenceSync", Form("Run %s: hdaqTimeDifferenceSync; #Delta T_{Sync}, #mus", file.Data()), 10000, 0, 100000);
  auto hdaqTimeDifferenceVSMultiplicity2 = make_shared<TH2F>("hdaqTimeDifferenceVSMultiplicity2", Form("Run %s: hdaqTimeDifferenceVSMultiplicity2; #Delta T_{trig}, #mus; N hits", file.Data()), 120, 0, 12000, 129, 0, 129);  
  auto hdaqTimeDifferenceVSTime = make_shared<TH2F>("hdaqTimeDifferenceVSTime", Form("Run %s: hdaqTimeDifferenceVSTime; time, s; #Delta T_{trig}, #mus", file.Data()), 600, 0, 600, 120, 0, 12000);
  auto hNPeriodsBenweenSync = make_shared<TH1F>("hNPeriodsBenweenSync", Form("Run %s: hNPeriodsBenweenSync; N sync periods", file.Data()), 12+500, -1, 11+500);
  
  auto hClusterShiftBetweenLayers01 = make_shared<TH1F>("hClusterShiftBetweenLayers01", Form("Run %s: hClusterShiftBetweenLayers01", file.Data()), 200, -100, 100);
  auto hClusterShiftBetweenLayers02 = make_shared<TH1F>("hClusterShiftBetweenLayers02", Form("Run %s: hClusterShiftBetweenLayers02", file.Data()), 200, -100, 100);
  auto hClusterShiftBetweenLayers12 = make_shared<TH1F>("hClusterShiftBetweenLayers12", Form("Run %s: hClusterShiftBetweenLayers12", file.Data()), 200, -100, 100);
  
  vector<shared_ptr<TH1F>> hMaxQ, hMaxQTime, hProfile /*, hTriggerShiftByMaxQ*/;
  vector<shared_ptr<TH2F>> hPositionVSMaxQ, hPositionVSMaxQTime /*, hTriggerShiftByMaxQ*/;
  for(auto i = 0; i < nAPVLayers; i++){
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

  vector<shared_ptr<TH1F>> hClusterMaxQ, hClusterQ, hClusterPosition, hClusterSize;
  vector<shared_ptr<TH2F>> hClusterPositionVSSize, hClusterPositionVSMaxQ, hClusterPositionVSQ;
  for(auto i = 0; i < nAPVLayers; i++){
    dirs.at(i)->cd();
    hClusterPosition.push_back(make_shared<TH1F>(Form("l%d_clusterPosition", i), Form("Run %s: l%d_clusterPosition", file.Data(), i), 361, 0, 361));
    hClusterMaxQ.push_back(make_shared<TH1F>(Form("l%d_hClusterMaxQ", i), Form("Run %s: l%d_hClusterMaxQ", file.Data(), i), 4096, 0, 4096));
    hClusterQ.push_back(make_shared<TH1F>(Form("l%d_hClusterQ", i), Form("Run %s: l%d_hClusterQ", file.Data(), i), 10240, 0, 10240));
    hClusterSize.push_back(make_shared<TH1F>(Form("l%d_clusterSize", i), Form("Run %s: l%d_clusterSize", file.Data(), i), 50, 0, 50));
    hClusterPositionVSMaxQ.push_back(make_shared<TH2F>(Form("l%d_clusterPositionVSMaxQ", i), Form("Run %s: l%d_clusterPositionVSMaxQ", file.Data(), i), 361, 0, 361, 4096, 0, 4096));
    hClusterPositionVSQ.push_back(make_shared<TH2F>(Form("l%d_clusterPositionVSQ", i), Form("Run %s: l%d_clusterPositionVSQ", file.Data(), i), 361, 0, 361, 10240, 0, 10240));
    hClusterPositionVSSize.push_back(make_shared<TH2F>(Form("l%d_clusterPositionVSSize", i), Form("Run %s: l%d_clusterPositionVSSize", file.Data(), i), 361, 0, 361, 40, 0, 40));
    out->cd();
  }
  // auto hClusterPositionAll = make_shared<TH1F>(Form("hClusterPosition"), Form("Run %s: hClusterPosition", file.Data()), 361, 0, 361);
  // auto hClusterMaxQAll = make_shared<TH1F>(Form("hClusterQ"), Form("Run %s: hClusterQ", file.Data()), 4096, 0, 4096);
  // auto hClusterQAll = make_shared<TH1F>(Form("hClusterQ"), Form("Run %s: hClusterQ", file.Data()), 10240, 0, 10240);
  // auto hClusterSizeAll = make_shared<TH1F>(Form("hClusterSize"), Form("Run %s: hClusterSize", file.Data()), 50, 0, 50);
  // auto hClusterPositionVSMaxQAll = make_shared<TH2F>(Form("hClusterPositionVSMaxQ"), Form("Run %s: hClusterPositionVSMaxQ", file.Data()), 361, 0, 361, 4096, 0, 4096);
  // auto hClusterPositionVSQAll = make_shared<TH2F>(Form("hClusterPositionVSQ"), Form("Run %s: hClusterPositionVSQ", file.Data()), 361, 0, 361, 10240, 0, 10240);
  // auto hClusterPositionVSSizeAll = make_shared<TH2F>(Form("hClusterPositionVSSize"), Form("Run %s: hClusterPositionVSSize", file.Data()), 361, 0, 361, 40, 0, 40);

  vector<shared_ptr<TH1F>> hPedMeanVal, hPedStdevVal, hPedSigmaVal, hPed;
  for(auto i = 0; i < nAPVLayers; i++){
    dirs.at(i)->cd();
    hPedMeanVal.push_back(make_shared<TH1F>(Form("l%d_hPedMeanVal", i), Form("Run %s: l%d_hPedMeanVal", file.Data(), i), 361, 0, 361));
    hPedStdevVal.push_back(make_shared<TH1F>(Form("l%d_hPedStdevVal", i), Form("Run %s: l%d_hPedStdevVal", file.Data(), i), 361, 0, 361));
    hPedSigmaVal.push_back(make_shared<TH1F>(Form("l%d_hPedSigmaVal", i), Form("Run %s: l%d_profilehPedSigmaVal", file.Data(), i), 361, 0, 361));
    hPed.push_back(make_shared<TH1F>(Form("l%d_pedestals", i), Form("Run %s: l%d_pedestals", file.Data(), i), 361, 0, 361));
    out->cd();
  }

  vector<shared_ptr<TH2F>> hClusterPositionXY, hClusterPositionXY_all;
  for(auto &i: {0, 1, 2}){
    dirs.at(i)->cd();
    hClusterPositionXY.push_back(make_shared<TH2F>(Form("l%d_hClusterPositionXY", i), Form("Run %s: l%d_hClusterPositionXY", file.Data(), i), 361, 0, 361, 361, 0, 361));
    hClusterPositionXY_all.push_back(make_shared<TH2F>(Form("l%d_hClusterPositionXY_all", i), Form("Run %s: l%d_hClusterPositionXY_all", file.Data(), i), 361, 0, 361, 361, 0, 361));
    out->cd();
  }

  auto hapvPulserMultiplicity = make_shared<TH1F>("hapvPulserMultiplicity", Form("Run %s: hapvPulserMultiplicity", file.Data()), 129, 0, 129);

  unsigned long nEventsWHitsTwoLayers = 0, nEventsWHitsThreeLayers = 0;
  /* Cluster histograms */
  // =============================== TDO & distributions ===============================

  // printf("Chain tree: %p\n", fChain);
  if(isChain()){
    printf("Chain n events: %lld\n", fChain->GetEntries());
    Long64_t nentries = fChain->GetEntries();

    if(n > 0 && nentries > n)
      nentries = n;

    unsigned long long previousTimestamp = 0;
    unsigned long long previousTimestampSync = 0;
    set<unsigned int> channelsAPVPulser = {};
    unsigned long long firstEventTime = daqTimeSec * 1E6 + daqTimeMicroSec;
  
    // nentries = 2000;
    for (auto event = 0; event < nentries; event++){
      // for(auto event = 80330; event <  nentries; event++){
      Long64_t ientry = LoadTree(event);
      if (ientry < 0) break;
      GetEntry(event);
      // printf("Entry: %d: error - %d\n", event, error);
      if(error)
        continue;

      // /* Checking if there is any mapped channels */
      // bool notErr = false;
      // for (auto &r: *mmReadout)
      //   notErr = notErr || (r != 'E');
      // if(!notErr) continue;
    
      /* Filling entry histogram */
      hevts->Fill(evt);
      hdaqTimeSec->Fill(daqTimeSec);
      hdaqTimeMSec->Fill(daqTimeMicroSec);
      hsrsTrigger->Fill(srsTrigger);

      unsigned long long currentTimestamp = daqTimeSec * 1E6 + daqTimeMicroSec;

      printf("Event parameters: evt %lld, time: %d & %d, timestamp: %d, trigger: %d;", evt, daqTimeSec, daqTimeMicroSec, srsTimeStamp, srsTrigger);
      // printf(" Unique timestamp: %llu;", unique_srs_time_stamp(daqTimeSec, daqTimeMicroSec, srsTimeStamp));
      printf("\n");

      // printf("  Channels (%lu):", max_q->size());
      /* Per-channel */
      hits.clear();
      channelsAPVPulser.clear();
      for (int j = 0; j < max_q->size(); j++){
        // printf("Record inside entry: %d\n", j);
        auto chip = srsChip->at(j);
        auto chan = srsChan->at(j);
        if(chip == pulserAPV) channelsAPVPulser.emplace(chan);

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
        // for(unsigned int qBin = 0; qBin < raw_q->at(j).size(); qBin++){
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
      hapvPulserMultiplicity->Fill(channelsAPVPulser.size());

      /* Constructing clusters */
      constructClusters();
    
      bool clusterInRange0 = false, clusterInRange1 = false, clusterInRange2 = false;
      for(auto &c: clusters){
        if(c.getLayer() == 0 && c.center() >= 124 && c.center() <= 168-16)
          clusterInRange0 = true;
        else if(c.getLayer() == 1 && c.center() >= 134 && c.center() <= 184-16)
          clusterInRange1 = true;
        else if(c.getLayer() == 2 && c.center() >= 150 && c.center() <= 209-16 && (c.center() < 170 || c.center() >173))
          clusterInRange2 = true;
      }
      if(clusterInRange0 && clusterInRange1){
        nEventsWHitsTwoLayers++;
        if(clusterInRange2)
          nEventsWHitsThreeLayers++;
      }

      {
        vector<unsigned long> clustersY;
        for(unsigned long i = 0; i < clusters.size(); i++){
          if(clusters.at(i).getLayer() == 3)
            clustersY.push_back(i);
        }
        for(auto &c:clusters){
          if(c.getLayer() == 3)
            continue;
          vector<unsigned long> clustersY_selected;
          for(auto &i: clustersY){
            hClusterPositionXY_all.at(c.getLayer())->Fill(c.center(), clusters.at(i).center());
            if(c.maxQTime() == clusters.at(i).maxQTime())
              clustersY_selected.push_back(i);
          }
          if(clustersY_selected.size() != 1)
            continue;
          for(auto &i: clustersY_selected){
            hClusterPositionXY.at(c.getLayer())->Fill(c.center(), clusters.at(i).center());
          }
        }
      }

      for(auto &c: clusters){
        printf(" ");
        c.print();
      
        hClusterPosition.at(c.getLayer())->Fill(c.center());
        hClusterMaxQ.at(c.getLayer())->Fill(c.maxQ());
        hClusterQ.at(c.getLayer())->Fill(c.q());
        hClusterSize.at(c.getLayer())->Fill(c.nHits());

        hClusterPositionVSMaxQ.at(c.getLayer())->Fill(c.center(), c.maxQ());;
        hClusterPositionVSQ.at(c.getLayer())->Fill(c.center(), c.q());;
        hClusterPositionVSSize.at(c.getLayer())->Fill(c.center(), c.nHits());
      }

      if(clusters.size() == 3){
        if(clusters.at(0).getLayer() == 0 &&
           clusters.at(1).getLayer() == 1 &&
           clusters.at(2).getLayer() == 2 &&
           (clusters.at(2).center() <= 122 || clusters.at(2).center() > 238)){
          hClusterShiftBetweenLayers01->Fill(clusters.at(0).center() - clusters.at(1).center());
          hClusterShiftBetweenLayers02->Fill(clusters.at(0).center() - clusters.at(2).center());
          hClusterShiftBetweenLayers12->Fill(clusters.at(1).center() - clusters.at(2).center());
        }
      }
      // printf("::N events with hits in two   layers: %lu (of %lld events -> %.2f)\n", nEventsWHitsTwoLayers, event+1, static_cast<double>(nEventsWHitsTwoLayers) / static_cast<double>(event+1));  
      // printf("::N events with hits in three layers: %lu (of %lld events -> %.2f)\n", nEventsWHitsThreeLayers, event+1, static_cast<double>(nEventsWHitsThreeLayers) / static_cast<double>(event+1));  
    
      // clusterTree->Fill();
      if(channelsAPVPulser.size() >= 127){
        if(previousTimestampSync > 0){
          auto timestampSyncDiff = currentTimestamp - previousTimestampSync;
          hdaqTimeDifferenceSync->Fill(timestampSyncDiff);
          unsigned int syncPeriod = 10000; // #mus
          unsigned int maxDiff = 1000; // #mus
          auto nSignalsBetween = static_cast<int>(round(static_cast<double>(timestampSyncDiff) / static_cast<double>(syncPeriod)));
          if(abs(static_cast<long>(timestampSyncDiff) - static_cast<long>(syncPeriod * nSignalsBetween)) > maxDiff){
            printf("For event %lld: Time difference between sync is: %llu, what is about %d signal period and %ld difference\n",
                   ientry, timestampSyncDiff, nSignalsBetween, abs(static_cast<long>(timestampSyncDiff) - static_cast<long>(syncPeriod * nSignalsBetween))
              );
            hNPeriodsBenweenSync->Fill(-1);
          } else
            hNPeriodsBenweenSync->Fill(nSignalsBetween);
          
        }
        previousTimestampSync = currentTimestamp;
      }
      if(previousTimestamp > 0){
        hdaqTimeDifference->Fill(currentTimestamp - previousTimestamp);
        hdaqTimeDifferenceVSMultiplicity2->Fill(currentTimestamp - previousTimestamp, channelsAPVPulser.size());
        hdaqTimeDifferenceVSTime->Fill(static_cast<double>(currentTimestamp - firstEventTime) / 1E6, currentTimestamp - previousTimestamp);
      }
      previousTimestamp = currentTimestamp;
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
#endif
