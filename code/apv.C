#define apv_cxx
#include "apv.h"
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <limits>

void apv::Loop()
{
  printf("apv::Loop\n");
  if (fChainSignal == nullptr)
    return;
  fChainSignal->GetEntry(0);
  TFile *out = new TFile("../out/out_apv_" + file + ending, "RECREATE"); // PATH where to save out_*.root file

  // fast check plots
  auto hevts = make_shared<TH1F>("hevt", Form("Run %s: hevt", file.Data()), 2500, 0, 2500);
  auto hdaqTimeSec = make_shared<TH1F>("hdaqTimeSec", Form("Run %s: hdaqTimeSec", file.Data()), 2500, daqTimeSec, daqTimeSec + 2500);
  auto hdaqTimeMSec = make_shared<TH1F>("hdaqTimeMSec", Form("Run %s: hdaqTimeMSec", file.Data()), 1000, 0, 1000E3);
  auto hsrsTrigger = make_shared<TH1F>("hsrsTrigger", Form("Run %s: hsrsTrigger", file.Data()), 500, srsTrigger, srsTrigger + 5000);
  vector<shared_ptr<TH1F>> hmaxQ, hmaxQTime, hTriggerShiftByMaxQ, hHitSpot, hmaxQFullHistory, hmaxQTimeFullHistory;
  vector<shared_ptr<TH2F>> hChePerTrigger;
  for(auto &i: {0, 1, 2}){
    hmaxQ.push_back(make_shared<TH1F>(Form("maxQ%d", i), Form("Run %s: maxQ%d", file.Data(), i), 2500, 0, 2500));
    hmaxQTime.push_back(make_shared<TH1F>(Form("maxQTime%d", i), Form("Run %s: maxQTime%d", file.Data(), i), 30, 0, 30));
    hHitSpot.push_back(make_shared<TH1F>(Form("hHitSpotL%d", i), Form("Run %s: hHitSpotL%d", file.Data(), i), 360, 0, 360));
    hTriggerShiftByMaxQ.push_back(make_shared<TH1F>(Form("triggerShiftByMaxQ%d", i), Form("Run %s: triggerShiftByMaxQ%d", file.Data(), i), 192, 0, 192));
    hChePerTrigger.push_back(make_shared<TH2F>(Form("hChePerTrigger%d", i), Form("Run %s: hChePerTrigger%d", file.Data(), i), 500, srsTrigger, srsTrigger + 5000, 360, 0, 360));
    hmaxQFullHistory.push_back(make_shared<TH1F>(Form("maxQFullHistory%d", i), Form("Run %s: maxQFullHistory%d", file.Data(), i), 2500, 0, 2500));
    hmaxQTimeFullHistory.push_back(make_shared<TH1F>(Form("maxQTimeFullHistory%d", i), Form("Run %s: maxQTimeFullHistory%d", file.Data(), i), 192, 0, 192));
  }
  auto hTriggerShiftByMaxQAll = make_shared<TH1F>("triggerShiftByMaxQ", Form("Run %s: triggerShiftByMaxQ", file.Data()), 192, 0, 192);
  auto hmaxQFullHistoryAll = make_shared<TH1F>("maxQFullHistory", Form("Run %s: maxQFullHistory", file.Data()), 2500, 0, 2500);
  auto hmaxQTimeFullHistoryAll = make_shared<TH1F>("maxQTimeFullHistory", Form("Run %s: maxQTimeFullHistory", file.Data()), 192, 0, 192);

  auto hClasterShiftBetweenLayers01 = make_shared<TH1F>("hClasterShiftBetweenLayers01", Form("Run %s: hClasterShiftBetweenLayers01", file.Data()), 200, -100, 100);
  auto hClasterShiftBetweenLayers02 = make_shared<TH1F>("hClasterShiftBetweenLayers02", Form("Run %s: hClasterShiftBetweenLayers02", file.Data()), 200, -100, 100);
  auto hClasterShiftBetweenLayers12 = make_shared<TH1F>("hClasterShiftBetweenLayers12", Form("Run %s: hClasterShiftBetweenLayers12", file.Data()), 200, -100, 100);

  /* Claster histograms */
  vector<shared_ptr<TH1F>> clasterQ, clasterPos;
  for(auto &i: {0, 1, 2}){
    clasterQ.push_back(make_shared<TH1F>(Form("clasterQ%d", i), Form("Run %s: clasterQ%d", file.Data(), i), 10240, 0, 10240));
    clasterPos.push_back(make_shared<TH1F>(Form("clasterPos%d", i), Form("Run %s: clasterPos%d", file.Data(), i), 360, 0, 360));
  }
  auto clasterEnergy = make_shared<TH1F>("clasterQ", Form("Run %s: clasterQ", file.Data()), 10240, 0, 10240);
  
  // =============================== TDO & distributions ===============================

  // printf("Chain tree: %p\n", fChainSignal);
  printf("Chain n events: %lld\n", fChainSignal->GetEntries());
  Long64_t nentries = fChainSignal->GetEntries();

  // nentries = 200;
  for (Long64_t event = 0; event < nentries; event++){
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

    printf("Evevt parameters: evt %lld, time: %d & %d, timestamp: %d, trigger: %d;", evt, daqTimeSec, daqTimeMicroSec, srsTimeStamp, srsTrigger);
    printf(" Unique timestamp: %llu;", unique_srs_time_stamp(daqTimeSec, daqTimeMicroSec, srsTimeStamp));
    printf("\n");

    printf("  Channels (%lu):", max_q->size());
    /* Per-channel */
    vector<apvHit> hits;
    for (int j = 0; j < max_q->size(); j++){
    // printf("Record inside entry: %d\n", j);
      auto readout = mmReadout->at(j);
      if(readout == 'E') //non-mapped channel
        continue;
      auto layer = mmLayer->at(j);
      auto strip = mmStrip->at(j);
      hHitSpot.at(layer)->Fill(strip);
      auto maxQ = max_q->at(j);
      hmaxQ.at(layer)->Fill(maxQ);
      
      printf(" %d-%d (%d)", layer, strip, maxQ);

      auto maxTime = t_max_q->at(j);
      hmaxQTime.at(layer)->Fill(maxTime);
      
      hits.push_back({layer, strip, maxQ, maxTime, raw_q->at(j)});
      
      // printf("Raw q pointer: %p\n", raw_q);
      int maxADCBin = -1;
      // printf("raw_q->at(%d).size(): %d\n", j, raw_q->at(j).size());
      for(uint qBin = 0; qBin < raw_q->at(j).size(); qBin++){
        if(raw_q->at(j).at(qBin) == maxQ){
          hTriggerShiftByMaxQ.at(layer)->Fill(qBin);
          hTriggerShiftByMaxQAll->Fill(qBin);
        }
        if(maxADCBin < 0 || raw_q->at(j).at(qBin) > raw_q->at(j).at(maxADCBin))
          maxADCBin = qBin;
      }
      hChePerTrigger.at(layer)->Fill(srsTrigger, strip);
      auto maxADC = raw_q->at(j).at(maxADCBin);
      hmaxQFullHistory.at(layer)->Fill(maxADC);
      hmaxQFullHistoryAll->Fill(maxADC);
      hmaxQTimeFullHistory.at(layer)->Fill(maxADCBin);
      hmaxQTimeFullHistoryAll->Fill(maxADCBin);
    }
    printf("\n");

    /* Constructing clasters */
    clasters.clear();
    std::sort(hits.begin(), hits.end(), [](auto h1, auto h2){return h1.strip < h2.strip;});
    for(auto &hit: hits)
      addHitToClasters(hit);

    /* Remove small clasters or clasters with small energy*/
    clasters.erase(std::remove_if(clasters.begin(), 
                                  clasters.end(),
                                  [](auto c){return c.nHits() < 3;}),
                   clasters.end());
    clasters.erase(std::remove_if(clasters.begin(), 
                                  clasters.end(),
                                  [](auto c){return c.maxQ() < 300;}),
                   clasters.end());
    
    std::sort(clasters.begin(), clasters.end(), [](auto c1, auto c2){return (c1.layer < c2.layer) || (c1.center() < c2.center());});

    for(auto &c: clasters){
      c.print();
      clasterPos.at(c.layer)->Fill(c.center());
      clasterQ.at(c.layer)->Fill(c.q());
      clasterEnergy->Fill(c.q());
    }

    if(clasters.size() == 3){
      if(clasters.at(0).layer == 0 &&
         clasters.at(1).layer == 1 &&
         clasters.at(2).layer == 2){
        hClasterShiftBetweenLayers01->Fill(clasters.at(0).center() - clasters.at(1).center());
        hClasterShiftBetweenLayers02->Fill(clasters.at(0).center() - clasters.at(2).center());
        hClasterShiftBetweenLayers12->Fill(clasters.at(1).center() - clasters.at(2).center());
      }
    }

  }
  
  // straw31_vs_straw30_banana_bcid->Write();
  out->Write();
  out->Close();
}
