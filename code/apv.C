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

  // =============================== TDO & distributions ===============================

  // printf("Chain tree: %p\n", fChainSignal);
  printf("Chain n events: %lld\n", fChainSignal->GetEntries());
  Long64_t nentries = fChainSignal->GetEntries();

  // nentries = 200;
  int previousEvt = -1;
  int previousTimeStamp = -1;
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
    if(previousEvt == evt)
      printf("Possible second part of event! \n");
    hevts->Fill(evt);
    hdaqTimeSec->Fill(daqTimeSec);
    hdaqTimeMSec->Fill(daqTimeMicroSec);
    hsrsTrigger->Fill(srsTrigger);

    printf("Evevt parameters: evt %lld, time: %d & %d, timestamp: %d, trigger: %d;", evt, daqTimeSec, daqTimeMicroSec, srsTimeStamp, srsTrigger);
    printf(" Timestamp difference: %d;", srsTimeStamp-previousTimeStamp);
    printf(" Unique timestamp: %llu;", unique_srs_time_stamp(daqTimeSec, daqTimeMicroSec, srsTimeStamp));
    printf("\n");

    printf("  Channels (%lu):", max_q->size());
    /* Per-channel */
    for (int j = 0; j < max_q->size(); j++){
    // printf("Record inside entry: %d\n", j);
      auto readout = mmReadout->at(j);
      if(readout == 'E') //non-mapped channel
        continue;
      auto layer = mmLayer->at(j);
      auto strip = mmStrip->at(j);
      hHitSpot.at(layer)->Fill(strip);
      auto maxQ = max_q->at(j);
      printf(" %d-%d (%d)", layer, strip, maxQ);
      hmaxQ.at(layer)->Fill(maxQ);
      auto maxTime = t_max_q->at(j);
      hmaxQTime.at(layer)->Fill(maxTime);
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
    previousEvt = evt;
    previousTimeStamp = srsTimeStamp;
  }
  
  // straw31_vs_straw30_banana_bcid->Write();
  out->Write();
  out->Close();
}
