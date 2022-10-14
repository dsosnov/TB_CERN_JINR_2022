#define tiger_cxx
#include "tiger.h"

#include "TH1F.h"
#include "TH2F.h"

map<pair<int, int>, string> detectorNames = {
  {{0, 0}, "scintillator coincidence"},
  {{0, 4}, "50us pulser"},
  {{0, 5}, "10ms pulser"},
  {{6, 0}, "SHiP straw"},
};

constexpr bool ENERGY_CUTS = false;
void tiger::Loop(unsigned long n)
{
  printf("tiger::Loop()\n");

  if (fChain == 0)
    return;

  vector<int> detMin(nDetectorTypes, -1), detMax(nDetectorTypes, -1);
  set<pair<int,int>> mappedTigers;
  for(auto &s: channelMap){
    if(!mappedTigers.count({get<0>(s.first), get<1>(s.first)}))
      mappedTigers.emplace(make_pair(get<0>(s.first), get<1>(s.first)));
    if(s.second.first < 0 || s.second.first >= nDetectorTypes)
      continue;
    if(detMin.at(s.second.first) < 0 || detMin.at(s.second.first) > s.second.second)
      detMin.at(s.second.first) = s.second.second;
    if(detMax.at(s.second.first) < 0 || detMax.at(s.second.first) < s.second.second)
      detMax.at(s.second.first) = s.second.second;
  }
  printf("Scintillators (det 0): %d-%d\n", detMin.at(0), detMax.at(0));
  printf("Straws (det 1): %d-%d\n", detMin.at(1), detMax.at(1));
  for(int i = 0; i < 4; i++)
    printf("MM %d (det %d): %d-%d\n", i, i+2, detMin.at(i+2), detMax.at(i+2));
  printf("Additional straws (det 6): %d-%d\n", detMin.at(6), detMax.at(6));
  printf("Lemo (det 7): %d-%d\n", detMin.at(7), detMax.at(7));

  TFile *out = new TFile("../out/out_tiger" + ((runFolder == "") ? "" : runFolder + "-") + file + ending, "RECREATE"); // PATH where to save out_*.root file

  map<pair<int, int>, TH1F*> hTigerProfile, hTigerProfileMapped;
  out->mkdir("tiger_profiles")->cd();
  for(auto m: mappedTigers){
    hTigerProfile.emplace(m, new TH1F(Form("profile_gr%d_t%d", m.first, m.second), Form("%s: profile for gemroc %d tiger %d;channel", file.Data(), m.first, m.second), 64, 0, 64));
    hTigerProfileMapped.emplace(m, new TH1F(Form("profileMapped_gr%d_t%d", m.first, m.second), Form("%s: profile for mapped channels for gemroc %d tiger %d;channel", file.Data(), m.first, m.second), 64, 0, 64));
  }
  out->cd();

  auto straw_vs_sci = new TH1F("straw_vs_sci", Form("%s: straw vs scint;#Deltat, ns", file.Data()), 1000, -500, 500);
  vector<TH1F*> straw_vs_mm, mm_vs_sci;
  for(int i = 0; i < 4; i++){
    straw_vs_mm.push_back(new TH1F(Form("straw_vs_mm%d", i), Form("%s: straw vs microMegas %d;#Deltat, ns", file.Data(), i), 1000, -500, 500));
    mm_vs_sci.push_back(new TH1F(Form("mm%d_vs_sci", i), Form("%s: microMegas %d vs scint;#Deltat, ns", file.Data(), i), 1000, -500, 500));
  }

  vector<TH2F*> straw_vs_mm_spatial_corr;
  for(int i = 0; i < 4; i++){
    straw_vs_mm_spatial_corr.push_back(new TH2F(Form("straw_vs_mm%d_spatial_corr", i), Form("%s: microMegas %d vs straw spatial correaltion;straw ch;MM ch", file.Data(), i),
                                                detMax.at(1) - detMin.at(1) + 1, detMin.at(1), detMax.at(1) + 1, detMax.at(i+2) - detMin.at(i+2) + 1, detMin.at(i+2), detMax.at(i+2)));
  }

  vector<TH1F*> hprofile;
  vector<TH2F*> hchargeToT, hchargeSH;
  vector<TH2F*> htimeFine, hFullTime;
  vector<TH2F*> htCoarse, heCoarse, htFine, heFine;
  vector<TH2F*> htCoarse10bit;
  map<int, TH2F*> hDeltaTPrev;
  map<pair<int, int>, TH2F*> hDeltaTPrevPerCharge, hChargePerTime, heFinePerTime;
  map<pair<int, int>, TH1F*> hFullTimePerChannel;
  for(auto i = 0; i < nDetectorTypes; i++){
    out->mkdir(Form("det%d", i))->cd();
    hprofile.push_back(new TH1F(Form("profile_det%d", i), Form("%s: profile for detector %d;channel", file.Data(), i), detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1));
    hchargeToT.push_back(new TH2F(Form("charge_tot_det%d", i), Form("%s: charge (Time over Threshold mode) for detector %d;channel;charge",
                                                                    file.Data(), i), detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 1025, 0, 1025));
    hchargeSH.push_back(new TH2F(Form("charge_sh_det%d", i), Form("%s: charge (Sample and Hold mode) for detector %d;channel;charge",
                                                                  file.Data(), i), detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 1025, 0, 1025));
    htimeFine.push_back(new TH2F(Form("timeFine_det%d", i), Form("%s: timeFine for detector %d;channel;time, ns", file.Data(), i), detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 4096, 0, 409600));

    htCoarse.push_back(new TH2F(Form("tCoarse_det%d", i), Form("%s: tCoarse for detector %d;channel;tCoarse", file.Data(), i), detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 65536, 0, 65536));
    heCoarse.push_back(new TH2F(Form("eCoarse_det%d", i), Form("%s: eCoarse for detector %d;channel;eCoarse", file.Data(), i), detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 1024, 0, 1024));
    htFine.push_back(new TH2F(Form("tFine_det%d", i), Form("%s: eFine for detector %d;channel;tFine", file.Data(), i), detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 1024, 0, 1024));
    heFine.push_back(new TH2F(Form("eFine_det%d", i), Form("%s: eFine for detector %d;channel;eFine", file.Data(), i), detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 1024, 0, 1024));

    htCoarse10bit.push_back(new TH2F(Form("tCoarse10bit_det%d", i), Form("%s: last 10 bit of tCoarse for detector %d;channel;tCoarse %% 0x400",
                                                                         file.Data(), i), detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 1024, 0, 1024));

    hFullTime.push_back(new TH2F(Form("fullTime_det%d", i), Form("%s: fullTime for detector %d;channel;time, s", file.Data(), i), detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 6000, 0, 60));

    if(i == 0 || i == 7)
      hDeltaTPrev.emplace(i, new TH2F(Form("DeltaTPrev_det%d", i), Form("%s: #Delta T for detector %d;channel;#Delta T, #mus", file.Data(), i), detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 10000, 0, 100000));

    if(i == 0 || i == 7)
    for(auto j = detMin.at(i); j <= detMax.at(i); j++)
      hDeltaTPrevPerCharge.emplace(make_pair(i, j),
                                  new TH2F(Form("DeltaTPrevPerCharge_det%d_ch%d", i, j), Form("%s: #Delta T for detector %d, channel %d;#Delta T, #mus; charge", file.Data(), i, j), 5000, 0, 50000, 512, 0, 1024));
    
    for(auto j = detMin.at(i); j <= detMax.at(i); j++){
      hChargePerTime.emplace(make_pair(i, j),
                             new TH2F(Form("ChargePerTime_det%d_ch%d", i, j),
                                      Form("%s: Charge for detector %d, channel %d;full time, s; charge%s", file.Data(), i, j, ((energyMode == TigerEnergyMode::SampleAndHold) ? " = 1024 - eFine": "")),
                                      6000, 0, 60, 512, 0, 1024));
      hFullTimePerChannel.emplace(make_pair(i, j),
                                  new TH1F(Form("hFullTimePerChannel_det%d_ch%d", i, j), Form("%s: fullTime for detector %d channel %d;time, s", file.Data(), i, j), 6000, 0, 60));
      heFinePerTime.emplace(make_pair(i, j),
                                  new TH2F(Form("eFinePerTime_det%d_ch%d", i, j), Form("%s: eFine for detector %d, channel %d;full time, s; eFine", file.Data(), i, j), 6000, 0, 60, 512, 0, 1024));
    }
    
    out->cd();
  }
  map<pair<int, int>, tigerHitTL> prevHit;
  
  // auto straw_rt_dir = out->mkdir("straw_rt");
  // straw_rt_dir->cd();
  // map<int, TH2F*> straw_rt, straw_rt_0;
  // for(auto i = detMin.at(1); i <= detMax.at(1); i++){
  //   straw_rt.emplace(i,
  //                    new TH2F(Form("straw%d_rt", i),
  //                             Form("%s: straw %d v-shape sci ch 60;R, mm;T, ns", file.Data(), i),
  //                             32, -4, 4, 300, -100, 200));
  //   straw_rt_0.emplace(i,
  //                      new TH2F(Form("straw%d_rt_0", i),
  //                               Form("%s: straw %d v-shape sci ch 0;R, mm;T, ns", file.Data(), i),
  //                               32, -4, 4, 300, -100, 200));
  // }
  // out->cd();

  out->mkdir("straw_deltat_corr")->cd();
  map<int, TH1F*> straw_deltat, straw_deltat_0;
  for(auto i = detMin.at(1); i <= detMax.at(1); i++){
    straw_deltat.emplace(i,
                         new TH1F(Form("straw%d_vs_sci60", i),
                                  Form("%s: straw%d_vs_sci60;#Delta t", file.Data(), i), 1000, -500, 500));
    straw_deltat_0.emplace(i,
                           new TH1F(Form("straw%d_vs_sci0", i),
                                    Form("%s: straw%d_vs_sci0;#Delta t", file.Data(), i), 1000, -500, 500));
  }
  out->cd();

  out->mkdir("straw_banana")->cd();
  map<int, TH2F*> straw_banana, straw_banana_0 ;
  for(auto i = detMin.at(1); i < detMax.at(1); i++){
    straw_banana.emplace(i,
                         new TH2F(Form("straw%d-%d_banana", i, i+1),
                                  Form("%s: Time difference between straws %d, %d and sci60;T_{straw%d} - T_{scint}, [ns];T_{straw%d} - T_{scint}, [ns]", file.Data(), i, i+1, i, i+1),
                                  500, -250, 250, 500, -250, 250));
    straw_banana_0.emplace(i,
                           new TH2F(Form("straw%d-%d_banana_0", i, i+1),
                                    Form("%s: Time difference between straws %d, %d and sci0;T_{straw%d} - T_{scint}, [ns];T_{straw%d} - T_{scint}, [ns]", file.Data(), i, i+1, i, i+1),
                                    500, -250, 250, 500, -250, 250));
  }
  out->cd();

  out->mkdir("straw_vs_straw_deltat")->cd();
  map<int, TH1F*> straw_straw;
  for(auto i = detMin.at(1); i < detMax.at(1); i++){
    straw_straw.emplace(i,
                        new TH1F(Form("straw%d_vs_straw%d", i, i+1),
                                 Form("%s: straw%d_vs_straw%d;T_{straw%d} - T_{straw%d}", file.Data(), i, i+1, i, i+1), 1000, -500, 500));
  }
  out->cd();

  TH2F* DeltaTBetweenPulsers = new TH2F("DeltaTBetweenPulsers", Form("%s: DeltaTBetweenPulsers;time, s;T^{scint}_{50#mus} - T^{scint}_{10ms}, #mus", file.Data()), 6000, 0, 60, 10000, 0, 100000);
  
  
  Long64_t nentries = fChain->GetEntries();
  if(n > 0 && nentries > n)
    nentries = n;
  
  // =============================== CORRELATION FINDING ===============================
  Long64_t timeWindowNS = 1E3; // ns
  Long64_t firstHitInWindow = 0;
  tigerHitTL hitMain, hitSecondary, hitFirst;
  for (Long64_t jentry = 0; jentry < nentries; jentry++) // You can remove "/ 10" and use the whole dataset
  {
    if (!(jentry % 100000))
    {
      std::cout << "Entry " << jentry << "\t of \t" << nentries << "\n";
    }
    fChain->GetEntry(jentry);
    updateTigerHitTLCurrent(hitMain);
    if (!jentry){
      printf("First hit: ");
      hitMain.print();
      hitFirst = hitMain;
    }
    auto fchMapped = getMapped(hitMain);
    auto [fchD, fchM] = fchMapped;
    auto charge = hitMain.charge(energyMode);
    if(ENERGY_CUTS){
      if(fchD == 0){
        if(fchM == 4){
          if(energyMode == TigerEnergyMode::SampleAndHold)
            if(charge < 870) // 590
              continue;
        }
        else if(fchM == 5){
          if(energyMode == TigerEnergyMode::SampleAndHold)
            if(charge < 600)
              continue;
        }
      }
    }
    if(hTigerProfile.count({hitMain.gemrocID, hitMain.tigerID})){
      hTigerProfile.at({hitMain.gemrocID, hitMain.tigerID})->Fill(hitMain.channelID);
      if(fchD >=0)
        hTigerProfileMapped.at({hitMain.gemrocID, hitMain.tigerID})->Fill(hitMain.channelID);
    }
    if (fchD >=0 && fchD < nDetectorTypes){
      hprofile.at(fchD)->Fill(fchM);
      hchargeToT.at(fchD)->Fill(fchM, hitMain.chargeToT());
      hchargeSH.at(fchD)->Fill(fchM, hitMain.chargeSH());
      htimeFine.at(fchD)->Fill(fchM, hitMain.timeFine());
      htCoarse.at(fchD)->Fill(fchM, hitMain.tCoarse);
      heCoarse.at(fchD)->Fill(fchM, hitMain.eCoarse);
      htFine.at(fchD)->Fill(fchM, hitMain.tFine);
      heFine.at(fchD)->Fill(fchM, hitMain.eFine);
      htCoarse10bit.at(fchD)->Fill(fchM, hitMain.tCoarse%0x400);
      hFullTime.at(fchD)->Fill(fchM, timeDifferenceFineNS(hitMain, hitFirst) * 1E-9);
      hFullTimePerChannel.at(fchMapped)->Fill(timeDifferenceFineNS(hitMain, hitFirst) * 1E-9);
      if(!prevHit.count(fchMapped)){
        if(fchD == 0 || fchD == 7){
          printf("First hit in detector %d, channel %d (time since start: %lld us) : ", fchD, fchM, static_cast<long long>(timeDifferenceFineNS(hitMain, hitFirst) * 1E-3));
          hitMain.print();
        }
        prevHit.emplace(fchMapped, hitMain);
      }
      else{
        // if(fchD == 0)
        //   printf("D %d, ch %d, Diff to previous: %g\n", fchD, fchM, timeDifferenceFineNS(hitMain, prevHit.at(fchMapped)));
        if(hDeltaTPrev.count(fchD))
          hDeltaTPrev.at(fchD)->Fill(fchM, timeDifferenceFineNS(hitMain, prevHit.at(fchMapped)) * 1E-3);
        if(hDeltaTPrevPerCharge.count(fchMapped))
          hDeltaTPrevPerCharge.at(fchMapped)->Fill(timeDifferenceFineNS(hitMain, prevHit.at(fchMapped)) * 1E-3, hitMain.charge(energyMode));
        hChargePerTime.at(fchMapped)->Fill(timeDifferenceFineNS(hitMain, hitFirst) * 1E-9, hitMain.charge(energyMode));
        heFinePerTime.at(fchMapped)->Fill(timeDifferenceFineNS(hitMain, hitFirst) * 1E-9, hitMain.eFine);        
        prevHit[fchMapped] = hitMain;
      }
    }
    if(fchD == 0){ // Scintillator && clocks
      if(fchM == 4){
        if(prevHit.count({0, 5}))
          DeltaTBetweenPulsers->Fill(timeDifferenceFineNS(hitMain, hitFirst) * 1E-9, timeDifferenceFineNS(hitMain, prevHit.at({0, 5})) * 1E-3);
      }
    }
    // else if (fchD == 1){ // All straw ch
    //   /* Searching for secondary hit */
    //   /* TODO improve speed: maybe load hits to vector and work with vector in memory */
    //   pair<Long64_t, tigerHitTL> closestNeighbor = {-1, tigerHitTL()};
    //   pair<Long64_t, tigerHitTL> closestSci0 = {-1, tigerHitTL()};
    //   pair<Long64_t, tigerHitTL> closestSci60 = {-1, tigerHitTL()};
    //   map<pair<int, int>, pair<Long64_t, tigerHitTL>> closestMM; // key: mapped detector-strip
    //   for(Long64_t kentry = firstHitInWindow; kentry < nentries; kentry++){
    //     fChain->GetEntry(kentry);
    //     updateTigerHitTLCurrent(hitSecondary);
    //     /* Checking that second hit in maximum time window */
    //     if(timeDifferenceFineNS(hitMain, hitSecondary) > timeWindowNS){
    //       firstHitInWindow++;
    //       continue;
    //     } else if(timeDifferenceFineNS(hitSecondary, hitMain) > timeWindowNS){
    //       break;
    //     }
    //     auto [ffchD, ffchM] = getMapped(hitSecondary);
    //     auto timeDifferenceAbs = fabs(timeDifferenceFineNS(hitMain, hitSecondary));
    //     if (ffchD == 1 && ffchM == fchM + 1){
    //       if(closestNeighbor.first < 0 ||
    //          timeDifferenceAbs < fabs(timeDifferenceFineNS(hitMain, closestNeighbor.second)))
    //         closestNeighbor = {kentry, hitSecondary};
    //     } else if(ffchD >= 2 && ffchD <= 5){
    //       if(!closestMM.count({ffchD, ffchM}))
    //         closestMM.emplace(make_pair(ffchD, ffchM), make_pair(kentry, hitSecondary));
    //       else if(timeDifferenceAbs < fabs(timeDifferenceFineNS(hitMain, closestMM.at({ffchD, ffchM}).second)))
    //         closestMM[make_pair(ffchD, ffchM)] = make_pair(kentry, hitSecondary);
    //     } else if(ffchD == 0 && ffchM == 0){
    //       if(closestSci0.first < 0 ||
    //          timeDifferenceAbs < fabs(timeDifferenceFineNS(hitMain, closestSci0.second)))
    //         closestSci0 = {kentry, hitSecondary};
    //     } else if(ffchD == 0 && ffchM == 3){
    //       if(closestSci60.first < 0 ||
    //          timeDifferenceAbs < fabs(timeDifferenceFineNS(hitMain, closestSci60.second)))
    //         closestSci60 = {kentry, hitSecondary};
    //     } else if(ffchD == 0 && ffchM == 4){
    //     }
    //   }
    //   /* Filling histograms */
    //   if(closestSci0.first >= 0){
    //     straw_vs_sci->Fill(timeDifferenceFineNS(hitMain, closestSci0.second));
    //     straw_deltat_0.at(fchM)->Fill(timeDifferenceFineNS(hitMain, closestSci0.second));
    //   }
    //   if(closestSci60.first >= 0){
    //     straw_deltat.at(fchM)->Fill(timeDifferenceFineNS(hitMain, closestSci60.second));
    //   }
    //   if(closestNeighbor.first >= 0){
    //     straw_straw.at(fchM)->Fill(timeDifferenceFineNS(hitMain, closestNeighbor.second));
    //     if(closestSci0.first >= 0)
    //       straw_banana_0.at(fchM)->Fill(timeDifferenceFineNS(hitMain,
    //                                                          closestSci0.second),
    //                                     timeDifferenceFineNS(closestNeighbor.second,
    //                                                          closestSci0.second));
    //     if(closestSci60.first >= 0)
    //       straw_banana.at(fchM)->Fill(timeDifferenceFineNS(hitMain,
    //                                                        closestSci60.second),
    //                                   timeDifferenceFineNS(closestNeighbor.second,
    //                                                        closestSci60.second));
        
    //   }
    //   for(auto &mm: closestMM){
    //     straw_vs_mm.at(mm.first.first-2)->Fill(timeDifferenceFineNS(hitMain, mm.second.second));
    //     straw_vs_mm_spatial_corr.at(mm.first.first-2)->Fill(fchM, mm.first.second);
    //   }
    // }
    // else if (fchD >= 2 && fchD <= 5){ // All MM ch
    //   /* Searching for secondary hit */
    //   /* TODO improve speed: maybe load hits to vector and work with vector in memory */
    //   pair<Long64_t, tigerHitTL> closestSci0 = {-1, tigerHitTL()};
    //   for(Long64_t kentry = firstHitInWindow; kentry < nentries; kentry++){
    //     fChain->GetEntry(kentry);
    //     updateTigerHitTLCurrent(hitSecondary);
    //     /* Checking that second hit in maximum time window */
    //     if(timeDifferenceFineNS(hitMain, hitSecondary) > timeWindowNS){
    //       firstHitInWindow++;
    //       continue;
    //     } else if(timeDifferenceFineNS(hitSecondary, hitMain) > timeWindowNS){
    //       break;
    //     }
    //     auto [ffchD, ffchM] = getMapped(hitSecondary);
    //     auto timeDifferenceAbs = fabs(timeDifferenceFineNS(hitMain, hitSecondary));
    //     if (ffchD == 1 && ffchM == fchM + 1){
    //     } else if(ffchD >= 2 && ffchD <= 5){
    //     } else if(ffchD == 0 && ffchM == 0){
    //       if(closestSci0.first < 0 ||
    //          timeDifferenceAbs < fabs(timeDifferenceFineNS(hitMain, closestSci0.second)))
    //         closestSci0 = {kentry, hitSecondary};
    //     } else if(ffchD == 0 && ffchM == 3){
    //     } else if(ffchD == 0 && ffchM == 4){
    //     }
    //   }
    //   /* Filling histograms */
    //   if(closestSci0.first >= 0){
    //     mm_vs_sci.at(fchD)->Fill(timeDifferenceFineNS(hitMain, closestSci0.second));
    //   }
    // }
    // else if (fchD == 6){ // All straw ch
    //   /* Searching for secondary hit */
    //   /* TODO improve speed: maybe load hits to vector and work with vector in memory */
    //   pair<Long64_t, tigerHitTL> closestNeighbor = {-1, tigerHitTL()};
    //   pair<Long64_t, tigerHitTL> closestSci0 = {-1, tigerHitTL()};
    //   pair<Long64_t, tigerHitTL> closestSci60 = {-1, tigerHitTL()};
    //   map<pair<int, int>, pair<Long64_t, tigerHitTL>> closestMM; // key: mapped detector-strip
    //   for(Long64_t kentry = firstHitInWindow; kentry < nentries; kentry++){
    //     fChain->GetEntry(kentry);
    //     updateTigerHitTLCurrent(hitSecondary);
    //     /* Checking that second hit in maximum time window */
    //     if(timeDifferenceFineNS(hitMain, hitSecondary) > timeWindowNS){
    //       firstHitInWindow++;
    //       continue;
    //     } else if(timeDifferenceFineNS(hitSecondary, hitMain) > timeWindowNS){
    //       break;
    //     }
    //     auto [ffchD, ffchM] = getMapped(hitSecondary);
    //     auto timeDifferenceAbs = fabs(timeDifferenceFineNS(hitMain, hitSecondary));
    //     if(ffchD >= 2 && ffchD <= 5){
    //       if(!closestMM.count({ffchD, ffchM}))
    //         closestMM.emplace(make_pair(ffchD, ffchM), make_pair(kentry, hitSecondary));
    //       else if(timeDifferenceAbs < fabs(timeDifferenceFineNS(hitMain, closestMM.at({ffchD, ffchM}).second)))
    //         closestMM[make_pair(ffchD, ffchM)] = make_pair(kentry, hitSecondary);
    //     } else if(ffchD == 0 && ffchM == 0){
    //       if(closestSci0.first < 0 ||
    //          timeDifferenceAbs < fabs(timeDifferenceFineNS(hitMain, closestSci0.second)))
    //         closestSci0 = {kentry, hitSecondary};
    //     } else if(ffchD == 0 && ffchM == 3){
    //       if(closestSci60.first < 0 ||
    //          timeDifferenceAbs < fabs(timeDifferenceFineNS(hitMain, closestSci60.second)))
    //         closestSci60 = {kentry, hitSecondary};
    //     } else if(ffchD == 0 && ffchM == 4){
    //     }
    //   }
    //   /* Filling histograms */
    //   if(closestSci0.first >= 0){
    //     straw_vs_sci->Fill(timeDifferenceFineNS(hitMain, closestSci0.second));
    //     straw_deltat_0.at(fchM)->Fill(timeDifferenceFineNS(hitMain, closestSci0.second));
    //   }
    //   if(closestSci60.first >= 0){
    //     straw_deltat.at(fchM)->Fill(timeDifferenceFineNS(hitMain, closestSci60.second));
    //   }
    //   if(closestNeighbor.first >= 0){
    //     straw_straw.at(fchM)->Fill(timeDifferenceFineNS(hitMain, closestNeighbor.second));
    //     if(closestSci0.first >= 0)
    //       straw_banana_0.at(fchM)->Fill(timeDifferenceFineNS(hitMain,
    //                                                          closestSci0.second),
    //                                     timeDifferenceFineNS(closestNeighbor.second,
    //                                                          closestSci0.second));
    //     if(closestSci60.first >= 0)
    //       straw_banana.at(fchM)->Fill(timeDifferenceFineNS(hitMain,
    //                                                        closestSci60.second),
    //                                   timeDifferenceFineNS(closestNeighbor.second,
    //                                                        closestSci60.second));
        
    //   }
    //   for(auto &mm: closestMM){
    //     straw_vs_mm.at(mm.first.first-2)->Fill(timeDifferenceFineNS(hitMain, mm.second.second));
    //     straw_vs_mm_spatial_corr.at(mm.first.first-2)->Fill(fchM, mm.first.second);
    //   }
    // }
  }

  // auto straw_rt_normed_dir = out->mkdir("straw_rt_normed");
  // straw_rt_normed_dir->cd();
  // map<int, TH2F*> straw_rt_normed, straw_rt_0_normed ;
  // for(auto &h: straw_rt){
  //   auto hnew = static_cast<TH2F*>(h.second->Clone(Form("straw%d_rt_normed", h.first)));
  //   for(auto i = 1; i <= hnew->GetNbinsX(); i++){
  //     auto integ = hnew->Integral(i, i, 1, hnew->GetNbinsY());
  //     if(!integ) continue;
  //     for(auto j = 1; j <= hnew->GetNbinsY(); j++){
  //       auto c = hnew->GetBinContent(i, j);
  //       auto e = hnew->GetBinError(i, j);
  //       hnew->SetBinContent(i, j, static_cast<float>(c) / static_cast<float>(integ));
  //       hnew->SetBinError(i, j, static_cast<float>(e) / static_cast<float>(integ));
  //     }
  //   }
  //   straw_rt_normed.emplace(h.first, hnew);
  // }
  // for(auto &h: straw_rt_0){
  //   auto hnew = static_cast<TH2F*>(h.second->Clone(Form("straw%d_rt_0_normed", h.first)));
  //   for(auto i = 1; i <= hnew->GetNbinsX(); i++){
  //     auto integ = hnew->Integral(i, i, 1, hnew->GetNbinsY());
  //     if(!integ) continue;
  //     for(auto j = 1; j <= hnew->GetNbinsY(); j++){
  //       auto c = hnew->GetBinContent(i, j);
  //       auto e = hnew->GetBinError(i, j);
  //       hnew->SetBinContent(i, j, static_cast<float>(c) / static_cast<float>(integ));
  //       hnew->SetBinError(i, j, static_cast<float>(e) / static_cast<float>(integ));
  //     }
  //   }
  //   straw_rt_0_normed.emplace(h.first, hnew);
  // }
  // out->cd();


  // auto straw_vs_mm_spatial_corr_normed = static_cast<TH2F*>(straw_vs_mm_spatial_corr->Clone("straw_vs_mm_spatial_corr_normed"));
  // straw_vs_mm_spatial_corr_normed->SetTitle(Form("%s: microMegas vs straw spatial correaltion (normed);straw ch;MM ch", file.Data()));
  // for(auto i = 1; i <= straw_vs_mm_spatial_corr_normed->GetNbinsX(); i++){
  //   auto integ = straw_vs_mm_spatial_corr_normed->Integral(i, i, 1, straw_vs_mm_spatial_corr_normed->GetNbinsY());
  //   if(!integ) continue;
  //   for(auto j = 1; j <= straw_vs_mm_spatial_corr_normed->GetNbinsY(); j++){
  //     auto c = straw_vs_mm_spatial_corr_normed->GetBinContent(i, j);
  //     auto e = straw_vs_mm_spatial_corr_normed->GetBinError(i, j);
  //     straw_vs_mm_spatial_corr_normed->SetBinContent(i, j, static_cast<float>(c) / static_cast<float>(integ));
  //     straw_vs_mm_spatial_corr_normed->SetBinError(i, j, static_cast<float>(e) / static_cast<float>(integ));
  //   }
  // }

  for(auto &i: {0, 7})
    for(auto j = detMin.at(i); j <= detMax.at(i); j++){
      if(prevHit.count({i, j})){
        printf("Last hit in detector %d, channel %d (time since start: %lld us) : ", i, j, static_cast<long long>(timeDifferenceFineNS(prevHit.at({i,j}), hitFirst) * 1E-3));
        prevHit.at({i,j}).print();
      }
    }

  out->Write();
  out->Close();

}
