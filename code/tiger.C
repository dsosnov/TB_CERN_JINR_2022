#define tiger_cxx
#include "tiger.h"

#include "TH1F.h"
#include "TH2F.h"

map<pair<int, int>, string> detectorNames = {
  {{0, 0}, "scintillator channel"},
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
  {
    printf("Scintillators (det 0):");
    if(detMax.at(0) >= 0) printf(" %d-%d", detMin.at(0), detMax.at(0));
    printf("\n");
    printf("Straws (det 1):");
    if(detMax.at(1) >= 0) printf(" %d-%d", detMin.at(1), detMax.at(1));
    printf("\n");
    for(int i = 0; i < 4; i++){
      printf("MM %d (det %d):", i, i+2);
      if(detMax.at(i+2) >= 0) printf(" %d-%d", detMin.at(i+2), detMax.at(i+2));
      printf("\n");
    }
    printf("Additional straws (det 6):");
    if(detMax.at(6) >= 0) printf(" %d-%d", detMin.at(6), detMax.at(6));
    printf("\n");
    printf("Lemo (det 7):");
    if(detMax.at(7) >= 0) printf(" %d-%d", detMin.at(7), detMax.at(7));
    printf("\n");
  }

  TFile *out = new TFile("../out/out_tiger" + ((runFolder == "") ? "" : runFolder + "-") + file + ending, "RECREATE"); // PATH where to save out_*.root file

  map<pair<int, int>, TH1F*> hTigerProfile, hTigerProfileMapped;
  out->mkdir("tiger_profiles")->cd();
  for(auto &gr: {0})
    for(int t = 0; t < 8; t++){
      pair<int, int> m = {gr, t};
      hTigerProfile.emplace(m, new TH1F(Form("profile_gr%d_t%d", gr, t), Form("%s: profile for gemroc %d tiger %d;channel", file.Data(), gr, t), 64, 0, 64));
      hTigerProfileMapped.emplace(m, new TH1F(Form("profileMapped_gr%d_t%d", gr, t), Form("%s: profile for mapped channels for gemroc %d tiger %d;channel", file.Data(), gr, t), 64, 0, 64));
    }
  out->cd();
  map<pair<int, int>, TH2F*> hTigerChargeToT, hTigerChargeSH;
  map<pair<int, int>, TH2F*> hTigerTimeFine, hTigerFullTime;
  map<pair<int, int>, TH2F*> hTigertCoarse, hTigereCoarse, hTigertFine, hTigereFine;
  map<tuple<int, int, int>, TH2F*> hTigerChargePerTime, hTigereFinePerTime;
  map<tuple<int, int, int>, TH1F*> hTigerFullTimePerChannel;
  map<pair<int, int>, TH1F*> hTigerFullTime1D;
  for(auto &gr: {0}){
    auto grD = out->mkdir(Form("gemroc_%d", gr));
    grD->cd();
    for(int t = 0; t < 8; t++){
      grD->mkdir(Form("tiger_%d", t))->cd();
      pair<int, int> m = {gr, t};
      hTigerChargeToT.emplace(m, new TH2F(Form("charge_tot_gr%d_t%d", gr, t), Form("%s: charge (Time over Threshold mode) for gemroc %d tiger %d;channel;charge", file.Data(), gr, t),
                                          64, 0, 64, 1025, 0, 1025));
      hTigerChargeSH.emplace(m, new TH2F(Form("charge_sh_gr%d_t%d", gr, t), Form("%s: charge (Sample and Hold mode) for gemroc %d tiger %d;channel;charge", file.Data(), gr, t),
                                         64, 0, 64, 1025, 0, 1025));
      hTigerTimeFine.emplace(m, new TH2F(Form("timeFine_gr%d_t%d", gr, t), Form("%s: timeFine for gemroc %d tiger %d;channel;time, ns", file.Data(), gr, t), 64, 0, 64, 4096, 0, 409600));
      hTigertCoarse.emplace(m, new TH2F(Form("tCoarse_gr%d_t%d", gr, t), Form("%s: tCoarse for gemroc %d tiger %d;channel;tCoarse", file.Data(), gr, t), 64, 0, 64, 65536, 0, 65536));
      hTigereCoarse.emplace(m, new TH2F(Form("eCoarse_gr%d_t%d", gr, t), Form("%s: eCoarse for gemroc %d tiger %d;channel;eCoarse", file.Data(), gr, t), 64, 0, 64, 512, 0, 1024));
      hTigertFine.emplace(m, new TH2F(Form("tFine_gr%d_t%d", gr, t), Form("%s: tFine for gemroc %d tiger %d;tFine", file.Data(), gr, t), 64, 0, 64, 512, 0, 1024));
      hTigereFine.emplace(m, new TH2F(Form("eFine_gr%d_t%d", gr, t), Form("%s: eFine for gemroc %d tiger %d;channel;eFine", file.Data(), gr, t), 64, 0, 64, 512, 0, 1024));
      hTigerFullTime.emplace(m, new TH2F(Form("fullTime_gr%d_t%d", gr, t), Form("%s: fullTime for gemroc %d tiger %d;channel;time, s", file.Data(), gr, t), 64, 0, 64, 600, 0, 60));
      hTigerFullTime1D.emplace(m, new TH1F(Form("fullTime1D_gr%d_t%d", gr, t), Form("%s: fullTime for gemroc %d tiger %d;time, s", file.Data(), gr, t), 600, 0, 60));
      for(auto j = 0; j < 64; j++){
        hTigerChargePerTime.emplace(make_tuple(gr, t, j),
                                    new TH2F(Form("ChargePerTime_gr%d_t%d_ch%d", gr, t, j),
                                             Form("%s: Charge for gemroc %d tiger %d channel %d;full time, s; charge%s", file.Data(), gr, t, j, ((energyMode == TigerEnergyMode::SampleAndHold) ? " = 1024 - eFine": "")),
                                             600, 0, 60, 512, 0, 1024));
        hTigereFinePerTime.emplace(make_tuple(gr, t, j),
                                   new TH2F(Form("eFinePerTime_gr%d_t%d_ch%d", gr, t, j), Form("%s: eFine for gemroc %d tiger %d channel %d;full time, s; eFine", file.Data(), gr, t, j), 600, 0, 60, 512, 0, 1024));
        hTigerFullTimePerChannel.emplace(make_tuple(gr, t, j),
                                         new TH1F(Form("FullTimePerChannel_gr%d_t%d_ch%d", gr, t, j), Form("%s: fullTime for gemroc %d tiger %d channel %d;time, s", file.Data(), gr, t, j), 600, 0, 60));
      }
      grD->cd();
    }
    out->cd();
  }

  auto straw_vs_sci = new TH1F("straw_vs_sci", Form("%s: straw vs scint;#Deltat, ns", file.Data()), 1000, -500, 500);
  map<int, TH1F*> straw_vs_mm, mm_vs_sci, mm_vs_sciCoarse;
  for(int i = 0; i < 4; i++){
    straw_vs_mm.emplace(i+2, new TH1F(Form("straw_vs_mm%d", i), Form("%s: straw vs microMegas %d;#Deltat, ns", file.Data(), i), 1000, -500, 500));
    mm_vs_sci.emplace(i+2, new TH1F(Form("mm%d_vs_sci", i), Form("%s: microMegas %d vs scint;#Deltat, ns", file.Data(), i), 1000, -500, 500));
    mm_vs_sciCoarse.emplace(i+2, new TH1F(Form("mm%d_vs_sci_coarse", i), Form("%s: microMegas %d vs scint (coarse time);#DeltaT, ns", file.Data(), i), 160, -500, 500));
  }

  map<int, TH2F*> straw_vs_mm_spatial_corr;
  for(int i = 0; i < 4; i++){
    straw_vs_mm_spatial_corr.emplace(i+2, new TH2F(Form("straw_vs_mm%d_spatial_corr", i), Form("%s: microMegas %d vs straw spatial correaltion;straw ch;MM ch", file.Data(), i),
                                                   detMax.at(1) - detMin.at(1) + 1, detMin.at(1), detMax.at(1) + 1, detMax.at(i+2) - detMin.at(i+2) + 1, detMin.at(i+2), detMax.at(i+2)));
  }

  vector<TH1F*> hprofile;
  vector<TH2F*> hChargeToT, hChargeSH;
  vector<TH2F*> hTimeFine, hFullTime;
  vector<TH2F*> htCoarse, heCoarse, htFine, heFine;
  vector<TH2F*> htCoarse10bit;
  map<int, TH2F*> hDeltaTPrev;
  map<pair<int, int>, TH2F*> hDeltaTPrevPerCharge, hChargePerTime, heFinePerTime;
  map<pair<int, int>, TH1F*> hFullTimePerChannel;
  vector<TH1F*> hFullTime1D;
  for(auto i = 0; i < nDetectorTypes; i++){
    auto d = out->mkdir(Form("det%d", i));
    d->cd();
    hprofile.push_back(new TH1F(Form("profile_det%d", i), Form("%s: profile for detector %d;channel", file.Data(), i), detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1));
    hChargeToT.push_back(new TH2F(Form("charge_tot_det%d", i), Form("%s: charge (Time over Threshold mode) for detector %d;channel;charge", file.Data(), i),
                                  detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 1025, 0, 1025));
    hChargeSH.push_back(new TH2F(Form("charge_sh_det%d", i), Form("%s: charge (Sample and Hold mode) for detector %d;channel;charge", file.Data(), i),
                                 detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 1025, 0, 1025));
    hTimeFine.push_back(new TH2F(Form("timeFine_det%d", i), Form("%s: timeFine for detector %d;channel;time, ns", file.Data(), i), detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 4096, 0, 409600));
    htCoarse.push_back(new TH2F(Form("tCoarse_det%d", i), Form("%s: tCoarse for detector %d;channel;tCoarse", file.Data(), i), detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 65536, 0, 65536));
    heCoarse.push_back(new TH2F(Form("eCoarse_det%d", i), Form("%s: eCoarse for detector %d;channel;eCoarse", file.Data(), i), detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 512, 0, 1024));
    htFine.push_back(new TH2F(Form("tFine_det%d", i), Form("%s: tFine for detector %d;channel;tFine", file.Data(), i), detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 512, 0, 1024));
    heFine.push_back(new TH2F(Form("eFine_det%d", i), Form("%s: eFine for detector %d;channel;eFine", file.Data(), i), detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 512, 0, 1024));
    htCoarse10bit.push_back(new TH2F(Form("tCoarse10bit_det%d", i), Form("%s: last 10 bit of tCoarse for detector %d;channel;tCoarse %% 0x400", file.Data(), i),
                                     detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 512, 0, 1024));
    hFullTime.push_back(new TH2F(Form("fullTime_det%d", i), Form("%s: fullTime for detector %d;channel;time, s", file.Data(), i), detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 600, 0, 60));
    hFullTime1D.push_back(new TH1F(Form("fullTime1D_det%d", i), Form("%s: fullTime for detector %d;time, s", file.Data(), i), 600, 0, 60));

    if(detMax.at(i) >= 0){
      for(auto j = detMin.at(i); j <= detMax.at(i); j++){
        hChargePerTime.emplace(make_pair(i, j),
                               new TH2F(Form("ChargePerTime_det%d_ch%d", i, j),
                                        Form("%s: Charge for detector %d, channel %d;full time, s; charge%s", file.Data(), i, j, ((energyMode == TigerEnergyMode::SampleAndHold) ? " = 1024 - eFine": "")),
                                        600, 0, 60, 512, 0, 1024));
        heFinePerTime.emplace(make_pair(i, j),
                              new TH2F(Form("eFinePerTime_det%d_ch%d", i, j), Form("%s: eFine for detector %d, channel %d%s;full time, s; eFine", file.Data(), i, j,
                                                                                   (detectorNames.count({i,j}) ? (string() + " (" + detectorNames.at({i,j}) + ")").c_str() : "")),
                                       600, 0, 60, 512, 0, 1024));
        hFullTimePerChannel.emplace(make_pair(i, j),
                                    new TH1F(Form("FullTimePerChannel_det%d_ch%d", i, j), Form("%s: fullTime for detector %d channel %d;time, s", file.Data(), i, j), 600, 0, 60));
      }
    }
    if(i == 0 || i == 7){
      if(detMax.at(i) >= 0){
        d->mkdir("deltaT")->cd();
        hDeltaTPrev.emplace(i, new TH2F(Form("DeltaTPrev_det%d", i), Form("%s: #Delta T for detector %d;channel;#Delta T, #mus", file.Data(), i), detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 10000, 0, 100000));
        for(auto j = detMin.at(i); j <= detMax.at(i); j++){
          hDeltaTPrevPerCharge.emplace(make_pair(i, j),
                                       new TH2F(Form("DeltaTPrevPerCharge_det%d_ch%d", i, j), Form("%s: #Delta T for detector %d, channel %d;#Delta T, #mus; charge", file.Data(), i, j), 5000, 0, 50000, 512, 0, 1024));
        }
        d->cd();
      }
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

  TH2F* DeltaTBetweenPulsers = new TH2F("DeltaTBetweenPulsers", Form("%s: DeltaTBetweenPulsers;time, s;T^{scint}_{50#mus} - T^{scint}_{10ms}, #mus", file.Data()), 600, 0, 60, 10000, 0, 100000);
  map<int, TH2F*> heFineCorr, hNeighborsPerTime;
  map<pair<int,int>, TH2F*> heFinePerTimeCorr;
  out->mkdir("plots_corr_6")->cd();
  for(auto i = 2; i < 5; i++){
    hNeighborsPerTime.emplace(i, new TH2F(Form("hNeighborsPerTime_det%d", i), Form("%s: N neighbors per time for detector %d;full time, s; N neighbors", file.Data(), i),
                                          600, 0, 60, detMax.at(i) - detMin.at(i) + 1, 0, detMax.at(i) - detMin.at(i) + 1));
    heFineCorr.emplace(i, new TH2F(Form("eFine_det%d_corr", i), Form("%s: eFine for detector %d corellated with ship straw;channel;eFine", file.Data(), i),
                                     detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 512, 0, 1024));
    for(auto j = detMin.at(i); j <= detMax.at(i); j++){
      heFinePerTimeCorr.emplace(make_pair(i, j),
                                new TH2F(Form("eFinePerTime_corr_det%d_ch%d", i, j), Form("%s: eFine for detector %d, channel %d%s correllated with ship straw;full time, s; eFine", file.Data(), i, j,
                                                                                     (detectorNames.count({i,j}) ? (string() + " (" + detectorNames.at({i,j}) + ")").c_str() : "")),
                                         60, 0, 60, 512, 0, 1024));
    }
  }
  out->cd();

  map<int, TH1F*> hSciTimeToDet, hSciTimeToDetCoarse;
  for(int i = 1; i < 7; i++){
    hSciTimeToDet.emplace(i, new TH1F(Form("sci_vs_det%d", i), Form("%s: T_{scint} - T_{det %d};#Deltat, ns", file.Data(), i), 1000, -500, 500));
    hSciTimeToDetCoarse.emplace(i, new TH1F(Form("sci_vs_det%d_coarse", i), Form("%s: T_{scint} - T_{det %d} (coarse time);#DeltaT, ns", file.Data(), i), 160, -500, 500));
  }
  map<int, TH2F*> hShipRT;
  for(int i = 2; i <= 4; i++){
    hShipRT.emplace(i, new TH2F(Form("ship_rt_vs_mm%d", i-2), Form("%s: RT for SHiP straw and MM %d;strip;#DeltaT, ns", file.Data(), i-2),
                                detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i), 2000, -1000, 1000));
  }

  
  
  Long64_t nentries = fChain->GetEntries();
  if(n > 0 && nentries > n)
    nentries = n;
  
  // =============================== CORRELATION FINDING ===============================
  Long64_t timeWindowNS = 1E4; // ns
  Long64_t firstHitInWindow = 0;
  tigerHitTL hitMain, hitSecondary, hitFirst;
  for (Long64_t jentry = 0; jentry < nentries; jentry++) // You can remove "/ 10" and use the whole dataset
  {
    if (!(jentry % 100000)){
      std::cout << "Entry " << jentry << "\t of \t" << nentries << "\n";
      // if(jentry>0) break;
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
    auto timeSinceStart = timeDifferenceFineNS(hitMain, hitFirst)  *  1E-9;
    if(ENERGY_CUTS){
      if(fchD == 0 && fchM == 0){
        if(energyMode == TigerEnergyMode::SampleAndHold)
          if(charge < 640 && charge > 590) // May be, 674 and 590
            continue;
      } else if(fchD == 0 && fchM == 4){
        if(energyMode == TigerEnergyMode::SampleAndHold)
          if(charge < 870) // 590
            continue;
      } else if(fchD == 0 && fchM == 5){
        if(energyMode == TigerEnergyMode::SampleAndHold)
          if(charge < 600)
            continue;
      }
    }
    { // per-tiger histograms
      hTigerProfile.at({hitMain.gemrocID, hitMain.tigerID})->Fill(hitMain.channelID);
      if(fchD >=0)
        hTigerProfileMapped.at({hitMain.gemrocID, hitMain.tigerID})->Fill(hitMain.channelID);
      hTigerChargeToT.at({hitMain.gemrocID, hitMain.tigerID})->Fill(hitMain.channelID, hitMain.chargeToT());
      hTigerChargeSH.at({hitMain.gemrocID, hitMain.tigerID})->Fill(hitMain.channelID, hitMain.chargeSH());
      hTigerTimeFine.at({hitMain.gemrocID, hitMain.tigerID})->Fill(hitMain.channelID, hitMain.timeFine());
      hTigertCoarse.at({hitMain.gemrocID, hitMain.tigerID})->Fill(hitMain.channelID, hitMain.tCoarse);
      hTigereCoarse.at({hitMain.gemrocID, hitMain.tigerID})->Fill(hitMain.channelID, hitMain.eCoarse);
      hTigertFine.at({hitMain.gemrocID, hitMain.tigerID})->Fill(hitMain.channelID, hitMain.tFine);
      hTigereFine.at({hitMain.gemrocID, hitMain.tigerID})->Fill(hitMain.channelID, hitMain.eFine);
      hTigerFullTime.at({hitMain.gemrocID, hitMain.tigerID})->Fill(hitMain.channelID, timeSinceStart);
      hTigerFullTime1D.at({hitMain.gemrocID, hitMain.tigerID})->Fill(timeSinceStart);
      hTigerChargePerTime.at({hitMain.gemrocID, hitMain.tigerID, hitMain.channelID})->Fill(timeSinceStart, hitMain.charge(energyMode));
      hTigereFinePerTime.at({hitMain.gemrocID, hitMain.tigerID, hitMain.channelID})->Fill(timeSinceStart, hitMain.eFine);        
      hTigerFullTimePerChannel.at({hitMain.gemrocID, hitMain.tigerID, hitMain.channelID})->Fill(timeSinceStart);
    }
    if (fchD < 0) continue; // unmapped channels
    if (fchD < nDetectorTypes){ // per-detector histograms
      hprofile.at(fchD)->Fill(fchM);
      hChargeToT.at(fchD)->Fill(fchM, hitMain.chargeToT());
      hChargeSH.at(fchD)->Fill(fchM, hitMain.chargeSH());
      hTimeFine.at(fchD)->Fill(fchM, hitMain.timeFine());
      htCoarse.at(fchD)->Fill(fchM, hitMain.tCoarse);
      heCoarse.at(fchD)->Fill(fchM, hitMain.eCoarse);
      htFine.at(fchD)->Fill(fchM, hitMain.tFine);
      heFine.at(fchD)->Fill(fchM, hitMain.eFine);
      htCoarse10bit.at(fchD)->Fill(fchM, hitMain.tCoarse%0x400);
      hFullTime.at(fchD)->Fill(fchM, timeSinceStart);
      hFullTime1D.at(fchD)->Fill(timeSinceStart);
      hFullTimePerChannel.at(fchMapped)->Fill(timeSinceStart);
      hChargePerTime.at(fchMapped)->Fill(timeSinceStart, hitMain.charge(energyMode));
      heFinePerTime.at(fchMapped)->Fill(timeSinceStart, hitMain.eFine);        
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
        prevHit[fchMapped] = hitMain;
      }
    }

    if(fchD == 0 && fchM == 0){ // Scintillator
      map<int, map<int, tigerHitTL>> closestHits;
      for(Long64_t kentry = firstHitInWindow; kentry < nentries; kentry++){
        if(kentry == jentry) continue;
        fChain->GetEntry(kentry);
        updateTigerHitTLCurrent(hitSecondary);
        /* Checking that second hit in maximum time window */
        auto timeDifference = timeDifferenceFineNS(hitMain, hitSecondary);
        auto timeDifferenceAbs = fabs(timeDifference);
        if(timeDifferenceAbs > timeWindowNS){
          if(timeDifference > 0){
            firstHitInWindow++;
            continue;
          } else
            break;
        }
        auto [ffchD, ffchM] = getMapped(hitSecondary);
        if(!closestHits.count(ffchD)) closestHits[ffchD] = {};

        // Searct closest hits for all other channels
        if(!closestHits.at(ffchD).count(ffchM) || timeDifferenceAbs < fabs(timeDifferenceFineNS(hitMain, closestHits.at(ffchD).at(ffchM))))
          closestHits.at(ffchD)[ffchM] = hitSecondary;

      }
      if(closestHits.count(6) && closestHits.at(6).count(0) &&
         fabs(timeDifferenceFineNS(hitMain, closestHits.at(6).at(0))) < 1E3){
        for(auto i = 0; i < 4; i++){
          auto idet = i+2;
          if(!closestHits.count(idet))
            continue;
          for(auto &h: closestHits.at(idet)){
            if(fabs(timeDifferenceFineNS(hitMain, h.second)) > 1E3)
              continue;
            hShipRT.at(idet)->Fill(h.first, timeDifferenceFineNS(hitMain, closestHits.at(6).at(0)));
          }
        }
      }
      for(auto i = 1; i < 7; i++){
        if(!closestHits.count(i))
          continue;
        if(hSciTimeToDet.count(i)){
          for(auto &h: closestHits.at(i))
            hSciTimeToDet.at(i)->Fill(timeDifferenceFineNS(hitMain, h.second));
        }
        if(hSciTimeToDet.count(i)){
          for(auto &h: closestHits.at(i))
            hSciTimeToDetCoarse.at(i)->Fill(timeDifferenceCoarsePS(hitMain, h.second)/1E3);
        }
      }
    }
    else if(fchD == 0 && fchM == 4){ // 50us clocks
      if(prevHit.count({0, 5}))
        DeltaTBetweenPulsers->Fill(timeSinceStart, timeDifferenceFineNS(hitMain, prevHit.at({0, 5})) * 1E-3);      
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
    else if (fchD >= 2 && fchD <= 5){ // All MM ch
      /* Searching for secondary hit */
      /* TODO improve speed: maybe load hits to vector and work with vector in memory */
      optional<pair<Long64_t, tigerHitTL>> closestSci0 = nullopt;
      optional<pair<Long64_t, tigerHitTL>> closestShipStraw = nullopt;
      set<int> neighbors;
      for(Long64_t kentry = firstHitInWindow; kentry < nentries; kentry++){
        if(kentry == jentry) continue;
        fChain->GetEntry(kentry);
        updateTigerHitTLCurrent(hitSecondary);
        /* Checking that second hit in maximum time window */
        auto timeDifference = timeDifferenceFineNS(hitMain, hitSecondary);
        auto timeDifferenceAbs = fabs(timeDifference);
        if(timeDifferenceAbs > timeWindowNS){
          if(timeDifference > 0){
            firstHitInWindow++;
            continue;
          } else
            break;
        }
        auto [ffchD, ffchM] = getMapped(hitSecondary);
        if (ffchD == 1 && ffchM == fchM + 1){
        } else if(ffchD == fchD){ // same MM
          if(timeDifferenceAbs < 500)
            neighbors.emplace(ffchM);
        } else if(ffchD >= 2 && ffchD <= 5){ // other MM
        } else if(ffchD == 0 && ffchM == 0){
          if(!closestSci0 || timeDifferenceAbs < fabs(timeDifferenceFineNS(hitMain, closestSci0.value().second)))
            closestSci0 = {kentry, hitSecondary};
        } else if(ffchD == 0 && ffchM == 3){
        } else if(ffchD == 0 && ffchM == 4){
        } else if(ffchD == 0 && ffchM == 5){
        } else if(ffchD == 6 && ffchM == 0){ // SHiP straw
          if(!closestShipStraw || timeDifferenceAbs < fabs(timeDifferenceFineNS(hitMain, closestShipStraw.value().second)))
            closestShipStraw = {kentry, hitSecondary};
        }
      }
      /* Filling histograms */
      if(closestSci0){
        mm_vs_sci.at(fchD)->Fill(timeDifferenceFineNS(hitMain, closestSci0.value().second));
        mm_vs_sciCoarse.at(fchD)->Fill(timeDifferenceCoarsePS(hitMain, closestSci0.value().second)/1E3);
      }
      if(closestShipStraw && fabs(timeDifferenceFineNS(hitMain, closestShipStraw.value().second)) < 500){
        heFineCorr.at(fchD)->Fill(fchM, hitMain.eFine);
        heFinePerTimeCorr.at(fchMapped)->Fill(timeSinceStart, hitMain.eFine);
      }
      if(neighbors.size()){
        int firstInClaster, lastInClaster;
        for(auto j = fchM+1; j <= detMax.at(fchD); j++){
          if(!neighbors.count(j) && !neighbors.count(j+1)){
            lastInClaster = j-1;
            break;
          }
        }
        for(auto j = fchM-1; j >= detMin.at(fchD); j--){
          if(!neighbors.count(j) && !neighbors.count(j-1)){
            firstInClaster = j+1;
            break;
          }
        }
        hNeighborsPerTime.at(fchD)->Fill(timeSinceStart, lastInClaster - firstInClaster);
      }
    }
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
    //     straw_vs_mm.at(mm.first.first)->Fill(timeDifferenceFineNS(hitMain, mm.second.second));
    //     straw_vs_mm_spatial_corr.at(mm.first.first)->Fill(fchM, mm.first.second);
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

  printf("Finish\n");
}
