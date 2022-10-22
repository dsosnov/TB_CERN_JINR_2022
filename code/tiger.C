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

// output: ns
int maxTimeDiff(int det1 = -1, int det2 = -1){
  auto detMin = det1 > det2 ? det2 : det1;
  auto detMax = detMin == det2 ? det1 : det2;
  int maxDifference = 1000;
  if(detMin == -1 && detMax == -1) // maximal value
    maxDifference = 1000;
  else if(detMin == 0 && detMax == 0) // scint-scint
    maxDifference = 500;
  else if(detMin == 0 && detMax > 0) // scint - anything
    maxDifference = 1000;
  else if(detMin == 1 && detMax == 1) // straw-straw
    maxDifference = 500;
  else if(detMin == 1 && detMax >= 2 && detMax <= 5) // straw-MM
    maxDifference = 500;
  else if(detMin >= 2 && detMin <= 5 && detMax >= 2 && detMax <= 5) // MM-MM
    maxDifference = 500;
  else // anything else
    maxDifference = 1000;
  return maxDifference;
}

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

  auto out = make_shared<TFile>("../out/out_tiger_" + ((runFolder == "") ? "" : runFolder + "-") + file + ending, "RECREATE"); // PATH where to save out_*.root file

  map<pair<int, int>, shared_ptr<TH1F>> hTigerProfile, hTigerProfileMapped;
  out->mkdir("tiger_profiles")->cd();
  for(auto &gr: {0})
    for(int t = 0; t < 8; t++){
      pair<int, int> m = {gr, t};
      hTigerProfile.emplace(m, make_shared<TH1F>(Form("profile_gr%d_t%d", gr, t), Form("%s: profile for gemroc %d tiger %d;channel", file.Data(), gr, t), 64, 0, 64));
      hTigerProfileMapped.emplace(m, make_shared<TH1F>(Form("profileMapped_gr%d_t%d", gr, t), Form("%s: profile for mapped channels for gemroc %d tiger %d;channel", file.Data(), gr, t), 64, 0, 64));
    }
  out->cd();
  map<pair<int, int>, shared_ptr<TH2F>> hTigerChargeToT, hTigerChargeSH;
  map<pair<int, int>, shared_ptr<TH2F>> hTigerTimeFine, hTigerFullTime;
  map<pair<int, int>, shared_ptr<TH2F>> hTigertCoarse, hTigereCoarse, hTigertFine, hTigereFine;
  map<tuple<int, int, int>, shared_ptr<TH2F>> hTigerChargePerTime, hTigereFinePerTime;
  map<tuple<int, int, int>, shared_ptr<TH1F>> hTigerFullTimePerChannel;
  map<pair<int, int>, shared_ptr<TH1F>> hTigerFullTime1D;
  for(auto &gr: {0}){
    auto grD = out->mkdir(Form("gemroc_%d", gr));
    grD->cd();
    for(int t = 0; t < 8; t++){
      grD->mkdir(Form("tiger_%d", t))->cd();
      pair<int, int> m = {gr, t};
      hTigerChargeToT.emplace(m, make_shared<TH2F>(Form("charge_tot_gr%d_t%d", gr, t), Form("%s: charge (Time over Threshold mode) for gemroc %d tiger %d;channel;charge", file.Data(), gr, t),
                                                   64, 0, 64, 1025, 0, 1025));
      hTigerChargeSH.emplace(m, make_shared<TH2F>(Form("charge_sh_gr%d_t%d", gr, t), Form("%s: charge (Sample and Hold mode) for gemroc %d tiger %d;channel;charge", file.Data(), gr, t),
                                                  64, 0, 64, 1025, 0, 1025));
      hTigerTimeFine.emplace(m, make_shared<TH2F>(Form("timeFine_gr%d_t%d", gr, t), Form("%s: timeFine for gemroc %d tiger %d;channel;time, ns", file.Data(), gr, t), 64, 0, 64, 4096, 0, 409600));
      hTigertCoarse.emplace(m, make_shared<TH2F>(Form("tCoarse_gr%d_t%d", gr, t), Form("%s: tCoarse for gemroc %d tiger %d;channel;tCoarse", file.Data(), gr, t), 64, 0, 64, 65536, 0, 65536));
      hTigereCoarse.emplace(m, make_shared<TH2F>(Form("eCoarse_gr%d_t%d", gr, t), Form("%s: eCoarse for gemroc %d tiger %d;channel;eCoarse", file.Data(), gr, t), 64, 0, 64, 512, 0, 1024));
      hTigertFine.emplace(m, make_shared<TH2F>(Form("tFine_gr%d_t%d", gr, t), Form("%s: tFine for gemroc %d tiger %d;tFine", file.Data(), gr, t), 64, 0, 64, 512, 0, 1024));
      hTigereFine.emplace(m, make_shared<TH2F>(Form("eFine_gr%d_t%d", gr, t), Form("%s: eFine for gemroc %d tiger %d;channel;eFine", file.Data(), gr, t), 64, 0, 64, 512, 0, 1024));
      hTigerFullTime.emplace(m, make_shared<TH2F>(Form("fullTime_gr%d_t%d", gr, t), Form("%s: fullTime for gemroc %d tiger %d;channel;time, s", file.Data(), gr, t), 64, 0, 64, 600, 0, 60));
      hTigerFullTime1D.emplace(m, make_shared<TH1F>(Form("fullTime1D_gr%d_t%d", gr, t), Form("%s: fullTime for gemroc %d tiger %d;time, s", file.Data(), gr, t), 600, 0, 60));
      for(auto j = 0; j < 64; j++){
        hTigerChargePerTime.emplace(make_tuple(gr, t, j),
                                    make_shared<TH2F>(Form("ChargePerTime_gr%d_t%d_ch%d", gr, t, j),
                                                      Form("%s: Charge for gemroc %d tiger %d channel %d;full time, s; charge%s", file.Data(), gr, t, j, ((energyMode == TigerEnergyMode::SampleAndHold) ? " = 1024 - eFine": "")),
                                                      600, 0, 60, 512, 0, 1024));
        hTigereFinePerTime.emplace(make_tuple(gr, t, j),
                                   make_shared<TH2F>(Form("eFinePerTime_gr%d_t%d_ch%d", gr, t, j), Form("%s: eFine for gemroc %d tiger %d channel %d;full time, s; eFine", file.Data(), gr, t, j), 600, 0, 60, 512, 0, 1024));
        hTigerFullTimePerChannel.emplace(make_tuple(gr, t, j),
                                         make_shared<TH1F>(Form("FullTimePerChannel_gr%d_t%d_ch%d", gr, t, j), Form("%s: fullTime for gemroc %d tiger %d channel %d;time, s", file.Data(), gr, t, j), 600, 0, 60));
      }
      grD->cd();
    }
    out->cd();
  }

  auto straw_vs_sci = make_shared<TH1F>("straw_vs_sci", Form("%s: straw vs scint;#Deltat, ns", file.Data()), 1000, -500, 500);
  map<int, shared_ptr<TH1F>> straw_vs_mm, mm_vs_sci, mm_vs_sciCoarse;
  for(int i = 0; i <= 4; i++){
    straw_vs_mm.emplace(i+2, make_shared<TH1F>(Form("straw_vs_mm%d", i), Form("%s: straw vs microMegas %d;#Deltat, ns", file.Data(), i), 1000, -500, 500));
    mm_vs_sci.emplace(i+2, make_shared<TH1F>(Form("mm%d_vs_sci", i), Form("%s: microMegas %d vs scint;#Deltat, ns", file.Data(), i), 1000, -500, 500));
    mm_vs_sciCoarse.emplace(i+2, make_shared<TH1F>(Form("mm%d_vs_sci_coarse", i), Form("%s: microMegas %d vs scint (coarse time);#DeltaT, ns", file.Data(), i), 160, -500, 500));
  }

  map<int, shared_ptr<TH2F>> straw_vs_mm_spatial_corr;
  for(int i = 0; i <= 4; i++){
    if(i+2 == mmLayerY) continue;
    straw_vs_mm_spatial_corr.emplace(i+2, make_shared<TH2F>(Form("straw_vs_mm%d_spatial_corr", i), Form("%s: microMegas %d vs straw spatial correaltion;straw ch;MM ch", file.Data(), i),
                                                            detMax.at(1) - detMin.at(1) + 1, detMin.at(1), detMax.at(1) + 1, detMax.at(i+2) - detMin.at(i+2) + 1, detMin.at(i+2), detMax.at(i+2)));
  }

  vector<shared_ptr<TH1F>> hprofile;
  vector<shared_ptr<TH2F>> hChargeToT, hChargeSH;
  vector<shared_ptr<TH2F>> hTimeFine, hFullTime;
  vector<shared_ptr<TH2F>> htCoarse, heCoarse, htFine, heFine;
  vector<shared_ptr<TH2F>> htCoarse10bit;
  map<int, shared_ptr<TH2F>> hDeltaTPrev;
  map<pair<int, int>, shared_ptr<TH2F>> hDeltaTPrevPerCharge, hChargePerTime, heFinePerTime;
  map<pair<int, int>, shared_ptr<TH1F>> hFullTimePerChannel;
  vector<shared_ptr<TH1F>> hFullTime1D;
  for(auto i = 0; i < nDetectorTypes; i++){
    auto d = out->mkdir(Form("det%d", i));
    d->cd();
    hprofile.push_back(make_shared<TH1F>(Form("profile_det%d", i), Form("%s: profile for detector %d;channel", file.Data(), i), detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1));
    hChargeToT.push_back(make_shared<TH2F>(Form("charge_tot_det%d", i), Form("%s: charge (Time over Threshold mode) for detector %d;channel;charge", file.Data(), i),
                                           detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 1025, 0, 1025));
    hChargeSH.push_back(make_shared<TH2F>(Form("charge_sh_det%d", i), Form("%s: charge (Sample and Hold mode) for detector %d;channel;charge", file.Data(), i),
                                          detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 1025, 0, 1025));
    hTimeFine.push_back(make_shared<TH2F>(Form("timeFine_det%d", i), Form("%s: timeFine for detector %d;channel;time, ns", file.Data(), i), detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 4096, 0, 409600));
    htCoarse.push_back(make_shared<TH2F>(Form("tCoarse_det%d", i), Form("%s: tCoarse for detector %d;channel;tCoarse", file.Data(), i), detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 65536, 0, 65536));
    heCoarse.push_back(make_shared<TH2F>(Form("eCoarse_det%d", i), Form("%s: eCoarse for detector %d;channel;eCoarse", file.Data(), i), detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 512, 0, 1024));
    htFine.push_back(make_shared<TH2F>(Form("tFine_det%d", i), Form("%s: tFine for detector %d;channel;tFine", file.Data(), i), detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 512, 0, 1024));
    heFine.push_back(make_shared<TH2F>(Form("eFine_det%d", i), Form("%s: eFine for detector %d;channel;eFine", file.Data(), i), detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 512, 0, 1024));
    htCoarse10bit.push_back(make_shared<TH2F>(Form("tCoarse10bit_det%d", i), Form("%s: last 10 bit of tCoarse for detector %d;channel;tCoarse %% 0x400", file.Data(), i),
                                              detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 512, 0, 1024));
    hFullTime.push_back(make_shared<TH2F>(Form("fullTime_det%d", i), Form("%s: fullTime for detector %d;channel;time, s", file.Data(), i), detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 600, 0, 60));
    hFullTime1D.push_back(make_shared<TH1F>(Form("fullTime1D_det%d", i), Form("%s: fullTime for detector %d;time, s", file.Data(), i), 600, 0, 60));

    if(detMax.at(i) >= 0){
      for(auto j = detMin.at(i); j <= detMax.at(i); j++){
        hChargePerTime.emplace(make_pair(i, j),
                               make_shared<TH2F>(Form("ChargePerTime_det%d_ch%d", i, j),
                                                 Form("%s: Charge for detector %d, channel %d;full time, s; charge%s", file.Data(), i, j, ((energyMode == TigerEnergyMode::SampleAndHold) ? " = 1024 - eFine": "")),
                                                 600, 0, 60, 512, 0, 1024));
        heFinePerTime.emplace(make_pair(i, j),
                              make_shared<TH2F>(Form("eFinePerTime_det%d_ch%d", i, j), Form("%s: eFine for detector %d, channel %d%s;full time, s; eFine", file.Data(), i, j,
                                                                                            (detectorNames.count({i,j}) ? (string() + " (" + detectorNames.at({i,j}) + ")").c_str() : "")),
                                                600, 0, 60, 512, 0, 1024));
        hFullTimePerChannel.emplace(make_pair(i, j),
                                    make_shared<TH1F>(Form("FullTimePerChannel_det%d_ch%d", i, j), Form("%s: fullTime for detector %d channel %d;time, s", file.Data(), i, j), 600, 0, 60));
      }
    }
    if(i == 0 || i == 7){
      if(detMax.at(i) >= 0){
        d->mkdir("deltaT")->cd();
        hDeltaTPrev.emplace(i, make_shared<TH2F>(Form("DeltaTPrev_det%d", i), Form("%s: #Delta T for detector %d;channel;#Delta T, #mus", file.Data(), i), detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 10000, 0, 100000));
        for(auto j = detMin.at(i); j <= detMax.at(i); j++){
          hDeltaTPrevPerCharge.emplace(make_pair(i, j),
                                       make_shared<TH2F>(Form("DeltaTPrevPerCharge_det%d_ch%d", i, j), Form("%s: #Delta T for detector %d, channel %d;#Delta T, #mus; charge", file.Data(), i, j), 5000, 0, 50000, 512, 0, 1024));
        }
        d->cd();
      }
    }    
    out->cd();
  }
  map<pair<int, int>, tigerHitTL> prevHit;
  
  // auto straw_rt_dir = out->mkdir("straw_rt");
  // straw_rt_dir->cd();
  // map<int, shared_ptr<TH2F>> straw_rt, straw_rt_0;
  // for(auto i = detMin.at(1); i <= detMax.at(1); i++){
  //   straw_rt.emplace(i,
  //                    make_shared<TH2F>(Form("straw%d_rt", i),
  //                             Form("%s: straw %d v-shape sci ch 60;R, mm;T, ns", file.Data(), i),
  //                             32, -4, 4, 300, -100, 200));
  //   straw_rt_0.emplace(i,
  //                      make_shared<TH2F>(Form("straw%d_rt_0", i),
  //                               Form("%s: straw %d v-shape sci ch 0;R, mm;T, ns", file.Data(), i),
  //                               32, -4, 4, 300, -100, 200));
  // }
  // out->cd();

  out->mkdir("straw_deltat_corr")->cd();
  map<int, shared_ptr<TH1F>> straw_deltat, straw_deltat_0;
  for(auto i = detMin.at(1); i <= detMax.at(1); i++){
    straw_deltat.emplace(i,
                         make_shared<TH1F>(Form("straw%d_vs_sci60", i),
                                           Form("%s: straw%d_vs_sci60;#Delta t", file.Data(), i), 1000, -500, 500));
    straw_deltat_0.emplace(i,
                           make_shared<TH1F>(Form("straw%d_vs_sci0", i),
                                             Form("%s: straw%d_vs_sci0;#Delta t", file.Data(), i), 1000, -500, 500));
  }
  out->cd();

  out->mkdir("straw_banana")->cd();
  map<int, shared_ptr<TH2F>> straw_banana, straw_banana_0 ;
  for(auto i = detMin.at(1); i < detMax.at(1); i++){
    straw_banana.emplace(i,
                         make_shared<TH2F>(Form("straw%d-%d_banana", i, i+1),
                                           Form("%s: Time difference between straws %d, %d and sci60;T_{straw%d} - T_{scint}, [ns];T_{straw%d} - T_{scint}, [ns]", file.Data(), i, i+1, i, i+1),
                                           500, -250, 250, 500, -250, 250));
    straw_banana_0.emplace(i,
                           make_shared<TH2F>(Form("straw%d-%d_banana_0", i, i+1),
                                             Form("%s: Time difference between straws %d, %d and sci0;T_{straw%d} - T_{scint}, [ns];T_{straw%d} - T_{scint}, [ns]", file.Data(), i, i+1, i, i+1),
                                             500, -250, 250, 500, -250, 250));
  }
  out->cd();

  out->mkdir("straw_vs_straw_deltat")->cd();
  map<int, shared_ptr<TH1F>> straw_straw;
  for(auto i = detMin.at(1); i < detMax.at(1); i++){
    straw_straw.emplace(i,
                        make_shared<TH1F>(Form("straw%d_vs_straw%d", i, i+1),
                                          Form("%s: straw%d_vs_straw%d;T_{straw%d} - T_{straw%d}", file.Data(), i, i+1, i, i+1), 1000, -500, 500));
  }
  out->cd();

  shared_ptr<TH2F> DeltaTBetweenPulsers = make_shared<TH2F>("DeltaTBetweenPulsers", Form("%s: DeltaTBetweenPulsers;time, s;T^{scint}_{50#mus} - T^{scint}_{10ms}, #mus", file.Data()), 600, 0, 60, 10000, 0, 100000);
  map<int, shared_ptr<TH2F>> heFineCorr, hNeighborsPerTime;
  map<pair<int,int>, shared_ptr<TH2F>> heFinePerTimeCorr;
  if(detMax.at(6) > 0){
    out->mkdir("plots_corr_6")->cd();
    for(auto i = 2; i <= 5; i++){
      hNeighborsPerTime.emplace(i, make_shared<TH2F>(Form("hNeighborsPerTime_det%d", i), Form("%s: N neighbors per time for detector %d;full time, s; N neighbors", file.Data(), i),
                                                     600, 0, 60, detMax.at(i) - detMin.at(i) + 1, 0, detMax.at(i) - detMin.at(i) + 1));
      heFineCorr.emplace(i, make_shared<TH2F>(Form("eFine_det%d_corr", i), Form("%s: eFine for detector %d corellated with ship straw;channel;eFine", file.Data(), i),
                                              detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 512, 0, 1024));
      for(auto j = detMin.at(i); j <= detMax.at(i); j++){
        heFinePerTimeCorr.emplace(make_pair(i, j),
                                  make_shared<TH2F>(Form("eFinePerTime_corr_det%d_ch%d", i, j), Form("%s: eFine for detector %d, channel %d%s correllated with ship straw;full time, s; eFine", file.Data(), i, j,
                                                                                                     (detectorNames.count({i,j}) ? (string() + " (" + detectorNames.at({i,j}) + ")").c_str() : "")),
                                                    60, 0, 60, 512, 0, 1024));
      }
    }
    out->cd();
  }
  
  map<int, shared_ptr<TH1F>> hSciTimeToDet, hSciTimeToDetCoarse;
  map<int, shared_ptr<TH2F>> hSciTimeToDetCoarsePerTime;
  for(int i = 1; i <= nDetectorTypes; i++){
    hSciTimeToDet.emplace(i, make_shared<TH1F>(Form("sci_vs_det%d", i), Form("%s: T_{scint} - T_{det %d};#Deltat, ns", file.Data(), i), 1000, -500, 500));
    hSciTimeToDetCoarse.emplace(i, make_shared<TH1F>(Form("sci_vs_det%d_coarse", i), Form("%s: T_{scint} - T_{det %d} (coarse time);#DeltaT, ns", file.Data(), i), 160, -500, 500));
    hSciTimeToDetCoarsePerTime.emplace(i, make_shared<TH2F>(Form("sci_vs_det%d_coarse_per_time", i), Form("%s: T_{scint} - T_{det %d} (coarse time);time, s; #DeltaT, ns", file.Data(), i), 600, 0, 60, 160, -500, 500));
  }
  map<int, shared_ptr<TH2F>> hShipRT;
  if(detMax.at(6) > 0){
    for(int i = 2; i <= 4; i++){
      hShipRT.emplace(i, make_shared<TH2F>(Form("ship_rt_vs_mm%d", i-2), Form("%s: RT for SHiP straw and MM %d;strip;#DeltaT, ns", file.Data(), i-2),
                                           detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i), 2000, -1000, 1000));
    }
  }

  map<pair<int,int>, shared_ptr<TH2F>> mmCorrellations;
  out->mkdir("mm_corellations")->cd();
  for(auto i = 2; i <= 5; i++){
    for(auto j = 2; j <= 5; j++){
      if(i == j)
        continue;
      mmCorrellations[{i, j}] = make_shared<TH2F>(Form("mmCorrellations_det%d-%d", i, j),
                                                  Form("%s: N corellations between detectors %d and %d;strip (det %d); strip (det %d)", file.Data(), i, j, i, j),
                                                  detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1,
                                                  detMax.at(j) - detMin.at(j) + 1, detMin.at(j), detMax.at(j) + 1);
    }
  }
  out->cd();

  map<int, shared_ptr<TH2F>> hStrawRT, hStrawRTCoarse;
  if(detMax.at(1) > 0){
    for(int i = 2; i <= 5; i++){
      if(i == mmLayerY) continue;
      hStrawRT.emplace(i, make_shared<TH2F>(Form("straw_rt_vs_mm%d", i-2), Form("%s: RT for straws and MM %d;strip;#DeltaT, ns", file.Data(), i-2),
                                            detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i), 2000, -1000, 1000));
      hStrawRTCoarse.emplace(i, make_shared<TH2F>(Form("straw_rt_vs_mm%d_coarse", i-2), Form("%s: RT for straws and MM %d (coarse time);strip;#DeltaT (coarse), ns", file.Data(), i-2),
                                                  detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i), 320, -1000, 1000));
    }
  }
  
  
  Long64_t nentries = fChain->GetEntries();
  if(n > 0 && nentries > n)
    nentries = n;
  
  // =============================== CORRELATION FINDING ===============================
  Long64_t timeWindowNS = maxTimeDiff(); // ns
  Long64_t firstHitInWindow = 0;
  tigerHitTL *hitMain, *hitSecondary, hitFirst;
  for (Long64_t jentry = 0; jentry < nentries; jentry++)
  {
    if (!(jentry % 100000)){
      std::cout << "Entry " << jentry << "\t of \t" << nentries << "\n";
      // if(jentry>0) break;
      freeHitMap(firstHitInWindow);
    }
    hitMain = getHitFromTree(jentry, true);
    if (!jentry){
      printf("First hit: ");
      hitMain->print();
      hitFirst = *hitMain;
    }
    auto fchMapped = getMapped(hitMain);
    auto [fchD, fchM] = fchMapped;
    auto charge = hitMain->charge(energyMode);
    auto timeSinceStart = timeDifferenceFineNS(hitMain, &hitFirst)  *  1E-9;
    { // per-tiger histograms
      hTigerProfile.at({hitMain->gemrocID, hitMain->tigerID})->Fill(hitMain->channelID);
      if(fchD >=0)
        hTigerProfileMapped.at({hitMain->gemrocID, hitMain->tigerID})->Fill(hitMain->channelID);
      hTigerChargeToT.at({hitMain->gemrocID, hitMain->tigerID})->Fill(hitMain->channelID, hitMain->chargeToT());
      hTigerChargeSH.at({hitMain->gemrocID, hitMain->tigerID})->Fill(hitMain->channelID, hitMain->chargeSH());
      hTigerTimeFine.at({hitMain->gemrocID, hitMain->tigerID})->Fill(hitMain->channelID, hitMain->timeFine());
      hTigertCoarse.at({hitMain->gemrocID, hitMain->tigerID})->Fill(hitMain->channelID, hitMain->tCoarse);
      hTigereCoarse.at({hitMain->gemrocID, hitMain->tigerID})->Fill(hitMain->channelID, hitMain->eCoarse);
      hTigertFine.at({hitMain->gemrocID, hitMain->tigerID})->Fill(hitMain->channelID, hitMain->tFine);
      hTigereFine.at({hitMain->gemrocID, hitMain->tigerID})->Fill(hitMain->channelID, hitMain->eFine);
      hTigerFullTime.at({hitMain->gemrocID, hitMain->tigerID})->Fill(hitMain->channelID, timeSinceStart);
      hTigerFullTime1D.at({hitMain->gemrocID, hitMain->tigerID})->Fill(timeSinceStart);
      hTigerChargePerTime.at({hitMain->gemrocID, hitMain->tigerID, hitMain->channelID})->Fill(timeSinceStart, hitMain->charge(energyMode));
      hTigereFinePerTime.at({hitMain->gemrocID, hitMain->tigerID, hitMain->channelID})->Fill(timeSinceStart, hitMain->eFine);        
      hTigerFullTimePerChannel.at({hitMain->gemrocID, hitMain->tigerID, hitMain->channelID})->Fill(timeSinceStart);
    }
    if (fchD < 0) continue; // unmapped channels
    if(!isGoodHit(jentry))
      continue;
    if (fchD < nDetectorTypes){ // per-detector histograms
      hprofile.at(fchD)->Fill(fchM);
      hChargeToT.at(fchD)->Fill(fchM, hitMain->chargeToT());
      hChargeSH.at(fchD)->Fill(fchM, hitMain->chargeSH());
      hTimeFine.at(fchD)->Fill(fchM, hitMain->timeFine());
      htCoarse.at(fchD)->Fill(fchM, hitMain->tCoarse);
      heCoarse.at(fchD)->Fill(fchM, hitMain->eCoarse);
      htFine.at(fchD)->Fill(fchM, hitMain->tFine);
      heFine.at(fchD)->Fill(fchM, hitMain->eFine);
      htCoarse10bit.at(fchD)->Fill(fchM, hitMain->tCoarse%0x400);
      hFullTime.at(fchD)->Fill(fchM, timeSinceStart);
      hFullTime1D.at(fchD)->Fill(timeSinceStart);
      hFullTimePerChannel.at(fchMapped)->Fill(timeSinceStart);
      hChargePerTime.at(fchMapped)->Fill(timeSinceStart, hitMain->charge(energyMode));
      heFinePerTime.at(fchMapped)->Fill(timeSinceStart, hitMain->eFine);        
      if(!prevHit.count(fchMapped)){
        if(fchD == 0 || fchD == 7){
          printf("First hit in detector %d, channel %d (time since start: %lld us) : ", fchD, fchM, static_cast<long long>(timeDifferenceFineNS(hitMain, &hitFirst) * 1E-3));
          hitMain->print();
        }
        prevHit.emplace(fchMapped, *hitMain);
      }
      else{
        // if(fchD == 0)
        //   printf("D %d, ch %d, Diff to previous: %g\n", fchD, fchM, timeDifferenceFineNS(*hitMain, prevHit.at(fchMapped)));
        if(hDeltaTPrev.count(fchD))
          hDeltaTPrev.at(fchD)->Fill(fchM, timeDifferenceFineNS(hitMain, &prevHit.at(fchMapped)) * 1E-3);
        if(hDeltaTPrevPerCharge.count(fchMapped))
          hDeltaTPrevPerCharge.at(fchMapped)->Fill(timeDifferenceFineNS(hitMain, &prevHit.at(fchMapped)) * 1E-3, hitMain->charge(energyMode));
        prevHit[fchMapped] = *hitMain;
      }
    }

    map<int, map<int, tigerHitTL*>> closestHitsInLayer;

    if(fchD == 0 && fchM == 0){ // Scintillator
      closestHitsInLayer.clear();
      for(Long64_t kentry = firstHitInWindow; kentry < nentries; kentry++){
        if(kentry == jentry) continue;
        hitSecondary = getHitFromTree(kentry);
        if(!hitSecondary) continue;
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
        if(!closestHitsInLayer.count(ffchD)) closestHitsInLayer[ffchD] = {};

        // Search closest hits for all other channels
        if(!closestHitsInLayer.at(ffchD).count(ffchM) || timeDifferenceAbs < fabs(timeDifferenceFineNS(hitMain, closestHitsInLayer.at(ffchD).at(ffchM))))
          (closestHitsInLayer.at(ffchD))[ffchM] = hitSecondary;

      }
      // RT for Ship straw
      if(closestHitsInLayer.count(6) && closestHitsInLayer.at(6).count(0) &&
         fabs(timeDifferenceFineNS(hitMain, closestHitsInLayer.at(6).at(0))) < maxTimeDiff(fchD, 6)){
        for(auto i = 0; i <= 4; i++){
          auto idet = i+2;
          if(idet == mmLayerY) continue;
          if(!closestHitsInLayer.count(idet))
            continue;
          for(auto &h: closestHitsInLayer.at(idet)){
            if(fabs(timeDifferenceFineNS(hitMain, h.second)) > maxTimeDiff(fchD, idet))
              continue;
            hShipRT.at(idet)->Fill(h.first, timeDifferenceFineNS(hitMain, closestHitsInLayer.at(6).at(0)));
          }
        }
      }
      // RT for normal straws
      if(closestHitsInLayer.count(1)){
        for(auto &straw: closestHitsInLayer.at(1)){
          if(fabs(timeDifferenceFineNS(hitMain, straw.second)) > maxTimeDiff(fchD, 1))
            continue;
          for(auto i = 0; i <= 4; i++){
            auto idet = i+2;
            if(idet == mmLayerY) continue;
            if(!closestHitsInLayer.count(idet)) continue;
            for(auto &h: closestHitsInLayer.at(idet)){
              if(fabs(timeDifferenceFineNS(hitMain, h.second)) > maxTimeDiff(fchD, idet))
                continue;
              hStrawRT.at(idet)->Fill(h.first, timeDifferenceFineNS(hitMain, straw.second));
              hStrawRTCoarse.at(idet)->Fill(h.first, timeDifferenceCoarsePS(hitMain, straw.second) / 1E3);
            }
          }
        }
      }
      for(auto i = 1; i < 7; i++){
        if(!closestHitsInLayer.count(i))
          continue;
        if(hSciTimeToDet.count(i)){
          for(auto &h: closestHitsInLayer.at(i))
            hSciTimeToDet.at(i)->Fill(timeDifferenceFineNS(hitMain, h.second));
        }
        if(hSciTimeToDet.count(i)){
          for(auto &h: closestHitsInLayer.at(i)){
            hSciTimeToDetCoarse.at(i)->Fill(timeDifferenceCoarsePS(hitMain, h.second)/1E3);
            hSciTimeToDetCoarsePerTime.at(i)->Fill(timeSinceStart, timeDifferenceCoarsePS(hitMain, h.second)/1E3);
          }
        }
      }
    }
    else if(fchD == 0 && fchM == 4){ // 50us clocks
      if(prevHit.count({0, 5}))
        DeltaTBetweenPulsers->Fill(timeSinceStart, timeDifferenceFineNS(hitMain, &prevHit.at({0, 5})) * 1E-3);      
    }
    else if (fchD == 1){ // All straw ch
    }
    else if (fchD >= 2 && fchD <= 5){ // All MM ch
      /* Searching for secondary hit */
      closestHitsInLayer.clear();
      for(Long64_t kentry = firstHitInWindow; kentry < nentries; kentry++){
        if(kentry == jentry) continue;
        hitSecondary = getHitFromTree(kentry);
        if(!hitSecondary) continue;
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
        if(!closestHitsInLayer.count(ffchD)) closestHitsInLayer[ffchD] = {};
        if (ffchD == 1){ // Straw
        } else if(ffchD == fchD && ffchM == fchM){ // same channel
        } else if(ffchD == fchD){ // same MM
          if(!closestHitsInLayer.at(ffchD).count(ffchM) || timeDifferenceAbs < fabs(timeDifferenceFineNS(hitMain, closestHitsInLayer.at(ffchD).at(ffchM))))
            if(timeDifferenceAbs < maxTimeDiff(fchD, ffchD))
              (closestHitsInLayer.at(ffchD))[ffchM] = hitSecondary;
        } else if(ffchD >= 2 && ffchD <= 5){ // other MM
          if(!closestHitsInLayer.at(ffchD).count(ffchM) || timeDifferenceAbs < fabs(timeDifferenceFineNS(hitMain, closestHitsInLayer.at(ffchD).at(ffchM))))
            if(timeDifferenceAbs < maxTimeDiff(fchD, ffchD))
              (closestHitsInLayer.at(ffchD))[ffchM] = hitSecondary;
        } else if(ffchD == 0 && ffchM == 0){ // Scint. coinc.
          if(!closestHitsInLayer.at(ffchD).count(ffchM) || timeDifferenceAbs < fabs(timeDifferenceFineNS(hitMain, closestHitsInLayer.at(ffchD).at(ffchM))))
            if(timeDifferenceAbs < maxTimeDiff(fchD, ffchD))
              (closestHitsInLayer.at(ffchD))[ffchM] = hitSecondary;
        } else if(ffchD == 0 && ffchM == 3){
        } else if(ffchD == 0 && ffchM == 4){
        } else if(ffchD == 0 && ffchM == 5){
        } else if(ffchD == 6 && ffchM == 0){ // SHiP straw
          if(!closestHitsInLayer.at(ffchD).count(ffchM) || timeDifferenceAbs < fabs(timeDifferenceFineNS(hitMain, closestHitsInLayer.at(ffchD).at(ffchM))))
            if(timeDifferenceAbs < maxTimeDiff(fchD, ffchD))
              (closestHitsInLayer.at(ffchD))[ffchM] = hitSecondary;
        }
      }
      /* Filling histograms */
      if(closestHitsInLayer.count(0) && closestHitsInLayer.at(0).count(0)){
        mm_vs_sci.at(fchD)->Fill(timeDifferenceFineNS(hitMain, closestHitsInLayer.at(0).at(0)));
        mm_vs_sciCoarse.at(fchD)->Fill(timeDifferenceCoarsePS(hitMain, closestHitsInLayer.at(0).at(0))/1E3);
      }
      if(closestHitsInLayer.count(6) && closestHitsInLayer.at(6).count(0)){
        heFineCorr.at(fchD)->Fill(fchM, hitMain->eFine);
        heFinePerTimeCorr.at(fchMapped)->Fill(timeSinceStart, hitMain->eFine);
      }
      // calculating number of hits in the same 
      if(closestHitsInLayer.count(fchM) && closestHitsInLayer.at(fchM).size()){
        int firstInClaster, lastInClaster;
        for(auto j = fchM+1; j <= detMax.at(fchD); j++){
          if(!closestHitsInLayer.at(fchM).count(j) && !closestHitsInLayer.at(fchM).count(j+1)){
            lastInClaster = j-1;
            break;
          }
        }
        for(auto j = fchM-1; j >= detMin.at(fchD); j--){
          if(!closestHitsInLayer.at(fchM).count(j) && !closestHitsInLayer.at(fchM).count(j-1)){
            firstInClaster = j+1;
            break;
          }
        }
        hNeighborsPerTime.at(fchD)->Fill(timeSinceStart, lastInClaster - firstInClaster);
      }

      for(auto i = 2; i <= 5; i++){
        if(fchD == i) continue;
        if(!mmCorrellations.count({fchD, i})) continue;
        if(!closestHitsInLayer.count(i)) continue;
        for(auto &v: closestHitsInLayer.at(i))
          mmCorrellations.at({fchD, i})->Fill(fchM, v.first);
      }
    }
    // else if (fchD == 6){ // All straw ch
    // }
  }

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
