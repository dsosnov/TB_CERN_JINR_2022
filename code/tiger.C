#define tiger_cxx
#include "tiger.h"

#include "TH1F.h"
#include "TH2F.h"

void tiger::Loop(unsigned long n)
{
  printf("tiger::Loop()\n");

  if (fChain == 0)
    return;

  vector<int> detMin(nDetectorTypes, -1), detMax(nDetectorTypes, -1);
  for(auto &s: channelMap){
    if(s.second.first < 0 || s.second.first >= nDetectorTypes)
      continue;
    printf("For channel %d -- %d\n", s.second.first, s.second.second);
    if(detMin.at(s.second.first) < 0 || detMin.at(s.second.first) > s.second.second)
      detMin.at(s.second.first) = s.second.second;
    if(detMax.at(s.second.first) < 0 || detMax.at(s.second.first) < s.second.second)
      detMax.at(s.second.first) = s.second.second;
  }
  printf("Straws (det 1): %d-%d\n", detMin.at(1), detMax.at(1));
  for(int i = 0; i < 4; i++)
    printf("MM %d (det %d): %d-%d\n", i, i+2, detMin.at(i+2), detMax.at(i+2));
  printf("Additional straws (det 6): %d-%d\n", detMin.at(6), detMax.at(6));
  printf("Lemo (det 7): %d-%d\n", detMin.at(7), detMax.at(7));

  TFile *out = new TFile("../out/out_tiger" + file + ending, "RECREATE"); // PATH where to save out_*.root file

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
  vector<TH2F*> hpdo;
  for(auto i = 0; i < nDetectorTypes; i++){
    out->mkdir(Form("det%d", i))->cd();
    hprofile.push_back(new TH1F(Form("profile_det%d", i), Form("%s: profile for detector %d;ch", file.Data(), i), detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1));
    hpdo.push_back(new TH2F(Form("pdo_det%d", i), Form("%s: pdo for detector %d;ch;pdo", file.Data(), i), detMax.at(i) - detMin.at(i) + 1, detMin.at(i), detMax.at(i) + 1, 1025, 0, 1025));
    out->cd();
  }
  
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
  
  Long64_t nentries = fChain->GetEntries();
  if(n > 0 && nentries > n)
    nentries = n;
  
  // =============================== CORRELATION FINDING ===============================
  Long64_t timeWindowNS = 1E3; // ns
  Long64_t firstHitInWindow = 0;
  tigerHitTL hitMain, hitSecondary;
  for (Long64_t jentry = 0; jentry < nentries; jentry++) // You can remove "/ 10" and use the whole dataset
  {
    if (jentry % 10000 == 0)
    {
      std::cout << "Entry " << jentry << "\t of \t" << nentries << "\n";
    }
    fChain->GetEntry(jentry);
    updateTigerHitTLCurrent(hitMain);
    auto [fchD, fchM] = getMapped(hitMain);
    if (fchD >=0 && fchD < nDetectorTypes){
      hprofile.at(fchD)->Fill(fchM);
      hpdo.at(fchD)->Fill(fchM, hitMain.charge());
    }
    if (fchD == 1){ // All straw ch
      /* Searching for secondary hit */
      /* TODO improve speed: maybe load hits to vector and work with vector in memory */
      pair<Long64_t, tigerHitTL> closestNeighbor = {-1, tigerHitTL()};
      pair<Long64_t, tigerHitTL> closestSci0 = {-1, tigerHitTL()};
      pair<Long64_t, tigerHitTL> closestSci60 = {-1, tigerHitTL()};
      map<pair<int, int>, pair<Long64_t, tigerHitTL>> closestMM; // key: mapped detector-strip
      for(Long64_t kentry = firstHitInWindow; kentry < nentries; kentry++){
        fChain->GetEntry(kentry);
        updateTigerHitTLCurrent(hitSecondary);
        /* Checking that second hit in maximum time window */
        if(timeDifferenceFineNS(hitMain, hitSecondary) > timeWindowNS){
          firstHitInWindow++;
          continue;
        } else if(timeDifferenceFineNS(hitSecondary, hitMain) > timeWindowNS){
          break;
        }
        auto [ffchD, ffchM] = getMapped(hitSecondary);
        auto timeDifferenceAbs = fabs(timeDifferenceFineNS(hitMain, hitSecondary));
        if (ffchD == 1 && ffchM == fchM + 1){
          if(closestNeighbor.first < 0 ||
             timeDifferenceAbs < fabs(timeDifferenceFineNS(hitMain, closestNeighbor.second)))
            closestNeighbor = {kentry, hitSecondary};
        } else if(ffchD >= 2 && ffchD <= 5){
          bool save = false;
          if(!closestMM.count({ffchD, ffchM}))
            save = true;
          else if(timeDifferenceAbs < fabs(timeDifferenceFineNS(hitMain, closestMM.at({ffchD, ffchM}).second)))
            save = true;
          if(save)
            closestMM.emplace(make_pair(ffchD, ffchM), make_pair(kentry, hitSecondary));
        } else if(ffchD == 0 && ffchM == 0){
          if(closestSci0.first < 0 ||
             timeDifferenceAbs < fabs(timeDifferenceFineNS(hitMain, closestSci0.second)))
            closestSci0 = {kentry, hitSecondary};
        } else if(ffchD == 0 && ffchM == 3){
          if(closestSci60.first < 0 ||
             timeDifferenceAbs < fabs(timeDifferenceFineNS(hitMain, closestSci60.second)))
            closestSci60 = {kentry, hitSecondary};
        } else if(ffchD == 0 && ffchM == 4){
        }
      }
      /* Filling histograms */
      if(closestSci0.first >= 0){
        straw_vs_sci->Fill(timeDifferenceFineNS(hitMain, closestSci0.second));
        straw_deltat_0.at(fchM)->Fill(timeDifferenceFineNS(hitMain, closestSci0.second));
      }
      if(closestSci60.first >= 0){
        straw_deltat.at(fchM)->Fill(timeDifferenceFineNS(hitMain, closestSci60.second));
      }
      if(closestNeighbor.first >= 0){
        straw_straw.at(fchM)->Fill(timeDifferenceFineNS(hitMain, closestNeighbor.second));
        if(closestSci0.first >= 0)
          straw_banana_0.at(fchM)->Fill(timeDifferenceFineNS(hitMain,
                                                             closestSci0.second),
                                        timeDifferenceFineNS(closestNeighbor.second,
                                                             closestSci0.second));
        if(closestSci60.first >= 0)
          straw_banana.at(fchM)->Fill(timeDifferenceFineNS(hitMain,
                                                           closestSci60.second),
                                      timeDifferenceFineNS(closestNeighbor.second,
                                                           closestSci60.second));
        
      }
      for(auto &mm: closestMM){
        straw_vs_mm.at(mm.first.first-2)->Fill(timeDifferenceFineNS(hitMain, mm.second.second));
        straw_vs_mm_spatial_corr.at(mm.first.first-2)->Fill(fchM, mm.first.second);
      }
    }
    if (fchD >= 2 && fchD <= 5){ // All MM ch
      /* Searching for secondary hit */
      /* TODO improve speed: maybe load hits to vector and work with vector in memory */
      pair<Long64_t, tigerHitTL> closestSci0 = {-1, tigerHitTL()};
      for(Long64_t kentry = firstHitInWindow; kentry < nentries; kentry++){
        fChain->GetEntry(kentry);
        updateTigerHitTLCurrent(hitSecondary);
        /* Checking that second hit in maximum time window */
        if(timeDifferenceFineNS(hitMain, hitSecondary) > timeWindowNS){
          firstHitInWindow++;
          continue;
        } else if(timeDifferenceFineNS(hitSecondary, hitMain) > timeWindowNS){
          break;
        }
        auto [ffchD, ffchM] = getMapped(hitSecondary);
        auto timeDifferenceAbs = fabs(timeDifferenceFineNS(hitMain, hitSecondary));
        if (ffchD == 1 && ffchM == fchM + 1){
        } else if(ffchD >= 2 && ffchD <= 5){
        } else if(ffchD == 0 && ffchM == 0){
          if(closestSci0.first < 0 ||
             timeDifferenceAbs < fabs(timeDifferenceFineNS(hitMain, closestSci0.second)))
            closestSci0 = {kentry, hitSecondary};
        } else if(ffchD == 0 && ffchM == 3){
        } else if(ffchD == 0 && ffchM == 4){
        }
      }
      /* Filling histograms */
      if(closestSci0.first >= 0){
        mm_vs_sci.at(fchD)->Fill(timeDifferenceFineNS(hitMain, closestSci0.second));
      }
    }
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

  out->Write();
  out->Close();

}
