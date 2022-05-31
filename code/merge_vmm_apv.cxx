#include "apv.C"
#include "evBuilder.C"

vector<map<unsigned long, unsigned long>> tryToMerge(vector<pair<unsigned long, analysisGeneral::mm2CenterHitParameters>> &hits_vmm,
                                                     vector<pair<unsigned long, analysisGeneral::mm2CenterHitParameters>> &hits_apv,
                                                     unsigned long i_vmm, unsigned long i_apv,
                                                     long long nHitsAfterTrigVMM, bool ignoreCountersVMM, bool ignoreCountersAPV,
                                                     bool verbose = false){
  if(verbose)
    printf("tryToMerge({%lu}, {%lu}, %lld, %d, %d)\n", hits_vmm.size() - i_vmm, hits_apv.size() - i_apv, nHitsAfterTrigVMM, ignoreCountersVMM, ignoreCountersAPV);
  double limitsTime = 1E6 / 100; // MSec
  int limitsStrips = 5; // 10;
  double limitTrigCount = 0.1;  // relative
  unsigned int trigTimeDiff = 25;
  double pdoDifference = 0.25;  //

  // WARNING: temporary
  bool checkCounters = !(ignoreCountersVMM || ignoreCountersAPV);

  vector<map<unsigned long, unsigned long>> options = {};
  
  if(hits_vmm.empty() || hits_apv.empty())
    return options;
  if((i_vmm >= hits_vmm.size()) || (i_apv >= hits_apv.size()))
    return options;

  auto hit_vmm = hits_vmm.at(i_vmm);
  if(verbose){
    printf("  vmm (%lu): ", i_vmm); hit_vmm.second.print();
  }

  unsigned long long nHitsAfterTrigAPV = 0;

  for(auto i = i_apv; i < hits_apv.size(); i++){
    auto hit_apv = hits_apv.at(i);
    bool accepted = true;
    if(verbose){
      printf("  apv (%lu): ", i); hit_apv.second.print();
    }

    if(hit_apv.second.timeFull() - hit_vmm.second.timeFull() > limitsTime)
      break;
    else if(hit_vmm.second.timeFull() - hit_apv.second.timeFull() > limitsTime){
      if(verbose){
        printf("    ! time difference: %lld\n", hit_vmm.second.timeFull() - hit_apv.second.timeFull());
      }
      accepted = false;
    }
    else if(abs(static_cast<long long>(hit_apv.second.stripX) - static_cast<long long>(hit_vmm.second.stripX)) > limitsStrips){
      if(verbose){
        printf("    ! strip-strip: %lld\n", abs(static_cast<long long>(hit_apv.second.stripX) - static_cast<long long>(hit_vmm.second.stripX)));
      }
      accepted = false;
    }
    // else if(checkCounters && (abs(hit_apv.second.nHitsToPrev - (hit_vmm.second.nHitsToPrev + nHitsAfterTrigVMM)) > hit_vmm.second.nHitsToPrev * limitTrigCount)){
    //   if(verbose){
    //     printf("    ! Counters: %lld\n", abs(hit_apv.second.nHitsToPrev - (hit_vmm.second.nHitsToPrev + nHitsAfterTrigVMM)));
    //   }
    //   accepted = false;
    // }

    // else if(fabs(hit_apv.second.time - hit_vmm.second.time) > trigTimeDiff){
    //   if(verbose){
    //     printf("    ! Trig time diff: %f\n", fabs(hit_apv.second.time - hit_vmm.second.time));
    //   }
    //   accepted = false;
    // }
    
    // else if(!hit_apv.second.approximated && !hit_vmm.second.approximated && fabs(hit_apv.second.pdoRelative - hit_vmm.second.pdoRelative) > pdoDifference)
    //   accepted = false;
        
    if(!accepted){
      nHitsAfterTrigAPV += hit_apv.second.nHitsToPrev;
      continue;
    }

    // hit_apv.second.print();

    auto tryOption = tryToMerge(hits_vmm, hits_apv, i_vmm+1, i+1, 0, false, false);
    if(tryOption.empty()){
      tryOption.push_back({{hit_vmm.first, hit_apv.first}});
    } else {
      for(auto &opt: tryOption)
        opt.emplace(hit_vmm.first, hit_apv.first);
    }
    std::copy(tryOption.begin(), tryOption.end(), std::back_inserter(options));
  }

  /* No such hit in APV data */
  auto option_pass = tryToMerge(hits_vmm, hits_apv, i_vmm+1, i_apv, nHitsAfterTrigVMM + hit_vmm.second.nHitsToPrev, ignoreCountersVMM, ignoreCountersAPV);
  std::copy(option_pass.begin(), option_pass.end(), std::back_inserter(options));

  return options;
  
}

void merge_vmm_apv(){
  pair<string, string> run_pair = {"run_0258", "run166"};

  auto out = TFile::Open("merge_vmm_apv.root", "recreate");

  // unsigned long long from = 1653145627, to = 1653145628;
  // unsigned long long from = 1653145640, to = 1653145640;
  unsigned long long from = 0, to = 0;

  auto hTrigTimeAPV = make_shared<TH1F>("hTrigTimeAPV", "hTrigTimeAPV", 27, 0, 27*25);
  auto hTrigTimeVMM = make_shared<TH1F>("hTrigTimeVMM", "hTrigTimeVMM", 27, 0, 27*25);

  auto apvan = new apv(run_pair.second);
  auto hits_apv = apvan->GetCentralHits(from, to);
  auto vmman = new evBuilder(run_pair.first, "g1_p25_s100-0&60", "map-20220518");
  // vmman->Loop();
  auto hits_vmm = vmman->GetCentralHits(from, to);
  for(auto &hh: {&hits_vmm, &hits_apv}){
    double pdoMax = 0.0;
    for(auto &h: *hh)
      if((!h.second.approximated) && h.second.pdoRelative > pdoMax)
        pdoMax = h.second.pdoRelative;
    for(auto &h: *hh)
      h.second.pdoRelative /= pdoMax;          
  }
  for(auto &h: hits_vmm)
    h.second.time += 325;

  printf("APV hits (%lu)\n", hits_apv.size());
  for(auto &h: hits_apv){
    hTrigTimeAPV->Fill(h.second.time);
    h.second.print();
  }
  printf("VMM hits (%lu):\n", hits_vmm.size());
  for(auto &h: hits_vmm){
    hTrigTimeVMM->Fill(h.second.time);
    h.second.print();
  }

  hTrigTimeVMM->SaveAs("hTrigTimeVMM.root");
  hTrigTimeAPV->SaveAs("hTrigTimeAPV.root");

  return;

  vector<pair<unsigned long, analysisGeneral::mm2CenterHitParameters>> hits_vmm_v;
  hits_vmm_v.assign(hits_vmm.begin(), hits_vmm.end());
  vector<pair<unsigned long, analysisGeneral::mm2CenterHitParameters>> hits_apv_v;
  hits_apv_v.assign(hits_apv.begin(), hits_apv.end());

  auto options = tryToMerge(hits_vmm_v, hits_apv_v, 0, 0, 0, true, true);

  std::sort(options.begin(), options.end(), [](auto &a, auto &b){return a.size() > b.size();});

  printf("Possible options size: %lu\n", options.size());

  if(options.size() > 50)
    options.resize(50);

  auto dir_timediff = out->mkdir("daq_time_diff");
  auto dir_maxqtimediff = out->mkdir("maxq_time_diff");
  auto dir_stripdiff = out->mkdir("strip_diff");
  auto dir_pdodiff = out->mkdir("pdo_diff");
  auto dir_pdodiff_existed = out->mkdir("pdo_diff_existed");
  
  for(ulong i = 0; i < options.size(); i++){
    auto option = options.at(i);
    dir_timediff->cd();
    auto h_timediff = new TH1F(Form("h_timediff_%lu", i), Form("h_timediff_%lu (%lu hits)", i, option.size()), 2E5, -1E6, 1E6);
    dir_maxqtimediff->cd();
    auto h_maxqtimediff = new TH1F(Form("h_maxqtimediff_%lu", i), Form("h_maxqtimediff_%lu (%lu hits)", i, option.size()), 54, -675, 675);
    dir_stripdiff->cd();
    auto h_stripdiff = new TH1F(Form("h_stripdiff_%lu", i), Form("h_stripdiff_%lu (%lu hits)", i, option.size()), 40, -20, 20);
    dir_pdodiff->cd();
    auto h_pdodiff = new TH1F(Form("h_pdodiff_%lu", i), Form("h_pdodiff_%lu (%lu hits)", i, option.size()), 400, -2, 2);
    dir_pdodiff_existed->cd();
    auto h_pdodiff_existed = new TH1F(Form("h_pdodiff_existed_%lu", i), Form("h_pdodiff_existed_%lu (%lu hits)", i, option.size()), 400, -2, 2);
    printf("Option -- %lu pairs\n", option.size());

    for(auto &o: option){
      auto hit_vmm = hits_vmm.at(o.first);
      auto hit_apv = hits_apv.at(o.second);

      h_timediff->Fill(hit_vmm.timeFull() - hit_apv.timeFull());
      h_maxqtimediff->Fill(hit_vmm.time - hit_apv.time);
      h_stripdiff->Fill(hit_vmm.stripX - hit_apv.stripX);
      h_pdodiff->Fill(hit_vmm.pdoRelative - hit_apv.pdoRelative);
      if(!hit_vmm.approximated && !hit_apv.approximated)
        h_pdodiff_existed->Fill(hit_vmm.pdoRelative - hit_apv.pdoRelative);
    }    
  }

  out->Write();
  out->Close();

}
