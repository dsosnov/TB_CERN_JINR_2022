#include "apv.C"
#include "evBuilder.C"

vector<map<unsigned long, unsigned long>> tryToMerge(vector<pair<unsigned long, analysisGeneral::mm2CenterHitParameters>> &hits_vmm,
                                                     vector<pair<unsigned long, analysisGeneral::mm2CenterHitParameters>> &hits_apv,
                                                     unsigned long i_vmm, unsigned long i_apv,
                                                     long long nHitsAfterTrigVMM, bool ignoreCountersVMM, bool ignoreCountersAPV,
                                                     bool verbose = false){
  if(verbose)
    printf("tryToMerge({%lu}, {%lu}, %lld, %d, %d)\n", hits_vmm.size() - i_vmm, hits_apv.size() - i_apv, nHitsAfterTrigVMM, ignoreCountersVMM, ignoreCountersAPV);
  double limitsTime = 1E4; // MSec
  int limitsStrips = 5; // 10;
  // double limitTrigCount = 0.1;  // relative
  // unsigned int trigTimeDiff = 25;
  // double pdoDifference = 0.25;  //

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

    auto tryOption = tryToMerge(hits_vmm, hits_apv, i_vmm+1, i+1, 0, false, false, verbose);
    if(tryOption.empty()){
      tryOption.push_back({{hit_vmm.first, hit_apv.first}});
    } else {
      for(auto &opt: tryOption)
        opt.emplace(hit_vmm.first, hit_apv.first);
    }
    std::copy(tryOption.begin(), tryOption.end(), std::back_inserter(options));
  }

  /* No such hit in APV data */
  auto option_pass = tryToMerge(hits_vmm, hits_apv, i_vmm+1, i_apv, nHitsAfterTrigVMM + hit_vmm.second.nHitsToPrev, ignoreCountersVMM, ignoreCountersAPV, verbose);
  std::copy(option_pass.begin(), option_pass.end(), std::back_inserter(options));

  return options;
  
}

// void printfBrief(pair<unsigned long, analysisGeneral::mm2CenterHitParameters> hit, bool revert = false){
//   if(revert)
//     printf("%c %3d - %.2f - %.3f - %7lld - %2llu | %5lu",
//            (hit.second.approximated ? '*' : ' ' ),
//            hit.second.stripX, hit.second.time, hit.second.pdoRelative,
//            hit.second.timeFull() % int(1E7), hit.second.nHitsToPrev, hit.first);
//   else
//     printf("%5lu | %2llu - %7lld - %.3f - %.2f - %3d %c",
//            hit.first, hit.second.nHitsToPrev, hit.second.timeFull() % int(1E7),
//            hit.second.pdoRelative, hit.second.time, hit.second.stripX,
//            (hit.second.approximated ? '*' : ' ' ));
// }

void merge_vmm_apv(){
  bool REMOVE_APPROX = false;
  // pair<string, string> run_pair = {"run_0258", "run166"};
  pair<string, string> run_pair = {"run_0814", "run407"};

  auto out = TFile::Open("merge_vmm_apv.root", "recreate");

  // unsigned long long from = 1653145627, to = 1653145628;
  // unsigned long long from = 1653145640, to = 1653145640;
  unsigned long long from = 0, to = 0;

  auto hTrigTimeAPV = make_shared<TH1F>("hTrigTimeAPV", "hTrigTimeAPV", 27, 0, 27*25);
  auto hTrigTimeVMM = make_shared<TH1F>("hTrigTimeVMM", "hTrigTimeVMM", 27, 0, 27*25);

  auto apvan = new apv(run_pair.second);
  apvan->useSyncSignal();
  auto hits_apv = apvan->GetCentralHits(from, to);
  // auto vmman = new evBuilder(run_pair.first, "g1_p25_s100-0&60", "map-20220518");
  auto vmman = new evBuilder(run_pair.first, "g1_p25_s100-0&60", "map-20220605");
  vmman->useSyncSignal();
  auto hits_vmm = vmman->GetCentralHits(from, to);
  for(auto &h: hits_vmm)
    h.second.time += 325;

  vector<pair<unsigned long, analysisGeneral::mm2CenterHitParameters>> hits_vmm_v;
  hits_vmm_v.assign(hits_vmm.begin(), hits_vmm.end());
  vector<pair<unsigned long, analysisGeneral::mm2CenterHitParameters>> hits_apv_v;
  hits_apv_v.assign(hits_apv.begin(), hits_apv.end());

  /* Removing approximate hits */
  if(REMOVE_APPROX)
    hits_apv_v.erase(std::remove_if(hits_apv_v.begin(), 
                                    hits_apv_v.end(),
                                    [](auto c){return c.second.approx;}),
                     hits_apv_v.end());

  printf("APV hits (%lu) -- VMM hits (%lu)\n", hits_apv.size(), hits_vmm.size());
  for(unsigned long i = 0; i < hits_apv_v.size() || i < hits_vmm_v.size(); i++){
    if(i < hits_apv_v.size()){
      // printfBrief(hits_apv_v.at(i), false);
      printf("%5lu | ", hits_vmm_v.at(i).first);
      hits_apv_v.at(i).second.printfBrief();
    } else
      printf("%5c | %2c - %7c - %6c - %6c - %3c %c", ' ', ' ', ' ', ' ', ' ', ' ', ' ');
    printf(" | %*c %4lu %*c | ", 10, ' ', i, 10, ' ');
    if(i < hits_vmm_v.size()){
      // printfBrief(hits_vmm_v.at(i), true);
      hits_vmm_v.at(i).second.printfBrief(true);
      printf(" | %5lu", hits_vmm_v.at(i).first);
    } else
      printf("%c %3c - %6c - %6c - %7c - %2c | %5c", ' ', ' ', ' ', ' ', ' ', ' ', ' ');
    printf("\n");
  }
  printf("\n\n");

  for(auto &h: hits_apv){
    hTrigTimeAPV->Fill(h.second.time);
    // h.second.print();
  }
  for(auto &h: hits_vmm){
    hTrigTimeVMM->Fill(h.second.time);
    // h.second.print();
  }

  hTrigTimeVMM->SaveAs("hTrigTimeVMM.root");
  hTrigTimeAPV->SaveAs("hTrigTimeAPV.root");

  return; // TODO delete

  vector<map<unsigned long, unsigned long>> options;
    options = tryToMerge(hits_vmm_v, hits_apv_v, 0, 0, 0, true, true, false);
  // if(hits_vmm_v.size() < hits_apv_v.size())
  //   options = tryToMerge(hits_vmm_v, hits_apv_v, 0, 0, 0, true, true, false);
  // else{
  //   auto options_rev = tryToMerge(hits_apv_v, hits_vmm_v, 0, 0, 0, true, true, false);
  //   for(auto &o: options_rev){
  //     map<unsigned long, unsigned long> option;
  //     for(auto &p: o)
  //       option.emplace(p.second, p.first);
  //     options.push_back(option);
  //   }
  // }
  
  std::sort(options.begin(), options.end(), [&hits_apv, &hits_vmm](auto &a, auto &b){
    if(a.size() != b.size())
      return a.size() > b.size();
    else{
      double meanA = 0.0, meanB = 0.0;
      for(auto &o: a)
        meanA += (hits_vmm.at(o.first).stripX - hits_apv.at(o.second).stripX);
      for(auto &o: b)
        meanB += (hits_vmm.at(o.first).stripX - hits_apv.at(o.second).stripX);
      meanA /= a.size();
      meanB /= b.size();
      return meanA < meanB;
    }
    });

  printf("Possible options size: %lu\n", options.size());
  auto maximum_size = options.begin()->size();

  if(options.size() > 150)
    options.resize(150);

  auto dir_timediff = out->mkdir("daq_time_diff");
  auto dir_maxqtimediff = out->mkdir("maxq_time_diff");
  auto dir_stripdiff = out->mkdir("strip_diff");
  auto dir_pdodiff = out->mkdir("pdo_diff");
  auto dir_pdodiff_existed = out->mkdir("pdo_diff_existed");
  
  for(unsigned long i = 0; i < options.size(); i++){
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

    auto h_pdo2d = new TH2F(Form("h_pdo2d_%lu", i), Form("h_pdo2d_%lu (%lu hits); vmm pdo; apv pdo", i, option.size()), 400, -2, 2, 400, -2, 2);
    dir_pdodiff_existed->cd();
    for(auto &o: option){
      auto hit_vmm = hits_vmm.at(o.first);
      auto hit_apv = hits_apv.at(o.second);

      h_timediff->Fill(hit_vmm.timeFull() - hit_apv.timeFull());
      h_maxqtimediff->Fill(hit_vmm.time - hit_apv.time);
      h_stripdiff->Fill(hit_vmm.stripX - hit_apv.stripX);
      h_pdodiff->Fill(hit_vmm.pdoRelative - hit_apv.pdoRelative);
      h_pdo2d->Fill(hit_vmm.pdoRelative, hit_apv.pdoRelative);
      if(!hit_vmm.approx && !hit_apv.approx)
        h_pdodiff_existed->Fill(hit_vmm.pdoRelative - hit_apv.pdoRelative);
    }    
  }

  options.erase(std::remove_if(options.begin(), 
                               options.end(),
                               [maximum_size](auto c){return c.size() != maximum_size;}),
                options.end());

  for(unsigned long i = 0; i < options.size(); i++){
    printf("Sution %lu:\n", i);
    int j = 0;
    for(auto &o: options.at(i)){
      auto event_vmm = o.first;
      auto event_apv = o.second;
      printf("%5lu | ", event_apv);
      hits_apv.at(event_apv).printfBrief();
      // printfBrief({event_apv, hits_apv.at(event_apv)}, false);
      printf(" | %*c %4d %*c | ", 10, ' ', j++, 10, ' ');
      hits_vmm.at(event_vmm).printfBrief();
      printf(" | %5lu", event_vmm);
      // printfBrief({event_vmm, hits_vmm.at(event_vmm)}, true);
      printf("\n");
    }
      printf("\n");    
  }


  out->Write();
  out->Close();

}
