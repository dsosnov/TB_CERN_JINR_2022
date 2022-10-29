#include "tiger_fullTime.cxx"

#include <memory>
using std::shared_ptr, std::make_shared;

void tiger_create_calibration_tfine(TChain *chain){
  chain->SetBranchStatus("gemrocID",1);
  chain->SetBranchStatus("tigerID",1);
  chain->SetBranchStatus("channelID",1);
  chain->SetBranchStatus("tFine",1);
  chain->SetBranchStatus("tacID",1);

  auto gemroc = 0;
  chain->Draw("(tigerID * 64 + channelID) * 4 + tacID:tFine >> h_tfine(1024, 0, 1024, 2048, 0, 2048)", Form("gemrocID == %d", gemroc));
  auto h = static_cast<TH2F*>(gDirectory->Get("h_tfine"));
  

  auto fout = fopen("../out/tiger_tfine_calibration.txt", "w");
  for(auto t = 0; t < 8; t++){
    printf("Tiger: %d\n", t);        
    for(auto ch = 0; ch < 64; ch++){
      for(auto tac = 0; tac < 4; tac++){
        auto bin = (t * 64 + ch) * 4 + tac + 1;
        auto p = h->ProjectionX("", bin, bin);
        auto max = p->GetMaximum();
        auto firstBin = p->FindFirstBinAbove(max/2.0);
        auto lastBin = p->FindLastBinAbove(max/2.0);
        auto minTFine = static_cast<int>(p->GetXaxis()->GetBinLowEdge(firstBin));
        auto maxTFine = static_cast<int>(p->GetXaxis()->GetBinUpEdge(lastBin));
        if(minTFine >= 0 && maxTFine >= 0){
          fprintf(fout, "%d %d %d %d %d %d\n", gemroc, t, ch, tac, minTFine, maxTFine);        
          printf("%d %d %d %d %d %d\n", gemroc, t, ch, tac, minTFine, maxTFine);
        }
      }
    }
  }
  fclose(fout);
  // h->SaveAs(Form("h_%d%d%d%d.root", gemroc, t, ch, tac));
}


// for eFine:
// 1. Define time between spills
// tigerTL->Draw("tiger_fullTimeCoarse_auto(tCoarse,frameCount,frameCountLoops) >> h(60, 0, 60)", "tiger_fullTimeCoarse_auto(tCoarse,frameCount,frameCountLoops) < 60")
// void tiger_create_calibration_efine(TChain *chain){
//   chain->SetBranchStatus("gemrocID",1);
//   chain->SetBranchStatus("tigerID",1);
//   chain->SetBranchStatus("channelID",1);
//   chain->SetBranchStatus("eFine",1);
//   chain->SetBranchStatus("tCoarse",1);
//   chain->SetBranchStatus("frameCount",1);
//   chain->SetBranchStatus("frameCountLoops",1);

//   auto gemroc = 0;
//   // that can be done in interactive session
//   chain->Draw("tigerID * 64 + channelID:eFine >> h_efine(1024, 0, 1024, 512, 0, 512)", Form("gemrocID == %d && tiger_fullTimeCoarse_auto(tCoarse,frameCount,frameCountLoops) > 8 && tiger_fullTimeCoarse_auto(tCoarse,frameCount,frameCountLoops) < 40", gemroc));
//   auto h = static_cast<TH2F*>(gDirectory->Get("h_efine"));
  // auto fout = fopen("../out/tiger_efine_calibration.txt", "w");
  // for(auto t = 0; t < 8; t++){
  //   printf("Tiger: %d\n", t);        
  //   for(auto ch = 0; ch < 64; ch++){
  //     auto bin = t * 64 + ch + 1;
  //     auto p = h->ProjectionX("", bin, bin);
  //     if(p->GetEntries() < 1000) continue;
  //     auto max = p->GetMaximum();
  //     auto firstBin = p->FindFirstBinAbove(max/2.0);
  //     auto lastBin = p->FindLastBinAbove(max/2.0);
  //     auto minTFine = static_cast<int>(p->GetXaxis()->GetBinLowEdge(firstBin));
  //     auto maxTFine = static_cast<int>(p->GetXaxis()->GetBinUpEdge(lastBin));
  //     if(minTFine >= 0 && maxTFine >= 0){
  //       fprintf(fout, "%d %d %d %d %d\n", gemroc, t, ch, maxTFine, minTFine);        
  //       printf("%d %d %d %d %d\n", gemroc, t, ch, minTFine, maxTFine);
  //     }
  //   }
  // }
  // fclose(fout);
//   // h->SaveAs(Form("h_%d%d%d%d.root", gemroc, t, ch, tac));
// }

  
void tiger_create_calibration(string path, bool directory = false, bool createTFine = true, bool createEFine = true){  
  auto chain = new TChain("tigerTL");
  if(directory)
    chain->Add((path + "/*.root").c_str());
  else
    chain->Add(path.c_str());
  
  chain->SetBranchStatus("*",0);
  chain->SetBranchStatus("gemrocID",1);
  chain->SetBranchStatus("tigerID",1);
  chain->SetBranchStatus("channelID",1);
  chain->SetBranchStatus("tFine",1);
  chain->SetBranchStatus("eFine",1);
  chain->SetBranchStatus("tacID",1);
  Char_t   gemrocID;        //                "B" == Char_t   ==  int8_t
  Short_t  tigerID;         //  8 bit data -- "S" == Short_t  == int16_t
  Char_t   channelID;       //  6 bit data -- "B" == Char_t   ==  int8_t
  Short_t  tFine;           // 10 bit data -- "S" == Short_t  == int16_t
  Short_t  eFine;           // 10 bit data -- "S" == Short_t  == int16_t
  Char_t   tacID;           //  2 bit data -- "B" == Char_t   ==  int8_t
  chain->SetBranchAddress("gemrocID", &gemrocID);
  chain->SetBranchAddress("tigerID", &tigerID);
  chain->SetBranchAddress("channelID", &channelID);
  chain->SetBranchAddress("tacID", &tacID);
  chain->SetBranchAddress("tFine", &tFine);
  chain->SetBranchAddress("eFine", &eFine);

  if(createTFine && !createEFine){
    tiger_create_calibration_tfine(chain);
    return;
  // } else if(!createTFine && createEFine){
  //   tiger_create_calibration_efine(chain);
  //   return;
  }


  map<tuple<Char_t, Short_t, Char_t>, map<Char_t, shared_ptr<TH1D>>> tFines;
  map<tuple<Char_t, Short_t, Char_t>, Short_t> eFineMax, eFineMin;
  auto nEntries = chain->GetEntries();
  for(auto i = 0; i < nEntries; i++){
    if(!(i%10000000))
      printf("Entry %d of %lld\n", i, nEntries);
    chain->GetEntry(i);

    if(!tFines.count({gemrocID, tigerID, channelID})){
      for(auto j = 0; j < 4; j++)
        tFines[{gemrocID, tigerID, channelID}].emplace(j, make_shared<TH1D>(Form("tfine_%d%d%d%d",gemrocID, tigerID, channelID, j), "", 1024, 0, 1024));
    }
    tFines.at({gemrocID, tigerID, channelID}).at(tacID)->Fill(tFine);
    if(!eFineMax.count({gemrocID, tigerID, channelID}) || eFineMax.at({gemrocID, tigerID, channelID}) < eFine)
      eFineMax[{gemrocID, tigerID, channelID}] = eFine;
    if(!eFineMin.count({gemrocID, tigerID, channelID}) || eFineMin.at({gemrocID, tigerID, channelID}) > eFine)
      eFineMin[{gemrocID, tigerID, channelID}] = eFine;
  }

  if(createTFine){
    auto fout = fopen("../out/tiger_tfine_calibration.txt", "w");
    fprintf(fout, "# Generated from %s\n", path.c_str());
    for(auto &tfm: tFines){
      for(auto &tf: tfm.second){
        auto max = tf.second->GetMaximum();
        auto firstBin = tf.second->FindFirstBinAbove(max/2.0);
        auto lastBin = tf.second->FindLastBinAbove(max/2.0);
        auto minTFine = static_cast<int>(tf.second->GetXaxis()->GetBinLowEdge(firstBin));
        auto maxTFine = static_cast<int>(tf.second->GetXaxis()->GetBinUpEdge(lastBin));
        fprintf(fout, "%d %d %d %d %d %d\n", get<0>(tfm.first), get<1>(tfm.first), get<2>(tfm.first), tf.first, minTFine, maxTFine);
      }
    }
    fclose(fout);
  }
  if(createEFine){
    auto fout = fopen("../out/tiger_efine_calibration.txt", "w");
    fprintf(fout, "# Generated from %s\n", path.c_str());
    for(auto &p: eFineMin){
      fprintf(fout, "%d %d %d %d %d\n", get<0>(p.first), get<1>(p.first), get<2>(p.first), p.second, eFineMax.at(p.first));
    }
    fclose(fout);
  }
}
