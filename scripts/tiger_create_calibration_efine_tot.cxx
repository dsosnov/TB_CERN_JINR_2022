#include "tiger_fullTime.cxx"

#include <memory>
using std::shared_ptr, std::make_shared;

void tiger_create_calibration_efine_tot(TChain *chain){
  chain->SetBranchStatus("gemrocID",1);
  chain->SetBranchStatus("tigerID",1);
  chain->SetBranchStatus("channelID",1);
  chain->SetBranchStatus("eFine",1);

  auto gemroc = 0;
  chain->Draw("tigerID * 64 + channelID:eFine >> h_efine(1024, 0, 1024, 512, 0, 512)", Form("gemrocID == %d", gemroc));
  auto h = static_cast<TH2F*>(gDirectory->Get("h_efine"));  

  auto fout = fopen("../out/tiger_efine_calibration_ToT.txt", "w");
  for(auto t = 0; t < 8; t++){
    printf("Tiger: %d\n", t);
    for(auto ch = 0; ch < 64; ch++){
      auto bin = (t * 64 + ch) + 1;
      auto p = h->ProjectionX("", bin, bin);
      auto max = p->GetMaximum();
      auto firstBin = p->FindFirstBinAbove(max/2.0);
      auto lastBin = p->FindLastBinAbove(max/2.0);
      auto minTFine = static_cast<int>(p->GetXaxis()->GetBinLowEdge(firstBin));
      auto maxTFine = static_cast<int>(p->GetXaxis()->GetBinUpEdge(lastBin));
      if(minTFine >= 0 && maxTFine >= 0){
        fprintf(fout, "%d %d %d %d %d\n", gemroc, t, ch, minTFine, maxTFine);
        printf("%d %d %d %d %d\n", gemroc, t, ch, minTFine, maxTFine);
      }
    }
  }
  fclose(fout);
  // h->SaveAs(Form("h_%d%d%d%d.root", gemroc, t, ch, tac));
}



void tiger_create_calibration_efine_tot(string path, bool directory = false){
  auto chain = new TChain("tigerTL");
  if(directory)
    chain->Add((path + "/*.root").c_str());
  else
    chain->Add(path.c_str());

  tiger_create_calibration_efine_tot(chain);
}
