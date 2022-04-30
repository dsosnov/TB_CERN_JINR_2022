#define calibration_cxx
#include "calibration.h"
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>

void calibration::Loop()
{
   // TDO distribution for every Ch
   TH1D *h[64];
   char name[10];
   char title[20];
   for (Int_t i = 0; i < 64; i++)
   {
      sprintf(name, "h_ch%d", i);
      sprintf(title, "TDO of ch%d", i);
      h[i] = new TH1D(name, title, 128, 0, 256);
   }

   // =============================== TDO & distributions ===============================

   if (fChain == 0)
      return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry = 0; jentry < nentries; jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0)
         break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      for (int j = 0; j < channel->at(0).size(); j++)
      {
         int chtemp = channel->at(0).at(j);
         if (chtemp == 0) 
            std::cout << "!!!" << std::endl;
         h[chtemp]->Fill(tdo->at(0).at(j));
      }
   }
   // ===================================================================================

   vector<array<int, 2>> limits; // vector of TDO limits for every Ch [limits.size() == 64]!

   TFile *out = new TFile("../out/calibration_25_100" + ending, "RECREATE"); // PATH where to save out_*.root file
   TDirectory *tdo_dir = out->mkdir("TDO");
   tdo_dir->cd();

   // ================================== LIMITS SEARCH ==================================
   ofstream out_txt;
   out_txt.open("../out/calibration_25_100.txt");

   auto *gausFitF = new TF1("gausFitF", "gaus", 0, 256); //
   for (Int_t i = 0; i < 64; i++)
   {
      // gausFitF->SetParameter(0, 0);
      // gausFitF->SetParameter(1, 0);
      // gausFitF->SetParameter(2, 0);
      // h[i]->RebinX(4);
      // h[i]->Fit(gausFitF, "Q");                                                   // fit the TDO distribution with Gauss
      // int bin1 = gausFitF->GetParameter(1) - 2 * gausFitF->GetParameter(2); // 2 sigma to left edge
      // int bin2 = gausFitF->GetParameter(1) + 2 * gausFitF->GetParameter(2); // 2 sigma to right edge

      int bin1 = h[i]->GetBinCenter(h[i]->FindFirstBinAbove(5));
      int bin2 = h[i]->GetBinCenter(h[i]->FindLastBinAbove(5));

      limits.push_back({bin1, bin2});
      h[i]->Write();

      std::cout << "CH " << i << "\t L TDO: " << bin1 << "\t R TDO: " << bin2 << "\n";

      out_txt << bin1 << " " << bin2 << "\n";
   }
   // ===================================================================================

   tdo_dir->cd("../");
   out->Close();
}

int main(){
  auto c = new calibration();
  c->Loop();
  return 0;
}