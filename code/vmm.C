#define vmm_cxx
#include "vmm.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TGraph.h>

void vmm::Loop()
{
   auto *spills = new TH1D("spills", "Num of spills; t, sec", 300, 0., 900.);
   auto *bcid_sci = new TH1D("bcid_sci", "bcid_sci; BCID", 4200, 0, 4200);
   auto *tdo_sci = new TH1D("tdo_sci", "tdo_sci; TDO", 500, 0, 500);

   auto *sci_vs_straw = new TH1D("sci_vs_straw", "sci_vs_straw; #Delta t, ns", 500, -500, 500);

   auto *straw_check = new TH1D("straw_check", "straw_check; TDO - F(PDO), ns", 500, -500, 500);
   auto *sci_check = new TH1D("sci_check", "sci_check; TDO - F(PDO), ns", 500, -500, 500);

   auto *p0H = new TGraph();
   auto *p1H = new TGraph();
   auto *hwhm = new TH2D("ch check", "check; ch; #Delta t, ns", 64, 0, 64, 500, -500, 500);


   TH2D *h[64];
   char name[10];
   char title[20];
   for (Int_t i = 0; i < 64; i++) 
   {
      sprintf(name, "h_ch%d", i);
      sprintf(title, "TDO of ch%d", i);
      h[i] = new TH2D(name, title, 64, 0, 1024, 128, 0, 256);
   }

   int sum = 0;
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

      sum += channel->at(0).size();
      for (int j = 0; j < channel->at(0).size(); j++)
      {
         h[channel->at(0).at(j)]->Fill(pdo->at(0).at(j), tdo->at(0).at(j));
         if (channel->at(0).at(j) == 0)
         {
            bcid_sci->Fill(bcid->at(0).at(j));
            tdo_sci->Fill(tdo->at(0).at(j));
         }
         else
         {
            continue;
         }
      }
   }


   std::cout << "nom of hits " << sum << std::endl;

   TFile *out = new TFile("../out/out_" + file + ending, "RECREATE");
   spills->Write("spills");
   bcid_sci->Write("bcid_sci");
   tdo_sci->Write("tdc_sci");
   TDirectory *cdth2 = out->mkdir("TDO vs PDO");
   cdth2->cd();
   for (Int_t i = 0; i < 64; i++) 
   {
      h[i]->Write();
   }

   auto *fitfunc = new TF1("fitfunc", "[p0] + [p1] * pow(x, -1.25)", 100, 1024);
   auto *fitfunchw = new TF1("fitfunchw", "[p0]", 100, 700);

   TDirectory *cdprof = out->mkdir("profiles");
   cdprof->cd();
   TDirectory *cdhwgraph = out->mkdir("shifted TDO");

   vector<array<double, 2>> fit_params;

   for (Int_t i = 0; i < 64; i++) 
   {
      fitfunc->SetParameter(0, 0);
      fitfunc->SetParameter(1, 0);
      char profilename[10];
      sprintf(profilename, "profile_ch%d", i);
      auto *currentBin = new TH1D();
      auto *profile = h[i]->ProfileX(profilename, 1, h[i]->GetNbinsX(),"");
      profile->Fit(fitfunc);
      p0H->AddPoint(i, fitfunc->GetParameter(0));
      p1H->AddPoint(i, fitfunc->GetParameter(1));
      fit_params.push_back({fitfunc->GetParameter(0), fitfunc->GetParameter(1)});
      cdprof->cd();
      profile->Write(profilename);
   }


      
   cdhwgraph->cd();
   TH1D *hwhmhist[64];
   for (Int_t i = 0; i < 64; i++) 
   {
      char hwhmname[20];
      sprintf(hwhmname, "shifted_TDO_ch%d", i);
      hwhmhist[i] = new TH1D(hwhmname, hwhmname, 200, -100, 100);
   }


   nbytes = 0;
   for (Long64_t jentry = 0; jentry < nentries; jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0)
         break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;

      for (int j = 0; j < channel->at(0).size(); j++)
      {
         int fch = channel->at(0).at(j);
         int fpdo = pdo->at(0).at(j);
         int ftdo = tdo->at(0).at(j);
         fitfunc->SetParameter(0, fit_params[fch][0]);
         fitfunc->SetParameter(1, fit_params[fch][1]);
         hwhmhist[fch]->Fill(ftdo - fitfunc->Eval(fpdo));
      }
   }

   vector <array<double, 3>> full_params;
   for (Int_t i = 0; i < 64; i++) 
   {
      double bin1 = hwhmhist[i]->GetBinCenter(hwhmhist[i]->FindFirstBinAbove(hwhmhist[i]->GetMaximum()/5));
      double bin2 = hwhmhist[i]->GetBinCenter(hwhmhist[i]->FindLastBinAbove(hwhmhist[i]->GetMaximum()/5));
      double fhwhm = (bin2 - bin1) / 2;

      full_params.push_back({fit_params[i][0], fit_params[i][1], fhwhm});

      hwhmhist[i]->Write();
   }

   cdprof->cd("../");

   nbytes = 0;
   double sciT = 0;
   for (Long64_t jentry = 0; jentry < nentries; jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0)
         break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      
      for (int j = 0; j < channel->at(0).size(); j++)
      {
         int trsh = threshold->at(0).at(j);
         if (!trsh)
         continue;
         int fch = channel->at(0).at(j);
         if (fch == 35 || fch == 63)
         continue;
         int fpdo = pdo->at(0).at(j);
         if (fpdo < 100 || fpdo > 800)
         continue;
         int ftdo = tdo->at(0).at(j);
         int fbcid = bcid->at(0).at(j);
         double time = fbcid * 25.0 - (ftdo - (fitfunc->Eval(fpdo) - full_params[fch][2])) * 25.0 / (full_params[fch][2] * 2.0);
         fitfunc->SetParameter(0, full_params[fch][0]);
         fitfunc->SetParameter(1, full_params[fch][1]);
         hwhmhist[fch]->Fill(ftdo - fitfunc->Eval(fpdo));

         if (fch == 30)
         {
            straw_check->Fill((ftdo - (fitfunc->Eval(fpdo) - full_params[fch][2])) * 25.0 / (full_params[fch][2] * 2.0));
         }

         if (fch == 0) 
         {
            sci_check->Fill((ftdo - 96) * 25.0 / (160.0 - 96.0));
            sciT = fbcid * 25.0 - (ftdo - 96) * 25.0 / (160.0 - 96.0);
            // for (int k = 0; k < channel->at(0).size(); k++)
            // {
            //    if (j == k)
            //    continue;
               
            //    int ffch = channel->at(0).at(k);
            //    int ffpdo = pdo->at(0).at(k);
            //    int fftdo = tdo->at(0).at(k);
            //    int ffbcid = bcid->at(0).at(k);
            //    int ftime = ffbcid * 25 - (fftdo - (fitfunc->Eval(ffpdo) - full_params[ffch][2])) * 25 / (full_params[ffch][2] * 2);

            //    if (ffch == 30) 
            //    {
            //       sci_vs_straw->Fill(sciT - ftime);
            //    }
            // }
         }
         if (fch == 30) 
         {
            sci_vs_straw->Fill(sciT - time);
         }

         if (fch != 0) 
         {
            hwhm->Fill(fch, sciT - time);
         }

      }
   }



   p0H->Write("p0H");
   p1H->Write("p1H");
   hwhm->Write("hwhm");
   sci_vs_straw->Write("sci_vs_straw");
   straw_check->Write("straw_check");
   sci_check->Write("sci_check");
   out->Close();
}
