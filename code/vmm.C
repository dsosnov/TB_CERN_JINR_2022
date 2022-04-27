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

   auto *sci_vs_straw_all = new TH1D("sci_vs_straw_all", "sci_vs_straw ALL; #Delta t, ns", 500, -1000, 1000);
   auto *sci_vs_straw = new TH1D("sci_vs_straw", "sci_vs_straw; #Delta t, ns", 500, -500, 500);
   auto *straw_vs_straw_all = new TH1D("straw_vs_straw_all", "straw ch 31 vs straw ch 30 ALL; #Delta t, ns", 500, -1000, 1000);
   auto *straw_vs_straw = new TH1D("straw_vs_straw", "straw ch 31 vs straw ch 30; #Delta t, ns", 500, -500, 500);
   
   auto *straw_vs_sci = new TH1D("straw_vs_sci", "sci_vs_straw; #Delta t, ns", 500, -500, 500);


   auto *straw_check = new TH1D("straw_check", "straw_check; TDO - F(PDO), ns", 500, -500, 500);
   auto *sci_check = new TH1D("sci_check", "sci_check; TDO - F(PDO), ns", 500, -500, 500);

   auto *hwhm = new TH2D("ch check", "check; ch; #Delta t, ns", 64, 0, 64, 500, -500, 500);

   auto *threedetcorr = new TH2D("3detcorr", "#Delta t (straw_{ch30} - sci) vs #Delta t (straw_{ch31} - sci) with 1 #mu s time window between ch30 & ch31; #Delta t (straw_{ch31} - sci), ns; #Delta t (straw_{ch30} - sci), ns", 500, -500, 500, 500, -500, 500);

   auto *scistud = new TH2D("scistud", "sci TDO converted to ns vs #Delta t (straw_{ch30} - sci); #Delta t (straw_{ch30} - sci), ns; sci TDO in ns", 200, -100, 100, 500, -500, 500);

   auto *deltabcid = new TH1D("deltabcid", "straw BCID_{ch30} - BCID_{ch31}; #Delta BCID", 200, -100, 100);
   auto *deltabcid_2d = new TH2D("deltabcid_2d", "straw BCID_{ch30} vs BCID_{ch31};BCID;BCID", 1000, 0, 4096, 1000, 0, 4096);
   auto *deltabcid_sci = new TH1D("deltabcid_sci", "straw BCID_{ch30} - BCID_{sci}; #Delta BCID", 200, -100, 100);

   TH1D *h[64];
   char name[10];
   char title[20];
   for (Int_t i = 0; i < 64; i++)
   {
      sprintf(name, "h_ch%d", i);
      sprintf(title, "TDO of ch%d", i);
      h[i] = new TH1D(name, title, 128, 0, 256);
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
         h[channel->at(0).at(j)]->Fill(tdo->at(0).at(j));
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
   vector<array<int, 2>> limits;

   TFile *out = new TFile("../out/out_" + file + ending, "RECREATE");
   spills->Write("spills");
   bcid_sci->Write("bcid_sci");
   tdo_sci->Write("tdc_sci");
   TDirectory *cdth2 = out->mkdir("TDO");
   cdth2->cd();

   auto *g1 = new TF1("m1","gaus",0,256);
   for (Int_t i = 0; i < 64; i++)
   {
      g1->SetParameter(0, 0);
      g1->SetParameter(1, 0);
      g1->SetParameter(2, 0);
      h[i]->RebinX(4);
      h[i]->Fit(g1, "Q");
      int bin1 = g1->GetParameter(1) - 2 * g1->GetParameter(2);
      int bin2 = g1->GetParameter(1) + 2 * g1->GetParameter(2);

      // int bin1 = h[i]->GetBinCenter(h[i]->FindFirstBinAbove(h[i]->GetMaximum() / 10));
      // int bin2 = h[i]->GetBinCenter(h[i]->FindLastBinAbove(h[i]->GetMaximum() / 10));
      limits.push_back({bin1, bin2});
      h[i]->Write();

      std::cout << "CH " << i << "\t L TDO: " << bin1 << "\t R TDO: " << bin2 << "\n";
   }

   cdth2->cd("../");

   nbytes = 0;
   nb = 0;
   double sciT = 0, scidaq = 0;
   double t30 = 0, straw30daq = 0;
   double t31 = 0, straw31daq = 0;
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

         if (fch == 0)
         {
            // t30 = fbcid * 25.0 - (ftdo - limits[fch][0]) * 25.0 / (limits[fch][1] - limits[fch][0]);
            t30 = fbcid * 25.0 - (ftdo - 80) * 25.0 / (152 - 80);
            straw30daq = daq_timestamp_s->at(0) * 1e9 + daq_timestamp_ns->at(0);
            Long64_t mbytes = 0, mb = 0;
            double minT = 1e6;
            double minTsci = 1e6;
            for (Long64_t kentry = jentry - 20; kentry < jentry + 20; kentry++)
            {
               Long64_t iientry = LoadTree(kentry);
               if (iientry < 0)
                  continue;
               mb = fChain->GetEntry(kentry);
               mbytes += mb;

               for (int k = 0; k < channel->at(0).size(); k++)
               {
                  double ffdaq = daq_timestamp_s->at(0) * 1e9 + daq_timestamp_ns->at(0);
                  if (abs(ffdaq - straw30daq) > 1e6)
                     continue;

                  int ffpdo = pdo->at(0).at(k);
                  if (ffpdo < 100 || ffpdo > 800)
                     continue;

                  int ffch = channel->at(0).at(k);
                  if (ffch != 1)
                     continue;

                  int fftdo = tdo->at(0).at(k);
                  int ffbcid = bcid->at(0).at(k);
                  // double fft = ffbcid * 25.0 - (fftdo - limits[ffch][0]) * 25.0 / (limits[ffch][1] - limits[ffch][0]);
                  double fft = ffbcid * 25.0 - (fftdo - 104) * 25.0 / (168 - 104);
                  straw_check->Fill((fftdo - limits[ffch][0]) * 25.0 / (limits[ffch][1] - limits[ffch][0]));
                  straw_vs_straw_all->Fill(t30 - fft);
                  deltabcid->Fill(ffbcid - fbcid);
                  if (abs(t30 - fft) < minT)
                  {
                     minT = abs(t30 - fft);
                     straw31daq = ffdaq;
                     t31 = fft;
                  }
               }
            }

            mbytes = 0, mb = 0;
            for (Long64_t kentry = jentry - 20; kentry < jentry + 20; kentry++)
            {
               Long64_t iientry = LoadTree(kentry);
               if (iientry < 0)
                  continue;
               mb = fChain->GetEntry(kentry);
               mbytes += mb;

               for (int k = 0; k < channel->at(0).size(); k++)
               {
                  double ffdaq = daq_timestamp_s->at(0) * 1e9 + daq_timestamp_ns->at(0);
                  if (abs(ffdaq - straw30daq) > 1e6)
                     continue;

                  int ffpdo = pdo->at(0).at(k);
                  if (ffpdo < 100 || ffpdo > 800)
                     continue;

                  int ffch = channel->at(0).at(k);
                  if (ffch != 2)
                     continue;

                  int fftdo = tdo->at(0).at(k);
                  int ffbcid = bcid->at(0).at(k);
                  double fft = ffbcid * 25.0 - (fftdo - 96) * 25 / (160 - 96);
                  
                 
                  scistud->Fill(t30 - fft, (fftdo - 96) * 25.0 / (160.0 - 96.0));
                  sci_vs_straw_all->Fill(t30 - fft);
                  deltabcid_2d->Fill(ffbcid, fbcid);
                  deltabcid_sci->Fill(ffbcid - fbcid);

                  if (abs(t30 - fft) < minTsci)
                  {
                     minTsci = abs(t30 - fft);
                     scidaq = ffdaq;
                     sciT = fft;
                  }
               }
            }

            if (sciT != 0 && t31 != 0 && abs(t31 - sciT) < 1e6 && abs(t31 - t30) < 200)
            {
               sci_vs_straw->Fill(t31 - sciT);
               straw_vs_straw->Fill(t30 - t31);
               straw_vs_sci->Fill(t30 - sciT);
               threedetcorr->Fill(t31 - sciT, t30 - sciT);
            }
         }
         else
         {
            if (fch == 0)
            {
                sci_check->Fill((ftdo - 96) * 25.0 / (160.0 - 96.0));
            }
            continue;
         }
      }
   }

   hwhm->Write("hwhm");
   straw_check->Write("straw_check");
   sci_check->Write("sci_check");
   sci_vs_straw_all->Write("sci_vs_straw_all");
   sci_vs_straw->Write("sci_vs_straw");
   straw_vs_straw_all->Write("straw_vs_straw_all");
   straw_vs_straw->Write("straw_vs_straw");
   straw_vs_sci->Write("straw_vs_sci");
   threedetcorr->Write("threedetcorr");
   scistud->Write("scistud");
   deltabcid->Write("deltabcid");
   deltabcid_sci->Write("deltabcid_sci");
   deltabcid_2d->Write("deltabcid_2d");
   out->Close();
}
