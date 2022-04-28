#define vmm_cxx
#include "vmm.h"
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>

void vmm::Loop()
{
   // fast check plots
   auto *tdo_sci0 = new TH1D("tdo_sci0", "tdo_sci0; TDO", 128, 0, 256);
   auto *tdo_sci1 = new TH1D("tdo_sci1", "tdo_sci1; TDO", 128, 0, 256);
   auto *tdo_sci2 = new TH1D("tdo_sci2", "tdo_sci2; TDO", 128, 0, 256);
   auto *tdo_straw31 = new TH1D("tdo_straw31", "tdo_straw31; TDO", 128, 0, 256); // without 350 PDO counts cut
   auto *tdo_vs_pdo_straw31 = new TH2D("tdo_vs_pdo_straw31", "tdo_vs_pdo_straw31; TDO", 256, 0, 1024, 128, 0, 256);

   // correlation plots
   auto *straw31_vs_sci0 = new TH1D("straw31_vs_sci0", "straw31_vs_sci0; #Delta t, ns", 500, -500, 500);
   auto *straw31_vs_sci0_all = new TH1D("straw31_vs_sci0_all", "straw31_vs_sci0_all; #Delta t, ns", 500, -500, 500);
   auto *straw31_vs_sci1 = new TH1D("straw31_vs_sci1", "straw31_vs_sci1; #Delta t, ns", 500, -500, 500);
   auto *straw31_vs_sci2 = new TH1D("straw31_vs_sci2", "straw31_vs_sci2; #Delta t, ns", 500, -500, 500);

   auto *sci0_vs_sci1 = new TH1D("sci0_vs_sci1", "sci0_vs_sci1; #Delta t, ns", 500, -500, 500);
   auto *sci0_vs_sci2 = new TH1D("sci0_vs_sci2", "sci0_vs_sci2; #Delta t, ns", 500, -500, 500);
   auto *sci1_vs_sci2 = new TH1D("sci1_vs_sci2", "sci1_vs_sci2; #Delta t, ns", 500, -500, 500);

   auto *straw31_vs_straw30 = new TH1D("straw31_vs_straw30", "straw31_vs_straw30; #Delta t, ns", 500, -500, 500);
   auto *straw31_vs_straw30_all = new TH1D("straw31_vs_straw30_all", "straw31_vs_straw30_all; #Delta t, ns", 500, -500, 500);

   auto *straw31_vs_straw30_banana = new TH2D("straw31_vs_straw30_banana",
                                              "straw31_vs_straw30_banana; straw 31 #Deltat, ns; straw 30 #Deltat, ns", 500, -250, 250, 500, -250, 250);
   auto *straw31_vs_straw30_banana_all = new TH2D("straw31_vs_straw30_banana_all",
                                              "straw31_vs_straw30_banana_all; straw 31 #Deltat, ns; straw 30 #Deltat, ns", 500, -250, 250, 500, -250, 250);

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
         if (pdo->at(0).at(j) > 350) // !!! <------ setup for TDO distributions for every Ch
         {
            h[chtemp]->Fill(tdo->at(0).at(j));
         }
         if (chtemp == 0)
         {
            tdo_sci0->Fill(tdo->at(0).at(j));
         }
         else if (chtemp == 1)
         {
            tdo_sci1->Fill(tdo->at(0).at(j));
         }
         else if (chtemp == 2)
         {
            tdo_sci2->Fill(tdo->at(0).at(j));
         }
         else if (chtemp == 31)
         {
            tdo_straw31->Fill(tdo->at(0).at(j));
            tdo_vs_pdo_straw31->Fill(pdo->at(0).at(j), tdo->at(0).at(j));
         }
         else
         {
            continue;
         }
      }
   }
   // ===================================================================================

   vector<array<int, 2>> limits; // vector of TDO limits for every Ch [limits.size() == 64]!

   TFile *out = new TFile("../out/out_" + file + ending, "RECREATE"); // PATH where to save out_*.root file
   TDirectory *tdo_dir = out->mkdir("TDO");
   tdo_dir->cd();

   // ================================== LIMITS SEARCH ================================== //or get from file
   // auto *gausFitF = new TF1("gausFitF", "gaus", 0, 256); //
   ifstream myfile ("../out/calibration_25_100.txt");
   int bin1 = 0;
   int bin2 = 0;
   int i = 0;
   while (myfile >> bin1 >> bin2)
   {
      limits.push_back({bin1, bin2});
      std::cout << "CH " << i << "\t L TDO: " << bin1 << "\t R TDO: " << bin2 << "\n";
      i++;
   }
   for (Int_t i = 0; i < 64; i++)
   {
      // gausFitF->SetParameter(0, 0);
      // gausFitF->SetParameter(1, 0);
      // gausFitF->SetParameter(2, 0);
      // h[i]->RebinX(4);
      // h[i]->Fit(gausFitF, "Q");                                                   // fit the TDO distribution with Gauss
      // int bin1 = gausFitF->GetParameter(1) - 2 * gausFitF->GetParameter(2); // 2 sigma to left edge
      // int bin2 = gausFitF->GetParameter(1) + 2 * gausFitF->GetParameter(2); // 2 sigma to right edge

      // int bin1 = h[i]->GetBinCenter(h[i]->FindFirstBinAbove(h[i]->GetMaximum() / 10));
      // int bin2 = h[i]->GetBinCenter(h[i]->FindLastBinAbove(h[i]->GetMaximum() / 10));

      // limits.push_back({bin1, bin2});
      h[i]->Write();

      // std::cout << "CH " << i << "\t L TDO: " << bin1 << "\t R TDO: " << bin2 << "\n";
   }
   // ===================================================================================
   tdo_dir->cd("../");
   
   nbytes = 0;
   nb = 0;

   // =============================== CORRELATION FINDING ===============================
   for (Long64_t jentry = 0; jentry < nentries; jentry++) // You can remove "/ 10" and use the whole dataset
   {
      if (jentry % 10000 == 0)
      {
         std::cout << "Entry " << jentry << "\t of \t" << nentries << "\n";
      }
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0)
         break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;

      double t31 = 0;

      for (int j = 0; j < channel->at(0).size(); j++)
      {
         int fch = channel->at(0).at(j);
         if (fch == 35 || fch == 63)
            continue; // remove 'bad' ch for future tasks

         if (fch == 31)
         {
            int fpdo = pdo->at(0).at(j);
            int ftdo = tdo->at(0).at(j);
            int fbcid = bcid->at(0).at(j);

            t31 = fbcid * 25.0 - (ftdo - limits[fch][0]) * 25.0 / (limits[fch][1] - limits[fch][0]); // 'auto' limits
            // t31 = fbcid * 25.0 - (fftdo - X) * 25.0 / (Y - X); //'hand' limits

            Long64_t mbytes = 0, mb = 0;
            double minTsci0 = 1e3;
            double minTsci1 = 1e3;
            double minTsci2 = 1e3;
            double minT30 = 1e3;

            double sciT_ch0 = 0, sciT_ch1 = 0, sciT_ch2 = 0;
            double t30 = 0;

            // ========================         LOOP OVER 40 events around         ========================
            //                           jentry to find correlation with straw 30

            for (Long64_t kentry = jentry - 20; kentry < jentry + 20; kentry++)
            {
               Long64_t iientry = LoadTree(kentry);
               if (iientry < 0)
                  continue;
               mb = fChain->GetEntry(kentry);
               mbytes += mb;

               for (int k = 0; k < channel->at(0).size(); k++)
               {
                  int ffch = channel->at(0).at(k);
                  if (ffch != 30)
                     continue;

                  int ffpdo = pdo->at(0).at(k);
                  int fftdo = tdo->at(0).at(k);
                  int ffbcid = bcid->at(0).at(k);
                  double fft = ffbcid * 25.0 - (fftdo - limits[ffch][0]) * 25.0 / (limits[ffch][1] - limits[ffch][0]); // 'auto' limits
                  // double fft = ffbcid * 25.0 - (fftdo - X) * 25.0 / (Y - X); //'hand' limits

                  straw31_vs_straw30_all->Fill(t31 - fft);
                  if (abs(t31 - fft) < minT30)
                  {
                     minT30 = abs(t31 - fft);
                     t30 = fft;
                  }
               }
            }
            if (t30 != 0)
            {
               straw31_vs_straw30->Fill(t30 - t31);
            }
            // ============================ end of straw 30 correlation finding ===========================

            // ========================         LOOP OVER 40 events around         ========================
            //                           jentry to find correlation with sci 0

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
                  int ffch = channel->at(0).at(k);
                  if (ffch != 0)
                     continue;

                  int ffpdo = pdo->at(0).at(k);
                  int fftdo = tdo->at(0).at(k);
                  int ffbcid = bcid->at(0).at(k);
                  double fft = ffbcid * 25.0 - (fftdo - 96) * 25.0 / (156 - 96); //'hand' limits
                  // double fft = ffbcid * 25.0 - (fftdo - limits[ffch][0]) * 25.0 / (limits[ffch][1] - limits[ffch][0]); // 'auto' limits

                  straw31_vs_sci0_all->Fill(t31 - fft);

                  if (t30 != 0)
                  {
                     straw31_vs_straw30_banana_all->Fill(t31 - fft, t30 - fft);
                  }

                  if (abs(t31 - fft) < minTsci0)
                  {
                     minTsci0 = abs(t31 - fft);
                     sciT_ch0 = fft;
                  }
               }
            }

            if (sciT_ch0 != 0)
            {
               straw31_vs_sci0->Fill(t31 - sciT_ch0);
            }

            if (t30 != 0 && sciT_ch0 != 0)
            {
               straw31_vs_straw30_banana->Fill(t31 - sciT_ch0, t30 - sciT_ch0);
            }

            // ============================= end of sci 0 correlation finding =============================

            // ========================         LOOP OVER 40 events around         ========================
            //                           jentry to find correlation with sci 1

            // mbytes = 0, mb = 0;
            // for (Long64_t kentry = jentry - 20; kentry < jentry + 20; kentry++)
            // {
            //    Long64_t iientry = LoadTree(kentry);
            //    if (iientry < 0)
            //       continue;
            //    mb = fChain->GetEntry(kentry);
            //    mbytes += mb;

            //    for (int k = 0; k < channel->at(0).size(); k++)
            //    {
            //       int ffch = channel->at(0).at(k);
            //       if (ffch != 1)
            //          continue;

            //       int ffpdo = pdo->at(0).at(k);
            //       int fftdo = tdo->at(0).at(k);
            //       int ffbcid = bcid->at(0).at(k);
            //       double fft = ffbcid * 25.0 - (fftdo - 64) * 25.0 / (85 - 64); //'hand' limits
            //       // double fft = ffbcid * 25.0 - (fftdo - limits[ffch][0]) * 25.0 / (limits[ffch][1] - limits[ffch][0]); // 'auto' limits

            //       if (abs(t31 - fft) < minTsci1)
            //       {
            //          minTsci1 = abs(t31 - fft);
            //          sciT_ch1 = fft;
            //       }
            //    }
            // }

            // if (sciT_ch1 != 0)
            // {
            //    straw31_vs_sci1->Fill(t31 - sciT_ch1);
            // }

            // if (sciT_ch1 != 0 && sciT_ch0 != 0)
            // {
            //    sci0_vs_sci1->Fill(sciT_ch0 - sciT_ch1);
            // }

            // ============================= end of sci 1 correlation finding =============================

            // ========================         LOOP OVER 40 events around         ========================
            //                           jentry to find correlation with sci 2

            // mbytes = 0, mb = 0;
            // for (Long64_t kentry = jentry - 20; kentry < jentry + 20; kentry++)
            // {
            //    Long64_t iientry = LoadTree(kentry);
            //    if (iientry < 0)
            //       continue;
            //    mb = fChain->GetEntry(kentry);
            //    mbytes += mb;

            //    for (int k = 0; k < channel->at(0).size(); k++)
            //    {
            //       int ffch = channel->at(0).at(k);
            //       if (ffch != 2)
            //          continue;

            //       int ffpdo = pdo->at(0).at(k);
            //       int fftdo = tdo->at(0).at(k);
            //       int ffbcid = bcid->at(0).at(k);
            //       double fft = ffbcid * 25.0 - (fftdo - 56) * 25.0 / (75 - 56); //'hand' limits
            //       // double fft = ffbcid * 25.0 - (fftdo - limits[ffch][0]) * 25.0 / (limits[ffch][1] - limits[ffch][0]); // 'auto' limits

            //       if (abs(t31 - fft) < minTsci2)
            //       {
            //          minTsci2 = abs(t31 - fft);
            //          sciT_ch2 = fft;
            //       }
            //    }
            // }

            // if (sciT_ch2 != 0)
            // {
            //    straw31_vs_sci2->Fill(t31 - sciT_ch2);
            // }

            // if (sciT_ch2 != 0 && sciT_ch1 != 0 && sciT_ch0 != 0)
            // {
            //    sci1_vs_sci2->Fill(sciT_ch1 - sciT_ch2);
            //    sci0_vs_sci2->Fill(sciT_ch0 - sciT_ch2);
            // }

            // ============================= end of sci 2 correlation finding =============================
         }
         else
         {
            continue;
         }
      }
   }

   // sci0_vs_sci1->Fit("gaus", "Q", "", sci0_vs_sci1->GetBinCenter(sci0_vs_sci1->FindFirstBinAbove(sci0_vs_sci1->GetMaximum() / 10)), 
   //                                    sci0_vs_sci1->GetBinCenter(sci0_vs_sci1->FindLastBinAbove(sci0_vs_sci1->GetMaximum() / 10)));

   // sci0_vs_sci2->Fit("gaus", "Q", "", sci0_vs_sci2->GetBinCenter(sci0_vs_sci2->FindFirstBinAbove(sci0_vs_sci2->GetMaximum() / 10)), 
   //                                    sci0_vs_sci2->GetBinCenter(sci0_vs_sci2->FindLastBinAbove(sci0_vs_sci2->GetMaximum() / 10)));

   // sci1_vs_sci2->Fit("gaus", "Q", "", sci1_vs_sci2->GetBinCenter(sci1_vs_sci2->FindFirstBinAbove(sci1_vs_sci2->GetMaximum() / 10)), 
   //                                    sci1_vs_sci2->GetBinCenter(sci1_vs_sci2->FindLastBinAbove(sci1_vs_sci2->GetMaximum() / 10)));


   tdo_sci0->Write();
   tdo_sci1->Write();
   tdo_sci2->Write();
   tdo_straw31->Write();
   tdo_vs_pdo_straw31->Write();
   straw31_vs_sci0->Write();
   straw31_vs_sci1->Write();
   straw31_vs_sci2->Write();
   sci0_vs_sci1->Write();
   sci0_vs_sci2->Write();
   sci1_vs_sci2->Write();
   straw31_vs_straw30->Write();
   straw31_vs_straw30_banana->Write();
   straw31_vs_sci0_all->Write();
   straw31_vs_straw30_all->Write();
   straw31_vs_straw30_banana_all->Write();
   out->Close();
}
