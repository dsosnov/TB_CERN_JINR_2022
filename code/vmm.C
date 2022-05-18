#define vmm_cxx
#include "vmm.h"
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>

void vmm::Loop()
{

  int strawMin = -1, strawMax = -1, mmMin = -1, mmMax = -1;
  for(auto &s: channelMap){
    if(s.second.first == 1){
      if(strawMin < 0 || strawMin > s.second.second)
        strawMin = s.second.second;
      if(strawMax < 0 || strawMax < s.second.second)
        strawMax = s.second.second;
    } else if(s.second.first == 4){
      if(mmMin < 0 || mmMin > s.second.second)
        mmMin = s.second.second;
      if(mmMax < 0 || mmMax < s.second.second)
        mmMax = s.second.second;
    }
  }
  printf("Straws: %d-%d, MM: %d-%d\n", strawMin, strawMax, mmMin, mmMax);

   // fast check plots
   auto tdo_sci0 = new TH1D("tdo_sci0", Form("%s: cid_sci0; TDO", file.Data()), 128, 0, 256);
   auto bcid_sci0 = new TH1D("bcid_sci0", Form("%s: bcid_sci0; TDO", file.Data()), 4096, 0, 4096);
   auto tdo_sci1 = new TH1D("tdo_sci1", Form("%s: cid_sci0; TDO", file.Data()), 128, 0, 256);
   auto tdo_sci2 = new TH1D("tdo_sci2", Form("%s: cid_sci0; TDO", file.Data()), 128, 0, 256);
   auto tdo_straw31 = new TH1D("tdo_straw31", Form("%s: tdocid_sci0; TDO", file.Data()), 128, 0, 256); // without 350 PDO counts cut
   auto bcid_straw31 = new TH1D("bcid_straw31", Form("%s: bcid_cid_sci0; TDO", file.Data()), 4096, 0, 4096);
   auto bcid_straw30 = new TH1D("bcid_straw30", Form("%s: bcid_cid_sci0; TDO", file.Data()), 4096, 0, 4096);
   auto tdo_vs_pdo_straw31 = new TH2D("tdo_vs_pdo_straw31", Form("%s: tdo_vs_pdocid_sci0; TDO", file.Data()), 256, 0, 1024, 128, 0, 256);

   // correlation plots
   auto straw31_vs_sci0 = new TH1D("straw31_vs_sci0", Form("%s: straw31_vs_sci0;cid_sci0; TDO", file.Data()), 1000, -500, 500);
   auto straw31_vs_sci0_all = new TH1D("straw31_vs_sci0_all", Form("%s: straw31_vs_sci0_all;cid_sci0; TDO", file.Data()), 1000, -500, 500);
   auto straw31_vs_sci1 = new TH1D("straw31_vs_sci1", Form("%s: straw31_vs_sci1;cid_sci0; TDO", file.Data()), 1000, -500, 500);
   auto straw31_vs_sci2 = new TH1D("straw31_vs_sci2", Form("%s: straw31_vs_sci2;cid_sci0; TDO", file.Data()), 1000, -500, 500);

   auto sci0_vs_sci1 = new TH1D("sci0_vs_sci1", Form("%s: sci0_vs_sci1;cid_sci0; TDO", file.Data()), 250, -50, 50);
   auto sci0_vs_sci2 = new TH1D("sci0_vs_sci2", Form("%s: sci0_vs_sci2;cid_sci0; TDO", file.Data()), 250, -50, 50);
   auto sci1_vs_sci2 = new TH1D("sci1_vs_sci2", Form("%s: sci1_vs_sci2;cid_sci0; TDO", file.Data()), 250, -50, 50);

   auto sci0_vs_mean_triple = new TH1D("sci0_vs_mean_triple", Form("%s: sci0_vs_mean_triple; #Delta t;", file.Data()), 250, -50, 50);
   auto sci1_vs_mean_triple = new TH1D("sci1_vs_mean_triple", Form("%s: sci1_vs_mean_triple; #Delta t;", file.Data()), 250, -50, 50);
   auto sci2_vs_mean_triple = new TH1D("sci2_vs_mean_triple", Form("%s: sci2_vs_mean_triple; #Delta t;", file.Data()), 250, -50, 50);

   auto straw31_vs_straw30 = new TH1D("straw31_vs_straw30", Form("%s: straw31_vs_straw30;cid_sci0; TDO", file.Data()), 1000, -500, 500);
   auto straw31_vs_straw30_all = new TH1D("straw31_vs_straw30_all", Form("%s: straw31_vs_straw30_all;cid_sci0; TDO", file.Data()), 1000, -500, 500);

   auto straw31_vs_straw30_banana_ch0 = new TH2D("straw31_vs_straw30_banana_ch0",
                                              Form("%s: straw31_vs_straw30_banana_ch0; straw 31 #Deltat, ns; straw 3cid_sci0; TDO", file.Data()), 500, -250, 250, 500, -250, 250);

   auto straw31_vs_straw30_banana_ch1 = new TH2D("straw31_vs_straw30_banana_ch1",
                                              Form("%s: straw31_vs_straw30_banana_ch1; straw 31 #Deltat, ns; straw 3cid_sci0; TDO", file.Data()), 500, -250, 250, 500, -250, 250);

   auto straw31_vs_straw30_banana_ch2 = new TH2D("straw31_vs_straw30_banana_ch2",
                                              Form("%s: straw31_vs_straw30_banana_ch2; straw 31 #Deltat, ns; straw 3cid_sci0; TDO", file.Data()), 500, -250, 250, 500, -250, 250);

   auto straw31_vs_straw30_banana_all = new TH2D("straw31_vs_straw30_banana_all",
                                                  Form("%s: straw31_vs_straw30_banana_all; straw 31 #Deltat, ns; straw 3cid_sci0; TDO", file.Data()), 500, -250, 250, 500, -250, 250);

   auto straw31_vs_straw30_banana_bcid = new TH2D("straw31_vs_straw30_banana_bcid",
                                                   Form("%s: straw31_vs_straw30_banana_bcid; straw 31 #Delta BCID; straw 3cid_sci0; TDO", file.Data()), 500, -250, 250, 500, -250, 250);
   // TDO distribution for every Ch
   TH1D *h[64];
   for (Int_t i = 0; i < 64; i++)
   {
     h[i] = new TH1D(Form("h_ch%d", i), Form("%s: TDO of ch%d", file.Data(), i), 128, 0, 256);
   }

   // BCID check
   TH1D *h_bcid_31[121];
   TH1D *h_bcid_30[121];
   int index = 0;

   for (int i = -5; i < 6; i++)
   {
      for (int j = -5; j < 6; j++)
      {
        h_bcid_31[index] = new TH1D(Form("ch_31_h_bcid%d", index), Form("%s: #Delta BCID_{30} = %d and #Delta BCID_{31} = %d; PDO_{31}", file.Data(), i, j), 256, 0, 1024);
        h_bcid_30[index] = new TH1D(Form("ch_30_h_bcid%d", index), Form("%s: #Delta BCID_{30} = %d and #Delta BCID_{31} = %d; PDO_{30}", file.Data(), i, j), 256, 0, 1024);
        index++;
      }
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


   TFile *out = new TFile("../out/out_" + file + ending, "RECREATE"); // PATH where to save out_*.root file
   TDirectory *tdo_dir = out->mkdir("TDO");
   tdo_dir->cd();

   // auto gausFitF = new TF1("gausFitF", "gaus", 0, 256);
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
      // h[i]->Write();

      // std::cout << "CH " << i << "\t L TDO: " << bin1 << "\t R TDO: " << bin2 << "\n";
   }
   // ===================================================================================
   tdo_dir->cd("../");
   // out->Close();
   // return 0;
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
      int straw_bcid_ch31 = 0;
      int straw_pdo_ch31 = 0;

      int ch31_max_pdo = 0;
      int ch30_max_pdo = 0;
      int ch0_max_pdo = 0;
      int ch1_max_pdo = 0;
      int ch2_max_pdo = 0;

      for (int j = 0; j < channel->at(0).size(); j++)
      {
         int fch = channel->at(0).at(j);
         if (fch == 35 || fch == 63)
            continue; // remove 'bad' ch for future tasks

         if (fch == 26)
         {
            int fpdoUC = pdo->at(0).at(j); // Uncorrected PDO, used at time calibration
            int fpdo = correctPDO(fch, fpdoUC);
            // if (fpdo < 100 || fpdo > 900)
            //    continue;
            int ftdo = tdo->at(0).at(j);
            int fbcid = grayDecoded->at(0).at(j);
            // if (fbcid < 40)
            //          continue;

            straw_bcid_ch31 = fbcid;
            straw_pdo_ch31 = fpdo;
            t31 = getTime(fch, fbcid, ftdo, fpdoUC); // 'auto' limits
            // t31 = getTimeByHand(fbcid, ftdo, Y, Y); //'hand' limits

            Long64_t mbytes = 0, mb = 0;
            double minTsci0 = 1e3, minTsci0_ch30 = 1e3;
            double minTsci1 = 1e3, minTsci1_ch30 = 1e3;;
            double minTsci2 = 1e3, minTsci2_ch30 = 1e3;;
            double minT30 = 1e3;

            double sciT_ch0 = 0, sciT_ch1 = 0, sciT_ch2 = 0;
            double t30 = 0;

            int sci_bcid_ch0 = 0, straw_bcid_ch30 = 0;
            int sci_pdo_ch0 = 0, straw_pdo_ch30 = 0;

            // ========================         LOOP OVER 40 events around         ========================
            //                           jentry to find correlation with straw 30

            for (Long64_t kentry = jentry - 1; kentry < jentry + 1; kentry++)
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

                  int ffpdoUC = pdo->at(0).at(k); // Uncorrected PDO, used at time calibration
                  int ffpdo = correctPDO(ffch, ffpdoUC);
                  // if (ffpdo < 100 || ffpdo > 900)
                  //    continue;
                  int fftdo = tdo->at(0).at(k);
                  int ffbcid = grayDecoded->at(0).at(k);
                  // if (ffbcid < 40)
                  //    continue;
                  double fft = getTime(ffch, ffbcid, fftdo, ffpdoUC); // 'auto' limits
                  // double fft = getTimeByHand(ffbcid, fftdo, X, Y); //'hand' limits

                  straw31_vs_straw30_all->Fill(t31 - fft);
                  if (fabs(t31 - fft) < minT30)
                  {
                     minT30 = fabs(t31 - fft);
                     t30 = fft;
                     straw_bcid_ch30 = ffbcid;
                     straw_pdo_ch30 = ffpdo;
                  }
               }
            }
            // ============================ end of straw 30 correlation finding ===========================

            // ========================         LOOP OVER 40 events around         ========================
            //                           jentry to find correlation with sci 0

            mbytes = 0, mb = 0;
            for (Long64_t kentry = jentry - 1; kentry < jentry + 1; kentry++)
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

                  int ffpdoUC = pdo->at(0).at(k); // Uncorrected PDO, used at time calibration
                  int ffpdo = correctPDO(ffch, ffpdoUC);
                  int fftdo = tdo->at(0).at(k);
                  int ffbcid = grayDecoded->at(0).at(k);
                  // if (ffbcid < 40)
                  //    continue;
                  double fft = getTimeByHand(ffbcid, fftdo, 88, 140); //'hand' limits
                  // double fft = getTime(ffch, ffbcid, fftdo, ffpdoUC); // 'auto' limits

                  straw31_vs_sci0_all->Fill(t31 - fft);

                  if (t30 != 0)
                  {
                     straw31_vs_straw30_banana_all->Fill(t31 - fft, t30 - fft);
                  }

                  if (fabs(t31 - fft) < minTsci0 && fabs(t30 - fft) < minTsci0_ch30)
                  {
                     minTsci0 = fabs(t31 - fft);
                     minTsci0_ch30 = fabs(t30 - fft);
                     sciT_ch0 = fft;
                     sci_bcid_ch0 = ffbcid;
                  }
               }
            }

            if (sciT_ch0 != 0)
            {
               straw31_vs_sci0->Fill(t31 - sciT_ch0);
               bcid_sci0->Fill(sci_bcid_ch0);
               bcid_straw31->Fill(straw_bcid_ch31);
               bcid_straw30->Fill(straw_bcid_ch30);
            }

            if (t30 != 0 && sciT_ch0 != 0)
            {
               straw31_vs_straw30_banana_ch0->Fill(t31 - sciT_ch0, t30 - sciT_ch0);
               straw31_vs_straw30_banana_bcid->Fill(straw_bcid_ch31 - sci_bcid_ch0, straw_bcid_ch30 - sci_bcid_ch0);
               straw31_vs_sci0->Fill(t31 - sciT_ch0);
               bcid_sci0->Fill(sci_bcid_ch0);
               bcid_straw31->Fill(straw_bcid_ch31);
               bcid_straw30->Fill(straw_bcid_ch30);
               straw31_vs_straw30->Fill(t30 - t31);

               int findex = 0;
               for (int ii = -5; ii < 6; ii++)
               {
                  for (int jj = -5; jj < 6; jj++)
                  {
                     if (straw_bcid_ch30 - sci_bcid_ch0 == ii && straw_bcid_ch31 - sci_bcid_ch0 == jj)
                     {
                        h_bcid_31[findex]->Fill(straw_pdo_ch31);
                        h_bcid_30[findex]->Fill(straw_pdo_ch30);
                     }
                     findex++;
                  }
               }
            }

            // ============================= end of sci 0 correlation finding =============================

            // ========================         LOOP OVER 40 events around         ========================
            //                           jentry to find correlation with sci 1

            mbytes = 0, mb = 0;
            for (Long64_t kentry = jentry - 1; kentry < jentry + 1; kentry++)
            {
               Long64_t iientry = LoadTree(kentry);
               if (iientry < 0)
                  continue;
               mb = fChain->GetEntry(kentry);
               mbytes += mb;

               for (int k = 0; k < channel->at(0).size(); k++)
               {
                  int ffch = channel->at(0).at(k);
                  if (ffch != 1)
                     continue;

                  int ffpdoUC = pdo->at(0).at(k); // Uncorrected PDO, used at time calibration
                  int ffpdo = correctPDO(ffch, ffpdoUC);
                  int fftdo = tdo->at(0).at(k);
                  int ffbcid = grayDecoded->at(0).at(k);
                  double fft = getTimeByHand(ffbcid, fftdo, 110, 160); //'hand' limits
                  // double fft = getTime(ffch, ffbcid, fftdo, ffpdoUC); // 'auto' limits

                  if (fabs(t31 - fft) < minTsci1 && fabs(t30 - fft) < minTsci1_ch30)
                  {
                     minTsci1 = fabs(t31 - fft);
                     minTsci1_ch30 = fabs(t30 - fft);
                     sciT_ch1 = fft;
                  }
               }
            }

            if (sciT_ch1 != 0 && t30 != 0)
            {
               straw31_vs_sci1->Fill(t31 - sciT_ch1);
               straw31_vs_straw30_banana_ch1->Fill(t31 - sciT_ch1, t30 - sciT_ch1);
            }

            if (sciT_ch1 != 0 && sciT_ch0 != 0)
            {
               sci0_vs_sci1->Fill(sciT_ch0 - sciT_ch1);
            }

            // ============================= end of sci 1 correlation finding =============================

            // ========================         LOOP OVER 40 events around         ========================
            //                           jentry to find correlation with sci 2

            mbytes = 0, mb = 0;
            for (Long64_t kentry = jentry - 1; kentry < jentry + 1; kentry++)
            {
               Long64_t iientry = LoadTree(kentry);
               if (iientry < 0)
                  continue;
               mb = fChain->GetEntry(kentry);
               mbytes += mb;

               for (int k = 0; k < channel->at(0).size(); k++)
               {
                  int ffch = channel->at(0).at(k);
                  if (ffch != 2)
                     continue;

                  int ffpdoUC = pdo->at(0).at(k); // Uncorrected PDO, used at time calibration
                  int ffpdo = correctPDO(ffch, ffpdoUC);
                  int fftdo = tdo->at(0).at(k);
                  int ffbcid = grayDecoded->at(0).at(k);
                  double fft = getTimeByHand(ffbcid, fftdo, 96, 148); //'hand' limits
                  // double fft = getTime(ffch, ffbcid, fftdo, ffpdoUC); // 'auto' limits

                  if (fabs(t31 - fft) < minTsci2 && fabs(t30 - fft) < minTsci2_ch30)
                  {
                     minTsci2 = fabs(t31 - fft);
                     minTsci2_ch30 = fabs(t30 - fft);
                     sciT_ch2 = fft;
                  }
               }
            }

            if (sciT_ch2 != 0 && t30 != 0)
            {
               straw31_vs_sci2->Fill(t31 - sciT_ch2);
               straw31_vs_straw30_banana_ch2->Fill(t31 - sciT_ch2, t30 - sciT_ch2);
            }

            if (sciT_ch2 != 0)
            {
              if(sciT_ch1 != 0)
                sci1_vs_sci2->Fill(sciT_ch1 - sciT_ch2);
              if(sciT_ch0 != 0)
                sci0_vs_sci2->Fill(sciT_ch0 - sciT_ch2);
            }
            if(sciT_ch0 && sciT_ch1 && sciT_ch2)
            {
              auto meanTime = (sciT_ch0 + sciT_ch1 + sciT_ch2) / 3.0;
              sci0_vs_mean_triple->Fill(sciT_ch0 - meanTime);
              sci1_vs_mean_triple->Fill(sciT_ch1 - meanTime);
              sci2_vs_mean_triple->Fill(sciT_ch2 - meanTime);
            }

            // ============================= end of sci 2 correlation finding =============================
         }
         else
         {
            continue;
         }
      }
   }

   TDirectory *pdo_dir = out->mkdir("PDO");
   pdo_dir->cd();
   for (int i = 0; i < 121; i++)
   {
      h_bcid_31[i]->Write();
      h_bcid_30[i]->Write();
   }

   pdo_dir->cd("../");

   // sci0_vs_sci1->Fit("gaus", "Q", "", sci0_vs_sci1->GetBinCenter(sci0_vs_sci1->FindFirstBinAbove(sci0_vs_sci1->GetMaximum() / 10)),
   //                                    sci0_vs_sci1->GetBinCenter(sci0_vs_sci1->FindLastBinAbove(sci0_vs_sci1->GetMaximum() / 10)));

   // sci0_vs_sci2->Fit("gaus", "Q", "", sci0_vs_sci2->GetBinCenter(sci0_vs_sci2->FindFirstBinAbove(sci0_vs_sci2->GetMaximum() / 10)),
   //                                    sci0_vs_sci2->GetBinCenter(sci0_vs_sci2->FindLastBinAbove(sci0_vs_sci2->GetMaximum() / 10)));

   // sci1_vs_sci2->Fit("gaus", "Q", "", sci1_vs_sci2->GetBinCenter(sci1_vs_sci2->FindFirstBinAbove(sci1_vs_sci2->GetMaximum() / 10)),
   //                                    sci1_vs_sci2->GetBinCenter(sci1_vs_sci2->FindLastBinAbove(sci1_vs_sci2->GetMaximum() / 10)));

   tdo_sci0->Write();
   bcid_sci0->Write();
   tdo_sci1->Write();
   tdo_sci2->Write();
   tdo_straw31->Write();
   tdo_vs_pdo_straw31->Write();
   bcid_straw31->Write();
   bcid_straw30->Write();
   straw31_vs_sci0->Write();
   straw31_vs_sci1->Write();
   straw31_vs_sci2->Write();
   sci0_vs_sci1->Write();
   sci0_vs_sci2->Write();
   sci1_vs_sci2->Write();
   sci0_vs_mean_triple->Write();
   sci1_vs_mean_triple->Write();
   sci2_vs_mean_triple->Write();
   straw31_vs_straw30->Write();
   straw31_vs_straw30_banana_ch0->Write();
   straw31_vs_straw30_banana_ch1->Write();
   straw31_vs_straw30_banana_ch2->Write();
   straw31_vs_sci0_all->Write();
   straw31_vs_straw30_all->Write();
   straw31_vs_straw30_banana_all->Write();
   straw31_vs_straw30_banana_bcid->Write();
   out->Close();
}
