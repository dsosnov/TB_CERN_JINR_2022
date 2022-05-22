#define evBuilder_cxx
#include "evBuilder.h"
#include <TH1.h>

void evBuilder::threePlotDrawF(TH1D *h1, TH1D *h2, TH1D *h3, TString fileEnding)
{
    h1->SetLineColor(kGreen - 2);
    h2->SetLineColor(kMagenta);
    h3->SetLineColor(kBlack);
    h1->SetStats(0);
    h2->SetStats(0);
    h3->SetStats(0);


    TCanvas *three_plots = new TCanvas("3 det correlation", "3 det correlation", 1400, 900);
    three_plots->cd();

    h2->Draw();
    h1->Draw("SAME");
    h3->Draw("SAME");

    h1->Fit("gaus", "", "", -200, 200); // TMM - TScint
    h2->Fit("gaus", "", "", -200, 200); // TStraw - TScint
    h3->Fit("gaus", "", "", -200, 200); // TStraw - TMM

    TF1 *g1 = (TF1*)h1->GetListOfFunctions()->FindObject("gaus");
    TF1 *g2 = (TF1*)h2->GetListOfFunctions()->FindObject("gaus");
    TF1 *g3 = (TF1*)h3->GetListOfFunctions()->FindObject("gaus");

    double m1 = g1->GetParameter(1);
    double s1 = g1->GetParameter(2);

    double m2 = g2->GetParameter(1);
    double s2 = g2->GetParameter(2);

    double m3 = g3->GetParameter(1);
    double s3 = g3->GetParameter(2);


    gStyle->SetOptFit(1111);
    // draw fit parameters as legends:
    Char_t ndf[80];
    Char_t sigma[80];
    Char_t mean[80];
    Char_t constant[80];
    auto legend = new TLegend(0.65, 0.9, 1.0, 0.75, "TMM - TScint");
    sprintf(ndf, "#chi^{2}/NDF = %.2f / %.2i", h1->GetFunction("gaus")->GetChisquare(), h1->GetFunction("gaus")->GetNDF());
    legend->AddEntry(h1, ndf);
    sprintf(sigma, "#sigma = %.1f", h1->GetFunction("gaus")->GetParameter(2));
    legend->AddEntry(h1, sigma);
    sprintf(mean, "Mean = %.1f", h1->GetFunction("gaus")->GetParameter(1));
    legend->AddEntry(h1, mean);
    sprintf(constant, "Events under peak: %.f", h1->GetFunction("gaus")->Integral(m1 - 4 * s1, m1 + 4 * s1) / h1->GetBinWidth(1));
    legend->AddEntry(h1, constant);
    legend->Draw("same");

    Char_t ndf1[80];
    Char_t sigma1[80];
    Char_t mean1[80];
    Char_t constant1[80];
    auto legend1 = new TLegend(0.65, 0.75, 1.0, 0.60, "TStraw - TScint");
    sprintf(ndf1, "#chi^{2}/NDF = %.2f / %.2i", h2->GetFunction("gaus")->GetChisquare(), h2->GetFunction("gaus")->GetNDF());
    legend1->AddEntry(h2, ndf1);
    sprintf(sigma1, "#sigma = %.1f", h2->GetFunction("gaus")->GetParameter(2));
    legend1->AddEntry(h2, sigma1);
    sprintf(mean1, "Mean = %.1f", h2->GetFunction("gaus")->GetParameter(1));
    legend1->AddEntry(h2, mean1);
    sprintf(constant1, "Events under peak: %.f", h2->GetFunction("gaus")->Integral(m2 - 4 * s2, m2 + 4 * s2) / h2->GetBinWidth(1));
    legend1->AddEntry(h2, constant1);
    legend1->Draw("same");

    Char_t ndf2[80];
    Char_t sigma2[80];
    Char_t mean2[80];
    Char_t constant2[80];
    auto legend2 = new TLegend(0.65, 0.60, 1.0, 0.45, "TStraw - TMM");
    sprintf(ndf2, "#chi^{2}/NDF = %.2f / %.2i", h3->GetFunction("gaus")->GetChisquare(), h3->GetFunction("gaus")->GetNDF());
    legend2->AddEntry(h3, ndf2);
    sprintf(sigma2, "#sigma = %.1f", h3->GetFunction("gaus")->GetParameter(2));
    legend2->AddEntry(h3, sigma2);
    sprintf(mean2, "Mean = %.1f", h3->GetFunction("gaus")->GetParameter(1));
    legend2->AddEntry(h3, mean2);
    sprintf(constant2, "Events under peak: %.f", h3->GetFunction("gaus")->Integral(m3 - 4 * s3, m3 + 4 * s3) / h3->GetBinWidth(1));
    legend2->AddEntry(h3, constant2);
    legend2->Draw("same");

    three_plots->SaveAs("../out/3plots_" + file + fileEnding + ".pdf");
}

void evBuilder::Loop()
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

   if (fChain == 0)
        return;

    // map<int, float> firstStripForStraw = {
    //   {24, 21 + 150 + 24/1.0},
    //   {25, 21 + 150 + 24/2.0},
    //   {26, 21 + 150         }, //21 + (mmMin - 4)
    //   {27, 21 + 150 - 24/2.0},
    //   {28, 21 + 150 + 24/1.0},
    //   {29, 21 + 150 + 24/1.5},
    // };
    map<int, float> strawCenterMM = {
      {24, 207}, // 21 + 150 + 24/1.0 + 12
      {25, 198}, // 21 + 150 + 24/2.0 + 15
      {26, 182}, // 21 + 150 + 12 - 1 // 21 + (mmMin - 4)
      {27, 174}, // 21 + 150 - 24/2.0 + 15
      {28, 207}, // 21 + 150 + 24/1.0 + 12
      {29, 199}, // 21 + 150 + 24/1.5 + 12
    };

    TFile *out = new TFile("../out/out_" + file + ending, "RECREATE"); // PATH where to save out_*.root file

    auto straw_vs_sci = new TH1D("straw_vs_sci", Form("%s: straw vs scint;#Deltat, ns", file.Data()), 1000, -500, 500);
    auto straw_vs_mm = new TH1D("straw_vs_mm", Form("%s: straw vs microMegas;#Deltat, ns", file.Data()), 1000, -500, 500);
    auto hits_in_cluster = new TH1D("hits_in_cluster", Form("%s: hits in microMegas per cluster;N", file.Data()), 100, 0, 100);
    auto mm_vs_sci = new TH1D("mm_vs_sci", Form("%s: microMegas vs scint;#Deltat, ns", file.Data()), 1000, -500, 500);

    auto sci0_vs_sci60 = new TH1D("sci0_vs_sci60", Form("%s: sci ch 0 vs sci ch 60;#Deltat, ns", file.Data()), 1000, -500, 500);

    auto straw_vs_sci_3det_corr = new TH1D("straw_vs_sci_3det_corr", Form("%s: 3-det correlations;#Deltat, ns", file.Data()), 1000, -500, 500);
    auto straw_vs_mm_3det_corr = new TH1D("straw_vs_mm_3det_corr", Form("%s: 3-det correlations;#Deltat, ns", file.Data()), 1000, -500, 500);
    auto mm_vs_sci_3det_corr = new TH1D("mm_vs_sci_3det_corr", Form("%s: 3-det correlations;#Deltat, ns", file.Data()), 1000, -500, 500);
    auto straw_vs_sci_3det_corr_0 = new TH1D("straw_vs_sci0_3det_corr", Form("%s: 3-det correlations, sci0;#Deltat, ns", file.Data()), 1000, -500, 500);
    auto straw_vs_mm_3det_corr_0 = new TH1D("straw_vs_mm_3det_corr_vs_sci0", Form("%s: 3-det correlations, sci0;#Deltat, ns", file.Data()), 1000, -500, 500);
    auto mm_vs_sci_3det_corr_0 = new TH1D("mm_vs_sci0_3det_corr", Form("%s: 3-det correlations, sci0;#Deltat, ns", file.Data()), 1000, -500, 500);

    auto straw_vs_mm_spatial_corr = new TH2D("straw_vs_mm_spatial_corr", Form("%s: microMegas vs straw spatial correaltion;straw ch;MM ch", file.Data()),
                                             strawMax - strawMin + 1, strawMin, strawMax + 1, mmMax - mmMin + 1, mmMin, mmMax);

    auto mm_corr0 = new TH1D("mm_corr0", Form("%s: mm_corr0;MM ch", file.Data()), mmMax - mmMin + 1, mmMin, mmMax);
    auto mm_corr0_straw = new TH1D("mm_corr0_straw", Form("%s: mm_corr0_straw;MM ch", file.Data()), mmMax - mmMin + 1, mmMin, mmMax);
    auto mmSingle_corr0 = new TH1D("mmSingle_corr0", Form("%s: mmSingle_corr0;MM ch", file.Data()), mmMax - mmMin + 1, mmMin, mmMax);
    auto mmSingle_corr0_straw = new TH1D("mmSingle_corr0_straw", Form("%s: mmSingle_corr0_straw;MM ch", file.Data()), mmMax - mmMin + 1, mmMin, mmMax);

    auto straw_rt_dir = out->mkdir("straw_rt");
    straw_rt_dir->cd();
    map<int, TH2D*> straw_rt, straw_rt_0 ;
    for(auto i = strawMin; i <= strawMax; i++){
      straw_rt.emplace(i,
                       new TH2D(Form("straw%d_rt", i),
                                Form("%s: straw %d v-shape sci ch 60;R, mm;T, ns", file.Data(), i),
                                32, -4, 4, 300, -100, 200));
      straw_rt_0.emplace(i,
                         new TH2D(Form("straw%d_rt_0", i),
                                  Form("%s: straw %d v-shape sci ch 0;R, mm;T, ns", file.Data(), i),
                                  32, -4, 4, 300, -100, 200));
    }
    out->cd();

    auto straw_pdo_dir = out->mkdir("straw_pdo_corr");
    straw_pdo_dir->cd();
    map<int, TH1D*> straw_pdo, straw_pdo_0 ;
    for(auto i = strawMin; i <= strawMax; i++){
      straw_pdo.emplace(i,
                        new TH1D(Form("straw%d_pdo_corr", i),
                                 Form("%s: pdo for straw %d corellated with sci ch 60;pdo", file.Data(), i), 64, 0, 1024));
      straw_pdo_0.emplace(i,
                          new TH1D(Form("straw%d_pdo_corr_0", i),
                                   Form("%s: pdo for straw %d corellated with sci ch 0;pdo", file.Data(), i), 64, 0, 1024));
    }
    out->cd();

    auto straw_deltat_dir = out->mkdir("straw_deltat_corr");
    straw_deltat_dir->cd();
    map<int, TH1D*> straw_deltat, straw_deltat_0 ;
    for(auto i = strawMin; i <= strawMax; i++){
      straw_deltat.emplace(i,
                        new TH1D(Form("straw%d_vs_sci60", i),
                                 Form("%s: straw%d_vs_sci60;#Delta t", file.Data(), i), 1000, -500, 500));
      straw_deltat_0.emplace(i,
                        new TH1D(Form("straw%d_vs_sci0", i),
                                 Form("%s: straw%d_vs_sci0;#Delta t", file.Data(), i), 1000, -500, 500));
    }
    out->cd();

    auto straw_pdo_uc_dir = out->mkdir("straw_pdo_uncorr");
    straw_pdo_uc_dir->cd();
    map<int, TH1D*> straw_pdo_ucl_ucr, straw_pdo_ucl_cr;
    for(auto i = strawMin; i <= strawMax; i++){
      straw_pdo_ucl_ucr.emplace(i,
                          new TH1D(Form("straw%d_pdo_uncorrellated_uncorrected", i),
                                   Form("%s: pdo for straw %d - uncorellated, no pdo correction;pdo", file.Data(), i), 64, 0, 1024));
      straw_pdo_ucl_cr.emplace(i,
                        new TH1D(Form("straw%d_pdo_uncorrellated_corrected", i),
                                 Form("%s: pdo for straw %d - uncorellated, pdo correction;pdo", file.Data(), i), 64, 0, 1024));
    }
    out->cd();

    auto straw_banana_dir = out->mkdir("straw_banana");
    straw_banana_dir->cd();
    map<int, TH2D*> straw_banana, straw_banana_0 ;
    for(auto i = strawMin; i < strawMax; i++){
      straw_banana.emplace(i,
                       new TH2D(Form("straw%d-%d_banana", i, i+1),
                                Form("%s: Time difference between straws %d, %d and sci60;T_{straw%d} - T_{scint}, [ns];T_{straw%d} - T_{scint}, [ns]", file.Data(), i, i+1, i, i+1),
                                500, -250, 250, 500, -250, 250));
      straw_banana_0.emplace(i,
                       new TH2D(Form("straw%d-%d_banana_0", i, i+1),
                                Form("%s: Time difference between straws %d, %d and sci0;T_{straw%d} - T_{scint}, [ns];T_{straw%d} - T_{scint}, [ns]", file.Data(), i, i+1, i, i+1),
                                  500, -250, 250, 500, -250, 250));
    }
    out->cd();

    auto straw_vs_straw_deltat_dir = out->mkdir("straw_vs_straw_deltat");
    straw_vs_straw_deltat_dir->cd();
    map<int, TH1D*> straw_straw;
    for(auto i = strawMin; i < strawMax; i++){
      straw_straw.emplace(i,
                          new TH1D(Form("straw%d_vs_straw%d", i, i+1),
                                   Form("%s: straw%d_vs_straw%d;T_{straw%d} - T_{straw%d}", file.Data(), i, i+1, i, i+1), 1000, -500, 500));
    }
    out->cd();
    
    unsigned int nLoopEntriesAround = 0;
    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;

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

        double t_srtraw = 0;
        int straw_bcid_ch_srtraw = 0;
        int straw_pdo_ch_srtraw = 0;


        for (int j = 0; j < channel->at(0).size(); j++)
        {
            int fch = channel->at(0).at(j);
            int fchD = getMappedDetector(fch);
            int fchM = getMappedChannel(fch);

            if (fchD == 1) // All straw ch
            {
                int fpdoUC = pdo->at(0).at(j); // Uncorrected PDO, used at time calibration
                int fpdo = correctPDO(fch, fpdoUC);
                int ftdo = tdo->at(0).at(j);
                int fbcid = grayDecoded->at(0).at(j);

                straw_bcid_ch_srtraw = fbcid;
                straw_pdo_ch_srtraw = fpdo;
                t_srtraw = getTime(fch, fbcid, ftdo, fpdoUC); // 'auto' limits
                // t_srtraw = getTimeByHand(fbcid, ftdo, Y, Y); //'hand' limits

                straw_pdo_ucl_cr.at(fchM)->Fill(fpdo);
                straw_pdo_ucl_ucr.at(fchM)->Fill(fpdoUC);

                Long64_t mbytes = 0, mb = 0;
                double t30 = 0;
                double minTsci0 = 1e3;
                double sciT_ch0 = 0;
                int sci_bcid_ch0 = 0;
                double minTsci60 = 1e3;
                double sciT_ch60 = 0;
                int sci_bcid_ch60 = 0;
                vector<array<double, 3> > MmCluster;
                double neighborStrawTime = 0;
                double neighborMinStrawTime = 1E3;

                // ========================         LOOP OVER nLoopEntriesAround  events around         ========================
                //                              jentry to find correlation with MM
                MmCluster.clear();
                mbytes = 0, mb = 0;

                for (Long64_t kentry = jentry - nLoopEntriesAround; kentry <= jentry + nLoopEntriesAround; kentry++)
                {
                    Long64_t iientry = LoadTree(kentry);
                    if (iientry < 0)
                        continue;
                    mb = fChain->GetEntry(kentry);
                    mbytes += mb;

                    for (int k = 0; k < channel->at(0).size(); k++)
                    {
                      int ffch = channel->at(0).at(k);
                      int ffchD = getMappedDetector(ffch);
                      int ffchM = getMappedChannel(ffch);
                      if(ffchD == 4) // All MM channels
                      { 
                        int ffpdoUC = pdo->at(0).at(k); // Uncorrected PDO, used at time calibration
                        int ffpdo = correctPDO(ffch, ffpdoUC);
                        int fftdo = tdo->at(0).at(k);
                        int ffbcid = grayDecoded->at(0).at(k);
                        // double fft = getTimeByHand(ffbcid, fftdo, 110, 160); //'hand' limits
                        double fft = getTime(ffch, ffbcid, fftdo, ffpdoUC); // 'auto' limits

                        if (fabs(t_srtraw - fft) < 500)
                        {
                          straw_vs_mm ->Fill(t_srtraw - fft);
                          array<double, 3> mM_hit = {{ffchM * 1.0, ffpdo * 1.0, fft}};
                          MmCluster.push_back(mM_hit);
                          mm_corr0_straw->Fill(ffchM);
                        }
                      }
                      else if (ffchD == 0 && ffchM == 0) // Sci 0
                      {
                        int ffpdoUC = pdo->at(0).at(k); // Uncorrected PDO, used at time calibration
                        int ffpdo = correctPDO(ffch, ffpdoUC);
                        int fftdo = tdo->at(0).at(k);
                        int ffbcid = grayDecoded->at(0).at(k);
                        // if (ffbcid < 40)
                        //    continue;
                        double fft = getTimeByHand(ffbcid, fftdo, 88, 140); //'hand' limits
                        // double fft = getTime(ffch, ffbcid, fftdo, ffpdoUC); // 'auto' limits

                        // straw_vs_sci->Fill(t_srtraw - fft);

                        if (fabs(t_srtraw - fft) < minTsci0)
                        {
                          minTsci0 = fabs(t_srtraw - fft);
                          sciT_ch0 = fft;
                          sci_bcid_ch0 = ffbcid;
                        }

                      }
                      else if (ffchD == 0 && ffchM == 3) // triple sci coinsidence
                      {
                        int ffpdoUC = pdo->at(0).at(k); // Uncorrected PDO, used at time calibration
                        int ffpdo = correctPDO(ffch, ffpdoUC);
                        int fftdo = tdo->at(0).at(k);
                        int ffbcid = grayDecoded->at(0).at(k);
                        // if (ffbcid < 40)
                        //    continue;
                        // double fft = getTimeByHand(ffbcid, fftdo, 88, 140); //'hand' limits
                        double fft = getTime(ffch, ffbcid, fftdo, ffpdoUC); // 'auto' limits

                        // straw_vs_sci->Fill(t_srtraw - fft);

                        if (fabs(t_srtraw - fft) < minTsci60)
                        {
                          minTsci60 = fabs(t_srtraw - fft);
                          sciT_ch60 = fft;
                          sci_bcid_ch60 = ffbcid;
                        }

                      }
                      else if (ffchD == 1 && ffchM == fchM + 1) // Next straw channel
                      {
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
                        if (fabs(t_srtraw - fft) < neighborMinStrawTime)
                        {
                          neighborMinStrawTime = fabs(t_srtraw - fft);
                          neighborStrawTime = fft;
                        }
                      }
                      else
                      {
                        continue;
                      }
                    }
                }

                double meanT = 0;
                double meanCh = 0;
                double sum1 = 0;
                double sum2 = 0;
                double w_sum = 0;

                double minT_straw_mm = 600;
                int mmCh_min = 0;

                if (MmCluster.size() != 0)
                {
                    for (size_t l = 0; l < MmCluster.size(); l++)
                    {
                        if (abs(MmCluster.at(l)[2] - t_srtraw) < minT_straw_mm)
                        {
                            minT_straw_mm = abs(MmCluster.at(l)[2] - t_srtraw);
                            mmCh_min = MmCluster.at(l)[0];
                        }
                    }

                    for (size_t l = 0; l < MmCluster.size(); l++)
                    {
                        if (abs(MmCluster.at(l)[0] - mmCh_min) > 5)
                        {
                            MmCluster.erase(MmCluster.begin()+l);
                        }
                    }
    

                    for (size_t l = 0; l < MmCluster.size(); l++)
                    {
                        sum1 += MmCluster.at(l)[2] * MmCluster.at(l)[1] / 1024.0;
                        sum2 += MmCluster.at(l)[0] * MmCluster.at(l)[1] / 1024.0;
                        w_sum += MmCluster.at(l)[1] / 1024.0;
                    }
                    meanT = sum1 / w_sum;
                    meanCh = sum2 / w_sum;
                    straw_vs_mm_spatial_corr->Fill(fchM, meanCh);
                    // straw_vs_mm ->Fill(t_srtraw - meanT);
                    hits_in_cluster->Fill( MmCluster.size());
                    if(MmCluster.size() == 1)
                      mmSingle_corr0_straw->Fill(meanCh);
                }

                // ============================= end of sci MM correlation finding ============================

                if (sciT_ch0 != 0 && sciT_ch60 != 0) // WARNING! this is not real scintillator corellation!
                {
                    sci0_vs_sci60->Fill(sciT_ch0 - sciT_ch60);
                }

                if (sciT_ch0 != 0)
                {
                    if (meanT != 0)
                    {
                        mm_vs_sci->Fill(meanT - sciT_ch0);
                    }
                    straw_vs_sci->Fill(t_srtraw - sciT_ch0);
                    straw_deltat_0.at(fchM)->Fill(t_srtraw - sciT_ch0);
                    straw_pdo_0.at(fchM)->Fill(fpdo);
                }
                if (sciT_ch60 != 0)
                {
                    straw_deltat.at(fchM)->Fill(t_srtraw - sciT_ch60);
                    straw_pdo.at(fchM)->Fill(fpdo);
                }
                if (sciT_ch0 != 0 && meanT != 0)
                {
                    straw_rt_0.at(fchM)->Fill((meanCh - strawCenterMM.at(fchM)) * 0.25, 100 + t_srtraw - sciT_ch0);
                    // if (strawCh == 3 && (meanCh > 21 && meanCh < 47))
                    // {
                    //     straw26_rt_0->Fill((meanCh - 21) * 0.25, 100 + t_srtraw - sciT_ch0);
                    // }
                    straw_vs_sci_3det_corr_0->Fill(t_srtraw - sciT_ch0);
                    straw_vs_mm_3det_corr_0->Fill(t_srtraw - meanT);
                    mm_vs_sci_3det_corr_0->Fill(meanT - sciT_ch0);
                }
                if (sciT_ch60 != 0 && meanT != 0)
                {
                    straw_vs_sci_3det_corr->Fill(t_srtraw - sciT_ch60);
                    straw_vs_mm_3det_corr->Fill(t_srtraw - meanT);
                    mm_vs_sci_3det_corr->Fill(meanT - sciT_ch60);
                    straw_rt.at(fchM)->Fill((meanCh - strawCenterMM.at(fchM)) * 0.25, 100 + t_srtraw - sciT_ch60);
                    // if (strawCh == 3 && (meanCh > 21 && meanCh < 47))
                    // {
                    //     straw26_rt->Fill((meanCh - 21) * 0.25, 100 + t_srtraw - sciT_ch60);
                    // }
                }
                if(neighborStrawTime != 0)
                {
                  // straw_banana.at(fchM)->Fill()
                  straw_straw.at(fchM)->Fill(t_srtraw - neighborStrawTime);
                  if (sciT_ch0 != 0)
                    straw_banana_0.at(fchM)->Fill(t_srtraw - sciT_ch0, neighborStrawTime - sciT_ch0);
                  if (sciT_ch60 != 0)
                    straw_banana.at(fchM)->Fill(t_srtraw - sciT_ch60, neighborStrawTime - sciT_ch60);
                }

                // ============================= end of sci 0 correlation finding =============================

            }
            else if (fchD == 0 && fchM == 0) // Sci0
            {
                int fpdoUC = pdo->at(0).at(j); // Uncorrected PDO, used at time calibration
                int fpdo = correctPDO(fch, fpdoUC);
                int ftdo = tdo->at(0).at(j);
                int fbcid = grayDecoded->at(0).at(j);

                straw_bcid_ch_srtraw = fbcid;
                straw_pdo_ch_srtraw = fpdo;
                t_srtraw = getTime(fch, fbcid, ftdo, fpdoUC); // 'auto' limits
                // t_srtraw = getTimeByHand(fbcid, ftdo, Y, Y); //'hand' limits

                Long64_t mbytes = 0, mb = 0;
                double t30 = 0;
                double minTsci0 = 1e3;
                double sciT_ch0 = 0;
                int sci_bcid_ch0 = 0;
                double minTsci60 = 1e3;
                double sciT_ch60 = 0;
                int sci_bcid_ch60 = 0;
                vector<array<double, 3> > MmCluster;

                // ========================         LOOP OVER nLoopEntriesAround  events around         ========================
                //                              jentry to find correlation with MM
                MmCluster.clear();
                mbytes = 0, mb = 0;

                for (Long64_t kentry = jentry - nLoopEntriesAround; kentry <= jentry + nLoopEntriesAround; kentry++)
                {
                    Long64_t iientry = LoadTree(kentry);
                    if (iientry < 0)
                        continue;
                    mb = fChain->GetEntry(kentry);
                    mbytes += mb;

                    for (int k = 0; k < channel->at(0).size(); k++)
                    {
                        int ffch = channel->at(0).at(k);
                        int ffchD = getMappedDetector(ffch);
                        int ffchM = getMappedChannel(ffch);
                        if (ffchD != 4)
                            continue;

                        int ffpdoUC = pdo->at(0).at(k); // Uncorrected PDO, used at time calibration
                        int ffpdo = correctPDO(ffch, ffpdoUC);
                        int fftdo = tdo->at(0).at(k);
                        int ffbcid = grayDecoded->at(0).at(k);
                        // double fft = getTimeByHand(ffbcid, fftdo, 110, 160); //'hand' limits
                        double fft = getTime(ffch, ffbcid, fftdo, ffpdoUC); // 'auto' limits

                        if (fabs(t_srtraw - fft) < 500)
                        {
                            array<double, 3> mM_hit = {{ffchM * 1.0, ffpdo * 1.0, fft}};
                            MmCluster.push_back(mM_hit);
                            mm_corr0->Fill(ffchM);
                        }
                    }
                }

                double meanT = 0;
                double meanCh = 0;
                double sum1 = 0;
                double sum2 = 0;
                double w_sum = 0;

                double minT_straw_mm = 600;
                int mmCh_min = 0;

                if (MmCluster.size() != 0)
                {
                    for (size_t l = 0; l < MmCluster.size(); l++)
                    {
                        if (abs(MmCluster.at(l)[2] - t_srtraw) < minT_straw_mm)
                        {
                            minT_straw_mm = abs(MmCluster.at(l)[2] - t_srtraw);
                            mmCh_min = MmCluster.at(l)[0];
                        }
                    }

                    for (size_t l = 0; l < MmCluster.size(); l++)
                    {
                        if (abs(MmCluster.at(l)[0] - mmCh_min) > 5)
                        {
                            MmCluster.erase(MmCluster.begin()+l);
                        }
                    }
    

                    for (size_t l = 0; l < MmCluster.size(); l++)
                    {
                        sum1 += MmCluster.at(l)[2] * MmCluster.at(l)[1] / 1024.0;
                        sum2 += MmCluster.at(l)[0] * MmCluster.at(l)[1] / 1024.0;
                        w_sum += MmCluster.at(l)[1] / 1024.0;
                    }
                    meanT = sum1 / w_sum;
                    meanCh = sum2 / w_sum;
                    if(MmCluster.size() == 1)
                      mmSingle_corr0->Fill(meanCh);
                }

            }
            else
            {
                continue;
            }
        }
    }

    auto straw_rt_normed_dir = out->mkdir("straw_rt_normed");
    straw_rt_normed_dir->cd();
    map<int, TH2D*> straw_rt_normed, straw_rt_0_normed ;
    for(auto &h: straw_rt){
      auto hnew = static_cast<TH2D*>(h.second->Clone(Form("straw%d_rt_normed", h.first)));
      for(auto i = 1; i <= hnew->GetNbinsX(); i++){
        auto integ = hnew->Integral(i, i, 1, hnew->GetNbinsY());
        {
          auto xshift = hnew->GetXaxis()->GetBinLowEdge(i);
          auto meanCh = xshift * 4.0 + strawCenterMM.at(h.first);
          auto bin = mmSingle_corr0_straw->GetXaxis()->FindBin(meanCh);
          printf("For straw %d, bin %d: integral - %g, in single_corr_straw - %g (meanch %g) -> %g ---- %g\n",
                 h.first, i, integ,
                 mmSingle_corr0_straw->GetBinContent(bin), meanCh,
                 1.0 * integ / mmSingle_corr0_straw->GetBinContent(bin),
                 1.0 * integ / mm_corr0_straw->GetBinContent(bin));
        }
        if(!integ) continue;
        for(auto j = 1; j <= hnew->GetNbinsY(); j++){
          auto c = hnew->GetBinContent(i, j);
          auto e = hnew->GetBinError(i, j);
          hnew->SetBinContent(i, j, static_cast<float>(c) / static_cast<float>(integ));
          hnew->SetBinError(i, j, static_cast<float>(e) / static_cast<float>(integ));
        }
      }
      straw_rt_normed.emplace(h.first, hnew);
    }
    for(auto &h: straw_rt_0){
      auto hnew = static_cast<TH2D*>(h.second->Clone(Form("straw%d_rt_0_normed", h.first)));
      for(auto i = 1; i <= hnew->GetNbinsX(); i++){
        auto integ = hnew->Integral(i, i, 1, hnew->GetNbinsY());
        if(!integ) continue;
        for(auto j = 1; j <= hnew->GetNbinsY(); j++){
          auto c = hnew->GetBinContent(i, j);
          auto e = hnew->GetBinError(i, j);
          hnew->SetBinContent(i, j, static_cast<float>(c) / static_cast<float>(integ));
          hnew->SetBinError(i, j, static_cast<float>(e) / static_cast<float>(integ));
        }
      }
      straw_rt_0_normed.emplace(h.first, hnew);
    }
    out->cd();

    auto straw_rt_normed_proj_dir = out->mkdir("straw_rt_normed_proj");
    straw_rt_normed_proj_dir->cd();
    for(auto &m: {straw_rt_normed, straw_rt_0_normed})
      for(auto &h: m){
        auto hist = h.second->ProjectionY(Form("%s_projectiony", h.second->GetName()));
        hist->SetTitle(Form("%s - projectionY", h.second->GetTitle()));
      }
    out->cd();

    auto straw_rt_proj_dir = out->mkdir("straw_rt_proj");
    straw_rt_proj_dir->cd();
    vector<TH1D*> straw_rt_proj;
    for(auto &m: {straw_rt, straw_rt_0})
      for(auto &h: m){
        auto hist = h.second->ProjectionY(Form("%s_projectiony", h.second->GetName()));
        hist->SetTitle(Form("%s - projectionY", h.second->GetTitle()));
      }
    out->cd();

    auto straw_rt_normed_projsum_dir = out->mkdir("straw_rt_normed_proj-sum");
    straw_rt_normed_projsum_dir->cd();
    for(auto &m: {straw_rt_normed, straw_rt_0_normed})
      for(auto &h: m){
        auto hOld = h.second;
        auto hist = new TH1D(Form("%s_projection_sum", hOld->GetName()), Form("%s: %s_projection_sum", file.Data(), hOld->GetName()),
                             hOld->GetNbinsY(), hOld->GetYaxis()->GetBinLowEdge(1), hOld->GetYaxis()->GetBinUpEdge(hOld->GetNbinsY()));
        double err;
        for(auto i = 1; i < hOld->GetNbinsY(); i++){
          auto integral = hOld->IntegralAndError(1, hOld->GetNbinsX(), i, i, err);
          hist->SetBinContent(i, integral);
          hist->SetBinError(i, err);
        }
        // hist->Scale(1.0/hist->Integral());
      }
    out->cd();

    auto straw_rt_normed_alternative_dir = out->mkdir("straw_rt_normed_alternative");
    straw_rt_normed_alternative_dir->cd();
    map<int, TH2D*> straw_rt_0_normed_alt ;
    for(auto &h: straw_rt_0){
      auto hnew = static_cast<TH2D*>(h.second->Clone(Form("straw%d_rt_0_normed_alt", h.first)));
      for(auto i = 1; i <= hnew->GetNbinsX(); i++){
        auto xshift = hnew->GetXaxis()->GetBinLowEdge(i);
        auto meanCh = xshift * 4.0 + strawCenterMM.at(h.first);
        auto bin = mmSingle_corr0->GetXaxis()->FindBin(meanCh);
        auto integ = mmSingle_corr0->GetBinContent(bin);
        if(!integ) continue;
        for(auto j = 1; j <= hnew->GetNbinsY(); j++){
          auto c = hnew->GetBinContent(i, j);
          auto e = hnew->GetBinError(i, j);
          hnew->SetBinContent(i, j, static_cast<float>(c) / static_cast<float>(integ));
          hnew->SetBinError(i, j, static_cast<float>(e) / static_cast<float>(integ));
        }
      }
      straw_rt_0_normed_alt.emplace(h.first, hnew);
    }
    out->cd();


    
    threePlotDrawF(mm_vs_sci_3det_corr, straw_vs_sci_3det_corr, straw_vs_mm_3det_corr);
    threePlotDrawF(mm_vs_sci_3det_corr_0, straw_vs_sci_3det_corr_0, straw_vs_mm_3det_corr_0, "_0");

    out->Write();
    out->Close();
}
