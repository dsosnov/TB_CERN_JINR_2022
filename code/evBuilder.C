#define evBuilder_cxx
#include "evBuilder.h"
#include <TH1.h>


void evBuilder::threePlotDrawF(TH1D *h1, TH1D *h2, TH1D *h3, TString fileEnding)
{
  if(!(h1->GetEntries()) || !(h2->GetEntries()) || !(h3->GetEntries()))
    return;
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

vector<analysisGeneral::mm2CenterHitParameters> evBuilder::GetCentralHits(unsigned long long fromSec, unsigned long long toSec){
  printf("evBuilder::GetCentralHits(%llu, %llu)\n", fromSec, toSec);

  vector<analysisGeneral::mm2CenterHitParameters> outputData = {};
  unsigned long long hitsToPrev = 0;

  if (fChain == 0)
    return outputData;

  unsigned int pdoThr = 100;
  unsigned int nLoopEntriesAround = 0;
  Long64_t nentries = fChain->GetEntries();

  Long64_t nbytes = 0, nb = 0, mbytes = 0, mb = 0;

  // =============================== CORRELATION FINDING ===============================
  for (Long64_t jentry = 0; jentry < nentries; jentry++) // You can remove "/ 10" and use the whole dataset
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
      break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;


    if(fromSec > 0 && daq_timestamp_s->at(0) < fromSec)
      continue;
    else if(toSec > 0 && toSec >= fromSec && daq_timestamp_s->at(0) > toSec)
      break;

    bool isTrigger = false;
    int trigTime = 0;
    for (int j = 0; j < channel->at(0).size(); j++)
    {
      int fch = channel->at(0).at(j);
      int fchD = getMappedDetector(fch);
      int fchM = getMappedChannel(fch);
      // printf("VERSO: ch %d -> %d %d\n", fch, fchD, fchM);
      
      int fpdoUC = pdo->at(0).at(j); // Uncorrected PDO, used at time calibration
      int fpdo = correctPDO(fch, fpdoUC);
      int ftdo = tdo->at(0).at(j);
      int fbcid = grayDecoded->at(0).at(j);

      // isTrigger = isTrigger || (ffchD == 0 && ffchM == 3);
      if(fchD == 0 && fchM == 3){
        isTrigger = true;
        trigTime = getTime(fch, fbcid, ftdo, fpdoUC);
        break;
      }
    }
    if(!isTrigger){
      hitsToPrev++;
      continue;
    }

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
                      
        int ffpdoUC = pdo->at(0).at(k); // Uncorrected PDO, used at time calibration
        int ffpdo = correctPDO(ffch, ffpdoUC);
        int fftdo = tdo->at(0).at(k);
        int ffbcid = grayDecoded->at(0).at(k);
        // double fft = getTimeByHand(ffbcid, fftdo, 110, 160); //'hand' limits
        double fft = getTime(ffch, ffbcid, fftdo, ffpdoUC); // 'auto' limits

        if(ffpdo < pdoThr) continue;

        if(ffchD == 4) // All MM channels
        { 
          if (fabs(trigTime - fft) < 500)
          {
            MmCluster.push_back({ffchM * 1.0, ffpdo * 1.0, fft});
          }
        }
      }
    }

    double minT_straw_mm = 600;
    auto [meanT, meanCh, maxPdo] = getClusterParameters(trigTime, minT_straw_mm);

    if(meanCh == 0){
        hitsToPrev++;
        continue;
    }

    analysisGeneral::mm2CenterHitParameters hit;
    hit.timeSec = daq_timestamp_s->at(0);
    hit.timeMSec = daq_timestamp_ns->at(0) / 1000;
    hit.strip = meanCh;
    hit.pdo = static_cast<int>(maxPdo);
    hit.pdoRelative = maxPdo / 1024.0;
    hit.nHitsToPrev = hitsToPrev;
    hit.time = trigTime - meanT;
    hitsToPrev = 0;

    outputData.push_back(hit);
  }
  
  // for (auto event = 0; event < nentries; event++){
  //   Long64_t ientry = LoadTree(event);
  //   if (ientry < 0) break;
  //   GetEntry(event);
  //   // printf("Entry: %d: error - %d\n", event, error);
  //   if(error)
  //     continue;
      
  //   if(fromSec > 0 && daqTimeSec < fromSec)
  //     continue;
  //   if(toSec > 0 && toSec >= fromSec && daqTimeSec > toSec)
  //     break;

  //   /* Checking if there is any mapped channels */
  //   bool notErr = false;
  //   for (auto &r: *mmReadout)
  //     notErr = notErr || (r != 'E');
  //   if(!notErr) continue;
    
  //   unsigned long long currentTimestamp = daqTimeSec * 1E6 + daqTimeMicroSec;
  //   previousTimestamp = currentTimestamp;

  //   hits.clear();
  //   for (int j = 0; j < max_q->size(); j++){
  //     // printf("Record inside entry: %d\n", j);
  //     auto readout = mmReadout->at(j);
  //     if(readout == 'E') //non-mapped channel
  //       continue;
  //     auto layer = mmLayer->at(j);
  //     auto strip = mmStrip->at(j);
  //     auto maxQ = max_q->at(j);
      
  //     auto maxTime = t_max_q->at(j);
      
  //     hits.push_back({layer, strip, maxQ, maxTime, raw_q->at(j)});      
  //   }

  //   /* Constructing clasters */
  //   constructClasters();
    
  //   auto tracks = constructTracks(clasters);
  //   bool trackIn2Center = false;
  //   apvTrack *highestTrack = nullptr;
  //   for(auto &t: tracks){
  //     auto x2 = get<2>(getHitsForTrack(t));
  //     if(x2 > 153 && x2 < 210){
  //       trackIn2Center = true;
  //       if(highestTrack == nullptr || highestTrack->maxQ() < t.maxQ())
  //         highestTrack = &t;
  //     }
  //   }
  //   if(!trackIn2Center){
  //     hitsToPrev++;
  //     continue;
  //   }

  //   analysisGeneral::mm2CenterHitParameters hit;
  //   hit.timeSec = daqTimeSec;
  //   hit.timeMSec = daqTimeMicroSec;
  //   hit.strip = get<2>(getHitsForTrack(*highestTrack));
  //   hit.pdo = highestTrack->maxQ();
  //   hit.pdoRelative = static_cast<double>(hit.pdo) / 2048;
  //   hit.nHitsToPrev = hitsToPrev;
  //   hitsToPrev = 0;

  //   outputData.push_back(hit);
  // }
  return outputData;
}


void evBuilder::LoopSecond(unsigned long long sec){
  printf("evBuilder::LoopSecond(%llu)\n", sec);

  if (fChain == 0)
        return;

  unsigned int pdoThr = 100;
  unsigned int nLoopEntriesAround = 0;
  Long64_t nentries = fChain->GetEntries();

  Long64_t nbytes = 0, nb = 0, mbytes = 0, mb = 0;

  // =============================== CORRELATION FINDING ===============================
  for (Long64_t jentry = 0; jentry < nentries; jentry++) // You can remove "/ 10" and use the whole dataset
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
      break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;


    if(daq_timestamp_s->at(0) < sec)
      continue;
    else if(daq_timestamp_s->at(0) > sec)
      break;

    bool isTrigger = false;
    int trigTime = 0;
    for (int j = 0; j < channel->at(0).size(); j++)
    {
      int fch = channel->at(0).at(j);
      int fchD = getMappedDetector(fch);
      int fchM = getMappedChannel(fch);
      // printf("VERSO: ch %d -> %d %d\n", fch, fchD, fchM);
      
      int fpdoUC = pdo->at(0).at(j); // Uncorrected PDO, used at time calibration
      int fpdo = correctPDO(fch, fpdoUC);
      int ftdo = tdo->at(0).at(j);
      int fbcid = grayDecoded->at(0).at(j);

      // isTrigger = isTrigger || (ffchD == 0 && ffchM == 3);
      if(fchD == 0 && fchM == 3){
        isTrigger = true;
        trigTime = getTime(fch, fbcid, ftdo, fpdoUC);
        break;
      }
    }
    if(!isTrigger)
      continue;

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
                      
        int ffpdoUC = pdo->at(0).at(k); // Uncorrected PDO, used at time calibration
        int ffpdo = correctPDO(ffch, ffpdoUC);
        int fftdo = tdo->at(0).at(k);
        int ffbcid = grayDecoded->at(0).at(k);
        // double fft = getTimeByHand(ffbcid, fftdo, 110, 160); //'hand' limits
        double fft = getTime(ffch, ffbcid, fftdo, ffpdoUC); // 'auto' limits

        if(ffpdo < pdoThr) continue;

        if(ffchD == 4) // All MM channels
        { 
          if (fabs(trigTime - fft) < 500)
          {
            MmCluster.push_back({ffchM * 1.0, ffpdo * 1.0, fft});
          }
        }
      }
    }

    double minT_straw_mm = 600;
    auto [meanT, meanCh, maxPdo] = getClusterParameters(trigTime, minT_straw_mm);

    if(meanCh != 0){
      printf("VERSO: time %d - %d, Event ID %d:\n", daq_timestamp_s->at(0), daq_timestamp_ns->at(0), triggerCounter->at(0));
      printf("  MM hits to %g (time: meanT %.2g)\n", meanCh, meanT);
    }
  }
}

tuple<double, double, double> evBuilder::getClusterParameters(double t_srtraw, double minT_straw_mm, int workType){
  double meanT = 0;
  double meanCh = 0;
  double sum1 = 0;
  double sum2 = 0;
  double w_sum = 0;
  double maxPdo = 0;

  // double minT_straw_mm = 600;
  int mmCh_min = 0;

  double maximalPdoHit = -1;
  
  if (MmCluster.size() != 0)
  {
    switch(workType){
      case 0:
        for (size_t l = 0; l < MmCluster.size(); l++)
        {
          if (abs(MmCluster.at(l).time - t_srtraw) < minT_straw_mm)
          {
            minT_straw_mm = abs(MmCluster.at(l).time - t_srtraw);
            mmCh_min = MmCluster.at(l).channel;
          }
        }
        for (size_t l = 0; l < MmCluster.size(); l++)
        {
          if (abs(MmCluster.at(l).channel - mmCh_min) > 5)
          {
            MmCluster.erase(MmCluster.begin()+l);
          }
        }
        for (size_t l = 0; l < MmCluster.size(); l++)
        {
          sum1 += MmCluster.at(l).time * MmCluster.at(l).pdo / 1024.0;
          sum2 += MmCluster.at(l).channel * MmCluster.at(l).pdo / 1024.0;
          w_sum += MmCluster.at(l).pdo / 1024.0;
          if(MmCluster.at(l).pdo > maxPdo)
            maxPdo = MmCluster.at(l).pdo;
        }
        meanT = sum1 / w_sum;
        meanCh = sum2 / w_sum;
        break;
      case 1:
        for (size_t l = 0; l < MmCluster.size(); l++)
          if(maximalPdoHit < 0 || MmCluster.at(l).pdo > MmCluster.at(maximalPdoHit).pdo)
            maximalPdoHit = l;
        meanT = MmCluster.at(maximalPdoHit).time;
        meanCh = MmCluster.at(maximalPdoHit).channel;
        maxPdo = MmCluster.at(maximalPdoHit).pdo;
        break;
    }
  }
  return {meanT, meanCh, maxPdo};
}

void evBuilder::Loop()
{
  printf("evBuilder::Loop()\n");

  if (fChain == 0)
        return;

   int strawMin = -1, strawMax = -1, mmMin = -1, mmMax = -1;
   int addstrawMin = -1, addstrawMax = -1;
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
     } else if(s.second.first == 6){
       if(addstrawMin < 0 || addstrawMin > s.second.second)
         addstrawMin = s.second.second;
       if(addstrawMax < 0 || addstrawMax < s.second.second)
         addstrawMax = s.second.second;
     }
   }
    printf("Straws: %d-%d, MM: %d-%d, additional straws: %d-%d\n", strawMin, strawMax, mmMin, mmMax, addstrawMin, addstrawMax);

    TFile *out = new TFile("../out/out_" + file + ending, "RECREATE"); // PATH where to save out_*.root file

    auto straw_vs_sci = new TH1D("straw_vs_sci", Form("%s: straw vs scint;#Deltat, ns", file.Data()), 1000, -500, 500);
    auto straw_vs_mm = new TH1D("straw_vs_mm", Form("%s: straw vs microMegas;#Deltat, ns", file.Data()), 1000, -500, 500);
    auto straw_vs_mm_cluster = new TH1D("straw_vs_mm_cluster", Form("%s: straw vs microMegas cluster time;#Deltat, ns", file.Data()), 1000, -500, 500);
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
    auto straw_add_vs_mm_spatial_corr = new TH2D("straw_add_vs_mm_spatial_corr", Form("%s: microMegas vs additional straw spatial correaltion;straw ch;MM ch", file.Data()),
                                             addstrawMax - addstrawMin + 1, addstrawMin, addstrawMax + 1, mmMax - mmMin + 1, mmMin, mmMax);
    for(auto i = addstrawMin; i <= addstrawMax; i++)
      if(addStrawType.count(i))
        straw_add_vs_mm_spatial_corr->GetXaxis()->SetBinLabel(i+1, addStrawType.at(i).c_str());

    auto straw_rt_dir = out->mkdir("straw_rt");
    straw_rt_dir->cd();
    map<int, TH2D*> straw_rt, straw_rt_0, straw_rt_add, straw_rt_add_0 ;
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
    if(addstrawMin >= 0 && addstrawMax >= 0){
      for(auto i = addstrawMin; i <= addstrawMax; i++){
        straw_rt_add.emplace(i, new TH2D(Form("straw_%s_rt", addStrawType.at(i).c_str()),
                                         Form("%s: %s straw v-shape sci ch 60;R, mm;T, ns", file.Data(), addStrawType.at(i).c_str()),
                                         80, -10, 10, 1000, -100, 900));
        straw_rt_add_0.emplace(i, new TH2D(Form("straw_%s_rt_0", addStrawType.at(i).c_str()),
                                           Form("%s: %s straw v-shape sci ch 0;R, mm;T, ns", file.Data(), addStrawType.at(i).c_str()),
                                           80, -10, 10, 1000, -100, 900));
      }
    }
    out->cd();

    auto straw_pdo_dir = out->mkdir("straw_pdo_corr");
    straw_pdo_dir->cd();
    map<int, TH1D*> straw_pdo, straw_pdo_0, straw_pdo_add, straw_pdo_add_0 ;
    for(auto i = strawMin; i <= strawMax; i++){
      straw_pdo.emplace(i,
                        new TH1D(Form("straw%d_pdo_corr", i),
                                 Form("%s: pdo for straw %d corellated with sci ch 60;pdo", file.Data(), i), 64, 0, 1024));
      straw_pdo_0.emplace(i,
                          new TH1D(Form("straw%d_pdo_corr_0", i),
                                   Form("%s: pdo for additional straw %d corellated with sci ch 0;pdo", file.Data(), i), 64, 0, 1024));
    }
    if(addstrawMin >= 0 && addstrawMax >= 0){
      for(auto i = addstrawMin; i <= addstrawMax; i++){
        straw_pdo_add.emplace(i,
                              new TH1D(Form("straw_%s_pdo_corr", addStrawType.at(i).c_str()),
                                       Form("%s: pdo for %s straw corellated with sci ch 60;pdo", file.Data(), addStrawType.at(i).c_str()), 64, 0, 1024));
        straw_pdo_add_0.emplace(i,
                                new TH1D(Form("straw_%s_pdo_corr_0", addStrawType.at(i).c_str()),
                                         Form("%s: pdo for %s straw corellated with sci ch 0;pdo", file.Data(), addStrawType.at(i).c_str()), 64, 0, 1024));
      }
    }
    out->cd();

    auto straw_deltat_dir = out->mkdir("straw_deltat_corr");
    straw_deltat_dir->cd();
    map<int, TH1D*> straw_deltat, straw_deltat_0, straw_deltat_add, straw_deltat_add_0 ;
    for(auto i = strawMin; i <= strawMax; i++){
      straw_deltat.emplace(i,
                        new TH1D(Form("straw%d_vs_sci60", i),
                                 Form("%s: straw%d_vs_sci60;#Delta t", file.Data(), i), 1000, -500, 500));
      straw_deltat_0.emplace(i,
                        new TH1D(Form("straw%d_vs_sci0", i),
                                 Form("%s: straw%d_vs_sci0;#Delta t", file.Data(), i), 1000, -500, 500));
    }
    if(addstrawMin >= 0 && addstrawMax >= 0){
      for(auto i = addstrawMin; i <= addstrawMax; i++){
        straw_deltat_add.emplace(i,
                                 new TH1D(Form("straw_%s_vs_sci60", addStrawType.at(i).c_str()),
                                          Form("%s: %s straw_vs_sci60;#Delta t", file.Data(), addStrawType.at(i).c_str()), 1500, -500, 1000));
        straw_deltat_add_0.emplace(i,
                                   new TH1D(Form("straw_%s_vs_sci0", addStrawType.at(i).c_str()),
                                            Form("%s: %s straw_vs_sci0;#Delta t", file.Data(), addStrawType.at(i).c_str()), 1500, -500, 1000));
      }
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

    unsigned int pdoThr = 100;
    unsigned int nLoopEntriesAround = 0;
    Long64_t nentries = fChain->GetEntries();

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

            int fpdoUC = pdo->at(0).at(j); // Uncorrected PDO, used at time calibration
            int fpdo = correctPDO(fch, fpdoUC);
            int ftdo = tdo->at(0).at(j);
            int fbcid = grayDecoded->at(0).at(j);

            if(fpdo < pdoThr) continue;

            if (fchD == 1) // All straw ch
            {
                straw_bcid_ch_srtraw = fbcid;
                straw_pdo_ch_srtraw = fpdo;
                t_srtraw = getTime(fch, fbcid, ftdo, fpdoUC); // 'auto' limits
                // t_srtraw = getTimeByHand(fbcid, ftdo, Y, Y); //'hand' limits

                straw_pdo_ucl_cr.at(fchM)->Fill(fpdo);
                straw_pdo_ucl_ucr.at(fchM)->Fill(fpdoUC);

                Long64_t mbytes = 0, mb = 0;
                double t30 = 0;
                double minTsci0 = 1E3;
                double sciT_ch0 = 0;
                double minTsci60 = 1E3;
                double sciT_ch60 = 0;
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
                      if(ffch == fch) continue;
                      int ffchD = getMappedDetector(ffch);
                      int ffchM = getMappedChannel(ffch);
                      
                      int ffpdoUC = pdo->at(0).at(k); // Uncorrected PDO, used at time calibration
                      int ffpdo = correctPDO(ffch, ffpdoUC);
                      int fftdo = tdo->at(0).at(k);
                      int ffbcid = grayDecoded->at(0).at(k);
                      // double fft = getTimeByHand(ffbcid, fftdo, 110, 160); //'hand' limits
                      double fft = getTime(ffch, ffbcid, fftdo, ffpdoUC); // 'auto' limits

                      if(ffpdo < pdoThr) continue;

                      if(ffchD == 4) // All MM channels
                      { 
                        if (fabs(t_srtraw - fft) < 500)
                        {
                          straw_vs_mm ->Fill(t_srtraw - fft);
                          MmCluster.push_back({ffchM * 1.0, ffpdo * 1.0, fft});
                        }
                      }
                      else if (ffchD == 0 && ffchM == 0) // Sci 0
                      {
                        fft = getTimeByHand(ffbcid, fftdo, 88, 140); //'hand' limits
                        if (fabs(t_srtraw - fft) < minTsci0)
                        {
                          minTsci0 = fabs(t_srtraw - fft);
                          sciT_ch0 = fft;
                        }
                      }
                      else if (ffchD == 0 && ffchM == 3) // triple sci coinsidence
                      {
                        if (fabs(t_srtraw - fft) < minTsci60)
                        {
                          minTsci60 = fabs(t_srtraw - fft);
                          sciT_ch60 = fft;
                        }
                      }
                      else if (ffchD == 1 && ffchM == fchM + 1) // Next straw channel
                      {
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

                double minT_straw_mm = 600;
                auto [meanT, meanCh, maxPdo] = getClusterParameters(t_srtraw, minT_straw_mm);
                if (MmCluster.size() != 0){
                    straw_vs_mm_spatial_corr->Fill(fchM, meanCh);
                    straw_vs_mm_cluster ->Fill(t_srtraw - meanT);
                    hits_in_cluster->Fill( MmCluster.size());
                }

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
                  straw_rt_0.at(fchM)->Fill((meanCh - strawCenterMM.at({fchD, fchM})) * 0.25, 100 + t_srtraw - sciT_ch0);
                    straw_vs_sci_3det_corr_0->Fill(t_srtraw - sciT_ch0);
                    straw_vs_mm_3det_corr_0->Fill(t_srtraw - meanT);
                    mm_vs_sci_3det_corr_0->Fill(meanT - sciT_ch0);
                }
                if (sciT_ch60 != 0 && meanT != 0)
                {
                    straw_vs_sci_3det_corr->Fill(t_srtraw - sciT_ch60);
                    straw_vs_mm_3det_corr->Fill(t_srtraw - meanT);
                    mm_vs_sci_3det_corr->Fill(meanT - sciT_ch60);
                    straw_rt.at(fchM)->Fill((meanCh - strawCenterMM.at({fchD, fchM})) * 0.25, 100 + t_srtraw - sciT_ch60);
                }
                if(neighborStrawTime != 0)
                {
                  straw_straw.at(fchM)->Fill(t_srtraw - neighborStrawTime);
                  if (sciT_ch0 != 0)
                    straw_banana_0.at(fchM)->Fill(t_srtraw - sciT_ch0, neighborStrawTime - sciT_ch0);
                  if (sciT_ch60 != 0)
                    straw_banana.at(fchM)->Fill(t_srtraw - sciT_ch60, neighborStrawTime - sciT_ch60);
                }

                // ============================= end of sci 0 correlation finding =============================

            }
            else if (fchD == 6) // Additional straws
            {
                straw_bcid_ch_srtraw = fbcid;
                straw_pdo_ch_srtraw = fpdo;
                t_srtraw = getTime(fch, fbcid, ftdo, fpdoUC); // 'auto' limits
                // t_srtraw = getTimeByHand(fbcid, ftdo, Y, Y); //'hand' limits

                                // TODO straw_pdo_ucl_cr.at(fchM)->Fill(fpdo);
                                // TODO straw_pdo_ucl_ucr.at(fchM)->Fill(fpdoUC);

                Long64_t mbytes = 0, mb = 0;
                double t30 = 0;
                double minTsci0 = 1E3;
                double sciT_ch0 = 0;
                double minTsci60 = 1E3;
                double sciT_ch60 = 0;
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
                      if(ffch == fch) continue;
                      int ffchD = getMappedDetector(ffch);
                      int ffchM = getMappedChannel(ffch);
                      
                      int ffpdoUC = pdo->at(0).at(k); // Uncorrected PDO, used at time calibration
                      int ffpdo = correctPDO(ffch, ffpdoUC);
                      int fftdo = tdo->at(0).at(k);
                      int ffbcid = grayDecoded->at(0).at(k);
                      // double fft = getTimeByHand(ffbcid, fftdo, 110, 160); //'hand' limits
                      double fft = getTime(ffch, ffbcid, fftdo, ffpdoUC); // 'auto' limits

                      if(ffpdo < pdoThr) continue;

                      if(ffchD == 4) // All MM channels
                      { 
                        if (fabs(t_srtraw - fft) < 500)
                        {
                                // TODO straw_vs_mm ->Fill(t_srtraw - fft); // 
                          MmCluster.push_back({ffchM * 1.0, ffpdo * 1.0, fft});
                        }
                      }
                      else if (ffchD == 0 && ffchM == 0) // Sci 0
                      {
                        fft = getTimeByHand(ffbcid, fftdo, 88, 140); //'hand' limits
                        if (fabs(t_srtraw - fft) < minTsci0)
                        {
                          minTsci0 = fabs(t_srtraw - fft);
                          sciT_ch0 = fft;
                        }
                      }
                      else if (ffchD == 0 && ffchM == 3) // triple sci coinsidence
                      {
                        if (fabs(t_srtraw - fft) < minTsci60)
                        {
                          minTsci60 = fabs(t_srtraw - fft);
                          sciT_ch60 = fft;
                        }
                      }
                      else if (ffchD == 1 && ffchM == fchM + 1) // Next straw channel
                      {
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

                double minT_straw_mm = 600;
                auto [meanT, meanCh, maxPdo] = getClusterParameters(t_srtraw, minT_straw_mm);
                if (MmCluster.size() != 0){
                    straw_add_vs_mm_spatial_corr->Fill(fchM, meanCh);
                }

                if (sciT_ch0 != 0)
                {
                  straw_deltat_add_0.at(fchM)->Fill(t_srtraw - sciT_ch0);
                  straw_pdo_add_0.at(fchM)->Fill(fpdo);
                }
                if (sciT_ch60 != 0)
                {
                    straw_deltat_add.at(fchM)->Fill(t_srtraw - sciT_ch60);
                    straw_pdo_add.at(fchM)->Fill(fpdo);
                }
                if (sciT_ch0 != 0 && meanT != 0)
                {
                  straw_rt_add_0.at(fchM)->Fill((meanCh - strawCenterMM.at({fchD, fchM})) * 0.25, 100 + t_srtraw - sciT_ch0);
                }
                if (sciT_ch60 != 0 && meanT != 0)
                {
                    straw_rt_add.at(fchM)->Fill((meanCh - strawCenterMM.at({fchD, fchM})) * 0.25, 100 + t_srtraw - sciT_ch60);
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
    map<int, TH2D*> straw_rt_add_normed, straw_rt_add_0_normed ;
    for(auto &h: straw_rt_add){
      auto hnew = static_cast<TH2D*>(h.second->Clone(Form("%s_normed", h.second->GetName())));
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
      straw_rt_add_normed.emplace(h.first, hnew);
    }
    for(auto &h: straw_rt_add_0){
      auto hnew = static_cast<TH2D*>(h.second->Clone(Form("%s_normed", h.second->GetName())));
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
      straw_rt_add_0_normed.emplace(h.first, hnew);
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

    auto straw_vs_mm_spatial_corr_normed = static_cast<TH2D*>(straw_vs_mm_spatial_corr->Clone("straw_vs_mm_spatial_corr_normed"));
    straw_vs_mm_spatial_corr_normed->SetTitle(Form("%s: microMegas vs straw spatial correaltion (normed);straw ch;MM ch", file.Data()));
    for(auto i = 1; i <= straw_vs_mm_spatial_corr_normed->GetNbinsX(); i++){
      auto integ = straw_vs_mm_spatial_corr_normed->Integral(i, i, 1, straw_vs_mm_spatial_corr_normed->GetNbinsY());
      if(!integ) continue;
      for(auto j = 1; j <= straw_vs_mm_spatial_corr_normed->GetNbinsY(); j++){
        auto c = straw_vs_mm_spatial_corr_normed->GetBinContent(i, j);
        auto e = straw_vs_mm_spatial_corr_normed->GetBinError(i, j);
        straw_vs_mm_spatial_corr_normed->SetBinContent(i, j, static_cast<float>(c) / static_cast<float>(integ));
        straw_vs_mm_spatial_corr_normed->SetBinError(i, j, static_cast<float>(e) / static_cast<float>(integ));
      }
    }
    
    auto straw_add_vs_mm_spatial_corr_normed = static_cast<TH2D*>(straw_add_vs_mm_spatial_corr->Clone("straw_add_vs_mm_spatial_corr_normed"));
    straw_add_vs_mm_spatial_corr_normed->SetTitle(Form("%s: microMegas vs additional straw spatial correaltion (normed);straw ch;MM ch", file.Data()));
    for(auto i = 1; i <= straw_add_vs_mm_spatial_corr_normed->GetNbinsX(); i++){
      auto integ = straw_add_vs_mm_spatial_corr_normed->Integral(i, i, 1, straw_add_vs_mm_spatial_corr_normed->GetNbinsY());
      if(!integ) continue;
      for(auto j = 1; j <= straw_add_vs_mm_spatial_corr_normed->GetNbinsY(); j++){
        auto c = straw_add_vs_mm_spatial_corr_normed->GetBinContent(i, j);
        auto e = straw_add_vs_mm_spatial_corr_normed->GetBinError(i, j);
        straw_add_vs_mm_spatial_corr_normed->SetBinContent(i, j, static_cast<float>(c) / static_cast<float>(integ));
        straw_add_vs_mm_spatial_corr_normed->SetBinError(i, j, static_cast<float>(e) / static_cast<float>(integ));
      }
    }

    threePlotDrawF(mm_vs_sci_3det_corr, straw_vs_sci_3det_corr, straw_vs_mm_3det_corr);
    threePlotDrawF(mm_vs_sci_3det_corr_0, straw_vs_sci_3det_corr_0, straw_vs_mm_3det_corr_0, "_0");

    {
      for(auto &h: straw_banana_0){
        auto banana = new TCanvas(Form("banana_%d", h.first), Form("banana_%d", h.first), 1000, 900);
        banana->cd();
        gStyle->SetOptStat(0);
        h.second->Draw("COLZ");
        banana->SaveAs(Form("../out/banana_%s_straws%d-%d.pdf", file.Data(), h.first, h.first+1));
        banana->SaveAs(Form("../out/banana_%s_straws%d-%d.png", file.Data(), h.first, h.first+1));
      }
    }

    out->Write();
    out->Close();

}
