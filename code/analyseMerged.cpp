
#include "evBuilder.C"
#include "apv.C"

#include <optional>
#include <utility>
#include <algorithm>


using std::optional, std::nullopt;
using std::pair, std::make_pair;

constexpr unsigned int layerDoubleReadout = 2; // 0

optional<pair<double, double>> getMeanPosition(map<int, int> hitsPerLayer, int layer, bool vmmHits = false){
  long long sum = 0, sumWeights = 0, nHits = 0;
  optional<double> min = nullopt, max = nullopt;
  for(auto &h: hitsPerLayer){
    long long strip = h.first;
    long long pdo = h.second;
    if(!vmmHits && (h.first <= 118 || h.first >= 172))
      continue; // TODO check why
    // shifts from Stefano
    switch(layer){
      case 0:
        strip = h.first;
        break;
      case 1:
        strip = h.first * (1 - 2.29e-3) - 2.412 / 0.25;
        break;        
      case 2:
        strip = h.first * (1 - 8e-3) - 8.46 / 0.25;
        break;
      default:
        strip = h.first;
        break;
    }
    if(!min || min.value() < strip) min = {strip};
    if(!max || max.value() > strip) max = {strip};
    if(layer == 0)
      printf("Strip in layer 0: %lld\n", strip);
    sum += strip * pdo;
    sumWeights += pdo;
    nHits++;
  }
  optional<pair<double, double>> center = nullopt;
  if(nHits){
    double centerV = static_cast<double>(sum) / static_cast<double>(sumWeights);
    center.emplace(make_pair(centerV,
                             std::max(fabs(centerV-min.value()), fabs(centerV-max.value()))));
  }
  return center;
}

/*
 * positions:
 * L0 - L1: 285
 * L1 - L2: 345
 * L2 - Straw: 523
 * return: position in mm, from Layer 0
 */
int getLayerPosition(int layer){
  int y = 0;
  switch(layer){
    case 3:
      y += 523;
    case 2:
      y += 345;
    case 1:
      y += 285;
    case 0:
      y += 0;
  }
  return y;
}

pair<double, double> getEstimatedTrack(map<int, pair<double, double>> positions){
  int n = 0;
  double sumWXY = 0, sumWX = 0, sumWY = 0, sumWX2 = 0, sumY2 = 0, sumW = 0;
  double y, x, w, yE; // x: horizontal, y - vertical coordinate
  for(auto &pos: positions){
    x = getLayerPosition(pos.first);
    y = pos.second.first;
    yE = pos.second.second;
    w = 1.0 / yE;
    sumW += w;
    sumWXY += w * x * y;
    sumWX += w * x;
    sumWY += w * y;
    sumWX2 += w * x * x;
    // sumY2 += y * y;
    n++;
  }
  double d = sumW * sumWX2 - sumWX * sumWX;
  double a0 = (sumW*sumWXY - sumWX*sumWY) / d;
  double b0 = (sumWX2 * sumWY - sumWX*sumWXY)/d;
  double sigmaa = sqrt(sumW/d);
  double sigmab = sqrt(sumWX2/d);
  // printf("%d layers -- a0: %g, b0: %g; a1: %g, b1: %g;\n", n, a0, b0, a1, b1);
  for(auto &pos: positions){
    x = getLayerPosition(pos.first);
    y = pos.second.first;
    yE = pos.second.second;
    printf("X = %g; y = %g; y_calc = %g\t\tyE = %g\n", x, y, a0 * x + b0, yE);
  }
  return {a0, b0};
}

double estimatePositionInLayer(pair<double, double> trackAB, int layer){
  double x = getLayerPosition(layer);
  double y = trackAB.first * x + trackAB.second;
  return y;
}

TH2F* renormToUnityByY(TH2F* histIn){
  TString newName = TString(histIn->GetName()) + TString("_normY");
  auto histOut = static_cast<TH2F*>(histIn->Clone(newName));
  // histOut->SetTitle(Form("%s: microMegas vs additional straw spatial correaltion (normed);straw ch;MM ch", file.Data()));
  for(auto i = 1; i <= histOut->GetNbinsX(); i++){
    auto integ = histOut->Integral(i, i, 1, histOut->GetNbinsY());
    if(!integ) continue;
    for(auto j = 1; j <= histOut->GetNbinsY(); j++){
      auto c = histOut->GetBinContent(i, j);
      auto e = histOut->GetBinError(i, j);
      histOut->SetBinContent(i, j, static_cast<float>(c) / static_cast<float>(integ));
      histOut->SetBinError(i, j, static_cast<float>(e) / static_cast<float>(integ));
    }
  }
  return histOut;
}
TH2F* renormToUnityByX(TH2F* histIn){
  TString newName = TString(histIn->GetName()) + TString("_normX");
  auto histOut = static_cast<TH2F*>(histIn->Clone(newName));
  // histOut->SetTitle(Form("%s: microMegas vs additional straw spatial correaltion (normed);straw ch;MM ch", file.Data()));
  for(auto j = 1; j <= histOut->GetNbinsY(); j++){
    auto integ = histOut->Integral(1, histOut->GetNbinsX(), j, j);
    if(!integ) continue;
    for(auto i = 1; i <= histOut->GetNbinsX(); i++){
      auto c = histOut->GetBinContent(i, j);
      auto e = histOut->GetBinError(i, j);
      histOut->SetBinContent(i, j, static_cast<float>(c) / static_cast<float>(integ));
      histOut->SetBinError(i, j, static_cast<float>(e) / static_cast<float>(integ));
    }
  }
  return histOut;
}

// output: mm
double stripToCoord(double strip){
  return strip * 0.25;
}

bool drLayerFromVMM = true;
void analyseMerged(string runVMM = "0832", string runAPV = "423", bool tight = false, bool timefix = true, string partialfileEnd = ""){
  string tightText = tight ? "_tight" : "";
  string fixTimeText = timefix ? "_timefix" : "";
  // string mergedFileName = "../out/runMerged_run_0832_run423.root";
  // TString mergedFileName = "../out/runMerged_run_0832_run423_timefix_all-100.root";
  // TString mergedFileName = "../out/runMerged_run_0832_run423_timefix_5-10.root";
  auto mergedFileName = TString("../out/runMerged_run_"+runVMM+"_run"+runAPV+tightText+fixTimeText+partialfileEnd+".root");
  // auto f = TFile::Open(mergedFileName, "read");

  auto tVMM = new TChain("vmm");
  tVMM->AddFile(mergedFileName);
  // auto tVMM = static_cast<TTree*>(f->Get("vmm"));
  auto vmman = new evBuilder(tVMM, "g1_p25_s100-0&60", "map-20220605");
  vmman->useSyncSignal();
  auto tAPV = new TChain("apv_raw");
  tAPV->AddFile(mergedFileName);
  // auto tAPV = static_cast<TTree*>(f->Get("apv_raw"));
  auto apvan = new apv(tAPV, nullptr);
  apvan->useSyncSignal();

  auto fOut = TFile::Open(Form("_analysed_%s_%s%s.root", runVMM.c_str(), runAPV.c_str(), (tightText+fixTimeText+partialfileEnd).c_str()), "recreate");

  auto hL0Straw = new TH2F("hL0Straw", "hL0Straw; straw; L0 strip", 6, 24, 30, 206, 54, 260);
  auto hL1Straw = new TH2F("hL1Straw", "hL2Straw; straw; L1 strip", 6, 24, 30, 206, 54, 260);
  auto hL2Straw = new TH2F("hL2Straw", "hL2Straw; straw; L2 strip", 6, 24, 30, 206, 54, 260);
  auto estimatedStraw = new TH2F("estimatedStraw", "Estimated track position; straw; estimated position, mm", 6, 24, 30, 500, 0, 50);
  // // 1 in MM strips
  // auto estimatedStraw1 = new TH2F("estimatedStraw1", "Estimated track position; straw; estimated position, mm", 6, 24, 30, 206, 54, 260);

  map<pair<int,int>, TH2F*> dtHists;
  for(auto &straw: {25, 26, 27, 28, 29}){
    for(auto &l: {0, 1, 2}){
      dtHists.emplace(make_pair(l, straw), new TH2F(Form("hL%dStraw%dDT",l, straw), Form("hL%dStraw%dDT; MM strip; dT, ns",l, straw), 206, 54, 260, 200, -100, 100));
    }
  }
  map<int, TH2F*> dtHistsProj, dtHistsProj1;
  for(auto &straw: {25, 26, 27, 28, 29}){
    dtHistsProj.emplace(straw, new TH2F(Form("hProjecteddStraw%dDT",straw), Form("hProjecteddStraw%dDT; position, mm; dT, ns",straw), 500, 0, 50, 200, -100, 100));
    dtHistsProj1.emplace(straw, new TH2F(Form("hProjecteddStraw%dDT_1",straw), Form("hProjecteddStraw%dDT_1; position, strip; dT, ns",straw), 206, 54, 260, 200, -100, 100));
  }

  int nGood = 0;
  auto nEvents = apvan->GetEntries();
  for(auto event = 0; event < nEvents; event++){
    auto dataAPV = apvan->GetCentralHits2ROnlyData(event);
    auto dataVMM = vmman->getHits(event);
    map<int, pair<double,double>> positions;
    if(!getMeanPosition(dataAPV.hitsPerLayer.at(2), 2))
      continue;
    for(auto i = 0; i < 3; i++){
      auto layerData = dataAPV.hitsPerLayer.at(i);
      auto meanPos = getMeanPosition(layerData, i);
      if(!meanPos)
        continue;
      if(i < 2 || !drLayerFromVMM)
        positions.emplace(i, meanPos.value());
      // printf("Event %d, position for layer %d: %g\n", event, i, meanPos.value());
      TH2F* hist = nullptr;
      switch(i){
        case 0:
          hist = hL0Straw;
          break;
        case 1:
          hist = hL1Straw;
          break;
        case 2:
          hist = hL2Straw;
          break;
      }
      for(auto h: dataVMM){
        if(h.detector != 1)
          continue;
        if(hist != nullptr){
          for(auto &hAPV: layerData)
            hist->Fill(h.strip, hAPV.first);
        }
        for(auto &hAPV: layerData){
          // printf("Filling %d, %d\n",i, h.strip);
          dtHists.at(make_pair(i, h.strip))->Fill(hAPV.first, h.timeToScint);
        }
      }
    }
    if(positions.size() < 2)
      continue;
    // printf("Event %d\t", event);

    if(drLayerFromVMM){
      // Adding MM layer from VMM data
      map<int, int> hitsMMVMM;
      for(auto h: dataVMM){
        if(h.detector != 4) continue;
        // printf("MM VMM strip %d\n", h.strip);
        hitsMMVMM.emplace(h.strip, h.pdo);
      }
      // printf("Number of MM VMM strips %lu\n", hitsMMVMM.size());
      auto meanPosMMVMM = getMeanPosition(hitsMMVMM, 2, true);
      if(!meanPosMMVMM)
        continue;
      positions.emplace(2, meanPosMMVMM.value());
    }

    auto trackParam = getEstimatedTrack(positions);
    auto estimatedCoord = estimatePositionInLayer(trackParam, 3);
    printf("Event %d -- estimated: %g\n", event, estimatedCoord);
    for(auto h: dataVMM){
      if(h.detector != 1)
        continue;
      estimatedStraw->Fill(h.strip, stripToCoord(estimatedCoord));
      dtHistsProj.at(h.strip)->Fill(stripToCoord(estimatedCoord), h.timeToScint);
      dtHistsProj1.at(h.strip)->Fill(estimatedCoord, h.timeToScint);
    }
    nGood++;
  }

  for(auto &h: {hL0Straw, hL1Straw, hL2Straw}){
    renormToUnityByY(h);
    renormToUnityByX(h);
  }

  fOut->Write();
  fOut->Close();
  
  printf("nEvents: %u - %u | %d\n", apvan->GetEntries(), vmman->GetEntries(), nGood);
}
