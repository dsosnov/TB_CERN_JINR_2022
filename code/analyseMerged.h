#pragma once

#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"

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

