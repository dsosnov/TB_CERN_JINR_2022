#pragma once

template<class T>
void normByX(T *hIn, T *hOut){
  for(auto i = 1; i <= hOut->GetNbinsX(); i++){
    auto integ = hOut->Integral(i, i, 1, hOut->GetNbinsY());
    if(!integ) continue;
    for(auto j = 1; j <= hOut->GetNbinsY(); j++){
      auto c = hOut->GetBinContent(i, j);
      auto e = hOut->GetBinError(i, j);
      hOut->SetBinContent(i, j, static_cast<float>(c) / static_cast<float>(integ));
      hOut->SetBinError(i, j, static_cast<float>(e) / static_cast<float>(integ));
    }
  }
}

template<class T>
T* normByX(T *hIn){
  auto hOut = static_cast<T*>(hIn->Clone(Form("%s_norm", hIn->GetName())));
  normByX(hIn, hOut);
  return hOut;
}

template<class T>
void normByY(T *hIn, T *hOut){
  for(auto i = 1; i <= hOut->GetNbinsY(); i++){
    auto integ = hOut->Integral(1, hOut->GetNbinsY(), i, i);
    if(!integ) continue;
    for(auto j = 1; j <= hOut->GetNbinsX(); j++){
      auto c = hOut->GetBinContent(j, i);
      auto e = hOut->GetBinError(j, i);
      hOut->SetBinContent(j, i, static_cast<float>(c) / static_cast<float>(integ));
      hOut->SetBinError(j, i, static_cast<float>(e) / static_cast<float>(integ));
    }
  }
}

template<class T>
T* normByY(T *hIn){
  auto hOut = static_cast<T*>(hIn->Clone(Form("%s_norm", hIn->GetName())));
  normByY(hIn, hOut);
  return hOut;
}
