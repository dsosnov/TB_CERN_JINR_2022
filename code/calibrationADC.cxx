#include "calibrationADC.h"

int calibrationADC(){
  // string comment = "g1_p25_s100";
  // string callibpath = "calibration/t@t_g1_p25_s100";
  // vector<string> calibrationADC = {"run_0095.root", "run_0096.root", "run_0097.root", "run_0098.root", "run_0099.root"};

  string comment = "t@t_g3_p25_s100";
  string callibpath = "calibration/t@t_g3_p25_s100";
  vector<string> calibrationADC = {"run_0000.root", "run_0100.root", "run_0101.root", "run_0102.root", "run_0103.root", "run_0104.root"};

  auto outfile = TFile::Open(Form("calibration_pdo_%s.root", comment.c_str()), "recreate");

  vector<uint> channelsExclude = {/*0, 1, 2*/}; //0
  vector<uint> channelsDelete = {}; //35, 63
  uint nChannels = 65;

  map<string,shared_ptr<TChain>> calibrationADCf;
  for(auto &i: calibrationADC){
    calibrationADCf.emplace(i, make_shared<TChain>("vmm"));
    calibrationADCf.at(i)->Add((callibpath+"/"+i).c_str());
  }

  vector<singleTest> vectorData;
  for(auto &fileChain: calibrationADCf){
    auto calData = findCalibrationForChannels(fileChain.second, fileChain.first, nChannels, channelsExclude, channelsDelete, 120);
    vectorData.push_back(calData);
  }

  auto pdoCorrections = fitPDO(vectorData);

  for(auto &e: channelsExclude)
    pdoCorrections.push_back({int(e), 0, 0, 1, 0, 0, 0});
  std::sort(pdoCorrections.begin(), pdoCorrections.end(), [](fitPDOChannelResult r1, fitPDOChannelResult r2){return r1.nChannel < r2.nChannel;});

  auto f = fopen(Form("calibration_pdo_%s.txt", comment.c_str()), "w");
  for(auto &correction: pdoCorrections){
    correction.print();
    correction.print(true, f);
  }
  fclose(f);
  
  for(auto &fileChain: calibrationADCf){
    applyPDOCorrections(fileChain.second, fileChain.first, pdoCorrections);
  }

  outfile->Write();
  outfile->Close();
  
  return 0;
}
