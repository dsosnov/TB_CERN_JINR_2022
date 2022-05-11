#include "calibrationPDO.h"

int calibrationPDO(){
  string outdir = "../out";

  // string comment = "g1_p25_s100";
  // string callibpath = "../data-calibration/t@t_g1_p25_s100";
  // vector<string> calibrationPDO = {"run_0095.root", "run_0096.root", "run_0097.root", "run_0098.root", "run_0099.root"};

  string comment = "t@t_g3_p25_s100";
  string callibpath = "../data-calibration/t@t_g3_p25_s100";
  vector<string> calibrationPDO = {"run_0000.root", "run_0100.root", "run_0101.root", "run_0102.root", "run_0103.root", "run_0104.root"};

  auto outfile = TFile::Open(Form("%s/calibration_pdo_%s.root", outdir.c_str(), comment.c_str()), "recreate");

  vector<uint> channelsExclude = {/*0, 1, 2*/}; //0
  vector<uint> channelsDelete = {}; //35, 63
  uint nChannels = 65;

  map<string,shared_ptr<TChain>> calibrationPDOf;
  for(auto &i: calibrationPDO){
    calibrationPDOf.emplace(i, make_shared<TChain>("vmm"));
    calibrationPDOf.at(i)->Add((callibpath+"/"+i).c_str());
  }

  vector<singleTest> vectorData;
  for(auto &fileChain: calibrationPDOf){
    auto calData = findCalibrationForChannels(fileChain.second, fileChain.first, nChannels, channelsExclude, channelsDelete, 120);
    vectorData.push_back(calData);
  }

  auto pdoCorrections = fitPDO(vectorData);

  for(auto &e: channelsExclude)
    pdoCorrections.push_back({int(e), 0, 0, 1, 0, 0, 0});
  std::sort(pdoCorrections.begin(), pdoCorrections.end(), [](fitPDOChannelResult r1, fitPDOChannelResult r2){return r1.nChannel < r2.nChannel;});

  auto f = fopen(Form("%s/calibration_pdo_%s.txt", outdir.c_str(), comment.c_str()), "w");
  for(auto &correction: pdoCorrections){
    correction.print();
    correction.print(true, f);
  }
  fclose(f);
  
  for(auto &fileChain: calibrationPDOf){
    applyPDOCorrections(fileChain.second, fileChain.first, pdoCorrections);
  }

  outfile->Write();
  outfile->Close();
  
  return 0;
}
