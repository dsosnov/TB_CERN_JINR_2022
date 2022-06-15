#include "calibrationPDO.h"

int calibrationPDO(){
  string outdir = "../out";

  // string comment = "t@t_g1_p25_s100_sci0&60";
  // string callibpath = "../data-calibration/t@t_g1_p25_s100_sci0&60";
  // vector<string> calibrationPDO = {"run_0243.root", "run_0244.root", "run_0245.root", "run_0246.root", "run_0247.root"};
  // auto thr=40;

  string comment = "t@t_g3_p25_s100_sci0&60";
  string callibpath = "../data-calibration/t@t_g3_p25_s100_sci0&60";
  vector<string> calibrationPDO = {"run_0248.root", "run_0249.root", "run_0250.root", "run_0251.root", "run_0252.root"};
  auto thr=120;

  auto outfile = TFile::Open(Form("%s/calibration_pdo_%s.root", outdir.c_str(), comment.c_str()), "recreate");

  vector<unsigned int> channelsExclude = {/*0, 1, 2*/}; //0
  vector<unsigned int> channelsDelete = {}; //35, 63
  unsigned int nChannels = 65;

  map<string,shared_ptr<TChain>> calibrationPDOf;
  for(auto &i: calibrationPDO){
    calibrationPDOf.emplace(i, make_shared<TChain>("vmm"));
    calibrationPDOf.at(i)->Add((callibpath+"/"+i).c_str());
  }

  vector<singleTest> vectorData;
  for(auto &fileChain: calibrationPDOf){
    auto calData = findCalibrationForChannels(fileChain.second, fileChain.first, nChannels, channelsExclude, channelsDelete, thr);
    vectorData.push_back(calData);
  }
  
  outfile->Write();
  auto pdoCorrections = fitPDO(vectorData);

  for(auto &e: channelsExclude)
    pdoCorrections.push_back({int(e), 0, 0, 1, 0, 0, 0});
  std::sort(pdoCorrections.begin(), pdoCorrections.end(), [](fitPDOChannelResult r1, fitPDOChannelResult r2){return r1.nChannel < r2.nChannel;});

  auto f = fopen(Form("%s/calibration_pdo_%s.txt", outdir.c_str(), comment.c_str()), "w");
  auto prevCh = -1;
  for(auto &correction: pdoCorrections){
    for(auto i = prevCh + 1; i < correction.nChannel; i++)
      fitPDOChannelResult({i, 0, 0, 1.0, 0, 0, 0}).print(true, f);
    correction.print();
    correction.print(true, f);
    prevCh = correction.nChannel;
  }
  fclose(f);
  
  for(auto &fileChain: calibrationPDOf){
    applyPDOCorrections(fileChain.second, fileChain.first, pdoCorrections);
  }

  outfile->Write();
  outfile->Close();
  
  return 0;
}
