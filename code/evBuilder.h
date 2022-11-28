#ifndef evBuilder_h
#define evBuilder_h

#include "vmm.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>

#include <iostream>

#include <vector> 
#include <map>
#include <utility>
#include <string>
#include <memory>
#include <set>

#include <algorithm> // max_element, sort
#include <functional> // std::function

using std::vector;
using std::map;
using std::string;
using std::shared_ptr, std::make_shared;
using std::pair, std::make_pair;
using std::tuple, std::get;
using std::set;
using std::array;

// Header file for the classes stored in the TTree if any.
//#include "c++/v1/vector"

class evBuilder : public vmm {
public :
   evBuilder(TString, TString runType_ = "g1_p25_s100-0&60", TString mapFile_ = "map-20220523.txt");
   evBuilder(vector<TString> filenames, TString runType_ = "g1_p25_s100-0&60", TString mapFile_ = "map-20220523.txt");
   evBuilder(TChain *tree = nullptr, TString runType_ = "g1_p25_s100-0&60", TString mapFile_ = "map-20220523.txt");
   virtual ~evBuilder();
   virtual void     Init() override;
   virtual void     Loop(unsigned long n = 0) override;
   virtual map<unsigned long, mm2CenterHitParameters> GetCentralHits(unsigned long long fromSec = 0, unsigned long long toSec = 0, bool saveOnly = false) override;
   virtual mm2CenterHitParameters GetCentralHitsData(unsigned long event) override;

   map<pair<int, int>, float> strawCenterMM = {
     {{1,24}, 156}, // 213 - 21 - 24/1.0 - 12
     {{1,25}, 165}, // 213 - 21 - 24/2.0 - 15
     {{1,26}, 181}, // 213 - 21 - 12 + 1
     {{1,27}, 189}, // 213 - 21 + 24/2.0 - 15
     {{1,28}, 204}, // 213 - 21 + 24/1.0 - 12
     {{1,29}, 216}, // 213 - 21 + 24/1.5 - 12
     {{6, 0}, 198}, // SHiP Straw // 170 before June 1st, 198 after
     {{6, 1}, 180}, // Netron Straw
    };
   map<int, string> addStrawType = {
     {0, "SHiP"},
     {1, "Neutron"},
   };

   void threePlotDrawF(TH1D *h1, TH1D *h2, TH1D *h3, TString fileEnding = "");

   struct mmHit {
     double channel;
     double pdo;
     double time;
   };
   vector<mmHit> MmCluster;
   tuple<double, double, double> getClusterParameters(double t_srtraw, double minT_straw_mm, int workType = 0);
  long long findFirstGoodPulser(unsigned long long fromSec = 0, unsigned long long toSec = 0);

  virtual vector<hitParam> getHits(unsigned long) override;

  unsigned int mmDoubleReadout;
  std::function<bool(int)> pulserPdoAccepted;
};

evBuilder::evBuilder(TString filename, TString runType_, TString mapFile_) : vmm(filename, runType_, mapFile_)
{
}
evBuilder::evBuilder(vector<TString> filenames, TString runType_, TString mapFile_): vmm(filenames, runType_, mapFile_)
{
}

evBuilder::evBuilder(TChain *tree, TString runType_, TString mapFile_) : vmm(tree, runType_, mapFile_)
{
}

evBuilder::~evBuilder()
{
}

void evBuilder::Init(){
  switch(GetTestBeam()){
    case analysisGeneral::TestBeams::TB22_October:
    case analysisGeneral::TestBeams::TB22_August:
    case analysisGeneral::TestBeams::TB22_July:
      mmDoubleReadout = 2;
      pulserPdoAccepted = [](auto pdo){return pdo == 1012;};
      break;
    case analysisGeneral::TestBeams::TB22_April:
      mmDoubleReadout = 4;
      pulserPdoAccepted = [](auto pdo){return pdo == 948 || pdo == 965;};
      break;
    default:
      pulserPdoAccepted = [](auto pdo){return pdo > 650;};
      break;
  };
  vmm::Init();
}


#endif
#ifndef evBuilder_cxx
void evBuilder::Loop(unsigned long n) {};
void evBuilder::threePlotDrawF(TH1D *h1, TH1D *h2, TH1D *h3, TString fileEnding) {};
map<unsigned long, analysisGeneral::mm2CenterHitParameters> evBuilder::GetCentralHits(unsigned long long fromSec,
                                                                                      unsigned long long toSec,
                                                                                      bool saveOnly) {
  return {};
};
long long evBuilder::findFirstGoodPulser(unsigned long long fromSec, unsigned long long toSec){return -1;}
#endif // #ifdef evBuilder_cxx
