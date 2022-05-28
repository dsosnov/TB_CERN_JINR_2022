#ifndef evBuilder_h
#define evBuilder_h

#include "vmm.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
//#include "c++/v1/vector"

class evBuilder : public vmm {
public :
   evBuilder(TString, TString runType_ = "g1_p25_s100-0&60", TString mapFile_ = "map-20220523.txt");
   evBuilder(vector<TString> filenames, TString runType_ = "g1_p25_s100-0&60", TString mapFile_ = "map-20220523.txt");
   evBuilder(TChain *tree = nullptr);
   virtual ~evBuilder();
   // virtual void     Init() override;
   virtual void     Loop() override;
   virtual void     LoopSecond(unsigned long long sec) override;
   virtual vector<mm2CenterHitParameters> GetCentralHits(unsigned long long fromSec = 0, unsigned long long toSec = 0) override;

   map<pair<int, int>, float> strawCenterMM = {
     {{1,24}, 156}, // 213 - 21 - 24/1.0 - 12
     {{1,25}, 165}, // 213 - 21 - 24/2.0 - 15
     {{1,26}, 181}, // 213 - 21 - 12 + 1
     {{1,27}, 189}, // 213 - 21 + 24/2.0 - 15
     {{1,28}, 204}, // 213 - 21 + 24/1.0 - 12
     {{1,29}, 216}, // 213 - 21 + 24/1.5 - 12
     {{6, 0}, 170}, // SHiP Straw
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
};

evBuilder::evBuilder(TString filename, TString runType_, TString mapFile_) : vmm(filename, runType_, mapFile_)
{
}
evBuilder::evBuilder(vector<TString> filenames, TString runType_, TString mapFile_): vmm(filenames, runType_, mapFile_)
{
}

evBuilder::evBuilder(TChain *tree) : vmm(tree)
{
}

evBuilder::~evBuilder()
{
}

#endif
#ifndef evBuilder_cxx
void evBuilder::Loop() {};
void evBuilder::LoopSecond(unsigned long long sec) {};
void evBuilder::threePlotDrawF(TH1D *h1, TH1D *h2, TH1D *h3, TString fileEnding) {};
vector<analysisGeneral::mm2CenterHitParameters> evBuilder::GetCentralHits(unsigned long long fromSec = 0,
                                                                          unsigned long long toSec = 0) {
  return {};
};
#endif // #ifdef evBuilder_cxx
