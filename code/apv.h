// #pragma once
#ifndef apv_h
#define apv_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <vector> 
#include <map>
#include <string>
#include <memory>

#include <numeric> // accumulate
#include <algorithm> // max_element, sort

using std::vector;
using std::map;
using std::string;
using std::shared_ptr, std::make_shared;

struct apvHit{
  int layer = 0;
  int strip = 0;
  short max_q;
  int t_max_q;
  vector<short> raw_q;
};

class apvClaster{
private:
  int layer;
  vector<apvHit> hits;
  unsigned long sizeOnLastUpdate;
  // float center_;
  // int width_;
  int maxq_;
  long qsum_;
public:
  apvClaster(int clasterLayer = 0):
    layer(clasterLayer), hits({}), sizeOnLastUpdate(-1){
  }
  apvClaster(apvHit hit):
    layer(hit.layer), hits({}), sizeOnLastUpdate(-1){
    addHit(hit);
  }
  ~apvClaster(){
    hits.clear();
  }
  int getLayer(){return layer;}
  vector<apvHit> getHits(){return hits;}
  bool addHit(apvHit hit){
    if(hit.layer != layer)
      return false;
    hits.push_back(hit);
    sortHits();
    return true;
  }
  bool addHitAdjacent(apvHit hit){
    if(hit.layer != layer)
      return false;
    if(hit.strip != hits.at(0).strip - 1 && hit.strip != hits.back().strip + 1)
      return false;
    addHit(hit);
    return true;
  }
  float center(){
    sortHits();
    return (hits.at(0).strip + hits.back().strip)/2.0;
  }
  int width(){
    sortHits();
    return hits.back().strip - hits.at(0).strip + 1;
  }
  int firstStrip(){ sortHits(); return hits.at(0).strip; }
  int lastStrip(){ sortHits(); return hits.back().strip; }
  int maxQ(){
    if(!sizeOnLastUpdate && sizeOnLastUpdate == nHits())
      return maxq_;
    maxq_ = -1;
    for(auto &hit: hits)
      if(hit.max_q > maxq_)
        maxq_ = hit.max_q;
    return maxq_;
  }
  long q(){
    if(!sizeOnLastUpdate && sizeOnLastUpdate == nHits())
      return qsum_;    
    qsum_ = 0;
    for(auto &hit: hits)
      qsum_ += hit.max_q;
    return qsum_;
  }
  unsigned long nHits(){ return hits.size(); }
  void print(){
    printf("Claster: %lu hits on layer %d, with center %.2f, width %d and Q %ld (maximal: %d)\n", nHits(), layer, center(), width(), q(), maxQ()); 
  }
  void sortHits(){
    if(!sizeOnLastUpdate && sizeOnLastUpdate == nHits())
      return;
    std::sort(hits.begin(), hits.end(), [](const apvHit h1, const apvHit h2){return (h1.strip < h2.strip);});
    sizeOnLastUpdate = nHits();
  }
  bool merge(apvClaster claster){
    if(claster.layer != layer)
      return false;
    if(claster.lastStrip() != hits.at(0).strip - 1 && claster.firstStrip() != hits.back().strip + 1)
      return false;
    hits.insert(hits.end(), claster.hits.begin(), claster.hits.end());
    sortHits();
    return true;
  }
};

class apv {
public :
  TString folder = "../data-apv/";
  TString file = "run16";
  TString ending = ".root";

  TTree          *fChainSignal;   //!pointer to the analyzed TTree or TChain
  TTree          *fChainPedestal;   //!pointer to the analyzed TTree or TChain
  std::map<string, TString> configs;
  int           fCurrentData; //!current Tree number in a TChain
  int           fCurrentPedestal; //!current Tree number in a TChain

  // Event variables, signal
  unsigned long long evt; // event id, is one of: DaqTime, SrsTimeStamp, SrsTriggerNumber
  uint error;
  int daqTimeSec; // DaqTime (stamped by the daq)
  int daqTimeMicroSec; // DaqTime (stamped by the daq)
  int srsTimeStamp; // srs time stamp (counter of clock cycles)
  uint srsTrigger; // trigger number!
  vector<uint> *srsFec;
  vector<uint> *srsChip;
  vector<uint> *srsChan;
  vector<string> *mmChamber;
  vector<int> *mmLayer;
  vector<char> *mmReadout;
  vector<int> *mmStrip;
  vector<vector<short>> *raw_q;
  vector<short> *max_q;
  vector<int> *t_max_q;
  /*Pedestal variable*/
  unsigned long long evtPed;
  uint errorPed;
  int daqTimeSecPed;
  int daqTimeMicroSecPed;
  int srsTimeStampPed;
  uint srsTriggerPed;
  vector<uint> *srsFecPed;
  vector<uint> *srsChipPed;
  vector<uint> *srsChanPed;
  vector<string> *mmChamberPed;
  vector<int> *mmLayerPed;
  vector<char> *mmReadoutPed;
  vector<int> *mmStripPed;
  vector<double> *ped_meanPed;
  vector<double> *ped_stdevPed;
  vector<double> *ped_sigmaPed;

  // List of branches
  // Signal
  TBranch *b_evt;
  TBranch *b_error;
  TBranch *b_daqTimeSec;
  TBranch *b_daqTimeMicroSec;
  TBranch *b_srsTimeStamp;
  TBranch *b_srsTrigger;
  TBranch *b_srsFec;
  TBranch *b_srsChip;
  TBranch *b_srsChan;
  TBranch *b_mmChamber;
  TBranch *b_mmLayer;
  TBranch *b_mmReadout;
  TBranch *b_mmStrip;
  TBranch *b_raw_q;
  TBranch *b_max_q;
  TBranch *b_t_max_q;

  apv(TString);
  apv(TTree *tree = nullptr, TTree *treePed = nullptr);
  virtual ~apv();
  virtual int    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree, TTree *treePed);
  virtual void     Loop();

  static unsigned long long unique_srs_time_stamp(int, int, int);

  vector<apvClaster> clasters;
  bool addHitToClasters(int layer, int strip, short max_q, int t_max_q, vector<short> raw_q);
  bool addHitToClasters(apvHit hit);

};

#endif

#ifdef apv_cxx
apv::apv(TString filename) : file(filename),
                             fChainSignal(nullptr), fChainPedestal(nullptr),
                             srsFec(nullptr), srsChip(nullptr), srsChan(nullptr), mmChamber(nullptr),
                             mmLayer(nullptr), mmReadout(nullptr), mmStrip(nullptr),
                             raw_q(nullptr), max_q(nullptr), t_max_q(nullptr),
                             srsFecPed(nullptr), srsChipPed(nullptr), srsChanPed(nullptr), mmChamberPed(nullptr),
                             mmLayerPed(nullptr), mmReadoutPed(nullptr), mmStripPed(nullptr),
                             ped_meanPed(nullptr), ped_stdevPed(nullptr), ped_sigmaPed(nullptr)
{
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(folder + file + ending);
  if(!f || !f->IsOpen()) {
    f = new TFile(folder + file + ending);
  }
  if(!f->IsOpen()) {
    std::cout << "Problem with opening data file" << std::endl;
    exit(1);
  }
  TTree* tree = nullptr; 
  f->GetObject("apv_raw",tree);
  TTree* treePed = nullptr; 
  f->GetObject("apv_raw_ped",treePed);
  Init(tree, treePed);
}

apv::apv(TTree *tree, TTree *treePed) : fChainSignal(nullptr), fChainPedestal(nullptr),
                             srsFec(nullptr), srsChip(nullptr), srsChan(nullptr), mmChamber(nullptr),
                             mmLayer(nullptr), mmReadout(nullptr), mmStrip(nullptr),
                             raw_q(nullptr), max_q(nullptr), t_max_q(nullptr),
                             srsFecPed(nullptr), srsChipPed(nullptr), srsChanPed(nullptr), mmChamberPed(nullptr),
                             mmLayerPed(nullptr), mmReadoutPed(nullptr), mmStripPed(nullptr),
                             ped_meanPed(nullptr), ped_stdevPed(nullptr), ped_sigmaPed(nullptr)
{
  if(tree == nullptr) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(folder + file + ending);
    if(!f || !f->IsOpen()) {
      f = new TFile(folder + file + ending);
    }
    if(!f->IsOpen()) {
      std::cout << "Problem with opening data file" << std::endl;
      exit(1);
    }
    f->GetObject("apv_raw",tree);
    if(treePed == nullptr)
      f->GetObject("apv_raw_ped",treePed);
  }
  Init(tree, treePed);
}

apv::~apv(){
  if(fChainSignal)
    delete fChainSignal->GetCurrentFile();
  if(fChainPedestal)
    delete fChainPedestal->GetCurrentFile();
}

int apv::GetEntry(Long64_t entry){
  // Read contents of entry.
  if(!fChainSignal) return 0;
  return fChainSignal->GetEntry(entry);
}
Long64_t apv::LoadTree(Long64_t entry){
  // Set the environment to read one entry
  if(!fChainSignal) return -5;
  Long64_t centry = fChainSignal->LoadTree(entry);
  if(centry < 0) return centry;
  if(fChainSignal->GetTreeNumber() != fCurrentData) {
    fCurrentData = fChainSignal->GetTreeNumber();
  }
  return centry;
}

void apv::Init(TTree *tree, TTree *treePed){
  fChainSignal = tree;
  fChainPedestal = treePed;
  printf("Init:: File: %s, tree %p, treePed %p\n", file.Data(), tree, treePed);
  // Signal
  tree->SetBranchAddress("evt", &evt, &b_evt);
  tree->SetBranchAddress("error", &error, &b_error);
  tree->SetBranchAddress("daqTimeSec", &daqTimeSec, &b_daqTimeSec);
  tree->SetBranchAddress("daqTimeMicroSec", &daqTimeMicroSec, &b_daqTimeMicroSec);
  tree->SetBranchAddress("srsTimeStamp", &srsTimeStamp, &b_srsTimeStamp);
  tree->SetBranchAddress("srsTrigger", &srsTrigger, &b_srsTrigger);
  tree->SetBranchAddress("srsFec", &srsFec, &b_srsFec);
  tree->SetBranchAddress("srsChip", &srsChip, &b_srsChip);
  tree->SetBranchAddress("srsChan", &srsChan, &b_srsChan);
  tree->SetBranchAddress("mmChamber", &mmChamber, &b_mmChamber);
  tree->SetBranchAddress("mmLayer", &mmLayer, &b_mmLayer);
  tree->SetBranchAddress("mmReadout", &mmReadout, &b_mmReadout);
  tree->SetBranchAddress("mmStrip", &mmStrip, &b_mmStrip);
  tree->SetBranchAddress("raw_q", &raw_q, &b_raw_q);
  tree->SetBranchAddress("max_q", &max_q, &b_max_q);
  tree->SetBranchAddress("t_max_q", &t_max_q, &b_t_max_q);
  // Pedestal
  if(treePed){
    // treePed->Print();
    treePed->SetBranchAddress("evt", &evtPed);
    treePed->SetBranchAddress("error", &errorPed);
    treePed->SetBranchAddress("daqTimeSec", &daqTimeSecPed);
    treePed->SetBranchAddress("daqTimeMicroSec", &daqTimeMicroSecPed);
    treePed->SetBranchAddress("srsTimeStamp", &srsTimeStampPed);
    treePed->SetBranchAddress("srsTrigger", &srsTriggerPed);
    treePed->SetBranchAddress("srsFec", &srsFecPed);
    treePed->SetBranchAddress("srsChip", &srsChipPed);
    treePed->SetBranchAddress("srsChan", &srsChanPed);
    treePed->SetBranchAddress("mmChamber", &mmChamberPed);
    treePed->SetBranchAddress("mmLayer", &mmLayerPed);
    treePed->SetBranchAddress("mmReadout", &mmReadoutPed);
    treePed->SetBranchAddress("mmStrip", &mmStripPed);
    treePed->SetBranchAddress("ped_mean", &ped_meanPed);
    treePed->SetBranchAddress("ped_stdev", &ped_stdevPed);
    treePed->SetBranchAddress("ped_sigma", &ped_sigmaPed);
    treePed->GetEntry(0);
    // printf("PedEntries: %lld\n", treePed->GetEntries());
  }
}

//daq epoch timestamp in upper 32 bits, 8 bits for tenths of second, followed by srs time stamp of 25 ns clock cycles counter, bizzare
unsigned long long apv::unique_srs_time_stamp(int daq_time_stamp_seconds, int daq_time_stamp_microseconds, int m_srs_time_stamp){
  auto b1 = static_cast<unsigned long long>(daq_time_stamp_seconds) << 32;
  auto b2 = static_cast<unsigned long long>(daq_time_stamp_microseconds % 100000) << 24;
  // return (((unsigned long long)(daq_time_stamp_seconds)) << 32) | (static_cast<unsigned long long>(daq_time_stamp_microseconds % 100000)) << 24 | m_srs_time_stamp);
  return (b1 | b2 | m_srs_time_stamp);
}

bool apv::addHitToClasters(int layer, int strip, short max_q, int t_max_q, vector<short> raw_q){
  return addHitToClasters({layer, strip, max_q, t_max_q, raw_q});
}

bool apv::addHitToClasters(apvHit hit){
  for(auto &c: clasters)
    if(c.addHitAdjacent(hit))
      return true;
  clasters.push_back(apvClaster(hit));
  return true;
}



#endif // #ifdef apv_cxx
