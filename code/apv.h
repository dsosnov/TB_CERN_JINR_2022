// #pragma once
#ifndef apv_h
#define apv_h

#include "analysisGeneral.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <vector> 
#include <map>
#include <string>
#include <memory>

#include <algorithm> // max_element, sort

#include "apv_claster.h"

using std::vector;
using std::map;
using std::string;
using std::shared_ptr, std::make_shared;

class apv : public analysisGeneral {
public :

  TChain          *fChainPedestal;   //!pointer to the analyzed TTree or TChain
  std::map<string, TString> configs;
  int           fCurrentPedestal; //!current Tree number in a TChain

  bool isChain(){ return fChain != nullptr; }
  bool isChainPed(){ return fChainPedestal != nullptr; }

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
  apv(vector<TString>);
  apv(TChain *tree = nullptr, TChain *treePed = nullptr);
  virtual ~apv();
  virtual void     Init() override;
  virtual void     Loop() override;
  virtual void     LoopSecond(unsigned long long sec);

  static unsigned long long unique_srs_time_stamp(int, int, int);

  vector<apvClaster> clasters;
  bool addHitToClasters(int layer, int strip, short max_q, int t_max_q, vector<short> raw_q);
  bool addHitToClasters(apvHit hit);

  TTree* clasterTree;
  TBranch* clasterBranch;

  map<int, tuple<double, double, double>> shiftBetweenLayels = {
    {0, {0.255, -13.83, 3.516}}, // layer 0 - layer 1: distance between layers [m], mean strip shift (claster_l0 - claster_l1), rms of strip shift -- from run76
    {1, {0.345, -18.23, 7.712}}, // layer 1 - layer 2, distance between layers [m], mean strip shift (claster_l0 - claster_l1), rms of strip shift -- from run76
    // Distance to straw: 523
  };
  tuple<double,double,double> getHitsForTrack(pair<double, double> track);
  vector<pair<double, double>> constructTracks(vector<apvClaster> clasters);
    
};

tuple<double,double,double> apv::getHitsForTrack(pair<double, double> track){
  auto [x0, b] = track;
  auto x1 = x0 - get<1>(shiftBetweenLayels.at(0)) + b * get<0>(shiftBetweenLayels.at(0));
  auto x2 = x0 - (get<1>(shiftBetweenLayels.at(0)) + get<1>(shiftBetweenLayels.at(1))) +
    b * (get<0>(shiftBetweenLayels.at(0)) + get<0>(shiftBetweenLayels.at(1)));
  return {x0, x1, x2};
}

vector<pair<double, double>> apv::constructTracks(vector<apvClaster> clasters){
  vector<pair<double, double>> tracks = {};
  vector<apvClaster> clL0, clL1, clL2;
  std::copy_if(clasters.begin(), clasters.end(), std::back_inserter(clL0), [](auto c) { return c.getLayer() == 0; });
  std::copy_if(clasters.begin(), clasters.end(), std::back_inserter(clL1), [](auto c) { return c.getLayer() == 1; });
  std::copy_if(clasters.begin(), clasters.end(), std::back_inserter(clL2), [](auto c) { return c.getLayer() == 2; });
  if(!clL0.size() || !clL1.size())
    return tracks;

  
  // x0, b, index of clasters
  vector<tuple<double, double, vector<int>>> possibleTracks; // rough
  double bMax = 200; // 120 strips per 0.6m
  for(auto i0 = 0; i0 < clL0.size(); i0++){
    for(auto i1 = 0; i1 < clL1.size(); i1++){
      auto x0 = clL0.at(i0).center();
      auto x1 = clL1.at(i1).center();
      auto b = (x1 + get<1>(shiftBetweenLayels.at(0)) - x0) / get<0>(shiftBetweenLayels.at(0));
      auto possibleX2 = get<2>(getHitsForTrack({x0, b}));
      vector<int> hits = {i0, i1};
      for(auto i2 = 0; i2 < clL2.size(); i2++){
        if(fabs(clL2.at(i2).center() - possibleX2) < (get<2>(shiftBetweenLayels.at(0)) + get<2>(shiftBetweenLayels.at(1)))){
          hits.push_back(i2);
          break;
        }
      }
      if(fabs(b) < bMax)
        possibleTracks.push_back({x0, b, hits});
    }
  }
  std::sort(possibleTracks.begin(), possibleTracks.end(),
            [](auto t1, auto t2){return (get<2>(t1).size() != get<2>(t2).size()) ? get<2>(t1).size() > get<2>(t2).size() : fabs(get<1>(t1)) < fabs(get<1>(t2));});

  set<pair<int,int>> usedHits;
  for(auto &t: possibleTracks){
    auto [x0, b, hits] = t;
    bool h0Used = usedHits.count({0, hits.at(0)});
    bool h1Used = usedHits.count({1, hits.at(1)});

    if(h0Used || h1Used) continue;
    usedHits.emplace(make_pair(0, hits.at(0)));
    usedHits.emplace(make_pair(1, hits.at(1)));

    if(hits.size() > 2){
      bool h2Used = usedHits.count({2, hits.at(2)});
    usedHits.emplace(make_pair(2, hits.at(2)));
    }
    tracks.push_back({x0, b});
  }
  return tracks;
}

#endif

#ifdef apv_cxx
apv::apv(TString filename) : fChainPedestal(nullptr),
                             srsFec(nullptr), srsChip(nullptr), srsChan(nullptr), mmChamber(nullptr),
                             mmLayer(nullptr), mmReadout(nullptr), mmStrip(nullptr),
                             raw_q(nullptr), max_q(nullptr), t_max_q(nullptr),
                             srsFecPed(nullptr), srsChipPed(nullptr), srsChanPed(nullptr), mmChamberPed(nullptr),
                             mmLayerPed(nullptr), mmReadoutPed(nullptr), mmStripPed(nullptr),
                             ped_meanPed(nullptr), ped_stdevPed(nullptr), ped_sigmaPed(nullptr),
                             clasterTree(nullptr)
{
  file = filename;
  folder = "../data-apv/";
  fChain = GetTree(filename, "apv_raw");
  fChainPedestal = GetTree(filename, "apv_raw_ped");
  Init();
}

apv::apv(vector<TString> filenames)
{
  file = filenames.at(0);
  folder = "../data-apv/";
  fChain = GetTree(filenames.at(0), "apv_raw");
  fChainPedestal = GetTree(filenames.at(0), "apv_raw_ped");
  for(auto i = 1; i < filenames.size(); i++)
    fChain->Add(folder + filenames.at(i) + ending);
  Init();
}


apv::apv(TChain *tree, TChain *treePed) : analysisGeneral(tree), fChainPedestal(treePed),
                                        srsFec(nullptr), srsChip(nullptr), srsChan(nullptr), mmChamber(nullptr),
                                        mmLayer(nullptr), mmReadout(nullptr), mmStrip(nullptr),
                                        raw_q(nullptr), max_q(nullptr), t_max_q(nullptr),
                                        srsFecPed(nullptr), srsChipPed(nullptr), srsChanPed(nullptr), mmChamberPed(nullptr),
                                        mmLayerPed(nullptr), mmReadoutPed(nullptr), mmStripPed(nullptr),
                                        ped_meanPed(nullptr), ped_stdevPed(nullptr), ped_sigmaPed(nullptr),
                                        clasterTree(nullptr)
{
  folder = "../data-apv/";
  fChain = (tree == nullptr) ? GetTree("", "apv_raw") : tree;
  fChainPedestal = (treePed == nullptr) ? GetTree("", "apv_raw_ped") : treePed;
  Init();
}

apv::~apv(){
  if(fChainPedestal)
    delete fChainPedestal->GetCurrentFile();
}

void apv::Init(){
  printf("Init:: File: %s, tree %p, treePed %p\n", file.Data(), fChain, fChainPedestal);
  // Signal
  if(fChain){
    fChain->SetBranchAddress("evt", &evt, &b_evt);
    fChain->SetBranchAddress("error", &error, &b_error);
    fChain->SetBranchAddress("daqTimeSec", &daqTimeSec, &b_daqTimeSec);
    fChain->SetBranchAddress("daqTimeMicroSec", &daqTimeMicroSec, &b_daqTimeMicroSec);
    fChain->SetBranchAddress("srsTimeStamp", &srsTimeStamp, &b_srsTimeStamp);
    fChain->SetBranchAddress("srsTrigger", &srsTrigger, &b_srsTrigger);
    fChain->SetBranchAddress("srsFec", &srsFec, &b_srsFec);
    fChain->SetBranchAddress("srsChip", &srsChip, &b_srsChip);
    fChain->SetBranchAddress("srsChan", &srsChan, &b_srsChan);
    fChain->SetBranchAddress("mmChamber", &mmChamber, &b_mmChamber);
    fChain->SetBranchAddress("mmLayer", &mmLayer, &b_mmLayer);
    fChain->SetBranchAddress("mmReadout", &mmReadout, &b_mmReadout);
    fChain->SetBranchAddress("mmStrip", &mmStrip, &b_mmStrip);
    fChain->SetBranchAddress("raw_q", &raw_q, &b_raw_q);
    fChain->SetBranchAddress("max_q", &max_q, &b_max_q);
    fChain->SetBranchAddress("t_max_q", &t_max_q, &b_t_max_q);
  }
  // Pedestal
  if(fChainPedestal){
    // treePed->Print();
    fChainPedestal->SetBranchAddress("evt", &evtPed);
    fChainPedestal->SetBranchAddress("error", &errorPed);
    fChainPedestal->SetBranchAddress("daqTimeSec", &daqTimeSecPed);
    fChainPedestal->SetBranchAddress("daqTimeMicroSec", &daqTimeMicroSecPed);
    fChainPedestal->SetBranchAddress("srsTimeStamp", &srsTimeStampPed);
    fChainPedestal->SetBranchAddress("srsTrigger", &srsTriggerPed);
    fChainPedestal->SetBranchAddress("srsFec", &srsFecPed);
    fChainPedestal->SetBranchAddress("srsChip", &srsChipPed);
    fChainPedestal->SetBranchAddress("srsChan", &srsChanPed);
    fChainPedestal->SetBranchAddress("mmChamber", &mmChamberPed);
    fChainPedestal->SetBranchAddress("mmLayer", &mmLayerPed);
    fChainPedestal->SetBranchAddress("mmReadout", &mmReadoutPed);
    fChainPedestal->SetBranchAddress("mmStrip", &mmStripPed);
    fChainPedestal->SetBranchAddress("ped_mean", &ped_meanPed);
    fChainPedestal->SetBranchAddress("ped_stdev", &ped_stdevPed);
    fChainPedestal->SetBranchAddress("ped_sigma", &ped_sigmaPed);
    fChainPedestal->GetEntry(0);
    // printf("PedEntries: %lld\n", treePed->GetEntries());
  }

  if(clasterTree == nullptr)
    clasterTree = new TTree("clasters", "clasters");
  else
    clasterTree->Reset();
  clasterBranch = clasterTree->Branch("clasters", &clasters);
  clasterTree->AddFriend(fChain, "signals");
  
  if(fChainPedestal)
    clasterTree->AddFriend(fChainPedestal, "pedestals");
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
