// #pragma once
#ifndef apv_h
#define apv_h

#include "analysisGeneral.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <vector> 
#include <map>
#include <utility>
#include <string>
#include <memory>

#include <algorithm> // max_element, sort

#include "apv_cluster.h"

using std::vector;
using std::map;
using std::string;
using std::shared_ptr, std::make_shared;
using std::pair;

class apv : public analysisGeneral {
public :

  TChain          *fChainPedestal;   //!pointer to the analyzed TTree or TChain
  std::map<string, TString> configs;
  int           fCurrentPedestal; //!current Tree number in a TChain

  bool isChain(){ return fChain != nullptr; }
  bool isChainPed(){ return fChainPedestal != nullptr; }

  // Event variables, signal
  unsigned long long evt; // event id, is one of: DaqTime, SrsTimeStamp, SrsTriggerNumber
  unsigned int error;
  int daqTimeSec; // DaqTime (stamped by the daq)
  int daqTimeMicroSec; // DaqTime (stamped by the daq)
  int srsTimeStamp; // srs time stamp (counter of clock cycles)
  unsigned int srsTrigger; // trigger number!
  vector<unsigned int> *srsFec;
  vector<unsigned int> *srsChip;
  vector<unsigned int> *srsChan;
  vector<string> *mmChamber;
  vector<int> *mmLayer;
  vector<char> *mmReadout;
  vector<int> *mmStrip;
  vector<vector<short>> *raw_q;
  vector<short> *max_q;
  vector<int> *t_max_q;
  /*Pedestal variable*/
  unsigned long long evtPed;
  unsigned int errorPed;
  int daqTimeSecPed;
  int daqTimeMicroSecPed;
  int srsTimeStampPed;
  unsigned int srsTriggerPed;
  vector<unsigned int> *srsFecPed;
  vector<unsigned int> *srsChipPed;
  vector<unsigned int> *srsChanPed;
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
  virtual void     Loop(unsigned long n = 0) override;
  virtual map<unsigned long, analysisGeneral::mm2CenterHitParameters> GetCentralHits(unsigned long long fromSec = 0, unsigned long long toSec = 0) override;

  struct doubleReadoutHits{
    bool sync;
    unsigned long timeSec;
    unsigned int timeMSec;
    long long timeFull() const {
      return timeMSec + timeSec * 1E6;
    }
    vector<apvHit> hits;
    vector<apvHit> hitsSync;
  };
  map<unsigned long, doubleReadoutHits> GetCentralHits2ROnly(unsigned long long fromSec = 0, unsigned long long toSec = 0);

  static unsigned long long unique_srs_time_stamp(int, int, int);

  vector<apvHit> hits;

  vector<apvCluster> clusters;
  void constructClusters();
  bool addHitToClusters(int layer, int strip, short max_q, int t_max_q, vector<short> raw_q);
  bool addHitToClusters(apvHit hit);

  TTree* clusterTree;
  TBranch* clusterBranch;

  map<int, tuple<double, double, double>> shiftBetweenLayels = {
    {0, {0.255, -13.83, 3.516}}, // layer 0 - layer 1: distance between layers [m], mean strip shift (cluster_l0 - cluster_l1), rms of strip shift -- from run76
    {1, {0.345, -18.23, 7.712}}, // layer 1 - layer 2, distance between layers [m], mean strip shift (cluster_l0 - cluster_l1), rms of strip shift -- from run76
    // Distance to straw: 523
  };
  tuple<double,double,double> getHitsForTrack(apvTrack track);
  vector<apvTrack> constructTracks(vector<apvCluster> clusters);

  unsigned int nAPVs = 5;
};

tuple<double,double,double> apv::getHitsForTrack(apvTrack track){
  auto x0 = track.intersect();
  auto b = track.slope();
  auto x1 = x0 - get<1>(shiftBetweenLayels.at(0)) + b * get<0>(shiftBetweenLayels.at(0));
  auto x2 = x0 - (get<1>(shiftBetweenLayels.at(0)) + get<1>(shiftBetweenLayels.at(1))) +
    b * (get<0>(shiftBetweenLayels.at(0)) + get<0>(shiftBetweenLayels.at(1)));
  return {x0, x1, x2};
}

vector<apvTrack> apv::constructTracks(vector<apvCluster> clusters){
  vector<apvTrack> tracks = {};
  vector<apvCluster> clL0, clL1, clL2;
  std::copy_if(clusters.begin(), clusters.end(), std::back_inserter(clL0), [](auto c) { return c.getLayer() == 0; });
  std::copy_if(clusters.begin(), clusters.end(), std::back_inserter(clL1), [](auto c) { return c.getLayer() == 1; });
  std::copy_if(clusters.begin(), clusters.end(), std::back_inserter(clL2), [](auto c) { return c.getLayer() == 2; });
  if(!clL0.size() || !clL1.size())
    return tracks;
  
  vector<apvTrack> possibleTracks; // rough
  double bMax = 200; // 120 strips per 0.6m
  for(auto &c0: clL0){
    for(auto &c1: clL1){
      auto x0 = c0.center();
      auto x1 = c1.center();
      auto b = (x1 + get<1>(shiftBetweenLayels.at(0)) - x0) / get<0>(shiftBetweenLayels.at(0));
      auto possibleX2 = get<2>(getHitsForTrack(apvTrack(x0, b)));
      vector<apvCluster> hits = {c0, c1};
      for(auto &c2: clL2){
        if(fabs(c2.center() - possibleX2) < (get<2>(shiftBetweenLayels.at(0)) + get<2>(shiftBetweenLayels.at(1)))){
          hits.push_back(c2);
          break;
        }
      }
      if(fabs(b) < bMax)
        possibleTracks.push_back(apvTrack(x0, b, hits));
    }
  }
  std::sort(possibleTracks.begin(), possibleTracks.end(),
            [](auto t1, auto t2){return (t1.nClusters() != t2.nClusters()) ? t1.nClusters() > t2.nClusters() : fabs(t1.slope()) < fabs(t2.slope());});

  set<apvCluster> usedClusters;
  for(auto &t: possibleTracks){
    bool used = false;
    for(auto &c: t.getClusters())
      if(usedClusters.count(c))
        used = true;
    // auto used = std::any_of(t.getClusters().begin(), t.getClusters().end(), [&usedClusters](auto c){return usedClusters.count(c);});
    if(used) continue;
    for(auto &c: t.getClusters())
      usedClusters.emplace(c);
    tracks.push_back(t);
  }
  return tracks;
}

apv::apv(TString filename) : fChainPedestal(nullptr),
                             srsFec(nullptr), srsChip(nullptr), srsChan(nullptr), mmChamber(nullptr),
                             mmLayer(nullptr), mmReadout(nullptr), mmStrip(nullptr),
                             raw_q(nullptr), max_q(nullptr), t_max_q(nullptr),
                             srsFecPed(nullptr), srsChipPed(nullptr), srsChanPed(nullptr), mmChamberPed(nullptr),
                             mmLayerPed(nullptr), mmReadoutPed(nullptr), mmStripPed(nullptr),
                             ped_meanPed(nullptr), ped_stdevPed(nullptr), ped_sigmaPed(nullptr),
                             clusterTree(nullptr)
{
  file = filename;
  folder = "../data-apv/";
  fChain = GetTree(filename, "apv_raw");
  fChainPedestal = GetTree(filename, "apv_raw_ped");
  Init();
}

apv::apv(vector<TString> filenames): fChainPedestal(nullptr),
                                     srsFec(nullptr), srsChip(nullptr), srsChan(nullptr), mmChamber(nullptr),
                                     mmLayer(nullptr), mmReadout(nullptr), mmStrip(nullptr),
                                     raw_q(nullptr), max_q(nullptr), t_max_q(nullptr),
                                     srsFecPed(nullptr), srsChipPed(nullptr), srsChanPed(nullptr), mmChamberPed(nullptr),
                                     mmLayerPed(nullptr), mmReadoutPed(nullptr), mmStripPed(nullptr),
                                     ped_meanPed(nullptr), ped_stdevPed(nullptr), ped_sigmaPed(nullptr),
                                     clusterTree(nullptr)
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
                                        clusterTree(nullptr)
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
  fCurrent = -1;
  // Signal
  if(fChain){
    fChain->SetMakeClass(1);
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
    fChainPedestal->SetMakeClass(1);
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

  if(clusterTree == nullptr)
    clusterTree = new TTree("clusters", "clusters");
  else
    clusterTree->Reset();
  // Temporary disabled. TODO: enable
  // clusterBranch = clusterTree->Branch("clusters", &clusters);
  // clusterTree->AddFriend(fChain, "signals");
  // clusterTree->SetMakeClass(1);
  
  if(fChainPedestal)
    clusterTree->AddFriend(fChainPedestal, "pedestals");
}

//daq epoch timestamp in upper 32 bits, 8 bits for tenths of second, followed by srs time stamp of 25 ns clock cycles counter, bizzare
unsigned long long apv::unique_srs_time_stamp(int daq_time_stamp_seconds, int daq_time_stamp_microseconds, int m_srs_time_stamp){
  auto b1 = static_cast<unsigned long long>(daq_time_stamp_seconds) << 32;
  auto b2 = static_cast<unsigned long long>(daq_time_stamp_microseconds % 100000) << 24;
  // return (((unsigned long long)(daq_time_stamp_seconds)) << 32) | (static_cast<unsigned long long>(daq_time_stamp_microseconds % 100000)) << 24 | m_srs_time_stamp);
  return (b1 | b2 | m_srs_time_stamp);
}

bool apv::addHitToClusters(int layer, int strip, short max_q, int t_max_q, vector<short> raw_q){
  return addHitToClusters({layer, strip, max_q, t_max_q, raw_q});
}

bool apv::addHitToClusters(apvHit hit){
  for(auto &c: clusters)
    if(c.addHitAdjacent(hit))
      return true;
  clusters.push_back(apvCluster(hit));
  return true;
}

void apv::constructClusters(){
  clusters.clear();
  std::sort(hits.begin(), hits.end(), [](auto h1, auto h2){return h1.strip < h2.strip;});
  for(auto &hit: hits)
    addHitToClusters(hit);

  /* Remove small clusters or clusters with small energy for the first layers only */
  clusters.erase(std::remove_if(clusters.begin(), 
                                clusters.end(),
                                [](auto c){return (c.getLayer() != 2 && c.nHits() < 3) || (c.nHits() > 90);}),
                 clusters.end());
  clusters.erase(std::remove_if(clusters.begin(), 
                                clusters.end(),
                                [](auto c){
                                  auto qthr = (c.getLayer() == 2) ? 100 : 300;
                                  return c.maxQ() < qthr;}),
                 clusters.end());

  std::stable_sort(clusters.begin(), clusters.end(), [](auto c1, auto c2){return (c1.getLayer() < c2.getLayer()) || (c1.center() < c2.center());});
}



#endif

#ifndef apv_cxx
void apv::Loop(unsigned long n) {};
map<unsigned long, analysisGeneral::mm2CenterHitParameters> apv::GetCentralHits(unsigned long long fromSec,
                                                                                unsigned long long toSec) {
  return {};
};
map<unsigned long, apv::doubleReadoutHits> apv::GetCentralHits2ROnly(unsigned long long fromSec,
                                                                 unsigned long long toSec){
  return {};
}
#endif // #ifdef apv_cxx
