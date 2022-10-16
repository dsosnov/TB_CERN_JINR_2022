#ifndef vmm_h
#define vmm_h

#include "analysisGeneral.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <vector>
#include <map>
#include <utility>
#include <string>
#include <memory>
#include <tuple>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>

#include "tigerTree.h"

using std::vector;
using std::map;
using std::string;
using std::shared_ptr, std::make_shared;
using std::pair, std::make_pair;
using std::tuple, std::make_tuple, std::get;
using std::set;
using std::ifstream;

class tiger : public analysisGeneral {
public :
   TString runFolder = "";
   TString mapFile = "map-tiger-empty.txt";
   enum TigerEnergyMode : bool {SampleAndHold = 0, TimeOverThreshold = 1};
   TigerEnergyMode energyMode = TigerEnergyMode::SampleAndHold;

   Char_t   gemrocID;        //                "B" == Char_t   ==  int8_t
   Short_t  tigerID;         //  8 bit data -- "S" == Short_t  == int16_t
   Char_t   chipID;          //  2 bit data -- "B" == Char_t   ==  int8_t
   Char_t   channelID;       //  6 bit data -- "B" == Char_t   ==  int8_t
   /* 4 TAC per shaper for event de-randomization */
   Char_t   tacID;           //  2 bit data -- "B" == Char_t   ==  int8_t
   Int_t    tCoarse;         // 16 bit data -- "I" == Int_t    == int32_t
   Short_t  eCoarse;         // 10 bit data -- "S" == Short_t  == int16_t
   Short_t  tFine;           // 10 bit data -- "S" == Short_t  == int16_t
   Short_t  eFine;           // 10 bit data -- "S" == Short_t  == int16_t
   Int_t    frameCount;      // 16 bit data -- "I" == Int_t    == int32_t
   Int_t    seu;             // 16 bit data -- "I" == Int_t    == int32_t
   Long64_t frameCountLoops; //                "L" == Long64_t == int64_t
   Int_t    counterWord;     // 24 bit data -- "I" == Int_t    == int32_t

   tiger(TString, TString runFolder_ = "", TString mapFile_ = "map-tiger-empty.txt", short energyMode_ = 0);
   tiger(vector<TString>, TString runFolder_ = "", TString mapFile_ = "map-tiger-empty.txt", short energyMode_ = 0);
   tiger(TChain *tree = nullptr, TString mapFile_ = "map-tiger-empty.txt", short energyMode_ = 0);
   virtual ~tiger();
   virtual void     Init() override;
   virtual void     Loop(unsigned long n = 0) override;

   map<tuple<int,int,int>, pair<int, int>> channelMap;
   void addMap(TString filename, bool verbose = false);
   pair<int,int> getMapped(const tuple<int,int,int> channel) const;
   pair<int,int> getMapped(const int gemroc, const int tiger, const int channel) const{
     return getMapped(make_tuple(gemroc, tiger, channel));
   }
   pair<int,int> getMapped(const tigerHitTL hit) const{
     return getMapped(hit.gemrocID, hit.tigerID, hit.channelID);
   }
   int getMappedDetector(const tuple<int,int,int> channel) const{
     return getMapped(channel).first;
   }
   int getMappedDetector(const int gemroc, const int tiger, const int channel) const{
     return getMappedDetector(make_tuple(gemroc, tiger, channel));
   }
   int getMappedDetector(const tigerHitTL hit) const{
     return getMappedDetector(hit.gemrocID, hit.tigerID, hit.channelID);
   }
   int getMappedChannel(const tuple<int,int,int> channel) const{
     return getMapped(channel).second;
   }
   int getMappedChannel(const int gemroc, const int tiger, const int channel) const{
     return getMappedChannel(make_tuple(gemroc, tiger, channel));
   }
   int getMappedChannel(const tigerHitTL hit) const{
     return getMappedChannel(hit.gemrocID, hit.tigerID, hit.channelID);
   }

   tigerHitTL getTigerHitTLCurrent() const;
   void updateTigerHitTLCurrent(tigerHitTL &hit) const;

  // map<pair<int, int>, float> strawCenterMM = {
  //   {{1,24}, 156}, // 213 - 21 - 24/1.0 - 12
  //   {{1,25}, 165}, // 213 - 21 - 24/2.0 - 15
  //   {{1,26}, 181}, // 213 - 21 - 12 + 1
  //   {{1,27}, 189}, // 213 - 21 + 24/2.0 - 15
  //   {{1,28}, 204}, // 213 - 21 + 24/1.0 - 12
  //   {{1,29}, 216}, // 213 - 21 + 24/1.5 - 12
  //   {{6, 0}, 198}, // SHiP Straw // 170 before June 1st, 198 after
  //   {{6, 1}, 180}, // Netron Straw
  // };

  int nDetectorTypes = 8;
};

#endif

tiger::tiger(TString filename, TString runFolder_, TString mapFile_, short energyMode_) : runFolder(runFolder_), mapFile(mapFile_), energyMode(static_cast<TigerEnergyMode>(energyMode_))
{
  file = filename;
  folder = "../data/tiger/" + runFolder + "/";
  fChain = GetTree(filename, "tigerTL");
  Init();
}

tiger::tiger(vector<TString> filenames, TString runFolder_, TString mapFile_, short energyMode_) : runFolder(runFolder_), mapFile(mapFile_), energyMode(static_cast<TigerEnergyMode>(energyMode_))
{
  file = filenames.at(0);
  folder = "../data/tiger/" + runFolder + "/";
  fChain = GetTree(filenames.at(0), "tigerTL");
  for(auto i = 1; i < filenames.size(); i++)
    fChain->Add(folder + filenames.at(i) + ending);
  Init();
}

tiger::tiger(TChain *tree, TString mapFile_, short energyMode_) : analysisGeneral(tree), mapFile(mapFile_), energyMode(static_cast<TigerEnergyMode>(energyMode_))
{
  folder = "../data/tiger/";
  fChain = (tree == nullptr) ? GetTree("", "tigerTL") : tree;
  Init();
}

tiger::~tiger()
{
}

void tiger::addMap(TString filename, bool verbose){
   ifstream infile(Form("../configs/%s", filename.Data()));
   std::string line;
   int gr, t, ch, d, dch;
   while (std::getline(infile, line))
   {
     std::istringstream iss(line);
     if(iss.str().substr(0, 1) == string("#")) // in c++20 there is starts_with("#")
       continue;
     if (!(iss >> gr >> t >> ch >> d >> dch))
       break; // error
     if(verbose)
       printf("Map: %d: %d - %d\n", ch, d, dch);
     channelMap.emplace(make_tuple(gr,t, ch), make_pair(d, dch));
   }
}
pair<int,int> tiger::getMapped(const tuple<int,int,int> channel) const{
  if(!channelMap.count(channel))
    return {-1, -1};
  return channelMap.at(channel);
}

void tiger::updateTigerHitTLCurrent(tigerHitTL &hit) const{
  hit.gemrocID = gemrocID;
  hit.tigerID = tigerID;
  hit.channelID = channelID;

  hit.chipID = chipID;
  hit.tacID = tacID;
  hit.tCoarse = tCoarse;
  hit.eCoarse = eCoarse;
  hit.tFine = tFine;
  hit.eFine = eFine;

  hit.frameCount = frameCount;
  hit.seu = seu;
  hit.frameCountLoops = frameCountLoops;

  hit.counterWord = counterWord;
}
tigerHitTL tiger::getTigerHitTLCurrent() const{
  tigerHitTL hit;
  updateTigerHitTLCurrent(hit);
  return hit;
}

void tiger::Init()
{
   // Set branch addresses and branch pointers
   if (!fChain) return;
   printf("tiger::Init()\n");
   fCurrent = -1;
   fChain->SetMakeClass(1);
   fChain->SetBranchAddress("gemrocID", &gemrocID);
   fChain->SetBranchAddress("tigerID", &tigerID);
   fChain->SetBranchAddress("chipID", &chipID);
   fChain->SetBranchAddress("channelID", &channelID);
   /* 4 TAC per shaper for event de-randomization */
   fChain->SetBranchAddress("tacID", &tacID);
   fChain->SetBranchAddress("tCoarse", &tCoarse);
   fChain->SetBranchAddress("eCoarse", &eCoarse);
   fChain->SetBranchAddress("tFine", &tFine);
   fChain->SetBranchAddress("eFine", &eFine);
   fChain->SetBranchAddress("frameCount", &frameCount);
   fChain->SetBranchAddress("seu", &seu);
   fChain->SetBranchAddress("frameCountLoops", &frameCountLoops);
   fChain->SetBranchAddress("counterWord", &counterWord);

   if(!mapFile.EndsWith(".txt"))
     mapFile.Append(".txt");
   addMap(mapFile.Data());
}

#ifndef tiger_cxx
void tiger::Loop(unsigned long n){}
#endif // #ifdef tiger_cxx
