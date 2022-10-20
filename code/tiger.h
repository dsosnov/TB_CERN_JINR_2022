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
#include <optional>
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
using std::optional, std::nullopt;
using std::set;
using std::ifstream;

class tiger : public analysisGeneral {
public :
   TString runFolder = "";
   TString mapFile = "map-tiger-empty.txt";
   TString efineFile = "";
   TString tfineFile = "";
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

   tiger(TString, TString runFolder_, TString mapFile_, TString eFineFile_, TString tFineFile_, short energyMode_ = 0);
   tiger(TString, TString runFolder_ = "", TString mapFile_ = "map-tiger-empty.txt", TString calibration = "", short energyMode_ = 0);
   tiger(vector<TString>, TString runFolder_ = "", TString mapFile_ = "map-tiger-empty.txt", TString eFineFile_ = "", short energyMode_ = 0);
   tiger(TChain *tree = nullptr, TString mapFile_ = "map-tiger-empty.txt", TString eFineFile_ = "", short energyMode_ = 0);
   virtual ~tiger();
   virtual void     Init() override;
   virtual void     Loop(unsigned long n = 0) override;

   map<tuple<int,int,int>, pair<int, int>> channelMap;
   void addMap(TString filename, bool verbose = false);
   pair<int,int> getMapped(const tuple<int,int,int> channel) const;
   pair<int,int> getMapped(const int gemroc, const int tiger, const int channel) const{
     return getMapped(make_tuple(gemroc, tiger, channel));
   }
   pair<int,int> getMapped(const tigerHitTL* hit) const{
     return getMapped(hit->gemrocID, hit->tigerID, hit->channelID);
   }
   pair<int,int> getMapped(const tigerHitTL hit) const{
     return getMapped(&hit);
   }
   int getMappedDetector(const tuple<int,int,int> channel) const{
     return getMapped(channel).first;
   }
   int getMappedDetector(const int gemroc, const int tiger, const int channel) const{
     return getMappedDetector(make_tuple(gemroc, tiger, channel));
   }
   int getMappedDetector(const tigerHitTL* hit) const{
     return getMappedDetector(hit->gemrocID, hit->tigerID, hit->channelID);
   }
   int getMappedDetector(const tigerHitTL hit) const{
     return getMappedDetector(&hit);
   }
   int getMappedChannel(const tuple<int,int,int> channel) const{
     return getMapped(channel).second;
   }
   int getMappedChannel(const int gemroc, const int tiger, const int channel) const{
     return getMappedChannel(make_tuple(gemroc, tiger, channel));
   }
   int getMappedChannel(const tigerHitTL* hit) const{
     return getMappedChannel(hit->gemrocID, hit->tigerID, hit->channelID);
   }
   int getMappedChannel(const tigerHitTL hit) const{
     return getMappedChannel(&hit);
   }

   map<tuple<int,int,int>, pair<int, int>> eFineMap;
   void addCalibrationEFine(TString filename, bool verbose = false);
   map<tuple<int,int,int>, pair<int, int>> tFineMap;
   void addCalibrationTFine(TString filename, bool verbose = false);

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
  map<long long, tigerHitTL> hitsMap;
  tigerHitTL* getHitFromTree(long long);
  void freeHitMap(long long);
  int mmLayerY = 5;
};

#endif

tiger::tiger(TString filename, TString runFolder_, TString mapFile_, TString eFineFile_, TString tFineFile_, short energyMode_) : runFolder(runFolder_), mapFile(mapFile_), efineFile(eFineFile_), tfineFile(tFineFile_), energyMode(static_cast<TigerEnergyMode>(energyMode_))
{
  file = filename;
  folder = "../data/tiger/" + runFolder + "/";
  fChain = GetTree(filename, "tigerTL");
  Init();
}
tiger::tiger(TString filename, TString runFolder_, TString mapFile_, TString calibration, short energyMode_) : runFolder(runFolder_), mapFile(mapFile_), efineFile(calibration), tfineFile(calibration),
                                                                                                               energyMode(static_cast<TigerEnergyMode>(energyMode_)){
  file = filename;
  folder = "../data/tiger/" + runFolder + "/";
  fChain = GetTree(filename, "tigerTL");
  Init();
}

tiger::tiger(vector<TString> filenames, TString runFolder_, TString mapFile_, TString eFineFile_, short energyMode_) : runFolder(runFolder_), mapFile(mapFile_), efineFile(eFineFile_), energyMode(static_cast<TigerEnergyMode>(energyMode_))
{
  file = filenames.at(0);
  folder = "../data/tiger/" + runFolder + "/";
  fChain = GetTree(filenames.at(0), "tigerTL");
  for(auto i = 1; i < filenames.size(); i++)
    fChain->Add(folder + filenames.at(i) + ending);
  Init();
}

tiger::tiger(TChain *tree, TString mapFile_, TString eFineFile_, short energyMode_) : analysisGeneral(tree), mapFile(mapFile_), efineFile(eFineFile_), energyMode(static_cast<TigerEnergyMode>(energyMode_))
{
  folder = "../data/tiger/";
  fChain = (tree == nullptr) ? GetTree("", "tigerTL") : tree;
  Init();
}

tiger::~tiger()
{
}

void tiger::addMap(TString filename, bool verbose){
   if(filename == "") return;
   ifstream infile(Form("../configs/%s", filename.Data()));
   if(infile.fail()){
      if(!filename.BeginsWith("map-tiger-")){
         auto fn2 = TString("map-tiger-") + filename;
         printf("No map file %s found. Try to find file \"%s\"\n", filename.Data(), fn2.Data());
         addMap(fn2, verbose);
      }
      return;
   }
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

void tiger::addCalibrationEFine(TString filename, bool verbose){
   if(filename == "") return;
   ifstream infile(Form("../configs/%s", filename.Data()));
   if(infile.fail()){
      if(!filename.BeginsWith("tiger_efine_calibration-")){
         auto fn2 = TString("tiger_efine_calibration-") + filename;
         printf("No eFine calibration file %s found. Try to find file \"%s\"\n", filename.Data(), fn2.Data());
         addCalibrationEFine(fn2, verbose);
      }
      return;
   }
   std::string line;
   int gr, t, ch, min, max;
   while (std::getline(infile, line))
   {
     std::istringstream iss(line);
     if(iss.str().substr(0, 1) == string("#")) // in c++20 there is starts_with("#")
       continue;
     if (!(iss >> gr >> t >> ch >> min >> max))
       break; // error
     if(verbose)
       printf("efine calibration: %d %d: %d: %d - %i\n", gr, t, ch, min, max);
     eFineMap.emplace(make_tuple(gr,t, ch), make_pair(min, max));
   }
}
void tiger::addCalibrationTFine(TString filename, bool verbose){
   if(filename == "") return;
   ifstream infile(Form("../configs/%s", filename.Data()));
   if(infile.fail()){
      if(!filename.BeginsWith("tiger_tfine_calibration-")){
         auto fn2 = TString("tiger_tfine_calibration-") + filename;
         printf("No tFine calibration file %s found. Try to find file \"%s\"\n", filename.Data(), fn2.Data());
         addCalibrationEFine(fn2, verbose);
      }
      return;
   }
   std::string line;
   int gr, t, ch, min, max;
   while (std::getline(infile, line))
   {
     std::istringstream iss(line);
     if(iss.str().substr(0, 1) == string("#")) // in c++20 there is starts_with("#")
       continue;
     if (!(iss >> gr >> t >> ch >> min >> max))
       break; // error
     if(verbose)
       printf("efine calibration: %d %d: %d: %d - %i\n", gr, t, ch, min, max);
     tFineMap.emplace(make_tuple(gr,t, ch), make_pair(min, max));
   }
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
  hit.tFineLimits = (tFineMap.count({gemrocID, tigerID, channelID})) ?
     tFineMap.at({gemrocID, tigerID, channelID}) : make_pair(0, 1023);
}
tigerHitTL tiger::getTigerHitTLCurrent() const{
  tigerHitTL hit;
  updateTigerHitTLCurrent(hit);
  return hit;
}
tigerHitTL* tiger::getHitFromTree(long long entry){
  if(!hitsMap.count(entry)){
    fChain->GetEntry(entry);
    hitsMap.emplace(entry, getTigerHitTLCurrent());
  }
  return &hitsMap.at(entry);
}
void tiger::freeHitMap(long long minEntry){
  if(!hitsMap.size()) return;
  auto it = hitsMap.lower_bound(minEntry);
  hitsMap.erase(hitsMap.begin(), --it);
  // auto it = hitsMap.begin();
  // auto itEnd = hitsMap.end();
  // for(; it != itEnd; ) {
  //   if (it->first < minEntry) {
  //     it = hitsMap.erase(it);
  //   } else {
  //     ++it;
  //   }
  // }
// For those on C++20 there are built-in std::erase_if functions for map and unordered_map:
// const auto count = std::erase_if(data, [](const auto& item) { auto const& [key, value] = item; return (key & 1) == 1;});
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

   if(mapFile!="" && !mapFile.EndsWith(".txt"))
     mapFile.Append(".txt");
   addMap(mapFile.Data());
   if(efineFile!="" && !efineFile.EndsWith(".txt"))
     efineFile.Append(".txt");
   addCalibrationEFine(efineFile.Data());
   if(tfineFile!="" && !tfineFile.EndsWith(".txt"))
     tfineFile.Append(".txt");
   addCalibrationTFine(tfineFile.Data());
}

#ifndef tiger_cxx
void tiger::Loop(unsigned long n){}
#endif // #ifdef tiger_cxx
