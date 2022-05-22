#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2I.h"
#include "TFitResult.h"

#include <vector>
#include <string>
#include <memory>
#include <utility>

#include <tuple>
#include <map>

#include <algorithm>
#include <functional>

using std::vector;
using std::string;
using std::shared_ptr, std::make_shared;
using std::pair;
using std::tuple, std::make_tuple;
using std::map;

/* fec-chip-channel -> multilayer, layer, readout, strip */
// map<int, int> afterMap;
map<int, int> afterMap;

// void getMap(string mapfile = "micromegas-map-l2-center-new.txt", int chip = 9, int fec = 1){
//   ifstream infile(mapfile.c_str());
//   std::string line;
//   int channel;
//   int ml, l, s;
//   char r;
//   while (std::getline(infile, line))
//   {
//     std::istringstream iss(line);
//     if(iss.str().substr(0, 1) == string("#")) // in c++20 there is starts_with("#")
//       continue;
//     if (!(iss >> channel >> ml >> l >> r >> s))
//       break; // error
//     printf("%d, %d, %c, %d\n", ml, l, r, s);
//     channelsMap.emplace(make_tuple(fec, chip, channel), make_tuple(ml, l, r, s));
//   }
// }
// void getMap(string mapfile = "micromegas-map-l2-center.txt", string mapfileBefore = "micromegas-map-l2-center-new.txt"){
//   map<int, int> beforeMap, afterMap;
//   {
//     ifstream infile(mapfileBefore.c_str());
//     std::string line;
//     int channel;
//     int ml, l, s;
//     char r;
//     while (std::getline(infile, line))
//     {
//       std::istringstream iss(line);
//       if(iss.str().substr(0, 1) == string("#")) // in c++20 there is starts_with("#")
//         continue;
//       if (!(iss >> channel >> ml >> l >> r >> s))
//         break; // error
//       // printf("%d, %d, %c, %d\n", ml, l, r, s);
//       beforeMap.emplace(channel, s);
//     }
//   }
//   {
//     ifstream infile(mapfile.c_str());
//     std::string line;
//     int channel;
//     int ml, l, s;
//     char r;
//     while (std::getline(infile, line))
//     {
//       std::istringstream iss(line);
//       if(iss.str().substr(0, 1) == string("#")) // in c++20 there is starts_with("#")
//         continue;
//       if (!(iss >> channel >> ml >> l >> r >> s))
//         break; // error
//       // printf("%d, %d, %c, %d\n", ml, l, r, s);
//       afterMap.emplace(channel, s);
//     }
//   }
//   for(auto &k: beforeMap)
//     afterMap.emplace(k.second, afterMap.at(k.first));
// }
void getMap(string mapfile = "micromegas-map-l2-center-new.txt"){
  ifstream infile(mapfile.c_str());
  std::string line;
  int channel;
  int ml, l, s;
  char r;
  while (std::getline(infile, line))
  {
    std::istringstream iss(line);
    if(iss.str().substr(0, 1) == string("#")) // in c++20 there is starts_with("#")
      continue;
    if (!(iss >> channel >> ml >> l >> r >> s))
      break; // error
    // printf("%d, %d, %c, %d\n", ml, l, r, s);
    afterMap.emplace(channel, s);
  }
}

// tuple<int, char, int> getMapped(int chip, int fec, int ch){
//   auto srspos = make_tuple(chip, fec, ch);
//   auto outpos = make_tuple(0, '-', 0);
//   if(channelsMap.count(srspos)){
//     int sl = get<1>(channelsMap.at(srspos));
//     char sr = (get<3>(channelsMap.at(srspos)) > 0) ? get<2>(channelsMap.at(srspos)) : 'E';
//     int sp = get<2>(channelsMap.at(srspos));
//     outpos = make_tuple(sl, sr, sp);
//   }
//   return outpos;
// }

int remapAPV(string filenameIn = "../data-apv/run103.root"){
  getMap();
  string filenameOut = "_test_out.root"; //string() + "_remapped" + filenameIn;

  auto fIn = TFile::Open(filenameIn.c_str(), "read");
  
  auto tapvIn = static_cast<TTree*>(fIn->Get("apv_raw"));
  tapvIn->SetBranchStatus("*", true);
  tapvIn->SetBranchStatus("mmStrip", false);
  tapvIn->SetBranchStatus("mmReadout", false);
  tapvIn->SetBranchStatus("mmChamber", false);
  tapvIn->SetBranchStatus("mmLayer", false);
  auto tapvPedIn = static_cast<TTree*>(fIn->Get("apv_raw_ped"));
  tapvPedIn->SetBranchStatus("*", true);
  tapvPedIn->SetBranchStatus("mmStrip", false);
  tapvPedIn->SetBranchStatus("mmReadout", false);
  tapvPedIn->SetBranchStatus("mmChamber", false);
  tapvPedIn->SetBranchStatus("mmLayer", false);

  auto fOut = TFile::Open(filenameOut.c_str(), "recreate");
  
  auto tapvOut = tapvIn->CloneTree();
  tapvIn->SetBranchStatus("*", true);
  auto tapvPedOut = tapvPedIn->CloneTree();
  tapvPedIn->SetBranchStatus("*", true);

  static auto mmChamberIn = new vector<string>();
  static auto mmReadoutIn = new vector<char>();
  static auto mmStripIn = new vector<int>();
  tapvPedIn->SetBranchAddress("mmStrip", &mmStripIn);
  tapvPedIn->SetBranchAddress("mmReadout", &mmReadoutIn);
  tapvPedIn->SetBranchAddress("mmChamber", &mmChamberIn);
  static auto mmLayerIn = new vector<int>();
  static auto srsFecIn = new vector<uint>();
  static auto srsChipIn = new vector<uint>();
  static auto srsChanIn = new vector<uint>();
  tapvPedIn->SetBranchAddress("mmLayer", &mmLayerIn);
  tapvPedIn->SetBranchAddress("srsFec", &srsFecIn);
  tapvPedIn->SetBranchAddress("srsChip", &srsChipIn);
  tapvPedIn->SetBranchAddress("srsChan", &srsChanIn);

  static auto mmChamberOut = new vector<string>();
  static auto mmReadoutOut = new vector<char>();
  static auto mmStripOut = new vector<int>();
  static auto mmLayerOut = new vector<int>();
  TBranch *b_mmStripOut, *b_mmReadoutOut, *b_mmChamberOut, *b_mmLayerOut;
  b_mmStripOut = tapvPedOut->Branch("mmStrip", &mmStripOut);
  b_mmReadoutOut = tapvPedOut->Branch("mmReadout", &mmReadoutOut);
  b_mmChamberOut = tapvPedOut->Branch("mmChamber", &mmChamberOut);
  b_mmLayerOut = tapvPedOut->Branch("mmLayer", &mmLayerOut);

  {
    mmReadoutOut->clear();
    mmStripOut->clear();
    mmChamberOut->clear();
    mmLayerOut->clear();
    tapvPedIn->GetEntry(0);
    for(auto i = 0; i < mmStripIn->size(); i++){
      if(srsFecIn->at(i) == 1 && srsChipIn->at(i) == 9){
        // printf("chan %d:\n", srsChanIn->at(i)+2);
        auto newstrip = afterMap.at(srsChanIn->at(i)+2);
        if(newstrip > 0){
          mmLayerOut->push_back(2);
          mmReadoutOut->push_back('X');
          mmStripOut->push_back(newstrip);
          mmChamberOut->push_back(string("Micromegas"));
        } else {
          mmLayerOut->push_back(0);
          mmReadoutOut->push_back('E');
          mmStripOut->push_back(0);
          mmChamberOut->push_back(string("unmapped"));
        }
      } else {
        mmStripOut->push_back(mmStripIn->at(i));
        mmReadoutOut->push_back(mmReadoutIn->at(i));
        mmLayerOut->push_back(mmLayerIn->at(i));
        mmChamberOut->push_back(mmChamberIn->at(i));
      }
    }
    b_mmStripOut->Fill();
    b_mmReadoutOut->Fill();
    b_mmChamberOut->Fill();
    b_mmLayerOut->Fill();
  }

  auto nEntries = tapvOut->GetEntries();
  {
    tapvIn->SetBranchAddress("mmStrip", &mmStripIn);
    tapvIn->SetBranchAddress("mmReadout", &mmReadoutIn);
    tapvIn->SetBranchAddress("mmChamber", &mmChamberIn);
    tapvIn->SetBranchAddress("mmLayer", &mmLayerIn);
    tapvIn->SetBranchAddress("srsFec", &srsFecIn);
    tapvIn->SetBranchAddress("srsChip", &srsChipIn);
    tapvIn->SetBranchAddress("srsChan", &srsChanIn);
    b_mmChamberOut = tapvOut->Branch("mmChamber", &mmChamberOut);
    b_mmLayerOut = tapvOut->Branch("mmLayer", &mmLayerOut);
    b_mmStripOut = tapvOut->Branch("mmStrip", &mmStripOut);
    b_mmReadoutOut = tapvOut->Branch("mmReadout", &mmReadoutOut);
    b_mmStripOut->Reset();
    b_mmReadoutOut->Reset();
    for(auto j = 0; j < nEntries; j++){
      mmReadoutOut->clear();
      mmStripOut->clear();
      mmChamberOut->clear();
      mmLayerOut->clear();
      tapvIn->GetEntry(j);
      if(!(j%(int)1E3)) printf("Entry %d\n", j);
      for(auto i = 0; i < mmStripIn->size(); i++){
        if(srsFecIn->at(i) == 1 && srsChipIn->at(i) == 9){
          // printf("chan %d:\n", srsChanIn->at(i)+2);
          auto newstrip = afterMap.at(srsChanIn->at(i)+2);
          if(newstrip > 0){
            mmLayerOut->push_back(2);
            mmReadoutOut->push_back('X');
            mmStripOut->push_back(newstrip);
            mmChamberOut->push_back(string("Micromegas"));
          } else {
            mmLayerOut->push_back(0);
            mmReadoutOut->push_back('E');
            mmStripOut->push_back(0);
            mmChamberOut->push_back(string("unmapped"));
          }
        } else {
          mmStripOut->push_back(mmStripIn->at(i));
          mmReadoutOut->push_back(mmReadoutIn->at(i));
          mmLayerOut->push_back(mmLayerIn->at(i));
          mmChamberOut->push_back(mmChamberIn->at(i));
        }
      }
      b_mmStripOut->Fill();
      b_mmReadoutOut->Fill();
      b_mmChamberOut->Fill();
      b_mmLayerOut->Fill();
    }
  }
  
  tapvOut->Print();
  tapvPedOut->Print();
  tapvOut->Write("", TObject::kOverwrite);
  tapvPedOut->Write("", TObject::kOverwrite);

  
  auto dir = fIn->GetDirectory("config");
  auto dirO = fOut->mkdir("config");
  for(const auto&& k: *dir->GetListOfKeys()){
    auto key = (TKey*) k;
    key->Print();
    auto obj = key->ReadObj();
    // obj->Print();
    // obj->Write();
    dirO->WriteTObject(obj, key->GetName());
  }

  
  fOut->Write();
  fOut->Close();

  return 0;
}
