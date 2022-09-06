#include <cstdint>

#include "Rtypes.h"
#include "TString.h"
#include "TObjString.h"
#include "TFile.h"
#include "TSystemDirectory.h"
#include "TRegexp.h"
#include "TTree.h"

#include <string>
#include <sstream>
#include <vector>
#include <map>
using std::vector;
using std::string;
using std::map;

#include <fstream>
using std::ifstream;

// With enabled padding (by default) every type
// may starts with new byte.
// Solution: disabling padding (in commented block) or
// or moving all to single uint64_t
// /* pragma pack desables padding */
// #pragma pack(push, 1)
// struct TLFrameWord{
//   uint16_t seuCount: 15;
//   uint16_t frameCount: 16;
//   uint32_t reserved: 25;
//   uint8_t tiger: 3;
//   uint8_t key: 5; // should be 0x4
// };
// #pragma pack(pop)
struct TLFrameWord{
  uint64_t seuCount: 15,
    frameCount: 16,
    reserved: 25,
    tiger: 3,
    key: 5; // should be 0x4
};
struct TLEventWord{
  uint64_t eFine: 10,
    tFine: 10,
    eCoarse: 10,
    tCoarse: 16,
    tac: 2,
    channel: 6,
    reserved: 2,
    tiger: 3,
    key: 5; // should be 0x0
};
struct TLCountWord{
  uint64_t counter: 24,
    channel: 6,
    reserved: 26,
    tiger: 3,
    key: 5; // should be 0x8
};
enum class TLHitType{unknown, FrameWord, EventWord, CountWord};
TLHitType getTypeTL(const int64_t &hit){
  if(reinterpret_cast<const TLFrameWord*>(&hit)->key == 0x4)
    return TLHitType::FrameWord;
  else if(reinterpret_cast<const TLEventWord*>(&hit)->key == 0x0)
    return TLHitType::EventWord;
  else if(reinterpret_cast<const TLCountWord*>(&hit)->key == 0x8)
    return TLHitType::CountWord;
  return TLHitType::unknown;
}

/* pragma pack desables padding */
#pragma pack(push, 1)
struct TMHeader{
  uint16_t l1Timestamp: 16;  
  uint8_t countHits: 8;

  // in python converter:
  uint32_t l1LocalCount: 32;
  uint8_t reserved: 2;
  // in tiger_data_format.pdf:
  // uint64_t l1LocalCount: 34;

  uint8_t status: 3;
  uint8_t key: 3; // should be 0x6
};
struct TMTrailer{
  uint32_t lastCountWordData: 18;
  uint8_t lastCountWordCh: 6;
  uint8_t l1LocalCount: 3;
  uint8_t tiger: 3;
  uint8_t reserved: 2;
  uint8_t gemroc: 5;
  uint32_t l1LocalFramenum: 24;
  uint8_t key: 3; // should be 0x7
};
struct TMData{
  // in tiger_data_format.pdf:
  // uint64_t rawData: 56;
  // in python converter:
  uint16_t eFine: 10;
  uint16_t tFine: 10;
  uint16_t eCoarse: 10;
  uint8_t reserved: 2;
  uint16_t tCoarse: 16;
  uint8_t tac: 2;
  uint8_t channel: 6;
  uint8_t lastTigerFrameNumber: 3;
  uint8_t tiger: 3;
  uint8_t key: 2; // should be 0x0
};
struct TMUDPCounter{
  uint32_t reserved: 28;
  uint32_t udpFrameCount: 24;
  uint8_t gemroc: 5;
  uint32_t headerStatus: 3;
  uint8_t key: 4; // should be 0x4
};
#pragma pack(pop)
enum class TMHitType{unknown, Header, Trailer, Data, UDPCounter};
TMHitType getTypeTM(const int64_t &hit){
  if(reinterpret_cast<const TMTrailer*>(&hit)->key == 0x7)
    return TMHitType::Trailer;
  else if(reinterpret_cast<const TMHeader*>(&hit)->key == 0x6)
    return TMHitType::Header;
  else if(reinterpret_cast<const TMData*>(&hit)->key == 0x0)
    return TMHitType::Data;
  else if(reinterpret_cast<const TMUDPCounter*>(&hit)->key == 0x4)
    return TMHitType::UDPCounter;
  return TMHitType::unknown;
}

void convertTL(ifstream* fIn, Char_t gemroc){
  Char_t gemrocID = gemroc; // -- "B" = Char_t == int8_t
  Short_t tigerID; // 8 bit -- "S" == Short_t == int16_t
  Char_t chipID = 0; // 2 bit -- "B" = Char_t == int8_t
  Char_t channelID; // 6 bit -- "B" = Char_t == int8_t
  Char_t tacID; // 2 bit -- 4 TAC per shaper for event de-randomization -- "B" = Char_t == int8_t
  Int_t tCoarse; // 16 bit data -- "I" = Int_t == int32_t
  Short_t eCoarse; // 10 bit data -- "S" == Short_t == int16_t
  Short_t tFine; // 10 bit data -- "S" == Short_t == int16_t
  Short_t eFine; // 10 bit data -- "S" == Short_t == int16_t
  Int_t frameCount = 0; // 16 bit data -- "I" = Int_t == int32_t
  Int_t seu = 0; // 16 bit data -- "I" = Int_t == int32_t
  Long64_t frameCountLoops = 0; // -- "L" = Long64_t == int64_t
  Int_t counterWord = 0; // 24 bit data -- "I" = Int_t == int32_t

  auto tree = new TTree("tigerTL", "tigerTL");
  tree->Branch("gemrocID", &gemrocID, "gemrocID/B");
  tree->Branch("tigerID", &tigerID, "tigerID/S");
  tree->Branch("chipID", &chipID, "chipID/B");
  tree->Branch("channelID", &channelID, "channelID/B");
  tree->Branch("tacID", &tacID, "tacID/B");
  tree->Branch("tCoarse", &tCoarse, "tCoarse/I");
  tree->Branch("eCoarse", &eCoarse, "eCoarse/S");
  tree->Branch("tFine", &tFine, "tFine/S");
  tree->Branch("eFine", &eFine, "eFine/S");
  tree->Branch("frameCount", &frameCount, "frameCount/I");
  tree->Branch("seu", &seu, "seu/I");
  tree->Branch("frameCountLoops", &frameCountLoops, "frameCountLoops/L");
  tree->Branch("counterWord", &counterWord, "counterWord/I");

  vector<Long64_t> FCRollOvers(256, 0);
  vector<Int_t> lastFrameCount(256, 0);
  vector<Int_t> lastSEU(256, 0);
  vector<Int_t> lastCW(256*64, 0);
  int64_t frame = 0;
  auto frameFW = reinterpret_cast<TLFrameWord*>(&frame); // (TLFrameWord*)(&frame);
  auto frameEW = reinterpret_cast<TLEventWord*>(&frame); // (TLEventWord*)(&frame);
  auto frameCW = reinterpret_cast<TLCountWord*>(&frame); // (TLCountWord*)(&frame);
  auto frameType = TLHitType::unknown;
  while (!fIn->eof()){
    fIn->read(reinterpret_cast<char*>(&frame), sizeof(frame));    
    frameType = getTypeTL(frame);
    switch(frameType){
      case TLHitType::EventWord: { // EW
        tigerID = frameEW->tiger;
        channelID = frameEW->channel;
        tacID = frameEW->tac;
        tCoarse = frameEW->tCoarse;
        eCoarse = frameEW->eCoarse;
        tFine = frameEW->tFine;
        eFine = frameEW->eFine;
        frameCount = lastFrameCount.at(tigerID);
        seu = lastSEU.at(tigerID);
        frameCountLoops = FCRollOvers.at(tigerID);
        counterWord = lastCW.at(tigerID*64 + channelID);
        tree->Fill();
        // printf("TIGER %01X: EW: ChID: %02X tacID: %01X Tcoarse: %04X Ecoarse: %03X Tfine: %03X Efine: %03X \n",
        //        tigerID, channelID, tacID, tCoarse, eCoarse, tFine, eFine);
        break;
      }
      case TLHitType::FrameWord: { // HB
        tigerID = frameFW->tiger;
        frameCount = frameFW->frameCount;
        seu = frameFW->seuCount;
        if(frameCount < lastFrameCount.at(tigerID))
          FCRollOvers.at(tigerID)++;
        lastSEU.at(tigerID) = seu;
        lastFrameCount.at(tigerID) = frameCount;
        // printf("TIGER %01X: HB: Framecount: %08X SEUcount: %08X\n", tigerID, frameCount, seu);
        break;
      }
      case TLHitType::CountWord: { // CW
        tigerID = frameCW->tiger;
        channelID = frameCW->channel;
        counterWord = frameCW->counter;
        lastCW.at(tigerID*64 + channelID) = counterWord;
        // printf("TIGER %01X: CW: ChID: %02X CounterWord: %016X\n", tigerID, channelID, counterWord);
        break;
      }
      case TLHitType::unknown: {
        break;
      }
    }
  }
  tree->Write();
  delete tree;
}

void tiger_tree_converter_tl(string folderName){
  
  auto dir = new TSystemDirectory(folderName.c_str(), folderName.c_str());
  auto files = dir->GetListOfFiles();
  if(!files){
    delete dir;
    return;
  }

  vector<string> datoutFiles;

  TString name;
  for(auto f: *files){
    name = f->GetName();
    if (name.Index(TRegexp(Form("SubRUN_*_GEMROC_*_TL.dat"),kTRUE)) == kNPOS)
        continue;
    datoutFiles.push_back(name.Data());
    printf("%s\n", name.Data());
  }

  for(auto &fn: datoutFiles){
    int gemroc = static_cast<TObjString*>(TString(fn).Tokenize("_")->At(3))->String().Atoi();
    auto fIn = new ifstream((folderName+"/"+fn).c_str(), std::ios::binary);
    if(!fIn->is_open()){
      delete fIn;
      continue;
    }
    auto tfn = TString(fn.c_str());
    auto n = tfn.Length();
    tfn.Replace(n-4, 4, ".root");
    printf("fout: %s\n", tfn.Data());
    auto fOut = TFile::Open(TString(Form("%s/",folderName.c_str()))+tfn, "recreate");
    convertTL(fIn, gemroc);
    fIn->close();
    delete fIn;
    fOut->Close();
    delete fOut;
  }
  
}

void tiger_tree_converter_bin(string folderName){
  printf("trigger_tree_converter_bin\n");
  tiger_tree_converter_tl(folderName);
}
