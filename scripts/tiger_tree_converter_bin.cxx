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
#include <optional>
using std::vector;
using std::string;
using std::map;
using std::optional, std::nullopt;

#include <fstream>
using std::ifstream;

constexpr bool DEBUG_PRINT = false;

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
// union works predictably due to equal size of structs
union TLDataFrame{
  uint64_t unknown;
  struct {
    uint64_t seuCount: 15,
      frameCount: 16,
      reserved: 25,
      tiger: 3,
      key: 5; // should be 0x4
  } frameWord;
  struct {
    uint64_t eFine: 10,
      tFine: 10,
      eCoarse: 10,
      tCoarse: 16,
      tac: 2,
      channel: 6,
      reserved: 2,
      tiger: 3,
      key: 5; // should be 0x0
  } eventWord;
  struct {
    uint64_t counter: 24,
      channel: 6,
      reserved: 26,
      tiger: 3,
      key: 5; // should be 0x8
  } countWord;
};
enum class TLDataFrameType{unknown, FrameWord, EventWord, CountWord};
TLDataFrameType getTypeTL(const TLDataFrame &hit){
  if(hit.frameWord.key == 0x4)
    return TLDataFrameType::FrameWord;
  else if(hit.eventWord.key == 0x0)
    return TLDataFrameType::EventWord;
  else if(hit.countWord.key == 0x8)
    return TLDataFrameType::CountWord;
  return TLDataFrameType::unknown;
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

  map<Int_t, Long64_t> FCRollOvers;
  map<Int_t, optional<Int_t>> lastFrameCount;
  map<Int_t, Int_t> lastSEU;
  map<pair<Int_t, Int_t>, Int_t> lastCW;
  for(auto i = 0; i < 256; i++){
    FCRollOvers.emplace(i, 0);
    lastFrameCount.emplace(i, nullopt);
    lastSEU.emplace(i, 0);
    for(auto j = 0; j < 64; j++)
      lastCW.emplace(make_pair(i, j), 0);
  }
  TLDataFrame frame = {};
  for(fIn->read(reinterpret_cast<char*>(&frame), sizeof(frame));
      !fIn->eof();
      fIn->read(reinterpret_cast<char*>(&frame), sizeof(frame))){
    switch(getTypeTL(frame)){
      case TLDataFrameType::EventWord: { // EW
        tigerID = frame.eventWord.tiger;
        if(!lastFrameCount.at(tigerID))
          continue;
        channelID = frame.eventWord.channel;
        tacID = frame.eventWord.tac;
        tCoarse = frame.eventWord.tCoarse;
        eCoarse = frame.eventWord.eCoarse;
        tFine = frame.eventWord.tFine;
        eFine = frame.eventWord.eFine;
        frameCount = lastFrameCount.at(tigerID).value();
        seu = lastSEU.at(tigerID);
        frameCountLoops = FCRollOvers.at(tigerID);
        counterWord = lastCW.at({tigerID, channelID});
        tree->Fill();
        if(DEBUG_PRINT)
          printf("TIGER %01X: EW: ChID: %02X tacID: %01X Tcoarse: %04X Ecoarse: %03X Tfine: %03X Efine: %03X \n",
                 tigerID, channelID, tacID, tCoarse, eCoarse, tFine, eFine);
        break;
      }
      case TLDataFrameType::FrameWord: { // HB
        tigerID = frame.frameWord.tiger;
        frameCount = frame.frameWord.frameCount;
        seu = frame.frameWord.seuCount;
        if(lastFrameCount.at(tigerID) && frameCount < lastFrameCount.at(tigerID).value())
          FCRollOvers[tigerID]++;
        lastSEU[tigerID] = seu;
        lastFrameCount[tigerID] = frameCount;
        if(DEBUG_PRINT)
          printf("TIGER %01X: HB: Framecount: %08X SEUcount: %08X\n", tigerID, frameCount, seu);
        break;
      }
      case TLDataFrameType::CountWord: { // CW
        tigerID = frame.countWord.tiger;
        channelID = frame.countWord.channel;
        counterWord = frame.countWord.counter;
        lastCW[{tigerID, channelID}] = counterWord;
        if(DEBUG_PRINT)
          printf("TIGER %01X: CW: ChID: %02X CounterWord: %016X\n", tigerID, channelID, counterWord);
        break;
      }
      case TLDataFrameType::unknown: {
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

/* pragma pack desables padding */
#pragma pack(push, 1)
union TMDataFrame{
  uint64_t unknown;
  struct {
    uint16_t l1Timestamp: 16;
    uint8_t countHits: 8; // What happens if all 64*8 = 512 channels will be activated?

    // in python converter:
    uint32_t l1LocalCount: 32;
    uint8_t reserved: 2;
    // in tiger_data_format.pdf:
    // uint64_t l1LocalCount: 34;

    uint8_t status: 3;
    uint8_t key: 3; // should be 0x6
  } header;
  struct {
    uint32_t lastCountWordData: 18;
    uint8_t lastCountWordCh: 6;
    uint8_t l1LocalCount: 3;
    uint8_t tiger: 3;
    uint8_t reserved: 2;
    uint8_t gemroc: 5;
    uint32_t l1LocalFramenum: 24;
    uint8_t key: 3; // should be 0x7
  } trailer;
  struct {
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
  } data;
  struct {
    uint32_t reserved: 28;
    uint32_t udpFrameCount: 24;
    uint8_t gemroc: 5;
    uint32_t headerStatus: 3;
    uint8_t key: 4; // should be 0x4
  } udpCounter;
};
#pragma pack(pop)
enum class TMDataFrameType{unknown, Header, Trailer, Data, UDPCounter};
TMDataFrameType getTypeTM(const TMDataFrame &hit){
  if(hit.trailer.key == 0x7)
    return TMDataFrameType::Trailer;
  else if(hit.header.key == 0x6)
    return TMDataFrameType::Header;
  else if(hit.data.key == 0x0)
    return TMDataFrameType::Data;
  else if(hit.udpCounter.key == 0x4)
    return TMDataFrameType::UDPCounter;
  return TMDataFrameType::unknown;
}

string printBinary(const TMDataFrame &hit){
  string s = "";
  for(auto i = 63; i >= 0; i--){
    s += static_cast<bool>((hit >> i) & 0x1) ? "1" : "0";
  }
  return s;
}

void convertTM(ifstream* fIn, Char_t gemroc){
  static vector<Char_t> gemrocID = gemroc; // -- "B" = Char_t == int8_t
  static vector<Short_t> tigerID; // 8 bit -- "S" == Short_t == int16_t
  static vector<Char_t> channelID; // 6 bit -- "B" = Char_t == int8_t
  static vector<Char_t> tacID; // 2 bit -- 4 TAC per shaper for event de-randomization -- "B" = Char_t == int8_t
  static vector<Int_t> tCoarse; // 16 bit data -- "I" = Int_t == int32_t
  static vector<Short_t> eCoarse; // 10 bit data -- "S" == Short_t == int16_t
  static vector<Short_t> tFine; // 10 bit data -- "S" == Short_t == int16_t
  static vector<Short_t> eFine; // 10 bit data -- "S" == Short_t == int16_t
  static vector<Char_t> lastFrameNum; // 3 bit -- "B" = Char_t == int8_t
  // Int_t frameCount = 0; // 16 bit data -- "I" = Int_t == int32_t
  // Int_t seu = 0; // 16 bit data -- "I" = Int_t == int32_t
  // Long64_t frameCountLoops = 0; // -- "L" = Long64_t == int64_t
  // Int_t counterWord = 0; // 24 bit data -- "I" = Int_t == int32_t
  Long64_t l1LocalCount = 0; // 32 bit data -- "L" = Long64_t == int64_t
  Int_t l1Timestamp = 0; // 16 bit data -- "I" = Int_t == int32_t
  Char_t status = 0; // 3 bit data -- "B" = Char_t == int8_t
  int hitsLeft = 0;

  auto tree = new TTree("tigerTM", "tigerTM");
  tree->Branch("gemrocID", &gemrocID);
  tree->Branch("tigerID", &tigerID);
  tree->Branch("channelID", &channelID);
  tree->Branch("tacID", &tacID);
  tree->Branch("tCoarse", &tCoarse);
  tree->Branch("eCoarse", &eCoarse);
  tree->Branch("tFine", &tFine);
  tree->Branch("eFine", &eFine);
  tree->Branch("lastFrameNum", &lastFrameNum);
  // tree->Branch("frameCount", &frameCount, "frameCount/I");
  // tree->Branch("seu", &seu, "seu/I");
  // tree->Branch("frameCountLoops", &frameCountLoops, "frameCountLoops/L");
  tree->Branch("l1LocalCount", &l1LocalCount, "l1LocalCount/L");
  tree->Branch("l1Timestamp", &l1Timestamp, "l1Timestamp/I");
  tree->Branch("status", &status, "status/B");
  
  Int_t prevl1Timestamp = 0;

  struct {bool header = false, trailer = false, udp = false, errors = false;} hitStatus;
  TMDataFrame frame = {};
  for(fIn->read(reinterpret_cast<char*>(&frame), sizeof(frame));
      !fIn->eof();
      fIn->read(reinterpret_cast<char*>(&frame), sizeof(frame))){
    switch(getTypeTM(frame)){
      case TMDataFrameType::Header: { // HEADER
        if(hitStatus.header)
          printf("Event %d was not finished normally!\n", l1LocalCount);
        {
          gemrocID.clear();
          tigerID.clear();
          channelID.clear();
          tacID.clear();
          tCoarse.clear();
          eCoarse.clear();
          tFine.clear();
          eFine.clear();
          lastFrameNum.clear();
        }
        hitStatus.header = true;
        hitStatus.trailer = false;
        hitStatus.udp = false;
        hitStatus.errors = false;
        l1LocalCount = frame.Header.l1LocalCount;
        l1Timestamp = frame.Header.l1Timestamp;
        status = frame.Header.status;
        hitsLeft = frame.Header.countHits;
        auto timeDiff = l1Timestamp - prevl1Timestamp;
        timediff += (timediff > 0) ? 0 : (1<<16);
        if(DEBUG_PRINT)
          printf("%s  HEADER :  STATUS BIT[2:0]: %01X: LOCAL L1 COUNT: %08X HitCount: %02X LOCAL L1 TIMESTAMP: %04X; Diff w.r.t. previous L1_TS: %04f us\n",
                 printBinary(hit).c_str(),
                 status, l1LocalCount, hitsLeft, l1Timestamp, timediff * 25.0 / 1E3);
        prevl1Timestamp = l1Timestamp;
        break;
      }
        //Trailer, Data, UDPCounter
      case TMDataFrameType::Trailer: { // TRAILER -- will ignore for now
        if(!hitStatus.header || hitStatus.udp)
          continue;
        if(hitStatus.trailer)
          hitStatus.errors = true;
        hitStatus.trailer = true;
        if(DEBUG_PRINT)
          printf("%s  TRAILER: LOCAL L1  FRAMENUM [23:0]: %06X: GEMROC_ID: %02X TIGER_ID: %01X LOCAL L1 COUNT[2:0]: %01X LAST COUNT WORD FROM TIGER:CH_ID[5:0]: %02X LAST COUNT WORD FROM TIGER: DATA[17:0]: %05X \n",
                 printBinary(hit).c_str(),
                 frame.Trailer.l1LocalFramenum, frame.Trailer.gemroc, frame.Trailer.tiger, frame.Trailer.l1LocalCount, frame.Trailer.lastCountWordCh, frame.Trailer.lastCountWordData);
        break;
      }
      case TMDataFrameType::Data: { // DATA
        if(!hitStatus.header || hitStatus.udp) // TODO maybe also check for trailer?
          continue;
        if(hitsLeft <= 0)
          hitStatus.errors = true;

        eFine.push_back(frame.data.eFine);
        tFine.push_back(frame.data.tFine);
        eCoarse.push_back(frame.data.eCoarse);
        tCoarse.push_back(frame.data.tCoarse);
        tacID.push_back(frame.data.tac);
        channelID.push_back(frame.data.channel);
        tigerID.push_back(frame.data.tiger);
        lastFrameNum.push_back(frame.data.lastTigerFrameNumber);

        if(DEBUG_PRINT)
          printf("%s  DATA   : TIGER: %01X L1_TS - TIGERCOARSE_TS: %d LAST TIGER FRAME NUM[2:0]: %01X TIGER DATA: ChID [base10]: %d tacID: %01X Tcoarse: %04X Ecoarse: %03X Tfine: %03X Efine: %d \n",
                 printBinary(hit).c_str(),
                 tigerID.back(), l1LocalCount - tCoarse.back(), lastFrameNum.back(), channelID.back(), tacID.back(), tCoarse.back(), eCoarse.back(), tFine.back(), eFine.back());
        
        hitsLeft--;
        break;
      }
      case TMDataFrameType::UDPCounter: { // UDP_SEQNO
        if(hitStatus.udp)
          hitStatus.errors = true;
        hitStatus.udp = true;
        if(frame.udp.udpFrameCount != l1LocalCount + 1)
          printf("UDP Countel :: For event %d UDP frame count not the %d\n", l1LocalCount, l1LocalCount+1);
        if(DEBUG_PRINT)
          printf("%s  UDP_SEQNO: GEMROC_ID: %02X UDP_SEQNO_U48: %012X  STATUS BIT[5:3]:{}\n",
                 printBinary(hit).c_str(),
                 frame.udp.gemroc, frame.udp.udpFrameCount, frame.udp.headerStatus);
        break;
      }
      case TLDataFrameType::unknown: {
        break;
      }
    }
    if(hitStatus.header && hitStatus.udp){
      if(hitStatus.errors)
        printf("For event %d errors was found!\n", l1LocalCount);
      else if(!hitStatus.trailer)
        printf("For event %d trailes is missed!\n", l1LocalCount);
      else{
        tree->Fill();
        
        gemrocID.clear();
        tigerID.clear();
        channelID.clear();
        tacID.clear();
        tCoarse.clear();
        eCoarse.clear();
        tFine.clear();
        eFine.clear();
        lastFrameNum.clear();

        hitStatus.header = false;
        hitStatus.trailer = false;
        hitStatus.udp = false;
        hitStatus.errors = false;
      }
    }
  }
  tree->Write();
  delete tree;
}

void tiger_tree_converter_bin(string folderName){
  printf("trigger_tree_converter_bin\n");
  tiger_tree_converter_tl(folderName);
}
