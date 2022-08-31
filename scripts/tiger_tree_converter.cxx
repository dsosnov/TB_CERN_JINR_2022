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

#include <iostream>
#include <fstream>

using std::vector;
using std::string;
using std::map;

using std::ifstream;

void split(const string &s, char delimeter, vector<TString> &result){
    std::istringstream iss(s);
    string item;
    while (std::getline(iss, item, delimeter)){
      result.push_back(TString(item));
    }
}

void convertTL(ifstream* fIn, TFile* fOut, Char_t gemroc){
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
  string line;
  vector<TString> lineTokens;
  while (std::getline(*fIn, line)){
    lineTokens.clear();
    split(line, ' ', lineTokens);
    tigerID = TString::BaseConvert(lineTokens.at(1), 16, 10).Atoi();
    if(lineTokens.at(2) == "EW:"){
      // tiger - #1, channel - #4, tac - #6, tcoarse #8, ecoarse #10, tfine #12, efile #14
      channelID = TString::BaseConvert(lineTokens.at(4), 16, 10).Atoi();
      tacID = TString::BaseConvert(lineTokens.at(6), 16, 10).Atoi();
      tCoarse = TString::BaseConvert(lineTokens.at(8), 16, 10).Atoi();
      eCoarse = TString::BaseConvert(lineTokens.at(10), 16, 10).Atoi();
      tFine = TString::BaseConvert(lineTokens.at(12), 16, 10).Atoi();
      eFine = TString::BaseConvert(lineTokens.at(14), 16, 10).Atoi();
      frameCount = lastFrameCount.at(tigerID);
      seu = lastSEU.at(tigerID);
      frameCountLoops = FCRollOvers.at(tigerID);
      counterWord = lastCW.at(tigerID*64 + channelID);
      tree->Fill();
      // printf("Event Word: %s\n", line.c_str());
    } else if(lineTokens.at(2) == "HB:"){
      frameCount = TString::BaseConvert(lineTokens.at(4), 16, 10).Atoi();
      seu = TString::BaseConvert(lineTokens.at(6), 16, 10).Atoi();
      if(frameCount < lastFrameCount.at(tigerID))
        FCRollOvers.at(tigerID)++;
      lastSEU.at(tigerID) = seu;
      lastFrameCount.at(tigerID) = frameCount;
      // printf("Frame Word: %s\n", line.c_str());
      // tiger - #1, framecount - #4, seu - #6
    } else if(lineTokens.at(2) == "CW:"){ // CounterWord
      channelID = TString::BaseConvert(lineTokens.at(4), 16, 10).Atoi();
      counterWord = TString::BaseConvert(lineTokens.at(6), 16, 10).Atoi();
      lastCW.at(tigerID*64 + channelID) = counterWord;
      // printf("Counter Word: %s\n", line.c_str());
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
    if (name.Index(TRegexp(Form("SubRUN_*_GEMROC_*_TL.datout.txt"),kTRUE)) == kNPOS)
        continue;
    datoutFiles.push_back(name.Data());
    printf("%s\n", name.Data());
  }

  for(auto &fn: datoutFiles){
    int gemroc = static_cast<TObjString*>(TString(fn).Tokenize("_")->At(3))->String().Atoi();
    auto fIn = new ifstream((folderName+"/"+fn).c_str());
    if(!fIn->is_open()){
      delete fIn;
      continue;
    }
    auto tfn = TString(fn.c_str());
    auto n = tfn.Length();
    tfn.Replace(n-11, 11, ".root");
    printf("fout: %s\n", tfn.Data());
    auto fOut = TFile::Open(TString(Form("%s/",folderName.c_str()))+tfn, "recreate");
    convertTL(fIn, fOut, gemroc);
    fIn->close();
    delete fIn;
    fOut->Close();
    delete fOut;
  }
  
}

void tiger_tree_converter(string folderName){
  printf("trigger_tree_converter\n");
  tiger_tree_converter_tl(folderName);
}
