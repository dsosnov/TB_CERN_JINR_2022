#pragma once

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <vector>
#include <string>
#include <map>

using std::vector;
using std::string;
using std::map;

class analysisGeneral {
public :
   TString folder = "../data/";
   TString file = "run_0057";
   TString ending = ".root";

   TChain          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   analysisGeneral(TChain *tree=0);
   analysisGeneral(TString);
   analysisGeneral(vector<TString>);
   virtual ~analysisGeneral();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     InitChain(TTree *tree) {};
   virtual void     Init() {};
   virtual void     Loop(unsigned long n = 0) {};
   virtual void     LoopSecond(unsigned long long) {};
   virtual TChain* GetTree(TString filename = "", TString treeName = "vmm");

   bool syncSignal = false;
   void useSyncSignal(bool use = true) {syncSignal = use;}

   struct mm2CenterHitParameters{
     bool approx, sync;
     unsigned long timeSec;
     unsigned int timeMSec;
     unsigned int stripX, stripY;
     unsigned int pdo;
     double pdoRelative;
     long long nHitsToPrev;
     float time;
     double timeSinceSync;
     long long timeFull() const {
       return timeMSec + timeSec * 1E6;
     }
     string getSignalTypeText() const {
       string signalTypeText = ""; // = "        ";
       if(approx)
         signalTypeText += string(signalTypeText.empty() ? "" : " ") + "a";
       else if(sync)
         signalTypeText += string(signalTypeText.empty() ? "" : " ") + "s";
       if(!signalTypeText.empty())
         signalTypeText = string() + "["+ signalTypeText + "]";
       signalTypeText += string(" ", 5 - signalTypeText.size());
       return signalTypeText;
     }
     void print() const {
       printf("Hit to straw %d with relative pdo %.3f and time %.2f at daq time %lld. %s Previous synchrosignal was %.2f us ago.\n",
              stripX, pdoRelative, time, timeFull(), getSignalTypeText().c_str(), timeSinceSync);
     }
     void printfBrief(bool revert = false) const {
       if(revert)
         printf("%s %3d - %.2f - %.3f - %7lld - %.2g",
                getSignalTypeText().c_str(),
                stripX, time, pdoRelative,
                timeFull() % int(1E7), timeSinceSync);
       else
         printf("%.2g - %7lld - %.3f - %.2f - %3d %s",
                timeSinceSync, timeFull() % int(1E7),
                pdoRelative, time, stripX,
                getSignalTypeText().c_str());
     }

   };
  virtual map<unsigned long, mm2CenterHitParameters> GetCentralHits(unsigned long long fromSec = 0, unsigned long long toSec = 0) {return {};};
};

TChain* analysisGeneral::GetTree(TString filename, TString treeName){
  auto chain = new TChain(treeName);

  if(filename == TString(""))
    filename = file;
  chain->Add(folder + filename + ending);
  
  return chain;
}

analysisGeneral::analysisGeneral(TString filename) : file(filename), fChain(nullptr), fCurrent(-1)
{
  fChain = GetTree(filename);
  Init();
}

analysisGeneral::analysisGeneral(vector<TString> filenames) : file(filenames.at(0)), fChain(nullptr), fCurrent(-1)
{
  fChain = GetTree(filenames.at(0));
  for(auto i = 1; i < filenames.size(); i++)
    fChain->Add(folder + filenames.at(i) + ending);
  Init();
}

analysisGeneral::analysisGeneral(TChain *tree) : fChain(nullptr), fCurrent(-1)
{
   if (tree == nullptr) {
     fChain = GetTree();
   }
   if(tree)
      fChain = tree;
   Init();
}

analysisGeneral::~analysisGeneral()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t analysisGeneral::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t analysisGeneral::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
   }
   return centry;
}
