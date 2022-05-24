#pragma once

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

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
   virtual void     Loop() {};
   virtual void     LoopOnce(unsigned long long entry) {};
   virtual TChain* GetTree(TString filename = "", TString treeName = "vmm");

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
