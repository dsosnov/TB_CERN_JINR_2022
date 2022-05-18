#pragma once

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class analysisGeneral {
public :
   TString folder = "../data/";
   TString file = "run_0057";
   TString ending = ".root";

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   analysisGeneral(TString);
   analysisGeneral(TTree *tree=0);
   virtual ~analysisGeneral();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     InitChain(TTree *tree) {};
   virtual void     Init() {};
   virtual void     Loop() {};
   virtual TTree* GetTree(TString filename = "", TString treeName = "vmm");

};

TTree* analysisGeneral::GetTree(TString filename, TString treeName){
  if(filename == TString(""))
    filename = file;
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(folder + filename + ending);
  if (!f || !f->IsOpen()) {
    f = new TFile(folder + file + ending);
  }
  if(!f->IsOpen()) {
    std::cout << "Problem with opening data file" << std::endl;
    exit(1);
  }
  TTree* tree = nullptr; 
  f->GetObject(treeName.Data(), tree);
  return tree;
}

analysisGeneral::analysisGeneral(TString filename) : file(filename), fChain(nullptr), fCurrent(-1)
{
  auto tree = GetTree();
  if(tree)
    fChain = tree;
   Init();
}

analysisGeneral::analysisGeneral(TTree *tree) : fChain(nullptr), fCurrent(-1)
{
   if (tree == nullptr) {
     tree = GetTree();
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
