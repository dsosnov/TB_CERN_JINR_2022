#ifndef vmm_h
#define vmm_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <vector>
#include <fstream>
#include <iostream>
#include <map>
using std::vector;
using std::array;
using std::ifstream, std::ofstream;
using std::map;

// Header file for the classes stored in the TTree if any.
// #include "c++/v1/vector"

class vmm {
public :
   TString folder = "../data/";
   TString file = "run_0058";
   TString ending = ".root";

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   Int_t           eventFAFA;
   vector<int>     *triggerTimeStamp;
   vector<int>     *triggerCounter;
   vector<int>     *boardId;
   vector<int>     *chip;
   vector<int>     *eventSize;
   vector<int>     *daq_timestamp_s;
   vector<int>     *daq_timestamp_ns;
   vector<vector<int> > *tdo;
   vector<vector<int> > *pdo;
   vector<vector<int> > *flag;
   vector<vector<int> > *threshold;
   vector<vector<int> > *bcid;
   vector<vector<int> > *relbcid;
   vector<vector<int> > *overflow;
   vector<vector<int> > *orbitCount;
   vector<vector<int> > *grayDecoded;
   vector<vector<int> > *channel;
   vector<vector<int> > *febChannel;
   vector<vector<int> > *mappedChannel;
   vector<int>     *art_valid;
   vector<int>     *art;
   vector<int>     *art_trigger;

   // List of branches
   TBranch        *b_eventFAFA;   //!
   TBranch        *b_triggerTimeStamp;   //!
   TBranch        *b_triggerCounter;   //!
   TBranch        *b_boardId;   //!
   TBranch        *b_chip;   //!
   TBranch        *b_eventSize;   //!
   TBranch        *b_daq_timestamp_s;   //!
   TBranch        *b_daq_timestamp_ns;   //!
   TBranch        *b_tdo;   //!
   TBranch        *b_pdo;   //!
   TBranch        *b_flag;   //!
   TBranch        *b_threshold;   //!
   TBranch        *b_bcid;   //!
   TBranch        *b_relbcid;   //!
   TBranch        *b_overflow;   //!
   TBranch        *b_orbitCount;   //!
   TBranch        *b_grayDecoded;   //!
   TBranch        *b_channel;   //!
   TBranch        *b_febChannel;   //!
   TBranch        *b_mappedChannel;   //!
   TBranch        *b_art_valid;   //!
   TBranch        *b_art;   //!
   TBranch        *b_art_trigger;   //!

   vmm(TTree *tree=0);
   virtual ~vmm();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();

};

#endif

#ifdef vmm_cxx
vmm::vmm(TTree *tree) : fChain(0) 
{
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(folder + file + ending);
      if (!f || !f->IsOpen()) {
         f = new TFile(folder + file + ending);
      }
      f->GetObject("vmm",tree);

   }
   Init(tree);
}

vmm::~vmm()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t vmm::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t vmm::LoadTree(Long64_t entry)
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

void vmm::Init(TTree *tree)
{
   triggerTimeStamp = 0;
   triggerCounter = 0;
   boardId = 0;
   chip = 0;
   eventSize = 0;
   daq_timestamp_s = 0;
   daq_timestamp_ns = 0;
   tdo = 0;
   pdo = 0;
   flag = 0;
   threshold = 0;
   bcid = 0;
   relbcid = 0;
   overflow = 0;
   orbitCount = 0;
   grayDecoded = 0;
   channel = 0;
   febChannel = 0;
   mappedChannel = 0;
   art_valid = 0;
   art = 0;
   art_trigger = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventFAFA", &eventFAFA, &b_eventFAFA);
   fChain->SetBranchAddress("triggerTimeStamp", &triggerTimeStamp, &b_triggerTimeStamp);
   fChain->SetBranchAddress("triggerCounter", &triggerCounter, &b_triggerCounter);
   fChain->SetBranchAddress("boardId", &boardId, &b_boardId);
   fChain->SetBranchAddress("chip", &chip, &b_chip);
   fChain->SetBranchAddress("eventSize", &eventSize, &b_eventSize);
   fChain->SetBranchAddress("daq_timestamp_s", &daq_timestamp_s, &b_daq_timestamp_s);
   fChain->SetBranchAddress("daq_timestamp_ns", &daq_timestamp_ns, &b_daq_timestamp_ns);
   fChain->SetBranchAddress("tdo", &tdo, &b_tdo);
   fChain->SetBranchAddress("pdo", &pdo, &b_pdo);
   fChain->SetBranchAddress("flag", &flag, &b_flag);
   fChain->SetBranchAddress("threshold", &threshold, &b_threshold);
   fChain->SetBranchAddress("bcid", &bcid, &b_bcid);
   fChain->SetBranchAddress("relbcid", &relbcid, &b_relbcid);
   fChain->SetBranchAddress("overflow", &overflow, &b_overflow);
   fChain->SetBranchAddress("orbitCount", &orbitCount, &b_orbitCount);
   fChain->SetBranchAddress("grayDecoded", &grayDecoded, &b_grayDecoded);
   fChain->SetBranchAddress("channel", &channel, &b_channel);
   fChain->SetBranchAddress("febChannel", &febChannel, &b_febChannel);
   fChain->SetBranchAddress("mappedChannel", &mappedChannel, &b_mappedChannel);
   fChain->SetBranchAddress("art_valid", &art_valid, &b_art_valid);
   fChain->SetBranchAddress("art", &art, &b_art);
   fChain->SetBranchAddress("art_trigger", &art_trigger, &b_art_trigger);
}

#endif // #ifdef vmm_cxx
