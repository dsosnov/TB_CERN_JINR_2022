#ifndef vmm_h
#define vmm_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
//#include "c++/v1/vector"

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

   vmm(TString);
   vmm(TTree *tree=0);
   virtual ~vmm();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();

   map<unsigned int, vector<array<int, 2>>> TDOlimits;
   vector<array<float, 2>> pdoCorrection;

  void addLimits(int minLimit, TString filename, bool verbose = false);
   array<int, 2> getLimits(int channel, int pdo);
   int getLimitLow(int channel, int pdo);
   int getLimitUp(int channel, int pdo);
   double getTime(int channel, int bcid, int tdo, int pdo);
   static double getTimeByHand(int bcid, int tdo, int lowLimit, int upLimit);

   void addPDOCorrection(TString filename, bool verbose = false);
   int correctPDO(int channel, int pdoIn);
};

#endif

#ifdef vmm_cxx
vmm::vmm(TString filename) : file(filename)
{
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(folder + file + ending);
   if (!f || !f->IsOpen()) {
      f = new TFile(folder + file + ending);
   }
   if(!f->IsOpen()) {
      std::cout << "Problem with opening data file" << std::endl;
      exit(1);
   }
   TTree* tree = nullptr; 
   f->GetObject("vmm",tree);
   Init(tree);
}

vmm::vmm(TTree *tree) : fChain(0) 
{
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(folder + file + ending);
      if (!f || !f->IsOpen()) {
         f = new TFile(folder + file + ending);
      }
      if(!f->IsOpen()) {
         std::cout << "Problem with opening data file" << std::endl;
         exit(1);
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

void vmm::addLimits(int minLimit, TString filename, bool verbose){
   vector<array<int, 2>> limitsCurrent;
   ifstream myfile(Form("../out/%s", filename.Data()));
   int bin1 = 0;
   int bin2 = 0;
   int i = 0;
   while (myfile >> bin1 >> bin2)
   {
      limitsCurrent.push_back({bin1, bin2});
      if(verbose)
        std::cout << "MinPDO "<< minLimit << "\tCH " << i << "\t L TDO: " << bin1 << "\t R TDO: " << bin2 << "\n";
      i++;
   }
   if(!limitsCurrent.size()){
      std::cout << "Problem with TDO calibration file " << filename << std::endl;
      exit(1);
   }
   TDOlimits.emplace(minLimit, limitsCurrent);
}
array<int, 2> vmm::getLimits(int channel, int pdo){
   int maxPDOLower = -1;
   int minPDO = -1;
   for(auto &l: TDOlimits){
     if((maxPDOLower < 0 || l.first > maxPDOLower) && l.first < pdo){
        maxPDOLower = l.first;
     }
     if(minPDO < 0 || l.first < minPDO)
        minPDO = l.first;
   }
   maxPDOLower = (maxPDOLower < 0) ? minPDO : maxPDOLower;
   return TDOlimits.at(maxPDOLower)[channel];
}
int vmm::getLimitLow(int channel, int pdo){
   return getLimits(channel, pdo)[0];
}
int vmm::getLimitUp(int channel, int pdo){
   return getLimits(channel, pdo)[1];
}
double vmm::getTime(int channel, int bcid, int tdo, int pdo){
   auto TDOlimits = getLimits(channel, pdo);
   auto ll = getLimitLow(channel, pdo),
        ul = getLimitUp(channel, pdo);
   if(ll < 0 || ul < 0){
      ul = 256;
      ll = 0;
   }
   return getTimeByHand(bcid, tdo, ll, ul);
}
double vmm::getTimeByHand(int bcid, int tdo, int lowLimit, int upLimit){
   return bcid * 25.0 - (tdo - lowLimit) * 25.0 / (upLimit - lowLimit);
}

void vmm::addPDOCorrection(TString filename, bool verbose){
   ifstream myfile(Form("../out/%s", filename.Data()));
   float p0, p1;
   int i = 0;
   uint sizeBefore = pdoCorrection.size();
   while (myfile >> p0 >> p1)
   {
      pdoCorrection.push_back({p0, p1});
      if(verbose)
        std::cout << "PDO Corrections: CH " << i << "\t [0]: " << p0 << "\t [1]: " << p1 << "\n";
      i++;
   }
   if(sizeBefore == pdoCorrection.size()){
      std::cout << "Problem with PDO calibration file " << filename << std::endl;
      exit(1);
   }
}
int vmm::correctPDO(int channel, int pdoIn){
  auto p0 = pdoCorrection.at(channel)[0];
  auto p1 = pdoCorrection.at(channel)[1];
  return static_cast<int>(p0 + pdoIn * p1);
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

   // ================================== LIMITS SEARCH ================================== //or get from file
   // addLimits(0, "calibration_25_100");

   addLimits(100, "calibration_25_100_pdo100.txt");
   addLimits(150, "calibration_25_100_pdo150.txt");
   addLimits(200, "calibration_25_100_pdo200.txt");
   addLimits(250, "calibration_25_100_pdo250.txt");
   addLimits(300, "calibration_25_100_pdo300.txt");

   addPDOCorrection("calibration_pdo_t@t_g1_p25_s100.txt");
}

#endif // #ifdef vmm_cxx
