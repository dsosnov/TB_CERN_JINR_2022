#ifndef vmm_h
#define vmm_h

#include "analysisGeneral.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "vector"
using std::vector;

// Header file for the classes stored in the TTree if any.
//#include "c++/v1/vector"

class vmm : public analysisGeneral {
public :
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

   vmm(TString filename) : analysisGeneral(filename),
                           triggerTimeStamp(nullptr),
                           triggerCounter(nullptr),
                           boardId(nullptr),
                           chip(nullptr),
                           eventSize(nullptr),
                           daq_timestamp_s(nullptr),
                           daq_timestamp_ns(nullptr),
                           tdo(nullptr),
                           pdo(nullptr),
                           flag(nullptr),
                           threshold(nullptr),
                           bcid(nullptr),
                           relbcid(nullptr),
                           overflow(nullptr),
                           orbitCount(nullptr),
                           grayDecoded(nullptr),
                           channel(nullptr),
                           febChannel(nullptr),
                           mappedChannel(nullptr),
                           art_valid(nullptr),
                           art(nullptr),
                           art_trigger(nullptr)
    {
      printf("vmm constructor\n");
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(folder + file + ending);
      if (!f || !f->IsOpen()) {
        f = new TFile(folder + file + ending);
      }
      if(!f->IsOpen()) {
        std::cout << "Problem with opening data file" << std::endl;
        exit(1);
      }
      TTree* tree = nullptr; 
      f->GetObject("analysisGeneral",tree);
      if(tree)
        fChain = tree;
      Init();
    };
   vmm(TTree *tree = 0) : triggerTimeStamp(nullptr),
                          triggerCounter(nullptr),
                          boardId(nullptr),
                          chip(nullptr),
                          eventSize(nullptr),
                          daq_timestamp_s(nullptr),
                          daq_timestamp_ns(nullptr),
                          tdo(nullptr),
                          pdo(nullptr),
                          flag(nullptr),
                          threshold(nullptr),
                          bcid(nullptr),
                          relbcid(nullptr),
                          overflow(nullptr),
                          orbitCount(nullptr),
                          grayDecoded(nullptr),
                          channel(nullptr),
                          febChannel(nullptr),
                          mappedChannel(nullptr),
                          art_valid(nullptr),
                          art(nullptr),
                          art_trigger(nullptr),
                          analysisGeneral(tree)
    {
      printf("vmm constructor\n");
      Init();
    };
   virtual ~vmm() {};

   void Init() override;
   virtual void Loop() override;

   map<unsigned int, vector<array<int, 2>>> TDOlimits;
   vector<array<float, 2>> pdoCorrection;

   virtual void addLimits(int minLimit, TString filename);
   array<int, 2> getLimits(int channel, int pdo);
   virtual int getLimitLow(int channel, int pdo);
   virtual int getLimitUp(int channel, int pdo);
   virtual double getTime(int channel, int bcid, int tdo, int pdo);
   double getTimeByHand(int bcid, int tdo, int lowLimit, int upLimit);

   virtual void addPDOCorrection(TString filename);
   virtual int correctPDO(int channel, int pdoIn);

   map<int, pair<int, int>> channelMap;
   virtual void addMap(TString filename);
   virtual pair<int,int> getMapped(int channel);
   virtual int getMappedDetector(int channel);
   virtual int getMappedChannel(int channel);
};

void vmm::addLimits(int minLimit, TString filename){
   vector<array<int, 2>> limitsCurrent;
   ifstream myfile(Form("../out/%s", filename.Data()));
   int bin1 = 0;
   int bin2 = 0;
   int i = 0;
   while (myfile >> bin1 >> bin2)
   {
      limitsCurrent.push_back({bin1, bin2});
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

void vmm::addPDOCorrection(TString filename){
   ifstream myfile(Form("../out/%s", filename.Data()));
   float p0, p1;
   int i = 0;
   uint sizeBefore = pdoCorrection.size();
   while (myfile >> p0 >> p1)
   {
      pdoCorrection.push_back({p0, p1});
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

void vmm::addMap(TString filename){
   ifstream infile(Form("../out/%s", filename.Data()));
   std::string line;
   int ch, d, dch;
   while (std::getline(infile, line))
   {
     std::istringstream iss(line);
     if(iss.str().substr(0, 1) == string("#")) // in c++20 there is starts_with("#")
       continue;
     if (!(iss >> ch >> d >> dch))
       break; // error
     channelMap.emplace(ch, make_pair(d, dch));
   }
}
pair<int,int> vmm::getMapped(int channel){
  if(!channelMap.count(channel))
    return {-1, -1};
  return channelMap.at(channel);
}
int vmm::getMappedDetector(int channel){
  return getMapped(channel).first;
}
int vmm::getMappedChannel(int channel){
  return getMapped(channel).second;
}

void vmm::Init()
{
   printf("vmm::Init\n");

   if (fChain==nullptr){
     printf("fchain == nullptr\n");
     return;
   }
   fCurrent = -1;
   // fChain->SetMakeClass(1);

   // Set branch addresses and branch pointers
   // fChain->SetBranchAddress("eventFAFA", &eventFAFA, &b_eventFAFA);
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

   // addPDOCorrection("calibration_pdo_t@t_g3_p25_s100.txt");
   addPDOCorrection("calibration_pdo_t@t_g1_p25_s100.txt");

   addMap("map.txt");
}

#endif
