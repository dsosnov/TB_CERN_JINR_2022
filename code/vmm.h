#ifndef vmm_h
#define vmm_h

#include "analysisGeneral.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

#include <vector> 
#include <map>
#include <utility>
#include <string>
#include <memory>
#include <set>

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm> // max_element, sort

using std::vector;
using std::map;
using std::string;
using std::shared_ptr, std::make_shared;
using std::pair, std::make_pair;
using std::tuple, std::get;
using std::set;
using std::array;
using std::begin, std::end;
using std::ifstream;

// Header file for the classes stored in the TTree if any.
//#include "c++/v1/vector"

class vmm : public analysisGeneral {
public :
   TString runType = "g1_p25_s100";
   TString mapFile = "map-20220518.txt";

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

   vmm(TString, TString runType_ = "g1_p25_s100", TString mapFile_ = "map-strawOnly.txt");
   vmm(vector<TString>, TString runType_ = "g1_p25_s100", TString mapFile_ = "map-strawOnly.txt");
   vmm(TChain *tree = nullptr, TString runType_ = "g1_p25_s100", TString mapFile_ = "map-strawOnly.txt");
   virtual ~vmm();
   virtual void     Init() override;
   virtual void     Loop(unsigned long n = 0) override;

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

   map<int, pair<int, int>> channelMap;
   void addMap(TString filename, bool verbose = false);
   pair<int,int> getMapped(int channel);
   int getMappedDetector(int channel);
   int getMappedChannel(int channel);
};

#endif

vmm::vmm(TString filename, TString runType_, TString mapFile_) : runType(runType_), mapFile(mapFile_),
                                                                 triggerTimeStamp(nullptr), triggerCounter(nullptr),
                                                                 boardId(nullptr), chip(nullptr), eventSize(nullptr),
                                                                 daq_timestamp_s(nullptr), daq_timestamp_ns(nullptr),
                                                                 tdo(nullptr), pdo(nullptr), flag(nullptr), threshold(nullptr),
                                                                 bcid(nullptr), relbcid(nullptr), overflow(nullptr), orbitCount(nullptr), grayDecoded(nullptr),
                                                                 channel(nullptr), febChannel(nullptr), mappedChannel(nullptr),
                                                                 art_valid(nullptr), art(nullptr), art_trigger(nullptr)
{
  file = filename;
  fChain = GetTree(filename);
  Init();
}

vmm::vmm(vector<TString> filenames, TString runType_, TString mapFile_) : runType(runType_), mapFile(mapFile_),
                                                                          triggerTimeStamp(nullptr), triggerCounter(nullptr),
                                                                          boardId(nullptr), chip(nullptr), eventSize(nullptr),
                                                                          daq_timestamp_s(nullptr), daq_timestamp_ns(nullptr),
                                                                          tdo(nullptr), pdo(nullptr), flag(nullptr), threshold(nullptr),
                                                                          bcid(nullptr), relbcid(nullptr), overflow(nullptr), orbitCount(nullptr), grayDecoded(nullptr),
                                                                          channel(nullptr), febChannel(nullptr), mappedChannel(nullptr),
                                                                          art_valid(nullptr), art(nullptr), art_trigger(nullptr)
{
  file = filenames.at(0);
  fChain = GetTree(filenames.at(0));
  for(auto i = 1; i < filenames.size(); i++)
    fChain->Add(folder + filenames.at(i) + ending);
  Init();
}

vmm::vmm(TChain *tree, TString runType_, TString mapFile_) : analysisGeneral(tree),  runType(runType_), mapFile(mapFile_),
                                                             triggerTimeStamp(nullptr), triggerCounter(nullptr),
                                                             boardId(nullptr), chip(nullptr), eventSize(nullptr),
                                                             daq_timestamp_s(nullptr), daq_timestamp_ns(nullptr),
                                                             tdo(nullptr), pdo(nullptr), flag(nullptr), threshold(nullptr),
                                                             bcid(nullptr), relbcid(nullptr), overflow(nullptr), orbitCount(nullptr), grayDecoded(nullptr),
                                                             channel(nullptr), febChannel(nullptr), mappedChannel(nullptr),
                                                             art_valid(nullptr), art(nullptr), art_trigger(nullptr)
{
  fChain = (tree == nullptr) ? GetTree("") : tree;
  Init();
}

vmm::~vmm()
{
}

void vmm::addLimits(int minLimit, TString filename, bool verbose){
   vector<array<int, 2>> limitsCurrent;
   ifstream myfile(Form("../configs/%s", filename.Data()));
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
   ifstream myfile(Form("../configs/%s", filename.Data()));
   float p0, p1;
   int i = 0;
   unsigned int sizeBefore = pdoCorrection.size();
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

void vmm::addMap(TString filename, bool verbose){
   ifstream infile(Form("../configs/%s", filename.Data()));
   std::string line;
   int ch, d, dch;
   while (std::getline(infile, line))
   {
     std::istringstream iss(line);
     if(iss.str().substr(0, 1) == string("#")) // in c++20 there is starts_with("#")
       continue;
     if (!(iss >> ch >> d >> dch))
       break; // error
     if(verbose)
       printf("Map: %d: %d - %d\n", ch, d, dch);
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
   // Set branch addresses and branch pointers
   if (!fChain) return;
   printf("vmm::Init()\n");
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

   // addLimits(1, "calibration_tdo_t@t_g1_p25_s100.txt");
   if(runType.BeginsWith("g1_p25_s100")){
     addLimits(100, "calibration_tdo_t@t_g1_p25_s100_pdo100.txt");
     addLimits(150, "calibration_tdo_t@t_g1_p25_s100_pdo150.txt");
     addLimits(200, "calibration_tdo_t@t_g1_p25_s100_pdo200.txt");
     addLimits(250, "calibration_tdo_t@t_g1_p25_s100_pdo250.txt");
     addLimits(300, "calibration_tdo_t@t_g1_p25_s100_pdo300.txt");
   } else if(runType.BeginsWith("g3_p25_s100")){
     addLimits(150, "calibration_tdo_t@t_g3_p25_s100_pdo150.txt");
     addLimits(200, "calibration_tdo_t@t_g3_p25_s100_pdo200.txt");
     addLimits(280, "calibration_tdo_t@t_g3_p25_s100_pdo280.txt");
     addLimits(420, "calibration_tdo_t@t_g3_p25_s100_pdo420.txt");
     addLimits(570, "calibration_tdo_t@t_g3_p25_s100_pdo570.txt");
     addLimits(710, "calibration_tdo_t@t_g3_p25_s100_pdo710.txt");
   }
   if(runType == "g1_p25_s100"){
     addPDOCorrection("calibration_pdo_t@t_g1_p25_s100.txt");
     // addPDOCorrection("calibration_pdo_t@t_g1_p25_s100_sci0&60.txt");
   } else if (runType == "g3_p25_s100"){
     addPDOCorrection("calibration_pdo_t@t_g3_p25_s100.txt");
     // addPDOCorrection("calibration_pdo_t@t_g3_p25_s100_sci0&60.txt");
   } else if(runType == "g1_p25_s100-0&60"){
     addPDOCorrection("calibration_pdo_t@t_g1_p25_s100_sci0&60.txt");
   } else if (runType == "g3_p25_s100-0&60"){
     addPDOCorrection("calibration_pdo_t@t_g3_p25_s100_sci0&60.txt");
   }

   if(!mapFile.EndsWith(".txt"))
     mapFile.Append(".txt");
   addMap(mapFile.Data());
}

#ifndef vmm_cxx
void vmm::Loop(unsigned long n){}
#endif // #ifdef vmm_cxx
