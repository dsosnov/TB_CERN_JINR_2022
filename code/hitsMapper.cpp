#include "apv.C"
#include "evBuilder.C"
#include <iostream>
#include <fstream>

#include <vector>
#include <utility>
#include <tuple>
#include <map>
#include <set>
#include <optional>

using std::cout, std::endl;
using std::ifstream, std::ofstream;

using std::vector;
using std::pair, std::make_pair;
using std::tuple, std::get;
using std::map;
using std::set;
using std::optional, std::nullopt;

optional<int> calculateVMMNPulsers(int bcidDiff, int maxDiff = 2, int maxNPulsers = 9)
{
    bcidDiff += (bcidDiff > 0) ? 0 : 4096;
    int predicted;
    for(auto i = 1; i <= maxNPulsers; i++)
    {
      predicted = (2000 * i) % 4096;
      if(abs(predicted - bcidDiff) <= maxDiff)
        return i;
    }
    return nullopt;
}

// IMPORTANT change: fixSRSTime
// return: time, us
double apvTimeTimeSRSTimestamp(int diff, bool fixSRSTime = false){
    return diff * 25.0 * (fixSRSTime ? 400000.0/400037.0 : 1.0) / 1000.0;  // TODO add dependency on testbeam
}

int apv_time_from_SRS(int srs1, int srs2, bool fixSRSTime = false)
{
    int diff = 0;
    int srsT_period = 16777215;
    if (srs2 - srs1 >= 0)
    {
        diff = srs2 - srs1;
    }
    else
    {
        diff = srs2 + srsT_period - srs1;
    }

    // Since standard time between SRSTimeStamps is 400038 (400037), but not 400000, correct srsTimeStamp clock time
    // return round(diff * 25.0 * (fixSRSTime ? 400000.0/400037.0 : 1.0) / 1000.0);
    return round(apvTimeTimeSRSTimestamp(diff, fixSRSTime));
}

int vmmRemoveFirstNPuslers(vector<pair<unsigned long, analysisGeneral::mm2CenterHitParameters>> &hits_vmm_v, int nRemoved = 637){
    int firstSyncBcid = 0;
    int index = -nRemoved; // 650; 665

    unsigned long j = 0;
    for (j = 0; j < hits_vmm_v.size() && index < 0; j++)
    {
        if (hits_vmm_v.at(j).second.sync)
            index++;
    }
    hits_vmm_v.erase(hits_vmm_v.begin(), hits_vmm_v.begin() + j-1); // remove first j elements: from 0 to j-1
    if (index == 0) // check if arways?
    {
        firstSyncBcid = hits_vmm_v.at(0).second.bcid;
        hits_vmm_v.erase(hits_vmm_v.begin(), hits_vmm_v.begin()); // remove first elements
    }
    else if(index < 0)
    {
        cout << "Removed all VMM events due to incorrect index value" << endl;
    }
    return firstSyncBcid;
}
pair<int, int> vmmRemoveFirstNPuslers(TTree* hits_vmm_t, pair<unsigned long, analysisGeneral::mm2CenterHitParameters> *hits_vmm_event, int nRemoved = 637){
    int firstSyncBcid = 0;
    int index = -nRemoved; // 650; 665

    unsigned long j = 0;
    for (j = 0; j < hits_vmm_t->GetEntries() && index < 0; j++)
    {
        hits_vmm_t->GetEntry(j);
        if (hits_vmm_event->second.sync)
            index++;
    }
    if (index == 0) // check if arways?
    {
        hits_vmm_t->GetEntry(j);
        firstSyncBcid = hits_vmm_event->second.bcid;
        j++;
        printf("First VMM sync-signal: %lu\n", hits_vmm_event->first);
    }
    else if(index < 0)
    {
        cout << "Removed all VMM events due to incorrect index value" << endl;
    }
    return make_pair(j, firstSyncBcid);
}

bool freeMemory(map<long long, vector<pair<unsigned long, analysisGeneral::mm2CenterHitParameters>>> &hits_vmm_events_map, long long firstElement)
{
    bool released = false;
    vector<long long> keysToErase;
    for(auto &i: hits_vmm_events_map)
        if(i.first < firstElement)
            keysToErase.push_back(i.first);
    if(!keysToErase.empty())
        released = true;
    for(auto &i: keysToErase)
        hits_vmm_events_map.erase(i);
    return released;
}

long long loadNextVMM(long long firstElement, map<long long, vector<pair<unsigned long, analysisGeneral::mm2CenterHitParameters>>> &hits_vmm_events_map,
                 TTree* hits_vmm_t, pair<unsigned long, analysisGeneral::mm2CenterHitParameters> *hits_vmm_event,
                 long long nElements = 100000){
    // printf("loadNextVMM(%lld, %p, %p, %p, %lld)\n", firstElement, &hits_vmm_events_map, hits_vmm_t, hits_vmm_event, nElements);
    if(hits_vmm_events_map.count(firstElement))
        return hits_vmm_events_map.at(firstElement).size();
  
    vector<pair<unsigned long, analysisGeneral::mm2CenterHitParameters>> elementVector;
    hits_vmm_t->GetEntry(firstElement);
    if(!(hits_vmm_event->second.sync)){ // Not a sync event, it seems a problem!
        printf("First loaded event not the sync!\n");
        hits_vmm_t->GetEntry(firstElement-1);
        if(hits_vmm_event->second.sync)
            printf("First before First loaded is sync!\n");
        hits_vmm_events_map.emplace(firstElement, elementVector);
        return 0;
    };
    // hits_vmm_event->second.print();
    long long lastSyncIndex = 0;
    auto maxEntries = hits_vmm_t->GetEntries() - firstElement;
    for(auto i = 0; i < maxEntries && lastSyncIndex < nElements; i++){
        hits_vmm_t->GetEntry(firstElement + i);
        elementVector.push_back(make_pair(hits_vmm_event->first, hits_vmm_event->second));
        if(hits_vmm_event->second.sync)
            lastSyncIndex = i;
    }
    // remove last sync and all after
    elementVector.erase(elementVector.begin() + lastSyncIndex, elementVector.end());
    hits_vmm_events_map.emplace(firstElement, elementVector);
    return elementVector.size(); // next start index
}
long long loadNextVMM(long long firstElement,
                      map<long long, vector<pair<unsigned long, analysisGeneral::mm2CenterHitParameters>>> &hits_vmm_events_map,
                      evBuilder* vmman,
                      long long nElements = 100000){
    if(hits_vmm_events_map.count(firstElement))
        return hits_vmm_events_map.at(firstElement).size();
    // printf("loadNextVMM(%lld, %p, %lld)\n", firstElement, &hits_vmm_events_map, nElements);
  
    vector<pair<unsigned long, analysisGeneral::mm2CenterHitParameters>> elementVector;
    auto hit_vmm = vmman->GetCentralHitsData(firstElement);
    if(!(hit_vmm.sync)){ // Not a sync event, it seems a problem!
        printf("First loaded event not the sync!\n");
        auto hit_vmm_prev = vmman->GetCentralHitsData(firstElement-1);
        if(hit_vmm_prev.sync)
            printf("First before First loaded is sync!\n");
        hits_vmm_events_map.emplace(firstElement, elementVector);
        return 0;
    };
    // hits_vmm_event->second.print();
    long long lastSyncIndex = 0;
    auto maxEntries = vmman->GetEntries() - firstElement;
    for(auto i = 0; i < maxEntries && lastSyncIndex < nElements; i++){
        hit_vmm = vmman->GetCentralHitsData(firstElement + i);
        elementVector.push_back(make_pair(firstElement+i, hit_vmm));
        if(hit_vmm.sync)
            lastSyncIndex = i;
    }
    // remove last sync and all after
    elementVector.erase(elementVector.begin() + lastSyncIndex, elementVector.end());
    hits_vmm_events_map.emplace(firstElement, elementVector);
    return elementVector.size(); // next start index
}

double weightedMean(vector<pair<int, int>> hits){
    long long sum = 0, sumWeights = 0, nHits = 0;
    for(auto &h: hits){
        long long strip = h.first;
        long long pdo = h.second;
        sum += strip * pdo;
        sumWeights += pdo;
        nHits++;
    }
    double center = -1;
    if(nHits)
        center = static_cast<double>(sum) / static_cast<double>(sumWeights);
    return center;
}

/*
 * positions (april measurements):
 * L0 - L1: 285
 * L1 - L2: 345
 * L2 - Straw: 523
 *
 * Some distances (July measurements):
 * 183 - mm0-straw (layer -1 now)
 * 434 - mm0-mm1 - mm layer 3 (y) was connected to mm layer 1
 * 325 - mm3-mm2
 *
 * return: position in mm, from Layer 0
 */
int getLayerPosition(int layer, bool aprilTB = false){ // TODO do normally
  int y = 0;
  if(aprilTB){
    switch(layer){
      case -1:
        y = -183;
        break;
      case 2:
        y += 325;
      case 3:
      case 1:
        y += 434;
      case 0:
        y += 0;
        break;
      default:
        y += 0;
        break;
    }
  }else{
    switch(layer){
      case 3:
        y += 523;
      case 2:
        y += 345;
      case 1:
        y += 285;
      case 0:
        y += 0;
    }
  }
  return y;
}

pair<double, double> getEstimatedTrack(map<int, double> positions){
  int n = 0;
  double sumXY = 0, sumX = 0, sumY = 0, sumX2 = 0, sumY2 = 0;
  int y; double x; // x: horizontal, y - vertical coordinate
  for(auto &pos: positions){
    y = getLayerPosition(pos.first);
    x = pos.second;
    sumXY += x * y;
    sumX += x;
    sumY += y;
    sumX2 += x * x;
    // sumY2 += y * y;
    n++;
  }
  double a0 = (n * sumXY - sumX * sumY)/(n * sumX2 - sumX * sumX);
  double b0 = (sumY - a0*sumX)/n;
  // double b1 = (sumX2 * sumY - sumXY * sumX) / (n * sumX2 - sumX * sumX);
  // double a1 = (sumY - n*b1) / sumX;
  // printf("%d layers -- a0: %g, b0: %g; a1: %g, b1: %g;\n", n, a0, b0, a1, b1);
  return {a0, b0};
}

double estimatePositionInLayer(pair<double, double> trackAB, int layer){
  double y = getLayerPosition(layer);
  double x = (y - trackAB.second) / trackAB.first;
  return x;
}

map<pair<string, string>, pair<int, int>> firstPulserMap = {
    {{"run_0832", "run423"}, {180059, 7042}}, // First APV vs VMM dt = 23209318 or 2320 pulses
    {{"run_0832_cut", "run423_cut"}, {180059, 7042}}, // First APV vs VMM dt = 23209318 or 2320 pulses
    {{"run_0832_cut10m", "run423_cut10m"}, {180059, 7042}}, // First APV vs VMM dt = 23209318 or 2320 pulses
    {{"run_0831", "run422"}, {151740, 7675}}, // First APV vs VMM dt = 32988 or 3 pulses
    {{"run_0830", "run421"}, {239786, 7633}}, // First APV vs VMM dt = 17748110 or 1774 pulses
    {{"run_0828", "run420"}, {138863, 3669}}, // First APV vs VMM dt = 32934 or 3 pulses
    {{"run_0827", "run419"}, {254990, 7638}}, // First APV vs VMM dt = 18284390 or 1828 pulses -- Strange file (in firstPulserMapAlternative)
    {{"run_0826", "run418"}, {80142, 3928}}, // First APV vs VMM dt = 32524 or 3 pulses
    // {{"run_0825", "run417"} //No vmm pulses
    {{"run_0824", "run416"}, {166310, 3429}}, // First APV vs VMM dt = 1477467 or 147 pulses -- Strange file (in firstPulserMapAlternative)
    {{"run_0821", "run412"}, {245637, 7662}}, // First APV vs VMM dt = 23988401 or 2398 pulses
    {{"run_0818", "run411"}, {141763, 3709}}, // First APV vs VMM dt = 26874351 or 2687 pulses
    // {{"run_0817", "run410"}, {53049, 126}}, // First APV vs VMM dt = 33004 or 3 pulses -- during acces
    // {{"run_0816", "run409"} // No apv pulses; No vmm pulses
    {{"run_0815", "run408"}, {238541, 7439}}, // First APV vs VMM dt = 25239911 or 2523 pulses
    {{"run_0814", "run407"}, {227778, 7572}}, // First APV vs VMM dt = 33027 or 3 pulses ; gaps between first n pulses in APV
    // {{"run_0813", "run406"} // No apv pulses
    // {{"run_0812", "run405"} // No apv pulses
    // {{"run_0809", "run399"} // No apv pulses
    // {{"run_0808", "run398"} // No apv pulses; No vmm pulses
    // {{"run_0807", "run397"} // No apv pulses
    {{"run_0805", "run395"}, {202993, 7800}}, // First APV vs VMM dt = 32394 or 3 pulses
    {{"run_0802", "run394"}, {218605, 9372}}, // First APV vs VMM dt = 32976 or 3 pulses
    {{"run_0801", "run393"}, {229541, 11715}}, // First APV vs VMM dt = 33058 or 3 pulses
    // {{"run_0800", "run392"} // No apv pulses
    // {{"run_0797", "run388"} // No apv pulses
    {{"run_0794", "run385"}, {330200, 11716}}, // First APV vs VMM dt = 41550608 or 4155 pulses
    {{"run_0792", "run383"}, {265274, 11227}}, // First APV vs VMM dt = 18083779 or 1808 pulses
    // {{"run_0791", "run382"} // First APV vs VMM dt = -5158615 or -515 pulses; *(VMM is earlier than APV)
    // {{"run_0789", "run380"} // First APV vs VMM dt = -6584644 or -658 pulses *(VMM is earlier than APV)
    {{"run_0788", "run379"}, {350758, 18140}}, // First APV vs VMM dt = 15829384 or 1582 pulses
};
// Alternative: first VMM pulser, usiailly about 33ms before first APV
map<pair<string, string>, pair<int, int>> firstPulserMapAlternative = {
    {{"run_0832", "run423"}, {179396, 7042}}, // First APV vs VMM dt = 23209318 or 2320 pulses
    {{"run_0832_cut", "run423_cut"}, {179396, 7042}}, // First APV vs VMM dt = 23209318 or 2320 pulses
    {{"run_0832_cut10m", "run423_cut10m"}, {179396, 7042}}, // First APV vs VMM dt = 23209318 or 2320 pulses
    {{"run_0831", "run422"}, {151064, 7675}}, // First APV vs VMM dt = 32988 or 3 pulses
    {{"run_0830", "run421"}, {239117, 7633}}, // First APV vs VMM dt = 17748110 or 1774 pulses
    {{"run_0828", "run420"}, {138204, 3669}}, // First APV vs VMM dt = 32934 or 3 pulses
    // {{"run_0827", "run419"}, {, 7638}}, // First APV vs VMM dt = 18284390 or 1828 pulses // VMM started 30 sec before APV and counts normally
    {{"run_0826", "run418"}, {79477, 3928}}, // First APV vs VMM dt = 32524 or 3 pulses
    // {{"run_0825", "run417"} //No vmm pulses
    // {{"run_0824", "run416"}, {, 3429}}, // First APV vs VMM dt = 1477467 or 147 pulses // VMM started 0.7 sec before APV and counts normally
    {{"run_0821", "run412"}, {245098, 7662}}, // First APV vs VMM dt = 23988401 or 2398 pulses
    {{"run_0818", "run411"}, {141090, 3709}}, // First APV vs VMM dt = 26874351 or 2687 pulses
    // {{"run_0817", "run410"}, {53049, 126}}, // First APV vs VMM dt = 33004 or 3 pulses -- during acces
    // {{"run_0816", "run409"} // No apv pulses; No vmm pulses
    {{"run_0815", "run408"}, {237884, 7439}}, // First APV vs VMM dt = 25239911 or 2523 pulses
    {{"run_0814", "run407"}, {193345, 7572}}, // First APV vs VMM dt = 33027 or 3 pulses ; gaps between first n pulses in APV
    // {{"run_0813", "run406"} // No apv pulses
    // {{"run_0812", "run405"} // No apv pulses
    // {{"run_0809", "run399"} // No apv pulses
    // {{"run_0808", "run398"} // No apv pulses; No vmm pulses
    // {{"run_0807", "run397"} // No apv pulses
    {{"run_0805", "run395"}, {202356, 7800}}, // First APV vs VMM dt = 32394 or 3 pulses
    {{"run_0802", "run394"}, {217939, 9372}}, // First APV vs VMM dt = 32976 or 3 pulses
    {{"run_0801", "run393"}, {228878, 11715}}, // First APV vs VMM dt = 33058 or 3 pulses
    // {{"run_0800", "run392"} // No apv pulses
    // {{"run_0797", "run388"} // No apv pulses
    {{"run_0794", "run385"}, {329538, 11716}}, // First APV vs VMM dt = 41550608 or 4155 pulses
    {{"run_0792", "run383"}, {264602, 11227}}, // First APV vs VMM dt = 18083779 or 1808 pulses
    // {{"run_0791", "run382"} // First APV vs VMM dt = -5158615 or -515 pulses; *(VMM is earlier than APV)
    // {{"run_0789", "run380"} // First APV vs VMM dt = -6584644 or -658 pulses *(VMM is earlier than APV)
    {{"run_0788", "run379"}, {350095, 18140}}, // First APV vs VMM dt = 15829384 or 1582 pulses
};

constexpr bool PRINT_TO_FILE = true;
constexpr bool DEBUG_PRINT = false;
constexpr bool findBestVMM = true;
constexpr bool printMerged = false;
constexpr bool saveBackup = false;
constexpr bool saveTemporaryParameters = true;
constexpr bool reloopVMMFile = false;

// AlternativeStart -- first VMM selected as first good VMM pulser in case the difference between first good APV and VMM is about 33ms
void hitsMapper(bool tight = false, bool fixSRSTime = false, int nAll = 1, int n = 0, string runVMM="0832", string runAPV = "423", bool alternativeStart = false, int mergeTimeWindow = 1000)
{
    pair<string, string> run_pair = {Form("run_%s", runVMM.c_str()), Form("run%s", runAPV.c_str())};
    if(!firstPulserMap.count(run_pair))
    {
        printf("No such event in firstPulserMap! Exiting...\n");
        exit(0);
    }
    // pair<string, string> run_pair = {"run_0832_cut", "run423_cut"};
    // pair<int,int> firstSelectedEntries = {180038 + 21, 7042};
    pair<int,int> firstSelectedEntries = alternativeStart ? firstPulserMapAlternative.at(run_pair) : firstPulserMap.at(run_pair);

    auto apvan = new apv(run_pair.second);
    apvan->useSyncSignal();
    auto vmman = new evBuilder(run_pair.first, "g1_p25_s100-0&60", "map-20220605");
    vmman->useSyncSignal();

    cout << "Num of events: APV -- " << apvan->GetEntries() << "; VMM -- " << vmman->GetEntries() << endl;

    if(nAll < 1) nAll = 1;  
    int nEntriesPerRun = apvan->GetEntries() / nAll;
    int firstEntry = n * nEntriesPerRun;
    int lastEntry = (n == nAll - 1) ? apvan->GetEntries() - 1 : (n+1) * nEntriesPerRun - 1;    
    string numberingText = (nAll == 1) ? "" : Form("_%d-%d", n, nAll);
    printf("APV entries for analysis: [%d, %d]\n", firstEntry, lastEntry);
    
    long long startT_apv = 0;
    long long finishT_apv = 0;
    long long startT_pulse_apv = 0;
    long long T_apv = 0;
    long long T_apv_sincePulse = 0;
    long long T_apv_sincePulseSRSTS = 0; // in timestamps
    long long T_apv_pulse_prev = 0;
    long long T_vmm = 0;
    long long T_vmm_pulse_prev = 0;
    long long nPeriodsAPV = 0;
    // if time between two syncs larger then multiple pulser periods, merging fails due to calculation N periods for APV from pulser time only
    long long nPeriodsAPV_corrected = 0;

    int prevSRS = -1;
    int prev_pulse_SRS = -1;
    long long prev_pulse_time_s, prev_pulse_time_ms = 0;
    long long pulser_T = 0;
    long long apv_hit_T = 0;

    string tightText = tight ? "_tight" : "";
    string fixTimeText = fixSRSTime ? "_timefix" : "";
    string alternativeText = alternativeStart ? "_alternative" : "";
    string dupText = findBestVMM ? "" : "_firstVMM";

    ofstream out_APV;
    if(PRINT_TO_FILE)
        out_APV.open(TString("../out/APV_hits_maped_"+run_pair.first+"_"+run_pair.second+tightText+fixTimeText+dupText+alternativeText+numberingText+".txt").Data());

    // ofstream out_VMM;
    // out_VMM.open(TString("../out/VMM_hits_"+run_pair.first+"_"+run_pair.second+"_after"+tightText+fixTimeText+dupText+alternativeText+numberingText+".txt").Data());

    ofstream out_VMM_hits;
    if(PRINT_TO_FILE)
        out_VMM_hits.open(TString("../out/VMM_hits_UNmaped_"+run_pair.first+"_"+run_pair.second+tightText+fixTimeText+dupText+alternativeText+numberingText+".txt").Data());

    int numOfMapped = 0;

    auto out = TFile::Open(TString("../out/mapped_"+run_pair.first+"_"+run_pair.second+tightText+fixTimeText+dupText+alternativeText+numberingText+".root"), "recreate");

    static long long eventNumAPV = -1, eventNumVMM = -1;
    static long long deltaT;
    auto mappedEventNums = new TTree("mappedEvents", "");
    // mappedEventNums->AutoSave("1000");
    mappedEventNums->Branch("apv", &eventNumAPV);
    mappedEventNums->Branch("vmm", &eventNumVMM);
    mappedEventNums->Branch("deltaT", &deltaT);
    TFile* mappedEventBackupFile = nullptr;
    TString mapBackupFileName = TString("../out/mappedEvents_"+run_pair.first+"_"+run_pair.second+tightText+fixTimeText+dupText+alternativeText+numberingText+"_bak.root");
    // map<long long, pair<long long, int>> mappedEventMap = {};

    auto stripsVMM = make_shared<TH1F>("stripsVMM", "stripsVMM", 360, 0, 360);
    auto mappedHitsPdo = make_shared<TH1F>("mappedHitsPdo", "mappedHitsPdo", 2000, 0, 2000);
    auto mappedHitsPdo_apv = make_shared<TH1F>("mappedHitsPdo_apv", "mappedHitsPdo_apv", 2000, 0, 2000);
    auto hitsPdo = make_shared<TH2F>("hitsPdo", "hitsPdo", 360, 0, 360, 2000, 0, 2000);
    auto hitsPdo_apv = make_shared<TH2F>("hitsPdo_apv", "hitsPdo_apv", 360, 0, 360, 2000, 0, 2000);
    auto hitsPdoAll = make_shared<TH1F>("hitsPdoAll", "hitsPdoAll", 2000, 0, 2000);
    auto hitsPdoAll_apv = make_shared<TH1F>("hitsPdoAll_apv", "hitsPdoAll_apv", 2000, 0, 2000);
    auto layer2_vs_layer0 = make_shared<TH2F>("layer2_vs_layer0", "layer2_vs_layer0; layer2 strip №; #Delta strips", 200, 0, 200, 100, 0, 100);
    auto tbh_vmm = make_shared<TH1F>("tbh_vmm", "#Delta t between events with hits VMM in us", 2000, 0, 100000);
    auto tbh_vmm_mapped = make_shared<TH1F>("tbh_vmm_mapped", "#Delta t between events with mapped hits VMM in us", 2000, 0, 100000);
    auto tbh_vmm_umapped = make_shared<TH1F>("tbh_vmm_umapped", "#Delta t between events with unmapped hits VMM in us", 2000, 0, 100000);
    auto tbh_apv_mapped = make_shared<TH1F>("tbh_apv_mapped", "#Delta t between events with mapped hits APV in us", 2000, 0, 100000);
    auto tbh_apv_umapped = make_shared<TH1F>("tbh_apv_umapped", "#Delta t between events with unmapped hits APV in us", 2000, 0, 100000);
    auto tb_vmm_apv = make_shared<TH1F>("tb_vmm_apv", "#Delta t between events with mapped hits VMM and APV in us", 2000, 0, 100000);
    auto tb_u_vmm_apv = make_shared<TH1F>("tb_u_vmm_apv", "#Delta t between events with unmapped hits VMM and APV in us", 2000, 0, 100000);
    auto tbh_apv = make_shared<TH1F>("tbh_apv", "#Delta t between events with hits APV in us", 2000, 0, 20000);
    auto tracksPerTime = make_shared<TH1F>("tracksPerTime", "tracksPerTime; APV DAQ time, s", 400, 0, 400);
    auto vmmHits2DOPerTime = make_shared<TH1F>("vmmHits2DOPerTime", "vmmHits2DOPerTime; APV DAQ time, s", 400, 0, 400);
    auto vmmChannelMerged = make_shared<TH1F>("vmmChannelMerged", "vmmChannelMerged; channel", 56, 154, 210);
    auto vmmChannelAll = make_shared<TH1F>("vmmChannelAll", "vmmChannelAll; channel", 56, 154, 210);
    auto vmmChannelMultiplicityMerged = make_shared<TH1F>("vmmChannelMultiplicityMerged", "vmmChannelMultiplicityMerged; N channels", 64, 0, 64);
    auto vmmChannelMultiplicityAll = make_shared<TH1F>("vmmChannelMultiplicityAll", "vmmChannelMultiplicityAll; N channels", 64, 0, 64);
    auto vmmTrigScintAll = make_shared<TH1F>("vmmTrigScintAll", "vmmTrigScintAll; time, s", 400, 0, 400);
    auto vmmType3All = make_shared<TH1F>("vmmType3All", "vmmType3All; time, s", 400, 0, 400);
    auto vmmChannelTimeAll = make_shared<TH2F>("vmmChannelTimeAll", "vmmChannelTimeAll; time, s; channel", 400, 0, 400, 56, 154, 210);

    auto hDACTimeDiff = make_shared<TH1F>("hDACTimeDiff", "hDACTimeDiff; #Delta DAC time, #mus", 20000, -10000, 10000);
    auto hDACTimeDiffPerTime = make_shared<TH2F>("hDACTimeDiffPerTime", "hDACTimeDiffPerTime; time, s; #Delta DAC time, #mus", 4000, 0, 400, 2000, -10000, 10000);
    auto hnMerged = make_shared<TH1F>("hnMerged", "hnMerged; time, s;", 4000, 0, 400);
    auto hTimeDiffPerTime = make_shared<TH2F>("hTimeDiffPerTime", "hTimeDiffPerTime; time, s; #Delta T, #mus", 1000, 0, 1000, 2000, -1000, 1000);
    auto hTimeDiff2D = make_shared<TH2F>("hTimeDiff2D", "hTimeDiff2D; #Delta T (estimated), #mus; #Delta T (DAC), #mus", 2000, -1000, 1000, 2000, -1000, 1000);
    // dac time difference
    // N merged per time
    auto hAPVEstimatedTime = make_shared<TH1F>("hAPVEstimatedTime", "hAPVEstimatedTime; time, s; T_{DAC} - T_{estimated} , #mus", 4000, 0, 400);
    auto hAPVEstimatedTimeInSpill = make_shared<TH1F>("hAPVEstimatedTimeInSpill", "hAPVEstimatedTimeInSpill; time, s; (T_{pulse}-T_{pulse-1}) - #Sum(T_{estimated}) , #mus", 4000, 0, 400);
    auto hbcidDiffIgnored = make_shared<TH1F>("hbcidDiffIgnored", "hbcidDiffIgnored; BCID diff", 4096, 0, 4096);

    auto hVMMEstimatedTimeBetweenPulsers = make_shared<TH2F>("hVMMEstimatedTimeBetweenPulsers", "hVMMEstimatedTimeBetweenPulsers; time, s; T_{pulse}-T_{pulse-1} , #mus", 4000, 0, 400, 2000, -10000, 10000);
    
    std::function<int(int)> thrPerStrip = [](int strip){
        int thr = 50;
        switch(strip){
            case 161:
            case 194:
            case 199:
                thr = 1025; // effectively, disabled
                break;
            default:
                thr = 50;
                break;
        }
        return thr;
    };


    // out_VMM.close();

    // std::cout << "\t ---> TOTAL N of VMM events " << hits_vmm_t->GetEntries() << "\n";

    unsigned long startIndexAPV = firstSelectedEntries.second;
    long long vectorPositionInTree = firstSelectedEntries.first;
    unsigned long startIndex = 0;
    int prevSyncBcid = vmman->GetCentralHitsData(vectorPositionInTree).bcid;
    T_vmm_pulse_prev = vmman->GetCentralHitsData(vectorPositionInTree).timeFull();
    int prevPrevSyncBcid = 0;
    long long nPeriods = 0;
    long long nPeriodsAddAPV = 0;
    long long pulseTime = 0;
    long long T_vmm_pulse_first = T_vmm_pulse_prev;

    /* Try to load variables from previous runs */
    if(saveTemporaryParameters)
    {
        TString tmpFileName = TString("../out/beforeLastPulserParameters_"+run_pair.first+"_"+run_pair.second+tightText+fixTimeText+alternativeText+".root");
        auto file = TFile::Open(tmpFileName, "read");
        if(file != nullptr)
        {
            auto tree = static_cast<TTree*>(file->Get("beforeLastPulserParameters"));
            unsigned long apvN;
            long long bestAPVn, bestAPVni = -1;;
            tree->SetBranchAddress("apvN", &apvN);
            for(auto i = 0; i < tree->GetEntries(); i++)
            {
                tree->GetEntry(i);
                if(apvN > firstEntry) continue;
                if(bestAPVni < 0 || bestAPVn < apvN)
                {
                    bestAPVn = apvN;
                    bestAPVni = i;
                }
            }
            if(bestAPVni >= 0)
            {
                printf("Found parameters from previous running for event %lld\n", bestAPVn);
                tree->SetBranchAddress("T_apv", &T_apv);
                tree->SetBranchAddress("T_apv_sincePulseSRSTS", &T_apv_sincePulseSRSTS);
                tree->SetBranchAddress("T_apv_sincePulse", &T_apv_sincePulse);
                tree->SetBranchAddress("T_apv_pulse_prev", &T_apv_pulse_prev);
                tree->SetBranchAddress("nPeriodsAPV", &nPeriodsAPV);
                tree->SetBranchAddress("nPeriodsAPV_corrected", &nPeriodsAPV_corrected);
                tree->SetBranchAddress("prevSRS", &prevSRS);
                tree->SetBranchAddress("prev_pulse_SRS", &prev_pulse_SRS);
                tree->SetBranchAddress("vectorPositionInTree", &vectorPositionInTree);
                tree->SetBranchAddress("index", &startIndex);
                tree->SetBranchAddress("prevSyncBcid", &prevSyncBcid);
                tree->SetBranchAddress("prevPrevSyncBcid", &prevPrevSyncBcid);
                tree->SetBranchAddress("nPeriods", &nPeriods);
                tree->SetBranchAddress("pulseTime", &pulseTime);
                tree->SetBranchAddress("T_vmm_pulse_prev", &T_vmm_pulse_prev);
                tree->GetEntry(bestAPVni);
                tree->ResetBranchAddresses();
                startIndexAPV = bestAPVn;
            }
            delete tree;
            file->Close();
            delete file;
       }
    }
    
    tuple<long long, unsigned long, int, int, long long, long long, long long> beforeLastPulserParameters = {vectorPositionInTree, startIndex, prevSyncBcid, prevPrevSyncBcid, nPeriods, pulseTime, T_vmm_pulse_prev};
    map<long long, vector<pair<unsigned long, analysisGeneral::mm2CenterHitParameters>>> hits_vmm_events_map;
    auto nLoaded = loadNextVMM(get<0>(beforeLastPulserParameters), hits_vmm_events_map, vmman);
    // printf("nLoaded: %lld\n", nLoaded);
    if(nLoaded <= 1)
        return;
    auto beforeLastPulserParametersCurrent = beforeLastPulserParameters;
    tuple<long long, long long, tuple<long long, unsigned long, int, int, long long, long long, long long>, unsigned long, int> bestHit;
    map<long long, int> pairedVMM = {};
    analysisGeneral::mm2CenterHitParameters* currEvent;

    int mappedHitsVMM = 0;
    int UNmappedHitsVMM = 0;

    int bad;
    long long currentEventsMap;
    bool hitMapped;

    vector<pair<int, int>> vmm_hits_vec;
    vector<pair<int, int>> apv_hits_vec;
    vector<pair<int, int>> apv_hits_vec_l0;
    vector<pair<int, int>> apv_hits_vec_l1;
    int diff_hit;
    long long dt_apv_vmm;
    int strip, pdo;
    int diff, diffDiff;
    optional<int> nPeriodsAdd;
    long long diffTvmm;

    int maxAPVPulserCountDifference = tight ? 0 : 1;
    bool ignoreNextSync = true;

    bool printed = false;
    unsigned long i;
    for (i = startIndexAPV; i < apvan->GetEntries(); i++)
    {
        auto hit_apv = apvan->GetCentralHits2ROnlyData(i);
        T_apv = hit_apv.timeSec*1E6 + hit_apv.timeMSec;
        if (prevSRS == -1)
        {
            startT_apv = T_apv;
            printf("startT_apv: %lld\n", T_apv);
        }
        else
        {
            T_apv_sincePulseSRSTS += (hit_apv.timeSrs - prevSRS) + ((prevSRS > hit_apv.timeSrs)?16777216:0);
            // T_apv_sincePulse += apv_time_from_SRS(prevSRS, hit_apv.timeSrs, fixSRSTime);
            T_apv_sincePulse = static_cast<long long>(apvTimeTimeSRSTimestamp(T_apv_sincePulseSRSTS, fixSRSTime));
        }
        prevSRS = hit_apv.timeSrs;

        if (hit_apv.sync)
        {
            if (prev_pulse_SRS == -1)
            {
                startT_pulse_apv = T_apv;
                printf("First APV sync-signal: %lu\n", i);
            }
            else
            {
                hAPVEstimatedTimeInSpill->Fill((T_apv - startT_apv)*1.0/1E6, round((T_apv - T_apv_pulse_prev) * 1.0 / 1e4)*1e4 - T_apv_sincePulse);
                nPeriodsAPV += round((T_apv - T_apv_pulse_prev) * 1.0 / 1e4);
                pulser_T = nPeriodsAPV * 1e4;
            }
            // printf("T_apv: %lld;\tT_apv_prev: %lld;\tdiff: %g (%g)\n", T_apv, T_apv_pulse_prev, (T_apv - T_apv_pulse_prev) * 1.0 / 1e4, round((T_apv - T_apv_pulse_prev) * 1.0 / 1e4));
            // printf("N periods APV: %lld\n", nPeriodsAPV);
            // if(T_apv_sincePulse > 1e6)
            //     printf("SRS Timestamps since last pulser: %lld (%%400037 = %lld) = (%lld ms if 25ns, %lld ms if smaller), %lld found -- theoretically, %g ns per srsTimeStamp\n",
            //            T_apv_sincePulseSRSTS,
            //            T_apv_sincePulseSRSTS%400037,
            //            static_cast<long long>(apvTimeTimeSRSTimestamp(T_apv_sincePulseSRSTS, false)),
            //            static_cast<long long>(apvTimeTimeSRSTimestamp(T_apv_sincePulseSRSTS, true)),
            //            T_apv_sincePulse,
            //            1000.0 * T_apv_sincePulse / T_apv_sincePulseSRSTS
            //         );
            //     // printf("SRS Timestamps since last pulser: %lld = %lld microseconds (%lld)\n", T_apv_sincePulseSRSTS, static_cast<long long>(apvTimeTimeSRSTimestamp(T_apv_sincePulseSRSTS, fixSRSTime)), T_apv_sincePulse);
            T_apv_pulse_prev = T_apv;
            prev_pulse_SRS = hit_apv.timeSrs;
            T_apv_sincePulse = 0;
            T_apv_sincePulseSRSTS = 0;
            // std::cout << "Period " << nPeriodsAPV << "--- is sync! N = " << hit_apv.hitsPerLayer.at(0).size() << "\n";
            // apvTimeTimeSRSTimestamp(int diff, bool fixSRSTime = false)
        }
        nPeriodsAPV_corrected = nPeriodsAPV + T_apv_sincePulse / 10000;
        // if(nPeriodsAPV_corrected != nPeriodsAPV)
        //     printf("nPeriodsAPV_corrected != nPeriodsAPV: time %lld, previous pulser: %lld, difference by DAC: %lld, by srsTimeStamp: %lld\n", T_apv, T_apv_pulse_prev, T_apv - T_apv_pulse_prev, T_apv_sincePulse);

        if (i < firstEntry)
            continue;
        else if (i == firstEntry)
            printf("Merging started...\n");
        else if (i == lastEntry)
            printf("Finish time: %lld\n", T_apv);
        else if(i > lastEntry)
            break;

        hAPVEstimatedTime->Fill((T_apv - startT_apv)*1.0/1E6, (T_apv - (T_apv_pulse_prev + T_apv_sincePulse)));

        if(!printed)
        {
            printed = true;
            printf("Start time: %lld\n", T_apv);
        }
        
        if(!hit_apv.hits.size())
            continue;
        
        if (prev_pulse_SRS != -1)
        {
            // std::cout << "\t Total:" << prev_pulse_SRS << "\t" << hit_apv.timeSrs << "\n";

            // apv_hit_T = pulser_T + apv_time_from_SRS(prev_pulse_SRS, hit_apv.timeSrs, fixSRSTime);
            if (hit_apv.hitsPerLayer.at(0).size() != 0)
            {
                apv_hits_vec.clear();
                apv_hits_vec_l0.clear();
                apv_hits_vec_l1.clear();
                // std::cout << "---> with hits:" << prev_pulse_SRS << "\t" << hit_apv.timeSrs << "\n";
                // out_APV << "------- APV Period " << nPeriodsAPV << " -------- \n";
                for (auto &h : hit_apv.hitsPerLayer.at(0))
                {
                    strip = h.first; // TODO add dependency on testbeam
                    pdo = h.second;
                    if (strip > 118 && strip < 172)
                        apv_hits_vec_l0.push_back(make_pair(strip, pdo));
                }
                for (auto &h : hit_apv.hitsPerLayer.at(1))
                {
                    strip = h.first * (1 - 2.29e-3) - 2.412 / 0.25; // TODO add dependency on testbeam
                    pdo = h.second;
                    if (strip > 118 && strip < 172)
                        apv_hits_vec_l1.push_back(make_pair(strip, pdo));
                }
                for (auto &h : hit_apv.hitsPerLayer.at(2))
                {
                    strip = h.first * (1 - 8e-3) - 8.46 / 0.25;// TODO add dependency on testbeam
                    pdo = h.second;
                    if (strip > 118 && strip < 172)
                        apv_hits_vec.push_back(make_pair(strip, pdo));
                }

                if (apv_hits_vec_l0.size() && apv_hits_vec_l1.size())
                    tracksPerTime->Fill((T_apv - startT_apv)/1E6);
                
                if (!apv_hits_vec.size() || !apv_hits_vec_l0.size() || !apv_hits_vec_l1.size())
                    continue;


                map<int, double> means = {{0, weightedMean(apv_hits_vec_l0)}, {1, weightedMean(apv_hits_vec_l1)}}; // TODO add dependency on testbeam
                auto tr = getEstimatedTrack(means);
                int propogated = static_cast<int>(round(estimatePositionInLayer(tr, 2)));

                vmm_hits_vec.clear();

                bad = 0;
                prevSyncBcid = get<2>(beforeLastPulserParameters);
                prevPrevSyncBcid = get<3>(beforeLastPulserParameters);
                nPeriods = get<4>(beforeLastPulserParameters);

                pulseTime = get<5>(beforeLastPulserParameters);
                hitMapped = false;
                beforeLastPulserParametersCurrent = beforeLastPulserParameters;
                bestHit = {-1, 0, beforeLastPulserParametersCurrent, 0, 0};

                vectorPositionInTree = get<0>(beforeLastPulserParameters); // TODO
                currentEventsMap = hits_vmm_events_map.at(vectorPositionInTree).size();

                T_vmm_pulse_prev = get<6>(beforeLastPulserParameters);

                if(DEBUG_PRINT)
                    printf("For APV event %lu search starting from VMM event %lu\n", i, get<1>(beforeLastPulserParameters));
                for (unsigned long j = get<1>(beforeLastPulserParameters); j <= currentEventsMap; j++)
                {
                    if(j == currentEventsMap)
                    {
                        vectorPositionInTree = get<0>(beforeLastPulserParameters) + hits_vmm_events_map.at(vectorPositionInTree).size();
                        j = 0;
                    }
                    if(j == 0)
                    {
                        auto nLoaded = loadNextVMM(vectorPositionInTree, hits_vmm_events_map, vmman);
                        if(nLoaded <= 1) break;
                        currentEventsMap = hits_vmm_events_map.at(vectorPositionInTree).size();
                    }
                    currEvent = &(hits_vmm_events_map.at(vectorPositionInTree).at(j).second);
                    if (currEvent->sync)
                    {
                        T_vmm = currEvent->timeSec*1E6 + currEvent->timeMSec;
                        if(ignoreNextSync)
                        {
                            ignoreNextSync = false;
                        }
                        else
                        {
                            beforeLastPulserParametersCurrent = {vectorPositionInTree, j, prevSyncBcid, prevPrevSyncBcid, nPeriods, pulseTime, T_vmm_pulse_prev};
                            diff = (currEvent->bcid - prevSyncBcid > 0) ? currEvent->bcid - prevSyncBcid : currEvent->bcid + 4096 - prevSyncBcid;
                            diffDiff = (currEvent->bcid - prevPrevSyncBcid > 0) ? currEvent->bcid - prevPrevSyncBcid: currEvent->bcid + 4096 - prevPrevSyncBcid;

                            nPeriodsAdd = calculateVMMNPulsers(diff, 3, 4); // TODO add dependency on testbeam
                            
                            diffTvmm = T_vmm-T_vmm_pulse_prev;

                            if(!nPeriodsAdd)
                            {
                                continue;
                            }
                            nPeriods += nPeriodsAdd.value();
                        }
                        prevPrevSyncBcid = prevSyncBcid;
                        prevSyncBcid = currEvent->bcid;
                        pulseTime = nPeriods * 50;
                        T_vmm_pulse_prev = T_vmm;
                    }

                    if(DEBUG_PRINT)
                        printf("vmm %lu | T APV: %lld; T VMM: %lld; difference: %lld; nPeriodsAPV: %lld; nPeriods APV, corrected: %lld; nPeriodsVMM: %lld; difference: %lld (%lld); time since pulser: %lld\n",
                           j, T_apv, currEvent->timeFull(), T_apv - currEvent->timeFull(),
                           nPeriodsAPV, nPeriodsAPV_corrected, nPeriods, nPeriodsAPV_corrected*200 - nPeriods,
                           static_cast<long long>((nPeriodsAPV_corrected*200 - nPeriods)*50),
                           T_apv_sincePulse%10000);
                    
                    if (nPeriods / 200 > nPeriodsAPV_corrected + maxAPVPulserCountDifference)
                    {
                        break;
                    }
                    else if (nPeriods / 200 < nPeriodsAPV_corrected - maxAPVPulserCountDifference)
                    {
                        // printf("%lld / 200 = %lld < %lld - %d\n", nPeriods, nPeriods / 200, nPeriodsAPV_corrected, maxAPVPulserCountDifference);
                        if(beforeLastPulserParametersCurrent != beforeLastPulserParameters)
                        {
                            beforeLastPulserParameters = beforeLastPulserParametersCurrent;
                            if(DEBUG_PRINT)
                                printf("beforeLastPulserParameters updated\n");
                        }
                        if(i == firstEntry && j == 0)
                        {
                          freeMemory(hits_vmm_events_map, get<0>(beforeLastPulserParameters));
                        }
                        continue;
                    }
                    else if (currEvent->hitsX.size() == 0)
                        continue;
                    else if (!currEvent->trigger) // IMPORTANT change: passing event without triple scint coincsidence
                        continue;
 
                    diff_hit = (currEvent->bcid - prevSyncBcid >= 0) ? currEvent->bcid - prevSyncBcid : currEvent->bcid + 4096 - prevSyncBcid;
                    dt_apv_vmm = T_apv_sincePulse - static_cast<long long>((nPeriods - nPeriodsAPV*200) * 50) - static_cast<long long>(round(diff_hit * 25.0 / 1000.0));
                    // printf("dt_apv_vmm = %lld - (%lld - %lld*200) * 50 - %d*25/1000 = %lld\n", T_apv_sincePulse, nPeriods, nPeriodsAPV, diff_hit, dt_apv_vmm);
                    // int dt_apv_vmm = T_apv - startT_pulse_apv - pulseTime;

                    // std::cout << nPeriods / 200 << " \t " << hits_vmm_event->second.hitsX.size() << "\n";
                    // out_APV << "------- VMM Period " << nPeriods / 200 << "  (" << nPeriods % 200 << ") -------- dT = " << dt_apv_vmm << "\n";
                    // printf("APV %lu, VMM %lld -- %lld\n", i, vectorPositionInTree + j, dt_apv_vmm);
                    if (abs(dt_apv_vmm) > mergeTimeWindow)
                        continue;

                    for (auto it = currEvent->hitsX.begin(); it != currEvent->hitsX.end(); ++it)
                    {
                        strip = it->first * (1 - 8e-3) - 8.46 / 0.25; // TODO add dependency on testbeam
                        pdo = it->second;
                        if (pdo < thrPerStrip(it->first))
                            continue;
                        vmm_hits_vec.push_back(make_pair(strip, pdo));
                        // out_APV << "\t Strip: " << it->first << "\n";
                        // out_APV << "\t PDO: " << it->second << "\n";
                    }
                    if (vmm_hits_vec.size() == 0)
                        continue;

                    hitMapped = false;

                    for (unsigned long l = 0; l < vmm_hits_vec.size(); l++)
                    {
                        if (abs(propogated - vmm_hits_vec.at(l).first) < 5){
                            hitMapped = true;
                            break;
                        }
                    }

                    // for (unsigned long k = 0; k < apv_hits_vec.size(); k++)
                    // {
                    //     for (unsigned long l = 0; l < vmm_hits_vec.size(); l++)
                    //     {
                    //         if (vmm_hits_vec.size() != 0 && abs(apv_hits_vec.at(k).first - vmm_hits_vec.at(l).first) < 5){
                    //             hitMapped = true;
                    //             break;
                    //         }
                    //     }
                    //     if (hitMapped)
                    //         break;
                    // }

                    if (hitMapped)
                    {
                        if(get<0>(bestHit) < 0 || abs(get<1>(bestHit)) > abs(dt_apv_vmm))
                        {
                            int hitsMappedInEvent = 0;
                            for (unsigned long l = 0; l < vmm_hits_vec.size(); l++)
                                hitsMappedInEvent++;
                            bestHit = {hits_vmm_events_map.at(vectorPositionInTree).at(j).first, dt_apv_vmm, beforeLastPulserParametersCurrent, j, hitsMappedInEvent};
                        }
                        {
                            if(DEBUG_PRINT)
                            {
                            std::cout << "APV event: " << i << " (" << i << ")" << "\t";
                            std::cout << "VMM event: " << j << " (" << hits_vmm_events_map.at(vectorPositionInTree).at(j).first << ")" << "\t Total of mapped " << numOfMapped << "\n";
                            printf("dt_apv_vmm = %lld - %lld - %lld - round(%d * 25.0 / 1000.0) = %lld\n", T_apv, startT_pulse_apv, pulseTime, diff_hit, dt_apv_vmm);
                            }
                            // printf("Event parameters: %lld, %lu, %d, %d, %d, %lld\n", get<0>(beforeLastPulserParameters), get<1>(beforeLastPulserParameters), get<2>(beforeLastPulserParameters),
                            //        get<3>(beforeLastPulserParameters), get<4>(beforeLastPulserParameters), get<5>(beforeLastPulserParameters));
                            // std::cout << "!!! Mapped hits: " << numOfMapped << " !!! \n";
                            if(PRINT_TO_FILE)
                            {
                                if (T_apv_sincePulse < 10e3)
                                {
                                    out_APV << "------- APV Double ReadOut Period shiiiish " << nPeriodsAPV << " -------- T = " << T_apv - startT_pulse_apv << " [us] \n";
                                }
                                else
                                {
                                    out_APV << "------- APV Double ReadOut Period " << nPeriodsAPV << " -------- T = " << T_apv - startT_pulse_apv << " [us] \n";
                                }
                            }
                            for (int k = 0; k < apv_hits_vec.size(); k++)
                            {
                                if(PRINT_TO_FILE)
                                {
                                    out_APV << k << "\t Strip: " << apv_hits_vec.at(k).first << "\n";
                                    out_APV << "  \t PDO: " << apv_hits_vec.at(k).second << "\n";
                                }
                                mappedHitsPdo_apv->Fill(apv_hits_vec.at(k).second);
                            }
                            if(PRINT_TO_FILE)
                            {
                                out_APV << "------- APV Layer 0 -------- \n";
                                for (int k = 0; k < apv_hits_vec_l0.size(); k++)
                                {
                                    out_APV << k << "\t Strip: " << apv_hits_vec_l0.at(k).first << "\n";
                                    out_APV << "  \t PDO: " << apv_hits_vec_l0.at(k).second << "\n";
                                }
                                out_APV << "------- APV Layer 1 -------- \n";
                                for (int k = 0; k < apv_hits_vec_l1.size(); k++)
                                {
                                    out_APV << k << "\t Strip: " << apv_hits_vec_l1.at(k).first << "\n";
                                    out_APV << "  \t PDO: " << apv_hits_vec_l1.at(k).second << "\n";
                                }

                                out_APV << "------- VMM Period " << nPeriods / 200 << "  (" << nPeriods % 200 << ") -------- T = " << pulseTime + diff_hit * 25 / 1000 << "[us] (dT = " << dt_apv_vmm << " [us]) \n";
                            }

                            for (int l = 0; l < vmm_hits_vec.size(); l++)
                            {
                                if(PRINT_TO_FILE)
                                {
                                    out_APV << "\t Strip: " << vmm_hits_vec.at(l).first << "\n";
                                    out_APV << "\t PDO: " << vmm_hits_vec.at(l).second << "\n";
                                }
                                mappedHitsPdo->Fill(vmm_hits_vec.at(l).second);
                                // mappedHitsVMM++;
                            }
                            // tbh_apv_mapped->Fill(T_apv - startT_pulse_apv - prevT_apv);
                            // tbh_vmm_mapped->Fill(pulseTime + diff_hit * 25 / 1000 - prevT_vmm);
                            // tb_vmm_apv->Fill(abs(dt_apv_vmm));

                            // prevT_apv = T_apv - startT_pulse_apv;
                            // prevT_vmm = pulseTime + diff_hit * 25 / 1000;
                        }
                        if(findBestVMM)
                        {
                            hitMapped = false;
                            continue;
                        }
                        else
                        {
                            break;
                        }
                    }
                    else
                    {
                        if (vmm_hits_vec.size() != 0)
                        {
                            if(PRINT_TO_FILE)
                                out_VMM_hits << "------- VMM Period " << nPeriods / 200 << "  (" << nPeriods % 200 << ") -------- T = " << pulseTime + diff_hit * 25 / 1000 << "[us] (dT = " << dt_apv_vmm << " [us]) \n";

                            for (int l = 0; l < vmm_hits_vec.size(); l++)
                            {
                                if(PRINT_TO_FILE)
                                {
                                    out_VMM_hits << "\t Strip: " << vmm_hits_vec.at(l).first << "\n";
                                    out_VMM_hits << "\t PDO: " << vmm_hits_vec.at(l).second << "\n";
                                }
                                UNmappedHitsVMM++;
                            }
                        }
                    }

                }
                if (get<0>(bestHit) >= 0)
                {
                    numOfMapped++;
                    // beforeLastPulserParameters = get<2>(bestHit);
                    if(printMerged)
                    {
                        printf("APV event: %9lu (%9lu)\t", i, i);
                        printf("VMM event: %9lu (%9lld)\t", get<3>(bestHit), get<0>(bestHit));
                        printf("Time difference: %9lld\t", get<1>(bestHit));
                        printf("Total of mapped: %d\n", numOfMapped);
                    }
                    eventNumAPV = i;
                    eventNumVMM = get<0>(bestHit);
                    deltaT = get<1>(bestHit);
                    mappedEventNums->Fill();


                    {
                        // auto vmmEvent = vmman->GetCentralHitsData(eventNumVMM);
                        auto vmmEvent = hits_vmm_events_map.at(get<0>(get<2>(bestHit))).at(get<3>(bestHit)).second;
                        vmmChannelMultiplicityMerged->Fill(vmmEvent.hitsX.size());
                        auto dacTimeDiff = T_apv - static_cast<long long>(vmmEvent.timeSec*1E6 + vmmEvent.timeMSec);
                        hDACTimeDiff->Fill(dacTimeDiff);
                        hDACTimeDiffPerTime->Fill((T_apv-startT_apv)*1.0/1E6, dacTimeDiff);
                        hTimeDiffPerTime->Fill((T_apv-startT_apv)*1.0/1E6,deltaT);
                        hTimeDiff2D->Fill(deltaT,dacTimeDiff);
                        hnMerged->Fill((T_apv-startT_apv)*1.0/1E6);
                        for (auto &h: vmmEvent.hitsX)
                        {
                            pdo = h.second;
                            if (pdo < thrPerStrip(h.first)) continue;
                            vmmChannelMerged->Fill(h.first);
                        }
                    }

                    // clear memory -- remove unused vectors with VMM events
                    freeMemory(hits_vmm_events_map, get<0>(beforeLastPulserParameters));
                    if(!(numOfMapped %1000))
                    {
                        if(saveBackup)
                        {
                            mappedEventNums->SaveAs(mapBackupFileName);
                            delete mappedEventNums;
                            if(mappedEventBackupFile != nullptr )
                                mappedEventBackupFile->Close();
                            mappedEventBackupFile = TFile::Open(mapBackupFileName);
                            mappedEventNums = static_cast<TTree*>(mappedEventBackupFile->Get("mappedEvents"));
                            // mappedEventNums->AutoSave("1000");
                            mappedEventNums->SetBranchAddress("apv", &eventNumAPV);
                            mappedEventNums->SetBranchAddress("vmm", &eventNumVMM);
                            mappedEventNums->SetBranchAddress("deltaT", &deltaT);
                        }
                        if(!printMerged)
                        {
                            printf("APV event: %9lu;\t", i);
                            printf("Total of mapped: %d\n", numOfMapped);
                        }
                    }
                    mappedHitsVMM+=get<4>(bestHit);
                }
            }
        }
    }
    std::cout << "Total number of mapped: " << numOfMapped << "\n";
    std::cout << "N of MAPPED hits in VMM " << mappedHitsVMM << "\n";
    out_VMM_hits << "N of unMAPPED hits in VMM " << UNmappedHitsVMM << "\n";

    out_APV.close();
    out_VMM_hits.close();

    if(reloopVMMFile)
    {
        long long prevPulse = -1;
        for (i = firstSelectedEntries.first; i < vmman->GetEntries(); i++)
        {
            auto hit_vmm = vmman->GetCentralHitsData(i);
            auto time = hit_vmm.timeSec*1E6 + hit_vmm.timeMSec;
            vmmChannelMultiplicityAll->Fill(hit_vmm.hitsX.size());
            if(time < startT_apv)
                continue;
            else if (time > T_apv)
                break;
            if(!hit_vmm.hitsX.size())
                continue;
            auto timeToAPV = (time - startT_apv)/1E6;
            bool haveMMSignals = false;
            for (auto &h: hit_vmm.hitsX){
                if (h.second < thrPerStrip(h.first)) continue;
                haveMMSignals = true;
                vmmChannelAll->Fill(h.first);
                vmmChannelTimeAll->Fill(timeToAPV, h.first);
            }
            if(hit_vmm.sync){
                if(prevPulse >= 0)
                    hVMMEstimatedTimeBetweenPulsers->Fill(timeToAPV, time - prevPulse);
                prevPulse = time;
            }
            if(!haveMMSignals)
                continue;
            vmmHits2DOPerTime->Fill(timeToAPV);
            // hit_vmm.signal = false;
            if(hit_vmm.sync || hit_vmm.trigger){
                vmmTrigScintAll->Fill(timeToAPV);
            } else {
                vmmType3All->Fill(timeToAPV);
            }

            // vmm_hits_vec.clear();
            // for (auto it = hit_vmm.hitsX.begin(); it != hit_vmm.hitsX.end(); ++it)
            // {
            //     strip = it->first * (1 - 8e-3) - 8.46 / 0.25;
            //     pdo = it->second;
            //     if (pdo < 50)
            //         continue;
            //     vmm_hits_vec.push_back(make_pair(strip, pdo));
            // }
            // if(!
        
        }
    }

    hTimeDiff2D->Scale(1.0 / mappedEventNums->GetEntries());

    out->Write();
    out->Close();
    delete apvan;
    delete vmman;

    if(saveTemporaryParameters)
    {
        TString tmpFileName = TString("../out/beforeLastPulserParameters_out_"+run_pair.first+"_"+run_pair.second+tightText+fixTimeText+dupText+alternativeText+numberingText+".root");
        auto file = TFile::Open(tmpFileName, "recreate");
        auto tree = new TTree("beforeLastPulserParameters", "");
        tree->SetDirectory(file);
        tree->Branch("apvN", &i);
        tree->Branch("T_apv", &T_apv);
        tree->Branch("T_apv_sincePulseSRSTS", &T_apv_sincePulseSRSTS);
        tree->Branch("T_apv_sincePulse", &T_apv_sincePulse);
        tree->Branch("T_apv_pulse_prev", &T_apv_pulse_prev);
        tree->Branch("nPeriodsAPV", &nPeriodsAPV);
        tree->Branch("nPeriodsAPV_corrected", &nPeriodsAPV_corrected);
        tree->Branch("prevSRS", &prevSRS);
        tree->Branch("prev_pulse_SRS", &prev_pulse_SRS);
        long long vectorPositionInTree;
        unsigned long index;
        int prevSyncBcid;
        int prevPrevSyncBcid;
        long long nPeriods;
        long long pulseTime;
        tree->Branch("vectorPositionInTree", &vectorPositionInTree);
        tree->Branch("index", &index);
        tree->Branch("prevSyncBcid", &prevSyncBcid);
        tree->Branch("prevPrevSyncBcid", &prevPrevSyncBcid);
        tree->Branch("nPeriods", &nPeriods);
        tree->Branch("pulseTime", &pulseTime);
        tree->Branch("T_vmm_pulse_prev", &T_vmm_pulse_prev);
        vectorPositionInTree = get<0>(beforeLastPulserParameters); 
        index = get<1>(beforeLastPulserParameters);
        prevSyncBcid = get<2>(beforeLastPulserParameters);
        prevPrevSyncBcid = get<3>(beforeLastPulserParameters);
        nPeriods = get<4>(beforeLastPulserParameters);
        pulseTime = get<5>(beforeLastPulserParameters);
        T_vmm_pulse_prev = get<6>(beforeLastPulserParameters);
        tree->Fill();
        tree->Write();
        file->Close();
    }

}

int main(int argc, char** argv)
{
    hitsMapper();
    return 0;
}
