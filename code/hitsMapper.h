#pragma once

#include "analysisGeneral.h"
#include "apv.C"
#include "evBuilder.C"
#include <iostream>
#include <fstream>

#include <vector>
#include <utility>
#include <tuple>
#include <map>
#include <set>

using std::cout, std::endl;
using std::ifstream, std::ofstream;

using std::vector;
using std::pair, std::make_pair;
using std::tuple, std::get;
using std::map;
using std::set;

int calculateVMMNPulsers(int bcidDiff, int maxDiff = 2, int maxNPulsers = 9)
{
    bcidDiff += (bcidDiff > 0) ? 0 : 4096;
    int predicted;
    for(auto i = 1; i <= maxNPulsers; i++)
    {
      predicted = (2000 * i) % 4096;
      if(abs(predicted - bcidDiff) <= maxDiff)
        return i;
    }
    return -1;
}

// Long should be enough (4 bytes), but to be sure we can use long long (8 bytes)
int calculateAPVNPulsers(long long TimeStampDiff, int maxDiff = 1, int maxNPulsers = 1000) // 1000 puslesr = 10s
{
    static constexpr long long maxCounter = static_cast<long long>(1<<24);
    static constexpr double pulserLength = 400036.0;
    TimeStampDiff += (TimeStampDiff > 0) ? 0 : maxCounter;
    long long predicted;
    for(auto i = 1; i <= maxNPulsers; i++)
    {
        predicted = static_cast<long long>(round(pulserLength * i)) % maxCounter;
        if(abs(predicted - TimeStampDiff) <= maxDiff)
            return i;
    }
    return -1;
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

double weightedMean(vector<pair<int, int>> hits){ // TODO correct to pair of position and error
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

constexpr int layerStrawOdd = -1;
constexpr int layerStrawEven = -2;
constexpr int layerDR = -3;
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
 * Layers layerStrawOdd and layerStrawEven are uned for straws (odd and even-numbered)
 * Layer layerDR is for auto-selection double read-out layer
 * return: position in mm, from Layer 0
 */
int getLayerPosition(int layer){
  static const auto tb = analysisGeneral::GetTestBeam();
  int y = 0;

  switch(tb){
    case analysisGeneral::TestBeams::TB22_April:
      switch(layer){
        case layerStrawEven:
        case layerStrawOdd:
        case 3:
          y += 523;
        case layerDR:
        case 2:
          y += 345;
        case 1:
          y += 285;
        case 0:
          y += 0;
      }
      break;
    case analysisGeneral::TestBeams::TB22_July:
      switch(layer){
        case layerStrawEven:
        case layerStrawOdd:
          y = -183;
          break;
        case 2:
          y += 325;
        case 3:
          y += 32;
        case 1:
          y += 434;
        case layerDR:
        case 0:
          y += 0;
          break;
        default:
          y += 0;
          break;
      }
      break;
    case analysisGeneral::TestBeams::TB22_August:
      break;
    case analysisGeneral::TestBeams::TB22_October:
      break;
    default:
      break;
  }
  return y;
}

double correctAlignment(int strip, int layer){
  static const auto tb = analysisGeneral::GetTestBeam();
  double stripD = 0;
  switch(tb){
    case analysisGeneral::TestBeams::TB22_October:
      break;
    case analysisGeneral::TestBeams::TB22_August:
      break;
    case analysisGeneral::TestBeams::TB22_July:
      switch(layer){
        case layerDR:
        case 0:
          stripD = static_cast<double>(strip) * 1.09 - 6.16;
          break;
        case 1:
          stripD = static_cast<double>(strip);
          break;
        case 2:
          stripD = static_cast<double>(strip) * 1.01 + 3.77;
          break;
        default:
          stripD = static_cast<double>(strip);
          break;
      }
      break;
    case analysisGeneral::TestBeams::TB22_April:
      // shifts from Stefano
      switch(layer){
        case 0:
          stripD = static_cast<double>(strip);
          break;
        case 1:
          stripD = static_cast<double>(strip) * (1 - 2.29e-3) - 2.412 / 0.25;
          break;
        case layerDR:
        case 2:
          stripD = static_cast<double>(strip) * (1 - 8e-3) - 8.46 / 0.25;
          break;
        default:
          stripD = static_cast<double>(strip);
          break;
      }
      break;
    default:
      break;
  }
  return stripD;
}

pair<double, double> getEstimatedTrack(map<int, pair<double, double>> positions){
  int n = 0;
  double sumWXY = 0, sumWX = 0, sumWY = 0, sumWX2 = 0, sumY2 = 0, sumW = 0;
  double y, x, w, yE; // x: horizontal, y - vertical coordinate
  for(auto &pos: positions){
    x = getLayerPosition(pos.first);
    y = pos.second.first;
    yE = pos.second.second;
    w = 1.0 / yE;
    sumW += w;
    sumWXY += w * x * y;
    sumWX += w * x;
    sumWY += w * y;
    sumWX2 += w * x * x;
    // sumY2 += y * y;
    n++;
  }
  double d = sumW * sumWX2 - sumWX * sumWX;
  double a0 = (sumW*sumWXY - sumWX*sumWY) / d;
  double b0 = (sumWX2 * sumWY - sumWX*sumWXY)/d;
  double sigmaa = sqrt(sumW/d);
  double sigmab = sqrt(sumWX2/d);
  return {a0, b0};
}
pair<double, double> getEstimatedTrack(map<int, double> positions){
  map<int, pair<double, double>> positionsW;  
  for(auto &pos: positions){
    positionsW.emplace(pos.first, make_pair(pos.second, 1.0));
  }
  return getEstimatedTrack(positionsW);
}

double estimatePositionInLayer(pair<double, double> trackAB, int layer){
  double x = getLayerPosition(layer);
  double y = trackAB.first * x + trackAB.second;
  return y;
}
