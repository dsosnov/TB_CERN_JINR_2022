#pragma once

#include <optional>
#include <utility>
#include <algorithm>
#include <map>
#include <vector>

using std::vector;
using std::map;
using std::optional, std::nullopt;
using std::pair, std::make_pair;
constexpr unsigned int layerDoubleReadout = 2; // 0

// From analyseMerged

vector<map<double, int>> splitByDistance(map<double, int> hitsPerLayer, double maxDistance = 5){
  vector<map<double, int>> out;
  optional<double> prev = nullopt;
  for(auto &h: hitsPerLayer){ // hits are sorted by the rules of map structure
    if(!out.size() || !prev || fabs(prev.value() - h.first) > maxDistance)
      out.push_back({});
    out.back().emplace(h.first, h.second);
    prev = h.first;
  }
  return out;
}

double correctAlignment(int strip, int layer){
    // shifts from Stefano
  double stripD = 0;
   switch(layer){
      case 0:
        stripD = static_cast<double>(strip);
        break;
      case 1:
        stripD = static_cast<double>(strip) * (1 - 2.29e-3) - 2.412 / 0.25;
        break;        
      case 2:
        stripD = static_cast<double>(strip) * (1 - 8e-3) - 8.46 / 0.25;
        break;
      default:
        stripD = static_cast<double>(strip);
        break;
    }
   return stripD;
}

void filterClusters(vector<map<double, int>> &clusters, int layer, bool vmmHits){
  /* Remove small clusters or clusters with small energy for the first layers only */
  clusters.erase(std::remove_if(clusters.begin(), clusters.end(),
                                [layer, vmmHits](auto c){
                                  if(c.size() > 90) return true;
                                  if(!vmmHits){
                                    if(layer != layerDoubleReadout && c.size() < 3) return true;
                                    else if(layer == layerDoubleReadout && c.size() < 2) return true;
                                  }
                                  return false;}),
                 clusters.end());
  clusters.erase(std::remove_if(clusters.begin(), clusters.end(),
                                [layer, vmmHits](auto c){
                                  auto qthr = 50;
                                  // if(vmmHits)
                                  //   qthr = 50;
                                  // else if(layer == layerDoubleReadout)
                                  //   qthr = 100;          
                                  for(auto &h: c)
                                    if(h.second > qthr)
                                      return false;
                                  return true;
                                }),
                 clusters.end());
}

vector<pair<double, double>> getMeanClusterPositions(map<int, int> hitsPerLayer, int layer, bool vmmHits = false, double singleHitAccuracy = 5){
  map<double, int> hitsPerLayerShifted;
  for(auto &h: hitsPerLayer){
    double strip = correctAlignment(h.first, layer);
    long long pdo = h.second;
    if(!vmmHits && (h.first <= 118 || h.first >= 172)) continue; // TODO check why
    hitsPerLayerShifted.emplace(strip, pdo);
  }

  auto clusters = splitByDistance(hitsPerLayerShifted);
  filterClusters(clusters, layer, vmmHits);
  std::stable_sort(clusters.begin(), clusters.end(), [](auto c1, auto c2){return (c1.begin()->first < c2.begin()->first);});

  vector<pair<double, double>> out;
  for(auto &c: clusters){
    long long sum = 0, sumWeights = 0, nHits = 0;
    optional<double> min = nullopt, max = nullopt;
    for(auto &h: c){
      long long strip = h.first;
      long long pdo = h.second;
      if(!min || min.value() < strip) min = {strip};
      if(!max || max.value() > strip) max = {strip};
      sum += strip * pdo;
      sumWeights += pdo;
      nHits++;
    }
    if(!nHits)
      continue;
    double centerV = static_cast<double>(sum) / static_cast<double>(sumWeights);
    double maxDiff = (nHits == 1) ? singleHitAccuracy : std::max(fabs(centerV-min.value()), fabs(centerV-max.value()));
    out.push_back(make_pair(centerV, maxDiff));
  }
  return out;
}

optional<pair<double, double>> getMeanPosition(map<int, int> hitsPerLayer, int layer, bool vmmHits = false){
  double sum = 0;
  long long sumWeights = 0, nHits = 0;
  optional<double> min = nullopt, max = nullopt;
  for(auto &h: hitsPerLayer){
    if(!vmmHits && (h.first <= 118 || h.first >= 172)) continue; // TODO check why
    double strip = correctAlignment(h.first, layer);
    long long pdo = h.second;
    if(!min || min.value() < strip) min = {strip};
    if(!max || max.value() > strip) max = {strip};
    // if(layer == 0)
    sum += strip * pdo;
    sumWeights += pdo;
    nHits++;
  }
  optional<pair<double, double>> center = nullopt;
  if(nHits){
    double centerV = static_cast<double>(sum) / static_cast<double>(sumWeights);
    center.emplace(make_pair(centerV,
                             (nHits == 1) ? 1.0 : std::max(fabs(centerV-min.value()), fabs(centerV-max.value()))));
  }
  return center;
}

/*
 * positions:
 * L0 - L1: 285
 * L1 - L2: 345
 * L2 - Straw: 523
 * return: position in mm, from Layer 0
 */
int getLayerPosition(int layer){
  int y = 0;
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
  return y;
}

double estimatePositionInLayer(pair<double, double> trackAB, int layer){
  double x = getLayerPosition(layer);
  double y = trackAB.first * x + trackAB.second;
  return y;
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

double getDistanceToTrack(map<int, pair<double, double>> positions, pair<double, double> track){
  double distance = 0;
  for(auto &pos: positions){
    auto distanceToEstimated = pos.second.first - estimatePositionInLayer(track, pos.first);
    distance += distanceToEstimated * distanceToEstimated;
  }
  distance = sqrt(distance);
  return distance;
}

void printCurrentMap(vector<map<int, pair<double, double>>> currentMap){
  for(auto &m: currentMap){
    printf("Map:");
    for(auto &v: m){
      printf("\tLayer %d -- %g (%g)", v.first, v.second.first, v.second.second);
    }
    printf("\n");
  }
}
vector<map<int, pair<double, double>>> getAllClusterCombinations(map<int, vector<pair<double, double>>> positions){
  vector<map<int, pair<double, double>>> positionsMap = {{}};
  vector<map<int, pair<double, double>>> currentMap = {};
  for(auto &pos: positions){
    for(auto &pair: pos.second){
      for(auto m: positionsMap){
        m.emplace(pos.first, pair);
        currentMap.push_back(m);
      }
    }
    positionsMap.clear();
    positionsMap.insert(positionsMap.end(), currentMap.begin(), currentMap.end());
    currentMap.clear();
  }
  return positionsMap;
}
vector<pair<double, double>> getEstimatedTracks(vector<map<int, pair<double, double>>> positions){
  vector<pair<double, double>> out = {};
  std::transform(positions.begin(), positions.end(), std::back_inserter(out), [](auto p){return getEstimatedTrack(p);});
  return out;
}
vector<pair<double, double>> getEstimatedTracks(map<int, vector<pair<double, double>>> positions){
  vector<map<int, pair<double, double>>> positionsMap = getAllClusterCombinations(positions);
  auto out = getEstimatedTracks(positionsMap);
  return out;
}


// From hitsMapper

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
double apvTimeFromSRSTimestamp(int diff, bool fixSRSTime = false){
    return diff * 25.0 * (fixSRSTime ? 400000.0/400037.0 : 1.0) / 1000.0;
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
    return round(apvTimeFromSRSTimestamp(diff, fixSRSTime));
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
