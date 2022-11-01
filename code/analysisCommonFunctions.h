#pragma once

#include "analysisGeneral.h"

#include <vector>
#include <string>
#include <map>
#include <optional>

using std::vector;
using std::string;
using std::map;
using std::optional, std::nullopt;

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

/*
 * positions (april measurements):
 * L0 - L1: 285
 * L1 - L2: 345
 * L2 - L3 ( & Straw): 523
 *
 * Some distances (July measurements):
 * 183 - mm0-straw (layer -1 now)
 * 434 - mm0-mm1 - mm layer 3 (y) was connected to mm layer 1
 * 325 - mm3-mm2
 *
 * return: position in mm, from Layer 0
 * Layer -1 -> straw, odd(?)
 * Layer -2 -> straw, even(?)
 */
int getLayerPosition(int layer, analysisGeneral::TestBeams tb){ // TODO do normally
  int y = 0;
  switch(tb){
    case analysisGeneral::TestBeams::TB22_October:
      switch(layer){
        case 2:
          y += 195;
        case -2:
          y += 7;
        case -1:
          y += 233;
        case 3:
          y += 32;
        case 1:
          y += 332;
        case 0:
          y += 0;
      }
      break;
    case analysisGeneral::TestBeams::TB22_August:
      break;
    case analysisGeneral::TestBeams::TB22_July:
      switch(layer){
        case -2:
          y += 7;
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
      break;
    case analysisGeneral::TestBeams::TB22_April:
      switch(layer){
        case -2:
          y += 7;
        case -1:
        case 3:
          y += 523;
        case 2:
          y += 345;
        case 1:
          y += 285;
        case 0:
          y += 0;
      }
      break;
    default:
      break;
  }
  return y;
}

double correctAlignment(int strip, int layer, analysisGeneral::TestBeams tb){
  double stripD = static_cast<double>(strip);
  switch(tb){
    case analysisGeneral::TestBeams::TB22_October:
      break;
    case analysisGeneral::TestBeams::TB22_August:
      break;
    case analysisGeneral::TestBeams::TB22_July:
      switch(layer){
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
vector<pair<double, double>> getMeanClusterPositions(vector<pair<int, int>> hitsPerLayer, int layer, bool vmmHits = false, double singleHitAccuracy = 5){
  map<int, int> hitsMap(hitsPerLayer.begin(), hitsPerLayer.end());
  return getMeanClusterPositions(hitsMap, layer, vmmHits, singleHitAccuracy);
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

double getDistanceToTrack(map<int, pair<double, double>> positions, pair<double, double> track){
  double distance = 0;
  for(auto &pos: positions){
    auto distanceToEstimated = pos.second.first - estimatePositionInLayer(track, pos.first);
    distance += distanceToEstimated * distanceToEstimated;
  }
  distance = sqrt(distance);
  return distance;
}
