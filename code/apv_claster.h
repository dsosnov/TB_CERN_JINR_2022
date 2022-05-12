#pragma once

#include <vector> 
#include <map>
#include <string>
#include <memory>

#include <numeric> // accumulate
#include <algorithm> // max_element, sort

struct apvHit{
  int layer = 0;
  int strip = 0;
  short max_q;
  int t_max_q;
  vector<short> raw_q;
};

class apvClaster{
private:
  int layer;
  vector<apvHit> hits;
  unsigned long sizeOnLastUpdate;
  // float center_;
  // int width_;
  int maxq_;
  long qsum_;
public:
  apvClaster(int clasterLayer = 0):
    layer(clasterLayer), hits({}), sizeOnLastUpdate(-1){
  }
  apvClaster(apvHit hit):
    layer(hit.layer), hits({}), sizeOnLastUpdate(-1){
    addHit(hit);
  }
  ~apvClaster(){
    hits.clear();
  }
  int getLayer(){return layer;}
  vector<apvHit> getHits(){return hits;}
  bool addHit(apvHit hit){
    if(hit.layer != layer)
      return false;
    hits.push_back(hit);
    sortHits();
    return true;
  }
  bool addHitAdjacent(apvHit hit){
    if(hit.layer != layer)
      return false;
    if(hit.strip != hits.at(0).strip - 1 && hit.strip != hits.back().strip + 1)
      return false;
    addHit(hit);
    return true;
  }
  float center(){
    sortHits();
    return (hits.at(0).strip + hits.back().strip)/2.0;
  }
  int width(){
    sortHits();
    return hits.back().strip - hits.at(0).strip + 1;
  }
  int firstStrip(){ sortHits(); return hits.at(0).strip; }
  int lastStrip(){ sortHits(); return hits.back().strip; }
  int maxQ(){
    if(!sizeOnLastUpdate && sizeOnLastUpdate == nHits())
      return maxq_;
    maxq_ = -1;
    for(auto &hit: hits)
      if(hit.max_q > maxq_)
        maxq_ = hit.max_q;
    return maxq_;
  }
  long q(){
    if(!sizeOnLastUpdate && sizeOnLastUpdate == nHits())
      return qsum_;    
    qsum_ = 0;
    for(auto &hit: hits)
      qsum_ += hit.max_q;
    return qsum_;
  }
  unsigned long nHits(){ return hits.size(); }
  void print(){
    printf("Claster: %lu hits on layer %d, with center %.2f, width %d and Q %ld (maximal: %d)\n", nHits(), layer, center(), width(), q(), maxQ()); 
  }
  void sortHits(){
    if(!sizeOnLastUpdate && sizeOnLastUpdate == nHits())
      return;
    std::sort(hits.begin(), hits.end(), [](const apvHit h1, const apvHit h2){return (h1.strip < h2.strip);});
    sizeOnLastUpdate = nHits();
  }
  bool merge(apvClaster claster){
    if(claster.layer != layer)
      return false;
    if(claster.lastStrip() != hits.at(0).strip - 1 && claster.firstStrip() != hits.back().strip + 1)
      return false;
    hits.insert(hits.end(), claster.hits.begin(), claster.hits.end());
    sortHits();
    return true;
  }
};

