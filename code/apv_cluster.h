#pragma once

#include <TObject.h>

#include <vector> 
#include <map>
#include <string>
#include <memory>
#include <set> 

#include <numeric> // accumulate
#include <algorithm> // max_element, sort

using std::vector;

struct apvHit{
  int layer = 0;
  int strip = 0;
  short max_q = 0;
  int t_max_q = 0;
  vector<short> raw_q = {};
  void print(bool verbose = false) const {
    printf("Hit on layer %d, strip %d with max Q %d at time bin %d", layer, strip, max_q, t_max_q);
    if(verbose){
      printf("; Q values: [");
      for(auto i = 0; i < raw_q.size(); i++){
        printf("%d", raw_q.at(i));
        if(i < raw_q.size()-1)
          printf(", ");        
      }
      printf("]");
    }
    printf("\n");
  }
  friend bool operator==( apvHit const& a, apvHit const& b ){
    return (a.layer == b.layer) &&
      (a.strip == b.strip) &&
      (a.max_q == b.max_q) &&
      (a.t_max_q == b.t_max_q) &&
      (a.raw_q == b.raw_q);
  }
  friend bool operator!=( apvHit const& a, apvHit const& b ){
    return !(a == b);
  }
  friend bool operator<( apvHit const& a, apvHit const& b ){
    if(a.layer != b.layer)
      return a.layer < b.layer;
    else if(a.strip != b.strip)
      return a.strip < b.strip;
    else if(a.max_q != b.max_q)
      return a.max_q < b.max_q;
    else if(a.t_max_q != b.t_max_q)
      return a.t_max_q < b.t_max_q;
    else
      return false;
  }
};

class apvCluster: public TObject{
private:
  int layer;
  vector<apvHit> hits;
  unsigned long sizeOnLastUpdate;
  float center_;
  // int width_;
  int maxq_;
  long qsum_;
  int maxqtime_;
  float meanqtime_;
  void update(){
    if(sizeOnLastUpdate && sizeOnLastUpdate == nHits())
      return;
    sortHits();
    calcQ();
    calcMaxQ();
    calcCenter();
    calcMaxQTime();
    calcMeanQTime();
    sizeOnLastUpdate = nHits();
  }
  void sortHits(){
    if(sizeOnLastUpdate && sizeOnLastUpdate == nHits())
      return;
    std::sort(hits.begin(), hits.end(), [](const apvHit h1, const apvHit h2){return (h1.strip < h2.strip);});
  }
  void calcCenter(){
    if(sizeOnLastUpdate && sizeOnLastUpdate == nHits())
      return;
    long long sum = 0, sumw = 0;
    for(auto &hit: hits){
      sum += hit.strip * hit.max_q;
      sumw += hit.max_q;
    }
    center_ = static_cast<float>(sum) / static_cast<float>(sumw);
  }
  void calcMaxQ(){
    if(sizeOnLastUpdate && sizeOnLastUpdate == nHits())
      return;
    maxq_ = -1;
    for(auto &hit: hits)
      if(hit.max_q > maxq_)
        maxq_ = hit.max_q;
  }
  void calcQ(){
    if(sizeOnLastUpdate && sizeOnLastUpdate == nHits())
      return;
    qsum_ = 0;
    for(auto &hit: hits)
      qsum_ += hit.max_q;
  }
  void calcMaxQTime(){
    if(sizeOnLastUpdate && sizeOnLastUpdate == nHits())
      return;
    maxqtime_ = -1;
    for(auto &hit: hits){
      if(hit.max_q == maxq_)
        maxqtime_ = hit.t_max_q;
    }
  }
  void calcMeanQTime(){
    if(sizeOnLastUpdate && sizeOnLastUpdate == nHits())
      return;
    meanqtime_ = -1;
    long long sum = 0, sumw = 0;
    for(auto &hit: hits){
      sum += hit.t_max_q * hit.max_q;
      sumw += hit.max_q;
    }
    meanqtime_ = static_cast<float>(sum) / static_cast<float>(sumw);
  }
public:
  apvCluster(int clusterLayer = 0):
    layer(clusterLayer), hits({}), sizeOnLastUpdate(0){
    update();
  }
  apvCluster(apvHit hit):
    layer(hit.layer), hits({}), sizeOnLastUpdate(0){
    addHit(hit);
    update();
  }
  ~apvCluster(){
    hits.clear();
    update();
  }
  int getLayer() const {return layer;}
  vector<apvHit> getHits() const {return hits;}
  bool addHit(apvHit hit){
    if(hit.layer != layer)
      return false;
    hits.push_back(hit);
    update();
    return true;
  }
  bool addHitAdjacent(apvHit hit, int neighborAccepted = 3){
    if(hit.layer != layer)
      return false;
    if(hit.strip < hits.at(0).strip - neighborAccepted || hit.strip >= hits.back().strip + neighborAccepted)
      return false;
    addHit(hit);
    return true;
  }
  unsigned long nHits() const { return hits.size(); }
  void print(bool verbose = false) const {
    printf("Cluster: %lu hits on layer %d, with center %.2f, width %d and Q %ld (maximal: %d)\n", nHits(), layer, center(), width(), q(), maxQ());
    if(verbose)
      for(auto &h: hits)
        h.print(true);
  }
  bool merge(apvCluster cluster){
    if(cluster.layer != layer)
      return false;
    if(cluster.lastStrip() != hits.at(0).strip - 1 && cluster.firstStrip() != hits.back().strip + 1)
      return false;
    hits.insert(hits.end(), cluster.hits.begin(), cluster.hits.end());
    update();
    return true;
  }

  float center() const { return center_; }
  int width() const { return hits.back().strip - hits.at(0).strip + 1; }
  int firstStrip() const { return hits.at(0).strip; }
  int lastStrip() const { return hits.back().strip; }
  int maxQ() const { return maxq_; }
  long q() const { return qsum_; }
  int maxQTime() const { return maxqtime_; }
  float meanQTime() const { return meanqtime_; }

  friend
  auto operator<( apvCluster const& a, apvCluster const& b ){
    if(a.getLayer() != b.getLayer())
      return a.getLayer() < b.getLayer();
    else if(a.hits.size() != b.hits.size())
      return a.hits.size() < b.hits.size();
    else{
      auto hitsA = a.getHits();
      auto hitsB = b.getHits();
      sort(hitsA.begin(), hitsA.end());
      sort(hitsB.begin(), hitsB.end());
      for(auto i = 0; i < hitsA.size(); i++){
        if(hitsA.at(i) != hitsB.at(i)){
          return hitsA.at(i) < hitsB.at(i);
        }
      }
      return false;
    }
  }
  friend
  auto operator==( apvCluster const& a, apvCluster const& b ){
    if(a.getLayer() != b.getLayer())
      return false;
    else
      return a.hits == b.hits;
  }

};

class apvTrack : public TObject{
private:
  double x0, b;
  std::set<apvCluster> clusters;
  std::set<apvCluster> clustersY; // TODO
public:
  apvTrack(double intersect = 0, double slope = 0): x0(intersect), b(slope), clusters({}) {};
  template <class T>
  apvTrack(double intersect, double slope, T clusters_): x0(intersect), b(slope){
    for(auto &c: clusters_)
      addCluster(c);
  };
  void addCluster(apvCluster cluster){clusters.emplace(cluster);}
  unsigned long nClusters(){return clusters.size();}
    std::set<apvCluster> getClusters(){return clusters;}
  double intersect() const {return x0;}
  double slope() const {return b;}
  void setIntersect(double intersect){x0 = intersect;}
  void setSlope(double slope){b = slope;}
  int maxQ() const {
    int maxq = 0;
    for(auto &c: clusters){
      if(c.maxQ() > maxq)
        maxq = c.maxQ();
    }
    return maxq;
  }
  bool isX2() const {
    for(auto &c: clusters)
      if(c.getLayer() == 2)
        return true;
    return false;
  }
  apvCluster* getX2Cluster() const {
    if(!isX2())
      return nullptr;
    for(auto &&c: clusters){
      if(c.getLayer() == 2)
        return const_cast<apvCluster*>(&c);
    }
    return nullptr;
  }
  float meanQTime() const {
    float meanqtime = 0;
    int nhits = 0;
    for(auto &c: clusters){
      if(!isX2() || (isX2() && c.getLayer() == 2)){
        meanqtime += c.meanQTime();
        nhits++;
      }
    }
    return meanqtime / nhits;
  }
  float maxQTime() const {
    float maxqtime = 0;
    int nhits = 0;
    for(auto &c: clusters){
      if(!isX2() || (isX2() && c.getLayer() == 2)){
        maxqtime += c.maxQTime();
        nhits++;
      }
    }
    return maxqtime / nhits;
  }
};
