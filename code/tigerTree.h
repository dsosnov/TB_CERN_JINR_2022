#pragma once

#include "Rtypes.h"
#include <cstdio>
#include <utility>
#include <tuple>

// For simplisity of calculation use signed-sized types with increased size.
struct tigerHit {
  Char_t   gemrocID;        //                "B" == Char_t   ==  int8_t
  Short_t  tigerID;         //  8 bit data -- "S" == Short_t  == int16_t
  Char_t   channelID;       //  6 bit data -- "B" == Char_t   ==  int8_t
  virtual void print(bool hex = false) const;
  virtual double stepSize_fC(double VcaspVth) const;
  virtual double convertEFineSHTofC(Short_t eFine, double saturationValue, Short_t vcasp, Short_t thrEFine, Short_t maximumEFine) const;
};
void tigerHit::print(bool hex) const {
  printf("TL hit: ");
  if(hex)
    printf("[ROC %2X, TIGER %2X, ch %2X] ", gemrocID, tigerID, channelID);
  else
    printf("[ROC %3d, TIGER %3d, ch %2d] ", gemrocID, tigerID, channelID);
  printf("\n");
};
double inline tigerHit::stepSize_fC(double VcaspVth) const {
  double gain = 12.25;
  return (VcaspVth * -0.621 + 39.224) / gain;
};
double inline tigerHit::convertEFineSHTofC(Short_t eFine, double saturationValue, Short_t vcasp, Short_t thrEFine, Short_t maximumEFine) const {
  if(eFine > 1007)
    return convertEFineSHTofC(eFine - 1024, saturationValue, vcasp, thrEFine, maximumEFine);
  auto thrStep = stepSize_fC(vcasp);
  return saturationValue - (static_cast<double>(eFine) / static_cast<double>(maximumEFine)) * (thrStep * thrEFine);
};

struct tigerHitTL : public tigerHit {
  Char_t   chipID;          //  2 bit data -- "B" == Char_t   ==  int8_t
  /* 4 TAC per shaper for event de-randomization */
  Char_t   tacID;           //  2 bit data -- "B" == Char_t   ==  int8_t
  Int_t    tCoarse;         // 16 bit data -- "I" == Int_t    == int32_t
  Short_t  eCoarse;         // 10 bit data -- "S" == Short_t  == int16_t
  Short_t  tFine;           // 10 bit data -- "S" == Short_t  == int16_t
  Short_t  eFine;           // 10 bit data -- "S" == Short_t  == int16_t

  Int_t    frameCount;      // 16 bit data -- "I" == Int_t    == int32_t
  Int_t    seu;             // 16 bit data -- "I" == Int_t    == int32_t
  Long64_t frameCountLoops; //                "L" == Long64_t == int64_t

  Int_t    counterWord;     // 24 bit data -- "I" == Int_t    == int32_t

  std::pair<Short_t, Short_t> tFineLimits = {0, 1023};
  std::pair<Short_t, Short_t> eFineLimits = {0, 1023};
  std::tuple<Double_t, Short_t, Short_t, Short_t> eFineCalibrationSH = {45.0, 55, 63, 1007};
  double tFineCorrected() const;
  double eFineCorrected() const;

  double timeFine() const; // ns
  double chargeToT(const bool fine = true) const;
  double chargeSH() const;
  double charge(const bool ToTMode = false) const;
  void print(bool hex = false) const override;

  friend bool isLater(const tigerHitTL hit1, const tigerHitTL hit2);
  friend bool isLater(const tigerHitTL *hit1, const tigerHitTL *hit2);
  friend Long64_t timeDifferenceCoarsePS(const tigerHitTL hit1, const tigerHitTL hit2);
  friend Long64_t timeDifferenceCoarsePS(const tigerHitTL *hit1, const tigerHitTL *hit2);
  friend double timeDifferenceFineNS(const tigerHitTL hit1, const tigerHitTL hit2);
  friend double timeDifferenceFineNS(const tigerHitTL *hit1, const tigerHitTL *hit2);
};

inline double tigerHitTL::tFineCorrected() const {
  if(tFine < tFineLimits.first){
    // printf("Hit tFine below limits! ");
    // print();
    return 0;
  }else if(tFine > tFineLimits.second){
    // printf("Hit tFine above limits! ");
    // print();
    return 1.0;
  }else{
    return static_cast<double>(tFine - tFineLimits.first) / static_cast<double>(tFineLimits.second - tFineLimits.first + 1);
  }
}

inline double tigerHitTL::eFineCorrected() const {
  if(eFine < eFineLimits.first){
    printf("Hit eFine below limits! ");
    print();
    return 0;
  }else if(eFine > eFineLimits.second){
    printf("Hit eFine above limits! ");
    print();
    return 1.0;
  }else{
    return static_cast<double>(eFine - eFineLimits.first) / static_cast<double>(eFineLimits.second - eFineLimits.first + 1);
  }
}

double tigerHitTL::timeFine() const { // ns
  double tCounts = double(tCoarse) - tFineCorrected();
  return tCounts * 6.25;
}

double tigerHitTL::chargeToT(const bool fine) const {
  double ediff = eCoarse - tCoarse%0x400;
  ediff += (ediff >= 0) ? 0 : 1024;
  if(fine)
    ediff -= eFineCorrected() - tFineCorrected();
  return ediff;    
}

double tigerHitTL::chargeSH() const { // in fC
  return convertEFineSHTofC(eFine,
                   std::get<0>(eFineCalibrationSH),
                   std::get<1>(eFineCalibrationSH),
                   std::get<2>(eFineCalibrationSH),
                   std::get<3>(eFineCalibrationSH));
}

double tigerHitTL::charge(const bool ToTMode) const {
  if(ToTMode)
    return chargeToT(true);
  else
    return chargeSH();
}


void tigerHitTL::print(bool hex) const {
  printf("TL hit: ");

  if(hex)
    printf("[ROC %2X, TIGER %2X, ch %2X] ", gemrocID, tigerID, channelID);
  else
    printf("[ROC %3d, TIGER %3d, ch %2d] ", gemrocID, tigerID, channelID);

  if(hex){
    printf("tCoarse: %4X, tFine: %3X, eCoarse: %3X, eFine: %3X; ", tCoarse, tFine, eCoarse, eFine);
    printf("FrameCount: %4X, FC loops: %lld", frameCount, frameCountLoops);
  }else{
    printf("tCoarse: %5d, tFine: %4d, eCoarse: %4d, eFine: %4d; ", tCoarse, tFine, eCoarse, eFine);
    printf("FrameCount: %5d, FC loops: %lld", frameCount, frameCountLoops);
  }
  printf("\n");
}


bool isLater(const tigerHitTL *hit1, const tigerHitTL *hit2){
  bool later = false;
  if(hit1->frameCountLoops != hit2->frameCountLoops){
    later = hit1->frameCountLoops > hit2->frameCountLoops;
  }else if(hit1->frameCount != hit2->frameCount){
    later = hit1->frameCount > hit2->frameCount;
  }else if(std::abs(hit1->tCoarse - hit2->tCoarse) > (1<<15)){ // If difference between hits > 2^15, then there is roll-over between
    later = hit1->tCoarse < hit2->tCoarse;
  }else if(hit1->tCoarse != hit2->tCoarse){
    later = hit1->tCoarse > hit2->tCoarse;
  }else{
    later = hit1->tFineCorrected() < hit2->tFineCorrected();
  }
  return later;
}
bool isLater(const tigerHitTL hit1, const tigerHitTL hit2){
  return isLater(&hit1, &hit2);
}

/*The differance should be not larger then ~106 days */
Long64_t timeDifferenceCoarsePS(const tigerHitTL *hit1, const tigerHitTL *hit2){
  auto later = isLater(hit1, hit2);
  auto hitFirst = later ? hit2: hit1;
  auto hitLast = later ? hit1: hit2;

  // each framecount loop is 2**16 frameCounts.
  Long64_t FCLoopDiff = hitLast->frameCountLoops - hitFirst->frameCountLoops;
  if(FCLoopDiff < 0){
    fprintf(stderr, "FrameCountLoop roll-over should never happens (3E12 years) !\n");
    FCLoopDiff += 9223372036854775807;  // max value of Long64_t = 2^63-1
    FCLoopDiff += 1;
  }

  // each two framecount is whole tCoarse roll-over loop.
  Long64_t FCDiff = FCLoopDiff * 65536;
  FCDiff += hitLast->frameCount - hitFirst->frameCount;
  // Roll-overs: int(floor((FC2_{Last} - FC_{First})/2))
  Long64_t tCoarseRollOvers = (FCDiff - (FCDiff%2 ? 1 : 0)) / 2;
  Long64_t tCoarseDiff = hitLast->tCoarse - hitFirst->tCoarse;
  if(FCDiff%2){ // odd FCDiff
    if(tCoarseDiff < 0)
      tCoarseRollOvers++;
  } else { // even FCDiff
    if(std::abs(tCoarseDiff) > (1<<15)){
      if(tCoarseDiff < 0)
        tCoarseDiff += (1<<16);
      else
        tCoarseDiff -= (1<<16);
    }
  }
  Long64_t clk_periods = tCoarseDiff + tCoarseRollOvers * (1<<16);
  Long64_t absTime = clk_periods * 6250;

  return later ? absTime : -absTime; // picoseconds
}
Long64_t timeDifferenceCoarsePS(const tigerHitTL hit1, const tigerHitTL hit2){
  return timeDifferenceCoarsePS(&hit1, &hit2);
}

double timeDifferenceFineNS(const tigerHitTL *hit1, const tigerHitTL *hit2){
  auto later = isLater(hit1, hit2);
  auto hitFirst = later ? hit2: hit1;
  auto hitLast = later ? hit1: hit2;
  double fineDiff = hitFirst->tFineCorrected() - hitLast->tFineCorrected();

  double timeDiffAbs = double(timeDifferenceCoarsePS(*hitLast, *hitFirst)) / 1E3 + fineDiff * 6.25;
  
  return later ? timeDiffAbs : -timeDiffAbs; // nanoseconds
}
double timeDifferenceFineNS(const tigerHitTL hit1, const tigerHitTL hit2){
  return timeDifferenceFineNS(&hit1, &hit2);
}
