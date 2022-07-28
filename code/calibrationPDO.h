#pragma once

#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2I.h"
#include "TFitResult.h"

#include <vector>
#include <string>
#include <memory>
#include <utility>

#include <tuple>
#include <map>

#include <algorithm>
#include <functional>

using std::vector;
using std::string;
using std::shared_ptr, std::make_shared;
using std::pair;
using std::tuple;
using std::map;

// https://ned.ipac.caltech.edu/level5/Leo/Stats4_5.html
pair<double, double> weightedMean(vector<pair<double, double>> values){
  double sumValues = 0, sumWeights = 0;
  for(auto &v: values){
    auto value = v.first;
    if(!value) continue;
    auto sigma = v.second;
    auto w = 1 / (sigma * sigma);
    sumValues += value * w;
    sumWeights += w;
  }
  double mu, sigma;
  if(sumWeights){
    mu = sumValues / sumWeights;
    sigma = sqrt( 1.0 / sumWeights);
  } else {
    mu = sumValues / values.size();
    sigma = 0;
  }
  return {mu, sigma};
}

struct singleChannelSingleTest{
  int channelNum;
  double pdoMean, pdoMin, pdoMax;
  double pdoMaxDiff() const {return std::max(pdoMax - pdoMean,pdoMean - pdoMin);}
  unsigned long nHits;
};
struct singleTest{
  string testName;
  double mean, meanE;
  vector<singleChannelSingleTest> channels;
  TH2I* pdoVsChannel;
  TH2I* pdoVsChannelWOTHR;
};

singleTest findCalibrationForChannels(const shared_ptr<TChain> chain,
                                      const string calibrationFileName,
                                      const unsigned int nChannels = 64,
                                      const vector<unsigned int> channelsExclude = {},
                                      const vector<unsigned int> channelsDelete = {},
                                      const double pdoThr = 120){
  static auto eventFAFAIn = new int;
  static auto triggerTimeStampIn = new vector<int>;
  static auto triggerCounterIn = new vector<int>;
  static auto boardIdIn = new vector<int>;
  static auto chipIn = new vector<int>;
  static auto eventSizeIn = new vector<int>;
  static auto daq_timestamp_sIn = new vector<int>;
  static auto daq_timestamp_nsIn = new vector<int>;
  static auto tdoIn = new vector<vector<int>>;
  static auto pdoIn = new vector<vector<int>>;
  static auto flagIn = new vector<vector<int>>;
  static auto thresholdIn = new vector<vector<int>>;
  static auto bcidIn = new vector<vector<int>>;
  static auto relbcidIn = new vector<vector<int>>;
  static auto overflowIn = new vector<vector<int>>;
  static auto orbitCountIn = new vector<vector<int>>;
  static auto grayDecodedIn = new vector<vector<int>>;
  static auto channelIn = new vector<vector<int>>;
  static auto febChannelIn = new vector<vector<int>>;
  static auto mappedChannelIn = new vector<vector<int>>;
  static auto art_validIn = new vector<int>;
  static auto artIn = new vector<int>;
  static auto art_triggerIn = new vector<int>;
  {
    chain->SetBranchAddress("eventFAFA", &eventFAFAIn);
    chain->SetBranchAddress("triggerTimeStamp", &triggerTimeStampIn);
    chain->SetBranchAddress("triggerCounter", &triggerCounterIn);
    chain->SetBranchAddress("chip", &chipIn);
    chain->SetBranchAddress("boardId", &boardIdIn);
    chain->SetBranchAddress("eventSize", &eventSizeIn);
    chain->SetBranchAddress("daq_timestamp_s", &daq_timestamp_sIn);
    chain->SetBranchAddress("daq_timestamp_ns", &daq_timestamp_nsIn);
    chain->SetBranchAddress("tdo", &tdoIn);
    chain->SetBranchAddress("pdo", &pdoIn);
    chain->SetBranchAddress("flag", &flagIn);
    chain->SetBranchAddress("threshold", &thresholdIn);
    chain->SetBranchAddress("bcid", &bcidIn);
    chain->SetBranchAddress("relbcid", &relbcidIn);
    chain->SetBranchAddress("overflow", &overflowIn);
    chain->SetBranchAddress("orbitCount", &orbitCountIn);
    chain->SetBranchAddress("grayDecoded", &grayDecodedIn);
    chain->SetBranchAddress("channel", &channelIn);
    chain->SetBranchAddress("febChannel", &febChannelIn);
    chain->SetBranchAddress("mappedChannel", &mappedChannelIn);
    chain->SetBranchAddress("art_valid", &art_validIn);
    chain->SetBranchAddress("art", &artIn);
    chain->SetBranchAddress("art_trigger", &art_triggerIn);
  }

  if(!chain->GetEntries()){
    printf("No events in file! Exiting...\n");
    exit(0);
  }

  vector<singleChannelSingleTest> channels;
  for(auto i = 0; i < nChannels; i++){
    channels.push_back({i, 0, 0, 0, 0});
  }
  auto pdoPerChannelWOTHR = new TH2I(Form("pdo_per_channel_%s_wo_thr", calibrationFileName.c_str()), Form("pdo_per_channel_%s_wo_thr", calibrationFileName.c_str()), 64, 0, 64, 1024, 0, 1024);
  auto pdoPerChannel = new TH2I(Form("pdo_per_channel_%s", calibrationFileName.c_str()), Form("pdo_per_channel_%s", calibrationFileName.c_str()), 64, 0, 64, 1024, 0, 1024);
  // fill pdo information for channel
  for(auto i = 0; i < chain->GetEntries(); i++){
    chain->GetEntry(i);
    for(unsigned long triggerN = 0; triggerN < triggerCounterIn->size(); triggerN++){
      for(unsigned long j = 0; j < channelIn->at(triggerN).size(); j++){
        auto channel = channelIn->at(triggerN).at(j);
        if(std::find(channelsExclude.begin(), channelsExclude.end(), channel) != channelsExclude.end())
          continue;
        if(std::find(channelsDelete.begin(), channelsDelete.end(), channel) != channelsDelete.end())
          continue;
        pdoPerChannelWOTHR->Fill(channel, pdoIn->at(triggerN).at(j));
        if(pdoIn->at(triggerN).at(j) < pdoThr)
          continue;
        pdoPerChannel->Fill(channel, pdoIn->at(triggerN).at(j));

        if(!channels.at(channel).nHits || pdoIn->at(triggerN).at(j) < channels.at(channel).pdoMin)
          channels.at(channel).pdoMin = pdoIn->at(triggerN).at(j);
        if(!channels.at(channel).nHits || pdoIn->at(triggerN).at(j) > channels.at(channel).pdoMax)
          channels.at(channel).pdoMax = pdoIn->at(triggerN).at(j);
        channels.at(channel).pdoMean += pdoIn->at(triggerN).at(j);
        channels.at(channel).nHits++;
      }
    }
  }
  for(auto &c: channels) if(c.nHits) c.pdoMean /= c.nHits;

  std::function<vector<singleChannelSingleTest>(vector<singleChannelSingleTest>)> filterOut = [](vector<singleChannelSingleTest> channels){
    vector<singleChannelSingleTest> channelsFiltered;
    std::copy_if (channels.begin(), channels.end(), std::back_inserter(channelsFiltered), [](singleChannelSingleTest c){return c.nHits;} );
    return channelsFiltered;
  };
  channels = filterOut(channels);

  vector<pair<double, double>> pdoMean;
  for(auto i = 0; i < channels.size(); i++){
    pdoMean.push_back({channels.at(i).pdoMean, channels.at(i).pdoMaxDiff()});
  }
  auto meanValue = weightedMean(pdoMean);

  return {calibrationFileName, meanValue.first, meanValue.second, channels, pdoPerChannel, pdoPerChannelWOTHR};
}

// Data for correcting pdo to the mean one
struct fitPDOChannelResult{
  int nChannel;
  // [p0] + [p1]*x
  double p0, p0E;
  double p1, p1E;
  double chi2;
  unsigned int ndf;
  double correctPDO(const double pdo) const {return p0 + p1 * pdo;}
  long correctPDOI(const double pdo) const {return static_cast<long>(correctPDO(pdo));}
  void print(bool compact = false, FILE* file = stdout) const {
    if(compact)
      fprintf(file, "%.2g\t%.2g\n", p0, p1);
    else
      fprintf(file, "fitPDOChannelResult channel %d: p0 = %.2g pm %.2g; p1 = %.2g pm %g [chi2/ndf: %g/%u]\n", nChannel, p0, p0E, p1, p1E, chi2, ndf);
  }
};
vector<fitPDOChannelResult> fitPDO(const vector<singleTest> valuesForFit){
  map<int, TH1F*> fitData;
  for(auto &test: valuesForFit){
    for(auto &channel: test.channels){
      if(!fitData.count(channel.channelNum))
        fitData.emplace(channel.channelNum,
                       new TH1F(Form("h_calibrate_ch%i", channel.channelNum), Form("h_calibrate_ch%i", channel.channelNum), 1024, 0, 1024));
      fitData.at(channel.channelNum)->SetBinContent(test.mean, channel.pdoMean);
      fitData.at(channel.channelNum)->SetBinError(test.mean, channel.pdoMaxDiff());
    }
  }

  
  vector<fitPDOChannelResult> out;
  for(auto &fitPair: fitData){
    // printf("Channel %d fitting ...\n", fitPair.first);
    auto h = fitPair.second;
    if(h->GetEntries() < 2){
      out.push_back({fitPair.first, 0, 0, 1.0, 0, 0, 0});
      continue;
    }
    auto fitR = h->Fit("pol1", "S");
    auto b = fitR->Parameter(0);
    auto bE = fitR->ParError(0);
    auto a = fitR->Parameter(1);
    auto aE = fitR->ParError(1);
    auto chi2 = fitR->Chi2();
    auto ndf = fitR->Ndf();
    out.push_back({fitPair.first, -b / a, bE / a, 1.0 / a, aE / a / a, chi2, ndf});
    // break;
  }
    
  return out;
}

void applyPDOCorrections(const shared_ptr<TChain> chain,
                         const string calibrationFileName,
                         vector<fitPDOChannelResult> corrections){
  static auto triggerCounterIn = new vector<int>;
  static auto pdoIn = new vector<vector<int>>;
  static auto channelIn = new vector<vector<int>>;
  {
    chain->SetBranchAddress("triggerCounter", &triggerCounterIn);
    chain->SetBranchAddress("pdo", &pdoIn);
    chain->SetBranchAddress("channel", &channelIn);
  }

  if(!chain->GetEntries()){
    printf("No events in file! Exiting...\n");
    exit(0);
  }

  auto pdoPerChannel = new TH2F(Form("corrected_pdo_per_channel_%s", calibrationFileName.c_str()),
                                Form("corrected_pdo_per_channel_%s", calibrationFileName.c_str()), 64, 0, 64, 1024, 0, 1024);
  for(auto i = 0; i < chain->GetEntries(); i++){
    chain->GetEntry(i);
    for(unsigned long triggerN = 0; triggerN < triggerCounterIn->size(); triggerN++){
      for(unsigned long j = 0; j < channelIn->at(triggerN).size(); j++){
        auto channel = channelIn->at(triggerN).at(j);
        auto energy = pdoIn->at(triggerN).at(j);
        for(auto &c: corrections)
          if(c.nChannel == channel){
            // printf("PDO in channel %d was: %i", channel, energy);
            energy = c.correctPDOI(energy);
            // printf(" -> becomes %i\n", energy);
            // c.print();
            break;
          }
        pdoPerChannel->Fill(channel, energy);
      }
    }
  }
  
}
