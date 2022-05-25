void printNoise(string file="run_0358.root"){
  auto chain = new TChain("vmm");
  chain->Add(file.c_str());

  // static auto eventFAFAIn = new int;
  // static auto triggerTimeStampIn = new vector<int>;
  // static auto triggerCounterIn = new vector<int>;
  // static auto boardIdIn = new vector<int>;
  // static auto chipIn = new vector<int>;
  // static auto eventSizeIn = new vector<int>;
  static auto daq_timestamp_sIn = new vector<int>;
  static auto daq_timestamp_nsIn = new vector<int>;
  // static auto tdoIn = new vector<vector<int>>;
  // static auto pdoIn = new vector<vector<int>>;
  // static auto flagIn = new vector<vector<int>>;
  // static auto thresholdIn = new vector<vector<int>>;
  // static auto bcidIn = new vector<vector<int>>;
  // static auto relbcidIn = new vector<vector<int>>;
  // static auto overflowIn = new vector<vector<int>>;
  // static auto orbitCountIn = new vector<vector<int>>;
  // static auto grayDecodedIn = new vector<vector<int>>;
  static auto channelIn = new vector<vector<int>>;
  // static auto febChannelIn = new vector<vector<int>>;
  // static auto mappedChannelIn = new vector<vector<int>>;
  // static auto art_validIn = new vector<int>;
  // static auto artIn = new vector<int>;
  // static auto art_triggerIn = new vector<int>;
  {
    // chain->SetBranchAddress("eventFAFA", &eventFAFAIn);
    // chain->SetBranchAddress("triggerTimeStamp", &triggerTimeStampIn);
    // chain->SetBranchAddress("triggerCounter", &triggerCounterIn);
    // chain->SetBranchAddress("chip", &chipIn);
    // chain->SetBranchAddress("boardId", &boardIdIn);
    // chain->SetBranchAddress("eventSize", &eventSizeIn);
    chain->SetBranchAddress("daq_timestamp_s", &daq_timestamp_sIn);
    chain->SetBranchAddress("daq_timestamp_ns", &daq_timestamp_nsIn);
    // chain->SetBranchAddress("tdo", &tdoIn);
    // chain->SetBranchAddress("pdo", &pdoIn);
    // chain->SetBranchAddress("flag", &flagIn);
    // chain->SetBranchAddress("threshold", &thresholdIn);
    // chain->SetBranchAddress("bcid", &bcidIn);
    // chain->SetBranchAddress("relbcid", &relbcidIn);
    // chain->SetBranchAddress("overflow", &overflowIn);
    // chain->SetBranchAddress("orbitCount", &orbitCountIn);
    // chain->SetBranchAddress("grayDecoded", &grayDecodedIn);
    chain->SetBranchAddress("channel", &channelIn);
    // chain->SetBranchAddress("febChannel", &febChannelIn);
    // chain->SetBranchAddress("mappedChannel", &mappedChannelIn);
    // chain->SetBranchAddress("art_valid", &art_validIn);
    // chain->SetBranchAddress("art", &artIn);
    // chain->SetBranchAddress("art_trigger", &art_triggerIn);
  }

  chain->Draw("channel >> hchannel(64, 0, 64)", "channel != 0 && channel != 60");
  auto hchannel = static_cast<TH1F*>(gDirectory->Get("hchannel"));

  unsigned long long daqtimeDiff = 0;
  chain->GetEntry(0);
  daqtimeDiff = daq_timestamp_sIn->at(0) * 1E9 + daq_timestamp_nsIn->at(0);
  chain->GetEntry(chain->GetEntries()-1);
  auto nhits = daq_timestamp_sIn->size();
  daqtimeDiff = daq_timestamp_sIn->at(nhits-1) * 1E9 + daq_timestamp_nsIn->at(nhits-1) - daqtimeDiff;

  double daqtimesecdiff = static_cast<double>(daqtimeDiff) / 1E9;
  
  hchannel->Scale(1.0 / daqtimesecdiff);

  printf("File %s\n", file.c_str());
  printf("Time difference %g sec:\n", daqtimesecdiff);  
  printf("Noise levels :\n");
  for(auto i = 1; i <= hchannel->GetNbinsX(); i++){
    int channel = static_cast<int>(TMath::Floor(hchannel->GetXaxis()->GetBinLowEdge(i)));
    printf(" channel %d - %.2g Hz (± %.2g)\n", channel,
           hchannel->GetBinContent(i),
           hchannel->GetBinError(i));
  }
  printf("\n");
}
