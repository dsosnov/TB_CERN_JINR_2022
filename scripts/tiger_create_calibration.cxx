void tiger_create_calibration(string path, bool directory = false){
  auto chain = new TChain("tigerTL");
  if(directory)
    chain->Add((path + "/*.root").c_str());
  else
    chain->Add(path.c_str());

  chain->SetBranchStatus("*",0);
  chain->SetBranchStatus("gemrocID",1);
  chain->SetBranchStatus("tigerID",1);
  chain->SetBranchStatus("channelID",1);
  chain->SetBranchStatus("tFine",1);
  chain->SetBranchStatus("eFine",1);
  Char_t   gemrocID;        //                "B" == Char_t   ==  int8_t
  Short_t  tigerID;         //  8 bit data -- "S" == Short_t  == int16_t
  Char_t   channelID;       //  6 bit data -- "B" == Char_t   ==  int8_t
  Short_t  tFine;           // 10 bit data -- "S" == Short_t  == int16_t
  Short_t  eFine;           // 10 bit data -- "S" == Short_t  == int16_t
  chain->SetBranchAddress("gemrocID", &gemrocID);
  chain->SetBranchAddress("tigerID", &tigerID);
  chain->SetBranchAddress("channelID", &channelID);
  chain->SetBranchAddress("tFine", &tFine);
  chain->SetBranchAddress("eFine", &eFine);

  map<tuple<Char_t, Short_t, Char_t>, Short_t> tFineMax, tFineMin;
  map<tuple<Char_t, Short_t, Char_t>, Short_t> eFineMax, eFineMin;
  auto nEntries = chain->GetEntries();
  for(auto i = 0; i < nEntries; i++){
    if(!(i%10000000))
      printf("Entry %d of %lld\n", i, nEntries);
    chain->GetEntry(i);
    if(!tFineMax.count({gemrocID, tigerID, channelID}) || tFineMax.at({gemrocID, tigerID, channelID}) < tFine)
      tFineMax[{gemrocID, tigerID, channelID}] = tFine;
    if(!tFineMin.count({gemrocID, tigerID, channelID}) || tFineMin.at({gemrocID, tigerID, channelID}) > tFine)
      tFineMin[{gemrocID, tigerID, channelID}] = tFine;
    if(!eFineMax.count({gemrocID, tigerID, channelID}) || eFineMax.at({gemrocID, tigerID, channelID}) < eFine)
      eFineMax[{gemrocID, tigerID, channelID}] = eFine;
    if(!eFineMin.count({gemrocID, tigerID, channelID}) || eFineMin.at({gemrocID, tigerID, channelID}) > eFine)
      eFineMin[{gemrocID, tigerID, channelID}] = eFine;
  }

  {
    auto fout = fopen("../out/tiger_tfine_calibration.txt", "w");
    fprintf(fout, "# Generated from %s\n", path.c_str());
    for(auto &p: tFineMin){
      fprintf(fout, "%d %d %d %d %d\n", get<0>(p.first), get<1>(p.first), get<2>(p.first), p.second, tFineMax.at(p.first));
    }
    fclose(fout);
  }
  {
    auto fout = fopen("../out/tiger_efine_calibration.txt", "w");
    fprintf(fout, "# Generated from %s\n", path.c_str());
    for(auto &p: eFineMin){
      fprintf(fout, "%d %d %d %d %d\n", get<0>(p.first), get<1>(p.first), get<2>(p.first), p.second, eFineMax.at(p.first));
    }
    fclose(fout);
  }
}
