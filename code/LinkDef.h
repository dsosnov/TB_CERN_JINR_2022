#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class vector<float>+;
#pragma link C++ class vector<vector<float>>+;
#pragma link C++ class vector<double>+;
#pragma link C++ class vector<vector<double>>+;
#pragma link C++ class vector<int>+;
#pragma link C++ class vector<vector<int>>+;
#pragma link C++ class vector<short>+;
#pragma link C++ class vector<vector<short>>+;
#pragma link C++ class vector<unsigned int>+;
#pragma link C++ class vector<vector<unsigned int>>+;
#pragma link C++ class vector<char>+;
#pragma link C++ class vector<vector<char>>+;
#pragma link C++ class vector<string>+;
#pragma link C++ class vector<vector<string>>+;

#pragma link C++ class vector<apvHit>+;
#pragma link C++ class vector<vector<apvHit>>+;
#pragma link C++ class vector<apvCluster>+;
#pragma link C++ class vector<vector<apvCluster>>+;
#pragma link C++ class vector<apvTrack>+;
#pragma link C++ class vector<vector<apvTrack>>+;

#pragma link C++ class map<int, int>+;
#pragma link C++ class map<int, map<int, int>>+;

#pragma link C++ class analysisGeneral::mm2CenterHitParameters+;
#pragma link C++ class apv::doubleReadoutHits+;

#pragma link C++ class map<int, pair<double, double>>+;
#pragma link C++ class map<int, double>+;

#pragma link C++ class mmCluster+;
#pragma link C++ class vector<mmCluster>+;

#endif
