#include "tiger.C"

void process_tiger(string filename, string directory, string mapFile, string calibFile){
  auto t = new tiger(filename.c_str(), directory.c_str(), mapFile.c_str(), calibFile.c_str());
  t->Loop();
  delete t;
}
