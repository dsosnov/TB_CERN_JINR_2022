#include "tiger.C"
#include <vector>

void process_tiger(string filename, string directory, string mapFile, string calibFile, std::vector<short> eModes){
  auto t = new tiger(filename.c_str(), directory.c_str(), mapFile.c_str(), calibFile.c_str(), eModes);
  t->Loop();
  delete t;
}
