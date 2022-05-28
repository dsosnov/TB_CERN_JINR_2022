#include "apv.C"
#include "evBuilder.C"

void merge_vmm_apv(){
  pair<string, string> run_pair = {"run_0258", "run166"};

  // unsigned long long from = 1653145627, to = 1653145627;
  unsigned long long from = 0, to = 0;

  auto apvan = new apv(run_pair.second);
  auto hits_apv = apvan->GetCentralHits(from, to);
  auto vmman = new evBuilder(run_pair.first, "g1_p25_s100-0&60", "map-20220518");
  vmman->Loop();
  auto hits_vmm = vmman->GetCentralHits(from, to);
  {
    double pdoMax = 0.0;
    for(auto &h: hits_vmm)
      if(h.pdoRelative > pdoMax)
        pdoMax = h.pdoRelative;
    for(auto &h: hits_vmm)
      h.pdoRelative /= pdoMax;          
  }
  {
    double pdoMax = 0.0;
    for(auto &h: hits_apv)
      if(h.pdoRelative > pdoMax)
        pdoMax = h.pdoRelative;
    for(auto &h: hits_apv)
      h.pdoRelative /= pdoMax;          
  }

  printf("APV hits (%lu)\n", hits_apv.size());
  for(auto &h: hits_apv) h.print();
  printf("VMM hits (%lu):\n", hits_vmm.size());
  for(auto &h: hits_vmm) h.print();
}
