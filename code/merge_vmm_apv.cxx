#include "apv.C"
#include "evBuilder.C"

void merge_vmm_apv(){
  pair<string, string> run_pair = {"run_0258", "run166"};
  string run_apv = "run166";
  string run_vmm = "run_0258";

  auto apvan = new apv(run_pair.second);
  auto hits_apv = apvan->GetCentralHits(1653145627, 1653145627);
  auto vmman = new evBuilder(run_pair.first, "g1_p25_s100-0&60", "map-20220518");
  auto hits_vmm = vmman->GetCentralHits(1653145627, 1653145627);

  printf("APV hits:\n");
  for(auto &h: hits_apv) h.print();
  printf("VMM hits:\n");
  for(auto &h: hits_vmm) h.print();
}
