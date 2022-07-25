// .L ../scripts/tiger_pdo.cxx
// tiger_pdo(tCoarse,eCoarse,tFine,eFine)
// the same as
// (eCoarse-tCoarse%1024) - (eFine / 1024.0 - tFine / 1024.0) + (eCoarse<tCoarse%1024)*1024.0
double tiger_pdo(Int_t tCoarse, Short_t eCoarse, Short_t tFine, Short_t eFine){
  double ediff = eCoarse - tCoarse%0x400;
  ediff += (ediff >= 0) ? 0 : 1024;
  ediff -= (eFine / 1024.0 - tFine / 1024.0);
  return ediff;
}
