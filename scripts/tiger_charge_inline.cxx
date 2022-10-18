// (eCoarse>=tCoarse&0x3FF)*(eCoarse-tCoarse&0x3FF)+(eCoarse<tCoarse&0x3FF)*(eCoarse+0x400-tCoarse&0x3FF)
// eCoarse>=(tCoarse&0x3FF) + (eCoarse<tCoarse&0x3FF)*0x400
// tCoarse&0x3ff is the same as tCoarse%1024
// double tiger_charge_inline(){
  // if(eCoarse >= tCoarse%0x400)
  //   return eCoarse - tCoarse%0x400;
  // else
  //   return eCoarse + 0x400 - tCoarse%0x400;
// }

// (eCoarse-tCoarse%1024) - (eFine / 1024.0 - tFine / 1024.0) + (eCoarse<tCoarse%1024)*1024.0
double tiger_charge_inline(){
  double ediff = eCoarse - tCoarse%0x400;
  ediff += (ediff >= 0) ? 0 : 1024;
  ediff -= (eFine / 1024.0 - tFine / 1024.0);
  return ediff;
}
