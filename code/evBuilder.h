#ifndef evBuilder_h
#define evBuilder_h

#include "vmm.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class evBuilder : public vmm {
public :
   evBuilder(TString filename) : vmm(filename) {} ;
   evBuilder(TTree *tree = 0) : vmm(tree) {};
   virtual ~evBuilder() {};

   virtual void     Loop() override;

   void threePlotDrawF(TH1D *h1, TH1D *h2, TH1D *h3);
};

#endif
