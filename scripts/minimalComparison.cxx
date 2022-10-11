#include <vector>
#include <string>
#include <map>
#include <utility>

#include <memory>

#include <cstdio>

#include "TCanvas.h"
#include "TH1I.h"
#include "THStack.h"
#include "TFile.h"
#include "TH2D.h"

using std::string;
using std::vector;
using std::map;
using std::pair;
using std::shared_ptr, std::make_shared;

void compareBananas(TPad *pad, string filename1, string comment1, string filename2, string comment2, string filename3 = "", string comment3 = ""){
  auto b1 = static_cast<TH2D*>(TFile::Open(filename1.c_str(),"read")->Get("straw31_vs_straw30_banana_ch0"));
  auto b2 = static_cast<TH2D*>(TFile::Open(filename2.c_str(),"read")->Get("straw31_vs_straw30_banana_ch0"));
  auto b3 = static_cast<TH2D*>(TFile::Open(filename3.c_str(),"read")->Get("straw31_vs_straw30_banana_ch0"));
  
  pad->Divide(3);
  auto t1 = pad->cd(1);
  b1->SetTitle(comment1.c_str());
  b1->Draw("colz");
  auto t2 = pad->cd(2);
  b2->SetTitle(comment2.c_str());
  b2->Draw("colz");
  auto t3 = pad->cd(3);
  b3->SetTitle(comment3.c_str());
  b3->Draw("colz");
  pad->cd(0);
  // for(auto &b: {b1, b2, b3}){
  //   b->GetXaxis()->SetRangeUser(-60, 20);
  //   b->GetYaxis()->SetRangeUser(-60, 20);
  // }
  b1->GetXaxis()->SetRangeUser(-70, 10);
  b1->GetYaxis()->SetRangeUser(-70, 10);
  b2->GetXaxis()->SetRangeUser(-80, 00);
  b2->GetYaxis()->SetRangeUser(-80, 00);
  b3->GetXaxis()->SetRangeUser(-60, 20);
  b3->GetYaxis()->SetRangeUser(-60, 20);
  pad->Modified();
  pad->Update();

}

void strawvssci(TPad *pad, string title, string filename1, string comment1, string filename2, string comment2, string filename3 = "", string comment3 = ""){
  auto b1 = static_cast<TH1D*>(TFile::Open(filename1.c_str(),"read")->Get("straw31_vs_sci0"));
  auto b2 = static_cast<TH1D*>(TFile::Open(filename2.c_str(),"read")->Get("straw31_vs_sci0"));
  auto b3 = static_cast<TH1D*>(TFile::Open(filename3.c_str(),"read")->Get("straw31_vs_sci0"));
  b1->SetTitle(comment1.c_str());
  b1->SetLineWidth(2);
  b1->SetLineColor(kSpring - 7);
  b2->SetTitle(comment2.c_str());
  b2->SetLineWidth(2);
  b2->SetLineColor(4);
  b3->SetTitle(comment3.c_str());
  b3->SetLineWidth(2);
  b3->SetLineColor(2);

  for(auto &b: {b1,b2,b3})
    b->Scale(1.0 / b->Integral());
  
  auto stack = new THStack("",(title+";#Deltat;").c_str()); // "Run 007: Straw 31 vs Sci0"
  stack->Add(b1);
  stack->Add(b2);
  stack->Add(b3);
  stack->Draw("nostack h");
  pad->BuildLegend(0.6, 0.7, 0.95, 0.9);
  pad->Modified();
  pad->Update();
}

void strawvsstraw(TPad *pad, string title, string filename1, string comment1, string filename2, string comment2, string filename3 = "", string comment3 = ""){
  auto b1 = static_cast<TH1D*>(TFile::Open(filename1.c_str(),"read")->Get("straw31_vs_straw30"));
  auto b2 = static_cast<TH1D*>(TFile::Open(filename2.c_str(),"read")->Get("straw31_vs_straw30"));
  auto b3 = static_cast<TH1D*>(TFile::Open(filename3.c_str(),"read")->Get("straw31_vs_straw30"));
  b1->SetTitle(comment1.c_str());
  b1->SetLineWidth(2);
  b1->SetLineColor(kSpring - 7);
  b2->SetTitle(comment2.c_str());
  b2->SetLineWidth(2);
  b2->SetLineColor(4);
  b3->SetTitle(comment3.c_str());
  b3->SetLineWidth(2);
  b3->SetLineColor(2);

  for(auto &b: {b1,b2,b3})
    b->Scale(1.0 / b->Integral());
  
  auto stack = new THStack("",(title+";#Deltat;").c_str()); // "Run 007: Straw 31 vs Sci0"
  stack->Add(b1);
  stack->Add(b2);
  stack->Add(b3);
  stack->Draw("nostack h");
  pad->BuildLegend(0.6, 0.7, 0.95, 0.9);
  pad->Modified();
  pad->Update();
}

void scivsmean(TPad *pad, string filename1, string comment1, string filename2, string comment2, string filename3 = "", string comment3 = ""){
  auto b10 = static_cast<TH2D*>(TFile::Open(filename1.c_str(),"read")->Get("sci0_vs_mean_triple"));
  auto b11 = static_cast<TH2D*>(TFile::Open(filename1.c_str(),"read")->Get("sci1_vs_mean_triple"));
  auto b12 = static_cast<TH2D*>(TFile::Open(filename1.c_str(),"read")->Get("sci2_vs_mean_triple"));
  auto b20 = static_cast<TH2D*>(TFile::Open(filename2.c_str(),"read")->Get("sci0_vs_mean_triple"));
  auto b21 = static_cast<TH2D*>(TFile::Open(filename2.c_str(),"read")->Get("sci1_vs_mean_triple"));
  auto b22 = static_cast<TH2D*>(TFile::Open(filename2.c_str(),"read")->Get("sci2_vs_mean_triple"));
  auto b30 = static_cast<TH2D*>(TFile::Open(filename3.c_str(),"read")->Get("sci0_vs_mean_triple"));
  auto b31 = static_cast<TH2D*>(TFile::Open(filename3.c_str(),"read")->Get("sci1_vs_mean_triple"));
  auto b32 = static_cast<TH2D*>(TFile::Open(filename3.c_str(),"read")->Get("sci2_vs_mean_triple"));

  for(auto &b: {b10, b20, b30, b11, b21, b31, b12, b22, b32}){
    b->SetLineWidth(2);
    b->Scale(1.0 / b->Integral());
  }
  for(auto &b: {b10, b20, b30}){
    b->SetLineColor(kSpring - 7);
  }
  for(auto &b: {b11, b21, b31}){
    b->SetLineColor(4);
  }
  for(auto &b: {b12, b22, b32}){
    b->SetLineColor(2);
  }
  
  auto s1 = new THStack("s1",Form("%s: Sci vs mean (triple);#Deltat;", comment1.c_str()));
  s1->Add(b10); s1->Add(b11); s1->Add(b12);
  auto s2 = new THStack("s2",Form("%s: Sci vs mean (triple);#Deltat;", comment2.c_str()));
  s2->Add(b20); s2->Add(b21); s2->Add(b22);
  auto s3 = new THStack("s3",Form("%s: Sci vs mean (triple);#Deltat;", comment3.c_str()));
  s3->Add(b30); s3->Add(b31); s3->Add(b32);
  
  pad->Divide(3);
  auto t1 = pad->cd(1);
  s1->Draw("nostack h");
  t1->BuildLegend(0.6, 0.7, 0.95, 0.85);
  auto t2 = pad->cd(2);
  s2->Draw("nostack h");
  auto t3 = pad->cd(3);
  s3->Draw("nostack h");
  pad->cd(0);
  
  pad->Modified();
  pad->Update();
}

void scivssci(TPad *pad, string filename1, string comment1, string filename2, string comment2, string filename3 = "", string comment3 = ""){
  auto b10 = static_cast<TH2D*>(TFile::Open(filename1.c_str(),"read")->Get("sci0_vs_sci1"));
  auto b11 = static_cast<TH2D*>(TFile::Open(filename1.c_str(),"read")->Get("sci0_vs_sci2"));
  auto b12 = static_cast<TH2D*>(TFile::Open(filename1.c_str(),"read")->Get("sci1_vs_sci2"));
  auto b20 = static_cast<TH2D*>(TFile::Open(filename2.c_str(),"read")->Get("sci0_vs_sci1"));
  auto b21 = static_cast<TH2D*>(TFile::Open(filename2.c_str(),"read")->Get("sci0_vs_sci2"));
  auto b22 = static_cast<TH2D*>(TFile::Open(filename2.c_str(),"read")->Get("sci1_vs_sci2"));
  auto b30 = static_cast<TH2D*>(TFile::Open(filename3.c_str(),"read")->Get("sci0_vs_sci1"));
  auto b31 = static_cast<TH2D*>(TFile::Open(filename3.c_str(),"read")->Get("sci0_vs_sci2"));
  auto b32 = static_cast<TH2D*>(TFile::Open(filename3.c_str(),"read")->Get("sci1_vs_sci2"));

  for(auto &b: {b10, b20, b30, b11, b21, b31, b12, b22, b32}){
    b->SetLineWidth(2);
    b->Scale(1.0 / b->Integral());
  }
  for(auto &b: {b10, b11, b12}){
    b->SetLineColor(kSpring - 7);
    b->SetTitle(comment1.c_str());
  }
  for(auto &b: {b20, b21, b22}){
    b->SetLineColor(4);
    b->SetTitle(comment2.c_str());
  }
  for(auto &b: {b30, b31, b32}){
    b->SetLineColor(2);
    b->SetTitle(comment3.c_str());
  }
  
  auto s1 = new THStack("s1","sci0 vs sci1;#Deltat;");
  s1->Add(b10); s1->Add(b20); s1->Add(b30);
  auto s2 = new THStack("s2","sci0 vs sci2;#Deltat;");
  s2->Add(b11); s2->Add(b21); s2->Add(b31);
  auto s3 = new THStack("s3","sci1 vs sci2;#Deltat;");
  s3->Add(b12); s3->Add(b22); s3->Add(b32);
  
  pad->Divide(3);
  auto t1 = pad->cd(1);
  s1->Draw("nostack h");
  t1->BuildLegend(0.6, 0.7, 0.95, 0.85);
  auto t2 = pad->cd(2);
  s2->Draw("nostack h");
  auto t3 = pad->cd(3);
  s3->Draw("nostack h");
  pad->cd(0);
  
  pad->Modified();
  pad->Update();
}

void scivssci2(TPad *pad, string filename1, string comment1, string filename2, string comment2, string filename3 = "", string comment3 = ""){
  auto b10 = static_cast<TH2D*>(TFile::Open(filename1.c_str(),"read")->Get("sci0_vs_sci1"));
  auto b11 = static_cast<TH2D*>(TFile::Open(filename1.c_str(),"read")->Get("sci0_vs_sci2"));
  auto b12 = static_cast<TH2D*>(TFile::Open(filename1.c_str(),"read")->Get("sci1_vs_sci2"));
  auto b20 = static_cast<TH2D*>(TFile::Open(filename2.c_str(),"read")->Get("sci0_vs_sci1"));
  auto b21 = static_cast<TH2D*>(TFile::Open(filename2.c_str(),"read")->Get("sci0_vs_sci2"));
  auto b22 = static_cast<TH2D*>(TFile::Open(filename2.c_str(),"read")->Get("sci1_vs_sci2"));
  auto b30 = static_cast<TH2D*>(TFile::Open(filename3.c_str(),"read")->Get("sci0_vs_sci1"));
  auto b31 = static_cast<TH2D*>(TFile::Open(filename3.c_str(),"read")->Get("sci0_vs_sci2"));
  auto b32 = static_cast<TH2D*>(TFile::Open(filename3.c_str(),"read")->Get("sci1_vs_sci2"));
  
  for(auto &b: {b10, b20, b30, b11, b21, b31, b12, b22, b32}){
    b->SetLineWidth(2);
    b->Scale(1.0 / b->Integral());
  }
  for(auto &b: {b10, b20, b30}){
    b->SetLineColor(kSpring - 7);
  }
  for(auto &b: {b11, b21, b31}){
    b->SetLineColor(4);
  }
  for(auto &b: {b12, b22, b32}){
    b->SetLineColor(2);
  }

  auto s1 = new THStack("s1",Form("%s: sci vs sci1;#Deltat;", comment1.c_str()));
  s1->Add(b10); s1->Add(b11); s1->Add(b12);
  auto s2 = new THStack("s2",Form("%s: sci vs sci;#Deltat;", comment2.c_str()));
  s2->Add(b20); s2->Add(b21); s2->Add(b22);
  auto s3 = new THStack("s3",Form("%s: sci vs sci;#Deltat;", comment3.c_str()));
  s3->Add(b30); s3->Add(b31); s3->Add(b32);
  
  pad->Divide(3);
  auto t1 = pad->cd(1);
  s1->Draw("nostack h");
  t1->BuildLegend(0.6, 0.7, 0.95, 0.85);
  auto t2 = pad->cd(2);
  s2->Draw("nostack h");
  auto t3 = pad->cd(3);
  s3->Draw("nostack h");
  pad->cd(0);
  
  pad->Modified();
  pad->Update();
}

void minimalComparison(string filename1, string comment1, string filename2, string comment2, string filename3 = "", string comment3 = ""){
  // auto t1 = new TPad("t1", "", 0.0, 1.0, 0.5, 1.0); t1->Draw(); t1->cd();
  vector<string> ext = {"png", "pdf"};
  string outdir = "../out";
  {
    auto c = make_shared<TCanvas>("cmin1", "cmin1", 1200*2, 1200);
    c->Divide(2);
    auto t1 = static_cast<TPad*>(c->cd(1));
    strawvssci(t1, "straw31_vs_sci0", filename1, comment1, filename2, comment2, filename3, comment3);
    auto t2 = static_cast<TPad*>(c->cd(2));
    strawvsstraw(t2, "straw31_vs_straw30", filename1, comment1, filename2, comment2, filename3, comment3);
    c->cd();
    c->Modified();
    c->Update();
      c->SaveAs((string() + outdir + "/" + "minimalComparison" + "1" + "." + "png").c_str());
    // for(auto &e: ext)
    //   c->SaveAs((string() + outdir + "/" + "minimalComparison" + "1" + "." + e).c_str());
  }
  {
    auto c = make_shared<TCanvas>("cmin2", "cmin2", 1200*3, 1200);
    compareBananas(c.get(), filename1, comment1, filename2, comment2, filename3, comment3);
    c->Modified();
    c->Update();
    for(auto &e: ext)
      c->SaveAs((string() + outdir + "/" + "minimalComparison" + "2" + "." + e).c_str());
  }
  // {
  //   auto c = make_shared<TCanvas>("cmin3", "cmin3", 1200*3, 1200*3);
  //   c->Divide(1,3);
  //   auto t1 = static_cast<TPad*>(c->cd(1));
  //   scivsmean(t1, filename1, comment1, filename2, comment2, filename3, comment3);
  //   auto t2 = static_cast<TPad*>(c->cd(2));
  //   scivssci(t2, filename1, comment1, filename2, comment2, filename3, comment3);
  //   auto t3 = static_cast<TPad*>(c->cd(3));
  //   scivssci2(t3, filename1, comment1, filename2, comment2, filename3, comment3);
  //   c->cd();
  //   c->Modified();
  //   c->Update();
  //   for(auto &e: ext)
  //     c->SaveAs((string() + outdir + "/" + "minimalComparison" + "3" + "." + e).c_str());
  //   // c->SaveAs((string() + outdir + "/" + "minimalComparison" + "3" + "." + "png").c_str());
  // }
    {
    auto c = make_shared<TCanvas>("cmin3-1", "cmin3-1", 1200*3, 1200);
    scivsmean(c.get(), filename1, comment1, filename2, comment2, filename3, comment3);
    c->cd();
    c->Modified();
    c->Update();
    for(auto &e: ext)
      c->SaveAs((string() + outdir + "/" + "minimalComparison" + "3-1" + "." + e).c_str());
  }
    {
    auto c = make_shared<TCanvas>("cmin3-2", "cmin3-2", 1200*3, 1200);
    scivssci(c.get(), filename1, comment1, filename2, comment2, filename3, comment3);
    c->cd();
    c->Modified();
    c->Update();
    for(auto &e: ext)
      c->SaveAs((string() + outdir + "/" + "minimalComparison" + "3-2" + "." + e).c_str());
  }
    {
    auto c = make_shared<TCanvas>("cmin3-3", "cmin3-3", 1200*3, 1200);
    scivssci2(c.get(), filename1, comment1, filename2, comment2, filename3, comment3);
    c->cd();
    c->Modified();
    c->Update();
    for(auto &e: ext)
      c->SaveAs((string() + outdir + "/" + "minimalComparison" + "3-3" + "." + e).c_str());
  }

}
