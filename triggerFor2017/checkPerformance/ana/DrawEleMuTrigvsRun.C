/***************************************
drawing trigger efficiency vs run period
for dilepton and trilepton
*****************************************/

#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"

#include <iostream>

using namespace std;

void GraphStyle1(TGraphAsymmErrors *g1) {
  g1->SetMarkerStyle(20);
  g1->SetMarkerColor(2);
  g1->SetLineColor(2);
}

void GraphStyle2(TGraphAsymmErrors *g1) {
  g1->SetMarkerStyle(20);
  g1->SetMarkerColor(kGreen+2);
  g1->SetLineColor(kGreen+2);
}

void GraphStyle3(TGraphAsymmErrors *g1) {
  g1->SetMarkerStyle(21);
  g1->SetMarkerColor(4);
  g1->SetLineColor(4);
}

void GraphStyle4(TGraphAsymmErrors *g1) {
  g1->SetMarkerStyle(21);
  g1->SetMarkerColor(kOrange);
  g1->SetLineColor(kOrange);
}

void GraphStyle5(TGraphAsymmErrors *g1) {
  g1->SetMarkerStyle(22);
  g1->SetMarkerColor(4);
  g1->SetLineColor(4);
}

void GraphStyle6(TGraphAsymmErrors *g1) {
  g1->SetMarkerStyle(22);
  g1->SetMarkerColor(kOrange);
  g1->SetLineColor(kOrange);
}

void DrawEleMuTrigvsRun(TString leg = "ele") {

  gStyle->SetOptStat(1111);

  cout << "starting program" << endl;
  TString charname = "output/DiEleMu_TriEleMu";
  TString suf;
  if (leg.Contains("ele") ) suf = "_eleLeg_EleEndcap";
  else suf = "_muonLeg";

  cout << "read input file" << endl;
  //TFile *frunD = new TFile(charname + "_2017D_PromptReco_v1_eleLeg_TagMuMatchMuLeg_EleBarrel.root", "READ");
  TFile *frunD = new TFile(charname + "_2017D_PromptReco_v1" + suf + ".root", "READ");
  TFile *frunE = new TFile(charname + "_2017E_PromptReco_v1" + suf + ".root", "READ");
  TFile *frunF = new TFile(charname + "_2017F_PromptReco_v1" + suf + ".root", "READ");

  //file of Ele23_Ele12 trigger
  TFile *frunD_diele = new TFile("output/DiEle_Ele23_Ele12_2017D_PromptReco_v1.root", "READ");
  TFile *frunE_diele = new TFile("output/DiEle_Ele23_Ele12_2017E_PromptReco_v1.root", "READ");
  TFile *frunF_diele = new TFile("output/DiEle_Ele23_Ele12_2017F_PromptReco_v1.root", "READ");


  cout << "define tgraph" << endl;

  //-----for Mu23_Ele12-------//
  //ele leg
  TGraphAsymmErrors *gElePt_Mu23Ele12_2017D;
  TGraphAsymmErrors *gElePt_Mu23Ele12_2017E;
  TGraphAsymmErrors *gElePt_Mu23Ele12_2017F;

  TGraphAsymmErrors *gEleEta_Mu23Ele12_2017D;
  TGraphAsymmErrors *gEleEta_Mu23Ele12_2017E;
  TGraphAsymmErrors *gEleEta_Mu23Ele12_2017F;

  TGraphAsymmErrors *gElePhi_Mu23Ele12_2017D;
  TGraphAsymmErrors *gElePhi_Mu23Ele12_2017E;
  TGraphAsymmErrors *gElePhi_Mu23Ele12_2017F;

  TGraphAsymmErrors *gEleVtx_Mu23Ele12_2017D;
  TGraphAsymmErrors *gEleVtx_Mu23Ele12_2017E;
  TGraphAsymmErrors *gEleVtx_Mu23Ele12_2017F;

  //muon leg
  TGraphAsymmErrors *gMuPt_Mu23Ele12_2017D;
  TGraphAsymmErrors *gMuPt_Mu23Ele12_2017E;
  TGraphAsymmErrors *gMuPt_Mu23Ele12_2017F;

  TGraphAsymmErrors *gMuEta_Mu23Ele12_2017D;
  TGraphAsymmErrors *gMuEta_Mu23Ele12_2017E;
  TGraphAsymmErrors *gMuEta_Mu23Ele12_2017F;

  TGraphAsymmErrors *gMuPhi_Mu23Ele12_2017D;
  TGraphAsymmErrors *gMuPhi_Mu23Ele12_2017E;
  TGraphAsymmErrors *gMuPhi_Mu23Ele12_2017F;

  TGraphAsymmErrors *gMuVtx_Mu23Ele12_2017D;
  TGraphAsymmErrors *gMuVtx_Mu23Ele12_2017E;
  TGraphAsymmErrors *gMuVtx_Mu23Ele12_2017F;

  //----for Mu12_Ele23----//
  //ele leg
  TGraphAsymmErrors *gElePt_Mu12Ele23_2017D;
  TGraphAsymmErrors *gElePt_Mu12Ele23_2017E;
  TGraphAsymmErrors *gElePt_Mu12Ele23_2017F;

  TGraphAsymmErrors *gEleEta_Mu12Ele23_2017D;
  TGraphAsymmErrors *gEleEta_Mu12Ele23_2017E;
  TGraphAsymmErrors *gEleEta_Mu12Ele23_2017F;

  TGraphAsymmErrors *gElePhi_Mu12Ele23_2017D;
  TGraphAsymmErrors *gElePhi_Mu12Ele23_2017E;
  TGraphAsymmErrors *gElePhi_Mu12Ele23_2017F;

  TGraphAsymmErrors *gEleVtx_Mu12Ele23_2017D;
  TGraphAsymmErrors *gEleVtx_Mu12Ele23_2017E;
  TGraphAsymmErrors *gEleVtx_Mu12Ele23_2017F;

  //muon leg
  TGraphAsymmErrors *gMuPt_Mu12Ele23_2017D;
  TGraphAsymmErrors *gMuPt_Mu12Ele23_2017E;
  TGraphAsymmErrors *gMuPt_Mu12Ele23_2017F;

  TGraphAsymmErrors *gMuEta_Mu12Ele23_2017D;
  TGraphAsymmErrors *gMuEta_Mu12Ele23_2017E;
  TGraphAsymmErrors *gMuEta_Mu12Ele23_2017F;

  TGraphAsymmErrors *gMuPhi_Mu12Ele23_2017D;
  TGraphAsymmErrors *gMuPhi_Mu12Ele23_2017E;
  TGraphAsymmErrors *gMuPhi_Mu12Ele23_2017F;

  TGraphAsymmErrors *gMuVtx_Mu12Ele23_2017D;
  TGraphAsymmErrors *gMuVtx_Mu12Ele23_2017E;
  TGraphAsymmErrors *gMuVtx_Mu12Ele23_2017F;

  //efficiency of Ele23 leg from Ele23_Ele12
  TGraphAsymmErrors *gEle_EE_leg1_2017D = (TGraphAsymmErrors*) frunD_diele->Get("gEle_EE_leg1");
  TGraphAsymmErrors *gEle_EE_leg1_2017E = (TGraphAsymmErrors*) frunE_diele->Get("gEle_EE_leg1");
  TGraphAsymmErrors *gEle_EE_leg1_2017F = (TGraphAsymmErrors*) frunF_diele->Get("gEle_EE_leg1");

  TGraphAsymmErrors *gEle_EE_leg2_2017D = (TGraphAsymmErrors*) frunD_diele->Get("gEle_EE_leg2");
  TGraphAsymmErrors *gEle_EE_leg2_2017E = (TGraphAsymmErrors*) frunE_diele->Get("gEle_EE_leg2");
  TGraphAsymmErrors *gEle_EE_leg2_2017F = (TGraphAsymmErrors*) frunF_diele->Get("gEle_EE_leg2");

  GraphStyle4(gEle_EE_leg1_2017D);
  GraphStyle4(gEle_EE_leg1_2017E);
  GraphStyle4(gEle_EE_leg1_2017F);

  GraphStyle4(gEle_EE_leg2_2017D);
  GraphStyle4(gEle_EE_leg2_2017E);
  GraphStyle4(gEle_EE_leg2_2017F);

  cout << "read graph from input file" << endl;
  //barrel
  //for Mu23_Ele12
  //ele leg
  gElePt_Mu23Ele12_2017D = (TGraphAsymmErrors*) frunD->Get("gelePt_Mu23Ele12");
  gElePt_Mu23Ele12_2017E = (TGraphAsymmErrors*) frunE->Get("gelePt_Mu23Ele12");
  gElePt_Mu23Ele12_2017F = (TGraphAsymmErrors*) frunF->Get("gelePt_Mu23Ele12");

  gEleEta_Mu23Ele12_2017D = (TGraphAsymmErrors*) frunD->Get("geleEta_Mu23Ele12");
  gEleEta_Mu23Ele12_2017E = (TGraphAsymmErrors*) frunE->Get("geleEta_Mu23Ele12");
  gEleEta_Mu23Ele12_2017F = (TGraphAsymmErrors*) frunF->Get("geleEta_Mu23Ele12");

  gElePhi_Mu23Ele12_2017D = (TGraphAsymmErrors*) frunD->Get("gelePhi_Mu23Ele12");
  gElePhi_Mu23Ele12_2017E = (TGraphAsymmErrors*) frunE->Get("gelePhi_Mu23Ele12");
  gElePhi_Mu23Ele12_2017F = (TGraphAsymmErrors*) frunF->Get("gelePhi_Mu23Ele12");

  gEleVtx_Mu23Ele12_2017D = (TGraphAsymmErrors*) frunD->Get("gnVtx_Mu23Ele12_eleleg");
  gEleVtx_Mu23Ele12_2017E = (TGraphAsymmErrors*) frunE->Get("gnVtx_Mu23Ele12_eleleg");
  gEleVtx_Mu23Ele12_2017F = (TGraphAsymmErrors*) frunF->Get("gnVtx_Mu23Ele12_eleleg");

  //muon leg
  gMuPt_Mu23Ele12_2017D = (TGraphAsymmErrors*) frunD->Get("gmuPt_Mu23Ele12");
  gMuPt_Mu23Ele12_2017E = (TGraphAsymmErrors*) frunE->Get("gmuPt_Mu23Ele12");
  gMuPt_Mu23Ele12_2017F = (TGraphAsymmErrors*) frunF->Get("gmuPt_Mu23Ele12");

  gMuEta_Mu23Ele12_2017D = (TGraphAsymmErrors*) frunD->Get("gmuEta_Mu23Ele12");
  gMuEta_Mu23Ele12_2017E = (TGraphAsymmErrors*) frunE->Get("gmuEta_Mu23Ele12");
  gMuEta_Mu23Ele12_2017F = (TGraphAsymmErrors*) frunF->Get("gmuEta_Mu23Ele12");

  gMuPhi_Mu23Ele12_2017D = (TGraphAsymmErrors*) frunD->Get("gmuPhi_Mu23Ele12");
  gMuPhi_Mu23Ele12_2017E = (TGraphAsymmErrors*) frunE->Get("gmuPhi_Mu23Ele12");
  gMuPhi_Mu23Ele12_2017F = (TGraphAsymmErrors*) frunF->Get("gmuPhi_Mu23Ele12");

  gMuVtx_Mu23Ele12_2017D = (TGraphAsymmErrors*) frunD->Get("gnVtx_Mu23Ele12_muleg");
  gMuVtx_Mu23Ele12_2017E = (TGraphAsymmErrors*) frunE->Get("gnVtx_Mu23Ele12_muleg");
  gMuVtx_Mu23Ele12_2017F = (TGraphAsymmErrors*) frunF->Get("gnVtx_Mu23Ele12_muleg");

  //------for Mu12_Ele23--------//
  // ele leg
  gElePt_Mu12Ele23_2017D = (TGraphAsymmErrors*) frunD->Get("gelePt_Ele23Mu12");
  gElePt_Mu12Ele23_2017E = (TGraphAsymmErrors*) frunE->Get("gelePt_Ele23Mu12");
  gElePt_Mu12Ele23_2017F = (TGraphAsymmErrors*) frunF->Get("gelePt_Ele23Mu12");

  gEleEta_Mu12Ele23_2017D = (TGraphAsymmErrors*) frunD->Get("geleEta_Ele23Mu12");
  gEleEta_Mu12Ele23_2017E = (TGraphAsymmErrors*) frunE->Get("geleEta_Ele23Mu12");
  gEleEta_Mu12Ele23_2017F = (TGraphAsymmErrors*) frunF->Get("geleEta_Ele23Mu12");

  gElePhi_Mu12Ele23_2017D = (TGraphAsymmErrors*) frunD->Get("gelePhi_Ele23Mu12");
  gElePhi_Mu12Ele23_2017E = (TGraphAsymmErrors*) frunE->Get("gelePhi_Ele23Mu12");
  gElePhi_Mu12Ele23_2017F = (TGraphAsymmErrors*) frunF->Get("gelePhi_Ele23Mu12");

  gEleVtx_Mu12Ele23_2017D = (TGraphAsymmErrors*) frunD->Get("gnVtx_Ele2323Mu12_eleleg");
  gEleVtx_Mu12Ele23_2017E = (TGraphAsymmErrors*) frunE->Get("gnVtx_Ele2323Mu12_eleleg");
  gEleVtx_Mu12Ele23_2017F = (TGraphAsymmErrors*) frunF->Get("gnVtx_Ele2323Mu12_eleleg");

  //muon leg
  gMuPt_Mu12Ele23_2017D = (TGraphAsymmErrors*) frunD->Get("gmuPt_Ele23Mu12");
  gMuPt_Mu12Ele23_2017E = (TGraphAsymmErrors*) frunE->Get("gmuPt_Ele23Mu12");
  gMuPt_Mu12Ele23_2017F = (TGraphAsymmErrors*) frunF->Get("gmuPt_Ele23Mu12");

  gMuEta_Mu12Ele23_2017D = (TGraphAsymmErrors*) frunD->Get("gmuEta_Ele23Mu12");
  gMuEta_Mu12Ele23_2017E = (TGraphAsymmErrors*) frunE->Get("gmuEta_Ele23Mu12");
  gMuEta_Mu12Ele23_2017F = (TGraphAsymmErrors*) frunF->Get("gmuEta_Ele23Mu12");

  gMuPhi_Mu12Ele23_2017D = (TGraphAsymmErrors*) frunD->Get("gmuPhi_Ele23Mu12");
  gMuPhi_Mu12Ele23_2017E = (TGraphAsymmErrors*) frunE->Get("gmuPhi_Ele23Mu12");
  gMuPhi_Mu12Ele23_2017F = (TGraphAsymmErrors*) frunF->Get("gmuPhi_Ele23Mu12");

  gMuVtx_Mu12Ele23_2017D = (TGraphAsymmErrors*) frunD->Get("gnVtx_Ele2323Mu12_muleg");
  gMuVtx_Mu12Ele23_2017E = (TGraphAsymmErrors*) frunE->Get("gnVtx_Ele2323Mu12_muleg");
  gMuVtx_Mu12Ele23_2017F = (TGraphAsymmErrors*) frunF->Get("gnVtx_Ele2323Mu12_muleg");

  cout << "graph style" << endl;
  //style for graph
  //for Mu23_Ele12
  GraphStyle1(gElePt_Mu23Ele12_2017D);
  GraphStyle2(gElePt_Mu23Ele12_2017E);
  GraphStyle3(gElePt_Mu23Ele12_2017F);
  
  GraphStyle1(gEleEta_Mu23Ele12_2017D);
  GraphStyle2(gEleEta_Mu23Ele12_2017E);
  GraphStyle3(gEleEta_Mu23Ele12_2017F);

  GraphStyle1(gElePhi_Mu23Ele12_2017D);
  GraphStyle2(gElePhi_Mu23Ele12_2017E);
  GraphStyle3(gElePhi_Mu23Ele12_2017F);

  GraphStyle1(gEleVtx_Mu23Ele12_2017D);
  GraphStyle2(gEleVtx_Mu23Ele12_2017E);
  GraphStyle3(gEleVtx_Mu23Ele12_2017F);

  GraphStyle1(gMuPt_Mu23Ele12_2017D);
  GraphStyle2(gMuPt_Mu23Ele12_2017E);
  GraphStyle3(gMuPt_Mu23Ele12_2017F);
  
  GraphStyle1(gMuEta_Mu23Ele12_2017D);
  GraphStyle2(gMuEta_Mu23Ele12_2017E);
  GraphStyle3(gMuEta_Mu23Ele12_2017F);

  GraphStyle1(gMuPhi_Mu23Ele12_2017D);
  GraphStyle2(gMuPhi_Mu23Ele12_2017E);
  GraphStyle3(gMuPhi_Mu23Ele12_2017F);

  GraphStyle1(gMuVtx_Mu23Ele12_2017D);
  GraphStyle2(gMuVtx_Mu23Ele12_2017E);
  GraphStyle3(gMuVtx_Mu23Ele12_2017F);

  //for Mu12_Ele23
  GraphStyle1(gElePt_Mu12Ele23_2017D);
  GraphStyle1(gElePt_Mu12Ele23_2017E);
  GraphStyle1(gElePt_Mu12Ele23_2017F);
  
  GraphStyle1(gEleEta_Mu12Ele23_2017D);
  GraphStyle2(gEleEta_Mu12Ele23_2017E);
  GraphStyle3(gEleEta_Mu12Ele23_2017F);

  GraphStyle1(gElePhi_Mu12Ele23_2017D);
  GraphStyle2(gElePhi_Mu12Ele23_2017E);
  GraphStyle3(gElePhi_Mu12Ele23_2017F);

  GraphStyle1(gEleVtx_Mu12Ele23_2017D);
  GraphStyle2(gEleVtx_Mu12Ele23_2017E);
  GraphStyle3(gEleVtx_Mu12Ele23_2017F);
  
  GraphStyle1(gMuPt_Mu12Ele23_2017D);
  GraphStyle2(gMuPt_Mu12Ele23_2017E);
  GraphStyle3(gMuPt_Mu12Ele23_2017F);
  
  GraphStyle1(gMuEta_Mu12Ele23_2017D);
  GraphStyle2(gMuEta_Mu12Ele23_2017E);
  GraphStyle3(gMuEta_Mu12Ele23_2017F);

  GraphStyle1(gMuPhi_Mu12Ele23_2017D);
  GraphStyle2(gMuPhi_Mu12Ele23_2017E);
  GraphStyle3(gMuPhi_Mu12Ele23_2017F);

  GraphStyle1(gMuVtx_Mu12Ele23_2017D);
  GraphStyle2(gMuVtx_Mu12Ele23_2017E);
  GraphStyle3(gMuVtx_Mu12Ele23_2017F);
  
  cout << "done graph style" << endl;

  TLatex tx;
  tx.SetNDC(kTRUE);
  tx.SetTextFont(52);
  tx.SetTextSize(0.045);

  TLine *l1 = new TLine(0., 1., 500., 1.);
  l1->SetLineStyle(6);
  l1->SetLineColor(kGreen);
  l1->SetLineWidth(2);

  TLine *l2 = new TLine(-2.5, 1., 2.5, 1.);
  l2->SetLineStyle(6);
  l2->SetLineColor(kGreen);
  l2->SetLineWidth(2);

  cout << "starting canvas for plotting" << endl;
  cout << "graph of Mu23_Ele12 path" << endl;

  /*
  TCanvas* cpt_ele = new TCanvas("cpt_ele", "cpt_ele", 650, 650);
  cpt_ele->cd();
  gElePt_Mu23Ele12_2017F->GetYaxis()->SetTitle("trigger efficiency");
  gElePt_Mu23Ele12_2017F->GetXaxis()->SetTitle("p_{T}^{e} (GeV/c)");
  //gElePt_Mu23Ele12_2017F->GetXaxis()->SetLimits(0, 200.);
  gElePt_Mu23Ele12_2017F->GetYaxis()->SetRangeUser(0., 1.2);
  gElePt_Mu23Ele12_2017F->Draw("APZ");
  gElePt_Mu23Ele12_2017D->Draw("PZ");
  gElePt_Mu23Ele12_2017E->Draw("PZ");
  TLegend *lg = new TLegend(0.35, 0.26, 0.75, 0.40);
  lg->AddEntry(gElePt_Mu23Ele12_2017D, "Run2017D_PromptReco_v1", "pl");
  lg->AddEntry(gElePt_Mu23Ele12_2017E, "Run2017E_PromptReco_v1", "pl");
  lg->AddEntry(gElePt_Mu23Ele12_2017F, "Run2017F_PromptReco_v1", "pl");
  //lg->AddEntry(gElePt_EB_2017B_v2_leg1, "Run2017B_PromptReco_v2", "pl");
  //lg->AddEntry(gElePt_EB_2017C_v1_leg1, "Run2017C_PromptReco_v1", "pl");
  //lg->AddEntry(gElePt_EB_2017C_v2_leg1, "Run2017C_PromptReco_v2", "pl");
  lg->SetLineColor(0);
  lg->SetFillColor(0);
  lg->SetShadowColor(0);
  lg->SetTextSize(0.035);
  lg->SetTextFont(42);
  lg->Draw();
  l1->Draw();
  tx.DrawLatex(0.45, 0.6 ,"EleLeg (HLT_Mu23_Ele12)"); 

  TCanvas* ceta_ele = new TCanvas("ceta_ele", "ceta_ele", 650, 650);
  ceta_ele->cd();
  gEleEta_Mu23Ele12_2017F->GetYaxis()->SetTitle("trigger efficiency");
  gEleEta_Mu23Ele12_2017F->GetXaxis()->SetTitle("#eta^{e}");
  //gEleEta_Mu23Ele12_2017F->GetXaxis()->SetLimits(0, 200.);
  gEleEta_Mu23Ele12_2017F->GetYaxis()->SetRangeUser(0., 1.2);
  gEleEta_Mu23Ele12_2017F->Draw("APZ");
  gEleEta_Mu23Ele12_2017D->Draw("PZ");
  gEleEta_Mu23Ele12_2017E->Draw("PZ");
  lg->Draw();
  l2->Draw();
  tx.DrawLatex(0.45, 0.6 ,"EleLeg (HLT_Mu23_Ele12)"); 

  TCanvas* cphi_ele = new TCanvas("cphi_ele", "cphi_ele", 650, 650);
  cphi_ele->cd();
  gElePhi_Mu23Ele12_2017F->GetYaxis()->SetTitle("trigger efficiency");
  gElePhi_Mu23Ele12_2017F->GetXaxis()->SetTitle("#phi^{e}");
  //gElePhi_Mu23Ele12_2017F->GetXaxis()->SetLimits(0, 200.);
  gElePhi_Mu23Ele12_2017F->GetYaxis()->SetRangeUser(0., 1.2);
  gElePhi_Mu23Ele12_2017F->Draw("APZ");
  gElePhi_Mu23Ele12_2017D->Draw("PZ");
  gElePhi_Mu23Ele12_2017E->Draw("PZ");
  lg->Draw();
  l2->Draw();
  tx.DrawLatex(0.45, 0.6 ,"EleLeg (HLT_Mu23_Ele12)"); 

  TCanvas* cvtx_ele = new TCanvas("cvtx_ele", "cvtx_ele", 650, 650);
  cvtx_ele->cd();
  gEleVtx_Mu23Ele12_2017F->GetYaxis()->SetTitle("trigger efficiency");
  gEleVtx_Mu23Ele12_2017F->GetXaxis()->SetTitle("nVtx");
  //gEleVtx_Mu23Ele12_2017F->GetXaxis()->SetLimits(0, 200.);
  gEleVtx_Mu23Ele12_2017F->GetYaxis()->SetRangeUser(0., 1.2);
  gEleVtx_Mu23Ele12_2017F->Draw("APZ");
  gEleVtx_Mu23Ele12_2017D->Draw("PZ");
  gEleVtx_Mu23Ele12_2017E->Draw("PZ");
  lg->Draw();
  l2->Draw();
  tx.DrawLatex(0.45, 0.6 ,"EleLeg (HLT_Mu23_Ele12)"); 


  TCanvas* cpt_mu = new TCanvas("cpt_mu", "cpt_mu", 650, 650);
  cpt_mu->cd();
  gMuPt_Mu23Ele12_2017F->GetYaxis()->SetTitle("trigger efficiency");
  gMuPt_Mu23Ele12_2017F->GetXaxis()->SetTitle("p_{T}^{e} (GeV/c)");
  //gMuPt_Mu23Ele12_2017F->GetXaxis()->SetLimits(0, 200.);
  gMuPt_Mu23Ele12_2017F->GetYaxis()->SetRangeUser(0., 1.2);
  gMuPt_Mu23Ele12_2017F->Draw("APZ");
  gMuPt_Mu23Ele12_2017D->Draw("PZ");
  gMuPt_Mu23Ele12_2017E->Draw("PZ");
  lg->Draw();
  l1->Draw();
  tx.DrawLatex(0.45, 0.6 ,"MuLeg (HLT_Mu23_Ele12)"); 


  TCanvas* ceta_mu = new TCanvas("ceta_mu", "ceta_mu", 650, 650);
  ceta_mu->cd();
  gMuEta_Mu23Ele12_2017F->GetYaxis()->SetTitle("trigger efficiency");
  gMuEta_Mu23Ele12_2017F->GetXaxis()->SetTitle("#eta^{#mu}");
  //gMuEta_Mu23Ele12_2017F->GetXaxis()->SetLimits(0, 200.);
  gMuEta_Mu23Ele12_2017F->GetYaxis()->SetRangeUser(0., 1.2);
  gMuEta_Mu23Ele12_2017F->Draw("APZ");
  gMuEta_Mu23Ele12_2017D->Draw("PZ");
  gMuEta_Mu23Ele12_2017E->Draw("PZ");
  lg->Draw();
  l2->Draw();
  tx.DrawLatex(0.45, 0.6 ,"MuLeg (HLT_Mu23_Ele12)"); 

  TCanvas* cphi_mu = new TCanvas("cphi_mu", "cphi_mu", 650, 650);
  cphi_mu->cd();
  gMuPhi_Mu23Ele12_2017F->GetYaxis()->SetTitle("trigger efficiency");
  gMuPhi_Mu23Ele12_2017F->GetXaxis()->SetTitle("#phi^{#mu}");
  //gMuPhi_Mu23Ele12_2017F->GetXaxis()->SetLimits(0, 200.);
  gMuPhi_Mu23Ele12_2017F->GetYaxis()->SetRangeUser(0., 1.2);
  gMuPhi_Mu23Ele12_2017F->Draw("APZ");
  gMuPhi_Mu23Ele12_2017D->Draw("PZ");
  gMuPhi_Mu23Ele12_2017E->Draw("PZ");
  lg->Draw();
  l2->Draw();
  tx.DrawLatex(0.45, 0.6 ,"MuLeg (HLT_Mu23_Ele12)"); 

  TCanvas* cvtx_mu = new TCanvas("cvtx_mu", "cvtx_mu", 650, 650);
  cvtx_mu->cd();
  gMuVtx_Mu23Ele12_2017F->GetYaxis()->SetTitle("trigger efficiency");
  gMuVtx_Mu23Ele12_2017F->GetXaxis()->SetTitle("nVtx");
  //gMuVtx_Mu23Ele12_2017F->GetXaxis()->SetLimits(0, 200.);
  gMuVtx_Mu23Ele12_2017F->GetYaxis()->SetRangeUser(0., 1.2);
  gMuVtx_Mu23Ele12_2017F->Draw("APZ");
  gMuVtx_Mu23Ele12_2017D->Draw("PZ");
  gMuVtx_Mu23Ele12_2017E->Draw("PZ");
  lg->Draw();
  l2->Draw();
  tx.DrawLatex(0.45, 0.6 ,"MuLeg (HLT_Mu23_Ele12)"); 


  cout << "plotting for Mu12_Ele23 path" << endl;
  TCanvas* cpt_ele_mu12ele23 = new TCanvas("cpt_ele_mu12ele23", "cpt_ele_mu12ele23", 650, 650);
  cpt_ele_mu12ele23->cd();
  gElePt_Mu12Ele23_2017F->GetYaxis()->SetTitle("trigger efficiency");
  gElePt_Mu12Ele23_2017F->GetXaxis()->SetTitle("p_{T}^{e} (GeV/c)");
  //gElePt_Mu12Ele23_2017F->GetXaxis()->SetLimits(0, 200.);
  gElePt_Mu12Ele23_2017F->GetYaxis()->SetRangeUser(0., 1.2);
  gElePt_Mu12Ele23_2017F->Draw("APZ");
  gElePt_Mu12Ele23_2017D->Draw("PZ");
  gElePt_Mu12Ele23_2017E->Draw("PZ");
  TLegend *lg1 = new TLegend(0.35, 0.26, 0.75, 0.40);
  lg1->AddEntry(gElePt_Mu12Ele23_2017D, "Run2017D_PromptReco_v1", "pl");
  lg1->AddEntry(gElePt_Mu12Ele23_2017E, "Run2017E_PromptReco_v1", "pl");
  lg1->AddEntry(gElePt_Mu12Ele23_2017F, "Run2017F_PromptReco_v1", "pl");
  //lg->AddEntry(gElePt_EB_2017B_v2_leg1, "Run2017B_PromptReco_v2", "pl");
  //lg->AddEntry(gElePt_EB_2017C_v1_leg1, "Run2017C_PromptReco_v1", "pl");
  //lg->AddEntry(gElePt_EB_2017C_v2_leg1, "Run2017C_PromptReco_v2", "pl");
  lg1->SetLineColor(0);
  lg1->SetFillColor(0);
  lg1->SetShadowColor(0);
  lg1->SetTextSize(0.035);
  lg1->SetTextFont(42);
  lg1->Draw();
  l1->Draw();
  tx.DrawLatex(0.45, 0.6 ,"EleLeg (HLT_Mu12_Ele23)"); 

  TCanvas* ceta_ele_mu12ele23 = new TCanvas("ceta_ele_mu12ele23", "ceta_ele_mu12ele23", 650, 650);
  ceta_ele_mu12ele23->cd();
  gEleEta_Mu12Ele23_2017F->GetYaxis()->SetTitle("trigger efficiency");
  gEleEta_Mu12Ele23_2017F->GetXaxis()->SetTitle("#eta^{e}");
  //gEleEta_Mu12Ele23_2017F->GetXaxis()->SetLimits(0, 200.);
  gEleEta_Mu12Ele23_2017F->GetYaxis()->SetRangeUser(0., 1.2);
  gEleEta_Mu12Ele23_2017F->Draw("APZ");
  gEleEta_Mu12Ele23_2017D->Draw("PZ");
  gEleEta_Mu12Ele23_2017E->Draw("PZ");
  lg1->Draw();
  l2->Draw();
  tx.DrawLatex(0.45, 0.6 ,"EleLeg (HLT_Mu12_Ele23)"); 

  TCanvas* cphi_ele_mu12ele23 = new TCanvas("cphi_ele_mu12ele23", "cphi_ele_mu12ele23", 650, 650);
  cphi_ele_mu12ele23->cd();
  gElePhi_Mu12Ele23_2017F->GetYaxis()->SetTitle("trigger efficiency");
  gElePhi_Mu12Ele23_2017F->GetXaxis()->SetTitle("#phi^{e}");
  gElePhi_Mu12Ele23_2017F->GetYaxis()->SetRangeUser(0., 1.2);
  gElePhi_Mu12Ele23_2017F->Draw("APZ");
  gElePhi_Mu12Ele23_2017D->Draw("PZ");
  gElePhi_Mu12Ele23_2017E->Draw("PZ");
  lg1->Draw();
  l2->Draw();
  tx.DrawLatex(0.45, 0.6 ,"EleLeg (HLT_Mu12_Ele23)"); 

  TCanvas* cvtx_ele_mu12ele23 = new TCanvas("cvtx_ele_mu12ele23", "cvtx_ele_mu12ele23", 650, 650);
  cvtx_ele_mu12ele23->cd();
  gEleVtx_Mu12Ele23_2017F->GetYaxis()->SetTitle("trigger efficiency");
  gEleVtx_Mu12Ele23_2017F->GetXaxis()->SetTitle("nVtx");
  gEleVtx_Mu12Ele23_2017F->GetYaxis()->SetRangeUser(0., 1.2);
  gEleVtx_Mu12Ele23_2017F->Draw("APZ");
  gEleVtx_Mu12Ele23_2017D->Draw("PZ");
  gEleVtx_Mu12Ele23_2017E->Draw("PZ");
  lg1->Draw();
  l2->Draw();
  tx.DrawLatex(0.45, 0.6 ,"EleLeg (HLT_Mu12_Ele23)"); 


  TCanvas* cpt_mu_mu12ele23 = new TCanvas("cpt_mu_mu12ele23", "cpt_mu_mu12ele23", 650, 650);
  cpt_mu_mu12ele23->cd();
  gMuPt_Mu12Ele23_2017F->GetYaxis()->SetTitle("trigger efficiency");
  gMuPt_Mu12Ele23_2017F->GetXaxis()->SetTitle("p_{T}^{e} (GeV/c)");
  //gMuPt_Mu12Ele23_2017F->GetXaxis()->SetLimits(0, 200.);
  gMuPt_Mu12Ele23_2017F->GetYaxis()->SetRangeUser(0., 1.2);
  gMuPt_Mu12Ele23_2017F->Draw("APZ");
  gMuPt_Mu12Ele23_2017D->Draw("PZ");
  gMuPt_Mu12Ele23_2017E->Draw("PZ");
  lg1->Draw();
  l1->Draw();
  tx.DrawLatex(0.45, 0.6 ,"MuLeg (HLT_Mu12_Ele23)"); 

  TCanvas* ceta_mu_mu12ele23 = new TCanvas("ceta_mu_mu12ele23", "ceta_mu_mu12ele23", 650, 650);
  ceta_mu_mu12ele23->cd();
  gMuEta_Mu12Ele23_2017F->GetYaxis()->SetTitle("trigger efficiency");
  gMuEta_Mu12Ele23_2017F->GetXaxis()->SetTitle("#eta^{#mu}");
  //gMuEta_Mu12Ele23_2017F->GetXaxis()->SetLimits(0, 200.);
  gMuEta_Mu12Ele23_2017F->GetYaxis()->SetRangeUser(0., 1.2);
  gMuEta_Mu12Ele23_2017F->Draw("APZ");
  gMuEta_Mu12Ele23_2017D->Draw("PZ");
  gMuEta_Mu12Ele23_2017E->Draw("PZ");
  lg1->Draw();
  l2->Draw();
  tx.DrawLatex(0.45, 0.6 ,"MuLeg (HLT_Mu12_Ele23)"); 

  TCanvas* cphi_mu_mu12ele23 = new TCanvas("cphi_mu_mu12ele23", "cphi_mu_mu12ele23", 650, 650);
  cphi_mu_mu12ele23->cd();
  gMuPhi_Mu12Ele23_2017F->GetYaxis()->SetTitle("trigger efficiency");
  gMuPhi_Mu12Ele23_2017F->GetXaxis()->SetTitle("#phi^{#mu}");
  gMuPhi_Mu12Ele23_2017F->GetYaxis()->SetRangeUser(0., 1.2);
  gMuPhi_Mu12Ele23_2017F->Draw("APZ");
  gMuPhi_Mu12Ele23_2017D->Draw("PZ");
  gMuPhi_Mu12Ele23_2017E->Draw("PZ");
  lg1->Draw();
  l2->Draw();
  tx.DrawLatex(0.45, 0.6 ,"MuLeg (HLT_Mu12_Ele23)"); 


  TCanvas* cvtx_mu_mu12ele23 = new TCanvas("cvtx_mu_mu12ele23", "cvtx_mu_mu12ele23", 650, 650);
  cvtx_mu_mu12ele23->cd();
  gMuVtx_Mu12Ele23_2017F->GetYaxis()->SetTitle("trigger efficiency");
  gMuVtx_Mu12Ele23_2017F->GetXaxis()->SetTitle("nVtx");
  gMuVtx_Mu12Ele23_2017F->GetYaxis()->SetRangeUser(0., 1.2);
  gMuVtx_Mu12Ele23_2017F->Draw("APZ");
  gMuVtx_Mu12Ele23_2017D->Draw("PZ");
  gMuVtx_Mu12Ele23_2017E->Draw("PZ");
  lg1->Draw();
  l2->Draw();
  tx.DrawLatex(0.45, 0.6 ,"MuLeg (HLT_Mu12_Ele23)"); 


  //saving plots
  if (leg.Contains("ele")) {
    cpt_ele->Print("plots/dilepton/Eff_EleLeg_vsPt_HLT_Mu23Ele12_RunDtoF.pdf");
    ceta_ele->Print("plots/dilepton/Eff_EleLeg_vsEta_HLT_Mu23Ele12_RunDtoF.pdf");
    cphi_ele->Print("plots/dilepton/Eff_EleLeg_vsPhi_HLT_Mu23Ele12_RunDtoF.pdf");
    cvtx_ele->Print("plots/dilepton/Eff_EleLeg_vsVtx_HLT_Mu23Ele12_RunDtoF.pdf");
    cpt_ele_mu12ele23->Print("plots/dilepton/Eff_EleLeg_vsPt_HLT_Mu12Ele23_RunDtoF.pdf");
    ceta_ele_mu12ele23->Print("plots/dilepton/Eff_EleLeg_vsEta_HLT_Mu12Ele23_RunDtoF.pdf");
    cphi_ele_mu12ele23->Print("plots/dilepton/Eff_EleLeg_vsPhi_HLT_Mu12Ele23_RunDtoF.pdf");
    cvtx_ele_mu12ele23->Print("plots/dilepton/Eff_EleLeg_vsVtx_HLT_Mu12Ele23_RunDtoF.pdf");
  }
  else {
    cpt_mu->Print("plots/dilepton/Eff_muLeg_vsPt_HLT_Mu23Ele12_RunDtoF.pdf");
    ceta_mu->Print("plots/dilepton/Eff_muLeg_vsEta_HLT_Mu23Ele12_RunDtoF.pdf");
    cphi_mu->Print("plots/dilepton/Eff_muLeg_vsPhi_HLT_Mu23Ele12_RunDtoF.pdf");
    cvtx_mu->Print("plots/dilepton/Eff_muLeg_vsVtx_HLT_Mu23Ele12_RunDtoF.pdf");
    cpt_mu_mu12ele23->Print("plots/dilepton/Eff_muLeg_vsPt_HLT_Mu12Ele23_RunDtoF.pdf");
    ceta_mu_mu12ele23->Print("plots/dilepton/Eff_muLeg_vsEta_HLT_Mu12Ele23_RunDtoF.pdf");
    cphi_mu_mu12ele23->Print("plots/dilepton/Eff_muLeg_vsPhi_HLT_Mu12Ele23_RunDtoF.pdf");
    cvtx_mu_mu12ele23->Print("plots/dilepton/Eff_muLeg_vsVtx_HLT_Mu12Ele23_RunDtoF.pdf");
  }

  */

  //for comparing with HLT_Ele23_Ele12

  TCanvas* c1 = new TCanvas("c1", "c1", 650, 650);
  c1->cd();
  gElePt_Mu12Ele23_2017D->GetYaxis()->SetTitle("trigger efficiency");
  gElePt_Mu12Ele23_2017D->GetXaxis()->SetTitle("p_{T}^{e} (GeV/c)");
  gElePt_Mu12Ele23_2017D->GetXaxis()->SetLimits(0, 200.);
  gElePt_Mu12Ele23_2017D->GetYaxis()->SetRangeUser(0., 1.2);
  gElePt_Mu12Ele23_2017D->Draw("APZ");
  gEle_EE_leg1_2017D->Draw("PZ");
  TLegend *lg = new TLegend(0.35, 0.26, 0.75, 0.40);
  lg->AddEntry(gElePt_Mu12Ele23_2017D, "HLT_Ele23_Mu12", "pl");
  lg->AddEntry(gEle_EE_leg1_2017D, "HLT_Ele23_Ele12", "pl");
  lg->SetLineColor(0);
  lg->SetFillColor(0);
  lg->SetShadowColor(0);
  lg->SetTextSize(0.035);
  lg->SetTextFont(42);
  lg->Draw();
  l1->Draw();
  tx.DrawLatex(0.45, 0.66 ,"Run2017D");
  tx.DrawLatex(0.45, 0.6 ,"EleLeg_Pt23"); 

  TCanvas* c2 = new TCanvas("c2", "c2", 650, 650);
  c2->cd();
  gElePt_Mu12Ele23_2017E->GetYaxis()->SetTitle("trigger efficiency");
  gElePt_Mu12Ele23_2017E->GetXaxis()->SetTitle("p_{T}^{e} (GeV/c)");
  gElePt_Mu12Ele23_2017E->GetXaxis()->SetLimits(0, 200.);
  gElePt_Mu12Ele23_2017E->GetYaxis()->SetRangeUser(0., 1.2);
  gElePt_Mu12Ele23_2017E->Draw("APZ");
  gEle_EE_leg1_2017E->Draw("PZ");
  lg->Draw();
  l1->Draw();
  tx.DrawLatex(0.45, 0.66 ,"Run2017E");
  tx.DrawLatex(0.45, 0.6 ,"EleLeg_Pt23"); 

  TCanvas* c3 = new TCanvas("c3", "c3", 650, 650);
  c3->cd();
  gElePt_Mu12Ele23_2017F->GetYaxis()->SetTitle("trigger efficiency");
  gElePt_Mu12Ele23_2017F->GetXaxis()->SetTitle("p_{T}^{e} (GeV/c)");
  gElePt_Mu12Ele23_2017F->GetXaxis()->SetLimits(0, 200.);
  gElePt_Mu12Ele23_2017F->GetYaxis()->SetRangeUser(0., 1.2);
  gElePt_Mu12Ele23_2017F->Draw("APZ");
  gEle_EE_leg1_2017F->Draw("PZ");
  lg->Draw();
  l1->Draw();
  tx.DrawLatex(0.45, 0.66 ,"Run2017F");
  tx.DrawLatex(0.45, 0.6 ,"EleLeg_Pt23"); 


  TCanvas* c4 = new TCanvas("c4", "c4", 650, 650);
  c4->cd();
  gElePt_Mu23Ele12_2017D->GetYaxis()->SetTitle("trigger efficiency");
  gElePt_Mu23Ele12_2017D->GetXaxis()->SetTitle("p_{T}^{e} (GeV/c)");
  gElePt_Mu23Ele12_2017D->GetXaxis()->SetLimits(0, 200.);
  gElePt_Mu23Ele12_2017D->GetYaxis()->SetRangeUser(0., 1.2);
  gElePt_Mu23Ele12_2017D->Draw("APZ");
  gEle_EE_leg2_2017D->Draw("PZ");
  lg->Draw();
  l1->Draw();
  tx.DrawLatex(0.45, 0.66 ,"Run2017D");
  tx.DrawLatex(0.45, 0.6 ,"EleLeg_Pt12"); 

  TCanvas* c5 = new TCanvas("c5", "c5", 650, 650);
  c5->cd();
  gElePt_Mu23Ele12_2017E->GetYaxis()->SetTitle("trigger efficiency");
  gElePt_Mu23Ele12_2017E->GetXaxis()->SetTitle("p_{T}^{e} (GeV/c)");
  gElePt_Mu23Ele12_2017E->GetXaxis()->SetLimits(0, 200.);
  gElePt_Mu23Ele12_2017E->GetYaxis()->SetRangeUser(0., 1.2);
  gElePt_Mu23Ele12_2017E->Draw("APZ");
  gEle_EE_leg2_2017E->Draw("PZ");
  lg->Draw();
  l1->Draw();
  tx.DrawLatex(0.45, 0.66 ,"Run2017E");
  tx.DrawLatex(0.45, 0.6 ,"EleLeg_Pt12"); 

  TCanvas* c6 = new TCanvas("c6", "c6", 650, 650);
  c6->cd();
  gElePt_Mu23Ele12_2017F->GetYaxis()->SetTitle("trigger efficiency");
  gElePt_Mu23Ele12_2017F->GetXaxis()->SetTitle("p_{T}^{e} (GeV/c)");
  gElePt_Mu23Ele12_2017F->GetXaxis()->SetLimits(0, 200.);
  gElePt_Mu23Ele12_2017F->GetYaxis()->SetRangeUser(0., 1.2);
  gElePt_Mu23Ele12_2017F->Draw("APZ");
  gEle_EE_leg2_2017F->Draw("PZ");
  lg->Draw();
  l1->Draw();
  tx.DrawLatex(0.45, 0.66 ,"Run2017F");
  tx.DrawLatex(0.45, 0.6 ,"EleLeg_Pt12"); 


  c1->SaveAs("plots/dilepton/Eff_Ele23_Endcap_DiEle23Ele12_vs_Ele23Mu12_RunD.pdf");
  c2->SaveAs("plots/dilepton/Eff_Ele23_Endcap_DiEle23Ele12_vs_Ele23Mu12_RunE.pdf");
  c3->SaveAs("plots/dilepton/Eff_Ele23_Endcap_DiEle23Ele12_vs_Ele23Mu12_RunF.pdf");


}

