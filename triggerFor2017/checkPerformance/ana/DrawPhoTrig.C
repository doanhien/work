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

void DrawPhoTrig(TString mode = "tripho35") {

  gStyle->SetOptStat(1111);
  cout << "starting program" << endl;
  TString charname = "";
  if (mode.Contains("tripho20") ) charname = "output/TriPho20";
  else if (mode.Contains("tripho30") )charname = "output/TriPho303010";
  else if (mode.Contains("tripho35") )charname = "output/TriPho35355";

  cout << "read input file" << endl;

  TFile *frunD = new TFile(charname + "_2017D_PromptReco_v1_v3.1.root", "READ");
  TFile *frunE = new TFile(charname + "_2017E_PromptReco_v1_v3.1.root", "READ");
  TFile *frunF = new TFile(charname + "_2017F_PromptReco_v1_v3.1.root", "READ");

  cout << "define tgrap" << endl;
  //pt in barrel
  TGraphAsymmErrors *gElePt_EB_2017D_leg1;
  TGraphAsymmErrors *gElePt_EB_2017E_leg1;
  TGraphAsymmErrors *gElePt_EB_2017F_leg1;
  //TGraphAsymmErrors *gElePt_EB_2017B_v2_leg1;
  //TGraphAsymmErrors *gElePt_EB_2017C_v1_leg1;
  //TGraphAsymmErrors *gElePt_EB_2017C_v2_leg1;

  TGraphAsymmErrors *gElePt_EB_2017D_leg2;
  TGraphAsymmErrors *gElePt_EB_2017E_leg2;
  TGraphAsymmErrors *gElePt_EB_2017F_leg2;
  //TGraphAsymmErrors *gElePt_EB_2017B_v2_leg2;
  //TGraphAsymmErrors *gElePt_EB_2017C_v1_leg2;
  //TGraphAsymmErrors *gElePt_EB_2017C_v2_leg2;

  TGraphAsymmErrors *gElePt_EB_2017D_leg3;
  TGraphAsymmErrors *gElePt_EB_2017E_leg3;
  TGraphAsymmErrors *gElePt_EB_2017F_leg3;
  //TGraphAsymmErrors *gElePt_EB_2017B_v2_leg3;
  //TGraphAsymmErrors *gElePt_EB_2017C_v1_leg3;
  //TGraphAsymmErrors *gElePt_EB_2017C_v2_leg3;

  //pt in endcaps
  TGraphAsymmErrors *gElePt_EE_2017D_leg1;
  TGraphAsymmErrors *gElePt_EE_2017E_leg1;
  TGraphAsymmErrors *gElePt_EE_2017F_leg1;
  //TGraphAsymmErrors *gElePt_EE_2017B_v2_leg1;
  //TGraphAsymmErrors *gElePt_EE_2017C_v1_leg1;
  //TGraphAsymmErrors *gElePt_EE_2017C_v2_leg1;

  TGraphAsymmErrors *gElePt_EE_2017D_leg2;
  TGraphAsymmErrors *gElePt_EE_2017E_leg2;
  TGraphAsymmErrors *gElePt_EE_2017F_leg2;
  //TGraphAsymmErrors *gElePt_EE_2017B_v2_leg2;
  //TGraphAsymmErrors *gElePt_EE_2017C_v1_leg2;
  //TGraphAsymmErrors *gElePt_EE_2017C_v2_leg2;

  TGraphAsymmErrors *gElePt_EE_2017D_leg3;
  TGraphAsymmErrors *gElePt_EE_2017E_leg3;
  TGraphAsymmErrors *gElePt_EE_2017F_leg3;
  //TGraphAsymmErrors *gElePt_EE_2017B_v2_leg3;
  //TGraphAsymmErrors *gElePt_EE_2017C_v1_leg3;
  //TGraphAsymmErrors *gElePt_EE_2017C_v2_leg3;

  //eta
  /*
  TGraphAsymmErrors *gEleEta_2017D_leg1;
  TGraphAsymmErrors *gEleEta_2017E_leg1;
  TGraphAsymmErrors *gEleEta_2017F_leg1;
  TGraphAsymmErrors *gEleEta_2017B_v2_leg1;
  TGraphAsymmErrors *gEleEta_2017C_v1_leg1;
  TGraphAsymmErrors *gEleEta_2017C_v2_leg1;

  TGraphAsymmErrors *gEleEta_2017D_leg2;
  TGraphAsymmErrors *gEleEta_2017E_leg2;
  TGraphAsymmErrors *gEleEta_2017F_leg2;
  TGraphAsymmErrors *gEleEta_2017B_v2_leg2;
  TGraphAsymmErrors *gEleEta_2017C_v1_leg2;
  TGraphAsymmErrors *gEleEta_2017C_v2_leg2;

  TGraphAsymmErrors *gEleEta_2017D_leg3;
  TGraphAsymmErrors *gEleEta_2017E_leg3;
  TGraphAsymmErrors *gEleEta_2017F_leg3;
  TGraphAsymmErrors *gEleEta_2017B_v2_leg3;
  TGraphAsymmErrors *gEleEta_2017C_v1_leg3;
  TGraphAsymmErrors *gEleEta_2017C_v2_leg3;

  //phi
  TGraphAsymmErrors *gElePhi_2017D_leg1;
  TGraphAsymmErrors *gElePhi_2017E_leg1;
  TGraphAsymmErrors *gElePhi_2017F_leg1;
  TGraphAsymmErrors *gElePhi_2017B_v2_leg1;
  TGraphAsymmErrors *gElePhi_2017C_v1_leg1;
  TGraphAsymmErrors *gElePhi_2017C_v2_leg1;

  TGraphAsymmErrors *gElePhi_2017D_leg2;
  TGraphAsymmErrors *gElePhi_2017E_leg2;
  TGraphAsymmErrors *gElePhi_2017F_leg2;
  TGraphAsymmErrors *gElePhi_2017B_v2_leg2;
  TGraphAsymmErrors *gElePhi_2017C_v1_leg2;
  TGraphAsymmErrors *gElePhi_2017C_v2_leg2;

  TGraphAsymmErrors *gElePhi_2017D_leg3;
  TGraphAsymmErrors *gElePhi_2017E_leg3;
  TGraphAsymmErrors *gElePhi_2017F_leg3;
  TGraphAsymmErrors *gElePhi_2017B_v2_leg3;
  TGraphAsymmErrors *gElePhi_2017C_v1_leg3;
  TGraphAsymmErrors *gElePhi_2017C_v2_leg3;
  */
  cout << "read graph from input file" << endl;


  //barrel
  gElePt_EB_2017D_leg1 = (TGraphAsymmErrors*) frunD->Get("gPhoPt_EB_leg1");
  gElePt_EB_2017E_leg1 = (TGraphAsymmErrors*) frunE->Get("gPhoPt_EB_leg1");
  gElePt_EB_2017F_leg1 = (TGraphAsymmErrors*) frunF->Get("gPhoPt_EB_leg1");

  gElePt_EB_2017D_leg2 = (TGraphAsymmErrors*) frunD->Get("gPhoPt_EB_leg2");
  gElePt_EB_2017E_leg2 = (TGraphAsymmErrors*) frunE->Get("gPhoPt_EB_leg2");
  gElePt_EB_2017F_leg2 = (TGraphAsymmErrors*) frunF->Get("gPhoPt_EB_leg2");

  //for TriPho20
  //gElePt_EB_2017D_leg2 = (TGraphAsymmErrors*) frunD->Get("gPhoPt_EB");
  //gElePt_EB_2017E_leg2 = (TGraphAsymmErrors*) frunE->Get("gPhoPt_EB");
  //gElePt_EB_2017F_leg2 = (TGraphAsymmErrors*) frunF->Get("gPhoPt_EB");

  gElePt_EB_2017D_leg3 = (TGraphAsymmErrors*) frunD->Get("gPhoPt_EB_leg3");
  gElePt_EB_2017E_leg3 = (TGraphAsymmErrors*) frunE->Get("gPhoPt_EB_leg3");
  gElePt_EB_2017F_leg3 = (TGraphAsymmErrors*) frunF->Get("gPhoPt_EB_leg3");

  //endcaps
  gElePt_EE_2017D_leg1 = (TGraphAsymmErrors*) frunD->Get("gPhoPt_EE_leg1");
  gElePt_EE_2017E_leg1 = (TGraphAsymmErrors*) frunE->Get("gPhoPt_EE_leg1");
  gElePt_EE_2017F_leg1 = (TGraphAsymmErrors*) frunF->Get("gPhoPt_EE_leg1");

  //gElePt_EE_2017D_leg2 = (TGraphAsymmErrors*) frunD->Get("gPhoPt_EE");
  //gElePt_EE_2017E_leg2 = (TGraphAsymmErrors*) frunE->Get("gPhoPt_EE");
  //gElePt_EE_2017F_leg2 = (TGraphAsymmErrors*) frunF->Get("gPhoPt_EE");

  gElePt_EE_2017D_leg2 = (TGraphAsymmErrors*) frunD->Get("gPhoPt_EE_leg2");
  gElePt_EE_2017E_leg2 = (TGraphAsymmErrors*) frunE->Get("gPhoPt_EE_leg2");
  gElePt_EE_2017F_leg2 = (TGraphAsymmErrors*) frunF->Get("gPhoPt_EE_leg2");

  gElePt_EE_2017D_leg3 = (TGraphAsymmErrors*) frunD->Get("gPhoPt_EE_leg3");
  gElePt_EE_2017E_leg3 = (TGraphAsymmErrors*) frunE->Get("gPhoPt_EE_leg3");
  gElePt_EE_2017F_leg3 = (TGraphAsymmErrors*) frunF->Get("gPhoPt_EE_leg3");


  //R9
  gElePtR9_EB_2017D_leg1 = (TGraphAsymmErrors*) frunD->Get("gPhoPtR9_EB_leg1");
  gElePtR9_EB_2017E_leg1 = (TGraphAsymmErrors*) frunE->Get("gPhoPtR9_EB_leg1");
  gElePtR9_EB_2017F_leg1 = (TGraphAsymmErrors*) frunF->Get("gPhoPtR9_EB_leg1");

  //gElePtR9_EB_2017D_leg2 = (TGraphAsymmErrors*) frunD->Get("gPhoPtR9_EB");
  //gElePtR9_EB_2017E_leg2 = (TGraphAsymmErrors*) frunE->Get("gPhoPtR9_EB");
  //gElePtR9_EB_2017F_leg2 = (TGraphAsymmErrors*) frunF->Get("gPhoPtR9_EB");

  gElePtR9_EB_2017D_leg2 = (TGraphAsymmErrors*) frunD->Get("gPhoPtR9_EB_leg2");
  gElePtR9_EB_2017E_leg2 = (TGraphAsymmErrors*) frunE->Get("gPhoPtR9_EB_leg2");
  gElePtR9_EB_2017F_leg2 = (TGraphAsymmErrors*) frunF->Get("gPhoPtR9_EB_leg2");
  gElePtR9_EB_2017D_leg3 = (TGraphAsymmErrors*) frunD->Get("gPhoPtR9_EB_leg3");
  gElePtR9_EB_2017E_leg3 = (TGraphAsymmErrors*) frunE->Get("gPhoPtR9_EB_leg3");
  gElePtR9_EB_2017F_leg3 = (TGraphAsymmErrors*) frunF->Get("gPhoPtR9_EB_leg3");

  gElePtR9_EE_2017D_leg1 = (TGraphAsymmErrors*) frunD->Get("gPhoPtR9_EE_leg1");
  gElePtR9_EE_2017E_leg1 = (TGraphAsymmErrors*) frunE->Get("gPhoPtR9_EE_leg1");
  gElePtR9_EE_2017F_leg1 = (TGraphAsymmErrors*) frunF->Get("gPhoPtR9_EE_leg1");

  //gElePtR9_EE_2017D_leg2 = (TGraphAsymmErrors*) frunD->Get("gPhoPtR9_EE");
  //gElePtR9_EE_2017E_leg2 = (TGraphAsymmErrors*) frunE->Get("gPhoPtR9_EE");
  //gElePtR9_EE_2017F_leg2 = (TGraphAsymmErrors*) frunF->Get("gPhoPtR9_EE");

  gElePtR9_EE_2017D_leg2 = (TGraphAsymmErrors*) frunD->Get("gPhoPtR9_EE_leg2");
  gElePtR9_EE_2017E_leg2 = (TGraphAsymmErrors*) frunE->Get("gPhoPtR9_EE_leg2");
  gElePtR9_EE_2017F_leg2 = (TGraphAsymmErrors*) frunF->Get("gPhoPtR9_EE_leg2");
  gElePtR9_EE_2017D_leg3 = (TGraphAsymmErrors*) frunD->Get("gPhoPtR9_EE_leg3");
  gElePtR9_EE_2017E_leg3 = (TGraphAsymmErrors*) frunE->Get("gPhoPtR9_EE_leg3");
  gElePtR9_EE_2017F_leg3 = (TGraphAsymmErrors*) frunF->Get("gPhoPtR9_EE_leg3");


  //eta
  /*gEleEta_2017D_leg1 = (TGraphAsymmErrors*) frunD->Get("gEleEta_leg1");
  gEleEta_2017E_leg1 = (TGraphAsymmErrors*) frunE->Get("gEleEta_leg1");
  gEleEta_2017F_leg1 = (TGraphAsymmErrors*) frunF->Get("gEleEta_leg1");
  //gEleEta_2017B_v2_leg1 = (TGraphAsymmErrors*) fB2->Get("gEleEta_leg1");
  //gEleEta_2017C_v1_leg1 = (TGraphAsymmErrors*) fC1->Get("gEleEta_leg1");
  //gEleEta_2017C_v2_leg1 = (TGraphAsymmErrors*) fC2->Get("gEleEta_leg1");

  gEleEta_2017D_leg2 = (TGraphAsymmErrors*) frunD->Get("gEleEta_leg2");
  gEleEta_2017E_leg2 = (TGraphAsymmErrors*) frunE->Get("gEleEta_leg2");
  gEleEta_2017F_leg2 = (TGraphAsymmErrors*) frunF->Get("gEleEta_leg2");
  //gEleEta_2017B_v2_leg2 = (TGraphAsymmErrors*) fB2->Get("gEleEta_leg2");
  //gEleEta_2017C_v1_leg2 = (TGraphAsymmErrors*) fC1->Get("gEleEta_leg2");
  //gEleEta_2017C_v2_leg2 = (TGraphAsymmErrors*) fC2->Get("gEleEta_leg2");

  //phi
  gElePhi_2017D_leg1 = (TGraphAsymmErrors*) frunD->Get("gElePhi_leg1");
  gElePhi_2017E_leg1 = (TGraphAsymmErrors*) frunE->Get("gElePhi_leg1");
  gElePhi_2017F_leg1 = (TGraphAsymmErrors*) frunF->Get("gElePhi_leg1");
  //gElePhi_2017B_v2_leg1 = (TGraphAsymmErrors*) fB2->Get("gElePhi_leg1");
  //gElePhi_2017C_v1_leg1 = (TGraphAsymmErrors*) fC1->Get("gElePhi_leg1");
  //gElePhi_2017C_v2_leg1 = (TGraphAsymmErrors*) fC2->Get("gElePhi_leg1");

  gElePhi_2017D_leg2 = (TGraphAsymmErrors*) frunD->Get("gElePhi_leg2");
  gElePhi_2017E_leg2 = (TGraphAsymmErrors*) frunE->Get("gElePhi_leg2");
  gElePhi_2017F_leg2 = (TGraphAsymmErrors*) frunF->Get("gElePhi_leg2");
  //gElePhi_2017B_v2_leg2 = (TGraphAsymmErrors*) fB2->Get("gElePhi_leg2");
  //gElePhi_2017C_v1_leg2 = (TGraphAsymmErrors*) fC1->Get("gElePhi_leg2");
  //gElePhi_2017C_v2_leg2 = (TGraphAsymmErrors*) fC2->Get("gElePhi_leg2");
  */

  cout << "graph style" << endl;
  //style for graph


  GraphStyle1(gElePt_EB_2017D_leg1);
  GraphStyle2(gElePt_EB_2017E_leg1);
  GraphStyle3(gElePt_EB_2017F_leg1);
  
  GraphStyle1(gElePt_EB_2017D_leg2);
  GraphStyle2(gElePt_EB_2017E_leg2);
  GraphStyle3(gElePt_EB_2017F_leg2);

  GraphStyle1(gElePt_EB_2017D_leg3);
  GraphStyle2(gElePt_EB_2017E_leg3);
  GraphStyle3(gElePt_EB_2017F_leg3);
  
  //for endcaps
  GraphStyle1(gElePt_EE_2017D_leg1);
  GraphStyle2(gElePt_EE_2017E_leg1);
  GraphStyle3(gElePt_EE_2017F_leg1);

  GraphStyle1(gElePt_EE_2017D_leg2);
  GraphStyle2(gElePt_EE_2017E_leg2);
  GraphStyle3(gElePt_EE_2017F_leg2);
  
  GraphStyle1(gElePt_EE_2017D_leg3);
  GraphStyle2(gElePt_EE_2017E_leg3);
  GraphStyle3(gElePt_EE_2017F_leg3);


  //-----------with R9
  GraphStyle1(gElePtR9_EB_2017D_leg1);
  GraphStyle2(gElePtR9_EB_2017E_leg1);
  GraphStyle3(gElePtR9_EB_2017F_leg1);
  
  GraphStyle1(gElePtR9_EB_2017D_leg2);
  GraphStyle2(gElePtR9_EB_2017E_leg2);
  GraphStyle3(gElePtR9_EB_2017F_leg2);

  GraphStyle1(gElePtR9_EB_2017D_leg3);
  GraphStyle2(gElePtR9_EB_2017E_leg3);
  GraphStyle3(gElePtR9_EB_2017F_leg3);

  //for endcaps
  GraphStyle1(gElePtR9_EE_2017D_leg1);
  GraphStyle2(gElePtR9_EE_2017E_leg1);
  GraphStyle3(gElePtR9_EE_2017F_leg1);

  GraphStyle1(gElePtR9_EE_2017D_leg2);
  GraphStyle2(gElePtR9_EE_2017E_leg2);
  GraphStyle3(gElePtR9_EE_2017F_leg2);

  GraphStyle1(gElePtR9_EE_2017D_leg3);
  GraphStyle2(gElePtR9_EE_2017E_leg3);
  GraphStyle3(gElePtR9_EE_2017F_leg3);


  cout << "done graph style" << endl;

  TLatex tx;
  tx.SetNDC(kTRUE);
  tx.SetTextFont(52);
  tx.SetTextSize(0.045);

  TLine *l1 = new TLine(0., 1., 200., 1.);
  l1->SetLineStyle(6);
  l1->SetLineColor(kGreen);
  l1->SetLineWidth(2);

  cout << "starting canvas for plotting" << endl;

  //pt barrel
  TCanvas* cpt_eb_leg1 = new TCanvas("cpt_eb_leg1", "cpt_eb_leg1", 650, 650);
  cpt_eb_leg1->cd();
  gElePt_EB_2017F_leg1->GetYaxis()->SetTitle("trigger efficiency");
  gElePt_EB_2017F_leg1->GetXaxis()->SetTitle("p_{T}^{e} (GeV/c)");
  gElePt_EB_2017F_leg1->GetXaxis()->SetLimits(0, 200.);
  gElePt_EB_2017F_leg1->GetYaxis()->SetRangeUser(0., 1.2);
  gElePt_EB_2017F_leg1->Draw("APZ");
  gElePt_EB_2017D_leg1->Draw("PZ");
  gElePt_EB_2017E_leg1->Draw("PZ");
  TLegend *lg = new TLegend(0.35, 0.26, 0.75, 0.40);
  lg->AddEntry(gElePt_EB_2017D_leg1, "Run2017D_PromptReco_v1", "pl");
  lg->AddEntry(gElePt_EB_2017E_leg1, "Run2017E_PromptReco_v1", "pl");
  lg->AddEntry(gElePt_EB_2017F_leg1, "Run2017F_PromptReco_v1", "pl");
  lg->SetLineColor(0);
  lg->SetFillColor(0);
  lg->SetShadowColor(0);
  lg->SetTextSize(0.035);
  lg->SetTextFont(42);
  lg->Draw();
  l1->Draw();
  //if (mode.Contains("diele") ) tx.DrawLatex(0.45, 0.6 ,"Barrel Leg1_Ele23");
  //else tx.DrawLatex(0.45, 0.6 ,"Barrel Leg1_Ele16");
  //else tx.DrawLatex(0.45, 0.6 ,"Leg1_Ele16");

  TCanvas* cpt_eb_leg2 = new TCanvas("cpt_eb_leg2", "cpt_eb_leg2", 650, 650);
  cpt_eb_leg2->cd();
  gElePt_EB_2017F_leg2->GetYaxis()->SetTitle("trigger efficiency");
  gElePt_EB_2017F_leg2->GetXaxis()->SetTitle("p_{T}^{e} (GeV/c)");
  gElePt_EB_2017F_leg2->GetXaxis()->SetLimits(0, 200.);
  gElePt_EB_2017F_leg2->GetYaxis()->SetRangeUser(0., 1.2);
  gElePt_EB_2017F_leg2->Draw("APZ");
  gElePt_EB_2017D_leg2->Draw("PZ");
  gElePt_EB_2017E_leg2->Draw("PZ");
  lg->Draw();
  l1->Draw();
  //tx.DrawLatex(0.45, 0.6 ,"Leg2_Ele12");


  TCanvas* cpt_eb_leg3 = new TCanvas("cpt_eb_leg3", "cpt_eb_leg3", 650, 650);
  cpt_eb_leg3->cd();
  gElePt_EB_2017F_leg3->GetYaxis()->SetTitle("trigger efficiency");
  gElePt_EB_2017F_leg3->GetXaxis()->SetTitle("p_{T}^{e} (GeV/c)");
  gElePt_EB_2017F_leg3->GetXaxis()->SetLimits(0, 200.);
  gElePt_EB_2017F_leg3->GetYaxis()->SetRangeUser(0., 1.2);
  gElePt_EB_2017F_leg3->Draw("APZ");
  gElePt_EB_2017D_leg3->Draw("PZ");
  gElePt_EB_2017E_leg3->Draw("PZ");
  lg->Draw();
  l1->Draw();


  TCanvas* cpt_ee_leg1 = new TCanvas("cpt_ee_leg1", "cpt_ee_leg1", 650, 650);
  cpt_ee_leg1->cd();
  gElePt_EE_2017D_leg1->GetYaxis()->SetTitle("trigger efficiency");
  gElePt_EE_2017D_leg1->GetXaxis()->SetTitle("p_{T}^{e} (GeV/c)");
  gElePt_EE_2017D_leg1->GetXaxis()->SetLimits(0, 200.);
  gElePt_EE_2017D_leg1->GetYaxis()->SetRangeUser(0., 1.2);
  gElePt_EE_2017D_leg1->Draw("APZ");
  gElePt_EE_2017E_leg1->Draw("PZ");
  gElePt_EE_2017F_leg1->Draw("PZ");
  lg->Draw();
  l1->Draw();
  //if (mode.Contains("diele") ) tx.DrawLatex(0.45, 0.6 ,"Endcaps Leg1_Ele23");
  //else tx.DrawLatex(0.45, 0.6 ,"Endcaps Leg1_Ele16");

  TCanvas* cpt_ee_leg2 = new TCanvas("cpt_ee_leg2", "cpt_ee_leg2", 650, 650);
  cpt_ee_leg2->cd();
  gElePt_EE_2017F_leg2->GetYaxis()->SetTitle("trigger efficiency");
  gElePt_EE_2017F_leg2->GetXaxis()->SetTitle("p_{T}^{e} (GeV/c)");
  gElePt_EE_2017F_leg2->GetXaxis()->SetLimits(0, 200.);
  gElePt_EE_2017F_leg2->GetYaxis()->SetRangeUser(0., 1.2);
  gElePt_EE_2017F_leg2->Draw("APZ");
  gElePt_EE_2017D_leg2->Draw("PZ");
  gElePt_EE_2017E_leg2->Draw("PZ");
  lg->Draw();
  l1->Draw();
  //tx.DrawLatex(0.45, 0.6 ,"Endcaps Leg2_Ele12");


  TCanvas* cpt_ee_leg3 = new TCanvas("cpt_ee_leg3", "cpt_ee_leg3", 650, 650);
  cpt_ee_leg3->cd();
  gElePt_EE_2017F_leg3->GetYaxis()->SetTitle("trigger efficiency");
  gElePt_EE_2017F_leg3->GetXaxis()->SetTitle("p_{T}^{e} (GeV/c)");
  gElePt_EE_2017F_leg3->GetXaxis()->SetLimits(0, 200.);
  gElePt_EE_2017F_leg3->GetYaxis()->SetRangeUser(0., 1.2);
  gElePt_EE_2017F_leg3->Draw("APZ");
  gElePt_EE_2017D_leg3->Draw("PZ");
  gElePt_EE_2017E_leg3->Draw("PZ");
  lg->Draw();
  l1->Draw();


  //---------------with R9
  TCanvas* cptr9_eb_leg1 = new TCanvas("cptr9_eb_leg1", "cptr9_eb_leg1", 650, 650);
  cptr9_eb_leg1->cd();
  gElePtR9_EB_2017F_leg1->GetYaxis()->SetTitle("trigger efficiency");
  gElePtR9_EB_2017F_leg1->GetXaxis()->SetTitle("p_{T}^{e} (GeV/c)");
  gElePtR9_EB_2017F_leg1->GetXaxis()->SetLimits(0, 200.);
  gElePtR9_EB_2017F_leg1->GetYaxis()->SetRangeUser(0., 1.2);
  gElePtR9_EB_2017F_leg1->Draw("APZ");
  gElePtR9_EB_2017D_leg1->Draw("PZ");
  gElePtR9_EB_2017E_leg1->Draw("PZ");
  TLegend *lg1 = new TLegend(0.35, 0.26, 0.75, 0.40);
  lg1->AddEntry(gElePtR9_EB_2017D_leg1, "Run2017D_PromptReco_v1", "pl");
  lg1->AddEntry(gElePtR9_EB_2017E_leg1, "Run2017E_PromptReco_v1", "pl");
  lg1->AddEntry(gElePtR9_EB_2017F_leg1, "Run2017F_PromptReco_v1", "pl");
  lg1->SetLineColor(0);
  lg1->SetFillColor(0);
  lg1->SetShadowColor(0);
  lg1->SetTextSize(0.035);
  lg1->SetTextFont(42);
  lg1->Draw();
  l1->Draw();
  //if (mode.Contains("diele") ) tx.DrawLatex(0.45, 0.6 ,"Barrel Leg1_Ele23");
  //else tx.DrawLatex(0.45, 0.6 ,"Barrel Leg1_Ele16");
  //else tx.DrawLatex(0.45, 0.6 ,"Leg1_Ele16");

  TCanvas* cptr9_eb_leg2 = new TCanvas("cptr9_eb_leg2", "cptr9_eb_leg2", 650, 650);
  cptr9_eb_leg2->cd();
  gElePtR9_EB_2017F_leg2->GetYaxis()->SetTitle("trigger efficiency");
  gElePtR9_EB_2017F_leg2->GetXaxis()->SetTitle("p_{T}^{e} (GeV/c)");
  gElePtR9_EB_2017F_leg2->GetXaxis()->SetLimits(0, 200.);
  gElePtR9_EB_2017F_leg2->GetYaxis()->SetRangeUser(0., 1.2);
  gElePtR9_EB_2017F_leg2->Draw("APZ");
  gElePtR9_EB_2017D_leg2->Draw("PZ");
  gElePtR9_EB_2017E_leg2->Draw("PZ");
  lg1->Draw();
  l1->Draw();
  //tx.DrawLatex(0.45, 0.6 ,"Leg2_Ele12");


  TCanvas* cptr9_eb_leg3 = new TCanvas("cptr9_eb_leg3", "cptr9_eb_leg3", 650, 650);
  cptr9_eb_leg3->cd();
  gElePtR9_EB_2017F_leg3->GetYaxis()->SetTitle("trigger efficiency");
  gElePtR9_EB_2017F_leg3->GetXaxis()->SetTitle("p_{T}^{e} (GeV/c)");
  gElePtR9_EB_2017F_leg3->GetXaxis()->SetLimits(0, 200.);
  gElePtR9_EB_2017F_leg3->GetYaxis()->SetRangeUser(0., 1.2);
  gElePtR9_EB_2017F_leg3->Draw("APZ");
  gElePtR9_EB_2017D_leg3->Draw("PZ");
  gElePtR9_EB_2017E_leg3->Draw("PZ");
  lg1->Draw();
  l1->Draw();


  TCanvas* cptR9_ee_leg1 = new TCanvas("cptR9_ee_leg1", "cptR9_ee_leg1", 650, 650);
  cptR9_ee_leg1->cd();
  gElePtR9_EE_2017D_leg1->GetYaxis()->SetTitle("trigger efficiency");
  gElePtR9_EE_2017D_leg1->GetXaxis()->SetTitle("p_{T}^{e} (GeV/c)");
  gElePtR9_EE_2017D_leg1->GetXaxis()->SetLimits(0, 200.);
  gElePtR9_EE_2017D_leg1->GetYaxis()->SetRangeUser(0., 1.2);
  gElePtR9_EE_2017D_leg1->Draw("APZ");
  gElePtR9_EE_2017E_leg1->Draw("PZ");
  gElePtR9_EE_2017F_leg1->Draw("PZ");
  lg1->Draw();
  l1->Draw();
  //if (mode.Contains("diele") ) tx.DrawLatex(0.45, 0.6 ,"Endcaps Leg1_Ele23");
  //else tx.DrawLatex(0.45, 0.6 ,"Endcaps Leg1_Ele16");

  TCanvas* cptR9_ee_leg2 = new TCanvas("cptR9_ee_leg2", "cptR9_ee_leg2", 650, 650);
  cptR9_ee_leg2->cd();
  gElePtR9_EE_2017F_leg2->GetYaxis()->SetTitle("trigger efficiency");
  gElePtR9_EE_2017F_leg2->GetXaxis()->SetTitle("p_{T}^{e} (GeV/c)");
  gElePtR9_EE_2017F_leg2->GetXaxis()->SetLimits(0, 200.);
  gElePtR9_EE_2017F_leg2->GetYaxis()->SetRangeUser(0., 1.2);
  gElePtR9_EE_2017F_leg2->Draw("APZ");
  gElePtR9_EE_2017D_leg2->Draw("PZ");
  gElePtR9_EE_2017E_leg2->Draw("PZ");
  lg1->Draw();
  l1->Draw();
  //tx.DrawLatex(0.45, 0.6 ,"Endcaps Leg2_Ele12");


  TCanvas* cptR9_ee_leg3 = new TCanvas("cptR9_ee_leg3", "cptR9_ee_leg3", 650, 650);
  cptR9_ee_leg3->cd();
  gElePtR9_EE_2017F_leg3->GetYaxis()->SetTitle("trigger efficiency");
  gElePtR9_EE_2017F_leg3->GetXaxis()->SetTitle("p_{T}^{e} (GeV/c)");
  gElePtR9_EE_2017F_leg3->GetXaxis()->SetLimits(0, 200.);
  gElePtR9_EE_2017F_leg3->GetYaxis()->SetRangeUser(0., 1.2);
  gElePtR9_EE_2017F_leg3->Draw("APZ");
  gElePtR9_EE_2017D_leg3->Draw("PZ");
  gElePtR9_EE_2017E_leg3->Draw("PZ");
  lg1->Draw();
  l1->Draw();



  TString cname;
  if (mode.Contains("tripho20")) cname = "TriPho20";
  else if (mode.Contains("tripho30")) cname = "TriPho303010";
  else cname = "TriPho35355";

  cpt_eb_leg1->SaveAs("plots/" + cname + "_leg1_pt_EB_RunDtoF.pdf");
  cpt_eb_leg2->SaveAs("plots/" + cname + "_leg2_pt_EB_RunDtoF.pdf");
  cpt_eb_leg3->SaveAs("plots/" + cname + "_leg3_pt_EB_RunDtoF.pdf");

  cpt_ee_leg1->SaveAs("plots/" +cname+"_leg1_pt_EE_RunDtoF.pdf");
  cpt_ee_leg2->SaveAs("plots/"+cname+"_leg2_pt_EE_RunDtoF.pdf");
  cpt_ee_leg3->SaveAs("plots/"+cname+"_leg3_pt_EE_RunDtoF.pdf");

  cptr9_eb_leg1->SaveAs("plots/" + cname + "_R9_leg1_pt_EB_RunDtoF.pdf");
  cptr9_eb_leg2->SaveAs("plots/" + cname + "_R9_leg2_pt_EB_RunDtoF.pdf");
  cptr9_eb_leg3->SaveAs("plots/" + cname + "_R9_leg3_pt_EB_RunDtoF.pdf");

  cptR9_ee_leg1->SaveAs("plots/" +cname + "_R9_leg1_pt_EE_RunDtoF.pdf");
  cptR9_ee_leg2->SaveAs("plots/" +cname + "_R9_leg2_pt_EE_RunDtoF.pdf");
  cptR9_ee_leg3->SaveAs("plots/" +cname + "_R9_leg3_pt_EE_RunDtoF.pdf");




}

int main(int argc, char **argv)
{
  if(argc == 2) {
    DrawPhoTrig(argv[1]);
    return 0;
  }
  else {
    return 1;
  }
}
