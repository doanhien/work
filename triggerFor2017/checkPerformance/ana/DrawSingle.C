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

void DrawSingle() {

  gStyle->SetOptStat(1111);
  cout << "starting program" << endl;
  TString charname = "Ele35_WPTight";

  cout << "read input file" << endl;


  TFile *fD = new TFile(charname + "_2017D_PromptReco_v1.root", "READ");
  TFile *fE = new TFile(charname + "_2017E_PromptReco_v1.root", "READ");
  TFile *fF = new TFile(charname + "_2017F_PromptReco_v1.root", "READ");

  cout << "define tgrap" << endl;
  //pt in barrel
  TGraphAsymmErrors *gElePt_EB_2017D;
  TGraphAsymmErrors *gElePt_EB_2017E;
  TGraphAsymmErrors *gElePt_EB_2017F;

  //pt in endcap
  TGraphAsymmErrors *gElePt_EE_2017D;
  TGraphAsymmErrors *gElePt_EE_2017E;
  TGraphAsymmErrors *gElePt_EE_2017F;

  cout << "read graph from input file" << endl;
  //barrel
  gElePt_EB_2017D = (TGraphAsymmErrors*) fD->Get("gEle_EB");
  gElePt_EB_2017E = (TGraphAsymmErrors*) fE->Get("gEle_EB");
  gElePt_EB_2017F = (TGraphAsymmErrors*) fF->Get("gEle_EB");

  //endcaps
  gElePt_EE_2017D = (TGraphAsymmErrors*) fD->Get("gEle_EE");
  gElePt_EE_2017E = (TGraphAsymmErrors*) fE->Get("gEle_EE");
  gElePt_EE_2017F = (TGraphAsymmErrors*) fF->Get("gEle_EE");


  cout << "graph style" << endl;
  //style for graph
  GraphStyle1(gElePt_EB_2017D);
  GraphStyle2(gElePt_EB_2017E);
  GraphStyle3(gElePt_EB_2017F);
  
  //for endcaps
  GraphStyle1(gElePt_EE_2017D);
  GraphStyle2(gElePt_EE_2017E);
  GraphStyle3(gElePt_EE_2017F);

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
  TCanvas* cpt_eb = new TCanvas("cpt_eb", "cpt_eb", 650, 650);
  cpt_eb->cd();
  gElePt_EB_2017F->GetYaxis()->SetTitle("trigger efficiency");
  gElePt_EB_2017F->GetXaxis()->SetTitle("p_{T}^{e} (GeV/c)");
  gElePt_EB_2017F->GetXaxis()->SetLimits(0, 200.);
  gElePt_EB_2017F->GetYaxis()->SetRangeUser(0., 1.2);
  gElePt_EB_2017F->Draw("APZ");
  gElePt_EB_2017D->Draw("PZ");
  gElePt_EB_2017E->Draw("PZ");
  TLegend *lg = new TLegend(0.35, 0.26, 0.75, 0.40);
  lg->AddEntry(gElePt_EB_2017D, "Run2017D_PromptReco_v1", "pl");
  lg->AddEntry(gElePt_EB_2017E, "Run2017E_PromptReco_v1", "pl");
  lg->AddEntry(gElePt_EB_2017F, "Run2017F_PromptReco_v1", "pl");
  lg->SetLineColor(0);
  lg->SetFillColor(0);
  lg->SetShadowColor(0);
  lg->SetTextSize(0.035);
  lg->SetTextFont(42);
  lg->Draw();
  l1->Draw();
  tx.DrawLatex(0.45, 0.6 ,"Barrel");


  TCanvas* cpt_ee = new TCanvas("cpt_ee", "cpt_ee", 650, 650);
  cpt_ee->cd();
  gElePt_EE_2017D->GetYaxis()->SetTitle("trigger efficiency");
  gElePt_EE_2017D->GetXaxis()->SetTitle("p_{T}^{e} (GeV/c)");
  gElePt_EE_2017D->GetXaxis()->SetLimits(0, 200.);
  gElePt_EE_2017D->GetYaxis()->SetRangeUser(0., 1.2);
  gElePt_EE_2017D->Draw("APZ");
  gElePt_EE_2017E->Draw("PZ");
  gElePt_EE_2017F->Draw("PZ");
  lg->Draw();
  l1->Draw();
  tx.DrawLatex(0.45, 0.6 ,"Endcaps");

  cpt_eb->SaveAs("plots/SingleEle35_pt_EB_eff_RunDtoF.pdf");
  cpt_ee->SaveAs("plots/SingleEle35_pt_EE_eff_RunDtoF.pdf");

}
