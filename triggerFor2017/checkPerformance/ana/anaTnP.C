/******************************
script for trigger efficiency of each leg
from di-lepton (e-mu) and tri-lepton (e-mu)
by matching trigger filter, method: cut and cout
******************************/

#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH3D.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include "untuplizer.h"
#include "treeTnP.h"
#include "ElectronSelection.h"
#include "MuonSelection.h"
//#include "cat.h"
#include "puweicalc.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>


void anaTnP(TString inputfile = "HiForest.root", TString outputfile = "XZg.root", bool mc = true)
{

  outputfile = "minitrees/";

  if (inputfile.Contains("SingleMu") ) outputfile += "SingleMu_";
  else if (inputfile.Contains("SingleEle") ) outputfile += "SingleEle_";
  else if (inputfile.Contains("Commissioning_") ) outputfile += "Commissioning_";
  else if (inputfile.Contains("JetHT_") ) outputfile += "JetHT_";
  else if (inputfile.Contains("DYJetsToLL_")) outputfile += "DYJetsToLL_";
  else if (inputfile.Contains("DoubleEG/")) outputfile += "DoubleEG_";
  else if (inputfile.Contains("MuonEG/")) outputfile += "MuonEG_";

  if (inputfile.Contains("2017A_PromptReco_v2_DCSJSON")) outputfile += "2017A_PromptReco_v2";
  if (inputfile.Contains("2017A_PromptReco_v3_DCSJSON")) outputfile += "2017A_PromptReco_v3";
  if (inputfile.Contains("2017B_PromptReco_v1_DCSJSON")) outputfile += "2017B_PromptReco_v1";
  if (inputfile.Contains("2017B_PromptReco_v2_DCSJSON")) outputfile += "2017B_PromptReco_v2";
  if (inputfile.Contains("2017C_PromptReco_v1_DCSJSON")) outputfile += "2017C_PromptReco_v1";
  if (inputfile.Contains("2017C_PromptReco_v2_DCSJSON")) outputfile += "2017C_PromptReco_v2";
  if (inputfile.Contains("2017D_PromptReco_v1")) outputfile += "2017D_PromptReco_v1";
  if (inputfile.Contains("2017E_PromptReco_v1")) outputfile += "2017E_PromptReco_v1";
  if (inputfile.Contains("2017F_PromptReco_v1")) outputfile += "2017F_PromptReco_v1";


  //if (mc)  outputfile = outputfile + "13TeV";
  //if (!mc) outputfile = outputfile + "13TeV_2017";

  if ( inputfile.Contains("eraA") ) outputfile += "_eraA";
  if ( inputfile.Contains("eraB") ) outputfile += "_eraB";
  //outputfile += "_StablePixel_TnP.root";
  //outputfile += "_DiEleMu_TriEleMu.root"; 
  outputfile += "_triLepTrig.root";
  
  TreeReader data(inputfile, "ggNtuplizer/EventTree");

  TFile *fo = TFile::Open(outputfile, "RECREATE");

  TTree *outtree = new TTree("outtree", "passing Id tree");

  inittree(outtree);

  Long64_t nTotal(0) ;
  Long64_t nPassTag(0) ;
  Long64_t nFill(0) ;
  Long64_t ntot(0);

  PUWeightCalculator puCalcGJ;
  if ( data.HasMC() ) {
    puCalcGJ.Init("external/PU_histo_13TeV_GoldenJSON_summer16_69000nb.root");
    //puCalcGJ.Init("external/PU_histo_13TeV_GoldenJSON_spring16_69000nb.root");
  }

  Long64_t TotalEntries = data.GetEntriesFast();  
  for ( Long64_t ev = 0; ev < TotalEntries; ev++) {   
    //for ( Long64_t ev = 0; ev < 30000000; ev++) {
    if (ev % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data.GetEntriesFast());

    data.GetEntry(ev);
    readggtree(data);

    nTotal++ ;

    if ( nEle < 1 && nMu < 1 ) continue ;

    vz_ = vz;
    vx_ = vx;
    vy_ = vy;
    nVtx_ = nVtx;
    isPVGood_ = isPVGood;
    run_ = run;
    event_ = event;
    lumi_ = lumis;

    //pfmet
    pfMet = pfMET;
    pfMetPhi = pfMETPhi;

    //PU weight for MC
    if (data.HasMC()) {
      float* puTrue = data.GetPtrFloat("puTrue");
      puweigj = (float) puCalcGJ.GetWeight(run, puTrue[1]); // in-time PU
      genWeight_ = (genWeight > 0) ? 1. : -1.;
      ntot += genWeight_;
    }

    //ele trigger
    //trig_Ele27_Eta2p1_WPTight = hlt>>1&1;
    trig_Ele27_WPTight = hlt>>4&1;
    trig_Ele23_Ele12_DZ = hlt>>5&1;
    trig_Ele23_Ele12 = hlt>>40&1;
    trig_DiSC30 = hlt>>30&1;

    //muon trigger
    trig_Mu17_Mu8_DZ = hlt >>19&1;
    trig_Mu17_TkMu8_DZ = hlt >>20&1;
    trig_TkMu17_TkMu8_DZ = hlt >>22&1;
    trig_Mu17_Mu8 = hlt >>18&1;
    trig_Mu17_TkMu8 = hlt >>21&1;
    trig_Mu17 = hlt>>29&1;
    trig_IsoMu27 = hlt>>19&1;
    trig_IsoTkMu24 = hlt>>20&1;

    //ele-mu trigger
    trig_Mu23_Ele12 = hlt>>45&1;
    trig_Mu23_Ele12_DZ = hlt>>46&1;
    trig_Mu12_Ele23 = hlt>>23&1;
    trig_Mu12_Ele23_DZ = hlt>>24&1;

    trig_Mu8_DiEle12 = hlt>>27&1;
    trig_Mu8_DiEle12_DZ = hlt>>52&1;
    trig_DiMu9_Ele9 = hlt>>28&1;
    trig_DiMu9_Ele9_DZ = hlt>>51&1;

    //trig_Pho30 = hltPho>>1&1;

    //--- electron loop ---//
    vector<int> acc_tag_ele, acc_probe_ele;
    vector<float> ele_probe_pt;
    int nZ = 0;

    for (int i=0; i<nEle; i++) {
      if (elePt[i] < 5.) continue;
      if (fabs(eleSCEta[i]) > 2.5) continue;
      if (fabs(eleSCEta[i]) > 1.4442 && fabs(eleSCEta[i]) < 1.566) continue;
      if (mc) {
        bool match = eleMatcher(data, i);
        if (!match) continue;
      }

      //if ( (eleID[i] >>3&1)==1 && elePt[i]>30 && fabs(eleSCEta[i])<2.1 && (eleFiredSingleTrgs[i] >>39&1)==1) {
      //tag ele matched to filter of HLT_Ele35
      if ( (eleID[i] >>3&1)==1 && (eleFiredSingleTrgs[i] >>39&1)==1) 
	acc_tag_ele.push_back(i);

      if ( (eleID[i] >>2&1)==1) {
	//if ( pass_cutbased_80X(i,1) ) {
	  acc_probe_ele.push_back(i);
	  ele_probe_pt.push_back(elePt[i]);
      }
    }

      //if (acc_tag_ele.size() < 1) continue;

    //-------------------------------------------------//
    // ---- muon loop --- //
    vector<int> acc_tag_mu, acc_probe_mu;
    vector<float> muon_tag_pt, muon_probe_pt;

    for (int i=0; i<nMu; i++) {
      if (muPt[i] < 5.) continue;
      if (fabs(muEta[i]) > 2.4) continue;

      if (mc) {
	bool match = muMatcher(data, i);
	if (!match) continue;
      }

      //tag muon matched to filter of IsoMu27, tight muon for tag
      if (muPt[i]>26 && (muIDbit[i]>>2 & 1)==1 && (muPFChIso[i] + TMath::Max(0., muPFPhoIso[i]+muPFNeuIso[i]-0.5 * muPFPUIso[i]))/muPt[i] < 0.15 && (muFiredTrgs[i]>>31&1)==1) { //for e-mu, die-mi
	//if (muPt[i]>26 && (muIDbit[i]>>2 & 1)==1 && (muPFChIso[i] + TMath::Max(0., muPFPhoIso[i]+muPFNeuIso[i]-0.5 * muPFPUIso[i]))/muPt[i] < 0.15) {
	acc_tag_mu.push_back(i);
	muon_tag_pt.push_back(muPt[i]);
      }
      //loose muon for probe
      if ( (muIDbit[i] >> 1 & 1) ==1 && (muPFChIso[i] + TMath::Max(0., muPFPhoIso[i]+muPFNeuIso[i]-0.5 * muPFPUIso[i]))/muPt[i] < 0.25) {
	acc_probe_mu.push_back(i);
	muon_probe_pt.push_back(muPt[i]);
      }
    }


    matching_EleLeg_Mu23Ele12 = 0;
    matching_EleLeg_Mu12Ele23 = 0;
    matching_EleLeg_Mu8DiEle12 = 0;
    matching_EleLeg_DiMu9Ele9 = 0;

    matching_MuLeg_Mu23Ele12 = 0;
    matching_MuLeg_Mu12Ele23 = 0;
    matching_MuLeg_Mu8DiEle12 = 0;
    matching_MuLeg_DiMu9Ele9 = 0;
    match_Mu9Ele9DZFilter = 0;

    //initial value
    Tag_elePt = -99.;
    Tag_eleEta = -99.;
    Tag_elePhi = -99.;
    Probe_elePt1 = - 99.;
    Probe_eleEta1 = -99.;
    Probe_elePhi1 = -99.;
    Probe_elePt2 = -99.;
    Probe_eleEta2 = -99.;
    Probe_elePhi2 = -99.;

    Tag_muPt = -99.;
    Tag_muEta = -99.;
    Tag_muPhi = -99.;
    Tag_matching_MuLeg_Mu23Ele12 = 0;
    Tag_matching_MuLeg_Mu12Ele23 = 0;
    Probe_muPt1 = - 99.;
    Probe_muEta1 = -99.;
    Probe_muPhi1 = -99.;
    Probe_muPt2 = -99.;
    Probe_muEta2 = -99.;
    Probe_muPhi2 = -99.;

    int siz_tag_ele = (int) acc_tag_ele.size();
    //use highest pt for tag
    if ( siz_tag_ele > 0) {
      Tag_elePt = elePt[acc_tag_ele[0]];
      Tag_eleEta = eleEta[acc_tag_ele[0]];
      Tag_elePhi = elePhi[acc_tag_ele[0]];
    }

    nProbe_ele = (int) acc_probe_ele.size();
    int ind_probe_ele[nProbe_ele];
    //cout << "n probe of electron" << nProbe_ele << endl;
    if ( nProbe_ele == 1) {
      Probe_elePt1 = elePt[acc_probe_ele[0]];
      Probe_eleEta1 = eleEta[acc_probe_ele[0]];
      Probe_elePhi1 = elePhi[acc_probe_ele[0]];
      matching_EleLeg_Mu23Ele12 = eleFiredSingleTrgs[acc_probe_ele[0]] >>24&1;
      matching_EleLeg_Mu12Ele23 = eleFiredSingleTrgs[acc_probe_ele[0]] >>38&1;
      matching_EleLeg_DiMu9Ele9 = eleFiredSingleTrgs[acc_probe_ele[0]] >>23&1;
      matching_EleLeg_Mu8DiEle12 = eleFiredSingleTrgs[acc_probe_ele[0]] >>25&1;
      match_Mu9Ele9DZFilter = eleFiredSingleTrgs[acc_probe_ele[0]] >>36&1;
    }
    else if (nProbe_ele > 1) {
      TMath::Sort(nProbe_ele, &ele_probe_pt.front(), ind_probe_ele); //sort in order of pt

      Probe_elePt1 = elePt[ind_probe_ele[0]];
      Probe_eleEta1 = eleEta[ind_probe_ele[0]];
      Probe_elePhi1 = elePhi[ind_probe_ele[0]];
      matching_EleLeg_Mu23Ele12 = eleFiredSingleTrgs[ind_probe_ele[0]] >>24&1;
      matching_EleLeg_Mu12Ele23 = eleFiredSingleTrgs[ind_probe_ele[0]] >>38&1;
      matching_EleLeg_DiMu9Ele9 = eleFiredSingleTrgs[ind_probe_ele[0]] >>23&1;
      match_Mu9Ele9DZFilter     = eleFiredSingleTrgs[ind_probe_ele[0]] >>36&1;
      matching_EleLeg_Mu8DiEle12 = eleFiredSingleTrgs[ind_probe_ele[0]] >>25&1;

      Probe_elePt2 = elePt[ind_probe_ele[1]];
      Probe_eleEta2 = eleEta[ind_probe_ele[1]];
      Probe_elePhi2 = elePhi[ind_probe_ele[1]];
      //matching_EleLeg_Mu8DiEle12 = eleFiredSingleTrgs[ind_probe_ele[1]] >>25&1;
    }

    int siz_tag = (int) acc_tag_mu.size();
    //use highest pt for tag
    if ( siz_tag > 0) {
      Tag_muPt = muPt[acc_tag_mu[0]];
      Tag_muEta = muEta[acc_tag_mu[0]];
      Tag_muPhi = muPhi[acc_tag_mu[0]];
      Tag_matching_MuLeg_Mu23Ele12 = muFiredTrgs[acc_tag_mu[0]] >>4&1;
      Tag_matching_MuLeg_Mu12Ele23 = muFiredTrgs[acc_tag_mu[0]] >>25&1;
    }

    //cout << "probe muon" << endl;
    nProbe_mu = (int) acc_probe_mu.size();
    int ind_probe_mu[nProbe_mu];
    //first two highest pt of probe
    if ( nProbe_mu == 1) {
      Probe_muPt1 = muPt[acc_probe_mu[0]];
      Probe_muEta1 = muEta[acc_probe_mu[0]];
      Probe_muPhi1 = muPhi[acc_probe_mu[0]];
      matching_MuLeg_Mu23Ele12 = muFiredTrgs[acc_probe_mu[0]] >>4&1;
      matching_MuLeg_Mu12Ele23 = muFiredTrgs[acc_probe_mu[0]] >>25&1;
      matching_MuLeg_Mu8DiEle12 = muFiredTrgs[acc_probe_mu[0]] >>5&1;
      matching_MuLeg_DiMu9Ele9 = muFiredTrgs[acc_probe_mu[0]] >>3&1;
    }
    else if (nProbe_mu > 1) {
      TMath::Sort(nProbe_mu, &muon_probe_pt.front(), ind_probe_mu); //sort in order of pt
      //cout << "get muon probe 1 and 2" << endl;

      Probe_muPt1 = muPt[ind_probe_mu[0]];
      Probe_muEta1 = muEta[ind_probe_mu[0]];
      Probe_muPhi1 = muPhi[ind_probe_mu[0]];
      matching_MuLeg_Mu23Ele12 = muFiredTrgs[ind_probe_mu[0]] >>4&1;
      matching_MuLeg_Mu12Ele23 = muFiredTrgs[ind_probe_mu[0]] >>25&1;
      matching_MuLeg_Mu8DiEle12 = muFiredTrgs[ind_probe_mu[0]] >>5&1;
      matching_MuLeg_DiMu9Ele9 = muFiredTrgs[ind_probe_mu[0]] >>3&1;

      Probe_muPt2 = muPt[ind_probe_mu[1]];
      Probe_muEta2 = muEta[ind_probe_mu[1]];
      Probe_muPhi2 = muPhi[ind_probe_mu[1]];
      //matching_MuLeg_DiMu9Ele9 = muFiredTrgs[ind_probe_mu[1]] >>3&1;
    }

    //cout << "filling tree" << endl;
    outtree->Fill();


  }// end of loop
  cout << "total events = " << nTotal << endl ;
  //cout << "pas probes     = " << nPassProbe << endl ;
  cout << "filled       = " << nFill << endl ;
  cout << "endloop " << endl ;
  
  fo->Write();
  fo->Close();
  delete fo;

}
