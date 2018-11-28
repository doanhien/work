/************************
script for trigger efficiency in event
count events passing given trigger from events passing selection
events: ele passing Id+Iso, muon passing Id+Iso
************************/

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
#include "tree.h"
#include "ElectronSelection.h"
#include "MuonSelection.h"
//#include "cat.h"
#include "puweicalc.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>


void xAna(TString inputfile = "HiForest.root", TString outputfile = "XZg.root", bool mc = true, bool ele = true)
{


  outputfile = "minitrees/";

  if (inputfile.Contains("SingleMu") ) outputfile += "SingleMu_";
  if (inputfile.Contains("Commissioning_") ) outputfile += "Commissioning_"; 
  if (inputfile.Contains("JetHT_") ) outputfile += "JetHT_";
  if (inputfile.Contains("DYJetsToLL_")) outputfile += "DYJetsToLL_";
  if (inputfile.Contains("DoubleEG/")) outputfile += "DoubleEG_";
  if (inputfile.Contains("MuonEG/")) outputfile += "MuonEG_";
  if (inputfile.Contains("SingleEle")) outputfile += "SingleEle_";
  if (inputfile.Contains("SingleMuon/")) outputfile += "SingleMuon_";
  if (inputfile.Contains("MET_Run2017A_PromptReco_v2")) outputfile += "Met_Run2017A_PromptReco_v2";
  if (inputfile.Contains("MET_Run2017A_PromptReco_v3")) outputfile += "Met_Run2017A_PromptReco_v3";
  if (inputfile.Contains("MET_Run2017B_PromptReco_v1")) outputfile += "Met_Run2017B_PromptReco_v1";
  if (inputfile.Contains("MET_Run2017B_PromptReco_v2")) outputfile += "Met_Run2017B_PromptReco_v2";
  if (inputfile.Contains("MET_Run2017C_PromptReco_v1")) outputfile += "Met_Run2017C_PromptReco_v1";
  if (inputfile.Contains("MET_Run2017C_PromptReco_v2")) outputfile += "Met_Run2017C_PromptReco_v2";

  if (inputfile.Contains("_2017A_PromptReco_v2") ) outputfile += "Run2017A_PromptReco_v2_DCSJSON";
  if (inputfile.Contains("_2017A_PromptReco_v3") ) outputfile += "Run2017A_PromptReco_v3_DCSJSON";
  if (inputfile.Contains("_2017B_PromptReco_v1") ) outputfile += "Run2017B_PromptReco_v1_DCSJSON";
  if (inputfile.Contains("_2017B_PromptReco_v2") ) outputfile += "Run2017B_PromptReco_v2_DCSJSON";
  if (inputfile.Contains("_2017C_PromptReco_v1") ) outputfile += "Run2017C_PromptReco_v1_DCSJSON";
  if (inputfile.Contains("_2017C_PromptReco_v2") ) outputfile += "Run2017C_PromptReco_v2_DCSJSON";

  //if (mc)  outputfile = outputfile + "13TeV";
  //if (!mc) outputfile = outputfile + "13TeV_2017";

  //outputfile += "_eraB_StablePixel.root";
  outputfile += ".root";

  TreeReader data(inputfile, "ggNtuplizer/EventTree");

  TFile *fo = TFile::Open(outputfile, "RECREATE");

  TTree *outtree = new TTree("outtree", "output tree");

  inittree(outtree);

  Long64_t nTotal(0) ;
  Long64_t nPassHLT(0) ;
  Long64_t nPassTag(0) ;
  Long64_t nFill(0) ;
  Long64_t ntot(0);


  PUWeightCalculator puCalcGJ;
  if ( data.HasMC() ) {
    puCalcGJ.Init("external/PU_histo_13TeV_GoldenJSON_summer16_69000nb.root");
    //puCalcGJ.Init("external/PU_histo_13TeV_GoldenJSON_spring16_69000nb.root");
  }

  TFile *fpu = new TFile("PUweight.root", "READ");
  TH1F *hpu = (TH1F*) fpu->Get("puwei");

  Long64_t TotalEntries = data.GetEntriesFast();  
  for ( Long64_t ev = 0; ev < TotalEntries; ev++) {   
    //for ( Long64_t ev = 0; ev < 10000000; ev++) {
    if (ev % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data.GetEntriesFast());

    data.GetEntry(ev);
    readggtree(data);

    nTotal++ ;

    if ( nEle < 1 && nMu < 1) continue ;

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
      //float* puTrue = data.GetPtrFloat("puTrue");
      //puweigj = (float) puCalcGJ.GetWeight(run, puTrue[1]); // in-time PU

      //pu from early data
      int ibin = hpu->FindBin(nVtx);
      puweigj = hpu->GetBinContent(ibin);

      genWeight_ = (genWeight > 0) ? 1. : -1.;
      ntot += genWeight_;
    }

    //electron trigger for 92X
    trig_Ele27_WPTight = hlt>>4&1;
    trig_Ele23_Ele12_DZ = hlt>>5&1;
    trig_Ele23_Ele12 = hlt>>40&1;
    trig_Ele16_Ele12_Ele8 = hlt>>6;

    //HLT_DiSC30_18_EIso_AND_HE_Mass70_v
    trig_DiSC30 = hlt>>30&1;

    //muon trigger
    trig_Mu17_Mu8_DZ = hlt >>19&1;
    trig_Mu17_TkMu8_DZ = hlt >>20&1;
    trig_TkMu17_TkMu8_DZ = hlt >>22&1;
    trig_Mu17_Mu8 = hlt >>18&1;
    trig_Mu17_TkMu8 = hlt >>21&1;
    trig_Mu17 = hlt>>29&1;
    trig_IsoMu24 = hlt>>10&1;

    //mu trigger 
    //trig_IsoMu24 = hlt>>19&1;

    // ele - mu trigger
    trig_Mu23_Ele12_DZ = hlt>>26&1;
    trig_Mu12_Ele23_DZ = hlt>>24&1;

    trig_Mu8_DiEle12 = hlt>>27&1;
    trig_Mu8_DiEle12_DZ = hlt>>52&1;
    trig_DiMu9_Ele9 = hlt>>28&1;
    trig_DiMu9_Ele9_DZ = hlt>>51&1;

    //met trigger
    trig_PFMet120_PFMHT120 = hltJet>>27&1;

    TLorentzVector pair,ele[10], muon[5];
    ne = 0, nmu = 0;

    for (int i=0; i<nEle; i++) {
      if (elePt[i] < 7.) continue;
      if (fabs(eleSCEta[i]) > 2.5) continue;
      //if (fabs(eleSCEta[i]) > 1.4442 && fabs(eleSCEta[i]) < 1.566) continue;
      if (mc) {
        bool match = eleMatcher(data, i);
        if (!match) continue;
      }
      //if ( !pass_cutbased_80X(i,1)) continue;  //1- WP Medium, 2- Tight
      if ( (eleID[i]>>2&1) !=1) continue;

      passSieie[ne] = 1;
      passdEta[ne] = 1;
      passdPhi[ne] = 1;
      passHoE[ne] = 1;
      passEoPInv[ne] = 1;
      passConvVeto[ne] = 1;
      passMissHits[ne] = 1;

      if (fabs(eleSCEta[i]) < 1.5) {
	if (eleSigmaIEtaIEta_Full5x5[i] > 0.00998) passSieie[ne] = 0;
	if (fabs(eledEtaseedAtVtx[i]) > 0.00311) passdEta[ne] = 0;
	if (fabs(eledPhiAtVtx[i]) > 0.103) passdPhi[ne] = 0;
	if (eleHoverE[i] > 0.253) passHoE[ne] = 0;
	if (eleEoverPInv[i] > 0.134) passEoPInv[ne] = 0;
	if (eleConvVeto[i] != 1) passConvVeto[ne] = 0;
	if (eleMissHits[i] > 1) passMissHits[ne] = 0;
      }
      else {
	if (eleSigmaIEtaIEta_Full5x5[i] > 0.0298) passSieie[ne] = 0;
	if (fabs(eledEtaseedAtVtx[i]) > 0.00609) passdEta[ne] = 0;
	if (fabs(eledPhiAtVtx[i]) > 0.045) passdPhi[ne] = 0;
	if (eleHoverE[i] > 0.0878) passHoE[ne] = 0;
	if (eleEoverPInv[i] > 0.13) passEoPInv[ne] = 0;
	if (eleConvVeto[i] != 1) passConvVeto[ne] = 0;
	if (eleMissHits[i] > 1) passMissHits[ne] = 0;
      }

      elePt_[ne] = elePt[i];
      eleEta_[ne] = eleEta[i];
      elePhi_[ne] = elePhi[i];
      eleEoP_[ne] = eleEoverP[i];
      eleHoE_[ne] = eleHoverE[i];
      eleEoPInv_[ne] = eleEoverPInv[i];
      eledEta_[ne] = eledEtaAtVtx[i];
      eledPhi_[ne] = eledPhiAtVtx[i];
      eleSigIetaIeta_[ne] = eleSigmaIEtaIEta_Full5x5[i];
      eleConvVeto_[ne] = eleConvVeto[i];
      eleMissHits_[ne] = eleMissHits[i];
      eleCh_[ne] = eleCharge[i];
      eleR9_[ne] = eleR9[i];

      eleChIso_[ne] = elePFChIso[i];
      elePhoIso_[ne] = elePFPhoIso[i];
      eleNeuIso_[ne] = elePFNeuIso[i];
      eleIso_[ne] = (elePFChIso[i] + elePFPhoIso[i] + elePFNeuIso[i]) / elePt[i];

      eleCh_[ne] = eleCharge[i];
      ele[ne].SetPtEtaPhiM(elePt[i],eleEta[i],elePhi[i], 0.511*0.001);

      //matching trigger
      matchTrigEle27_[ne] = eleFiredSingleTrgs[i] >>1&1;
      matchEle23Leg_[ne] = eleFiredDoubleTrgs[i] >>2&1;
      matchEle12Leg_[ne] = eleFiredDoubleTrgs[i] >>3&1;
      matchDZ_Ele23_Ele12_[ne] = eleFiredDoubleTrgs[i] >>4&1;
      matchEle16Leg_TriEle_[ne] = eleFiredDoubleTrgs[i] >>9&1;
      matchEle12Leg_TriEle_[ne] = eleFiredDoubleTrgs[i] >>10&1;
      matchEle8Leg_TriEle_[ne] = eleFiredDoubleTrgs[i] >>11&1;


      ne++;
    } //endl of electron loop

    for (int i=0; i<nMu; i++) {

      if (muPt[i] < 5.) continue;
      if (fabs(muEta[i]) > 2.4) continue;

      if (mc) {
        bool match = muMatcher(data, i);
        if (!match) continue;
      }

      if ((muIDbit[i] >> 1 & 1) == 0) continue;  // >>2- tight id, >>1-medium id
      if ((muPFChIso[i] + TMath::Max(0., muPFPhoIso[i] + muPFNeuIso[i] -0.5 * muPFPUIso[i]))/muPt[i] > 0.25) continue; //loose iso: 0.25; tight iso:0.15

      muPt_[nmu] = muPt[i];
      muEta_[nmu] = muEta[i];
      muPhi_[nmu] = muPhi[i];

      muChIso_[nmu] = muPFChIso[i];
      muPhoIso_[nmu] = muPFPhoIso[i];
      muNeuIso_[nmu] = muPFNeuIso[i];
      muIso_[nmu] = (muPFChIso[i] + muPFPhoIso[i] + muPFNeuIso[i]) / muPt[i];

      //matching trigger
      //MatchTrig24_[nmu] = muFiredSingleTrgs[i] >>11&1;

      muCh_[nmu] = muCharge[i];
      muon[nmu].SetPtEtaPhiM(muPt[i],muEta[i],muPhi[i], 105.6*0.001);
      nmu++;
    }

    if (ne<1 && nmu<1) continue;

    Ztype = 0;
    if (ne > 1) {
      pair = ele[0] + ele[1];
      Ztype = 11;
      Zmass = pair.M();
      Zpt = pair.Pt();
      Zy = pair.Rapidity();
      Zeta = pair.Eta();
      Zphi = pair.Phi();
      Zcharge = eleCh_[0] + eleCh_[1];
      dPhidiLep = deltaPhi(ele[0].Phi(),ele[1].Phi());
      dRdiLep = deltaR(ele[0].Eta(), ele[0].Phi(),ele[1].Eta(), ele[1].Phi());
    }
    else if (nmu > 1) {
      pair = muon[0] + muon[1];
      Ztype =13;
      Zmass = pair.M();
      Zpt = pair.Pt();
      Zy = pair.Rapidity();
      Zeta = pair.Eta();
      Zphi = pair.Phi();
      Zcharge = muCh_[0] + muCh_[1];
      dPhidiLep = deltaPhi(muon[0].Phi(),muon[1].Phi());
      dRdiLep = deltaR(muon[0].Eta(), muon[0].Phi(),muon[1].Eta(), muon[1].Phi());
    }
    else {
      Ztype = -1;
      Zmass = -99.;
      Zpt = -99.;
      Zy = -99.;
      Zeta = -99.;
      Zphi  = -99;
      Zcharge = -99;
    }

    //if (Zmass < 50.) continue;
    outtree->Fill();

  }// end of loop
  cout << "total events = " << nTotal << endl ;
  cout << "pass HLT     = " << nPassHLT << endl ;
  cout << "pas tags     = " << nPassTag << endl ;
  cout << "filled       = " << nFill << endl ;
  cout << "endloop " << endl ;
  
  fo->Write();
  fo->Close();
  delete fo;

}
