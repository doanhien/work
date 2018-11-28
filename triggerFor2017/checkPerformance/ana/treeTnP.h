#ifndef tree_h
#define tree_h

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TMath.h>
#include <TLorentzVector.h>

#include <iostream>

float deltaPhi(float phi1, float phi2) {

  float dPhi = phi1 - phi2;
  if (dPhi > TMath::Pi()) dPhi -= 2*TMath::Pi();
  if (dPhi <= -TMath::Pi()) dPhi += 2* TMath::Pi();

  return dPhi;

}

float deltaR(float eta1, float phi1, float eta2, float phi2) {

  float dEta = eta1 - eta2;
  float dPhi = phi1 - phi2;
  if (dPhi > TMath::Pi()) dPhi -= 2*TMath::Pi();
  if (dPhi <= -TMath::Pi()) dPhi += 2* TMath::Pi();

  return sqrt(pow(dEta,2) + pow (dPhi,2));

}

const int maxPar = 50;

Int_t run_;
Int_t lumi_;
Long64_t event_;
Int_t nVtx_;
Float_t vx_;
Float_t vy_;
Float_t vz_;
Float_t genWeight_;
Float_t puweigj;
Int_t   isPVGood_;

Int_t     trig_Ele27_WPTight;
Int_t     trig_Ele35_WPTight;
Int_t     trig_Ele38_WPTight;
Int_t     trig_Ele40_WPTight;
Int_t     trig_DoubleEle33;
Int_t     trig_Ele23_Ele12;
Int_t     trig_Ele23_Ele12_DZ;
Int_t     trig_DiSC30;

Int_t     trig_Mu23_Ele12;
Int_t     trig_Mu23_Ele12_DZ;
Int_t     trig_Mu12_Ele23;
Int_t     trig_Mu12_Ele23_DZ ;
Int_t     trig_Mu8_DiEle12;
Int_t     trig_Mu8_DiEle12_DZ;
Int_t     trig_DiMu9_Ele9;
Int_t     trig_DiMu9_Ele9_DZ;

Int_t     trig_Mu17;
Int_t     trig_IsoMu27;
Int_t     trig_IsoTkMu24;
Int_t     trig_Mu17_Mu8;
Int_t     trig_Mu17_TkMu8;
Int_t     trig_TkMu17_TkMu8;
Int_t     trig_Mu17_Mu8_DZ;
Int_t     trig_Mu17_TkMu8_DZ;
Int_t     trig_TkMu17_TkMu8_DZ;

Int_t     trig_doublePho60;
Int_t     trig_Ele17_Ele12;

Float_t   pfMet;
Float_t   pfMetPhi;
Int_t     ne;
Int_t     nTag;
Int_t     TagIndex;
Float_t   Tag_Pt;
Float_t   Tag_Eta;
Float_t   Tag_Phi;
Float_t   Tag_SCEta;
Float_t   Tag_sieie;
Float_t   Tag_chIso;
Float_t   Tag_neuIso;
Float_t   Tag_phoIso;
Float_t   Tag_hoe;
Float_t   Tag_eop;
Float_t   Tag_ch;
Float_t   Tag_dEtaSeed;
Float_t   Tag_dEtaVtx;
Float_t   Tag_dPhiVtx;
Float_t   Tag_IDMVA;
Int_t     ele1_MatchTrig27;
Float_t   mW_Tag_Met;

Int_t     nProbe;
Int_t     ProbeIndex;
Int_t     nPassId;
Int_t     nFailId;
Int_t     nPassIdIso;
Int_t     nFailIdIso;
Float_t   Probe_Pt;
Float_t   Probe_Eta;
Float_t   Probe_Phi;
Float_t   Probe_SCEta;
Float_t   Probe_sieie;
Float_t   Probe_chIso;
Float_t   Probe_neuIso;
Float_t   Probe_phoIso;
Float_t   Probe_hoe;
Float_t   Probe_eop;
Float_t   Probe_ch;
Float_t   Probe_dEtaSeed;
Float_t   Probe_dEtaVtx;
Float_t   Probe_dPhiVtx;
Float_t   Probe_miniIso;
Int_t   LooseID;
Int_t   MediumID;
Int_t   TightID;
Int_t   LooseIDOnly;
Int_t   passingMVA;
Int_t   passingHLT;

Float_t Tag_elePt;
Float_t Tag_eleEta;
Float_t Tag_elePhi;

Int_t  nProbe_ele;
Float_t Probe_elePt1;
Float_t Probe_elePt2;
Float_t Probe_eleEta1;
Float_t Probe_eleEta2;
Float_t Probe_elePhi1;
Float_t Probe_elePhi2;
Int_t   matching_EleLeg_Mu23Ele12;
Int_t   matching_EleLeg_Mu12Ele23;
Int_t   matching_EleLeg_Mu8DiEle12;
Int_t   matching_EleLeg_DiMu9Ele9;
Int_t   match_Mu9Ele9DZFilter;

Float_t Tag_muPt;
Float_t Tag_muEta;
Float_t Tag_muPhi;
Int_t   Tag_matching_MuLeg_Mu23Ele12;
Int_t   Tag_matching_MuLeg_Mu12Ele23;

Int_t nProbe_mu;
Float_t Probe_muPt1;
Float_t Probe_muPt2;
Float_t Probe_muEta1;
Float_t Probe_muEta2;
Float_t Probe_muPhi1;
Float_t Probe_muPhi2;
Int_t   matching_MuLeg_Mu23Ele12;
Int_t   matching_MuLeg_Mu12Ele23;
Int_t   matching_MuLeg_Mu8DiEle12;
Int_t   matching_MuLeg_DiMu9Ele9;


void inittree(TTree* tree) {

  tree->Branch("run",         &run_);
  tree->Branch("lumi",        &lumi_);
  tree->Branch("event",       &event_);
  tree->Branch("nVtx",        &nVtx_);
  tree->Branch("vx",          &vx_);
  tree->Branch("vy",          &vy_);
  tree->Branch("vz",          &vz_);
  tree->Branch("isPVGood",    &isPVGood_);
  tree->Branch("genWeight",   &genWeight_);
  tree->Branch("puweigj",     &puweigj);

  //electron trigger
  tree->Branch("trig_Ele27_WPTight",      &trig_Ele27_WPTight);
  tree->Branch("trig_Ele35_WPTight",      &trig_Ele35_WPTight);
  tree->Branch("trig_Ele38_WPTight",      &trig_Ele38_WPTight);
  tree->Branch("trig_Ele40_WPTight",      &trig_Ele40_WPTight);
  tree->Branch("trig_DiSC30",             &trig_DiSC30);

  tree->Branch("trig_Ele23_Ele12",        &trig_Ele23_Ele12);
  tree->Branch("trig_Ele23_Ele12_DZ",     &trig_Ele23_Ele12_DZ);
  tree->Branch("trig_DoubleEle33",        &trig_DoubleEle33);

  tree->Branch("trig_Mu12_Ele23",         &trig_Mu12_Ele23);
  tree->Branch("trig_Mu12_Ele23_DZ",      &trig_Mu12_Ele23_DZ);
  tree->Branch("trig_Mu23_Ele12",        &trig_Mu23_Ele12);
  tree->Branch("trig_Mu23_Ele12_DZ",      &trig_Mu23_Ele12_DZ);
  tree->Branch("trig_Mu8_DiEle12",        &trig_Mu8_DiEle12);
  tree->Branch("trig_Mu8_DiEle12_DZ",     &trig_Mu8_DiEle12_DZ);
  tree->Branch("trig_DiMu9_Ele9",         &trig_DiMu9_Ele9);
  tree->Branch("trig_DiMu9_Ele9_DZ",      &trig_DiMu9_Ele9_DZ);

  //muon trigger
  tree->Branch("trig_Mu17",               &trig_Mu17);
  tree->Branch("trig_IsoMu27",            &trig_IsoMu27);
  tree->Branch("trig_IsoTkMu24",          &trig_IsoTkMu24);
  tree->Branch("trig_Mu17_Mu8",           &trig_Mu17_Mu8);
  tree->Branch("trig_Mu17_TkMu8",         &trig_Mu17_TkMu8);
  tree->Branch("trig_TkMu17_TkMu8",       &trig_TkMu17_TkMu8);
  tree->Branch("trig_Mu17_Mu8_DZ",        &trig_Mu17_Mu8_DZ);
  tree->Branch("trig_Mu17_TkMu8_DZ",      &trig_Mu17_TkMu8_DZ);
  tree->Branch("trig_TkMu17_TkMu8_DZ",    &trig_TkMu17_TkMu8_DZ);

  //ele-muon trigger
  tree->Branch("trig_Mu12_Ele23_DZ",      &trig_Mu12_Ele23_DZ);
  tree->Branch("trig_Mu23_Ele12_DZ",      &trig_Mu23_Ele12_DZ);

  tree->Branch("pfMet",        &pfMet);
  tree->Branch("pfMetPhi",     &pfMetPhi);
  tree->Branch("ne",           &ne);
  tree->Branch("nTag",         &nTag);
  tree->Branch("TagIndex",     &TagIndex);
  tree->Branch("Tag_Pt",       &Tag_Pt);
  tree->Branch("Tag_Eta",      &Tag_Eta);
  tree->Branch("Tag_Phi",      &Tag_Phi);
  tree->Branch("Tag_SCEta",    &Tag_SCEta);
  tree->Branch("Tag_sieie",    &Tag_sieie);
  tree->Branch("Tag_chIso",    &Tag_chIso);
  tree->Branch("Tag_neuIso",   &Tag_neuIso);
  tree->Branch("Tag_phoIso",   &Tag_phoIso);
  tree->Branch("Tag_hoe",      &Tag_hoe);
  tree->Branch("Tag_eop",      &Tag_eop);
  tree->Branch("Tag_ch",       &Tag_ch);
  tree->Branch("Tag_dEtaSeed",     &Tag_dEtaSeed);
  tree->Branch("Tag_dEtaVtx",      &Tag_dEtaVtx);
  tree->Branch("Tag_dPhiVtx",      &Tag_dPhiVtx);
  tree->Branch("Tag_IDMVA",        &Tag_IDMVA);
  tree->Branch("ele1_MatchTrig27", &ele1_MatchTrig27);
  tree->Branch("mW_Tag_Met",       &mW_Tag_Met);

  tree->Branch("nProbe",         &nProbe);
  tree->Branch("nPassId",        &nPassId);
  tree->Branch("nFailId",        &nFailId);
  tree->Branch("nPassIdIso",     &nPassIdIso);
  tree->Branch("nFailIdIso",     &nFailIdIso);
  tree->Branch("ProbeIndex",     &ProbeIndex);
  tree->Branch("Probe_Pt",       &Probe_Pt);
  tree->Branch("Probe_Eta",      &Probe_Eta);
  tree->Branch("Probe_Phi",      &Probe_Phi);
  tree->Branch("Probe_SCEta",    &Probe_SCEta);
  tree->Branch("Probe_sieie",    &Probe_sieie);
  tree->Branch("Probe_chIso",    &Probe_chIso);
  tree->Branch("Probe_neuIso",   &Probe_neuIso);
  tree->Branch("Probe_phoIso",   &Probe_phoIso);
  tree->Branch("Probe_hoe",      &Probe_hoe);
  tree->Branch("Probe_eop",      &Probe_eop);
  tree->Branch("Probe_dEtaSeed", &Probe_dEtaSeed);
  tree->Branch("Probe_dEtaVtx",  &Probe_dEtaVtx);
  tree->Branch("Probe_dPhiVtx",  &Probe_dPhiVtx);
  tree->Branch("Probe_ch",       &Probe_ch);
  tree->Branch("Probe_miniIso",  &Probe_miniIso);
  tree->Branch("LooseID",        &LooseID);
  tree->Branch("MediumID",       &MediumID);
  tree->Branch("TightID",        &TightID);
  tree->Branch("LooseIDOnly",    &LooseIDOnly);
  tree->Branch("passingMVA",     &passingMVA);
  tree->Branch("passingHLT",     &passingHLT);

  tree->Branch("Tag_elePt",        &Tag_elePt);
  tree->Branch("Tag_eleEta",       &Tag_eleEta);
  tree->Branch("Tag_elePhi",       &Tag_elePhi);

  tree->Branch("nProbe_ele",       &nProbe_ele);
  tree->Branch("Probe_elePt1",     &Probe_elePt1);
  tree->Branch("Probe_eleEta1",    &Probe_eleEta1);
  tree->Branch("Probe_elePhi1",    &Probe_elePhi1);
  tree->Branch("Probe_elePt2",     &Probe_elePt2);
  tree->Branch("Probe_eleEta2",    &Probe_eleEta2);
  tree->Branch("Probe_elePhi2",    &Probe_elePhi2);
  tree->Branch("matching_EleLeg_Mu23Ele12",     &matching_EleLeg_Mu23Ele12);
  tree->Branch("matching_EleLeg_Mu12Ele23",     &matching_EleLeg_Mu12Ele23);
  tree->Branch("matching_EleLeg_Mu8DiEle12",    &matching_EleLeg_Mu8DiEle12);
  tree->Branch("matching_EleLeg_DiMu9Ele9",     &matching_EleLeg_DiMu9Ele9);
  tree->Branch("match_Mu9Ele9DZFilter",         &match_Mu9Ele9DZFilter);

  tree->Branch("Tag_muPt",        &Tag_muPt);
  tree->Branch("Tag_muEta",       &Tag_muEta);
  tree->Branch("Tag_muPhi",       &Tag_muPhi);
  tree->Branch("Tag_matching_MuLeg_Mu12Ele23",  &Tag_matching_MuLeg_Mu12Ele23);
  tree->Branch("Tag_matching_MuLeg_Mu23Ele12",  &Tag_matching_MuLeg_Mu23Ele12);

  tree->Branch("nProbe_mu",       &nProbe_mu);
  tree->Branch("Probe_muPt1",     &Probe_muPt1);
  tree->Branch("Probe_muEta1",    &Probe_muEta1);
  tree->Branch("Probe_muPhi1",    &Probe_muPhi1);
  tree->Branch("Probe_muPt2",     &Probe_muPt2);
  tree->Branch("Probe_muEta2",    &Probe_muEta2);
  tree->Branch("Probe_muPhi2",    &Probe_muPhi2);
  tree->Branch("matching_MuLeg_Mu23Ele12",     &matching_MuLeg_Mu23Ele12);
  tree->Branch("matching_MuLeg_Mu12Ele23",     &matching_MuLeg_Mu12Ele23);
  tree->Branch("matching_MuLeg_Mu8DiEle12",    &matching_MuLeg_Mu8DiEle12);
  tree->Branch("matching_MuLeg_DiMu9Ele9",     &matching_MuLeg_DiMu9Ele9);

}

Int_t     run;
Long64_t  event;
Int_t     lumis;
Bool_t    isData;
ULong64_t hlt;
ULong64_t hltPho;
Int_t nVtx;
Float_t vz;
Float_t vx; 
Float_t vy;
Int_t   isPVGood;
Float_t genWeight;
Float_t rho;
Float_t pfMET;
Float_t pfMETPhi;

Int_t nEle ;
Int_t* eleCharge;
Int_t* eleChargeConsistent;
Float_t *eleEn;
Float_t* elePt ;
Float_t* eleEta ;
Float_t* elePhi ;
Float_t* eleHoverE ;
Float_t* eleEoverP ;
Float_t* eleEoverPInv ;
Float_t* eleSigmaIEtaIEta_Full5x5 ;
Float_t* eleSigmaIPhiIPhi;
Int_t* eleMissHits ;
Float_t* eleD0 ;
Float_t* eleDz ;
Float_t* eledEtaAtVtx ;
Float_t* eledPhiAtVtx ;
Float_t* eledEtaseedAtVtx ;
Int_t* eleConvVeto;
Int_t* eleEcalDrivenSeed;
Float_t* eleE1x5Full5x5;
Float_t* eleE2x5Full5x5;
Float_t* eleE5x5Full5x5;
Float_t* eleSCEta ;
Float_t* eleSCPhi;
Float_t* eleSCEtaWidth;
Float_t* eleSCPhiWidth;
Float_t* eleSCEn ;
Float_t* elePFChIso ;
Float_t* elePFPhoIso ;
Float_t* elePFNeuIso ;
Float_t* elePFPUIso ;
Float_t* elePFMiniIso ;
Short_t* eleID;
Float_t* eleIDMVA;
Float_t* eleIDMVANonTrg;
//Int_t*   eleFiredSingleTrgs;
Long64_t*   eleFiredSingleTrgs;
Long64_t*   eleFiredDoubleTrgs;
Float_t* elePFClusEcalIso;
Float_t* elePFClusHcalIso;
Float_t* eleDr03TkSumPt;

Int_t nPho;
Float_t* phoEt;
Float_t* phoSCEta;
Float_t* phoSCPhi;
Float_t* phoEta;
Float_t* phoPhi;
Float_t* phoCalibEt;
Float_t* phoHoverE;
Float_t* phoSigmaIEtaIEtaFull5x5;
Float_t* phoR9;
Float_t* phoPFChIso;
Float_t* phoPFNeuIso;
Float_t* phoPFPhoIso;
Int_t* phohasPixelSeed;
Int_t* phoEleVeto;
Short_t* phoIDbit;
Float_t* phoIDMVA;


Int_t    nMu;
Float_t* muPt;
Float_t* muEta;
Float_t* muPhi;
Int_t*   muCh;
Float_t* muChi2NDF;
Int_t*   muMuonHits;
Int_t*   muStations;
Int_t*   muTrkLayers;
Int_t*   muPixelHits;
Float_t* muInnerD0;
Float_t* muInnerDz;
Float_t* muPFChIso;
Float_t* muPFPhoIso;
Float_t* muPFNeuIso;
Float_t* muPFPUIso;
Int_t*   muCharge;
Short_t* muIDbit;
Long64_t* muFiredTrgs;

//gen level
Int_t nMC;
Int_t* mcPID;
Int_t* mcMomPID;
Float_t* mcPt;
Float_t* mcPhi;
Float_t* mcEta;
UShort_t* mcStatusFlag;


void readggtree(TreeReader &data) {

  run = data.GetInt("run");
  event = data.GetLong64("event");
  lumis = data.GetInt("lumis");
  hlt = data.GetLong64("HLTEleMuX");
  hltPho = data.GetLong64("HLTPho");
  nVtx = data.GetInt("nVtx");
  vx = data.GetFloat("vtx");
  vy = data.GetFloat("vty");
  vz = data.GetFloat("vtz");
  isPVGood = data.GetBool("isPVGood");
  rho = data.GetFloat("rho");
  pfMET = data.GetFloat("pfMET");
  pfMETPhi = data.GetFloat("pfMETPhi");
  nEle = data.GetInt("nEle");
  eleCharge = data.GetPtrInt("eleCharge");
  //elePt = data.GetPtrFloat("eleCalibPt");
  elePt = data.GetPtrFloat("elePt");
  //eleEn = data.GetPtrFloat("eleCalibEn");
  eleEn = data.GetPtrFloat("eleEn");
  eleEta = data.GetPtrFloat("eleEta");
  elePhi = data.GetPtrFloat("elePhi");
  eleHoverE = data.GetPtrFloat("eleHoverE");
  eleEoverP = data.GetPtrFloat("eleEoverP");
  eleEoverPInv = data.GetPtrFloat("eleEoverPInv");
  eleSigmaIEtaIEta_Full5x5 = data.GetPtrFloat("eleSigmaIEtaIEtaFull5x5");
  eleMissHits = data.GetPtrInt("eleMissHits");
  eleD0 = data.GetPtrFloat("eleD0");
  eleDz = data.GetPtrFloat("eleDz");
  eledEtaAtVtx = data.GetPtrFloat("eledEtaAtVtx");
  eledPhiAtVtx = data.GetPtrFloat("eledPhiAtVtx");
  eledEtaseedAtVtx = data.GetPtrFloat("eledEtaseedAtVtx");
  eleConvVeto = data.GetPtrInt("eleConvVeto");
  eleEcalDrivenSeed = data.GetPtrInt("eleEcalDrivenSeed");
  eleE1x5Full5x5 = data.GetPtrFloat("eleE1x5Full5x5");
  eleE2x5Full5x5 = data.GetPtrFloat("eleE2x5Full5x5");
  eleE5x5Full5x5 = data.GetPtrFloat("eleE5x5Full5x5");
  eleSCEta = data.GetPtrFloat("eleSCEta");
  eleSCPhi = data.GetPtrFloat("eleSCPhi");
  eleSCEtaWidth = data.GetPtrFloat("eleSCEtaWidth");
  eleSCPhiWidth = data.GetPtrFloat("eleSCPhiWidth");
  eleSCEn = data.GetPtrFloat("eleSCEn");
  elePFChIso = data.GetPtrFloat("elePFChIso");
  elePFPhoIso = data.GetPtrFloat("elePFPhoIso");
  elePFNeuIso = data.GetPtrFloat("elePFNeuIso");
  elePFMiniIso = data.GetPtrFloat("elePFMiniIso");
  eleID = data.GetPtrShort("eleIDbit");
  eleIDMVA = data.GetPtrFloat("eleIDMVA");
  //eleIDMVANonTrg = data.GetPtrFloat("eleIDMVANonTrg");
  //eleFiredSingleTrgs = data.GetPtrInt("eleFiredSingleTrgs");
  eleFiredSingleTrgs = data.GetPtrLong64("eleFiredSingleTrgs");
  eleFiredDoubleTrgs = data.GetPtrLong64("eleFiredDoubleTrgs");
  elePFClusEcalIso = data.GetPtrFloat("elePFClusEcalIso");
  elePFClusHcalIso = data.GetPtrFloat("elePFClusHcalIso");
  eleDr03TkSumPt = data.GetPtrFloat("eleDr03TkSumPt");

  nPho = data.GetInt("nPho");
  phoEt = data.GetPtrFloat("phoEt");
  phoEta = data.GetPtrFloat("phoEta");
  phoPhi = data.GetPtrFloat("phoPhi");
  phoSCEta = data.GetPtrFloat("phoSCEta");
  phoSCPhi = data.GetPtrFloat("phoSCPhi");
  phoCalibEt = data.GetPtrFloat("phoCalibEt");
  //phoCalibEt = data.GetPtrFloat("phoEt");
  phoHoverE = data.GetPtrFloat("phoHoverE");
  phoSigmaIEtaIEtaFull5x5 = data.GetPtrFloat("phoSigmaIEtaIEtaFull5x5");
  phoR9 = data.GetPtrFloat("phoR9");
  phoPFChIso = data.GetPtrFloat("phoPFChIso");
  phoPFNeuIso = data.GetPtrFloat("phoPFNeuIso");
  phoPFPhoIso = data.GetPtrFloat("phoPFPhoIso");
  phohasPixelSeed = data.GetPtrInt("phohasPixelSeed");
  phoEleVeto = data.GetPtrInt("phoEleVeto");
  phoIDbit = data.GetPtrShort("phoIDbit");
  phoIDMVA = data.GetPtrFloat("phoIDMVA");

  nMu = data.GetInt("nMu");
  muPt = data.GetPtrFloat("muPt");
  muEta = data.GetPtrFloat("muEta");
  muPhi = data.GetPtrFloat("muPhi");
  muCh = data.GetPtrInt("muCharge");
  muChi2NDF = data.GetPtrFloat("muChi2NDF");
  muMuonHits = data.GetPtrInt("muMuonHits");
  muStations = data.GetPtrInt("muStations");
  muTrkLayers = data.GetPtrInt("muTrkLayers");
  muPixelHits = data.GetPtrInt("muPixelHits");
  muInnerD0 = data.GetPtrFloat("muInnerD0");
  muInnerDz = data.GetPtrFloat("muInnerDz");
  muPFChIso = data.GetPtrFloat("muPFChIso");
  muPFNeuIso = data.GetPtrFloat("muPFNeuIso");
  muPFPhoIso = data.GetPtrFloat("muPFPhoIso");
  muPFPUIso = data.GetPtrFloat("muPFPUIso");
  muCharge = data.GetPtrInt("muCharge");
  muIDbit = data.GetPtrShort("muIDbit");
  muFiredTrgs = data.GetPtrLong64("muFiredTrgs");

  if (data.HasMC()) {
  nMC = data.GetInt("nMC");
  mcPID = data.GetPtrInt("mcPID");
  mcMomPID = data.GetPtrInt("mcMomPID");
  mcPt = data.GetPtrFloat("mcPt");
  mcPhi = data.GetPtrFloat("mcPhi");
  mcEta = data.GetPtrFloat("mcEta");
  mcStatusFlag = (UShort_t*) data.GetPtrShort("mcStatusFlag");
  genWeight = data.GetFloat("genWeight");
}

}


#endif
