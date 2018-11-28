#ifndef MAKETNPTREE_h
#define MAKETNPTREE_h

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

Int_t     trig_doublePho60;
Int_t     trig_Ele27_WPTight;
Int_t     trig_Ele35_WPTight;
Int_t     trig_Ele17_Ele12;
Int_t     trig_Ele23_Ele12;
Int_t     trig_DoubleEle33;

Int_t     ne;
Int_t     nTag;
Int_t     TagIndex;
Float_t   tag_elePt;
Float_t   tag_eleEta;
Float_t   tag_elePhi;
Float_t   tag_eleSCEta;
Int_t     tag_MatchFilter_Ele35;
Int_t     tag_MatchFilter_Leg1_Ele23;
Int_t     tag_MatchFilter_Leg1_Ele16_Ele12_Ele8;
Int_t     tag_MatchFilter_Leg2_Ele16_Ele12_Ele8;
Int_t     nProbe;
Int_t     ProbeIndex;
Float_t   probe_elePt;
Float_t   probe_eleEta;
Float_t   probe_elePhi;
Float_t   probe_eleSCEta;
Int_t     probe_eleID;
Int_t     probe_MatchFilter_Ele35;
Int_t     probe_MatchFilter_Leg1_Ele23;
Int_t     probe_MatchFilter_Leg2_Ele12;
Int_t     probe_MatchFilter_Leg1_Ele16_Ele12_Ele8;
Int_t     probe_MatchFilter_Leg2_Ele16_Ele12_Ele8;
Int_t     probe_MatchFilter_Leg3_Ele16_Ele12_Ele8;

Float_t   Zm;
Float_t   Zeta;
Float_t   Zy;
Int_t     npair;

Int_t     tag_ele_match_Leg1_DiPho70;
Int_t     tag_ele_match_Leg1_TriPho20_CaloId;
Int_t     tag_ele_match_Leg1_TriPho20_CaloR9Id;
Int_t     tag_ele_match_Leg1_TriPho303010_CaloId;
Int_t     tag_ele_match_Leg1_TriPho303010_CaloR9Id;
Int_t     tag_ele_match_Leg1_TriPho35355_CaloR9Id;
Int_t     nProbe_matchPho;
Float_t   probe_phoPt;
Float_t   probe_phoEta;
Float_t   probe_phoPhi;
Float_t   probe_phoSCEta;
Int_t     probe_MatchFilter_Leg1_DiPho70;
Int_t     probe_MatchFilter_Leg2_DiPho70;
Int_t     probe_MatchFilter_Leg1_TriPho20_CaloId;
Int_t     probe_MatchFilter_Leg2_TriPho20_CaloId;
Int_t     probe_MatchFilter_Leg3_TriPho20_CaloId;
Int_t     probe_MatchFilter_Leg1_TriPho20_CaloR9Id;
Int_t     probe_MatchFilter_Leg2_TriPho20_CaloR9Id;
Int_t     probe_MatchFilter_Leg3_TriPho20_CaloR9Id;
Int_t     probe_MatchFilter_Leg1_TriPho303010_CaloId;
Int_t     probe_MatchFilter_Leg2_TriPho303010_CaloId;
Int_t     probe_MatchFilter_Leg3_TriPho303010_CaloId;
Int_t     probe_MatchFilter_TriPho303010_CaloId;
Int_t     probe_MatchFilter_Leg1_TriPho303010_CaloR9Id;
Int_t     probe_MatchFilter_Leg2_TriPho303010_CaloR9Id;
Int_t     probe_MatchFilter_Leg3_TriPho303010_CaloR9Id;
Int_t     probe_MatchFilter_TriPho303010_CaloR9Id;
Int_t     probe_MatchFilter_Leg1_TriPho35355_CaloR9Id;
Int_t     probe_MatchFilter_Leg2_TriPho35355_CaloR9Id;
Int_t     probe_MatchFilter_Leg3_TriPho35355_CaloR9Id;
Int_t     probe_MatchFilter_TriPho35355_CaloR9Id;

Float_t   Zm_matchPho;
Float_t   Zeta_matchPho;
Float_t   Zy_matchPho;
Int_t     npaireg;


void inittree(TTree* tree) {

  tree->Branch("run",         &run_);
  tree->Branch("lumi",        &lumi_);
  tree->Branch("event",       &event_);
  tree->Branch("nVtx",        &nVtx_);
  tree->Branch("vx",          &vx_);
  tree->Branch("vy",          &vy_);
  tree->Branch("vz",          &vz_);
  tree->Branch("genWeight",   &genWeight_);
  tree->Branch("puweigj",     &puweigj);

  tree->Branch("trig_doublePho60",        &trig_doublePho60);
  tree->Branch("trig_Ele27_WPTight",      &trig_Ele27_WPTight);
  tree->Branch("trig_Ele35_WPTight",      &trig_Ele35_WPTight);
  tree->Branch("trig_Ele17_Ele12",        &trig_Ele17_Ele12);
  tree->Branch("trig_Ele23_Ele12",        &trig_Ele23_Ele12);
  tree->Branch("trig_DoubleEle33",        &trig_DoubleEle33);

  tree->Branch("ne",              &ne);
  tree->Branch("nTag",            &nTag);
  tree->Branch("TagIndex",        &TagIndex);
  tree->Branch("tag_elePt",       &tag_elePt);
  tree->Branch("tag_eleEta",      &tag_eleEta);
  tree->Branch("tag_elePhi",      &tag_elePhi);
  tree->Branch("tag_eleSCEta",    &tag_eleSCEta);
  tree->Branch("tag_MatchFilter_Ele35",  &tag_MatchFilter_Ele35);
  tree->Branch("tag_MatchFilter_Leg1_Ele23",  &tag_MatchFilter_Leg1_Ele23);
  tree->Branch("tag_MatchFilter_Leg1_Ele16_Ele12_Ele8",  &tag_MatchFilter_Leg1_Ele16_Ele12_Ele8);
  tree->Branch("tag_MatchFilter_Leg2_Ele16_Ele12_Ele8",  &tag_MatchFilter_Leg2_Ele16_Ele12_Ele8);

  tree->Branch("nProbe",          &nProbe);
  tree->Branch("ProbeIndex",      &ProbeIndex);
  tree->Branch("probe_elePt",     &probe_elePt);
  tree->Branch("probe_eleEta",    &probe_eleEta);
  tree->Branch("probe_elePhi",    &probe_elePhi);
  tree->Branch("probe_eleSCEta",  &probe_eleSCEta);
  tree->Branch("probe_eleID",     &probe_eleID);
  tree->Branch("probe_MatchFilter_Ele35",  &probe_MatchFilter_Ele35);
  tree->Branch("probe_MatchFilter_Leg1_Ele23",  &probe_MatchFilter_Leg1_Ele23);
  tree->Branch("probe_MatchFilter_Leg2_Ele12",  &probe_MatchFilter_Leg2_Ele12);
  tree->Branch("probe_MatchFilter_Leg1_Ele16_Ele12_Ele8",  &probe_MatchFilter_Leg1_Ele16_Ele12_Ele8);
  tree->Branch("probe_MatchFilter_Leg2_Ele16_Ele12_Ele8",  &probe_MatchFilter_Leg2_Ele16_Ele12_Ele8);
  tree->Branch("probe_MatchFilter_Leg3_Ele16_Ele12_Ele8",  &probe_MatchFilter_Leg3_Ele16_Ele12_Ele8);

  tree->Branch("Zm",           &Zm);
  tree->Branch("Zeta",         &Zeta);
  tree->Branch("Zy",           &Zy);
  tree->Branch("npair",        &npair);

  tree->Branch("tag_ele_match_Leg1_DiPho70",               &tag_ele_match_Leg1_DiPho70);
  tree->Branch("tag_ele_match_Leg1_TriPho20_CaloId",       &tag_ele_match_Leg1_TriPho20_CaloId);
  tree->Branch("tag_ele_match_Leg1_TriPho20_CaloR9Id",     &tag_ele_match_Leg1_TriPho20_CaloR9Id);
  tree->Branch("tag_ele_match_Leg1_TriPho303010_CaloId",   &tag_ele_match_Leg1_TriPho303010_CaloId);
  tree->Branch("tag_ele_match_Leg1_TriPho303010_CaloR9Id", &tag_ele_match_Leg1_TriPho303010_CaloR9Id);
  tree->Branch("tag_ele_match_Leg1_TriPho35355_CaloR9Id",  &tag_ele_match_Leg1_TriPho35355_CaloR9Id);
  tree->Branch("nProbe_matchPho",       &nProbe_matchPho);
  tree->Branch("probe_phoPt",           &probe_phoPt);
  tree->Branch("probe_phoEta",          &probe_phoEta);
  tree->Branch("probe_phoPhi",          &probe_phoPhi);
  tree->Branch("probe_phoSCEta",        &probe_phoSCEta);
  tree->Branch("probe_MatchFilter_Leg1_DiPho70", &probe_MatchFilter_Leg1_DiPho70);
  tree->Branch("probe_MatchFilter_Leg2_DiPho70", &probe_MatchFilter_Leg2_DiPho70);
  tree->Branch("probe_MatchFilter_Leg1_TriPho20_CaloId",  &probe_MatchFilter_Leg1_TriPho20_CaloId);
  tree->Branch("probe_MatchFilter_Leg2_TriPho20_CaloId",  &probe_MatchFilter_Leg2_TriPho20_CaloId);
  tree->Branch("probe_MatchFilter_Leg3_TriPho20_CaloId",  &probe_MatchFilter_Leg3_TriPho20_CaloId);

  tree->Branch("probe_MatchFilter_Leg1_TriPho20_CaloR9Id",  &probe_MatchFilter_Leg1_TriPho20_CaloR9Id);
  tree->Branch("probe_MatchFilter_Leg2_TriPho20_CaloR9Id",  &probe_MatchFilter_Leg2_TriPho20_CaloR9Id);
  tree->Branch("probe_MatchFilter_Leg3_TriPho20_CaloR9Id",  &probe_MatchFilter_Leg3_TriPho20_CaloR9Id);

  tree->Branch("probe_MatchFilter_Leg1_TriPho303010_CaloId",  &probe_MatchFilter_Leg1_TriPho303010_CaloId);
  tree->Branch("probe_MatchFilter_Leg2_TriPho303010_CaloId",  &probe_MatchFilter_Leg2_TriPho303010_CaloId);
  tree->Branch("probe_MatchFilter_Leg3_TriPho303010_CaloId",  &probe_MatchFilter_Leg3_TriPho303010_CaloId);
  tree->Branch("probe_MatchFilter_TriPho303010_CaloId",       &probe_MatchFilter_TriPho303010_CaloId);

  tree->Branch("probe_MatchFilter_Leg1_TriPho303010_CaloR9Id",  &probe_MatchFilter_Leg1_TriPho303010_CaloR9Id);
  tree->Branch("probe_MatchFilter_Leg2_TriPho303010_CaloR9Id",  &probe_MatchFilter_Leg2_TriPho303010_CaloR9Id);
  tree->Branch("probe_MatchFilter_Leg3_TriPho303010_CaloR9Id",  &probe_MatchFilter_Leg3_TriPho303010_CaloR9Id);
  tree->Branch("probe_MatchFilter_TriPho303010_CaloR9Id",       &probe_MatchFilter_TriPho303010_CaloR9Id);

  tree->Branch("probe_MatchFilter_Leg1_TriPho35355_CaloR9Id",  &probe_MatchFilter_Leg1_TriPho35355_CaloR9Id);
  tree->Branch("probe_MatchFilter_Leg2_TriPho35355_CaloR9Id",  &probe_MatchFilter_Leg2_TriPho35355_CaloR9Id);
  tree->Branch("probe_MatchFilter_Leg3_TriPho35355_CaloR9Id",  &probe_MatchFilter_Leg3_TriPho35355_CaloR9Id);
  tree->Branch("probe_MatchFilter_TriPho35355_CaloR9Id",       &probe_MatchFilter_TriPho35355_CaloR9Id);

  tree->Branch("Zm_matchPho",           &Zm_matchPho);
  tree->Branch("Zeta_matchPho",         &Zeta_matchPho);
  tree->Branch("Zy_matchPho",           &Zy_matchPho);
  tree->Branch("npaireg",               &npaireg);

}


//for muon
Int_t     trig_IsoMu27;
Int_t     trig_IsoTkMu24;

Int_t     nmu;
Float_t   tag_muPt;
Float_t   tag_muEta;
Float_t   tag_muPhi;
Float_t   probe_muPt;
Float_t   probe_muEta;
Float_t   probe_muPhi;
Int_t     probe_MatchFilter_Leg1_TriMu;
Int_t     probe_MatchFilter_Leg2_TriMu;
Int_t     probe_MatchFilter_Leg3_TriMu;
Int_t     probe_MatchFilter_DoubleMu207;
Int_t     probe_MatchFilter_Pho23;

void initmutree(TTree* tree) {

  tree->Branch("run",         &run_);
  tree->Branch("lumi",        &lumi_);
  tree->Branch("event",       &event_);
  tree->Branch("nVtx",        &nVtx_);
  tree->Branch("vx",          &vx_);
  tree->Branch("vy",          &vy_);
  tree->Branch("vz",          &vz_);
  tree->Branch("genWeight",   &genWeight_);
  tree->Branch("puweigj",     &puweigj);

  tree->Branch("nmu",             &nmu);
  tree->Branch("nTag",            &nTag);
  //tree->Branch("TagIndex",        &TagIndex);
  tree->Branch("tag_muPt",        &tag_muPt);
  tree->Branch("tag_muEta",       &tag_muEta);
  tree->Branch("tag_muPhi",       &tag_muPhi);
  tree->Branch("probe_muPt",      &probe_muPt);
  tree->Branch("probe_muEta",     &probe_muEta);
  tree->Branch("probe_muPhi",     &probe_muPhi);
  tree->Branch("probe_MatchFilter_Leg1_TriMu",  &probe_MatchFilter_Leg1_TriMu);
  tree->Branch("probe_MatchFilter_Leg2_TriMu",  &probe_MatchFilter_Leg2_TriMu);
  tree->Branch("probe_MatchFilter_Leg3_TriMu",  &probe_MatchFilter_Leg3_TriMu);
  tree->Branch("probe_MatchFilter_DoubleMu207", &probe_MatchFilter_DoubleMu207);

  tree->Branch("probe_phoPt",           &probe_phoPt);
  tree->Branch("probe_phoEta",          &probe_phoEta);
  tree->Branch("probe_phoPhi",          &probe_phoPhi);
  tree->Branch("probe_phoSCEta",        &probe_phoSCEta);
  tree->Branch("probe_MatchFilter_Pho23", &probe_MatchFilter_Pho23);

  tree->Branch("Zm",           &Zm);
  tree->Branch("Zeta",         &Zeta);
  tree->Branch("Zy",           &Zy);
  tree->Branch("npair",        &npair);
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
Float_t genWeight;
Float_t rho;

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
ULong64_t*   eleFiredSingleTrgs;
Long64_t*   eleFiredDoubleTrgs;
Float_t* elePFClusEcalIso;
Float_t* elePFClusHcalIso;
Float_t* eleDr03TkSumPt;

Int_t    nMu;
Float_t* muPt;
Float_t* muEta;
Float_t* muPhi;
Int_t*   muCh;
Short_t* muIDbit;
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
Long64_t*   muFiredTrgs;

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
Long64_t*   phoFiredSingleTrgs;
Long64_t*   phoFiredDoubleTrgs;
Long64_t*   phoFiredTripleTrgs;


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
  rho = data.GetFloat("rho");

  nEle = data.GetInt("nEle");
  eleCharge = data.GetPtrInt("eleCharge");
  elePt = data.GetPtrFloat("eleCalibPt");
  //elePt = data.GetPtrFloat("elePt");
  eleEn = data.GetPtrFloat("eleCalibEn");
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
  eleFiredSingleTrgs = (ULong64_t*) data.GetPtrLong64("eleFiredSingleTrgs");
  eleFiredDoubleTrgs = data.GetPtrLong64("eleFiredDoubleTrgs");
  elePFClusEcalIso = data.GetPtrFloat("elePFClusEcalIso");
  elePFClusHcalIso = data.GetPtrFloat("elePFClusHcalIso");
  eleDr03TkSumPt = data.GetPtrFloat("eleDr03TkSumPt");

  nMu = data.GetInt("nMu");
  muPt = data.GetPtrFloat("muPt");
  muEta = data.GetPtrFloat("muEta");
  muPhi = data.GetPtrFloat("muPhi");
  muCh = data.GetPtrInt("muCharge");
  muIDbit = data.GetPtrShort("muIDbit");
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
  muFiredTrgs = data.GetPtrLong64("muFiredTrgs");

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
  phoFiredSingleTrgs = data.GetPtrLong64("phoFiredSingleTrgs");
  phoFiredDoubleTrgs = data.GetPtrLong64("phoFiredDoubleTrgs");
  phoFiredTripleTrgs = data.GetPtrLong64("phoFiredTripleTrgs");


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
