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
Int_t   isPVGood_;
Int_t   nPU;
Float_t  puTrue;

Float_t   pfMet;
Float_t   pfMetPhi;

Int_t     trig_Ele27_WPTight;
Int_t     trig_Ele35_WPTight;
Int_t     trig_Ele38_WPTight;
Int_t     trig_Ele40_WPTight;
Int_t     trig_DoubleEle33;
Int_t     trig_Ele23_Ele12;
Int_t     trig_Ele23_Ele12_DZ;
Int_t     trig_DiSC30;
Int_t     trig_Ele16_Ele12_Ele8;

Int_t     trig_Mu23_Ele12_DZ;
Int_t     trig_Mu12_Ele23_DZ ;
Int_t     trig_Mu8_DiEle12;
Int_t     trig_Mu8_DiEle12_DZ;
Int_t     trig_DiMu9_Ele9;
Int_t     trig_DiMu9_Ele9_DZ;

Int_t     trig_Mu17;
Int_t     trig_IsoMu24;
Int_t     trig_Mu17_Mu8;
Int_t     trig_Mu17_TkMu8;
Int_t     trig_TkMu17_TkMu8;
Int_t     trig_Mu17_Mu8_DZ;
Int_t     trig_Mu17_TkMu8_DZ;
Int_t     trig_TkMu17_TkMu8_DZ;

Int_t     trig_PFMet120_PFMHT120;

Int_t     ne;
Float_t   eleChIso_[maxPar];
Float_t   elePhoIso_[maxPar];
Float_t   eleNeuIso_[maxPar];
Float_t   eleMiniIso_[maxPar];
Float_t   eleIso_[maxPar];
Float_t   elePt_[maxPar];
Float_t   eleEta_[maxPar];
Float_t   elePhi_[maxPar];
Float_t   eleEoP_[maxPar];
Float_t   eleHoE_[maxPar];
Float_t   eleEoPInv_[maxPar];
Float_t   eleSigIetaIeta_[maxPar];
Float_t   eledEta_[maxPar];
Float_t   eledPhi_[maxPar];
Float_t   eled0_[maxPar];
Float_t   eledz_[maxPar];
Float_t   eleR9_[maxPar];
Float_t   eleSCEt_[maxPar];
Float_t   eleSCEn_[maxPar];
Float_t   eleSCEta_[maxPar];
Float_t   eleSCPhi_[maxPar];
Float_t   eleSCEtaWidth_[maxPar];
Float_t   eleSCPhiWidth_[maxPar];
Int_t     eleConvVeto_[maxPar];
Int_t     eleMissHits_[maxPar];
Int_t     eleCh_[maxPar];
Int_t     eleIdLoose[maxPar];
Int_t     eleIdMedium[maxPar];
Int_t     eleIdTight[maxPar];
Int_t     eleIdHeep[maxPar];
Int_t     eleIdIChep[maxPar];
Int_t     passSieie[maxPar];
Int_t     passdEta[maxPar];
Int_t     passdPhi[maxPar];
Int_t     passHoE[maxPar];
Int_t     passEoPInv[maxPar];
Int_t     passConvVeto[maxPar];
Int_t     passMissHits[maxPar];

Int_t     matchTrigEle27_[maxPar];
Int_t     matchTrigEle35_[maxPar];
Int_t     matchTrigEle38_[maxPar];
Int_t     matchTrigEle40_[maxPar];
Int_t     matchEle23Leg_[maxPar];
Int_t     matchEle12Leg_[maxPar];
Int_t     matchDZ_Ele23_Ele12_[maxPar];
Int_t     matchEle16Leg_TriEle_[maxPar];
Int_t     matchEle12Leg_TriEle_[maxPar];
Int_t     matchEle8Leg_TriEle_[maxPar];

Int_t     nmu;
Float_t   muEta_[maxPar];
Float_t   muPhi_[maxPar];
Float_t   muPt_[maxPar];
Float_t   muChIso_[maxPar];
Float_t   muPhoIso_[maxPar];
Float_t   muNeuIso_[maxPar];
Float_t   muIso_[maxPar];
Int_t     muCh_[maxPar];
Bool_t    muMatch_[maxPar];
Int_t     matchTrigMu17_[maxPar];
Int_t     matchTrigMu8_[maxPar];
Int_t     matchTrigTkMu17_[maxPar];
Int_t     matchTrigTkMu8_[maxPar];
Int_t     matchTrigMu17_Mu8_DZ_[maxPar];
Int_t     matchTrigMu17_TkMu8_DZ_[maxPar];
Int_t     matchTrigTkMu17_TkMu8_DZ_[maxPar];

Int_t     Ztype;
Float_t   Zmass;
Float_t   Zpt;
Float_t   Zeta;
Float_t   Zy;
Float_t   Zphi;
Int_t     Zcharge;
Float_t   dPhidiLep;
Float_t   dRdiLep;

Int_t npho;
Float_t phoEt_[maxPar];
Float_t phoEta_[maxPar];
Float_t phoPhi_[maxPar];
Float_t phoSCEta_[maxPar];
Float_t phoR9_[maxPar];
Float_t phoSigmaIEtaIEta_[maxPar];
Float_t phoHoverE_[maxPar];
Float_t phoChIso_[maxPar];
Float_t dRPhoLep1_[maxPar];
Float_t dRPhoLep2_[maxPar];
Int_t   phoPixSeed_[maxPar];
Int_t   phoEleVeto_[maxPar];

Float_t puweigj;
Float_t puweisj;

Int_t ngen;
Float_t   genMuPt[maxPar];
Float_t   genMuPhi[maxPar];
Float_t   genMuEta[maxPar];
Int_t     genMuCh[maxPar];
Float_t   genElePt[maxPar];
Float_t   genElePhi[maxPar];
Float_t   genEleEta[maxPar];
Int_t     genEleCh[maxPar];
Float_t   genPhoEt;
Float_t   genPhoPhi;
Float_t   genPhoEta;
Float_t   gendRPhoLep1;
Float_t   gendRPhoLep2;
Float_t   genZm;
Float_t   genZpt;
Float_t   genZy;
Int_t     genZch;
Float_t   genXm;
Float_t   genXpt;


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
  tree->Branch("isPVGood",    &isPVGood_);

  tree->Branch("pfMet",        &pfMet);
  tree->Branch("pfMetPhi",     &pfMetPhi);

  //ele trigger
  tree->Branch("trig_Ele27_WPTight",      &trig_Ele27_WPTight);
  tree->Branch("trig_Ele35_WPTight",      &trig_Ele35_WPTight);
  tree->Branch("trig_Ele38_WPTight",      &trig_Ele38_WPTight);
  tree->Branch("trig_Ele40_WPTight",      &trig_Ele40_WPTight);
  tree->Branch("trig_DiSC30",             &trig_DiSC30);

  tree->Branch("trig_Ele23_Ele12",        &trig_Ele23_Ele12);
  tree->Branch("trig_Ele23_Ele12_DZ",     &trig_Ele23_Ele12_DZ);
  tree->Branch("trig_DoubleEle33",        &trig_DoubleEle33);
  tree->Branch("trig_Ele16_Ele12_Ele8",     &trig_Ele16_Ele12_Ele8);

  tree->Branch("trig_Mu12_Ele23_DZ",      &trig_Mu12_Ele23_DZ);
  tree->Branch("trig_Mu23_Ele12_DZ",      &trig_Mu23_Ele12_DZ);
  tree->Branch("trig_Mu8_DiEle12",        &trig_Mu8_DiEle12);
  tree->Branch("trig_Mu8_DiEle12_DZ",     &trig_Mu8_DiEle12_DZ);
  tree->Branch("trig_DiMu9_Ele9",         &trig_DiMu9_Ele9);
  tree->Branch("trig_DiMu9_Ele9_DZ",      &trig_DiMu9_Ele9_DZ);

  //muon trigger
  tree->Branch("trig_Mu17",           &trig_Mu17);
  tree->Branch("trig_IsoMu24",           &trig_IsoMu24);
  tree->Branch("trig_Mu17_Mu8",           &trig_Mu17_Mu8);
  tree->Branch("trig_Mu17_TkMu8",         &trig_Mu17_TkMu8);
  tree->Branch("trig_TkMu17_TkMu8",       &trig_TkMu17_TkMu8);
  tree->Branch("trig_Mu17_Mu8_DZ",        &trig_Mu17_Mu8_DZ);
  tree->Branch("trig_Mu17_TkMu8_DZ",      &trig_Mu17_TkMu8_DZ);
  tree->Branch("trig_TkMu17_TkMu8_DZ",    &trig_TkMu17_TkMu8_DZ);

  //ele-mu trigger
  tree->Branch("trig_Mu12_Ele23_DZ",      &trig_Mu12_Ele23_DZ);
  tree->Branch("trig_Mu23_Ele12_DZ",      &trig_Mu23_Ele12_DZ);

  //met
  tree->Branch("trig_PFMet120_PFMHT120", &trig_PFMet120_PFMHT120);

  tree->Branch("ne",          &ne);
  tree->Branch("eleChIso",    &eleChIso_,   "eleChIso[ne]/F");
  tree->Branch("elePhoIso",   &elePhoIso_,  "elePhoIso[ne]/F");
  tree->Branch("eleNeuIso",   &eleNeuIso_,  "eleNeuIso[ne]/F");
  tree->Branch("eleMiniIso",  &eleMiniIso_, "eleMiniIso[ne]/F");
  tree->Branch("eleIso",      &eleIso_,     "eleIso[ne]/F");
  tree->Branch("elePt",       &elePt_,      "elePt[ne]/F");
  tree->Branch("eleEta",      &eleEta_,     "eleEta[ne]/F");
  tree->Branch("elePhi",      &elePhi_,     "elePhi[ne]/F");
  tree->Branch("eleEoP",      &eleEoP_,     "eleEoP[ne]/F");
  tree->Branch("eleEoPInv",   &eleEoPInv_,  "eleEoPInv[ne]/F");
  tree->Branch("eleHoE",      &eleHoE_,     "eleHoE[ne]/F");
  tree->Branch("eled0",       &eled0_,      "eled0[ne]/F");
  tree->Branch("eledz",       &eledz_,      "eledz[ne]/F");
  tree->Branch("eleR9",       &eleR9_,      "eleR9[ne]/F");
  tree->Branch("eledEta",     &eledEta_,    "eledEta[ne]/F");
  tree->Branch("eledPhi",     &eledPhi_,    "eledPhi[ne]/F");
  tree->Branch("eleSigIetaIeta",  &eleSigIetaIeta_,  "eleSigIetaIeta[ne]/F");
  tree->Branch("eleSCEt",     &eleSCEt_,    "eleSCEt[ne]/F");
  tree->Branch("eleSCEn",     &eleSCEn_,   "eleSCEn_[ne]/F");
  tree->Branch("eleSCEta",    &eleSCEta_,  "eleSCEta[ne]/F" );
  tree->Branch("eleSCPhi",    &eleSCPhi_,  "eleSCPhi[ne]/F" );
  tree->Branch("eleSCEtaWidth",    &eleSCEtaWidth_,  "eleSCEtaWidth[ne]/F" );
  tree->Branch("eleSCPhiWidth",    &eleSCPhiWidth_,  "eleSCPhiWidth[ne]/F" );
  tree->Branch("eleConvVeto",      &eleConvVeto_,    "eleConvVeto[ne]/F");
  tree->Branch("eleMissHits",      &eleMissHits_,    "eleMissHits[ne]/F");
  tree->Branch("eleCh",            &eleCh_,          "eleCh[ne]/I");
  tree->Branch("eleIdLoose",       &eleIdLoose,      "eleIdLoose[ne]/I");
  tree->Branch("eleIdMedium",      &eleIdMedium,     "eleIdMedium[ne]/I");
  tree->Branch("eleIdTight",       &eleIdTight,      "eleIdTight[ne]/I");
  tree->Branch("eleIdHeep",        &eleIdHeep,       "eleIdHeep[ne]/I");
  tree->Branch("eleIdIChep",       &eleIdIChep,      "eleIdIChep[ne]/I");
  tree->Branch("passSieie",        &passSieie,       "passSieie[ne]/I");
  tree->Branch("passdEta",         &passdEta,        "passdEta[ne]/I");
  tree->Branch("passdPhi",         &passdPhi,        "passdPhi[ne]/I");
  tree->Branch("passHoE",          &passHoE,         "passHoE[ne]/I");
  tree->Branch("passEoPInv",       &passEoPInv,      "passEoPInv[ne]/I");
  tree->Branch("passConvVeto",     &passConvVeto,    "passConvVeto[ne]/I");
  tree->Branch("passMissHits",     &passMissHits,    "passMissHits[ne]/I");

  //ele trigger matching to filter
  tree->Branch("matchTrigEle27",        &matchTrigEle27_,      "matchTrigEle27[ne]/I");
  tree->Branch("matchTrigEle35",        &matchTrigEle35_,      "matchTrigEle35[ne]/I");
  tree->Branch("matchTrigEle38",        &matchTrigEle38_,      "matchTrigEle38[ne]/I");
  tree->Branch("matchTrigEle40",        &matchTrigEle40_,      "matchTrigEle40[ne]/I");
  tree->Branch("matchEle23Leg",         &matchEle23Leg_,       "matchEle23Leg[ne]/I");
  tree->Branch("matchEle12Leg",         &matchEle12Leg_,       "matchEle12Leg[ne]/I");
  tree->Branch("matchDZ_Ele23_Ele12_",  &matchDZ_Ele23_Ele12_, "matchDZ_Ele23_Ele12[ne]/I");
  tree->Branch("matchEle16Leg_TriEle",  &matchEle16Leg_TriEle_,    "matchEle16Leg_TriEle[ne]/I");
  tree->Branch("matchEle12Leg_TriEle",  &matchEle12Leg_TriEle_,    "matchEle12Leg_TriEle[ne]/I");
  tree->Branch("matchEle8Leg_TriEle",   &matchEle8Leg_TriEle_,     "matchEle8Leg_TriEle[ne]/I");

  tree->Branch("nmu",         &nmu);
  tree->Branch("muPt",        &muPt_,      "muPt[nmu]/F");
  tree->Branch("muEta",       &muEta_,     "muEta[nmu]/F");
  tree->Branch("muPhi",       &muPhi_,     "muPhi[nmu]/F");
  tree->Branch("muMatch",     &muMatch_,   "muMatch[nmu]/O");
  tree->Branch("muChIso",     &muChIso_,   "muChIso[nmu]/F");
  tree->Branch("muPhoIso",    &muPhoIso_,  "muPhoIso[nmu]/F");
  tree->Branch("muNeuIso",    &muNeuIso_,  "muNeuIso[nmu]/F");
  tree->Branch("muIso",       &muIso_,     "muIso[nmu]/F");
  tree->Branch("muCh",        &muCh_,      "muCh[nmu]/I");
  tree->Branch("matchTrigMu17",               &matchTrigMu17_,          "matchTrigMu17[nmu]/I");
  tree->Branch("matchTrigMu17",               &matchTrigMu17_,          "matchTrigMu8[nmu]/I");
  tree->Branch("matchTrigTkMu17",             &matchTrigTkMu17_,        "matchTrigTkMu17[nmu]/I");
  tree->Branch("matchTrigTkMu8",              &matchTrigTkMu8_,         "matchTrigTkMu8[nmu]/I");
  tree->Branch("matchTrigMu17_Mu8_DZ",        &matchTrigMu17_Mu8_DZ_,   "matchTrigMu17_Mu8_DZ[nmu]I");
  tree->Branch("matchTrigMu17_TkMu8_DZ",      &matchTrigMu17_TkMu8_DZ_, "matchTrigMu17_TkMu8_DZ[nmu]I");
  tree->Branch("matchTrigTkMu17_TkMu8_DZ",    &matchTrigTkMu17_TkMu8_DZ_, "matchTrigTkMu17_TkMu8_DZ[nmu]I");

  tree->Branch("Ztype",       &Ztype);
  tree->Branch("Zmass",       &Zmass);
  tree->Branch("Zpt",         &Zpt);
  tree->Branch("Zeta",        &Zeta);
  tree->Branch("Zphi",        &Zphi);
  tree->Branch("Zy",          &Zy);
  tree->Branch("Zcharge",     &Zcharge);
  tree->Branch("dPhidiLep",   &dPhidiLep);
  tree->Branch("dRdiLep",     &dRdiLep);

  tree->Branch("ngen",        &ngen);
  tree->Branch("genMuPt",     &genMuPt,    "genMuPt[ngen]/F");
  tree->Branch("genMuPhi",    &genMuPhi,   "genMuPhi[ngen]/F");
  tree->Branch("genMuEta",    &genMuEta,   "genMuEta[ngen]/F");
  tree->Branch("genMuCh",     &genMuCh,    "genMuCh[ngen]/I");
  tree->Branch("genElePt",    &genElePt,   "genElePt[ngen]/F");
  tree->Branch("genElePhi",   &genElePhi,  "genElePhi[ngen]/F");
  tree->Branch("genEleEta",   &genEleEta,  "genEleEta[ngen]/F");
  tree->Branch("genEleCh",    &genEleCh,   "genEleCh[ngen]/I");
  tree->Branch("genPhoEt",    &genPhoEt,   "genPhoEt/F");
  tree->Branch("genPhoPhi",   &genPhoPhi,  "genPhoPhi/F");
  tree->Branch("genPhoEta",   &genPhoEta,  "genPhoEta/F");
  tree->Branch("gendRPhoLep1",&gendRPhoLep1);
  tree->Branch("gendRPhoLep2",&gendRPhoLep2);
  tree->Branch("genZm",       &genZm);
  tree->Branch("genZpt",      &genZpt);
  tree->Branch("genZy",       &genZy);
  tree->Branch("genZch",      &genZch);
  tree->Branch("genXm",       &genXm);
  tree->Branch("genXpt",      &genXpt);
}

Int_t     run;
Long64_t  event;
Int_t     lumis;
Bool_t    isData;
ULong64_t hlt;
ULong64_t hltPho;
ULong64_t hltJet;
Int_t nVtx;
Float_t vz;
Float_t vx; 
Float_t vy;
Float_t genWeight;
Float_t rho;
Int_t   isPVGood;
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
Float_t* eleR9;
Float_t* eleSCEn ;
Float_t* elePFChIso ;
Float_t* elePFPhoIso ;
Float_t* elePFNeuIso ;
Float_t* elePFPUIso ;
Float_t* elePFMiniIso ;
Short_t* eleID;
Float_t* eleIDMVA;
Float_t* eleIDMVAHZZ;
Float_t* elePFClusEcalIso;
Float_t* elePFClusHcalIso;
Float_t* eleDr03TkSumPt;
Int_t*   eleFiredSingleTrgs;
Int_t*   eleFiredDoubleTrgs;


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
Int_t* muFiredTrgs;

Int_t nPho;
Float_t* phoEt;
Float_t* phoSCEta;
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
  hltJet = data.GetLong64("HLTJet");
  nVtx = data.GetInt("nVtx");
  vx = data.GetFloat("vtx");
  vy = data.GetFloat("vty");
  vz = data.GetFloat("vtz");
  rho = data.GetFloat("rho");
  isPVGood = data.GetBool("isPVGood");
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
  eleR9 = data.GetPtrFloat("eleR9");
  eleSCEn = data.GetPtrFloat("eleSCEn");
  elePFChIso = data.GetPtrFloat("elePFChIso");
  elePFPhoIso = data.GetPtrFloat("elePFPhoIso");
  elePFNeuIso = data.GetPtrFloat("elePFNeuIso");
  elePFMiniIso = data.GetPtrFloat("elePFMiniIso");
  eleID = data.GetPtrShort("eleIDbit");
  eleIDMVA = data.GetPtrFloat("eleIDMVA");
  eleIDMVAHZZ = data.GetPtrFloat("eleIDMVAHZZ");
  //eleIDMVANonTrg = data.GetPtrFloat("eleIDMVANonTrg");
  elePFClusEcalIso = data.GetPtrFloat("elePFClusEcalIso");
  elePFClusHcalIso = data.GetPtrFloat("elePFClusHcalIso");
  eleDr03TkSumPt = data.GetPtrFloat("eleDr03TkSumPt");
  eleFiredSingleTrgs = data.GetPtrInt("eleFiredSingleTrgs");
  eleFiredDoubleTrgs = data.GetPtrInt("eleFiredDoubleTrgs");

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
  muFiredTrgs = data.GetPtrInt("muFiredTrgs");

  nPho = data.GetInt("nPho");
  phoEt = data.GetPtrFloat("phoEt");
  phoEta = data.GetPtrFloat("phoEta");
  phoSCEta = data.GetPtrFloat("phoSCEta");
  phoPhi = data.GetPtrFloat("phoPhi");
  //phoCalibEt = data.GetPtrFloat("phoCalibEt");
  phoCalibEt = data.GetPtrFloat("phoEt");
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
