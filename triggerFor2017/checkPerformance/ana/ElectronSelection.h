#ifndef ElectronSelection_h__
#define ElectronSelection_h__

bool pass_cutbased_80X (int iEle, int iWP) {  // iWP: 0 -loose, 1 - medium, 2- tight

  bool pass = true;

  Float_t sIeIeCut_EB[3]     = {0.011,   0.00998,  0.00998};
  Float_t dEtaInCut_EB[3]    = {0.00477, 0.00311,  0.00308 };
  Float_t dPhiInCut_EB[3]    = {0.222,   0.103,    0.0816};
  Float_t HoECut_EB[3]       = {0.298,   0.253,    0.0414};
  Float_t diffEpCut_EB[3]    = {0.241,   0.134,    0.0129 };
  Int_t ConvVeto_EB[3]       = {1,       1,        1};
  Int_t missHitsCut_EB[3]    = {1,       1,        1};
  Double_t isoCut_EB[3]      = {0.0994,  0.0695,   0.0588};

  Float_t sIeIeCut_EE[3]     = {0.0314,  0.0298,  0.0292};
  Float_t dEtaInCut_EE[3]    = {0.00868, 0.00609, 0.00605};
  Float_t dPhiInCut_EE[3]    = {0.213,   0.045,   0.0394 };
  Float_t HoECut_EE[3]       = {0.101,   0.0878,  0.0641};
  Float_t diffEpCut_EE[3]    = {0.14,    0.13,    0.0129 };
  Int_t ConvVeto_EE[3]       = {1,       1,       1};
  Int_t missHitsCut_EE[3]    = {1,       1,       1};
  Double_t isoCut_EE[3]      = {0.107,   0.0821,  0.0571};

  //effective area
  Double_t EA[7] = {0.1703, 0.1715, 0.1213, 0.1230, 0.1635, 0.1937, 0.2393};
  int index = -1;
  if ( fabs(eleSCEta[iEle])<1.0) index = 0;
  if ( fabs(eleSCEta[iEle])>1.0 && fabs(eleSCEta[iEle])<1.479) index = 1;
  if ( fabs(eleSCEta[iEle])>1.479 && fabs(eleSCEta[iEle])<2.0) index = 2;
  if ( fabs(eleSCEta[iEle])>2.0 && fabs(eleSCEta[iEle])<2.2) index = 3;
  if ( fabs(eleSCEta[iEle])>2.2 && fabs(eleSCEta[iEle])<2.3) index = 4;
  if ( fabs(eleSCEta[iEle])>2.3 && fabs(eleSCEta[iEle])<2.4) index = 5;
  if ( fabs(eleSCEta[iEle])>2.4 && fabs(eleSCEta[iEle])<5.0) index = 6;

  Double_t eleiso = elePFPhoIso[iEle] + elePFNeuIso[iEle]-rho*EA[index];

  if (elePt[iEle] < 5.) pass = false;
  if (fabs(eleSCEta[iEle]) > 2.5) pass = false;
  //if (fabs(eleSCEta[iEle]) > 1.4442 && fabs(eleSCEta[iEle]) < 1.566) pass = false;

  if (fabs(eleSCEta[iEle]) < 1.4442) {
    if (eleSigmaIEtaIEta_Full5x5[iEle] > sIeIeCut_EB[iWP]) pass = false;
    if (fabs(eledEtaseedAtVtx[iEle]) > dEtaInCut_EB[iWP]) pass = false;
    if (fabs(eledPhiAtVtx[iEle]) > dPhiInCut_EB[iWP]) pass = false;
    if (eleHoverE[iEle] > HoECut_EB[iWP]) pass = false;
    if (eleEoverPInv[iEle] > diffEpCut_EB[iWP]) pass = false;
    if (eleConvVeto[iEle] != ConvVeto_EB[iWP]) pass = false;
    if (eleMissHits[iEle] > missHitsCut_EB[iWP]) pass = false;
    if ( (elePFChIso[iEle] + TMath::Max(eleiso, 0.0))/ elePt[iEle] > isoCut_EB[iWP]) pass = false; 

  } else {
    if (eleSigmaIEtaIEta_Full5x5[iEle] > sIeIeCut_EE[iWP]) pass = false;
    if (fabs(eledEtaseedAtVtx[iEle]) > dEtaInCut_EE[iWP]) pass = false;
    if (fabs(eledPhiAtVtx[iEle]) > dPhiInCut_EE[iWP]) pass = false;
    if (eleHoverE[iEle] > HoECut_EE[iWP]) pass = false;
    if (eleEoverPInv[iEle] > diffEpCut_EE[iWP]) pass = false;
    if (eleConvVeto[iEle] != ConvVeto_EE[iWP]) pass = false;                                                                                                   
    if (eleMissHits[iEle] > missHitsCut_EE[iWP]) pass = false;
    if ( (elePFChIso[iEle] + TMath::Max(eleiso, 0.0))/ elePt[iEle] > isoCut_EE[iWP]) pass = false;

  }

  //if (elePFMiniIso[iEle] > 0.1) pass = false;

  return pass;
}

bool select_eleMVA (int iEle) {

  bool pass = true;

  if ( abs(eleSCEta[iEle]) < 0.8 && eleIDMVA[iEle] < 0.972153 ) pass = false;
  else if ( (abs(eleSCEta[iEle]) > 0.8 && abs(eleSCEta[iEle]) < 1.479) && eleIDMVA[iEle] < 0.922126) pass = false;
  else if ( abs(eleSCEta[iEle]) > 1.479 && eleIDMVA[iEle] < 0.610764) pass = false;

  return pass;
}

Bool_t eleMatcher(TreeReader& data, Int_t iEle) {

  if ( !data.HasMC()) return false ;
  Int_t    nMC        = data.GetInt("nMC");
  //Int_t    nEle       = data.GetInt("nEle");
  Int_t*   mcPID      = data.GetPtrInt("mcPID");
  Int_t*   mcMomPID   = data.GetPtrInt("mcMomPID");
  Int_t*   mcGMomPID  = data.GetPtrInt("mcGMomPID");
  Float_t* mcPt       = data.GetPtrFloat("mcPt");
  Float_t* eleEta     = data.GetPtrFloat("eleEta");
  Float_t* elePhi     = data.GetPtrFloat("elePhi");
  Float_t* mcEta      = data.GetPtrFloat("mcEta");
  Float_t* mcPhi      = data.GetPtrFloat("mcPhi");
  UShort_t* mcStatusFlag = (UShort_t*) data.GetPtrShort("mcStatusFlag");

  //return false;

  for (int iMC = 0; iMC < nMC; ++iMC)
    {
      //if (mcPt[iMC]< 10.) continue ;
      if (fabs(mcEta[iMC]) > 2.5) continue;
      if ( fabs( mcPID[iMC]) != 11 ) continue ;
      //if ( mcMomPID[iMC]  != 23 && mcMomPID[iMC]  != 24 ) continue ;
      //if ( !(mcStatusFlag[iMC]>>0&1) && !(mcStatusFlag[iMC]>>1&1) ) continue;
      if ( !(mcStatusFlag[iMC]>>0&1) ) continue;
      if ( deltaR(mcEta[iMC],mcPhi[iMC],eleEta[iEle],elePhi[iEle]) < 0.1 ) return true ;
    }
  return false;

}


Double_t getGenZMass(TreeReader& data)
{
  if ( !data.HasMC()) return Double_t(0) ;
  Int_t    nMC        = data.GetInt("nMC");
  Int_t*   mcPID      = data.GetPtrInt("mcPID");
  Float_t* mcMass     = data.GetPtrFloat("mcMass");
  for (int iMC = 0; iMC < nMC; ++iMC)
    if ( mcPID[iMC] == 23 ) return Double_t(mcMass[iMC]);
  return Double_t(0);
}

Int_t eleGenIndexFinder(TreeReader& data, Int_t iEle)
{
  if ( !data.HasMC()) return Int_t(-1) ;
  Int_t    nMC        = data.GetInt("nMC");
  Int_t    nEle       = data.GetInt("nEle");
  Int_t*   mcPID      = data.GetPtrInt("mcPID");
  Int_t*   mcMomPID   = data.GetPtrInt("mcMomPID");
  Int_t*   mcGMomPID  = data.GetPtrInt("mcGMomPID");
  Float_t* mcPt       = data.GetPtrFloat("mcPt");
  Float_t* eleEta     = data.GetPtrFloat("eleEta");
  Float_t* elePhi     = data.GetPtrFloat("elePhi");
  Float_t* mcEta      = data.GetPtrFloat("mcEta");
  Float_t* mcPhi      = data.GetPtrFloat("mcPhi");
  for (int iMC = 0; iMC < nMC; ++iMC)
    {
      if ( fabs( mcPID[iMC]) != 11 ) continue ;
      if (       mcMomPID[iMC]  != 23 && mcGMomPID[iMC]  != 23 ) continue ;
      if ( deltaR(mcEta[iMC],mcPhi[iMC],eleEta[iEle],elePhi[iEle]) < 0.2 ) return iMC ;
    }
  return Int_t(-1) ;
}

#endif
