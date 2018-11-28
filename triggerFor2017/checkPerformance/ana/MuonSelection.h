#ifndef MuonSelection_h__
#define MuonSelection_h__

Bool_t muMatcher(TreeReader& data, Int_t iMu) {

  if ( !data.HasMC()) return false ;
  Int_t    nMC        = data.GetInt("nMC");
  Int_t*   mcPID      = data.GetPtrInt("mcPID");
  Int_t*   mcMomPID   = data.GetPtrInt("mcMomPID");
  Int_t*   mcGMomPID  = data.GetPtrInt("mcGMomPID");
  Float_t* mcPt       = data.GetPtrFloat("mcPt");


  for (int iMC = 0; iMC < nMC; ++iMC)
    {
      if (mcPt[iMC]< 10.) continue ;
      if (fabs(mcEta[iMC]) > 2.5) continue;
      if ( fabs( mcPID[iMC]) != 13 ) continue ;
      //if (       mcMomPID[iMC]  != 23 && mcMomPID[iMC]  != 24 ) continue ;                                                                                                   
      if ( !( (mcStatusFlag[iMC]>>0&1)==1 || (mcStatusFlag[iMC]>>1&1)==1) ) continue;
      if ( deltaR(mcEta[iMC],mcPhi[iMC],muEta[iMu],muPhi[iMu]) < 0.1 ) return true ;
    }
  return false;
}


#endif
