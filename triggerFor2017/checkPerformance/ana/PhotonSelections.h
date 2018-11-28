#include <algorithm>

//#include "TMath.h"
#include <TMVA/Reader.h>

bool select_photonXZg( int iPho) {

  bool pass = false;

  if (phoCalibEt[iPho] > 20. && fabs(phoSCEta[iPho]) < 2.5) {
    if (fabs(phoSCEta[iPho]) < 1.479
	//&& phohasPixelSeed[iPho]==0
	//&& phoEleVeto[iPho] ==1
	&& phoSigmaIEtaIEtaFull5x5[iPho] < 0.0102
	&& phoHoverE[iPho] < 0.05
	&& phoPFChIso[iPho] < 2.5)
      pass = true;
    else if (fabs(phoSCEta[iPho]) > 1.57
	     //&& phohasPixelSeed[iPho]==0
	     //&& phoEleVeto[iPho] ==1
	     && phoSigmaIEtaIEtaFull5x5[iPho] < 0.0274
	     && phoHoverE[iPho] < 0.05 
	     && phoPFChIso[iPho] < 2.5 )
      pass = true;
  }

  return pass;

}

bool select_photon_80X (int iPho) {

  bool pass = false;

  int index;
  if (fabs(phoSCEta[iPho]) < 1.0) index = 0;
  else if (fabs(phoSCEta[iPho])>1.0 && fabs(phoSCEta[iPho])<1.479) index = 1;
  else if (fabs(phoSCEta[iPho])>1.479 && fabs(phoSCEta[iPho])<2.0) index = 2;
  else if (fabs(phoSCEta[iPho])>2.0 && fabs(phoSCEta[iPho])<2.2) index = 3;
  else if (fabs(phoSCEta[iPho])>2.2 && fabs(phoSCEta[iPho])<2.3) index = 4;
  else if (fabs(phoSCEta[iPho])>2.3 && fabs(phoSCEta[iPho])<2.4) index = 5;
  else index = 6;

  Double_t EACh[7] = {0.0360, 0.0377, 0.0306, 0.0283, 0.0254, 0.0217, 0.0167};
  Double_t EANeu[7] = {0.0597, 0.0807, 0.0629, 0.0197, 0.0184, 0.0284, 0.0591};
  Double_t EAPho[7] = {0.1210, 0.1107, 0.0699, 0.1056, 0.1457, 0.1719, 0.1998};

  Double_t rhoChIso = TMath::Max(phoPFChIso[iPho] - rho*EACh[index], 0.0);
  Double_t rhoNeuIso = TMath::Max((phoPFNeuIso[iPho] - rho*EANeu[index]), 0.0);
  Double_t rhoPhoIso = TMath::Max(phoPFPhoIso[iPho] - rho*EAPho[index], 0.);

  if (phoCalibEt[iPho] > 20. && fabs(phoSCEta[iPho]) < 2.5) {
    if (fabs(phoSCEta[iPho]) < 1.479
	//&& phohasPixelSeed[iPho]==0
	//&& phoEleVeto[iPho] == 1  
	&& phoHoverE[iPho] < 0.0597
	&& phoSigmaIEtaIEtaFull5x5[iPho] < 0.01031
	&& rhoChIso < 1.295
	&& rhoNeuIso < 10.910+0.0148*phoCalibEt[iPho]+0.000017*pow(phoCalibEt[iPho],2) 
	&& rhoPhoIso < 3.630+0.0047*phoCalibEt[iPho])
      pass = true;
    else if ( fabs(phoSCEta[iPho]) > 1.57
	      //&& phohasPixelSeed[iPho]==0
	      //&& phoEleVeto[iPho] == 1
	      && phoHoverE[iPho] < 0.0481
	      && phoSigmaIEtaIEtaFull5x5[iPho] < 0.03013
	      && rhoChIso < 1.01 
	      && rhoNeuIso < 5.931+0.0163*phoCalibEt[iPho]+0.000014*pow(phoCalibEt[iPho],2)
	      && rhoPhoIso < 6.641+0.0034*phoCalibEt[iPho])
      pass = true;

  }

  return pass;

}



bool phoMatcher (TreeReader &data, int iPho) {

  bool match = false;
  if ( !data.HasMC()) match= false ;

  Int_t nMC        = data.GetInt("nMC");
  Int_t *mcIndex   = data.GetPtrInt("mcIndex");
  Int_t *mcPID     = data.GetPtrInt("mcPID");
  Int_t *mcMomPID  = data.GetPtrInt("mcMomPID");
  float *mcMomMass = data.GetPtrFloat("mcMomMass");
  float *mcEta     = data.GetPtrFloat("mcEta");
  float *mcPhi     = data.GetPtrFloat("mcPhi");
  float *mcMass    = data.GetPtrFloat("mcMass");
  UShort_t* mcStatusFlag = (UShort_t*) data.GetPtrShort("mcStatusFlag");

  for(int imc=0; imc<nMC; imc++){
    if (mcPID[imc]!=22 || ( (mcStatusFlag[imc]>>0&1) != 1 && (mcStatusFlag[imc]>>1&1) != 1) ) continue; //for signal and Zg
    //if (mcPID[imc]!=22 || ( (mcStatusFlag[imc]>>0&1) != 0 || (mcStatusFlag[imc]>>1&1) != 0) ) continue;  //for Z+jets
    //if (abs(mcPID[imc])!=11 || ( (mcStatusFlag[imc]>>0&1) != 1 && (mcStatusFlag[imc]>>1&1) != 1 )) continue; //for WZ and ZZ
    if ( deltaR(mcEta[imc],mcPhi[imc],phoEta[iPho],phoPhi[iPho]) < 0.1 ) match = true ;
  }

  return match;
}
 
