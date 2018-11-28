/**********************************
TnP method for tri-leptons( 3e or 3mu) and tri-photons trigger
for each leg by matching filter using cut and cout number of Z
with lepton/photon passing selection
**********************************/

#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TF1.h>
#include <TMath.h>

#include "untuplizer.h"
#include "makeTnPtree.h"
#include "ElectronSelection.h"
#include "MuonSelection.h"
#include "PhotonSelections.h"
#include "puweicalc.h"


void makeTnPtree(TString inputfile = "HiForest.root", TString outputfile = "XZg.root", bool mc = true) {

  outputfile = "minitrees/";
  if (inputfile.Contains("SingleEle") ) outputfile += "SingleEle";
  if (inputfile.Contains("SingleMuon") ) outputfile += "SingleMu";
  if (mc)  outputfile = outputfile + "DYJets_amcnlo_summer16";
  //else outputfile = outputfile + "data";

  if (inputfile.Contains("_Run2016B") ) outputfile += "_Run2016B";
  if (inputfile.Contains("_Run2016C") ) outputfile += "_Run2016C";
  if (inputfile.Contains("_Run2016D") ) outputfile += "_Run2016D";
  if (inputfile.Contains("_Run2016E") ) outputfile += "_Run2016E";
  if (inputfile.Contains("_Run2016F_SepRereco1") ) outputfile += "_Run2016F_SepRereco1";
  if (inputfile.Contains("_Run2016F_SepRereco2") ) outputfile += "_Run2016F_SepRereco2";
  if (inputfile.Contains("_Run2016G") ) outputfile += "_Run2016G";
  if (inputfile.Contains("run2016H/PromptReco_v2") ) outputfile += "_Run2016H_PRv2";
  if (inputfile.Contains("run2016H/PromptReco_v3") ) outputfile += "_Run2016H_PRv3";

  if (inputfile.Contains("2017A_PromptReco_v2") ) outputfile += "_Run2017A_PromptReco_v2";
  if (inputfile.Contains("2017A_PromptReco_v3") ) outputfile += "_Run2017A_PromptReco_v3";
  if (inputfile.Contains("2017B_PromptReco_v1") ) outputfile += "_Run2017B_PromptReco_v1";
  if (inputfile.Contains("2017B_PromptReco_v2") ) outputfile += "_Run2017B_PromptReco_v2";
  if (inputfile.Contains("2017C_PromptReco_v1") ) outputfile += "_Run2017C_PromptReco_v1";
  if (inputfile.Contains("2017C_PromptReco_v2") ) outputfile += "_Run2017C_PromptReco_v2";
  if (inputfile.Contains("2017C_PromptReco_v3") ) outputfile += "_Run2017C_PromptReco_v3";
  if (inputfile.Contains("2017D_PromptReco_v1") ) outputfile += "_Run2017D_PromptReco_v1";
  if (inputfile.Contains("2017E_PromptReco_v1") ) outputfile += "_Run2017E_PromptReco_v1";
  if (inputfile.Contains("2017F_PromptReco_v1") ) outputfile += "_Run2017F_PromptReco_v1";


  outputfile += "_TnP_HEEPID.root";

  TreeReader data1(inputfile, "ggNtuplizer/EventTree");

  TFile *fo = TFile::Open(outputfile, "RECREATE");

  TTree *eletree = new TTree("eletree", "output tree");
  TTree *photree = new TTree("photree", "output tree for pho");
  TTree *mutree = new TTree("mutree", "output tree");
  TTree *muphotree = new TTree("muphotree", "output tree of dimu and mupho");

  TTree *dieletree = new TTree("dieletree", "output tree of di-ele trigger");
  TTree *trieletree = new TTree("trieletree", "output tree of tri-ele trigger");

  TTree *pho20Calotree = new TTree("pho20Calotree", "output tree for tripho20 Calo");


  inittree(eletree);
  inittree(photree);
  initmutree(mutree);
  initmutree(muphotree);

  inittree(dieletree);
  inittree(trieletree);
  inittree(pho20Calotree);

  TH1D *hntotweight = new TH1D("hntotweight", "", 2, 1, 3);

  Long64_t ev1 = 0;
  Long64_t ev2 = -1;

  if (ev2 < 0) ev2 = data1.GetEntriesFast();
  if (ev2 > data1.GetEntriesFast()) ev2 = data1.GetEntriesFast();

  PUWeightCalculator puCalcGJ;

  if ( data1.HasMC() ) {
    puCalcGJ.Init("external/PU_histo_13TeV_GoldenJSON_summer16_69000nb.root");
  }

  double ntot = 0;

  //for (Long64_t ev = ev1; ev < ev2; ++ev) {
  for (Long64_t ev = ev1; ev < 20000000; ++ev) {
    if (ev % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data1.GetEntriesFast());

    data1.GetEntry(ev);
    readggtree(data1);

    vz_ = vz;
    vx_ = vx;
    vy_ = vy;
    nVtx_ = nVtx;
    run_ = run;
    event_ = event;
    lumi_ = lumis;

    trig_Ele27_WPTight = hlt>>4&1;
    trig_Ele35_WPTight = hlt>>3&1;
    trig_Ele23_Ele12 = hlt>>5&1;
    trig_DoubleEle33 = hlt>>11&1;

    trig_IsoMu27 = hlt>>19&1;
    trig_IsoTkMu24 = hlt>>20&1;

    trig_doublePho60 = hltPho >> 22&1;

    //PU reweighting for MC
    if (data1.HasMC()) {
      float* puTrue = data1.GetPtrFloat("puTrue");
      puweigj = (float) puCalcGJ.GetWeight(run, puTrue[1]); // in-time PU
      genWeight_ = (genWeight > 0) ? 1. : -1.;
      ntot += genWeight_;
    }


    //if (nEle < 2) continue;

    vector<int> acc_tag_ele, acc_probe_ele;
    vector<int> acc_tag_ele_leg2_diele, acc_tag_ele_leg2_triele;
    vector<int> acc_tag_ele_leg2_tripho20Calo, acc_tag_ele_leg2_tripho20CaloR9;
    vector<int> acc_tag_ele_leg2_tripho30Calo, acc_tag_ele_leg2_tripho30CaloR9;
    vector<int> acc_tag_ele_leg2_tripho35CaloR9;

    for (int i=0; i<nEle; i++) {
      if (elePt[i] < 10.) continue;
      if (fabs(eleSCEta[i]) > 2.5) continue;
      if (fabs(eleSCEta[i]) > 1.4442 && fabs(eleSCEta[i]) < 1.566) continue;
      if (mc) {
        bool match = eleMatcher(data1, i);
        if (!match) continue;
      }

      //tight id for electron
      //if ( (eleID[i] >>3&1)==1 && (eleFiredSingleTrgs[i] >>12&1)==1) acc_tag_ele.push_back(i); // tag --WP_Ele27
      if ( (eleID[i] >>3&1)==1 && (eleFiredSingleTrgs[i] >>39&1)==1) acc_tag_ele.push_back(i); // tag --WP_Ele35
      if ( (eleID[i] >>3&1)==1 && (eleFiredSingleTrgs[i] >>39&1)==1 && (eleFiredDoubleTrgs[i] >>2&1)==1) acc_tag_ele_leg2_diele.push_back(i); 
      if ( (eleID[i] >>3&1)==1 && (eleFiredSingleTrgs[i] >>12&1)==1 && (eleFiredDoubleTrgs[i] >>9&1)==1) acc_tag_ele_leg2_triele.push_back(i);
      //if ( (eleID[i] >>3&1)==1 && (eleFiredDoubleTrgs[i] >>8&1)==1) acc_tag_ele_leg2_triele.push_back(i);

      //tag ele matching to first filter of triple photon trigger
      //if ( (eleID[i] >>3&1)==1 && (eleFiredSingleTrgs[i] >>12&1)==1 && (phoFiredTripleTrgs[i]>>0&1)==1) acc_tag_ele_leg2_tripho20Calo.push_back(i);
      //if ( (eleID[i] >>3&1)==1 && (eleFiredSingleTrgs[i] >>12&1)==1 && (phoFiredTripleTrgs[i]>>2&1)==1) acc_tag_ele_leg2_tripho20CaloR9.push_back(i);

      //probe ele
      if ( (eleID[i] >>4&1)==1) acc_probe_ele.push_back(i); //probe
      //if (eleIDMVA[i] > 0.2) acc_probe_ele.push_back(i);

    }

    //cout << "number of tag = " << acc_tag_ele.size() << "\t number of probe = " << acc_probe_ele.size() << endl;
    if (acc_tag_ele.size() < 1) continue;
    if (acc_probe_ele.size() < 1) continue;

    int nZ = 0, nZ_matchPho = 0;
    ne = 0;
    npair = 0;
    npaireg = 0;

    //loop on tag ele matching to singleEle trigger
    for (vector<int>::iterator itag = acc_tag_ele.begin() ; itag != acc_tag_ele.end(); ++itag) {
      //initial values
      tag_elePt = -99.;
      tag_eleEta = -99.;
      tag_eleSCEta = -99.;
      tag_elePhi = -99.;
      tag_MatchFilter_Ele35= 0;
      tag_MatchFilter_Leg1_Ele16_Ele12_Ele8 = 0;
      tag_MatchFilter_Leg2_Ele16_Ele12_Ele8 = 0;

      ProbeIndex = -99.;
      probe_elePt = -99.;
      probe_eleEta = -99.;
      probe_eleSCEta= -99.;
      probe_elePhi = -99.;
      probe_MatchFilter_Ele35= 0;
      probe_MatchFilter_Leg1_Ele23 = 0;
      probe_MatchFilter_Leg2_Ele12 = 0;
      probe_MatchFilter_Leg1_Ele16_Ele12_Ele8 = 0;
      probe_MatchFilter_Leg2_Ele16_Ele12_Ele8 = 0;
      probe_MatchFilter_Leg3_Ele16_Ele12_Ele8 = 0;

      Zm = -99.;
      Zeta = -99.;
      //Zpt = -99.;

      if ( elePt[*itag] < 30.) continue;
      nTag++;
      TagIndex = *itag;

      //cout << "passing tag selection" << endl;
      TLorentzVector vTag ;
      vTag.SetPtEtaPhiM(elePt[*itag], eleEta[*itag], elePhi[*itag], 0.511*0.001);

      vector<float> Zmass, ZEta, ZPt;
      vector<int> index;
      nProbe = 0;

      for (unsigned int iele = 0; iele < acc_probe_ele.size(); iele++) {
	if ( *itag >= acc_probe_ele[iele]) continue; //leading electron is tag
	if ( elePt[acc_probe_ele[iele]] < 10.) continue;

	TLorentzVector vProbe ,vZ ;
	vProbe.SetPtEtaPhiM(elePt[acc_probe_ele[iele]], eleEta[acc_probe_ele[iele]], elePhi[acc_probe_ele[iele]],0.511*0.001);
	vZ = vTag + vProbe;
	if ( vZ.M() < 60 || vZ.M() > 120 ) continue ;

	Zmass.push_back(vZ.M());
        ZEta.push_back(vZ.Eta());
        ZPt.push_back(vZ.Pt());
        index.push_back(acc_probe_ele[iele]);
      }

      nZ = Zmass.size();
      if (nZ < 1) continue;
      float minZ = 999.;
      int iProbe = -1;
      for (int i = 0; i < nZ; i++) {
        float deltaM = fabs(91.19-Zmass[i]);
        if (minZ > deltaM) {
          minZ = deltaM;
          iProbe = index[i];
          Zm = Zmass[i];
          Zeta = ZEta[i];
          //Zpt = ZPt[i];
        }
      }

      if (iProbe < 0) continue;
      if ( *itag == iProbe) continue;

      tag_elePt = elePt[*itag];
      tag_eleEta = eleEta[*itag];
      tag_eleSCEta= eleSCEta[*itag];
      tag_elePhi = elePhi[*itag];
      tag_MatchFilter_Leg1_Ele23 = eleFiredDoubleTrgs[*itag] >> 2&1;
      tag_MatchFilter_Leg1_Ele16_Ele12_Ele8 = eleFiredDoubleTrgs[*itag] >> 9&1;
      tag_MatchFilter_Leg2_Ele16_Ele12_Ele8 = eleFiredDoubleTrgs[*itag] >> 10&1;
	
      ProbeIndex = iProbe;
      probe_elePt = elePt[iProbe];
      probe_eleEta = eleEta[iProbe];
      probe_eleSCEta= eleSCEta[iProbe];
      probe_elePhi = elePhi[iProbe];
      probe_eleID = eleID[iProbe]>>3&1;
      probe_MatchFilter_Ele35 = eleFiredSingleTrgs[iProbe] >> 39&1;
      probe_MatchFilter_Leg1_Ele23 = eleFiredDoubleTrgs[iProbe] >> 2&1;
      probe_MatchFilter_Leg2_Ele12 = eleFiredDoubleTrgs[iProbe] >> 3&1;
      probe_MatchFilter_Leg1_Ele16_Ele12_Ele8 = eleFiredDoubleTrgs[iProbe] >> 9&1;
      probe_MatchFilter_Leg2_Ele16_Ele12_Ele8 = eleFiredDoubleTrgs[iProbe] >> 10&1;
      probe_MatchFilter_Leg3_Ele16_Ele12_Ele8 = eleFiredDoubleTrgs[iProbe] >> 11&1;

      //probe_ele_hltDiEG30CaloIdLV2EtUnseededFilter = ;
      //probe_ele_hltDiEG30CaloIdLV2R9IdVLEtUnseededFilter = ;
      //probe_ele_hltDiEG35CaloIdLV2R9IdVLEtUnseededFilter = ;

      npair++;
      //cout << "filling tnp tree for e-e" << endl;
      eletree->Fill();
    }

    
    /*---------loop on tag ele matching first leg filter of di-ele trigger---------*/
    for (vector<int>::iterator itag = acc_tag_ele_leg2_diele.begin(); itag != acc_tag_ele_leg2_diele.end(); ++itag) {
      //initial values
      tag_elePt = -99.;
      tag_eleEta = -99.;
      tag_eleSCEta = -99.;
      tag_elePhi = -99.;
      //ele1_MatchTrig27 = -99;
      tag_MatchFilter_Ele35 = 0;
      tag_MatchFilter_Leg1_Ele16_Ele12_Ele8 = 0;
      tag_MatchFilter_Leg2_Ele16_Ele12_Ele8 = 0;

      ProbeIndex = -99.;
      probe_elePt = -99.;
      probe_eleEta = -99.;
      probe_eleSCEta= -99.;
      probe_elePhi = -99.;
      probe_MatchFilter_Ele35 = 0;
      probe_MatchFilter_Leg1_Ele23 = 0;
      probe_MatchFilter_Leg2_Ele12 = 0;
      probe_MatchFilter_Leg1_Ele16_Ele12_Ele8 = 0;
      probe_MatchFilter_Leg2_Ele16_Ele12_Ele8 = 0;
      probe_MatchFilter_Leg3_Ele16_Ele12_Ele8 = 0;

      Zm = -99.;
      Zeta = -99.;

      if ( elePt[*itag] < 35.) continue;
      nTag++;

      //cout << "passing tag selection" << endl;
      TLorentzVector vTag ;
      vTag.SetPtEtaPhiM(elePt[*itag], eleEta[*itag], elePhi[*itag], 0.511*0.001);

      vector<float> Zmass, ZEta, ZPt;
      vector<int> index;
      nProbe = 0;

      for (vector<int>::iterator iele = acc_probe_ele.begin(); iele != acc_probe_ele.end(); ++iele) {
	//if ( *itag >= acc_probe_ele[iele]) continue; //leading electron is tag
	if (*itag >= *iele) continue;
	if ( elePt[*iele] < 10.) continue;

	TLorentzVector vProbe ,vZ ;
	vProbe.SetPtEtaPhiM(elePt[*iele], eleEta[*iele], elePhi[*iele],0.511*0.001);
	vZ = vTag + vProbe;
	if ( vZ.M() < 60 || vZ.M() > 120 ) continue ;

	Zmass.push_back(vZ.M());
        ZEta.push_back(vZ.Eta());
        ZPt.push_back(vZ.Pt());
        index.push_back(*iele);
      }

      nZ = Zmass.size();
      if (nZ < 1) continue;
      float minZ = 999.;
      int iProbe = -1;
      for (int i = 0; i < nZ; i++) {
        float deltaM = fabs(91.19-Zmass[i]);
        if (minZ > deltaM) {
          minZ = deltaM;
          iProbe = index[i];
          Zm = Zmass[i];
          Zeta = ZEta[i];
          //Zpt = ZPt[i];
        }
      }

      if (iProbe < 0) continue;
      //if ( *itag == iProbe) continue;

      tag_elePt = elePt[*itag];
      tag_eleEta = eleEta[*itag];
      tag_eleSCEta= eleSCEta[*itag];
      tag_elePhi = elePhi[*itag];
      tag_MatchFilter_Ele35 = eleFiredSingleTrgs[*itag] >> 39&1;
      tag_MatchFilter_Leg1_Ele23 = eleFiredDoubleTrgs[*itag] >> 2&1;
      tag_MatchFilter_Leg1_Ele16_Ele12_Ele8 = eleFiredDoubleTrgs[*itag] >> 9&1;
      tag_MatchFilter_Leg2_Ele16_Ele12_Ele8 = eleFiredDoubleTrgs[*itag] >> 10&1;
	
      ProbeIndex = iProbe;
      probe_elePt = elePt[iProbe];
      probe_eleEta = eleEta[iProbe];
      probe_eleSCEta= eleSCEta[iProbe];
      probe_elePhi = elePhi[iProbe];
      probe_MatchFilter_Ele35 = eleFiredSingleTrgs[iProbe] >> 39&1;
      probe_MatchFilter_Leg1_Ele23 = eleFiredDoubleTrgs[iProbe] >> 2&1;
      probe_MatchFilter_Leg2_Ele12 = eleFiredDoubleTrgs[iProbe] >> 3&1;
      probe_MatchFilter_Leg1_Ele16_Ele12_Ele8 = eleFiredDoubleTrgs[iProbe] >> 9&1;
      probe_MatchFilter_Leg2_Ele16_Ele12_Ele8 = eleFiredDoubleTrgs[iProbe] >> 10&1;
      probe_MatchFilter_Leg3_Ele16_Ele12_Ele8 = eleFiredDoubleTrgs[iProbe] >> 11&1;

      npair++;
      dieletree->Fill();

    } //end of di-ele tree
    

    /*---------loop on tag ele matching first leg filter of tri-ele trigger---------*/
    for (vector<int>::iterator itag = acc_tag_ele_leg2_triele.begin(); itag != acc_tag_ele_leg2_triele.end(); ++itag) {
      //initial values
      tag_elePt = -99.;
      tag_eleEta = -99.;
      tag_eleSCEta = -99.;
      tag_elePhi = -99.;

      ProbeIndex = -99.;
      probe_elePt = -99.;
      probe_eleEta = -99.;
      probe_eleSCEta= -99.;
      probe_elePhi = -99.;
      probe_MatchFilter_Leg1_Ele23 = 0;
      probe_MatchFilter_Leg2_Ele12 = 0;
      probe_MatchFilter_Leg1_Ele16_Ele12_Ele8 = 0;
      probe_MatchFilter_Leg2_Ele16_Ele12_Ele8 = 0;
      probe_MatchFilter_Leg3_Ele16_Ele12_Ele8 = 0;

      Zm = -99.;
      Zeta = -99.;

      if ( elePt[*itag] < 30.) continue;
      nTag++;

      //cout << "passing tag selection" << endl;
      TLorentzVector vTag ;
      vTag.SetPtEtaPhiM(elePt[*itag], eleEta[*itag], elePhi[*itag], 0.511*0.001);

      vector<float> Zmass, ZEta, ZPt;
      vector<int> index;
      nProbe = 0;

      for (vector<int>::iterator iele = acc_probe_ele.begin(); iele != acc_probe_ele.end(); ++iele) {
	//if ( *itag >= acc_probe_ele[iele]) continue; //leading electron is tag
	if (*itag >= *iele) continue;
	if ( elePt[*iele] < 10.) continue;

	TLorentzVector vProbe ,vZ ;
	vProbe.SetPtEtaPhiM(elePt[*iele], eleEta[*iele], elePhi[*iele],0.511*0.001);
	vZ = vTag + vProbe;
	if ( vZ.M() < 60 || vZ.M() > 120 ) continue ;

	Zmass.push_back(vZ.M());
        ZEta.push_back(vZ.Eta());
        ZPt.push_back(vZ.Pt());
        index.push_back(*iele);
      }

      nZ = Zmass.size();
      if (nZ < 1) continue;
      float minZ = 999.;
      int iProbe = -1;
      for (int i = 0; i < nZ; i++) {
        float deltaM = fabs(91.19-Zmass[i]);
        if (minZ > deltaM) {
          minZ = deltaM;
          iProbe = index[i];
          Zm = Zmass[i];
          Zeta = ZEta[i];
          //Zpt = ZPt[i];
        }
      }

      if (iProbe < 0) continue;

      tag_elePt = elePt[*itag];
      tag_eleEta = eleEta[*itag];
      tag_eleSCEta= eleSCEta[*itag];
      tag_elePhi = elePhi[*itag];
	
      ProbeIndex = iProbe;
      probe_elePt = elePt[iProbe];
      probe_eleEta = eleEta[iProbe];
      probe_eleSCEta= eleSCEta[iProbe];
      probe_elePhi = elePhi[iProbe];
      probe_MatchFilter_Leg1_Ele23 = eleFiredDoubleTrgs[iProbe] >> 2&1;
      probe_MatchFilter_Leg2_Ele12 = eleFiredDoubleTrgs[iProbe] >> 3&1;
      probe_MatchFilter_Leg1_Ele16_Ele12_Ele8 = eleFiredDoubleTrgs[iProbe] >> 9&1;
      probe_MatchFilter_Leg2_Ele16_Ele12_Ele8 = eleFiredDoubleTrgs[iProbe] >> 10&1;
      probe_MatchFilter_Leg3_Ele16_Ele12_Ele8 = eleFiredDoubleTrgs[iProbe] >> 11&1;

      npair++;
      trieletree->Fill();

    } //end for tri-ele tree

    /*----------- for photon trigger, leading ele matching filter of SingleEle trigger ------------*/
    for (vector<int>::iterator itag = acc_tag_ele.begin() ; itag != acc_tag_ele.end(); ++itag) {

      //initial values
      ProbeIndex = -99.;
      tag_elePt = -99.;
      tag_eleEta = -99.;
      tag_eleSCEta = -99.;
      tag_elePhi = -99.;
      tag_ele_match_Leg1_DiPho70 = 0;
      tag_ele_match_Leg1_TriPho20_CaloId = 0;
      tag_ele_match_Leg1_TriPho20_CaloR9Id = 0;
      tag_ele_match_Leg1_TriPho303010_CaloId = 0;
      tag_ele_match_Leg1_TriPho303010_CaloR9Id = 0;
      tag_ele_match_Leg1_TriPho35355_CaloR9Id = 0;
      
      probe_phoPt = -99.;
      probe_phoEta = -99.;
      probe_phoSCEta = -99.;
      probe_phoPhi = -99.;
      probe_MatchFilter_Leg1_DiPho70 = 0;
      probe_MatchFilter_Leg2_DiPho70 = 0;
      probe_MatchFilter_Leg1_TriPho20_CaloId = 0;
      probe_MatchFilter_Leg2_TriPho20_CaloId = 0;
      probe_MatchFilter_Leg3_TriPho20_CaloId = 0;

      probe_MatchFilter_Leg1_TriPho20_CaloR9Id = 0;
      probe_MatchFilter_Leg2_TriPho20_CaloR9Id = 0;
      probe_MatchFilter_Leg3_TriPho20_CaloR9Id = 0;

      probe_MatchFilter_Leg1_TriPho303010_CaloId = 0;
      probe_MatchFilter_Leg2_TriPho303010_CaloId = 0;
      probe_MatchFilter_Leg3_TriPho303010_CaloId = 0;
      probe_MatchFilter_TriPho303010_CaloId = 0;

      probe_MatchFilter_Leg1_TriPho303010_CaloR9Id = 0;
      probe_MatchFilter_Leg2_TriPho303010_CaloR9Id = 0;
      probe_MatchFilter_Leg3_TriPho303010_CaloR9Id = 0;
      probe_MatchFilter_TriPho303010_CaloR9Id = 0;

      probe_MatchFilter_Leg1_TriPho35355_CaloR9Id = 0;
      probe_MatchFilter_Leg2_TriPho35355_CaloR9Id = 0;
      probe_MatchFilter_Leg3_TriPho35355_CaloR9Id = 0;
      probe_MatchFilter_TriPho35355_CaloR9Id = 0;


      Zm_matchPho = -99.;
      Zeta_matchPho = -99.;
      Zy_matchPho = -99.;
      TagIndex = -99;
   
      if ( elePt[*itag] < 30.) continue;
      nTag++;
      TagIndex = *itag;
      
      TLorentzVector vTag ;
      vTag.SetPtEtaPhiM(elePt[*itag], eleEta[*itag], elePhi[*itag], 0.511*0.001);
   
      vector<float> Zmass, ZEta, ZPt;
      vector<int> index;
      
      //for (int iele = 0; iele < nEle; iele++) {
      //if (*itag >= iele) continue;
      //if (elePt[iele] < 5.) continue;
      //if (fabs(eleSCEta[iele]) > 2.5) continue;
      //if (fabs(eleSCEta[iele]) > 1.4442 && fabs(eleSCEta[iele]) < 1.566) continue;
	
	if (nPho < 1) continue;
	for (int ipho = 0; ipho < nPho; ipho++) {
	  if (phoCalibEt[ipho] < 5.) continue;
	  if (fabs(phoSCEta[ipho]) > 2.5) continue;
	  if (fabs(phoSCEta[ipho]) > 1.4442 && fabs(phoSCEta[ipho]) < 1.566) continue;
	  if ( phoIDMVA[ipho] < 0.2) continue;   
	  //if ( fabs(phoSCEta[ipho] - eleSCEta[iele]) >= 0.05 || fabs(deltaPhi(phoSCPhi[ipho], eleSCPhi[iele]) >= 0.05) ) continue;
	  //if (phoEleVeto[ipho] !=1 ) continue;
	  
	  TLorentzVector vProbe_matchPho ,vZ ;
	  vProbe_matchPho.SetPtEtaPhiM(phoCalibEt[ipho] , phoEta[ipho] , phoPhi[ipho] , 0. );
	  vZ = vTag + vProbe_matchPho ;
	  if (vZ.M()< 60 || vZ.M() > 120 ) continue ;
	  
	  Zmass.push_back(vZ.M());
	  ZEta.push_back(vZ.Eta());
	  ZPt.push_back(vZ.Pt());
	  index.push_back(ipho);
	}
	//}
      
      nZ = Zmass.size();
      if (nZ < 1) continue;
      float minZ = 999.;
      int iProbe = -1;
      for (int i = 0; i < nZ; i++) {
	float deltaM = fabs(91.19-Zmass[i]);
	if (minZ > deltaM) {
	  minZ = deltaM;
	  iProbe = index[i];
	  Zm_matchPho = Zmass[i];
	  Zeta_matchPho = ZEta[i];
	  //Zy_matchPho = ZRapidity[i];
	}
      }
      
      if (iProbe < 0) continue;
      //veto photon matched to leading electron
      if ( fabs(phoSCEta[iProbe] - eleSCEta[*itag]) <= 0.05 || fabs(deltaPhi(phoSCPhi[iProbe], eleSCPhi[*itag]) <= 0.05) ) continue;

      tag_elePt = elePt[*itag];
      tag_eleEta = eleEta[*itag];
      tag_eleSCEta= eleSCEta[*itag];
      tag_elePhi = elePhi[*itag];
      //matching trigger for tag ele
      tag_ele_match_Leg1_DiPho70 = phoFiredDoubleTrgs[*itag]>>30&1;
      tag_ele_match_Leg1_TriPho20_CaloId = phoFiredTripleTrgs[*itag]>>0&1;
      tag_ele_match_Leg1_TriPho20_CaloR9Id = phoFiredTripleTrgs[*itag]>>2&1;
      tag_ele_match_Leg1_TriPho303010_CaloId = phoFiredTripleTrgs[*itag]>>4&1;
      tag_ele_match_Leg1_TriPho303010_CaloR9Id = phoFiredTripleTrgs[*itag]>>8&1;
      tag_ele_match_Leg1_TriPho35355_CaloR9Id = phoFiredTripleTrgs[*itag]>>12&1;
      
      probe_phoPt = phoCalibEt[iProbe];
      probe_phoEta = phoEta[iProbe];
      probe_phoSCEta = phoSCEta[iProbe];
      probe_phoPhi = phoPhi[iProbe];
      probe_MatchFilter_Leg1_DiPho70 = phoFiredDoubleTrgs[iProbe]>>30&1;
      probe_MatchFilter_Leg2_DiPho70 = phoFiredDoubleTrgs[iProbe]>>31&1;

      probe_MatchFilter_Leg1_TriPho20_CaloId = phoFiredTripleTrgs[iProbe]>>0&1;
      probe_MatchFilter_Leg2_TriPho20_CaloId = phoFiredTripleTrgs[iProbe]>>1&1;

      probe_MatchFilter_Leg1_TriPho20_CaloR9Id = phoFiredTripleTrgs[iProbe]>>2&1;
      probe_MatchFilter_Leg2_TriPho20_CaloR9Id = phoFiredTripleTrgs[iProbe]>>3&1;

      probe_MatchFilter_Leg1_TriPho303010_CaloId = phoFiredTripleTrgs[iProbe]>>4&1;
      probe_MatchFilter_Leg2_TriPho303010_CaloId = phoFiredTripleTrgs[iProbe]>>5&1;
      probe_MatchFilter_Leg3_TriPho303010_CaloId = phoFiredTripleTrgs[iProbe]>>6&1;
      probe_MatchFilter_TriPho303010_CaloId = phoFiredTripleTrgs[iProbe]>>7&1;

      probe_MatchFilter_Leg1_TriPho303010_CaloR9Id = phoFiredTripleTrgs[iProbe]>>8&1;
      probe_MatchFilter_Leg2_TriPho303010_CaloR9Id = phoFiredTripleTrgs[iProbe]>>9&1;
      probe_MatchFilter_Leg3_TriPho303010_CaloR9Id = phoFiredTripleTrgs[iProbe]>>10&1;
      probe_MatchFilter_TriPho303010_CaloR9Id = phoFiredTripleTrgs[iProbe]>>11&1;

      probe_MatchFilter_Leg1_TriPho35355_CaloR9Id = phoFiredTripleTrgs[iProbe]>>12&1;
      probe_MatchFilter_Leg2_TriPho35355_CaloR9Id = phoFiredTripleTrgs[iProbe]>>13&1;
      probe_MatchFilter_Leg3_TriPho35355_CaloR9Id = phoFiredTripleTrgs[iProbe]>>14&1;
      probe_MatchFilter_TriPho35355_CaloR9Id = phoFiredTripleTrgs[iProbe]>>15&1;

      npaireg++;
      
      photree->Fill();
    }

    //cout << "loop on tag electron matching to filter of pho" << endl;
    /*-------------------------------------------------------------------------------------------- */
    /*----------- for photon trigger, leading ele match first filter of TriPho trigger ------------*/
    /*for (vector<int>::iterator itag = acc_tag_ele_leg2_tripho20Calo.begin() ; itag != acc_tag_ele_leg2_tripho20Calo.end(); ++itag) {

      //initial values
      ProbeIndex = -99.;
      tag_elePt = -99.;
      tag_eleEta = -99.;
      tag_eleSCEta = -99.;
      tag_elePhi = -99.;
      
      probe_phoPt = -99.;
      probe_phoEta = -99.;
      probe_phoSCEta = -99.;
      probe_phoPhi = -99.;
      probe_MatchFilter_Leg1_DiPho70 = 0;
      probe_MatchFilter_Leg2_DiPho70 = 0;
      probe_MatchFilter_Leg1_TriPho20_CaloId = 0;
      probe_MatchFilter_Leg2_TriPho20_CaloId = 0;
      probe_MatchFilter_Leg3_TriPho20_CaloId = 0;

      probe_MatchFilter_Leg1_TriPho20_CaloR9Id = 0;
      probe_MatchFilter_Leg2_TriPho20_CaloR9Id = 0;
      probe_MatchFilter_Leg3_TriPho20_CaloR9Id = 0;

      probe_MatchFilter_Leg1_TriPho303010_CaloId = 0;
      probe_MatchFilter_Leg2_TriPho303010_CaloId = 0;
      probe_MatchFilter_Leg3_TriPho303010_CaloId = 0;
      probe_MatchFilter_TriPho303010_CaloId = 0;

      probe_MatchFilter_Leg1_TriPho303010_CaloR9Id = 0;
      probe_MatchFilter_Leg2_TriPho303010_CaloR9Id = 0;
      probe_MatchFilter_Leg3_TriPho303010_CaloR9Id = 0;
      probe_MatchFilter_TriPho303010_CaloR9Id = 0;

      probe_MatchFilter_Leg1_TriPho35355_CaloR9Id = 0;
      probe_MatchFilter_Leg2_TriPho35355_CaloR9Id = 0;
      probe_MatchFilter_Leg3_TriPho35355_CaloR9Id = 0;
      probe_MatchFilter_TriPho35355_CaloR9Id = 0;


      Zm_matchPho = -99.;
      Zeta_matchPho = -99.;
      Zy_matchPho = -99.;
      TagIndex = -99;
   
      if ( elePt[*itag] < 30.) continue;
      nTag++;
      TagIndex = *itag;
      
      TLorentzVector vTag ;
      vTag.SetPtEtaPhiM(elePt[*itag], eleEta[*itag], elePhi[*itag], 0.511*0.001);
   
      vector<float> Zmass, ZEta, ZPt;
      vector<int> index;
      
      for (int iele = 0; iele < nEle; iele++) {
	if (*itag >= iele) continue;
	if (elePt[iele] < 5.) continue;
	if (fabs(eleSCEta[iele]) > 2.5) continue;
	if (fabs(eleSCEta[iele]) > 1.4442 && fabs(eleSCEta[iele]) < 1.566) continue;
	
	if (nPho < 1) continue;
	for (int ipho = 0; ipho < nPho; ipho++) {
	  if (phoCalibEt[ipho] < 5.) continue;
	  if (fabs(phoSCEta[ipho]) > 2.5) continue;
	  if (fabs(phoSCEta[ipho]) > 1.4442 && fabs(phoSCEta[ipho]) < 1.566) continue;
	  if ( phoIDMVA[ipho] < 0.2) continue;   
	  if ( fabs(phoSCEta[ipho] - eleSCEta[iele]) >= 0.05 || fabs(deltaPhi(phoSCPhi[ipho], eleSCPhi[iele]) >= 0.05) ) continue;
	  if (phoEleVeto[ipho] !=1 ) continue;
	  
	  TLorentzVector vProbe_matchPho ,vZ ;
	  vProbe_matchPho.SetPtEtaPhiM(phoCalibEt[ipho] , phoEta[ipho] , phoPhi[ipho] , 0. );
	  vZ = vTag + vProbe_matchPho ;
	  if (vZ.M()< 60 || vZ.M() > 120 ) continue ;
	  
	  Zmass.push_back(vZ.M());
	  ZEta.push_back(vZ.Eta());
	  ZPt.push_back(vZ.Pt());
	  index.push_back(ipho);
	}
      }
      
      nZ = Zmass.size();
      if (nZ < 1) continue;
      float minZ = 999.;
      int iProbe = -1;
      for (int i = 0; i < nZ; i++) {
	float deltaM = fabs(91.19-Zmass[i]);
	if (minZ > deltaM) {
	  minZ = deltaM;
	  iProbe = index[i];
	  Zm_matchPho = Zmass[i];
	  Zeta_matchPho = ZEta[i];
	  //Zy_matchPho = ZRapidity[i];
	}
      }
      
      if (iProbe < 0) continue;
      //veto photon matched to leading electron
      if ( fabs(phoSCEta[iProbe] - eleSCEta[*itag]) <= 0.05 || fabs(deltaPhi(phoSCPhi[iProbe], eleSCPhi[*itag]) <= 0.05) ) continue;

      tag_elePt = elePt[*itag];
      tag_eleEta = eleEta[*itag];
      tag_eleSCEta= eleSCEta[*itag];
      tag_elePhi = elePhi[*itag];
      
      probe_phoPt = phoCalibEt[iProbe];
      probe_phoEta = phoEta[iProbe];
      probe_phoSCEta = phoSCEta[iProbe];
      probe_phoPhi = phoPhi[iProbe];
      probe_MatchFilter_Leg1_DiPho70 = phoFiredDoubleTrgs[iProbe]>>30&1;
      probe_MatchFilter_Leg2_DiPho70 = phoFiredDoubleTrgs[iProbe]>>31&1;

      probe_MatchFilter_Leg1_TriPho20_CaloId = phoFiredTripleTrgs[iProbe]>>0&1;
      probe_MatchFilter_Leg2_TriPho20_CaloId = phoFiredTripleTrgs[iProbe]>>1&1;

      probe_MatchFilter_Leg1_TriPho20_CaloR9Id = phoFiredTripleTrgs[iProbe]>>2&1;
      probe_MatchFilter_Leg2_TriPho20_CaloR9Id = phoFiredTripleTrgs[iProbe]>>3&1;

      probe_MatchFilter_Leg1_TriPho303010_CaloId = phoFiredTripleTrgs[iProbe]>>4&1;
      probe_MatchFilter_Leg2_TriPho303010_CaloId = phoFiredTripleTrgs[iProbe]>>5&1;
      probe_MatchFilter_Leg2_TriPho303010_CaloId = phoFiredTripleTrgs[iProbe]>>6&1;
      probe_MatchFilter_TriPho303010_CaloId = phoFiredTripleTrgs[iProbe]>>7&1;

      probe_MatchFilter_Leg1_TriPho303010_CaloR9Id = phoFiredTripleTrgs[iProbe]>>8&1;
      probe_MatchFilter_Leg2_TriPho303010_CaloR9Id = phoFiredTripleTrgs[iProbe]>>9&1;
      probe_MatchFilter_Leg2_TriPho303010_CaloR9Id = phoFiredTripleTrgs[iProbe]>>10&1;
      probe_MatchFilter_TriPho303010_CaloR9Id = phoFiredTripleTrgs[iProbe]>>11&1;

      probe_MatchFilter_Leg1_TriPho35355_CaloR9Id = phoFiredTripleTrgs[iProbe]>>12&1;
      probe_MatchFilter_Leg2_TriPho35355_CaloR9Id = phoFiredTripleTrgs[iProbe]>>13&1;
      probe_MatchFilter_Leg2_TriPho35355_CaloR9Id = phoFiredTripleTrgs[iProbe]>>14&1;
      probe_MatchFilter_TriPho35355_CaloR9Id = phoFiredTripleTrgs[iProbe]>>15&1;

      npaireg++;
      
      pho20Calotree->Fill();
    }
    */

    //muon loop
    vector<int> acc_tag_mu, acc_probe_mu;

    for (int i=0; i<nMu; i++) {
      if (muPt[i] < 10.) continue;
      if (fabs(muEta[i]) > 2.4) continue;
      if (mc) {
        bool match = muMatcher(data1, i);
        if (!match) continue;
      }

      //tag muon: tight id + tight iso
      if ((muIDbit[i] >> 2 & 1) == 1) {
	if ((muPFChIso[i] + TMath::Max(0., muPFPhoIso[i] + muPFNeuIso[i] -0.5 * muPFPUIso[i]))/muPt[i] < 0.15) {
	  //if ( (muFiredTrgs[i] >> 17&1) || (muFiredTrgs[i] >> 19&1) ) acc_tag_mu.push_back(i);
	  //if ( (muFiredTrgs[i] >> 0&1) || (muFiredTrgs[i] >> 1&1) ) acc_tag_mu.push_back(i);
	  if ( (muFiredTrgs[i] >> 29&1) ) acc_tag_mu.push_back(i);
	  acc_tag_mu.push_back(i); 
	}
      }

      //probe muon: medium id + medium iso
      if ((muIDbit[i] >> 1 & 1) == 1) {
	if ((muPFChIso[i] + TMath::Max(0., muPFPhoIso[i] + muPFNeuIso[i] -0.5 * muPFPUIso[i]))/muPt[i] < 0.25) {
	  acc_probe_mu.push_back(i);
	}
      }
    } // end loop on muon

    //cout << "number of muon tag: " << acc_tag_mu.size() << "\t probe: " << acc_probe_mu.size() << endl;
    if (acc_tag_mu.size() < 1) continue;
    if (acc_probe_mu.size() < 1) continue;

    nZ = 0;
    nmu = 0;
    npair = 0;

    for (unsigned int itag = 0 ; itag < acc_tag_mu.size(); itag++) {

      //initial values
      TagIndex = -99;
      tag_muPt = -99.;
      tag_muEta = -99.;
      tag_muPhi = -99.;

      ProbeIndex = -99.;
      probe_muPt = -99.;
      probe_muEta = -99.;
      probe_muPhi = -99.;
      probe_MatchFilter_Leg1_TriMu = 0;
      probe_MatchFilter_Leg2_TriMu = 0;
      probe_MatchFilter_Leg3_TriMu = 0;
      probe_MatchFilter_DoubleMu207 = 0;
      Zm = -99.;
      Zeta = -99.;
      //Zpt = -99.;


      if ( muPt[acc_tag_mu[itag]] < 27.) continue;
      nTag++;
      TagIndex = acc_tag_mu[itag];

      //cout << "passing tag selection" << endl;
      TLorentzVector vTag ;
      vTag.SetPtEtaPhiM(muPt[acc_tag_mu[itag]], muEta[acc_tag_mu[itag]], muPhi[acc_tag_mu[itag]], 0.511*0.001);

      vector<float> Zmass, ZEta, ZPt;
      vector<int> index;
      nProbe = 0;

      for (unsigned int imu = 0; imu < acc_probe_mu.size(); imu++) {
	if ( acc_tag_mu[itag] >= acc_probe_mu[imu]) continue; //leading muon is tag
	if ( muPt[acc_probe_mu[imu]] < 10.) continue;

	TLorentzVector vProbe ,vZ ;
	vProbe.SetPtEtaPhiM(muPt[acc_probe_mu[imu]], muEta[acc_probe_mu[imu]], muPhi[acc_probe_mu[imu]],0.511*0.001);
	vZ = vTag + vProbe;
	if ( vZ.M() < 60 || vZ.M() > 120 ) continue ;

	Zmass.push_back(vZ.M());
        ZEta.push_back(vZ.Eta());
        ZPt.push_back(vZ.Pt());
        index.push_back(acc_probe_mu[imu]);
      }

      nZ = Zmass.size();
      if (nZ < 1) continue;
      float minZ = 999.;
      int iProbe = -1;
      for (int i = 0; i < nZ; i++) {
        float deltaM = fabs(91.19-Zmass[i]);
        if (minZ > deltaM) {
          minZ = deltaM;
          iProbe = index[i];
          Zm = Zmass[i];
          Zeta = ZEta[i];
          //Zpt = ZPt[i];
        }
      }

      if (iProbe < 0) continue;
      if ( acc_tag_mu[itag] == iProbe) continue;

      tag_muPt = muPt[acc_tag_mu[itag]];
      tag_muEta = muEta[acc_tag_mu[itag]];
      tag_muPhi = muPhi[acc_tag_mu[itag]];
      //mu1_MatchTrig27 = muFiredSingleTrgs[acc_tag_mu[itag]] >>12&1;
	
      probe_muPt = muPt[iProbe];
      probe_muEta = muEta[iProbe];
      probe_muPhi = muPhi[iProbe];
      probe_MatchFilter_Leg1_TriMu = muFiredTrgs[iProbe] >> 2&1;
      probe_MatchFilter_Leg2_TriMu = muFiredTrgs[iProbe] >> 3&1;
      probe_MatchFilter_Leg3_TriMu = muFiredTrgs[iProbe] >> 3&1;
      probe_MatchFilter_DoubleMu207 = muFiredTrgs[iProbe] >> 30&1;

      npair++;
      //cout << "filling tnp tree for e-e" << endl;
      mutree->Fill();
    }

    //mu-pho tree
    for (vector<int>::iterator itag = acc_tag_mu.begin() ; itag != acc_tag_mu.end(); ++itag) {

      //initial values
      tag_muPt = -99.;
      tag_muEta = -99.;
      tag_muPhi = -99.;

      probe_muPt = -99.;
      probe_muEta = -99.;
      probe_muPhi = -99.;
      probe_MatchFilter_Leg1_TriMu = 0;
      probe_MatchFilter_Leg2_TriMu = 0;
      probe_MatchFilter_Leg3_TriMu = 0;
      probe_MatchFilter_DoubleMu207 = 0;
      Zm = -99.;
      Zeta = -99.;
      //Zpt = -99.;

      probe_phoPt = -99.;
      probe_phoEta = -99.;
      probe_phoSCEta = -99.;
      probe_phoPhi = -99.;
      probe_MatchFilter_Pho23 = -99.;


      if ( muPt[*itag] < 27.) continue;
      nTag++;
      TagIndex = *itag;

      //cout << "passing tag selection" << endl;
      TLorentzVector vTag ;
      vTag.SetPtEtaPhiM(muPt[*itag], muEta[*itag], muPhi[*itag], 0.511*0.001);

      vector<float> Zmass, ZEta, ZPt;
      vector<int> index;
      int idxPho = -1;

      for (vector<int>::iterator imu = acc_probe_mu.end(); imu != acc_probe_mu.end(); ++imu) {
	if ( *itag >= *imu) continue; //leading muon is tag
	if ( muPt[*imu] < 10.) continue;

	TLorentzVector vProbe ,vZ , vPho;
	vProbe.SetPtEtaPhiM(muPt[*imu], muEta[*imu], muPhi[*imu],0.511*0.001);
	int npho = 0;
	for (int ipho = 0; ipho < nPho; ipho++) {
	  if (phoCalibEt[ipho] < 15.) continue;
	  if (fabs(phoSCEta[ipho]) > 2.5) continue;
	  if (fabs(phoSCEta[ipho]) > 1.4442 && fabs(phoSCEta[ipho]) < 1.566) continue;
	  if ( phoIDMVA[ipho] < 0.2) continue;   
	  if (phoEleVeto[ipho] !=1 ) continue;
	  if (deltaR(phoEta[ipho], phoPhi[ipho], muEta[*imu], muPhi[*imu]) > 0.3) continue; //choose FSR pho
	  if (deltaR(phoEta[ipho], phoPhi[ipho], muEta[*itag], muPhi[*itag]) > 0.3) continue; //choose FSR pho
	  vPho.SetPtEtaPhiM(phoCalibEt[ipho], phoEta[ipho], phoPhi[ipho], 0.);
	  idxPho = ipho;
	  npho++;
	}
	if (npho > 1) continue;  //select only one pho
	vZ = vTag + vProbe + vPho;
	if ( vZ.M() < 60 || vZ.M() > 120 ) continue ;

	Zmass.push_back(vZ.M());
        ZEta.push_back(vZ.Eta());
        ZPt.push_back(vZ.Pt());
        index.push_back(*imu);
      }

      nZ = Zmass.size();
      if (nZ < 1) continue;
      float minZ = 999.;
      int iProbe = -1;
      for (int i = 0; i < nZ; i++) {
        float deltaM = fabs(91.19-Zmass[i]);
        if (minZ > deltaM) {
          minZ = deltaM;
          iProbe = index[i];
          Zm = Zmass[i];
          Zeta = ZEta[i];
          //Zpt = ZPt[i];
        }
      }

      if (iProbe < 0) continue;
      if ( *itag == iProbe) continue;

      tag_muPt = muPt[*itag];
      tag_muEta = muEta[*itag];
      tag_muPhi = muPhi[*itag];
      //mu1_MatchTrig27 = muFiredSingleTrgs[*itag] >>12&1;
	
      probe_muPt = muPt[iProbe];
      probe_muEta = muEta[iProbe];
      probe_muPhi = muPhi[iProbe];
      probe_MatchFilter_Leg1_TriMu = muFiredTrgs[iProbe] >> 2&1;
      probe_MatchFilter_Leg2_TriMu = muFiredTrgs[iProbe] >> 3&1;
      probe_MatchFilter_Leg3_TriMu = muFiredTrgs[iProbe] >> 3&1;
      probe_MatchFilter_DoubleMu207 = muFiredTrgs[iProbe] >> 30&1;

      probe_phoPt = phoCalibEt[idxPho];
      probe_phoEta = phoEta[idxPho];
      probe_phoSCEta = phoSCEta[idxPho];
      probe_phoPhi = phoPhi[idxPho];
      probe_MatchFilter_Pho23 = phoFiredSingleTrgs[idxPho] >> 15&1;

      npair++;
      //cout << "filling tnp tree for e-e" << endl;
      muphotree->Fill();
    }

  }

  fo->Write();
  fo->Close();
  
}

