{
  gSystem->SetBuildDir("tmpdir", kTRUE);
  gROOT->LoadMacro("makeTnPtree.C+");

  runData(0);
  //runMC(1);

}

void runData(bool mc = 0) {

  TString inpaths;

  //inpaths = "/home/tdoan/analysis/triggerFor2017/checkPerformance/samples/run2017_eraB/DCS_JSON/DoubleEG/ggtree_data_DoubleEG*.root"; 
  //inpaths = "/home/tdoan/analysis/triggerFor2017/checkPerformance/samples/run2017_eraB/StablePixel/SingleMuon/ggtree_data_DoubleEG*.root";  
  //inpaths = "/home/tdoan/analysis/triggerFor2017/checkPerformance/samples/run2017_eraB/StablePixel/SingleElectron/ggtree_data_DoubleEG*.root";
  //inpaths = "root://eoscms.cern.ch//eos/cms/store/group/phys_smp/ggNtuples/13TeV/data/V08_00_26_01/job_SingleEle_Run2016D_FebReminiAOD/ggtree_data_*.root";
  //inpaths = "/home/tdoan/samples/data/2017/merged/SingleEle_Run2017B_PromptReco_v1_Aug2.root";
  //inpaths = "/home/tdoan/analysis/triggerFor2017/checkPerformance/samples/run2016H/PromptReco_v2/ggtree_data_*.root";
  //inpaths = "/home/tdoan/samples/data/2017/SingleEle_2017A_PromptReco_v2_DCSJSON/ggtree_data_*.root";
  //inpaths = "/home/tdoan/samples/data/2017/SingleEle_2017A_PromptReco_v3_DCSJSON/ggtree_data_*.root";
  //inpaths = "/home/tdoan/samples/data/2017/SingleEle_2017B_PromptReco_v1_DCSJSON/ggtree_data_*.root";
  //inpaths = "/home/tdoan/samples/data/2017/SingleEle_2017B_PromptReco_v2_DCSJSON/ggtree_data_*.root";
  //inpaths = "/home/tdoan/samples/data/2017/SingleEle_2017C_PromptReco_v1_DCSJSON/ggtree_data_*.root";
  //inpaths = "/home/tdoan/samples/data/2017/SingleEle_2017C_PromptReco_v2_DCSJSON/ggtree_data_*.root";

  //inpaths = "/home/tdoan/samples/data/2017/SingleEle_2017D_PromptReco_v1_DCSJSON/ggtree_data_*.root";
  //inpaths = "/home/tdoan/samples/data/2017/SingleMuon_2017A_PromptReco_v2_DCSJSON/ggtree_data_*.root";
  //inpaths = "/home/tdoan/samples/data/2017/SingleMuon_2017A_PromptReco_v3_DCSJSON/ggtree_data_*.root";

  //inpaths = "/data1/ggNtuples/V929/SingleElectron_Run2017D_PromptReco_v1_Nov22/ggtree_data_*.root";
  inpaths = "/data1/ggNtuples/V929/SingleElectron_Run2017E_PromptReco_v1_Nov22/ggtree_data_*.root";
  //inpaths = "/data1/ggNtuples/V929/SingleElectron_Run2017F_PromptReco_v1_Nov22/ggtree_data_*.root";


  cout << "Process sample: " << inpaths << endl;
  makeTnPtree(inpaths, "out.root",  mc);
}

void runMC(bool mc = 1) {

  TString inpaths;
  //inpaths = "/data6/ggNtuples/V08_00_24_00/summer16_Zg_aMCatNLO/ggtree_mc*.root";
  //inpaths = "/data6/ggNtuples/V08_00_24_00/summer16_Zg_pt130/ggtree_mc*.root";
  //inpaths = "/data6/ggNtuples/V08_00_24_00/summer16_DYJetsToLL_m50_MG/ggtree_mc*.root";
  //inpaths = "/data6/ggNtuples/V08_00_24_00/spring16_DYJetsToLL_m50_aMCatNLO/ggtree_mc*.root";
  //inpaths = "/data6/ggNtuples/V08_00_24_00/spring16_DYJetsToLL_m50_aMCatNLO_MiniIso0p05/ggtree_mc*.root";
  inpaths = "/data6/ggNtuples/V08_00_26_01/job_summer16_DYJets_M50_aMCatNLO/ggtree_mc*.root";

  cout << "Process sample: " << inpaths << endl;
  makeTnPtree(inpaths, "out.root",  mc);

}
