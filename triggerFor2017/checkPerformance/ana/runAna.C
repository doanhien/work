{
  gSystem->SetBuildDir("tmpdir", kTRUE);
  gROOT->LoadMacro("xAna.C+");

  runData(0, 1);
  //runMC(1, 0);

}

void runData(bool mc = 0, bool ele = 0) {

  TString inpaths;

  //inpaths = "/home/tdoan/analysis/triggerFor2017/checkPerformance/samples/run2017_eraB/DCS_JSON/DoubleEG/ggtree_data_DoubleEG*.root";
  //inpaths = "/home/tdoan/analysis/triggerFor2017/checkPerformance/samples/run2017_eraB/StablePixel/SingleMuon/ggtree_data_DoubleEG*.root";
  //inpaths = "/home/tdoan/analysis/triggerFor2017/checkPerformance/samples/run2017_eraB/StablePixel/SingleElectron/ggtree_data_DoubleEG*.root";
  //inpaths = "/home/tdoan/analysis/triggerFor2017/checkPerformance/samples/run2017_eraB/MuonEG/ggtree_data_DoubleEG*.root";
  //inpaths = "root://eoscms.cern.ch//eos/cms/store/user/tdoan//eos/cms/store/user/tdoan/job_JetHLT_Run2017A_296173_296174/JetHT/crab_job_JetHLT_Run2017A_296173_296174/170611_105930/0000/ggtree_data_*.root";
  //inpaths = "/home/tdoan/samples/data/2017/merged/SingleEle_Run2017B_PromptReco_v1_Aug2.root";
  //inpaths = "/home/tdoan/samples/data/2017/MET_Run2017A_PromptReco_v2/ggtree_data_*.root";
  //inpaths = "/home/tdoan/samples/data/2017/MET_Run2017A_PromptReco_v3/ggtree_data_*.root";
  //inpaths = "/home/tdoan/samples/data/2017/MET_Run2017C_PromptReco_v2/ggtree_data_*.root";
  //inpaths = "/home/tdoan/samples/data/2017/MET_Run2017C_PromptReco_v1/ggtree_data_*.root";
  //inpaths = "/home/tdoan/samples/data/2017/MET_Run2017B_PromptReco_v1/ggtree_data_*.root";
  //inpaths = "/home/tdoan/samples/data/2017/MET_Run2017B_PromptReco_v2/ggtree_data_*.root";

  //inpaths = "/home/tdoan/samples/data/2017/SingleEle_2017A_PromptReco_v2_DCSJSON/ggtree_data_*.root";
  inpaths = "/home/tdoan/samples/data/2017/SingleEle_2017A_PromptReco_v3_DCSJSON/ggtree_data_*.root";
  //inpaths = "/home/tdoan/samples/data/2017/SingleEle_2017B_PromptReco_v1_DCSJSON/ggtree_data_*.root";
  //inpaths = "/home/tdoan/samples/data/2017/SingleEle_2017B_PromptReco_v2_DCSJSON/ggtree_data_*.root";
  //inpaths = "/home/tdoan/samples/data/2017/SingleEle_2017C_PromptReco_v1_DCSJSON/ggtree_data_*.root";
  //inpaths = "/home/tdoan/samples/data/2017/SingleEle_2017C_PromptReco_v2_DCSJSON/ggtree_data_*.root"; 


  cout << "Process sample: " << inpaths << endl;
  xAna(inpaths, "out.root",  mc, ele);
}

void runMC(bool mc = 1, bool ele = 1) {

  TString inpaths;
  //inpaths = "/data6/ggNtuples/V08_00_26_01/job_summer16_DYJets_M50_aMCatNLO/ggtree_mc*.root";
  inpaths = "root://eoscms.cern.ch//eos/cms/store/user/tdoan/job_spring17_DYJetsToLL/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_job_spring17_DYJetsToLL/170611_144132/0000/ggtree_mc*.root";

  cout << "Process sample: " << inpaths << endl;
  xAna(inpaths, "out.root",  mc, ele);

}
