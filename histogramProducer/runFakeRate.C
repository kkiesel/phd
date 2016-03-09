// Shared libray should be produced with compiling
// https://github.com/cms-susy-photon-rwth-1b/TreeWriter
// The lib should be in the CMSSW lib path
// Alternative: Compile https://github.com/cms-susy-photon-rwth-1b/tools/tree/master/templates

// start with root -l -q run.C

int runFakeRate( string infile="/user/kiesel/nTuples/v09/SingleElectron_Run2015C_25ns-05Oct2015-v1_nTuple.root") {

  ifstream f(infile.c_str());
  if(!f.good()) {
    cout << "No file found, just compile." << endl;
    auto out = TSelector::GetSelector("FakeRateEstimator.cc+");
    if(!out) return 1;
  }

  gSystem->Load("pluginTreeWriterTreeWriterAuto.so");

  TChain ch( "TreeWriter/eventTree" );
  ch.AddFile( infile.c_str() );

  double start_time = time(NULL);
  ch.Process( "FakeRateEstimator.cc+" );
  cout << " in " << 1.*( time(NULL) - start_time)/60 << " min" << endl;

  return 0;
}
