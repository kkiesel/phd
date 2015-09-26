// Shared libray should be produced with compiling
// https://github.com/cms-susy-photon-rwth-1b/TreeWriter
// The lib should be in the CMSSW lib path

//R__LOAD_LIBRARY(TreeParticles)

// start with root -l -q run.C
// Using ROOT6

void run( string infile="/user/kiesel/nTuples/V02/SinglePhoton_V02_nTuple.root") {

  ifstream f(infile.c_str());
  if(!f.good()) {
    cout << "No file found, just compile." << endl;
    TSelector::GetSelector("HistogramProducer.cc+");
    return;
  }

  gSystem->Load("pluginTreeWriterTreeWriterAuto.so");

  TChain ch( "TreeWriter/eventTree" );
  ch.AddFile( infile.c_str() );

  double start_time = time(NULL);
  ch.Process( "HistogramProducer.cc+" );
  cout << "Job needed " << 1.*( time(NULL) - start_time)/60 << " min real time." << endl;

}
