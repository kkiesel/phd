// Shared libray should be produced with compiling
// https://github.com/cms-susy-photon-rwth-1b/TreeWriter
// The lib should be in the CMSSW lib path

//R__LOAD_LIBRARY(TreeParticles)

// start with root -l -q run.C
// Using ROOT6
void run( string infile="/user/kiesel/nTuples/2015-07-20/WJetsToLNu_HT-400To600.root") {
  gSystem->Load("pluginTreeWriterTreeWriterAuto.so");

  TChain ch( "TreeWriter/eventTree" );
  ch.AddFile( infile.c_str() );

  double start_time = time(NULL);
  ch.Process( "HistogramProducer.cc+" );
  cout << "Job needed " << 1.*( time(NULL) - start_time)/360 << " min real time." << endl;

}
