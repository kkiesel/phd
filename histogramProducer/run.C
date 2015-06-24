// Shared libray should be produced with compiling
// https://github.com/cms-susy-photon-rwth-1b/TreeWriter
// The lib should be in the CMSSW lib path

//R__LOAD_LIBRARY(TreeParticles)

// start with root -l -q run.C
// Using ROOT6
//void run( string infile="photon_ntuple_mva_mini_38.root") {
void run( string infile="/user/kiesel/nTuples/QCD_HT_1000ToInf_2_nTuple.root") {
  gSystem->Load("pluginTreeWriterTreeWriterAuto.so");

  TChain ch( "TreeWriter/eventTree" );
  ch.AddFile( infile.c_str() );
  //ch.AddFile( "/user/kiesel/nTuples/GJets_HT-400to600_nTuple.root" );
  //ch.AddFile( "/user/kiesel/nTuples/GJets_HT-600toInf_nTuple.root" );
  //ch.AddFile( "/user/kiesel/nTuples/ZJetsToNuNu_HT-400to600_nTuple.root" );
  //ch.AddFile( "/user/kiesel/nTuples/ZJetsToNuNu_HT-600toInf_nTuple.root" );

  double start_time = time(NULL);
  ch.Process( "HistogramProducer.cc+" );
  cout << "Job needed " << 1.*( time(NULL) - start_time)/360 << " min real time." << endl;

}
