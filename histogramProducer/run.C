// Shared libray should be produced with compiling
// https://github.com/cms-susy-photon-rwth-1b/TreeWriter
// The lib should be in the CMSSW lib path

R__LOAD_LIBRARY(TreeParticles)

// start with root -l -q run.C
// Using ROOT6
void run( string infile="photon_ntuple_mva_mini_38.root") {

  TChain ch( "TreeWriter/eventTree" );
  ch.AddFile( infile.c_str() );
  //ch.AddFile( "/user/kiesel/nTuples/GJets_HT-400to600_nTuple.root" );
  //ch.AddFile( "/user/kiesel/nTuples/GJets_HT-600toInf_nTuple.root" );
  //ch.AddFile( "/user/kiesel/nTuples/ZJetsToNuNu_HT-400to600_nTuple.root" );
  //ch.AddFile( "/user/kiesel/nTuples/ZJetsToNuNu_HT-600toInf_nTuple.root" );
  ch.Process( "HistogramProducer.cc+" );

}
