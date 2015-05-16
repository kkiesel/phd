
// start with root -l -q run.C
// Using ROOT6
void run() {

  // Shared libray should be produced with compiling
  // https://github.com/cms-susy-photon-rwth-1b/TreeWriter
  // The lib should be in the CMSSW lib path
  gSystem->Load("pluginTreeWriterTreeWriterAuto.so");

  TChain ch( "TreeWriter/eventTree" );
  ch.AddFile( "/user/kiesel/nTuples/GJets_HT-400to600_nTuple.root" );
  //ch.AddFile( "photon_ntuple_mva_mini_38.root" );
  ch.Process( "TreeReader.cc+" );

}
