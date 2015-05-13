
void run() {

  gSystem->Load("pluginTreeWriterTreeWriterAuto.so");

  TChain ch( "TreeWriter/eventTree" );
  ch.AddFile( "/user/kiesel/nTuples/GJets_HT-400to600_nTuple.root" );
  //ch.AddFile( "photon_ntuple_mva_mini_38.root" );
  ch.Process( "TreeReader.cc+" );

}
