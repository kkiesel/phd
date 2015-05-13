#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TSelector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"


#include "TreeParticles.hpp"


/*
struct BaseHistograms {
  TH1F h_met( "", ";E^{miss}_{T} (GeV)", 100, 0, 1000 );
  TH1F h_mt_g_met( "", ";m_{T}(p_{T},E^{miss}_{T}) (GeV)", 100, 0, 1000 );
  TH1F h_metAndL( "", ";E^{miss}_{T} #vec{+} l (GeV)", 100, 0, 1000 );

  TH1F h_ht( "", ";H_{T}", 200, 0, 2000 );
  TH1F h_ht_g( "", ";#gammaH_{T}", 200, 0, 2000 );
  TH1F h_st( "", ";S_{T}", 200, 0, 2000 );

  TH1F h_g_pt( "", ";p_{T} (GeV)", 100, 0, 1500 );
  TH1F h_g_eta( "", ";#eta", 100, -3, 3 );

  TH1F h_j1_pt( "", ";p_{T}^{1.jet} (GeV)", 100, 0, 1500 );
  TH1F h_j1_eta( "", ";#eta^{1.jet}", 100, -3, 3 );
  TH1F h_j2_pt( "", ";p_{T}^{2.jet} (GeV)", 100, 0, 1500 );
  TH1F h_j2_eta( "", ";#eta^{2.jet}", 100, -3, 3 );
  TH1F h_j3_pt( "", ";p_{T}^{3.jet} (GeV)", 100, 0, 1500 );
  TH1F h_j3_eta( "", ";#eta^{3.jet}", 100, -3, 3 );

  TH1F h_bj1_pt( "", ";p_{T}^{1.b-jet} (GeV)", 100, 0, 1500 );
  TH1F h_bj1_eta( "", ";#eta^{1.b-jet}", 100, -3, 3 );
  TH1F h_bj2_pt( "", ";p_{T}^{2.b-jet} (GeV)", 100, 0, 1500 );
  TH1F h_bj2_eta( "", ";#eta^{2.b-jet}", 100, -3, 3 );
  TH1F h_bj3_pt( "", ";p_{T}^{3.b-jet} (GeV)", 100, 0, 1500 );
  TH1F h_bj3_eta( "", ";#eta^{3.b-jet}", 100, -3, 3 );


  TH1F h_dphi_met_g( "", ";#Delta#phi(E_{T}^{miss},#gamma)", 100, -3.5, 3.5 );
  TH1F h_dphi_met_j1( "", ";#Delta#phi(E_{T}^{miss},1.jet)", 100, -3.5, 3.5 );
  TH1F h_dphi_met_j2( "", ";#Delta#phi(E_{T}^{miss},2.jet)", 100, -3.5, 3.5 );

  // multiplicities
  TH1F h_n_vertex( "", ";Vertices", 61, -0.5, 60.5 );
  TH1F h_n_photon( "", ";Photons", 4, -0.5, 3.5 );
  TH1F h_n_jet( "", ";Jets", 13, -0.5, 12.5 );
  TH1F h_n_bjet( "", ";B Jets", 13, -0.5, 12.5 );
  TH1F h_n_electron( "", ";Electrons", 4, -0.5, 3.5 );
  TH1F h_n_muon( "", ";Muons", 4, -0.5, 3.5 );

  // jet matching
  TH2F h2_match_jet_photon( "", ";#Delta R;p_{T}^{jet}/p_{T}^{#gamma}", 100, 0, 0.5, 100, 0, 4 );
  TH2F h2_match_jet_electron( "", ";#Delta R;p_{T}^{jet}/p_{T}^{e}", 100, 0, 0.5, 100, 0, 4 );
  TH2F h2_match_jet_muon( "", ";#Delta R;p_{T}^{jet}/p_{T}^{#mu}", 100, 0, 0.5, 100, 0, 4 );

  // profile
  TProfile profile_g_pt_met( "", ";E^{miss}_{T} (GeV);#LTp_{T}#GT (GeV)", 100, 0, 1000 );
  TProfile profile_ht_met( "", ";E^{miss}_{T} (GeV);#LTH_{T}#GT (GeV)", 100, 0, 1000 );
  TProfile profile_jets_met( "", ";E^{miss}_{T} (GeV);#LTJets#GT", 100, 0, 1000 );

};
*/


class TreeReader : public TSelector {

 public:

  TreeReader();
  virtual ~TreeReader() { }

  virtual void    Init(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual Bool_t  Process(Long64_t entry);
  virtual void    Terminate();
  virtual Int_t   Version() const { return 2; }

  //ClassDef(TreeReader,0);
 private:
  TTreeReader fReader;
  TTreeReaderArray<tree::Photon> photons;
  TTreeReaderArray<tree::Jet> jets;
  TTreeReaderArray<tree::Electron> electrons;
  TTreeReaderArray<tree::Muon> muons;
  TTreeReaderArray<tree::Particle> genJets;

  TTreeReaderValue<tree::MET> met;
  TTreeReaderValue<Int_t> nGoodVertices;
  TTreeReaderValue<Float_t> w;

  TH1F h_met;
  TH1F h_ht;


};

TreeReader::TreeReader():
  photons( fReader, "photons" ),
  jets( fReader, "jets" ),
  electrons( fReader, "electrons" ),
  muons( fReader, "muons" ),
  genJets( fReader, "genJets" ),
  met( fReader, "met" ),
  nGoodVertices( fReader, "nGoodVertices" ),
  w( fReader, "pu_weight" )
{
}

void TreeReader::Init(TTree *tree)
{
  fReader.SetTree(tree);
}

void TreeReader::SlaveBegin(TTree *tree)
{
  h_met = TH1F( "h_met", ";E^{miss}_{T} (GeV)", 100, 0, 1000 );
  h_ht = TH1F( "h_ht", ";H_{T} (GeV)", 200, 0, 2000 );
}

Bool_t TreeReader::Process(Long64_t entry)
{
  fReader.SetEntry(entry);
  h_met.Fill( met->p.Pt() );
  float ht = 0;
  for( auto jet : jets )
    ht += jet.p.Pt();
  h_ht.Fill( ht );

/*
  if( entry > 3 ) return true;
  printf( "n photons = %i\n", (int)photons.GetSize() );
  for( auto photon : photons ) {
    printf( "pt = %f\n", photon.p.Pt() );
  }
*/

  return kTRUE;
}

void TreeReader::Terminate()
{
  TFile file("out.root", "recreate");
  h_met.Write();
  h_ht.Write();
}
