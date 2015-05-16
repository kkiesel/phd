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

class BaseHistograms;

class TreeReader : public TSelector {

 public:

  TreeReader();
  virtual ~TreeReader() { }

  virtual void    Init(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual Bool_t  Process(Long64_t entry);
  virtual void    Terminate();
  virtual Int_t   Version() const { return 2; }

  void resetSelection();

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

  vector<tree::Photon*> selPhotons;
  vector<tree::Jet*> selJets;
  vector<tree::Jet*> selBJets;
  vector<tree::Electron*> selElectrons;
  vector<tree::Muon*> selMuons;

  BaseHistograms* baseHistograms;
};

class BaseHistograms {
 public:

  map<string,TH1F> h;
  map<string,TH2F> h2;
  map<string,TProfile> prof;

  BaseHistograms() {
    h["h_met"] = TH1F( "", ";E^{miss}_{T} (GeV)", 100, 0, 1000 );
    h["h_mt_g_met"] = TH1F( "", ";m_{T}(p_{T},E^{miss}_{T}) (GeV)", 100, 0, 1000 );
    h["h_metAndL"] = TH1F( "", ";E^{miss}_{T} #vec{+} l (GeV)", 100, 0, 1000 );

    h["h_ht"] = TH1F( "", ";H_{T}", 200, 0, 2000 );
    h["h_ht_g"] = TH1F( "", ";#gammaH_{T}", 200, 0, 2000 );
    h["h_st"] = TH1F( "", ";S_{T}", 200, 0, 2000 );

    h["h_g_pt"] = TH1F( "", ";p_{T} (GeV)", 100, 0, 1500 );
    h["h_g_eta"] = TH1F( "", ";#eta", 100, -3, 3 );

    h["h_j1_pt"] = TH1F( "", ";p_{T}^{1.jet} (GeV)", 100, 0, 1500 );
    h["h_j1_eta"] = TH1F( "", ";#eta^{1.jet}", 100, -3, 3 );
    h["h_j2_pt"] = TH1F( "", ";p_{T}^{2.jet} (GeV)", 100, 0, 1500 );
    h["h_j2_eta"] = TH1F( "", ";#eta^{2.jet}", 100, -3, 3 );
    h["h_j3_pt"] = TH1F( "", ";p_{T}^{3.jet} (GeV)", 100, 0, 1500 );
    h["h_j3_eta"] = TH1F( "", ";#eta^{3.jet}", 100, -3, 3 );

    h["h_bj1_pt"] = TH1F( "", ";p_{T}^{1.b-jet} (GeV)", 100, 0, 1500 );
    h["h_bj1_eta"] = TH1F( "", ";#eta^{1.b-jet}", 100, -3, 3 );
    h["h_bj2_pt"] = TH1F( "", ";p_{T}^{2.b-jet} (GeV)", 100, 0, 1500 );
    h["h_bj2_eta"] = TH1F( "", ";#eta^{2.b-jet}", 100, -3, 3 );
    h["h_bj3_pt"] = TH1F( "", ";p_{T}^{3.b-jet} (GeV)", 100, 0, 1500 );
    h["h_bj3_eta"] = TH1F( "", ";#eta^{3.b-jet}", 100, -3, 3 );


    h["h_dphi_met_g"] = TH1F( "", ";#Delta#phi(E_{T}^{miss},#gamma)", 100, -3.5, 3.5 );
    h["h_dphi_met_j1"] = TH1F( "", ";#Delta#phi(E_{T}^{miss},1.jet)", 100, -3.5, 3.5 );
    h["h_dphi_met_j2"] = TH1F( "", ";#Delta#phi(E_{T}^{miss},2.jet)", 100, -3.5, 3.5 );

    // multiplicities
    h["h_n_vertex"] = TH1F( "", ";Vertices", 61, -0.5, 60.5 );
    h["h_n_photon"] = TH1F( "", ";Photons", 6, -0.5, 5.5 );
    h["h_n_jet"] = TH1F( "", ";Jets", 13, -0.5, 12.5 );
    h["h_n_bjet"] = TH1F( "", ";B Jets", 13, -0.5, 12.5 );
    h["h_n_electron"] = TH1F( "", ";Electrons", 4, -0.5, 3.5 );
    h["h_n_muon"] = TH1F( "", ";Muons", 4, -0.5, 3.5 );

    // matching
    h2["h2_match_jet_photon"] = TH2F( "", ";#Delta R;p_{T}^{jet}/p_{T}^{#gamma}", 100, 0, 0.5, 100, 0, 4 );
    h2["h2_match_jet_electron"] = TH2F( "", ";#Delta R;p_{T}^{jet}/p_{T}^{e}", 100, 0, 0.5, 100, 0, 4 );
    h2["h2_match_jet_muon"] = TH2F( "", ";#Delta R;p_{T}^{jet}/p_{T}^{#mu}", 100, 0, 0.5, 100, 0, 4 );

    // profile
    prof["profile_g_pt_met"] = TProfile( "", ";E^{miss}_{T} (GeV);#LTp_{T}#GT (GeV)", 100, 0, 1000 );
    prof["profile_ht_met"] = TProfile( "", ";E^{miss}_{T} (GeV);#LTH_{T}#GT (GeV)", 100, 0, 1000 );
    prof["profile_jets_met"] = TProfile( "", ";E^{miss}_{T} (GeV);#LTJets#GT", 100, 0, 1000 );
  }

  void fill( TreeReader& tr ) {
    float w = *(tr.w);

    auto metAndL = tr.met->p;
//    printf("met = %f\n", tr.met->p.Pt() );
//    printf("n ele = %i\n", (int)(tr.selElectrons.size()) );
//    if( tr.selElectrons.size() > 0 )
//      printf("lepton pt = = %f\n", tr.selElectrons.at(0)->p.Pt() );
//    for( auto l : tr.selElectrons )
//      metAndL = metAndL + l->p;
//    for( auto& l : tr.selMuons )
//      metAndL += l->p;

    float ht = 0;
    for( auto& j : tr.selJets )
      ht += j->p.Pt();

    float ht_g = ht;
    for( auto& p : tr.selPhotons )
      ht_g += p->p.Pt();

    float st = ht_g + tr.met->p.Pt();




    h["h_met"].Fill( tr.met->p.Pt(), w );
    h["h_metAndL"].Fill( metAndL.Pt(), w );

    h["h_ht"].Fill( ht, w );
    h["h_ht_g"].Fill( ht_g, w );
    h["h_st"].Fill( st, w );

    if( tr.selPhotons.size() > 0 ) {
      h["h_mt_g_met"].Fill( (tr.selPhotons[0]->p + tr.met->p).Pt(), w );
      h["h_g_pt"].Fill( tr.met->p.Pt(), w );
      h["h_g_eta"].Fill( tr.met->p.Pt(), w );
    }

    if( tr.selJets.size() > 2 ) {
      h["h_j3_pt"].Fill( tr.selJets[2]->p.Pt(), w );
      h["h_j3_eta"].Fill( tr.selJets[2]->p.Eta(), w );
    }
    if( tr.selJets.size() > 1 ) {
      h["h_j2_pt"].Fill( tr.selJets[1]->p.Pt(), w );
      h["h_j2_eta"].Fill( tr.selJets[1]->p.Eta(), w );
    }
    if( tr.selJets.size() > 0 ) {
      h["h_j1_pt"].Fill( tr.selJets[0]->p.Pt(), w );
      h["h_j1_eta"].Fill( tr.selJets[0]->p.Eta(), w );
    }

    if( tr.selBJets.size() > 2 ) {
      h["h_bj3_pt"].Fill( tr.selBJets[2]->p.Pt(), w );
      h["h_bj3_eta"].Fill( tr.selBJets[2]->p.Eta(), w );
    }
    if( tr.selBJets.size() > 1 ) {
      h["h_bj2_pt"].Fill( tr.selBJets[1]->p.Pt(), w );
      h["h_bj2_eta"].Fill( tr.selBJets[1]->p.Eta(), w );
    }
    if( tr.selBJets.size() > 0 ) {
      h["h_bj1_pt"].Fill( tr.selBJets[0]->p.Pt(), w );
      h["h_bj1_eta"].Fill( tr.selBJets[0]->p.Eta(), w );
    }

    if( tr.selPhotons.size() > 0 )
      h["h_dphi_met_g"].Fill( tr.met->p.DeltaPhi( tr.selPhotons[0]->p ), w );
    if( tr.selJets.size() > 0 )
      h["h_dphi_met_j1"].Fill( tr.met->p.DeltaPhi( tr.selJets[0]->p ), w );
    if( tr.selPhotons.size() > 1 )
      h["h_dphi_met_j2"].Fill( tr.met->p.DeltaPhi( tr.selJets[1]->p ), w );

    h["h_n_vertex"].Fill( *tr.nGoodVertices, w );
    h["h_n_photon"].Fill( tr.selPhotons.size(), w );
    h["h_n_jet"].Fill( tr.selJets.size(), w );
    h["h_n_bjet"].Fill( tr.selBJets.size(), w );
    h["h_n_electron"].Fill( tr.selElectrons.size(), w );
    h["h_n_muon"].Fill( tr.selMuons.size(), w );

    // matching
    for( auto& jet : tr.selJets ) {
      for( auto& p : tr.selPhotons )
        h2["h2_match_jet_photon"].Fill( jet->p.DeltaR( p->p ), tr.met->p.Pt(), w );
      for( auto& p : tr.selElectrons )
        h2["h2_match_jet_electron"].Fill( tr.met->p.Pt(), w );
      for( auto& p : tr.selMuons )
        h2["h2_match_jet_muon"].Fill( tr.met->p.Pt(), w );
    }

    // profile
    if( tr.selPhotons.size() > 0 )
      prof["profile_g_pt_met"].Fill( tr.met->p.Pt(), tr.selPhotons[0]->p.Pt(), w );
    prof["profile_ht_met"].Fill( tr.met->p.Pt(), ht, w );
    prof["profile_jets_met"].Fill( tr.met->p.Pt(), tr.selJets.size(), w );

  } // end fill

  void save( const string& namePrefix="" ) {
    for( auto& mapIt : h )
      mapIt.second.Write( (mapIt.first+namePrefix).c_str(), TObject::kWriteDelete );
    for( auto& mapIt : h2 )
      mapIt.second.Write( (mapIt.first+namePrefix).c_str(), TObject::kWriteDelete );
    for( auto& mapIt : prof )
      mapIt.second.Write( (mapIt.first+namePrefix).c_str(), TObject::kWriteDelete );
  }

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
  baseHistograms = new BaseHistograms();
}

Bool_t TreeReader::Process(Long64_t entry)
{
  if( entry > 3 ) return true;
  if( entry > 20000 ) return true;
  if(!( entry%1000 )) printf( "%lli / %lli\n", entry, fReader.GetEntries(false) );
  fReader.SetEntry(entry);

  for( auto& photon : photons ) {
    selPhotons.push_back( &photon );
  }
  for( auto& jet : jets ) {
    selJets.push_back( &jet );
  }
  for( auto& jet : jets ) {
    selBJets.push_back( &jet );
  }
  for( auto& mu : muons ) {
    selMuons.push_back( &mu );
  }
  printf("original e size %hi\n", electrons.GetSize() );
  for( auto& el : electrons ) {
    printf("filling e\n" );
    selElectrons.push_back( &el );
  }
  printf("sel  e size %hi\n", selElectrons.size() );

  baseHistograms->fill( *this );
  resetSelection();

  return kTRUE;
}

void TreeReader::Terminate()
{
  TFile file("out.root", "update");
  baseHistograms->save();
}

void TreeReader::resetSelection() {
  selPhotons.clear();
  selJets.clear();
  selBJets.clear();
  selElectrons.clear();
  selMuons.clear();
}
