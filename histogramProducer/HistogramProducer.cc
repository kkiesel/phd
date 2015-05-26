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
 private:
  static map<string,int> bTaggingWorkingPoints8TeV;
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
  BaseHistograms* looseCuts;
  BaseHistograms* tightPhoton;
};

map<string,int> TreeReader::bTaggingWorkingPoints8TeV = {
  { "CSVL", 0.244 },
  { "CSVM", 0.679 },
  { "CSVT", 0.989 }
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

    h["h_ht"] = TH1F( "", ";H_{T}", 200, 0, 2500 );
    h["h_ht_g"] = TH1F( "", ";#gammaH_{T}", 200, 0, 2500 );
    h["h_st"] = TH1F( "", ";S_{T}", 200, 0, 2500 );

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
    h["h_n_jet"] = TH1F( "", ";Jets", 21, -0.5, 20.5 );
    h["h_n_bjet"] = TH1F( "", ";B Jets", 21, -0.5, 20.5 );
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

    // isolation dependencies
    prof["profile_cIso_met"] = TProfile( "", ";E^{miss}_{T} (GeV);#LTI_{#pi}#GT (GeV)", 100, 0, 1000 );
    prof["profile_nIso_met"] = TProfile( "", ";E^{miss}_{T} (GeV);#LTI_{n}#GT (GeV)", 100, 0, 1000 );
    prof["profile_pIso_met"] = TProfile( "", ";E^{miss}_{T} (GeV);#LTI_{#gamma}#GT (GeV)", 100, 0, 1000 );
    prof["profile_sie_met"] = TProfile( "", ";E^{miss}_{T} (GeV);#LT#sigma_{i#etai#eta}#GT", 100, 0, 1000 );
    prof["profile_hOverE_met"] = TProfile( "", ";E^{miss}_{T} (GeV);#LTH/E#GT", 100, 0, 1000 );
    prof["profile_mva_met"] = TProfile( "", ";E^{miss}_{T} (GeV);#LTy_{MVA}#GT", 100, 0, 1000 );

    prof["profile_cIso_met_gloose"] = TProfile( "", ";E^{miss}_{T} (GeV);#LTI_{#pi}#GT (GeV)", 100, 0, 1000 );
    prof["profile_nIso_met_gloose"] = TProfile( "", ";E^{miss}_{T} (GeV);#LTI_{n}#GT (GeV)", 100, 0, 1000 );
    prof["profile_pIso_met_gloose"] = TProfile( "", ";E^{miss}_{T} (GeV);#LTI_{#gamma}#GT (GeV)", 100, 0, 1000 );
    prof["profile_sie_met_gloose"] = TProfile( "", ";E^{miss}_{T} (GeV);#LT#sigma_{i#etai#eta}#GT", 100, 0, 1000 );
    prof["profile_hOverE_met_gloose"] = TProfile( "", ";E^{miss}_{T} (GeV);#LTH/E#GT", 100, 0, 1000 );
    prof["profile_mva_met_gloose"] = TProfile( "", ";E^{miss}_{T} (GeV);#LTy_{MVA}#GT", 100, 0, 1000 );

  }

  void fill( TreeReader& tr ) {
    float w = *(tr.w);

    auto metAndL = tr.met->p;
    for( auto l : tr.selElectrons )
      metAndL = metAndL + l->p;
    for( auto& l : tr.selMuons )
      metAndL += l->p;

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
      h["h_g_pt"].Fill( tr.selPhotons[0]->p.Pt(), w );
      h["h_g_eta"].Fill( tr.selPhotons[0]->p.Eta(), w );
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
        h2["h2_match_jet_photon"].Fill( jet->p.DeltaR( p->p ), jet->p.Pt()/p->p.Pt(), w );
      for( auto& p : tr.selElectrons )
        h2["h2_match_jet_electron"].Fill( jet->p.DeltaR( p->p ), jet->p.Pt()/p->p.Pt(), w );
      for( auto& p : tr.selMuons )
        h2["h2_match_jet_muon"].Fill( jet->p.DeltaR( p->p ), jet->p.Pt()/p->p.Pt(), w );
    }

    // profile
    if( tr.selPhotons.size() > 0 )
      prof["profile_g_pt_met"].Fill( tr.met->p.Pt(), tr.selPhotons[0]->p.Pt(), w );
    prof["profile_ht_met"].Fill( tr.met->p.Pt(), ht, w );
    prof["profile_jets_met"].Fill( tr.met->p.Pt(), tr.selJets.size(), w );

    // isolations
    for( auto& p : tr.selPhotons ) {
      prof["profile_cIso_met"].Fill( tr.met->p.Pt(), p->isoChargedHadronsWithEA, w );
      prof["profile_nIso_met"].Fill( tr.met->p.Pt(), p->isoNeutralHadronsWithEA, w );
      prof["profile_pIso_met"].Fill( tr.met->p.Pt(), p->isoPhotonsWithEA, w );
      prof["profile_sie_met"].Fill( tr.met->p.Pt(), p->full5x5_sigmaIetaIeta, w );
      prof["profile_hOverE_met"].Fill( tr.met->p.Pt(), p->hOverE, w );
      prof["profile_mva_met"].Fill( tr.met->p.Pt(), p->mvaValue, w );
    }
    for( auto& p : tr.selPhotons ) {
      if( p->p.Pt() < 50 || abs(p->p.Eta()) > 1.4442 ) continue;
      if(
        p->hOverE < 0.028
        && p->full5x5_sigmaIetaIeta < 0.0107
        && p->isoNeutralHadronsWithEA < 7.23 + exp(0.0028*p->p.Pt()+0.5408)
        && p->isoPhotonsWithEA < 2.11 + 0.0014*p->p.Pt()
      )
      prof["profile_cIso_met_gloose"].Fill( tr.met->p.Pt(), p->isoChargedHadronsWithEA, w );
      if(
        p->hOverE < 0.028
        && p->full5x5_sigmaIetaIeta < 0.0107
        && p->isoChargedHadronsWithEA < 2.67
        && p->isoPhotonsWithEA < 2.11 + 0.0014*p->p.Pt()
      )
      prof["profile_nIso_met_gloose"].Fill( tr.met->p.Pt(), p->isoNeutralHadronsWithEA, w );
      if(
        p->hOverE < 0.028
        && p->full5x5_sigmaIetaIeta < 0.0107
        && p->isoChargedHadronsWithEA < 2.67
        && p->isoNeutralHadronsWithEA < 7.23 + exp(0.0028*p->p.Pt()+0.5408)
      )
      prof["profile_pIso_met_gloose"].Fill( tr.met->p.Pt(), p->isoPhotonsWithEA, w );
      if(
        p->hOverE < 0.028
        && p->isoChargedHadronsWithEA < 2.67
        && p->isoNeutralHadronsWithEA < 7.23 + exp(0.0028*p->p.Pt()+0.5408)
        && p->isoPhotonsWithEA < 2.11 + 0.0014*p->p.Pt()
      )
      prof["profile_sie_met_gloose"].Fill( tr.met->p.Pt(), p->full5x5_sigmaIetaIeta, w );
      if(
        p->full5x5_sigmaIetaIeta < 0.0107
        && p->isoChargedHadronsWithEA < 2.67
        && p->isoNeutralHadronsWithEA < 7.23 + exp(0.0028*p->p.Pt()+0.5408)
        && p->isoPhotonsWithEA < 2.11 + 0.0014*p->p.Pt()
      )
      prof["profile_hOverE_met_gloose"].Fill( tr.met->p.Pt(), p->hOverE, w );
      if(
        p->hOverE < 0.028
        && p->full5x5_sigmaIetaIeta < 0.0107
        && p->isoChargedHadronsWithEA < 2.67
        && p->isoNeutralHadronsWithEA < 7.23 + exp(0.0028*p->p.Pt()+0.5408)
        && p->isoPhotonsWithEA < 2.11 + 0.0014*p->p.Pt()
      )
      prof["profile_mva_met_gloose"].Fill( tr.met->p.Pt(), p->mvaValue, w );
    }


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
  looseCuts = new BaseHistograms();
  tightPhoton = new BaseHistograms();
}

Bool_t TreeReader::Process(Long64_t entry)
{
  resetSelection();
//  if( entry > 3 ) return true;
//  if( entry > 20000 ) return true;
  if(!( entry%10000 )) printf( "\r%lli / %lli", entry, fReader.GetEntries(false) );
  fReader.SetEntry(entry);

  if( photons.GetSize() ) {
  for( auto& photon : photons ) {
    selPhotons.push_back( &photon );
  }
  }
  if( jets.GetSize() ) {
  for( auto& jet : jets ) {
    selJets.push_back( &jet );
  }
  }
  if( jets.GetSize() ) {
  for( auto& jet : jets ) {
    selBJets.push_back( &jet );
  }
  }
  if( muons.GetSize() ) {
  for( auto& mu : muons ) {
    selMuons.push_back( &mu );
  }
  }
  if( electrons.GetSize() ) {
  for( auto& el : electrons ) {
    selElectrons.push_back( &el );
  }
  }
  baseHistograms->fill( *this );

  // New selection ============================================================
  resetSelection();

  if( photons.GetSize() ) {
  for( auto& photon : photons ) {
    if( !photon.isLoose || photon.p.Pt() < 50 ) continue;
    selPhotons.push_back( &photon );
  }
  }
  if( jets.GetSize() ) {
  for( auto& jet : jets ) {
    if( !jet.isLoose || jet.p.Pt() < 40 || abs(jet.p.Eta()) > 3 ) continue;
    selJets.push_back( &jet );
  }
  }
  if( jets.GetSize() ) {
  for( auto& jet : jets ) {
    if( !jet.isLoose || jet.p.Pt() < 40 || abs(jet.p.Eta()) > 2.5 || jet.bDiscriminator < bTaggingWorkingPoints8TeV["CSVL"] ) continue;
    selBJets.push_back( &jet );
  }
  }
  if( muons.GetSize() ) {
  for( auto& mu : muons ) {
    if( mu.p.Pt() < 15 ) continue;
    selMuons.push_back( &mu );
  }
  }
  if( electrons.GetSize() ) {
  for( auto& el : electrons ) {
    if( !el.isLoose || el.p.Pt() < 15 ) continue;
    selElectrons.push_back( &el );
  }
  }
  if( selPhotons.size() )
    looseCuts->fill( *this );

  // New photon id ============================================================
  selPhotons.clear();

  if( photons.GetSize() ) {
  for( auto& photon : photons ) {
    if( !photon.isTight || photon.p.Pt() < 50 ) continue;
    selPhotons.push_back( &photon );
  }
  }
  if( selPhotons.size() )
    tightPhoton->fill( *this );



  return kTRUE;
}

void TreeReader::Terminate()
{
  cout << endl;
  string originalName = fReader.GetTree()->GetCurrentFile()->GetName();
  auto startPos = originalName.rfind("/");
  auto endPos = originalName.find("nTuple.root");
  string outputName = "out.root";
  if( endPos != string::npos ) {
    outputName = originalName.substr( startPos+1, endPos-startPos-1 ) + "hists.root";
  }
  cout << "Writing to output file " << outputName << endl;
  TFile file( outputName.c_str(), "RECREATE");
  baseHistograms->save();
  looseCuts->save("_loose");
  tightPhoton->save("_tightPhoton");
}

void TreeReader::resetSelection() {
  selPhotons.clear();
  selJets.clear();
  selBJets.clear();
  selElectrons.clear();
  selMuons.clear();
}
