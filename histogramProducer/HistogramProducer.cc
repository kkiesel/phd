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


///////////////////////////////////////////////////////////////////////////////
// User Functions
///////////////////////////////////////////////////////////////////////////////

map<string,int> bTaggingWorkingPoints8TeV = {
  { "CSVL", 0.244 },
  { "CSVM", 0.679 },
  { "CSVT", 0.989 }
};

pair<TVector3,TVector3> megajets( const vector<TVector3>& jets ) {
  // code from https://twiki.cern.ch/twiki/bin/view/CMSPublic/RazorLikelihoodHowTo

  TVector3 j1, j2;
  int N_comb = (int) pow( 2, jets.size() );

  double M_min = numeric_limits<double>::max();
  int j_count;
  for(int i = 1; i < N_comb-1; i++){
    TVector3 j_temp1, j_temp2;
    int itemp = i;
    j_count = N_comb/2;
    int count = 0;
    while(j_count > 0){
      if(itemp/j_count == 1){
        j_temp1 += jets[count];
      } else {
        j_temp2 += jets[count];
      }
      itemp -= j_count*(itemp/j_count);
      j_count /= 2;
      count++;
    }
    double M_temp = j_temp1.Mag2()+j_temp2.Mag2();
    // smallest mass
    if(M_temp < M_min){
      M_min = M_temp;
      j1 = j_temp1;
      j2 = j_temp2;
    }
  }

  if(j2.Pt() > j1.Pt()) return pair<TVector3,TVector3>(j2,j1);
  else                  return pair<TVector3,TVector3>(j1,j2);
}


pair<float,float> razorVariables( const pair<TVector3,TVector3>& megajets, const TVector3& met ) {
  float mr2 = pow( megajets.first.Mag() + megajets.second.Mag(), 2 ) - pow( megajets.first.Z() + megajets.second.Z(), 2 );
  float mrt2 = 0.5 * ( met.Mag()*( megajets.first.Pt() + megajets.second.Pt() ) - met.XYvector() * ( megajets.first.XYvector() + megajets.second.XYvector() ) );
  float r2 = mrt2 / mr2;
  return pair<float,float>(sqrt(mr2),r2);
}

string getOutputFilename( string inputFileName ) {

  // Converts "/path/to/ntuple/QCD_nTuple.root" to "QCD_hists.root"

  auto startPos = inputFileName.rfind("/");
  auto endPos = inputFileName.find("nTuple.root");
  string outputName = "out.root";
  if( endPos != string::npos ) {
    outputName = inputFileName.substr( startPos+1, endPos-startPos-1 ) + "hists.root";
  }
  return outputName;

}



///////////////////////////////////////////////////////////////////////////////
// User Functions
///////////////////////////////////////////////////////////////////////////////



class BaseHistograms;

class HistogramProducer : public TSelector {
 public:

  HistogramProducer();
  virtual ~HistogramProducer() { }

  virtual void    Init(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual Bool_t  Process(Long64_t entry);
  virtual void    Terminate();
  virtual Int_t   Version() const { return 2; }

  void resetSelection();
  void defaultSelection();

  TTreeReader fReader;
  TTreeReaderArray<tree::Photon> photons;
  TTreeReaderArray<tree::Jet> jets;
  TTreeReaderArray<tree::Electron> electrons;
  TTreeReaderArray<tree::Muon> muons;
  TTreeReaderArray<tree::Particle> genJets;

  TTreeReaderValue<tree::MET> met;
  TTreeReaderValue<Int_t> nGoodVertices;
  TTreeReaderValue<Float_t> w;

  vector<tree::Photon*> selPhotons;
  vector<tree::Jet*> selJets;
  vector<tree::Jet*> selBJets;
  vector<tree::Electron*> selElectrons;
  vector<tree::Muon*> selMuons;

  map<string,BaseHistograms*> hMap;
  map<string,TH2F> h2;
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

    h2["h2_razorPlane"] = TH2F( "", ";M_{R} (GeV); R^{2}", 100, 0, 2000, 100, 0, .5 );

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

  void fill( HistogramProducer& tr ) {
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


HistogramProducer::HistogramProducer():
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

void HistogramProducer::Init(TTree *tree)
{
  fReader.SetTree(tree);
}

void HistogramProducer::SlaveBegin(TTree *tree)
{
  std::vector<string> strs = { "base", "loose", "loose_met200", "looseElectron", "loose_genPhoton", "tightPhoton" };
  for( auto& v : strs )
    hMap[v] = new BaseHistograms();
  h2["h2_razorPlane"] = TH2F( "", ";M_{R} (GeV); R^{2}", 100, 0, 2000, 100, 0, .5 );
}

void HistogramProducer::defaultSelection()
{
  for( auto& mu : muons ) {
    if( mu.p.Pt() < 15 ) continue;
    selMuons.push_back( &mu );
  }
  for( auto& el : electrons ) {
    if( !el.isLoose || el.p.Pt() < 15 ) continue;
    selElectrons.push_back( &el );
  }
  for( auto& jet : jets ) {
    if( !jet.isLoose || jet.p.Pt() < 40 || abs(jet.p.Eta()) > 3 ) continue;

    for( auto& p: selPhotons   ) { if( p->p.DeltaR( jet.p ) < .4 ) goto closeToOther; }
    for( auto& p: selElectrons ) { if( p->p.DeltaR( jet.p ) < .4 ) goto closeToOther; }
    for( auto& p: selMuons     ) { if( p->p.DeltaR( jet.p ) < .4 ) goto closeToOther; }

    selJets.push_back( &jet );
    if( jet.bDiscriminator > bTaggingWorkingPoints8TeV["CSVL"] )
      selBJets.push_back( &jet );
    closeToOther:;
  }
}

Bool_t HistogramProducer::Process(Long64_t entry)
{
  resetSelection();
//  if( entry > 3 ) return true;
  if(!( entry%10000 )) printf( "\r%lli / %lli", entry, fReader.GetEntries(false) );
  fReader.SetEntry(entry);

  for( auto& photon : photons ) {
    selPhotons.push_back( &photon );
  }
  for( auto& jet : jets ) {
    selJets.push_back( &jet );
    selBJets.push_back( &jet );
  }
  for( auto& mu : muons ) {
    selMuons.push_back( &mu );
  }
  for( auto& el : electrons ) {
    selElectrons.push_back( &el );
  }
  hMap["base"]->fill( *this );


  float ht = 0;
  for( auto& jet : jets ){
    if( jet.p.Pt() > 40 && abs(jet.p.Eta()) < 3 ) {
      ht += jet.p.Pt();
    }
  }

  // New selection ============================================================
  resetSelection();

  for( auto& photon : photons ) {
    if( !photon.isLoose || photon.p.Pt() < 100 || abs( photon.p.Eta() ) > 1.4442 || photon.hasPixelSeed ) continue;
    selPhotons.push_back( &photon );
  }
  defaultSelection();

  if( selPhotons.size() && ht > 600 ) {
    hMap["loose"]->fill( *this );
    if( met->p.Pt() > 200 )
      hMap["loose_met200"]->fill( *this );
    if( selPhotons[0]->isTrueAlternative )
      hMap["loose_genPhoton"]->fill( *this );
  }

  // extra stuff for razor
  if( selPhotons.size() && ht > 600 ) {
    vector<TVector3> js;
    for( auto& p : selJets ) js.push_back( p->p );
    for( auto& p : selPhotons ) js.push_back( p->p );

    auto razorPair = razorVariables( megajets( js ), met->p );
    h2["h2_razorPlane"].Fill( razorPair.first, razorPair.second, *w );
  }


  // electron sample
  resetSelection();

  for( auto& photon : photons ) {
    if( !photon.isLoose || photon.p.Pt() < 100 || abs( photon.p.Eta() ) > 1.4442 || !photon.hasPixelSeed ) continue;
    selPhotons.push_back( &photon );
  }
  defaultSelection();

  if( selPhotons.size() && ht > 600 ) {
     hMap["looseElectron"]->fill( *this );
  }


  // New photon id ============================================================
  selPhotons.clear();

  for( auto& photon : photons ) {
    if( !photon.isTight || photon.p.Pt() < 100 || abs( photon.p.Eta() ) > 1.4442 || photon.hasPixelSeed ) continue;
    selPhotons.push_back( &photon );
  }
  defaultSelection();

  if( selPhotons.size() && ht > 600 ) {
    hMap["tightPhoton"]->fill( *this );
  }

  return kTRUE;
}

void HistogramProducer::Terminate()
{
  cout << endl;

  auto outputName = getOutputFilename( fReader.GetTree()->GetCurrentFile()->GetName() );
  cout << "Writing to output file " << outputName << endl;

  // save all defined histograms to file
  TFile file( outputName.c_str(), "RECREATE");
  for( auto& hMapIt : hMap )
    hMapIt.second->save( string("_")+hMapIt.first );

  for( auto& mapIt : h2 )
    mapIt.second.Write( mapIt.first.c_str(), TObject::kWriteDelete );

  fReader.GetTree()->GetCurrentFile()->Get("TreeWriter/hCutFlow")->Write();

}

void HistogramProducer::resetSelection() {
  selPhotons.clear();
  selJets.clear();
  selBJets.clear();
  selElectrons.clear();
  selMuons.clear();
}
