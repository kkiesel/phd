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
#include "TEfficiency.h"

//#include "TreeParticles.hpp"
#include "../../CMSSW/treewriter/CMSSW_7_4_5/src/TreeWriter/TreeWriter/plugins/TreeParticles.hpp"


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
  auto endPos = inputFileName.find(".root");
  string outputName = "out.root";
  if( endPos != string::npos ) {
    outputName = inputFileName.substr( startPos+1, endPos-startPos-1 ) + "_hists.root";
  }
  return outputName;

}



///////////////////////////////////////////////////////////////////////////////
// End User Functions
///////////////////////////////////////////////////////////////////////////////




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

  void initSelection( string const& s );
  void fillSelection( string const& s );

  void initObjects( string const& s );
  void fillObjects( string const& s );

  void initBkgEst( string const& s );
  void fillBkgEst( string const& s );

  TTreeReader fReader;
  TTreeReaderArray<tree::Photon> photons;
  TTreeReaderArray<tree::Jet> jets;
  TTreeReaderArray<tree::Electron> electrons;
  TTreeReaderArray<tree::Muon> muons;
  TTreeReaderArray<tree::Particle> genJets;
  TTreeReaderArray<tree::GenParticle> genParticles;

  TTreeReaderValue<tree::MET> met;
  TTreeReaderValue<Int_t> nGoodVertices;
  TTreeReaderValue<Float_t> w;
  TTreeReaderValue<Bool_t> signalTriggerFired;
  TTreeReaderValue<Bool_t> crossTriggerPhotonFired;
  TTreeReaderValue<Bool_t> crossTriggerHtFired;

  vector<tree::Photon*> selPhotons;
  vector<tree::Jet*> selJets;
  vector<tree::Jet*> selBJets;
  vector<tree::Electron*> selElectrons;
  vector<tree::Muon*> selMuons;

  pair<float,float> mrr2;

  map<string,TH1F> h;
  map<string,TH2F> h2;
  map<string,TProfile> prof;
  map<string,TEfficiency> eff;

};

void HistogramProducer::initSelection( string const& s ) {
  h["h_met__"+s] = TH1F( "", ";E^{miss}_{T} (GeV)", 100, 0, 1000 );
  h["h_mt_g_met__"+s] = TH1F( "", ";m_{T}(p_{T},E^{miss}_{T}) (GeV)", 100, 0, 1000 );

  h["h_ht__"+s] = TH1F( "", ";H_{T}", 200, 0, 2500 );
  h["h_ht_g__"+s] = TH1F( "", ";HE_{T}", 200, 0, 2500 );
  h["h_st__"+s] = TH1F( "", ";S_{T}", 200, 0, 2500 );

  h["h_g_pt__"+s] = TH1F( "", ";p_{T} (GeV)", 100, 0, 1500 );
  h["h_g_eta__"+s] = TH1F( "", ";#eta", 100, -3, 3 );

  // jet
  h["h_j1_pt__"+s] = TH1F( "", ";p_{T}^{1.jet} (GeV)", 100, 0, 1500 );
  h["h_j1_eta__"+s] = TH1F( "", ";#eta^{1.jet}", 100, -3, 3 );
  h["h_j2_pt__"+s] = TH1F( "", ";p_{T}^{2.jet} (GeV)", 100, 0, 1500 );
  h["h_j2_eta__"+s] = TH1F( "", ";#eta^{2.jet}", 100, -3, 3 );
  h["h_j3_pt__"+s] = TH1F( "", ";p_{T}^{3.jet} (GeV)", 100, 0, 1500 );
  h["h_j3_eta__"+s] = TH1F( "", ";#eta^{3.jet}", 100, -3, 3 );

  // b-jets
  h["h_bj1_pt__"+s] = TH1F( "", ";p_{T}^{1.b-jet} (GeV)", 100, 0, 1500 );
  h["h_bj1_eta__"+s] = TH1F( "", ";#eta^{1.b-jet}", 100, -3, 3 );
  h["h_bj2_pt__"+s] = TH1F( "", ";p_{T}^{2.b-jet} (GeV)", 100, 0, 1500 );
  h["h_bj2_eta__"+s] = TH1F( "", ";#eta^{2.b-jet}", 100, -3, 3 );
  h["h_bj3_pt__"+s] = TH1F( "", ";p_{T}^{3.b-jet} (GeV)", 100, 0, 1500 );
  h["h_bj3_eta__"+s] = TH1F( "", ";#eta^{3.b-jet}", 100, -3, 3 );

  // angles
  h["h_dphi_met_g__"+s] = TH1F( "", ";#Delta#phi(E_{T}^{miss},#gamma)", 100, -3.5, 3.5 );
  h["h_dphi_met_j1__"+s] = TH1F( "", ";#Delta#phi(E_{T}^{miss},1.jet)", 100, -3.5, 3.5 );
  h["h_dphi_met_j2__"+s] = TH1F( "", ";#Delta#phi(E_{T}^{miss},2.jet)", 100, -3.5, 3.5 );
  h["h_dphi_met_j3__"+s] = TH1F( "", ";#Delta#phi(E_{T}^{miss},3.jet)", 100, -3.5, 3.5 );
  h["h_dphi_met_recoil__"+s] = TH1F( "", ";#Delta#phi(E_{T}^{miss},#Sigma jet)", 100, -3.5, 3.5 );

  // multiplicities
  h["h_n_vertex__"+s] = TH1F( "", ";Vertices", 61, -0.5, 60.5 );
  h["h_n_photon__"+s] = TH1F( "", ";Photons", 6, -0.5, 5.5 );
  h["h_n_jet__"+s] = TH1F( "", ";Jets", 21, -0.5, 20.5 );
  h["h_n_bjet__"+s] = TH1F( "", ";B Jets", 21, -0.5, 20.5 );
  h["h_n_electron__"+s] = TH1F( "", ";Electrons", 4, -0.5, 3.5 );
  h["h_n_muon__"+s] = TH1F( "", ";Muons", 4, -0.5, 3.5 );

  h2["h2_razorPlane__"+s] = TH2F( "", ";M_{R} (GeV); R^{2}", 100, 0, 2000, 100, 0, .5 );

}

void HistogramProducer::initObjects( string const& s ) {

  eff["h_j_looseID__"+s] = TEfficiency( "", ";loose jet ID", 2, 0., 2. );

  eff["h_gCol_sigmaIetaIeta__"+s] = TEfficiency( "", ";#sigma_{i#etai#eta}", 100, 0, 0.07 );
  eff["h_gCol_hOverE__"+s] = TEfficiency( "", ";H/E", 100, 0, 0.15 );
  eff["h_gCol_r9__"+s] = TEfficiency( "", ";r9", 110, 0, 1.1 );
  eff["h_gCol_cIso__"+s] = TEfficiency( "", ";I_{#pi}", 200, 0, 150 );
  eff["h_gCol_mva__"+s] = TEfficiency( "", ";y_{mva}", 120, -1.1, 1.1 );

  eff["h_g_sigmaIetaIeta__"+s] = TEfficiency( "", ";#sigma_{i#etai#eta}", 100, 0, 0.07 );
  eff["h_g_hOverE__"+s] = TEfficiency( "", ";H/E", 100, 0, 0.15 );
  eff["h_g_r9__"+s] = TEfficiency( "", ";r9", 110, 0, 1.1 );
  eff["h_g_cIso__"+s] = TEfficiency( "", ";I_{#pi}", 200, 0, 150 );
  eff["h_g_nIso__"+s] = TEfficiency( "", ";I_{n}", 200, 0, 150 );
  eff["h_g_pIso__"+s] = TEfficiency( "", ";I_{#gamma}", 200, 0, 150 );


  // matching
  h2["h2_match_jet_photon__"+s] = TH2F( "", ";#Delta R;p_{T}^{jet}/p_{T}^{#gamma}", 100, 0, 0.5, 100, 0, 4 );
  h2["h2_match_jet_electron__"+s] = TH2F( "", ";#Delta R;p_{T}^{jet}/p_{T}^{e}", 100, 0, 0.5, 100, 0, 4 );
  h2["h2_match_jet_muon__"+s] = TH2F( "", ";#Delta R;p_{T}^{jet}/p_{T}^{#mu}", 100, 0, 0.5, 100, 0, 4 );

  h2["h2_match_photon_genElectron__"+s] = TH2F( "", ";#Delta R;p_{T}/p_{T}^{gen e}", 100, 0, 0.5, 100, 0, 4 );
  h2["h2_match_photon_genPhoton__"+s] = TH2F( "", ";#Delta R;p_{T}/p_{T}^{gen #gamma}", 100, 0, 0.5, 100, 0, 4 );

  eff["eff_hlt_pt__"+s] = TEfficiency( "", ";p_{T};#varepsilon", 200, 0, 2000 );
  eff["eff_hlt_ht__"+s] = TEfficiency( "", ";H_{T};#varepsilon", 200, 0, 2000 );

}

void HistogramProducer::initBkgEst( string const& s ) {
  h["h_metAndL__"+s] = TH1F( "", ";E^{miss}_{T} #vec{+} l (GeV)", 100, 0, 1000 );

  // profile
  prof["profile_g_pt_met__"+s] = TProfile( "", ";E^{miss}_{T} (GeV);#LTp_{T}#GT (GeV)", 100, 0, 1000 );
  prof["profile_jets_met__"+s] = TProfile( "", ";E^{miss}_{T} (GeV);#LTJets#GT", 100, 0, 1000 );

  // isolation dependencies
  prof["profile_cIso_met__"+s] = TProfile( "", ";E^{miss}_{T} (GeV);#LTI_{#pi}#GT (GeV)", 100, 0, 1000 );
  prof["profile_nIso_met__"+s] = TProfile( "", ";E^{miss}_{T} (GeV);#LTI_{n}#GT (GeV)", 100, 0, 1000 );
  prof["profile_pIso_met__"+s] = TProfile( "", ";E^{miss}_{T} (GeV);#LTI_{#gamma}#GT (GeV)", 100, 0, 1000 );
  prof["profile_sie_met__"+s] = TProfile( "", ";E^{miss}_{T} (GeV);#LT#sigma_{i#etai#eta}#GT", 100, 0, 1000 );
  prof["profile_hOverE_met__"+s] = TProfile( "", ";E^{miss}_{T} (GeV);#LTH/E#GT", 100, 0, 1000 );
  prof["profile_mva_met__"+s] = TProfile( "", ";E^{miss}_{T} (GeV);#LTy_{MVA}#GT", 100, 0, 1000 );

  prof["profile_cIso_met_gloose__"+s] = TProfile( "", ";E^{miss}_{T} (GeV);#LTI_{#pi}#GT (GeV)", 100, 0, 1000 );
  prof["profile_nIso_met_gloose__"+s] = TProfile( "", ";E^{miss}_{T} (GeV);#LTI_{n}#GT (GeV)", 100, 0, 1000 );
  prof["profile_pIso_met_gloose__"+s] = TProfile( "", ";E^{miss}_{T} (GeV);#LTI_{#gamma}#GT (GeV)", 100, 0, 1000 );
  prof["profile_sie_met_gloose__"+s] = TProfile( "", ";E^{miss}_{T} (GeV);#LT#sigma_{i#etai#eta}#GT", 100, 0, 1000 );
  prof["profile_hOverE_met_gloose__"+s] = TProfile( "", ";E^{miss}_{T} (GeV);#LTH/E#GT", 100, 0, 1000 );
  prof["profile_mva_met_gloose__"+s] = TProfile( "", ";E^{miss}_{T} (GeV);#LTy_{MVA}#GT", 100, 0, 1000 );
}


void HistogramProducer::fillSelection( string const& s ) {

  float ht = 0;
  for( auto& j : selJets )
    ht += j->p.Pt();

  float ht_g = ht;
  for( auto& p : selPhotons )
    ht_g += p->p.Pt();

  float st = ht_g + met->p.Pt();

  TVector3 recoil(0,0,0);
  for( auto& p : selJets )
    recoil += p->p;

  h["h_met__"+s].Fill( met->p.Pt(), *w );

  h["h_ht__"+s].Fill( ht, *w );
  h["h_ht_g__"+s].Fill( ht_g, *w );
  h["h_st__"+s].Fill( st, *w );

  if( selPhotons.size() > 0 ) {
    h["h_mt_g_met__"+s].Fill( (selPhotons[0]->p + met->p).Pt(), *w );
    h["h_g_pt__"+s].Fill( selPhotons[0]->p.Pt(), *w );
    h["h_g_eta__"+s].Fill( selPhotons[0]->p.Eta(), *w );
  }

  if( selJets.size() > 2 ) {
    h["h_j3_pt__"+s].Fill( selJets[2]->p.Pt(), *w );
    h["h_j3_eta__"+s].Fill( selJets[2]->p.Eta(), *w );
  }
  if( selJets.size() > 1 ) {
    h["h_j2_pt__"+s].Fill( selJets[1]->p.Pt(), *w );
    h["h_j2_eta__"+s].Fill( selJets[1]->p.Eta(), *w );
  }
  if( selJets.size() > 0 ) {
    h["h_j1_pt__"+s].Fill( selJets[0]->p.Pt(), *w );
    h["h_j1_eta__"+s].Fill( selJets[0]->p.Eta(), *w );
  }

  if( selBJets.size() > 2 ) {
    h["h_bj3_pt__"+s].Fill( selBJets[2]->p.Pt(), *w );
    h["h_bj3_eta__"+s].Fill( selBJets[2]->p.Eta(), *w );
  }
  if( selBJets.size() > 1 ) {
    h["h_bj2_pt__"+s].Fill( selBJets[1]->p.Pt(), *w );
    h["h_bj2_eta__"+s].Fill( selBJets[1]->p.Eta(), *w );
  }
  if( selBJets.size() > 0 ) {
    h["h_bj1_pt__"+s].Fill( selBJets[0]->p.Pt(), *w );
    h["h_bj1_eta__"+s].Fill( selBJets[0]->p.Eta(), *w );
  }

  if( selPhotons.size() > 0 )
    h["h_dphi_met_g__"+s].Fill( met->p.DeltaPhi( selPhotons[0]->p ), *w );
  if( selJets.size() > 0 )
    h["h_dphi_met_j1__"+s].Fill( met->p.DeltaPhi( selJets[0]->p ), *w );
  if( selJets.size() > 1 )
    h["h_dphi_met_j2__"+s].Fill( met->p.DeltaPhi( selJets[1]->p ), *w );
  if( selJets.size() > 2 )
    h["h_dphi_met_j3__"+s].Fill( met->p.DeltaPhi( selJets[2]->p ), *w );
  h["h_dphi_met_recoil__"+s].Fill( met->p.DeltaPhi( recoil ), *w );

  h["h_n_vertex__"+s].Fill( *nGoodVertices, *w );
  h["h_n_photon__"+s].Fill( selPhotons.size(), *w );
  h["h_n_jet__"+s].Fill( selJets.size(), *w );
  h["h_n_bjet__"+s].Fill( selBJets.size(), *w );
  h["h_n_electron__"+s].Fill( selElectrons.size(), *w );
  h["h_n_muon__"+s].Fill( selMuons.size(), *w );

  if( mrr2.first > 1e-7 ) {
    h2["h2_razorPlane__"+s].Fill( mrr2.first, mrr2.second, *w );
  }

} // end fill

void HistogramProducer::fillObjects( string const& s ) {

  // matching
  for( auto& jet : selJets ) {
    for( auto& p : selPhotons )
      h2["h2_match_jet_photon__"+s].Fill( jet->p.DeltaR( p->p ), jet->p.Pt()/p->p.Pt(), *w );
    for( auto& p : selElectrons )
      h2["h2_match_jet_electron__"+s].Fill( jet->p.DeltaR( p->p ), jet->p.Pt()/p->p.Pt(), *w );
    for( auto& p : selMuons )
      h2["h2_match_jet_muon__"+s].Fill( jet->p.DeltaR( p->p ), jet->p.Pt()/p->p.Pt(), *w );

    h["h_j_looseID__"+s].Fill( float(jet->isLoose) );
  }

  for( auto& photon : selPhotons ) {
    for( auto& p : genParticles ) {
      if( abs(p.pdgId) == 22 )
        h2["h2_match_photon_genPhoton__"+s].Fill( photon->p.DeltaR( p.p ), photon->p.Pt()/p.p.Pt(), *w );
      if( abs(p.pdgId) == 11 )
        h2["h2_match_photon_genElectron__"+s].Fill( photon->p.DeltaR( p.p ), photon->p.Pt()/p.p.Pt(), *w );
    }

    eff["h_gCol_sigmaIetaIeta__"+s].FillWeighted( photon->isTrue, *w, photon->sigmaIetaIeta );
    eff["h_gCol_hOverE__"+s].FillWeighted( photon->isTrue, *w, photon->hOverE );
    eff["h_gCol_r9__"+s].FillWeighted( photon->isTrue, *w, photon->r9 );
    eff["h_gCol_cIso__"+s].FillWeighted( photon->isTrue, *w, photon->isoChargedHadronsEA );
    eff["h_gCol_mva__"+s].FillWeighted( photon->isTrue, *w, photon->mvaValue );

    if( photon->p.Pt() < 90 || abs(photon->p.Eta()) > 1.4442 ) continue;
    if(
      photon->hOverE < 0.028
      && photon->sigmaIetaIeta < 0.0107
      && photon->isoNeutralHadronsEA < 7.23 + exp(0.0028*photon->p.Pt()+0.5408)
      && photon->isoPhotonsEA < 2.11 + 0.0014*photon->p.Pt()
    )
    eff["h_g_cIso__"+s].FillWeighted( photon->isTrue, *w, photon->isoChargedHadronsEA );
    if(
      photon->hOverE < 0.028
      && photon->sigmaIetaIeta < 0.0107
      && photon->isoChargedHadronsEA < 2.67
      && photon->isoPhotonsEA < 2.11 + 0.0014*photon->p.Pt()
    )
    eff["h_g_nIso__"+s].FillWeighted( photon->isTrue, *w, photon->isoNeutralHadronsEA );
    if(
      photon->hOverE < 0.028
      && photon->sigmaIetaIeta < 0.0107
      && photon->isoChargedHadronsEA < 2.67
      && photon->isoNeutralHadronsEA < 7.23 + exp(0.0028*photon->p.Pt()+0.5408)
    )
    eff["h_g_nIso__"+s].FillWeighted( photon->isTrue, *w, photon->isoPhotonsEA );
    if(
      photon->hOverE < 0.028
      && photon->isoChargedHadronsEA < 2.67
      && photon->isoNeutralHadronsEA < 7.23 + exp(0.0028*photon->p.Pt()+0.5408)
      && photon->isoPhotonsEA < 2.11 + 0.0014*photon->p.Pt()
    )
    eff["h_g_sigmaIetaIeta__"+s].FillWeighted( photon->isTrue, *w, photon->sigmaIetaIeta );
    if(
      photon->sigmaIetaIeta < 0.0107
      && photon->isoChargedHadronsEA < 2.67
      && photon->isoNeutralHadronsEA < 7.23 + exp(0.0028*photon->p.Pt()+0.5408)
      && photon->isoPhotonsEA < 2.11 + 0.0014*photon->p.Pt()
    )
    eff["h_g_hOverE__"+s].FillWeighted( photon->isTrue, *w, photon->hOverE );

  }


  float ht = 0;
  for( auto& j : selJets ) {
    if( j->p.Pt() > 40 && fabs(j->p.Eta()) < 3. ) {
      ht += j->p.Pt();
    }
  }

  // trigger efficiencies
  if( *crossTriggerHtFired && ht > 700 )
    eff["eff_hlt_pt__"+s].FillWeighted( *signalTriggerFired, *w, selPhotons[0]->p.Pt() );
  if( *crossTriggerPhotonFired && selPhotons[0]->p.Pt() > 100 )
    eff["eff_hlt_ht__"+s].FillWeighted( *signalTriggerFired, *w, ht );

}

void HistogramProducer::fillBkgEst( string const& s ) {
  auto metAndL = met->p;
  for( auto l : selElectrons )
    metAndL = metAndL + l->p;
  for( auto& l : selMuons )
    metAndL += l->p;

  h["h_metAndL__"+s].Fill( metAndL.Pt(), *w );

  // profile
  if( selPhotons.size() > 0 )
    prof["profile_g_pt_met__"+s].Fill( met->p.Pt(), selPhotons[0]->p.Pt(), *w );
  prof["profile_jets_met__"+s].Fill( met->p.Pt(), selJets.size(), *w );

  // isolations
  for( auto& p : selPhotons ) {
    prof["profile_cIso_met__"+s].Fill( met->p.Pt(), p->isoChargedHadronsEA, *w );
    prof["profile_nIso_met__"+s].Fill( met->p.Pt(), p->isoNeutralHadronsEA, *w );
    prof["profile_pIso_met__"+s].Fill( met->p.Pt(), p->isoPhotonsEA, *w );
    prof["profile_sie_met__"+s].Fill( met->p.Pt(), p->sigmaIetaIeta, *w );
    prof["profile_hOverE_met__"+s].Fill( met->p.Pt(), p->hOverE, *w );
    prof["profile_mva_met__"+s].Fill( met->p.Pt(), p->mvaValue, *w );
  }
  for( auto& p : selPhotons ) {
    if( p->p.Pt() < 50 || abs(p->p.Eta()) > 1.4442 ) continue;
    if(
      p->hOverE < 0.028
      && p->sigmaIetaIeta < 0.0107
      && p->isoNeutralHadronsEA < 7.23 + exp(0.0028*p->p.Pt()+0.5408)
      && p->isoPhotonsEA < 2.11 + 0.0014*p->p.Pt()
    )
    prof["profile_cIso_met_gloose__"+s].Fill( met->p.Pt(), p->isoChargedHadronsEA, *w );
    if(
      p->hOverE < 0.028
      && p->sigmaIetaIeta < 0.0107
      && p->isoChargedHadronsEA < 2.67
      && p->isoPhotonsEA < 2.11 + 0.0014*p->p.Pt()
    )
    prof["profile_nIso_met_gloose__"+s].Fill( met->p.Pt(), p->isoNeutralHadronsEA, *w );
    if(
      p->hOverE < 0.028
      && p->sigmaIetaIeta < 0.0107
      && p->isoChargedHadronsEA < 2.67
      && p->isoNeutralHadronsEA < 7.23 + exp(0.0028*p->p.Pt()+0.5408)
    )
    prof["profile_pIso_met_gloose__"+s].Fill( met->p.Pt(), p->isoPhotonsEA, *w );
    if(
      p->hOverE < 0.028
      && p->isoChargedHadronsEA < 2.67
      && p->isoNeutralHadronsEA < 7.23 + exp(0.0028*p->p.Pt()+0.5408)
      && p->isoPhotonsEA < 2.11 + 0.0014*p->p.Pt()
    )
    prof["profile_sie_met_gloose__"+s].Fill( met->p.Pt(), p->sigmaIetaIeta, *w );
    if(
      p->sigmaIetaIeta < 0.0107
      && p->isoChargedHadronsEA < 2.67
      && p->isoNeutralHadronsEA < 7.23 + exp(0.0028*p->p.Pt()+0.5408)
      && p->isoPhotonsEA < 2.11 + 0.0014*p->p.Pt()
    )
    prof["profile_hOverE_met_gloose__"+s].Fill( met->p.Pt(), p->hOverE, *w );
    if(
      p->hOverE < 0.028
      && p->sigmaIetaIeta < 0.0107
      && p->isoChargedHadronsEA < 2.67
      && p->isoNeutralHadronsEA < 7.23 + exp(0.0028*p->p.Pt()+0.5408)
      && p->isoPhotonsEA < 2.11 + 0.0014*p->p.Pt()
    )
    prof["profile_mva_met_gloose__"+s].Fill( met->p.Pt(), p->mvaValue, *w );

  }

}

HistogramProducer::HistogramProducer():
  photons( fReader, "photons" ),
  jets( fReader, "jets" ),
  electrons( fReader, "electrons" ),
  muons( fReader, "muons" ),
  genJets( fReader, "genJets" ),
  genParticles( fReader, "genParticles" ),
  met( fReader, "met" ),
  nGoodVertices( fReader, "nGoodVertices" ),
  w( fReader, "pu_weight" ),
  signalTriggerFired( fReader, "HLT_Photon90_CaloIdL_PFHT500_v" ),
  crossTriggerPhotonFired( fReader, "HLT_Photon90_v" ),
  crossTriggerHtFired( fReader, "HLT_PFHT600_v" )
{
}

void HistogramProducer::Init(TTree *tree)
{
  fReader.SetTree(tree);
}

void HistogramProducer::SlaveBegin(TTree *tree)
{
  vector<string> strs = { "tr", "tr_met200", "tr_met200_dphi", "tr_met200_dphi_j2", "tr_met200_dphi_j2_genElectron", "tr_met200_dphi_j2_genPhoton" };
  for( auto& v : strs ) initSelection(v);
  initObjects("base");

  // after all initializations
  for( auto& it : eff ) it.second.SetUseWeightedEvents();
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
  //if( entry > 3 ) return true;
  if(!( entry%10000 )) printf( "\r%lli / %lli", entry, fReader.GetEntries(false) );
  fReader.SetEntry(entry);

  // Base selection, without cuts
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
  fillObjects("base");

  resetSelection();
  // New selection
  float ht = 0;
  for( auto& jet : jets ){
    if( jet.p.Pt() > 40 && abs(jet.p.Eta()) < 3 ) {
      ht += jet.p.Pt();
    }
  }

  for( auto& photon : photons ) {
    if( !photon.isLoose || photon.p.Pt() < 100 || abs( photon.p.Eta() ) > 1.4442 || photon.hasPixelSeed ) continue;
    selPhotons.push_back( &photon );
  }
  defaultSelection();

  if( selPhotons.size() && ht > 600 ) {
    // Calculate razor variables
    if( true ) { // attention: computing intensive
      vector<TVector3> js;
      for( auto& p : selJets ) js.push_back( p->p );
      for( auto& p : selPhotons ) js.push_back( p->p );
      mrr2 = razorVariables( megajets( js ), met->p );
    }

    fillSelection("tr");
    if( met->p.Pt() > 200 ) {
      fillSelection("tr_met200");
      // angle between met,j1 and met,j2 must be larger than 0.4
      if( !selJets.size() || ( fabs(met->p.DeltaPhi( selJets[0]->p )) > 0.4 && ( selJets.size() == 1 || fabs(met->p.DeltaPhi( selJets[1]->p )) > 0.4 ) ) ) {
        fillSelection("tr_met200_dphi");
        if( selJets.size() > 1 ) {
          fillSelection("tr_met200_dphi_j2");
          if( selPhotons[0]->isTrueAlternative ) {
            fillSelection("tr_met200_dphi_j2_genPhoton");
          }
          // find generated electron nearyby
          for( auto& p : genParticles ) {
            if( abs(p.pdgId) == 11 && p.p.DeltaR( selPhotons[0]->p ) < 0.1 ) {
              fillSelection("tr_met200_dphi_j2_genElectron");
              break;
            }
          }
        }
      }
    }
  }


  resetSelection();
  /*
  // electron sample

  for( auto& photon : photons ) {
    if( !photon.isLoose || photon.p.Pt() < 100 || abs( photon.p.Eta() ) > 1.4442 || !photon.hasPixelSeed ) continue;
    selPhotons.push_back( &photon );
  }
  defaultSelection();

  if( selPhotons.size() && ht > 600 ) {
     fillSelection("looseElectron");
  }

  // jet sample
  resetSelection();
  for( auto& photon : photons ) {
    if( photon.p.Pt() < 100 || abs( photon.p.Eta() ) > 1.4442 || photon.hasPixelSeed ) continue;
    if( photon.hOverE < 0.028
        && photon.sigmaIetaIeta < 0.0107
    //    && photon.isoChargedHadronsEA < 2.67
        && photon.isoChargedHadronsEA > 2.67 // more chargediso
        && photon.isoNeutralHadronsEA < 7.23 + exp(0.0028*photon.p.Pt()+0.5408)
        && photon.isoPhotonsEA < 2.11 + 0.0014*photon.p.Pt()
      )
      selPhotons.push_back( &photon );
  }
  defaultSelection();

  if( selPhotons.size() && ht > 600 ) {
     fillSelection("looseJet");
  }
*/
  return kTRUE;
}

void HistogramProducer::Terminate()
{
  cout << endl;

  auto outputName = getOutputFilename( fReader.GetTree()->GetCurrentFile()->GetName() );
  cout << "Writing to output file " << outputName << endl;

  // save all defined histograms to file
  TFile file( outputName.c_str(), "RECREATE");

  for( auto& mapIt : h )
    mapIt.second.Write( mapIt.first.c_str(), TObject::kWriteDelete );
  for( auto& mapIt : h2 )
    mapIt.second.Write( mapIt.first.c_str(), TObject::kWriteDelete );
  for( auto& mapIt : prof )
    mapIt.second.Write( mapIt.first.c_str(), TObject::kWriteDelete );
  for( auto& mapIt : eff )
    mapIt.second.Write( mapIt.first.c_str(), TObject::kWriteDelete );


  fReader.GetTree()->GetCurrentFile()->Get("TreeWriter/hCutFlow")->Write();

}

void HistogramProducer::resetSelection() {
  mrr2.first = 0;
  mrr2.second = 0;
  selPhotons.clear();
  selJets.clear();
  selBJets.clear();
  selElectrons.clear();
  selMuons.clear();
}
