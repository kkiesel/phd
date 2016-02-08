#include<regex>

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

#include "TreeParticles.hpp"
#include "UserFunctions.h"
#include "Weighter.h"
#include "CutFlow.h"


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
  const TVector3* matchedJet( const tree::Particle& p );

  void initSelection( string const& s );
  void fillSelection( string const& s );

  void initObjects( string const& s );
  void fillObjects( string const& s );

  TTreeReader fReader;
  TTreeReaderArray<tree::Photon> photons;
  TTreeReaderArray<tree::Jet> jets;
  TTreeReaderArray<tree::Electron> electrons;
  TTreeReaderArray<tree::Muon> muons;
  TTreeReaderArray<tree::Particle> genJets;
  TTreeReaderArray<tree::GenParticle> genParticles;

  TTreeReaderValue<tree::MET> met;
  TTreeReaderValue<Float_t> pu_weight;
  TTreeReaderValue<Char_t> mc_weight;

  TTreeReaderValue<Int_t> nGoodVertices;
  TTreeReaderValue<Int_t> genLeptonsFromW;
  TTreeReaderValue<Float_t> genHt;
  TTreeReaderValue<Float_t> dR_recoGenJet;

  TTreeReaderValue<ULong64_t> eventNo;
  TTreeReaderValue<UInt_t> runNo;
  TTreeReaderValue<UInt_t> lumNo;

  TTreeReaderValue<Bool_t> hlt_photon90_ht500;
  TTreeReaderValue<Bool_t> hlt_photon90;
  TTreeReaderValue<Bool_t> hlt_photon175;
  TTreeReaderValue<Bool_t> hlt_ht600;

  vector<tree::Photon*> selPhotons;
  vector<tree::Jet*> selJets;
  vector<tree::Jet*> selBJets;
  vector<tree::Electron*> selElectrons;
  vector<tree::Muon*> selMuons;

  vector<tree::Photon> fakeJets;

  float selW=1.; // weight
  float selHt=0;

  pair<float,float> mrr2;

  map<string,TH1F> h;
  map<string,TH2F> h2;
  map<string,TProfile> prof;
  map<string,TEfficiency> eff;

  map<int,pair<int,int>> rawEff_vs_run;

  bool isData;

  //Weighter qcdWeighter;
  Weighter htWeighter;
  Weighter qcdNFakeJetWeighter;

  CutFlowPhoton looseCutFlowPhoton;

};

void HistogramProducer::initSelection( string const& s ) {

  // event
  h["h_met__"+s] = TH1F( "", ";E^{miss}_{T} (GeV)", 100, 0, 1000 );
  h["h_metStar__"+s] = TH1F( "", ";E^{miss}_{T} (GeV)", 100, 0, 1000 );
  h["h_met_dn__"+s] = TH1F( "", ";E^{miss}_{T} (GeV)", 100, 0, 1000 );
  h["h_met_dn10__"+s] = TH1F( "", ";E^{miss}_{T} (GeV)", 100, 0, 1000 );
  h["h_met_dn09__"+s] = TH1F( "", ";E^{miss}_{T} (GeV)", 100, 0, 1000 );
  h["h_met_dn08__"+s] = TH1F( "", ";E^{miss}_{T} (GeV)", 100, 0, 1000 );
  h["h_met_dn07__"+s] = TH1F( "", ";E^{miss}_{T} (GeV)", 100, 0, 1000 );
  h["h_mt_g_met__"+s] = TH1F( "", ";m_{T}(p_{T},E^{miss}_{T}) (GeV)", 100, 0, 1000 );

  h["h_tremht__"+s] = TH1F( "", ";EMH_{T}^{trigger-like} (GeV)", 250, 0, 2500 );
  h["h_emht__"+s] = TH1F( "", ";EMH_{T} (GeV)", 250, 0, 2500 );
  h["h_emhtStar__"+s] = TH1F( "", ";EMH_{T}* (GeV)", 250, 0, 2500 );
  h["h_ht__"+s] = TH1F( "", ";H_{T} (GeV)", 250, 0, 2500 );
  h["h_st__"+s] = TH1F( "", ";S_{T} (GeV)", 250, 0, 2500 );
  h["h_recoilt__"+s] = TH1F( "", ";#vec{H}_{T} (GeV)", 150, 0, 1500 );
  h["h_meg__"+s] = TH1F( "", ";m_{ee} (GeV)", 600, 0, 600 );

  // photon
  h["h_g_pt__"+s] = TH1F( "", ";p_{T} (GeV)", 150, 0, 1500 );
  h["h_g_ptStar__"+s] = TH1F( "", ";p_{T}* (GeV)", 150, 0, 1500 );
  h["h_g_eta__"+s] = TH1F( "", ";|#eta|", 1500, 0, 1.5 );
  h["h_g_phi__"+s] = TH1F( "", ";#phi", 200, -3.2, 3.2 );
  h["h_g_sie__"+s] = TH1F( "", ";#sigma_{i#etai#eta}", 400, 0, 0.02 );
  h["h_g_hoe__"+s] = TH1F( "", ";H/E", 150, 0, 0.15 );
  h["h_g_cIso__"+s] = TH1F( "", ";I_{#pi} (GeV)", 100, 0, 10 );
  h["h_g_nIso__"+s] = TH1F( "", ";I_{n} (GeV)", 100, 0, 20 );
  h["h_g_pIso__"+s] = TH1F( "", ";I_{#gamma} (GeV)", 100, 0, 20 );

  // jet
  h["h_j1_pt__"+s] = TH1F( "", ";p_{T}^{1.jet} (GeV)", 150, 0, 1500 );
  h["h_j1_eta__"+s] = TH1F( "", ";|#eta^{1.jet}|", 150, 0, 3 );
  h["h_j2_pt__"+s] = TH1F( "", ";p_{T}^{2.jet} (GeV)", 150, 0, 1500 );
  h["h_j2_eta__"+s] = TH1F( "", ";|#eta^{2.jet}|", 150, 0, 3 );
  h["h_j3_pt__"+s] = TH1F( "", ";p_{T}^{3.jet} (GeV)", 150, 0, 1500 );
  h["h_j3_eta__"+s] = TH1F( "", ";|#eta^{3.jet}|", 150, 0, 3 );
  h["h_j1_eb_pt__"+s] = TH1F( "", ";p_{T}^{1.jet} (GeV)", 150, 0, 1500 );
  h["h_j1_ee_pt__"+s] = TH1F( "", ";p_{T}^{1.jet} (GeV)", 150, 0, 1500 );

  h["h_j1_chf__"+s] = TH1F( "", ";1.jet charged hadron fraction", 100, 0, 1 );
  h["h_j1_nhf__"+s] = TH1F( "", ";1.jet neutral hadron fraction", 100, 0, 1 );
  h["h_j1_cef__"+s] = TH1F( "", ";1.jet charged electomagnetic fraction", 100, 0, 1 );
  h["h_j1_nef__"+s] = TH1F( "", ";1.jet neutral electomagnetic fraction", 100, 0, 1 );
  h["h_j1_nch__"+s] = TH1F( "", ";1.jet charged multiplicity", 41, -0.5, 40.5 );
  h["h_j1_nconstituents__"+s] = TH1F( "", ";1.jet constituent multiplicity", 71, -0.5, 70.5 );

  // b-jets
  h["h_bj1_pt__"+s] = TH1F( "", ";p_{T}^{1.b-jet} (GeV)", 150, 0, 1500 );
  h["h_bj1_eta__"+s] = TH1F( "", ";|#eta^{1.b-jet}|", 150, 0, 3 );
  h["h_bj2_pt__"+s] = TH1F( "", ";p_{T}^{2.b-jet} (GeV)", 150, 0, 1500 );
  h["h_bj2_eta__"+s] = TH1F( "", ";|#eta^{2.b-jet}|", 150, 0, 3 );
  h["h_bj3_pt__"+s] = TH1F( "", ";p_{T}^{3.b-jet} (GeV)", 150, 0, 1500 );
  h["h_bj3_eta__"+s] = TH1F( "", ";|#eta^{3.b-jet}|", 150, 0, 3 );

  // angles
  h["h_dphi_met_g__"+s] = TH1F( "", ";|#Delta#phi(E_{T}^{miss},#gamma)|", 70, 0, 3.5 );
  h["h_dphi_met_j1__"+s] = TH1F( "", ";|#Delta#phi(E_{T}^{miss},1.jet)|", 70, 0, 3.5 );
  h["h_dphi_met_j2__"+s] = TH1F( "", ";|#Delta#phi(E_{T}^{miss},2.jet)|", 70, 0, 3.5 );
  h["h_dphi_met_j3__"+s] = TH1F( "", ";|#Delta#phi(E_{T}^{miss},3.jet)|", 70, 0, 3.5 );
  h["h_dphi_met_j13_min__"+s] = TH1F( "", ";min_{i}|#Delta#phi(E_{T}^{miss},i.jet)|", 70, 0, 3.5 );
  h["h_dphi_met_recoil__"+s] = TH1F( "", ";|#Delta#phi(E_{T}^{miss},#Sigma jet)|", 70, 0, 3.5 );
  h["h_dphi_g_j1__"+s] = TH1F( "", ";|#Delta#phi(#gamma,1.jet)|", 70, 0, 3.5 );
  h["h_dphi_g_j2__"+s] = TH1F( "", ";|#Delta#phi(#gamma,2.jet)|", 70, 0, 3.5 );

  // multiplicities
  h["h_n_vertex__"+s] = TH1F( "", ";vertex multiplicity", 61, -0.5, 60.5 );
  h["h_n_photon__"+s] = TH1F( "", ";photon multiplicity", 6, -0.5, 5.5 );
  h["h_n_jet__"+s] = TH1F( "", ";jet multiplicity", 21, -0.5, 20.5 );
  h["h_n_bjet__"+s] = TH1F( "", ";b-jet multiplicity", 21, -0.5, 20.5 );
  h["h_n_electron__"+s] = TH1F( "", ";electron multiplicity", 4, -0.5, 3.5 );
  h["h_n_muon__"+s] = TH1F( "", ";muon multiplicity", 4, -0.5, 3.5 );

  h["h_n_fakeJet__"+s] = TH1F( "", ";photon-like jet multiplicity", 11, -0.5, 10.5 );
  h["h_n_photonAndfakeJet__"+s] = TH1F( "", ";photon and photon-like jet multiplicity", 11, -0.5, 10.5 );

  h2["h2_razorPlane__"+s] = TH2F( "", ";M_{R} (GeV); R^{2}", 100, 0, 2000, 100, 0, .5 );
  h2["h2_g_pt_emht__"+s] = TH2F( "", ";p_{T} (GeV); EMH_{T} (GeV)", 100, 0, 1000, 150, 500, 2000 );
  h2["h2_n_fakeJets_photonPosition__"+s] = TH2F( "", ";jet fake multiplicity;photon position", 9, -0.5, 8.5, 9, -0.5, 8.5 );

}

void HistogramProducer::initObjects( string const& s ) {

  eff["eff_gCol_sigmaIetaIeta__"+s] = TEfficiency( "", ";#sigma_{i#etai#eta}", 100, 0, 0.02 );
  eff["eff_gCol_hOverE__"+s] = TEfficiency( "", ";H/E", 100, 0, 0.15 );
  eff["eff_gCol_r9__"+s] = TEfficiency( "", ";r9", 110, 0, 1.1 );
  eff["eff_gCol_cIso__"+s] = TEfficiency( "", ";I_{#pi}", 200, 0, 150 );
  eff["eff_gCol_mva__"+s] = TEfficiency( "", ";y_{mva}", 120, -1.1, 1.1 );

  eff["eff_g_sigmaIetaIeta__"+s] = TEfficiency( "", ";#sigma_{i#etai#eta}", 100, 0, 0.02 );
  eff["eff_g_hOverE__"+s] = TEfficiency( "", ";H/E", 100, 0, 0.15 );
  eff["eff_g_r9__"+s] = TEfficiency( "", ";r9", 110, 0, 1.1 );
  eff["eff_g_cIso__"+s] = TEfficiency( "", ";I_{#pi}", 100, 0, 10 );
  eff["eff_g_nIso__"+s] = TEfficiency( "", ";I_{n}", 200, 0, 20 );
  eff["eff_g_pIso__"+s] = TEfficiency( "", ";I_{#gamma}", 200, 0, 20 );


  // matching
  h2["h2_match_jet_photon__"+s] = TH2F( "", ";#Delta R;p_{T}^{jet}/p_{T}^{#gamma}", 100, 0, 0.5, 100, 0, 4 );
  h2["h2_match_jet_electron__"+s] = TH2F( "", ";#Delta R;p_{T}^{jet}/p_{T}^{e}", 100, 0, 0.5, 100, 0, 4 );
  h2["h2_match_jet_muon__"+s] = TH2F( "", ";#Delta R;p_{T}^{jet}/p_{T}^{#mu}", 100, 0, 0.5, 100, 0, 4 );
  h2["h2_match_photon_pt_jet_pt_"+s] = TH2F( "", ";p_{T}^{photon} (GeV);p_{T}^{jet} (GeV)", 100, 0, 1000, 100, 0, 1000 );

  h2["h2_match_photon_genElectron__"+s] = TH2F( "", ";#Delta R;p_{T}/p_{T}^{gen e}", 100, 0, 0.5, 100, 0, 4 );
  h2["h2_match_photon_genPhoton__"+s] = TH2F( "", ";#Delta R;p_{T}/p_{T}^{gen #gamma}", 100, 0, 0.5, 100, 0, 4 );

  eff["eff_hlt_pt__offlineHT650__"+s] = TEfficiency( "", ";p_{T} (GeV);#varepsilon_{EMH_{T}^{trigger-like}>650GeV}", 250, 0, 1000 );
  eff["eff_hlt_pt__"+s] = TEfficiency( "", ";p_{T} (GeV);#varepsilon", 250, 0, 1000 );
  eff["eff_hlt_ht__offlinePT100__"+s] = TEfficiency( "", ";EMH_{T}^{trigger-like} (GeV);#varepsilon_{p_{T}>100GeV}", 200, 0, 2000 );
  eff["eff_hlt_ht__"+s] = TEfficiency( "", ";EMH_{T}^{trigger-like} (GeV);#varepsilon", 200, 0, 2000 );
  eff["eff_hlt_nVertex__"+s] = TEfficiency( "", ";Vertex multiplicity", 41, -0.5, 40.5 );
  eff["eff_hlt_sie__"+s] = TEfficiency( "", ";#sigma_{i#etai#eta}", 400, 0, 0.02 );
  eff["eff_hlt_hoe__"+s] = TEfficiency( "", ";H/E", 100, 0, 0.15 );
  eff["eff_hlt_cIso__"+s] = TEfficiency( "", ";I_{#pi} (GeV)", 100, 0, 10 );
  eff["eff_hlt_nIso__"+s] = TEfficiency( "", ";I_{n} (GeV)", 100, 0, 20 );
  eff["eff_hlt_pIso__"+s] = TEfficiency( "", ";I_{#gamma} (GeV)", 100, 0, 20 );
  eff["eff_hlt_n_jet_uncleaned__"+s] = TEfficiency( "", ";uncleaned jet multiplicity", 15, -0.5, 14.5 );

  h["h_genHt"] = TH1F( "", ";H_{T}^{gen} (GeV)", 6000, 0, 3000 );

}

void HistogramProducer::fillSelection( string const& s ) {

  float ht = 0;
  TVector3 recoil(0,0,0);
  for( auto& j: selJets ){
    ht += j->p.Pt();
    recoil += j->p;
  }
  float ht_g = ht;
  for( auto& p : selPhotons ) ht_g += p->p.Pt();
  float st = ht_g + met->p.Pt();
  float ht_gStar = ht;
  for( auto& p : selPhotons ) {
    ht_gStar += matchedJet(*p)->Pt();
  }

  h.at("h_tremht__"+s).Fill( selHt, selW );
  h.at("h_emht__"+s).Fill( ht_g, selW );
  h.at("h_emhtStar__"+s).Fill( ht_gStar, selW );
  h.at("h_ht__"+s).Fill( ht, selW );
  h.at("h_st__"+s).Fill( st, selW );
  h.at("h_recoilt__"+s).Fill( recoil.Pt(), selW );

  h.at("h_met__"+s).Fill( met->p.Pt(), selW );
  h.at("h_metStar__"+s).Fill( (met->p-(*matchedJet(*selPhotons.at(0))-selPhotons.at(0)->p)).Pt(), selW );
  h.at("h_met_dn__"+s).Fill( (met->p-0.1*selPhotons.at(0)->p).Pt(), selW );
  h.at("h_met_dn10__"+s).Fill( (met->p-0.1*selPhotons.at(0)->p).Pt(), selW );
  h.at("h_met_dn09__"+s).Fill( (met->p-0.09*selPhotons.at(0)->p).Pt(), selW );
  h.at("h_met_dn08__"+s).Fill( (met->p-0.08*selPhotons.at(0)->p).Pt(), selW );
  h.at("h_met_dn07__"+s).Fill( (met->p-0.07*selPhotons.at(0)->p).Pt(), selW );

  if( selElectrons.size() && selPhotons.size() ){
    h.at("h_meg__"+s).Fill( (selPhotons.at(0)->p+selElectrons.at(0)->p).Mag() );
  }

  if( selPhotons.size() > 0 ) {
    h.at("h_mt_g_met__"+s).Fill( (selPhotons.at(0)->p + met->p).Pt(), selW );
    h.at("h_g_pt__"+s).Fill( selPhotons.at(0)->p.Pt(), selW );
    h.at("h_g_ptStar__"+s).Fill( matchedJet(*selPhotons.at(0))->Pt(), selW );
    h.at("h_g_eta__"+s).Fill( fabs(selPhotons.at(0)->p.Eta()), selW );
    h.at("h_g_phi__"+s).Fill( selPhotons.at(0)->p.Phi(), selW );
    h.at("h_g_sie__"+s).Fill( selPhotons.at(0)->sigmaIetaIeta );
    h.at("h_g_hoe__"+s).Fill( selPhotons.at(0)->hOverE );
    h.at("h_g_cIso__"+s).Fill( selPhotons.at(0)->isoChargedHadronsEA );
    h.at("h_g_nIso__"+s).Fill( selPhotons.at(0)->isoNeutralHadronsEA );
    h.at("h_g_pIso__"+s).Fill( selPhotons.at(0)->isoPhotonsEA );
    h2.at("h2_g_pt_emht__"+s).Fill( selPhotons.at(0)->p.Pt(), selHt, selW );
    unsigned photonPosition=0;
    for(;photonPosition<fakeJets.size() && fakeJets.at(photonPosition).p.Pt() > selPhotons.at(0)->p.Pt();photonPosition++);
    h2.at("h2_n_fakeJets_photonPosition__"+s).Fill( fakeJets.size(), photonPosition, selW );
  }

  if( selJets.size() > 2 ) {
    h.at("h_j3_pt__"+s).Fill( selJets.at(2)->p.Pt(), selW );
    h.at("h_j3_eta__"+s).Fill( fabs(selJets.at(2)->p.Eta()), selW );
  }
  if( selJets.size() > 1 ) {
    h.at("h_j2_pt__"+s).Fill( selJets.at(1)->p.Pt(), selW );
    h.at("h_j2_eta__"+s).Fill( fabs(selJets.at(1)->p.Eta()), selW );
  }
  if( selJets.size() > 0 ) {
    h.at("h_j1_pt__"+s).Fill( selJets.at(0)->p.Pt(), selW );
    h.at("h_j1_eta__"+s).Fill( fabs(selJets.at(0)->p.Eta()), selW );
    if( fabs(selJets.at(0)->p.Eta())<1.45 ) h.at("h_j1_eb_pt__"+s).Fill( selJets.at(0)->p.Pt(), selW );
    else                                 h.at("h_j1_ee_pt__"+s).Fill( selJets.at(0)->p.Pt(), selW );
    h.at("h_j1_chf__"+s).Fill( selJets.at(0)->chf, selW );
    h.at("h_j1_nhf__"+s).Fill( selJets.at(0)->nhf, selW );
    h.at("h_j1_cef__"+s).Fill( selJets.at(0)->cef, selW );
    h.at("h_j1_nef__"+s).Fill( selJets.at(0)->nef, selW );
    h.at("h_j1_nch__"+s).Fill( selJets.at(0)->nch, selW );
    h.at("h_j1_nconstituents__"+s).Fill( selJets.at(0)->nconstituents, selW );
  }

  if( selBJets.size() > 2 ) {
    h.at("h_bj3_pt__"+s).Fill( selBJets.at(2)->p.Pt(), selW );
    h.at("h_bj3_eta__"+s).Fill( fabs(selBJets.at(2)->p.Eta()), selW );
  }
  if( selBJets.size() > 1 ) {
    h.at("h_bj2_pt__"+s).Fill( selBJets.at(1)->p.Pt(), selW );
    h.at("h_bj2_eta__"+s).Fill( fabs(selBJets.at(1)->p.Eta()), selW );
  }
  if( selBJets.size() > 0 ) {
    h.at("h_bj1_pt__"+s).Fill( selBJets.at(0)->p.Pt(), selW );
    h.at("h_bj1_eta__"+s).Fill( fabs(selBJets.at(0)->p.Eta()), selW );
  }

  if( selPhotons.size() > 0 )
    h.at("h_dphi_met_g__"+s).Fill( fabs(met->p.DeltaPhi( selPhotons.at(0)->p )), selW );

  vector<float> minDeltaPhiMetJet;
  if( selJets.size() > 0 ) {
    float dphi = fabs(met->p.DeltaPhi( selJets.at(0)->p ));
    h.at("h_dphi_met_j1__"+s).Fill( dphi, selW );
    if( selPhotons.size() > 0 ) h.at("h_dphi_g_j1__"+s).Fill( fabs(selPhotons.at(0)->p.DeltaPhi(selJets.at(0)->p)), selW );
    minDeltaPhiMetJet.push_back( dphi );
  }
  if( selJets.size() > 1 ) {
    float dphi = fabs(met->p.DeltaPhi( selJets.at(1)->p ));
    h.at("h_dphi_met_j2__"+s).Fill( dphi, selW );
    if( selPhotons.size() > 0 ) h.at("h_dphi_g_j2__"+s).Fill( fabs(selPhotons.at(0)->p.DeltaPhi(selJets.at(1)->p)), selW );
    minDeltaPhiMetJet.push_back( dphi );
  }
  if( selJets.size() > 2 ) {
    float dphi = fabs(met->p.DeltaPhi( selJets.at(2)->p ));
    h.at("h_dphi_met_j3__"+s).Fill( dphi, selW );
    minDeltaPhiMetJet.push_back( dphi );
  }
  if( minDeltaPhiMetJet.size() ) {
    h.at("h_dphi_met_j13_min__"+s).Fill( *min_element( minDeltaPhiMetJet.begin(), minDeltaPhiMetJet.end() ), selW);
  }

  h.at("h_dphi_met_recoil__"+s).Fill( fabs(met->p.DeltaPhi( recoil )), selW );

  h.at("h_n_vertex__"+s).Fill( *nGoodVertices, selW );
  h.at("h_n_photon__"+s).Fill( selPhotons.size(), selW );
  h.at("h_n_jet__"+s).Fill( selJets.size(), selW );
  h.at("h_n_bjet__"+s).Fill( selBJets.size(), selW );
  h.at("h_n_electron__"+s).Fill( selElectrons.size(), selW );
  h.at("h_n_muon__"+s).Fill( selMuons.size(), selW );

  h.at("h_n_fakeJet__"+s).Fill( fakeJets.size(), selW );
  unsigned nFakeJets = fakeJets.size();
  if( selPhotons.size() && selPhotons.at(0)->r9 > 0) nFakeJets++; // only for real photons, else it is contained in fakeJets
  h.at("h_n_photonAndfakeJet__"+s).Fill( nFakeJets, selW );

  if( mrr2.first > 1e-7 ) {
    h2.at("h2_razorPlane__"+s).Fill( mrr2.first, mrr2.second, selW );
  }

} // end fill

void HistogramProducer::fillObjects( string const& s ) {
  h.at("h_genHt").Fill( *genHt, selW );

  // matching
  for( auto& jet : selJets ) {
    for( auto& p : selPhotons ) {
      h2.at("h2_match_jet_photon__"+s).Fill( jet->p.DeltaR( p->p ), jet->p.Pt()/p->p.Pt(), selW );
      if( p->isLoose && p->p.Pt()>100 && fabs(p->p.Eta())<photonsEtaMaxBarrel && jet->p.DeltaR(p->p)<.3 ) h2.at("h2_match_photon_pt_jet_pt_"+s).Fill( p->p.Pt(), jet->p.Pt(), selW );
    }
    for( auto& p : selElectrons )
      h2.at("h2_match_jet_electron__"+s).Fill( jet->p.DeltaR( p->p ), jet->p.Pt()/p->p.Pt(), selW );
    for( auto& p : selMuons )
      h2.at("h2_match_jet_muon__"+s).Fill( jet->p.DeltaR( p->p ), jet->p.Pt()/p->p.Pt(), selW );

  }

  for( auto& photon : selPhotons ) {
    for( auto& p : genParticles ) {
      if( fabs(p.pdgId) == 22 )
        h2.at("h2_match_photon_genPhoton__"+s).Fill( photon->p.DeltaR( p.p ), photon->p.Pt()/p.p.Pt(), selW );
      if( fabs(p.pdgId) == 11 )
        h2.at("h2_match_photon_genElectron__"+s).Fill( photon->p.DeltaR( p.p ), photon->p.Pt()/p.p.Pt(), selW );
    }

    eff.at("eff_gCol_sigmaIetaIeta__"+s).FillWeighted( photon->isTrue, selW, photon->sigmaIetaIeta );
    eff.at("eff_gCol_hOverE__"+s).FillWeighted( photon->isTrue, selW, photon->hOverE );
    eff.at("eff_gCol_r9__"+s).FillWeighted( photon->isTrue, selW, photon->r9 );
    eff.at("eff_gCol_cIso__"+s).FillWeighted( photon->isTrue, selW, photon->isoChargedHadronsEA );
    eff.at("eff_gCol_mva__"+s).FillWeighted( photon->isTrue, selW, photon->mvaValue );

    if( photon->p.Pt() > 100 && fabs(photon->p.Eta()) < photonsEtaMaxBarrel ) {
      if(
        looseCutFlowPhoton.check( *photon )
        ) eff.at("eff_g_r9__"+s).FillWeighted( photon->isTrue, selW, photon->r9 );
      if(
        looseCutFlowPhoton.passHoe()
        && looseCutFlowPhoton.passSie()
        && looseCutFlowPhoton.passNIso()
        && looseCutFlowPhoton.passPIso()
        ) eff.at("eff_g_cIso__"+s).FillWeighted( photon->isTrue, selW, photon->isoChargedHadronsEA );
      if(
        looseCutFlowPhoton.passHoe()
        && looseCutFlowPhoton.passSie()
        && looseCutFlowPhoton.passCIso()
        && looseCutFlowPhoton.passPIso()
        ) eff.at("eff_g_nIso__"+s).FillWeighted( photon->isTrue, selW, photon->isoNeutralHadronsEA );
      if(
        looseCutFlowPhoton.passHoe()
        && looseCutFlowPhoton.passSie()
        && looseCutFlowPhoton.passCIso()
        && looseCutFlowPhoton.passNIso()
        ) eff.at("eff_g_pIso__"+s).FillWeighted( photon->isTrue, selW, photon->isoPhotonsEA );
      if(
        looseCutFlowPhoton.passHoe()
        && looseCutFlowPhoton.passIso()
      ) eff.at("eff_g_sigmaIetaIeta__"+s).FillWeighted( photon->isTrue, selW, photon->sigmaIetaIeta );
      if(
        looseCutFlowPhoton.passSie()
        && looseCutFlowPhoton.passIso()
      ) eff.at("eff_g_hOverE__"+s).FillWeighted( photon->isTrue, selW, photon->hOverE );
    }
  }

  // trigger efficiencies
  tree::Photon* thisPhoton=0;
  for( auto& photon : selPhotons ) {
    if( photon->p.Pt() > 15 && fabs(photon->p.Eta()) < photonsEtaMaxBarrel && photon->isLoose && !photon->hasPixelSeed ) {
      if( !thisPhoton || thisPhoton->p.Pt() < photon->p.Pt() ) {
        thisPhoton = photon;
      }
    }
  }
  if( thisPhoton ) {

    float ht = 0;
    for( auto& jet : selJets ) {
      if( jet->p.Pt() > 40 && fabs(jet->p.Eta()) < 3 ) {
        ht += jet->p.Pt();
      }
    }

    if( *hlt_ht600 ) {
      eff.at("eff_hlt_pt__"+s).Fill( *hlt_photon90_ht500, thisPhoton->p.Pt() );
      if( ht > 650 ) eff.at("eff_hlt_pt__offlineHT650__"+s).Fill( *hlt_photon90_ht500, thisPhoton->p.Pt() );
    }
    if( *hlt_photon90 ) {
      eff.at("eff_hlt_ht__"+s).Fill( *hlt_photon90_ht500, ht );
      if( thisPhoton->p.Pt() > 100 ) eff.at("eff_hlt_ht__offlinePT100__"+s).Fill( *hlt_photon90_ht500, ht );
    }
    if( thisPhoton->p.Pt()>100 && *hlt_ht600 ) {
      eff.at("eff_hlt_nVertex__"+s).Fill( *hlt_photon90_ht500, *nGoodVertices );
      eff.at("eff_hlt_hoe__"+s).Fill( *hlt_photon90_ht500, thisPhoton->hOverE );
      eff.at("eff_hlt_sie__"+s).Fill( *hlt_photon90_ht500, thisPhoton->sigmaIetaIeta );
      eff.at("eff_hlt_cIso__"+s).Fill( *hlt_photon90_ht500, thisPhoton->isoChargedHadronsEA );
      eff.at("eff_hlt_nIso__"+s).Fill( *hlt_photon90_ht500, thisPhoton->isoNeutralHadronsEA );
      eff.at("eff_hlt_pIso__"+s).Fill( *hlt_photon90_ht500, thisPhoton->isoPhotonsEA );
      eff.at("eff_hlt_n_jet_uncleaned__"+s).Fill( *hlt_photon90_ht500, selJets.size() );
      if( !rawEff_vs_run.count( *runNo ) ) rawEff_vs_run[*runNo] = make_pair(0,0);
      if( *hlt_photon90_ht500 ) rawEff_vs_run.at(*runNo).first += 1;
      rawEff_vs_run.at(*runNo).second += 1;
    }

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
  pu_weight( fReader, "pu_weight" ),
  mc_weight( fReader, "mc_weight" ),
  genLeptonsFromW( fReader, "genLeptonsFromW" ),
  genHt( fReader, "genHt" ),
  dR_recoGenJet( fReader, "dR_recoGenJet" ),
  runNo( fReader, "runNo" ),
  hlt_photon90_ht500( fReader, "HLT_Photon90_CaloIdL_PFHT500_v" ),
  hlt_photon90( fReader, "HLT_Photon90_v" ),
  hlt_photon175( fReader, "HLT_Photon175_v" ),
  hlt_ht600( fReader, "HLT_PFHT600_v" ),
  //qcdWeighter( "../plotter/weights_unweighted.root", "weight__data_g_pt__tr_jControl"),
  htWeighter( "../plotter/weights_unweighted.root", "weight_mc_data_ht__tr"),
  qcdNFakeJetWeighter( "../plotter/toreador_weight_nJet.root", "fakeJetsWeight"),
  looseCutFlowPhoton( {{"sigmaIetaIeta_eb",0.0102}, {"cIso_eb",3.32}, {"nIso1_eb",1.92}, {"nIso2_eb",0.014}, {"nIso3_eb",0.000019}, {"pIso1_eb",0.81}, {"pIso2_eb",0.0053},
    {"sigmaIetaIeta_ee",0.0274}, {"cIso_ee",1.97}, {"nIso1_ee",11.86}, {"nIso2_ee",0.0139}, {"nIso3_ee",0.000025}, {"pIso1_ee",0.83}, {"pIso2_ee",0.0034} })
{
}

void HistogramProducer::Init(TTree *tree)
{
  fReader.SetTree(tree);
  string inputName = fReader.GetTree()->GetCurrentFile()->GetName();
  isData = inputName.find("SinglePhoton_Run") != string::npos
    || inputName.find("JetHT_Run") != string::npos
    || inputName.find("MET_Run") != string::npos;
}

void HistogramProducer::SlaveBegin(TTree *tree)
{
  initObjects("base");

  vector<string> strs = {
    "tr",
    "tr_genElectron",
    "tr_met100_genElectron",
    "tr_eControl",
    "tr_eControl_genElectron",
    "tr_eControl_met100_genElectron",
    "tr_met100_eControl",
    "tr_jControlPhotonSB",
    "tr_jControlPhotonAll",
    "tr_jControlLeadingJet",
    "tr_jControlTrailingJet",
    "tr_jControlRandomJet",
    "tr_jControlTrailingJet_nJetWeighted",
    "tr_htWeighted",
  };
  for( auto& v : strs ) initSelection(v);

  // after all initializations
  for( auto& it : eff ) it.second.SetUseWeightedEvents();
}

void HistogramProducer::defaultSelection()
{
  for( auto& mu : muons ) {
    if( mu.p.Pt() < 15 ) continue;
    if( indexOfMatchedParticle<tree::Photon*>( mu, selPhotons ) >= 0 ) continue;
    selMuons.push_back( &mu );
  }
  for( auto& el : electrons ) {
    if( !el.isLoose || el.p.Pt() < 15 ) continue;
    if( indexOfMatchedParticle<tree::Photon*>( el, selPhotons ) >= 0 ) continue;
    selElectrons.push_back( &el );
  }
  for( auto& jet : jets ) {
    if( !jet.isLoose || jet.p.Pt() < 40 || fabs(jet.p.Eta()) > 3 ) continue;
    if( indexOfMatchedParticle<tree::Photon*>( jet, selPhotons, .3 ) >= 0 ) continue;
    if( indexOfMatchedParticle<tree::Electron*>( jet, selElectrons, .3 ) >= 0 ) continue;
    if( indexOfMatchedParticle<tree::Muon*>( jet, selMuons, .3 ) >= 0 ) continue;
    if( jet.p.Pt() > 100 && fabs(jet.p.Eta()) < photonsEtaMaxBarrel ) {
      tree::Photon photon;
      photon.p = jet.p;
      fakeJets.push_back( photon );
    }
    // the jets are not cleared by the fakeJets yet
    selJets.push_back( &jet );
    if( jet.bDiscriminator > bTaggingWorkingPoints.at(CSVv2M) )
      selBJets.push_back( &jet );
  }

}

const TVector3* HistogramProducer::matchedJet( const tree::Particle& p ) {
  vector<tree::Jet*> allJets;
  for( auto& j : jets ) allJets.push_back(&j);
  auto i = indexOfMatchedParticle<tree::Jet*>( p, allJets, .1 );
  return i>=0 ? &(jets[i].p) : &(p.p);
}

Bool_t HistogramProducer::Process(Long64_t entry)
{
  resetSelection();
  fReader.SetEntry(entry);

  // set weight
  selW = *mc_weight * *pu_weight;
  float originalW = selW;

  // https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2552/1/1/1.html
  // The signal trigger effiency in this run is low.
  // Perhaps this has something to do with the bad beam spot in this and other runs
  if( isData && *runNo == 259637 ) return true;

  selHt = 0;
  for( auto& jet : jets ){
    if( jet.p.Pt()>40 && fabs(jet.p.Eta())<3 ) {
      selHt += jet.p.Pt();
    }
  }

  /////////////////////////////////////////////////////////////////////////////
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
  /////////////////////////////////////////////////////////////////////////////
  // electron sample

  for( auto& photon : photons ) {
    if( photon.isLoose && photon.hasPixelSeed && photon.p.Pt() > 100 && fabs(photon.p.Eta()) < photonsEtaMaxBarrel ) {
      selPhotons.push_back( &photon );
    }
  }
  defaultSelection();

  if( selPhotons.size() && selHt > 700 && (*hlt_photon90_ht500 || !isData) ) {
    fillSelection("tr_eControl");
    for( auto& p : genParticles ) {
      if( fabs(p.pdgId)==11 && selPhotons.at(0)->p.DeltaR(p.p)<0.2 ) {
        fillSelection("tr_eControl_genElectron");
        if( met->p.Pt()>100 ) fillSelection("tr_eControl_met100_genElectron");
        break; // fill only once
      }
    }
    if(met->p.Pt()>100) fillSelection("tr_met100_eControl");
  }

  resetSelection();
  /////////////////////////////////////////////////////////////////////////////
  // jet sample
  for( auto& photon : photons ) {
    if( photon.p.Pt() > 100 && fabs(photon.p.Eta()) < photonsEtaMaxBarrel && !photon.isLoose ) {
      selPhotons.push_back( &photon );
    }
  }
  defaultSelection();
  if( selPhotons.size() && selHt > 700 && (*hlt_photon90_ht500 || !isData) ) {
     fillSelection("tr_jControlPhotonAll");
  }

  resetSelection();
  /////////////////////////////////////////////////////////////////////////////
  // jet sample
  for( auto& photon : photons ) {
    looseCutFlowPhoton.check(photon);
    if( photon.p.Pt() > 100 && fabs(photon.p.Eta()) < photonsEtaMaxBarrel  && !photon.hasPixelSeed
        && !photon.isLoose
        && looseCutFlowPhoton.passHoe()
        && looseCutFlowPhoton.passNIso()
        && looseCutFlowPhoton.passPIso()
        && photon.isoChargedHadronsEA < 15
        && photon.sigmaIetaIeta < 0.014
        && (looseCutFlowPhoton.passCIso() ^ looseCutFlowPhoton.passSie() )
    ) {
      selPhotons.push_back( &photon );
    }
  }
  defaultSelection();
  if( selPhotons.size() && selHt > 700 && (*hlt_photon90_ht500 || !isData) ) {
     fillSelection("tr_jControlPhotonSB");
  }

  resetSelection();
  /////////////////////////////////////////////////////////////////////////////
  // Main (signal) selection
  for( auto& photon : photons ) {
    if( photon.isLoose && !photon.hasPixelSeed && photon.p.Pt() > 100 && fabs(photon.p.Eta()) < photonsEtaMaxBarrel ) {
       selPhotons.push_back( &photon );
    }
  }
  defaultSelection();

  // Calculate razor variables
  if( false ) { // attention: computing intensive
    vector<TVector3> js;
    for( auto& p : selJets ) js.push_back( p->p );
    for( auto& p : selPhotons ) js.push_back( p->p );
    for( auto& p : selElectrons ) js.push_back( p->p );
    for( auto& p : selMuons ) js.push_back( p->p );
    mrr2 = razorVariables( megajets( js ), met->p );
  }

  if( selPhotons.size() && selHt > 700 && (*hlt_photon90_ht500 || !isData) ) {
    fillSelection("tr");
    for( auto& p : genParticles ) {
      if( fabs(p.pdgId)==11 && selPhotons.at(0)->p.DeltaR(p.p)<0.2 ) {
        fillSelection("tr_genElectron");
        if( met->p.Pt()>100 ) fillSelection("tr_met100_genElectron");
        break; // fill only once
      }
    }

    float ht=0;
    for( auto& j: selJets ) ht += j->p.Pt();
    if(!isData) selW = originalW*htWeighter.getWeight( ht );
    fillSelection("tr_htWeighted");
    selW = originalW;
  }

  if( fakeJets.size() && !selPhotons.size() && selHt > 700 && (*hlt_ht600 || !isData) ) {

    auto saveSelJets = selJets;
    selJets.clear();

    selPhotons.push_back( &fakeJets.at(0) );
    for( auto& j : saveSelJets ) { if (j->p.DeltaR(selPhotons.at(0)->p)>0.3) selJets.push_back( j ); }
    fillSelection("tr_jControlLeadingJet");
    selPhotons.clear();
    selJets.clear();

    selPhotons.push_back( &fakeJets.at(fakeJets.size()-1) );
    for( auto& j : saveSelJets ) { if (j->p.DeltaR(selPhotons.at(0)->p)>0.3) selJets.push_back( j ); }
    fillSelection("tr_jControlTrailingJet");
    selW = originalW*qcdNFakeJetWeighter.getWeight( fakeJets.size() );
    fillSelection("tr_jControlTrailingJet_nJetWeighted");
    selW = originalW;
    selPhotons.clear();
    selJets.clear();

    selPhotons.push_back( &fakeJets.at(rand()%fakeJets.size()) );
    for( auto& j : saveSelJets ) { if (j->p.DeltaR(selPhotons.at(0)->p)>0.3) selJets.push_back( j ); }
    fillSelection("tr_jControlRandomJet");
    selPhotons.clear();
    selJets.clear();

  }
  resetSelection();

  return kTRUE;
}

void HistogramProducer::Terminate()
{

  auto outputName = getOutputFilename( fReader.GetTree()->GetCurrentFile()->GetName() );
  cout << "Created " << outputName;

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

  file.WriteObject( &rawEff_vs_run, "rawEff_vs_run" );

  fReader.GetTree()->GetCurrentFile()->Get("TreeWriter/hCutFlow")->Write();

}

void HistogramProducer::resetSelection() {
  selPhotons.clear();
  selJets.clear();
  selBJets.clear();
  selElectrons.clear();
  selMuons.clear();
  fakeJets.clear();
}
