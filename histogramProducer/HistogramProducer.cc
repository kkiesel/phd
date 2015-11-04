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

// https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation74X
map<string,float> bTaggingWorkingPoints = {
  { "CSVv2L", 0.605 },
  { "CSVv2M", 0.89 },
  { "CSVv2T", 0.97}
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
  auto endPos = inputFileName.find("_nTuple.root");
  string outputName = "out.root";
  if( endPos != string::npos ) {
    outputName = inputFileName.substr( startPos+1, endPos-startPos-1 ) + "_hists.root";
  }
  return outputName;

}

template <typename VectorClass>
int indexOfMatchedParticle( const tree::Particle& tag, const std::vector<VectorClass>& particles, float deltaR=.1, float relPt=-1 ) {
  int match=-1;
  for( int i=0; i<(int)particles.size(); ++i ) {
    if(    ( deltaR<0 || particles.at(i)->p.DeltaR( tag.p ) < deltaR )
        && ( relPt<0  || fabs(particles.at(i)->p.Pt()-tag.p.Pt())/tag.p.Pt() < relPt ) ) {
      return match;
    }
  }
  return match;
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

  TTreeReaderValue<Bool_t> signalTrigger;
  TTreeReaderValue<Bool_t> crossTriggerPhoton;
  TTreeReaderValue<Bool_t> crossTriggerHt;

  vector<tree::Photon*> selPhotons;
  vector<tree::Jet*> selJets;
  vector<tree::Jet*> selBJets;
  vector<tree::Electron*> selElectrons;
  vector<tree::Muon*> selMuons;

  float selW=1.; // weight
  float selHt=0;

  pair<float,float> mrr2;

  map<string,TH1F> h;
  map<string,TH2F> h2;
  map<string,TProfile> prof;
  map<string,TEfficiency> eff;

  bool isData;

};

void HistogramProducer::initSelection( string const& s ) {
  h["h_met__"+s] = TH1F( "", ";E^{miss}_{T} (GeV)", 100, 0, 1000 );
  h["h_mt_g_met__"+s] = TH1F( "", ";m_{T}(p_{T},E^{miss}_{T}) (GeV)", 100, 0, 1000 );

  h["h_ht__"+s] = TH1F( "", ";H_{T}", 200, 0, 2500 );
  h["h_ht_g__"+s] = TH1F( "", ";HE_{T}", 200, 0, 2500 );
  h["h_st__"+s] = TH1F( "", ";S_{T}", 200, 0, 2500 );

  h["h_g_pt__"+s] = TH1F( "", ";p_{T} (GeV)", 100, 0, 1500 );
  h["h_g_eta__"+s] = TH1F( "", ";|#eta|", 100, 0, 3 );

  // jet
  h["h_j1_pt__"+s] = TH1F( "", ";p_{T}^{1.jet} (GeV)", 100, 0, 1500 );
  h["h_j1_eta__"+s] = TH1F( "", ";|#eta^{1.jet}|", 100, 0, 3 );
  h["h_j2_pt__"+s] = TH1F( "", ";p_{T}^{2.jet} (GeV)", 100, 0, 1500 );
  h["h_j2_eta__"+s] = TH1F( "", ";|#eta^{2.jet}|", 100, 0, 3 );
  h["h_j3_pt__"+s] = TH1F( "", ";p_{T}^{3.jet} (GeV)", 100, 0, 1500 );
  h["h_j3_eta__"+s] = TH1F( "", ";|#eta^{3.jet}|", 100, 0, 3 );

  // b-jets
  h["h_bj1_pt__"+s] = TH1F( "", ";p_{T}^{1.b-jet} (GeV)", 100, 0, 1500 );
  h["h_bj1_eta__"+s] = TH1F( "", ";|#eta^{1.b-jet}|", 100, 0, 3 );
  h["h_bj2_pt__"+s] = TH1F( "", ";p_{T}^{2.b-jet} (GeV)", 100, 0, 1500 );
  h["h_bj2_eta__"+s] = TH1F( "", ";|#eta^{2.b-jet}|", 100, 0, 3 );
  h["h_bj3_pt__"+s] = TH1F( "", ";p_{T}^{3.b-jet} (GeV)", 100, 0, 1500 );
  h["h_bj3_eta__"+s] = TH1F( "", ";|#eta^{3.b-jet}|", 100, 0, 3 );

  // angles
  h["h_dphi_met_g__"+s] = TH1F( "", ";|#Delta#phi(E_{T}^{miss},#gamma)|", 100, 0, 3.5 );
  h["h_dphi_met_j1__"+s] = TH1F( "", ";|#Delta#phi(E_{T}^{miss},1.jet)|", 100, 0, 3.5 );
  h["h_dphi_met_j2__"+s] = TH1F( "", ";|#Delta#phi(E_{T}^{miss},2.jet)|", 100, 0, 3.5 );
  h["h_dphi_met_j3__"+s] = TH1F( "", ";|#Delta#phi(E_{T}^{miss},3.jet)|", 100, 0, 3.5 );
  h["h_dphi_met_j13_min__"+s] = TH1F( "", ";min_{i}|#Delta#phi(E_{T}^{miss},i.jet)|", 100, 0, 3.5 );
  h["h_dphi_met_recoil__"+s] = TH1F( "", ";|#Delta#phi(E_{T}^{miss},#Sigma jet)|", 100, 0, 3.5 );

  // multiplicities
  h["h_n_vertex__"+s] = TH1F( "", ";Vertex multiplicity", 61, -0.5, 60.5 );
  h["h_n_photon__"+s] = TH1F( "", ";Photon multiplicity", 6, -0.5, 5.5 );
  h["h_n_jet__"+s] = TH1F( "", ";jet multiplicity", 21, -0.5, 20.5 );
  h["h_n_bjet__"+s] = TH1F( "", ";b-jet multiplicity", 21, -0.5, 20.5 );
  h["h_n_electron__"+s] = TH1F( "", ";electron multiplicity", 4, -0.5, 3.5 );
  h["h_n_muon__"+s] = TH1F( "", ";muon multiplicity", 4, -0.5, 3.5 );

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

  eff["eff_hlt_pt__"+s] = TEfficiency( "", ";p_{T} (GeV);#varepsilon", 200, 0, 2000 );
  eff["eff_hlt_ht__"+s] = TEfficiency( "", ";H_{T} (GeV);#varepsilon", 200, 0, 2000 );
  eff["eff_hlt_nVertex__"+s] = TEfficiency( "", ";Vertex multiplicity", 41, -0.5, 40.5 );
  eff["eff_hlt_sie__"+s] = TEfficiency( "", ";#sigma_{i#etai#eta}", 41, -0.5, 40.5 );
  eff["eff_hlt_hoe__"+s] = TEfficiency( "", ";H/E", 41, -0.5, 40.5 );
  eff["eff_hlt_cIso__"+s] = TEfficiency( "", ";I_{#pi} (GeV)", 41, -0.5, 40.5 );
  eff["eff_hlt_nIso__"+s] = TEfficiency( "", ";I_{n} (GeV)", 41, -0.5, 40.5 );
  eff["eff_hlt_pIso__"+s] = TEfficiency( "", ";I_{#gamma} (GeV)", 41, -0.5, 40.5 );

}

void HistogramProducer::fillSelection( string const& s ) {

  float ht_g = selHt;
  for( auto& p : selPhotons )
    ht_g += p->p.Pt();

  float st = ht_g + met->p.Pt();

  TVector3 recoil(0,0,0);
  for( auto& p : selJets )
    recoil += p->p;

  h["h_met__"+s].Fill( met->p.Pt(), selW );

  h["h_ht__"+s].Fill( selHt, selW );
  h["h_ht_g__"+s].Fill( ht_g, selW );
  h["h_st__"+s].Fill( st, selW );

  if( selPhotons.size() > 0 ) {
    h["h_mt_g_met__"+s].Fill( (selPhotons[0]->p + met->p).Pt(), selW );
    h["h_g_pt__"+s].Fill( selPhotons[0]->p.Pt(), selW );
    h["h_g_eta__"+s].Fill( fabs(selPhotons[0]->p.Eta()), selW );
  }

  if( selJets.size() > 2 ) {
    h["h_j3_pt__"+s].Fill( selJets[2]->p.Pt(), selW );
    h["h_j3_eta__"+s].Fill( fabs(selJets[2]->p.Eta()), selW );
  }
  if( selJets.size() > 1 ) {
    h["h_j2_pt__"+s].Fill( selJets[1]->p.Pt(), selW );
    h["h_j2_eta__"+s].Fill( fabs(selJets[1]->p.Eta()), selW );
  }
  if( selJets.size() > 0 ) {
    h["h_j1_pt__"+s].Fill( selJets[0]->p.Pt(), selW );
    h["h_j1_eta__"+s].Fill( fabs(selJets[0]->p.Eta()), selW );
  }

  if( selBJets.size() > 2 ) {
    h["h_bj3_pt__"+s].Fill( selBJets[2]->p.Pt(), selW );
    h["h_bj3_eta__"+s].Fill( fabs(selBJets[2]->p.Eta()), selW );
  }
  if( selBJets.size() > 1 ) {
    h["h_bj2_pt__"+s].Fill( selBJets[1]->p.Pt(), selW );
    h["h_bj2_eta__"+s].Fill( fabs(selBJets[1]->p.Eta()), selW );
  }
  if( selBJets.size() > 0 ) {
    h["h_bj1_pt__"+s].Fill( selBJets[0]->p.Pt(), selW );
    h["h_bj1_eta__"+s].Fill( fabs(selBJets[0]->p.Eta()), selW );
  }

  if( selPhotons.size() > 0 )
    h["h_dphi_met_g__"+s].Fill( fabs(met->p.DeltaPhi( selPhotons[0]->p )), selW );

  vector<float> minDeltaPhiMetJet;
  if( selJets.size() > 0 ) {
    float dphi = fabs(met->p.DeltaPhi( selJets[0]->p ));
    h["h_dphi_met_j1__"+s].Fill( dphi, selW );
    minDeltaPhiMetJet.push_back( dphi );
  }
  if( selJets.size() > 1 ) {
    float dphi = fabs(met->p.DeltaPhi( selJets[1]->p ));
    h["h_dphi_met_j2__"+s].Fill( dphi, selW );
    minDeltaPhiMetJet.push_back( dphi );
  }
  if( selJets.size() > 2 ) {
    float dphi = fabs(met->p.DeltaPhi( selJets[2]->p ));
    h["h_dphi_met_j3__"+s].Fill( dphi, selW );
    minDeltaPhiMetJet.push_back( dphi );
  }
  if( minDeltaPhiMetJet.size() ) {
    h["h_dphi_met_j13_min__"+s].Fill( *min_element( minDeltaPhiMetJet.begin(), minDeltaPhiMetJet.end() ), selW);
  }

  h["h_dphi_met_recoil__"+s].Fill( fabs(met->p.DeltaPhi( recoil )), selW );

  h["h_n_vertex__"+s].Fill( *nGoodVertices, selW );
  h["h_n_photon__"+s].Fill( selPhotons.size(), selW );
  h["h_n_jet__"+s].Fill( selJets.size(), selW );
  h["h_n_bjet__"+s].Fill( selBJets.size(), selW );
  h["h_n_electron__"+s].Fill( selElectrons.size(), selW );
  h["h_n_muon__"+s].Fill( selMuons.size(), selW );

  if( mrr2.first > 1e-7 ) {
    h2["h2_razorPlane__"+s].Fill( mrr2.first, mrr2.second, selW );
  }

} // end fill

void HistogramProducer::fillObjects( string const& s ) {

  // matching
  for( auto& jet : selJets ) {
    for( auto& p : selPhotons )
      h2["h2_match_jet_photon__"+s].Fill( jet->p.DeltaR( p->p ), jet->p.Pt()/p->p.Pt(), selW );
    for( auto& p : selElectrons )
      h2["h2_match_jet_electron__"+s].Fill( jet->p.DeltaR( p->p ), jet->p.Pt()/p->p.Pt(), selW );
    for( auto& p : selMuons )
      h2["h2_match_jet_muon__"+s].Fill( jet->p.DeltaR( p->p ), jet->p.Pt()/p->p.Pt(), selW );

    h["h_j_looseID__"+s].Fill( float(jet->isLoose) );
  }

  for( auto& photon : selPhotons ) {
    for( auto& p : genParticles ) {
      if( abs(p.pdgId) == 22 )
        h2["h2_match_photon_genPhoton__"+s].Fill( photon->p.DeltaR( p.p ), photon->p.Pt()/p.p.Pt(), selW );
      if( abs(p.pdgId) == 11 )
        h2["h2_match_photon_genElectron__"+s].Fill( photon->p.DeltaR( p.p ), photon->p.Pt()/p.p.Pt(), selW );
    }

    eff["h_gCol_sigmaIetaIeta__"+s].FillWeighted( photon->isTrue, selW, photon->sigmaIetaIeta );
    eff["h_gCol_hOverE__"+s].FillWeighted( photon->isTrue, selW, photon->hOverE );
    eff["h_gCol_r9__"+s].FillWeighted( photon->isTrue, selW, photon->r9 );
    eff["h_gCol_cIso__"+s].FillWeighted( photon->isTrue, selW, photon->isoChargedHadronsEA );
    eff["h_gCol_mva__"+s].FillWeighted( photon->isTrue, selW, photon->mvaValue );

/*    if( photon->p.Pt() < 90 || abs(photon->p.Eta()) > 1.4442 ) continue;
    if(
      photon->hOverE < 0.028
      && photon->sigmaIetaIeta < 0.0107
      && photon->isoNeutralHadronsEA < 7.23 + exp(0.0028*photon->p.Pt()+0.5408)
      && photon->isoPhotonsEA < 2.11 + 0.0014*photon->p.Pt()
    )
    eff["h_g_cIso__"+s].FillWeighted( photon->isTrue, selW, photon->isoChargedHadronsEA );
    if(
      photon->hOverE < 0.028
      && photon->sigmaIetaIeta < 0.0107
      && photon->isoChargedHadronsEA < 2.67
      && photon->isoPhotonsEA < 2.11 + 0.0014*photon->p.Pt()
    )
    eff["h_g_nIso__"+s].FillWeighted( photon->isTrue, selW, photon->isoNeutralHadronsEA );
    if(
      photon->hOverE < 0.028
      && photon->sigmaIetaIeta < 0.0107
      && photon->isoChargedHadronsEA < 2.67
      && photon->isoNeutralHadronsEA < 7.23 + exp(0.0028*photon->p.Pt()+0.5408)
    )
    eff["h_g_nIso__"+s].FillWeighted( photon->isTrue, selW, photon->isoPhotonsEA );
    if(
      photon->hOverE < 0.028
      && photon->isoChargedHadronsEA < 2.67
      && photon->isoNeutralHadronsEA < 7.23 + exp(0.0028*photon->p.Pt()+0.5408)
      && photon->isoPhotonsEA < 2.11 + 0.0014*photon->p.Pt()
    )
    eff["h_g_sigmaIetaIeta__"+s].FillWeighted( photon->isTrue, selW, photon->sigmaIetaIeta );
    if(
      photon->sigmaIetaIeta < 0.0107
      && photon->isoChargedHadronsEA < 2.67
      && photon->isoNeutralHadronsEA < 7.23 + exp(0.0028*photon->p.Pt()+0.5408)
      && photon->isoPhotonsEA < 2.11 + 0.0014*photon->p.Pt()
    )
    eff["h_g_hOverE__"+s].FillWeighted( photon->isTrue, selW, photon->hOverE );

  */
  }

  // trigger efficiencies
  tree::Photon* thisPhoton=0;
  for( auto& photon : selPhotons ) {
    if( photon->p.Pt() > 15 && fabs(photon->p.Eta()) < 1.4442 && photon->isLoose ) {
      thisPhoton = photon;
      break; // take leading photon
    }
  }
  if( thisPhoton ) {

    float ht = 0;
    for( auto& jet : selJets ) {
      if( jet->p.Pt() > 40 && fabs(jet->p.Eta()) < 3 ) {
        ht += jet->p.Pt();
      }
    }

    if( *crossTriggerHt && ht > 650 ) {
      eff["eff_hlt_pt__"+s].Fill( *signalTrigger, thisPhoton->p.Pt() );
    }
    if( *crossTriggerPhoton && thisPhoton->p.Pt() > 100 ) {
      eff["eff_hlt_ht__"+s].Fill( *signalTrigger, ht );
    }
    if( thisPhoton->p.Pt()>100 && ht>600 && *crossTriggerHt ) {
      eff["eff_hlt_nVertex__"+s].Fill( *signalTrigger, *nGoodVertices );
      eff["eff_hlt_hoe__"+s].Fill( *signalTrigger, thisPhoton->hOverE );
      eff["eff_hlt_sie__"+s].Fill( *signalTrigger, thisPhoton->sigmaIetaIeta );
      eff["eff_hlt_cIso__"+s].Fill( *signalTrigger, thisPhoton->isoChargedHadronsEA );
      eff["eff_hlt_nIso__"+s].Fill( *signalTrigger, thisPhoton->isoNeutralHadronsEA );
      eff["eff_hlt_pIso__"+s].Fill( *signalTrigger, thisPhoton->isoPhotonsEA );
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
  signalTrigger( fReader, "HLT_Photon90_CaloIdL_PFHT500_v" ),
  crossTriggerPhoton( fReader, "HLT_Photon90_v" ),
  crossTriggerHt( fReader, "HLT_PFHT600_v" )
{
}

void HistogramProducer::Init(TTree *tree)
{
  fReader.SetTree(tree);
  string inputName = fReader.GetTree()->GetCurrentFile()->GetName();
  isData = inputName.find("SinglePhoton") != string::npos
    || inputName.find("JetHT") != string::npos
    || inputName.find("MET") != string::npos;
}

void HistogramProducer::SlaveBegin(TTree *tree)
{
  initObjects("base");
  h["h_genHt"] = TH1F( "", ";H_{T}^{gen} (GeV)", 6000, 0, 3000 );

  vector<string> strs = { "trBit", "tr", "tr_met200", "tr_genElectron", "tr_genPhoton", "tr_eControl", "tr_jControl" };
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
    if( !jet.isLoose || jet.p.Pt() < 40 || abs(jet.p.Eta()) > 3 ) continue;
    if( indexOfMatchedParticle<tree::Photon*>( jet, selPhotons, .3 ) >= 0 ) continue;
    if( indexOfMatchedParticle<tree::Electron*>( jet, selElectrons, .3 ) >= 0 ) continue;
    if( indexOfMatchedParticle<tree::Muon*>( jet, selMuons, .3 ) >= 0 ) continue;

    selJets.push_back( &jet );
    if( jet.bDiscriminator > bTaggingWorkingPoints["CSVv2M"] )
      selBJets.push_back( &jet );
  }
}

Bool_t HistogramProducer::Process(Long64_t entry)
{
  resetSelection();
  //if( entry > 3 ) return true;
  //if(!( entry%10000 )) printf( "\r%lli / %lli", entry, fReader.GetEntries(false) );
  fReader.SetEntry(entry);

  // set weight
  selW = *mc_weight * *pu_weight;

  h["h_genHt"].Fill( *genHt, selW );

  selHt = 0;
  for( auto& jet : jets ){
    if( jet.p.Pt()>40 && fabs(jet.p.Eta())<3) {
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
  // New selection
  for( auto& photon : photons ) {
    if( photon.p.Pt() > 100 && fabs(photon.p.Eta()) < 1.4442  && !photon.hasPixelSeed && photon.isLoose ) {
       selPhotons.push_back( &photon );
    }
  }
  defaultSelection();

  // Calculate razor variables
  if( true ) { // attention: computing intensive
    vector<TVector3> js;
    for( auto& p : selJets ) js.push_back( p->p );
    for( auto& p : selPhotons ) js.push_back( p->p );
    for( auto& p : selElectrons ) js.push_back( p->p );
    for( auto& p : selMuons ) js.push_back( p->p );
    mrr2 = razorVariables( megajets( js ), met->p );
  }


  if( *signalTrigger ) {
    fillSelection("trBit");
  }

  if( selPhotons.size() && selHt > 600 && (*signalTrigger || !isData) ) {
    fillSelection("tr");
    for( auto& p : genParticles ) {
        if( abs(p.pdgId) == 11 && p.p.DeltaR( selPhotons[0]->p ) < 0.1 ) {
          fillSelection("tr_genElectron");
          break;
        }
    }
    if( selPhotons[0]->isTrue ) {
      fillSelection("tr_genPhoton");
    }
    if( met->p.Pt() > 200 ) {
      fillSelection("tr_met200");
    }

  }


  resetSelection();
  /////////////////////////////////////////////////////////////////////////////
  // electron sample

  for( auto& photon : photons ) {
    if( photon.p.Pt() > 100 && fabs(photon.p.Eta()) < 1.4442  && photon.hasPixelSeed && photon.isLoose ) {
      selPhotons.push_back( &photon );
    }
  }
  defaultSelection();

  if( selPhotons.size() && selHt > 600 && (*signalTrigger || !isData) ) {
    fillSelection("tr_eControl");
  }

  resetSelection();
  /////////////////////////////////////////////////////////////////////////////
  // jet sample
  for( auto& photon : photons ) {
    if( photon.p.Pt() > 100 && fabs(photon.p.Eta()) < 1.4442  && !photon.hasPixelSeed
        && !photon.isLoose
        && photon.hOverE < 0.05
        && photon.sigmaIetaIeta < 0.0103
    ) {
      selPhotons.push_back( &photon );
    }
  }
  defaultSelection();
  if( selPhotons.size() && selHt > 600 && (*signalTrigger || !isData) ) {
     fillSelection("tr_jControl");
  }

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


  fReader.GetTree()->GetCurrentFile()->Get("TreeWriter/hCutFlow")->Write();

}

void HistogramProducer::resetSelection() {
  selPhotons.clear();
  selJets.clear();
  selBJets.clear();
  selElectrons.clear();
  selMuons.clear();
}
