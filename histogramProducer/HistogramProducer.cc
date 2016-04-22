#include "regex"
#include "time.h"

#include "TROOT.h"
#include "TFile.h"
#include "TSelector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH2F.h"
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
  void fillSelection( string const& s );
  void fillSelectedPhotons( const tree::Particle& p );
  tree::Jet* matchedJet( const tree::Particle& p );

  void initTriggerStudies();
  void fillTriggerStudies();

  int genMatch( const tree::Particle& p ) ;


  TTreeReader fReader;
  TTreeReaderValue<std::vector<tree::Photon>> photons;
  TTreeReaderValue<std::vector<tree::Jet>> jets;
  TTreeReaderValue<std::vector<tree::Electron>> electrons;
  TTreeReaderValue<std::vector<tree::Muon>> muons;
  TTreeReaderValue<std::vector<tree::Particle>> genJets;
  TTreeReaderValue<std::vector<tree::GenParticle>> genParticles;
  TTreeReaderValue<std::vector<tree::IntermediateGenParticle>> intermediateGenParticles;
  TTreeReaderValue<tree::MET> met;
  TTreeReaderValue<Float_t> pu_weight;
  TTreeReaderValue<Char_t> mc_weight;
  TTreeReaderValue<Int_t> nGoodVertices;
  TTreeReaderValue<Float_t> genHt;
  TTreeReaderValue<ULong64_t> eventNo;
  TTreeReaderValue<UInt_t> runNo;
  TTreeReaderValue<UInt_t> lumNo;
  TTreeReaderValue<Bool_t> hlt_photon90_ht500;
  TTreeReaderValue<Bool_t> hlt_photon90;
  TTreeReaderValue<Bool_t> hlt_ht600;
  TTreeReaderValue<Int_t> hlt_ht600_pre;

  vector<tree::Photon*> selPhotons;
  vector<tree::Jet*> selJets;
  vector<tree::Jet*> selBJets;
  vector<tree::Jet*> selHEJets;
  vector<tree::Electron*> selElectrons;
  vector<tree::Muon*> selMuons;

  vector<tree::Photon> artificialPhotons;

  float selW=1.; // weight

  map<string,map<string,TH1F>> hMapMap;
  map<string,map<string,TH2F>> hMapMap2;
  map<string,TEfficiency> effMap;
  map<int,pair<int,int>> rawEff_vs_run;

  bool isData;
  bool zToMet;
  int gluinoMass;

  JetSelector jetSelector;
  Weighter nJetWeighter;
  CutFlowPhoton looseCutFlowPhoton;

  double startTime;
};

int HistogramProducer::genMatch( const tree::Particle& p ) {

  // match to daghters of massive particles
  for(auto const& genP :*intermediateGenParticles ) {
    for( auto const & d : genP.daughters ) {
      auto id = fabs(d.pdgId);
      auto dr = p.p.DeltaR(d.p);
      auto dpt = p.p.Pt()/d.p.Pt();
      if( dr < 0.5 ) return id;
    }
  }

  for( auto const& genP : *genParticles ) {
    auto id = fabs(genP.pdgId);
    auto dr = p.p.DeltaR(genP.p);
    auto dpt = p.p.Pt()/genP.p.Pt();
    if( dr<0.15 && fabs(dpt-1) < 0.15) {
      if( genP.isPrompt ) return id;
      else return -id;
    }
  }

  return 0;
}

void HistogramProducer::initTriggerStudies() {
  effMap["eff_hlt_pt"] = TEfficiency( "", ";p_{T} (GeV);#varepsilon", 250, 0, 1000 );
  effMap["eff_hlt_eta"] = TEfficiency( "", ";|#eta|;#varepsilon", 15, 0, 1.5 );
  effMap["eff_hlt_ht"] = TEfficiency( "", ";EMH_{T} (GeV);#varepsilon", 200, 0, 2000 );
  effMap["eff_hlt_ht_ptMin"] = TEfficiency( "", ";EMH_{T} (GeV);#varepsilon", 200, 0, 2000 );
  effMap["eff_hlt_ht_etaMax"] = TEfficiency( "", ";EMH_{T} (GeV);#varepsilon", 310, 0, 3.1 );
  effMap["eff_hlt_ht_ct"] = TEfficiency( "", ";EMH_{T} (GeV);#varepsilon", 200, 0, 2000 );
  effMap["eff_hlt_ht_ct_preScaled"] = TEfficiency( "", ";EMH_{T} (GeV);#varepsilon (prescaled)", 200, 0, 2000 );
  effMap["eff_hlt_ht_ct2"] = TEfficiency( "", ";EMH_{T} (GeV);#varepsilon", 200, 0, 2000 );
  effMap["eff_hlt_ht_ct2_preScaled"] = TEfficiency( "", ";EMH_{T} (GeV);#varepsilon (prescaled)", 200, 0, 2000 );
  effMap["eff_hlt_nVertex"] = TEfficiency( "", ";Vertex multiplicity", 41, -0.5, 40.5 );
  effMap["eff_hlt_sie"] = TEfficiency( "", ";#sigma_{i#etai#eta}", 400, 0, 0.02 );
  effMap["eff_hlt_hoe"] = TEfficiency( "", ";H/E", 100, 0, 0.15 );
  effMap["eff_hlt_r9"] = TEfficiency( "", ";r9", 110, 0, 1.1 );
  effMap["eff_hlt_cIso"] = TEfficiency( "", ";I_{#pi} (GeV)", 100, 0, 10 );
  effMap["eff_hlt_nIso"] = TEfficiency( "", ";I_{n} (GeV)", 100, 0, 20 );
  effMap["eff_hlt_pIso"] = TEfficiency( "", ";I_{#gamma} (GeV)", 100, 0, 20 );
  effMap["eff_hlt_nJet"] = TEfficiency( "", ";uncleaned jet multiplicity", 15, -0.5, 14.5 );
  effMap["eff_hlt_met"] = TEfficiency( "", ";E_{T}^{miss} (GeV)", 150, 0, 150 );
  effMap["eff_hlt_met_ct"] = TEfficiency( "", ";E_{T}^{miss} (GeV)", 150, 0, 150 );
  effMap["eff_hlt_met_ct2"] = TEfficiency( "", ";E_{T}^{miss} (GeV)", 150, 0, 150 );
}

void HistogramProducer::fillTriggerStudies() {

  tree::Photon* thisPhoton=0;
  for( auto& photon : *photons ) {
    if( photon.p.Pt() > 15 && fabs(photon.p.Eta()) < photonsEtaMaxBarrel && photon.isLoose && !photon.hasPixelSeed ) {
      thisPhoton = &photon;
      break; // take leading(first) photon
    }
  }

  float ht = 0;
  if( thisPhoton ) ht += thisPhoton->p.Pt();
  tree::Jet *minPtJet=0;
  tree::Jet *maxEtaJet=0;
  for( auto& jet : *jets ) {
    if( jet.p.Pt() > 40 && fabs(jet.p.Eta()) < 3 ) {
      if( !minPtJet || minPtJet->p.Pt() > jet.p.Pt() ) minPtJet = &jet;
      if( !maxEtaJet || fabs(minPtJet->p.Eta()) < fabs(jet.p.Eta()) ) maxEtaJet = &jet;
      if( !thisPhoton || jet.p.DeltaR(thisPhoton->p) > 0.3 ) {
        ht += jet.p.Pt();
      }
    }
  }

  if( thisPhoton ) {

    // get trigger efficiencies for each run
    if( thisPhoton->p.Pt()>100  && ht>700 && *hlt_ht600) {
      if( !rawEff_vs_run.count( *runNo ) ) rawEff_vs_run[*runNo] = make_pair(0,0);
      if( *hlt_photon90_ht500 ) rawEff_vs_run.at(*runNo).first += 1;
      rawEff_vs_run.at(*runNo).second += 1;
    }

    // main trigger plots
    if( thisPhoton->p.Pt()>100 && *hlt_photon90 ) {
      effMap.at("eff_hlt_ht").Fill( *hlt_photon90_ht500, ht );
      effMap.at("eff_hlt_ht_ptMin").Fill( *hlt_photon90_ht500, minPtJet->p.Pt() );
      if( ht>600 ) effMap.at("eff_hlt_ht_etaMax").Fill( *hlt_photon90_ht500, fabs(maxEtaJet->p.Eta()) );
    }
    if( thisPhoton->p.Pt()>100 && *hlt_photon90_ht500 ) {
      effMap.at("eff_hlt_ht_ct").Fill( *hlt_ht600, ht );
      effMap.at("eff_hlt_ht_ct_preScaled").FillWeighted( *hlt_ht600, *hlt_ht600_pre, ht );
      if( ht > 700 ) effMap.at("eff_hlt_met_ct").Fill( *hlt_ht600, met->p.Pt() );
    }
    if( thisPhoton->p.Pt()>100 && *hlt_photon90 ) {
      effMap.at("eff_hlt_ht_ct2").Fill( *hlt_ht600, ht );
      effMap.at("eff_hlt_ht_ct2_preScaled").FillWeighted( *hlt_ht600, *hlt_ht600_pre, ht );
      if( ht > 700 ) effMap.at("eff_hlt_met_ct2").Fill( *hlt_ht600, met->p.Pt() );
    }
    if( *hlt_ht600 && ht > 700 ) {
      effMap.at("eff_hlt_pt").Fill( *hlt_photon90_ht500, thisPhoton->p.Pt() );
      effMap.at("eff_hlt_eta").Fill( *hlt_photon90_ht500, fabs(thisPhoton->p.Eta()) );
    }
  }

  if( *hlt_ht600 && ht > 700 ) {
    for( auto& photon : *photons ) {
      if( photon.p.Pt() > 100 && fabs(photon.p.Eta()) < photonsEtaMaxBarrel ) {
        if( looseCutFlowPhoton.check( photon ) ) {
          effMap.at("eff_hlt_r9").Fill( *hlt_photon90_ht500, photon.r9 );
          effMap.at("eff_hlt_nVertex").Fill( *hlt_photon90_ht500, *nGoodVertices );
          effMap.at("eff_hlt_nJet").Fill( *hlt_photon90_ht500, jets->size() );
          effMap.at("eff_hlt_met").Fill( *hlt_photon90_ht500, met->p.Pt() );
        }
        if(
          looseCutFlowPhoton.passHoe()
          && looseCutFlowPhoton.passSie()
          && looseCutFlowPhoton.passNIso()
          && looseCutFlowPhoton.passPIso()
          ) effMap.at("eff_hlt_cIso").Fill( *hlt_photon90_ht500, photon.cIso );
        if(
          looseCutFlowPhoton.passHoe()
          && looseCutFlowPhoton.passSie()
          && looseCutFlowPhoton.passCIso()
          && looseCutFlowPhoton.passPIso()
          ) effMap.at("eff_hlt_nIso").Fill( *hlt_photon90_ht500, photon.nIso );
        if(
          looseCutFlowPhoton.passHoe()
          && looseCutFlowPhoton.passSie()
          && looseCutFlowPhoton.passCIso()
          && looseCutFlowPhoton.passNIso()
          ) effMap.at("eff_hlt_pIso").Fill( *hlt_photon90_ht500, photon.pIso );
        if(
          looseCutFlowPhoton.passHoe()
          && looseCutFlowPhoton.passIso()
        ) effMap.at("eff_hlt_sie").Fill( *hlt_photon90_ht500, photon.sigmaIetaIeta );
        if(
          looseCutFlowPhoton.passSie()
          && looseCutFlowPhoton.passIso()
        ) effMap.at("eff_hlt_hoe").Fill( *hlt_photon90_ht500, photon.hOverE );
      }
    } // photon loop
  } // ht cross trigger
}

map<string,TH2F> initHistograms2(){
  map<string,TH2F> hMap;

  hMap["n_heJets_vs_photonPosition"] = TH2F("","",10, -0.5, 9.5, 10, -0.5, 9.5 );
  hMap["g_eta_vs_g_phi"] = TH2F("","",100, -1.5, 1.5, 100, -3.1, 3.1 );
  hMap["met_vs_emht"] = TH2F("", ";E_{T}^{miss} (GeV);EMH_{T} (GeV)", 500, 0, 5000, 500, 0, 5000 );
  hMap["metRaw_vs_emht"] = TH2F("", ";uncorrected E_{T}^{miss} (GeV);EMH_{T} (GeV)", 500, 0, 5000, 500, 0, 5000 );

  return hMap;
}

map<string,TH1F> initHistograms(){
  map<string,TH1F> hMap;

  hMap["met"] = TH1F( "", ";E^{miss}_{T} (GeV)", 200, 0, 2000 );
  hMap["metStar"] = TH1F( "", ";E^{miss}_{T}* (GeV)", 200, 0, 2000 );
  hMap["metStar2"] = TH1F( "", ";E^{miss}_{T}* (GeV)", 200, 0, 2000 );
  hMap["metUp"] = TH1F( "", ";E^{miss}_{T} up (GeV)", 200, 0, 2000 );
  hMap["metDn"] = TH1F( "", ";E^{miss}_{T} down (GeV)", 200, 0, 2000 );
  hMap["metUpJec"] = TH1F( "", ";E^{miss}_{T} up (GeV)", 200, 0, 2000 );
  hMap["metDnJec"] = TH1F( "", ";E^{miss}_{T} down (GeV)", 200, 0, 2000 );
  hMap["metPar"] = TH1F( "", ";E^{miss}_{T} #parallel (GeV)", 400, -2000, 2000 );
  hMap["metPer"] = TH1F( "", ";E^{miss #perp  }_{T} (GeV)", 200, 0, 2000 );
  hMap["metRaw"] = TH1F( "", ";uncorrected E^{miss}_{T} (GeV)", 200, 0, 2000 );
  hMap["metParRaw"] = TH1F( "", ";uncorrected E^{miss}_{T} #parallel (GeV)", 400, -2000, 2000 );
  hMap["metPerRaw"] = TH1F( "", ";uncorrected E^{miss #perp  }_{T} (GeV)", 200, 0, 2000 );
  hMap["mt_g_met"] = TH1F( "", ";m_{T}(#gamma,E^{miss}_{T}) (GeV)", 150, 0, 1500 );

  hMap["tremht"] = TH1F( "", ";EMH_{T}^{trigger-like} (GeV)", 300, 0, 3000 );
  hMap["emht"] = TH1F( "", ";EMH_{T} (GeV)", 300, 0, 3000 );
  hMap["emhtStar"] = TH1F( "", ";EMH_{T}* (GeV)", 300, 0, 3000 );
  hMap["ht"] = TH1F( "", ";H_{T} (GeV)", 300, 0, 3000 );
  hMap["st"] = TH1F( "", ";S_{T} (GeV)", 200, 0, 4000 );
  hMap["emrecoilt"] = TH1F( "", ";#vec{EMH}_{T} (GeV)", 150, 0, 1500 );
  hMap["recoilt"] = TH1F( "", ";#vec{H}_{T} (GeV)", 150, 0, 1500 );

  // photon
  hMap["g_pt"] = TH1F( "", ";p_{T} (GeV)", 150, 0, 1500 );
  hMap["g_ptStar"] = TH1F( "", ";p_{T}* (GeV)", 150, 0, 1500 );
  hMap["g_eta"] = TH1F( "", ";|#eta|", 1500, 0, 1.5 );

  // he-jet
  hMap["hej1_pt"] = TH1F( "", ";p_{T}^{1.jet} (GeV)", 150, 0, 1500 );
  hMap["hej1_eta"] = TH1F( "", ";|#eta^{1.jet}|", 150, 0, 3 );
  hMap["hej2_pt"] = TH1F( "", ";p_{T}^{2.jet} (GeV)", 150, 0, 1500 );
  hMap["hej2_eta"] = TH1F( "", ";|#eta^{2.jet}|", 150, 0, 3 );
  hMap["hej3_pt"] = TH1F( "", ";p_{T}^{3.jet} (GeV)", 150, 0, 1500 );
  hMap["hej3_eta"] = TH1F( "", ";|#eta^{3.jet}|", 150, 0, 3 );

  // angles
  hMap["dphi_met_g"] = TH1F( "", ";|#Delta#phi(E_{T}^{miss},#gamma)|", 70, 0, 3.5 );
  hMap["dphi_met_hej1"] = TH1F( "", ";|#Delta#phi(E_{T}^{miss},1.jet)|", 70, 0, 3.5 );
  hMap["dphi_met_hej2"] = TH1F( "", ";|#Delta#phi(E_{T}^{miss},2.jet)|", 70, 0, 3.5 );
  hMap["dphi_met_hej3"] = TH1F( "", ";|#Delta#phi(E_{T}^{miss},3.jet)|", 70, 0, 3.5 );
  hMap["dphi_met_recoil"] = TH1F( "", ";|#Delta#phi(E_{T}^{miss},#Sigma jet)|", 70, 0, 3.5 );
  hMap["dphi_g_hej1"] = TH1F( "", ";|#Delta#phi(#gamma,1.jet)|", 70, 0, 3.5 );
  hMap["dphi_g_hej2"] = TH1F( "", ";|#Delta#phi(#gamma,2.jet)|", 70, 0, 3.5 );

  // multiplicities
  hMap["n_vertex"] = TH1F( "", ";vertex multiplicity", 61, -0.5, 60.5 );
  hMap["n_vertex_unw"] = TH1F( "", ";vertex multiplicity", 61, -0.5, 60.5 );
  hMap["n_photon"] = TH1F( "", ";photon multiplicity", 6, -0.5, 5.5 );
  hMap["n_jet"] = TH1F( "", ";jet multiplicity", 16, -0.5, 15.5 );
  hMap["n_bjet"] = TH1F( "", ";b-jet multiplicity", 16, -0.5, 15.5 );
  hMap["n_electron"] = TH1F( "", ";electron multiplicity", 4, -0.5, 3.5 );
  hMap["n_muon"] = TH1F( "", ";muon multiplicity", 4, -0.5, 3.5 );
  hMap["n_heJet"] = TH1F( "", ";photon-like jet multiplicity", 11, -0.5, 10.5 );

  hMap["genMatch"] = TH1F( "", ";pdg id for gen match", 24, -0.5, 23.5 );
  hMap["genHt"] = TH1F( "", ";H_{T}^{gen}", 3000, 0, 3000 );

  return hMap;
}

void HistogramProducer::fillSelection( string const& s ) {

  if( !hMapMap.count(s) ) {
    hMapMap[s] = initHistograms();
    hMapMap2[s] = initHistograms2();
  }

  // calculate variables
  float ht = 0;
  TVector3 recoil(0,0,0);
  for( auto& j: selJets ){
    ht += j->p.Pt();
    recoil += j->p;
  }
  TVector3 emrecoil = recoil;
  float ht_g = ht;
  float ht_gStar = ht;
  for( auto& p : selPhotons ) {
    emrecoil += p->p;
    ht_g += p->p.Pt();
    auto mJet = matchedJet(*p);
    if(mJet) ht_gStar += mJet->p.Pt();
    else ht_gStar = p->p.Pt();
  }
  float st = ht_g + met->p.Pt();

  float tremht = 0;
  for( auto& jet : *jets ){
    auto pt = jet.p.Pt();
    if( pt>40 && fabs(jet.p.Eta())<3 ) {
      tremht += pt;
    }
  }

  hMapMap.at(s).at("tremht").Fill( tremht, selW );
  hMapMap.at(s).at("emht").Fill( ht_g, selW );
  hMapMap.at(s).at("emhtStar").Fill( ht_gStar, selW );
  hMapMap.at(s).at("ht").Fill( ht, selW );
  hMapMap.at(s).at("st").Fill( st, selW );
  hMapMap.at(s).at("recoilt").Fill( recoil.Pt(), selW );
  hMapMap.at(s).at("emrecoilt").Fill( emrecoil.Pt(), selW );

  hMapMap.at(s).at("met").Fill( met->p.Pt(), selW );
  hMapMap.at(s).at("metUp").Fill( met->p.Pt()+met->uncertainty, selW );
  hMapMap.at(s).at("metDn").Fill( met->p.Pt()-met->uncertainty, selW );
  hMapMap.at(s).at("metRaw").Fill( met->p_raw.Pt(), selW );

  if( selPhotons.size() > 0 ) {
    hMapMap.at(s).at("genMatch").Fill( genMatch(*selPhotons.at(0)), selW );
    auto mJet = matchedJet(*selPhotons.at(0));
    if(mJet){
      hMapMap.at(s).at("metStar").Fill( (met->p+selPhotons.at(0)->p-mJet->p).Pt(), selW );
      hMapMap.at(s).at("metStar2").Fill( (met->p-selPhotons.at(0)->p+mJet->p).Pt(), selW );
      hMapMap.at(s).at("metUpJec").Fill( (met->p+(mJet->p*mJet->uncert)).Pt(), selW );
      hMapMap.at(s).at("metDnJec").Fill( (met->p-(mJet->p*mJet->uncert)).Pt(), selW );
      hMapMap.at(s).at("g_ptStar").Fill( mJet->p.Pt(), selW );
    }

    hMapMap.at(s).at("mt_g_met").Fill( (selPhotons.at(0)->p + met->p).Pt(), selW );
    hMapMap.at(s).at("g_pt").Fill( selPhotons.at(0)->p.Pt(), selW );
    hMapMap.at(s).at("g_eta").Fill( fabs(selPhotons.at(0)->p.Eta()), selW );
    float dphi_met_g = fabs(met->p.DeltaPhi( selPhotons.at(0)->p ));
    hMapMap.at(s).at("dphi_met_g").Fill( dphi_met_g, selW );
    hMapMap.at(s).at("metPar").Fill( met->p.Pt()*cos(dphi_met_g), selW );
    hMapMap.at(s).at("metPer").Fill( fabs(met->p.Pt()*sin(dphi_met_g)), selW );
    hMapMap.at(s).at("metParRaw").Fill( met->p_raw.Pt()*cos(met->p_raw.DeltaPhi(selPhotons.at(0)->p)), selW );
    hMapMap.at(s).at("metPerRaw").Fill( fabs(met->p_raw.Pt()*sin(met->p_raw.DeltaPhi(selPhotons.at(0)->p))), selW );
    unsigned photonPosition=0;
    for(;photonPosition<selHEJets.size() && selHEJets.at(photonPosition)->p.Pt() > selPhotons.at(0)->p.Pt();photonPosition++);
    hMapMap2.at(s).at("n_heJets_vs_photonPosition").Fill( selHEJets.size(), photonPosition, selW );
    hMapMap2.at(s).at("g_eta_vs_g_phi").Fill( selPhotons.at(0)->p.Eta(), selPhotons.at(0)->p.Phi(), selW );
  }

  if( selHEJets.size() > 2 ) {
    hMapMap.at(s).at("hej3_pt").Fill( selHEJets.at(2)->p.Pt(), selW );
    hMapMap.at(s).at("hej3_eta").Fill( fabs(selHEJets.at(2)->p.Eta()), selW );
    hMapMap.at(s).at("dphi_met_hej3").Fill( fabs(met->p.DeltaPhi( selHEJets.at(2)->p )), selW );
  }
  if( selHEJets.size() > 1 ) {
    hMapMap.at(s).at("hej2_pt").Fill( selHEJets.at(1)->p.Pt(), selW );
    hMapMap.at(s).at("hej2_eta").Fill( fabs(selHEJets.at(1)->p.Eta()), selW );
    hMapMap.at(s).at("dphi_met_hej2").Fill( fabs(met->p.DeltaPhi( selHEJets.at(1)->p )), selW );
    if( selPhotons.size() > 0 ) hMapMap.at(s).at("dphi_g_hej2").Fill( fabs(selPhotons.at(0)->p.DeltaPhi(selHEJets.at(1)->p)), selW );
  }
  if( selHEJets.size() > 0 ) {
    hMapMap.at(s).at("hej1_pt").Fill( selHEJets.at(0)->p.Pt(), selW );
    hMapMap.at(s).at("hej1_eta").Fill( fabs(selHEJets.at(0)->p.Eta()), selW );
    hMapMap.at(s).at("dphi_met_hej1").Fill( fabs(met->p.DeltaPhi( selHEJets.at(0)->p )), selW );
    if( selPhotons.size() > 0 ) hMapMap.at(s).at("dphi_g_hej1").Fill( fabs(selPhotons.at(0)->p.DeltaPhi(selHEJets.at(0)->p)), selW );
  }

  hMapMap.at(s).at("dphi_met_recoil").Fill( fabs(met->p.DeltaPhi( recoil )), selW );

  hMapMap.at(s).at("n_vertex").Fill( *nGoodVertices, selW );
  hMapMap.at(s).at("n_vertex_unw").Fill( *nGoodVertices );
  hMapMap.at(s).at("n_photon").Fill( selPhotons.size(), selW );
  hMapMap.at(s).at("n_jet").Fill( selJets.size(), selW );
  hMapMap.at(s).at("n_bjet").Fill( selBJets.size(), selW );
  hMapMap.at(s).at("n_electron").Fill( selElectrons.size(), selW );
  hMapMap.at(s).at("n_muon").Fill( selMuons.size(), selW );
  hMapMap.at(s).at("n_heJet").Fill( selHEJets.size(), selW );
  hMapMap.at(s).at("genHt").Fill( *genHt, selW );

  hMapMap2.at(s).at("met_vs_emht").Fill( met->p.Pt(), ht_g, selW );
  hMapMap2.at(s).at("metRaw_vs_emht").Fill( met->p_raw.Pt(), ht_g, selW );

} // end fill


HistogramProducer::HistogramProducer():
  photons( fReader, "photons" ),
  jets( fReader, "jets" ),
  electrons( fReader, "electrons" ),
  muons( fReader, "muons" ),
  genJets( fReader, "genJets" ),
  genParticles( fReader, "genParticles" ),
  intermediateGenParticles( fReader, "intermediateGenParticles" ),
  met( fReader, "met" ),
  nGoodVertices( fReader, "nGoodVertices" ),
  pu_weight( fReader, "pu_weight" ),
  mc_weight( fReader, "mc_weight" ),
  genHt( fReader, "genHt" ),
  runNo( fReader, "runNo" ),
  hlt_photon90_ht500( fReader, "HLT_Photon90_CaloIdL_PFHT500_v" ),
  hlt_photon90( fReader, "HLT_Photon90_v" ),
  hlt_ht600( fReader, "HLT_PFHT600_v" ),
  hlt_ht600_pre( fReader, "HLT_PFHT600_v_pre" ),
  jetSelector("../plotter/gammaPosition.root", "gqcd"),
  nJetWeighter("../plotter/weights.root", "weight_n_heJet"),
  looseCutFlowPhoton( {{"sigmaIetaIeta_eb",0.0102}, {"cIso_eb",3.32}, {"nIso1_eb",1.92}, {"nIso2_eb",0.014}, {"nIso3_eb",0.000019}, {"pIso1_eb",0.81}, {"pIso2_eb",0.0053},
    {"sigmaIetaIeta_ee",0.0274}, {"cIso_ee",1.97}, {"nIso1_ee",11.86}, {"nIso2_ee",0.0139}, {"nIso3_ee",0.000025}, {"pIso1_ee",0.83}, {"pIso2_ee",0.0034} }),
  gluinoMass(0),
  startTime(time(NULL))
{
}

void HistogramProducer::Init(TTree *tree)
{
  fReader.SetTree(tree);
  string inputName = fReader.GetTree()->GetCurrentFile()->GetName();
  isData = inputName.find("Run201") != string::npos;
  zToMet = inputName.find("ZGTo2LG") != string::npos;

  std::smatch sm;
  if( regex_match( inputName, sm, regex(".*/T5.*_(\\d+)_(\\d+).root") ) ) {
      gluinoMass = stoi(sm[1]);
  }
  initTriggerStudies();
}

void HistogramProducer::SlaveBegin(TTree *tree)
{
}

void HistogramProducer::defaultSelection()
{
  for( auto& mu : *muons ) {
    if( mu.p.Pt() < 15 ) continue;
    if( indexOfMatchedParticle<tree::Photon*>( mu, selPhotons, .3 ) >= 0 ) continue;
    selMuons.push_back( &mu );
  }
  for( auto& el : *electrons ) {
    if( !el.isLoose || el.p.Pt() < 15 ) continue;
    if( indexOfMatchedParticle<tree::Photon*>( el, selPhotons, .3 ) >= 0 ) continue;
    selElectrons.push_back( &el );
  }
  for( auto& jet : *jets ) {
    if( !jet.isLoose
      || jet.hasPhotonMatch || jet.hasElectronMatch || jet.hasMuonMatch
      || jet.p.Pt() < 40 || fabs(jet.p.Eta()) > 3 ) continue;

    selJets.push_back( &jet );
    if( jet.p.Pt() > 100 && fabs(jet.p.Eta()) < photonsEtaMaxBarrel ) {
      selHEJets.push_back( &jet );
    }
    if( jet.bDiscriminator > bTaggingWorkingPoints.at(CSVv2M) && fabs(jet.p.Eta())<2.5 )
      selBJets.push_back( &jet );
  }
}

void HistogramProducer::fillSelectedPhotons( const tree::Particle& p ) {
  tree::Photon photon;
  photon.p = p.p;
  artificialPhotons.push_back( photon );
  selPhotons.push_back( &artificialPhotons.back() );
}

tree::Jet* HistogramProducer::matchedJet( const tree::Particle& p ) {
  for( auto& jet : *jets ) {
    if( jet.p.DeltaR(p.p)<.2 ) {
      return &jet;
    }
  }
  return 0;
}


Bool_t HistogramProducer::Process(Long64_t entry)
{
  resetSelection();
  fReader.SetEntry(entry);

  float zToMetPt = -1;
  if( zToMet && intermediateGenParticles->size() ) {
    auto z = intermediateGenParticles->at(0);
    if( z.daughters.size() > 2 ) {
      return kTRUE;
    }
    met->p += z.p;
    for( auto d : z.daughters ) {
      for( auto p = photons->begin(); p!=photons->end(); ) {
        if (p->p.DeltaR(d.p)<0.1) p = photons->erase(p);
        else p++;
      }
      for( auto p = jets->begin(); p!=jets->end(); ) {
        if (p->p.DeltaR(d.p)<0.1) p = jets->erase(p);
        else p++;
      }
      for( auto p = electrons->begin(); p!=electrons->end(); ) {
        if (p->p.DeltaR(d.p)<0.1) p = electrons->erase(p);
        else p++;
      }
      for( auto p = muons->begin(); p!=muons->end(); ) {
        if (p->p.DeltaR(d.p)<0.1) p = muons->erase(p);
        else p++;
      }
    }
    for( auto p : *genParticles ) {
      if( p.pdgId==22 ) {
        zToMetPt = p.p.Pt();
        break;
      }
    }
  }
  if (zToMet && zToMetPt<0) return kTRUE;


  if(isData) fillTriggerStudies();

  // set weight
  selW = *mc_weight * *pu_weight;
  float originalW = selW;


  // For FastSim CMSSW7X, there are events with large pt jets (pt>sqrt(2))
  // It was recommended for the 2015 analysis to ignore events with objets > 2 times gluino mass
  // TODO: remove these events from the acceptance, change nGen (should be small influence)
  // see here: https://hypernews.cern.ch/HyperNews/CMS/get/met/432.html
  if( gluinoMass > 0 && jets->size() && jets->at(0).p.Pt() > 2*gluinoMass ) return kTRUE;


  /////////////////////////////////////////////////////////////////////////////
  // signal sample
  /////////////////////////////////////////////////////////////////////////////

  for( auto& photon : *photons ) {
    if( photon.isLoose && !photon.hasPixelSeed && photon.p.Pt() > 100 && fabs(photon.p.Eta()) < photonsEtaMaxBarrel ) {
      selPhotons.push_back( &photon );
    }
  }
  defaultSelection();

  float myHt=0;
  for( auto& p : selPhotons ) myHt += p->p.Pt();
  for( auto& p : selJets ) myHt += p->p.Pt();

  if( selPhotons.size() && myHt > 700 && (*hlt_photon90_ht500 || !isData) ) {
    fillSelection("tr");
    if(zToMet&&zToMetPt>130) fillSelection("tr_130pt");
    if(zToMet&&zToMetPt<130) fillSelection("tr_0pt130");
    fillSelection("tr");
    if(myHt>2500) fillSelection("tr_highHt");
    else fillSelection("tr_lowHt");
    auto gMatch = genMatch(*selPhotons.at(0));
    if( gMatch == 11 ) fillSelection("tr_genE");
    else if (gMatch<11 || gMatch>16) fillSelection("tr_noGenLep");
    if(met->p.Pt()>300) fillSelection("tr_highMet");
    else if(met->p.Pt()>70) fillSelection("tr_mediumMet");
    else fillSelection("tr_lowMet");
    if(selHEJets.size()==0) fillSelection("tr_he0");
    if(selHEJets.size()==1) fillSelection("tr_he1");
    if(selHEJets.size()==2) fillSelection("tr_he2");
    if(selHEJets.size()==3) fillSelection("tr_he3");
  }

  if( !selPhotons.size() && selHEJets.size() && myHt > 700 && (*hlt_ht600 || !isData) ) {
    selW *= *hlt_ht600_pre;
    auto indexToRemove = jetSelector.getJetN(selHEJets.size()-1); // -1, since one jet is erased from list (gets photon proxy)
    fillSelectedPhotons( *selHEJets.at(indexToRemove) );
    // clean jets and other he objects
    selHEJets.erase( selHEJets.begin()+indexToRemove );
    for( unsigned i=0; i<selJets.size();){
      if( selPhotons.at(0)->p.DeltaR(selJets.at(i)->p)<0.3 ) selJets.erase(selJets.begin()+i);
      else i++;
    }
    fillSelection("tr_jControl");
    if(myHt>2500) fillSelection("tr_jControl_highHt");
    else fillSelection("tr_jControl_lowHt");
    if(selHEJets.size()==0) fillSelection("tr_jControl_he0");
    if(selHEJets.size()==1) fillSelection("tr_jControl_he1");
    if(selHEJets.size()==2) fillSelection("tr_jControl_he2");
    if(selHEJets.size()==3) fillSelection("tr_jControl_he3");
    float wnJet=nJetWeighter.getWeight(selHEJets.size());
    selW = originalW * wnJet * *hlt_ht600_pre;
    fillSelection("tr_jControl_wnjet");
    if(met->p.Pt()>300) fillSelection("tr_jControl_highMet");
    else if(met->p.Pt()>70) fillSelection("tr_jControl_mediumMet");
    else fillSelection("tr_jControl_lowMet");
    selW = originalW;
  }

  resetSelection();
  /////////////////////////////////////////////////////////////////////////////
  // electron sample
  /////////////////////////////////////////////////////////////////////////////

  for( auto& photon : *photons ) {
    if( photon.isLoose && photon.hasPixelSeed && photon.p.Pt() > 100 && fabs(photon.p.Eta()) < photonsEtaMaxBarrel ) {
      selPhotons.push_back( &photon );
    }
  }
  defaultSelection();

  myHt=0;
  selW = originalW;
  for( auto& p : selPhotons ) myHt += p->p.Pt();
  for( auto& p : selJets ) myHt += p->p.Pt();
  if( selPhotons.size() && myHt > 700 && (*hlt_photon90_ht500 || !isData) ) {
    fillSelection("tr_eControl");
    if( genMatch(*selPhotons.at(0)) == 11 ) fillSelection("tr_eControl_genE");
  }


  return kTRUE;
}

void HistogramProducer::Terminate()
{

  auto outputName = getOutputFilename( fReader.GetTree()->GetCurrentFile()->GetName() );

  // save all defined histograms to file
  TFile file( outputName.c_str(), "RECREATE");

  for( auto& hMapIt : hMapMap ) {
    file.mkdir(hMapIt.first.c_str());
    file.cd(hMapIt.first.c_str());
    for( auto& h : hMapIt.second ) {
      h.second.Write( h.first.c_str(), TObject::kWriteDelete );
    }
    file.cd();
  }
  for( auto& hMapIt : hMapMap2 ) {
    if( !file.cd(hMapIt.first.c_str()) ) {
      file.mkdir(hMapIt.first.c_str());
      file.cd(hMapIt.first.c_str());
    }
    for( auto& h : hMapIt.second ) {
      h.second.Write( h.first.c_str(), TObject::kWriteDelete );
    }
    file.cd();
  }

  // trigger studies
  file.mkdir("triggerStudies");
  file.cd("triggerStudies");
  for( auto& effIt : effMap ) effIt.second.Write( effIt.first.c_str(), TObject::kWriteDelete );
  file.GetDirectory("triggerStudies")->WriteObject( &rawEff_vs_run, "rawEff_vs_run" );
  file.cd();

  fReader.GetTree()->GetCurrentFile()->Get("TreeWriter/hCutFlow")->Write("hCutFlow");

  file.Close();
  cout << "Created " << outputName << " in " << (time(NULL) - startTime)/60 << " min" << endl;
}

void HistogramProducer::resetSelection() {
  selPhotons.clear();
  selJets.clear();
  selBJets.clear();
  selHEJets.clear();
  selElectrons.clear();
  selMuons.clear();
  artificialPhotons.clear();
}
