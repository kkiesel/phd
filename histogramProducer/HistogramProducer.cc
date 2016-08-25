#include "regex"
#include "time.h"

#include "TROOT.h"
#include "TFile.h"
#include "TSelector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TEfficiency.h"
#include "TRandom2.h"
#include "TTree.h"

#include "TreeParticles.hpp"
#include "UserFunctions.h"
#include "Weighter.h"
#include "CutFlow.h"


class HistogramProducer : public TSelector {
 public:

  HistogramProducer();
  virtual ~HistogramProducer() { }

  virtual void Init(TTree *tree);
  virtual void SlaveBegin(TTree *tree);
  virtual Bool_t Process(Long64_t entry);
  virtual void Terminate();
  virtual Int_t Version() const { return 2; }

  void resetSelection();
  void defaultSelection();
  void fillSelection(string const& s);
  void fillSelectedPhotons(const tree::Particle& p);
  tree::Jet* matchedJet(const tree::Particle& p);

  void initTriggerStudies();
  void fillTriggerStudies();

  void initUncut();
  void fillUncut();

  int genMatch(const tree::Particle& p);

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
  TTreeReaderValue<Bool_t> hlt_ht800;
  TTreeReaderValue<Int_t> hlt_ht600_pre;
//  TTreeReaderValue<std::string> modelName;

  vector<tree::Photon*> selPhotons;
  vector<tree::Jet*> selJets;
  vector<tree::Jet*> selBJets;
  vector<tree::Jet*> selHEJets;
  vector<tree::Electron*> selElectrons;
  vector<tree::Muon*> selMuons;

  vector<tree::Photon> artificialPhotons;

  float selW = 1.; // weight
  float sampleW = 1.;

  map<string,map<string,TH1F>> h1Maps;
  map<string,map<string,TH2F>> h2Maps;
  map<string,map<string,TH3F>> h3Maps;
  map<string,TEfficiency> effMap;
  map<string,TTree*> treeMap;
  map<int,pair<int,int>> rawEff_vs_run;

  bool isData;
  bool zToMet;
  int gluinoMass;

  JetSelector jetSelector;
  Weighter nJetWeighter;
  Weighter emhtWeighter;
  CutFlowPhoton looseCutFlowPhoton;

  double startTime;
  TRandom2 rand;
};

int HistogramProducer::genMatch(const tree::Particle& p) {

  // match to daghters of massive particles
  for (auto const& genP : *intermediateGenParticles) {
    for (auto const & d : genP.daughters) {
      auto id = fabs(d.pdgId);
      auto dr = p.p.DeltaR(d.p);
      auto dpt = p.p.Pt()/d.p.Pt();
      if (dr < 0.5) return id;
    }
  }

  for (auto const& genP : *genParticles) {
    auto id = fabs(genP.pdgId);
    auto dr = p.p.DeltaR(genP.p);
    auto dpt = p.p.Pt()/genP.p.Pt();
    if (dr < 0.15 && fabs(dpt-1) < 0.15) {
      if (genP.isPrompt) return id;
      else return -id;
    }
  }

  return 0;
}

void HistogramProducer::initTriggerStudies() {
  effMap["eff_hlt_pt"] = TEfficiency("", ";#it{p}_{T} (GeV);#varepsilon", 250, 0, 1000);
  effMap["eff_hlt_eta"] = TEfficiency("", ";|#eta|;#varepsilon", 15, 0, 1.5);
  effMap["eff_hlt_ht"] = TEfficiency("", ";#it{EMH}_{T} (GeV);#varepsilon", 200, 0, 2000);
  effMap["eff_hlt_ht_ptMin"] = TEfficiency("", ";#it{EMH}_{T} (GeV);#varepsilon", 200, 0, 2000);
  effMap["eff_hlt_ht_etaMax"] = TEfficiency("", ";#it{EMH}_{T} (GeV);#varepsilon", 310, 0, 3.1);
  effMap["eff_hlt_ht_ct"] = TEfficiency("", ";#it{EMH}_{T} (GeV);#varepsilon", 200, 0, 2000);
  effMap["eff_hlt_ht_ct_preScaled"] = TEfficiency("", ";#it{EMH}_{T} (GeV);#varepsilon (prescaled)", 200, 0, 2000);
  effMap["eff_hlt_ht_ct2"] = TEfficiency("", ";#it{EMH}_{T} (GeV);#varepsilon", 200, 0, 2000);
  effMap["eff_hlt_ht_ct2_preScaled"] = TEfficiency("", ";#it{EMH}_{T} (GeV);#varepsilon (prescaled)", 200, 0, 2000);
  effMap["eff_hlt_nVertex"] = TEfficiency("", ";vertex multiplicity", 41, -0.5, 40.5);
  effMap["eff_hlt_sie"] = TEfficiency("", ";#sigma_{i#etai#eta}", 400, 0, 0.02);
  effMap["eff_hlt_hoe"] = TEfficiency("", ";H/E", 100, 0, 0.15);
  effMap["eff_hlt_r9"] = TEfficiency("", ";r9", 110, 0, 1.1);
  effMap["eff_hlt_cIso"] = TEfficiency("", ";I_{#pi} (GeV)", 100, 0, 10);
  effMap["eff_hlt_nIso"] = TEfficiency("", ";I_{n} (GeV)", 100, 0, 20);
  effMap["eff_hlt_pIso"] = TEfficiency("", ";I_{#gamma} (GeV)", 100, 0, 20);
  effMap["eff_hlt_nJet"] = TEfficiency("", ";uncleaned jet multiplicity", 15, -0.5, 14.5);
  effMap["eff_hlt_met"] = TEfficiency("", ";#it{E}_{T}^{miss} (GeV)", 150, 0, 150);
  effMap["eff_hlt_met_ct"] = TEfficiency("", ";#it{E}_{T}^{miss} (GeV)", 150, 0, 150);
  effMap["eff_hlt_met_ct2"] = TEfficiency("", ";#it{E}_{T}^{miss} (GeV)", 150, 0, 150);
}

void HistogramProducer::fillTriggerStudies() {

  tree::Photon* thisPhoton = 0;
  for (auto& photon : *photons) {
    if (photon.p.Pt() > 15 && fabs(photon.p.Eta()) < photonsEtaMaxBarrel && photon.isLoose && !photon.hasPixelSeed) {
      thisPhoton = &photon;
      break; // take leading(first) photon
    }
  }

  float ht = 0;
  if (thisPhoton) ht += thisPhoton->p.Pt();
  tree::Jet *minPtJet = 0;
  tree::Jet *maxEtaJet = 0;
  for (auto& jet : *jets) {
    if (jet.p.Pt() > 40 && fabs(jet.p.Eta()) < 3) {
      if (!minPtJet || minPtJet->p.Pt() > jet.p.Pt()) minPtJet = &jet;
      if (!maxEtaJet || fabs(minPtJet->p.Eta()) < fabs(jet.p.Eta())) maxEtaJet = &jet;
      if (!thisPhoton || jet.p.DeltaR(thisPhoton->p) > 0.3) {
        ht += jet.p.Pt();
      }
    }
  }

  if (thisPhoton) {

    // get trigger efficiencies for each run
    if (thisPhoton->p.Pt() > 100  && ht > 700 && *hlt_ht600) {
      if (!rawEff_vs_run.count(*runNo)) rawEff_vs_run[*runNo] = make_pair(0,0);
      if (*hlt_photon90_ht500) rawEff_vs_run.at(*runNo).first += 1;
      rawEff_vs_run.at(*runNo).second += 1;
    }

    // main trigger plots
    if (thisPhoton->p.Pt() > 100 && *hlt_photon90) {
      effMap.at("eff_hlt_ht").Fill(*hlt_photon90_ht500, ht);
      effMap.at("eff_hlt_ht_ptMin").Fill(*hlt_photon90_ht500, minPtJet->p.Pt());
      if (ht>600) effMap.at("eff_hlt_ht_etaMax").Fill(*hlt_photon90_ht500, fabs(maxEtaJet->p.Eta()));
    }
    if (thisPhoton->p.Pt() > 100 && *hlt_photon90_ht500) {
      effMap.at("eff_hlt_ht_ct").Fill(*hlt_ht600, ht);
      effMap.at("eff_hlt_ht_ct_preScaled").FillWeighted(*hlt_ht600, *hlt_ht600_pre, ht);
      if (ht > 700) effMap.at("eff_hlt_met_ct").Fill(*hlt_ht600, met->p.Pt());
    }
    if (thisPhoton->p.Pt() > 100 && *hlt_photon90) {
      effMap.at("eff_hlt_ht_ct2").Fill(*hlt_ht600, ht);
      effMap.at("eff_hlt_ht_ct2_preScaled").FillWeighted(*hlt_ht600, *hlt_ht600_pre, ht);
      if (ht > 700) effMap.at("eff_hlt_met_ct2").Fill(*hlt_ht600, met->p.Pt());
    }
    if (ht > 700 && *hlt_ht600) {
      effMap.at("eff_hlt_pt").Fill(*hlt_photon90_ht500, thisPhoton->p.Pt());
      effMap.at("eff_hlt_eta").Fill(*hlt_photon90_ht500, fabs(thisPhoton->p.Eta()));
    }
  }

  if (ht > 700 && *hlt_ht600) {
    for (auto& photon : *photons) {
      if (photon.p.Pt() > 100 && fabs(photon.p.Eta()) < photonsEtaMaxBarrel) {
        if (looseCutFlowPhoton.check(photon)) {
          effMap.at("eff_hlt_r9").Fill(*hlt_photon90_ht500, photon.r9);
          effMap.at("eff_hlt_nVertex").Fill(*hlt_photon90_ht500, *nGoodVertices);
          effMap.at("eff_hlt_nJet").Fill(*hlt_photon90_ht500, jets->size());
          effMap.at("eff_hlt_met").Fill(*hlt_photon90_ht500, met->p.Pt());
        }
        if (
          looseCutFlowPhoton.passHoe()
          && looseCutFlowPhoton.passSie()
          && looseCutFlowPhoton.passNIso()
          && looseCutFlowPhoton.passPIso()
         ) effMap.at("eff_hlt_cIso").Fill(*hlt_photon90_ht500, photon.cIso);
        if (
          looseCutFlowPhoton.passHoe()
          && looseCutFlowPhoton.passSie()
          && looseCutFlowPhoton.passCIso()
          && looseCutFlowPhoton.passPIso()
         ) effMap.at("eff_hlt_nIso").Fill(*hlt_photon90_ht500, photon.nIso);
        if (
          looseCutFlowPhoton.passHoe()
          && looseCutFlowPhoton.passSie()
          && looseCutFlowPhoton.passCIso()
          && looseCutFlowPhoton.passNIso()
         ) effMap.at("eff_hlt_pIso").Fill(*hlt_photon90_ht500, photon.pIso);
        if (
          looseCutFlowPhoton.passHoe()
          && looseCutFlowPhoton.passIso()
       ) effMap.at("eff_hlt_sie").Fill(*hlt_photon90_ht500, photon.sigmaIetaIeta);
        if (
          looseCutFlowPhoton.passSie()
          && looseCutFlowPhoton.passIso()
       ) effMap.at("eff_hlt_hoe").Fill(*hlt_photon90_ht500, photon.hOverE);
      }
    } // photon loop
  } // ht cross trigger
}

void HistogramProducer::initUncut() {
  map<string,TH1F> h;
  h1Maps["uncut"] = h;

  map<string,TH2F> h2;
  h2["dr_vs_relpt"] = TH2F("", ";#DeltaR;#it{p}_{T}^{jet}/#it{p}_{T}^{#gamma}", 100, 0, 0.5, 300, 0, 3);
  h2["dr_vs_relpt_genG"] = TH2F("", ";#DeltaR;#it{p}_{T}^{jet}/#it{p}_{T}^{#gamma}", 100, 0, 0.5, 300, 0, 3);
  h2Maps["uncut"] = h2;
}

void HistogramProducer::fillUncut() {
  h1Maps["uncut"]["genHt"].Fill(*genHt);

  for (auto& g : *photons) {
    if (g.isLoose && !g.hasPixelSeed && g.p.Pt() > 100 && fabs(g.p.Eta()) < photonsEtaMaxBarrel) {

      // search gen match
      bool genMatch = false;
      for (auto& gen : *genParticles) {
        if (gen.pdgId == 22 && gen.isPrompt && gen.fromHardProcess && gen.p.DeltaR(g.p) < 0.3) {
          genMatch = true;
          break;
        }
      }

      for (auto& j : *jets) {
        auto dr = g.p.DeltaR(j.p);
        auto relPt = j.p.Pt()/g.p.Pt();
        h2Maps["uncut"]["dr_vs_relpt"].Fill(dr, relPt);
        if (genMatch) {
          h2Maps["uncut"]["dr_vs_relpt_genG"].Fill(dr, relPt);
        }
      }
    }
  }
}

map<string,TH2F> initHistograms2() {
  map<string,TH2F> hMap;

  hMap["n_jets_vs_photonPosition"] = TH2F("",";jet multiplicity;#gamma position",10, -0.5, 9.5, 10, -0.5, 9.5);
  hMap["g_eta_vs_g_phi"] = TH2F("",";|#eta|;|#phi|", 260, 2.6, 2.6, 100, -3.1, 3.1);
  hMap["met_vs_emht"] = TH2F("", ";#it{E}_{T}^{miss} (GeV);#it{EMH}_{T} (GeV)", 300, 0, 3000, 450, 500, 5000);
  hMap["metPar_vs_emht"] = TH2F("", ";#it{E}_{T}^{miss} #parallel (GeV);#it{EMH}_{T} (GeV)", 600, -3000, 3000, 450, 500, 5000);
  hMap["metPer_vs_emht"] = TH2F("", ";#it{E}_{T}^{miss #perp  } (GeV);#it{EMH}_{T} (GeV)", 300, 0, 3000, 450, 500, 5000);
  hMap["metRaw_vs_emht"] = TH2F("", ";uncorrected #it{E}_{T}^{miss} (GeV);#it{EMH}_{T} (GeV)", 300, 0, 3000, 450, 500, 5000);
  hMap["metParRaw_vs_emht"] = TH2F("", ";uncorrected #it{E}_{T}^{miss} #parallel (GeV);#it{EMH}_{T} (GeV)", 600, -3000, 3000, 450, 500, 5000);
  hMap["metPerRaw_vs_emht"] = TH2F("", ";uncorrected #it{E}_{T}^{miss #perp  } (GeV);#it{EMH}_{T} (GeV)", 300, 0, 3000, 450, 500, 5000);
  hMap["met_vs_n_jet"] = TH2F("", ";#it{E}_{T}^{miss} (GeV);N_{jet}", 300, 0, 3000, 15, -.5, 14.5);
  hMap["met_vs_n_obj"] = TH2F("", ";#it{E}_{T}^{miss} (GeV);N_{jet}", 300, 0, 3000, 15, -.5, 14.5);
  hMap["metRaw_vs_n_jet"] = TH2F("", ";uncorrected #it{E}_{T}^{miss} (GeV);N_{jet}", 300, 0, 3000, 15, -.5, 14.5);
  hMap["metRaw_vs_n_obj"] = TH2F("", ";uncorrected #it{E}_{T}^{miss} (GeV);N_{jet}", 300, 0, 3000, 15, -.5, 14.5);
  hMap["met_vs_g_pt"] = TH2F("", ";#it{E}_{T}^{miss} (GeV);#it{p}_{T} (GeV)", 300, 0, 3000, 100, 0, 1000);
  hMap["metRaw_vs_g_pt"] = TH2F("", ";#it{E}_{T}^{miss} (GeV);#it{p}_{T} (GeV)", 300, 0, 3000, 100, 0, 1000);
  hMap["met_vs_g_e"] = TH2F("", ";#it{E}_{T}^{miss} (GeV);#it{E} (GeV)", 300, 0, 3000, 100, 0, 1000);
  hMap["met_vs_j1_pt"] = TH2F("", ";#it{E}_{T}^{miss} (GeV);#it{p}_{T}^{jet} (GeV)", 300, 0, 3000, 200, 0, 2000);
  hMap["met_vs_mht"] = TH2F("", ";#it{E}_{T}^{miss} (GeV);|#vec{H}_{T}| (GeV)", 300, 0, 3000, 1500, 0, 1500);
  hMap["memht_vs_emht"] = TH2F("", ";#it{EMH}_{T}^{miss} (GeV);#it{EMH}_{T} (GeV)", 300, 0, 3000, 450, 500, 5000);
  hMap["mht_vs_emht"] = TH2F("", ";#it{H}_{T}^{miss} (GeV);#it{EMH}_{T} (GeV)", 300, 0, 3000, 450, 500, 5000);
  hMap["razor"] = TH2F("",";M_{R} (GeV);R^{2}", 500, 0, 5000, 200, 0, 2);
  hMap["st_vs_n_jet"] = TH2F("",";S_{T} (GeV);jet multiplicity", 300, 0, 3000, 10, -.5, 9.5);

  return hMap;
}

map<string,TH1F> initHistograms() {
  map<string,TH1F> hMap;

  hMap["met"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 200, 0, 2000);
  hMap["metUp"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 200, 0, 2000);
  hMap["metDn"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 200, 0, 2000);
  hMap["metPar"] = TH1F("", ";#it{E}_{T}^{miss} #parallel (GeV)", 400, -2000, 2000);
  hMap["metParUp"] = TH1F("", ";#it{E}_{T}^{miss} #parallel (GeV)", 400, -2000, 2000);
  hMap["metParDn"] = TH1F("", ";#it{E}_{T}^{miss} #parallel (GeV)", 400, -2000, 2000);
  hMap["metPer"] = TH1F("", ";#it{E}_{T}^{miss #perp  } (GeV)", 200, 0, 2000);
  hMap["metRaw"] = TH1F("", ";uncorrected #it{E}_{T}^{miss} (GeV)", 200, 0, 2000);
  hMap["metParRaw"] = TH1F("", ";uncorrected #it{E}_{T}^{miss} #parallel (GeV)", 400, -2000, 2000);
  hMap["metPerRaw"] = TH1F("", ";uncorrected #it{E}_{T}^{miss #perp  } (GeV)", 200, 0, 2000);
  hMap["mt_g_met"] = TH1F("", ";m(#gamma,#it{E}_{T}^{miss}_{T}) (GeV)", 150, 0, 1500);

  hMap["metSig"] = TH1F("", ";#it{S}", 3000, 0, 3000);
  hMap["tremht"] = TH1F("", ";#it{EMH}_{T}^{trigger-like} (GeV)", 300, 0, 3000);
  hMap["emht"] = TH1F("", ";#it{EMH}_{T} (GeV)", 300, 0, 3000);
  hMap["emhtNoLep"] = TH1F("", ";lepton cleaned #it{EMH}_{T} (GeV)", 300, 0, 3000);
  hMap["emhtStar"] = TH1F("", ";#it{EMH}_{T}* (GeV)", 300, 0, 3000);
  hMap["ht"] = TH1F("", ";#it{H}_{T} (GeV)", 300, 0, 3000);
  hMap["st"] = TH1F("", ";#it{S}_{T} (GeV)", 200, 0, 4000);
  hMap["emrecoilt"] = TH1F("", ";#vec{#it{EMH}}_{T} (GeV)", 150, 0, 1500);
  hMap["recoilt"] = TH1F("", ";#vec{#it{H}}_{T} (GeV)", 150, 0, 1500);

  // photon
  hMap["g_pt"] = TH1F("", ";#it{p}_{T} (GeV)", 150, 0, 1500);
  hMap["g_e"] = TH1F("", ";#it{E} (GeV)", 150, 0, 1500);
  hMap["g_ptStar"] = TH1F("", ";#it{p}_{T}* (GeV)", 150, 0, 1500);
  hMap["g_eta"] = TH1F("", ";|#eta|", 2600, 0, 2.6);

  // he-jet
  hMap["j1_pt"] = TH1F("", ";#it{p}_{T}^{1.jet} (GeV)", 150, 0, 1500);
  hMap["j1_eta"] = TH1F("", ";|#eta^{1.jet}|", 150, 0, 3);
  hMap["j2_pt"] = TH1F("", ";#it{p}_{T}^{2.jet} (GeV)", 150, 0, 1500);
  hMap["j2_eta"] = TH1F("", ";|#eta^{2.jet}|", 150, 0, 3);
  hMap["j3_pt"] = TH1F("", ";#it{p}_{T}^{3.jet} (GeV)", 150, 0, 1500);
  hMap["j3_eta"] = TH1F("", ";|#eta^{3.jet}|", 150, 0, 3);

  // angles
  hMap["dphi_met_g"] = TH1F("", ";|#Delta#phi(#it{E}_{T}^{miss},#gamma)|", 70, 0, 3.5);
  hMap["dphi_met_j1"] = TH1F("", ";|#Delta#phi(#it{E}_{T}^{miss},1.jet)|", 70, 0, 3.5);
  hMap["dphi_met_j2"] = TH1F("", ";|#Delta#phi(#it{E}_{T}^{miss},2.jet)|", 70, 0, 3.5);
  hMap["dphi_met_j3"] = TH1F("", ";|#Delta#phi(#it{E}_{T}^{miss},3.jet)|", 70, 0, 3.5);
  hMap["dphi_met_recoil"] = TH1F("", ";|#Delta#phi(#it{E}_{T}^{miss},#vec{#it{H}}_{T})|", 70, 0, 3.5);
  hMap["dphi_g_j1"] = TH1F("", ";|#Delta#phi(#gamma,1.jet)|", 70, 0, 3.5);
  hMap["dphi_g_j2"] = TH1F("", ";|#Delta#phi(#gamma,2.jet)|", 70, 0, 3.5);

  // multiplicities
  hMap["n_vertex"] = TH1F("", ";vertex multiplicity", 61, -0.5, 60.5);
  hMap["n_vertex_unw"] = TH1F("", ";vertex multiplicity", 61, -0.5, 60.5);
  hMap["n_photon"] = TH1F("", ";photon multiplicity", 6, -0.5, 5.5);
  hMap["n_jet"] = TH1F("", ";jet multiplicity", 16, -0.5, 15.5);
  hMap["n_bjet"] = TH1F("", ";b-jet multiplicity", 16, -0.5, 15.5);
  hMap["n_electron"] = TH1F("", ";electron multiplicity", 4, -0.5, 3.5);
  hMap["n_muon"] = TH1F("", ";muon multiplicity", 4, -0.5, 3.5);
  hMap["n_heJet"] = TH1F("", ";photon-like jet multiplicity", 11, -0.5, 10.5);

  hMap["genMatch"] = TH1F("", ";pdg id for gen match", 24, -0.5, 23.5);
  hMap["genHt"] = TH1F("", ";#it{H}_{T}^{gen}", 3000, 0, 3000);

  return hMap;
}

void HistogramProducer::fillSelection(string const& s) {

  float tree_m, tree_mRaw, tree_w, tree_emht, tree_pt;
  UInt_t tree_njet, tree_genMatch;
  if (!h1Maps.count(s)) {
    h1Maps[s] = initHistograms();
    h2Maps[s] = initHistograms2();
    treeMap[s] = new TTree("simpleTree", "");
    treeMap[s]->Branch("met", &tree_m);
    treeMap[s]->Branch("metRaw", &tree_mRaw);
    treeMap[s]->Branch("emht", &tree_emht);
    treeMap[s]->Branch("weight", &tree_w);
    treeMap[s]->Branch("pt", &tree_pt, "pt");
    treeMap[s]->Branch("njet", &tree_njet, "njet/i");
    treeMap[s]->Branch("genMatch", &tree_genMatch, "genMatch/i");
  }
  auto m1 = &h1Maps[s];
  auto m2 = &h2Maps[s];

  tree_m = met->p.Pt();
  tree_mRaw = met->p.Pt();
  tree_w = selW * sampleW;
  tree_genMatch = true;
  tree_pt = 0;
  if (!isData and jets->size()) {
    tree_genMatch = false;
    auto& p = jets->at(0).p;
    for (auto& gj : *genJets) {
      if (gj.p.DeltaR(p)<0.5 and p.Pt()/gj.p.Pt()<3) {
        tree_genMatch = true;
        break;
      }
    }
  }
  if (selPhotons.size()) {
    tree_pt = selPhotons.at(0)->p.Pt();
  }

  // calculate variables
  float ht = 0;
  float emhtNoLep = 0;
  TVector3 recoil(0,0,0);
  for (auto& j : selJets) {
    ht += j->p.Pt();
    recoil += j->p;
    if (!j->hasPhotonMatch && !j->hasElectronMatch && !j->hasMuonMatch) {
      emhtNoLep += j->p.Pt();
    }
  }

  TVector3 emrecoil = recoil;
  float emht = ht;
  float emhtStar = ht;
  for (auto& p : selPhotons) {
    emrecoil += p->p;
    emht += p->p.Pt();
    emhtNoLep += p->p.Pt();
    auto mJet = matchedJet(*p);
    if (mJet) emhtStar += mJet->p.Pt();
    else emhtStar = p->p.Pt();
  }
  tree_emht = emht;
  tree_njet = selJets.size();
  treeMap[s]->Fill();
  float st = emht + met->p.Pt();

  float tremht = 0;
  for (auto& jet : *jets) {
    auto pt = jet.p.Pt();
    if ( pt > 40 && fabs(jet.p.Eta()) < 3) {
      tremht += pt;
    }
  }

  m1->at("metSig").Fill(met->sig, selW);
  m1->at("emhtNoLep").Fill(emhtNoLep, selW);
  m1->at("tremht").Fill(tremht, selW);
  m1->at("emht").Fill(emht, selW);
  m1->at("emhtStar").Fill(emhtStar, selW);
  m1->at("ht").Fill(ht, selW);
  m1->at("st").Fill(st, selW);
  m1->at("recoilt").Fill(recoil.Pt(), selW);
  m1->at("emrecoilt").Fill(emrecoil.Pt(), selW);

  m1->at("met").Fill(met->p.Pt(), selW);
  m1->at("metRaw").Fill(met->p.Pt(), selW);

  if (selPhotons.size() > 0) {
    auto g = selPhotons.at(0);
    m1->at("genMatch").Fill(genMatch(*g), selW);
    m1->at("metUp").Fill((met->p+g->sigmaPt/g->p.Pt()*g->p).Pt(), selW);
    m1->at("metDn").Fill((met->p-g->sigmaPt/g->p.Pt()*g->p).Pt(), selW);
    auto mJet = matchedJet(*g);
    if (mJet) {
      m1->at("g_ptStar").Fill(mJet->p.Pt(), selW);
    }

    m1->at("mt_g_met").Fill((g->p + met->p).Pt(), selW);
    m1->at("g_pt").Fill(g->p.Pt(), selW);
    m1->at("g_e").Fill(g->p.Mag(), selW);
    m1->at("g_eta").Fill(fabs(g->p.Eta()), selW);
    float dphi_met_g = fabs(met->p.DeltaPhi(g->p));
    m1->at("dphi_met_g").Fill(dphi_met_g, selW);
    m1->at("metPar").Fill(met->p.Pt()*cos(dphi_met_g), selW);
    m1->at("metParUp").Fill(met->p.Pt()*cos(dphi_met_g)+g->sigmaPt, selW);
    m1->at("metParDn").Fill(met->p.Pt()*cos(dphi_met_g)-g->sigmaPt, selW);
    m1->at("metPer").Fill(fabs(met->p.Pt()*sin(dphi_met_g)), selW);
    m1->at("metParRaw").Fill(met->p.Pt()*cos(met->p.DeltaPhi(g->p)), selW);
    m1->at("metPerRaw").Fill(fabs(met->p.Pt()*sin(met->p.DeltaPhi(g->p))), selW);
    m2->at("metPar_vs_emht").Fill(met->p.Pt()*cos(dphi_met_g), emht, selW);
    m2->at("metPer_vs_emht").Fill(met->p.Pt()*sin(dphi_met_g), emht, selW);
    m2->at("metParRaw_vs_emht").Fill(met->p.Pt()*cos(dphi_met_g), emht, selW);
    m2->at("metPerRaw_vs_emht").Fill(met->p.Pt()*sin(dphi_met_g), emht, selW);
    unsigned photonPosition=0;
    for (;photonPosition<selJets.size() && selJets.at(photonPosition)->p.Pt() > g->p.Pt();photonPosition++);
    m2->at("n_jets_vs_photonPosition").Fill(selJets.size(), photonPosition, selW);
    m2->at("g_eta_vs_g_phi").Fill(g->p.Eta(), g->p.Phi(), selW);
  }

  if (selJets.size() > 2) {
    m1->at("j3_pt").Fill(selJets.at(2)->p.Pt(), selW);
    m1->at("j3_eta").Fill(fabs(selJets.at(2)->p.Eta()), selW);
    m1->at("dphi_met_j3").Fill(fabs(met->p.DeltaPhi(selJets.at(2)->p)), selW);
  }
  if (selJets.size() > 1) {
    m1->at("j2_pt").Fill(selJets.at(1)->p.Pt(), selW);
    m1->at("j2_eta").Fill(fabs(selJets.at(1)->p.Eta()), selW);
    m1->at("dphi_met_j2").Fill(fabs(met->p.DeltaPhi(selJets.at(1)->p)), selW);
    if (selPhotons.size() > 0) m1->at("dphi_g_j2").Fill(fabs(selPhotons.at(0)->p.DeltaPhi(selJets.at(1)->p)), selW);
  }
  if (selJets.size() > 0) {
    m1->at("j1_pt").Fill(selJets.at(0)->p.Pt(), selW);
    m1->at("j1_eta").Fill(fabs(selJets.at(0)->p.Eta()), selW);
    m1->at("dphi_met_j1").Fill(fabs(met->p.DeltaPhi(selJets.at(0)->p)), selW);
    if (selPhotons.size() > 0) m1->at("dphi_g_j1").Fill(fabs(selPhotons.at(0)->p.DeltaPhi(selJets.at(0)->p)), selW);
  }

  m1->at("dphi_met_recoil").Fill(fabs(met->p.DeltaPhi(recoil)), selW);

  m1->at("n_vertex").Fill(*nGoodVertices, selW);
  m1->at("n_vertex_unw").Fill(*nGoodVertices);
  m1->at("n_photon").Fill(selPhotons.size(), selW);
  m1->at("n_jet").Fill(selJets.size(), selW);
  m1->at("n_bjet").Fill(selBJets.size(), selW);
  m1->at("n_electron").Fill(selElectrons.size(), selW);
  m1->at("n_muon").Fill(selMuons.size(), selW);
  m1->at("n_heJet").Fill(selJets.size(), selW);
  m1->at("genHt").Fill(*genHt, selW);

  m2->at("met_vs_n_jet").Fill(met->p.Pt(), selJets.size(), selW);
  m2->at("met_vs_n_obj").Fill(met->p.Pt(), selJets.size()+selPhotons.size(), selW);
  m2->at("metRaw_vs_n_jet").Fill(met->p.Pt(), selJets.size(), selW);
  m2->at("metRaw_vs_n_obj").Fill(met->p.Pt(), selJets.size()+selPhotons.size(), selW);
  if (selJets.size()) m2->at("met_vs_j1_pt").Fill(met->p.Pt(), selJets.at(0)->p.Pt() , selW);
  if (selPhotons.size()) {
    m2->at("met_vs_g_pt").Fill(met->p.Pt(), selPhotons.at(0)->p.Pt() , selW);
    m2->at("metRaw_vs_g_pt").Fill(met->p.Pt(), selPhotons.at(0)->p.Pt() , selW);
    m2->at("met_vs_g_e").Fill(met->p.Pt(), selPhotons.at(0)->p.Pt() , selW);
  } else {
    m2->at("metRaw_vs_g_pt").Fill(met->p.Pt(), 0., selW);
  }
  m2->at("met_vs_emht").Fill(met->p.Pt(), emht, selW);
  m2->at("metRaw_vs_emht").Fill(met->p.Pt(), emht, selW);
  m2->at("met_vs_mht").Fill(met->p.Pt(), recoil.Pt(), selW);
  m2->at("memht_vs_emht").Fill(emrecoil.Pt(), emht, selW);
  m2->at("mht_vs_emht").Fill(recoil.Pt(), emht, selW);

  // razor
  std::vector<TVector3> razorObjects;
  for (auto j: *jets) {
    razorObjects.push_back(j.p);
  }

  auto megaJets = megajets(razorObjects);
  auto razorVar = razorVariables(megaJets, met->p);
  m2->at("razor").Fill(razorVar.first, razorVar.second, selW);

  // st stuff
  m2->at("st_vs_n_jet").Fill(st, selJets.size(), selW);


} // end fill


HistogramProducer::HistogramProducer():
  photons(fReader, "photons"),
  jets(fReader, "jets"),
  electrons(fReader, "electrons"),
  muons(fReader, "muons"),
  genJets(fReader, "genJets"),
  genParticles(fReader, "genParticles"),
  intermediateGenParticles(fReader, "intermediateGenParticles"),
  met(fReader, "met"),
  nGoodVertices(fReader, "nGoodVertices"),
  pu_weight(fReader, "pu_weight"),
  mc_weight(fReader, "mc_weight"),
  genHt(fReader, "genHt"),
  runNo(fReader, "runNo"),
  hlt_photon90_ht500(fReader, "HLT_Photon90_CaloIdL_PFHT500_v"),
  hlt_photon90(fReader, "HLT_Photon90_v"),
  hlt_ht600(fReader, "HLT_PFHT600_v"),
  hlt_ht800(fReader, "HLT_PFHT800_v"),
  hlt_ht600_pre(fReader, "HLT_PFHT600_v_pre"),
//  modelName(fReader, "modelName"),
  jetSelector("../plotter/gammaPosition.root", "gqcd"),
  nJetWeighter("../plotter/weights.root", "weight_n_heJet"),
  emhtWeighter("../plotter/weights.root", "weight_emht_gqcd"),
  looseCutFlowPhoton({{"sigmaIetaIeta_eb",0.0102}, {"cIso_eb",3.32}, {"nIso1_eb",1.92}, {"nIso2_eb",0.014}, {"nIso3_eb",0.000019}, {"pIso1_eb",0.81}, {"pIso2_eb",0.0053},
    {"sigmaIetaIeta_ee",0.0274}, {"cIso_ee",1.97}, {"nIso1_ee",11.86}, {"nIso2_ee",0.0139}, {"nIso3_ee",0.000025}, {"pIso1_ee",0.83}, {"pIso2_ee",0.0034} }),
  gluinoMass(0),
  startTime(time(NULL)),
  rand()
{
}

void HistogramProducer::Init(TTree *tree)
{
  fReader.SetTree(tree);
  string inputName = fReader.GetTree()->GetCurrentFile()->GetName();
  isData = inputName.find("Run201") != string::npos;
  zToMet = inputName.find("ZGTo2LGmod") != string::npos;

  float lumi = 2.32e3; // pb^{-1}
  float nGen = ((TH1F*)fReader.GetTree()->GetCurrentFile()->Get("TreeWriter/hCutFlow"))->GetBinContent(2);
  sampleW = isData ? 1. : lumi * sampleCrossSection(inputName) / nGen;


  std::smatch sm;
  if (regex_match(inputName, sm, regex(".*/T5.*_(\\d+)_(\\d+).root"))) {
      gluinoMass = stoi(sm[1]);
  }
  initTriggerStudies();
  initUncut();
}

void HistogramProducer::SlaveBegin(TTree *tree)
{
}

void HistogramProducer::defaultSelection()
{
  for (auto& mu : *muons) {
    if (mu.p.Pt() < 15) continue;
    if (indexOfMatchedParticle<tree::Photon*>(mu, selPhotons, .3) >= 0) continue;
    selMuons.push_back(&mu);
  }
  for (auto& el : *electrons) {
    if (!el.isLoose || el.p.Pt() < 15) continue;
    if (indexOfMatchedParticle<tree::Photon*>(el, selPhotons, .3) >= 0) continue;
    selElectrons.push_back(&el);
  }
  for (auto& jet : *jets) {
    if (!jet.isLoose
//      || jet.hasPhotonMatch || jet.hasElectronMatch || jet.hasMuonMatch
      || indexOfMatchedParticle<tree::Photon*>(jet, selPhotons, .3) >= 0
      || jet.p.Pt() < 40 || fabs(jet.p.Eta()) > 3) continue;

    selJets.push_back(&jet);
    if (jet.bDiscriminator > bTaggingWorkingPoints.at(CSVv2M) && fabs(jet.p.Eta()) < 2.5)
      selBJets.push_back(&jet);
  }
}

void HistogramProducer::fillSelectedPhotons(const tree::Particle& p) {
  tree::Photon photon;
  photon.p = p.p;
  artificialPhotons.push_back(photon);
  selPhotons.push_back(&artificialPhotons.back());
}

tree::Jet* HistogramProducer::matchedJet(const tree::Particle& p) {
  for (auto& jet : *jets) {
    if (jet.p.DeltaR(p.p) < .2) {
      return &jet;
    }
  }
  return 0;
}


Bool_t HistogramProducer::Process(Long64_t entry)
{
  resetSelection();
  fReader.SetEntry(entry);
  fillUncut();
//  cout << *modelName << endl;

  float zToMetPt = -1;
  if (zToMet && intermediateGenParticles->size()) {
    auto z = intermediateGenParticles->at(0);
    if (z.daughters.size() > 2) {
      return kTRUE;
    }
    met->p += z.p;
    for (auto d : z.daughters) {
      for (auto p = photons->begin(); p!=photons->end();) {
        if (p->p.DeltaR(d.p) < 0.1) p = photons->erase(p);
        else p++;
      }
      for (auto p = jets->begin(); p!=jets->end();) {
        if (p->p.DeltaR(d.p) < 0.1) p = jets->erase(p);
        else p++;
      }
      for (auto p = electrons->begin(); p!=electrons->end();) {
        if (p->p.DeltaR(d.p) < 0.1) p = electrons->erase(p);
        else p++;
      }
      for (auto p = muons->begin(); p!=muons->end();) {
        if (p->p.DeltaR(d.p) < 0.1) p = muons->erase(p);
        else p++;
      }
    }
    for (auto p : *genParticles) {
      if (p.pdgId == 22) {
        zToMetPt = p.p.Pt();
        break;
      }
    }
  }
  if (zToMet && zToMetPt<0) return kTRUE;


  if (isData) fillTriggerStudies();

  // set weight
  selW = *mc_weight * *pu_weight;

  // For FastSim CMSSW7X, there are events with large pt jets (pt>sqrt(2))
  // It was recommended for the 2015 analysis to ignore events with objets > 2 times gluino mass
  // TODO: remove these events from the acceptance, change nGen (should be small influence)
  // see here: https://hypernews.cern.ch/HyperNews/CMS/get/met/432.html
  if (gluinoMass > 0 && jets->size() && jets->at(0).p.Pt() > 2*gluinoMass) return kTRUE;


  /////////////////////////////////////////////////////////////////////////////
  // signal sample
  /////////////////////////////////////////////////////////////////////////////

  for (auto& photon : *photons) {
    if (photon.isLoose && !photon.hasPixelSeed && photon.p.Pt() > 100 && fabs(photon.p.Eta()) < photonsEtaMaxBarrel) {
      selPhotons.push_back(&photon);
    }
  }
  defaultSelection();

  float myHt=0;
  for (auto& p : selPhotons) myHt += p->p.Pt();
  for (auto& p : selJets) myHt += p->p.Pt();

  if (selPhotons.size() && myHt > 700 && (*hlt_photon90_ht500 || !isData)) {
    fillSelection("tr");
    if (selPhotons.at(0)->isTrue and selPhotons.at(0)->isTrueAlternative) fillSelection("tr_true");
    if (selPhotons.at(0)->isTight) fillSelection("tr_tight");
    if (met->p.Pt() < 100) fillSelection("tr_0met100");
    else                    fillSelection("tr_100met");
    auto gMatch = genMatch(*selPhotons.at(0));
    if (gMatch == 11) fillSelection("tr_genE");
    else              fillSelection("tr_noGenE");
  }

  if (!selPhotons.size() && myHt > 700 && (*hlt_ht600 || !isData)) {
    fillSelection("tr_jControl");
  }

  resetSelection();
  /////////////////////////////////////////////////////////////////////////////
  // electron sample
  /////////////////////////////////////////////////////////////////////////////

  for (auto& photon : *photons) {
    if (photon.isLoose && photon.hasPixelSeed && photon.p.Pt() > 100 && fabs(photon.p.Eta()) < photonsEtaMaxBarrel) {
      selPhotons.push_back(&photon);
    }
  }
  defaultSelection();

  myHt=0;
  for (auto& p : selPhotons) myHt += p->p.Pt();
  for (auto& p : selJets) myHt += p->p.Pt();
  if (selPhotons.size() && myHt > 700 && (*hlt_photon90_ht500 || !isData)) {
    fillSelection("tr_eControl");
    if (genMatch(*selPhotons.at(0)) == 11) fillSelection("tr_eControl_genE");
  }

  resetSelection();
  /////////////////////////////////////////////////////////////////////////////
  // endcap selection
  /////////////////////////////////////////////////////////////////////////////

  for (auto& photon : *photons) {
    if (photon.isLoose && !photon.hasPixelSeed && photon.p.Pt() > 100 && fabs(photon.p.Eta()) > photonsEtaMaxBarrel) {
      selPhotons.push_back(&photon);
    }
  }
  defaultSelection();

  myHt=0;
  for (auto& p : selPhotons) myHt += p->p.Pt();
  for (auto& p : selJets) myHt += p->p.Pt();

  if (selPhotons.size() && myHt > 700 && (*hlt_photon90_ht500 || !isData)) {
    fillSelection("tr_ee");
  }

  resetSelection();
  /////////////////////////////////////////////////////////////////////////////
  // low pt selection
  /////////////////////////////////////////////////////////////////////////////

  for (auto& photon : *photons) {
    if (photon.isLoose && !photon.hasPixelSeed && photon.p.Pt() > 40 && fabs(photon.p.Eta()) < photonsEtaMaxBarrel) {
      selPhotons.push_back(&photon);
    }
  }
  defaultSelection();

  myHt=0;
  for (auto& p : selPhotons) myHt += p->p.Pt();
  for (auto& p : selJets) myHt += p->p.Pt();

  if (selPhotons.size() && myHt > 700 && (*hlt_ht600 || !isData)) {
    fillSelection("tr_lowPt");
  }

  resetSelection();
  /////////////////////////////////////////////////////////////////////////////
  // low pt ee selection
  /////////////////////////////////////////////////////////////////////////////

  for (auto& photon : *photons) {
    if (photon.isLoose && !photon.hasPixelSeed && photon.p.Pt() > 40 && fabs(photon.p.Eta()) > photonsEtaMaxBarrel) {
      selPhotons.push_back(&photon);
    }
  }
  defaultSelection();

  myHt=0;
  for (auto& p : selPhotons) myHt += p->p.Pt();
  for (auto& p : selJets) myHt += p->p.Pt();

  if (selPhotons.size() && myHt > 700 && (*hlt_ht600 || !isData)) {
    fillSelection("tr_lowPt_ee");
  }

  return kTRUE;
}

template<typename T>
void save2File(const map<string,map<string,T>>& hMaps, TFile& file)
{
  for (auto& hMapIt : hMaps) {
    if (!file.Get(hMapIt.first.c_str())) {
      file.mkdir(hMapIt.first.c_str());
    }
    file.cd(hMapIt.first.c_str());
    for (auto& h : hMapIt.second) {
      h.second.Write(h.first.c_str(), TObject::kWriteDelete);
    }
    file.cd();
  }
}

void HistogramProducer::Terminate()
{
  auto outputName = getOutputFilename(fReader.GetTree()->GetCurrentFile()->GetName());
  TFile file(outputName.c_str(), "RECREATE");
  save2File(h1Maps, file);
  save2File(h2Maps, file);
  for (auto& metsMapIt : treeMap) {
    file.cd(metsMapIt.first.c_str());
    metsMapIt.second->Write("simpleTree", TObject::kWriteDelete);
    file.cd();
  }

  file.mkdir("triggerStudies");
  file.cd("triggerStudies");
  for (auto& effIt : effMap) effIt.second.Write(effIt.first.c_str(), TObject::kWriteDelete);
  file.GetDirectory("triggerStudies")->WriteObject(&rawEff_vs_run, "rawEff_vs_run");
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
