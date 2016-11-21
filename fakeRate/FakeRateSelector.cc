#include <math.h>
#include <regex>
#include <time.h>

#include "TH2F.h"
#include "TH3F.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TSelector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TRandom2.h"

#include "TreeParticles.hpp"
#include "UserFunctions.h"

static const float triggerDRcut = .2;
static const float tpSeperationDRcut = .3;

float m(const TVector3& p1, const TVector3& p2) {
  TLorentzVector a(p1, p1.Mag());
  TLorentzVector b(p2, p2.Mag());
  return (a+b).M();
}

static const float zBosonMass = 91.1876;

class FakeRateSelector : public TSelector {
 public:

  FakeRateSelector();
  virtual ~FakeRateSelector() { }

  virtual void Init(TTree *tree);
  virtual void SlaveBegin(TTree *tree);
  virtual Bool_t Process(Long64_t entry);
  virtual void Terminate();
  virtual Int_t Version() const { return 2; }

  void fillSelection(tree::Electron* tag, tree::Photon* probe, const string& name="");
  void fillGenSelection(const tree::Photon& probe, const string& name="");
  TTreeReader fReader;
  TTreeReaderValue<std::vector<tree::Photon>> photons;
  TTreeReaderValue<std::vector<tree::Jet>> jets;
  TTreeReaderValue<std::vector<tree::Electron>> electrons;
  TTreeReaderValue<std::vector<tree::Muon>> muons;
  TTreeReaderValue<std::vector<tree::Particle>> triggerObjects;
  TTreeReaderValue<std::vector<tree::IntermediateGenParticle>> intermediateGenParticles;
  TTreeReaderValue<tree::MET> met;
  TTreeReaderValue<Int_t> nGoodVertices;
  TTreeReaderValue<Int_t> nTracksPV;
  TTreeReaderValue<Bool_t> hlt;

  double startTime;
  map<string,map<string,TEfficiency>> effMaps;
  map<string,TH1F> h1Map;
  bool isData;
};


FakeRateSelector::FakeRateSelector():
  photons(fReader, "photons"),
  jets(fReader, "jets"),
  electrons(fReader, "electrons"),
  muons(fReader, "muons"),
  triggerObjects(fReader, "triggerObjects"),
  intermediateGenParticles(fReader, "intermediateGenParticles"),
  met(fReader, "met"),
  nGoodVertices(fReader, "nGoodVertices"),
  nTracksPV(fReader, "nTracksPV"),
  hlt(fReader, "HLT_Ele27_eta2p1_WPLoose_Gsf_v"),
  startTime(time(NULL))
{
}

void FakeRateSelector::Init(TTree *tree)
{
  fReader.SetTree(tree);
  string inputName = fReader.GetTree()->GetCurrentFile()->GetName();
  isData = inputName.find("Run201") != string::npos;
}

void FakeRateSelector::SlaveBegin(TTree *tree)
{
  h1Map["n_tags"] = TH1F("", ";number of tags;", 4, -.5, 3.5);
  h1Map["n_probes"] = TH1F("", ";number of probes;", 4, -.5, 3.5);
  h1Map["n_muons"] = TH1F("", ";number of muons;", 4, -.5, 3.5);
}

map<string,TEfficiency> initGenEfficiencies() {
  map<string,TEfficiency> m;
  m["pt"] = TEfficiency("",";#it{p}_{T} (GeV);#varepsilon", 40, 0, 200);
  m["met"] = TEfficiency("",";#it{E}_{T}^{miss} (GeV);#varepsilon", 40, 0, 200);
  m["emht"] = TEfficiency("",";#it{EMH}_{T} (GeV);#varepsilon", 40, 0, 400);
  m["eta"] = TEfficiency("",";|#eta|;#varepsilon", 300, 0, 3);
  m["jets"] = TEfficiency("",";jet multiplicity;#varepsilon", 21, -0.5, 20.5);
  m["vtx"] = TEfficiency("",";vertex multiplicity;#varepsilon", 36, 0.5, 36.5);
  m["sie"] = TEfficiency("",";#sigma_{i#etai#eta};#varepsilon", 60, 0, 0.03);
  m["sip"] = TEfficiency("",";#sigma_{i#phii#phi};#varepsilon", 80, 0, 0.08);
  m["hoe"] = TEfficiency("",";H/E;#varepsilon", 100, 0, 0.05);
  m["r9"] = TEfficiency("",";r9;#varepsilon", 60, 0, 1.1);
  m["cIso"] = TEfficiency("",";I_{#pm} (GeV);#varepsilon", 35, 0, 3.5);
  m["nIso"] = TEfficiency("",";I_{0} (GeV);#varepsilon", 50, 0, 5);
  m["pIso"] = TEfficiency("",";I_{#gamma} (GeV);#varepsilon", 20, 0, 2);
  m["cIsoWorst"] = TEfficiency("",";worst I_{#pm} (GeV);#varepsilon", 80, 0, 40);
  m["nTracksPV"] = TEfficiency("",";PV track multiplicity;#varepsilon", 50, -.5, 49.5);
  m["pt_vs_vtx"] = TEfficiency("", ";#it{p}_{T} (GeV);vertex multiplicity", 40, 0, 200, 36, 0.5, 36.5);
  m["pt_vs_nTracksPV"] = TEfficiency("", ";#it{p}_{T} (GeV);PV track multiplicity", 40, 0, 200, 50, -0.5, 49.5);
  m["vtx_vs_nTracksPV"] = TEfficiency("", ";m (GeV);vertex multiplicity ;PV track multiplicity", 200, 40, 140, 36, .5, 36.5, 50, -0.5, 49.5);
  m["met_vs_vtx"] = TEfficiency("", ";#it{E}_{T}^{miss} (GeV);vertex multiplicity", 40, 0, 200, 36, 0.5, 36.5);
  m["met_vs_jets"] = TEfficiency("", ";#it{E}_{T}^{miss} (GeV);jet multiplicity", 40, 0, 200, 8, -0.5, 7.5);
  m["met_vs_pt"] = TEfficiency("", ";#it{E}_{T}^{miss} (GeV);#it{p}_{T} (GeV)", 40, 0, 200, 40, 0, 200);
  m["met_vs_jetPt"] = TEfficiency("", ";#it{E}_{T}^{miss} (GeV);#it{p}_{T}^{1.jet} (GeV)", 40, 0, 200, 40, 0, 200);
  m["eta_vs_phi"] = TEfficiency("", ";#eta;#phi", 250, -2.5, 2.5, 64, -3.2, 3.2);
  return m;
}

void FakeRateSelector::fillGenSelection(const tree::Photon& probe, const string& name)
{
  if (!effMaps.count(name)) {
    effMaps[name] = initGenEfficiencies();
  }
  auto effs = &effMaps.at(name);

  bool hasPixelSeed = probe.hasPixelSeed;
  float thisMet = met->p.Pt();
  unsigned nJet = 0;
  float emht = probe.p.Pt();;
  for (auto& j : *jets ) {
    if (j.p.Pt()>30 && j.p.Eta()<3 && probe.p.DeltaR(j.p)>0.4) {
      nJet++;
      emht += j.p.Pt();
    }
  }
  effs->at("pt").Fill(!hasPixelSeed, probe.p.Pt());
  effs->at("pt").Fill(!hasPixelSeed, probe.p.Pt());
  effs->at("eta").Fill(!hasPixelSeed, fabs(probe.p.Eta()));
  effs->at("vtx").Fill(!hasPixelSeed, *nGoodVertices);
  effs->at("emht").Fill(!hasPixelSeed, emht);
  effs->at("jets").Fill(!hasPixelSeed, nJet);
  effs->at("met").Fill(!hasPixelSeed, thisMet);
  effs->at("sie").Fill(!hasPixelSeed, probe.sigmaIetaIeta);
  effs->at("sip").Fill(!hasPixelSeed, probe.sigmaIphiIphi);
  effs->at("hoe").Fill(!hasPixelSeed, probe.hOverE);
  effs->at("r9").Fill(!hasPixelSeed, probe.r9);
  effs->at("cIso").Fill(!hasPixelSeed, probe.cIso);
  effs->at("nIso").Fill(!hasPixelSeed, probe.nIso);
  effs->at("pIso").Fill(!hasPixelSeed, probe.pIso);
  effs->at("cIsoWorst").Fill(!hasPixelSeed, probe.cIsoWorst);
  effs->at("nTracksPV").Fill(!hasPixelSeed, *nTracksPV);

  effs->at("pt_vs_vtx").Fill(!hasPixelSeed, probe.p.Pt(), *nGoodVertices);
  effs->at("pt_vs_nTracksPV").Fill(!hasPixelSeed, probe.p.Pt(), *nTracksPV);
  effs->at("vtx_vs_nTracksPV").Fill(!hasPixelSeed, *nGoodVertices, *nTracksPV);
  effs->at("met_vs_vtx").Fill(!hasPixelSeed, thisMet, *nGoodVertices);
  effs->at("met_vs_jets").Fill(!hasPixelSeed, thisMet, nJet);
  effs->at("met_vs_pt").Fill(!hasPixelSeed, thisMet, probe.p.Pt());
  //effs->at("met_vs_jetPt").Fill(!hasPixelSeed, thisMet, leadingJetPt);
  effs->at("eta_vs_phi").Fill(!hasPixelSeed, probe.p.Eta(), probe.p.Phi());
}

map<string,TEfficiency> initEfficiencies() {
  map<string,TEfficiency> m;
  m["pt"] = TEfficiency("",";m (GeV);#it{p}_{T} (GeV);#varepsilon", 200, 40, 140, 40, 0, 200);
  m["met"] = TEfficiency("",";m (GeV);#it{E}_{T}^{miss} (GeV);#varepsilon", 200, 40, 140, 40, 0, 200);
  m["emht"] = TEfficiency("",";m (GeV);#it{EMH}_{T} (GeV);#varepsilon", 200, 40, 140, 40, 0, 400);
  m["eta"] = TEfficiency("",";m (GeV);|#eta|;#varepsilon", 200, 40, 140, 300, 0, 3);
  m["jets"] = TEfficiency("",";m (GeV);jet multiplicity;#varepsilon", 200, 40, 140, 21, -0.5, 20.5);
  m["vtx"] = TEfficiency("",";m (GeV);vertex multiplicity;#varepsilon", 200, 40, 140, 36, 0.5, 36.5);
  m["sie"] = TEfficiency("",";m (GeV);#sigma_{i#etai#eta};#varepsilon", 200, 40, 140, 60, 0, 0.03);
  m["sip"] = TEfficiency("",";m (GeV);#sigma_{i#phii#phi};#varepsilon", 200, 40, 140, 80, 0, 0.08);
  m["hoe"] = TEfficiency("",";m (GeV);H/E;#varepsilon", 200, 40, 140, 100, 0, 0.05);
  m["r9"] = TEfficiency("",";m (GeV);r9;#varepsilon", 200, 40, 140, 60, 0, 1.1);
  m["cIso"] = TEfficiency("",";m (GeV);I_{#pm} (GeV);#varepsilon", 200, 40, 140, 50, 0, 5);
  m["nIso"] = TEfficiency("",";m (GeV);I_{0} (GeV);#varepsilon", 200, 40, 140, 50, 0, 5);
  m["pIso"] = TEfficiency("",";m (GeV);I_{#gamma} (GeV);#varepsilon", 200, 40, 140, 50, 0, 5);
  m["cIsoWorst"] = TEfficiency("",";m (GeV);worst I_{#pm} (GeV);#varepsilon", 200, 40, 140, 80, 0, 40);
  m["nTracksPV"] = TEfficiency("",";m (GeV);PV track multiplicity;#varepsilon", 200, 40, 140, 50, -.5, 49.5);
  m["pt_vs_vtx"] = TEfficiency("", ";m (GeV);#it{p}_{T} (GeV);vertex multiplicity", 200, 40, 140, 40, 0, 200, 36, 0.5, 36.5);
  m["pt_vs_nTracksPV"] = TEfficiency("", ";m (GeV);#it{p}_{T} (GeV);PV track multiplicity", 200, 40, 140, 40, 0, 200, 50, -0.5, 49.5);
  m["vtx_vs_nTracksPV"] = TEfficiency("", ";m (GeV);vertex multiplicity ;PV track multiplicity", 200, 40, 140, 36, .5, 36.5, 50, -0.5, 49.5);
  m["met_vs_vtx"] = TEfficiency("", ";m (GeV);#it{E}_{T}^{miss} (GeV);vertex multiplicity", 200, 40, 140, 40, 0, 200, 36, 0.5, 36.5);
  m["met_vs_jets"] = TEfficiency("", ";m (GeV);#it{E}_{T}^{miss} (GeV);jet multiplicity", 200, 40, 140, 40, 0, 200, 8, -0.5, 7.5);
  m["met_vs_pt"] = TEfficiency("", ";m (GeV);#it{E}_{T}^{miss} (GeV);#it{p}_{T} (GeV)", 200, 40, 140, 40, 0, 200, 40, 0, 200);
  m["met_vs_jetPt"] = TEfficiency("", ";m (GeV);#it{E}_{T}^{miss} (GeV);#it{p}_{T}^{1.jet} (GeV)", 200, 40, 140, 40, 0, 200, 40, 0, 200);
  m["eta_vs_phi"] = TEfficiency("", ";m (GeV);#eta;#phi", 200, 40, 140, 250, -2.5, 2.5, 64, -3.2, 3.2);
  return m;
}

void FakeRateSelector::fillSelection(tree::Electron* tag, tree::Photon* probe, const string& name)
{
  if (tag->p.DeltaR(probe->p)<tpSeperationDRcut) return;

  if (!effMaps.count(name)) {
    effMaps[name] = initEfficiencies();
  }
  auto effs = &effMaps.at(name);

  bool hasPixelSeed = probe->hasPixelSeed;
  float thisMet = met->p.Pt();
  float mll = m(tag->p, probe->p);
  float emht = tag->p.Pt() + probe->p.Pt();
  unsigned nJet = 0;
  float leadingJetPt = 0;
  for (auto& j : *jets ) {
    if (j.p.Pt()>30 && j.p.Eta()<3 && tag->p.DeltaR(j.p)>0.4 && probe->p.DeltaR(j.p)>0.4) {
      if (leadingJetPt<1) leadingJetPt = j.p.Pt();
      nJet++;
      emht += j.p.Pt();
    }
  }

  effs->at("pt").Fill(!hasPixelSeed, mll, probe->p.Pt());
  effs->at("pt").Fill(!hasPixelSeed, mll, probe->p.Pt());
  effs->at("eta").Fill(!hasPixelSeed, mll, fabs(probe->p.Eta()));
  effs->at("vtx").Fill(!hasPixelSeed, mll, *nGoodVertices);
  effs->at("emht").Fill(!hasPixelSeed, mll, emht);
  effs->at("jets").Fill(!hasPixelSeed, mll, nJet);
  effs->at("met").Fill(!hasPixelSeed, mll, thisMet);
  effs->at("sie").Fill(!hasPixelSeed, mll, probe->sigmaIetaIeta);
  effs->at("sip").Fill(!hasPixelSeed, mll, probe->sigmaIphiIphi);
  effs->at("hoe").Fill(!hasPixelSeed, mll, probe->hOverE);
  effs->at("r9").Fill(!hasPixelSeed, mll, probe->r9);
  effs->at("cIso").Fill(!hasPixelSeed, mll, probe->cIso);
  effs->at("nIso").Fill(!hasPixelSeed, mll, probe->nIso);
  effs->at("pIso").Fill(!hasPixelSeed, mll, probe->pIso);
  effs->at("cIsoWorst").Fill(!hasPixelSeed, mll, probe->cIsoWorst);
  effs->at("nTracksPV").Fill(!hasPixelSeed, mll, *nTracksPV);

  effs->at("pt_vs_vtx").Fill(!hasPixelSeed, mll, probe->p.Pt(), *nGoodVertices);
  effs->at("pt_vs_nTracksPV").Fill(!hasPixelSeed, mll, probe->p.Pt(), *nTracksPV);
  effs->at("vtx_vs_nTracksPV").Fill(!hasPixelSeed, mll, *nGoodVertices, *nTracksPV);
  effs->at("met_vs_vtx").Fill(!hasPixelSeed, mll, thisMet, *nGoodVertices);
  effs->at("met_vs_jets").Fill(!hasPixelSeed, mll, thisMet, nJet);
  effs->at("met_vs_pt").Fill(!hasPixelSeed, mll, thisMet, probe->p.Pt());
  effs->at("met_vs_jetPt").Fill(!hasPixelSeed, mll, thisMet, leadingJetPt);
  effs->at("eta_vs_phi").Fill(!hasPixelSeed, mll, probe->p.Eta(), probe->p.Phi());
}

Bool_t FakeRateSelector::Process(Long64_t entry)
{
  if (!(entry%int(2e6))) cout << 1.*entry / fReader.GetEntries(false) << endl;
  fReader.SetEntry(entry);
  if (!*hlt && isData) return true;


  ////////////////////////////////////////////
  // Gen selection
  ///////////////////////////////////////////
  if (!isData) {
    for (auto& pho : *photons) {
      auto eta = fabs(pho.p.Eta());
      auto pt = pho.p.Pt();
      if (!pho.isLoose) continue;
      if (photonsEtaMaxBarrel < eta and eta < photonsEtaMinEndcap) continue;
      if (fabs(genMatchWZDecay(pho, *intermediateGenParticles)) != 11) continue;
      fillGenSelection(pho, "all_gen");
      if (pt>40) fillGenSelection(pho, "40pt_gen");
      if (eta<photonsEtaMaxBarrel) {
        fillGenSelection(pho, "EB_gen");
        if (pt>40) {
          fillGenSelection(pho, "EB_40pt_gen");
          if (eta>0.1) fillGenSelection(pho, "EB_01eta_40pt_gen");
        }
        if (pt>100) fillGenSelection(pho, "EB_100pt_gen");
      }
    }
  }

  // search tag
  vector<tree::Electron*> selTags;
  for (auto& el : *electrons) {
    if (el.p.Pt() < 30 || fabs(el.p.Eta())>2.1 || !el.isTight) continue;
    for (const auto& tObj : *triggerObjects) {
      if (el.p.DeltaR(tObj.p)<triggerDRcut) {
        selTags.push_back(&el);
        break;
      }
    }
    if (not isData) selTags.push_back(&el); // trigger not simulated
  }
  if (!selTags.size()) return true;

  // search probe
  vector<tree::Photon*> selProbes;
  for (auto& pho : *photons) {
    if (!pho.isLoose) continue;
    if (pho.p.Pt()<25) continue;
    auto eta = fabs(pho.p.Eta());
    if (photonsEtaMaxBarrel < eta and eta < photonsEtaMinEndcap) continue;
    selProbes.push_back(&pho);
  }

  // search muons
  vector<tree::Photon> artificialPhotons;
  vector<tree::Photon*> selMuons;
  for (auto& mu : *muons) {
    if (mu.p.Pt() < 25) continue;
    auto eta = fabs(mu.p.Eta());
    if (photonsEtaMaxBarrel < eta and eta < photonsEtaMinEndcap) continue;
    tree::Photon x;
    x.p = mu.p;
    x.sigmaIetaIeta = -1;
    x.sigmaIphiIphi = -1;
    x.hOverE = -1;
    x.hasPixelSeed = false;
    x.passElectronVeto = false;
    x.r9 = -1;
    x.sigmaPt = -1;
    x.cIso = -1;
    x.nIso = -1;
    x.pIso = -1;
    x.cIsoWorst = -1;
    artificialPhotons.push_back(x);
    selMuons.push_back(&artificialPhotons.back());
  }
  h1Map.at("n_tags").Fill(selTags.size());
  h1Map.at("n_probes").Fill(selProbes.size());
  h1Map.at("n_muons").Fill(selMuons.size());

  for (auto& tag : selTags) {
    for (auto& probe : selProbes) {
      auto eta = abs(probe->p.Eta());
      fillSelection(tag, probe, "all");
      if (probe->p.Pt()>40) fillSelection(tag, probe, "40pt");
      if (eta<photonsEtaMaxBarrel) {
        fillSelection(tag, probe, "EB");
        if (probe->p.Pt()>40) fillSelection(tag, probe, "EB_40pt");
        if (probe->p.Pt()>40 && *nTracksPV>2 ) fillSelection(tag, probe, "EB_40pt_3nTracks");
        if (eta>0.1 && probe->p.Pt()>40) fillSelection(tag, probe, "EB_01eta_40pt");
      } else if ( eta>photonsEtaMinEndcap && eta<2.3) {
        if (probe->p.Pt()>40) fillSelection(tag, probe, "EE_40pt");
        if (probe->p.Pt()>40 && *nTracksPV>2) fillSelection(tag, probe, "EE_40pt_3nTracks");
      }
      if (fabs(genMatchWZDecay(*tag, *intermediateGenParticles)) == 11 && fabs(genMatchWZDecay(*probe, *intermediateGenParticles)) == 11) {
        fillSelection(tag, probe, "all_Zmatch");
        if (probe->p.Pt()>40) fillSelection(tag, probe, "40pt_Zmatch");
        if (eta<photonsEtaMaxBarrel) {
          fillSelection(tag, probe, "EB_Zmatch");
          if (probe->p.Pt()>40) fillSelection(tag, probe, "EB_40pt_Zmatch");
          if (probe->p.Pt()>40 && *nTracksPV>2 ) fillSelection(tag, probe, "EB_40pt_3nTracks_Zmatch");
          if (eta>0.1 && probe->p.Pt()>40) fillSelection(tag, probe, "EB_01eta_40pt_Zmatch");
        } else if (eta<2.3) {
          if (probe->p.Pt()>40) fillSelection(tag, probe, "EE_40pt_Zmatch");
          if (probe->p.Pt()>40 && *nTracksPV>2) fillSelection(tag, probe, "EE_40pt_3nTracks_Zmatch");
        }
      }
    }
    for (auto& probe : selMuons) {
      auto eta = abs(probe->p.Eta());
      fillSelection(tag, probe, "all_bkg");
      if (probe->p.Pt()>40) fillSelection(tag, probe, "40pt_bkg");
      if (eta<photonsEtaMaxBarrel) {
        fillSelection(tag, probe, "EB_bkg");
        if (probe->p.Pt()>40) fillSelection(tag, probe, "EB_40pt_bkg");
        if (probe->p.Pt()>40 && *nTracksPV>2) fillSelection(tag, probe, "EB_40pt_3nTracks_bkg");
        if (eta>0.1 && probe->p.Pt()>40) fillSelection(tag, probe, "EB_01eta_40pt_bkg");
      } else if (eta<2.3) {
        if (probe->p.Pt()>40) fillSelection(tag, probe, "EE_40pt_bkg");
        if (probe->p.Pt()>40 && *nTracksPV>2) fillSelection(tag, probe, "EE_40pt_3nTracks_bkg");
      }
    }
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

void FakeRateSelector::Terminate()
{
  auto outputName = getOutputFilename(fReader.GetTree()->GetCurrentFile()->GetName(), "fake");
  TFile file(outputName.c_str(), "RECREATE");
  for (auto& it : h1Map ) it.second.Write(it.first.c_str());
  save2File(effMaps, file);
  file.Close();
  cout << "Created " << outputName << " in " << (time(NULL) - startTime)/60 << " min" << endl;
}

