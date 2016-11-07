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

  map<string,TEfficiency> effs;
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
  hlt(fReader, "HLT_Ele27_eta2p1_WPLoose_Gsf_v")
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

void FakeRateSelector::fillGenSelection(const tree::Photon& probe, const string& name) {

  if (!effs.count("gen_pt"+name)) {
    effs["gen_pt"+name] = TEfficiency("",";#it{p}_{T} (GeV);#varepsilon", 40, 0, 200);
    effs["gen_met"+name] = TEfficiency("",";#it{E}_{T}^{miss} (GeV);#varepsilon", 40, 0, 200);
    effs["gen_emht"+name] = TEfficiency("",";#it{EMH}_{T} (GeV);#varepsilon", 40, 0, 400);
    effs["gen_eta"+name] = TEfficiency("",";|#eta|;#varepsilon", 300, 0, 3);
    effs["gen_jets"+name] = TEfficiency("",";jet multiplicity;#varepsilon", 21, -0.5, 20.5);
    effs["gen_vtx"+name] = TEfficiency("",";vertex multiplicity;#varepsilon", 36, 0.5, 36.5);
    effs["gen_sie"+name] = TEfficiency("",";#sigma_{i#etai#eta};#varepsilon", 60, 0, 0.03);
    effs["gen_sip"+name] = TEfficiency("",";#sigma_{i#phii#phi};#varepsilon", 80, 0, 0.08);
    effs["gen_hoe"+name] = TEfficiency("",";H/E;#varepsilon", 100, 0, 0.05);
    effs["gen_r9"+name] = TEfficiency("",";r9;#varepsilon", 60, 0, 1.1);
    effs["gen_cIso"+name] = TEfficiency("",";I_{#pm} (GeV);#varepsilon", 50, 0, 5);
    effs["gen_nIso"+name] = TEfficiency("",";I_{0} (GeV);#varepsilon", 50, 0, 5);
    effs["gen_pIso"+name] = TEfficiency("",";I_{#gamma} (GeV);#varepsilon", 50, 0, 5);
    effs["gen_nTracksPV"+name] = TEfficiency("",";PV track multiplicity;#varepsilon", 50, -.5, 49.5);
    effs["gen_cIsoWorst"+name] = TEfficiency("",";worst I_{#pm} (GeV);#varepsilon", 80, 0, 40);
    effs["gen_pt_vs_vtx"+name] = TEfficiency("", ";#it{p}_{T} (GeV);vertex multiplicity", 40, 0, 200, 36, 0.5, 36.5);
    effs["gen_met_vs_vtx"+name] = TEfficiency("", ";#it{E}_{T}^{miss} (GeV);vertex multiplicity", 40, 0, 200, 36, 0.5, 36.5);
  }

  bool hasPixelSeed = probe.hasPixelSeed;
  float thisMet = met->p.Pt();
  auto nVertex = *nGoodVertices;
  unsigned nJet = 0;
  float emht = probe.p.Pt();;
  for (auto& j : *jets ) {
    if (j.p.Pt()>30 && j.p.Eta()<3 && probe.p.DeltaR(j.p)>0.4) {
      nJet++;
      emht += j.p.Pt();
    }
  }

  effs.at("gen_pt"+name).Fill(!hasPixelSeed, probe.p.Pt());
  effs.at("gen_eta"+name).Fill(!hasPixelSeed, fabs(probe.p.Eta()));
  effs.at("gen_vtx"+name).Fill(!hasPixelSeed, nVertex);
  effs.at("gen_emht"+name).Fill(!hasPixelSeed, emht);
  effs.at("gen_jets"+name).Fill(!hasPixelSeed, nJet);
  effs.at("gen_met"+name).Fill(!hasPixelSeed, thisMet);
  effs.at("gen_sie"+name).Fill(!hasPixelSeed, probe.sigmaIetaIeta);
  effs.at("gen_sip"+name).Fill(!hasPixelSeed, probe.sigmaIphiIphi);
  effs.at("gen_hoe"+name).Fill(!hasPixelSeed, probe.hOverE);
  effs.at("gen_r9"+name).Fill(!hasPixelSeed, probe.r9);
  effs.at("gen_cIso"+name).Fill(!hasPixelSeed, probe.cIso);
  effs.at("gen_nIso"+name).Fill(!hasPixelSeed, probe.nIso);
  effs.at("gen_pIso"+name).Fill(!hasPixelSeed, probe.pIso);
  effs.at("gen_cIsoWorst"+name).Fill(!hasPixelSeed, probe.cIsoWorst);
  effs.at("gen_nTracksPV"+name).Fill(!hasPixelSeed, *nTracksPV);

  effs.at("gen_pt_vs_vtx"+name).Fill(!hasPixelSeed, probe.p.Pt(), nVertex);
  effs.at("gen_met_vs_vtx"+name).Fill(!hasPixelSeed, thisMet, nVertex);
}

void FakeRateSelector::fillSelection(tree::Electron* tag, tree::Photon* probe, const string& name) {
  if (tag->p.DeltaR(probe->p)<tpSeperationDRcut) return;

  if (!effs.count("pt__"+name)) {
    effs["pt"+name] = TEfficiency("",";m (GeV);#it{p}_{T} (GeV);#varepsilon", 200, 40, 140, 40, 0, 200);
    effs["met"+name] = TEfficiency("",";m (GeV);#it{E}_{T}^{miss} (GeV);#varepsilon", 200, 40, 140, 40, 0, 200);
    effs["emht"+name] = TEfficiency("",";m (GeV);#it{EMH}_{T} (GeV);#varepsilon", 200, 40, 140, 40, 0, 400);
    effs["eta"+name] = TEfficiency("",";m (GeV);|#eta|;#varepsilon", 200, 40, 140, 300, 0, 3);
    effs["jets"+name] = TEfficiency("",";m (GeV);jet multiplicity;#varepsilon", 200, 40, 140, 21, -0.5, 20.5);
    effs["vtx"+name] = TEfficiency("",";m (GeV);vertex multiplicity;#varepsilon", 200, 40, 140, 36, 0.5, 36.5);
    effs["sie"+name] = TEfficiency("",";m (GeV);#sigma_{i#etai#eta};#varepsilon", 200, 40, 140, 60, 0, 0.03);
    effs["sip"+name] = TEfficiency("",";m (GeV);#sigma_{i#phii#phi};#varepsilon", 200, 40, 140, 80, 0, 0.08);
    effs["hoe"+name] = TEfficiency("",";m (GeV);H/E;#varepsilon", 200, 40, 140, 100, 0, 0.05);
    effs["r9"+name] = TEfficiency("",";m (GeV);r9;#varepsilon", 200, 40, 140, 60, 0, 1.1);
    effs["cIso"+name] = TEfficiency("",";m (GeV);I_{#pm} (GeV);#varepsilon", 200, 40, 140, 50, 0, 5);
    effs["nIso"+name] = TEfficiency("",";m (GeV);I_{0} (GeV);#varepsilon", 200, 40, 140, 50, 0, 5);
    effs["pIso"+name] = TEfficiency("",";m (GeV);I_{#gamma} (GeV);#varepsilon", 200, 40, 140, 50, 0, 5);
    effs["cIsoWorst"+name] = TEfficiency("",";m (GeV);worst I_{#pm} (GeV);#varepsilon", 200, 40, 140, 80, 0, 40);
    effs["nTracksPV"+name] = TEfficiency("",";PV track multiplicity;#varepsilon", 50, -.5, 49.5);
    effs["pt_vs_vtx"+name] = TEfficiency("", ";m (GeV);#it{p}_{T} (GeV);vertex multiplicity", 200, 40, 140, 40, 0, 200, 36, 0.5, 36.5);
    effs["pt_vs_nTracksPV"+name] = TEfficiency("", ";m (GeV);#it{p}_{T} (GeV);PV track multiplicity", 200, 40, 140, 40, 0, 200, 50, -0.5, 49.5);
    effs["met_vs_vtx"+name] = TEfficiency("", ";m (GeV);#it{E}_{T}^{miss} (GeV);vertex multiplicity", 200, 40, 140, 40, 0, 200, 36, 0.5, 36.5);
    effs["met_vs_jets"+name] = TEfficiency("", ";m (GeV);#it{E}_{T}^{miss} (GeV);jet multiplicity", 200, 40, 140, 40, 0, 200, 8, -0.5, 7.5);
    effs["met_vs_pt"+name] = TEfficiency("", ";m (GeV);#it{E}_{T}^{miss} (GeV);#it{p}_{T} (GeV)", 200, 40, 140, 40, 0, 200, 40, 0, 200);
    effs["met_vs_jetPt"+name] = TEfficiency("", ";m (GeV);#it{E}_{T}^{miss} (GeV);#it{p}_{T}^{1.jet} (GeV)", 200, 40, 140, 40, 0, 200, 40, 0, 200);
  }

  bool hasPixelSeed = probe->hasPixelSeed;
  auto nVertex = *nGoodVertices;
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

  effs.at("pt"+name).Fill(!hasPixelSeed, mll, probe->p.Pt());
  effs.at("pt"+name).Fill(!hasPixelSeed, mll, probe->p.Pt());
  effs.at("eta"+name).Fill(!hasPixelSeed, mll, fabs(probe->p.Eta()));
  effs.at("vtx"+name).Fill(!hasPixelSeed, mll, nVertex);
  effs.at("emht"+name).Fill(!hasPixelSeed, mll, emht);
  effs.at("jets"+name).Fill(!hasPixelSeed, mll, nJet);
  effs.at("met"+name).Fill(!hasPixelSeed, mll, thisMet);
  effs.at("sie"+name).Fill(!hasPixelSeed, mll, probe->sigmaIetaIeta);
  effs.at("sip"+name).Fill(!hasPixelSeed, mll, probe->sigmaIphiIphi);
  effs.at("hoe"+name).Fill(!hasPixelSeed, mll, probe->hOverE);
  effs.at("r9"+name).Fill(!hasPixelSeed, mll, probe->r9);
  effs.at("cIso"+name).Fill(!hasPixelSeed, mll, probe->cIso);
  effs.at("nIso"+name).Fill(!hasPixelSeed, mll, probe->nIso);
  effs.at("pIso"+name).Fill(!hasPixelSeed, mll, probe->pIso);
  effs.at("cIsoWorst"+name).Fill(!hasPixelSeed, mll, probe->cIsoWorst);
  effs.at("nTracksPV"+name).Fill(!hasPixelSeed, mll, *nTracksPV);

  effs.at("pt_vs_vtx"+name).Fill(!hasPixelSeed, mll, probe->p.Pt(), nVertex);
  effs.at("pt_vs_nTracksPV"+name).Fill(!hasPixelSeed, mll, probe->p.Pt(), *nTracksPV);
  effs.at("met_vs_vtx"+name).Fill(!hasPixelSeed, mll, thisMet, nVertex);
  effs.at("met_vs_jets"+name).Fill(!hasPixelSeed, mll, thisMet, nJet);
  effs.at("met_vs_pt"+name).Fill(!hasPixelSeed, mll, thisMet, probe->p.Pt());
  effs.at("met_vs_jetPt"+name).Fill(!hasPixelSeed, mll, thisMet, leadingJetPt);
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
      fillGenSelection(pho);
      if (pt>40) fillGenSelection(pho, "_40pt");
      if (eta<photonsEtaMaxBarrel) {
        fillGenSelection(pho, "_EB");
        if (pt>40) {
          fillGenSelection(pho, "_EB_40pt");
          if (eta>0.1) fillGenSelection(pho, "_EB_01eta_40pt");
        }
        if (pt>100) fillGenSelection(pho, "_EB_100pt");
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
      fillSelection(tag, probe, "__");
      if (probe->p.Pt()>40) fillSelection(tag, probe, "__40pt");
      if (abs(probe->p.Eta())<photonsEtaMaxBarrel) fillSelection(tag, probe, "__EB");
      if (abs(probe->p.Eta())<photonsEtaMaxBarrel && probe->p.Pt()>40) fillSelection(tag, probe, "__EB_40pt");
      if (abs(probe->p.Eta())>0.1 && abs(probe->p.Eta())<photonsEtaMaxBarrel && probe->p.Pt()>40) fillSelection(tag, probe, "__EB_01eta_40pt");
      if (abs(probe->p.Eta())<photonsEtaMaxBarrel && probe->p.Pt()>100) fillSelection(tag, probe, "__EB_100pt");
      if (fabs(genMatchWZDecay(*probe, *intermediateGenParticles)) == 11) {
        fillSelection(tag, probe, "__gen");
        if (probe->p.Pt()>40) fillSelection(tag, probe, "__gen_40pt");
        if (abs(probe->p.Eta())<photonsEtaMaxBarrel) fillSelection(tag, probe, "__gen_EB");
        if (abs(probe->p.Eta())<photonsEtaMaxBarrel && probe->p.Pt()>40) fillSelection(tag, probe, "__gen_EB_40pt");
        if (abs(probe->p.Eta())>0.1 && abs(probe->p.Eta())<photonsEtaMaxBarrel && probe->p.Pt()>40) fillSelection(tag, probe, "__gen_EB_01eta_40pt");
        if (abs(probe->p.Eta())<photonsEtaMaxBarrel && probe->p.Pt()>100) fillSelection(tag, probe, "__gen_EB_100pt");
      }
    }
    for (auto& probe : selMuons) {
      fillSelection(tag, probe, "__bkg");
      if (probe->p.Pt()>40) fillSelection(tag, probe, "__40pt_bkg");
      if (abs(probe->p.Eta())<photonsEtaMaxBarrel) fillSelection(tag, probe, "__EB_bkg");
      if (abs(probe->p.Eta())<photonsEtaMaxBarrel && probe->p.Pt()>40) fillSelection(tag, probe, "__EB_40pt_bkg");
      if (abs(probe->p.Eta())>0.1 && abs(probe->p.Eta())<photonsEtaMaxBarrel && probe->p.Pt()>40) fillSelection(tag, probe, "__EB_01eta_40pt_bkg");
      if (abs(probe->p.Eta())<photonsEtaMaxBarrel && probe->p.Pt()>100) fillSelection(tag, probe, "__EB_100pt_bkg");
    }
  }

  return kTRUE;
}

void FakeRateSelector::Terminate()
{
  auto outputName = getOutputFilename(fReader.GetTree()->GetCurrentFile()->GetName(), "fake");
  TFile file(outputName.c_str(), "RECREATE");
  for (auto& it : effs ) it.second.Write(it.first.c_str());
  file.Close();
  cout << "Created " << outputName << endl;
}

