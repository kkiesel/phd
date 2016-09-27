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

class FourHists {
  // Creates four histograms for {signal,background}x{pixelseed, }
 public:
  FourHists() {}
  FourHists(const string& name, unsigned nBins, float min, float max, unsigned nMaxMixings_=10) {
    nMaxMixings = nMaxMixings_;
    eg_sig = TH2F( ("eg_sig_"+name).c_str(), "", 120, 60, 120, nBins, min, max);
    ee_sig = TH2F( ("ee_sig_"+name).c_str(), "", 120, 60, 120, nBins, min, max);
    eg_bkg = TH2F( ("eg_bkg_"+name).c_str(), "", 120, 60, 120, nBins, min, max);
    ee_bkg = TH2F( ("ee_bkg_"+name).c_str(), "", 120, 60, 120, nBins, min, max);
    binnedTags = std::vector<std::vector<TVector3>>(nBins+2, std::vector<TVector3>(0));
    binnedProbes = std::vector<std::vector<TVector3>>(nBins+2, std::vector<TVector3>(0));
  }
  ~FourHists() { }
  void Fill(const TVector3& tag, const TVector3& probe, float val, bool signal) {
    auto bin = eg_sig.GetYaxis()->FindFixBin(val);
    auto mll = m(tag, probe);
    ee_sig.Fill(mll, val);
    if (signal) eg_sig.Fill(mll, val);
    for (const auto& t : binnedTags.at(bin)) {
      if (t.DeltaR(probe)<tpSeperationDRcut) continue;
      mll = m(t, probe);
      ee_bkg.Fill(mll, val);
      if (signal) eg_bkg.Fill(mll, val);
    }
    for (const auto& p : binnedProbes.at(bin)) {
      if (tag.DeltaR(p)<tpSeperationDRcut) continue;
      mll = m(tag, p);
      ee_bkg.Fill(mll, val);
      if (signal) eg_bkg.Fill(mll, val);
    }
    if (binnedTags.at(bin).size()<nMaxMixings) {
      binnedTags.at(bin).push_back(tag);
      binnedProbes.at(bin).push_back(probe);
    }
  }
  void Write() {
    eg_sig.Write();
    ee_sig.Write();
    eg_bkg.Write();
    ee_bkg.Write();
  }
 private:
  TH2F eg_sig;
  TH2F eg_bkg;
  TH2F ee_sig;
  TH2F ee_bkg;
  vector<vector<TVector3>> binnedTags;
  vector<vector<TVector3>> binnedProbes;
  unsigned nMaxMixings;
};


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

  TTreeReader fReader;
  TTreeReaderValue<std::vector<tree::Photon>> photons;
  TTreeReaderValue<std::vector<tree::Jet>> jets;
  TTreeReaderValue<std::vector<tree::Electron>> electrons;
  TTreeReaderValue<std::vector<tree::Particle>> triggerObjects;
  TTreeReaderValue<tree::MET> met;
  TTreeReaderValue<tree::MET> metRaw;
  TTreeReaderValue<Float_t> pu_weight;
  TTreeReaderValue<Char_t> mc_weight;
  TTreeReaderValue<Int_t> nGoodVertices;
  TTreeReaderValue<ULong64_t> eventNo;
  TTreeReaderValue<UInt_t> runNo;
  TTreeReaderValue<UInt_t> lumNo;
  TTreeReaderValue<Bool_t> hlt;

  map<string,TEfficiency> effs;

  FourHists pt;
  FourHists eta;
  FourHists vtx;

  TTree outTree;
  TVector3 tag;
  TVector3 probe;
  UInt_t nJet, nVertex;
  float thisMet, emht, mll;
  Bool_t hasPixelSeed;

  double startTime;
};


FakeRateSelector::FakeRateSelector():
  photons(fReader, "photons"),
  jets(fReader, "jets"),
  electrons(fReader, "electrons"),
  triggerObjects(fReader, "triggerObjects"),
  met(fReader, "met"),
  metRaw(fReader, "met_raw"),
  nGoodVertices(fReader, "nGoodVertices"),
  pu_weight(fReader, "pu_weight"),
  runNo(fReader, "runNo"),
  hlt(fReader, "HLT_Ele27_eta2p1_WPLoose_Gsf_v"),
  outTree("fakeTree", ""),
  startTime(time(NULL))
{
}

void FakeRateSelector::Init(TTree *tree)
{
  fReader.SetTree(tree);
}

void FakeRateSelector::SlaveBegin(TTree *tree)
{
  outTree.Branch("tag", &tag);
  outTree.Branch("probe", &probe);
  outTree.Branch("mll", &mll);
  outTree.Branch("met", &thisMet);
  outTree.Branch("emht", &emht);
  outTree.Branch("nVertex", &nVertex, "nVertex/i");
  outTree.Branch("nJet", &nJet, "nJet/i");
  outTree.Branch("hasPixelSeed", &hasPixelSeed, "hasPixelSeed/O");

  pt = FourHists("pt", 12, 0, 120);
  eta = FourHists("eta", 300, 0, 3);
  vtx = FourHists("vtx", 36, 0.5, 36.5);

  effs["pt"] = TEfficiency("","", 120, 60, 120, 12, 0, 120);
}



Bool_t FakeRateSelector::Process(Long64_t entry)
{
  //if (entry>1e4) return true;
  if (!(entry%int(2e6))) cout << 1.*entry / fReader.GetEntries(false) << endl;
  fReader.SetEntry(entry);
  if (!*hlt) return true;

  // search tag
  vector<tree::Electron*> selTags;
  for (auto& el : *electrons) {
    if (el.p.Pt() < 30 || fabs(el.p.Eta())>2.1) continue;
    bool triggerMatch = false;
    for (const auto& tObj : *triggerObjects) {
      if (el.p.DeltaR(tObj.p)<triggerDRcut) {
        triggerMatch = true;
        break;
      }
    }
    if (triggerMatch) {
      selTags.push_back(&el);
    }
  }
  if (!selTags.size()) return true;

  // search probe
  vector<tree::Photon*> selProbes;
  for (auto& pho : *photons) {
    if (!pho.isLoose) continue;
    bool tagMatch = false;
    for (const auto& t : selTags) {
      if (pho.p.DeltaR(t->p)<tpSeperationDRcut) {
        tagMatch = true;
        break;
      }
    }
    if (!tagMatch) {
      selProbes.push_back(&pho);
    }
  }
  if (!selProbes.size()) return true;

  // select exactly one tag and exactly one probe
  if (selProbes.size()>1 and selTags.size()>1) {
    float diffZmass = 100;
    tree::Electron *a = 0;
    tree::Photon *b = 0;
    for (unsigned i=0; i<selTags.size(); i++) {
      for (unsigned j=0; j<selProbes.size(); j++) {
        float dz = std::abs( zBosonMass - m(selTags.at(i)->p,selProbes.at(j)->p) ); // maasss
        if (dz<diffZmass) {
          diffZmass = dz;
          a = selTags.at(i);
          b = selProbes.at(j);
        }
      }
    }
    selTags.clear();
    selProbes.clear();
    selTags.push_back(a);
    selProbes.push_back(b);
  }

  tag = selTags.at(0)->p;
  probe = selProbes.at(0)->p;

  hasPixelSeed = selProbes.at(0)->hasPixelSeed;
  nVertex = *nGoodVertices;
  thisMet = met->p.Pt();
  mll = m(tag, probe);
  emht = 0;
  nJet = 0;
  for (auto& j : *jets ) {
    if (j.p.Pt()>30 && j.p.Eta()<3 && tag.DeltaR(j.p)>0.4 && probe.DeltaR(j.p)>0.4) {
      nJet++;
      emht += j.p.Pt();
    }
  }
  outTree.Fill();
  if (abs(probe.Eta())<1.4442) pt.Fill(tag, probe, probe.Pt(), !hasPixelSeed);
  if (probe.Pt()>40) eta.Fill(tag, probe, abs(probe.Eta()), !hasPixelSeed);
  vtx.Fill(tag, probe, nVertex, !hasPixelSeed);
  if (abs(probe.Eta())<1.4442) effs.at("pt").Fill(!hasPixelSeed, mll, probe.Pt());

  return kTRUE;
}

void FakeRateSelector::Terminate()
{
  auto outputName = getOutputFilename(fReader.GetTree()->GetCurrentFile()->GetName(), "fake");
  TFile file(outputName.c_str(), "RECREATE");
  outTree.Write("fakeTree", TObject::kWriteDelete);
  pt.Write();
  eta.Write();
  vtx.Write();
  for (auto& it : effs ) it.second.Write(it.first.c_str());
  file.Close();
  cout << "Created " << outputName << " in " << (time(NULL) - startTime)/60 << " min" << endl;
}

