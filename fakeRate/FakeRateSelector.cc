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

  TTreeReader fReader;
  TTreeReaderValue<std::vector<tree::Photon>> photons;
  TTreeReaderValue<std::vector<tree::Jet>> jets;
  TTreeReaderValue<std::vector<tree::Electron>> electrons;
  TTreeReaderValue<std::vector<tree::Particle>> triggerObjects;
  TTreeReaderValue<tree::MET> met;
  TTreeReaderValue<Int_t> nGoodVertices;
  TTreeReaderValue<Bool_t> hlt;

  map<string,TEfficiency> effs;

  bool isData;
};


FakeRateSelector::FakeRateSelector():
  photons(fReader, "photons"),
  jets(fReader, "jets"),
  electrons(fReader, "electrons"),
  triggerObjects(fReader, "triggerObjects"),
  met(fReader, "met"),
  nGoodVertices(fReader, "nGoodVertices"),
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
  effs["pt"] = TEfficiency("",";m (GeV);#it{p}_{T} (GeV);#varepsilon", 120, 60, 120, 12, 0, 120);
  effs["met"] = TEfficiency("",";m (GeV);#it{E}_{T}^{miss} (GeV);#varepsilon", 120, 60, 120, 20, 0, 200);
  effs["emht"] = TEfficiency("",";m (GeV);#it{EMH}_{T} (GeV);#varepsilon", 120, 60, 120, 30, 0, 300);
  effs["eta"] = TEfficiency("",";m (GeV);|#eta|;#varepsilon", 120, 60, 120, 300, 0, 3);
  effs["jets"] = TEfficiency("",";m (GeV);jet multiplicity;#varepsilon", 120, 60, 120, 21, -0.5, 20.5);
  effs["vtx"] = TEfficiency("",";m (GeV);vertex multiplicity;#varepsilon", 120, 60, 120, 36, 0.5, 36.5);
}



Bool_t FakeRateSelector::Process(Long64_t entry)
{
  //if (entry>1e4) return true;
  if (!(entry%int(2e6))) cout << 1.*entry / fReader.GetEntries(false) << endl;
  fReader.SetEntry(entry);
  if (!*hlt && isData) return true;

  // search tag
  vector<tree::Electron*> selTags;
  for (auto& el : *electrons) {
    if (el.p.Pt() < 30 || fabs(el.p.Eta())>2.1) continue;
    bool triggerMatch = !isData;
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

  TVector3 tag = selTags.at(0)->p;
  TVector3 probe = selProbes.at(0)->p;

  bool hasPixelSeed = selProbes.at(0)->hasPixelSeed;
  auto nVertex = *nGoodVertices;
  float thisMet = met->p.Pt();
  float mll = m(tag, probe);
  float emht = 0;
  unsigned nJet = 0;
  for (auto& j : *jets ) {
    if (j.p.Pt()>30 && j.p.Eta()<3 && tag.DeltaR(j.p)>0.4 && probe.DeltaR(j.p)>0.4) {
      nJet++;
      emht += j.p.Pt();
    }
  }

  effs.at("vtx").Fill(!hasPixelSeed, mll, nVertex);
  if (abs(probe.Eta())<1.4442) effs.at("pt").Fill(!hasPixelSeed, mll, probe.Pt());
  if (probe.Pt()>40) effs.at("eta").Fill(!hasPixelSeed, mll, fabs(probe.Eta()));
  effs.at("emht").Fill(!hasPixelSeed, mll, emht);
  effs.at("jets").Fill(!hasPixelSeed, mll, nJet);
  effs.at("met").Fill(!hasPixelSeed, mll, thisMet);

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

