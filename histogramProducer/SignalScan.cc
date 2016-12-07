#include <math.h>
#include <regex>
#include <time.h>

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
#include "Resolution.h"


class SignalScan : public TSelector {
 public:

  SignalScan();
  virtual ~SignalScan() { }

  virtual void Init(TTree *tree);
  virtual void SlaveBegin(TTree *tree){};
  virtual Bool_t Process(Long64_t entry);
  virtual void Terminate();
  virtual Int_t Version() const { return 2; }

  TTreeReader fReader;
  TTreeReaderValue<std::vector<tree::Photon>> photons;
  TTreeReaderValue<std::vector<tree::Jet>> jets;
  TTreeReaderValue<std::vector<tree::Electron>> electrons;
  TTreeReaderValue<std::vector<tree::Muon>> muons;
  TTreeReaderValue<std::vector<tree::Particle>> genJets;
  TTreeReaderValue<std::vector<tree::GenParticle>> genParticles;
  TTreeReaderValue<std::vector<tree::IntermediateGenParticle>> intermediateGenParticles;
  TTreeReaderValue<tree::MET> met;
  TTreeReaderValue<tree::MET> metRaw;
  TTreeReaderValue<Float_t> pu_weight;
  TTreeReaderValue<Char_t> mc_weight;
  TTreeReaderValue<Int_t> nGoodVertices;
  TTreeReaderValue<Float_t> genHt;
  TTreeReaderValue<Float_t> rho;
  TTreeReaderValue<ULong64_t> eventNo;
  TTreeReaderValue<UInt_t> runNo;
  TTreeReaderValue<UInt_t> lumNo;
  TTreeReaderValue<Bool_t> hlt_photon90_ht600;
  TTreeReaderValue<Bool_t> hlt_photon90;
  TTreeReaderValue<Bool_t> hlt_ht600;
  TTreeReaderValue<Bool_t> hlt_ht800;
  TTreeReaderValue<Int_t> hlt_ht600_pre;

  // signal scan
  TTreeReaderValue<UShort_t> signal_nBinos;
  TTreeReaderValue<UShort_t> signal_m1;
  TTreeReaderValue<UShort_t> signal_m2;

  map<string,map<string,TH2F>> h2Maps;
  map<string,unsigned int> nEventMap;

  string inputName;

  vector<tree::Photon*> selPhotons;
  vector<tree::Jet*> selJets;
  vector<tree::Electron*> selElectrons;
  vector<tree::Muon*> selMuons;

  double startTime;
  ClassDef(SignalScan, 1)
};

SignalScan::SignalScan():
  photons(fReader, "photons"),
  jets(fReader, "jets"),
  electrons(fReader, "electrons"),
  muons(fReader, "muons"),
  genJets(fReader, "genJets"),
  genParticles(fReader, "genParticles"),
  intermediateGenParticles(fReader, "intermediateGenParticles"),
  met(fReader, "met"),
  metRaw(fReader, "met_raw"),
  nGoodVertices(fReader, "nGoodVertices"),
  pu_weight(fReader, "pu_weight"),
  mc_weight(fReader, "mc_weight"),
  genHt(fReader, "genHt"),
  rho(fReader, "rho"),
  runNo(fReader, "runNo"),
  hlt_photon90_ht600(fReader, "HLT_Photon90_CaloIdL_PFHT600_v"),
  hlt_photon90(fReader, "HLT_Photon90_v"),
  hlt_ht600(fReader, "HLT_PFHT600_v"),
  hlt_ht800(fReader, "HLT_PFHT800_v"),
  hlt_ht600_pre(fReader, "HLT_PFHT600_v_pre"),
  signal_nBinos(fReader, "signal_nBinos"),
  signal_m1(fReader, "signal_m1"),
  signal_m2(fReader, "signal_m2"),
  startTime(time(NULL))
{
}

void SignalScan::Init(TTree *tree)
{
  fReader.SetTree(tree);
}

map<string,TH2F> initHistograms() {
  map<string,TH2F> hMap;
  hMap["met_vs_emht"] = TH2F("", ";#it{E}_{T}^{miss} (GeV);#it{EMH}_{T} (GeV)", 300, 0, 3000, 450, 500, 5000);
  return hMap;
}


Bool_t SignalScan::Process(Long64_t entry)
{
  if (!(entry%int(2e6))) cout << 1.*entry / fReader.GetEntries(false) << endl;
  fReader.SetLocalEntry(entry);
  auto pn = getSignalPointName(*signal_nBinos, *signal_m1, *signal_m2);

  if (!nEventMap.count(pn)) {
    unsigned int nGen = 0;
    string cutFlowName = "TreeWriter/hCutFlow";
    inputName = fReader.GetTree()->GetCurrentFile()->GetName();
    if (inputName.find("T5") != string::npos) {
      cutFlowName += "T5Wg";
    } else if (inputName.find("T6") != string::npos) {
      cutFlowName += "T6Wg";
    } else {
      cout << "Could not find correct cutflowname for " << inputName << endl;
    }
    cutFlowName += "_"+to_string(*signal_m1);
    cutFlowName += "_"+to_string(*signal_m2);
    TH1F* cutFlow = (TH1F*)fReader.GetTree()->GetCurrentFile()->Get(cutFlowName.c_str());
    if (cutFlow) {
      nGen = cutFlow->GetBinContent(2);
    } else {
      cout << "Could not read cutFlow histogram " << cutFlowName << endl;
    }
    if (inputName.find("T6") != string::npos or inputName.find("T5") != string::npos) {
      switch (*signal_nBinos) {
        case 0: nGen /= 4; break;
        case 1: nGen /= 2; break;
        case 2: nGen /= 4; break;
        default: cout << "Do not know what to do with " << *signal_nBinos << " binos" << endl;
      }
    } else {
      cout << "do stuff for other scans" << endl;
    }
    nEventMap[pn] = nGen;
  }

  auto selW = *mc_weight * *pu_weight/nEventMap[pn];
  if (!h2Maps.count(pn)) {
    h2Maps[pn] = initHistograms();
  }
  auto m2 = &h2Maps[pn];

  selPhotons.clear();
  selElectrons.clear();
  selMuons.clear();
  selJets.clear();
  for (auto& photon : *photons) {
    if (photon.isLoose && !photon.hasPixelSeed && photon.p.Pt() > 100 && fabs(photon.p.Eta()) < photonsEtaMaxBarrel) {
      selPhotons.push_back(&photon);
    }
  }
  if (!selPhotons.size()) return kTRUE;
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
  }
  float myHt=0;
  for (auto& p : selPhotons) myHt += p->p.Pt();
  for (auto& p : selJets) myHt += p->p.Pt();

  if (myHt<700) return kTRUE;
  m2->at("met_vs_emht").Fill(met->p.Pt(), myHt, selW);

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

void SignalScan::Terminate()
{
  inputName = fReader.GetTree()->GetCurrentFile()->GetName();
  auto outputName = getOutputFilename(inputName, "signalScan");
  TFile file(outputName.c_str(), "RECREATE");
  save2File(h2Maps, file);
  file.Close();
  cout << "Created " << outputName << " in " << (time(NULL) - startTime)/60 << " min" << endl;
}

