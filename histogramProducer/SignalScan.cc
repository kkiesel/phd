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

  float getPhotonWeight(const tree::Photon& p);
  void fillSignalSelection(const string p, string const& s, float weight);

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
  TTreeReaderValue<tree::MET> met_JESu;
  TTreeReaderValue<tree::MET> met_JESd;
  TTreeReaderValue<tree::MET> met_JERu;
  TTreeReaderValue<tree::MET> met_JERd;
  TTreeReaderValue<Float_t> pu_weight;
  TTreeReaderValue<Char_t> mc_weight;
  TTreeReaderValue<Int_t> nGoodVertices;
  TTreeReaderValue<Float_t> rho;
  TTreeReaderValue<ULong64_t> eventNo;
  TTreeReaderValue<UInt_t> runNo;
  TTreeReaderValue<UInt_t> lumNo;
  TTreeReaderValue<Int_t> nTruePV;

  // signal scan
  TTreeReaderValue<UShort_t> signal_nBinos;
  TTreeReaderValue<UShort_t> signal_m1;
  TTreeReaderValue<UShort_t> signal_m2;

  map<string,map<string,map<string,TH1F>>> h1MapsMaps;
  map<string,unsigned int> nEventMap;

  string inputName;

  vector<tree::Photon*> selPhotons;
  vector<tree::Jet*> selJets;
  vector<tree::Electron*> selElectrons;
  vector<tree::Muon*> selMuons;

  map<string,Weighter> weighters;

  double startTime;
  ClassDef(SignalScan, 1)
};

SignalScan::SignalScan():
  photons(fReader, "photons"),
  jets(fReader, "jets"),
//  electrons(fReader, "electrons"),
//  muons(fReader, "muons"),
//  genJets(fReader, "genJets"),
//  genParticles(fReader, "genParticles"),
//  intermediateGenParticles(fReader, "intermediateGenParticles"),
 met(fReader, "met"),
//  metRaw(fReader, "met_raw"),
  met_JESu(fReader, "met_JESu"),
  met_JESd(fReader, "met_JESd"),
  met_JERu(fReader, "met_JERu"),
  met_JERd(fReader, "met_JERd"),
//  nGoodVertices(fReader, "nGoodVertices"),
  pu_weight(fReader, "pu_weight"),
  mc_weight(fReader, "mc_weight"),
//  rho(fReader, "rho"),
//  runNo(fReader, "runNo"),
  nTruePV(fReader, "true_nPV"),
  signal_nBinos(fReader, "signal_nBinos"),
  signal_m1(fReader, "signal_m1"),
  signal_m2(fReader, "signal_m2"),
  startTime(time(NULL))
{
  weighters["sf_photon_id_loose"] = Weighter("../plotter/data/dataMcScaleFactors_80X.root", "EGamma_SF2D");
  weighters["sf_photon_pixel"] = Weighter("../plotter/data/EleVeto_SFs_80X.root", "Scaling Factors_HasPix_InclusiveR9");
  string puUp = "pileupWeightUp_mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU";
  string puDn = "pileupWeightDown_mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU";
  weighters["puWeightUp"] = Weighter("../../CMSSW/treewriter/CMSSW_8_0_25/src/TreeWriter/PUreweighting/data/puWeights.root", puUp);
  weighters["puWeightDn"] = Weighter("../../CMSSW/treewriter/CMSSW_8_0_25/src/TreeWriter/PUreweighting/data/puWeights.root", puDn);
  weighters.at("sf_photon_id_loose").fillOverflow2d();
  weighters.at("sf_photon_pixel").fillOverflow2d();
}

void SignalScan::Init(TTree *tree)
{
  fReader.SetTree(tree);
}

map<string,TH1F> initSignalHistograms() {
  map<string,TH1F> hMap;
  hMap["met"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_puUp"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_puDn"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_jesUp"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_jesDn"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_jerUp"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_jerDn"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  return hMap;
}

void SignalScan::fillSignalSelection(const string p, string const& s, float weight)
{
  if (std::isnan(met->p.X()) and std::isnan(met->p.Y())) return;

  if (!h1MapsMaps.count(p)) {
    h1MapsMaps[p] = {{"signal_lowEMHT", initSignalHistograms()}, {"signal_highEMHT", initSignalHistograms()}};
  }
  auto m1 = &h1MapsMaps[p][s];
  auto _met = met->p.Pt();
  m1->at("met").Fill(_met, weight);
  m1->at("met_puUp").Fill(_met, weight*weighters.at("puWeightUp").getWeight(*nTruePV)/ *pu_weight);
  m1->at("met_puDn").Fill(_met, weight*weighters.at("puWeightDn").getWeight(*nTruePV)/ *pu_weight);
  m1->at("met_jesUp").Fill(met_JESu->p.Pt(), weight);
  m1->at("met_jesDn").Fill(met_JESd->p.Pt(), weight);
  m1->at("met_jerUp").Fill(met_JERu->p.Pt(), weight);
  m1->at("met_jerDn").Fill(met_JERd->p.Pt(), weight);
}


float SignalScan::getPhotonWeight(const tree::Photon& p) {
  float weight = 1.;
  auto pt = p.p.Pt();
  auto eta = p.p.Eta();
  // sf for id and electron veto
  weight = weighters.at("sf_photon_id_loose").getWeight(eta, pt) * weighters.at("sf_photon_pixel").getWeight(fabs(eta), pt);
  // trigger efficiency
  weight *= fabs(eta)<photonsEtaMaxBarrel ? 0.964 : 0.94;

  return weight;
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

  selPhotons.clear();
  selJets.clear();
  for (auto& photon : *photons) {
    if (photon.isLoose && !photon.hasPixelSeed && photon.p.Pt() > 100 && fabs(photon.p.Eta()) < photonsEtaMaxBarrel) {
      selPhotons.push_back(&photon);
    }
  }
  if (!selPhotons.size()) return kTRUE;

  for (auto& jet : *jets) {
    if (!jet.isLoose
//      || jet.hasPhotonMatch || jet.hasElectronMatch || jet.hasMuonMatch
      || indexOfMatchedParticle<tree::Photon*>(jet, selPhotons, .4) >= 0
      || jet.p.Pt() < 30 || fabs(jet.p.Eta()) > 3) continue;
    selJets.push_back(&jet);
  }
  float emht=0;
  for (auto& p : selPhotons) emht += p->p.Pt();
  for (auto& p : selJets) emht += p->p.Pt();

  if (emht<700) return kTRUE;
  selW *= getPhotonWeight(*selPhotons.at(0));
  if (emht<2000) fillSignalSelection(pn, "signal_lowEMHT", selW);
  else fillSignalSelection(pn, "signal_highEMHT", selW);

  return kTRUE;
}

void cdNewDir(TFile& file, const string& name) {
  if (!file.Get(name.c_str())) {
    file.mkdir(name.c_str());
  }
  file.cd(name.c_str());
}

template<typename T>
void save2File(const map<string,map<string,map<string,T>>>& hMapsMaps, TFile& file)
{
  for (auto& hMapMapIt: hMapsMaps) {
    file.mkdir(hMapMapIt.first.c_str());
    for (auto& hMapIt: hMapMapIt.second) {
      file.mkdir((hMapMapIt.first+"/"+hMapIt.first).c_str());
      file.cd((hMapMapIt.first+"/"+hMapIt.first).c_str());
      for (auto& h: hMapIt.second) {
        h.second.Write(h.first.c_str(), TObject::kWriteDelete);
      }
    }
  }
}

void SignalScan::Terminate()
{
  inputName = fReader.GetTree()->GetCurrentFile()->GetName();
  auto outputName = getOutputFilename(inputName, "signalScan");
  TFile file(outputName.c_str(), "RECREATE");
  save2File(h1MapsMaps, file);
  file.Close();
  cout << "Created " << outputName << " in " << (time(NULL) - startTime)/60 << " min" << endl;
}

