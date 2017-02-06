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
  void fillSelection(string const& s, float w, bool);
  void fillSignalSelection(string const& s, float w);
  void fillSelectedPhotons(const tree::Particle& p);
  tree::Jet* matchedJet(const tree::Particle& p);
  float getPhotonWeight(const tree::Photon& p);

  void initTriggerStudies();
  void fillTriggerStudies();

  void initUncut();
  void fillUncut();

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
  TTreeReaderValue<std::vector<Float_t>> pdf_weights;
  TTreeReaderValue<Int_t> nGoodVertices;
  TTreeReaderValue<Int_t> nTracksPV;
  TTreeReaderValue<Float_t> genHt;
  TTreeReaderValue<Float_t> rho;
  TTreeReaderValue<Int_t> nTruePV;
  TTreeReaderValue<ULong64_t> evtNo;
  TTreeReaderValue<UInt_t> runNo;
  TTreeReaderValue<UInt_t> lumNo;
  TTreeReaderValue<Bool_t> hlt_photon90_ht600;
  TTreeReaderValue<Bool_t> hlt_photon90;
  TTreeReaderValue<Bool_t> hlt_ht600;
  TTreeReaderValue<Bool_t> hlt_ht800;
  TTreeReaderValue<Bool_t> hlt_el27;
  TTreeReaderValue<Int_t> hlt_ht600_pre;

  // signal scan
  TTreeReaderValue<UShort_t> signal_nBinos;
  TTreeReaderValue<UShort_t> signal_m1;
  TTreeReaderValue<UShort_t> signal_m2;

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
  bool genHt600;
  bool noPromptPhotons;
  bool isSignal;
  TH1F cutFlow;
  string inputName;

  CutFlowPhoton looseCutFlowPhoton;
  map<string,Weighter> weighters;
  Resolution resolution;

  double startTime;
  TRandom2 rand;
  ClassDef(HistogramProducer, 1)
};

