#include "TSelector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TFile.h"
#include "TEfficiency.h"

#include "TreeParticles.hpp"
#include "UserFunctions.h"

float mass( const TVector3& p1, const TVector3& p2 ) {
  return sqrt(2*(p1.Mag()*p2.Mag()-(p1*p2)));
}

class FakeRateEstimator : public TSelector {
 public:

  FakeRateEstimator();
  virtual ~FakeRateEstimator() { }

  virtual void    Init(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual Bool_t  Process(Long64_t entry);
  virtual void    Terminate();
  virtual Int_t   Version() const { return 2; }

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
  TTreeReaderValue<Float_t> genHt;
  TTreeReaderValue<ULong64_t> eventNo;
  TTreeReaderValue<UInt_t> runNo;
  TTreeReaderValue<UInt_t> lumNo;
  TTreeReaderValue<Bool_t> hlt_photon90_ht500;
  TTreeReaderValue<Bool_t> hlt_photon90;
  TTreeReaderValue<Bool_t> hlt_photon175;
  TTreeReaderValue<Bool_t> hlt_ht600;

  tree::Electron* tag;
  tree::Photon* probe;
  map<string,TEfficiency> eff;

  float selW=1.; // weight

  bool isData;

};

FakeRateEstimator::FakeRateEstimator():
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
  genHt( fReader, "genHt" ),
  runNo( fReader, "runNo" ),
  hlt_photon90_ht500( fReader, "HLT_Photon90_CaloIdL_PFHT500_v" ),
  hlt_photon90( fReader, "HLT_Photon90_v" ),
  hlt_photon175( fReader, "HLT_Photon175_v" ),
  hlt_ht600( fReader, "HLT_PFHT600_v" )
{
}

void FakeRateEstimator::Init(TTree *tree)
{
  fReader.SetTree(tree);
  string inputName = fReader.GetTree()->GetCurrentFile()->GetName();
  isData = inputName.find("Run201") != string::npos;

  // x axis is always mee
  eff["pt"] = TEfficiency("","", 200, 50, 150, 150, 0, 150 );
  eff["eta"] = TEfficiency("","", 200, 50, 150, 150, 0, 1.5 );
  eff["njet"] = TEfficiency("","", 200, 50, 150, 8, -.5, 7.5 );
  eff["nvtx"] = TEfficiency("","", 200, 50, 150, 21, -.5, 20.5 );
  eff["ptVSeta"] = TEfficiency("","", 200, 50, 150, 150, 0, 150, 150, 0, 1.5 );

}

void FakeRateEstimator::SlaveBegin(TTree *tree)
{
}

Bool_t FakeRateEstimator::Process(Long64_t entry)
{
  fReader.SetEntry(entry);
  for( auto& electron : electrons ) {
    if( electron.isLoose
       //&& electron.p.DeltaR( /* trigger matched electron */ ) < .1
      ) {
      tag = &electron;
      break;
    }
  }
  for( auto& photon : photons ) {
    if( photon.isLoose && photon.p.DeltaR(tag->p)>.4 ) {
      probe = &photon;
      break;
    }
  }
  //int nJet = count_if(jets.begin(), jets.end(), [*probe](const tree::Jet& jet) { return jet.p.Pt()>40 && fabs(jet.p.Eta())<3 && jet.p.DeltaR(probe->p)>.3;});
  int nJet = 0;
  float m = mass( tag->p, probe->p );
  eff.at("pt").Fill( !probe->hasPixelSeed, m, probe->p.Pt() );
  eff.at("eta").Fill( !probe->hasPixelSeed, m, fabs(probe->p.Eta()) );
  eff.at("njet").Fill( !probe->hasPixelSeed, m, nJet );
  eff.at("nvtx").Fill( !probe->hasPixelSeed, m, *nGoodVertices );
  eff.at("ptVSeta").Fill( !probe->hasPixelSeed, m, probe->p.Pt(), fabs(probe->p.Eta()) );

  return kTRUE;
}

void FakeRateEstimator::Terminate()
{
  auto outputName = getOutputFilename( fReader.GetTree()->GetCurrentFile()->GetName(), "fakeRate" );
  cout << "Created " << outputName;
  TFile file( outputName.c_str(), "RECREATE");

  for( auto& it : eff ) {
    it.second.Write( it.first.c_str(), TObject::kWriteDelete );
  }

}

