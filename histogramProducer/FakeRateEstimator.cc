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
  TTreeReaderArray<tree::Particle> triggerObjects;
  TTreeReaderValue<tree::MET> met;
  TTreeReaderValue<Float_t> pu_weight;
  TTreeReaderValue<Char_t> mc_weight;
  TTreeReaderValue<Int_t> nGoodVertices;
  TTreeReaderValue<Int_t> nTracksPV;
  TTreeReaderValue<Float_t> genHt;
  TTreeReaderValue<ULong64_t> eventNo;
  TTreeReaderValue<UInt_t> runNo;
  TTreeReaderValue<UInt_t> lumNo;
  TTreeReaderValue<Bool_t> hlt_photon90_ht500;
  TTreeReaderValue<Bool_t> hlt_photon90;
  TTreeReaderValue<Bool_t> hlt_photon175;
  TTreeReaderValue<Bool_t> hlt_ht600;
  TTreeReaderValue<Bool_t> hlt_e22_wp75;
  TTreeReaderValue<Bool_t> hlt_e22_wpl;
  TTreeReaderValue<Bool_t> hlt_e27_wpl;

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
  triggerObjects( fReader, "triggerObjects" ),
  met( fReader, "met" ),
  nGoodVertices( fReader, "nGoodVertices" ),
  nTracksPV( fReader, "nTracksPV" ),
  pu_weight( fReader, "pu_weight" ),
  mc_weight( fReader, "mc_weight" ),
  genHt( fReader, "genHt" ),
  runNo( fReader, "runNo" ),
  hlt_photon90_ht500( fReader, "HLT_Photon90_CaloIdL_PFHT500_v" ),
  hlt_photon90( fReader, "HLT_Photon90_v" ),
  hlt_photon175( fReader, "HLT_Photon175_v" ),
  hlt_ht600( fReader, "HLT_PFHT600_v" ),
  hlt_e22_wp75( fReader, "HLT_Ele22_eta2p1_WP75_Gsf_v" ),
  hlt_e22_wpl( fReader, "HLT_Ele22_eta2p1_WPLoose_Gsf_v" ),
  hlt_e27_wpl( fReader, "HLT_Ele27_eta2p1_WPLoose_Gsf_v" )
{
}

void FakeRateEstimator::Init(TTree *tree)
{
  fReader.SetTree(tree);
  string inputName = fReader.GetTree()->GetCurrentFile()->GetName();
  isData = inputName.find("Run201") != string::npos;

  // x axis is always mee
  eff["pt"] = TEfficiency("",";m (GeV);p_{T} (GeV)", 200, 50, 150, 150, 0, 150 );
  eff["eta"] = TEfficiency("",";m (GeV);|#eta|", 200, 50, 150, 150, 0, 1.5 );
  eff["njet"] = TEfficiency("",";m (GeV);jet multiplicity", 200, 50, 150, 8, -.5, 7.5 );
  eff["nvtx"] = TEfficiency("",";m (GeV);vertex multiplicity", 200, 50, 150, 21, -.5, 20.5 );
  eff["ntrk"] = TEfficiency("",";m (GeV);track multiplicity", 200, 50, 150, 21, -.5, 20.5 );
  eff["ptVSeta"] = TEfficiency("",";m (GeV);p_{T} (GeV);|#eta|", 200, 50, 150, 150, 0, 150, 150, 0, 1.5 );

  eff["hlt_ptLoose_e22_75"] = TEfficiency("","",100, 0, 80 );
  eff["hlt_ptMedium_e22_75"] = TEfficiency("","",100, 0, 80 );
  eff["hlt_ptTight_e22_75"] = TEfficiency("","",100, 0, 80 );
  eff["hlt_ptLoose_e22_l"] = TEfficiency("","",100, 0, 80 );
  eff["hlt_ptMedium_e22_l"] = TEfficiency("","",100, 0, 80 );
  eff["hlt_ptTight_e22_l"] = TEfficiency("","",100, 0, 80 );
  eff["hlt_ptLoose_e27_l"] = TEfficiency("","",100, 0, 80 );
  eff["hlt_ptMedium_e27_l"] = TEfficiency("","",100, 0, 80 );
  eff["hlt_ptTight_e27_l"] = TEfficiency("","",100, 0, 80 );


}

void FakeRateEstimator::SlaveBegin(TTree *tree)
{
}

Bool_t FakeRateEstimator::Process(Long64_t entry)
{
  fReader.SetEntry(entry);
  tag=0;
  probe=0;
  for( auto& electron : electrons ) {
    if( triggerObjects.GetSize() && triggerObjects[0].p.DeltaR(electron.p) < .2 && fabs(electron.p.Eta())<2.1 )
    {
      if( electron.isLoose ) {
        eff["hlt_ptLoose_e22_75"].Fill( *hlt_e22_wp75, electron.p.Pt() );
        eff["hlt_ptLoose_e22_l"].Fill( *hlt_e22_wpl, electron.p.Pt() );
        eff["hlt_ptLoose_e27_l"].Fill( *hlt_e27_wpl, electron.p.Pt() );
      }
      if( electron.isMedium ) {
        eff["hlt_ptMedium_e22_75"].Fill( *hlt_e22_wp75, electron.p.Pt() );
        eff["hlt_ptMedium_e22_l"].Fill( *hlt_e22_wpl, electron.p.Pt() );
        eff["hlt_ptMedium_e27_l"].Fill( *hlt_e27_wpl, electron.p.Pt() );
      }
      if( electron.isTight ) {
        eff["hlt_ptTight_e22_75"].Fill( *hlt_e22_wp75, electron.p.Pt() );
        eff["hlt_ptTight_e22_l"].Fill( *hlt_e22_wpl, electron.p.Pt() );
        eff["hlt_ptTight_e27_l"].Fill( *hlt_e27_wpl, electron.p.Pt() );
      }

      if( electron.isLoose && electron.p.Pt()>40 ) {
        tag = &electron;
        break;
      }
    }
  }
  if(!tag) return kTRUE;

  for( auto& photon : photons ) {
    if( photon.isLoose && photon.p.DeltaR(tag->p)>.4 ) {
      probe = &photon;
      break;
    }
  }
  if(!probe) return kTRUE;

  //int nJet = count_if(jets.begin(), jets.end(), [*probe](const tree::Jet& jet) { return jet.p.Pt()>40 && fabs(jet.p.Eta())<3 && jet.p.DeltaR(probe->p)>.3;});
  int nJet = 0;
  float m = mass( tag->p, probe->p );
  eff.at("pt").Fill( !probe->hasPixelSeed, m, probe->p.Pt() );
  eff.at("eta").Fill( !probe->hasPixelSeed, m, fabs(probe->p.Eta()) );
  eff.at("njet").Fill( !probe->hasPixelSeed, m, nJet );
  eff.at("nvtx").Fill( !probe->hasPixelSeed, m, *nGoodVertices );
  eff.at("ntrk").Fill( !probe->hasPixelSeed, m, *nTracksPV );
  eff.at("ptVSeta").Fill( !probe->hasPixelSeed, m, probe->p.Pt(), fabs(probe->p.Eta()) );

  return kTRUE;
}

void FakeRateEstimator::Terminate()
{
  auto outputName = getOutputFilename( fReader.GetTree()->GetCurrentFile()->GetName(), "fakeRate" );
  TFile file( outputName.c_str(), "RECREATE");

  for( auto& it : eff ) {
    it.second.Write( it.first.c_str(), TObject::kWriteDelete );
  }

}

