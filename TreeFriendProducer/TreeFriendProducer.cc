#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TFile.h>

#include "TreeObjects.h"

class TreeFriendProducer : public TSelector {

 private:

  TTreeReader fReader;
  TTreeReaderArray<tree::Photon> photons;
  TTreeReaderArray<tree::Jet> jets;
  TTreeReaderArray<tree::Particle> genPhotons;
  TTreeReaderValue<float> met;

  TreeFriendProducer(TTree * = 0):
    photons( fReader, "photons" ),
    jets( fReader, "jets" ),
    genPhotons( fReader, "genPhotons" ),
    met( fReader, "met" ) {}
  virtual ~TreeFriendProducer() { }
  virtual void    Init(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual Bool_t  Process(Long64_t entry);
  virtual void    Terminate();
  virtual Int_t   Version() const { return 0; }

  TTree* outputTree;
  float htHLT;

};


void TreeFriendProducer::SlaveBegin(TTree *tree)
{
  outputTree = new TTree("treeFrined", tree->GetName() );
  outputTree->Branch( "htHLT", &htHLT );

  GetOutputList()->Add( outputTree );
}

Bool_t TreeFriendProducer::Process(Long64_t entry)
{
  fReader.SetEntry(entry);

  htHLT = 0;

  for( auto j : jets ) {
    if( std::abs( j.eta ) < 3 && j.pt > 40 ) {
      htHLT += j.pt;
    }
  }

  /*for( auto p : photons ) {
    float drMin = 100;
    for( auto j : jets ) {
      float dr = p.DeltaR( j );
      if( dr < drMin ) drMin = dr;
    }
    std::cout << drMin << std::endl;
  }

*/

  outputTree->Fill();

  return kTRUE;
}

void TreeFriendProducer::Init(TTree *tree)
{
  fReader.SetTree(tree);
}

void TreeFriendProducer::Terminate()
{
  auto outputfile = TFile::Open("output.root", "recreate" );
  outputTree->Write( 0, kOverwrite );
  outputfile->Close();
}
