#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TFile.h>

#include "TreeObjects.h"

class TreeFriendProducer : public TSelector {

 public:

  TreeFriendProducer(TTree * = 0):
    photons( fReader, "photons" ),
    jets( fReader, "jets" ),
    genPhotons( fReader, "genPhotons" ),
    met( fReader, "met" ) {
      std::cout << " constructed" << std::endl;
      }
  virtual ~TreeFriendProducer() { 
      std::cout << " destructed" << std::endl;
}
  virtual void    Init(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual Bool_t  Process(Long64_t entry);
  virtual void    Terminate();
  virtual Int_t   Version() const { return 0; }

 private:

  TTreeReader fReader;
  TTreeReaderArray<tree::Photon> photons;
  TTreeReaderArray<tree::Jet> jets;
  TTreeReaderArray<tree::Particle> genPhotons;
  TTreeReaderValue<float> met;

  TTree* outputTree;
  float htHLT;

};


void TreeFriendProducer::SlaveBegin(TTree *tree)
{
      std::cout << " slave begin" << std::endl;
  outputTree = new TTree("treeFrined", tree->GetName() );
  outputTree->Branch( "htHLT", &htHLT );

  GetOutputList()->Add( outputTree );
}

Bool_t TreeFriendProducer::Process(Long64_t entry)
{
      std::cout << " process" << std::endl;
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
  std::cout << htHLT << std::endl;

  return kTRUE;
}

void TreeFriendProducer::Init(TTree *tree)
{
      std::cout << " init" << std::endl;
  fReader.SetTree(tree);
}

void TreeFriendProducer::Terminate()
{
      std::cout << " terminated" << std::endl;
  auto outputfile = TFile::Open("output.root", "recreate" );
  outputTree->Write( 0, kOverwrite );
  outputfile->Close();
}
