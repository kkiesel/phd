#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TFile.h>

#include "TreeObjects.h"

class TreeFriendProducer : public TSelector {

public :

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
   virtual Int_t   Version() const { return 2; }

   TTree* outputTree;
   float ht4030;

};


