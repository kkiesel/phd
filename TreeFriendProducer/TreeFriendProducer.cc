#include "TreeFriendProducer.h"

void TreeFriendProducer::Init(TTree *tree)
{
   fReader.SetTree(tree);
}

void TreeFriendProducer::SlaveBegin(TTree *tree)
{
    outputTree = new TTree("treeFrined", tree->GetName() );
        outputTree->Branch( "ht4030", &ht4030 );

        GetOutputList()->Add( outputTree );
    }

    Bool_t TreeFriendProducer::Process(Long64_t entry)
    {
       fReader.SetEntry(entry);
       ht4030 = 0;

       for( auto j : jets ) {
           if( std::abs( j.eta) < 3 && j.pt> 40 ) {
               ht4030 += j.pt;
           }
       }

       outputTree->Fill();

       return kTRUE;
    }

    void TreeFriendProducer::Terminate()
    {
        auto outputfile = TFile::Open("output.root", "recreate" );
        outputTree->Write( 0, kOverwrite );
    outputfile->Close();
}
