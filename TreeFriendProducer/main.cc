#include <TFile.h>
#include <TTree.h>

int main( int argc, char** argv ) {
    auto file = TFile::Open( "TTJets_V03.31_tree.root" );
    auto tree = (TTree*)file->Get("photonTree");
    tree->Process( "TreeFriendProducer.cc+");
    return 0;
}

