#!/usr/bin/env python3

import ROOT

class BaseHistograms:
    def __init__( self, nameSuffix="" ):
        self.h_met = ROOT.TH1F( "", ";E^{miss}_T (GeV)", 100, 0, 1000 )

        self.h_ht = ROOT.TH1F( "", ";H_{T}", 200, 0, 2000 )
        self.h_ht_g = ROOT.TH1F( "", ";#gammaH_{T}", 200, 0, 2000 )
        self.h_st = ROOT.TH1F( "", ";S_{T}", 200, 0, 2000 )

        self.h_g_pt = ROOT.TH1F( "", ";p_{T} (GeV)", 100, 0, 1500 )
        self.h_g_eta = ROOT.TH1F( "", ";#eta", 100, -3, 3 )

        self.h_j1_pt = ROOT.TH1F( "", ";p_{T}^{1.jet} (GeV)", 100, 0, 1500 )
        self.h_j1_eta = ROOT.TH1F( "", ";#eta^{1.jet}", 100, -3, 3 )
        self.h_j2_pt = ROOT.TH1F( "", ";p_{T}^{2.jet} (GeV)", 100, 0, 1500 )
        self.h_j2_eta = ROOT.TH1F( "", ";#eta^{2.jet}", 100, -3, 3 )
        self.h_j3_pt = ROOT.TH1F( "", ";p_{T}^{3.jet} (GeV)", 100, 0, 1500 )
        self.h_j3_eta = ROOT.TH1F( "", ";#eta^{3.jet}", 100, -3, 3 )

        # multiplicities
        self.h_n_vertex = ROOT.TH1F( "", ";Vertices", 61, -0.5, 60.5 )
        self.h_n_photon = ROOT.TH1F( "", ";Photons", 4, -0.5, 3.5 )
        self.h_n_jet = ROOT.TH1F( "", ";Jets", 13, -0.5, 12.5 )
        self.h_n_electron = ROOT.TH1F( "", ";Muons", 4, -0.5, 3.5 )
        self.h_n_muon = ROOT.TH1F( "", ";Muons", 4, -0.5, 3.5 )

        for name, obj in self.__dict__.items():
            obj.SetName( name+nameSuffix )

    def fill( self, event ):

        self.h_met.Fill( event.met.p.Pt(), event.pu_weight )

        ht = 0
        ht_g = 0
        st = 0

        for jet in event.jets:
            ht += jet.p.Pt()

        st += event.met.p.Pt()

        self.h_ht.Fill( ht, event.pu_weight )
        self.h_ht_g.Fill( ht_g, event.pu_weight )
        self.h_st.Fill( st, event.pu_weight )

        self.h_g_pt.Fill( event., event.pu_weight )
        self.h_g_eta.Fill( event., event.pu_weight )

        self.h_j1_pt.Fill( event., event.pu_weight )
        self.h_j1_eta.Fill( event., event.pu_weight )
        self.h_j2_pt.Fill( event., event.pu_weight )
        self.h_j2_eta.Fill( event., event.pu_weight )
        self.h_j3_pt.Fill( event., event.pu_weight )
        self.h_j3_eta.Fill( event., event.pu_weight )

        self.h_n_vertex.Fill( event.nGoodVertices, event.pu_weight )
        self.h_n_photon.Fill( event.photons.size(), event.pu_weight )
        self.h_n_jet.Fill( event.jets.size(), event.pu_weight )
        self.h_n_electron.Fill( event.electrons.size(), event.pu_weight )
        self.h_n_muon.Fill( event.muons.size(), event.pu_weight )


def makeHistograms( filename ):
    ch = ROOT.TChain( "TreeWriter/eventTree" )
    ch.AddFile( filename )

    base = BaseHistograms()
    for e in ch:
        base.h_n_vertex.Fill( e.nGoodVertices, e.pu_weight )







if __name__ == "__main__":
    filenames = [ "photon_ntuple_mva_mini_38.root" ]

    for filename in filenames:
        makeHistograms( filename )


