#!/usr/bin/env python

import ROOT

ROOT.gSystem.Load( "pluginTreeWriterTreeWriterAuto.so" )

class BaseHistograms:
    def __init__( self, nameSuffix="" ):
        self.h_met = ROOT.TH1F( "", ";E^{miss}_{T} (GeV)", 100, 0, 1000 )
        self.h_mt_g_met = ROOT.TH1F( "", ";m_{T}(p_{T},E^{miss}_{T}) (GeV)", 100, 0, 1000 )
        self.h_metAndL = ROOT.TH1F( "", ";E^{miss}_{T} #vec{+} l (GeV)", 100, 0, 1000 )

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

        self.h_bj1_pt = ROOT.TH1F( "", ";p_{T}^{1.b-jet} (GeV)", 100, 0, 1500 )
        self.h_bj1_eta = ROOT.TH1F( "", ";#eta^{1.b-jet}", 100, -3, 3 )
        self.h_bj2_pt = ROOT.TH1F( "", ";p_{T}^{2.b-jet} (GeV)", 100, 0, 1500 )
        self.h_bj2_eta = ROOT.TH1F( "", ";#eta^{2.b-jet}", 100, -3, 3 )
        self.h_bj3_pt = ROOT.TH1F( "", ";p_{T}^{3.b-jet} (GeV)", 100, 0, 1500 )
        self.h_bj3_eta = ROOT.TH1F( "", ";#eta^{3.b-jet}", 100, -3, 3 )


        self.h_dphi_met_g = ROOT.TH1F( "", ";#Delta#phi(E_{T}^{miss},#gamma)", 100, -3.5, 3.5 )
        self.h_dphi_met_j1 = ROOT.TH1F( "", ";#Delta#phi(E_{T}^{miss},1.jet)", 100, -3.5, 3.5 )
        self.h_dphi_met_j2 = ROOT.TH1F( "", ";#Delta#phi(E_{T}^{miss},2.jet)", 100, -3.5, 3.5 )

        # multiplicities
        self.h_n_vertex = ROOT.TH1F( "", ";Vertices", 61, -0.5, 60.5 )
        self.h_n_photon = ROOT.TH1F( "", ";Photons", 4, -0.5, 3.5 )
        self.h_n_jet = ROOT.TH1F( "", ";Jets", 13, -0.5, 12.5 )
        self.h_n_bjet = ROOT.TH1F( "", ";B Jets", 13, -0.5, 12.5 )
        self.h_n_electron = ROOT.TH1F( "", ";Electrons", 4, -0.5, 3.5 )
        self.h_n_muon = ROOT.TH1F( "", ";Muons", 4, -0.5, 3.5 )

        # jet matching
        self.h2_match_jet_photon = ROOT.TH2F( "", ";#Delta R;p_{T}^{jet}/p_{T}^{#gamma}", 100, 0, 0.5, 100, 0, 4 )
        self.h2_match_jet_electron = ROOT.TH2F( "", ";#Delta R;p_{T}^{jet}/p_{T}^{e}", 100, 0, 0.5, 100, 0, 4 )
        self.h2_match_jet_muon = ROOT.TH2F( "", ";#Delta R;p_{T}^{jet}/p_{T}^{#mu}", 100, 0, 0.5, 100, 0, 4 )

        self.profile_g_pt_met = ROOT.TProfile( "", ";E^{miss}_{T} (GeV);#LTp_{T}#GT (GeV)", 100, 0, 1000 )
        self.profile_ht_met = ROOT.TProfile( "", ";E^{miss}_{T} (GeV);#LTH_{T}#GT (GeV)", 100, 0, 1000 )
        self.profile_jets_met = ROOT.TProfile( "", ";E^{miss}_{T} (GeV);#LTJets#GT", 100, 0, 1000 )

        for name, obj in self.__dict__.items():
            obj.SetName( name+nameSuffix )

    def fill( self, event ):

        ht = 0
        ht_g = 0
        st = 0
        bjets = []

        for jet in event.jets:
            if jet.bDiscriminator>0:
                bjets.append( jet )


        for jet in event.jets:
            ht += jet.p.Pt()

        st += event.met.p.Pt()

        metAndL = event.met.p
        for electron in event.electrons:
            metAndL = metAndL + electron.p
        for muon in event.muons:
            metAndL = metAndL + muon.p


        self.h_mt_g_met.Fill( ( event.photons.at(0).p + event.met.p ).Pt() )
        self.h_dphi_met_g.Fill( event.met.p.DeltaPhi( event.photons.at(0).p ) )
        self.h_dphi_met_j1.Fill( event.met.p.DeltaPhi( event.jets.at(0).p ) )
        self.h_dphi_met_j2.Fill( event.met.p.DeltaPhi( event.jets.at(1).p ) )

        self.h_metAndL.Fill( metAndL.Pt(), event.pu_weight )

        self.h_met.Fill( event.met.p.Pt(), event.pu_weight )
        self.h_ht.Fill( ht, event.pu_weight )
        self.h_ht_g.Fill( ht_g, event.pu_weight )
        self.h_st.Fill( st, event.pu_weight )

        self.h_g_pt.Fill( event.photons.at(0).p.Pt(), event.pu_weight )
        self.h_g_eta.Fill( event.photons.at(0).p.Eta(), event.pu_weight )

        if event.jets.size() > 0:
            self.h_j1_pt.Fill( event.jets.at(0).p.Pt(), event.pu_weight )
            self.h_j1_eta.Fill( event.jets.at(0).p.Eta(), event.pu_weight )
        if event.jets.size() > 1:
            self.h_j2_pt.Fill( event.jets.at(1).p.Pt(), event.pu_weight )
            self.h_j2_eta.Fill( event.jets.at(1).p.Eta(), event.pu_weight )
        if event.jets.size() > 2:
            self.h_j3_pt.Fill( event.jets.at(2).p.Pt(), event.pu_weight )
            self.h_j3_eta.Fill( event.jets.at(2).p.Eta(), event.pu_weight )

        if len(bjets) > 0:
            self.h_bj1_pt.Fill( bjets[0].p.Pt(), event.pu_weight )
            self.h_bj1_eta.Fill( bjets[0].p.Eta(), event.pu_weight )
        if len(bjets) > 1:
            self.h_bj2_pt.Fill( bjets[1].p.Pt(), event.pu_weight )
            self.h_bj2_eta.Fill( bjets[1].p.Eta(), event.pu_weight )
        if len(bjets)> 2:
            self.h_bj3_pt.Fill( bjets[2].p.Pt(), event.pu_weight )
            self.h_bj3_eta.Fill( bjets[2].p.Eta(), event.pu_weight )

        self.h_n_vertex.Fill( event.nGoodVertices, event.pu_weight )
        self.h_n_photon.Fill( event.photons.size(), event.pu_weight )
        self.h_n_jet.Fill( event.jets.size(), event.pu_weight )
        self.h_n_bjet.Fill( len(bjets), event.pu_weight )
        self.h_n_electron.Fill( event.electrons.size(), event.pu_weight )
        self.h_n_muon.Fill( event.muons.size(), event.pu_weight )

        # matching
        for jet in event.jets:
            for photon in event.photons:
                self.h2_match_jet_photon.Fill( photon.p.DeltaR( jet.p ), jet.p.Pt()/photon.p.Pt(), event.pu_weight )
            for electron in event.electrons:
                self.h2_match_jet_electron.Fill( electron.p.DeltaR( jet.p ), jet.p.Pt()/electron.p.Pt(), event.pu_weight )
            for muon in event.muons:
                self.h2_match_jet_muon.Fill( muon.p.DeltaR( jet.p ), jet.p.Pt()/muon.p.Pt(), event.pu_weight )


        # profiles
        self.profile_g_pt_met.Fill( event.met.p.Pt(), event.photons.at(0).p.Pt() )
        self.profile_g_pt_met.Fill( event.met.p.Pt(), ht )
        self.profile_jets_met.Fill( event.met.p.Pt(), event.jets.size() )

    def save( self, filename ):
        f = ROOT.TFile( filename, "update" )
        f.cd()
        for name, obj in self.__dict__.items():
            obj.Write( "", ROOT.TObject.kWriteDelete )
        f.Close()


def getOutputFilename( inputfilename ):
    tmpDir = "./"
    nameAppendix = "_histos"

    import os
    basename = os.path.basename( inputfilename )

    return basename.replace( ".root", nameAppendix+".root" )


def makeHistograms( filename ):
    chain = ROOT.TChain( "TreeWriter/eventTree" )
    chain.AddFile( filename )

    base = BaseHistograms()
    for event in chain:
        base.fill( event )
    base.save( getOutputFilename( filename ) )






if __name__ == "__main__":
    filenames = [ "photon_ntuple_mva_mini_38.root" ]
    filenames = [ "../TreeFriendProducer/photon_ntuple_mva_mini_38.root" ]
    filenames = [ "/user/kiesel/nTuples/photon_ntuple_mva_mini_tmp.root" ]

    for filename in filenames:
        makeHistograms( filename )


