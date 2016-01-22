#!/usr/bin/env python2
import ROOT
import argparse
import re

ROOT.gSystem.Load("../../tools/templates/TreeParticles.so")

def fill( p, h, appendix ):
    h["pt"+appendix].Fill( p.p.Pt() )
    h["eta"+appendix].Fill( abs(p.p.Eta()) )
    h["mvaValue"+appendix].Fill( p.mvaValue )
    h["r9"+appendix].Fill( p.r9 )
    h["sie"+appendix].Fill( p.sigmaIetaIeta )
    h["hoe"+appendix].Fill( p.hOverE )
    h["hasPixelSeed"+appendix].Fill( p.hasPixelSeed )
    h["passEVeto"+appendix].Fill( p.passElectronVeto )
    h["cIso"+appendix].Fill( p.isoChargedHadronsEA )
    h["nIso"+appendix].Fill( p.isoNeutralHadronsEA )
    h["pIso"+appendix].Fill( p.isoPhotonsEA )
    h["isLoose"+appendix].Fill( p.isLoose  )
    h["isMedium"+appendix].Fill( p.isMedium  )
    h["isTight"+appendix].Fill( p.isTight )

def isFake( photon ):
    return not photon.isLoose \
        and photon.hOverE<0.05 \
        and photon.isoNeutralHadronsEA < 1.06+0.014*photon.p.Pt()+0.000019*(photon.p.Pt())**2 \
        and photon.isoPhotonsEA < 0.28+0.0053*photon.p.Pt() \
        and not photon.hasPixelSeed \
        and ( \
        ( 0.012 < photon.sigmaIetaIeta and photon.sigmaIetaIeta <0.014 ) \
        != ( 1.37 < photon.isoChargedHadronsEA and photon.isoChargedHadronsEA < 15 ) )


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('inputFile' )
    args = parser.parse_args()


    ch = ROOT.TChain("TreeWriter/eventTree")
    ch.AddFile( args.inputFile )

    baseH = {
        "pt": ROOT.TH1F("", ";p_{T} (GeV)", 100, 0, 3000 ),
        "eta": ROOT.TH1F("", ";|#eta|", 100, 0, 2.6 ),
        "mvaValue": ROOT.TH1F("", ";y_{MVA}", 100, -1.01, 1.01 ),
        "r9": ROOT.TH1F("", ";r9", 100, 0, 1.01 ),
        "sie": ROOT.TH1F("", ";#sigma_{i#etai#eta}", 100, 0, 0.03 ),
        "hoe": ROOT.TH1F("", ";H/E", 100, 0, 0.05 ),
        "hasPixelSeed": ROOT.TH1F("", ";has pixel seed", 2, 0, 2 ),
        "passEVeto": ROOT.TH1F("", ";pass conversion safe electron veto", 2, 0, 2 ),
        "cIso": ROOT.TH1F("", ";I_{#pi}", 100, 0, 3 ),
        "nIso": ROOT.TH1F("", ";I_{n}", 100, 0, 7 ),
        "pIso": ROOT.TH1F("", ";I_{#gamma}", 100, 0, 5 ),
        "isLoose": ROOT.TH1F("", ";isLoose", 2, 0, 2 ),
        "isMedium": ROOT.TH1F("", ";isMedium", 2, 0, 2 ),
        "isTight": ROOT.TH1F("", ";isTight", 2, 0, 2 ),
    }

    selections = [
        "_loose_genPhoton", "_loose_eb_genPhoton", "_loose_ee_genPhoton",
        "_medium_genPhoton", "_medium_eb_genPhoton", "_medium_ee_genPhoton",
        "_tight_genPhoton", "_tight_eb_genPhoton", "_tight_ee_genPhoton",
        "_fake", "_fake_eb", "_fake_ee",
    ]
    h = {}
    for s in selections:
        for name, hist in baseH.iteritems():
            h[name+s] = hist.Clone()

    h["h3_pt_eta_nV_gen"] = ROOT.TH3F("","p_{T} (GeV);|#eta|;vertex multiplicity", 30, 0, 3000, 100, 0, 2.6, 36, -0.5, 35.5 )
    import copy
    h["h3_pt_eta_nV_loose"] = copy.copy( h["h3_pt_eta_nV_gen"] )
    h["h3_pt_eta_nV_medium"] = copy.copy( h["h3_pt_eta_nV_gen"] )
    h["h3_pt_eta_nV_tight"] = copy.copy( h["h3_pt_eta_nV_gen"] )
    h["h3_pt_eta_nV_fake"] = copy.copy( h["h3_pt_eta_nV_gen"] )


    for ievent, event in enumerate(ch):
        #if ievent > 1000: break

        for photon in event.photons:
            eta = abs(photon.p.Eta())

            if photon.isLoose and photon.isTrueAlternative:
                fill( photon, h, "_loose_genPhoton" )
                if eta<1.4442:
                    fill( photon, h, "_loose_eb_genPhoton" )
                elif 1.5 < eta and eta < 2.4:
                    fill( photon, h, "_loose_ee_genPhoton" )

            if photon.isMedium and photon.isTrueAlternative:
                fill( photon, h, "_medium_genPhoton" )
                if eta<1.4442:
                    fill( photon, h, "_medium_eb_genPhoton" )
                elif 1.5 < eta and eta < 2.4:
                    fill( photon, h, "_medium_ee_genPhoton" )

            if photon.isTight and photon.isTrueAlternative:
                fill( photon, h, "_tight_genPhoton" )
                if eta<1.4442:
                    fill( photon, h, "_tight_eb_genPhoton" )
                elif 1.5 < eta and eta < 2.4:
                    fill( photon, h, "_tight_ee_genPhoton" )

            if not photon.isLoose \
                and photon.hOverE<0.05 \
                and photon.isoNeutralHadronsEA < 1.06+0.014*photon.p.Pt()+0.000019*(photon.p.Pt())**2 \
                and photon.isoPhotonsEA < 0.28+0.0053*photon.p.Pt() \
                and not photon.hasPixelSeed \
                and ( \
                ( 0.012 < photon.sigmaIetaIeta and photon.sigmaIetaIeta <0.015 ) \
                != ( 1.37 < photon.isoChargedHadronsEA and photon.isoChargedHadronsEA < 15 ) ):

                fill( photon, h, "_fake" )
                if eta<1.4442:
                    fill( photon, h, "_fake" )
                elif 1.5 < eta and eta < 2.4:
                    fill( photon, h, "_fake" )

        for genP in event.genParticles:
            if genP.pdgId != 22: continue
            pt = genP.p.Pt()
            eta = abs(genP.p.Eta())
            nv = event.nGoodVertices
            h["h3_pt_eta_nV_gen"].Fill( pt, eta, nv )
            for photon in event.photons:
                if photon.p.DeltaR( genP.p ) > 0.1 or (photon.p.Pt()-pt)/pt > 0.1: continue

                if photon.isLoose: h["h3_pt_eta_nV_loose"].Fill( pt, eta, nv )
                if photon.isMedium: h["h3_pt_eta_nV_medium"].Fill( pt, eta, nv )
                if photon.isTight: h["h3_pt_eta_nV_tight"].Fill( pt, eta, nv )
                if isFake(photon): h["h3_pt_eta_nV_fake"].Fill( pt, eta, nv )



    outFileName = "out.root"
    m = re.match( ".*/([^/]*).root", args.inputFile )
    if m: outFileName = m.group(1)+"_hists.root"


    inFile = ROOT.TFile( args.inputFile )
    cutFlow = inFile.Get("TreeWriter/hCutFlow")
    cutFlow.SetDirectory(0)

    out = ROOT.TFile( outFileName, "recreate")
    out.cd()
    for name, hist in h.iteritems():
        hist.Write( name, ROOT.TObject.kWriteDelete )
    cutFlow.Write( "", ROOT.TObject.kWriteDelete )

    out.Close()
