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

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('inputFile' )
    args = parser.parse_args()


    ch = ROOT.TChain("TreeWriter/eventTree")
    ch.AddFile( args.inputFile )

    baseH = {
        "pt": ROOT.TH1F("", ";p_{T} (GeV)", 100, 15, 150 ),
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
        "isTight": ROOT.TH1F("", ";isTight", 2, 0, 2 )
    }

    selections = "_loose", "_loose_eb", "_loose_ee", "_loose_eb_genPhoton", "_loose_ee_genPhoton", "_loose_eb_genElectron", "_loose_ee_genElectron"

    h = {}
    for s in selections:
        for name, hist in baseH.iteritems():
            h[name+s] = hist.Clone()



    for ievent, event in enumerate(ch):
        #if ievent > 1000: break

        for photon in event.photons:

            if photon.isLoose:
                fill( photon, h, "_loose" )
                eta = abs(photon.p.Eta())

                # find out if generated electron
                genElectron = False
                for ele in event.genParticles:
                    if abs(ele.pdgId) == 11 \
                    and ele.p.DeltaR( photon.p ) < 0.2 \
                    and abs( (ele.p.Pt()-photon.p.Pt())/ele.p.Pt() ) < 0.1:
                        genElectron=True
                        break

                if eta<1.4442:
                    fill( photon, h, "_loose_eb" )
                    if photon.isTrueAlternative:
                        fill( photon, h, "_loose_eb_genPhoton" )
                    if genElectron:
                        fill( photon, h, "_loose_eb_genElectron" )

                elif 1.5 < eta and eta < 2.4:
                    fill( photon, h, "_loose_ee" )
                    if photon.isTrueAlternative:
                        fill( photon, h, "_loose_ee_genPhoton" )
                    if genElectron:
                        fill( photon, h, "_loose_ee_genElectron" )

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
