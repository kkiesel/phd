#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from include import *

def drawSameHistogram(sampleNames, name, bkg=[], additional=[], binning=None, binningName=""):
    can = ROOT.TCanvas()
    m = multiplot.Multiplot()

    dataIntegral = 0
    for d in additional:
        if d == data or d == dataHt:
            h = d.getHist(name)
            dataIntegral = h.Integral(0,-1)
    bkgHist = None
    for d in bkg:
        h = d.getHist(name)
        if bkgHist: bkgHist.Add(h)
        else: bkgHist = h
    bkgIntegral = bkgHist.Integral(0,-1)

    for d in bkg[-1::-1]:
        h = d.getHist(name)
        if not h: continue
        if not h.Integral(): continue
        if True:
            h.Scale(dataIntegral/bkgIntegral)
        if binning: h = aux.rebin( h, binning )

        aux.appendFlowBin( h )
        h.SetYTitle( aux.getYAxisTitle( h ) )
        m.addStack( h, d.label )

    dataHist = None
    for d in additional:
        h = d.getHist( name )
        if not h: continue
        if not h.Integral(): continue
        if binning: h = aux.rebin( h, binning )
        aux.appendFlowBin( h )

        if h.GetLineColor() == ROOT.kBlack: # data
            h.drawOption_ = "ep"
            h.SetMarkerStyle(20)
            # disable errors for data, so that ErrorOption is working
            # TODO: kPoisson also for rebinned and scaled histograms
            if not binning: h.Sumw2(False)
            h.SetBinErrorOption( ROOT.TH1.kPoisson )
            dataHist = h
        else:
            h.drawOption_ = "hist"
            h.SetLineWidth(3)

        m.add( h, d.label )

    m.sortStackByIntegral()
    if m.Draw():

        # ratio
        hsm = m.hists[0].GetStack().Last()
        if dataHist:
            r = ratio.Ratio( "Data/SM", dataHist, hsm )
            rMean = dataHist.Integral()/hsm.Integral()
            if rMean > 2 or rMean < 0.25: # for jcontrol with prescaled data
                r.draw(rMean/2,1.5*rMean)
            else:
                r.draw(.5,1.5)

        l = aux.Label(sim=data not in additional)

        if binningName: binningName = "_"+binningName
        name = name.replace("/","__")
        saveName = "sameHistograms_{}_{}{}".format(sampleNames, name, binningName )
        aux.save( saveName )


def drawSameHistograms( sampleNames="test", stack=[], additional=[] ):
    file = stack[0].files[0] if stack else additional[0].files[0]
    names = aux.getObjectNames( file, "tr", [ROOT.TH1F] )
    additionalHt = [ x for x in additional if x is not data ]
    if data in additional: additionalHt += [ dataHt ]

    if data in additional:
        if "genMatch" in names: names.remove("genMatch")
        if "genHt" in names: names.remove("genHt")

    for name in names:
        for binningName, binning in aux.getBinningsFromName( name ).iteritems():
            drawSameHistogram( sampleNames, "tr_dPhi3/"+name, stack, additional, binning, binningName )
            drawSameHistogram( sampleNames, "tr_jControl_dPhi3/"+name, stack, additional, binning, binningName )
            drawSameHistogram( sampleNames, "tr/"+name, stack, additional, binning, binningName )
            drawSameHistogram( sampleNames, "tr/"+name, stack, additional, binning, binningName )
            drawSameHistogram( sampleNames, "tr_tight/"+name, stack, additional, binning, binningName )
            drawSameHistogram( sampleNames, "tr_eControl/"+name, stack, additional, binning, binningName )
            drawSameHistogram( sampleNames, "tr_0met100/"+name, stack, additional, binning, binningName )
            #drawSameHistogram( sampleNames, "tr_100met/"+name, stack, additional, binning, binningName )
            drawSameHistogram( sampleNames, "tr_jControl/"+name, stack, additionalHt, binning, binningName )
            drawSameHistogram( sampleNames, "tr_jControl_noLep/"+name, stack, additionalHt, binning, binningName )
            drawSameHistogram( sampleNames, "tr_jControl_highHt_neutralEM9/"+name, stack, additionalHt, binning, binningName )
            drawSameHistogram( sampleNames, "tr_jControl_highHt_neutralEM7/"+name, stack, additionalHt, binning, binningName )

def transitions():
    drawSameHistogram( "gjets", "tr_jControl/genHt", [gjets40,gjets100,gjets200,gjets400,gjets600] )
    drawSameHistogram( "qcd", "tr_jControl/genHt", [qcd100,qcd200,qcd300,qcd500,qcd700,qcd1000,qcd1500,qcd2000] )
    drawSameHistogram( "wjets", "tr_jControl/genHt", [wjets100,wjets200,wjets400,wjets600,wjets800,wjets1200,wjets2500] )
    drawSameHistogram( "wjets_tr_jControlans1", "tr_jControl/genHt", [wjets100,wjets200], binning=range(90,410,1) )
    drawSameHistogram( "wjets_tr_jControlans2", "tr_jControl/genHt", [wjets200,wjets400], binning=range(190,610,1) )
    drawSameHistogram( "wjets_tr_jControlans3", "tr_jControl/genHt", [wjets400,wjets600], binning=range(390,810,1) )
    drawSameHistogram( "wjets_tr_jControlans4", "tr_jControl/genHt", [wjets600,wjets800], binning=range(590,1210,2) )
    drawSameHistogram( "wjets_tr_jControlans5", "tr_jControl/genHt", [wjets800,wjets1200], binning=range(790,2510,2) )
    drawSameHistogram( "wjets_tr_jControlans6", "tr_jControl/genHt", [wjets1200,wjets2500], binning=range(1190,3000,2) )
    drawSameHistogram( "wjets_tr_jControlans6", "tr_jControl/genHt", [wjets1200,wjets2500], binning=range(1200,3000,20) )
    drawSameHistogram( "wjets_tr_jControlans45", "tr_jControl/genHt", [wjets600,wjets800,wjets1200], binning=range(590,1510,5) )
    #drawSameHistogram( "wg_500", "h_g_pt__tr", [wg_pt500], [wg_mg], scaleToData=False )
    #drawSameHistogram( "_zg", "h_g_pt__tr", [zg_130], [znunu], scaleToData=False )

if __name__ == "__main__":
    #transitions()
    #drawSameHistograms( "gqcd_data", [gjets, qcd], additional=[data])
    #drawSameHistograms( "mc_data", [gjets, qcd, ttjets, ttg, wjets, wg_mg, zg, znunu], additional=[data])
    #drawSameHistograms( "mc", [gjets, qcd, ttjets, ttg, wjets, wg_mg, zg, znunu], additional=[])

