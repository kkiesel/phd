#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import ConfigParser
import ROOT
import math
import argparse
import re
from random import randint
from sys import maxint

# private libs
import ratio
import style
import multiplot

import auxiliary as aux

intLumi = 1280.23 # /pb https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2522.html
intLumi = 1264 # /pb only RunD
intLumi = 594.6-16 # /pb, only okt05 file

def getHistoFromDataset( dataset, name ):
    h0 = None
    for i in range( len(dataset.files) ):
        h = aux.getFromFile( dataset.files[i], name )
        if isinstance( h, ROOT.TH1 ):
            if dataset.xsecs[i]:
                h.Scale( intLumi * dataset.xsecs[i] / dataset.ngens[i] )
            h.SetLineColor( dataset.color )
            h.SetMarkerColor( dataset.color )

        if h0: h0.Add( h )
        else: h0 = h

    return h0

def compare( datasets, name, saveName ):
    m = multiplot.Multiplot()

    for d in datasets:
        h = getHistoFromDataset( d, name )
        if not h.Integral(): continue
        h.Scale( 1./h.Integral() )
        m.add( h, d.label )

    m.Draw()

    aux.save( "compare%s_%s"%(saveName,name) )

def drawH2( dataset, name, savename="test" ):
    x = style.style2d()
    c = ROOT.TCanvas()
    h = getHistoFromDataset( dataset, name )
    h.Draw("colz")
    l = aux.Label(sim=savename!="data")
    aux.save( "h2_%s_%s"%(savename,name) )
    style.defaultStyle()

def subtractH2( dataset_num, dataset_den, name, savename="test" ):
    x = style.style2d()
    x.SetPalette( 1 )
    c = ROOT.TCanvas()
    num = getHistoFromDataset( dataset_num, name )
    den = getHistoFromDataset( dataset_den, name )
    h = num.Clone()
    h.GetZaxis().SetTitle("( Data-Simulation ) / #sigma_{stat}               ")

    for xbin in range(h.GetNbinsX()+2):
        for ybin in range(h.GetNbinsY()+2):
            n = num.GetBinContent(xbin,ybin)
            e_n = num.GetBinError(xbin,ybin)
            d = den.GetBinContent(xbin,ybin)
            e_d = den.GetBinError(xbin,ybin)
            if e_n == 0: e_n = 1.4
            if e_d == 0: e_d = 1.4
            h.SetBinContent(xbin,ybin, (n-d)/math.sqrt(e_n**2 + e_d**2) )

    absMax = max( [ abs(h.GetMaximum()),abs(h.GetMinimum()) ] )
    h.SetMaximum( absMax )
    h.SetMinimum( -absMax )

    h.Draw("colz")
    l = aux.Label()
    aux.save( "h2subtract_%s_%s"%(savename,name) )
    style.defaultStyle()


def compareAll( saveName="test", *datasets ):
    names = aux.getObjectNames( datasets[0].files[0], "" )

    for name in names:
        if name.startswith("h_"):
            compare( datasets, name, saveName )

def drawSameHistogram( saveName, name, data, bkg, additional=[], binning=None ):

    can = ROOT.TCanvas()
    m = multiplot.Multiplot()

    """
    if "dphi" in name:
        m.leg = ROOT.TLegend(.5,.6,.83,.92)
        m.leg.SetFillColor( ROOT.kWhite )

    scale = 1.
    for d in additional:
        h_data = getHistoFromDataset( d, name )
        if True:
            sumBkg = None
            for d in bkg[-1::-1]:
                h = getHistoFromDataset( d, name )
                if sumBkg: sumBkg.Add( h )
                else: sumBkg = h
            if sumBkg.Integral(): scale = h_data.Integral() / sumBkg.Integral()
    """


    for d in bkg[-1::-1]:
        h = d.getHist( name )
        if not h.Integral(): continue
        if binning: h = aux.rebin( h, binning )

        aux.appendFlowBin( h )
        h.SetYTitle( aux.getYAxisTitle( h ) )
        m.addStack( h, d.label )

    dataHist = None
    for d in additional:
        h = d.getHist( name )
        if not h.Integral(): continue
        if binning: h = aux.rebin( h, binning )

        if h.GetLineColor() == ROOT.kBlack: # data
            h.drawOption_ = "ep"
            h.SetMarkerStyle(20)
            h.SetBinErrorOption( ROOT.TH1.kPoisson )
            # TODO: make sure kPoisson works
            dataHist = h
        else:
            h.drawOption_ = "hist"
            h.SetLineWidth(3)

        m.add( h, d.label )

    if m.Draw():

        # ratio
        hsm = m.hists[0].GetStack().Last()
        if dataHist:
            r = ratio.Ratio( "Data/SM", dataHist, hsm )
            r.draw(0,2)

        l = aux.Label()
        aux.save( "sameHistogram%s_%s"%(saveName,name) )
        can.SetLogy()
        aux.save( "sameHistogram%s_%s_log"%(saveName,name) )


def drawSameHistograms( saveName="test", data=None, bkg=[], additional=[] ):
    names = aux.getObjectNames( bkg[0].files[0] )

    #names = ["h_met__tr"] # test plot
    #names = ["h_g_pt__tr"] # test plot
    #names = ["h_dphi_met_g__tr_bit"] # test plot
    #names = [ "h_ht__tr"]

    for name in names:
        if not name.startswith("h_"): continue

        for binningName, binning in aux.getBinnigsFromName( name ).iteritems():
            drawSameHistogram( saveName+"_bin"+binningName, name, data, bkg, additional, binning )

def getProjections( h2, alongX=True ):
    hs = []
    label = h2.GetYaxis().GetTitle()

    for ybin in range( h2.GetNbinsY()+2 ):

        ylow = h2.GetYaxis().GetBinLowEdge(ybin)
        yhigh = h2.GetYaxis().GetBinUpEdge(ybin)
        name = "{} #leq {} < {}".format( ylow, label, yhigh )
        if ybin == 0: name = "{} < {}".format( label, yhigh )
        if ybin == h2.GetNbinsY()+1: name = "{} #leq {}".format( ylow, label )


        h = h2.ProjectionX( name, ybin, ybin )
        h.SetLineColor( ybin+2)
        if h.GetEntries():
            h.Scale( 1./h.GetEntries() )
            hs.append( h )

    return hs


def drawRazor( dataset ):
    h2 = getHistoFromDataset( dataset, "h2_razorPlane" )
    h2.Rebin2D( 1, 20 )
    razorFit = ROOT.TF2("razorFitFunc", "[0]*( [1]*(x[0]-[2])*(x[1]-[3]) - 1 ) * exp( -[1]*(x[0]-[2])*(x[1]-[3]) )", 0, 2000, 0, 0.5 )
    razorFit.SetParameters( h2.GetEntries(), 0.0005, 170, 0.00001 )
    razorFit.FixParameter( 2, 170 )
    fr = h2.Fit( "razorFitFunc" )
    h2.Draw("cont2")
    razorFit.Draw("same")
    aux.save( "razorPlane" )

    pX = getProjections( h2 )
    pX[0].Draw()
    for h in pX[1:]: h.Draw("same")

    leg = ROOT.TLegend( .7, .7, .95, .95 )
    for h in pX: leg.AddEntry( h, h.GetName(), "l" )
    leg.Draw()

    aux.save( "razorAlongX" )

def qcdClosure( dataset, samplename="" ):
    names = aux.getObjectNames( dataset.files[0], "", [ROOT.TH1F] )

    gSet = "tr_reco"
    cSet = "jControl"

    for name in names:
        if cSet not in name: continue


        hdir = getHistoFromDataset( dataset, name )
        hdir.SetLineColor(1)
        hdir.SetMarkerColor(1)
        hdir.SetMarkerStyle(20)
        hdir.drawOption_ = "p"

        hpre = getHistoFromDataset( dataset, name.replace( cSet, gSet ) )
        scale = hdir.Integral() / hpre.Integral()
        hpre.Scale( scale )
        hpre.drawOption_ = "hist"

        for h in hdir, hpre:
            h.SetYTitle( aux.getYAxisTitle( h ) )

        for binningName, binning in getBinnigsFromName( name ).iteritems():
            can = ROOT.TCanvas()
            m = multiplot.Multiplot()
            mod_dir = hdir
            mod_pre = hpre

            if binning:
                mod_dir = aux.rebin(mod_dir, binning)
                mod_pre = aux.rebin(mod_pre, binning)


            if name == "h_g_eta__jControl":
                m.maximum = 7000
                m.minimum = 2000
            m.add( mod_dir, "#gamma" )
            m.add( mod_pre, "#gamma-like" )

            m.Draw()

            l = aux.Label()
            aux.save( "qcdClosure_"+name+samplename+binningName )
            can.SetLogy()
            aux.save( "qcdClosure_"+name+samplename+binningName+"_log" )


def ewkClosure( dataset, samplename="" ):
    names = aux.getObjectNames( dataset.files[0], "", [ROOT.TH1F] )

    gSet = "tr_genElectron"
    cSet = "tr_eControl"

    for name in names:
        if gSet not in name: continue

        hdir = getHistoFromDataset( dataset, name )
        hdir.SetLineColor(1)
        hdir.SetMarkerColor(1)
        hdir.SetMarkerStyle(20)
        hdir.drawOption_ = "p"

        hpre = getHistoFromDataset( dataset, name.replace( gSet, cSet ) )
        scale = hdir.Integral() / hpre.Integral()
        hpre.Scale( scale )
        hpre.drawOption_ = "hist"
        for h in hdir, hpre:
            h.SetYTitle( aux.getYAxisTitle( h ) )

        for binningName, binning in aux.getBinnigsFromName( name ).iteritems():
            can = ROOT.TCanvas()
            m = multiplot.Multiplot()
            mod_dir = hdir
            mod_pre = hpre

            if binning:
                mod_dir = aux.rebin(mod_dir, binning)
                mod_pre = aux.rebin(mod_pre, binning)

            if name == "h_met__tr_reco_genElectron": mod_dir.SetMaximum( mod_dir.GetMaximum() *5 )
            m.add( mod_dir, "#gamma (gen e)" )
            m.add( mod_pre, "{:.2f}% #times #gamma_{{pixel}}".format(100*scale) )

            m.Draw()

            l = aux.Label()
            aux.save( "ewkClosure_"+name+samplename+binningName )
            can.SetLogy()
            aux.save( "ewkClosure_"+name+samplename+binningName+"_log" )


def drawROCs():

    for hname in ["h_gCol_mva__base", "h_g_cIso__base", "h_g_nIso__base", "h_g_pIso__base", "h_g_sigmaIetaIeta__base", "h_g_hOverE__base"]:
        c = ROOT.TCanvas()

        effSig = getHistoFromDataset( gjets, hname )
        effBkg = getHistoFromDataset( qcd, hname )

        a = aux.getROC( effSig.GetPassedHistogram(), effBkg.GetTotalHistogram(), "mva" in hname )
        if not a: continue
        a.Draw()
        a.GetXaxis().SetRangeUser(0.01,1)
        a.GetYaxis().SetRangeUser(0.01,1)
        l = aux.Label()
        aux.save("roc_"+hname)

def efficienciesDataMC( dataset_data, dataset_mc, savename="" ):
    names = aux.getObjectNames( dataset_data.files[0], "", [ROOT.TEfficiency] )

    for name in names:
        h = getHistoFromDataset( dataset_data, name )
        if h.UsesWeights(): h.SetStatisticOption( ROOT.TEfficiency.kFNormal )

        h_pass = h.GetPassedHistogram()
        h_tot = h.GetTotalHistogram()

        h_tot.SetTitle("")

        h_tot.SetLineColor(ROOT.kGray)
        h_tot.SetFillColor(ROOT.kGray)
        h_tot.SetMaximum( 1.3*h_pass.GetMaximum() )
        h_tot.Draw("hist")
        h_pass.SetLineColor(ROOT.kBlack)
        h_pass.Draw("same hist")



        """
        if "hlt" in name:
            eff = ROOT.TEfficiency( h_pass, h_tot )
        else:
            eff = h
        eff.Draw()

        if name == "eff_hlt_ht__base":
            cutValue = 600
        elif name == "eff_hlt_pt__base":
            cutValue = 100
        else:
            cutValue = None

        if cutValue:
            bin = h_pass.FindFixBin( cutValue )
            passed = int(h_pass.Integral( bin, -1 ))
            total = int(h_tot.Integral( bin, -1 ))
            conf = 0.682689492137
            e = 1.*passed/total
            e_up = ROOT.TEfficiency.ClopperPearson( total, passed, conf, True )
            e_dn = ROOT.TEfficiency.ClopperPearson( total, passed, conf, False )
            print "ε = {:.1%} + {:.1%} - {:.1%}".format(e, e_up-e,e-e_dn )

            # graphical representation
            l = ROOT.TLine()
            l.SetLineWidth(2)
            l.SetLineColor( ROOT.kRed )
            xmax = eff.CreateGraph().GetHistogram().GetXaxis().GetXmax()
            l.DrawLine( cutValue, e, xmax, e )

            l.SetLineStyle(2)
            ymin = eff.CreateGraph().GetHistogram().GetYaxis().GetXmin()
            ymax = eff.CreateGraph().GetHistogram().GetYaxis().GetXmax()
            l.DrawLine( cutValue, ymin, cutValue, ymax )


        """
        l = aux.Label()
        aux.save( "efficiencyDataMC_"+savename+name )


def efficiencies( dataset, savename="" ):
    names = aux.getObjectNames( dataset.files[0], "", [ROOT.TEfficiency] )

    for name in names:
        h = getHistoFromDataset( dataset, name )
        if h.UsesWeights(): h.SetStatisticOption( ROOT.TEfficiency.kFNormal )

        h_pas = h.GetPassedHistogram()
        h_tot = h.GetTotalHistogram()
        if "hlt" in name:
            eff = ROOT.TEfficiency( h_pas, h_tot )
        else:
            eff = h

        eff.Draw()
        ROOT.gPad.Update()
        #eff.GetPaintedGraph().GetYaxis().SetRangeUser(0.9, 1.01)

        if name == "eff_hlt_ht__base":
            cutValue = 600
        elif name == "eff_hlt_pt__base":
            cutValue = 100
        else:
            cutValue = None

        if cutValue:
            bin = h_pas.FindFixBin( cutValue )
            passed = int(h_pas.Integral( bin, -1 ))
            total = int(h_tot.Integral( bin, -1 ))
            conf = 0.682689492137
            e = 1.*passed/total
            e_up = ROOT.TEfficiency.ClopperPearson( total, passed, conf, True )
            e_dn = ROOT.TEfficiency.ClopperPearson( total, passed, conf, False )
            print "ε = {:.1%} + {:.1%} - {:.1%}".format(e, e_up-e,e-e_dn )

            # graphical representation
            l = ROOT.TLine()
            l.SetLineWidth(2)
            l.SetLineColor( ROOT.kRed )
            xmax = eff.CreateGraph().GetHistogram().GetXaxis().GetXmax()
            l.DrawLine( cutValue, e, xmax, e )
            l.DrawLine( cutValue, e_up, xmax, e_up )
            l.DrawLine( cutValue, e_dn, xmax, e_dn )

            # cut line
            l.SetLineStyle(2)
            ymin = eff.GetPaintedGraph().GetYaxis().GetXmin()
            ymax = eff.GetPaintedGraph().GetYaxis().GetXmax()
            l.DrawLine( cutValue, ymin, cutValue, ymax )


        l = aux.Label()
        aux.save( "efficiency_"+savename+name )

        h_tot.SetLineColor(2)
        h_tot.Draw("hist")
        h_pas.Draw("same e*")
        aux.save( "efficiency_"+savename+name+"_raw" )

def ewkIsrSamplesSplitting( dataset, isrDataset, saveName="test" ):
    names = aux.getObjectNames( dataset.files[0] )
    names = [ n for n in names if "_genPhoton" in n and n.startswith("h_") ]
    #names = [ "h_g_pt__tr_genPhoton" ]

    for name in names:
        can = ROOT.TCanvas()
        m = multiplot.Multiplot()

        h_tot = dataset.getHist( name.replace("_genPhoton","" ) )
        h_g = dataset.getHist( name )
        h_e = dataset.getHist( name.replace("_genPhoton","_genElectron" ) )

        h_isr_tot = isrDataset.getHist( name.replace("_genPhoton","" ) )
        h_isr_g = isrDataset.getHist( name )
        h_isr_e = isrDataset.getHist( name.replace("_genPhoton","_genElectron" ) )

        for h in h_tot, h_isr_tot:
            h.SetMarkerStyle(20)
            h.SetMarkerSize(0.5)
            h.drawOption_="ep"

        for h in h_g, h_isr_g:
            h.SetLineStyle(2)
            h.drawOption_="hist"

        for h in h_e, h_isr_e:
            h.SetLineStyle(3)
            h.drawOption_="hist"

        m.add( h_tot, dataset.label )
        m.add( h_g, "#gamma" )
        m.add( h_e, "e" )
        m.add( h_isr_tot, isrDataset.label )
        m.add( h_isr_g, "#gamma" )
        m.add( h_isr_e, "e" )

        if m.Draw():

            l = aux.Label()
            aux.save( "ewkIsrSampleSplitting_%s_%s"%(saveName,name) )
            can.SetLogy()
            aux.save( "ewkIsrSampleSplitting_%s_%s_log"%(saveName,name) )




def getRazorProjection( dataset, xaxis, cut ):
    h2 = getHistoFromDataset( dataset, "h2_razorPlane__tr" )
    if xaxis:
        h = h2.ProjectionX("_px", h2.GetYaxis().FindFixBin( cut ), -1)
    else:
        h = h2.ProjectionY("_py", h2.GetXaxis().FindFixBin( cut ), -1)
    return h


def compareRazorProjections( bkg, additional, xaxis=True, cut=-1, binning=[] ):

    can = ROOT.TCanvas()
    m = multiplot.Multiplot()

    for d in bkg[-1::-1]:
        h = getRazorProjection( d, xaxis,cut )
        if binning: h = aux.rebin( h, binning )
        aux.appendFlowBin( h )
        h.SetYTitle( aux.getYAxisTitle( h ) )
        if not h.Integral(): continue
        m.addStack( h, d.label )


    for d in additional:
        h = getRazorProjection( d, xaxis,cut )
        if not h.Integral(): continue
        if binning: h = aux.rebin( h, binning )

        if h.GetLineColor() == ROOT.kBlack: # data
            h.drawOption_ = "ep"
            h.SetMarkerStyle(20)
            h.SetBinErrorOption( ROOT.TH1.kPoisson )
        else:
            h.drawOption_ = "hist"
            h.SetLineWidth(3)

        m.add( h, d.label )

    if m.Draw():

        l = aux.Label()
        aux.save( "razorProjections" )
        can.SetLogy()
        aux.save( "razorProjections_log" )

def razorStudies():
    compareRazorProjections( [gjets,qcd,ttjets,wjets], [t5wg_1500_125], False )


def transitions():
    drawSameHistogram( "_gjets", "h_genHt", [], [gjets40,gjets100,gjets200,gjets400,gjets600] )
    drawSameHistogram( "_gjets_inc15", "h_genHt", [], [gjets_pt15,gjets40,gjets100], binning=range(0,120) )
    drawSameHistogram( "_qcd", "h_genHt", [], [qcd100,qcd200,qcd300,qcd500,qcd700,qcd1000,qcd1500,qcd2000] )
    drawSameHistogram( "_wjets", "h_genHt", [], [wjets100,wjets200,wjets400,wjets600,wjets800,wjets1200,wjets2500] )
    drawSameHistogram( "_wjets_trans1", "h_genHt", [], [wjets100,wjets200], binning=range(90,410,1) )
    drawSameHistogram( "_wjets_trans2", "h_genHt", [], [wjets200,wjets400], binning=range(190,610,1) )
    drawSameHistogram( "_wjets_trans3", "h_genHt", [], [wjets400,wjets600], binning=range(390,810,1) )
    drawSameHistogram( "_wjets_trans4", "h_genHt", [], [wjets600,wjets800], binning=range(590,1210,2) )
    drawSameHistogram( "_wjets_trans5", "h_genHt", [], [wjets800,wjets1200], binning=range(790,2510,2) )
    drawSameHistogram( "_wjets_trans6", "h_genHt", [], [wjets1200,wjets2500], binning=range(1190,3000,2) )
    drawSameHistogram( "_wjets_trans6", "h_genHt", [], [wjets1200,wjets2500], binning=range(1200,3000,20) )
    drawSameHistogram( "_wjets_trans45", "h_genHt", [], [wjets600,wjets800,wjets1200], binning=range(590,1510,5) )


def main():
    #transitions()
    #compareAll( "_all", gjets400, gjets600, znn400, znn600 )
    #compareAll( "_GjetsVsZnn", gjets, znn )
    #compareAll( "_allMC", gjets, znn, qcd, wjets )
    #drawSameHistograms( "_gqcd_data", bkg=[ gjets, qcd], additional=[data])
    #drawSameHistograms( "_gjet15_data", bkg=[gjets_pt15, qcd], additional=[data])
    #drawSameHistograms( "_mc_data", bkg=[gjets, qcd, ttjets, ttg, wjets, dy], additional=[data,t5wg_1500_125, t5wg_1500_1475 ])
    #drawSameHistograms( "_mc_data", bkg=[gjets, qcd, wjets, ttjets, ttg], additional=[data,t5wg_1500_125, t5gg_1000_200])
    #drawSameHistograms( "_mc", bkg=[gjets, qcd, wjets, ttjets, ttg], additional=[t5wg_1500_125, t5gg_1000_200])
    #drawSameHistograms( "_QCD", bkg=[ qcd2000, qcd1500, qcd1000, qcd700, qcd500, qcd300 ] )
    #drawRazor( ttjets )

    #ewkClosure( ttjets, "_tt" )
    #ewkClosure( wjets, "_w" )
    #ewkClosure( wjets+ttjets, "_ewk" )
    #qcdClosure( qcd+gjets, "_qcd-gjets" )

    #efficiencies( ttjets+qcd+gjets+wjets, "allMC_" )
    #efficiencies( qcd+gjets, "gqcd_" )
    #efficiencies( data, "singlePhoton_" )
    #efficiencies( dataHt, "jetHt_" )
    #efficienciesDataMC( dataHt, ttjets+qcd+gjets+wjets, "jetHt_mc_" )

    #drawROCs()
    ewkIsrSamplesSplitting( ttjets, ttg, "tt" )

    """
    for h2name in aux.getObjectNames( data.files[0], objects=[ROOT.TH2]):
        break
        drawH2( data, h2name, "data" )
        drawH2( qcd+gjets, h2name, "gqcd" )
        subtractH2( data, qcd+gjets, h2name, "data-gqcd" )
    """

    #razorStudies()

if __name__ == "__main__":
    from datasets import *
    main()

