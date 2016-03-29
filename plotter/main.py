#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys
if sys.version_info[:2] == (2,6):
    print "Initialize correct python version first!"
    sys.exit()

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

intLumi = 2.26e3 # https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2577.html
#intLumi = 2110.588 # https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2544.html

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

    # draw linear fit
    if name in ["h2_match_photon_pt_jet_pt_base"]:
        #h = aux.diagonalFlip( h )
        #h.Draw("colz")
        h.Fit("pol1")

    aux.save( "h2_%s_%s"%(savename,name) )
    c.SetLogz()
    aux.save( "h2_%s_%s_log"%(savename,name) )

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


def divideDatasetIntegrals( numerator, denominator, name ):
    numMerged = sum(numerator)
    h_num = numMerged.getHist( name )
    num = h_num.Integral(0,-1)
    denMerged = sum(denominator)
    h_den = denMerged.getHist( name )
    den = h_den.Integral(0,-1)
    return num/den if den else 1.

def drawSameHistogram( sampleNames, name, bkg=[], additional=[], binning=None, binningName="", scaleToData=False ):
    if name == "h_j_looseID__base": return

    can = ROOT.TCanvas()
    m = multiplot.Multiplot()

    scale = 1.
    if scaleToData: scale = divideDatasetIntegrals( [ i for i in additional if "Data" in i.label ], bkg, name )
    if "__trPhoton90" in name: scale = 1/90.
    if "__tr_jControl" in name and "Jet" in name: scale = 1/6.

    for d in bkg[-1::-1]:
        h = d.getHist( name )
        if not h.Integral(): continue
        h.Scale(scale)
        if binning: h = aux.rebin( h, binning )

        aux.appendFlowBin( h )
        h.SetYTitle( aux.getYAxisTitle( h ) )
        m.addStack( h, d.label )

    dataHist = None
    for d in additional:
        h = d.getHist( name )
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

    if name == "h_met__tr" and binningName == "3": m.minimum = 1e-3
    if name == "h_met__tr_eControl" and binningName == "3": m.minimum = 1e-1
    if name == "h_met__tr_jControlLeadingJet" and binningName == "1": m.minimum = 0.9

    if name == "h_ht__tr" and binningName == "3": m.minimum = 1e-3
    if name == "h_ht__tr_eControl" and binningName == "1": m.minimum = 1e-1
    if name == "h_ht__tr_jControlLeadingJet" and binningName == "1": m.minimum = 1e-2

    if name == "h_n_jet__tr" and binningName == "3": m.minimum = 1
    if name == "h_n_jet__tr_eControl" and binningName == "1": m.minimum = 1
    if name == "h_n_jet__tr_jControlLeadingJet" and binningName == "1": m.minimum = 1

    if name == "h_g_pt__tr" and binningName == "3": m.minimum = 1e-4
    if name == "h_j1_pt__tr" and binningName == "3": m.minimum = 1e-3

    if "h_met__tr" in name:
        m.leg.SetX1(0.45)
        m.leg.SetY1(0.44)
    else:
        m.leg.SetX1(1)
        m.leg.SetY1(1)
        m.leg.SetX2(1)
        m.leg.SetY2(1)

    if "jControl" in name: m.sortStackByIntegral()
    ROOT.gPad.SetLogy()
    if m.Draw():

        # ratio
        hsm = m.hists[0].GetStack().Last()
        if dataHist:
            r = ratio.Ratio( "Data/SM", dataHist, hsm )
            r.draw(0.5,1.5)

            if name  == "h_ht__tr" and binningName == "1":
                aux.writeWeight( r.ratio, name, sampleNames )

        info = ""
        l = aux.Label(info="#scale[0.7]{%s}"%info, sim=data not in additional)

        if binningName: binningName = "_"+binningName
        saveName = "sameHistograms_{}_{}{}".format(sampleNames, name, binningName )
        aux.save( saveName )
        can.SetLogy()
        aux.save( saveName+"_log" )


def drawSameHistograms( sampleNames="test", stack=[], additional=[] ):
    file = stack[0].files[0] if stack else additional[0].files[0]
    names = aux.getObjectNames( file )

    #names = [ "h_met__tr_jControlLeadingJet" ] # test plot
    names = ["h_met__tr_eControl"]



    for name in names:
        if not name.startswith("h_"): continue

        for binningName, binning in aux.getBinnigsFromName( name ).iteritems():
            drawSameHistogram( sampleNames, name, stack, additional, binning, binningName )

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
    x = style.style2d()
    c = ROOT.TCanvas()
    h2 = getHistoFromDataset( dataset, "h2_razorPlane__tr" )
    h2 = aux.rebin2d( h2, range(0, 1000,80), aux.drange(0,0.5, 10) )
    h2.Draw("colz")
    aux.save("razorPlane")
    style.defaultStyle()
    """
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
    """

def multiQcdClosure( dataset, controlDataset, name, samplename, binning, binningName, dirDir="tr", preDir="tr_jControl" ):
    if not controlDataset: controlDataset = dataset
    can = ROOT.TCanvas()
    m = multiplot.Multiplot()


    hdir = dataset.getHist( dirDir+"/"+name )
    dirInt = hdir.Integral(0,-1)
    if not hdir.Integral(): return
    if binning: hdir = aux.rebin( hdir, binning )
    aux.appendFlowBin( hdir )
    hdir.SetYTitle( aux.getYAxisTitle( hdir ) )
    hdir.SetLineColor(1)
    hdir.SetMarkerColor(1)
    hdir.SetMarkerStyle(20)
    hdir.SetMarkerSize(0.8)
    hdir.drawOption_ = "pe x0"
    m.add( hdir, "#gamma  ({})".format(aux.metricPrefix(dirInt)) )

    hpre = controlDataset.getHist( preDir+"/"+name )
    preInt = hpre.Integral(0,-1)
    if not preInt: return
    if binning: hpre = aux.rebin( hpre, binning )
    aux.appendFlowBin( hpre )
    hpre.Scale( dirInt/preInt )
    hpre.SetYTitle( aux.getYAxisTitle( hpre ) )
    hpre.SetLineColor(ROOT.kRed)
    hpre.drawOption_ = "hist e"
    m.add( hpre, "jet ({})".format(aux.metricPrefix(preInt)) )

    if name=="met":
        hpreUp = controlDataset.getHist( preDir+"/"+name+"Up" )
        if binning: hpreUp = aux.rebin( hpreUp, binning )
        aux.appendFlowBin( hpreUp )
        hpreUp.Scale( dirInt/preInt )
        hpreUp.SetLineColor(ROOT.kRed)
        hpreUp.SetLineStyle(3)
        hpreUp.drawOption_ = "hist"
        m.add( hpreUp, "jet +" )
        hpreDn = controlDataset.getHist( preDir+"/"+name+"Dn" )
        if binning: hpreDn = aux.rebin( hpreDn, binning )
        aux.appendFlowBin( hpreDn )
        hpreDn.Scale( dirInt/preInt )
        hpreDn.SetLineColor(ROOT.kRed)
        hpreDn.SetLineStyle(2)
        hpreDn.drawOption_ = "hist"
        m.add( hpreDn, "jet -" )

    m.Draw()
    l = aux.Label(info=dataset.label, sim=dataset is not data)
    r = ratio.Ratio("#gamma/jet",hdir,hpre)
    r.draw(.5,1.5)
    if name=="met":
        hPlus = hpreUp.Clone(aux.randomName())
        hPlus.Divide(hpre)
        hPlus.Draw("same")
        hMinus = hpreDn.Clone(aux.randomName())
        hMinus.Divide(hpre)
        hMinus.Draw("same")


    if binningName: binningName = "_"+binningName
    saveName = "multiQcdClosure_{}_{}_{}{}".format(samplename, dirDir, name, binningName )
    aux.save( saveName )
    can.SetLogy()
    aux.save( saveName+"_log" )

    if name == "n_heJet" and preDir=="tr_jControl":
        aux.write2File( r.ratio.Clone(), "weight_n_heJet", "weights.root" )

def selectionComparison( dataset, controlDataset, name, samplename, binning, binningName, dirDir="tr", preDir="tr_jControl" ):
    if not controlDataset: controlDataset = dataset
    can = ROOT.TCanvas()
    m = multiplot.Multiplot()

    if "jControl" in dirDir and "jControl" in preDir: m.leg.SetHeader("Jet control sample")
    else: m.leg.SetHeader("Photon sample")

    hdir = dataset.getHist( dirDir+"/"+name )
    dirInt = hdir.Integral(0,-1)
    if not hdir.Integral(): return
    if binning: hdir = aux.rebin( hdir, binning )
    aux.appendFlowBin( hdir )
    hdir.SetYTitle( aux.getYAxisTitle( hdir ) )
    hdir.SetLineColor(1)
    hdir.SetMarkerColor(1)
    hdir.SetMarkerStyle(20)
    hdir.SetMarkerSize(0.8)
    hdir.drawOption_ = "pe x0"
    m.add( hdir, "E_{{T}}^{{miss}}>70 ({:.1f}k)".format(dirInt/1e3) )

    hpre = controlDataset.getHist( preDir+"/"+name )
    preInt = hpre.Integral(0,-1)
    if not preInt: return
    if binning: hpre = aux.rebin( hpre, binning )
    aux.appendFlowBin( hpre )
    hpre.Scale( dirInt/preInt )
    hpre.SetYTitle( aux.getYAxisTitle( hpre ) )
    hpre.SetLineColor(ROOT.kRed)
    hpre.drawOption_ = "hist e"
    m.add( hpre, "inclusive ({:.1f}k)".format(preInt/1e3) )

    m.Draw()
    l = aux.Label(info=dataset.label, sim=dataset is not data)
    r = ratio.Ratio("E_{T}^{miss}>70/all",hdir,hpre)
    r.draw(.5,1.5)

    if binningName: binningName = "_"+binningName
    saveName = "multiQcdClosureCompareSelection_{}_{}_{}{}".format(samplename, dirDir, name, binningName )
    aux.save( saveName )
    #can.SetLogy()
    #aux.save( saveName+"_log" )

def selectionComparison2( dataset, controlDataset, name, samplename, binning, binningName, dirDirs ):
    can = ROOT.TCanvas()
    m = multiplot.Multiplot()

    integral=None

    for dir in dirDirs:
        hdir = dataset.getHist( dir+"/"+name )
        if not integral: integral= hdir.Integral(0,-1)
        if binning: hdir = aux.rebin( hdir, binning )
        aux.appendFlowBin( hdir )
        hdir.SetYTitle( aux.getYAxisTitle( hdir ) )
        hdir.Scale(integral/hdir.Integral(0,-1))
        legend=""
        if "jControl" in dir:
            hdir.drawOption_ = "hist"
            legend+="jet"
        else:
            hdir.drawOption_ = "pe x0"
            hdir.SetMarkerStyle(20)
            hdir.SetMarkerSize(0.8)
            legend+="#gamma  "

        if "highMet" in dir:
            hdir.SetMarkerColor(1)
            hdir.SetLineColor(1)
            legend+=" high E_{T}^{miss}"
        else:
            hdir.SetLineColor(2)
            hdir.SetMarkerColor(2)
        if "lowMet" in dir:
            legend+=" low E_{T}^{miss}"
        m.add( hdir, legend )

    m.Draw()
    l = aux.Label(info=dataset.label, sim=dataset is not data)
    #r = ratio.Ratio("E_{T}^{miss}>70/all",hdir,hpre)
    #r.draw(.5,1.5)

    if binningName: binningName = "_"+binningName
    saveName = "multiQcdClosureCompareSelection2_{}_{}_{}{}".format(samplename, dirDirs[0], name, binningName )
    aux.save( saveName )
    can.SetLogy()
    aux.save( saveName+"_log" )




def multiQcdClosures( dataset, samplename, controlDataset=None ):
    names = aux.getObjectNames( dataset.files[0], "tr", [ROOT.TH1F] )

    for name in names:
        for binningName, binning in aux.getBinnigsFromName( name ).iteritems():
            multiQcdClosure( dataset, controlDataset, name, samplename, binning, binningName )
            multiQcdClosure( dataset, controlDataset, name, samplename, binning, binningName, dirDir="tr_he0", preDir="tr_jControl_he0" )
            multiQcdClosure( dataset, controlDataset, name, samplename, binning, binningName, dirDir="tr_he1", preDir="tr_jControl_he1" )
            multiQcdClosure( dataset, controlDataset, name, samplename, binning, binningName, dirDir="tr_he2", preDir="tr_jControl_he2" )
            multiQcdClosure( dataset, controlDataset, name, samplename, binning, binningName, dirDir="tr_he3", preDir="tr_jControl_he3" )
            multiQcdClosure( dataset, controlDataset, name, samplename, binning, binningName+"_wnjet", preDir="tr_jControl_wnjet" )
            multiQcdClosure( dataset, controlDataset, name, samplename, binning, binningName, dirDir="tr_highMet", preDir="tr_jControl_highMet" )
            selectionComparison( dataset, controlDataset, name, samplename, binning, binningName, dirDir="tr_highMet", preDir="tr" )
            selectionComparison( dataset, controlDataset, name, samplename, binning, binningName, dirDir="tr_jControl_highMet", preDir="tr_jControl" )
            selectionComparison2( dataset, controlDataset, name, samplename, binning, binningName, dirDirs=["tr_highMet", "tr_lowMet", "tr_jControl_highMet","tr_jControl_lowMet"] )




def metCorrections( dataset, samplename="" ):

    for binningName, binning in aux.getBinnigsFromName( "met" ).iteritems():
        for dir in "tr", "tr_jControl_wnjet_wht", "tr_jControl","tr_jControl_wnjet":
            can = ROOT.TCanvas()
            m = multiplot.Multiplot()
            for name in "met","metPar","metPer","metRaw","metParRaw","metPerRaw":
                hdir = dataset.getHist( dir+"/"+name )
                if not hdir.Integral(): return
                if binning: hdir = aux.rebin( hdir, binning )
                aux.appendFlowBin( hdir )
                hdir.SetYTitle( aux.getYAxisTitle( hdir ) )
                hdir.SetXTitle( "E_{T}^{miss} (GeV)" )
                hdir.drawOption_ = "pe x0"
                if "Raw" in name:
                    hdir.SetLineColor(ROOT.kRed)
                    hdir.SetMarkerColor(ROOT.kRed)
                else:
                    hdir.SetLineColor(1)
                    hdir.SetMarkerColor(1)
                if "Par" in name:
                    hdir.drawOption_ = "hist"
                    hdir.SetLineStyle(2)
                elif "Per" in name:
                    hdir.drawOption_ = "hist"
                    hdir.SetLineStyle(3)
                else:
                    hdir.SetMarkerStyle(20)
                legendName = name.replace("met","E_{T}^{miss}").replace("Par","#parallel").replace("Per"," {}^{#perp  }").replace("Raw","(uncorr)")
                m.add( hdir, legendName )
            m.Draw()
            can.SetLogy()
            aux.save("metComparison_{}_{}_{}_log".format(samplename,dir,binningName))



def qcdClosure( dataset, samplename="", gSet="tr_ht700", cSet="tr_jControl2" ):
    names = aux.getObjectNames( dataset.files[0], "", [ROOT.TH1F] )

    for name in names:
        if not name.endswith(cSet): continue

        hdir = getHistoFromDataset( dataset, name.replace( cSet, gSet ) )
        hdir.SetLineColor(1)
        hdir.SetMarkerColor(1)
        hdir.SetMarkerStyle(20)
        hdir.drawOption_ = "p"

        hpre = getHistoFromDataset( dataset, name )
        if cSet=="tr_jControl2" and dataset==data: hpre = getHistoFromDataset( dataHt, name ) # temporary fix
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

            m.add( mod_dir, "#gamma" )
            m.add( mod_pre, "#gamma-like" )

            if m.Draw():

                r = ratio.Ratio( "#gamma/#gamma-like", mod_dir, mod_pre )
                r.draw(0.5,1.5)
                if name  == "h_g_pt__tr_jControl" and binningName == "2":
                    aux.write2File( r.ratio.Clone(), name.replace("h_","weight_{}_".format(samplename)), "weights.root" )
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

        #scale = 0.018/2
        #scale = 0.006
        scale_e = scale/3
        scale_e = 0

        hpre.Scale( scale )
        hpre_sys = hpre.Clone( aux.randomName() )
        hpre_sys.SetMarkerSize(0)
        hpre_sys.SetFillColor(hpre_sys.GetLineColor())
        hpre_sys.SetFillStyle(3446)
        hpre_sys.drawOption_ = "e2"

        hpre.drawOption_ = "hist"
        for h in hdir, hpre:
            h.SetYTitle( aux.getYAxisTitle( h ) )

        for binningName, binning in aux.getBinnigsFromName( name ).iteritems():
            can = ROOT.TCanvas()
            m = multiplot.Multiplot()
            mod_dir = hdir
            mod_pre = hpre
            mod_pre_sys = hpre_sys

            if binning:
                mod_dir = aux.rebin(mod_dir, binning)
                mod_pre = aux.rebin(mod_pre, binning)
                mod_pre_sys = aux.rebin(mod_pre_sys, binning)
            for bin in range(mod_pre_sys.GetNbinsX()+2): mod_pre_sys.SetBinError( bin, mod_pre_sys.GetBinContent(bin)*scale_e/scale )

            m.add( mod_dir, "e#rightarrow#gamma" )
            m.add( mod_pre, "Normalized #gamma-like e" )
            #m.add( mod_pre, "{:.2f}% #times #gamma-like e".format(100*scale) )
            if scale_e: m.add( mod_pre_sys, "" )

            m.leg.SetX1(0.5)
            m.leg.SetY1(0.8)
            m.Draw()

            l = aux.Label(sim=True,info=dataset.label)
            if binningName: binningName = "_"+binningName
            aux.save( "ewkClosure_{}_{}{}".format(samplename,name,binningName ) )
            can.SetLogy()
            aux.save( "ewkClosure_{}_{}{}_log".format(samplename,name,binningName ) )


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


def efficiency( dataset, name, savename="" ):
    c = ROOT.TCanvas()
    c.SetLogy(0)
    h = getHistoFromDataset( dataset, name )
    if h.UsesWeights(): h.SetStatisticOption( ROOT.TEfficiency.kFNormal )

    h_pas = h.GetPassedHistogram()
    h_tot = h.GetTotalHistogram()
    if "hlt" in name:
        newBins = []
        if name == "eff_hlt_pt": newBins = range(0,80,8) + range(80,108,4) + range(108,300,12) + [300,400,500, 1000]
        if name == "eff_hlt_ht": newBins = range(0,1001,40) + [1000,1500,2000]
        if newBins:
            h_pas = aux.rebin( h_pas, newBins, False )
            h_tot = aux.rebin( h_tot, newBins, False )
        eff = ROOT.TEfficiency( h_pas, h_tot )
    else:
        eff = h

    eff.Draw()
    ROOT.gPad.Update()
    eff.GetPaintedGraph().GetYaxis().SetRangeUser(0., 1.1)

    if "_pt" in name:
        cutValue = 100
    elif "_ht" in name:
        cutValue = 600
    else:
        cutValue = 0

    if cutValue != None:
        bin = h_pas.FindFixBin( cutValue )
        passed = int(h_pas.Integral( bin, -1 ))
        total = int(h_tot.Integral( bin, -1 ))
        if not total: return
        conf = ROOT.TEfficiency().GetConfidenceLevel()
        e = 1.*passed/total
        e_up = ROOT.TEfficiency.ClopperPearson( total, passed, conf, True )
        e_dn = ROOT.TEfficiency.ClopperPearson( total, passed, conf, False )
        eLabel = ROOT.TLatex( 0.7, .15, "#varepsilon = {:.1f}^{{#plus{:.1f}}}_{{#minus{:.1f}}}%".format(100*e, 100*(e_up-e),100*(e-e_dn) ) )
        eLabel.SetNDC()
        eLabel.Draw()

        # graphical representation
        l = ROOT.TLine()
        l.SetLineWidth(2)
        l.SetLineColor( ROOT.kRed )
        xmax = eff.CreateGraph().GetHistogram().GetXaxis().GetXmax()
        l.DrawLine( cutValue, e, xmax, e )
        l.DrawLine( cutValue, e_up, xmax, e_up )
        l.DrawLine( cutValue, e_dn, xmax, e_dn )

        if cutValue > eff.CreateGraph().GetHistogram().GetXaxis().GetXmin():
            # cut line
            l.SetLineStyle(2)
            ymin = eff.GetPaintedGraph().GetYaxis().GetXmin()
            ymax = eff.GetPaintedGraph().GetYaxis().GetXmax()
            l.DrawLine( cutValue, ymin, cutValue, ymax )


    l = aux.Label(sim="Data" not in dataset.label)
    l.lum.DrawLatexNDC( .1, l.lum.GetY(), dataset.label )
    name = "_"+name.split("/")[-1]
    aux.save( "efficiency_"+savename+name )

    h_tot.SetLineColor(2)
    h_tot.Draw("hist")
    h_pas.Draw("same e*")
    aux.save( "efficiency_"+savename+name+"_raw" )

def efficiencies( dataset, savename="" ):
    names = ["triggerStudies/"+x for x in aux.getObjectNames( dataset.files[0], "triggerStudies", [ROOT.TEfficiency] ) ]

    for name in names:
        efficiency( dataset, name, savename )

def ewkIsrSamplesSplitting( dataset, isrDataset, saveName="test" ):
    names = aux.getObjectNames( dataset.files[0] )
    names = [ n for n in names if "_genPhoton" in n and n.startswith("h_") ]
    #names = [ "h_g_pt__tr_genPhoton" ]

    for name in names:
        for binningName, binning in aux.getBinnigsFromName( name ).iteritems():
            can = ROOT.TCanvas()
            m = multiplot.Multiplot()

            h_tot = dataset.getHist( name.replace("_genPhoton","" ) )
            h_g = dataset.getHist( name )
            h_e = dataset.getHist( name.replace("_genPhoton","_genElectron" ) )

            h_isr_tot = isrDataset.getHist( name.replace("_genPhoton","" ) )
            h_isr_g = isrDataset.getHist( name )
            h_isr_e = isrDataset.getHist( name.replace("_genPhoton","_genElectron" ) )

            if binning:
                h_tot =aux.rebin(h_tot, binning)
                h_g =aux.rebin(h_g, binning)
                h_e =aux.rebin(h_e, binning)

                h_isr_tot =aux.rebin(h_isr_tot, binning)
                h_isr_g =aux.rebin(h_isr_g, binning)
                h_isr_e =aux.rebin(h_isr_e, binning)


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
            m.add( h_g, "#gamma^{gen} match" )
            m.add( h_e, "e^{gen} match" )
            m.add( h_isr_tot, isrDataset.label )
            m.add( h_isr_g, "#gamma^{gen} match" )
            m.add( h_isr_e, "e^{gen} match" )

            if m.Draw():

                l = aux.Label()
                aux.save( "ewkIsrSampleSplitting_%s_%s_%s"%(saveName,name,binningName) )
                can.SetLogy()
                aux.save( "ewkIsrSampleSplitting_%s_%s_%s_log"%(saveName,name,binningName) )




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
    drawSameHistogram( "_gjets", "h_genHt", [gjets40,gjets100,gjets200,gjets400,gjets600] )
    drawSameHistogram( "_gjets_inc15", "h_genHt", [gjets_pt15,gjets40,gjets100], binning=range(0,120) )
    drawSameHistogram( "_qcd", "h_genHt", [qcd100,qcd200,qcd300,qcd500,qcd700,qcd1000,qcd1500,qcd2000] )
    drawSameHistogram( "_wjets", "h_genHt", [wjets100,wjets200,wjets400,wjets600,wjets800,wjets1200,wjets2500] )
    drawSameHistogram( "_wjets_trans1", "h_genHt", [wjets100,wjets200], binning=range(90,410,1) )
    drawSameHistogram( "_wjets_trans2", "h_genHt", [wjets200,wjets400], binning=range(190,610,1) )
    drawSameHistogram( "_wjets_trans3", "h_genHt", [wjets400,wjets600], binning=range(390,810,1) )
    drawSameHistogram( "_wjets_trans4", "h_genHt", [wjets600,wjets800], binning=range(590,1210,2) )
    drawSameHistogram( "_wjets_trans5", "h_genHt", [wjets800,wjets1200], binning=range(790,2510,2) )
    drawSameHistogram( "_wjets_trans6", "h_genHt", [wjets1200,wjets2500], binning=range(1190,3000,2) )
    drawSameHistogram( "_wjets_trans6", "h_genHt", [wjets1200,wjets2500], binning=range(1200,3000,20) )
    drawSameHistogram( "_wjets_trans45", "h_genHt", [wjets600,wjets800,wjets1200], binning=range(590,1510,5) )
    drawSameHistogram( "_wg_500", "h_g_pt__tr", [wg_pt500], [wg_mg], scaleToData=False )
    drawSameHistogram( "_zg", "h_g_pt__tr", [zg_130], [znunu], scaleToData=False )


def drawH2s():
    names = aux.getObjectNames( data.files[0], objects=[ROOT.TH2])

    for h2name in names:
        drawH2( data, h2name, "data" )
        drawH2( qcd+gjets, h2name, "gqcd" )
        #subtractH2( data, qcd+gjets, h2name, "data-gqcd" )

def drawISRsplitting():
    for d in ttjets,wjets,znunu:
        d.color=ROOT.kBlack
    for d in ttg,wg_mc,wg_mg,zg_130:
        d.color=ROOT.kRed

    ewkIsrSamplesSplitting( ttjets, ttg, "tt" )
    ewkIsrSamplesSplitting( wjets, wg_mc, "w_mc" )
    ewkIsrSamplesSplitting( wjets, wg_mg, "w_mg" )
    ewkIsrSamplesSplitting( znunu, zg_130, "zg" )

def checkGJetsQcdNlo():
    can = ROOT.TCanvas()
    m = multiplot.Multiplot()

    mc = gjets+qcd+ttjets+ttg+wjets+wg_mg+zg_130+znunu

    name = "h_tremht__tr"
    cutName = "tr_jControlLeadingJet"
    binningName = "1"

    h_mc_sr = mc.getHist( name )
    h_mc_cr = mc.getHist( name.replace("__tr","__"+cutName) )

    h_da_sr = data.getHist( name )
    h_da_cr = dataHt.getHist( name.replace("__tr","__"+cutName) )

    for h in h_mc_cr,h_da_cr:
        h.SetLineColor(ROOT.kRed)
        h.SetMarkerColor(ROOT.kRed)
    for h in h_mc_sr,h_da_sr:
        h.SetLineColor(ROOT.kBlack)
        h.SetMarkerColor(ROOT.kBlack)

    for h in h_mc_sr,h_mc_cr:
        h.SetLineStyle( 2 )
    for h in h_mc_cr,h_mc_sr,h_da_cr,h_da_sr:
        h.drawOption_="hist"
        h.Scale( 1./h.Integral() )
        h = aux.rebin( h, aux.getBinnigsFromName( name )["1"] )
    h_mc_cr = aux.rebin( h_mc_cr, aux.getBinnigsFromName( name )["1"] )
    h_mc_sr = aux.rebin( h_mc_sr, aux.getBinnigsFromName( name )["1"] )
    h_da_cr = aux.rebin( h_da_cr, aux.getBinnigsFromName( name )["1"] )
    h_da_sr = aux.rebin( h_da_sr, aux.getBinnigsFromName( name )["1"] )

    m.add( h_mc_sr, "MC #gamma" )
    m.add( h_mc_cr, "MC jet" )
    m.add( h_da_sr, "Data #gamma" )
    m.add( h_da_cr, "Data jet" )

    m.Draw()
    can.SetLogy()
    aux.save("checkGjetsQcdNlo")

def photonPosition( dataset, savename, dir="tr", normalize=False ):
    can = ROOT.TCanvas()
    m = multiplot.Multiplot()

    h2 = dataset.getHist( dir+"/n_heJets_vs_photonPosition" )

    colors = [ROOT.kBlack,ROOT.kBlue,ROOT.kCyan,ROOT.kGreen+4,ROOT.kRed,ROOT.kMagenta]+range(100)

    #for xbin in range(1,h2.GetNbinsX()+1):
    for xbin in range(1,6):
        h = h2.ProjectionY(str(xbin)+aux.randomName(),xbin,xbin)
        h.SetTitle(";photon position;Events")
        h.GetXaxis().SetTitleOffset(0.9)
        if not h.Integral(0,-1): continue
        if normalize: h.Scale(1./h.Integral())
        if normalize: h.GetYaxis().SetTitle("Normalized events")
        h.SetLineColor(colors[xbin-1])
        h.drawOption_="hist e"
        m.add( h, str(xbin-1) )

    m.leg.SetHeader("number of additional high-energetic jets")
    m.leg.SetX2(.99)
    m.leg.SetX1(.4)
    m.leg.SetY1(.5)
    m.Draw()
    aux.save("photonPosition_"+savename+"_"+dir )
    if False: aux.write2File( h2, savename, "gammaPosition.root" )


def main():
    pass
    #transitions()
    #compareAll( "_all", gjets400, gjets600, znn400, znn600 )
    #compareAll( "_GjetsVsZnn", gjets, znn )
    #compareAll( "_allMC", gjets, znn, qcd, wjets )
    #drawSameHistograms( "gqcd_data", [gjets600,gjets400,gjets200, gjets100,gjets40,qcd], [data])
    #drawSameHistograms( "gqcd_dataHt", [gjets600,gjets400,gjets200, gjets100,gjets40,qcd], [dataHt])
    #drawSameHistograms( "emqcd_data", [emqcd], [data])
    #drawSameHistograms( "_gjet15_data", [gjets_pt15, qcd], additional=[data])
    #drawSameHistograms( "mc_data", [gjets, qcd, ttjets, ttg, wjets,wg_mg,zg_130,znunu,dy], additional=[data])
    #drawSameHistograms( "mc_data", [gjets, qcd, ttjets, ttg, wjets,wg_mg,znunu,zg_130], additional=[data,signal["T5gg_1400_200"], signal["T5gg_1400_1200"]])
    #drawSameHistograms( "mc_data", [gjets, qcd, ttjets, ttg, wjets,wg_mg,znunu,zg_130], additional=[data])
    #drawSameHistograms( "mc_dataHt", [gjets, qcd, ttjets, ttg, wjets,wg_mg,znunu,zg_130], additional=[dataHt])
    #drawSameHistograms( "_mc_data", bkg=[gjets, qcd, ttjets, ttg, wjets, dy,znunu], additional=[data])
    #drawSameHistograms( "mc", [gjets, qcd, wjets, ttjets, ttg,znunu], additional=[signal["T5gg_1400_1200"], signal["T5gg_1400_200"]])
    #drawSameHistograms( "_QCD", bkg=[ qcd2000, qcd1500, qcd1000, qcd700, qcd500, qcd300 ] )
    #drawRazor( data )

    #ewkClosure( ttjets, "tt" )
    #ewkClosure( wjets, "w" )
    #ewkClosure( wjets+ttjets, "ewk" )
    #qcdClosure( qcd+gjets, "_gqcd" )
    #qcdClosure( data, "_data" )
    #multiQcdClosures( qcd+gjets, "gqcd" )
    #multiQcdClosures( zg_130+znunu, "znunu" )
    #multiQcdClosures( ttg+ttjets, "tt" )
    #multiQcdClosures( wg_mg+wjets, "w" )
    #multiQcdClosures( data, "data", dataHt )
    #multiQcdClosures( signal["T5Wg_1550_100"], "signal" )
    #multiQcdClosures( gjets+qcd+ttjets+ttg+wjets+wg_mg+znunu+zg_130,"mc" )
    #gjets.label = "GJets Data"
    #drawSameHistogram( "gjets_qcd", "h_genHt", [qcd], [gjets], scaleToData=True, binning=range(0,3000,20))
    #checkGJetsQcdNlo()

    #efficiencies( ttjets+qcd+gjets+wjets, "allMC_" )
    #efficiencies( qcd+gjets, "gqcd_" )
    #efficiencies( gjets, "gjet_" )
    #efficiencies( data, "singlePhoton" )
    #efficiencies( dataHt, "jetHt" )
    #efficienciesDataMC( dataHt, ttjets+qcd+gjets+wjets, "jetHt_mc_" )
    #efficiency( data_2015D, "eff_hlt_ht__offlinePT100__base", "singlePhotonD_" )
    #efficiency( data_prompt, "eff_hlt_ht__offlinePT100__base", "singlePhotonPrompt_" )
    #efficiency( data_prompt, "eff_hlt_ht__offlinePT100__extra__base", "singlePhotonPrompt_" )
    #efficiency( data, "eff_hlt_ht__base", "singlePhoton_" )
    #efficiency( dataHt_2015D, "eff_hlt_pt__offlineHT650__base", "dataHtD_" )
    #efficiency( dataHt_prompt, "eff_hlt_pt__offlineHT650__base", "dataHtPrompt_" )
    #efficiency( dataHt_prompt, "eff_hlt_pt__offlineHT650__extra__base", "dataHtPrompt_" )

    #drawROCs()
    #drawH2s()
    #drawISRsplitting()
    #razorStudies()

    #photonPosition( gjets+qcd, "gqcd" )
    #photonPosition( ttg+ttjets, "gjtt" )
    #photonPosition( gjets+qcd, "gqcd_jet", "tr_jControl" )
    #photonPosition( gjets+qcd, "gqcd_jet", "tr_jControl_wnjet" )
    #photonPosition( data, "data" )

    #metCorrections( gjets+qcd, "gqcd" )


if __name__ == "__main__":
    from datasets import *
    main()

