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
from rwthColors import rwth

import auxiliary as aux

intLumi = 2.26e3 # https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2577.html

def compare( datasets, name, saveName ):
    m = multiplot.Multiplot()

    for d in datasets:
        d.getHist(name)
        if not h.Integral(): continue
        h.Scale( 1./h.Integral() )
        m.add( h, d.label )

    m.Draw()

    aux.save( "compare%s_%s"%(saveName,name) )

def drawH2( dataset, name, savename="test" ):
    x = style.style2d()
    c = ROOT.TCanvas()
    h = dataset.getHist( name )
    h.Draw("colz")
    l = aux.Label(sim=savename!="data")
    aux.save( "h2_%s_%s"%(savename,name) )
    c.SetLogz()
    aux.save( "h2_%s_%s_log"%(savename,name) )
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

    can = ROOT.TCanvas()
    m = multiplot.Multiplot()

    scale = 1.
    if scaleToData: scale = divideDatasetIntegrals( [ i for i in additional if "Data" in i.label ], bkg, name )

    for d in bkg[-1::-1]:
        h = d.getHist( name )
        if not h: continue
        if not h.Integral(): continue
        h.Scale(scale)
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
            r.draw(0.5,1.5)

            if name  == "h_ht__tr" and binningName == "1":
                aux.writeWeight( r.ratio, name, sampleNames )

        info = ""
        l = aux.Label(info="#scale[0.7]{%s}"%info, sim=data not in additional)

        if binningName: binningName = "_"+binningName
        name = name.replace("/","__")
        saveName = "sameHistograms_{}_{}{}".format(sampleNames, name, binningName )
        aux.save( saveName )


def drawSameHistograms( sampleNames="test", stack=[], additional=[] ):
    file = stack[0].files[0] if stack else additional[0].files[0]
    names = aux.getObjectNames( file, "tr", [ROOT.TH1F] )
    additionalHt = [ x for x in additional if x is not data ]
    if data in additional: additionalHt += [ dataHt ]

    #names = ["met"]

    for name in names:
        for binningName, binning in aux.getBinnigsFromName( name ).iteritems():
            drawSameHistogram( sampleNames, "tr/"+name, stack, additional, binning, binningName )
            drawSameHistogram( sampleNames, "tr_highMet/"+name, stack, additional, binning, binningName )
            drawSameHistogram( sampleNames, "tr_mediumMet/"+name, stack, additional, binning, binningName )
            drawSameHistogram( sampleNames, "tr_highHt/"+name, stack, additional, binning, binningName )
            drawSameHistogram( sampleNames, "tr_jControl_highHt/"+name, stack, additional, binning, binningName )
            drawSameHistogram( sampleNames, "tr_jControl/"+name, stack, additionalHt, binning, binningName )

def multiQcdClosure( dataset, controlDataset, name, samplename, binning, binningName, dirDir="tr", preDir="tr_jControl" ):
    if not controlDataset: controlDataset = dataset
    can = ROOT.TCanvas()
    m = multiplot.Multiplot()


    hdir = dataset.getHist( dirDir+"/"+name )
    dirInt,dirIntErr = aux.integralAndError(hdir)
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
    preInt,preIntErr = aux.integralAndError(hpre)
    if not preInt: return
    if binning: hpre = aux.rebin( hpre, binning )
    aux.appendFlowBin( hpre )
    hpre.Scale( dirInt/preInt )
    scaleErr = dirInt/preInt * math.sqrt( (dirIntErr/dirInt)**2 + (preIntErr/preInt)**2 )
    m.leg.SetHeader("scale = ({:3.1f}#pm{:3.1f})m".format(dirInt/preInt*1000,1000*scaleErr))
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




def multiQcdClosures( dataset, samplename, controlDataset=None ):
    names = aux.getObjectNames( dataset.files[0], "tr", [ROOT.TH1F] )

    names = ["met"]

    for name in names:
        for binningName, binning in aux.getBinnigsFromName( name ).iteritems():
            multiQcdClosure( dataset, controlDataset, name, samplename, binning, binningName )
            multiQcdClosure( dataset, controlDataset, name, samplename, binning, binningName, dirDir="tr_he0", preDir="tr_jControl_he0" )
            multiQcdClosure( dataset, controlDataset, name, samplename, binning, binningName, dirDir="tr_he1", preDir="tr_jControl_he1" )
            multiQcdClosure( dataset, controlDataset, name, samplename, binning, binningName, dirDir="tr_he2", preDir="tr_jControl_he2" )
            multiQcdClosure( dataset, controlDataset, name, samplename, binning, binningName, dirDir="tr_he3", preDir="tr_jControl_he3" )
            multiQcdClosure( dataset, controlDataset, name, samplename, binning, binningName+"_wnjet", preDir="tr_jControl_wnjet" )
            multiQcdClosure( dataset, controlDataset, name, samplename, binning, binningName, dirDir="tr_highMet", preDir="tr_jControl_highMet" )
            multiQcdClosure( dataset, controlDataset, name, samplename, binning, binningName, dirDir="tr_mediumMet", preDir="tr_jControl_mediumMet" )
            selectionComparison( dataset, controlDataset, name, samplename, binning, binningName, dirDir="tr_highMet", preDir="tr" )
            selectionComparison( dataset, controlDataset, name, samplename, binning, binningName, dirDir="tr_jControl_highMet", preDir="tr_jControl" )
            selectionComparison2( dataset, controlDataset, name, samplename, binning, binningName, dirDirs=["tr_highMet", "tr_lowMet", "tr_jControl_highMet","tr_jControl_lowMet"] )




def metCorrections( dataset, samplename="" ):

    for binningName, binning in aux.getBinnigsFromName( "met" ).iteritems():
        for dir in "tr", "tr_jControl_wnjet", "tr_jControl","tr_jControl_wnjet":
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
                    hdir.SetLineColor(rwth.red)
                    hdir.SetMarkerColor(rwth.red)
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
            aux.save("metComparison_{}_{}_{}".format(samplename,dir,binningName))




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
    h = dataset.getHist( name )
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
    aux.save( "efficiency_"+savename+name, log=False )

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

def drawISRsplitting():
    for d in ttjets,wjets,znunu:
        d.color=ROOT.kBlack
    for d in ttg,wg_mc,wg_mg,zg_130:
        d.color=ROOT.kRed

    ewkIsrSamplesSplitting( ttjets, ttg, "tt" )
    ewkIsrSamplesSplitting( wjets, wg_mc, "w_mc" )
    ewkIsrSamplesSplitting( wjets, wg_mg, "w_mg" )
    ewkIsrSamplesSplitting( znunu, zg_130, "zg" )


def photonPosition( dataset, savename, dir="tr", normalize=False ):
    can = ROOT.TCanvas()
    m = multiplot.Multiplot()

    h2 = dataset.getHist( dir+"/n_heJets_vs_photonPosition" )

    colors = [ rwth.black, rwth.red, rwth.blue, rwth.lila, rwth.petrol50, rwth.green ]+range(100)

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
    #drawSameHistograms( "mc_data", [gjets, qcd, ttjets, ttg, wjets,wg_mg,zg_130,znunu], additional=[signal["T5Wg_1550_100"],data])
    #drawSameHistograms( "mc_dataHt", [gjets, qcd, ttjets, ttg, wjets,wg_mg,znunu,zg_130], additional=[dataHt])
    #drawSameHistograms( "mc", [gjets, qcd, ttjets, ttg, wjets,wg_mg,zg_130,znunu], additional=[signal["T5Wg_1550_100"],signal["T5Wg_1550_1500"]])

    #ewkClosure( ttjets, "tt" )
    #ewkClosure( wjets, "w" )
    #ewkClosure( wjets+ttjets, "ewk" )

    #multiQcdClosures( gjets+qcd, "gqcd" )
    #multiQcdClosures( zg_130, "znunu",znunu )
    #multiQcdClosures( ttg, "tt", ttjets )
    #multiQcdClosures( wg_mg, "w", wjets )
    #multiQcdClosures( data, "data", dataHt )
    #multiQcdClosures( signal["T5Wg_1550_100"], "signal" )
    #multiQcdClosures( gjets+qcd+ttjets+ttg+wjets+wg_mg+znunu+zg_130,"mc" )

    #efficiencies( ttjets+qcd+gjets+wjets, "allMC_" )
    #efficiencies( qcd+gjets, "gqcd_" )
    #efficiencies( gjets, "gjet_" )
    #efficiencies( data, "singlePhoton" )
    #efficiencies( dataHt, "jetHt" )

    #drawROCs()
    #drawH2s()
    #drawISRsplitting()

    #photonPosition( gjets+qcd, "gqcd" )
    #photonPosition( ttg+ttjets, "gjtt" )
    #photonPosition( gjets+qcd, "gqcd_jet", "tr_jControl" )
    #photonPosition( gjets+qcd, "gqcd_jet", "tr_jControl_wnjet" )
    #photonPosition( data, "data" )

    #metCorrections( gjets+qcd, "gqcd" )


if __name__ == "__main__":
    from datasets import *
    #ROOT.gErrorIgnoreLevel = ROOT.kFatal
    main()

