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
intLumi = 2.32e3 # https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2611.html
# and here https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2015Analysis#ReReco_at_25_ns

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

            if sampleNames == "mc_data" and name  == "tr_jControl/emht" and binningName == "1":
                aux.write2File( r.ratio, "mc_data_tr_jControl__emht_1", "weights.root" )

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
    if data in additional:
        if "genMatch" in names: names.remove("genMatch")
        if "genHt" in names: names.remove("genHt")

    for name in names:
        for binningName, binning in aux.getBinnigsFromName( name ).iteritems():
            drawSameHistogram( sampleNames, "tr/"+name, stack, additional, binning, binningName )
            drawSameHistogram( sampleNames, "tr_eControl/"+name, stack, additional, binning, binningName )
            drawSameHistogram( sampleNames, "tr_0met100/"+name, stack, additional, binning, binningName )
            drawSameHistogram( sampleNames, "tr_jControl/"+name, stack, additionalHt, binning, binningName )
            drawSameHistogram( sampleNames, "tr_jControl_prescale/"+name, stack, additionalHt, binning, binningName )
            drawSameHistogram( sampleNames, "tr_ee/"+name, stack, additionalHt, binning, binningName )

def multiQcdClosure( dataset, controlDataset, name, samplename, binning, binningName, dirDir="tr", preDir="tr_jControl" ):
    if not controlDataset: controlDataset = dataset
    can = ROOT.TCanvas()
    m = multiplot.Multiplot()


    hdir = dataset.getHist( dirDir+"/"+name )
    dirInt,dirIntErr = aux.integralAndError(hdir)
    if not hdir.Integral(): return
    if binning: hdir = aux.rebin( hdir, binning )
    aux.appendFlowBin( hdir )
    if style.divideByBinWidth: hdir.Scale(1,"width")
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
    if style.divideByBinWidth: hpre.Scale(1,"width")
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
    saveName = "multiQcdClosure_{}_{}_vs_{}_{}{}".format(samplename, dirDir, preDir, name, binningName )
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
    saveName = "multiQcdClosureCompareSelection_{}_{}_{}{}".format(samplename, dirDir+preDir, name, binningName )
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

    #names = ["met","emht"]

    for name in names:
        for binningName, binning in aux.getBinnigsFromName( name ).iteritems():
            multiQcdClosure( dataset, controlDataset, name, samplename, binning, binningName )
            multiQcdClosure( dataset, controlDataset, name, samplename, binning, binningName, dirDir="tr", preDir="tr_ee" )
            multiQcdClosure( dataset, controlDataset, name, samplename, binning, binningName, dirDir="tr_jControl", preDir="tr_jControl_ee" )
            multiQcdClosure( dataset, controlDataset, name, samplename, binning, binningName, dirDir="tr_0met100", preDir="tr_jControl_0met100" )
            multiQcdClosure( dataset, controlDataset, name, samplename, binning, binningName, dirDir="tr_100met", preDir="tr_jControl_100met" )
            multiQcdClosure( dataset, controlDataset, name, samplename, binning, binningName, dirDir="tr_noGenLep")
            multiQcdClosure( dataset, controlDataset, name, samplename, binning, binningName, dirDir="tr_he0", preDir="tr_jControl_he0" )
            multiQcdClosure( dataset, controlDataset, name, samplename, binning, binningName, dirDir="tr_he1", preDir="tr_jControl_he1" )
            multiQcdClosure( dataset, controlDataset, name, samplename, binning, binningName, dirDir="tr_he2", preDir="tr_jControl_he2" )
            multiQcdClosure( dataset, controlDataset, name, samplename, binning, binningName, dirDir="tr_he3", preDir="tr_jControl_he3" )
            multiQcdClosure( dataset, controlDataset, name, samplename, binning, binningName+"_wnjet", preDir="tr_jControl_wnjet" )



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


def ewkClosure( dataset, controlDataset, name, samplename, binning, binningName, dirDir="tr_genE", preDir="tr_eControl" ):
    if not controlDataset: controlDataset = dataset

    f, ef = 1.07/100., 30/100.

    can = ROOT.TCanvas()
    m = multiplot.Multiplot()
    hdir = dataset.getHist( dirDir+"/"+name )
    if not hdir.Integral(): return
    if binning: hdir = aux.rebin( hdir, binning )
    aux.appendFlowBin( hdir )
    hdir.SetYTitle( aux.getYAxisTitle( hdir ) )
    hdir.SetLineColor(1)
    hdir.SetMarkerColor(1)
    hdir.SetMarkerStyle(20)
    hdir.SetMarkerSize(0.8)
    hdir.drawOption_ = "pe x0"
    m.add( hdir, "e#rightarrow#gamma" )

    hpre = controlDataset.getHist( preDir+"/"+name )
    if not hpre.Integral(): return
    if binning: hpre = aux.rebin( hpre, binning )
    aux.appendFlowBin( hpre )
    hpre.Scale(0.01)
    hpre.SetYTitle( aux.getYAxisTitle( hpre ) )
    hpre.SetLineColor(ROOT.kRed)
    hpre.drawOption_ = "hist e"
    m.add( hpre, "" )

    hpre_sys = hpre.Clone( aux.randomName() )
    hpre_sys.SetMarkerSize(0)
    hpre_sys.SetFillColor(hpre_sys.GetLineColor())
    hpre_sys.SetFillStyle(3446)
    hpre_sys.drawOption_ = "e2"

    for bin in range(hpre_sys.GetNbinsX()+2): hpre_sys.SetBinError( bin, hpre_sys.GetBinContent(bin)*ef )
    m.add( hpre_sys, "Prediction" )

    m.leg.SetX1(0.5)
    m.leg.SetY1(0.8)
    m.Draw()

    r = ratio.Ratio("e#rightarrow#gamma/Pred", hdir, hpre, hpre_sys)
    r.draw(.5,1.5)

    l = aux.Label(sim=True,info=dataset.label)
    if binningName: binningName = "_"+binningName
    aux.save( "ewkClosure_{}_{}{}".format(samplename,name,binningName ) )

def ewkClosures( dataset, samplename="", controlDataset=None ):
    names = aux.getObjectNames( dataset.files[0], "tr_eControl", [ROOT.TH1F] )

    #names = ["met"]

    for name in names:
        for binningName, binning in aux.getBinnigsFromName( name ).iteritems():
            ewkClosure( dataset, controlDataset, name, samplename, binning, binningName )


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


def efficiency( dataset, name, savename="", binning=None, binningName="" ):
    c = ROOT.TCanvas()
    c.SetLogy(0)
    eff = dataset.getHist( name )
    if eff.UsesWeights(): eff.SetStatisticOption( ROOT.TEfficiency.kFNormal )

    h_pas = eff.GetPassedHistogram()
    h_tot = eff.GetTotalHistogram()

    if name.endswith("_preScaled"):
        eff2 = dataset.getHist( name.replace("_preScaled", "") )
        h_tot = eff2.GetTotalHistogram()

    if binning:
        h_pas = aux.rebin( h_pas, binning, False )
        h_tot = aux.rebin( h_tot, binning, False )

    if name.endswith("_preScaled"):
        ratio = h_pas.Clone(aux.randomName())
        ratio.Divide(h_tot)
        gr = ROOT.TGraphAsymmErrors(ratio)
        gr.SetLineColor(1)
        gr.Draw("ap")
    else:
        x = h_pas.Clone(aux.randomName())
        y = h_tot.Clone(aux.randomName())
        eff = ROOT.TEfficiency(h_pas, h_tot)
        h_pas, h_tot = x, y
        eff.Draw()
        ROOT.gPad.Update()
        gr = eff.GetPaintedGraph()

    gr.GetYaxis().SetRangeUser(0., 1.1)

    if name.endswith("eff_hlt_pt"):
        cutValue = 100
    elif name.endswith("eff_hlt_ht") \
            or name.endswith("eff_hlt_ht_ct") \
            or name.endswith("eff_hlt_ht_ct_preScaled"):
        cutValue = 700
    else:
        cutValue = 0

    if cutValue or True:
        bin = h_pas.FindFixBin( cutValue )
        passed = int(h_pas.Integral( bin, -1 ))
        total = int(h_tot.Integral( bin, -1 ))
        if not total: return
        e = 1.*passed/total
        if passed<=total:
            conf = ROOT.TEfficiency().GetConfidenceLevel()
            e_up = ROOT.TEfficiency.ClopperPearson( total, passed, conf, True )
            e_dn = ROOT.TEfficiency.ClopperPearson( total, passed, conf, False )
        else:
            passed, epassed = aux.integralAndError(h_pas,bin,-1)
            total, etotal = aux.integralAndError(h_pas,bin,-1)
            ee = e * math.sqrt( (epassed/passed)**2 + (etotal/total)**2 )
            e_up = e + ee
            e_dn = e - ee
        eLabel = ROOT.TLatex( 0.7, .15, "#varepsilon = {:.1f}^{{#plus{:.1f}}}_{{#minus{:.1f}}}%".format(100*e, 100*(e_up-e),100*(e-e_dn) ) )
        eLabel.SetNDC()
        eLabel.Draw()

        # graphical representation
        l = ROOT.TLine()
        l.SetLineWidth(2)
        l.SetLineColor( ROOT.kRed )
        xmax = gr.GetHistogram().GetXaxis().GetXmax()
        l.DrawLine( cutValue, e, xmax, e )
        l.DrawLine( cutValue, e_up, xmax, e_up )
        l.DrawLine( cutValue, e_dn, xmax, e_dn )

        if cutValue > eff.CreateGraph().GetHistogram().GetXaxis().GetXmin():
            # cut line
            l.SetLineStyle(2)
            ymin = gr.GetYaxis().GetXmin()
            ymax = gr.GetYaxis().GetXmax()
            l.DrawLine( cutValue, ymin, cutValue, ymax )


    l = aux.Label(sim="Data" not in dataset.label)
    l.lum.DrawLatexNDC( .1, l.lum.GetY(), dataset.label )
    name = "_"+name.split("/")[-1]
    if binningName: binningName = "_"+binningName
    aux.save( "efficiency_"+savename+name+binningName, log=False )

    if False:
        h_tot.SetLineColor(2)
        h_tot.Draw("hist")
        h_pas.Draw("same e*")
        aux.save( "efficiency_"+savename+name+"_raw" )

def efficiencies( dataset, savename="" ):
    names = ["triggerStudies/"+x for x in aux.getObjectNames( dataset.files[0], "triggerStudies", [ROOT.TEfficiency] ) ]

    #names = ["triggerStudies/eff_hlt_ht_ct2_preScaled"]

    for name in names:
        efficiency( dataset, name, savename )
        if name.endswith("eff_hlt_pt"):
            efficiency( dataset, name, savename, range(0,80,8) + range(80,108,4) + range(108,300,12) + [300,400,500, 1000], "1" )
        if name.endswith("eff_hlt_ht"):
            efficiency( dataset, name, savename, range(0,1001,40) + range(1000,1500,2000), "1" )
        if name.endswith("eff_hlt_ht_ct") or name.endswith("eff_hlt_ht_ct_preScaled"):
            efficiency( dataset, name, savename, range(450,1001,10), "1" )
        if name.endswith("eff_hlt_ht_ct2") or name.endswith("eff_hlt_ht_ct2_preScaled"):
            efficiency( dataset, name, savename, range(500,1001,50), "1" )
        if "_met_" in name:
            efficiency( dataset, name, savename, range(0, 100, 5)+range(100,151,10), "1" )



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




def significanceMetHt():
    x = style.style2d()
    allDatasets = qcd+ttjets+ttg+wg_mg+znunu+zg_130+gjets
    sigDatasets = signal["T5Wg_1550_100"]

    bkgH = allDatasets.getHist("tr/met_vs_emht")
    sigH = sigDatasets.getHist("tr/met_vs_emht")
    for h in bkgH,sigH:
        h.Rebin2D(2,2)

    significanceH = sigH.Clone(aux.randomName())
    significanceH.SetTitle(";minimal E_{T}^{miss} (GeV);minimal H_{T} (GeV);Significance")
    significanceH.GetXaxis().SetRangeUser(0,1000)
    significanceH.GetZaxis().SetRangeUser(0,4)

    for xbin, ybin in aux.loopH2(sigH):
        s = sigH.Integral(xbin,-1,ybin,-1)
        be = ROOT.Double()
        b = bkgH.IntegralAndError(xbin,-1,ybin,-1,be)
        if not s or b <= 0: continue
        sig = s/math.sqrt(b+(b/2)**2+be**2)
        significanceH.SetBinContent(xbin,ybin,sig)

    significanceH.Draw("colz")
    aux.save("significanceMax2D_{}".format("test"), log=False)

def zToMet():
    names = aux.getObjectNames( zgll.files[0], "tr", [ROOT.TH1F] )

    names = ["met", "g_pt"]

    bfNNToLL = 1.980707905005249 # branching fraction Z(νν)/Z(ll)

    for name in names:
        binningSettings = aux.getBinnigsFromName( name )
        if name == "met": binningSettings["10"] = range(0,1000,20)
        if name == "met": binningSettings["11"] = range(0,1000,100)
        if name == "g_pt": binningSettings["10"] = range(0,1000,10)
        for binningName, binning in binningSettings.iteritems():

            c = ROOT.TCanvas()
            m = multiplot.Multiplot()
            h130 = aux.stdHist(zg_130, "tr/"+name, binning)
            m.addStack(h130, ">130")

            h0To130 = aux.stdHist(zgll, "tr_0pt130/"+name, binning)
            h0To130.Scale(bfNNToLL)
            m.addStack(h0To130, "<130 from Z#rightarrowll")

            m.Draw()
            aux.save("zToMet_{}_{}".format(name,binningName))

            c2 = ROOT.TCanvas()
            m2 = multiplot.Multiplot()
            h130 = aux.stdHist(zg_130, "tr/"+name, binning)
            h130.drawOption_ = "hist e"
            m.addStack(h130, ">130")
            m2.add(h130, ">130")

            h130pt = aux.stdHist(zgll, "tr_130pt/"+name, binning )
            h130pt.Scale(bfNNToLL)
            m2.add(h130pt, ">130 from Z#rightarrowll")

            m2.Draw()
            aux.save("zToMetCompare_{}_{}".format(name,binningName))
            del c

def htStuff():
    allDatasets = gjets, qcd, ttjets, ttg, wjets, wg_mg, znunu, zg_130
    signals = signal["T5Wg_1550_100"], signal["T5Wg_1550_1500"]

    names = ["tr/met_vs_emht","tr_jControl/met_vs_emht"]

    metBinnings = aux.getBinnigsFromName("met")
    emhtBinnings = aux.getBinnigsFromName("emht")

    for name in names:
        for binningName, binning in metBinnings.iteritems():
            for cut1, cut2 in [(0,1e6), (0,2000), (2000,1e6)]:
                c = ROOT.TCanvas()
                m = multiplot.Multiplot()
                for d in allDatasets:
                    h = aux.stdHist(d, name, binning=binning, xCut=False, cut1=cut1, cut2=cut2 )
                    h.SetXTitle("E_{T}^{miss} (GeV)")
                    m.addStack(h,d.label)
                for s in signals:
                    h = aux.stdHist(s, name, binning=binning, xCut=False, cut1=cut1, cut2=cut2 )
                    h.drawOption_="hist"
                    h.SetLineWidth(3)
                    m.add(h,s.label)
                m.sortStackByIntegral()
                m.Draw()
                info = "EMH_{T}/GeV"
                if cut1: info = str(cut1)+"<"+info
                if cut2<1e5: info = info+"<"+str(cut2)
                if "<" not in info and ">" not in info: info=""
                l = aux.Label(sim=True,info=info)
                aux.save("htStuff_met_{}_{}_{}ht{}".format(name.split("/")[0], binningName, int(cut1), int(cut2)))

        for binningName, binning in emhtBinnings.iteritems():
            for cut1, cut2 in [(0,1e6), (0,200), (200,1e6)]:
                c = ROOT.TCanvas()
                m = multiplot.Multiplot()
                for d in allDatasets:
                    h = aux.stdHist(d, name, binning=binning, xCut=True, cut1=cut1, cut2=cut2 )
                    h.SetXTitle("EMH_{T} (GeV)")
                    m.addStack(h,d.label)
                for s in signals:
                    h = aux.stdHist(s, name, binning=binning, xCut=True, cut1=cut1, cut2=cut2 )
                    h.drawOption_="hist"
                    h.SetLineWidth(3)
                    m.add(h,s.label)
                m.sortStackByIntegral()
                m.Draw()
                info = "E_{T}^{miss}/GeV"
                if cut1: info = str(cut1)+"<"+info
                if cut2<1e5: info = info+"<"+str(cut2)
                if "<" not in info and ">" not in info: info=""
                l = aux.Label(sim=True,info=info)
                aux.save("htStuff_emht_{}_{}_{}met{}".format(name.split("/")[0], binningName, int(cut1), int(cut2)))

def qcdPrediction3d(totalNum, totalDen, xCut=100, save=False):

    xCutBin = totalNum.GetXaxis().FindBin(xCut) - 1
    for h in totalNum, totalDen:
        h.GetXaxis().SetRange(0, xCutBin)

    num = totalNum.Project3D("zye")
    den = totalDen.Project3D("zye2") # different name
    num.Divide(den)

    if save or True:
        aux.write2File( num, "weight_emht_vs_pt_gqcd", "weights.root" )

    pre = totalDen.Clone(aux.randomName())
    preSys = totalDen.Clone(aux.randomName())
    for xbin, ybin, zbin in aux.loopH3(pre):
        w = num.GetBinContent(ybin, zbin)
        we = num.GetBinError(ybin, zbin)
        c = totalDen.GetBinContent(xbin, ybin, zbin)
        ce = totalDen.GetBinError(xbin, ybin, zbin)
        pre.SetBinContent(xbin, ybin, zbin, c*w)
        pre.SetBinError(xbin, ybin, zbin, ce*w )
        preSys.SetBinContent(xbin, ybin, zbin, c*w )
        preSys.SetBinError(xbin, ybin, zbin, c*we )
    return pre, preSys

def qcdClosure3d(dirSet, name, dirName="tr", preSet=None):
    if not preSet: preSet = dirSet

    metBinning = aux.getBinnigsFromName("met")["3"]
    emhtBinning = aux.getBinnigsFromName("emht")["2"]
    ptBinning = aux.getBinnigsFromName("g_pt")["1"]
    ptBinning = None

    hName = "metRaw_vs_emht_vs_njet"
    h3dir = dirSet.getHist(dirName+"/"+hName)
    h3pre = preSet.getHist("tr_jControl/"+hName)

    h3dir = aux.rebin3d(h3dir, metBinning, emhtBinning, ptBinning)
    h3pre = aux.rebin3d(h3pre, metBinning, emhtBinning, ptBinning)

    h3w, h3wSys = qcdPrediction3d(h3dir, h3pre, 100)

    # compute total weight for simple total scaling
    metCut = 100
    metCutBin = h3dir.GetXaxis().FindBin(metCut) - 1
    cuts = [
        ("xyz", 0, -1, 0, -1),
        ("xyz", 0, 2000, 0, -1),
        ("xyz", 2000, -1, 0, -1),
        ("xyz", 0, 800, 0, -1),
        ("yxz", 0, 100, 0, -1),
        ("yxz", 100, -1, 0, -1),
        ("zxy", 0, -1, 0, -1),
    ]
    for dir, cut1low, cut1high, cut2low, cut2high in cuts:
        cut1BinLow = aux.getAxis(h3dir, dir[1]).FindBin(cut1low)
        cut1BinHigh = aux.getAxis(h3dir, dir[1]).FindBin(cut1high-1e-6) if cut1high > 0 else cut1high
        cut2BinLow = aux.getAxis(h3dir, dir[2]).FindBin(cut2low)
        cut2BinHigh = aux.getAxis(h3dir, dir[2]).FindBin(cut2high-1e-6) if cut2high > 0 else cut2high
        aux.getAxis(h3dir, dir[0]).SetRange()
        aux.getAxis(h3dir, dir[1]).SetRange(cut1BinLow, cut1BinHigh)
        aux.getAxis(h3dir, dir[2]).SetRange(cut2BinLow, cut2BinHigh)
        aux.getAxis(h3w, dir[0]).SetRange()
        aux.getAxis(h3w, dir[1]).SetRange(cut1BinLow, cut1BinHigh)
        aux.getAxis(h3w, dir[2]).SetRange(cut2BinLow, cut2BinHigh)
        aux.getAxis(h3wSys, dir[0]).SetRange()
        aux.getAxis(h3wSys, dir[1]).SetRange(cut1BinLow, cut1BinHigh)
        aux.getAxis(h3wSys, dir[2]).SetRange(cut2BinLow, cut2BinHigh)
        h1dir = h3dir.Project3D(dir[0]+"e_dir")
        h1w = h3w.Project3D(dir[0]+"e_w")
        h1wSys = h3wSys.Project3D(dir[0]+"e_we")

        for h in h1dir, h1w, h1wSys:
            h.SetTitle("")
            aux.appendFlowBin(h)

        aux.drawOpt(h1dir, "data")
        h1w.drawOption_ = "hist e"
        for h in h1w, h1wSys:
            h.SetLineColor(ROOT.kRed)
            if dir == "y": h.SetTitleOffset(1)

        aux.drawOpt(h1wSys, "sys")

        c = ROOT.TCanvas()
        m = multiplot.Multiplot()
        m.add(h1dir, "#gamma")
        m.add(h1w, "Jet EMH_{T} weighted")
        m.add(h1wSys, "Jet weight uncert.")

        info = ""
        if cut1low: info += "{}<".format(cut1low)
        if cut1low or cut1high>0: info += aux.getAxis(h3dir,dir[1]).GetTitle().replace(" (GeV)","")
        if cut1high>0: info += "<{}".format(cut1high)
        if cut2low or cut2high>0: info += " "
        if cut2low: info += "{}<".format(cut2low)
        if cut2low or cut2high>0: info += aux.getAxis(h3dir,dir[2]).GetTitle().replace(" (GeV)","")
        if cut2high>0: info += "<{}".format(cut2high)
        m.leg.SetHeader(info)
        m.Draw()

        l = aux.Label(sim=True, info=dirSet.label)

        r = ratio.Ratio("#gamma/jet", h1dir, h1w, h1wSys)
        r.draw(0.5,1.5)

        cutName = dir[0] + "_"
        if cut1low: cutName += str(cut1low)
        if cut1low or cut1high>0: cutName += dir[1]
        if cut1high>0: cutName += str(cut1high)
        if cut2low: cutName += str(cut2low)
        if cut2low or cut2high>0: cutName += dir[2]
        if cut2high>0: cutName += str(cut2high)

        aux.save("qcdClosure_3d_{}_{}".format(name, cutName))

def qcdPrediction2dsimple(h2num, h2den, xCut=100, yCut=0):
    xCutBin = h2num.GetXaxis().FindBin(xCut) - 1
    h2denW = h2den.Clone(aux.randomName())
    h2denWsys = h2den.Clone(aux.randomName())
    ybin = h2num.GetYaxis().FindFixBin(yCut)
    h2denW.Scale(h2num.Integral(0,-1,ybin,-1)/h2den.Integral(0,-1,ybin,-1))
    h2denWsys.Scale(h2num.Integral(0,-1,ybin,-1)/h2den.Integral(0,-1,ybin,-1))
    return h2denW, h2denWsys

def qcdPrediction2d(h2num, h2den, xCut=100, save=False):
    xCutBin = h2num.GetXaxis().FindBin(xCut) - 1
    h1numY = h2num.ProjectionY(aux.randomName(), 0, xCutBin)
    h1denY = h2den.ProjectionY(aux.randomName(), 0, xCutBin)
    h1numY.Divide(h1denY)

    if save:
        aux.write2File( h1numY, "weight_emht_gqcd", "weights.root" )

    h2denW = h2den.Clone(aux.randomName())
    h2denWsys = h2den.Clone(aux.randomName())
    for xbin, ybin in aux.loopH2(h2den):
        w = h1numY.GetBinContent(ybin)
        we = h1numY.GetBinError(ybin)
        c = h2den.GetBinContent(xbin,ybin)
        h2denW.SetBinContent(xbin, ybin, c*w )
        h2denW.SetBinError(xbin, ybin, h2den.GetBinError(xbin, ybin) * w )
        h2denWsys.SetBinContent(xbin, ybin, c*w )
        h2denWsys.SetBinError(xbin, ybin, c*we )
    return h2denW, h2denWsys

def htRebinning(dSets, name, dirName="tr", predSets=None, simpleScale=False):
    if not predSets: predSets = dSets

    metBinning = aux.getBinnigsFromName("met")["3"]
    emhtBinning = aux.getBinnigsFromName("emht")["2"]

    hName = "metRaw_vs_emht"
    if "Raw" not in hName: name+="_typeI"
    h2GJet = dSets.getHist(dirName+"/"+hName)
    h2Qcd  = predSets.getHist("tr_jControl/"+hName)

    h2GJet = aux.rebin2d( h2GJet, metBinning, emhtBinning )
    h2Qcd = aux.rebin2d( h2Qcd, metBinning, emhtBinning )

    if simpleScale:
        h2QcdW, h2QcdWsys = qcdPrediction2dsimple(h2GJet, h2Qcd, 100)
    else:
        h2QcdW, h2QcdWsys = qcdPrediction2d(h2GJet, h2Qcd, 100)
    # compute total weight for simple total scaling
    metCut = 100
    metCutBin = h2GJet.GetXaxis().FindBin(metCut) - 1
    wTot = h2GJet.Integral(0,metCutBin,0,-1)/h2Qcd.Integral(0,metCutBin,0,-1)
    for dir, cut1, cut2 in [
        ("x", 0, 1e6),
        ("x", 0, 2000),
        ("x", 0, 800),
        ("x", 0, 1000),
        ("x", 1000, 2000),
        ("x", 2000, 1e6),
        ("x", 700, 1e6),
        ("x", 800, 1e6),
        ("x", 900, 1e6),
        ("x", 1000, 1e6),
        ("x", 1100, 1e6),
        ("x", 1200, 1e6),
        ("x", 1300, 1e6),
        ("x", 1400, 1e6),
        ("x", 1500, 1e6),
        ("x", 2000, 1e6),
        ("y", 0, 1e6),
        ("y", 0, 100),
        ("y", 100, 1e6),
        ]:
        if dir=="x":
            cut1Bin = h2GJet.GetYaxis().FindBin(cut1)
            cut2Bin = h2GJet.GetYaxis().FindBin(cut2-1e-6)
            h1GJetMet = h2GJet.ProjectionX(aux.randomName(), cut1Bin, cut2Bin)
            h1QcdMet = h2Qcd.ProjectionX(aux.randomName(), cut1Bin, cut2Bin)
            h1QcdMetW = h2QcdW.ProjectionX(aux.randomName(), cut1Bin, cut2Bin)
            h1QcdMetWsys = h2QcdWsys.ProjectionX(aux.randomName(), cut1Bin, cut2Bin)
        elif dir=="y":
            cut1Bin = h2GJet.GetXaxis().FindBin(cut1)
            cut2Bin = h2GJet.GetXaxis().FindBin(cut2-1e-6)
            h1GJetMet = h2GJet.ProjectionY(aux.randomName(), cut1Bin, cut2Bin)
            h1QcdMet = h2Qcd.ProjectionY(aux.randomName(), cut1Bin, cut2Bin)
            h1QcdMetW = h2QcdW.ProjectionY(aux.randomName(), cut1Bin, cut2Bin)
            h1QcdMetWsys = h2QcdWsys.ProjectionY(aux.randomName(), cut1Bin, cut2Bin)
        else:
            print "Do not know what to do with", dir

        h1QcdMet.Scale(wTot)

        for h in h1GJetMet, h1QcdMet, h1QcdMetW, h1QcdMetWsys:
            h.SetYTitle("Events/Bin")
            #h.Scale(1, "width")
            aux.appendFlowBin( h )

        aux.drawOpt(h1GJetMet, "data")

        h1QcdMet.SetLineColor(rwth.violett50)
        h1QcdMet.drawOption_ = "hist"
        h1QcdMetW.drawOption_ = "hist e"

        for h in h1QcdMetW, h1QcdMetWsys:
            h.Scale(h1GJetMet.Integral()/h.Integral())
            h.SetLineColor(ROOT.kRed)
            if dir == "y": h.SetTitleOffset(1)

        aux.drawOpt(h1QcdMetWsys, "sys")

        c = ROOT.TCanvas()
        m = multiplot.Multiplot()
        m.add(h1GJetMet, "#gamma")
        m.add(h1QcdMetW, "Jet EMH_{T} weighted")
        m.add(h1QcdMetWsys, "Jet weight uncert.")
        m.add(h1QcdMet, "Jet unweighted")

        info = "E_{T}^{miss}/GeV" if dir=="y" else "H_{T}/GeV"
        if cut1: info = str(cut1)+"<"+info
        if cut2<1e5: info = info+"<"+str(cut2)
        if "<" not in info and ">" not in info: info=""
        m.leg.SetHeader(info)
        m.Draw()

        l = aux.Label(sim=True, info=dSets.label)

        r = ratio.Ratio( "#gamma/jet", h1GJetMet, h1QcdMetW, h1QcdMetWsys )
        r.draw(0.5,1.5)

        simpleString = "_simple" if simpleScale else ""
        aux.save("htReweighting_{}{}_{}_{}to{}".format(name,simpleString,dir,int(cut1),int(cut2)), normal=False)

def finalPrediction(allSets, simpleScale=False):
    jetSet = dataHt if allSets == data else allSets

    hName = "metRaw_vs_emht"
    m = re.match("(.*)_vs_(.*)", hName)
    if m:
        axVars = {"x":m.group(1), "y": m.group(2)}
    else:
        print "Could not determine x and y variable from", hName

    metBinning = aux.getBinnigsFromName("met")["3"]
    emhtBinning = aux.getBinnigsFromName("emht")["2"]
    #emhtBinning = [0, 800, 1000, 2000 ]
    #emhtBinning = [0, 100, 120, 130, 140, 150, 200, 250, 300, 400, 500]
    #emhtBinning = None

    h2data = allSets.getHist("tr/"+hName)
    h2jetControl = jetSet.getHist("tr_jControl/"+hName)
    h2eControl = allSets.getHist("tr_eControl/"+hName)
    h2ttg = ttg.getHist("tr/"+hName)
    h2tg = tg.getHist("tr/"+hName)
    h2wg = wg_mg.getHist("tr/"+hName)
    h2zg = zg_130.getHist("tr/"+hName)
    h2dy = zgll.getHist("tr/"+hName)
    h2s1 = signal["T5Wg_1550_100"].getHist("tr/"+hName)
    h2s2 = signal["T5Wg_1550_1500"].getHist("tr/"+hName)

    h2data = aux.rebin2d(h2data, metBinning, emhtBinning)
    h2jetControl = aux.rebin2d(h2jetControl, metBinning, emhtBinning)
    h2eControl = aux.rebin2d(h2eControl, metBinning, emhtBinning)
    h2ttg = aux.rebin2d(h2ttg, metBinning, emhtBinning)
    h2tg = aux.rebin2d(h2tg, metBinning, emhtBinning)
    h2wg = aux.rebin2d(h2wg, metBinning, emhtBinning)
    h2zg = aux.rebin2d(h2zg, metBinning, emhtBinning)
    h2dy = aux.rebin2d(h2dy, metBinning, emhtBinning)
    h2s1 = aux.rebin2d(h2s1, metBinning, emhtBinning)
    h2s2 = aux.rebin2d(h2s2, metBinning, emhtBinning)

    # subtract mc from numerator (small influence on denominator)
    num = h2data.Clone(aux.randomName())
    num.Add(h2ttg, -1)
    num.Add(h2tg, -1)
    num.Add(h2wg, -1)
    num.Add(h2zg, -1)
    num.Add(h2dy, -1)
    num.Add(h2eControl, -0.0197)

    if simpleScale:
        h2QcdW, h2QcdWsys = qcdPrediction2dsimple(num, h2jetControl, 100)
    else:
        h2QcdW, h2QcdWsys = qcdPrediction2d(num, h2jetControl, 100)

    if allSets == data:
        h2eControl.Scale(0.0197) # fakeRate Data
    else:
        h2eControl.Scale(0.0107) # fakeRate MC
    for h in h2ttg, h2tg, h2wg, h2zg, h2dy, h2s1, h2s2:
        h.Scale(0.983) # trigger efficiency, uncertainty: 0.002

    cuts = [
        ("x", 0, 1e6),
        ("x", 0, 800),
        ("x", 0, 1000),
        ("x", 800, 1000),
        ("x", 0, 2000),
        ("x", 0, 1500),
        ("x", 700, 1e6),
        ("x", 800, 1e6),
        ("x", 900, 1e6),
        ("x", 1000, 1e6),
        ("x", 1100, 1e6),
        ("x", 1200, 1e6),
        ("x", 1300, 1e6),
        ("x", 1400, 1e6),
        ("x", 1500, 1e6),
        ("x", 2000, 1e6),
        ("y", 0, 1e6),
        ("y", 0, 100),
        ("y", 100, 1e6),
        ]
    #if emhtBinning: cuts = [ ("x", low, high) for low, high in zip(emhtBinning,emhtBinning[1:]+[1e6]) ]

    if axVars["y"] == "emht" and False:
        cuts = [("x", 0, 1e6), ("x", 0, 800), ("x", 0, 1000), ("x", 800, 1000), ("x", 1000, 1e6), ("x", 2000, 1e6), ("y", 0, 1e6)]
    elif axVars["y"] == "mht":
        cuts = [("x", 0, 10), ("x", 40, 100), ("x", 100, 200), ("x", 200, 300), ("x", 300, 400), ("x", 400, 1e6), ("y", 0, 1e6)]
    elif axVars["y"] == "gPt":
        cuts = [("x", 100, 300), ("x", 300, 500), ("x", 500, 1e6), ("y", 0, 1e6)]
    elif axVars["y"] == "jetPt":
        cuts = [("x", 100, 400), ("x", 400, 600), ("x", 600, 1e6), ("y", 0, 1e6)]
    elif axVars["y"] == "njet":
        cuts = [("x", 0, 2), ("x", 3, 4), ("x", 5, 6), ("x", 7, 1e6), ("y", 0, 1e6)]

    for dir, cut1, cut2 in cuts:
        parDir = "xy".replace(dir, "") # get the other axis
        cut1Bin = aux.getAxis(h2data,parDir).FindBin(cut1)
        cut2Bin = aux.getAxis(h2data,parDir).FindBin(cut2-1e-6)
        h1data = aux.getProjection(h2data, dir, cut1Bin, cut2Bin)
        h1QcdW = aux.getProjection(h2QcdW, dir, cut1Bin, cut2Bin)
        h1QcdWsys = aux.getProjection(h2QcdWsys, dir, cut1Bin, cut2Bin)
        h1e = aux.getProjection(h2eControl, dir, cut1Bin, cut2Bin)
        h1ttg = aux.getProjection(h2ttg, dir, cut1Bin, cut2Bin)
        h1tg = aux.getProjection(h2tg, dir, cut1Bin, cut2Bin)
        h1wg = aux.getProjection(h2wg, dir, cut1Bin, cut2Bin)
        h1zg = aux.getProjection(h2zg, dir, cut1Bin, cut2Bin)
        h1dy = aux.getProjection(h2dy, dir, cut1Bin, cut2Bin)
        h1s1 = aux.getProjection(h2s1, dir, cut1Bin, cut2Bin)
        h1s2 = aux.getProjection(h2s2, dir, cut1Bin, cut2Bin)

        for h in h1data, h1QcdW, h1QcdWsys, h1e, h1ttg, h1tg, h1wg, h1zg, h1dy, h1s1, h1s2:
            h.SetYTitle("Events/Bin")
            aux.appendFlowBin(h)
            if dir == "y": h.SetTitleOffset(1)
        h1QcdW_integral = h1QcdW.Integral()
        if h1QcdW_integral:
            h1QcdW.Scale(h1data.Integral()/h1QcdW_integral)
            h1QcdWsys.Scale(h1data.Integral()/h1QcdW_integral)

        if allSets == data:
            if dir=="y": maxi=2000e5
            if dir=="x": maxi=160
            for bin in range(h1data.FindBin(maxi),h1data.GetNbinsX()+2):
                h1data.SetBinContent(bin,0)

        h1eSys = aux.getSysHisto(h1e,.3)
        h1ttgSys = aux.getSysHisto(h1ttg,.3)
        h1tgSys = aux.getSysHisto(h1tg,.3)
        h1wgSys = aux.getSysHisto(h1wg,.3)
        h1zgSys = aux.getSysHisto(h1zg,.5)
        h1dySys = aux.getSysHisto(h1dy,.5)

        aux.drawOpt(h1data, "data")
        h1QcdW.SetLineColor(rwth.blue75)
        h1e.SetLineColor(rwth.green)
        h1ttg.SetLineColor(rwth.bordeaux)
        h1tg.SetLineColor(rwth.bordeaux75)
        h1wg.SetLineColor(rwth.red75)
        h1zg.SetLineColor(rwth.yellow)
        h1dy.SetLineColor(rwth.yellow50)
        h1s1.SetLineWidth(3)
        h1s1.drawOption_ = "hist"
        h1s2.SetLineWidth(3)
        h1s2.drawOption_ = "hist"
        h1s2.SetLineStyle(2)

        dataLabel = "(Pseudo)Data"
        if allSets == data: dataLabel = "Data"
        if allSets == gjets+qcd: dataLabel = "#gammaJet+QCD"
        c = ROOT.TCanvas()
        m = multiplot.Multiplot()
        m.add(h1data, dataLabel)
        if allSets == data or True:
            m.addStack(h1tg, "t#gamma")
            m.addStack(h1dy, "DY")
            m.addStack(h1e, "e#rightarrow#gamma")
            m.addStack(h1ttg, "t#bar{t}#gamma")
            m.addStack(h1zg, "Z#gamma")
            m.addStack(h1wg, "W#gamma")
            m.add(h1s1, signal["T5Wg_1550_100"].label)
            m.add(h1s2, signal["T5Wg_1550_1500"].label)
        m.addStack(h1QcdW, "#gammaJet+QCD pred.")

        # systematics
        sysStack = ROOT.THStack()
        if allSets == data or True:
            systHists = h1QcdWsys, h1eSys, h1ttgSys, h1wgSys, h1zgSys, h1dy
        else:
            systHists = [h1QcdWsys]
        for h in systHists:
            sysStack.Add(h)
        sysHist = sysStack.GetStack().Last()
        aux.drawOpt(sysHist, "sys")

        info = aux.getAxis(h2data, parDir).GetTitle().replace(" (GeV)","")
        if cut1: info = str(cut1)+"<"+info
        if cut2<1e5: info = info+"<"+str(cut2)
        if "<" not in info and ">" not in info: info=""
        m.leg.SetHeader(info)
        m.Draw()
        sysHist.Draw("same e2")

        l = aux.Label(sim=allSets is not data)

        r = ratio.Ratio("Data/SM", h1data, m.hists[0].GetStack().Last(), sysHist)
        r.draw(0.5,1.5)

        appendix = ""
        if allSets == data: appendix = "_data"
        if any([ x.startswith("T5") for x in allSets.names ]): appendix = "_sigCont"

        title = "finalPlot"# if allSets == data else "closure"
        simpleString = "_simple" if simpleScale else ""
        aux.save("{}{}{}_{}_{}{}{}".format(title,appendix,simpleString,axVars[dir],int(cut1),axVars[parDir],int(cut2)), normal=False)

def metInfluence( dataset, savename="test", dirs=["tr"] ):
    c = ROOT.TCanvas()
    m = multiplot.Multiplot()

    for dir in dirs:
        for iname, name in enumerate(["emht", "njet", "jetPt", "gPt", "charged2JetPt"]):
            name = dir+"/metRaw_vs_"+name
            h2 = dataset.getHist(name)
            metBinning = aux.getBinnigsFromName("met")["3"]

            h2 = aux.rebin2d(h2, metBinning)
            p = h2.ProfileX()
            p.SetLineColor(iname+1)
            if iname == 5: p.SetLineColor(ROOT.kGreen+4)
            label = h2.GetYaxis().GetTitle().replace(" (GeV)","")
            if "gPt" in name: label = label.replace("^{jet}","")
            if "charged" in name: label += "charged"
            if dir.endswith("tr_jControl"):
                p.SetLineStyle(2)
                label += " jet"
            if dir.endswith("tr_jControl_wemht"):
                p.SetLineStyle(3)
                label += " jet EMH_{T} weighted"
            if name.endswith("njet"):
                p.Scale(200)
                label += "#times 200"
            m.add(p, label)

    m.Draw()
    aux.save("metInfluence_"+savename, log=False)

def metDependency( dataset, savename="test", dirs=["tr"] ):
    for name in ["emht", "njet", "jetPt", "gPt", "charged2JetPt"]:
        c = ROOT.TCanvas()
        m = multiplot.Multiplot()
        for idir, dir in enumerate(dirs):
            h2 = dataset.getHist(dir+"/metRaw_vs_"+name)
            yBinning = None
            if name == "emht": yBinning = aux.getBinnigsFromName("emht")["1"]
            if name == "jetPt" or "charged" in name: yBinning = aux.getBinnigsFromName("j1_pt")["1"]
            if name == "gPt": yBinning = aux.getBinnigsFromName("g_pt")["1"]

            h2 = aux.rebin2d(h2, None, yBinning)
            p = h2.ProfileY()
            p.SetLineColor(idir+1)
            p.SetTitleOffset(1)
            p.SetYTitle("#LT{}#GT".format(h2.GetXaxis().GetTitle()))
            if name == "njet": p.SetMaximum(120)
            m.add(p, dir)
        m.Draw()
        aux.save("metDependency_{}_{}".format(savename,name), log=False)

def gammaFakeRatio():
    metBinning = aux.getBinnigsFromName("met")["1"]

    c = ROOT.TCanvas()
    hg = gjets.getHist("tr/met")
    hj = qcd.getHist("tr/met")
    hg = aux.rebin(hg, metBinning)
    hj = aux.rebin(hj, metBinning)

    hj.Add(hg)
    hg.Divide(hj)

    hg.SetLineColor(1)
    hg.SetYTitle("#gamma+Jet/(#gamma+Jet+QCD)")
    hg.GetYaxis().SetRangeUser(0,1)
    hg.Draw("hist e")
    aux.save("gammaFakeRatio",log=False)



def main():
    pass
    #transitions()
    #compareAll( "_all", gjets400, gjets600, znn400, znn600 )
    #compareAll( "_GjetsVsZnn", gjets, znn )
    #compareAll( "_allMC", gjets, znn, qcd, wjets )

    #drawSameHistograms( "gqcd_data", [gjets, qcd], additional=[data])
    #drawSameHistograms( "mc_data", [gjets, qcd, ttjets, ttg, wjets,wg_mg,zg_130,znunu], additional=[signal["T5Wg_1550_100"],data])
    #drawSameHistograms( "mc", [gjets, qcd, ttjets, ttg, wjets,wg_mg,zg_130,znunu], additional=[signal["T5Wg_1550_100"],signal["T5Wg_1550_1500"]])

    #ewkClosures( ttjets, "tt" )
    #ewkClosures( wjets, "w" )
    #ewkClosures( wjets+ttjets, "ewk" )

    multiQcdClosures( gjets+qcd, "gqcd" )
    #multiQcdClosures( zg_130, "znunu", znunu)
    #multiQcdClosures( ttg, "tt",ttjets )
    #multiQcdClosures( wg_mg, "w", wjets )
    #multiQcdClosures( data, "data", dataHt )
    #multiQcdClosures( signal["T5Wg_1550_100"], "signal" )
    #multiQcdClosures( gjets+qcd+znunu+wjets+ttjets+wg_mg+ttg+zg_130, "mc", gjets+qcd+wjets+ttjets+znunu+wg_mg+ttg+zg_130)

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

    #significanceMetHt()

    #zToMet()

    #htStuff()

    #htRebinning(gjets+qcd, "gqcd")
    #htRebinning(gjets+qcd, "gqcd", simpleScale=True)
    #htRebinning(gjets400 + gjets600, "gqcd2", "tr", qcd700 + qcd1000 + qcd1500 + qcd2000 )
    #htRebinning(gjets+qcd, "gqcdVSall", "tr", gjets+qcd+ttjets+ttg+wjets+wg_mg+zg_130)
    #dirSets = gjets+qcd+wjets+ttjets+znunu
    #dirSets.label = "#gamma/MutiJet,W,#bar{t}t,Z(#nu#nu)"
    #htRebinning(dirSets, "gqcdtwzVSall", "tr_noGenE", gjets+qcd+ttjets+wjets+znunu+ttg+wg_mg+zg_130)
    #finalPrediction(data)
    #finalPrediction(data, simpleScale=True)
    #finalPrediction(gjets+qcd, simpleScale=True)
    #finalPrediction(gjets+qcd)
    #finalPrediction(gjets+qcd+ttjets+wjets+znunu+ttg+wg_mg+zg_130, simpleScale=True)
    #finalPrediction(gjets+qcd+ttjets+wjets+znunu+ttg+wg_mg+zg_130+signal["T5Wg_1550_1500"])


    #metInfluence( gjets, "gjet" )
    #metInfluence( data, "data" )
    #metInfluence( qcd, "qcd", ["tr_jControl"] )
    #metInfluence( gjets+qcd, "gqcd", ["tr","tr_jControl","tr_jControl_wemht"] )
    #metDependency( data, "data", ["tr","tr_jControl","tr_jControl_wemht"] )

    #gammaFakeRatio()

    #qcdClosure3d(gjets, "test", dirName="tr", preSet=qcd)

if __name__ == "__main__":
    from datasets import *
    #ROOT.gErrorIgnoreLevel = ROOT.kFatal
    main()

