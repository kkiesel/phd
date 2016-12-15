#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from include import *

ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
ROOT.RooMsgService.instance().setSilentMode(True)

# settings
fnameData = "../fakeRate/SingleElectron_Run2016-PromptReco_fake.root"
fnameMC = "../fakeRate/DYJetsToLL_M-50_ext_fake.root"

selectionTex = {
    "all":"p_{T}>25",
    "40pt":"p_{T}>40",
    "EB":"EB",
    "EB_40pt":"EB, p_{T}>40",
    "EB_40pt_3nTracks":"EB, p_{T}>40, #geq3tr",
    "EB_01eta_40pt":"EB, p_{T}>40, |#eta|>0.1",
    "EB_100pt":"EB, p_{T}>100",
    "EE_40pt": "EE, p_{T}>40",
    "EE_40pt_3nTracks":"EE, p_{T}>40 #geq3tr",
}

def drawUnfitted(name, hData, hSig, hBkg, infoText=""):
    var = "m_{e#gamma}" if "num" in name else "m_{ee}"
    aux.drawOpt(hData, "data")
    hData.SetTitle(";{} (GeV);Pairs [abitrary normalization]".format(var))
    hData.Draw("ep")
    dataMax = hData.GetMaximum() if hSig.GetMaximum() else 1.
    if hSig.GetMaximum():
        hSig.SetLineColor(ROOT.kRed)
        hSig.Scale(hData.GetMaximum()/hSig.GetMaximum())
        hSig.Draw("hist same")
    if hBkg.Integral():
        hBkg.SetLineColor(ROOT.kGreen)
        hBkg.Scale(hData.Integral()/hBkg.Integral())
        hBkg.Draw("hist same")
    l = aux.Label(info=infoText, sim=False)
    leg = ROOT.TLegend(.6, .6, .93, .92)
    leg.SetFillStyle(0)
    leg.AddEntry(hData, "Data")
    leg.AddEntry(hSig, "Signal shape")
    leg.AddEntry(hBkg, "Background shape")
    leg.Draw()
    aux.save(name+"_raw", log=False)

def histToKdePdf(h, x, name):
    centers = aux.array.array( 'd', [h.GetBinCenter(i) for i in range(1, h.GetNbinsX()+1)] )
    contents = aux.array.array( 'd', [h.GetBinContent(i) for i in range(1, h.GetNbinsX()+1)] )
    kde = ROOT.TKDE(len(centers), centers, contents, h.GetBinLowEdge(1), h.GetBinLowEdge(h.GetNbinsX()+1), "", 1./aux.sqrt(12))
    func = kde.GetFunction()
    func.Draw()
    pdf = ROOT.RooTFnPdfBinding(name, "", func, ROOT.RooArgList(x))
    return pdf

def fitHist(name, hData, hSig, hBkg, infoText=""):
    drawUnfitted(name, hData, hSig, hBkg, infoText)
    if not hData.GetEntries() or not hSig.GetEntries() or not hBkg.GetEntries(): return 0, 0
    totalInt = hData.Integral(0,-1)
    var = "m_{e#gamma}" if "num" in name else "m_{ee}"
    x = ROOT.RooRealVar("x", var, 60, 120, "GeV")

    # datahist
    dh = ROOT.RooDataHist("dhData", "", ROOT.RooArgList(x), ROOT.RooFit.Import(hData))

    # background shape
    pdfBkg = histToKdePdf(hBkg, x, "pdfBkg")
    nBkg = ROOT.RooRealVar("nBkg", "number of events", totalInt/10, 0, 2*totalInt)
    bkg = ROOT.RooExtendPdf("bkg", "", pdfBkg, nBkg)

    # signal shape
    dhSig = ROOT.RooDataHist("dhSig", "", ROOT.RooArgList(x), ROOT.RooFit.Import(hSig))
    pdfSig = ROOT.RooHistPdf("histpdf1", "", ROOT.RooArgSet(x), dhSig, 2)
    #pdfSig = histToKdePdf(hSig, x, "pdfSig")
    meanSig = ROOT.RooRealVar("meanSig", "meanSig", 0, -2, 2)
    widthSig = ROOT.RooRealVar("widthSig", "widthSig", 1, 0, 3)
    smearSig = ROOT.RooGaussian("gausSig", "", x, meanSig, widthSig)
    smearedSig = ROOT.RooFFTConvPdf("smearedSig","", x, pdfSig, smearSig)
    nSig = ROOT.RooRealVar("nSig","number of events", totalInt, 0, 2*totalInt)
    signal = ROOT.RooExtendPdf("signal", "", smearedSig, nSig)

    # add signal and background shapes
    total = ROOT.RooAddPdf("total", "sig+bkg", ROOT.RooArgList(signal, bkg))

    # fit
    x.setBins(10000, "cache") # for the convolution
    x.setRange("peak", 60, 120) # range can be decreased to avoid boundary effects
    total.fitTo(dh, ROOT.RooFit.Range("peak"), ROOT.RooFit.PrintEvalErrors(-1))

    # draw
    frame = x.frame(ROOT.RooFit.Title(" "))
    dh.plotOn(frame)
    frame.SetYTitle(frame.GetYaxis().GetTitle().replace("Events", "Pairs"))
    total.plotOn(frame)
    total.plotOn(frame, ROOT.RooFit.Components("pdfBkg"), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.LineStyle(ROOT.kDashed))
    total.plotOn(frame, ROOT.RooFit.Components("smearedSig"), ROOT.RooFit.LineColor(ROOT.kRed))

    residuals = frame.residHist()
    frameRes = x.frame(ROOT.RooFit.Title(" "))
    frameRes.addPlotable(residuals, "P")

    c = ROOT.TCanvas()
    frame.Draw()
    #ratio.createBottomPad()
    #frameRes.Draw()
    l = aux.Label(info=infoText, sim=False)
    aux.save(name+"_fit", log=False)
    return nSig.getVal(), nSig.getError()

def draw1dEfficiency(savename, eff, selection="", xbins=None):
    c = ROOT.TCanvas()
    if xbins:
        num = aux.rebin(eff.GetCopyPassedHisto(), xbins, scale=False)
        den = aux.rebin(eff.GetCopyTotalHisto(), xbins, scale=False)
        eff = ROOT.TEfficiency(num, den)
    g = eff.CreateGraph()
    g.SetMaximum(0.05)
    g.SetMinimum(0)
    g.Draw("ap")
    l = aux.Label(sim=True, info="Gen-Match "+selectionTex[selection.replace("_gen","")])
    aux.save(savename, log=False)

def inclusiveFakeRate(savename, variable, selection, isData, intOnly=False):
    effMc = aux.getFromFile(fnameMC, "{}{}/{}".format(selection, "" if selection.endswith("_gen") else "_Zmatch", variable) )
    nDim = effMc.GetDimension()
    if nDim == 1:
        print "not implemented"
        return 0

    hMcNum = effMc.GetCopyPassedHisto()
    hMcDen = effMc.GetCopyTotalHisto()

    if isData:
        effData = aux.getFromFile(fnameData, "{}/{}".format(selection, variable))
        hDataNum = effData.GetCopyPassedHisto()
        hDataDen = effData.GetCopyTotalHisto()
        effBkg = aux.getFromFile(fnameData, "{}_bkg/{}".format(selection, variable))
        hBkg = effBkg.GetCopyTotalHisto()

    if nDim == 2:
        hMcNum_mll = hMcNum.ProjectionX(aux.randomName(), 0, -1)
        hMcDen_mll = hMcDen.ProjectionX(aux.randomName(), 0, -1)
        xMin = hMcNum_mll.GetXaxis().FindFixBin(80)
        xMax = hMcNum_mll.GetXaxis().FindFixBin(100)-1
        if isData:
            hDataNum_mll = hDataNum.ProjectionX(aux.randomName(), 0, -1)
            hDataDen_mll = hDataDen.ProjectionX(aux.randomName(), 0, -1)
            hBkg_mll = hBkg.ProjectionX(aux.randomName(), 0, -1)
            if not hBkg_mll.Integral():
                hBkg_mll = hBkg.ProjectionX(aux.randomName(), 0, -1)
            if intOnly:
                num = hDataNum_mll.Integral(xMin, xMax)
                den = hDataDen_mll.Integral(xMin, xMax)
            else:
                num, numE = fitHist("{}_bin{}_num".format(savename,"incl"), hDataNum_mll, hMcDen_mll, hBkg_mll, "")
                den, denE = fitHist("{}_bin{}_den".format(savename,"incl"), hDataDen_mll, hMcDen_mll, hBkg_mll, "")
        else:
            num = hMcNum_mll.Integral(xMin, xMax)
            den = hMcDen_mll.Integral(xMin, xMax)
        if den:
            return num/den*100
        else:
            return 0
    elif nDim == 3:
        print "not implemented"
        return 0
        style.style2d()
        c = ROOT.TCanvas()
        style.defaultStyle()

        hOut = hMcNum.Project3D("yz")
        hOut.Reset("ICESM")
        for xbin in range(1, hOut.GetNbinsX()+1):
            print 1.*xbin/(hOut.GetNbinsX()+1)
            for ybin in range(1, hOut.GetNbinsY()+1):
                num, den = 0, 0
                hMcNum_mll = hMcNum.ProjectionX(aux.randomName(), ybin, ybin, xbin, xbin)
                hMcDen_mll = hMcDen.ProjectionX(aux.randomName(), ybin, ybin, xbin, xbin)
                xMin = hMcNum_mll.GetXaxis().FindFixBin(80)
                xMax = hMcNum_mll.GetXaxis().FindFixBin(100)-1
                if isData:
                    hDataNum_mll = hDataNum.ProjectionX(aux.randomName(), ybin, ybin, xbin, xbin)
                    hDataDen_mll = hDataDen.ProjectionX(aux.randomName(), ybin, ybin, xbin, xbin)
                    hBkg_mll = hBkg.ProjectionX(aux.randomName(), ybin, ybin, xbin, xbin)
                    if not hBkg_mll.Integral():
                        hBkg_mll = hBkg.ProjectionX(aux.randomName())
                    if intOnly:
                        num = hDataNum_mll.Integral(xMin, xMax)
                        den = hDataDen_mll.Integral(xMin, xMax)
                    else:
                        infoText = "bin{}_{}".format(xbin,ybin)
                        num, numE = fitHist("{}_bin{}_{}_num".format(savename,xbin,ybin), hDataNum_mll, hMcDen_mll, hBkg_mll, infoText)
                        den, denE = fitHist("{}_bin{}_{}_den".format(savename,xbin,ybin), hDataDen_mll, hMcDen_mll, hBkg_mll, infoText)
                else:
                    num = hMcNum_mll.Integral(xMin, xMax)
                    den = hMcDen_mll.Integral(xMin, xMax)
                if den:
                    hOut.SetBinContent(xbin, ybin, 100*num/den)
                    hOut.SetBinError(xbin, ybin, 100*aux.sqrt(num)/den*aux.sqrt(1.+num/den))
        hOut.Draw("colz")
        hOut.SetZTitle("f_{e#rightarrow#gamma} (%)")
        hOut.SetMinimum(0)
        hOut.SetMaximum(5)
        hOut.SetTitle("")

        """
        m = multiplot.Multiplot()
        ROOT.gStyle.SetPalette(55)
        proj = aux.getProjections(h2Num, scale=False)
        style.defaultStyle()
        for ih, h in enumerate(proj):
            h.SetMinimum(0)
            h.SetMaximum(5)
            h.drawOption_ = "hist e"
            m.add(h, h.GetName())
        m.Draw()
        """
        l = aux.Label(drawAll=False, sim=not isData, info=selectionTex[selection])
        l.cms.SetX(0.01)
        l.cms.SetY(0.05)
        l.pub.SetX(0.01)
        l.pub.SetY(0.01)
        l.draw()
        aux.save(savename, log=False)


def binnedFakeRate(variable, selection, isData, xbins=None, ybins=None, intOnly=False, drawInclusive=False, saveToFile=False):
    binning = [None, xbins, ybins]
    savename = "fakeRate/"
    savename += "data" if isData else "sim"
    if intOnly: savename += "_intOnly"
    if drawInclusive: savename += "_errorBand"
    if xbins: savename += "_x{}".format(len(xbins))
    if ybins: savename += "_y{}".format(len(ybins))
    savename += "_{}_{}".format(selection, variable)
    effMc = aux.getFromFile(fnameMC, "{}{}/{}".format(selection, "" if selection.endswith("_gen") else "_Zmatch", variable) )
    nDim = effMc.GetDimension()
    if nDim == 1:
        draw1dEfficiency(savename+"_1d", effMc, selection, xbins)
        return

    hMcNum = effMc.GetCopyPassedHisto()
    hMcDen = effMc.GetCopyTotalHisto()
    hMcDen.Add(hMcNum, -1)
    hMcNum = aux.rebinX(hMcNum, *binning)
    hMcDen = aux.rebinX(hMcDen, *binning)

    if isData:
        effData = aux.getFromFile(fnameData, "{}/{}".format(selection, variable))
        hDataNum = effData.GetCopyPassedHisto()
        hDataDen = effData.GetCopyTotalHisto()
        hDataDen.Add(hDataNum, -1)
        effBkg = aux.getFromFile(fnameData, "{}_bkg/{}".format(selection, variable))
        hBkg = effBkg.GetCopyTotalHisto()

        hDataNum = aux.rebinX(hDataNum, *binning)
        hDataDen = aux.rebinX(hDataDen, *binning)
        hBkg = aux.rebinX(hBkg, *binning)

    if nDim == 2:
        aux.appendFlowBin2d(hMcNum, mergeX=False)
        aux.appendFlowBin2d(hMcDen, mergeX=False)
        if isData:
            aux.appendFlowBin2d(hDataNum, mergeX=False)
            aux.appendFlowBin2d(hDataDen, mergeX=False)
            aux.appendFlowBin2d(hBkg, mergeX=False)

        hOut = hMcNum.ProjectionY()
        hOut.Reset("ICESM")
        for bin in range(1, hOut.GetNbinsX()+1):
            #if bin not in [14]: continue
            ax = hMcNum.GetYaxis()
            yMin, yMax = ax.GetBinLowEdge(bin), ax.GetBinUpEdge(bin)
            if not yMin - int(yMin): yMin = int(yMin)
            if not yMax - int(yMax): yMax = int(yMax)
            infoText = "{} < {} < {}".format(yMin, ax.GetTitle(), yMax)
            hMcNum_mll = hMcNum.ProjectionX(aux.randomName(), bin, bin)
            hMcDen_mll = hMcDen.ProjectionX(aux.randomName(), bin, bin)
            xMin = hMcNum_mll.GetXaxis().FindFixBin(80)
            xMax = hMcNum_mll.GetXaxis().FindFixBin(100)-1
            if isData:
                hDataNum_mll = hDataNum.ProjectionX(aux.randomName(), bin, bin)
                hDataDen_mll = hDataDen.ProjectionX(aux.randomName(), bin, bin)
                hBkg_mll = hBkg.ProjectionX(aux.randomName(), bin, bin)
                if not hBkg_mll.Integral():
                    hBkg_mll = hBkg.ProjectionX(aux.randomName(), 0, -1)
                if intOnly:
                    num = hDataNum_mll.Integral(xMin, xMax)
                    den = hDataDen_mll.Integral(xMin, xMax)
                    numE, denE = aux.sqrt(num), aux.sqrt(den)
                else:
                    num, numE = fitHist("{}_bin{}_num".format(savename,bin), hDataNum_mll, hMcDen_mll, hBkg_mll, infoText)
                    den, denE = fitHist("{}_bin{}_den".format(savename,bin), hDataDen_mll, hMcDen_mll, hBkg_mll, infoText)
            else:
                num = hMcNum_mll.Integral(xMin, xMax)
                den = hMcDen_mll.Integral(xMin, xMax)
                numE, denE = aux.sqrt(num), aux.sqrt(den)
            if den:
                hOut.SetBinContent(bin, 100*num/den)
                hOut.SetBinError(bin, 100*aux.sqrt((numE/den)**2 + (denE*num/den**2)**2))

        hOut.SetTitle("")
        hOut.SetYTitle("N_{e#gamma}/N_{ee} (%)")
        hOut.SetTitleOffset(0.9,"y")
        aux.drawOpt(hOut, "data")
        hOut.SetMaximum(5)
        hOut.SetMinimum(0)
        if drawInclusive: incFakeRate = inclusiveFakeRate(savename, variable, selection, isData, intOnly)
        c = ROOT.TCanvas()
        style.defaultStyle()
        hOut.Draw("e hist")
        l = aux.Label(sim=not isData, info=selectionTex[selection])
        if drawInclusive:
            hError = ROOT.TH1F("","",1, hOut.GetBinLowEdge(0), hOut.GetBinLowEdge(hOut.GetNbinsX()+2))
            hError.SetBinContent(1, incFakeRate)
            hError.SetBinError(1, hError.GetBinContent(1)*.3)
            aux.drawOpt(hError, "sysUnc")
            hError.Draw("same"+hError.drawOption_)
            text = ROOT.TLatex()
            text.DrawLatexNDC(.6, .15, "N_{{e#gamma}}/N_{{ee}} = {:.2f}%".format(incFakeRate))
        if variable == "vtx":
            res = hOut.Fit("pol1", "QS")
            f = hOut.GetFunction("pol1")
            f.SetLineColor(ROOT.kRed)
            f.Draw("same")
            text = ROOT.TLatex()
            text.DrawLatexNDC(.2, .7, "f(x)={:.3f}+{:.4f}x".format(res.Parameter(0), res.Parameter(1)))
            text.DrawLatexNDC(.2, .6, "#chi^{{2}}/NDF = {:.2f}".format(res.Chi2()/res.Ndf()))
        if saveToFile:
            if not hOut.GetBinContent(hOut.GetNbinsX()+1): hOut.SetBinContent(hOut.GetNbinsX()+1,hOut.GetNbinsX()+1,  hOut.GetBinContent(hOut.GetNbinsX())) # fill overflow bin
            hOut.Scale(0.01) # from percent to normal
            aux.write2File(hOut, savename.replace("/","__"), "weights.root")
        aux.save(savename, log=False)
    elif nDim == 3:
        aux.appendFlowBin3d(hMcNum, mergeX=False)
        aux.appendFlowBin3d(hMcDen, mergeX=False)
        if isData:
            aux.appendFlowBin3d(hDataNum, mergeX=False)
            aux.appendFlowBin3d(hDataDen, mergeX=False)
            aux.appendFlowBin3d(hBkg, mergeX=False)
        style.style2d()
        c = ROOT.TCanvas()
        style.defaultStyle()

        hOut = hMcNum.Project3D("yz")
        hOut.Reset("ICESM")
        for xbin in range(1, hOut.GetNbinsX()+1):
            print 1.*xbin/(hOut.GetNbinsX()+1)
            for ybin in range(1, hOut.GetNbinsY()+1):
                num, den = 0, 0
                hMcNum_mll = hMcNum.ProjectionX(aux.randomName(), ybin, ybin, xbin, xbin)
                hMcDen_mll = hMcDen.ProjectionX(aux.randomName(), ybin, ybin, xbin, xbin)
                xMin = hMcNum_mll.GetXaxis().FindFixBin(80)
                xMax = hMcNum_mll.GetXaxis().FindFixBin(100)-1
                if isData:
                    hDataNum_mll = hDataNum.ProjectionX(aux.randomName(), ybin, ybin, xbin, xbin)
                    hDataDen_mll = hDataDen.ProjectionX(aux.randomName(), ybin, ybin, xbin, xbin)
                    hBkg_mll = hBkg.ProjectionX(aux.randomName(), ybin, ybin, xbin, xbin)
                    if not hBkg_mll.Integral():
                        hBkg_mll = hBkg.ProjectionX(aux.randomName())
                    if intOnly:
                        num = hDataNum_mll.Integral(xMin, xMax)
                        den = hDataDen_mll.Integral(xMin, xMax)
                        numE, denE = aux.sqrt(num), aux.sqrt(den)
                    else:
                        infoText = "bin{}_{}".format(xbin,ybin)
                        num, numE = fitHist("{}_bin{}_{}_num".format(savename,xbin,ybin), hDataNum_mll, hMcDen_mll, hBkg_mll, infoText)
                        den, denE = fitHist("{}_bin{}_{}_den".format(savename,xbin,ybin), hDataDen_mll, hMcDen_mll, hBkg_mll, infoText)
                else:
                    num = hMcNum_mll.Integral(xMin, xMax)
                    den = hMcDen_mll.Integral(xMin, xMax)
                    numE, denE = aux.sqrt(num), aux.sqrt(den)
                if den:
                    hOut.SetBinContent(xbin, ybin, 100*num/den)
                    hOut.SetBinError(xbin, ybin, 100*aux.sqrt((numE/den)**2 + (denE*num/den**2)**2))
        hOut.Draw("colz")
        hOut.SetZTitle("N_{e#gamma}/N_{ee} (%)")
        hOut.SetMinimum(0)
        hOut.SetMaximum(5)
        hOut.SetTitle("")
        l = aux.Label(drawAll=False, sim=not isData, info=selectionTex[selection])
        l.cms.SetX(0.01)
        l.cms.SetY(0.05)
        l.pub.SetX(0.01)
        l.pub.SetY(0.01)
        l.draw()
        aux.save(savename, log=False)

        if True:
            c = ROOT.TCanvas()
            m = multiplot.Multiplot()
            ROOT.gStyle.SetPalette(55)
            proj = aux.getProjections(hOut, scale=False)
            style.defaultStyle()
            for ih, h in enumerate(proj):
                h.SetMinimum(0)
                h.SetMaximum(5)
                h.drawOption_ = "hist e"
                m.add(h, h.GetName())
            m.Draw()
            l = aux.Label(sim=not isData, info=selectionTex[selection])
            aux.save(savename+"_x", log=False)

            m = multiplot.Multiplot()
            ROOT.gStyle.SetPalette(55)
            proj = aux.getProjections(hOut, "y", scale=False)
            style.defaultStyle()
            for ih, h in enumerate(proj):
                h.SetMinimum(0)
                h.SetMaximum(5)
                h.drawOption_ = "hist e"
                m.add(h, h.GetName())
            m.Draw()
            l = aux.Label(sim=not isData, info=selectionTex[selection])
            aux.save(savename+"_y", log=False)





binnings = {
    "jets": [-.5, .5, 1.5, 2.5],
    "jetPt": [0,5, 30, 40, 50, 80],
    "vtx": [0.5, 10.5, 15.5, 20.5, 25.5],
    #"pt": [15, 20, 25, 30, 40, 60, 90],
    "pt": range(30, 90, 5) + [90,100,120,150,200],
    "met": range(0,50,5)+[50,60,70,100,120],
    "cIsoWorst": range(0, 8) + [8, 10, 15, 20],
    "cIso": [0, .1, .5, 1., 1.5, 2.0, 2.5, 3.5],
    "pIso": [0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1., 1.5, 2],
    "hoe": list(aux.frange(0, 0.05, 0.005)),
    "nTracksPV": [-.5, .5, 1.5, 2.5, 3.5, 4.5, 25.5, 49.5],
}

if __name__ == "__main__":
    selections = "all" , "40pt", "EB", "EB_40pt" , "EB_01eta_40pt", "EB_100pt"
    variables = "pt", "met", "emht", "eta", "jets", "vtx", "sie", "sip", "hoe", "r9", "cIso", "nIso", "pIso", "cIsoWorst", "nTracksPV" \
       "pt_vs_vtx", "pt_vs_nTracksPV", "met_vs_vtx", "met_vs_jets", "met_vs_pt", "met_vs_jetPt", "eta_vs_phi", "cIso_vs_vtx"

    #binnedFakeRate(variable, selection, isData, xbins=None, ybins=None, intOnly=False, drawInclusive=False)
    # Main variables
    #binnedFakeRate("pt", "EB_40pt", True, drawInclusive=True)
    #binnedFakeRate("pt", "EB_40pt", True, xbins=binnings["pt"], drawInclusive=True)
    #binnedFakeRate("eta", "EB_40pt", True, drawInclusive=True)
    #binnedFakeRate("met", "EB_40pt", True, binnings["met"], drawInclusive=True)
    #binnedFakeRate("vtx", "EB_40pt", True, drawInclusive=True)
    #binnedFakeRate("cIso", "EB_40pt", True, xbins=binnings["cIso"], drawInclusive=True)
    #binnedFakeRate("emht", "EB_40pt", True, drawInclusive=True)

    # Main variables simulation
    #binnedFakeRate("pt", "EB_40pt", False, drawInclusive=True)
    #binnedFakeRate("pt", "EB_40pt", False, xbins=binnings["pt"], drawInclusive=True)
    #binnedFakeRate("eta", "EB_40pt", False, drawInclusive=True)
    #binnedFakeRate("met", "EB_40pt", False, binnings["met"], drawInclusive=True)
    #binnedFakeRate("vtx", "EB_40pt", False, drawInclusive=True)
    #binnedFakeRate("cIso", "EB_40pt", False, xbins=binnings["cIso"], drawInclusive=True)
    #binnedFakeRate("emht", "EB_40pt", False, drawInclusive=True)

    #binnedFakeRate("cIso_vs_vtx", "EB_40pt", True, binnings["cIso"], binnings["vtx"])
    #binnedFakeRate("cIso_vs_vtx", "EB_40pt_3nTracks", True, binnings["cIso"], binnings["vtx"])


    """
    for stat in True, False: # isData
        #for sel in "all", "EB_40pt", "40pt":
        #    binnedFakeRate("pt", sel, stat, binning=[None, binnings["pt"]])
        #    binnedFakeRate("eta", sel, stat)
        #binnedFakeRate("met", "EB_40pt", stat, binning=[None, binnings["met"]])
        #binnedFakeRate("vtx", "EB_40pt", stat)
        #binnedFakeRate("cIso", "EB_40pt", stat)
        #binnedFakeRate("emht", "EB_40pt", stat)
        #binnedFakeRate("jets", "EB_40pt", stat)
        pass
        binnedFakeRate("cIso", "EB_40pt", stat, intOnly=True)
        binnedFakeRate("nIso", "EB_40pt", stat, intOnly=True)
        binnedFakeRate("pIso", "EB_40pt", stat, [None, binnings["pIso"]], intOnly=True)
    """
    #binnedFakeRate("nTracksPV", "all", isData=False, intOnly=True)
    #binnedFakeRate("nTracksPV", "EB_40pt", isData=True, intOnly=False, binning=[None, aux.frange(-.5, 36,1)+[39.5, 49.5]])
    #binnedFakeRate("nTracksPV", "EB_40pt", isData=False, intOnly=True)
    #binnedFakeRate("nTracksPV", "EB_01eta_40pt", isData=False, intOnly=True)
    #binnedFakeRate("nTracksPV", "all", isData=False, intOnly=True)
    #binnedFakeRate("nTracksPV", "all", isData=True, intOnly=False)
    #binnedFakeRate("eta_vs_phi", "all", isData=False, intOnly=True)
    #binnedFakeRate("eta_vs_phi", "all", isData=True, intOnly=True)
    #binnedFakeRate("nTracksPV", "all", isData=False, intOnly=True)
    #binnedFakeRate("nTracksPV", "all", isData=True, intOnly=True, binning=[None, binnings["nTracksPV"]])
    #binnedFakeRate("pt_vs_nTracksPV", "all", isData=True, intOnly=True, binning=[None, binnings["pt"], binnings["nTracksPV"]])
    #binnedFakeRate("pt_vs_nTracksPV", "all", isData=True, intOnly=True)
    #binnedFakeRate("vtx_vs_nTracksPV", "EB_40pt", isData=True, intOnly=True)
    #binnedFakeRate("vtx_vs_nTracksPV", "EB_40pt", isData=False, intOnly=True)
    #binnedFakeRate("vtx_vs_nTracksPV", "EB_40pt", isData=True, intOnly=False)
    #binnedFakeRate("vtx", "EB_40pt", isData=False, intOnly=True)
    #binnedFakeRate("vtx", "EB_40pt", isData=True, intOnly=False)
    #binnedFakeRate("cIso", "EB_40pt_3nTracks", isData=True, intOnly=False)
    #binnedFakeRate("nTracksPV", "EB_40pt_3nTracks", isData=True, intOnly=True)
    #binnedFakeRate("vtx", "EB_40pt_3nTracks", isData=True, intOnly=True)
    #binnedFakeRate("vtx", "EE_40pt_3nTracks", isData=True, intOnly=True)
    #binnedFakeRate("nTracksPV", "EE_40pt", isData=True, intOnly=True)
    #binnedFakeRate("pt", "EB_40pt_3nTracks", isData=True, binning=[None, binnings["pt"]], intOnly=True)
    #binnedFakeRate("vtx", "EB_40pt_3nTracks", isData=True, intOnly=True)
    #binnedFakeRate("pt_vs_nTracksPV", "all", isData=True, intOnly=True)
    #binnedFakeRate("pt_vs_nTracksPV", "all", True, binning=[None,binnings["pt"],binnings["vtx"]],intOnly=True)
    #binnedFakeRate("cIso_vs_vtx", "all", True, binning=[None,binnings["cIso"],binnings["vtx"]],intOnly=True)
    #binnedFakeRate("cIso_vs_vtx", "EB_40pt_3nTracks", True, binning=[None,binnings["cIso"],binnings["vtx"]],intOnly=True)

    #binnedFakeRate("met", "EB_40pt_3nTracks", isData=True, xbins=binnings["met"], intOnly=False, drawInclusive=True)
    #binnedFakeRate("pt", "EB_40pt_3nTracks", isData=True, xbins=binnings["pt"], intOnly=False, drawInclusive=True)
    #binnedFakeRate("emht", "EB_40pt_3nTracks", isData=True, intOnly=False, drawInclusive=True)
    #binnedFakeRate("eta", "EB_40pt_3nTracks", isData=True, intOnly=False, drawInclusive=True)
    #binnedFakeRate("jets", "EB_40pt_3nTracks", isData=True, intOnly=False, drawInclusive=True)
    #binnedFakeRate("vtx", "EB_40pt_3nTracks", isData=True, intOnly=False, drawInclusive=True)
    #binnedFakeRate("sie", "EB_40pt_3nTracks", isData=True, intOnly=False, drawInclusive=True)
    #binnedFakeRate("cIso", "EB_40pt_3nTracks", isData=True, intOnly=False, drawInclusive=True)
    #binnedFakeRate("nTracksPV", "EB_40pt_3nTracks", isData=True, intOnly=False, drawInclusive=True)
    #binnedFakeRate("pt", "EB_gen", False, xbins=binnings["pt"])
    #binnedFakeRate("pt", "EB", False, xbins=binnings["pt"])
    #binnedFakeRate("pt", "EB_100pt", True, xbins=binnings["pt"], drawInclusive=True)
    #binnedFakeRate("pt", "EB", True, xbins=binnings["pt"], drawInclusive=True)
    #binnedFakeRate("eta", "40pt", True, saveToFile=True)
    #binnedFakeRate("eta", "40pt", False, saveToFile=True)
    #binnedFakeRate("pt", "EB_40pt", True, xbins=binnings["pt"], drawInclusive=True)
    #binnedFakeRate("pt", "EB_40pt", False, xbins=binnings["pt"], drawInclusive=True)
    #binnedFakeRate("pt", "EB_40pt", True, xbins=binnings["pt"], drawInclusive=True)

    #binnedFakeRate("met", "EB_40pt", True, binnings["met"], intOnly=False, drawInclusive=True)
    #binnedFakeRate("pt", "EB_40pt", True, binnings["pt"], intOnly=False, drawInclusive=True)
    #binnedFakeRate("emht", "EB_40pt", True, intOnly=False, drawInclusive=True)
    #binnedFakeRate("eta", "EB_40pt", True, intOnly=False, drawInclusive=True)
    #binnedFakeRate("jets", "EB_40pt", True, binnings["jets"], intOnly=False, drawInclusive=True)
    #binnedFakeRate("vtx", "EB_40pt", True, intOnly=False, drawInclusive=True)
    #binnedFakeRate("sie", "EB_40pt", True, intOnly=False, drawInclusive=True)
    #binnedFakeRate("cIso", "EB_40pt", True, intOnly=False, drawInclusive=True)
    #binnedFakeRate("nTracksPV", "EB_40pt", True, binnings["nTracksPV"], intOnly=False, drawInclusive=True)
    binnedFakeRate("pt", "EE_40pt", True, binnings["pt"], intOnly=False, drawInclusive=True)

