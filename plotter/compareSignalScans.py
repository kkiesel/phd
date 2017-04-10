#!/usr/bin/env python2

from include import *


def compare(myFile, joFile, oldFile=""):
    includeSUS15012 = False
    saveName = "compareLimits_"
    info = ""

    h2 = aux.getFromFile(myFile, "obs_hist")
    h2.Reset()
    h2.GetXaxis().SetRangeUser(1200,h2.GetXaxis().GetXmax())
    if "T6gg" in myFile:
        saveName += "T6gg"
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        info = "pp #rightarrow #tilde{q}#tilde{q}, #tilde{q} #rightarrow q%s, %s #rightarrow #gamma#tilde{G}"%(lsp_s,lsp_s)
        h2.SetTitle(";m_{#tilde{q}} (GeV);m_{%s} (GeV)"%lsp_s)
    elif "T6Wg" in myFile:
        saveName += "T6Wg"
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        info = "pp #rightarrow #tilde{q}#tilde{q}, #tilde{q} #rightarrow q%s, %s #rightarrow #gamma/W^{#pm}#tilde{G}"%(lsp_s,lsp_s)
        h2.SetTitle(";m_{#tilde{q}} (GeV);m_{%s} (GeV)"%lsp_s)
    elif "T5gg" in myFile:
        saveName += "T5gg"
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        info = "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow qq%s, %s #rightarrow #gamma#tilde{G}"%(lsp_s,lsp_s)
        h2.SetTitle(";m_{#tilde{g}} (GeV);m_{%s} (GeV)"%lsp_s)
    elif "T5Wg" in myFile:
        saveName += "T5Wg"
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        info = "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow qq%s, %s #rightarrow #gamma/W^{#pm}#tilde{G}"%(lsp_s,lsp_s)
        h2.SetTitle(";m_{#tilde{g}} (GeV);m_{%s} (GeV)"%lsp_s)

    obsGr_knut = aux.getFromFile(myFile, "obs")
    expGr_knut = aux.getFromFile(myFile, "exp")
    obsGr_knut.SetLineColor(ROOT.kRed)
    expGr_knut.SetLineColor(ROOT.kRed)

    obsGr_jo = aux.getFromFile(joFile, "gr_obsC_sm")
    expGr_jo = aux.getFromFile(joFile, "gr_expC_sm")
    obsGr_jo.SetLineColor(ROOT.kBlue)
    expGr_jo.SetLineColor(ROOT.kBlue)

    for gr in obsGr_knut, expGr_knut, obsGr_jo, expGr_jo:
        gr.SetLineWidth(2)

    if oldFile:
        obsGr_old = aux.getFromFile(oldFile, "obs")
        expGr_old = aux.getFromFile(oldFile, "exp")
        obsGr_old.SetLineColor(ROOT.kMagenta)
        expGr_old.SetLineColor(ROOT.kMagenta)
        expGr_old.SetLineStyle(2)
        for gr in obsGr_old, expGr_old:
            gr.SetLineWidth(2)

    expGr_jo.SetLineStyle(2)
    expGr_knut.SetLineStyle(2)

    h2.GetYaxis().SetLimits(0,3200)
    if "T6" in myFile: h2.GetYaxis().SetLimits(0,2800)
    h2.Draw("colz")
    obsGr_knut.Draw("same")
    expGr_knut.Draw("same")
    obsGr_jo.Draw("same")
    expGr_jo.Draw("same")
    if oldFile:
        expGr_old.Draw("same")
        obsGr_old.Draw("same")

    if "gg" in myFile and includeSUS15012:
        obsGr_old2 = ROOT.TGraph()
        expGr_old2 = ROOT.TGraph()
        obsGr_old2.SetLineColor(ROOT.kGreen+2)
        expGr_old2.SetLineColor(ROOT.kGreen+2)
        expGr_old2.SetLineStyle(2)
        obsGr_old2.SetLineWidth(2)
        expGr_old2.SetLineWidth(2)
        if "T5gg" in myFile:
            obsGr_old2.SetPoint(0, 1560, 100)
            obsGr_old2.SetPoint(1, 1610, 500)
            obsGr_old2.SetPoint(2, 1610, 1400)
            expGr_old2.SetPoint(0, 1510, 100)
            expGr_old2.SetPoint(1, 1560, 500)
            expGr_old2.SetPoint(2, 1560, 1400)
        if "T6gg" in myFile:
            obsGr_old2.SetPoint(0, 1350, 100)
            obsGr_old2.SetPoint(1, 1370, 500)
            obsGr_old2.SetPoint(2, 1370, 1200)
            expGr_old2.SetPoint(0, 1330, 100)
            expGr_old2.SetPoint(1, 1350, 500)
            expGr_old2.SetPoint(2, 1350, 1200)
        obsGr_old2.Draw("same")
        expGr_old2.Draw("same")


    xmin = h2.GetXaxis().GetBinLowEdge(h2.GetXaxis().GetFirst())
    xmax = h2.GetXaxis().GetXmax()
    line = ROOT.TLine(xmin,xmin,xmax,xmax)
    line.SetLineColor(ROOT.kGray+2)
    line.SetLineStyle(3)
    line.Draw()

    text = ROOT.TLatex()
    text.DrawLatexNDC(0.15, .95, "#font[61]{CMS} #scale[0.76]{Preliminary}" )
    text.DrawLatexNDC(0.2, .87, info )
    text.DrawLatexNDC( .59, .95, "%.1f fb^{-1} (%s TeV)"%(aux.intLumi/1000., 13) )

    obsGr = obsGr_knut.Clone(aux.randomName())
    expGr = expGr_knut.Clone(aux.randomName())
    obsGr.SetLineColor(ROOT.kBlack)
    expGr.SetLineColor(ROOT.kBlack)

    legBlack = ROOT.TLegend(.67,.74,.97,.84)
    legBlack.SetFillStyle(0)
    legBlack.AddEntry(obsGr, "Observed")
    legBlack.AddEntry(expGr, "Expected")
    legBlack.Draw()

    leg = ROOT.TLegend(.16,.74,.65,.84)
    if oldFile: leg = ROOT.TLegend(.16,.65,.69,.84)
    leg.SetFillStyle(0)
    leg.AddEntry(obsGr_knut, "SUS-16-047: high EMH_{T}")
    leg.AddEntry(obsGr_jo, "SUS-16-046: high S_{T}^{#gamma}")
    if oldFile: leg.AddEntry(obsGr_old, "SUS-16-023: high S_{T}^{#gamma} 2.3 fb^{-1}")
    if "gg" in myFile and includeSUS15012:
        leg.AddEntry(obsGr_old2, "SUS-15-012: 2#gamma 2.3fb^{-1}", "l")
    leg.Draw()

    aux.save("forChristian_"+saveName, log=False)




if __name__ == "__main__":
    compare("limitCalculations/T6gg_v11/saved_graphs1d_limit.root", "/home/home4/institut_1b/kiesel/other_photon_limits/CMS-SUS-16-046-T6gg.root")
    compare("limitCalculations/T6Wg_v11/saved_graphs1d_limit.root", "/home/home4/institut_1b/kiesel/other_photon_limits/CMS-SUS-16-046-T6wg.root")
    compare("limitCalculations/T5gg_v11/saved_graphs1d_limit.root", "/home/home4/institut_1b/kiesel/other_photon_limits/CMS-SUS-16-046-T5gg.root", "/home/home4/institut_1b/kiesel/other_photon_limits/CMS-PAS-SUS-16-023-T5gg.root")
    compare("limitCalculations/T5Wg_v11/saved_graphs1d_limit.root", "/home/home4/institut_1b/kiesel/other_photon_limits/CMS-SUS-16-046-T5wg.root")
