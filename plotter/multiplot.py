import ROOT

class Multiplot:
    def __init__( self ):
        self.hists = []
        self.histsToStack = []

        # todo: impmlement setter and getter
        self.minimum = None
        self.maximum = None

        self.leg = ROOT.TLegend(.7,.7,.93,.92)
        self.leg.SetFillColor( ROOT.kWhite )

    def add( self, h, label="" ):
        h.SetName(label)
        self.hists.append( h )

    def addStack( self, h, label="" ):
        h.SetName( label )
        self.histsToStack.append( h )

    def getMinimum( self ):
        return min( [ h.GetMinimum() for h in self.hists ] )

    def getMaximum( self ):
        return max( [ h.GetMaximum() for h in self.hists ] )

    def stackHists( self ):
        if not self.histsToStack:
            return
        stack = ROOT.THStack()
        stack.SetTitle( ";%s;%s"%(self.histsToStack[0].GetXaxis().GetTitle(),self.histsToStack[0].GetYaxis().GetTitle()) )
        stack.drawOption_ = ""
        for h in self.histsToStack:
            h.Sumw2( False )
            h.SetFillColor( h.GetLineColor() )
            h.SetLineColor( ROOT.kBlack )
            stack.Add( h )

        self.hists.append( stack )


    def Draw( self ):
        if not self.hists and not self.histsToStack:
            return False
        self.stackHists()

        minimum = self.getMinimum()
        maximum = self.getMaximum()

        if ROOT.gPad and ROOT.gPad.GetLogy():
            maximum = 2.5*maximum
            minimum = 0.3*minimum
        else:
            maximum = maximum + (maximum-minimum)*.1

        if self.maximum != None:
            maximum = self.maximum
        if self.minimum != None:
            minimum = self.minimum

        if minimum==0: minimum=0.1

        # fill data first in legend?

        # fill legend (in correct order)
        for h in self.histsToStack[-1::-1]:
            h.SetLineColor(0)
            self.leg.AddEntry( h, h.GetName(), "f" )
        for h in self.hists:
            if isinstance( h, ROOT.THStack ): continue
            self.leg.AddEntry( h, h.GetName(), "l" )

        # change the order for drawing
        self.hists.reverse()
        self.hists[0].SetMinimum( minimum )
        self.hists[0].SetMaximum( maximum )
        self.hists[0].Draw( self.hists[0].drawOption_ )

        for h in self.hists[1:]:
            h.Draw( "same %s"%h.drawOption_ )

        self.leg.Draw()

        return True
