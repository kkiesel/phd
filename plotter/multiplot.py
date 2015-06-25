import ROOT

class Multiplot:
    def __init__( self ):
        self.hists = []
        self.histsToStack = []

        # todo: impmlement setter and getter
        self.minimum = None
        self.maximum = None

        self.leg = ROOT.TLegend(.7,.7,.93,.92)

    def add( self, h, label="" ):
        # important histograms first, they will be drawn on top
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
            return
        self.stackHists()

        minimum = self.getMinimum()
        maximum = self.getMaximum()

        if ROOT.gPad and ROOT.gPad.GetLogy():
            maximum = 2.5*maximum
            minimum = 0.3*minimum
        else:
            maximum = maximum + (maximum-minimum)*.1
            minimum = max(0.1,minimum - (maximum-minimum)*.1)

        if self.maximum != None:
            maximum = self.maximum
        if self.minimum != None:
            minimum = self.minimum

        # fill legend (in correct order)
        for h in self.hists:
            if isinstance( h, ROOT.THStack ): continue
            # todo: special cases?
            self.leg.AddEntry( h, h.GetName(), "lpe" )
        for h in self.histsToStack[-1::-1]:
            self.leg.AddEntry( h, h.GetName(), "f" )

        # change the order for drawing
        self.hists.reverse()
        self.hists[0].SetMinimum( minimum )
        self.hists[0].SetMaximum( maximum )
        self.hists[0].Draw( self.hists[0].drawOption_ )

        for h in self.hists[1:]:
            h.Draw( "same %s"%h.drawOption_ )

        self.leg.Draw()
