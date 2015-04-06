import ROOT

# todo: replace with 'style.sethisttopmargin()'
def increaseMinMax( minimum, maximum ):
     if not ROOT.gPad and ROOT.gPad.GetLogy():
        maximum = 2.5*maximum
        minimum = 0.3*minimum
    else:
        maximum = maximum + (maximum-minimum)*.1
        minimum = minimum - (maximum-minimum)*.1
    return minimum, maximum


class Multiplot:
    def __init__( self ):
        self.hists = []
        self.histsToStack = []

        # todo: impmlement setter and getter
        self._minimum = None
        self._maximum = None

        self.leg = ROOT.TLegend(.7,.7,.93,.92)

    def add( self, h ):
        # important histograms first, they will be drawn on top
        self.hists.append( h )

    def addStack( self, h ):
        self.histsToStack.append( h )

    @property
    def minimum( self ):
        return min( [ h.GetMinimum(0) for h in self.hists ] ) if self._minimum == None else self._minimum

    @property
    def maximum( self ):
        return max( [ h.GetMaximum(0) for h in self.hists ] ) if self._maximum == None else self._maximum

    @minimum.setter
    def minimum( self, x ):
        self._minimum = x

    @maximum.setter
    def maximum( self, x ):
        self._maximum = x

    def stackHists( self ):
        if not self.histsToStack:
            return
        stack = ROOT.THStack( self.histsToStack[0] ) # get title, etc
        for h in histsToStack:
            h.SetFillColor( h.GetLineColor() )
            h.SetLineColor( ROOT.kBlack )
            stack.Add( h )

        self.hists.append( stack )


    def Draw( self ):
        self.stackHists()

        minimum, maximum = inceaseMinMax( self.minimum, self.maximum )

        # fill legend (in correct order)
        for h in self.hists:
            # todo: special cases?
            self.leg.AddEntry( h, h.GetTitle(), "lpe" )
        for h in self.histsToStack:
            self.leg.AddEntry( h, h.GetTitle(), "f" )

        # change the order for drawing
        h.reverse()
        hists[0].SetMinimum( minimum )
        hists[0].SetMaximum( maximum )
        hists[0].Draw()

        for h in hists[1:]:
            h.Draw("same")
