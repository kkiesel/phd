import ROOT
# Usage of the corporate design colors of RWTH Aachen University
# RGB values copied from here:
# http://www9.rwth-aachen.de/global/show_document.asp?id=aaaaaaaaaadrbpl#page=27

# definition of the colors:
colors = {
    "blue": (0,84,159),
    "blue75": (64,127,183),
    "blue50": (142,186,229),
    "blue25": (199,221,242),
    "blue10": (232,241,250),
    "black": (0,0,0),
    "black75": (100,101,103),
    "black50": (156,158,159),
    "black25": (207,209,210),
    "black10": (236,237,237),
    "magenta": (227,0,102),
    "magenta75": (233,96,136),
    "magenta50": (241,158,177),
    "magenta25": (249,210,218),
    "magenta10": (253,238,240),
    "yellow": (255,237,0),
    "yellow75": (255,240,85),
    "yellow50": (255,245,155),
    "yellow25": (255,250,209),
    "yellow10": (255,253,238),
}

# This class will contain all colors as members
class rwth:
    colors = []
    pass

for i, (name,(r,g,b)) in enumerate(colors.iteritems()):
    exec "rwth.{} = {:d}".format(name,1001+i)
    rwth.colors.append( ROOT.TColor(1001+i,r/255.,g/255.,b/255.) )

# Usage: histogram.SetLineColor( rwth.blue50 )
