import ConfigParser

def getlist( self, section, option, delimiter=" " ):
    str = self.get( section, option )
    return str.split( delimiter )

def getints( self, section, option, delimiter=" " ):
    return [ int(i) for i in self.getlist( section, option, delimiter ) ]

def getfloats( self, section, option, delimiter=" " ):
    return [ float(i) for i in self.getlist( section, option, delimiter ) ]

ConfigParser.SafeConfigParser.getlist = getlist
ConfigParser.SafeConfigParser.getints = getints
ConfigParser.SafeConfigParser.getfloats = getfloats

cfg = ConfigParser.SafeConfigParser()
cfg.readfp(open('rebin.cfg'))

