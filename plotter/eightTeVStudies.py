#! /usr/bin/env python2
# -*- coding: utf-8 -*-

# Python libraries
import ConfigParser


# Own includes
import style
import multiplot
import auxiliary as aux

xVarCfg = ConfigParser.SafeConfigParser()
xVarCfg.read( "xVar.cfg" )

datasetCfg = ConfigParser.SafeConfigParser()
datasetCfg.read( "dataset.cfg" )

selectionCfg = ConfigParser.SafeConfigParser()
selectionCfg.read( "selection.cfg" )

if __name__ == "__main__":

    for xVarCfgEntry in xVarCfg.sections():
        print xVarCfg.get( xVarCfgEntry, "title" )

