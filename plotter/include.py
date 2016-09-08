#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys
if sys.version_info[:2] == (2,6):
    print "Initialize correct python version first!"
    sys.exit()

import suppressor
with suppressor.suppress_stdout_stderr():
    import ROOT

# private libs
import ratio
import style
import multiplot
from rwthColors import rwth
import auxiliary as aux
from datasets import *

import math
import collections
