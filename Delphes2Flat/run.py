#!/usr/bin/env python

import sys
if len(sys.argv) < 3:
    print "%s INPUT.root OUTPUT_PREFIX"
    sys.exit(1)
inputFile = sys.argv[1]
prefix = sys.argv[2]

import os
from ROOT import *
delphesPath = "Delphes"
gSystem.AddIncludePath('-I"%s"' % delphesPath)
gSystem.AddDynamicPath(delphesPath)
gSystem.AddLinkedLibs('-L"%s"' % delphesPath)
gSystem.Load("libDelphes")

gROOT.ProcessLine(".L makeFlatTuple.C+")

outFile = "./output/%s%s" % (prefix, os.path.basename(inputFile))
gROOT.ProcessLine('makeFlatTuple("%s", "%s");' % (inputFile, outFile))

