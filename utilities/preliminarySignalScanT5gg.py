#!/usr/bin/env python

import re
import subprocess
import ROOT
import style
style.style2d()

out = subprocess.check_output('das_client.py --limit 0 --query="dataset=/SMS-T5gg/weinberg-mGluino-*_mNeutralino-*/USER instance=prod/phys03"', shell=True )
out = subprocess.check_output('das_client.py --limit 0 --query="dataset=/SMS-T5gg/kiesel-SMS-T5gg_mGluino-*_mNeutralino-*3d7be4403ea17498be45eb057fcb0278/USER instance=prod/phys03"', shell=True )

scan = ROOT.TH2F("", ";m(#tilde{g}) (GeV);m(#tilde{#chi}^{0}_{1}) (GeV)", 5, 1100, 2100, 9, 100, 1900 )

for line in out.split("\n"):
    m = re.match( ".*mGluino-(\d+)_mNeutralino-(\d+)[-_].*", line )
    if m:
        mG = int(m.group(1))
        mN = int(m.group(2))
        scan.Fill( mG, mN )


scan.Draw("colz")
ROOT.gPad.SaveAs("SMS-T5gg_miniaods.pdf")

