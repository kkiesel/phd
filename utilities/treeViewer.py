import ROOT

files = ["/home/home4/institut_1b/kiesel/nTuples/T5Wg_1000_0.root"]
files = ["root://xrootd-cms.infn.it///store/user/kiesel/13TeV/nTuples/SMS-T5gg_mGluino-1900-1950_mLSP-0to1800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/v11/160407_145105/0000/photonTree_%s.root"%i for i in range(1,25)]

ch = ROOT.TChain("TreeWriter/eventTree")

for f in files:
    ch.AddFile(f)

models = []

for e in ch:
    mn = e.modelName.data()
    if mn not in models:
        print "adding", mn
        models.append(mn)
