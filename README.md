# Software used for CMS-PAS-SUS-16-047 #

nTuples are assumed to be generated using the [TreeWriter]{https://github.com/cms-susy-photon-rwth-1b/TreeWriter}
software.

## Creating histograms ##

There are several classes producing histograms from the nTuples.
In classes inheriting from ROOT::TSelector, the event processing and filling of
histograms is defined. There are 'run.py' files, in which single files can be
processed, and 'process.py' files, in which several files can be processed,
either on the local PC or on the local condor cluster.

### histogramProducer/HistogramProducer.cc ##
Used for the background estimation

### histogramProducer/SignalScan.cc ###
Fills pTmiss histograms for several signal scans. They are also filled in the
HistogramProducer to check consistency between these two classes. (You have to
synchronize those by hand!)

### fakeRate/FakeRateSelector.cc ###
Selects ee and e#gamma events to calculate R(gumma/e)


## Plotting histograms ##

Plotting scripts are found in the "plotter" directory. Most important are:

### finalPredictions.py ###

### signalScan.py ###

### datasets.py ###
Here the datasets are defined, including the link to the files produced by the
histogramProducer

### main.py ###
Many different functions, eg creating the event yield table, covariances, old versions
of the background estimation, ...

### dataMCcomparisons.py ###

### triggerStudies.py ###
