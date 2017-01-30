#!/bin/zsh
rm ~/phd/fakeRate/SingleElectron_Run2016_fake.root
hadd ~/phd/fakeRate/SingleElectron_Run2016_fake.root ~/phd/fakeRate/SingleElectron_Run2016(B|C|D|E|F|G|H)-*_fake.root
rm ~/phd/fakeRate/DYJetsToLL_M-50_merged_fake.root
hadd ~/phd/fakeRate/DYJetsToLL_M-50_merged_fake.root ~/phd/fakeRate/DYJetsToLL_M-50_fake.root ~/phd/fakeRate/DYJetsToLL_M-50_ext_fake.root

