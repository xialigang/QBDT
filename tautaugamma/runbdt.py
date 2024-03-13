
import string
import os
from ROOT import *
import ROOT

import sys

#from array import array

print('sys.argv =', sys.argv)
if len(sys.argv)<2:
   print('No specify the trees dir')
   sys.exit(1)
if len(sys.argv)<3:
   print('Not specify the number of systematic sources')
   sys.exit(1)


treedir0 = str(sys.argv[1])
Nsyst0 = int(sys.argv[2])

# flag to switch on systematics or not into training
syston0 = 0
if len(sys.argv)<4:
   print('Not specify whether to turn on systematics or not')
   print('Do NOT consider systmatics by default')
else:
   syston0 = int(sys.argv[3])

if Nsyst0 == 0:
    syston0 = 0

ntrees0 = 3

if len(sys.argv) < 5:
   print('Not specifying number of trees')
   print('Build 100 trees by default')
else:
   ntrees0 = int(sys.argv[4])


if not os.path.isdir(treedir0):
   os.mkdir(treedir0)

gROOT.LoadMacro("AtlasStyle.C")
from ROOT import SetAtlasStyle
SetAtlasStyle()

gROOT.SetBatch(True)

from qbdtmodule import qbdtmodule
#from adabdtmodule import adabdtmodule
#from gradbdtmodule import gradbdtmodule


bdt = qbdtmodule()
treename0 = 'NOMINAL_train'

bdt.syston = syston0


bdt.treename = treename0

bdt.signals = {'h2atata': ( 'h2atata', 'root_dir/fhist_h2atata.root') }
bdt.backgrounds = {'atata' : ('atata', 'root_dir/fhist_atata.root') }
bdt.weightvar = 'fweight'
bdt.trainflagvar = 'trainflag'
bdt.trainflagcut = 0.5

#not relevant for the training
#but you can set it correctly
#so that it will give the correct significance somewhere ...
bdt.totNsig = 1. #21158*0.006
bdt.totNbkg = 1. #50176*0.803996;

# max number of nodes
bdt.maxnode = 10 #40 #20

# max number of trees
bdt.maxtree = ntrees0

# training results and testing results will be put here
bdt.treedir = treedir0

# a cut to stop the tree splitting
# set it to 0 by default (no cut at all and recommended)
bdt.mindQ = 0.0 #0.001

scaleNbins = 1.0


bdt.variables = {
        # input variable name : (name, xtitle, 'other purpose', binning, 'log (show log y-axis) or nomr1 (renomarlized to unit area)')
'pt_lep' : ('pt_lep', 'p_{T}(l) [GeV]', '',                         [int(scaleNbins*60), 20, 80], 'norm1'),
'pt_tau' : ('pt_tau', 'p_{T}(#tau_{had}) [GeV]', '',                      [int(scaleNbins*60), 20, 80], 'norm1'),
'pt_pho' : ('pt_pho', 'p_{T}(#gamma) [GeV]', '',                    [int(scaleNbins*70), 10, 80], 'norm1'),
'pt_met' : ('pt_met', 'E_{T}^{miss} [GeV]', '',                       [int(scaleNbins*60), 0, 60], 'norm1'),
'm_lephad' : ('m_lephad', 'm(l#tau_{had}) [GeV]', '',                     [int(scaleNbins*60), 10, 120], 'norm1'),
'm_lephadpho' : ('m_lephadpho', 'm(l#tau_{had}#gamma) [GeV]', '',         [int(scaleNbins*75), 50,150], 'norm1'),
'm_lephadmet' : ('m_lephadmet', 'm(l#tau_{had}#nu) [GeV]', '',            [int(scaleNbins*80), 20, 180], 'norm1'),
'm_lephadphomet' : ('m_lephadphomet', 'm(l#tau_{had}#gamma#nu) [GeV]', '',[int(scaleNbins*65), 70, 200], 'norm1'),

}



if Nsyst0 == 1:
   bdt.bkgsysts = {
      # 1 syst
      # systematics name: ( name, tree name, background sample name, root file name, root file directory, high/low variation,
      'tes' : ('tes', treename0, 'atata', 'fhist_atata_syst_tes.root', 'root_dir/', 'high'),
   }
elif Nsyst0 == 3:
   bdt.bkgsysts = {
      # 3 systs
      'tes' : ('tes', treename0, 'atata', 'fhist_atata_syst_tes.root', 'root_dir/', 'high', 0.31),
      'tauid' : ('tauid', treename0, 'atata', 'fhist_atata_syst_tauid.root', 'root_dir/', 'high', 0.16),
      'met' : ('met', treename0, 'atata', 'fhist_atata_syst_met.root', 'root_dir/', 'high', 0.40),
   }
elif Nsyst0 == 7:
   bdt.bkgsysts = {
      # 7 systs
      'tes' : ('tes', treename0, 'atata', 'fhist_atata_syst_tes.root', 'root_dir/', 'high', 0.40),
      'tauid' : ('tauid', treename0, 'atata', 'fhist_atata_syst_tauid.root', 'root_dir/', 'high', 0.42),
      'met' : ('met', treename0, 'atata', 'fhist_atata_syst_met.root', 'root_dir/', 'high', 0.57),
      'leppt' : ('leppt', treename0, 'atata', 'fhist_atata_syst_leppt.root', 'root_dir/', 'high', 0.53),
      'lepid' : ('lepid', treename0, 'atata', 'fhist_atata_syst_lepid.root', 'root_dir/', 'high', 0.53),
      'phopt' : ('phopt', treename0, 'atata', 'fhist_atata_syst_phopt.root', 'root_dir/', 'high', 0.71),
      'phoid' : ('phoid', treename0, 'atata', 'fhist_atata_syst_phoid.root', 'root_dir/', 'high', 0.44),
   }
else:
   bdt.syston = 0

# load files
bdt.load_files()
# building trees
bdt.build_trees()

