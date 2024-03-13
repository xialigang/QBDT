
import string
import os
from ROOT import *
import ROOT

#from array import array

import sys

if len(sys.argv)<2:
   print('Not specify the treedir')
   sys.exit(1)
if len(sys.argv)<3:
   print('Not specify the number of systematics')
   print('Will NOT calculate QBDT score for systematic samples')
if len(sys.argv)<4:
   print('Not specify the number of trees')
   print('Use the mininum of 100 and the number of trees under directory of ', str(sys.argv[1]))

print('sys.argv =', sys.argv)
treedir0 = str(sys.argv[1])

Nsyst = 0
if len(sys.argv)>2:
    Nsyst = int(sys.argv[2])

maxtree0 = 100
if len(sys.argv) > 3:
   maxtree0 = int(sys.argv[3])

gROOT.LoadMacro("AtlasStyle.C")
from ROOT import SetAtlasStyle
SetAtlasStyle()

gROOT.SetBatch(True)

from qbdtmodule import qbdtmodule
#from adabdtmodule import adabdtmodule
#from gradbdtmodule import gradbdtmodule


qbdt = 1
bdt = qbdtmodule()

realtest = 1

treename0 = 'NOMINAL_train'
if realtest:
   treename0 = 'NOMINAL'

bdt.treename = treename0

bdt.signals = {'h2atata': ( 'h2atata', 'root_dir/fhist_h2atata.root') }
bdt.backgrounds = {'atata' : ('atata', 'root_dir/fhist_atata.root') }
bdt.weightvar = 'fweight'
bdt.trainflagvar = 'trainflag'
bdt.trainflagcut = 0.5

#bdt.totNsig = 12464*0.00585756
#bdt.totNbkg = 4161*35.3655

bdt.totNsig = 12338*0.00585756
bdt.totNbkg = 4082*35.3655

bdt.treedir = treedir0

bdt.maxtree = 1
list_files = os.listdir(bdt.treedir)
ntree = -1
for a in list_files:
   if 'tree_' not in a:
      continue
   n = a.replace('tree_', '')
   n = n.replace('.txt', '')
   n = int(n)
   if n >= bdt.maxtree:
      bdt.maxtree = n+1

print('bdt.maxtree =', bdt.maxtree)

if bdt.maxtree < maxtree0:
   print('The number of trees', bdt.maxtree,'is <',maxtree0)
else:
   bdt.maxtree = maxtree0

print('bdt.maxtree =', bdt.maxtree)



bdt.mindQ = 0.0 #0.001
scaleNbins = 1
bdt.variables = {
'pt_lep' : ('pt_lep', 'p_{T}(l) [GeV]', '',                         [int(scaleNbins*60), 20, 80], 'norm1'),
'pt_tau' : ('pt_tau', 'p_{T}(#tau) [GeV]', '',                      [int(scaleNbins*60), 20, 80], 'norm1'),
'pt_pho' : ('pt_pho', 'p_{T}(#gamma) [GeV]', '',                    [int(scaleNbins*60), 20, 80], 'norm1'),
'pt_met' : ('pt_met', 'p_{T}(MET) [GeV]', '',                       [int(scaleNbins*60), 0, 60], 'norm1'),
'm_lephad' : ('m_lephad', 'm(l#tau) [GeV]', '',                     [int(scaleNbins*60), 10, 120], 'norm1'),
'm_lephadpho' : ('m_lephadpho', 'm(l#tau#gamma) [GeV]', '',         [int(scaleNbins*75), 50,150], 'norm1'),
'm_lephadmet' : ('m_lephadmet', 'm(l#tauMET) [GeV]', '',            [int(scaleNbins*80), 20, 180], 'norm1'),
'm_lephadphomet' : ('m_lephadphomet', 'm(l#tau#gammaMET) [GeV]', '',[int(scaleNbins*65), 70, 200], 'norm1'),
}

bdt.bkgsysts = {
      #name : ( name, treename, sample, filename, filepath )
      'tes' : ('tes', 'NOMINAL', 'atata', 'fhist_atata_syst_tes.root', 'root_dir/', 'high'),
      'tauid' : ('tauid', 'NOMINAL', 'atata', 'fhist_atata_syst_tauid.root', 'root_dir/', 'high'),
      'met' : ('met', 'NOMINAL', 'atata', 'fhist_atata_syst_met.root', 'root_dir/', 'high'),
      'leppt' : ('leppt', 'NOMINAL', 'atata', 'fhist_atata_syst_leppt.root', 'root_dir/', 'high'),
      'lepid' : ('lepid', 'NOMINAL', 'atata', 'fhist_atata_syst_lepid.root', 'root_dir/', 'high'),
      'phopt' : ('phopt', 'NOMINAL', 'atata', 'fhist_atata_syst_phopt.root', 'root_dir/', 'high'),
      'phoid' : ('phoid', 'NOMINAL', 'atata', 'fhist_atata_syst_phoid.root', 'root_dir/', 'high'),
      }
bdt.corrmatrix = [
      [1, ],
      ]

if Nsyst == 1:
   syst = ['tes' ]
elif Nsyst == 3:
   syst = ['tes', 'tauid', 'met' ]
elif Nsyst == 7:
   syst = ['tes', 'tauid', 'met', 'leppt', 'lepid', 'phopt', 'phoid']
else:
   syst = []

# load root files
bdt.load_files(realtest)
# show input variables
bdt.show_variables()
# generate QBDT score for nominal and systematic samples
# show QBDT score distribution for nominal samples
bdt.show_performance(realtest, syst)

