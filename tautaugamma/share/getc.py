
from ROOT import *
import ROOT
import math


gROOT.LoadMacro("AtlasStyle.C")
from ROOT import SetAtlasStyle
SetAtlasStyle()

gROOT.SetBatch(True)

def calsig(hs, hb):
    nbins = hs.GetNbinsX()
    #print 'nbins =', nbins
    Q = 0
    for i in range(nbins):
        ns = hs.GetBinContent(i+1)
        nb = hb.GetBinContent(i+1)
        if nb == 0:
            #print 'bin',i,' == 0'
            continue
        Q += 2*((ns+nb)*math.log(1+ns/nb)-ns)
    Z = math.sqrt(Q)
    return Z

def gethist(tree, treename, var, wcut):
    tree.Draw(var+'>>htmp(20,-1,1)',wcut)
    h1 = gROOT.FindObject('htmp')
    h1.SetName(treename)
    #print h1, h1.Integral()
    return h1

def getZ(c=1., bdtmode = -1, ntrees = 100):
    if bdtmode < 0:
       if bdtmode == -1:
          file0 = TFile('official_BDTG/TMVA.root', 'read')
       else:
          file0 = TFile('official_BDTGp/TMVA.root', 'read')
       tree = file0.Get('dataset/TestTree')
       var = 'tanh('+str(c)+'*atanh(BDTG))'
       wcut = 'fweight'
       hb = gethist(tree, 'hb', var, wcut+'*(classID==1)')
       hs = gethist(tree, 'hs', var, wcut+'*(classID==0)')
    else:
       treedir = 'trees'+str(bdtmode)
       tail = ''
       tail = '_'+str(ntrees)
       files = TFile(treedir+'/bdt_sig'+tail+'.root', 'read')
       trees = files.Get('bdt')
       fileb = TFile(treedir+'/bdt_bkg'+tail+'.root', 'read')
       treeb = fileb.Get('bdt')
       var = 'tanh('+str(c)+'*qbdt)'
       wcut = 'fweight*(trainflag != 1)'
       hb = gethist(treeb, 'hb', var, wcut)
       hs = gethist(trees, 'hs', var, wcut)
    #print hb, hb.Integral()
    #print hs, hs.Integral()
    Z = calsig(hs, hb)
    okflag = 0
    nlast = 0.
    nlast = hb.GetBinContent(hb.GetNbinsX())
    if nlast > 0:
        okflag = 1
    return okflag, nlast, Z

def getc(bdtmode=0, ntrees=100):
    print '************************************working on bdtmode =', bdtmode, 'ntrees =',ntrees
    Npx = 100
    nmax = 0
    for i in range(Npx+1):
        c = 1./Npx*i
        okflag, nlast, Z = getZ(c,bdtmode, ntrees)
        if okflag and nmax<10:
            print c, '  ', okflag, '    ',nlast, '  ',Z
            nmax +=1
        #print c, '  ', okflag, '    ',nlast, '  ',Z
    return


getc(-1)
getc(-2)
getc(0)
getc(1)
getc(3)
getc(7)
getc(7, 150)
getc(7, 200)
#getc(3, 100)
#getc(3, 150)
#getc(3, 200)
