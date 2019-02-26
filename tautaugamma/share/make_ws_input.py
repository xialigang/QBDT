
import os
from ROOT import *
import ROOT
import math
from array import array


gROOT.LoadMacro("AtlasStyle.C")
from ROOT import SetAtlasStyle
SetAtlasStyle()

gROOT.SetBatch(True)



channels = [
      'SR',
      ]
samples = [
      'data',
      'sig',
      'bkg',
      ]
systs = [
      'nominal',
      'tes',
      'tauid',
      'met',
      'leppt',
      'lepid',
      'phopt',
      'phoid',
      ]

#totNsig = 12464*0.00585756
#totNbkg = 4161*35.3655

def prepare_dirs(filename):
   file0 = TFile(filename, 'recreate')
   for chan in channels:
      file0.mkdir(chan)
      dir0 = file0.GetDirectory(chan)
      dir0.cd()
      for samp in samples:
         dir0.mkdir(samp)
         dir1 = dir0.GetDirectory(samp)
   file0.Write()
   return
#tmpsigchain.Draw(var+'>>hq_sig_train'+histopt,'trainflag==1')
def get_hist(chain, var, cut, histopt='(40, -1., 1.)', remove_zero = 0):
   #print chain
   chain.Draw(var+'>>htmp'+histopt, cut)
   hist = gROOT.FindObject('htmp')
   print hist, hist.Integral()
   if remove_zero:
      nbins = hist.GetNbinsX()
      for i in range(nbins):
         x = hist.GetBinContent(i+1)
         if x == 0:
            if i==0:
               hist.SetBinContent(i+1, hist.GetBinContent(2)/2.)
               hist.SetBinError(i+1, hist.GetBinError(2))
            elif i == nbins -1:
               hist.SetBinContent(i+1, hist.GetBinContent(nbins-1)/2.)
               hist.SetBinError(i+1, hist.GetBinError(nbins-1))
            else:
               hist.SetBinContent(i+1, (hist.GetBinContent(i) + hist.GetBinContent(i+2))/2.)
               hist.SetBinError(i+1, math.sqrt(pow(hist.GetBinError(i),2) + pow(hist.GetBinError(i+2),2))/2.)
      print 'after removing zero bins', hist.Integral()
   return hist
def add_hist(h1, filename, chan, samp, systname):
   file0 = TFile(filename, 'update')
   if samp == '':
      dir1 = file0.GetDirectory(chan)
      dir1.cd()
      h1.Write(systname)
      return
   if systname != 'nominal':
      systname = systname+'_high'
   dir2 = file0.GetDirectory(chan+'/'+samp)
   dir2.cd()
   h1.Write(systname)
   file0.Save()
   return
def get_chain(rootdir, chan, samp, syst, ntrees=100):
   chain  = TChain('bdt')
   tail = '.root'
   if ntrees > 0 :
      tail = '_'+str(ntrees)+'.root'
   if samp == 'data':
      chain.Add(rootdir+'/'+'bdt_bkg'+tail)
   elif samp == 'sig':
      chain.Add(rootdir+'/'+'bdt_sig'+tail)
   elif samp == 'bkg':
      if syst == 'nominal':
         chain.Add(rootdir+'/'+'bdt_bkg'+tail)
      else:
         chain.Add(rootdir+'/'+'bdt_bkg_'+syst+tail)
   else:
      print 'get_chain ERROR !!!'
      return None
   return chain

def binhist(h1, binning, histname):
   #h1.Scale(1./h1.Integral())
   ly = []
   ldy = []
   lx = array('f')
   nbins = len(binning) - 1
   #print 'nbins =', nbins
   for i in range(len(binning)):
      if i == 0:
         lx.append(h1.GetBinLowEdge(binning[0]))
      else:
         lx.append(h1.GetBinLowEdge(binning[0]) + binning[i]*h1.GetXaxis().GetBinWidth(1))
      if i == 0:
         continue
      ia = binning[i-1]
      if i > 1:
         ia += 1
      ib = binning[i]
      n = 0
      dn = 0
      for j in range(ia, ib+1):
         n += h1.GetBinContent(j)
         dn+= pow(h1.GetBinError(j), 2)
      dn = math.sqrt(dn)
      ly.append(n)
      ldy.append(dn)
   #h1new = TH1F(histname, '', nbins, lx)
   h1new = TH1F(histname, '', nbins, 0, nbins)
   ysum = 0.
   for i in range(nbins):
      ysum += ly[i]
      bw = h1new.GetBinWidth(i+1)
      bw = 1. # We do not need to divide the bin width for input
      h1new.SetBinContent(i+1, ly[i]/bw)
      h1new.SetBinError(i+1, ldy[i]/bw)
   #print i, lx[i], lx[i+1], ly[i]
   print 'ysum =', ysum
   return h1new

def makewsinput(rootdir, filename, var, cut, binning = [], flagbin = 0, ntrees=100):
   list_chains = []
   for chan in channels:
      for samp in samples:
         systs = [ 'nominal', 'tes',  'tauid','met','leppt','lepid','phopt','phoid',]
         if 'QBDT1' in filename:
            systs = ['nominal', 'tes']
         if 'QBDT3' in filename:
            systs = ['nominal', 'tes', 'tauid', 'met']
         for syst in systs:
            if 'official' in rootdir:
               if samp == 'sig':
                  file0 = TFile(rootdir+'/TMVApp_'+'h2atata.root', 'read')
                  h1 = file0.Get('MVA_BDTG')
                  h1.Scale(2.)
                  #h1.Scale(10.)
               elif samp == 'data':
                  if syst != 'nominal':
                     continue
                  hlumi = TH1F("lumiininvpb", "", 1, 0, 1)
                  hlumi.SetBinContent(1, 80000.)
                  add_hist(hlumi, filename, chan, '', 'lumiininvpb')
                  file0 = TFile(rootdir+'/TMVApp_'+'atata.root', 'read')
                  h1 = file0.Get('MVA_BDTG')
                  h1.Scale(2.)
               elif samp == 'bkg':
                  if syst == 'nominal':
                     file0 = TFile(rootdir+'/TMVApp_'+'atata.root', 'read')
                  else:
                     file0 = TFile(rootdir+'/TMVApp_'+'atata_syst_'+syst+'.root', 'read')
                  h1 = file0.Get('MVA_BDTG')
                  h1.Scale(2.)
               print h1, h1.Integral()
               if binning == [] or not flagbin:
                  add_hist(h1, filename, chan, samp, syst)
               else:
                  h1new = binhist(h1, binning, 'h_'+chan+'_'+samp+'_'+syst) 
                  add_hist(h1new, filename, chan, samp, syst)
               continue
                  
            chain = get_chain(rootdir, chan, samp, syst, ntrees)
            if chain not in list_chains:
               list_chains.append(chain)
            if samp == 'data' and syst !='nominal':
               continue
            if binning == [] or not flagbin:
               #h1 = get_hist(chain, var, cut, '(20, -1., 1.)')
               #h1 = get_hist(chain, var, cut, '(10, -1., 1.)')
               remove_zero = 0
               #if samp == 'bkg':
                  #remove_zero = 1
               h1 = get_hist(chain, var, cut, '(20, -1., 1.)', remove_zero)
            else:
               h1 = get_hist(chain, var, cut, '(100, -1., 1.)')
            h1.Scale(2.)
            #if samp == 'sig':
               #h1.Scale(10.0)
            if samp == 'data':
               hlumi = TH1F("lumiininvpb", "", 1, 0, 1)
               hlumi.SetBinContent(1, 80000.)
               add_hist(hlumi, filename, chan, '', 'lumiininvpb')
            if binning == [] or not flagbin:
               add_hist(h1, filename, chan, samp, syst)
            else:
               print 'before new binning', h1.Integral()
               h1new = binhist(h1, binning, 'h_'+chan+'_'+samp+'_'+syst) 
               print 'after new binning', h1new.Integral()
               add_hist(h1new, filename, chan, samp, syst)
   return


binBDTG = [1, 21, 40, 52, 64, 69, 74, 75, 76, 77, 82, 85, 88, 90, 92, 93, 94, 96, 97, 98, 100]
bin3    = [1, 23, 36, 48, 56, 58, 61, 66, 69, 70, 72, 74, 75, 77, 79, 80, 85, 94, 95, 99, 100]
bin7    = [1, 26, 37, 45, 49, 57, 61, 65, 70, 75, 80, 81, 83, 84, 86, 87, 88, 91, 93, 94, 100]
bin11   = [1, 36, 48, 56, 64, 66, 67, 71, 77, 78, 79, 80, 81, 86, 87, 88, 90, 93, 95, 96, 100]


rd = [
      'official_BDTG',
      'official_BDTGp',
      'trees0',
      'trees1',
      'trees3',
      'trees7',
      ]
fn = [
      'wsinput_BDTG.root',
      'wsinput_BDTGp.root',
      'wsinput_BDTGpp.root',
      'wsinput_QBDT0.root',
      'wsinput_QBDT1.root',
      'wsinput_QBDT3.root',
      'wsinput_QBDT7.root',
      'wsinput_QBDT0x.root',
      'wsinput_QBDT1x.root',
      'wsinput_QBDT3x.root',
      'wsinput_QBDT7x.root',
      'wsinput_QBDT0xx.root',
      'wsinput_QBDT1xx.root',
      'wsinput_QBDT3xx.root',
      'wsinput_QBDT7xx.root',
      ]
binnings = [ 
      binBDTG,
      bin3,
      bin7,
      bin11
      ]

for i in range(len(fn)):
   #if i not in [6,7,8,9,10,11,12,13]:
   #if i not in [16]:
   #if i <2:
      #continue
   filename = fn[i]
   #rootdir = rd[i]
   #if 'x' not in filename and 'BDTGpp' not in filename:
   if '7' not in filename:
      continue
   wstag = filename.replace('wsinput_', '')
   wstag = wstag.replace('.root', '')
   if 'BDTG' in filename:
      rootdir = 'official_'+wstag
   else:
      for a in [ 0, 1, 3, 7]:
         if str(a) in wstag:
            rootdir = 'trees'+str(a)
            break
   #if '3p' not in filename and '7p' not in filename:
      #continue
   ntrees = 100
   if 'xxx' in filename:
      ntrees = 250
   elif 'xx' in filename:
      ntrees = 200
   elif 'x' in filename:
      ntrees = 150
   newbin = binnings[0]
   c = '0.25'
   if '7' in filename:
      if ntrees == 100:
         c = '0.35'
      elif ntrees == 150:
         c = '0.30'
      elif ntrees == 200:
         c = '0.25'
   var = 'tanh('+c+'*qbdt)'
   cut = '(trainflag != 1)*(fweight)'
   print '******************************************making input for', filename
   prepare_dirs(filename)
   makewsinput(rootdir, filename, var, cut, newbin, 0, ntrees)
