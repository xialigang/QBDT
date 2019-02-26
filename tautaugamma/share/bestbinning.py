
import os
from ROOT import *
import ROOT
import math
from array import array

gROOT.LoadMacro("AtlasStyle.C")
from ROOT import SetAtlasStyle
SetAtlasStyle()

gROOT.SetBatch(True)


def calQ(s,b):
   if b==0:
      return -1
   Q = 2*((s+b)*math.log(1+s/b)-s)
   return Q


def get_bestsplit(hs, hb, ia, ib):
   nbins = hs.GetNbinsX()
   cut = -9999
   dQmax = -9999
   s = hs.Integral(ia,ib)
   b = hb.Integral(ia,ib)
   Q = calQ(s,b)
   Ql = -9999
   Qr = -9999
   #print 'ia =', ia, type(ia)
   #print 'ib =', ib, type(ib)
   for i in range(ia, ib):
      sL = hs.Integral(ia,i)
      bL = hb.Integral(ia,i)
      sR = s - sL
      bR = b - bL
      if bL == 0:
         continue
      if bR == 0:
         break
      QL = calQ(sL, bL)
      QR = calQ(sR, bR)
      dQ = QL + QR - Q
      if dQ > dQmax:
         dQmax = dQ
         cut = int(i)
         Ql = QL
         Qr = QR
   return cut, dQmax

def find_bestsplit(hs, hb, binning = [], precuts = [], flag_stop = 0, Nbins = 20):
   nbins = hs.GetNbinsX()
   ia = 1
   ib = nbins
   if len(binning)==0 or len(precuts) == 0:
      ia = 1
      ib = nbins
      binning.append(ia)
      binning.append(ib)
   if len(precuts) != 0:
      Nprecuts = len(precuts)
      dQmax = -9999
      for i in range(Nprecuts):
         node = precuts[i]
         dQ = node[-1]
         if dQ > dQmax:
            dQmax = dQ
            ia = node[0]
            ib = node[1]
   #print 'ia =', ia, type(ia)
   #print 'ib =', ib, type(ib)
   cut, dQmax = get_bestsplit(hs,hb,ia,ib)
   binning.append(cut)
   binning.sort()
   #print len(binning)-1, 'binning =', binning
   if len(binning) - 1 >= Nbins:
      flag_stop = 1
   node = [ia, ib, cut, dQmax]
   if node in precuts:
      precuts.remove(node)
   cutL, dQmaxL = get_bestsplit(hs, hb, ia, cut)
   cutR, dQmaxR = get_bestsplit(hs, hb, cut+1, ib)
   if flag_stop:
      return
   else:
      precuts.append([ia, cut, cutL, dQmaxL])
      precuts.append([cut+1, ib, cutR, dQmaxR])
      find_bestsplit(hs, hb, binning, precuts, flag_stop, Nbins)

def get_hist(chain, var, cut, histopt='(40, -1., 1.)', histname= 'htmp'):
   #print chain
   chain.Draw(var+'>>htmp'+histopt, cut)
   hist = gROOT.FindObject('htmp')
   hist.SetName(histname)
   #print hist, hist.Integral()
   return hist

def plotBDT(hsig, hbkg, plotname):
   #hsig.Scale(1./hsig.Integral())
   #hbkg.Scale(1./hbkg.Integral())
   Cs = TCanvas('Cs', '', 600, 600)
   hsig.Draw('PE')
   hsig.SetMarkerColor(ROOT.kBlue)
   hsig.SetLineColor(ROOT.kBlue)
   hbkg.Draw('PEsame')
   hbkg.SetMarkerColor(ROOT.kRed)
   hbkg.SetLineColor(ROOT.kRed)
   ymax = hsig.GetMaximum()
   if hbkg.GetMaximum() > ymax:
      ymax = hbkg.GetMaximum()
   hsig.GetYaxis().SetRangeUser(0, ymax*1.5)
   leg = TLegend(0.6, 0.75, 0.95, 0.9)
   leg.AddEntry(hsig, 'sig.', 'LPE')
   leg.AddEntry(hbkg, 'bkg.', 'LPE')
   leg.SetFillStyle(0)
   leg.SetBorderSize(0)
   leg.Draw()
   Cs.SaveAs('Cs_'+plotname+'.png')
   Cs.SaveAs('Cs_'+plotname+'.pdf')

def binhist(h1, binning, histname):
   h1.Scale(1./h1.Integral())
   ly = []
   ldy = []
   lx = array('f')
   nbins = len(binning) - 1
   print 'nbins =', nbins
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
   h1new = TH1F(histname, '', nbins, lx)
   ysum = 0.
   for i in range(nbins):
      ysum += ly[i]
      bw = h1new.GetBinWidth(i+1)
      h1new.SetBinContent(i+1, ly[i]/bw)
      h1new.SetBinError(i+1, ldy[i]/bw)
      #print i, lx[i], lx[i+1], ly[i]
   print 'ysum =', ysum
   return h1new
def work(treename, sigfilepath, bkgfilepath, var, cut, histopt='(100, -1, 1)', Nbins = 20, cuts = '', cutb = ''):
   if not os.path.isfile(sigfilepath):
      print sigfilepath, 'not found!!!'
      return
   else:
      print 'signal file path =', sigfilepath
   if not os.path.isfile(bkgfilepath):
      print bkgfilepath, 'not found!!!'
      return
   else:
      print 'background file path =', bkgfilepath
   print 'treename =', treename
   if cuts != '' or cutb != '':
      file0 = TFile(sigfilepath, 'read')
      print file0
      dir0 = file0.Get('dataset')
      sigchain = dir0.Get(treename)
      bkgchain = sigchain
   else:
      sigchain = TChain(treename)
      sigchain.Add(sigfilepath)
      bkgchain = TChain(treename)
      bkgchain.Add(bkgfilepath)
      cuts = cut
      cuts = cut

   hsig = get_hist(sigchain, var, cuts, histopt, 'hsig')
   hbkg = get_hist(bkgchain, var, cutb, histopt, 'hbkg')
   print hsig, hsig.Integral()
   print hbkg, hbkg.Integral()
   hsig.Scale(1./hsig.Integral())
   hbkg.Scale(1./hbkg.Integral())
   print hsig, hsig.Integral()
   print hbkg, hbkg.Integral()
   plotBDT(hsig, hbkg, 'BDT0_old')
   binning = []
   #def find_bestsplit(hs, hb, binning = [], precuts = [], flag_stop = 0, Nbins = 20):
   find_bestsplit(hsig, hbkg, binning, [], 0, Nbins)
   print 'binning =', binning
   hsignew = binhist(hsig, binning, 'hsignew')
   hbkgnew = binhist(hbkg, binning, 'hbkgnew')
   nbins = hsignew.GetNbinsX()
   print 'new nbins =', nbins
   #for i in range(nbins):
      #print i, hsignew.GetBinContent(i+1), hbkgnew.GetBinContent(i+1)
   plotBDT(hsignew, hbkgnew, 'BDT0_new')
   #plotBDT(hsignew, hbkgnew, 'BDT0_new1')
   return


treedir = 'trees3'
#treedir = 'trees7'
#treedir = 'trees11'
treedir = 'official_BDTG'
treedir = 'trees2'
treedir = 'trees0'
treedir = 'trees1'
treedir = 'trees2'
treedir = 'trees4'
if 'trees' in treedir:
   treename = 'bdt'
   sigfilepath = treedir +'/' + 'bdt_sig.root'
   bkgfilepath = treedir +'/' + 'bdt_bkg.root'
   var = 'tanh(0.25*qbdt)'
   #var = 'tanh(0.2*qbdt)'
   cut = '(trainflag !=1)*(fweight)'
   cuts = ''
   cutb = ''
else:
   treename = 'TestTree' 
   sigfilepath = treedir +'/' + 'TMVA.root'
   bkgfilepath = treedir +'/' + 'TMVA.root'
   #var = 'BDTG'
   var = 'tanh(0.6*atanh(BDTG))'
   cut = ''
   cuts = '(classID == 0)*(fweight)'
   cutb = '(classID == 1)*(fweight)'

histopt = '(100, -1, 1)'
Nbins = 20
work(treename, sigfilepath, bkgfilepath, var, cut, histopt, Nbins, cuts, cutb)





   
   


