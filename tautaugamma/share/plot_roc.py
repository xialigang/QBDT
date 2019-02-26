import os.path
import copy
import pickle
import logging
import ROOT
import math
import sys, getopt
import string
import os
from ROOT import *

import ROOT

from array import array

gROOT.LoadMacro("AtlasStyle.C")
from ROOT import SetAtlasStyle
SetAtlasStyle()

gROOT.SetBatch(True)


def getoverflow(hist):
   nbins = hist.GetNbinsX()
   a = hist.GetBinContent(nbins)
   da = hist.GetBinError(nbins)
   a1 = hist.GetBinContent(nbins+1)
   da1 = hist.GetBinError(nbins+1)
   hist.SetBinContent(nbins, a+a1)
   hist.SetBinError(nbins, math.sqrt(da*da+da1*da1))
   a = hist.GetBinContent(1)
   da = hist.GetBinError(1)
   a1 = hist.GetBinContent(0)
   da1 = hist.GetBinError(0)
   hist.SetBinContent(1, a+a1)
   hist.SetBinError(1, math.sqrt(da*da+da1*da1))
   return hist

def get_roc_graph(hq_sig, hq_bkg):
   sig_acc = array('d')
   bkg_sup = array('d')
   nbins = hq_sig.GetNbinsX()
   #print 'nbins =',nbins
   for i in range(nbins):
      s = hq_sig.Integral(1,i+1)/hq_sig.Integral(1,nbins)
      b = hq_bkg.Integral(1,i+1)/hq_bkg.Integral(1,nbins)
      #print i+1, s, b
      sig_acc.append(1-s)
      bkg_sup.append(b)
   #for i in range(len(sig_acc)):
      #print i, sig_acc[i], bkg_sup[i]
   g_roc = TGraph(len(sig_acc), sig_acc, bkg_sup)
   return g_roc

def get_rocs(filepath, ntrees=100):
   histopt = '(110, -1.1, 1.1)'
   tmpsigchain=TChain('bdt')
   tmpsigchain.Add(filepath+'/bdt_sig_'+str(ntrees)+'.root')
   var = 'qbdt'
   var = '2./(1.+exp(-1.*'+var+')) -1.'
   tmpsigchain.Draw(var+'>>hq_sig_train'+histopt,'trainflag==1')
   tmpsigchain.Draw(var+'>>hq_sig_test'+histopt,'trainflag!=1')
   hq_sig_train = gROOT.FindObject('hq_sig_train')
   hq_sig_test = gROOT.FindObject('hq_sig_test')
   hq_sig_train = getoverflow(hq_sig_train)
   hq_sig_test = getoverflow(hq_sig_test)
   #for i in range(hq_sig_train.GetNbinsX()):
      #print i, hq_sig_train.GetBinContent(i)
   tmpbkgchain=TChain('bdt')
   tmpbkgchain.Add(filepath+'/bdt_bkg_'+str(ntrees)+'.root')
   tmpbkgchain.Draw(var+'>>hq_bkg_train'+histopt,'trainflag==1')
   tmpbkgchain.Draw(var+'>>hq_bkg_test'+histopt,'trainflag!=1')
   hq_bkg_train = gROOT.FindObject('hq_bkg_train')
   hq_bkg_test = gROOT.FindObject('hq_bkg_test')
   hq_bkg_train = getoverflow(hq_bkg_train)
   hq_bkg_test = getoverflow(hq_bkg_test)
   #for i in range(hq_bkg_train.GetNbinsX()):
      #print i, hq_bkg_train.GetBinContent(i)
   g_roc_train = get_roc_graph(hq_sig_train, hq_bkg_train)
   g_roc_test = get_roc_graph(hq_sig_test, hq_bkg_test)
   return g_roc_train, g_roc_test
def get_rocs_official(filepath):
   histopt = '(110, -1.1, 1.1)'
   file0 = TFile(filepath+'/TMVA.root', 'read')
   TrainTree = file0.Get('dataset/TrainTree')
   TestTree = file0.Get('dataset/TestTree')
   TrainTree.Draw('tanh(0.65*atanh(BDTG))>>hq_sig_train'+histopt,'classID==0')
   TestTree.Draw('tanh(0.65*atanh(BDTG))>>hq_sig_test'+histopt,'classID==0')
   hq_sig_train = gROOT.FindObject('hq_sig_train')
   hq_sig_test = gROOT.FindObject('hq_sig_test')
   hq_sig_train = getoverflow(hq_sig_train)
   hq_sig_test = getoverflow(hq_sig_test)

   TrainTree.Draw('tanh(0.65*atanh(BDTG))>>hq_bkg_train'+histopt,'classID==1')
   TestTree.Draw('tanh(0.65*atanh(BDTG))>>hq_bkg_test'+histopt,'classID==1')
   hq_bkg_train = gROOT.FindObject('hq_bkg_train')
   hq_bkg_test = gROOT.FindObject('hq_bkg_test')
   hq_bkg_train = getoverflow(hq_bkg_train)
   hq_bkg_test = getoverflow(hq_bkg_test)
   
   g_roc_train = get_roc_graph(hq_sig_train, hq_bkg_train)
   g_roc_test = get_roc_graph(hq_sig_test, hq_bkg_test)
   return g_roc_train, g_roc_test

def get_graph_ratio(g1, g0, showflag = 0):
   n = 20
   dx = 1.0/n
   ax = array('d')
   ay = array('d')
   #print g0.GetN()
   for i in range(n):
      x = (i+0.5)*dx
      y0 = g0.Eval(x)
      y1 = g1.Eval(x)
      if y0 > 0:
         y = y1/y0
      else:
         y = 1.
      y = y1 - y0
      #if y0>0:
      #   y = (y1-y0)/y0
      #else:
      #   y = 0.
      ax.append(x)
      ay.append(y)
   g = TGraph(n, ax, ay)
   return g

def plot_rocs(list_filepath,a, highscore = 0, plotnametag =''):
   list_g_roc = []
   list_rg_roc= []
   for fp in list_filepath:
      if 'official_BDTG' in fp:
         g0,g1 = get_rocs_official(fp)
      else:
         ntrees = 100
         if 'trees7' in fp:
         #if 'trees7' in fp:
            ntrees = 200
         g0,g1 = get_rocs(fp,ntrees)
      #print 'checking ', fp
      #for i in range(g1.GetN()):
         #print fp, i, g1.GetX()[i], g1.GetY()[i]
      list_g_roc.append([g0, g1])
   g0 = list_g_roc[0][0]
   g1 = list_g_roc[0][1]
   for g_roc in list_g_roc:
      rg0 = get_graph_ratio(g_roc[0], g0)
      rg1 = get_graph_ratio(g_roc[1], g1)
      list_rg_roc.append([rg0, rg1])
   Cs_roc = TCanvas('Cs_roc', '', 600, 600)
   padsize = 0.5
   tsize = 0.06
   pad1 = TPad('pad1', '', 0.0, padsize, 1.0, 1.0)
   pad2 = TPad('pad2', '', 0.0, 0.0, 1.0, padsize)
   pad1.Draw()
   pad2.Draw()
   pad1.cd()
   pad1.SetGrid()
   pad1.Modified()
   pad1.SetBottomMargin(0)
   xmin = 0
   xmax = 1
   ymin = 1e-3
   ymax = 1
   if highscore == 1:
      xmin = 0.8
      ymax = 0.9
   elif highscore == 2:
      xmax = 0.8
      ymin = 0.8
   h2 = TH2F('h2', '', 100, xmin, xmax, 100, ymin, ymax)
   h2.Draw()
   i = 0
   for g_roc in list_g_roc:
      #g_roc[0].Draw('Lsame')
      g_roc[0].SetLineStyle(ROOT.kDashed)
      g_roc[0].SetLineWidth(2)
      g_roc[0].SetLineColor(i+1)
      g_roc[1].Draw('Lsame')
      g_roc[1].SetLineWidth(2)
      g_roc[1].SetLineStyle(i+1)
      if i>=2:
         #g_roc[1].SetLineColor(ROOT.kBlue)
         g_roc[1].SetLineColor(i)
      i += 1
   #h2.GetXaxis().SetTitle('Signal acceptance')
   h2.GetXaxis().SetNdivisions(505)
   h2.GetYaxis().SetTitle('Background rejection')
   h2.GetYaxis().SetTitleSize(0.8*tsize/(1-padsize))
   h2.GetYaxis().SetTitleOffset(0.7)
   h2.GetYaxis().SetLabelSize(0.6*tsize/(1-padsize))
   h2.GetYaxis().SetNdivisions(505)
   leg = TLegend(0.2, 0.1, 0.6, 0.6)
   Ncase = len(list_g_roc)
   for i in range(Ncase):
      #leg.AddEntry(list_g_roc[i][0], a[i]+' (train)', 'L')
      leg.AddEntry(list_g_roc[i][1], a[i]+' (test)', 'L')
   #leg.SetFillStyle(0)
   leg.SetBorderSize(0)
   leg.Draw()
   pad2.cd()
   pad2.SetTopMargin(0)
   pad2.SetBottomMargin(0.4)
   #rh2 = TH2F('rh2', '', 100, xmin, xmax, 100, 0.95, 1.05)
   rh2 = TH2F('rh2', '', 100, xmin, xmax, 100, -0.03, 0.03)
   
   rh2.Draw()
   i = 0
   for g_roc in list_rg_roc:
      #g_roc[0].Draw('Lsame')
      g_roc[0].SetLineStyle(ROOT.kDashed)
      g_roc[0].SetLineWidth(2)
      g_roc[0].SetLineColor(i+1)
      g_roc[1].Draw('Lsame')
      g_roc[1].SetLineWidth(2)
      g_roc[1].SetLineStyle(i+1)
      if i>=2:
         g_roc[1].SetLineColor(i)
         #g_roc[1].SetLineColor(ROOT.kBlue)
      #g_roc[1].SetLineColor(i+1)
      i += 1
   rh2.GetYaxis().SetTitle('Difference')
   rh2.GetYaxis().SetTitleSize(0.8*tsize/padsize)
   rh2.GetYaxis().SetTitleOffset(0.7)
   rh2.GetYaxis().SetLabelSize(0.6*tsize/padsize)
   rh2.GetYaxis().SetNdivisions(505)
   rh2.GetXaxis().SetTitle('Signal acceptance')
   rh2.GetXaxis().SetTitleOffset(1)
   rh2.GetXaxis().SetNdivisions(505)
   rh2.GetXaxis().SetTitleSize(tsize/padsize)
   rh2.GetXaxis().SetLabelSize(0.6*tsize/padsize)

   plotname = 'Cs_rocs_'+str(highscore)
   if plotnametag != '':
      plotname += '_'+plotnametag
   Cs_roc.SaveAs(plotname +'.png')
   Cs_roc.SaveAs(plotname +'.pdf')
   return

def plot_rocs0(list_filepath,a, highscore = 0, plotnametag =''):
   list_g_roc = []
   for fp in list_filepath:
      if 'official_BDTG' in fp:
         g0,g1 = get_rocs_official(fp)
      else:
         ntrees = 100
         if 'trees7' in fp:
            ntrees = 200
         g0,g1 = get_rocs(fp,ntrees)
      list_g_roc.append([g0, g1])
   Cs_roc = TCanvas('Cs_roc', '', 600, 600)
   Cs_roc.SetGrid()
   Cs_roc.Modified()
   xmin = 0
   xmax = 1
   ymin = 0
   ymax = 1
   if highscore == 1:
      xmin = 0.8
      ymax = 0.9
   elif highscore == 2:
      xmax = 0.8
      ymin = 0.8
   h2 = TH2F('h2', '', 100, xmin, xmax, 100, ymin, ymax)
   
   h2.Draw()
   i = 0
   for g_roc in list_g_roc:
      #g_roc[0].Draw('Lsame')
      g_roc[0].SetLineStyle(ROOT.kDashed)
      g_roc[0].SetLineWidth(2)
      g_roc[0].SetLineColor(i+1)
      g_roc[1].Draw('Lsame')
      g_roc[1].SetLineWidth(2)
      g_roc[1].SetLineColor(i+1)
      i += 1
   h2.GetXaxis().SetTitle('Signal acceptance')
   h2.GetXaxis().SetNdivisions(505)
   h2.GetYaxis().SetTitle('Background rejection')
   h2.GetYaxis().SetNdivisions(505)
   leg = TLegend(0.2, 0.2, 0.7, 0.6)
   Ncase = len(list_g_roc)
   for i in range(Ncase):
      #leg.AddEntry(list_g_roc[i][0], a[i]+' (train)', 'L')
      leg.AddEntry(list_g_roc[i][1], a[i]+' (test)', 'L')
   leg.SetFillStyle(0)
   leg.SetBorderSize(0)
   leg.Draw()
   plotname = 'Cs_rocs_'+str(highscore)
   if plotnametag != '':
      plotname += '_'+plotnametag
   Cs_roc.SaveAs(plotname +'.png')
   Cs_roc.SaveAs(plotname +'.pdf')
   return

def func(x):
   num = math.tanh(x/2.)
   den = (1+math.exp(x))*(1+math.exp(-x))
   return num/den

def showfunc():
   for i in range(100):
      x = -10 + 20./100*i
      print x, func(x)

#showfunc()
tag = 'onesyst'
#tag = 'threesysts'
#tag = 'sevensysts'
tag = 'all'

if 'one' in tag:
   a = ['GradBDT','QBDT0','QBDT1','QBDT1\'',]
   dir_filepath = ['official_BDTG', 'trees0',  'trees1', 'trees7']
elif 'three' in tag:
   a = ['GradBDT','QBDT0','QBDT3','QBDT3\'',]
   dir_filepath = ['official_BDTG', 'trees0',  'trees3', 'trees4']
elif 'seven' in tag:
   a = ['GradBDT','QBDT0','QBDT7','QBDT7\'',]
   dir_filepath = ['official_BDTG', 'trees0',  'trees5', 'trees6']

if 'one' in tag:
   a = ['GradBDT','QBDT0','QBDT1']
   dir_filepath = ['official_BDTGp', 'trees0',  'trees1']
if 'all' in tag:
   a = ['GradBDT','QBDT0','QBDT1', 'QBDT3', 'QBDT7']
   dir_filepath = ['official_BDTG', 'trees0',  'trees1', 'trees3', 'trees7']


plot_rocs(dir_filepath,  a, 0, tag)
#plot_rocs(dir_filepath,  a, 1, tag)
#plot_rocs(dir_filepath,  a, 2, tag)
