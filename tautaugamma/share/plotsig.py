


import string
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
   Q = 2*((s+b)*math.log(1+s/b)-s)
   return Q

def plotQ():
   h2 = TH2F('h2', '', 100,0,1,100,0,1)
   for i in range(h2.GetNbinsX()):
      for j in range(h2.GetYaxis().GetNbins()):
         b = h2.GetXaxis().GetBinCenter(i+1)
         s = h2.GetYaxis().GetBinCenter(j+1)
         Q = calQ(s,b)
         h2.SetBinContent(i+1,j+1,Q)
         print i,j,b,s,Q
   Cs = TCanvas('Cs', '', 600, 600)
   Cs.SetRightMargin(0.2)
   h2.Draw('colz')
   h2.GetXaxis().SetTitle('b')
   h2.GetYaxis().SetTitle('s')
   h2.GetZaxis().SetTitle('Q')
   Cs.SaveAs('Q_s_b.png')
   return


plotQ()
