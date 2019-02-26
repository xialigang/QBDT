from ROOT import *
import os
import os.path
import copy
import pickle
import logging
import ROOT
from array import array
import shutil
import math

class qbdtmodule(object):

    def __init__(self):
        self.plotformats = ['png', 'pdf', 'root']
        self.signals = {}
        self.backgrounds = {}
        self.trainingfrac = 0.5
        self.weightvar = ''
        self.trainflagvar = ''
        self.variables = {'var': ('varname', 'xtitle', 'cut', [100, 0, 1]) }
        self.order = {}
        self.bestsplit = {}
        self.lastbestsplitvar = ''
        self.signal_hists = {}
        self.background_hists = {}
        self.treename = 'nominal'
        self.tree_epsil = -1
        self.tree_alpha = -1
        self.initial_weight_s = 1.0
        self.initial_weight_b = 1.0
        self.sigchain = None
        self.bkgchain = None
        self.sigchainsysts = {}
        self.bkgchainsysts = {}
        self.sigsysts = {}
        self.bkgsysts = {}
        self.corrmatrix = []
        self.syston = 0
        self.depth = 0
        self.maxdepth = 2
        self.node = 0
        self.totNbkg = 1
        self.totNsig = 1
        self.maxnode = 6
        self.mindQ = 0.01
        self.mindQr = 0.1
        self.tree = []
        self.pretree = []
        self.maxtree = 2
        self.treeq = []
        self.alltreeq = []
        self.istart = -1
        self.itree = -1
        self.treedir = 'trees'
        self.lasttree = []
        self.oldtree = []
        self.oldtrees = []
        self.trainflagcut = 0.10 #0.33
    def load_files(self, test=0):
        self.test = test
        if self.treename == '':
            print 'treename is not specified!!!'
            return
        self.sigchain = TChain(self.treename)
        for sig in self.signals:
            print 'loading signal', self.signals[sig][0], 'file', self.signals[sig][1]
            if not os.path.isfile(self.signals[sig][1]):
                print self.signals[sig][1],'is not a file!!!'
                continue
            filename = self.signals[sig][1]
            print 'filename =', filename
            filename = filename[filename.find('/')+1:]
            print 'filename =', filename
            if not os.path.isfile(self.treedir+'/'+filename):
                shutil.copyfile(self.signals[sig][1], self.treedir+'/'+filename)
            if test:
                self.sigchain.Add(self.signals[sig][1])
            else:
                self.sigchain.Add(self.treedir+'/'+filename)
        print 'signal total entries :', self.sigchain.GetEntries()
        self.bkgchain = TChain(self.treename)
        for bkg in self.backgrounds:
            print 'loading background', self.backgrounds[bkg][0], 'file', self.backgrounds[bkg][1]
            if not os.path.isfile(self.backgrounds[bkg][1]):
                print self.backgrounds[bkg][1],'is not a file!!!'
                continue
            filename = self.backgrounds[bkg][1]
            filename = filename[filename.find('/')+1:]
            if not os.path.isfile(self.treedir+'/'+filename):
                shutil.copyfile(self.backgrounds[bkg][1], self.treedir+'/'+filename)
            if test:
                self.bkgchain.Add(self.backgrounds[bkg][1])
            else:
                self.bkgchain.Add(self.treedir+'/'+filename)
        print 'background total entries :', self.bkgchain.GetEntries()
        print 'now we load background systemaitc files'
        #NOTE: Here below I only consider one background. In the future, we may combine all backgrounds into one for convenience.
        for syst in self.bkgsysts:
            treename = self.bkgsysts[syst][1]
            sample = self.bkgsysts[syst][2]
            filename = self.bkgsysts[syst][3]
            filepath = self.bkgsysts[syst][4]
            bkgchain = TChain(treename)
            rootfilepath = filepath + '/' + filename
            print 'loading syst', syst, 'for background', sample,' with file', filepath,'/', filename
            if not os.path.isfile(rootfilepath):
                print rootfilepath, 'is not a file!!!'
                continue
            if not os.path.isfile(self.treedir+'/'+filename):
                shutil.copyfile(rootfilepath, self.treedir+'/'+filename)
            bkgchain.Add(self.treedir+'/'+filename)
            self.bkgchainsysts[syst] = bkgchain
            print 'systematics', syst, 'background total entries :', self.bkgchainsysts[syst].GetEntries()
        print 'load_files DONE!!!'
    def get_hist(self, s, varname, cut):
        chaintag = 'signal'
        histopt = self.variables[varname][3]
        #cut0 = varname +'>' + str(histopt[1]) + ' && '+varname+'<'+ str(histopt[2])
        cut0 = ''
        histopt0 = '('+str(histopt[0]) + ',' + str(histopt[1]) + ',' + str(histopt[2]) + ')'
        nbins = histopt[0]
        if cut == '':
            cut = self.weightvar + '*(1)'
        else:
            cut = self.weightvar + '*('+cut+')'
        if s:
            self.sigchain.Draw(varname+'>>htmp' + histopt0, cut)
        else:
            chaintag = 'background'
            self.bkgchain.Draw(varname+'>>htmp' + histopt0, cut)
        htmp = gROOT.FindObject('htmp')
        c_lastbin = htmp.GetBinContent(nbins)
        e_lastbin = htmp.GetBinError(nbins)
        c_overflow = htmp.GetBinContent(nbins+1)
        e_overflow = htmp.GetBinError(nbins+1)
        htmp.SetBinContent(nbins, c_lastbin + c_overflow)
        htmp.SetBinError(nbins, math.sqrt(pow(e_lastbin,2) + pow(e_overflow,2)))
        htmp.SetName('hist_'+varname+str(s))
        #print 'Get hist from', chaintag, 'for', varname, 'with the cut', cut
        return htmp
    def get_sigbkg_hist(self, tag, varname, cut, syst = 'nominal'):
        histopt = self.variables[varname][3]
        histopt0 = '('+str(histopt[0]) + ',' + str(histopt[1]) + ',' + str(histopt[2]) + ')'
        nbins = histopt[0]
        fullcut = ''
        Nevt = 0.0
        if cut == '':
            fullcut = self.weightvar + '*(1)'
        else:
            fullcut = self.weightvar + '*('+cut+')'
        #print 'last tree length is', len(self.lasttree)
        hsig = TH1F('h'+tag+'_'+syst+'_'+varname, '', histopt[0], histopt[1], histopt[2])
        hsig.Sumw2()
        itree = self.itree
        if tag == 'sig' and syst == 'nominal' :
            samplechain = self.sigchain
        elif tag == 'bkg'and syst == 'nominal':
            samplechain = self.bkgchain
        elif tag == 'bkg' and syst != 'nominal':
            samplechain = self.bkgchainsysts[syst]
        else:
            print tag,'with', syst, 'is not found!!!'
            return
        #print samplechain, samplechain.GetName()
        if len(self.lasttree) == 0:
            cut = self.weightvar + '*(1)'
            samplechain.Draw(varname+'>>htmp' + histopt0, fullcut)
            htmp = gROOT.FindObject('htmp')
            #print tag, varname, 'htmp=',htmp.Integral()
            hsig.Add(htmp)
        else:
            if cut == '':
                cut = '1'
            bdtweight = '1.'
            if itree>0:
                bdtweight = 'bdtweight'+str(itree-1)
            fullcut0= self.weightvar + '*('+cut +')'
            fullcut = self.weightvar + '*'+bdtweight+'*('+cut +')'
            samplechain.Draw(varname + '>>htmp0' + histopt0, fullcut0)
            samplechain.Draw(varname + '>>htmp' + histopt0, fullcut)
            htmp0 = gROOT.FindObject('htmp0')
            htmp = gROOT.FindObject('htmp')
            #print tag, varname, 'htmp0=',htmp0.Integral(), 'htmp=',htmp.Integral()
            #if htmp.Integral() > 0:
                #htmp.Scale(htmp0.Integral() / htmp.Integral())
            hsig.Add(htmp)
        hsig = self.getoverflow(hsig)
        if tag == 'sig':
            hsig.Scale(self.initial_weight_s)
        else:
            hsig.Scale(self.initial_weight_b)
        #print 'Get background hist for', varname, 'from last tree'
        #print tag, 'hist integral =', hsig.Integral()
        return hsig
    def get_hists(self):
        for var in self.variables:
            hsig = self.get_hist(1, self.variables[var][0], '')
            hbkg = self.get_hist(0, self.variables[var][0], '')
            self.signal_hists[var] = hsig
            self.background_hists[var] = hbkg
    def show_variables(self):
        self.get_hists()
        Nvar = len(self.signal_hists)
        for var in self.variables:
            varname = self.variables[var][0]
            log = 0
            norm1 = 0
            if 'log' in self.variables[var][4]:
                log = 1
            if 'norm1' in self.variables[var][4]:
                norm1 = 1
            Cs = TCanvas('Cs_'+varname , '', 600, 600)
            if log:
                Cs.SetLogy()
            hsig = self.signal_hists[var]
            hbkg = self.background_hists[var]
            if norm1:
                hsig.Scale(1./hsig.Integral())
                hbkg.Scale(1./hbkg.Integral())
            hsig.SetLineColor(ROOT.kBlue)
            hbkg.SetLineColor(ROOT.kRed)
            hsig.Draw('hist')
            hbkg.Draw('hist,same')
            hsig.GetXaxis().SetTitle(self.variables[var][1])
            hsig.GetXaxis().SetNdivisions(505)
            hsig.GetYaxis().SetTitle('Fraction / %i GeV' % hsig.GetBinWidth(1))
            hsig.GetYaxis().SetNdivisions(505)
            ymax = hbkg.GetMaximum()
            if ymax < hsig.GetMaximum():
                ymax = hsig.GetMaximum()
            if log:
                hsig.GetYaxis().SetRangeUser(1.0e-4, ymax*10)
            else:
                hsig.GetYaxis().SetRangeUser(0, ymax*1.5)
            leg = TLegend(0.5, 0.7, 0.95, 0.95)
            leg.AddEntry(hsig, 'Signal', 'L')
            leg.AddEntry(hbkg, 'Background', 'L')
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            leg.Draw()
            for a in self.plotformats:
                Cs.SaveAs(self.treedir+'/Cs_var_'+varname+'.' + a)
    def calQ(self, s, b, db):
        # this can help to reduce the overfit
        #b = b + 0.001
        #b = b + 0.01
        #b = b + 0.05
        if b == 0 and s == 0:
            return -2,-2
        elif b == 0 and s > 0:
            return -1,-1
        elif s == 0 and b > 0:
            return 0.0,0.0
        else:
            # no background uncertainty
            Q0 = 2.0*((s+b)*math.log(1.0+s/b) - s)
            # consider background uncertainty
            db2 = pow(db,2)
            if db > 0 and self.syston:
                Q0 = 2.*((s+b)*math.log((s+b)*(b+db2)/(b*b+(s+b)*db2))-pow(b/db,2)*math.log(1+s*db2/b/(b+db2)))
            q = 5.0
            c = math.exp(-Q0/q)
            Q1 = Q0*c
            Q1 = Q0
            return Q0,Q1


    def get_dN(self, hnom, list_hsyst, bin0 = -1, bin1 = -1, list_syst=[]):
        if not self.syston:
            return 0.
        Nsyst = len(list_hsyst)
        dN = 0.
        if bin0 <0 or bin1 < 0:
            bin0 = 1
            bin1 = hnom.GetNbinsX()
        nom = hnom.Integral(bin0, bin1)
        for i in range(Nsyst):
            systi = list_hsyst[i].Integral(bin0, bin1)
            Deltai = systi - nom
            #if nom == 0:
            #    print 'nom =',nom
            #if (nom + pow(Deltai,2)) == 0:
            #    print 'nom =',nom
            #    print 'Deltai =',Deltai
            if (nom + pow(Deltai,2)) > 0:
                dN += pow(Deltai,2)/(nom + pow(Deltai,2))
            for j in range(Nsyst):
                if i==j:
                    continue
                systj = list_hsyst[j].Integral(bin0, bin1)
                Deltaj = systj - nom
                if (nom + pow(Deltai,2)) > 0 and (nom + pow(Deltaj,2)) > 0:
                    dN -= pow(Deltai*Deltaj,2)/(nom + pow(Deltai,2))/(nom + pow(Deltaj,2))
        if dN >= 0:
            #print 'dN =', dN, '>0'
            dN = math.sqrt(dN*nom)
        else:
            #print 'dN =', dN, '<0'
            dN = 0.0
        return dN

    def get_bestsplit(self, var, cut, show_bins = 0):
        hsig = self.get_sigbkg_hist('sig', self.variables[var][0], cut, 'nominal')
        hbkg = self.get_sigbkg_hist('bkg', self.variables[var][0], cut, 'nominal')
        list_hbkg_syst = []
        list_syst = []
        for syst in self.bkgsysts:
            if not self.syston:
                continue
            hbkg_syst = self.get_sigbkg_hist('bkg', self.variables[var][0], cut, syst)
            hbkg_syst.SetName('hbkg_'+syst)
            list_hbkg_syst.append(hbkg_syst)
            list_syst.append(syst)
        dNb = self.get_dN(hbkg, list_hbkg_syst, -1, -1, list_syst)
        nbins = hsig.GetNbinsX()
        show_bins = 0
        if show_bins:
            print 'cut =', cut
            print 'var =', var
            for i in range(nbins):
                print 'bin', i+1, hsig.GetBinLowEdge(i+1), hsig.GetBinContent(i+1), hbkg.GetBinContent(i+1)
        Ns = hsig.Integral()
        Nb = hbkg.Integral()
        #print 'for cut',cut,'Ns =',Ns, 'Nb =', Nb
        if Ns==0 and Nb==0:
            print 'var =',var, 'for cut', cut,' no sig/bkg events!!!'
        purity = Ns/(Ns + Nb + pow(dNb,2))
        #Q0 = (Ns+Nb)*purity*(1-purity)
        Q0,Q1 = self.calQ(Ns, Nb, dNb)
        dQmax = -9999
        x0 = 0.0
        Ql = 0.0
        Qr = 0.0
        pl = -1
        pr = -1
        flag_no_events = 0
        for i in range(1,nbins-1):
            NsL = hsig.Integral(1,i)
            NbL = hbkg.Integral(1,i)
            dNbL = self.get_dN(hbkg, list_hbkg_syst, 1, i, list_syst)
            dNbR = self.get_dN(hbkg, list_hbkg_syst, i+1, nbins, list_syst)
            #print i, dNb, dNbL, dNbR
            QL0,QL1 = self.calQ(NsL, NbL, dNbL)
            QR0,QR1 = self.calQ(Ns-NsL, Nb-NbL, dNbR)
            #dQ =  - (Q0 - QL0 - QR0)
            dQ = - (Q1 - QL1 - QR1)
            if dQ > dQmax:
                x0 = hsig.GetBinLowEdge(i+1)
                dQmax = dQ
                nsl = NsL
                nbl = NbL
        #if dQmax == -9999:
            #print 'dQmax =', dQmax, ' !!!'
            #print 'purity =', purity, ' !!!'
        # if there is no events af
        self.bestsplit[var] = (x0, dQmax,  Q0, purity, Ns, Nb)
        #self.bestsplit[var] = (x0, dQmax,  Q1, purity, Ns, Nb)
        return
    def get_bestsplit_var(self, cut, show_bins):
        dQmax = -9999
        self.bestsplit = {}
        var0 = ''
        for var in self.variables:
            self.get_bestsplit(var, cut, show_bins)
            #if self.bestsplit[var][1] < 0:
                #print 'var =',var,'dQ =', self.bestsplit[var][1]
            if self.bestsplit[var][1] > dQmax :
                dQmax = self.bestsplit[var][1]
                var0 = var
        #print 'The best variable with the cut', cut, 'is', var0, 'with dQmax =', dQmax
        if var0 == '':
            #print 'var0 =', var0, ' with cut =', cut
            return 'error'
        else:
            return var0
    def show_bestsplit(self):
        for var in self.variables:
            self.get_bestsplit(var, '')
            print var, self.bestsplit[var]
    def get_prenode(self, cut):
        var = self.get_bestsplit_var(cut, 0)
        if var == 'error':
            print 'var =', var,' !!!'
            #return ['', 'error', '', 0, -1, -1, -1, -1]
            for var0 in self.variables:
                var = var0
                break
        cutvalue = self.bestsplit[var][0]
        dQmax = self.bestsplit[var][1]
        Q = self.bestsplit[var][2]
        purity = self.bestsplit[var][3]
        Ns = self.bestsplit[var][4]
        Nb = self.bestsplit[var][5]
        prenode = [cut, var, cutvalue, dQmax, Q, purity, Ns, Nb]
        return prenode
    def check_prenode(self, prenode):
        cut = prenode[0]
        var = prenode[1]
        cutvalue = prenode[2]
        dQmax = prenode[3]
        Q = prenode[4]
        purity = prenode[5]
        Ns = prenode[6]
        Nb = prenode[7]
        #print 'checking prenode', cut,Ns, Nb
        samecut = 0
        if var+'<'+str(cutvalue) in cut:
            samecut = 1
        if var+'>'+str(cutvalue) in cut:
            samecut = 1
        #if prenode in self.pretree:
            #self.pretree.remove(prenode)
        #if dQmax < self.mindQ or Q<0 or Q==0 or self.node > self.maxnode or samecut:
        dQmaxcut = 0
        #print 'number of nodes =', self.node, self.maxnode
        if dQmax<0 or Q<=self.mindQ or Q<=0 or Ns == 0 or Nb == 0 or self.node > self.maxnode or samecut or Ns==0 or Nb==0:
            tguess = abs(Q)
            self.tree.append([self.depth, self.node, var, cut, Q, tguess, dQmax, purity, Ns, Nb])
            self.node += 1
            return 0
        else:
            return 1
    def tree_split(self, cut, iprenode, show_prenode=0, show_bins=0):
        Nprenode = len(self.pretree)
        if Nprenode == 0:
            prenode = self.get_prenode(cut)
            if prenode[1] == 'error':
                print 'prenode[1] = error'
                return
        else:
            if iprenode < 0:
                print 'ERROR for iprenode =', iprenode
                return
            prenode = self.pretree[iprenode]
        #print 'prenode =', prenode
        #print 'self.pretree =', self.pretree
        var = prenode[1]
        if var == 'error':
            #print 'var = error!!!'
            return
        cutvalue = prenode[2]
        dQmax = prenode[3]
        Q = prenode[4]
        purity = prenode[5]
        Ns = prenode[6]
        Nb = prenode[7]

        samecut = 0
        if var+'<'+str(cutvalue) in cut:
            samecut = 1
        if var+'>'+str(cutvalue) in cut:
            samecut = 1
        if prenode in self.pretree:
            self.pretree.remove(prenode)
        if 0:
            print 'self.node =', self.node
            print 'var =',var
            print 'dQmax =', dQmax
            print 'Q =', Q
            print 'purity =', purity
        #if dQmax < self.mindQ or Q<0 or Q==0 or self.node > self.maxnode or samecut:
        dQmaxcut = 0
        if dQmax < self.mindQ:
            dQmaxcut = 1
        #print 'number of nodes =', self.node, self.maxnode, self.node>self.maxnode
        #if dQmaxcut or Q<=0 or self.node > self.maxnode or samecut:
        if dQmax<0 or Q<=self.mindQ or Q<=0 or Ns == 0 or Nb == 0 or self.node > self.maxnode or samecut or Ns==0 or Nb==0:
            tguess = abs(Q)
            self.tree.append([self.depth, self.node, var, cut, Q, tguess, dQmax, purity, Ns, Nb])
            self.node += 1
        else:
            #print 'splitting ...', cut, var, 'cut =',cutvalue,'with', Q, dQmax, Ns, Nb
            prenode_left = self.get_prenode(cut+'&&'+var+'<'+str(cutvalue))
            if self.check_prenode(prenode_left):
                self.pretree.append(prenode_left)
            prenode_right = self.get_prenode(cut+'&&'+var+'>'+str(cutvalue))
            if self.check_prenode(prenode_right):
                self.pretree.append(prenode_right)
            self.node += 2
        #print 'to find the next split....'
        Nprenode = len(self.pretree)
        if Nprenode == 0:
            print 'No prenode, finish tree splitting...'
            return
        #print 'Nprenode =', Nprenode
        #print 'self.pretree =',self.pretree
        dQmax = -9999
        iprenode = -1
        cut = ''
        for i in range(Nprenode):
            #dQ = abs(self.pretree[i][3]/self.pretree[i][6])
            dQ = self.pretree[i][3]
            if dQ>dQmax:
                iprenode = i
                dQmax = dQ
                cut = self.pretree[i][0]
        self.depth += 1
        #print 'next split is', iprenode, self.pretree[iprenode]
        self.tree_split(cut,iprenode)


    def store_tree(self, i):
        if not os.path.isdir(self.treedir):
            os.mkdir(self.treedir)
        with open(self.treedir + '/tree_' + str(i) + '.txt', 'w') as f0:
            for node in self.tree:
                f0.write(str(node) + '\n')
    def turn_to_list(self,line0):
        line0 = line0.replace('[', '')
        line0 = line0.replace(']', '')
        line0 = line0.split(',')
        line = []
        for i in range(len(line0)):
            a = line0[i]
            if i==0 or i==1:
                a = int(a)
            elif i==4 or i==5:
                a = float(a)
            else:
                a = a.replace('\'','')
                a = a.strip()
            line.append(a)
        return line
    def cal_istart(self):
        if not os.path.isdir(self.treedir):
            return
        files = os.listdir(self.treedir)
        if len(files) == 0:
            return
        istart = -1
        i = 0
        while 1:
           if 'tree_'+str(i)+'.txt' in files:
              istart = i
              i += 1
           else:
              break
        print 'istart =', istart
        if istart < 0:
            return
        self.istart = istart
    def load_tree(self, i):
        self.lasttree = []
        if not os.path.isfile(self.treedir+'/tree_'+str(i)+'.txt'):
            return
        with open(self.treedir+'/tree_'+str(i)+'.txt') as f0:
            while 1:
                line0 = f0.readline()
                if line0 == '':
                    break
                line0 = line0.strip()
                line = self.turn_to_list(line0)
                print line
                self.lasttree.append(line)
    def build_trees(self):
        self.lasttree = []
        self.istart = -1
        self.cal_istart()
        self.itree = self.istart + 1
        self.oldtrees = []
        for i in range(self.istart + 1):
            self.load_tree(i)
            if self.lasttree in self.oldtrees:
                print 'i =', i
                os.remove(self.treedir+'/tree_'+str(i)+'.txt')
            else:
                self.oldtrees.append(self.lasttree)
        #return
        #print 'oldtrees =', self.oldtrees
        self.initial_weight_s = 1.0
        self.initial_weight_b = 1.0
        for i in range(self.maxtree):
            print '\n building tree', i
            if i <= self.istart:
                continue
            cut = '%s<%s' % (self.trainflagvar, str(self.trainflagcut))
            #cut = '1'
            if i>-1:
                self.sigchain = TChain(self.treename)
                self.sigchain.Add(self.treedir+'/fhist_h2atata.root')
                self.bkgchain = TChain(self.treename)
                self.bkgchain.Add(self.treedir+'/fhist_atata.root')
                print 'sigchain entries=', self.sigchain.GetEntries()
                print 'bkgchain entries=', self.bkgchain.GetEntries()
                self.bkgchainsyst = {}
                for syst in self.bkgsysts:
                    if not self.syston:
                        continue
                    # 'tes' : ('tes', 'NOMINAL', 'atata', 'fhist_atata_syst_tes.root', 'root_dir/', 'low'),
                    treename = self.bkgsysts[syst][1]
                    filename = self.bkgsysts[syst][3]
                    bkgchain = TChain(treename)
                    bkgchain.Add(self.treedir+'/'+filename)
                    self.bkgchainsysts[syst] = bkgchain
                    print 'bkgchain entries=', bkgchain.GetEntries(), 'for syst', syst
            if i == 0 or self.initial_weight_s == 1.0 or self.initial_weight_b == 1.0: #set initial weight
                var = ''
                for var in self.variables:
                    break
                print 'to get initial weight, the variables is', var
                hsig = self.get_sigbkg_hist('sig', self.variables[var][0], cut)
                hbkg = self.get_sigbkg_hist('bkg', self.variables[var][0], cut)
                #print 'N(hsig) =',hsig.Integral()
                #print 'N(hbkg) =',hbkg.Integral()
                self.initial_weight_s = 1./hsig.Integral()
                self.initial_weight_b = 1./hbkg.Integral()
            self.tree = []
            self.pretree = []
            self.depth = 0
            self.node = 0
            show_prenode = 0
            if i>0:
                show_prenode = 0
            self.tree_split(cut, -1, 1, 0)
            self.tree_epsil = -1
            self.tree_alpha = -1
            sum_Ns = 0.0
            sum_Nb = 0.0
            sum_Q = 0.0
            sum_Ws = 0.0
            sum_Wb = 0.0
            #self.tree.append([self.depth, self.node, var, cut, Q, tguess, dQmax, purity, Ns, Nb])
            for node in self.tree:
                Q = node[4]
                Ns = node[8]
                Nb = node[9]
                sum_Ns += Ns
                sum_Nb += Nb
                if Q < 0:
                    continue
                sum_Q += node[4]
            for node in self.tree:
                Q = node[4]
                purity = node[7]
                Ns = node[8]
                Nb = node[9]
                #if Q < 0:
                    #continue
                #print 'Q*sum_Q =',Q,'*',sum_Q,'=',Q*sum_Q
                #sum_Ws += Ns*math.exp(-Q*sum_Q)
                #sum_Wb += Nb*math.exp(Q*sum_Q)
                sum_Ws += Ns*math.exp(-purity*sum_Q)
                sum_Wb += Nb*math.exp(purity*sum_Q)
            #self.tree_epsil = (sum_Ns+sum_Nb)/2.0
            #if self.tree_epsil < 0 or self.tree_epsil > 1:
            #    for node in self.tree:
            #        print 'Ns, Nb =', node[8], node[9]
            #    print 'tree_epsil =', self.tree_epsil, ' < 0 or > 1 !!! stop ...'
            #    break
            #if self.tree_epsil >=0.5:
            #    print 'tree_epsil =', self.tree_epsil, ' >=0.5 !!!! stop ...'
            #    break
            #self.tree_alpha = 0.5*math.log((1-self.tree_epsil)/self.tree_epsil)
            #weight = math.exp(self.tree_alpha)
            #weight_s = 1.0/(sum_Ns*weight + (1.0-sum_Ns)*1.0)
            #weight_b = 1.0/(sum_Nb*weight + (1.0-sum_Nb)*1.0)
            self.tree_epsil = sum_Q
            self.tree_alpha = sum_Q
            weight_s = sum_Ns / sum_Ws
            weight_b = sum_Nb / sum_Wb
            for node in self.tree:
                node.append(self.tree_epsil)
                node.append(self.tree_alpha)
            self.lasttree = self.tree
            if self.lasttree in self.oldtrees:
                print 'SAME tree is built! stop...'
                break
            self.store_tree(i)
            self.itree += 1
            self.oldtrees.append(self.lasttree)
            if len(self.pretree)>0:
                print 'ERROR in splitting tree',i
                print 'Here is the self.pretree.'
                for a in self.pretree:
                    print a
                return
            self.sigchain = None
            self.bkgchain = None
            print 'now adding bdtweight as a branch to the signal ntuple'
            self.add_bdtweight2ntuple('fhist_h2atata.root', self.treename, 'sig', i, self.lasttree, sum_Q, weight_s, self.tree_alpha)
            print 'now adding bdtweight as a branch to the background ntuple'
            self.add_bdtweight2ntuple('fhist_atata.root', self.treename, 'bkg', i, self.lasttree, sum_Q, weight_b, self.tree_alpha)
            print 'now adding bdtweight as a branch to the background systematics ntuples'
            if not self.syston:
                continue
            for syst in self.bkgsysts:
                print 'syst =',syst
                treename = self.bkgsysts[syst][1]
                self.add_bdtweight2ntuple('fhist_atata_syst_'+syst+'.root', treename, 'bkg', i, self.lasttree, sum_Q, weight_b, self.tree_alpha)
        return
    def read_treeq(self, i):
        treefilepath = self.treedir+'/tree_'+str(i)+'.txt'
        self.treeq = []
        with open(treefilepath, 'r') as f0:
            while 1:
                line0 = f0.readline()
                if line0 == '':
                    break
                line0 = line0.strip()
                line0 = line0.replace('[', '')
                line0 = line0.replace(']', '')
                line0 = line0.split(',')
                self.treeq.append([line0[3].replace('\'',''), float(line0[7]), float(line0[-1])])
    def read_alltreeq(self):
        self.alltreeq = []
        for i in range(self.maxtree):
            self.read_treeq(i)
            #print 'read tree',i,self.treeq
            self.alltreeq.append(self.treeq)
    def get_q(self, event):
        q = 0.0
        qq = 1.0
        list_regions = []
        list_q = []
        for i in range(self.maxtree):
            inode = 0
            flag = 0
            for node in self.alltreeq[i]:
            #self.tree.append([self.depth, self.node, var, cut, Q, tguess, dQmax, purity, Ns, Nb])
                nodecond = node[0].replace('&&', '&&event.')
                #print 'type of nodecond is', type(nodecond)
                #print 'nodecond =', nodecond
                nodecond = nodecond.replace('&&', ' and ')
                nodecond = nodecond.replace('%s<%s and' % (self.trainflagvar, str(self.trainflagcut)), '')
                #print 'fallintonode = %s' % nodecond
                exec('fallintonode = %s' % nodecond)

                if fallintonode:
                    flag += 1
                    #if nodecond in list_regions:
                       #continue
                    list_regions.append(nodecond)
                    #print 'tree',i, 'node', inode, '=', node
                    purity = node[1]
                    alpha = node[2]
                    #q += purity*alpha
                    q += (2.*purity-1.)*alpha
                    break
                inode += 1
            if flag !=1:
                print 'tree',i,'is not complete with flag =', flag
        #qq = 1.0 - qq
        #q = q/self.maxtree
        #q = 2./(1.+math.exp(-2*q)) - 1.
        return q
    def get_q_hists_old(self, chain, tag, syst ='nominal'):
        nbins = 100
        xmin = -1.1
        xmax = 1.1
        hq_train = TH1F('hq_train', '', nbins, xmin, xmax)
        hq_test = TH1F('hq_test', '', nbins, xmin, xmax)
        ievent = 0
        bdtfilename = self.treedir+'/bdt_'+tag
        if syst != 'nominal':
            bdtfilename += '_'+syst
        bdtfilename += '_'+str(self.maxtree)+'.root'
        print 'We are producing test samples for', tag, syst
        file0 = TFile(bdtfilename, 'recreate')
        tree0 = TTree('bdt', '')
        trflag = array('f', [0])
        qbdt = array('f', [0])
        fweight = array('f', [0])
        tree0.Branch('trainflag', trflag, 'trainflag/F')
        tree0.Branch('qbdt', qbdt, 'qbdt/F')
        tree0.Branch('fweight', fweight, 'fweight/F')
        for event in chain:
            exec('trainflag = %s.%s<%s' % ('event', self.trainflagvar, str(self.trainflagcut)))
            exec('weight = %s.%s' %('event', self.weightvar))
            q = self.get_q(event)
            #print 'trainflag =', trainflag
            if ievent % 2000 == 0:
                print ievent, 'q =',q, 'weight =', weight
            #testflag = 1
            #if ievent > 8000 and testflag:
                #break
            trflag[0] = trainflag
            qbdt[0] = q
            fweight[0] = weight
            tree0.Fill()
            if trainflag:
                hq_train.Fill(q, weight)
            else:
                hq_test.Fill(q, weight)
            ievent += 1
        tree0.Write()
        file0.Close()
        hq_train.SetName('hq_train_'+tag)
        hq_test.SetName('hq_test_'+tag)
        hq_train = self.getoverflow(hq_train)
        hq_test = self.getoverflow(hq_test)
        hq_train.Scale(1.0 / hq_train.Integral())
        hq_test.Scale(1.0 / hq_test.Integral())
        return hq_train, hq_test
    def getoverflow(self, hist):
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
    def get_q_hists(self, chain, tag):
        nbins = 22
        nbins = 100
        xmin = -1.1
        xmax = 1.1
        hq_train = TH1F('hq_train', '', nbins, xmin, xmax)
        hq_test = TH1F('hq_test', '', nbins, xmin, xmax)
        #var = 'qbdt'+str(self.maxtree-1)+'/'+str(self.maxtree)
        var = '0.25*qbdt'+str(self.maxtree-1)
        #var = 'pow(qbdt'+str(self.maxtree-1)+', 1)'
        #var = '(1-exp(-qbdt'+str(self.maxtree-1)+'))'
        #var = '2./(1.+exp(-1.*'+var+')) -1.'
        var = 'tanh('+var+')'
        traincut = '%s<%s' % (self.trainflagvar, str(self.trainflagcut))
        chain.Draw(var+'>>htrain('+str(nbins)+','+str(xmin)+','+str(xmax)+')', 'fweight*('+traincut+')')
        chain.Draw(var+'>>htest('+str(nbins)+','+str(xmin)+','+str(xmax)+')', 'fweight*('+traincut+')')
        hq_train = gROOT.FindObject('htrain')
        hq_test = gROOT.FindObject('htest')
        print hq_train,hq_test
        hq_train = self.getoverflow(hq_train)
        hq_test = self.getoverflow(hq_test)
        hq_train.SetName('hq_train_'+tag)
        hq_test.SetName('hq_test_'+tag)
        hq_train.Scale(1.0 / hq_train.Integral())
        hq_test.Scale(1.0 / hq_test.Integral())
        return hq_train, hq_test
    def get_roc_graph(self, hq_sig, hq_bkg):
        sig_acc = array('d')
        bkg_sup = array('d')
        nbins = hq_sig.GetNbinsX()
        for i in range(nbins):
           s = hq_sig.Integral(1,i+1)/hq_sig.Integral()
           b = hq_bkg.Integral(1,i+1)/hq_bkg.Integral()
           #print i+1, s, b
           sig_acc.append(1-s)
           bkg_sup.append(b)
        g_roc = TGraph(len(sig_acc), sig_acc, bkg_sup)
        return g_roc
    def show_roc(self, hq_sig_train, hq_bkg_train, hq_sig_test, hq_bkg_test):
        print 'hq_sig', hq_sig_train.Integral(), hq_sig_test.Integral()
        print 'hq_bkg', hq_bkg_train.Integral(), hq_bkg_test.Integral()
        g_roc_test = self.get_roc_graph(hq_sig_test, hq_bkg_test)
        g_roc_train = self.get_roc_graph(hq_sig_train, hq_bkg_train)
        Cs_roc = TCanvas('Cs_roc', '', 600, 600)
        Cs_roc.SetGrid()
        Cs_roc.Modified()
        h2 = TH2F('h2', '', 100, 0, 1, 100, 0, 1)
        h2.Draw()
        g_roc_test.Draw('Lsame')
        g_roc_test.SetLineWidth(2)
        g_roc_test.SetLineColor(ROOT.kBlue)
        g_roc_train.Draw('Lsame')
        g_roc_train.SetLineStyle(ROOT.kDashed)
        g_roc_train.SetLineWidth(2)
        h2.GetXaxis().SetTitle('Sig. Acc.')
        h2.GetYaxis().SetTitle('Bkg. Red.')
        leg = TLegend(0.2, 0.3, 0.5, 0.5)
        leg.AddEntry(g_roc_train, 'ROC (train)', 'L')
        leg.AddEntry(g_roc_test, 'ROC (test)', 'L')
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.Draw()
        for a in self.plotformats:
            Cs_roc.SaveAs(self.treedir+'/Cs_roc_'+str(self.maxtree)+'.' + a)
    def show_q_hist(self, log, hq_sig_train, hq_bkg_train, hq_sig_test, hq_bkg_test):
        #if not os.path.isfile(filepath):
        #   print filepath,'is not found!'
        #   return
        #file0 = TFile(filepath, 'read')
        #hq_sig_train =  file0.Get('hq_train_sig')
        #hq_bkg_train =  file0.Get('hq_train_bkg')
        #hq_sig_test =  file0.Get('hq_test_sig')
        #hq_bkg_test =  file0.Get('hq_test_bkg')
        print 'hq_sig', hq_sig_train.Integral(), hq_sig_test.Integral()
        print 'hq_bkg', hq_bkg_train.Integral(), hq_bkg_test.Integral()
        Cs_q = TCanvas('Cs_q', '', 600, 600)
        if log:
            Cs_q.SetLogy()
        hq_sig_train.Draw('hist')
        hq_bkg_train.Draw('hist,same')
        hq_sig_test.Draw('EP, same')
        hq_bkg_test.Draw('EP,same')
        hq_sig_train.SetLineColor(ROOT.kBlue)
        hq_bkg_train.SetLineColor(ROOT.kRed)
        hq_sig_train.SetFillColor(ROOT.kBlue)
        hq_bkg_train.SetFillColor(ROOT.kRed)

        hq_sig_train.SetFillStyle(3001)
        hq_bkg_train.SetFillStyle(3004)

        hq_sig_test.SetLineColor(ROOT.kBlue)
        hq_bkg_test.SetLineColor(ROOT.kRed)
        hq_sig_test.SetMarkerColor(ROOT.kBlue)
        hq_bkg_test.SetMarkerColor(ROOT.kRed)
        ymax = hq_bkg_train.GetMaximum()
        if ymax < hq_sig_train.GetMaximum():
            ymax = hq_sig_train.GetMaximum()
        if log:
            hq_sig_train.GetYaxis().SetRangeUser(1.0e-6, ymax*10)
        else:
            hq_sig_train.GetYaxis().SetRangeUser(0, ymax*1.2)
        hq_sig_train.GetXaxis().SetTitle('q')
        leg = TLegend(0.6, 0.6, 0.95, 0.95)
        leg.AddEntry(hq_sig_train, 'sig. (train)', 'LF')
        leg.AddEntry(hq_bkg_train, 'bkg. (train)', 'LF')
        leg.AddEntry(hq_sig_test, 'sig. (test)', 'LPE')
        leg.AddEntry(hq_bkg_test, 'bkg. (test)', 'LPE')
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.Draw()
        Cs_q.RedrawAxis()
        for a in self.plotformats:
            Cs_q.SaveAs(self.treedir+'/Cs_q_'+str(log) +'_' +str(self.maxtree)+'.'+ a)
    def rearrange_q_hists(self, hsig, hbkg):
        nbins = hbkg.GetNbinsX()
        list_nbkg = []
        list_nsig = []
        list_nevt = []
        for i in range(nbins):
            nbkg = hbkg.GetBinContent(i+1)
            nsig = hsig.GetBinContent(i+1)
            list_nbkg.append(hbkg.GetBinContent(i+1))
            list_nsig.append(hsig.GetBinContent(i+1))
            ruler = -1
            if nbkg == 0:
                ruler = 9999
            else:
                #ruler = nsig/nbkg
                ruler = 1./nbkg
            list_nevt.append((ruler, nbkg,nsig))
        #print 'list_nevt =', list_nevt
        #list_nevt.sort(reverse=True)
        list_nevt.sort()
        #print 'list_nevt =', list_nevt
        for i in range(nbins):
           hbkg.SetBinContent(i+1,  list_nevt[i][1])
           hsig.SetBinContent(i+1,  list_nevt[i][2])
        #return hsig, hbkg

    def getZ(self, hsig, hbkg):
        nbins = hsig.GetNbinsX()
        Z = 0.0
        for i in range(nbins):
            Ns = hsig.GetBinContent(i+1)*self.totNsig
            Nb = hbkg.GetBinContent(i+1)*self.totNbkg
            if Nb == 0:
                continue
            Z += 2*((Ns+Nb)*math.log(1.+Ns/Nb) - Ns)
        Z = math.sqrt(Z)
        return Z
    def show_performance(self, realtest=0, systs = ['nominal']):
        self.read_alltreeq()
        if realtest:
            if not os.path.exists(self.treedir+'/bdt_sig_'+str(self.maxtree)+'.root'):
                hq_sig_train, hq_sig_test = self.get_q_hists_old(self.sigchain, 'sig')
            if not os.path.exists(self.treedir+'/bdt_bkg_'+str(self.maxtree)+'.root'):
                hq_bkg_train, hq_bkg_test = self.get_q_hists_old(self.bkgchain, 'bkg')
            for syst in systs:
                if syst != 'nominal' and not os.path.exists(self.treedir+'/bdt_bkg_'+syst+'_'+str(self.maxtree)+'.root'):
                    htmp_train, htmp_test = self.get_q_hists_old(self.bkgchainsysts[syst], 'bkg', syst)
        else:
            hq_sig_train, hq_sig_test = self.get_q_hists(self.sigchain, 'sig')
            hq_bkg_train, hq_bkg_test = self.get_q_hists(self.bkgchain, 'bkg')
        if realtest:
            histopt = '(20, -1, 1)'
            tmpsigchain=TChain('bdt')
            tmpsigchain.Add(self.treedir+'/bdt_sig_'+str(self.maxtree)+'.root')
            var = 'qbdt'
            var = 'tanh(0.25*'+var+')'
            tmpsigchain.Draw(var+'>>hq_sig_train'+histopt,'(trainflag==1)*(fweight)')
            tmpsigchain.Draw(var+'>>hq_sig_test'+histopt,'(trainflag!=1)*(fweight)')
            hq_sig_train = gROOT.FindObject('hq_sig_train')
            hq_sig_test = gROOT.FindObject('hq_sig_test')
            hq_sig_train = self.getoverflow(hq_sig_train)
            hq_sig_test = self.getoverflow(hq_sig_test)
            tmpbkgchain=TChain('bdt')
            tmpbkgchain.Add(self.treedir+'/bdt_bkg_'+str(self.maxtree)+'.root')
            tmpbkgchain.Draw(var+'>>hq_bkg_train'+histopt,'(trainflag==1)*(fweight)')
            tmpbkgchain.Draw(var+'>>hq_bkg_test'+histopt,'(trainflag!=1)*(fweight)')
            hq_bkg_train = gROOT.FindObject('hq_bkg_train')
            hq_bkg_test = gROOT.FindObject('hq_bkg_test')
            hq_bkg_train = self.getoverflow(hq_bkg_train)
            hq_bkg_test = self.getoverflow(hq_bkg_test)

            hq_sig_train.Scale(1.0 / hq_sig_train.Integral())
            hq_sig_test.Scale(1.0 / hq_sig_test.Integral())
            hq_bkg_train.Scale(1.0 / hq_bkg_train.Integral())
            hq_bkg_test.Scale(1.0 / hq_bkg_test.Integral())
        Z_train = self.getZ(hq_sig_train, hq_bkg_train)
        Z_test = self.getZ(hq_sig_test, hq_bkg_test)
        print 'Z_train =', Z_train
        print 'Z_test =', Z_test
        self.show_q_hist(0, hq_sig_train, hq_bkg_train, hq_sig_test, hq_bkg_test)
        self.show_q_hist(1, hq_sig_train, hq_bkg_train, hq_sig_test, hq_bkg_test)
        self.show_roc( hq_sig_train, hq_bkg_train, hq_sig_test, hq_bkg_test)
        if realtest:
            del tmpsigchain
            del tmpbkgchain
    def add_bdtweight2ntuple(self, filename, treename,  tag, iweight, lasttree, sum_Q, weight_total, alpha):
        #file0 = TFile(self.train'bdtweight_'+tag+'.root', 'update')
        file0 = TFile(self.treedir+'/'+filename, 'update')
        ntuple = file0.Get(treename)
        bdtweight = array('f',[0])
        qbdt = array('f',[0])
        Bbdtweight = ntuple.Branch('bdtweight'+str(iweight), bdtweight, 'bdtweight'+str(iweight)+'/F')
        Bqbdt = ntuple.Branch('qbdt'+str(iweight), qbdt, 'qbdt'+str(iweight)+'/F')
        for event in ntuple:
            cuttrain = 'event.'+self.trainflagvar+'<'+str(self.trainflagcut)
            exec('cutflag = %s' % cuttrain)
            exec('weight = event.%s' % self.weightvar)
            bdtweight[0] = 1.0
            qbdt[0] = 0.0
            if not cutflag:
                Bbdtweight.Fill()
                Bqbdt.Fill()
                continue
            if iweight>0:
                exec('oldbdtweight = event.bdtweight%s' % str(iweight-1))
                exec('oldqbdt = event.qbdt%s' % str(iweight-1))
                bdtweight[0] = oldbdtweight
                qbdt[0] = oldqbdt
            for node in lasttree:
                cutnode = node[3]
                cutnode = cutnode.replace(self.trainflagvar, 'event.'+self.trainflagvar)
                cutnode = cutnode.replace('&&', '&&event.')
                cutnode = cutnode.replace('&&', ' and ')
                exec('cutflag = %s' % cutnode)
                if not cutflag:
                    continue

                Q = node[4]
                purity = node[7]
                if Q < 0:
                    if tag != 'sig':
                        print 'Q =', Q, '<0 for ', tag, ' !!!!'
                scoresign = 1
                if purity < 0.5:
                    scoresign = -1
                qbdt[0] += (2.*purity-1.)*alpha
                if tag == 'sig':
                    bdtweight[0] *= weight_total*math.exp(-purity*alpha)
                else:
                    bdtweight[0] *= weight_total*math.exp(purity*alpha)
                break
            Bbdtweight.Fill()
            Bqbdt.Fill()
        ntuple.Write()
        file0.Close()
        del file0
        print 'success in adding bdtweight branch!'
        return


