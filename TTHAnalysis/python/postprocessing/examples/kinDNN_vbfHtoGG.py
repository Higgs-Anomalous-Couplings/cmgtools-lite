import ROOT, itertools

from CMGTools.TTHAnalysis.postprocessing.framework.datamodel import Collection 
from CMGTools.TTHAnalysis.postprocessing.framework.eventloop import Module
from CMGTools.TTHAnalysis.tools.mvaTool import *

P="%s/src/CMGTools/TTHAnalysis/data/kinMVA/vbfh2gg/" % os.environ['CMSSW_BASE'];
        
class KinDNN:
    def __init__(self, algo, varset, path, rarity=False):
        self.catnames = ['BKG','SM','BSM']
        self.MVAs = {}
        if varset == "minimal":
            self._vars = [  
                MVAVar("dijet_dipho_dphi := dijet_dipho_dphi", func = lambda ev : ev.dijet_dipho_dphi if ev.dijet_Mjj > 0 else 0),
                MVAVar("dipho_lead_ptoM := dipho_lead_ptoM", func = lambda ev : ev.dipho_lead_ptoM if ev.dijet_Mjj > 0 else 0),
                MVAVar("dijet_dipho_pt := dijet_dipho_pt", func = lambda ev : ev.dijet_dipho_pt if ev.dijet_Mjj > 0 else 0),
                MVAVar("dijet_abs_dEta := dijet_abs_dEta", func = lambda ev : ev.dijet_abs_dEta if ev.dijet_Mjj > 0 else 0),
                MVAVar("dijet_Zep := dijet_Zep", func = lambda ev : ev.dijet_Zep if ev.dijet_Mjj > 0 else 0),
                MVAVar("dijet_minDRJetPho := dijet_minDRJetPho", func = lambda ev : ev.dijet_minDRJetPho if ev.dijet_Mjj > 0 else 0),
                MVAVar("dipho_sublead_ptoM := dipho_sublead_ptoM", func = lambda ev : ev.dipho_sublead_ptoM if ev.dijet_Mjj > 0 else 0),
                MVAVar("jet2_pt := jet2_pt", func = lambda ev : ev.jet2_pt if ev.dijet_Mjj > 0 else 0),
                MVAVar("CosPhi := CosPhi", func = lambda ev : ev.CosPhi if ev.dijet_Mjj > 0 else 0),
                MVAVar("dijet_dphi := dijet_dphi", func = lambda ev : ev.dijet_dphi if ev.dijet_Mjj > 0 else 0),
            ]
            training = path+("multiclass_bkgAndGGH_vbfsm_vbfbsm_%s.weights.xml" % (algo))
            self.MVAs[algo] = MVATool(algo, training, self._vars, rarity=rarity, nClasses=3) 
    def listBranches(self):
        return self.MVAs.keys()
    def __call__(self,event):
        out = {}
        for name, mva in self.MVAs.iteritems():
            x = mva(event)
            for i,cat in enumerate(self.catnames):
                out['%s_%s'%(name,cat)] = x[i]
        return out

class KinDNN_VbfH2GG(Module):
    def __init__(self, path):
        self.nodes = {'BKG' : 'pbkg', 'SM' : 'psm', 'BSM' : 'pbsm'}
        self.algo = "DNN"
        varset = "minimal"
        self.kinDNN = KinDNN(self.algo, varset, P)
        self.branches = ["vbf{algo}_{node}".format(algo=self.algo,node=n) for k,n in self.nodes.iteritems()]

    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        for B in self.branches:
            self.out.branch(B, "F")
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self,event):
        # prepare output
        ret = dict([(name,-999.) for name in self.branches])
        if getattr(event,'dijet_Mjj') > 0: # this is the only way to suppress events with <2 VBF jets... 
            out = self.kinDNN(event)
            for cat,var in self.nodes.iteritems():
                ret["vbf{algo}_{node}".format(algo=self.algo,node=self.nodes[cat])] = out["%s_%s"%(self.algo,cat)]
        else:
            for cat,var in self.nodes.iteritems():
                ret["vbf{algo}_{node}".format(algo=self.algo,node=self.nodes[cat])] = -999.

        for br in self.branches:
            self.out.fillBranch(br,ret[br])
        return True


if __name__ == '__main__':
    from sys import argv
    file = ROOT.TFile(argv[1])
    tree = file.Get("/tagsDumper/trees/vbfh_13TeV_VBFTag")
    tree.vectorTree = True
    #tree.AddFriend("sf/t",argv[2])
    class Tester(Module):
        def __init__(self, name):
            Module.__init__(self,name,None)
            self.sf = KinDNN_VbfH2GG('kinMVA','minimal','DNN',P)
        def analyze(self,ev):
            print "\nrun %6d lumi %4d event %d: mass %f" % (ev.run, ev.lumi, ev.event, ev.dipho_mass)
            print self.sf(ev)
    el = eventLoop([ Tester("tester") ])
    el.loop([tree], maxEvents = 50)

kinDNN = lambda : KinDNN_VbfH2GG(P)

