import ROOT, itertools

from CMGTools.TTHAnalysis.postprocessing.framework.datamodel import Collection 
from CMGTools.TTHAnalysis.postprocessing.framework.eventloop import Module

from JHUGenMELA.MELA.mela import Mela, SimpleParticle_t, SimpleParticleCollection_t, TVar
from JHUGenMELA.MELA.pythonmelautils import MultiDimensionalCppArray, SelfDParameter, SelfDCoupling

ROOT.gSystem.Load("libJHUGenMELAMELA.so")

constants = dict(c_0minus = 0.297979440554, c_0hplus = 0.271880048944, c_0hplusza = 0.130395173298, c_0minusza = 0.104503154335)

class EventVars2Gam(Module):
    def __init__(self, doSystJEC=True):
        self.namebranches = ["costheta1", "costheta2", "Phi1", "costhetastar", "Phi", "HJJpz","M2G","costheta1d","costheta2d","Phid","costhetastard","Phi1d"]
        self.namebranches += ["pg1","pg4","pg2","pg1g2","pg1g4","pg2za","pg4za","pg1g2za","pg1g4za","D0minus","D0hplus","DCP","Dint","D0minus_za","D0hplus_za","Dint_za","DCP_za" ]
        self.namebranches += ["q2V1", "q2V2","Dphijj"]
        self.namebranches += ["ptH","pxH","pyH","pzH","EH","rapH","rapHJJ",
                              "pxj1", "pyj1", "pzj1", "Ej1",
                              "pxj2", "pyj2", "pzj2", "Ej2",
                              "ptpho1","pxpho1","pypho1","pzpho1","Epho1",
                              "ptpho2","pxpho2","pypho2","pzpho2","Epho2",
        ]

        self.systsJEC = {0:""} #, 1:"_jecUp", -1:"_jecDown"} if doSystJEC else {0:""}
        self.branches = []
        for var in self.systsJEC: self.branches.extend([br+self.systsJEC[var] for br in self.namebranches])
        self.m = Mela(13, 125)

    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        for B in self.namebranches:
            self.out.branch(B, "F")
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self,event):
        allret = {}

        # prepare output
        ret = dict([(name,-999.) for name in self.namebranches])

        if getattr(event,'dijet_Mjj') > 0: # this is the only way to suppress events with <2 VBF jets... 

            # flashgg tree is a scalar tree with very bad indexing ... 
            # take few vars by hand
            photons = [ROOT.TLorentzVector() for i in xrange(2)]
            jets = [ROOT.TLorentzVector() for i in xrange(2)]
            daughters = SimpleParticleCollection_t()
            associated = SimpleParticleCollection_t()
            for i,index in enumerate(['lead','sublead']):
                photons[i].SetPtEtaPhiM(getattr(event,'dipho_%sPt'%index),
                                        getattr(event,'dipho_%sEta'%index),
                                        getattr(event,'dipho_%sPhi'%index),
                                        0)
                daughters.push_back(SimpleParticle_t(22,photons[i]))
                jets[i].SetPtEtaPhiM(getattr(event,'dijet_%sPt'%index), 
                                     getattr(event,'dijet_%sEta'%index),
                                     getattr(event,'dijet_%sPhi'%index),
                                     0)
                associated.push_back(SimpleParticle_t(1,jets[i]))
     
            self.m.setInputEvent(daughters,associated)
     
            # calculate production probabilities
            self.m.setProcess(TVar.HSMHiggs, TVar.JHUGen, TVar.JJVBF)
            self.m.ghz2 = 1
            ret["pg1"] = self.m.computeProdP(False)
     
            # pure pseudoscalar
            self.m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.JJVBF)
            self.m.ghz1 = 0
            self.m.ghz4 = 1
            ret["pg4"] = self.m.computeProdP(False)
     
            # scalar - pseudoscalar with interference
            self.m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.JJVBF)
            self.m.ghz1 = 1
            self.m.ghz4 = 1
            ret["pg1g4"] = self.m.computeProdP(False) - ret["pg1"] - ret["pg4"]
     
            # pure loop induced production (gg, aa, za)
            self.m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.JJVBF)
            self.m.ghz1 = 0
            self.m.ghz2 = 1
            ret["pg2"] = self.m.computeProdP(False)
     
            # pure loop induced production (gg, aa, za) with interference
            self.m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.JJVBF)
            self.m.ghz1 = 1
            self.m.ghz2 = 1
            ret["pg1g2"] = self.m.computeProdP(False) - ret["pg1"] - ret["pg2"]
     
            # Zgamma production: a2 induced
            self.m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.JJVBF)
            self.m.ghz1 = 0
            self.m.ghzgs2 = 1
            ret["pg2za"] = self.m.computeProdP(False)
     
            # Zgamma production: a2 induced with interference
            self.m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.JJVBF)
            self.m.ghz1 = 1
            self.m.ghzgs2 = 1
            ret["pg1g2za"] = self.m.computeProdP(False) - ret["pg1"] - ret["pg2za"]
     
            # Zgamma production: a4 induced
            self.m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.JJVBF)
            self.m.ghz1 = 0
            self.m.ghzgs4 = 1
            ret["pg4za"] = self.m.computeProdP(False)
     
            # Zgamma production: a4 induced with interference
            self.m.setProcess(TVar.SelfDefine_spin0, TVar.JHUGen, TVar.JJVBF)
            self.m.ghz1 = 1
            self.m.ghzgs4 = 1
            ret["pg1g4za"] = self.m.computeProdP(False) - ret["pg1"] - ret["pg2za"]
            
            ### DISCRIMINANTS ###
            ret["D0minus"] = ret["pg1"] / (ret["pg1"] + constants["c_0minus"]*constants["c_0minus"]*ret["pg4"]) if ret["pg1"]*ret["pg4"]>0 else -900
            ret["D0hplus"] = ret["pg1"] / (ret["pg1"] + constants["c_0hplus"]*constants["c_0hplus"]*ret["pg2"]) if ret["pg1"]*ret["pg2"]>0 else -900
            ret["DCP"] = ret["pg1g4"] / (2 * (ret["pg1"] * ret["pg4"]) ** 0.5) if ret["pg1"]*ret["pg4"]>0 else -900
            ret["Dint"] = ret["pg1g2"] / (2 * (ret["pg1"] * ret["pg2"]) ** 0.5) if ret["pg1"]*ret["pg2"]>0 else -900
            
            ret["D0minus_za"] = ret["pg1"] / (ret["pg1"] + constants["c_0minusza"]*constants["c_0minusza"]*ret["pg4za"]) if ret["pg1"]*ret["pg4za"]>0 else -900
            ret["D0hplus_za"] = ret["pg1"] / (ret["pg1"] + constants["c_0hplusza"]*constants["c_0hplusza"]*ret["pg2za"]) if ret["pg1"]*ret["pg4"]>0 else -900
            ret["DCP_za"] = ret["pg1g4za"] / (2 * (ret["pg1"] * ret["pg4za"]) ** 0.5) if ret["pg1"]*ret["pg4za"]>0 else -900
            ret["Dint_za"] = ret["pg1g2za"] / (2 * (ret["pg1"] * ret["pg2za"]) ** 0.5) if ret["pg1"]*ret["pg2za"]>0 else -900
     
            ### kinematic variables ###
            ret["q2V1"], ret["q2V2"], ret["costheta1"], ret["costheta2"], ret["Phi"], ret["costhetastar"], ret["Phi1"]= self.m.computeVBFAngles()
            ret["M2G"], _, _, ret["costheta1d"],ret["costheta2d"], ret["Phid"], ret["costhetastard"], ret["Phi1d"]= self.m.computeDecayAngles()
            ret["HJJpz"] = sum((particle.second for particle in itertools.chain(daughters, associated)), ROOT.TLorentzVector()).Pz()
     
            pj1 = associated[0].second
            pj2 = associated[1].second
            if pj1.Pt() > pj2.Pt() :
                ret["Dphijj"] = pj1.DeltaPhi(pj2)
            else:
                ret["Dphijj"] = pj2.DeltaPhi(pj1)
     
            pH = sum((particle.second for particle in daughters), ROOT.TLorentzVector())
            ret["ptH"] = pH.Pt()
            ret["pxH"] = pH.Px()
            ret["pyH"] = pH.Py()
            ret["pzH"] = pH.Pz()
            ret["EH"] = pH.E()
            ret["rapH"] = pH.Rapidity()
         
            # copy the kinematics to a scalar tree with humanly written names...
            ppho1 = daughters[0].second
            ret["ptpho1"] = ppho1.Pt()
            ret["pxpho1"] = ppho1.Px()
            ret["pypho1"] = ppho1.Py()
            ret["pzpho1"] = ppho1.Pz()
            ret["Epho1"] = ppho1.E()
            
            ppho2 = daughters[1].second
            ret["ptpho2"] = ppho2.Pt()
            ret["pxpho2"] = ppho2.Px()
            ret["pypho2"] = ppho2.Py()
            ret["pzpho2"] = ppho2.Pz()
            ret["Epho2"] = ppho2.E()
         
            pj1 = associated[0].second
            ret["pxj1"] = pj1.Px()
            ret["pyj1"] = pj1.Py()
            ret["pzj1"] = pj1.Pz()
            ret["Ej1"] = pj1.E()          
            pj2 = associated[1].second
            ret["pxj2"] = pj2.Px()
            ret["pyj2"] = pj2.Py()
            ret["pzj2"] = pj2.Pz()
            ret["Ej2"] = pj2.E()
            phjj = pH + pj1 + pj2
            ret["rapHJJ"] = phjj.Rapidity()

        #for br in self.namebranches:
        #     allret[br+self.label+self.systsJEC[var]] = ret[br]
	# return allret

        for br in self.namebranches:
            self.out.fillBranch(br,ret[br])
        return True

if __name__ == '__main__':
    from sys import argv
    file = ROOT.TFile(argv[1])
    tree = file.Get("/vbfTagDumper/trees/vbf_125_13TeV_GeneralDipho")
    tree.vectorTree = True
    #tree.AddFriend("sf/t",argv[2])
    class Tester(Module):
        def __init__(self, name):
            Module.__init__(self,name,None)
            self.sf = EventVars2Gam()
        def analyze(self,ev):
            print "\nrun %6d lumi %4d event %d: mass %f" % (ev.run, ev.lumi, ev.event, ev.dipho_mass)
            print self.sf(ev)
    el = eventLoop([ Tester("tester") ])
    el.loop([tree], maxEvents = 50)

        
recoMelaCentral = lambda : EventVars2Gam(doSystJEC=False)
