import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import sys
import copy

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import * #deltaR, matching etc..

from PhysicsTools.NanoAODTools.postprocessing.analysis.higgs.vhbb.applysmearing import SmearApplicator

class VHbbProducer(Module):
    def __init__(self, isMC, era, useCMVA=False,isVjets=False):
        self.era = era
        self.isMC = isMC
        self.useCMVA = useCMVA
        self.isVjets = isVjets
        self.smearapplicator = SmearApplicator(year=era, isdata=(not isMC), useV13=3, doscaling=False)
        pass
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Vtype",   "I");
        self.out.branch("V_pt",    "F");
        self.out.branch("V_eta",   "F");
        self.out.branch("V_phi",   "F");
        self.out.branch("V_mass",  "F");
        self.out.branch("V_mt",    "F");
        self.out.branch("TTtype",  "I");
        self.out.branch("TTmll",   "F");
        self.out.branch("TTdrll",  "F");
        self.out.branch("TT_lep1_pt",   "F");
        self.out.branch("TT_lep1_phi",  "F");
        self.out.branch("TT_lep1_eta",  "F");
        self.out.branch("TT_lep1_mass", "F");
        self.out.branch("TT_lep2_pt",   "F");
        self.out.branch("TT_lep2_phi",  "F");
        self.out.branch("TT_lep2_eta",  "F");
        self.out.branch("TT_lep2_mass", "F");
        self.out.branch("Jet_lepFilter",  "O", 1, "nJet");
        self.out.branch("vLidx",  "I", 2);
        self.out.branch("hJidx",  "I", 2);
        self.out.branch("hJidxCMVA",  "I", 2);
        self.out.branch("H_pt",  "F");
        self.out.branch("H_eta",  "F");
        self.out.branch("H_phi",  "F");
        self.out.branch("H_mass",  "F");
        self.out.branch("H_pt_JERUp",  "F");
        self.out.branch("H_eta_JERUp",  "F");
        self.out.branch("H_phi_JERUp",  "F");
        self.out.branch("H_mass_JERUp",  "F");
        self.out.branch("H_pt_JERDown",  "F");
        self.out.branch("H_eta_JERDown",  "F");
        self.out.branch("H_phi_JERDown",  "F");
        self.out.branch("H_mass_JERDown",  "F");
        self.out.branch("HCMVA_pt",  "F");
        self.out.branch("HCMVA_eta",  "F");
        self.out.branch("HCMVA_phi",  "F");
        self.out.branch("HCMVA_mass",  "F");
        self.out.branch("HFSR_pt",  "F");
        self.out.branch("HFSR_eta",  "F");
        self.out.branch("HFSR_phi",  "F");
        self.out.branch("HFSR_mass",  "F");
        self.out.branch("HFSR_pt_JERUp",  "F");
        self.out.branch("HFSR_eta_JERUp",  "F");
        self.out.branch("HFSR_phi_JERUp",  "F");
        self.out.branch("HFSR_mass_JERUp",  "F");
        self.out.branch("HFSR_pt_JERDown",  "F");
        self.out.branch("HFSR_eta_JERDown",  "F");
        self.out.branch("HFSR_phi_JERDown",  "F");
        self.out.branch("HFSR_mass_JERDown",  "F");
        self.out.branch("SA_Ht",  "F");
        self.out.branch("SA5",  "F");
        self.out.branch("Jet_Pt", "F", 1, "nJet");
        self.out.branch("Jet_PtReg", "F", 1, "nJet");
        self.out.branch("Jet_PtRegUp", "F", 1, "nJet");
        self.out.branch("Jet_PtRegDown", "F", 1, "nJet");
        self.out.branch("Jet_CvsL", "F", 1, "nJet")
        self.out.branch("Jet_CvsB", "F", 1, "nJet")
        self.out.branch("Jet_DeepFlavCvsL", "F", 1, "nJet")
        self.out.branch("Jet_DeepFlavCvsB", "F", 1, "nJet")
        self.out.branch("MET_Pt","F");
        self.out.branch("MET_Phi","F");
        self.out.branch("GenJetWithNeutrinos_pt","F",1, "nGenJet")
        self.out.branch("GenJetWithNeutrinos_eta","F",1, "nGenJet")
        self.out.branch("GenJetWithNeutrinos_phi","F",1, "nGenJet")
        self.out.branch("GenJetWithNeutrinos_mass","F",1, "nGenJet")

        ## for the boosted analysis
        self.out.branch("Pt_fjidx",  "I");        
        self.out.branch("Msd_fjidx",  "I");
        self.out.branch("Hbb_fjidx",  "I");
        
        self.out.branch("SAptfj_HT",  "F");
        self.out.branch("SAptfj5",  "F");
        self.out.branch("SAmfj_HT",  "F");
        self.out.branch("SAmfj5",  "F");
        self.out.branch("SAhbbfj_HT",  "F");
        self.out.branch("SAhbbfj5",  "F");
        
        self.out.branch("FatJet_lepFilter",  "O", 1, "nFatJet");
        self.out.branch("FatJet_Pt", "F", 1, "nFatJet");
        self.out.branch("FatJet_Msoftdrop", "F", 1, "nFatJet");
        
        self.out.branch("FatJet_FlavourComposition", "I", 1, "nFatJet"); #55 bb, #54 bc, #5 b, #4 c, #44 cc, #1 other
        
        self.out.branch("FatJet_HiggsProducts", "O", 1, "nFatJet");
        self.out.branch("FatJet_WProducts", "O", 1, "nFatJet");
        self.out.branch("FatJet_ZProducts", "O", 1, "nFatJet");
        
        ## Gen information
        self.out.branch("nGenStatus2bHad", "I");
        self.out.branch("GenBJ1_pt", "F");
        self.out.branch("GenBJ1_genjet_pt", "F");
        self.out.branch("GenBJ1_eta", "F");
        self.out.branch("GenBJ1_phi", "F");
        self.out.branch("GenBJ1_mass", "F");
        self.out.branch("GenBJ1_genjet_mass", "F");
        self.out.branch("GenBJ1_index", "I");
        self.out.branch("GenBJ1_genjet_index", "I");
        self.out.branch("GenBJ2_pt", "F");
        self.out.branch("GenBJ2_genjet_pt", "F");
        self.out.branch("GenBJ2_eta", "F");
        self.out.branch("GenBJ2_phi", "F");
        self.out.branch("GenBJ2_mass", "F");
        self.out.branch("GenBJ2_genjet_mass", "F");
        self.out.branch("GenBJ2_index", "I");
        self.out.branch("GenBJ2_genjet_index", "I");
        self.out.branch("GenBJJ_pt", "F");
        self.out.branch("GenBJJ_genjet_pt", "F");
        self.out.branch("GenBJJ_eta", "F");
        self.out.branch("GenBJJ_phi", "F");
        self.out.branch("GenBJJ_mass", "F");
        self.out.branch("GenBJJ_genjet_mass", "F");
        self.out.branch("GenBJJ_dPhi", "F");
        self.out.branch("GenBJJ_dR", "F");
        self.out.branch("GenBJJ_dEta", "F");
        self.out.branch("GenLepIndex1","I");
        self.out.branch("GenLepIndex2","I");
        self.out.branch("GenLep_GenBJ1_dR","F");
        self.out.branch("GenLep_GenBJ1_dEta","F");
        self.out.branch("GenLep_GenBJ1_dPhi","F");
        self.out.branch("GenTop1_mass","F");
        self.out.branch("GenLep_GenBJ2_dR","F");
        self.out.branch("GenLep_GenBJ2_dEta","F");
        self.out.branch("GenLep_GenBJ2_dPhi","F");
        self.out.branch("GenTop2_mass","F");
        self.out.branch("nGenTop","I")
        self.out.branch("nW","I")
        self.out.branch("nWlep","I")
        self.out.branch("nGenVBosons","I")
        self.out.branch("LeadGenVBoson_pt","F")
        self.out.branch("LeadGenVBoson_eta","F")
        self.out.branch("LeadGenVBoson_phi","F")
        self.out.branch("LeadGenVBoson_pdgId","F")
        self.out.branch("LeadNeutrinoFromGenVBoson_pt", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def matchSoftActivity(self,jets,saJets,dR=0.4) :
        matched=set()
        for saj in saJets:
            for j in jets :
                if deltaR(saj,j) < dR :
                    matched.add(saj)            
        return matched
    
    def matchSoftActivityFSR(self,jet1,jet2,saJets,dR=0.4) :
        matched=set()
        drjj = deltaR(jet1,jet2)
        sumDeltaRMin = drjj + 2*dR
        for saj in saJets:
            dr1 = deltaR(saj,jet1)
            dr2 = deltaR(saj,jet2)
            if ((dr1+dr2) < sumDeltaRMin):
                matched.add(saj)
        return matched
			
    def pt(self, jet, isMC, noReg=False, sysVar=0):
        ## the MC has JER smearing applied which has output branch Jet_pt_nom which should be compared 
        ## with data branch Jet_pt. This essentially aliases the two branches to one common jet pt variable.
        if noReg:
            if isMC:
                return jet.pt_nom
            else: #jet_pt_nom contains re-corrected jet pT if jecRecalibrator has been run
                if hasattr(jet, 'pt_nom'):
                    return jet.pt_nom
                else:
                    return jet.pt		 
        else:
            genJet=None
            if isMC and jet.genJetIdx >=0 and  jet.genJetIdx < len(self.genJetsWithNeutrinos) :
                genJet=self.genJetsWithNeutrinos[jet.genJetIdx]

            jet_pt = jet.pt_nom if hasattr(jet, 'pt_nom') else jet.pt

            # Not using rho at the moment (luckily?)
            # If initialized for data, does not do any smearing
            pt_variations = self.smearapplicator.get_smear(jet_pt, jet.bRegCorr, genJet.Pt() if genJet else 0)

            if sysVar==0: # nominal
                return pt_variations.nominal

            elif sysVar==1: # up
                return pt_variations.up

            elif sysVar==-1: # down
                return pt_variations.down
    
    def met(self, met, isMC):
        ## the MC has JER smearing applied which has output branch met_[pt/phi]_nom which should be compared 
        ## with data branch MET_[pt/phi]. This essentially aliases the two branches to one common variable.
        if isMC:
            return (met.pt_nom,met.phi_nom)
        else:
            if hasattr(met, 'pt_nom'):#pt_nom contains re-corrected pT if jecRecalibrator has been run
                return (met.pt_nom,met.phi_nom)
            else:
                return (met.pt,met.phi)
	 
    def msoftdrop(self, jet, isMC):
        ## the MC has JER smearing applied which has output branch Jet_pt_smeared which should be compared 
        ## with data branch Jet_pt. This essentially aliases the two branches to one common jet pt variable.
        if isMC:
            return jet.msoftdrop_nom
        else:
            return jet.msoftdrop
 
    def btag(self, jet):
        if (self.useCMVA):
            return jet.btagCMVA
        else:
            return jet.btagDeepB

    def cvsltag(self, jet):
        btagDeepL = 1.-(jet.btagDeepC+jet.btagDeepB)
        if jet.btagDeepB >= 0. and jet.btagDeepB < 1. and jet.btagDeepC >= 0. and btagDeepL >= 0.:
            return jet.btagDeepC/(1.-jet.btagDeepB)
        else:
            return -1

    def cvsbtag(self, jet):
        btagDeepL = 1.-(jet.btagDeepC+jet.btagDeepB)
        if jet.btagDeepB > 0. and jet.btagDeepC > 0. and btagDeepL >= 0.:
            return jet.btagDeepC/(jet.btagDeepC+jet.btagDeepB)
        else:
            return -1

    def deepflavcvsltag(self, jet):
        if not hasattr(jet, 'btagDeepFlavC'):           # only Nano V5 onwards
            return -99. 
        btagDeepFlavL = 1.-(jet.btagDeepFlavC+jet.btagDeepFlavB)
        if jet.btagDeepFlavB >= 0. and jet.btagDeepFlavB < 1. and jet.btagDeepFlavC >= 0. and btagDeepFlavL >= 0.:
            return jet.btagDeepFlavC/(1.-jet.btagDeepFlavB)
        else:
            return -1

    def deepflavcvsbtag(self, jet):
        if not hasattr(jet, 'btagDeepFlavC'):           # only Nano V5 onwards
            return -99.
        btagDeepFlavL = 1.-(jet.btagDeepFlavC+jet.btagDeepFlavB)
        if jet.btagDeepFlavB > 0. and jet.btagDeepFlavC > 0. and btagDeepFlavL >= 0.:
            return jet.btagDeepFlavC/(jet.btagDeepFlavC+jet.btagDeepFlavB)
        else:
            return -1

    def elid(self, el, wp):
        if (wp == "80"):
            return el.mvaFall17V2Iso_WP80
        elif (wp == "90"):
            return el.mvaFall17V2Iso_WP90

    def statusFlags_dict(self, bit):
        Dict={0 : "isPrompt", 1 : "isDecayedLeptonHadron", 2 : "isTauDecayProduct", 3 : "isPromptTauDecayProduct", 4 : "isDirectTauDecayProduct", 5 : "isDirectPromptTauDecayProduct", 6 : "isDirectHadronDecayProduct", 7 : "isHardProcess", 8 : "fromHardProcess", 9 : "isHardProcessTauDecayProduct", 10 : "isDirectHardProcessTauDecayProduct", 11 : "fromHardProcessBeforeFSR", 12 : "isFirstCopy", 13 : "isLastCopy", 14 : "isLastCopyBeforeFSR" }
        return Dict[bit] 

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        electrons = list(Collection(event, "Electron"))
        muons = list(Collection(event, "Muon"))
        jets = list(Collection(event, "Jet"))
        met = Object(event, "MET")
        sa = Collection(event, "SoftActivityJet")
        fatjets = list(Collection(event, "FatJet"))
        subjets = Collection(event, "SubJet")
        if self.isMC:
            genParticles = Collection(event, "GenPart")
            genJets = list(Collection(event, "GenJet"))
            self.genJetsWithNeutrinos = []
            genWNPts = [-99]*len(genJets)
            genWNEtas = [-99]*len(genJets)
            genWNPhis = [-99]*len(genJets)
            genWNMasses = [-99]*len(genJets)
            for iJet in xrange(len(genJets)):
                genJet = ROOT.TLorentzVector()
                genJet.SetPtEtaPhiM(genJets[iJet].pt,genJets[iJet].eta,genJets[iJet].phi,genJets[iJet].mass)
                self.genJetsWithNeutrinos.append(genJet)
                for genPart in genParticles:
                    if (genPart.status==1 and (abs(genPart.pdgId)==12 or abs(genPart.pdgId)==14 or abs(genPart.pdgId)==16)):
                        neutrino = ROOT.TLorentzVector() 
                        neutrino.SetPtEtaPhiM(genPart.pt,genPart.eta,genPart.phi,genPart.mass)
                        if (neutrino.DeltaR(self.genJetsWithNeutrinos[iJet])>0.4): continue 
                        # add back neutrino to genJet
                        self.genJetsWithNeutrinos[iJet] = self.genJetsWithNeutrinos[iJet] + neutrino
                genWNPts[iJet] = self.genJetsWithNeutrinos[iJet].Pt()
                genWNEtas[iJet] = self.genJetsWithNeutrinos[iJet].Eta()
                genWNPhis[iJet] = self.genJetsWithNeutrinos[iJet].Phi()
                genWNMasses[iJet] = self.genJetsWithNeutrinos[iJet].M()

            self.out.fillBranch("GenJetWithNeutrinos_pt",genWNPts)
            self.out.fillBranch("GenJetWithNeutrinos_eta",genWNEtas)
            self.out.fillBranch("GenJetWithNeutrinos_phi",genWNPhis)
            self.out.fillBranch("GenJetWithNeutrinos_mass",genWNMasses)

        metPt,metPhi = self.met(met,self.isMC)
        self.out.fillBranch("MET_Pt",metPt)
        self.out.fillBranch("MET_Phi",metPhi) 
      
        Vtype = -1

        wElectrons = [x for x in electrons if self.elid(x,"80") and x.pt > 25 and x.pfRelIso03_all < 0.12]      
        wMuons = [x for x in muons if x.pt > 25 and x.tightId >= 1 and x.pfRelIso04_all < 0.15 and abs(x.dxy) < 0.05 and abs(x.dz) < 0.2]
        zElectrons = [x for x in electrons if x.pt > 20 and self.elid(x,"90") and x.pfRelIso03_all < 0.15]
        zMuons = [x for x in muons if x.pt > 20 and x.pfRelIso04_all < 0.25 and abs(x.dxy) < 0.05 and abs(x.dz) < 0.2] # muons already preselected with looseId requirement

        zMuons.sort(key=lambda x:x.pt,reverse=True)
        zElectrons.sort(key=lambda x:x.pt,reverse=True)

        vLeptons = [] # decay products of V
        vLidx = [-1,-1] # indices in lepton collection of selected leptons
        if len(zMuons) >= 2:
            if zMuons[0].pt > 20:
                for i in xrange(1,len(zMuons)):
                    if zMuons[0].charge * zMuons[i].charge < 0:
                        Vtype = 0
                        vLeptons = [zMuons[0],zMuons[i]]
                        vLidx[0] = muons.index(zMuons[0])
                        vLidx[1] = muons.index(zMuons[1])
                        break
        elif len(zElectrons) >= 2:
            if zElectrons[0].pt > 20:
                for i in xrange(1,len(zElectrons)):
                    if zElectrons[0].charge * zElectrons[i].charge < 0:
                        Vtype = 1
                        vLeptons = [zElectrons[0],zElectrons[i]]
                        vLidx[0] = electrons.index(zElectrons[0])
                        vLidx[1] = electrons.index(zElectrons[1])
                        break
        elif len(wElectrons) + len(wMuons) == 1:
            if len(wMuons) == 1:
                Vtype = 2
                vLeptons = [wMuons[0]]
                vLidx[0] = muons.index(wMuons[0])
            if len(wElectrons) == 1:
                Vtype=3
                vLeptons = [wElectrons[0]]
                vLidx[0] = electrons.index(wElectrons[0])
        elif len(zElectrons) + len(zMuons) > 0:
            Vtype = 5
        else:
            Vtype = 4
            if metPt < 150:
                Vtype = -1
        self.out.fillBranch("Vtype",Vtype)
        self.out.fillBranch("vLidx",vLidx)

        ## add branches for some basic V kinematics
        V = ROOT.TLorentzVector()
        for vLepton in vLeptons:
            vLepton_4vec = ROOT.TLorentzVector()
            vLepton_4vec.SetPtEtaPhiM(vLepton.pt,vLepton.eta,vLepton.phi,vLepton.mass)
            V = V + vLepton_4vec
        if Vtype >=2 and Vtype<=4:
            met_4vec = ROOT.TLorentzVector()
            met_4vec.SetPtEtaPhiM(metPt,0.,metPhi,0.) # only use met vector to derive transverse quantities
            V = V + met_4vec
        self.out.fillBranch("V_pt",V.Pt())
        self.out.fillBranch("V_eta",V.Eta())
        self.out.fillBranch("V_phi",V.Phi())
        self.out.fillBranch("V_mass",V.M())
        self.out.fillBranch("V_mt",V.Mt())
        

        TTtype = -1
        TTmll = 0
        TTdrll = 99
        
        tightElectrons = [x for x in electrons if self.elid(x,"80") and x.pt > 25 and abs(x.eta)<2.4]      
        tightMuons =     [x for x in muons if x.pt > 25 and x.tightId >= 1 and x.pfRelIso04_all < 0.15]
        looseElectrons = [x for x in electrons if self.elid(x,"80") and x.pt > 15 and abs(x.eta)<2.4]      
        looseMuons =     [x for x in muons if x.pt > 15 and x.tightId >= 1 and x.pfRelIso04_all < 0.15]

        nLooseLep=len(looseElectrons)+len(looseMuons)
        nTightLep=len(tightElectrons)+len(tightMuons)

        ## TTtype lepton selection
        if nLooseLep==2 and nTightLep>0:
            mass_e=0.511e-6
            mass_m=105.658e-6
            lep1=ROOT.TLorentzVector()
            lep2=ROOT.TLorentzVector()
            if len(looseElectrons)==2:
                lep1.SetPtEtaPhiM(looseElectrons[0].pt,looseElectrons[0].eta,looseElectrons[0].phi,mass_e)
                lep2.SetPtEtaPhiM(looseElectrons[1].pt,looseElectrons[1].eta,looseElectrons[1].phi,mass_e)
                TTtype = 1
                self.out.fillBranch("TT_lep1_pt",   looseElectrons[0].pt)
                self.out.fillBranch("TT_lep1_phi",  looseElectrons[0].phi)
                self.out.fillBranch("TT_lep1_eta",  looseElectrons[0].eta)
                self.out.fillBranch("TT_lep1_mass", mass_e)
                self.out.fillBranch("TT_lep2_pt",   looseElectrons[1].pt)
                self.out.fillBranch("TT_lep2_phi",  looseElectrons[1].phi)
                self.out.fillBranch("TT_lep2_eta",  looseElectrons[1].eta)
                self.out.fillBranch("TT_lep2_mass", mass_e)
            elif len(looseMuons)==2:
                lep1.SetPtEtaPhiM(looseMuons[0].pt,looseMuons[0].eta,looseMuons[0].phi,mass_e)
                lep2.SetPtEtaPhiM(looseMuons[1].pt,looseMuons[1].eta,looseMuons[1].phi,mass_e)
                TTtype = 3
                self.out.fillBranch("TT_lep1_pt",   looseMuons[0].pt)
                self.out.fillBranch("TT_lep1_phi",  looseMuons[0].phi)
                self.out.fillBranch("TT_lep1_eta",  looseMuons[0].eta)
                self.out.fillBranch("TT_lep1_mass", mass_m)
                self.out.fillBranch("TT_lep2_pt",   looseMuons[1].pt)
                self.out.fillBranch("TT_lep2_phi",  looseMuons[1].phi)
                self.out.fillBranch("TT_lep2_eta",  looseMuons[1].eta)
                self.out.fillBranch("TT_lep2_mass", mass_m)
            else:
                lep1.SetPtEtaPhiM(looseMuons[0].pt,looseMuons[0].eta,looseMuons[0].phi,mass_e)
                lep2.SetPtEtaPhiM(looseElectrons[0].pt,looseElectrons[0].eta,looseElectrons[0].phi,mass_e)
                TTtype = 2
                self.out.fillBranch("TT_lep1_pt",   looseMuons[0].pt)
                self.out.fillBranch("TT_lep1_phi",  looseMuons[0].phi)
                self.out.fillBranch("TT_lep1_eta",  looseMuons[0].eta)
                self.out.fillBranch("TT_lep1_mass", mass_m)
                self.out.fillBranch("TT_lep2_pt",   looseElectrons[0].pt)
                self.out.fillBranch("TT_lep2_phi",  looseElectrons[0].phi)
                self.out.fillBranch("TT_lep2_eta",  looseElectrons[0].eta)
                self.out.fillBranch("TT_lep2_mass", mass_e)

            llp4=lep1+lep2
            TTmll=llp4.M()
            TTdrll=lep1.DeltaR(lep2)
        else:
            self.out.fillBranch("TT_lep1_pt",   -99);
            self.out.fillBranch("TT_lep1_phi",  -99);
            self.out.fillBranch("TT_lep1_eta",  -99);
            self.out.fillBranch("TT_lep1_mass", -99);
            self.out.fillBranch("TT_lep2_pt",   -99);
            self.out.fillBranch("TT_lep2_phi",  -99);
            self.out.fillBranch("TT_lep2_eta",  -99);
            self.out.fillBranch("TT_lep2_mass", -99);


        self.out.fillBranch("TTtype",TTtype)
        self.out.fillBranch("TTmll",TTmll)
        self.out.fillBranch("TTdrll",TTdrll)
        
        ## filter jets that overlap with any of the selected leptons
        allLeptons = zElectrons[:]
        allLeptons.extend(zMuons)
        allLeptons.extend(wElectrons)
        allLeptons.extend(wMuons)
        jetFilterFlags = [True]*len(jets)
        fatjetFilterFlags = [True]*len(fatjets)
        #for jet in jets:
        #    jet.jetFilter = True
        for lepton in allLeptons:
            jetInd = lepton.jetIdx
            if jetInd >= 0:
                jetFilterFlags[jetInd] = False
                #jets[jetInd].jetFilter = False
        self.out.fillBranch("Jet_lepFilter",jetFilterFlags)
 
        for fatjet in fatjets:
            fatjet.jetFilter = True
            for lepton in allLeptons:
               if deltaR(fatjet,lepton) < 0.8:
                  fatjetFilterFlags[fatjets.index(fatjet)] = False

        self.out.fillBranch("FatJet_lepFilter",fatjetFilterFlags)

        ## alias JER-smeared MC jet pT and data jet pT to the same
        ## branch name
        jetPts = [-99.]*len(jets)
        jetPtRegs = [-99.]*len(jets)
        jetPtRegsUp = [-99.]*len(jets)
        jetPtRegsDown = [-99.]*len(jets)
        jetCvsL = [-99.]*len(jets)
        jetCvsB = [-99.]*len(jets)
        jetDeepFlavCvsL = [-99.]*len(jets)
        jetDeepFlavCvsB = [-99.]*len(jets)
        for i in xrange(len(jets)):
            jetPts[i] = self.pt(jets[i],self.isMC, True)
            jetPtRegs[i] = self.pt(jets[i],self.isMC)
            jetPtRegsUp[i] = self.pt(jets[i],self.isMC,sysVar=1)
            jetPtRegsDown[i] = self.pt(jets[i],self.isMC,sysVar=-1)
            jetCvsL[i]       = self.cvsltag(jets[i])
            jetCvsB[i]       = self.cvsbtag(jets[i])
            jetDeepFlavCvsL[i]       = self.deepflavcvsltag(jets[i])
            jetDeepFlavCvsB[i]       = self.deepflavcvsbtag(jets[i])

        self.out.fillBranch("Jet_Pt",jetPts)
        self.out.fillBranch("Jet_PtReg",jetPtRegs)
        self.out.fillBranch("Jet_PtRegUp",jetPtRegsUp)
        self.out.fillBranch("Jet_PtRegDown",jetPtRegsDown)
        self.out.fillBranch("Jet_CvsL", jetCvsL)
        self.out.fillBranch("Jet_CvsB", jetCvsB)
        self.out.fillBranch("Jet_DeepFlavCvsL", jetDeepFlavCvsL)
        self.out.fillBranch("Jet_DeepFlavCvsB", jetDeepFlavCvsB)

        fatjetPts = [-99.]*len(fatjets)
        for i in xrange(len(fatjets)):
            fatjetPts[i] = self.pt(fatjets[i],self.isMC,True)

        fatjetMSD = [-99.]*len(fatjets)
        for i in xrange(len(fatjets)):
            fatjetMSD[i] = self.msoftdrop(fatjets[i],self.isMC)

        self.out.fillBranch("FatJet_Pt",fatjetPts)
        self.out.fillBranch("FatJet_Msoftdrop",fatjetMSD)

        ## Add explicit indices for selected H(bb) candidate jets
        jetsForHiggs = [x for x in jets if x.lepFilter and x.puId>0 and x.jetId>0 and self.pt(x,self.isMC)>20 and abs(x.eta)<2.5]
        if (len(jetsForHiggs) >= 2): 
            hJets = sorted(jetsForHiggs, key = lambda jet : self.btag(jet), reverse=True)[0:2]
            hJidx = [jets.index(x) for x in hJets]
            self.out.fillBranch("hJidx",hJidx)

            ## Save a few basic reco. H kinematics
            hj1 = ROOT.TLorentzVector()
            hj2 = ROOT.TLorentzVector()
            hj1.SetPtEtaPhiM(self.pt(jets[hJidx[0]],self.isMC),jets[hJidx[0]].eta,jets[hJidx[0]].phi,jets[hJidx[0]].mass)
            hj2.SetPtEtaPhiM(self.pt(jets[hJidx[1]],self.isMC),jets[hJidx[1]].eta,jets[hJidx[1]].phi,jets[hJidx[1]].mass)
            hbb = hj1 + hj2
            self.out.fillBranch("H_pt",hbb.Pt())
            self.out.fillBranch("H_phi",hbb.Phi())
            self.out.fillBranch("H_eta",hbb.Eta())
            self.out.fillBranch("H_mass",hbb.M())
            hj1up = ROOT.TLorentzVector()
            hj2up = ROOT.TLorentzVector()
            hj1up.SetPtEtaPhiM(self.pt(jets[hJidx[0]],self.isMC,sysVar=1),jets[hJidx[0]].eta,jets[hJidx[0]].phi,jets[hJidx[0]].mass)
            hj2up.SetPtEtaPhiM(self.pt(jets[hJidx[1]],self.isMC,sysVar=1),jets[hJidx[1]].eta,jets[hJidx[1]].phi,jets[hJidx[1]].mass)
            hbbup = hj1up + hj2up
            self.out.fillBranch("H_pt_JERUp",hbbup.Pt())
            self.out.fillBranch("H_phi_JERUp",hbbup.Phi())
            self.out.fillBranch("H_eta_JERUp",hbbup.Eta())
            self.out.fillBranch("H_mass_JERUp",hbbup.M())
            hj1down = ROOT.TLorentzVector()
            hj2down = ROOT.TLorentzVector()
            hj1down.SetPtEtaPhiM(self.pt(jets[hJidx[0]],self.isMC,sysVar=-1),jets[hJidx[0]].eta,jets[hJidx[0]].phi,jets[hJidx[0]].mass)
            hj2down.SetPtEtaPhiM(self.pt(jets[hJidx[1]],self.isMC,sysVar=-1),jets[hJidx[1]].eta,jets[hJidx[1]].phi,jets[hJidx[1]].mass)
            hbbdown = hj1down + hj2down
            self.out.fillBranch("H_pt_JERDown",hbbdown.Pt())
            self.out.fillBranch("H_phi_JERDown",hbbdown.Phi())
            self.out.fillBranch("H_eta_JERDown",hbbdown.Eta())
            self.out.fillBranch("H_mass_JERDown",hbbdown.M())
            
            ## calculate separately selected indices using CMVA, although keep in mind this is already the
            ## default for 2016
            hJetsCMVA = sorted(jetsForHiggs, key = lambda jet : jet.btagCMVA, reverse=True)[0:2]
            hJidxCMVA = [jets.index(x) for x in hJetsCMVA]
            self.out.fillBranch("hJidxCMVA",hJidxCMVA)
            
            ## Save a few basic reco. H kinematics (from CMVA)
            hj1cmva = ROOT.TLorentzVector()
            hj2cmva = ROOT.TLorentzVector()
            hj1cmva.SetPtEtaPhiM(self.pt(jets[hJidxCMVA[0]],self.isMC),jets[hJidxCMVA[0]].eta,jets[hJidxCMVA[0]].phi,jets[hJidxCMVA[0]].mass)
            hj2cmva.SetPtEtaPhiM(self.pt(jets[hJidxCMVA[1]],self.isMC),jets[hJidxCMVA[1]].eta,jets[hJidxCMVA[1]].phi,jets[hJidxCMVA[1]].mass)
            hbbcmva = hj1cmva + hj2cmva
            self.out.fillBranch("HCMVA_pt",hbbcmva.Pt())
            self.out.fillBranch("HCMVA_phi",hbbcmva.Phi())
            self.out.fillBranch("HCMVA_eta",hbbcmva.Eta())
            self.out.fillBranch("HCMVA_mass",hbbcmva.M())

            ## try to recover FSR
            jetsFromFSR = []
            for ijet in xrange(len(jets)):
                if ijet == hJidx[0] or ijet == hJidx[1]: continue
                jet = jets[ijet]
                if self.pt(jet,self.isMC,noReg=True)>20 and abs(jet.eta)<3.0 and jet.puId>0 and jet.jetId>0 and jet.lepFilter:
                   if min(deltaR(jet,jets[hJidx[0]]),deltaR(jet,jets[hJidx[1]])) < 0.8:
                       jetsFromFSR.append(jet)
            HFSR = hbb
            HFSRUp = hbbup
            HFSRDown = hbbdown
            for jet in jetsFromFSR:
                fsrJetToAdd = ROOT.TLorentzVector()
                fsrJetToAdd.SetPtEtaPhiM(jet.Pt,jet.eta,jet.phi,jet.mass)
                HFSR = HFSR + fsrJetToAdd
                HFSRUp = HFSRUp + fsrJetToAdd
                HFSRDown = HFSRDown + fsrJetToAdd
            self.out.fillBranch("HFSR_pt",HFSR.Pt())
            self.out.fillBranch("HFSR_phi",HFSR.Phi())
            self.out.fillBranch("HFSR_eta",HFSR.Eta())
            self.out.fillBranch("HFSR_mass",HFSR.M())
            self.out.fillBranch("HFSR_pt_JERUp",HFSRUp.Pt())
            self.out.fillBranch("HFSR_phi_JERUp",HFSRUp.Phi())
            self.out.fillBranch("HFSR_eta_JERUp",HFSRUp.Eta())
            self.out.fillBranch("HFSR_mass_JERUp",HFSRUp.M())
            self.out.fillBranch("HFSR_pt_JERDown",HFSRDown.Pt())
            self.out.fillBranch("HFSR_phi_JERDown",HFSRDown.Phi())
            self.out.fillBranch("HFSR_eta_JERDown",HFSRDown.Eta())
            self.out.fillBranch("HFSR_mass_JERDown",HFSRDown.M())
            


            ## Compute soft activity vetoing Higgs jets
            #find signal footprint
            toVeto = hJets + vLeptons
            matchedSAJets=self.matchSoftActivity(toVeto,sa)
            #matchedSAJets=self.matchSoftActivityFSR(hJets[0],hJets[1],sa)
            # update SA variables 
            softActivityJetHT=event.SoftActivityJetHT-sum([x.pt for x in matchedSAJets])
            self.out.fillBranch("SA_Ht",softActivityJetHT)

            matchedSAJetsPt5=[x for x in matchedSAJets if x.pt>5]
            softActivityJetNjets5=event.SoftActivityJetNjets5-len(matchedSAJetsPt5)
            self.out.fillBranch("SA5",softActivityJetNjets5)
            
        else:
            self.out.fillBranch("hJidx",[-1,-1])
            self.out.fillBranch("H_pt",-1)
            self.out.fillBranch("H_phi",-1)
            self.out.fillBranch("H_eta",-1)
            self.out.fillBranch("H_mass",-1)
            self.out.fillBranch("H_pt_JERUp",-1)
            self.out.fillBranch("H_phi_JERUp",-1)
            self.out.fillBranch("H_eta_JERUp",-1)
            self.out.fillBranch("H_mass_JERUp",-1)
            self.out.fillBranch("H_pt_JERDown",-1)
            self.out.fillBranch("H_phi_JERDown",-1)
            self.out.fillBranch("H_eta_JERDown",-1)
            self.out.fillBranch("H_mass_JERDown",-1)
            self.out.fillBranch("hJidxCMVA",[-1,-1])
            self.out.fillBranch("HCMVA_pt",-1)
            self.out.fillBranch("HCMVA_phi",-1)
            self.out.fillBranch("HCMVA_eta",-1)
            self.out.fillBranch("HCMVA_mass",-1)
            self.out.fillBranch("HFSR_pt",-1)
            self.out.fillBranch("HFSR_phi",-1)
            self.out.fillBranch("HFSR_eta",-1)
            self.out.fillBranch("HFSR_mass",-1)
            self.out.fillBranch("HFSR_pt_JERUp",-1)
            self.out.fillBranch("HFSR_phi_JERUp",-1)
            self.out.fillBranch("HFSR_eta_JERUp",-1)
            self.out.fillBranch("HFSR_mass_JERUp",-1)
            self.out.fillBranch("HFSR_pt_JERDown",-1)
            self.out.fillBranch("HFSR_phi_JERDown",-1)
            self.out.fillBranch("HFSR_eta_JERDown",-1)
            self.out.fillBranch("HFSR_mass_JERDown",-1)
            self.out.fillBranch("SA_Ht",-1)
            self.out.fillBranch("SA5",-1)
   
        ## indices for Hjets and soft activity (?)
        fatjetsForHiggs = [x for x in fatjets if x.lepFilter and x.jetId>0 and x.Pt>250 and x.Msoftdrop>40 and abs(x.eta)<2.5]
        if (len(fatjetsForHiggs) >= 1):

            jh = sorted(fatjetsForHiggs, key = lambda jet : jet.Pt, reverse=True)
            pt_idx = fatjets.index(jh[0])
            jh = sorted(fatjetsForHiggs, key = lambda jet : jet.Msoftdrop, reverse=True)
            msd_idx = fatjets.index(jh[0])
            jh = sorted(fatjetsForHiggs, key = lambda jet : jet.btagHbb, reverse=True)
            hbb_idx = fatjets.index(jh[0])
            self.out.fillBranch("Pt_fjidx",pt_idx)
            self.out.fillBranch("Msd_fjidx",msd_idx)
            self.out.fillBranch("Hbb_fjidx",hbb_idx)

            ## SA leading pt
            toVeto = [fatjets[pt_idx]]
            toVeto.extend(vLeptons)
            matchedSAJets=self.matchSoftActivity(toVeto,sa,0.8)
            matchedSAJetsPt5=[x for x in matchedSAJets if x.pt>5]
            softActivityJetHT=event.SoftActivityJetHT-sum([x.pt for x in matchedSAJets])
            self.out.fillBranch("SAptfj_HT",softActivityJetHT)
            softActivityJetNjets5=event.SoftActivityJetNjets5-len(matchedSAJetsPt5)
            self.out.fillBranch("SAptfj5",softActivityJetNjets5)

            ## SA leading mass
            toVeto = [fatjets[msd_idx]]
            toVeto.extend(vLeptons)
            matchedSAJets=self.matchSoftActivity(toVeto,sa,0.8)
            matchedSAJetsPt5=[x for x in matchedSAJets if x.pt>5]
            softActivityJetHT=event.SoftActivityJetHT-sum([x.pt for x in matchedSAJets])
            self.out.fillBranch("SAmfj_HT",softActivityJetHT)
            softActivityJetNjets5=event.SoftActivityJetNjets5-len(matchedSAJetsPt5)
            self.out.fillBranch("SAmfj5",softActivityJetNjets5)

            ## SA leading mass
            toVeto = [fatjets[hbb_idx]]
            toVeto.extend(vLeptons)
            matchedSAJets=self.matchSoftActivity(toVeto,sa,0.8)
            matchedSAJetsPt5=[x for x in matchedSAJets if x.pt>5]
            softActivityJetHT=event.SoftActivityJetHT-sum([x.pt for x in matchedSAJets])
            self.out.fillBranch("SAhbbfj_HT",softActivityJetHT)
            softActivityJetNjets5=event.SoftActivityJetNjets5-len(matchedSAJetsPt5)
            self.out.fillBranch("SAhbbfj5",softActivityJetNjets5)

        else:
            self.out.fillBranch("Pt_fjidx",-1)
            self.out.fillBranch("Msd_fjidx",-1)
            self.out.fillBranch("Hbb_fjidx",-1)
            self.out.fillBranch("SAptfj_HT",-1)
            self.out.fillBranch("SAptfj5",-1)
            self.out.fillBranch("SAmfj_HT",-1)
            self.out.fillBranch("SAmfj5",-1)
            self.out.fillBranch("SAhbbfj_HT",-1)
            self.out.fillBranch("SAhbbfj5",-1)

        FatJet_FlavourComposition=[1]*len(fatjets)
        FatJet_HiggsProducts=[False]*len(fatjets)
        FatJet_ZProducts=[False]*len(fatjets)
        FatJet_WProducts=[False]*len(fatjets)

        if self.isMC:
            ##First add gen information
            
            ##Store number of gen b hadrons with status 2:
            nGenStatus2bHad = 0
            for genP in genParticles:
                if genP.status==2:
                    if ( (abs(genP.pdgId)/100)%10==5 or (abs(genP.pdgId)/1000)%10==5):
                        nGenStatus2bHad=nGenStatus2bHad+1
            self.out.fillBranch("nGenStatus2bHad",nGenStatus2bHad)
            ##Info on b-quarks from top or Higgs decay
            nGenTop = 0;
            genb_id_1=-999;
            genb_id_2=-999;
            genBParts = [];
            genBPartsShallow = [];
            GenBJ1 = ROOT.TLorentzVector()
            GenBJ2 = ROOT.TLorentzVector()
            GenBJ1_genjet = ROOT.TLorentzVector()
            GenBJ2_genjet = ROOT.TLorentzVector()
            for genP in genParticles:
                if(abs(genP.pdgId)==6 and (genP.statusFlags&8192)==8192): nGenTop=nGenTop+1
                if(abs(genP.pdgId)==5 and genP.genPartIdxMother > -1 and (abs(genParticles[genP.genPartIdxMother].pdgId)==25 or abs(genParticles[genP.genPartIdxMother].pdgId)==6)):
                    genBParts.append(copy.copy(genP))
                    genBPartsShallow.append(genP)#Just to extract the IDs in the loop again
            genBParts.sort(key=lambda x:x.pt,reverse=True)
            self.out.fillBranch("nGenTop",nGenTop)
            if len(genBPartsShallow) > 1:
                for i in range(0,len(genParticles)):
                    if(genBPartsShallow[0]==genParticles[i]): genb_id_1=i
                    if(genBPartsShallow[1]==genParticles[i]): genb_id_2=i
            self.out.fillBranch("GenBJ1_index",genb_id_1)
            self.out.fillBranch("GenBJ2_index",genb_id_2)
            for part in genBParts:
                part.mass=4.2
            if len (genBParts) > 1:
                GenBJ1.SetPtEtaPhiM(genBParts[0].pt,genBParts[0].eta,genBParts[0].phi,genBParts[0].mass)
                GenBJ2.SetPtEtaPhiM(genBParts[1].pt,genBParts[1].eta,genBParts[1].phi,genBParts[1].mass)
                self.out.fillBranch("GenBJ1_pt",genBParts[0].pt)
                self.out.fillBranch("GenBJ1_eta",genBParts[0].eta)
                self.out.fillBranch("GenBJ1_phi",genBParts[0].phi)
                self.out.fillBranch("GenBJ1_mass",genBParts[0].mass)
                self.out.fillBranch("GenBJ2_pt",genBParts[1].pt)
                self.out.fillBranch("GenBJ2_eta",genBParts[1].eta)
                self.out.fillBranch("GenBJ2_phi",genBParts[1].phi)
                self.out.fillBranch("GenBJ2_mass",genBParts[1].mass)
                self.out.fillBranch("GenBJJ_pt",(GenBJ1+GenBJ2).Pt())
                self.out.fillBranch("GenBJJ_eta",(GenBJ1+GenBJ2).Eta())
                self.out.fillBranch("GenBJJ_phi",(GenBJ1+GenBJ2).Phi())
                self.out.fillBranch("GenBJJ_mass",(GenBJ1+GenBJ2).M())
                self.out.fillBranch("GenBJJ_dEta",abs(genBParts[0].eta-genBParts[1].eta))
                self.out.fillBranch("GenBJJ_dR", deltaR(genBParts[0],genBParts[1]))
                self.out.fillBranch("GenBJJ_dPhi", deltaPhi(genBParts[0],genBParts[1]))
            else:
                self.out.fillBranch("GenBJ1_pt",-999)
                self.out.fillBranch("GenBJ1_eta",-999)
                self.out.fillBranch("GenBJ1_phi",-999)
                self.out.fillBranch("GenBJ1_mass",-999)
                self.out.fillBranch("GenBJ2_pt",-999)
                self.out.fillBranch("GenBJ2_eta",-999)
                self.out.fillBranch("GenBJ2_phi",-999)
                self.out.fillBranch("GenBJ2_mass",-999)
                self.out.fillBranch("GenBJJ_pt",-999)
                self.out.fillBranch("GenBJJ_eta",-999)
                self.out.fillBranch("GenBJJ_phi",-999)
                self.out.fillBranch("GenBJJ_mass",-999)
                self.out.fillBranch("GenBJJ_dEta",-999)
                self.out.fillBranch("GenBJJ_dR",-999)
                self.out.fillBranch("GenBJJ_dPhi",-999)

            ##Closest gen jet to b-parton:
            genbj1_genjet_index=-999
            genbj2_genjet_index=-999
            if len(genBParts) > 1:
                copiedGenJets_firstb = copy.copy(genJets)
                copiedGenJets_secondb = copy.copy(genJets)
                copiedGenJets_firstb.sort(key=lambda x:deltaR(x,genBParts[0])) 
                copiedGenJets_secondb.sort(key=lambda x:deltaR(x,genBParts[1])) 
                sameJetIndex=False
                if len(copiedGenJets_firstb)>0:
                    if copiedGenJets_firstb[0]==copiedGenJets_secondb[0]: sameJetIndex=True
                if len(copiedGenJets_firstb) > 1:
                    for i in range(0,len(genJets)):
                        if copiedGenJets_firstb[0]==genJets[i]:genbj1_genjet_index=i
                        if sameJetIndex:
                            if copiedGenJets_secondb[1]==genJets[i]:genbj2_genjet_index=i
                        else :
                            if copiedGenJets_secondb[0]==genJets[i]:genbj2_genjet_index=i
                else:
                    for i in range(0,len(genJets)):
                        if copiedGenJets_firstb[0]==genJets[i]:genbj1_genjet_index=i
            self.out.fillBranch("GenBJ1_genjet_index",genbj1_genjet_index)
            self.out.fillBranch("GenBJ2_genjet_index",genbj2_genjet_index)
            if genbj2_genjet_index>-1:
                self.out.fillBranch("GenBJ1_genjet_pt",genJets[genbj1_genjet_index].pt)
                self.out.fillBranch("GenBJ1_genjet_mass",genJets[genbj1_genjet_index].mass)
                self.out.fillBranch("GenBJ2_genjet_pt",genJets[genbj2_genjet_index].pt)
                self.out.fillBranch("GenBJ2_genjet_mass",genJets[genbj2_genjet_index].mass)
                GenBJ1_genjet.SetPtEtaPhiM(genJets[genbj1_genjet_index].pt,genJets[genbj1_genjet_index].eta,genJets[genbj1_genjet_index].phi,genJets[genbj1_genjet_index].mass) 
                GenBJ2_genjet.SetPtEtaPhiM(genJets[genbj2_genjet_index].pt,genJets[genbj2_genjet_index].eta,genJets[genbj2_genjet_index].phi,genJets[genbj2_genjet_index].mass) 
                self.out.fillBranch("GenBJJ_genjet_pt",(GenBJ1_genjet+GenBJ2_genjet).Pt())
                self.out.fillBranch("GenBJJ_genjet_mass",(GenBJ1_genjet+GenBJ2_genjet).M())
            elif genbj1_genjet_index>-1:
                self.out.fillBranch("GenBJ1_genjet_pt",genJets[genbj1_genjet_index].pt)
                self.out.fillBranch("GenBJ1_genjet_mass",genJets[genbj1_genjet_index].mass)
                self.out.fillBranch("GenBJ2_genjet_pt", -999)
                self.out.fillBranch("GenBJ2_genjet_mass",-999)
            else:
                self.out.fillBranch("GenBJ1_genjet_pt", -999)
                self.out.fillBranch("GenBJ1_genjet_mass",-999)
                self.out.fillBranch("GenBJ2_genjet_pt", -999)
                self.out.fillBranch("GenBJ2_genjet_mass",-999)

               
         
            ##Information on closest gen lepton to b-jets
            genLeptons = [];
            genLep1 = ROOT.TLorentzVector()
            genLep2 = ROOT.TLorentzVector()
            genTop1 = ROOT.TLorentzVector()
            genTop2 = ROOT.TLorentzVector()
            genindex_1=-1
            genindex_2=-1
            if len(genBParts) > 1:
                for genP in genParticles:
                    if( (abs(genP.pdgId)==11 or abs(genP.pdgId)==13) and genP.status==1 and ( (genP.statusFlags&1)==1 or (genP.statusFlags&32)==32)):
                        genLeptons.append(genP)
             
                genLeptonsCp = copy.copy(genLeptons)
                genLeptons.sort(key=lambda x:deltaR(x,genBParts[0])) 
                genLeptonsCp.sort(key=lambda x:deltaR(x,genBParts[1])) 
                sameParticleIndex=False
                if len(genLeptons) > 0 :
                    if genLeptons[0]==genLeptonsCp[0]: sameParticleIndex=True
                if len(genLeptons) > 1:
                    for i in range(0,len(genParticles)):
                        if( genParticles[i]==genLeptons[0]): genindex_1=i
                        if sameParticleIndex:
                            if(genParticles[i]==genLeptonsCp[1]): genindex_2=i
                        else:
                            if(genParticles[i]==genLeptonsCp[0]): genindex_2=i
                elif len(genLeptons) > 0:
                    for i in range(0,len(genParticles)):
                        if( genParticles[i]==genLeptons[0]): genindex_1=i
            self.out.fillBranch("GenLepIndex1",genindex_1) 
            self.out.fillBranch("GenLepIndex2",genindex_2) 
            if genindex_1>-1:
                if (abs(genParticles[genindex_1].pdgId)==11): part_mass=0.000511; 
                if (abs(genParticles[genindex_1].pdgId)==13): part_mass=0.105; 
                genLep1.SetPtEtaPhiM(genParticles[genindex_1].pt,genParticles[genindex_1].eta,genParticles[genindex_1].phi,part_mass)
                genTop1 = genLep1+GenBJ1
                self.out.fillBranch("GenLep_GenBJ1_dR",deltaR(genParticles[genindex_1],genBParts[0]))
                self.out.fillBranch("GenLep_GenBJ1_dEta",abs(genParticles[genindex_1].eta-genBParts[0].eta))
                self.out.fillBranch("GenLep_GenBJ1_dPhi",deltaPhi(genParticles[genindex_1],genBParts[0]))
                self.out.fillBranch("GenTop1_mass",genTop1.M())
            else :
                self.out.fillBranch("GenLep_GenBJ1_dR",-999)
                self.out.fillBranch("GenLep_GenBJ1_dEta",-999)
                self.out.fillBranch("GenLep_GenBJ1_dPhi",-999)
                self.out.fillBranch("GenTop1_mass",-999)
            if genindex_2>-1:
                if (abs(genParticles[genindex_2].pdgId)==11): part_mass=0.000511; 
                if (abs(genParticles[genindex_2].pdgId)==13): part_mass=0.105; 
                genLep2.SetPtEtaPhiM(genParticles[genindex_2].pt,genParticles[genindex_2].eta,genParticles[genindex_2].phi,part_mass)
                genTop2 = genLep2+GenBJ2
                self.out.fillBranch("GenLep_GenBJ2_dR",deltaR(genParticles[genindex_2],genBParts[1]))
                self.out.fillBranch("GenLep_GenBJ2_dEta",abs(genParticles[genindex_2].eta-genBParts[1].eta))
                self.out.fillBranch("GenLep_GenBJ2_dPhi",deltaPhi(genParticles[genindex_2],genBParts[1]))
                self.out.fillBranch("GenTop2_mass",genTop2.M())
            else :
                self.out.fillBranch("GenLep_GenBJ2_dR",-999)
                self.out.fillBranch("GenLep_GenBJ2_dEta",-999)
                self.out.fillBranch("GenLep_GenBJ2_dPhi",-999)
                self.out.fillBranch("GenTop2_mass",-999)

            #Vector boson reconstruction
            vbosons = []
            hbosons = []
            wbosons = []
            wbosons_lep=[]
            vleptons=[]
            vneutrinos=[]
            nvbosons=-999
            vboson_pdgid=-999
            vboson_lead_pt=-999
            vboson_lead_eta=-999
            vboson_lead_phi=-999
            vneutrino_pt=-999
            for genP in genParticles:
                if ( abs(genP.pdgId)==24 and (genP.statusFlags&8192)==8192): 
                    vbosons.append(genP)
                    wbosons.append(genP)
                if ( abs(genP.pdgId)==23 and (genP.statusFlags&8192)==8192): vbosons.append(genP)
                if ( abs(genP.pdgId)==25 and (genP.statusFlags&8192)==8192): hbosons.append(genP)
            if len(wbosons) > 0:
                for genP in genParticles:
                    if( abs(genP.pdgId)>=11 and abs(genP.pdgId)<=16 and genP.genPartIdxMother>-1):
                        if( abs(genParticles[genP.genPartIdxMother].pdgId)==24 and (genParticles[genP.genPartIdxMother].statusFlags&8192)==8192):
                            wbosons_lep.append(genP.genPartIdxMother) 

            wbosons_lep = list(set(wbosons_lep)) #To keep only unique mother indices
            vbosons.sort(key=lambda x:x.pt, reverse=True)
            for genP in genParticles:
               if( ( (abs(genP.pdgId)==12 or abs(genP.pdgId)==14 or abs(genP.pdgId)==16) and genP.status==1 and (genP.statusFlags&1)==1 ) or ( abs(genP.pdgId)==16 and (genP.statusFlags&1)==1 and (genP.statusFlags&2)==2)):
                 vneutrinos.append(genP)
            vneutrinos.sort(key=lambda x:x.pt,reverse=True)
            if( len(vneutrinos) > 0):
                vneutrino_pt = vneutrinos[0].pt
            if( len(vbosons)==0 and self.isVjets):
                #If V+jets sample but no vector bosons in the event record, reconstruct them from the leptons in the event: either prompt final state leptons or prompt taus which satisfy isDecayedLeptonHadron.
                for genP in genParticles:
                   if( (abs(genP.pdgId)>=11 and abs(genP.pdgId)<=16 and genP.status==1 and (genP.statusFlags&1)==1 ) or ( abs(genP.pdgId)==15 and (genP.statusFlags&1)==1 and (genP.statusFlags&2)==2)):
                     vleptons.append(genP)
 
                vleptons.sort(key=lambda x:x.pt,reverse=True)
                if len(vleptons)>1:
                    lep1 = ROOT.TLorentzVector()
                    lep2 = ROOT.TLorentzVector()
                    lep1.SetPtEtaPhiM(vleptons[0].pt,vleptons[0].eta,vleptons[0].phi,vleptons[0].mass)
                    lep2.SetPtEtaPhiM(vleptons[1].pt,vleptons[1].eta,vleptons[1].phi,vleptons[1].mass)
                    nvbosons = 1
                    vboson_pdgid = 23
                    vboson_lead_pt = (lep1+lep2).Pt()
                    vboson_lead_eta = (lep1+lep2).Eta()
                    vboson_lead_phi = (lep1+lep2).Phi()
            
            if( len(vbosons) > 0):
                nvbosons = len(vbosons)
                vboson_pdgid = vbosons[0].pdgId
                vboson_lead_pt = vbosons[0].pt 
                vboson_lead_eta = vbosons[0].eta 
                vboson_lead_phi = vbosons[0].phi 
                #The following is a hack that makes sure the EWK/QCD reweighting does not get applied to diboson samples in which only one of the vector bosons is in the event record.
                if (len(vbosons) ==1 and not self.isVjets and len(hbosons)==0): nvbosons=2 

            self.out.fillBranch("nW", len(wbosons))
            self.out.fillBranch("nWlep", len(wbosons_lep))
            self.out.fillBranch("nGenVBosons",nvbosons)
            self.out.fillBranch("LeadGenVBoson_pt",vboson_lead_pt)
            self.out.fillBranch("LeadGenVBoson_eta",vboson_lead_eta)
            self.out.fillBranch("LeadGenVBoson_phi",vboson_lead_phi)
            self.out.fillBranch("LeadGenVBoson_pdgId", vboson_pdgid)
            self.out.fillBranch("LeadNeutrinoFromGenVBoson_pt", vneutrino_pt) #Basically just for leading neutrino in W(ln).

            #do flavour composition here
            for fatjet in fatjets:

                numBHadrons=0
                numCHadrons=0
                H_decay=False
                W_decay=False
                Z_decay=False

                for genP in genParticles:
                    if deltaR(fatjet,genP)<0.8:
                        id_at=max((abs(genP.pdgId)/1000) % 10,(abs(genP.pdgId)/100) % 10)

                        if (id_at==4 or id_at==5):
                            lastinChain=True
                            if genP.genPartIdxMother<0: id_mother=-1
                            else: id_mother=max((abs(genParticles[genP.genPartIdxMother].pdgId)/1000) % 10,(abs(genParticles[genP.genPartIdxMother].pdgId)/100) % 10)

                            if (id_mother==4 or id_mother==5):
                                lastinChain=False

                            if lastinChain:
                                if id_at==4: numCHadrons=numCHadrons+1
                                if id_at==5: numBHadrons=numBHadrons+1

                                #subjet matching?

                                #s1=fatjet.subJetIdx1
                                #s2=fatjet.subJetIdx2

                                #if s1>-1:
                                    #if deltaR(subjets[s1],genP)<0.3: print "it's in subjet 1"
                                #if s2>-1:
                                    #if deltaR(subjets[s2],genP)<0.3: print "it's in subjet 2"
                        if (abs(genP.pdgId)==24): W_decay=True
                        if (genP.pdgId==23): Z_decay=True
                        if (genP.pdgId==25): H_decay=True

                        #status flags: maybe we want ot check a few bits in the status before assigning True?

                        #if (genP.pdgId==25): 
                            #print "Higgs boson",genP.pdgId, genP.status, genP.statusFlags,genP.genPartIdxMother
                            #for bit in range(15):
                                #if ((genP.statusFlags>>bit)&1) : print self.statusFlags_dict(bit)

                if numBHadrons>=2:
                 FatJet_FlavourComposition[fatjets.index(fatjet)]=55
                elif numBHadrons==1 and numCHadrons>=1:
                 FatJet_FlavourComposition[fatjets.index(fatjet)]=54
                elif numBHadrons==1 and numCHadrons==0:
                 FatJet_FlavourComposition[fatjets.index(fatjet)]=5
                elif numBHadrons==0 and numCHadrons==2:
                 FatJet_FlavourComposition[fatjets.index(fatjet)]=44
                elif numBHadrons==0 and numCHadrons==1:
                 FatJet_FlavourComposition[fatjets.index(fatjet)]=4

                FatJet_HiggsProducts[fatjets.index(fatjet)]=H_decay
                FatJet_ZProducts[fatjets.index(fatjet)]=Z_decay
                FatJet_WProducts[fatjets.index(fatjet)]=W_decay

            self.out.fillBranch("FatJet_FlavourComposition",FatJet_FlavourComposition)
            self.out.fillBranch("FatJet_HiggsProducts",FatJet_HiggsProducts)
            self.out.fillBranch("FatJet_ZProducts",FatJet_ZProducts)
            self.out.fillBranch("FatJet_WProducts",FatJet_WProducts)
        else:
            # data event, fill with default values
            self.out.fillBranch("FatJet_FlavourComposition",FatJet_FlavourComposition)
            self.out.fillBranch("FatJet_HiggsProducts",FatJet_HiggsProducts)
            self.out.fillBranch("FatJet_ZProducts",FatJet_ZProducts)
            self.out.fillBranch("FatJet_WProducts",FatJet_WProducts)

        return True
                

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

vhbb2016 = lambda : VHbbProducer(True,"2016") 
vhbb2016_vjets = lambda : VHbbProducer(True,"2016",isVjets=True) 
vhbb2016_data = lambda : VHbbProducer(False,"2016") 
vhbb2017 = lambda : VHbbProducer(True,"2017") 
vhbb2017_vjets = lambda : VHbbProducer(True,"2017",isVjets=True) 
vhbb2017_data = lambda : VHbbProducer(False,"2017") 
vhbb2018 = lambda : VHbbProducer(True,"2018") 
vhbb2018_vjets = lambda : VHbbProducer(True,"2018",isVjets=True) 
vhbb2018_data = lambda : VHbbProducer(False,"2018") 
