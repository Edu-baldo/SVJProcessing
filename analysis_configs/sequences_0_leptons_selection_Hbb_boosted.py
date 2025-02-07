import numpy as np
import awkward as ak

from skimmer import skimmer_utils
from utils.variables_computation import event_variables as event_vars
from utils.data.triggers import primary_dataset_triggers
from utils.tree_maker.triggers import trigger_table
from analysis_configs import objects_definition_0_lepton_Hbb as obj
from utils import Decay_Vtype 

from utils.Logger import *


def apply_good_ak8_jet_filter(events):
    #analysis_jets = events.FatJet[obj.is_analysis_ak8_jet(events.FatJet)]
    #good_ak8_jets_filter = ak.all(obj.is_good_ak8_jet(analysis_jets), axis=1)
    ak8_jets = ak.zip(
        {
            "pt": events.FatJet_pt,
            "mass": events.FatJet_mass,
            "eta": events.FatJet_eta,
            "phi": events.FatJet_phi,
            "id": events.FatJet_jetId,
        },
        with_name="PtEtaPhiMLorentzVector",
    )
    analysis_jets = ak8_jets[obj.is_analysis_ak8_jet(ak8_jets)] #questo filtra un singolo jet
    good_ak8_jets_filter = ak.all(obj.is_good_ak8_jet(analysis_jets), axis=1)    # "ak.all"verifica che tutti i jet in ciascun evento soddisfano la condizione 

    events = events[good_ak8_jets_filter]
    return events


def apply_good_ak4_jet_filter(events):
    ak4_jets = ak.zip(
        {
            "pt": events.Jet_pt,
            "mass": events.Jet_mass,
            "eta": events.Jet_eta,
            "phi": events.Jet_phi,
            "id": events.Jet_jetId,
        },
        with_name="PtEtaPhiMLorentzVector",
    )
    analysis_jets = ak4_jets[obj.is_analysis_ak4_jet(ak4_jets)]
    good_ak4_jets_filter = ak.all(obj.is_good_ak4_jet(analysis_jets), axis=1)

    events = events[good_ak4_jets_filter]
    return events


def add_good_ak8_jet_branch(events):
    ak8_jets = ak.zip(
        {
            "pt": events.FatJet_pt,
            "mass": events.FatJet_mass,
            "eta": events.FatJet_eta,
            "phi": events.FatJet_phi,
            "id": events.FatJet_jetId,
        },
        with_name="PtEtaPhiMLorentzVector",
    )
    
    is_good_analysis_ak8_jet = (
        obj.is_analysis_ak8_jet(ak8_jets)
        & obj.is_good_ak8_jet(ak8_jets)
    )

    #add new branch to the events
    events["FatJet_isGood"] = is_good_analysis_ak8_jet

    return events

def add_good_ak4_jet_branch(events):
    ak4_jets = ak.zip(
        {
            "pt": events.Jet_pt,
            "mass": events.Jet_mass,
            "eta": events.Jet_eta,
            "phi": events.Jet_phi,
            "id": events.Jet_jetId,
        },
        with_name="PtEtaPhiMLorentzVector",
    )
    
    is_good_analysis_ak4_jet = (
        obj.is_analysis_ak4_jet(ak4_jets)
        & obj.is_good_ak4_jet(ak4_jets)
    )

    #add new branch to the events
    events["Jet_isGood"] = is_good_analysis_ak4_jet

    return events

def add_analysis_branches(events):

    # Event variables
    good_jets_ak8_lv = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
        pt=events.FatJet_pt[events.FatJet_isGood],
        eta=events.FatJet_eta[events.FatJet_isGood],
        phi=events.FatJet_phi[events.FatJet_isGood],
        mass=events.FatJet_mass[events.FatJet_isGood],
    )
    met = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
            pt=events.MET_pt,
            phi=events.MET_phi,
    )

    #add rt variable
    mt = event_vars.calculate_transverse_mass(good_jets_ak8_lv, met)
    rt = met.pt / mt
    events["RTFatJet"] = rt

    #add mt variable
    events["MT01FatJetMET"] = mt

    #add deltaeta the two leading jets
    events["DeltaEtaJ0J1FatJet"] = event_vars.calculate_delta_eta(good_jets_ak8_lv)
    
    #add minimum delta phi between the MET and the two leading jets
    events["DeltaPhiMinFatJetMET"] = event_vars.calculate_delta_phi_min(good_jets_ak8_lv, met)
    
    return events

def apply_good_electrons_filter(events):
    #creare una funzione dove inserissco sia gli elettroni che i muoni
    #la mia domanda è qual è la mia varibaile? sarà event.lepton? o qualcosa del genere?
    electrons = ak.zip(
        {
            "pt": events.Electron_pt,
            "eta": events.Electron_eta,
            "phi": events.Electron_phi,
            "mass": events.Electron_mass,
            "pfRelIso": events.Electron_pfRelIso03_all,
            "id": events.Electron_cutBased,
        },
        with_name="PtEtaPhiMLorentzVector",
    )
    
    analysis_electrons = electrons[obj.is_analysis_electron(electrons)]
    good_electrons_filter = ak.all(obj.is_veto_electron(analysis_electrons), axis=1)
    events = events[good_electrons_filter]

    return events 

def apply_good_muons_filter(events):
    muons = ak.zip(
        {
            "pt": events.Muon_pt,
            "eta": events.Muon_eta,
            "phi": events.Muon_phi,
            "mass": events.Muon_mass,
            "pfRelIso": events.Muon_pfRelIso03_all,
            "id": events.Muon_mediumId,
        },
        with_name="PtEtaPhiMLorentzVector",
    )

    analysis_muons = muons[obj.is_analysis_muon(muons)]
    good_muons_filter = ak.all(obj.is_veto_muon(analysis_muons), axis=1)
    events = events[good_muons_filter]

    return events

def add_good_electrons_branch(events):
    electrons = ak.zip(
        {
            "pt": events.Electron_pt,
            "eta": events.Electron_eta,
            "phi": events.Electron_phi,
            "mass": events.Electron_mass,
            "pfRelIso": events.Electron_pfRelIso03_all,
            "id": events.Electron_cutBased,
        },
        with_name="PtEtaPhiMLorentzVector",
    )
    
    is_good_analysis_electron = (
        obj.is_analysis_electron(electrons)
        & obj.is_veto_electron(electrons)
    )

    # Aggiungi il nuovo branch agli eventi
    events["Electron_isGood"] = is_good_analysis_electron

    return events

def add_good_muons_branch(events):
    muons = ak.zip(
        {
            "pt": events.Muon_pt,
            "eta": events.Muon_eta,
            "phi": events.Muon_phi,
            "mass": events.Muon_mass,
            "pfRelIso": events.Muon_pfRelIso03_all,
            "id": events.Muon_mediumId,
        },
        with_name="PtEtaPhiMLorentzVector",
    )
    
    is_good_analysis_muon = (
        obj.is_analysis_muon(muons)
        & obj.is_veto_muon(muons)
    )

    # Aggiungi il nuovo branch agli eventi
    events["Muon_isGood"] = is_good_analysis_muon

    return events

def remove_collections(events):
    events = events[[x for x in events.fields if x != "JetsAK15"]]
    events = events[[x for x in events.fields if x != "GenJetsAK15"]]
    # events = events[[x for x in events.fields if x != "run"]]
    # events = events[[x for x in events.fields if x != "luminosityBlock"]]
    # events = events[[x for x in events.fields if x != "event"]]
    # events = events[[x for x in events.fields if x != "ptoftheHiggsbosonasidentifiedinHTXS"]]
    # events = events[[x for x in events.fields if x != "HTXS_Higgs_pt"]]
    # events = events[[x for x in events.fields if x != "HTXS_Higgs_y"]]
    # events = events[[x for x in events.fields if x != "HTXS_stage1_1_cat_pTjet25GeV"]]
    # events = events[[x for x in events.fields if x != "HTXS_stage1_1_cat_pTjet30GeV"]]
    # events = events[[x for x in events.fields if x != "HTXS_stage1_1_fine_cat_pTjet25GeV"]]
    # events = events[[x for x in events.fields if x != "HTXS_stage1_1_fine_cat_pTjet30GeV"]]
    # events = events[[x for x in events.fields if x != "HTXS_stage_0_cat"]]
    # events = events[[x for x in events.fields if x != "HTXS_stage1_2_cat_pTjet25GeV"]]
    # events = events[[x for x in events.fields if x != "HTXS_stage1_2_cat_pTjet30GeV"]]
    # events = events[[x for x in events.fields if x != "HTXS_stage1_2_fine_cat_pTjet25GeV"]]
    # events = events[[x for x in events.fields if x != "HTXS_stage1_2_fine_cat_pTjet30GeV"]]
    # events = events[[x for x in events.fields if x != "HTXS_stage_0"]]
    # events = events[[x for x in events.fields if x != "HTXS_stage_1_pTjet25"]]
    # events = events[[x for x in events.fields if x != "HTXS_stage_1_pTjet30"]]
    # events = events[[x for x in events.fields if x != "HTXS_njets25"]]
    # events = events[[x for x in events.fields if x != "HTXS_njets30"]]
    # events = events[[x for x in events.fields if x != "nboostedTau"]]  
    # events = events[[x for x in events.fields if x != "btagWeight_CSVV2"]]
    # events = events[[x for x in events.fields if x != "btagWeight_DeepCSVB"]]
    # events = events[[x for x in events.fields if x != "CaloMET_phi"]]
    # events = events[[x for x in events.fields if x != "CaloMET_pt"]]
    # events = events[[x for x in events.fields if x != "CaloMET_sumEt"]]
    # events = events[[x for x in events.fields if x != "ChsMET_phi"]]
    # events = events[[x for x in events.fields if x != "ChsMET_pt"]]
    # events = events[[x for x in events.fields if x != "ChsMET_sumEt"]]
    # events = events[[x for x in events.fields if x != "nCorrT1METJet"]]
    # events = events[[x for x in events.fields if x != "nJetPFCands"]]
    # events = events[[x for x in events.fields if x != "nJetSVs"]]
    # events = events[[x for x in events.fields if x != "nFatJetPFCands"]]
    # events = events[[x for x in events.fields if x != "nFatJetSVs"]]
    # events = events[[x for x in events.fields if x != "nPFCands"]]
    # events = events[[x for x in events.fields if x != "DeepMETResolutionTune_phi"]]
    # events = events[[x for x in events.fields if x != "DeepMETResolutionTune_pt"]]
    # events = events[[x for x in events.fields if x != "DeepMETResponseTune_phi"]]
    # events = events[[x for x in events.fields if x != "DeepMETResponseTune_pt"]]
    # events = events[[x for x in events.fields if x != "DeepMETResponseTune_pt"]]
    # events = events[[x for x in events.fields if x != "nElectron"]]
    # events = events[[x for x in events.fields if x != "nFatJet"]]
    # events = events[[x for x in events.fields if x != "nFsrPhoton"]]
    # events = events[[x for x in events.fields if x != "nGenJetCands"]]
    # events = events[[x for x in events.fields if x != "nGenJetSVs"]]
    # events = events[[x for x in events.fields if x != "nGenVisnGenFatJetCandsTau"]]
    # events = events[[x for x in events.fields if x != "nGenVisTau"]]
    # events = events[[x for x in events.fields if x != "nGenVisTauCands"]]
    # events = events[[x for x in events.fields if x != "nGenVisTauCands_nGenVisTauCands"]]
    # events = events[[x for x in events.fields if x != "nGenFatJetCands"]]
    # events = events[[x for x in events.fields if x != "nGenFatJetSVs"]]
    # events = events[[x for x in events.fields if x != "nGenJetAK8"]]
    # events = events[[x for x in events.fields if x != "nGenJet"]]
    # events = events[[x for x in events.fields if x != "nGenCands"]]
    # events = events[[x for x in events.fields if x != "nGenPart"]]
    # events = events[[x for x in events.fields if x != "nSubGenJetAK8"]]
    # events = events[[x for x in events.fields if x != "Generator_binvar"]]
    # events = events[[x for x in events.fields if x != "Generator_scalePDF"]]
    # events = events[[x for x in events.fields if x != "Generator_weight"]]
    # events = events[[x for x in events.fields if x != "Generator_x1"]]
    # events = events[[x for x in events.fields if x != "Generator_x2"]]
    # events = events[[x for x in events.fields if x != "Generator_xpdf1"]]
    # events = events[[x for x in events.fields if x != "Generator_xpdf2"]]



    return events

def filter_isZnn(events): 

    # filtered_events, event_mask = Decay_Vtype.apply_vtype_mask(events, vtype_filter=4)
    # events["isZnn"] = filtered_events
    # events = events[filtered_events]
    # print("events_mask", event_mask)

    electrons = ak.zip(
        {
            "pt": events.Electron_pt,
            "eta": events.Electron_eta,
            "phi": events.Electron_phi,
            "mass": events.Electron_mass,
            "pfRelIso03_all": events.Electron_pfRelIso03_all,
            "id": events.Electron_cutBased,
            "charge": events.Electron_charge,
        },
        with_name="PtEtaPhiMLorentzVector",
    )

    muons = ak.zip(
        {
            "pt": events.Muon_pt,
            "eta": events.Muon_eta,
            "phi": events.Muon_phi,
            "mass": events.Muon_mass,
            "pfRelIso03_all": events.Muon_pfRelIso03_all,
            "tightId": events.Muon_tightId,
            "charge": events.Muon_charge,
            "dxy": events.Muon_dxy,
            "dz": events.Muon_dz,
        },
        with_name="PtEtaPhiMLorentzVector",
    )

    zElectrons = electrons[
        (electrons.pt > 20)
        & (electrons.id == 1) 
        & (electrons.pfRelIso03_all < 0.15)
    ]
    
    zMuons = muons[
        (muons.pt > 20)
        & (muons.pfRelIso03_all < 0.25)
        & (abs(muons.dxy) < 0.05)
        & (abs(muons.dz) < 0.2)
    ]

    wElectrons = electrons[
        (electrons.pt > 25)
        & (electrons.id == 2)
        & (electrons.pfRelIso03_all < 0.12)
    ]

    wMuons = muons[
        (muons.pt > 25)
        & (muons.tightId >= 1)
        & (muons.pfRelIso03_all < 0.15)
        & (abs(muons.dxy) < 0.05)
        & (abs(muons.dz) < 0.2)
    ]

    Znn_mask = (ak.num(zElectrons) == 0) & (ak.num(zMuons) == 0) & (events.MET_pt > 150)
    event_mask = Znn_mask
    events = events[event_mask] 
    print("events_mask", event_mask)

    return events

def calculate_mht_pt(events):

    # Scalar sum of pt of all jets
    jet_pt_sum = ak.sum(events.Jet_pt[events.Jet_isGood], axis=1)
    
    # Scalar sum of pt of all leptons (muons + electrons)
    lepton_pt_sum = (
        ak.sum(events.Muon_pt[events.Muon_isGood], axis=1) +
        ak.sum(events.Electron_pt[events.Electron_isGood], axis=1)
    )
    
    # Total MHT_pt
    mht_pt = jet_pt_sum + lepton_pt_sum
    return mht_pt