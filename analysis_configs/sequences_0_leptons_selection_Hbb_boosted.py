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
    return events

def filter_isZnn(events):

    filtered_events = Decay_Vtype.calculate_vtype(events, vtype_filter=4)
    events = events[filtered_events]

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