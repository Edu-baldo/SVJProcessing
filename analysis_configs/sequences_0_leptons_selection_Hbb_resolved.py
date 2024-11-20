import numpy as np
import awkward as ak

from skimmer import skimmer_utils
from utils.variables_computation import event_variables as event_vars
from utils.data.triggers import primary_dataset_triggers
from utils.tree_maker.triggers import trigger_table
from analysis_configs import objects_definition_0_lepton_Hbb as obj

from utils.Logger import *


def apply_good_ak4_jet_filter(events):
    #analysis_jets = events.FatJet[obj.is_analysis_ak8_jet(events.FatJet)]
    #good_ak8_jets_filter = ak.all(obj.is_good_ak8_jet(analysis_jets), axis=1)
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
    good_jets_ak4_lv = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
        pt=events.Jet_pt[events.Jet_isGood],
        eta=events.Jet_eta[events.Jet_isGood],
        phi=events.Jet_phi[events.Jet_isGood],
        mass=events.Jet_mass[events.Jet_isGood],
    )
    met = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
            pt=events.MET_pt,
            phi=events.MET_phi,
    )

    #add rt variable
    mt = event_vars.calculate_transverse_mass(good_jets_ak4_lv, met)
    rt = met.pt / mt
    events["RTJet"] = rt

    #add mt variable
    events["MT01JetMET"] = mt

    #add deltaeta the two leading jets
    events["DeltaEtaJ0J1Jet"] = event_vars.calculate_delta_eta(good_jets_ak4_lv)
    
    #add minimum delta phi between the MET and the two leading jets
    events["DeltaPhiMinJetMET"] = event_vars.calculate_delta_phi_min(good_jets_ak4_lv, met)
    
    return events


def remove_collections(events):
    events = events[[x for x in events.fields if x != "JetsAK15"]]
    events = events[[x for x in events.fields if x != "GenJetsAK15"]]
    return events
