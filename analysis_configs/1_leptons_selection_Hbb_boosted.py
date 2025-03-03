import os
import numpy as np
import awkward as ak

from skimmer import skimmer_utils
from utils.awkward_array_utilities import as_type
import analysis_configs.triggers_Hbb as trg
import utils.variables_computation.event_variables as event_vars
from analysis_configs.met_filters import met_filters_nanoaod as met_filters
from analysis_configs import sequences_1_leptons_selection_Hbb_boosted as sequences


def process(events, cut_flow, year, primary_dataset="", pn_tagger=False, **kwargs):
    """Hbb 0 leptons pre-selection boosted categorie."""

    # Trigger event selection
    triggers = getattr(trg, f"single_lepton_2017")
    events = skimmer_utils.apply_trigger_cut(events, triggers)
    skimmer_utils.update_cut_flow(cut_flow, "Trigger", events)

    # Good jets filters
    events = sequences.apply_good_ak4_jet_filter(events)
    skimmer_utils.update_cut_flow(cut_flow, "GoodJetsAK4", events)

    # Good Fatjet filters
    events = sequences.apply_good_ak8_jet_filter(events)
    skimmer_utils.update_cut_flow(cut_flow, "GoodJetsAK8", events)

    # Apply good electrons filter
    events = sequences.apply_good_electrons_filter(events)
    skimmer_utils.update_cut_flow(cut_flow, "Good electrons filter", events)

    # Apply good muons filter
    events = sequences.apply_good_muons_filter(events)
    skimmer_utils.update_cut_flow(cut_flow, "Good muons filter", events)

    # Adding good objects branches branch already so that it can be used
    # in the rest of the pre-selection
    events = sequences.add_good_ak4_jet_branch(events)
    events = sequences.add_good_ak8_jet_branch(events)
    events = sequences.add_good_electrons_branch(events)
    events = sequences.add_good_muons_branch(events)

    #Single lepton filter
    if len(events) != 0:
        # Filter events with exactly one muon or one electron
        single_lepton_mask = (ak.num(events.Muon_pt) == 1) | (ak.num(events.Electron_pt) == 1)
        events = events[single_lepton_mask]

        # # isWmunu
        # if ak.any(ak.num(events.Muon_pt) == 1): 
        #     filter = sequences.isWmunu(events)
        #     events["isWmunu"] = filter
        #     events = events[filter]
        #     skimmer_utils.update_cut_flow(cut_flow, "isWmunu_filter", events)

        # # isWenu
        # if ak.any(ak.num(events.Electron_pt) == 1): 
        #     filter = sequences.isWenu(events)
        #     events["isWenu"] = filter
        #     events = events[filter]
        #     skimmer_utils.update_cut_flow(cut_flow, "isWenu_filter", events)
    skimmer_utils.update_cut_flow(cut_flow, "Single_Lepton", events)

    # isWmunu
    if ak.any(ak.num(events.Muon_pt) == 1): 
        filter = sequences.isWmunu(events)
        events["isWmunu"] = filter
        events = events[filter]
    skimmer_utils.update_cut_flow(cut_flow, "isWmunu_filter", events)

    # isWenu
    if ak.any(ak.num(events.Electron_pt) == 1): 
        filter = sequences.isWenu(events)
        events["isWenu"] = filter
        events = events[filter]
    skimmer_utils.update_cut_flow(cut_flow, "isWenu_filter", events)

    # Requiring at least 1 good Jets at a specific event
    filter = ak.count(events.Jet_pt[events.Jet_isGood], axis=1) >= 1
    events = events[filter]
    skimmer_utils.update_cut_flow(cut_flow, "nJetsAK4Gt2", events)
    
    # Apply delta phi cut for leptons: dPhi(lep,MET) < 2
    if len(events) != 0:

        met = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
            pt=events.MET_pt,
            phi=events.MET_phi,
        )

        # For electrons
        if ak.any(events.Electron_pt):
            electrons = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
                pt=events.Electron_pt,
                eta=events.Electron_eta,
                phi=events.Electron_phi,
                mass=events.Electron_mass,
            )
            met_electrons = ak.broadcast_arrays(met, electrons)[0]
            delta_phi_electrons = abs(electrons.delta_phi(met_electrons))
            electron_filter = ak.any(delta_phi_electrons < 2, axis=1)
        else:
            electron_filter = False

        # For muons
        if ak.any(events.Muon_pt):
            muons = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
                pt=events.Muon_pt,
                eta=events.Muon_eta,
                phi=events.Muon_phi,
                mass=events.Muon_mass,
            )
            met_muons = ak.broadcast_arrays(met, muons)[0]
            delta_phi_muons = abs(muons.delta_phi(met_muons))
            muon_filter = ak.any(delta_phi_muons < 2, axis=1)
        else:
            muon_filter = False

        filter_deltaphi_lep_met = electron_filter | muon_filter
        filter_deltaphi_lep_met = as_type(filter_deltaphi_lep_met, bool)
        events = events[filter_deltaphi_lep_met]

    skimmer_utils.update_cut_flow(cut_flow, "DeltaPhi(lep, MET) selection", events)

    # Requiring at least 1 good FatJets at a specific event
    filter = ak.count(events.FatJet_pt[events.FatJet_isGood], axis=1) >= 1
    events = events[filter]
    skimmer_utils.update_cut_flow(cut_flow, "nJetsAK8Gt2", events)

    # events = sequences.add_analysis_branches(events)
    # events = sequences.remove_collections(events)

    return events, cut_flow

