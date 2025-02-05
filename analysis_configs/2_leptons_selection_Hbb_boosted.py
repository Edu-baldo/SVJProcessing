import os
import numpy as np
import awkward as ak

from skimmer import skimmer_utils
from utils.awkward_array_utilities import as_type
import analysis_configs.triggers_Hbb as trg
import utils.variables_computation.event_variables as event_vars
from analysis_configs.met_filters import met_filters_nanoaod as met_filters
from analysis_configs import sequences_2_leptons_selection_Hbb_boosted as sequences


def process(events, cut_flow, year, primary_dataset="", pn_tagger=False, **kwargs):
    """Hbb 0 leptons pre-selection boosted categorie."""

    # Trigger event selection
    triggers = getattr(trg, f"double_lepton_2017")
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

    #Double lepton filter
    # isWmumu
    if len(events.Muon_pt) == 2:  
        events = sequences.filter_isWmumu(events)
        skimmer_utils.update_cut_flow(cut_flow, "isWmumu", events)

    # isWee
    elif len(events.Electron_pt) == 2:  
        events = sequences.filter_isWemu(events)
        skimmer_utils.update_cut_flow(cut_flow, "isWee", events)

    # Requiring at least 1 good Jets at a specific event
    filter = ak.count(events.Jet_pt[events.Jet_isGood], axis=1) > 1
    events = events[filter]
    skimmer_utils.update_cut_flow(cut_flow, "nJetsAK4Gt2", events)
    
    # Requiring at least 1 good FatJets at a specific event
    filter = ak.count(events.FatJet_pt[events.FatJet_isGood], axis=1) >= 1
    events = events[filter]
    skimmer_utils.update_cut_flow(cut_flow, "nJetsAK8Gt2", events)

 
    events = sequences.add_analysis_branches(events)
    events = sequences.remove_collections(events)

    return events, cut_flow

