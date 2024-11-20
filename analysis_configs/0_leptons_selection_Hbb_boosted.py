import os
import numpy as np
import awkward as ak

from skimmer import skimmer_utils
from utils.awkward_array_utilities import as_type
import analysis_configs.triggers_Hbb as trg
import utils.variables_computation.event_variables as event_vars
from analysis_configs.met_filters import met_filters_nanoaod as met_filters
from analysis_configs import sequences_0_leptons_selection_Hbb_boosted as sequences


def process(events, cut_flow, year, primary_dataset="", pn_tagger=False, **kwargs):
    """Hbb 0 leptons pre-selection boosted categorie."""

    # Trigger event selection
    triggers = getattr(trg, f"no_lepton_2017")
    events = skimmer_utils.apply_trigger_cut(events, triggers)
    skimmer_utils.update_cut_flow(cut_flow, "Trigger", events)

    # Good jet filters
    events = sequences.apply_good_ak8_jet_filter(events)
    skimmer_utils.update_cut_flow(cut_flow, "GoodJetsAK8", events)


    # Adding JetsAK8_isGood branch already so that it can be used
    # in the rest of the pre-selection
    events = sequences.add_good_ak8_jet_branch(events)

    # Requiring at least 2 good FatJets at a specific event
    filter = ak.count(events.FatJet_pt[events.FatJet_isGood], axis=1) >= 1
    events = events[filter]
    skimmer_utils.update_cut_flow(cut_flow, "nJetsAK8Gt2", events)

    #leptons veto
    if len(events) != 0:
    # Apply good electrons filter
        events = sequences.apply_good_electrons_filter(events)
    skimmer_utils.update_cut_flow(cut_flow, "Good electrons filter", events)

    if len(events) != 0:
    # Apply good muons filter
        events = sequences.apply_good_muons_filter(events)
    skimmer_utils.update_cut_flow(cut_flow, "Good muons filter", events)

    # Apply pT miss selection
    if len(events) != 0:
        events = events[events.MET_pt > 250]  # Threshold for boosted category
    skimmer_utils.update_cut_flow(cut_flow, "pT miss selection", events)

    # MET filter event selection
    #events = skimmer_utils.apply_met_filters_cut(events, met_filters)
    #skimmer_utils.update_cut_flow(cut_flow, "METFilters", events)
    
    # Groomed jet mass selection
    if len(events) != 0:
        jets = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
            pt=events.FatJet_pt[events.FatJet_isGood],
            eta=events.FatJet_eta[events.FatJet_isGood],
            phi=events.FatJet_phi[events.FatJet_isGood],
            mass=events.FatJet_msoftdrop[events.FatJet_isGood],  # groomed jet mass
        )
        mass = jets.mass
        filter_mass = (mass > 90) & (mass < 150)
        filter_mass = as_type(filter_mass, bool)
        events = events[filter_mass]
    
    skimmer_utils.update_cut_flow(cut_flow, "Groomed Mass selection", events)


    # Delta phi cut
    if len(events) != 0:
    # jet_pT > 30 GeV
        jets_pt30 = events.jet_pt > 30
    
    # Delta Phi 
        delta_phi = np.abs(events.MET_phi - events.jet_phi[jets_pt30])
        delta_phi = np.where(delta_phi > np.pi, 2 * np.pi - delta_phi, delta_phi)
    
    # Delta Phi > 0.5
        filter_delta_phi = ak.all(delta_phi > 0.5, axis=1)
        events = events[filter_delta_phi]
    
    skimmer_utils.update_cut_flow(cut_flow, "Delta phi min cut", events)

 
    events = sequences.add_analysis_branches(events)
    events = sequences.remove_collections(events)

    return events, cut_flow

