import os
import awkward as ak

from skimmer import skimmer_utils
from utils.awkward_array_utilities import as_type
import analysis_configs.triggers_Hbb as trg
import utils.variables_computation.event_variables as event_vars
from analysis_configs.met_filters import met_filters_nanoaod as met_filters
from analysis_configs import sequences_0_leptons_selection_Hbb_resolved as sequences


def process(events, cut_flow, year, primary_dataset="", pn_tagger=False, **kwargs):
    """SVJ s-channel leptons pre-selection."""

    # Trigger event selection
    triggers = getattr(trg, f"no_lepton_2017")
    events = skimmer_utils.apply_trigger_cut(events, triggers)
    skimmer_utils.update_cut_flow(cut_flow, "Trigger", events)

    # Good jet filters
    events = sequences.apply_good_ak4_jet_filter(events)
    skimmer_utils.update_cut_flow(cut_flow, "GoodJetsAK84", events)


    # Adding JetsAK4_isGood branch already so that it can be used
    # in the rest of the pre-selection
    events = sequences.add_good_ak4_jet_branch(events)

    # Requiring at least 2 good Jets
    filter = ak.count(events.Jet_pt[events.Jet_isGood], axis=1) >= 2
    events = events[filter]
    skimmer_utils.update_cut_flow(cut_flow, "nJetsAK4Gt2", events)

    #apply RT filter (RT = MET over MT)
    if len(events) != 0:
        jets = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
            pt=events.Jet_pt[events.Jet_isGood],
            eta=events.Jet_eta[events.Jet_isGood],
            phi=events.Jet_phi[events.Jet_isGood],
            mass=events.Jet_mass[events.Jet_isGood],
        )
        met = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
            pt=events.MET_pt,
            phi=events.MET_phi,
        )
        mt = event_vars.calculate_transverse_mass(jets, met)
        rt = events.MET_pt / mt
        filter_rt = rt > 0.15
        filter_rt = as_type(filter_rt, bool)   #not needed
        events = events[filter_rt]
    
    skimmer_utils.update_cut_flow(cut_flow, "RT selection", events)


    #apply DeltaEta filter
    if len(events) != 0:
        jets = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
            pt=events.Jet_pt[events.Jet_isGood],
            eta=events.Jet_eta[events.Jet_isGood],
            phi=events.Jet_phi[events.Jet_isGood],
            mass=events.Jet_mass[events.Jet_isGood],
        )
        delta_eta = abs(event_vars.calculate_delta_eta(jets))
        filter_deltaeta = delta_eta < 2.5
        filter_deltaeta = as_type(filter_deltaeta, bool)
        events = events[filter_deltaeta]
    
    skimmer_utils.update_cut_flow(cut_flow, "DeltaEtaj0j1 selection", events)

    #apply MT selection
    if len(events) != 0:
        jets = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
            pt=events.Jet_pt[events.Jet_isGood],
            eta=events.Jet_eta[events.Jet_isGood],
            phi=events.Jet_phi[events.Jet_isGood],
            mass=events.Jet_mass[events.Jet_isGood],
        )
        met = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
            pt=events.MET_pt,
            phi=events.MET_phi,
        )
        mt = event_vars.calculate_transverse_mass(jets, met)
        filter_mt = mt > 1500
        filter_mt = as_type(filter_mt, bool)
        events = events[filter_mt]
    
    skimmer_utils.update_cut_flow(cut_flow, "MT selection", events)

    # MET filter event selection
    events = skimmer_utils.apply_met_filters_cut(events, met_filters)
    skimmer_utils.update_cut_flow(cut_flow, "METFilters", events)
    
    # Delta phi min cut
    if len(events) != 0:
        # If needed because the selection crashes due to the special ak type
        met = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
            pt=events.MET_pt,
            phi=events.MET_phi,
        )
        jets = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
            pt=events.Jet_pt[events.Jet_isGood],
            eta=events.Jet_eta[events.Jet_isGood],
            phi=events.Jet_phi[events.Jet_isGood],
            mass=events.Jet_mass[events.Jet_isGood],
        )

        met = ak.broadcast_arrays(met, jets)[0]
        delta_phi_min = ak.min(abs(jets.delta_phi(met)), axis=1)
        filter_deltaphi = delta_phi_min < 0.8
        # Needed otherwise type is not defined and skim cannot be written
        filter_deltaphi = as_type(filter_deltaphi, bool)
        events = events[filter_deltaphi]

    skimmer_utils.update_cut_flow(cut_flow, "DeltaPhiMin selection", events)
    
    events = sequences.add_analysis_branches(events)
    events = sequences.remove_collections(events)

    return events, cut_flow

