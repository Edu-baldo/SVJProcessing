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

    #veto all leptons
    if len(events) != 0:
        veto_electrons = ak.count(events.Electron_pt, axis=1) == 0
        veto_muons = ak.count(events.Muon_pt, axis=1) == 0
        events = events[veto_electrons & veto_muons]
    skimmer_utils.update_cut_flow(cut_flow, "Lepton veto", events)

    # Requiring at least 1 good FatJets at a specific event
    filter = ak.count(events.Jet_pt[events.Jet_isGood], axis=1) >= 1
    events = events[filter]
    skimmer_utils.update_cut_flow(cut_flow, "nJetsAK84Gt2", events)
    
    # Apply delta phi cut
    if len(events) != 0:
        #filter jets with pT > 30 GeV
        jet_filter = events.Jet_pt[events.Jet_isGood] > 30

        # If needed because the selection crashes due to the special ak type
        met = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
            pt=events.MET_pt,
            phi=events.MET_phi,
        )
        jets = skimmer_utils.make_pt_eta_phi_mass_lorentz_vector(
            pt=events.Jet_pt[jet_filter],
            eta=events.Jet_eta[jet_filter],
            phi=events.Jet_phi[jet_filter],
            mass=events.Jet_mass[jet_filter],
        )

        met = ak.broadcast_arrays(met, jets)[0]
        delta_phi = abs(jets.delta_phi(met))
        #require the Deltaphi between MET and any jet to be greater than 0.5
        filter_deltaphi = ak.any(delta_phi > 0.5, axis=1)
        filter_deltaphi = as_type(filter_deltaphi, bool)
        events = events[filter_deltaphi]

    skimmer_utils.update_cut_flow(cut_flow, "DeltaPhi selection", events)

    # Apply pT miss selection
    if len(events) != 0:
        events = events[events.MET_pt > 250]  # Threshold for boosted category
    skimmer_utils.update_cut_flow(cut_flow, "pT miss selection", events)

    # Requiring at least 2 good FatJets at a specific event
    filter = ak.count(events.FatJet_pt[events.FatJet_isGood], axis=1) >= 1
    events = events[filter]
    skimmer_utils.update_cut_flow(cut_flow, "nJetsAK8Gt2", events)

 
    events = sequences.add_analysis_branches(events)
    events = sequences.remove_collections(events)

    return events, cut_flow

