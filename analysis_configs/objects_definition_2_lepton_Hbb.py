import numpy as np 

def is_good_ak8_jet(jets_ak8):
    return jets_ak8.id == 6


def is_analysis_ak8_jet(jets_ak8):
    return (
        (jets_ak8.pt > 250)
        & (abs(jets_ak8.eta) < 2.5)
    )


def is_good_ak4_jet(jets):
    return jets.id == 6


def is_analysis_ak4_jet(jets):
    return (
        (jets.pt > 20)
        & (abs(jets.eta) < 2.5)
    )


def is_analysis_electron(electrons):
    return (
        (electrons.pt > 15)
        & (abs(electrons.eta) < 2.5)
    )


def is_analysis_muon(muons):
    return (
        (muons.pt > 15)
        & (abs(muons.eta) < 2.5)   #id for veto (no leptons) category
    )


def is_veto_electron(electrons):
    return (
        is_analysis_electron(electrons)
        & (abs(electrons.pfRelIso) < 0.1) #? iso cut value and which iso?  
    )


def is_veto_muon(muons):
    return (
        is_analysis_muon(muons)
        & (abs(muons.pfRelIso) < 0.4) #? iso cut value and which iso?
    )


def is_tag_electron(electrons):
    return (
        is_analysis_electron(electrons)
        & (electrons.mediumID == 1)
        & (electrons.iso < 0.1)  # mini-isolation, tight WP
    )


def is_tag_muon(muons):
    return (
        is_analysis_muon(muons)
        & (muons.mediumID == 1)
        & (muons.iso < 0.1)  # mini-isolation, tight WP
    )


def __is_isolated_electron(electrons):
    return electrons.pfRelIso < 0.15


def __is_isolated_muon(muons):
    return muons.pfRelIso < 0.15


def is_cleaning_electron(electrons):
    return (
        is_tag_electron(electrons)
        & __is_isolated_electron(electrons)
    )


def is_cleaning_muon(muons):
    return (
        is_tag_muon(muons)
        & __is_isolated_muon(muons)
    )

