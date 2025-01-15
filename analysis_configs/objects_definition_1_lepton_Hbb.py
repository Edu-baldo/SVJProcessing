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
        (jets.pt > 25)
        & (abs(jets.eta) < 2.5)
    )


def is_analysis_electron(electrons):
    return (
        (electrons.pt > 30)
        & (abs(electrons.eta) < 2.4)
        & (electrons.id == 2)
    )


def is_analysis_muon(muons):
    return (
        (muons.pt > 25)
        & (abs(muons.eta) < 2.4)   
        & (muons.id == 3)
    )


def is_veto_electron(electrons):
    return (
        is_analysis_electron(electrons)
        & (abs(electrons.pfRelIso) < 0.06) #? iso cut value and which iso?  
    )


def is_veto_muon(muons):
    return (
        is_analysis_muon(muons)
        & (abs(muons.pfRelIso) < 0.06) #? iso cut value and which iso?
    )
