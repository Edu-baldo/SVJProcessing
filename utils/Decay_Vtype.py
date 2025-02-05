import awkward as ak
from utils.Logger import *



import awkward as ak

def calculate_vtype(events, vtype_filter):

    # Calculate the Vtype branch for the event

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

    # Debug prints for the masks
    print(f"Number of zElectrons: {ak.num(zElectrons)}")
    print(f"Number of zMuons: {ak.num(zMuons)}")
    print(f"Number of wElectrons: {ak.num(wElectrons)}")
    print(f"Number of wMuons: {ak.num(wMuons)}")

    # # Define masks for each type of Vtype
    # if ak.all(ak.num(zMuons) >= 2):
    #     Zmm_mask = (ak.num(zMuons) >= 2) & (ak.num(zMuons[:, :2]) == 2) & (zMuons[:, 0].charge * zMuons[:, 1].charge < 0)
    # else:
    #     Zmm_mask = False

    # if ak.all(ak.num(zElectrons) >= 2):
    #     Zee_mask = (ak.num(zElectrons) >= 2) & (ak.num(zElectrons[:, :2]) == 2) & (zElectrons[:, 0].charge * zElectrons[:, 1].charge < 0)
    # else:
    #     Zee_mask = False

    # Wmu_mask = ak.all(ak.num(wMuons) == 1)
    # We_mask = ak.all(ak.num(wElectrons) == 1)
    # Znn_mask = ak.all((ak.num(zElectrons) == 0) & (ak.num(zMuons) == 0) & (events.MET_pt > 150))

    # Create a mask that selects only events with the desired Vtype
    if vtype_filter == 0:
        Zmm_mask = (ak.num(zMuons) >= 2) & (ak.num(zMuons[:, :2]) == 2) & (zMuons[:, 0].charge * zMuons[:, 1].charge < 0)
        event_mask = Zmm_mask
        print(f"Zmm_mask: {Zmm_mask}")
    elif vtype_filter == 1:
        Zee_mask = (ak.num(zElectrons) >= 2) & (ak.num(zElectrons[:, :2]) == 2) & (zElectrons[:, 0].charge * zElectrons[:, 1].charge < 0)
        event_mask = Zee_mask
        print(f"Zee_mask: {Zee_mask}")
    elif vtype_filter == 2:
        Wmu_mask = ak.all(ak.num(wMuons) == 1)
        event_mask = Wmu_mask
        print(f"Wmu_mask: {Wmu_mask}")
    elif vtype_filter == 3:
        We_mask = ak.all(ak.num(wElectrons) == 1)
        event_mask = We_mask
        print(f"We_mask: {We_mask}")
    elif vtype_filter == 4:
        Znn_mask = ak.all((ak.num(zElectrons) == 0) & (ak.num(zMuons) == 0) & (events.MET_pt > 150))
        event_mask = Znn_mask
        print(f"Znn_mask: {Znn_mask}")
    else:
        raise ValueError(f"Vtype {vtype_filter} is not valid.")

    print(f"Event mask: {event_mask}")
    print(f"Number of events before filtering: {len(events)}")

    filtered_events = events[event_mask]

    
    return filtered_events