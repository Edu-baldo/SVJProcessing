import os

from tests import helper


def __run_skimmer(input_files, output_file, flags, run_particle_net):
    n_workers = 4
    chunk_size = 10000
    executor = "futures"

    bash_command = f"python {os.environ['SVJ_PROCESSING_ROOT']}/skimmer/skim.py -i {input_files} -o {output_file} -e {executor} -n {n_workers} -c {chunk_size} {flags}"
    if run_particle_net:
        bash_command += " -pn_tagger"
    helper.test_command(bash_command)


def test_execution():
    params_list = [
        (
            [f"root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/signal_production_2Dscans/NTUPLE/Private2DUL18/SVJ_t-channel_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8/0_RA2AnalysisTree.root"],
            "-y 2018 -p analysis_configs.t_channel_pre_selection -pd dummy",
        ),
        (
            [f"root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/signal_production_2Dscans/NTUPLE/Private2DUL18/SVJ_t-channel_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8/0_RA2AnalysisTree.root"],
            "-y 2018 -p analysis_configs.t_channel_pre_selection -pd dummy --variation jer_up",
        ),
        (
            [f"root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/signal_production_2Dscans/NTUPLE/Private2DUL18/SVJ_t-channel_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8/0_RA2AnalysisTree.root"],
            "-y 2018 -p analysis_configs.t_channel_pre_selection -pd dummy --variation jec_up",
        ),
        (
            [f"root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/signal_production_2Dscans/NTUPLE/Private2DUL18/SVJ_t-channel_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1_13TeV-madgraphMLM-pythia8/0_RA2AnalysisTree.root"],
            "-y 2018 -p analysis_configs.t_channel_pre_selection -pd dummy --weight_variation scale pdf",
        ),
        (
            ["root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2018A-UL2018-v1/JetHT/0_RA2AnalysisTree.root"],
            "-y 2018 -p analysis_configs.t_channel_pre_selection -pd JetHT",
        ),
        (
            ["root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2018A-UL2018-v1/EGamma/0_RA2AnalysisTree.root"],
            "-y 2018 -p analysis_configs.t_channel_wnae_qcd_training_region -pd EGamma",
        ),
        (
            ["root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Run2018A-UL2018-v3/SingleMuon/4_RA2AnalysisTree.root"],
            "-y 2018 -p analysis_configs.t_channel_wnae_top_training_region -pd SingleMuon",
        ),
        (
            ["root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Summer20UL18/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/0_RA2AnalysisTree.root"],
            "-y 2018 -p analysis_configs.t_channel_lost_lepton_control_region -pd dummy",
        ),
    ]

    if "cmslpc" in os.environ["HOSTNAME"]:
        run_particle_net = True
    elif "t3ui" in os.environ["HOSTNAME"]:
        run_particle_net = False
    else:
        run_particle_net = False

    output_path = helper.get_temporary_directory()
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    for params in params_list:
        files_list, flags = params
        input_files = ",".join(files_list)
        output_name = helper.run_bash_command(f'echo {input_files}_$(date +"%Y%m%d-%H%M%S") | shasum | cut -d " " -f1')
        output_file = f"{output_path}/{output_name}.root"
            
        __run_skimmer(input_files, output_file, flags, run_particle_net)

        assert os.path.exists(output_file)


test_execution()

