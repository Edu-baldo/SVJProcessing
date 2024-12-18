import os

from tests import helper


def __run_skimmer(input_files, output_file, config, year, primary_dataset, run_particle_net):
    n_workers = 8
    chunk_size = 100000
    executor = "futures"

    bash_command = f"python {os.environ['SVJ_PROCESSING_ROOT']}/skimmer/skim.py -i {input_files} -o {output_file} -p {config} -y {year} -e {executor} -n {n_workers} -c {chunk_size} -pd {primary_dataset} -nano"
    if run_particle_net:
        bash_command += " -pn_tagger"
    helper.test_command(bash_command)


def test_execution():
    params_list = [
        (
            "2018",
            "analysis_configs.s_channel_leptons_mc_pre_selection",
            [f"root://storage01.lcg.cscs.ch:1096//pnfs/lcg.cscs.ch/cms/trivcat/store/user/cazzanig/darkshowers/samples/cms_offline/pfnano/PFNano_UL18_100k/PFNanoAOD/PFNanoAOD_SVJL_mMed-2000GeV_mDark-8GeV_rinv-0.3_alpha-peak_13TeV-pythia8_part-{i}.root" for i in range(2,9)],
            "dummy",
        ),
        
    ]

    if "cmslpc" in os.environ["HOSTNAME"]:
        run_particle_net = True
    elif "t3ui" in os.environ["HOSTNAME"]:
        run_particle_net = False

    output_path = helper.get_temporary_directory()
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    for params in params_list:
        year, config, files_list, primary_dataset = params
        input_files = ",".join(files_list)
        output_name = helper.run_bash_command(f'echo {input_files}_$(date +"%Y%m%d-%H%M%S") | shasum | cut -d " " -f1')
        output_file = f"{output_path}/{output_name}.root"
            
        __run_skimmer(input_files, output_file, config, year, primary_dataset, run_particle_net)

        assert os.path.exists(output_file)


test_execution()

