#!/bin/bash

dataset_directory=/pnfs/psi.ch/cms/trivcat/store/t3groups/ethz-susy/PFNanoVHbb/UL2017/cmssw
dataset_config=dataset_configs.analysis_selection_Hbb_datasets_paths

module=analysis_configs.analysis_selection_Hbb_datasets_paths
selection_name=analysis_selection_Hbb_datasets_paths

#module=analysis_configs.t_channel_wnae_qcd_training_region
#selection_name=t_channel_wnae_qcd_training_region

#module=analysis_configs.t_channel_wnae_top_training_region
#selection_name=t_channel_wnae_top_training_region

#module=analysis_configs.t_channel_lost_lepton_control_region
#selection_name=t_channel_lost_lepton_control_region

year=2017

dataset_names=(
    #
    # Signals
    #
    ggZH_HToBB_ZToLL
    ggZH_HToBB_ZToNuNu
    #
    # Backgrounds
    #
    #
    # Wjets
    #
    WJetsToLNu_Pt-100To250
    wjets/WJetsToLNu_Pt-250To400 
    WJetsToLNu_Pt-400To600
    WJetsToLNu_Pt-600ToInf
    #
    # Zjets
    #
    Z1JetsToNuNu_M-50_LHEFilterPtZ-150To250
    Z1JetsToNuNu_M-50_LHEFilterPtZ-50To150   
    Z2JetsToNuNu_M-50_LHEFilterPtZ-400ToInf
    Z1JetsToNuNu_M-50_LHEFilterPtZ-250To400
    Z2JetsToNuNu_M-50_LHEFilterPtZ-150To250
    Z2JetsToNuNu_M-50_LHEFilterPtZ-50To150
    Z1JetsToNuNu_M-50_LHEFilterPtZ-400ToInf
    Z2JetsToNuNu_M-50_LHEFilterPtZ-250To400
)

prepare_input_files_list() {

    local dataset_config=$1
    local dataset_directory=$2
    local module=$3
    local selection_name=$4
    local year=$5
    local dataset_name=$6

    echo ""
    echo "Preparing input files for dataset ${dataset_name} year ${year} and selection ${selection_name}"

    # Add the right path for each dataset
    if [[ ${dataset_name} == ggZH_HToBB* ]]; then
        dataset_path="${dataset_directory}/hbb/${dataset_name}"
    elif [[ ${dataset_name} == WJetsToLNu* ]]; then
        dataset_path="${dataset_directory}/wjets/${dataset_name}"
    elif [[ ${dataset_name} == Z1JetsToNuNu* || ${dataset_name} == Z2JetsToNuNu* ]]; then
        dataset_path="${dataset_directory}/zjets/${dataset_name}"
    else
        echo "Unknown dataset ${dataset_name}"
        return 1
    fi

    python list_dataset_files.py -d ${dataset_name} -y ${year} -c ${dataset_config} -o ${dataset_directory} -nano_scout
    python compute_unweighted_selection_efficiency.py -d ${dataset_name} -y ${year} -p ${module} -s ${selection_name} -i ${dataset_directory} -o ${dataset_directory} -n 6 -e futures -c 10000 -nano_scout
    python prepare_input_files_list.py -d ${dataset_name} -y ${year} -s ${selection_name} -i ${dataset_directory} -o ${dataset_directory} -m 50000
}

for dataset_name in ${dataset_names[@]}; do

    prepare_input_files_list ${dataset_config} ${dataset_directory} ${module} ${selection_name} ${year} ${dataset_name}

done
