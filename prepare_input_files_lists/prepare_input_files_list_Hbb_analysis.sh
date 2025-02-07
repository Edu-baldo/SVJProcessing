#!/bin/bash

dataset_directory=/work/ext-ebaldo/datasets_hbb/
dataset_config=dataset_configs.analysis_selection_Hbb_datasets_paths

module=analysis_configs.0_leptons_selection_Hbb_boosted
selection_name=0_leptons_selection_Hbb_boosted

# module=analysis_configs.1_leptons_selection_Hbb_boosted
# selection_name=1_leptons_selection_Hbb_boosted

# module=analysis_configs.2_leptons_selection_Hbb_boosted
# selection_name=2_leptons_selection_Hbb_boosted

year=2017

dataset_names=(
    # 
    # Signals
    #
    # ggZH_HToBB_ZToLL
    ggZH_HToBB_ZToNuNu
    #
    # Backgrounds
    #
    #
    # Wjets
    #
    #WJetsToLNu_Pt-100To250
    #WJetsToLNu_Pt-250To400 
    #WJetsToLNu_Pt-400To600
    #WJetsToLNu_Pt-600ToInf
    #
    # Zjets
    #
    # Z1JetsToNuNu_M-50_LHEFilterPtZ-150To250
    # Z1JetsToNuNu_M-50_LHEFilterPtZ-50To150   
    # Z2JetsToNuNu_M-50_LHEFilterPtZ-400ToInf
    # Z1JetsToNuNu_M-50_LHEFilterPtZ-250To400
    # Z2JetsToNuNu_M-50_LHEFilterPtZ-150To250
    # Z2JetsToNuNu_M-50_LHEFilterPtZ-50To150
    # Z1JetsToNuNu_M-50_LHEFilterPtZ-400ToInf
    # Z2JetsToNuNu_M-50_LHEFilterPtZ-250To400
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

    python list_dataset_files.py -d ${dataset_name} -y ${year} -c ${dataset_config} -o ${dataset_directory} -nano
    python compute_unweighted_selection_efficiency.py -d ${dataset_name} -y ${year} -p ${module} -s ${selection_name} -i ${dataset_directory} -o ${dataset_directory} -n 6 -e futures -c 10000 -nano
    python prepare_input_files_list.py -d ${dataset_name} -y ${year} -s ${selection_name} -i ${dataset_directory} -o ${dataset_directory} -m 50000
}

for dataset_name in ${dataset_names[@]}; do

    prepare_input_files_list ${dataset_config} ${dataset_directory} ${module} ${selection_name} ${year} ${dataset_name}

done 

