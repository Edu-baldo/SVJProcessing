#!/bin/bash

N_WORKERS=6

dataset_directory=/pnfs/psi.ch/cms/trivcat/store/t3groups/ethz-susy/PFNanoVHbb/UL2017/cmssw

selection_name=analysis_selection_Hbb_datasets_paths
#selection_name=t_channel_wnae_qcd_training_region
#selection_name=t_channel_wnae_top_training_region
#selection_name=t_channel_lost_lepton_control_region

year=2017

# Output directory for nominal samples - no variation of the uncertainties
output_directory=/work/ext-ebaldo/Hbb_analysis/output_selection_Hbb


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
    WJetsToLNu_Pt-250To400 
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


check_number_of_events() {

    local dataset_directory=$1
    local selection_name=$2
    local year=$3
    local dataset_name=$4
    local output_directory=$5

    # Path automatically built when preparing input files lists
    local files_list=${dataset_directory}/files_list/${year}/${dataset_name}.csv
    local output_directory=${output_directory}/${year}/${selection_name}/${dataset_name}

    python check_number_of_events.py -i ${files_list} -o ${output_directory} -n ${N_WORKERS}
}


for dataset_name in ${dataset_names[@]}; do

    check_number_of_events ${dataset_directory} ${selection_name} ${year} ${dataset_name} ${output_directory}

done

