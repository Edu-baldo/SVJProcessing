#!/bin/bash

MEMORY=10GB
TIME=12:00:00
PARTITION=standard
CORES=2
CHUNK_SIZE=1000
N_WORKERS=150
#EXECUTOR=dask/lpccondor    # HTCondor at LPC
EXECUTOR=dask/slurm         # slurm at PSI
FORCE_RECREATE=0   # 1 to recreate output file if it exists, 0 else
FIRST_FILE=0
LAST_FILE=-1  # Use -1 to skim all input files

dataset_directory=/work/cazzanig/datasets_mc_nohem/

module=analysis_configs.s_channel_leptons_mc_pre_selection_nohem
selection_name=s_channel_leptons_pre_selection


year=2017

output_directory=root://t3dcachedb03.psi.ch//pnfs/psi.ch/cms/trivcat/store/user/cazzanig/schannel_leptons_skims_nohem_mc_2018/SVJleptons_signal/


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


make_skims() {

    local dataset_directory=$1
    local module=$2
    local selection_name=$3
    local year=$4
    local dataset_name=$5
    local output_directory=$6
    local xsec=$7

    # Path automatically built when preparing input files lists
    local suffix_dir=${year}/${selection_name}/${dataset_name}
    local files_list_directory=${dataset_directory}/skim_input_files_list/${suffix_dir}
    local output_directory=${output_directory}/${suffix_dir}

    local output_redirector=$(echo ${output_directory} | cut -d/ -f 1-4)
    local output_dir=$(echo ${output_directory} | cut -d/ -f 4-)
    xrdfs ${output_redirector} ls ${output_dir} > /dev/null 2>&1
    if [ "$?" != "0" ]; then
        xrdfs ${output_redirector} mkdir -p ${output_dir}
    fi

    i_file=-1
    for files_list in $(ls ${files_list_directory} | sort -V); do
        ((i_file++))
        if [ ${i_file} -le ${LAST_FILE} ] || [ "${LAST_FILE}" == "-1" ]; then
            if [ ${i_file} -ge ${FIRST_FILE} ]; then

                local input_files=${files_list_directory}/${files_list}
                local output_file=${output_directory}/${files_list/.txt/.root}
                local output_file_name_tmp=$(echo ${ouput_file}_$(date +"%Y%m%d-%H%M%S") | shasum | cut -d " " -f1).root
                local output_file_tmp=/work/${USER}/tmp/${output_file_name_tmp}   

                echo ""
                echo "Making skim file ${output_file}"

                local output_redirector=$(echo ${output_file} | cut -d/ -f 1-4)
                local output_file_path=$(echo ${output_file} | cut -d/ -f 4-)
                xrdfs ${output_redirector} ls ${output_file_path} > /dev/null 2>&1
                if [ "$?" != "0" ] || [ "${FORCE_RECREATE}" == "1" ]; then
                    python skim.py -i ${input_files} -o ${output_file_tmp} -p ${module} -pd ${dataset_name} -y ${year} -e ${EXECUTOR} -n ${N_WORKERS} -c ${CHUNK_SIZE} --memory ${MEMORY} --cores ${CORES} --walltime ${TIME} --queue ${PARTITION} -xsec ${xsec}   -nano
                    xrdcp -f ${output_file_tmp} ${output_file}
                    echo ${output_file} has been saved.
                    rm ${output_file_tmp}
                else
                    echo ${output_file} already exists and FORCE_RECREATE is 0. Skipping.
                fi
            fi
        fi
    done
}


for dataset_name in ${dataset_names[@]}; do

    make_skims ${dataset_directory} ${module} ${selection_name} ${year} ${dataset_name} ${output_directory}

done
