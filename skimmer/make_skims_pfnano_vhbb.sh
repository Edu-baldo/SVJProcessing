#!/bin/bash

MEMORY=10GB
TIME=12:00:00
PARTITION=standard
CORES=2
CHUNK_SIZE=1000
N_WORKERS=150
#EXECUTOR=dask/lpccondor    # HTCondor at LPC
EXECUTOR=dask/slurm         # slurm at PSI
#EXECUTOR=futures        # run interactively
FORCE_RECREATE=0   # 1 to recreate output file if it exists, 0 else
FIRST_FILE=0
LAST_FILE=-1  # Use -1 to skim all input files

dataset_directory=/work/ext-ebaldo/datasets_hbb/

module=analysis_configs.0_leptons_selection_Hbb_boosted
selection_name=0_leptons_selection_Hbb_boosted


years=(
    #2016
    2017
    #2018
)

output_directory=/pnfs/psi.ch/cms/trivcat/store/user/ext-ebaldo/vh_bb_test_skims/
#root://t3dcachedb03.psi.ch/

dataset_names=(
    #
    # Signals
    #
    #ggZH_HToBB_ZToLL
    #ggZH_HToBB_ZToNuNu
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
    #Z1JetsToNuNu_M-50_LHEFilterPtZ-150To250
    Z1JetsToNuNu_M-50_LHEFilterPtZ-50To150   
    #Z2JetsToNuNu_M-50_LHEFilterPtZ-400ToInf
    #Z1JetsToNuNu_M-50_LHEFilterPtZ-250To400
    #Z2JetsToNuNu_M-50_LHEFilterPtZ-150To250
    #Z2JetsToNuNu_M-50_LHEFilterPtZ-50To150
    #Z1JetsToNuNu_M-50_LHEFilterPtZ-400ToInf
    #Z2JetsToNuNu_M-50_LHEFilterPtZ-250To400
)

variations=(
    nominal
#    # JEC/JER variations only for signal!
#    jec_up
#    jec_down
#    jer_up
#    jer_down
)


cross_sections=(
    #
    # Signals
    #
    #0.0062
    #0.0122
    #
    # Backgrounds
    #
    #
    # Wjets
    #
    #763.7
    #27.55
    #3.48
    #0.5415
    #
    # Zjets
    #
    #17.15
    583.4   
    #0.8318
    #1.972
    #29.28
    #308.2
    #0.216
    #4.96
)


make_skims() {

    local dataset_directory=$1
    local module=$2
    local selection_name=$3
    local year=$4
    local variation=$5
    local dataset_name=$6
    local output_directory=$7
    local xsec=$8

    # Path automatically built when preparing input files lists
    local files_list_directory=${dataset_directory}/skim_input_files_list/${year}/${selection_name}/${dataset_name}
    local output_directory=${output_directory}/${year}/${selection_name}/${variation}/${dataset_name}

    if [[ "${output_directory}"  == "root://"* ]]; then  # if the output has a redirector
        local output_redirector=$(echo ${output_directory} | cut -d/ -f 1-4)
        local output_dir=$(echo ${output_directory} | cut -d/ -f 4-)
        xrdfs ${output_redirector} ls ${output_dir} > /dev/null 2>&1
        if [ "$?" != "0" ]; then
            xrdfs ${output_redirector} mkdir -p ${output_dir}
        fi
    else
        local output_redirector="none"
        if [ ! -d ${output_directory} ]; then
            mkdir -p ${output_directory}
        fi
    fi

    i_file=-1
    for files_list in $(ls ${files_list_directory} | sort -V); do
        ((i_file++))
        if [ ${i_file} -le ${LAST_FILE} ] || [ "${LAST_FILE}" == "-1" ]; then
            if [ ${i_file} -ge ${FIRST_FILE} ]; then

                local input_files=${files_list_directory}/${files_list}
                local output_file=${output_directory}/${files_list/.txt/.root}
                local output_file_name_tmp=$(echo ${ouput_file}_$(date +"%Y%m%d-%H%M%S") | shasum | cut -d " " -f1).root
                local output_file_tmp=/scratch/${USER}/tmp/${output_file_name_tmp}

                echo ""
                echo "Making skim file ${output_file}"

                if [ "${output_redirector}" == "none" ]; then
                    ls ${output_file} > /dev/null 2>&1
                else
                    local output_file_path=$(echo ${output_file} | cut -d/ -f 4-)
                    xrdfs ${output_redirector} ls ${output_file_path} > /dev/null 2>&1
                fi
                if [ "$?" != "0" ] || [ "${FORCE_RECREATE}" == "1" ]; then
                    if [ "${variation}" == "nominal" ]; then
                        variation_flag=""
                    else
                        variation_flag="--variation ${variation}"
                    fi
                    if [[ ${dataset_name} == t-channel* ]]; then
                        weight_variation_flag="--weight_variation scale pdf"
                    else
                        weight_variation_flag=""
                    fi
                    python skim.py -i ${input_files} -o ${output_file_tmp} -p ${module} -pd ${dataset_name} -y ${year} -e ${EXECUTOR} -n ${N_WORKERS} -c ${CHUNK_SIZE} --memory ${MEMORY} --cores ${CORES} --walltime ${TIME} --queue ${PARTITION} -pn_tagger ${variation_flag} ${weight_variation_flag} -xsec ${xsec} -nano 
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


n_datasets=${#dataset_names[@]}


for ((i=0; i<$n_datasets; i++)); do
    dataset_name=${dataset_names[i]}
    cross_section=${cross_sections[i]}
  for year in ${years[@]}; do
    for variation in ${variations[@]}; do
      make_skims ${dataset_directory} ${module} ${selection_name} ${year} ${variation} ${dataset_name} ${output_directory} ${cross_section}
    done
  done
done