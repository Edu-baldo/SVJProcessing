###################################  README  ###################################
#
# Is called "dataset" a set of files corresponding to the same physics process.
# The object `datasets_info` describes the location of the different datasets.
# Its structure is the following:
#    * Keys are year
#    * Values are the "datasets_per_year_info"
# The structure of the `datasets_per_year_info` is the following:
#    * Keys are dataset names
#    * Values are the "dataset_info" defining which files belong to the dataset
# 
# The "dataset_info" has the following structure. It is a list of dict, which
# has 2 keys:
#    * "redirector": The XRootD redirector to the remote storage element
#    * "path": The path to the directory at which the files are located
#    * "regex": The regex to apply to select some files from that directory. 
#               The regex must be "" if no regex is applied.
#
################################################################################
#import pprint

year = "2017"

datasets_info = {year: {}}


signal_ggZH = [
    "ggZH_HToBB_ZToLL",
    "ggZH_HToBB_ZToNuNu",
]


wjets_bins = [
    "WJetsToLNu_Pt-100To250",
    "WJetsToLNu_Pt-250To400", 
    "WJetsToLNu_Pt-400To600",
    "WJetsToLNu_Pt-600ToInf",
]

zjets_bins = [
    "Z1JetsToNuNu_M-50_LHEFilterPtZ-150To250",  
    "Z1JetsToNuNu_M-50_LHEFilterPtZ-50To150",   
    "Z2JetsToNuNu_M-50_LHEFilterPtZ-400ToInf",
    "Z1JetsToNuNu_M-50_LHEFilterPtZ-250To400",
    "Z2JetsToNuNu_M-50_LHEFilterPtZ-150To250",
    "Z2JetsToNuNu_M-50_LHEFilterPtZ-50To150",
    "Z1JetsToNuNu_M-50_LHEFilterPtZ-400ToInf",
    "Z2JetsToNuNu_M-50_LHEFilterPtZ-250To400",
]


background_bins = wjets_bins + zjets_bins

#Signal ggZH
datasets_info[year].update({
    signal: [
        {
            #"redirector": "",      #"root://t3se01.psi.ch:1094//",
            "path": f"/pnfs/psi.ch/cms/trivcat/store/t3groups/ethz-susy/PFNanoVHbb/UL2017/cmssw/hbb/{signal}/",
            "regex": f"PFNANOAOD_{signal}_part",
        },
    ]
    for signal in signal_ggZH
})



for bin in background_bins:

    if "WJetsToLNu" in bin:
        datasets_info[year].update({
            bin: [
                {
                    #"redirector": "",     #"root://t3se01.psi.ch:1094//",
                    "path": f"/pnfs/psi.ch/cms/trivcat/store/t3groups/ethz-susy/PFNanoVHbb/UL2017/cmssw/wjets/{bin}/",
                    "regex": f"PFNANOAOD_{bin}_part",
                }
            ]
        })

    elif "JetsToNuNu" in bin:
        datasets_info[year].update({
            bin: [
                {
                    #"redirector": "",     #"root://t3se01.psi.ch:1094//",
                    "path": f"/pnfs/psi.ch/cms/trivcat/store/t3groups/ethz-susy/PFNanoVHbb/UL2017/cmssw/zjets/{bin}/",
                    "regex": f"PFNANOAOD_{bin}_part",
                }
            ]
        })
    else:
        print(f"Unknown background {bin}")
        exit(1)

      
      
