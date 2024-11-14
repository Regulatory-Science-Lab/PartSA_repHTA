This repository provides code to estimate the cost-effectiveness of tumour-agnostic therapy entrectinib, using partitioned survival analysis. Alongside this repository is a dashboard that enables adaptation of selected input parameters. No coding is required to use the dashboard, which is available [here](https://regulatorysciencelab.shinyapps.io/partsa_rephta/).

The script to run all stages of the analysis ([main.R](analysis/main.R)) is found in the [analysis](https://github.com/Regulatory-Science-Lab/PartSA_repHTA/tree/main/analysis) folder, along with scripts to execute the deterministic and probabilistic analyses. All function files are contained in the [R](https://github.com/Regulatory-Science-Lab/PartSA_repHTA/tree/main/R) folder. Input parameters, taken from publicly available Health Technology Assessment (HTA) reimbursement reviews,  and all files relating to setting up the parameter space and survival curves are found in [data-raw](https://github.com/Regulatory-Science-Lab/PartSA_repHTA/tree/main/data-raw). 

Before running this analysis, the packages 'pacman', 'devtools' and 'dampack' are required. Installation of these packages is found in the [Install_packages.R](https://github.com/Regulatory-Science-Lab/PartSA_repHTA/blob/main/analysis/Install_packages.R) file in the analysis folder. All other packages will be installed (if necessary) and loaded via pload in main.R. 