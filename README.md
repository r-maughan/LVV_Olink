# LVV_Olink

This Repo provides the code used to conduct the core analyses described in our recent publication: 
"Proteomic profiling of the large vessel vasculitis spectrum identifies shared signatures of innate immune activation and stromal remodelling" DOI https://doi.org/10.1101/2024.09.06.24313218


## Description

In this project, we analysed the levels of 184 inflammation and vascular related proteins in the plasma of 281 large vessel vasculitis patients and 124 control participants across two independent cohorts. The raw data is provided alongside R scripts which will carry out quality control (QC) analysis and filtering as well as differential abundance (DA) analysis. The QC and DA analyses are conducted on each cohort separately. No output data is provided, the code should be run in R within the directory structure of the Repo. This will then produce the various output tables and plots.


## Getting Started

### Dependencies

* Each script should check and automatically install the R packages required but the following packages are used: "OlinkAnalyze", "tidyverse", "FactoMineR", "factoextra", "gridExtra", "ggplot2","grid", "this.path", "tidyr",  "pheatmap",   "dplyr", "ggrepel"
* Written using R version 4.4.2, RStudio version 2024.09.0+375

### Running the scripts

* Download the full repo complete with all directories
* Run the QC and DA scripts within each Cohort's directory to conduct analysis and produce output files


## Authors

* Robert T. Maughan  [@R-maughan](https://github.com/r-maughan)  - _Primary_ 
* Erin Macdonald-Dunlop - _Initial work_
* Jimmy Peters - _Supervision_


## Version History

* 0.1 - Initial Release


## License

This project is licensed under the MIT License - see the LICENSE.md file for details


## Acknowledgments

Publication Authors: Robert T. Maughan, Erin MacDonald-Dunlop, Lubna Haroon-Rashid, Louise Sorensen, Natalie Chaddock, Shauna Masters, Andrew Porter, Marta Peverelli, Charis Pericleous, Andrew Hutchings, James Robinson, Taryn Youngstein, Raashid A. Luqmani, Justin C. Mason, Ann W. Morgan, James E. Peters

Authors and maintainers of packages cited above