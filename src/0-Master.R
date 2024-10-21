#!/usr/bin/env Rscript

# ======================================================================================= #
#                                                                                         #
#  Master R script - version 2.1 (October, 2024).                                         #
#                                                                                         #
#                                                                                         #
#  Purpose:                                                                               #
#  This is a 'master' script, designed to run all analysis scripts.                       #
#                                                                                         #
#                                                                                         #
#  Version:                                                                               #
#  2.1 supersedes all previous versions. (v2.0 -> v2.1)                                   #
#                                                                                         #
#  It was altered to be GitHub compatible, and the GitHub template created by Daniel      #
#  Morillo. Daniel Morillo's template is licensed under Creative Commons Attribution 4.0  #
#  International license. See: https://github.com/DaniMori/rproj-template  and:           #
#  https://creativecommons.org/licenses/by/4.0/                                           #
#                                                                                         #
#                                                                                         #
#  Related scripts:                                                                       #
#  This script is designed to work in conjunction with the TNT scripts created for        #
#  Young et al. (2024: ZJLS 200:547-617). See: https://doi.org/10.1093/zoolinnean/zlad165 #
#                                                                                         #
#                                                                                         #
#  Origin:                                                                                #
#  Version 1.0 of this script was based on Moon & Stubbs (2020: Communications Biology    #
#  3: 68). See: https://www.nature.com/articles/s42003-020-0779-6                         #
#                                                                                         #
#                                                                                         #
#  License and Rights Statement:                                                          #
#  This work is licensed under Creative Commons Attribution 4.0 International license     #
#  (CC BY 4.0). https://creativecommons.org/licenses/by/4.0/                              #
#  Users worldwide are free to copy, redistribute, remix, transform, and build upon the   #
#  material in any medium or format, even commercially. The licensor cannot revoke these  #
#  freedoms as long as you follow the license terms. You must give appropriate credit,    #
#  provide a link to the license, and indicate if changes were made. You may do so in any #
#  reasonable manner, but not in any way that suggests the licensor endorses you or your  #
#  use.                                                                                   #
#                                                                                         #
#                                                                                         #
#  Mark T. Young                                                                          #
#  marktyoung1984@gmail.com                                                               #
# ======================================================================================= #

# ================================================================ #
# 0. Installation of required packages and creation of directories #
# ================================================================ #


# -------------------------------------- #
# 0.1 Required packages and installation #
# -------------------------------------- #

 library(beepr)

### Packages on CRAN:
 list_of_packages <- c("ape", "beeper", "Claddis", "data.table", "doParallel", "geoscale", 
                       "geiger", "ggplot2", "graphic", "grid", "foreach", "magrittr",
                       "mailR", "maps", "Morpho", "nlme", "paleotree", "parallel", 
                       "phytools", "plot3D", "rr2", "strap", "TeachingDemos", 
                       "TreeTools", "vegan")

 new_packages <- list_of_packages[!(list_of_packages %in%
                                   installed.packages()[, "Package"])]


### Change CRAN repository as desired:
 if (length(new_packages)) {
   install.packages(new_packages, repos = "https://www.stats.bris.ac.uk/R/", dependencies = TRUE)
 }

 beep(1)

### IF you have never installed any of these packages, or after an R update, using "dependencies = TRUE" can help if something *odd* occurs ###
### i.e.: install.packages("ape", dependencies = TRUE) ###
### or: install.packages(list_of_packages, dependencies = TRUE) ###


 if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

 BiocManager::install("ggtree")


# --------------------------------------------------- #
# 0.2 Check for required directories and installation #
# --------------------------------------------------- #

### The repository file structure follows the GitHub template created by Daniel Morillo (see above for link to license and template)

### |
### |--- apps         (To store apps, e.g. in Shiny)
### |
### |--- dat          (To store input datasets)
### |
### |--- doc          (To store important documentation of the project)
### |    |
### |    |--- minutes (To store meeting minutes)
### |
### |--- notebooks    (Notebooks to explore data and test processes live here)
### |
### |--- output       (Processing outputs)
### |
### |--- R            (R functions created for this project)
### |
### |--- renv         (System library necessary for `renv` to work)
### |
### |--- src          (Source scripts that implement the main processes)
### |
### |--- www          (Project assets, e.g., images, bibliography files, etc.)



### check for directories:
 directories <- c("./apps",
                  "./dat", "./dat/Ages", "./dat/Trees", "./dat/NEXUS", "./dat/SinusData",
                  "./doc", "./doc/minutes",
                  "./notebooks",
                  "./output",
                  "./output/ParatympanicAnalyses", "./output/ParatympanicAnalyses/NonPhylo","./output/ParatympanicAnalyses/NonPhylo/Files","./output/ParatympanicAnalyses/NonPhylo/Figures",
		  "./output/ParatympanicAnalyses/EW_C0","./output/ParatympanicAnalyses/EW_C0/Files", "./output/ParatympanicAnalyses/EW_C0/Figures",
	  	  "./output/ParatympanicAnalyses/EW_C1","./output/ParatympanicAnalyses/EW_C1/Files", "./output/ParatympanicAnalyses/EW_C1/Figures",
		  "./output/ParatympanicAnalyses/EW_C3","./output/ParatympanicAnalyses/EW_C3/Files", "./output/ParatympanicAnalyses/EW_C3/Figures",
		  "./output/ParatympanicAnalyses/EW_C4","./output/ParatympanicAnalyses/EW_C4/Files", "./output/ParatympanicAnalyses/EW_C4/Figures",
	  	  "./output/ParatympanicAnalyses/EIWk15_C0","./output/ParatympanicAnalyses/EIWk15_C0/Files", "./output/ParatympanicAnalyses/EIWk15_C0/Figures",
	  	  "./output/ParatympanicAnalyses/EIWk15_C1","./output/ParatympanicAnalyses/EIWk15_C1/Files", "./output/ParatympanicAnalyses/EIWk15_C1/Figures",
		  "./output/ParatympanicAnalyses/EIWk15_C3","./output/ParatympanicAnalyses/EIWk15_C3/Files", "./output/ParatympanicAnalyses/EIWk15_C3/Figures",
		  "./output/ParatympanicAnalyses/EIWk15_C4","./output/ParatympanicAnalyses/EIWk15_C4/Files", "./output/ParatympanicAnalyses/EIWk15_C4/Figures",
                  "./output/ParanasalAnalyses", "./output/ParanasalAnalyses/NonPhylo", "./output/ParanasalAnalyses/NonPhylo/Files","./output/ParanasalAnalyses/NonPhylo/Figures",
		  "./output/ParanasalAnalyses/EW_C0","./output/ParanasalAnalyses/EW_C0/Files", "./output/ParanasalAnalyses/EW_C0/Figures",
	  	  "./output/ParanasalAnalyses/EW_C1","./output/ParanasalAnalyses/EW_C1/Files", "./output/ParanasalAnalyses/EW_C1/Figures",
		  "./output/ParanasalAnalyses/EW_C3","./output/ParanasalAnalyses/EW_C3/Files", "./output/ParanasalAnalyses/EW_C3/Figures",
		  "./output/ParanasalAnalyses/EW_C4","./output/ParanasalAnalyses/EW_C4/Files", "./output/ParanasalAnalyses/EW_C4/Figures",
	  	  "./output/ParanasalAnalyses/EIWk15_C0","./output/ParanasalAnalyses/EIWk15_C0/Files", "./output/ParanasalAnalyses/EIWk15_C0/Figures",
	  	  "./output/ParanasalAnalyses/EIWk15_C1","./output/ParanasalAnalyses/EIWk15_C1/Files", "./output/ParanasalAnalyses/EIWk15_C1/Figures",
		  "./output/ParanasalAnalyses/EIWk15_C3","./output/ParanasalAnalyses/EIWk15_C3/Files", "./output/ParanasalAnalyses/EIWk15_C3/Figures",
		  "./output/ParanasalAnalyses/EIWk15_C4","./output/ParanasalAnalyses/EIWk15_C4/Files", "./output/ParanasalAnalyses/EIWk15_C4/Figures",
                  "./R",
                  "./renv",
                  "./src",
                  "./www")


 missing_dir <- directories[!(directories %in% list.dirs())]

### Create new directories:
 for (dir in missing_dir) dir.create(dir)

 beep(1)


# ============================================================== #
# 1. PCoA, CVA, phylogenetic comparative methods. Sinus Analyses #
# ============================================================== #

# Uses cladistic dataset to run a PCoA and then CVAs on taxon-distance matrices. #
# Runs a series of multivariate statistics to test habitat differences.          #
# Creates contmaps and phenograms for each PCo axis.                             #
# The two scripts are identical, running for the two different sinus datasets.   #
# Within each script, the 8 different phylogenetic analyses are tested as        #
# sensitivity analyses (equal weights vs extended implied weights k =15) and     #
# the four different thalattosuchian positional hypotheses (C0, C1, C3, C4).     #

 source("src/1-ParatympanicSinusAnalysis.R")

 source("src/2-ParanasalSinusAnalysis.R")


# ========= #
# 0.3 Done! #
# ========= #

### Noise when script has finished (useful if you are running this in the background):
 library(beepr)

 beep(8)