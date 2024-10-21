#!/usr/bin/env Rscript

# ======================================================================================= #
#                                                                                         #
#  Paranasal sinus analysis script - version 2.1 (October, 2024).                         #
#                                                                                         #
#                                                                                         #
#  Purpose:                                                                               #
#  Runs all of the multivariate and phylogenetic comparative analyses for the snout       #
#  (paranasal) sinus dataset.                                                             #
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
#  It was also altered to be compatible with Claddis v 0.7.0                              #
#                                                                                         #
#                                                                                         #
#  Related scripts:                                                                       #
#  This script is designed to work in conjunction with the TNT scripts created for        #
#  Young et al. (2024: ZJLS 200:547-617). See: https://doi.org/10.1093/zoolinnean/zlad165 #
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

# ========================== #
# 0.1 Parallelisation set-up #
# ========================== #

### Required libraries:
 library(parallel)
 library(foreach)
 library(doParallel)

### Set-up (note, fork is used as I use Mac OS.)
 core_number <- detectCores() - 1
 cluster <- makeForkCluster(core_number)
 registerDoParallel(cluster)


# ========================= #
# 0.2 Colour palette set-up #
# ========================= #

 mycolours<-c("steelblue","coral","brown","violetred","darkolivegreen2","darkgoldenrod1","chocolate1","black","plum4","red")
 palette(mycolours)



# ============================== #
# 1. Input and format tree files #
# ============================== #

### Required libraries:
 library(ape)
 library(beepr)
 library(geoscale)
 library(TreeTools)
 library(magrittr)

###### NOTE: load in fully bifurcating i.e. binary phylogeny - needed for Claddis phylomorphospace and for Geiger analyses ######


### Loads in the random most parsimonious cladogram files from the TNT phylogenetic analyses, ladderises, and then converts them into nexus format using Ape. Plots them using Phytools. The %>% requires magrittr:
 EW_C0_phylogeny <- ReadTntTree("dat/Trees/xiw_equal_C0_randmpc.tre") %>% ladderize
 ape::write.tree(EW_C0_phylogeny, file='dat/Trees/xiw_equal_C0_randmpc.txt') 

 EW_C1_phylogeny <- ReadTntTree("dat/Trees/xiw_equal_C1_randmpc.tre") %>% ladderize
 ape::write.tree(EW_C1_phylogeny, file='dat/Trees/xiw_equal_C1_randmpc.txt') 

 EW_C3_phylogeny <- ReadTntTree("dat/Trees/xiw_equal_C3_randmpc.tre") %>% ladderize
 ape::write.tree(EW_C3_phylogeny, file='dat/Trees/xiw_equal_C3_randmpc.txt') 

 EW_C4_phylogeny <- ReadTntTree("dat/Trees/xiw_equal_C4_randmpc.tre") %>% ladderize
 ape::write.tree(EW_C4_phylogeny, file='dat/Trees/xiw_equal_C4_randmpc.txt') 

 EIWk15_C0_phylogeny <- ReadTntTree("dat/Trees/xiw_k15_C0_randmpc.tre") %>% ladderize
 ape::write.tree(EIWk15_C0_phylogeny, file='dat/Trees/xiw_k15_C0_randmpc.txt') 

 EIWk15_C1_phylogeny <- ReadTntTree("dat/Trees/xiw_k15_C1_randmpc.tre") %>% ladderize
 ape::write.tree(EIWk15_C1_phylogeny, file='dat/Trees/xiw_k15_C1_randmpc.txt') 

 EIWk15_C3_phylogeny <- ReadTntTree("dat/Trees/xiw_k15_C3_randmpc.tre") %>% ladderize
 ape::write.tree(EIWk15_C3_phylogeny, file='dat/Trees/xiw_k15_C3_randmpc.txt')

 EIWk15_C4_phylogeny <- ReadTntTree("dat/Trees/xiw_k15_C4_randmpc.tre") %>% ladderize
 ape::write.tree(EIWk15_C4_phylogeny, file='dat/Trees/xiw_k15_C4_randmpc.txt') 


### Prunes the paranasal-only dataset:
 paranasal.species<-c("Junggarsuchus_sloani","Protosuchus_haughtoni","cf_Hamadasuchus_rebouli","Campinasuchus_dinizi","Lohuecosuchus_megadontos","Gavialis_gangeticus","Tomistoma_schlegelii","Osteolaemus_tetraspis","Mecistops_cataphractus","Crocodylus_acutus","Crocodylus_johnstoni","Crocodylus_moreletii","Crocodylus_niloticus","Crocodylus_porosus","Crocodylus_rhombifer","Alligator_mississippiensis","Alligator_sinensis","Caiman_crocodilus","Paleosuchus_palpebrosus","Plagiophthalmosuchus_gracilirostris","Macrospondylus_bollensis","Pelagosaurus_typus","Eoneustes_gaudryi","Cricosaurus_schroederi","Cricosaurus_araucanensis","Thalattosuchus_superciliosus")


### Pruned phylogenies:
 EW_C0_phylogeny_SinPN_Pruned<-drop.tip(EW_C0_phylogeny,setdiff(EW_C0_phylogeny$tip.label,paranasal.species))
 plot(EW_C0_phylogeny_SinPN_Pruned, cex=0.5)

 EW_C1_phylogeny_SinPN_Pruned<-drop.tip(EW_C1_phylogeny,setdiff(EW_C1_phylogeny$tip.label,paranasal.species))
 plot(EW_C1_phylogeny_SinPN_Pruned, cex=0.5)

 EW_C3_phylogeny_SinPN_Pruned<-drop.tip(EW_C3_phylogeny,setdiff(EW_C3_phylogeny$tip.label,paranasal.species))
 plot(EW_C3_phylogeny_SinPN_Pruned, cex=0.5)

 EW_C4_phylogeny_SinPN_Pruned<-drop.tip(EW_C4_phylogeny,setdiff(EW_C4_phylogeny$tip.label,paranasal.species))
 plot(EW_C4_phylogeny_SinPN_Pruned, cex=0.5)

 EIWk15_C0_phylogeny_SinPN_Pruned<-drop.tip(EIWk15_C0_phylogeny,setdiff(EIWk15_C0_phylogeny$tip.label,paranasal.species))
 plot(EIWk15_C0_phylogeny_SinPN_Pruned, cex=0.5)

 EIWk15_C1_phylogeny_SinPN_Pruned<-drop.tip(EIWk15_C1_phylogeny,setdiff(EIWk15_C1_phylogeny$tip.label,paranasal.species))
 plot(EIWk15_C1_phylogeny_SinPN_Pruned, cex=0.5)

 EIWk15_C3_phylogeny_SinPN_Pruned<-drop.tip(EIWk15_C3_phylogeny,setdiff(EIWk15_C3_phylogeny$tip.label,paranasal.species))
 plot(EIWk15_C3_phylogeny_SinPN_Pruned, cex=0.5)

 EIWk15_C4_phylogeny_SinPN_Pruned<-drop.tip(EIWk15_C4_phylogeny,setdiff(EIWk15_C4_phylogeny$tip.label,paranasal.species))
 plot(EIWk15_C4_phylogeny_SinPN_Pruned, cex=0.5)


### checks if rooted and if binary:
 is.rooted(EW_C0_phylogeny_SinPN_Pruned)
 is.binary(EW_C0_phylogeny_SinPN_Pruned)
 is.rooted(EW_C1_phylogeny_SinPN_Pruned)
 is.binary(EW_C1_phylogeny_SinPN_Pruned)
 is.rooted(EW_C3_phylogeny_SinPN_Pruned)
 is.binary(EW_C3_phylogeny_SinPN_Pruned)
 is.rooted(EW_C4_phylogeny_SinPN_Pruned)
 is.binary(EW_C4_phylogeny_SinPN_Pruned)
 is.rooted(EIWk15_C0_phylogeny_SinPN_Pruned)
 is.binary(EIWk15_C0_phylogeny_SinPN_Pruned)
 is.rooted(EIWk15_C1_phylogeny_SinPN_Pruned)
 is.binary(EIWk15_C1_phylogeny_SinPN_Pruned)
 is.rooted(EIWk15_C3_phylogeny_SinPN_Pruned)
 is.binary(EIWk15_C3_phylogeny_SinPN_Pruned)
 is.rooted(EIWk15_C4_phylogeny_SinPN_Pruned)
 is.binary(EIWk15_C4_phylogeny_SinPN_Pruned)

 beep(1)


# ================================================================================= #
# 2. Creates time-dated phylogenies using Strap, 'equal method', and generates PDFs #
# ================================================================================= #

### Required libraries:
 library(beepr)
 library(strap)
 library(geiger)


### Loads in Age file:
 PNAges <- read.table("dat/Ages/SinusAgesPN.txt", header=T) 

### Use geiger to name check that the pruned phylogeny and age files match:
 name.check(EW_C0_phylogeny_SinPN_Pruned, PNAges)
 name.check(EW_C1_phylogeny_SinPN_Pruned, PNAges)
 name.check(EW_C3_phylogeny_SinPN_Pruned, PNAges)
 name.check(EW_C4_phylogeny_SinPN_Pruned, PNAges)
 name.check(EIWk15_C0_phylogeny_SinPN_Pruned, PNAges)
 name.check(EIWk15_C1_phylogeny_SinPN_Pruned, PNAges)
 name.check(EIWk15_C3_phylogeny_SinPN_Pruned, PNAges)
 name.check(EIWk15_C4_phylogeny_SinPN_Pruned, PNAges)


### C0 Equal weights phylogeny time-scaled with the equal method:
 EW_C0_phylogeny_SinPN_Pruned_EqualMethod <- DatePhylo(EW_C0_phylogeny_SinPN_Pruned, PNAges, method="equal", rlen=1)
 EW_C0_phylogeny_SinPN_Pruned_EqualMethod <- ladderize(EW_C0_phylogeny_SinPN_Pruned_EqualMethod)
 plot(EW_C0_phylogeny_SinPN_Pruned_EqualMethod)

 pdf("output/ParanasalAnalyses/EW_C0/Figures/EW_C0_EqualMethod_paranasalSinus_right.pdf", width=8, height=6) 
   plot(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, cex=0.7)
 dev.off()

 geoscalePhylo(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, PNAges, cex.age=0.6, cex.ts=0.6, cex.tip=0.5, quat.rm=TRUE, units = c("Period"), boxes="Period", direction="rightwards") 

 pdf("output/ParanasalAnalyses/EW_C0/Figures/EW_C0_EqualMethod_paranasalSinus_geoscale_right.pdf", width=17, height=9) 
   geoscalePhylo(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, PNAges, cex.age=0.6, cex.ts=0.6, cex.tip=0.5, quat.rm=TRUE, units = c("Period"), boxes="Period", direction="rightwards") 
 dev.off()


### C1 Equal weights phylogeny time-scaled with the equal method:
 EW_C1_phylogeny_SinPN_Pruned_EqualMethod <- DatePhylo(EW_C1_phylogeny_SinPN_Pruned, PNAges, method="equal", rlen=1)
 EW_C1_phylogeny_SinPN_Pruned_EqualMethod <- ladderize(EW_C1_phylogeny_SinPN_Pruned_EqualMethod)
 plot(EW_C1_phylogeny_SinPN_Pruned_EqualMethod)

 pdf("output/ParanasalAnalyses/EW_C1/Figures/EW_C1_EqualMethod_paranasalSinus_right.pdf", width=8, height=6) 
   plot(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, cex=0.7)
 dev.off()

 geoscalePhylo(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, PNAges, cex.age=0.6, cex.ts=0.6, cex.tip=0.5, quat.rm=TRUE, units = c("Period"), boxes="Period", direction="rightwards") 

 pdf("output/ParanasalAnalyses/EW_C1/Figures/EW_C1_EqualMethod_paranasalSinus_geoscale_right.pdf", width=17, height=9) 
   geoscalePhylo(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, PNAges, cex.age=0.6, cex.ts=0.6, cex.tip=0.5, quat.rm=TRUE, units = c("Period"), boxes="Period", direction="rightwards") 
 dev.off()


### C3 Equal weights phylogeny time-scaled with the equal method:
 EW_C3_phylogeny_SinPN_Pruned_EqualMethod <- DatePhylo(EW_C3_phylogeny_SinPN_Pruned, PNAges, method="equal", rlen=1)
 EW_C3_phylogeny_SinPN_Pruned_EqualMethod <- ladderize(EW_C3_phylogeny_SinPN_Pruned_EqualMethod)
 plot(EW_C3_phylogeny_SinPN_Pruned_EqualMethod)

 pdf("output/ParanasalAnalyses/EW_C3/Figures/EW_C3_EqualMethod_paranasalSinus_right.pdf", width=8, height=6) 
   plot(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, cex=0.7)
 dev.off()

 geoscalePhylo(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, PNAges, cex.age=0.6, cex.ts=0.6, cex.tip=0.5, quat.rm=TRUE, units = c("Period"), boxes="Period", direction="rightwards") 

 pdf("output/ParanasalAnalyses/EW_C3/Figures/EW_C3_EqualMethod_paranasalSinus_geoscale_right.pdf", width=17, height=9) 
   geoscalePhylo(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, PNAges, cex.age=0.6, cex.ts=0.6, cex.tip=0.5, quat.rm=TRUE, units = c("Period"), boxes="Period", direction="rightwards") 
 dev.off()


### C4 Equal weights phylogeny time-scaled with the equal method:
 EW_C4_phylogeny_SinPN_Pruned_EqualMethod <- DatePhylo(EW_C4_phylogeny_SinPN_Pruned, PNAges, method="equal", rlen=1)
 EW_C4_phylogeny_SinPN_Pruned_EqualMethod <- ladderize(EW_C4_phylogeny_SinPN_Pruned_EqualMethod)
 plot(EW_C4_phylogeny_SinPN_Pruned_EqualMethod)

 pdf("output/ParanasalAnalyses/EW_C4/Figures/EW_C4_EqualMethod_paranasalSinus_right.pdf", width=8, height=6) 
   plot(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, cex=0.7)
 dev.off()

 geoscalePhylo(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, PNAges, cex.age=0.6, cex.ts=0.6, cex.tip=0.5, quat.rm=TRUE, units = c("Period"), boxes="Period", direction="rightwards") 

 pdf("output/ParanasalAnalyses/EW_C4/Figures/EW_C4_EqualMethod_paranasalSinus_geoscale_right.pdf", width=17, height=9) 
   geoscalePhylo(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, PNAges, cex.age=0.6, cex.ts=0.6, cex.tip=0.5, quat.rm=TRUE, units = c("Period"), boxes="Period", direction="rightwards") 
 dev.off()


### C0 Extended implied weights (k=15) phylogeny time-scaled with the equal method:
 EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod <- DatePhylo(EIWk15_C0_phylogeny_SinPN_Pruned, PNAges, method="equal", rlen=1)
 EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod <- ladderize(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod)
 plot(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod)

 pdf("output/ParanasalAnalyses/EIWk15_C0/Figures/EIWk15_C0_EqualMethod_paranasalSinus_right.pdf", width=8, height=6) 
   plot(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, cex=0.7)
 dev.off()

 geoscalePhylo(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, PNAges, cex.age=0.6, cex.ts=0.6, cex.tip=0.5, quat.rm=TRUE, units = c("Period"), boxes="Period", direction="rightwards") 

 pdf("output/ParanasalAnalyses/EIWk15_C0/Figures/EIWk15_C0_EqualMethod_paranasalSinus_geoscale_right.pdf", width=17, height=9) 
   geoscalePhylo(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, PNAges, cex.age=0.6, cex.ts=0.6, cex.tip=0.5, quat.rm=TRUE, units = c("Period"), boxes="Period", direction="rightwards") 
 dev.off()


### C1 Extended implied weights (k=15) phylogeny time-scaled with the equal method:
 EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod <- DatePhylo(EIWk15_C1_phylogeny_SinPN_Pruned, PNAges, method="equal", rlen=1)
 EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod <- ladderize(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod)
 plot(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod)

 pdf("output/ParanasalAnalyses/EIWk15_C1/Figures/EIWk15_C1_EqualMethod_paranasalSinus_right.pdf", width=8, height=6) 
   plot(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, cex=0.7)
 dev.off()

 geoscalePhylo(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, PNAges, cex.age=0.6, cex.ts=0.6, cex.tip=0.5, quat.rm=TRUE, units = c("Period"), boxes="Period", direction="rightwards") 

 pdf("output/ParanasalAnalyses/EIWk15_C1/Figures/EIWk15_C1_EqualMethod_paranasalSinus_geoscale_right.pdf", width=17, height=9) 
   geoscalePhylo(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, PNAges, cex.age=0.6, cex.ts=0.6, cex.tip=0.5, quat.rm=TRUE, units = c("Period"), boxes="Period", direction="rightwards") 
 dev.off()


### C3 Extended implied weights (k=15) phylogeny time-scaled with the equal method:
 EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod <- DatePhylo(EIWk15_C3_phylogeny_SinPN_Pruned, PNAges, method="equal", rlen=1)
 EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod <- ladderize(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod)
 plot(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod)

 pdf("output/ParanasalAnalyses/EIWk15_C3/Figures/EIWk15_C3_EqualMethod_paranasalSinus_right.pdf", width=8, height=6) 
   plot(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, cex=0.7)
 dev.off()

 geoscalePhylo(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, PNAges, cex.age=0.6, cex.ts=0.6, cex.tip=0.5, quat.rm=TRUE, units = c("Period"), boxes="Period", direction="rightwards") 

 pdf("output/ParanasalAnalyses/EIWk15_C3/Figures/EIWk15_C3_EqualMethod_paranasalSinus_geoscale_right.pdf", width=17, height=9) 
   geoscalePhylo(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, PNAges, cex.age=0.6, cex.ts=0.6, cex.tip=0.5, quat.rm=TRUE, units = c("Period"), boxes="Period", direction="rightwards") 
 dev.off()


### C4 Extended implied weights (k=15) phylogeny time-scaled with the equal method:
 EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod <- DatePhylo(EIWk15_C4_phylogeny_SinPN_Pruned, PNAges, method="equal", rlen=1)
 EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod <- ladderize(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod)
 plot(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod)

 pdf("output/ParanasalAnalyses/EIWk15_C4/Figures/EIWk15_C4_EqualMethod_paranasalSinus_right.pdf", width=8, height=6) 
   plot(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, cex=0.7)
 dev.off()

 geoscalePhylo(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, PNAges, cex.age=0.6, cex.ts=0.6, cex.tip=0.5, quat.rm=TRUE, units = c("Period"), boxes="Period", direction="rightwards") 

 pdf("output/ParanasalAnalyses/EIWk15_C4/Figures/EIWk15_C4_EqualMethod_paranasalSinus_geoscale_right.pdf", width=17, height=9) 
   geoscalePhylo(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, PNAges, cex.age=0.6, cex.ts=0.6, cex.tip=0.5, quat.rm=TRUE, units = c("Period"), boxes="Period", direction="rightwards") 
 dev.off()

 beep(1)


# ======================================== #
# 3. The data for CVA and further analyses #
# ======================================== #

### Required libraries:
 library(beepr)
 library(data.table)
 library(geiger)
 library(Claddis)

### Paranasal-only csv data:
 CrocSinPNData<- read.csv("dat/SinusData/CrocSinPNData.csv",header = TRUE, row.names = 1)
 Habitat_SinSchwab<-CrocSinPNData$HabitatS2020
 Habitat_SinWilberg<-CrocSinPNData$HabitatW2019
 Habitat_SinSimple<-CrocSinPNData$HabitatSimple
 Terrestrial_Sin<-CrocSinPNData$Terrestrial
 Semiaquatic_Sin<-CrocSinPNData$Semiaquatic
 Pelagic_Sin<-CrocSinPNData$Pelagic
 Marine_Sin<-CrocSinPNData$Marine
 Freshwater_Sin<-CrocSinPNData$Freshwater
 Aquatic_Sin<-CrocSinPNData$Aquatic
 Clade_Sin<-CrocSinPNData$Clade
 Thalattosuchian_Sin<-CrocSinPNData$Thalattosuchian
 Metriorhynchid_Sin<-CrocSinPNData$Metriorhynchid
 Species_Sin<-CrocSinPNData$Species


### Use geiger to name check that the pruned paranasal-only phylogeny and data files match:
 name.check(EW_C0_phylogeny_SinPN_Pruned, CrocSinPNData)
 name.check(EW_C1_phylogeny_SinPN_Pruned, CrocSinPNData)
 name.check(EW_C3_phylogeny_SinPN_Pruned, CrocSinPNData)
 name.check(EW_C4_phylogeny_SinPN_Pruned, CrocSinPNData)
 name.check(EIWk15_C0_phylogeny_SinPN_Pruned, CrocSinPNData)
 name.check(EIWk15_C1_phylogeny_SinPN_Pruned, CrocSinPNData)
 name.check(EIWk15_C3_phylogeny_SinPN_Pruned, CrocSinPNData)
 name.check(EIWk15_C4_phylogeny_SinPN_Pruned, CrocSinPNData)


## Matrix files:
 CrocParanasalNexus<-read_nexus_matrix("dat/NEXUS/CrocParanasalMatrix.nex", equalize_weights = FALSE)

 beep(1)


# =========================================== #
# 4. Discrete Evolutionary Character analyses #
# =========================================== #

### Required libraries:
 library(beepr)
clusterCall(cluster, function () {
 library(Claddis)
})

 library(TeachingDemos)

### Time bins at the stage level:
 stage_fad <- as.vector(c(247.2,237.0,201.3,174.1,163.5,145.0,100.5,66.0,23.0,2.58))
 stage_lad <- as.vector(c(237.0,201.3,174.1,163.5,145.0,100.5,66.0,23.0,2.58,0.0))
 stage_matrix <- cbind(stage_fad, stage_lad)
 stage_bins <- matrix(stage_matrix, ncol = 2,  dimnames = list(LETTERS[1:10], c("fad", "lad")))
 class(stage_bins) <- "timeBins"


# ------------------------------------------------- #
# 4.1 Equal weights C0 paranasal character analysis #
# ------------------------------------------------- #

### Find edge labels:
 plot (EW_C0_phylogeny_SinPN_Pruned_EqualMethod)
 edgelabels()

 pdf("output/ParanasalAnalyses/EW_C0/Figures/EW_C0_EqualMethod_paranasalSinus_EdgeLabels.pdf", width=17, height=9) 
   plot (EW_C0_phylogeny_SinPN_Pruned_EqualMethod)
   edgelabels()
 dev.off()

### Every edge in cladogram rates analysis:
 txtStart("output/ParanasalAnalyses/EW_C0/Files/EW_C0_EqualMethod_paranasalSinus_BranchRates.txt")

 EW_C0_paranasalBranchRates <- test_rates(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, CrocParanasalNexus, time_bins = stage_bins, branch_partitions = lapply(X = as.list(x = 1:nrow(EW_C0_phylogeny_SinPN_Pruned_EqualMethod$edge)), as.list), alpha = 0.01)
 EW_C0_paranasalBranchRates$branch_test_results

 txtStop()


# ------------------------------------------------- #
# 4.2 Equal weights C1 paranasal character analysis #
# ------------------------------------------------- #

### Find edge labels:
 plot (EW_C1_phylogeny_SinPN_Pruned_EqualMethod)
 edgelabels()

 pdf("output/ParanasalAnalyses/EW_C1/Figures/EW_C1_EqualMethod_paranasalSinus_EdgeLabels.pdf", width=17, height=9) 
   plot (EW_C1_phylogeny_SinPN_Pruned_EqualMethod)
   edgelabels()
 dev.off()

### Every edge in cladogram rates analysis:
 txtStart("output/ParanasalAnalyses/EW_C1/Files/EW_C1_EqualMethod_paranasalSinus_BranchRates.txt")

 EW_C1_paranasalBranchRates <- test_rates(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, CrocParanasalNexus, time_bins = stage_bins, branch_partitions = lapply(X = as.list(x = 1:nrow(EW_C1_phylogeny_SinPN_Pruned_EqualMethod$edge)), as.list), alpha = 0.01)
 EW_C1_paranasalBranchRates$branch_test_results

 txtStop()


# ------------------------------------------------- #
# 4.3 Equal weights C3 paranasal character analysis #
# ------------------------------------------------- #

### Find edge labels:
 plot (EW_C3_phylogeny_SinPN_Pruned_EqualMethod)
 edgelabels()

 pdf("output/ParanasalAnalyses/EW_C3/Figures/EW_C3_EqualMethod_paranasalSinus_EdgeLabels.pdf", width=17, height=9) 
   plot (EW_C3_phylogeny_SinPN_Pruned_EqualMethod)
   edgelabels()
 dev.off()

### Every edge in cladogram rates analysis:
 txtStart("output/ParanasalAnalyses/EW_C3/Files/EW_C3_EqualMethod_paranasalSinus_BranchRates.txt")

 EW_C3_paranasalBranchRates <- test_rates(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, CrocParanasalNexus, time_bins = stage_bins, branch_partitions = lapply(X = as.list(x = 1:nrow(EW_C3_phylogeny_SinPN_Pruned_EqualMethod$edge)), as.list), alpha = 0.01)
 EW_C3_paranasalBranchRates$branch_test_results

 txtStop()


# ------------------------------------------------- #
# 4.4 Equal weights C4 paranasal character analysis #
# ------------------------------------------------- #

### Find edge labels:
 plot (EW_C4_phylogeny_SinPN_Pruned_EqualMethod)
 edgelabels()

 pdf("output/ParanasalAnalyses/EW_C4/Figures/EW_C4_EqualMethod_paranasalSinus_EdgeLabels.pdf", width=17, height=9) 
   plot (EW_C4_phylogeny_SinPN_Pruned_EqualMethod)
   edgelabels()
 dev.off()

### Every edge in cladogram rates analysis:
 txtStart("output/ParanasalAnalyses/EW_C4/Files/EW_C4_EqualMethod_paranasalSinus_BranchRates.txt")

 EW_C4_paranasalBranchRates <- test_rates(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, CrocParanasalNexus, time_bins = stage_bins, branch_partitions = lapply(X = as.list(x = 1:nrow(EW_C4_phylogeny_SinPN_Pruned_EqualMethod$edge)), as.list), alpha = 0.01)
 EW_C4_paranasalBranchRates$branch_test_results

 txtStop()


# ------------------------------------------------------------------- #
# 4.5 Extended implied weights (k=15) C0 paranasal character analysis #
# ------------------------------------------------------------------- #

### Find edge labels:
 plot (EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod)
 edgelabels()

 pdf("output/ParanasalAnalyses/EIWk15_C0/Figures/EIWk15_C0_EqualMethod_paranasalSinus_EdgeLabels.pdf", width=17, height=9) 
   plot (EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod)
   edgelabels()
 dev.off()

### Every edge in cladogram rates analysis:
 txtStart("output/ParanasalAnalyses/EIWk15_C0/Files/EIWk15_C0_EqualMethod_paranasalSinus_BranchRates.txt")

 EIWk15_C0_paranasalBranchRates <- test_rates(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, CrocParanasalNexus, time_bins = stage_bins, branch_partitions = lapply(X = as.list(x = 1:nrow(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod$edge)), as.list), alpha = 0.01)
 EIWk15_C0_paranasalBranchRates$branch_test_results

 txtStop()


# ------------------------------------------------------------------- #
# 4.6 Extended implied weights (k=15) C1 paranasal character analysis #
# ------------------------------------------------------------------- #

### Find edge labels:
 plot (EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod)
 edgelabels()

 pdf("output/ParanasalAnalyses/EIWk15_C1/Figures/EIWk15_C1_EqualMethod_paranasalSinus_EdgeLabels.pdf", width=17, height=9) 
   plot (EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod)
   edgelabels()
 dev.off()

### Every edge in cladogram rates analysis:
 txtStart("output/ParanasalAnalyses/EIWk15_C1/Files/EIWk15_C1_EqualMethod_paranasalSinus_BranchRates.txt")

 EIWk15_C1_paranasalBranchRates <- test_rates(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, CrocParanasalNexus, time_bins = stage_bins, branch_partitions = lapply(X = as.list(x = 1:nrow(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod$edge)), as.list), alpha = 0.01)
 EIWk15_C1_paranasalBranchRates$branch_test_results

 txtStop()


# ------------------------------------------------------------------- #
# 4.7 Extended implied weights (k=15) C3 paranasal character analysis #
# ------------------------------------------------------------------- #

### Find edge labels:
 plot (EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod)
 edgelabels()

 pdf("output/ParanasalAnalyses/EIWk15_C3/Figures/EIWk15_C3_EqualMethod_paranasalSinus_EdgeLabels.pdf", width=17, height=9) 
   plot (EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod)
   edgelabels()
 dev.off()

### Every edge in cladogram rates analysis:
 txtStart("output/ParanasalAnalyses/EIWk15_C3/Files/EIWk15_C3_EqualMethod_paranasalSinus_BranchRates.txt")

 EIWk15_C3_paranasalBranchRates <- test_rates(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, CrocParanasalNexus, time_bins = stage_bins, branch_partitions = lapply(X = as.list(x = 1:nrow(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod$edge)), as.list), alpha = 0.01)
 EIWk15_C3_paranasalBranchRates$branch_test_results

 txtStop()


# ------------------------------------------------------------------- #
# 4.8 Extended implied weights (k=15) C4 paranasal character analysis #
# ------------------------------------------------------------------- #

### Find edge labels:
 plot (EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod)
 edgelabels()

 pdf("output/ParanasalAnalyses/EIWk15_C4/Figures/EIWk15_C4_EqualMethod_paranasalSinus_EdgeLabels.pdf", width=17, height=9) 
   plot (EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod)
   edgelabels()
 dev.off()

### Every edge in cladogram rates analysis:
 txtStart("output/ParanasalAnalyses/EIWk15_C4/Files/EIWk15_C4_EqualMethod_paranasalSinus_BranchRates.txt")

 EIWk15_C4_paranasalBranchRates <- test_rates(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, CrocParanasalNexus, time_bins = stage_bins, branch_partitions = lapply(X = as.list(x = 1:nrow(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod$edge)), as.list), alpha = 0.01)
 EIWk15_C4_paranasalBranchRates$branch_test_results

 txtStop()

 beep(1)


# =========================================================== #
# 5. PCo analysis using Claddis, and generates PDFs and TIFFs #
# =========================================================== #

### Required libraries:
 library(beepr)
clusterCall(cluster, function () {
 library(phytools)
 library(Claddis)
})
 library(ggtree)
 library(ggplot2)
 library(TeachingDemos)

# --------------------------------- #
# 5.1 paranasal sinus PCoA analysis #
# --------------------------------- #

### Writes to file:
 txtStart("output/ParanasalAnalyses/NonPhylo/Files/Paranasal_PCoA.txt")

### Perform a PCoA ordination on the paranasal sinus dataset:
 SinPN_Pcoa_input<-ordinate_cladistic_matrix(cladistic_matrix = CrocParanasalNexus)

### Gives the eigenvalues, look at Reg_corr_eig and Cum_corr_eig for the % and cumulative % each axis explains:
 SinPN_Pcoa_input$values

### Gives the PCoA scores (###Check values to see how many axes contribute to 95%. Add that number in to the square brackets###):
 SinPN_Pcoa_scores<-SinPN_Pcoa_input$vectors[,1:18]

txtStop()


### Plot this as a simple bivarate morphospace:
 plot_morphospace(pcoa_input = SinPN_Pcoa_input)

### Add taxon names to morphospace:
 plot_morphospace(pcoa_input = SinPN_Pcoa_input, plot_taxon_names = TRUE)

### Define some simple taxon groups for the data as a named list:
SinPN_taxon_groups <- list('Sphenosuchia' = c("Junggarsuchus_sloani"),
  'Protosuchia' = c("Protosuchus_haughtoni"), 
  Notosuchia = c("cf_Hamadasuchus_rebouli","Campinasuchus_dinizi"),
Eusuchia = c("Lohuecosuchus_megadontos","Gavialis_gangeticus","Tomistoma_schlegelii","Osteolaemus_tetraspis","Mecistops_cataphractus","Crocodylus_acutus","Crocodylus_johnstoni","Crocodylus_moreletii","Crocodylus_niloticus","Crocodylus_porosus","Crocodylus_rhombifer","Alligator_mississippiensis","Alligator_sinensis","Caiman_crocodilus","Paleosuchus_palpebrosus"),
  BasalThalattosuchia = c("Plagiophthalmosuchus_gracilirostris","Macrospondylus_bollensis","Pelagosaurus_typus", "Eoneustes_gaudryi"),
  Metriorhynchidae = c("Cricosaurus_schroederi","Cricosaurus_araucanensis","Thalattosuchus_superciliosus"))


 class(SinPN_taxon_groups) <- "taxonGroups"


### Plot taxon groups including convex hulls:
 plot_morphospace(pcoa_input = SinPN_Pcoa_input, plot_taxon_names = TRUE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE)


### Generates PDFs of the PCoA ordination:
 pdf("output/ParanasalAnalyses/NonPhylo/Figures/Paranasal_PCoA_morphospace_no_taxnames_or_hulls.pdf", width=17, height=9) 
   plot_morphospace(pcoa_input = SinPN_Pcoa_input, plot_taxon_names = FALSE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = FALSE, plot_group_legend = FALSE)
 dev.off()

 pdf("output/ParanasalAnalyses/NonPhylo/Figures/Paranasal_PCoA_morphospace_no_taxnames.pdf", width=17, height=9) 
   plot_morphospace(pcoa_input = SinPN_Pcoa_input, plot_taxon_names = FALSE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE)
 dev.off()

 pdf("output/ParanasalAnalyses/NonPhylo/Figures/Paranasal_PCoA_morphospace.pdf", width=17, height=9) 
   plot_morphospace(pcoa_input = SinPN_Pcoa_input, plot_taxon_names = TRUE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE)
 dev.off()


### Generates high-resolution TIFF, using ggplot2, of the PCoA ordination:
 ggsave("output/ParanasalAnalyses/NonPhylo/Figures/Paranasal_PCoA_morphospace.tiff",
   plot_morphospace(pcoa_input = SinPN_Pcoa_input, plot_taxon_names = TRUE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)



# -------------------------------------------------------------------- #
# 5.2 paranasal sinus PCoA phylomorphospace analysis: Equal weights C0 #
# -------------------------------------------------------------------- #

### Make new ordination with tree included (enabling phylomorphospace):
 EW_C0_SinPN_phylomorphospace_pcoa_input <- ordinate_cladistic_matrix(cladistic_matrix = CrocParanasalNexus, time_tree = EW_C0_phylogeny_SinPN_Pruned_EqualMethod)

### Plot this as a simple bivarate phylomorphospace:
 plot_morphospace(pcoa_input = EW_C0_SinPN_phylomorphospace_pcoa_input)

### Add taxon names as well:
 plot_morphospace(pcoa_input = EW_C0_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE)

### Add taxon groups including convex hulls:
 plot_morphospace(pcoa_input = EW_C0_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE)


### Generates PDFs of the PCoA phylomorphospace ordination:
 pdf("output/ParanasalAnalyses/EW_C0/Figures/EW_C0_EqualMethod_paranasal_PCoA_phylomorphospace.pdf", width=17, height=9) 
   plot_morphospace(pcoa_input = EW_C0_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE)
 dev.off()

 pdf("output/ParanasalAnalyses/EW_C0/Figures/EW_C0_EqualMethod_paranasal_notaxnames_PCoA_phylomorphospace.pdf", width=17, height=9) 
   plot_morphospace(pcoa_input = EW_C0_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = FALSE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE)
 dev.off()

### Generates high-resolution TIFF, using ggplot2, of the PCoA phylomorphospace ordination:
 ggsave("output/ParanasalAnalyses/EW_C0/Figures/EW_C0_EqualMethod_paranasal_PCoA_phylomorphospace.tiff",
   plot_morphospace(pcoa_input = EW_C0_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)



# -------------------------------------------------------------------- #
# 5.3 paranasal sinus PCoA phylomorphospace analysis: Equal weights C1 #
# -------------------------------------------------------------------- #

### Make new ordination with tree included (enabling phylomorphospace):
 EW_C1_SinPN_phylomorphospace_pcoa_input <- ordinate_cladistic_matrix(cladistic_matrix = CrocParanasalNexus, time_tree = EW_C1_phylogeny_SinPN_Pruned_EqualMethod)

### Plot this as a simple bivarate phylomorphospace:
 plot_morphospace(pcoa_input = EW_C1_SinPN_phylomorphospace_pcoa_input)

### Add taxon names as well:
 plot_morphospace(pcoa_input = EW_C1_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE)

### Add taxon groups including convex hulls:
 plot_morphospace(pcoa_input = EW_C1_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE)


### Generates PDFs of the PCoA phylomorphospace ordination:
 pdf("output/ParanasalAnalyses/EW_C1/Figures/EW_C1_EqualMethod_paranasal_PCoA_phylomorphospace.pdf", width=17, height=9) 
   plot_morphospace(pcoa_input = EW_C1_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE)
 dev.off()

 pdf("output/ParanasalAnalyses/EW_C1/Figures/EW_C1_EqualMethod_paranasal_notaxnames_PCoA_phylomorphospace.pdf", width=17, height=9) 
   plot_morphospace(pcoa_input = EW_C1_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = FALSE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE)
 dev.off()

### Generates high-resolution TIFF, using ggplot2, of the PCoA phylomorphospace ordination:
 ggsave("output/ParanasalAnalyses/EW_C1/Figures/EW_C1_EqualMethod_paranasal_PCoA_phylomorphospace.tiff",
   plot_morphospace(pcoa_input = EW_C1_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)



# -------------------------------------------------------------------- #
# 5.4 paranasal sinus PCoA phylomorphospace analysis: Equal weights C3 #
# -------------------------------------------------------------------- #

### Make new ordination with tree included (enabling phylomorphospace):
 EW_C3_SinPN_phylomorphospace_pcoa_input <- ordinate_cladistic_matrix(cladistic_matrix = CrocParanasalNexus, time_tree = EW_C3_phylogeny_SinPN_Pruned_EqualMethod)

### Plot this as a simple bivarate phylomorphospace:
 plot_morphospace(pcoa_input = EW_C3_SinPN_phylomorphospace_pcoa_input)

### Add taxon names as well:
 plot_morphospace(pcoa_input = EW_C3_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE)

### Add taxon groups including convex hulls:
 plot_morphospace(pcoa_input = EW_C3_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE)


### Generates PDFs of the PCoA phylomorphospace ordination:
 pdf("output/ParanasalAnalyses/EW_C3/Figures/EW_C3_EqualMethod_paranasal_PCoA_phylomorphospace.pdf", width=17, height=9) 
   plot_morphospace(pcoa_input = EW_C3_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE)
 dev.off()

 pdf("output/ParanasalAnalyses/EW_C3/Figures/EW_C3_EqualMethod_paranasal_notaxnames_PCoA_phylomorphospace.pdf", width=17, height=9) 
   plot_morphospace(pcoa_input = EW_C3_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = FALSE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE)
 dev.off()

### Generates high-resolution TIFF, using ggplot2, of the PCoA phylomorphospace ordination:
 ggsave("output/ParanasalAnalyses/EW_C3/Figures/EW_C3_EqualMethod_paranasal_PCoA_phylomorphospace.tiff",
   plot_morphospace(pcoa_input = EW_C3_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)



# -------------------------------------------------------------------- #
# 5.5 paranasal sinus PCoA phylomorphospace analysis: Equal weights C4 #
# -------------------------------------------------------------------- #

### Make new ordination with tree included (enabling phylomorphospace):
 EW_C4_SinPN_phylomorphospace_pcoa_input <- ordinate_cladistic_matrix(cladistic_matrix = CrocParanasalNexus, time_tree = EW_C4_phylogeny_SinPN_Pruned_EqualMethod)

### Plot this as a simple bivarate phylomorphospace:
 plot_morphospace(pcoa_input = EW_C4_SinPN_phylomorphospace_pcoa_input)

### Add taxon names as well:
 plot_morphospace(pcoa_input = EW_C4_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE)

### Add taxon groups including convex hulls:
 plot_morphospace(pcoa_input = EW_C4_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE)


### Generates PDFs of the PCoA phylomorphospace ordination:
 pdf("output/ParanasalAnalyses/EW_C4/Figures/EW_C4_EqualMethod_paranasal_PCoA_phylomorphospace.pdf", width=17, height=9) 
   plot_morphospace(pcoa_input = EW_C4_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE)
 dev.off()

 pdf("output/ParanasalAnalyses/EW_C4/Figures/EW_C4_EqualMethod_paranasal_notaxnames_PCoA_phylomorphospace.pdf", width=17, height=9) 
   plot_morphospace(pcoa_input = EW_C4_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = FALSE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE)
 dev.off()

### Generates high-resolution TIFF, using ggplot2, of the PCoA phylomorphospace ordination:
 ggsave("output/ParanasalAnalyses/EW_C4/Figures/EW_C4_EqualMethod_paranasal_PCoA_phylomorphospace.tiff",
   plot_morphospace(pcoa_input = EW_C4_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)



# --------------------------------------------------------------- #
# 5.6 paranasal sinus PCoA phylomorphospace analysis: EIW k=15 C0 #
# --------------------------------------------------------------- #

### Make new ordination with tree included (enabling phylomorphospace):
 EIWk15_C0_SinPN_phylomorphospace_pcoa_input <- ordinate_cladistic_matrix(cladistic_matrix = CrocParanasalNexus, time_tree = EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod)

### Plot this as a simple bivarate phylomorphospace:
 plot_morphospace(pcoa_input = EIWk15_C0_SinPN_phylomorphospace_pcoa_input)

### Add taxon names as well:
 plot_morphospace(pcoa_input = EIWk15_C0_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE)

### Add taxon groups including convex hulls:
 plot_morphospace(pcoa_input = EIWk15_C0_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE)


### Generates PDFs of the PCoA phylomorphospace ordination:
 pdf("output/ParanasalAnalyses/EIWk15_C0/Figures/EIWk15_C0_EqualMethod_paranasal_PCoA_phylomorphospace.pdf", width=17, height=9) 
   plot_morphospace(pcoa_input = EIWk15_C0_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE)
 dev.off()

 pdf("output/ParanasalAnalyses/EIWk15_C0/Figures/EIWk15_C0_EqualMethod_paranasal_notaxnames_PCoA_phylomorphospace.pdf", width=17, height=9) 
   plot_morphospace(pcoa_input = EIWk15_C0_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = FALSE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE)
 dev.off()

### Generates high-resolution TIFF, using ggplot2, of the PCoA phylomorphospace ordination:
 ggsave("output/ParanasalAnalyses/EIWk15_C0/Figures/EIWk15_C0_EqualMethod_paranasal_PCoA_phylomorphospace.tiff",
   plot_morphospace(pcoa_input = EIWk15_C0_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)



# --------------------------------------------------------------- #
# 5.7 paranasal sinus PCoA phylomorphospace analysis: EIW k=15 C1 #
# --------------------------------------------------------------- #

### Make new ordination with tree included (enabling phylomorphospace):
 EIWk15_C1_SinPN_phylomorphospace_pcoa_input <- ordinate_cladistic_matrix(cladistic_matrix = CrocParanasalNexus, time_tree = EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod)

### Plot this as a simple bivarate phylomorphospace:
 plot_morphospace(pcoa_input = EIWk15_C1_SinPN_phylomorphospace_pcoa_input)

### Add taxon names as well:
 plot_morphospace(pcoa_input = EIWk15_C1_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE)

### Add taxon groups including convex hulls:
 plot_morphospace(pcoa_input = EIWk15_C1_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE)


### Generates PDFs of the PCoA phylomorphospace ordination:
 pdf("output/ParanasalAnalyses/EIWk15_C1/Figures/EIWk15_C1_EqualMethod_paranasal_PCoA_phylomorphospace.pdf", width=17, height=9) 
   plot_morphospace(pcoa_input = EIWk15_C1_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE)
 dev.off()

 pdf("output/ParanasalAnalyses/EIWk15_C1/Figures/EIWk15_C1_EqualMethod_paranasal_notaxnames_PCoA_phylomorphospace.pdf", width=17, height=9) 
   plot_morphospace(pcoa_input = EIWk15_C1_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = FALSE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE)
 dev.off()

### Generates high-resolution TIFF, using ggplot2, of the PCoA phylomorphospace ordination:
 ggsave("output/ParanasalAnalyses/EIWk15_C1/Figures/EIWk15_C1_EqualMethod_paranasal_PCoA_phylomorphospace.tiff",
   plot_morphospace(pcoa_input = EIWk15_C1_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)



# --------------------------------------------------------------- #
# 5.8 paranasal sinus PCoA phylomorphospace analysis: EIW k=15 C3 #
# --------------------------------------------------------------- #

### Make new ordination with tree included (enabling phylomorphospace):
 EIWk15_C3_SinPN_phylomorphospace_pcoa_input <- ordinate_cladistic_matrix(cladistic_matrix = CrocParanasalNexus, time_tree = EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod)

### Plot this as a simple bivarate phylomorphospace:
 plot_morphospace(pcoa_input = EIWk15_C3_SinPN_phylomorphospace_pcoa_input)

### Add taxon names as well:
 plot_morphospace(pcoa_input = EIWk15_C3_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE)

### Add taxon groups including convex hulls:
 plot_morphospace(pcoa_input = EIWk15_C3_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE)


### Generates PDFs of the PCoA phylomorphospace ordination:
 pdf("output/ParanasalAnalyses/EIWk15_C3/Figures/EIWk15_C3_EqualMethod_paranasal_PCoA_phylomorphospace.pdf", width=17, height=9) 
   plot_morphospace(pcoa_input = EIWk15_C3_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE)
 dev.off()

 pdf("output/ParanasalAnalyses/EIWk15_C3/Figures/EIWk15_C3_EqualMethod_paranasal_notaxnames_PCoA_phylomorphospace.pdf", width=17, height=9) 
   plot_morphospace(pcoa_input = EIWk15_C3_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = FALSE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE)
 dev.off()

### Generates high-resolution TIFF, using ggplot2, of the PCoA phylomorphospace ordination:
 ggsave("output/ParanasalAnalyses/EIWk15_C3/Figures/EIWk15_C3_EqualMethod_paranasal_PCoA_phylomorphospace.tiff",
   plot_morphospace(pcoa_input = EIWk15_C3_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)



# --------------------------------------------------------------- #
# 5.9 paranasal sinus PCoA phylomorphospace analysis: EIW k=15 C4 #
# --------------------------------------------------------------- #

### Make new ordination with tree included (enabling phylomorphospace):
 EIWk15_C4_SinPN_phylomorphospace_pcoa_input <- ordinate_cladistic_matrix(cladistic_matrix = CrocParanasalNexus, time_tree = EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod)

### Plot this as a simple bivarate phylomorphospace:
 plot_morphospace(pcoa_input = EIWk15_C4_SinPN_phylomorphospace_pcoa_input)

### Add taxon names as well:
 plot_morphospace(pcoa_input = EIWk15_C4_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE)

### Add taxon groups including convex hulls:
 plot_morphospace(pcoa_input = EIWk15_C4_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE)


### Generates PDFs of the PCoA phylomorphospace ordination:
 pdf("output/ParanasalAnalyses/EIWk15_C4/Figures/EIWk15_C4_EqualMethod_paranasal_PCoA_phylomorphospace.pdf", width=17, height=9) 
   plot_morphospace(pcoa_input = EIWk15_C4_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE)
 dev.off()

 pdf("output/ParanasalAnalyses/EIWk15_C4/Figures/EIWk15_C4_EqualMethod_paranasal_notaxnames_PCoA_phylomorphospace.pdf", width=17, height=9) 
   plot_morphospace(pcoa_input = EIWk15_C4_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = FALSE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE)
 dev.off()

### Generates high-resolution TIFF, using ggplot2, of the PCoA phylomorphospace ordination:
 ggsave("output/ParanasalAnalyses/EIWk15_C4/Figures/EIWk15_C4_EqualMethod_paranasal_PCoA_phylomorphospace.tiff",
   plot_morphospace(pcoa_input = EIWk15_C4_SinPN_phylomorphospace_pcoa_input, plot_taxon_names = TRUE, taxon_groups = SinPN_taxon_groups, plot_convex_hulls = TRUE, plot_group_legend = FALSE),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)

 beep(1)


# ======================================================================================= #
# 6. CV analysis using Morpho, using the morphospace scores, and generates PDFs and TIFFs #
# ======================================================================================= #

### Required libraries:
 library(beepr)
 library(Morpho)
 library(plot3D)
 library(TeachingDemos)


# ------------------- #
# 6.1 CVA input files #
# ------------------- #

### Paranasal-only PCo scores:
 PCo1_SinPN<-SinPN_Pcoa_scores[,1]
 PCo2_SinPN<-SinPN_Pcoa_scores[,2]
 PCo3_SinPN<-SinPN_Pcoa_scores[,3]
 PCos_SinPN<-SinPN_Pcoa_scores


# ------------------------------------ #
# 6.2 Paranasal-only sinus dataset CVA #
# ------------------------------------ #

### Writes to file:
 txtStart("output/ParanasalAnalyses/NonPhylo/Files/Paranasal_CVA.txt")


### CVA for PCo scores, Schwab 2020 habitat groups:
 SinPN_habitatSchwab_groups<-CrocSinPNData$HabitatS2020
 SinPN_habitatSchwab_groups<-as.factor(SinPN_habitatSchwab_groups)

 SinPN_cva_habitatSchwab<-CVA(PCos_SinPN, SinPN_habitatSchwab_groups, plot = T)
 SinPN_cva_habitatSchwab
 SinPN_cva_scores_habitatSchwab<-SinPN_cva_habitatSchwab$CVscores

 plot(SinPN_cva_habitatSchwab$CVscores,col=SinPN_habitatSchwab_groups,pch=16,cex=1.5)
   text(SinPN_cva_habitatSchwab$CVscores,rownames(PCos_SinPN),pos=3,cex=.6)
 legend(x="bottomleft", y=1, c("Pelagic", "Semiaquatic", "Terrestrial"), cex=.9,fill = mycolours)


### CVA for PCo scores, Wilberg 2019 habitat groups:
 SinPN_habitatWilberg_groups<-CrocSinPNData$HabitatW2019
 SinPN_habitatWilberg_groups <-as.factor(SinPN_habitatWilberg_groups)

 SinPN_cva_habitatWilberg<-CVA(PCos_SinPN, SinPN_habitatWilberg_groups, plot = T)
 SinPN_cva_habitatWilberg
 SinPN_cva_scores_habitatWilberg<-SinPN_cva_habitatWilberg$CVscores

 plot(SinPN_cva_habitatWilberg$CVscores,col=SinPN_habitatWilberg_groups,pch=16,cex=1.5)
   text(SinPN_cva_habitatWilberg$CVscores,rownames(PCos_SinPN),pos=3,cex=.6)
 legend(x="topright", y=1, c("Freshwater", "Marine", "Terrestrial"), cex=.9,fill = mycolours)


### CVA for PCo scores, Simple habitat groups:
 SinPN_habitatSimple_groups<-CrocSinPNData$HabitatSimple
 SinPN_habitatSimple_groups <-as.factor(SinPN_habitatSimple_groups)

 SinPN_cva_habitatSimple<-CVA(PCos_SinPN, SinPN_habitatSimple_groups, plot = T)
 SinPN_cva_habitatSimple
 SinPN_cva_scores_habitatSimple<-SinPN_cva_habitatSimple$CVscores

 plot(SinPN_cva_habitatSimple$CVscores,col=SinPN_habitatSimple_groups,pch=16,cex=1.5)
   text(SinPN_cva_habitatSimple$CVscores,rownames(PCos_SinPN),pos=3,cex=.6)
 legend(x="bottomright", y=1, c("Aquatic", "Terrestrial"), cex=.9,fill = mycolours)


### CVA for PCo scores, clade groups:
 SinPN_clade_groups<-CrocSinPNData$Clade
 SinPN_clade_groups <-as.factor(SinPN_clade_groups)

 SinPN_cva_clade<-CVA(PCos_SinPN, SinPN_clade_groups, plot = T)
 SinPN_cva_clade
 SinPN_cva_scores_clade<-SinPN_cva_clade$CVscores

 plot(SinPN_cva_clade$CVscores,col=SinPN_clade_groups,pch=16,cex=1.5)
   text(SinPN_cva_clade$CVscores,rownames(PCos_SinPN),pos=3,cex=.6)
 legend(x="bottomright", y=1, c("Neosuchia", "Notosuchia", "Protosuchia", "Sphenosuchia", "Thalattosuchia"), cex=.9,fill = mycolours)

 txtStop()


### Generates PDFS of the CVA ordinations:
 pdf("output/ParanasalAnalyses/NonPhylo/Figures/Paranasal_CVA_habitat_Schwab_text.pdf", width=17, height=9) 
   plot(SinPN_cva_habitatSchwab$CVscores,col=SinPN_habitatSchwab_groups,pch=16,cex=1.5)
   text(SinPN_cva_habitatSchwab$CVscores,rownames(PCos_SinPN),pos=3,cex=.6)
   legend(x="bottomleft", y=1, c("Pelagic", "Semiaquatic", "Terrestrial"), cex=.9,fill = mycolours)
 dev.off()

 pdf("output/ParanasalAnalyses/NonPhylo/Figures/Paranasal_CVA_habitat_Schwab_notext.pdf", width=17, height=9) 
   plot(SinPN_cva_habitatSchwab$CVscores,col=SinPN_habitatSchwab_groups,pch=16,cex=1.5)
   legend(x="bottomleft", y=1, c("Pelagic", "Semiaquatic", "Terrestrial"), cex=.9,fill = mycolours)
 dev.off()

 pdf("output/ParanasalAnalyses/NonPhylo/Figures/Paranasal_CVA_habitat_Wilberg_text.pdf", width=17, height=9) 
   plot(SinPN_cva_habitatWilberg$CVscores,col=SinPN_habitatWilberg_groups,pch=16,cex=1.5)
   text(SinPN_cva_habitatWilberg$CVscores,rownames(PCos_SinPN),pos=3,cex=.6)
   legend(x="topright", y=1, c("Freshwater", "Marine", "Terrestrial"), cex=.9,fill = mycolours)
 dev.off()

 pdf("output/ParanasalAnalyses/NonPhylo/Figures/Paranasal_CVA_habitat_Wilberg_notext.pdf", width=17, height=9) 
   plot(SinPN_cva_habitatWilberg$CVscores,col=SinPN_habitatWilberg_groups,pch=16,cex=1.5)
   legend(x="topright", y=1, c("Freshwater", "Marine", "Terrestrial"), cex=.9,fill = mycolours)
 dev.off()


 pdf("output/ParanasalAnalyses/NonPhylo/Figures/Paranasal_CVA_habitat_Simple_text.pdf", width=17, height=9) 
   plot(SinPN_cva_habitatSimple$CVscores,col=SinPN_habitatSimple_groups,pch=16,cex=1.5)
   text(SinPN_cva_habitatSimple$CVscores,rownames(PCos_SinPN),pos=3,cex=.6)
   legend(x="bottomright", y=1, c("Aquatic", "Terrestrial"), cex=.9,fill = mycolours)
 dev.off()

 pdf("output/ParanasalAnalyses/NonPhylo/Figures/Paranasal_CVA_habitat_Simple_notext.pdf", width=17, height=9) 
   plot(SinPN_cva_habitatSimple$CVscores,col=SinPN_habitatSimple_groups,pch=16,cex=1.5)
   legend(x="bottomright", y=1, c("Aquatic", "Terrestrial"), cex=.9,fill = mycolours)
 dev.off()

 pdf("output/ParanasalAnalyses/NonPhylo/Figures/Paranasal_CVA_clade_text.pdf", width=17, height=9) 
   plot(SinPN_cva_clade$CVscores,col=SinPN_clade_groups,pch=16,cex=1.5)
   text(SinPN_cva_clade$CVscores,rownames(PCos_SinPN),pos=3,cex=.6)
   legend(x="bottomright", y=1, c("Neosuchia", "Notosuchia", "Protosuchia", "Sphenosuchia", "Thalattosuchia"), cex=.9,fill = mycolours)
 dev.off()

 pdf("output/ParanasalAnalyses/NonPhylo/Figures/Paranasal_CVA_clade_notext.pdf", width=17, height=9) 
   plot(SinPN_cva_clade$CVscores,col=SinPN_clade_groups,pch=16,cex=1.5)
   legend(x="bottomright", y=1, c("Neosuchia", "Notosuchia", "Protosuchia", "Sphenosuchia", "Thalattosuchia"), cex=.9,fill = mycolours)
 dev.off()


### Generates high-resolution TIFFs, using ggplot2, of the CVA ordinations:
 library(ggtree)
 library(ggplot2)

 ggsave("output/ParanasalAnalyses/NonPhylo/Figures/Paranasal_CVA_habitat_Schwab.tiff",
   plot(SinPN_cva_habitatSchwab$CVscores,col=SinPN_habitatSchwab_groups,pch=16,cex=1.5),
   device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)
   legend(x="bottomleft", y=1, c("Pelagic", "Semiaquatic", "Terrestrial"), cex=.9,fill = mycolours)

ggsave("output/ParanasalAnalyses/NonPhylo/Figures/Paranasal_CVA_habitat_Wilberg.tiff",
 plot(SinPN_cva_habitatWilberg$CVscores,col=SinPN_habitatWilberg_groups,pch=16,cex=1.5),
   device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)
   legend(x="topright", y=1, c("Freshwater", "Marine", "Terrestrial"), cex=.9,fill = mycolours)

ggsave("output/ParanasalAnalyses/NonPhylo/Figures/Paranasal_CVA_habitat_Simple.tiff",
   plot(SinPN_cva_habitatSimple$CVscores,col=SinPN_habitatSimple_groups,pch=16,cex=1.5),
   device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)
   legend(x="bottomright", y=1, c("Aquatic", "Terrestrial"), cex=.9,fill = mycolours)

 ggsave("output/ParanasalAnalyses/NonPhylo/Figures/Paranasal_CVA_clade.tiff",
   plot(SinPN_cva_clade$CVscores,col=SinPN_clade_groups,pch=16,cex=1.5),
   device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)
   legend(x="bottomright", y=1, c("Neosuchia", "Notosuchia", "Protosuchia", "Sphenosuchia", "Thalattosuchia"), cex=.9,fill = mycolours)

 beep(1)


# ================================================================================ #
# 7. Creates contmaps and phenograms using Phytools, and generates PDFs and TIFFs  #
# ================================================================================ #

### Required libraries:
 library(beepr)
 library(phytools)
 library(ggtree)
 library(ggplot2)


# --------------------------------------- #
# 7.1a Equal weights C0 analysis contmaps #
# --------------------------------------- #

### contmap PCo1
EW_C0_SinPN_PCo1n<-setNames(PCo1_SinPN, rownames(CrocSinPNData))
range(EW_C0_SinPN_PCo1n)
fit<-fastAnc(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, EW_C0_SinPN_PCo1n,vars=TRUE,CI=TRUE)
fit
print(fit,printlen=10)
EW_C0_SinPN_PCo1_contmap<-contMap(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, EW_C0_SinPN_PCo1n,plot=FALSE, res=300)
EW_C0_SinPN_PCo1_contmap<-setMap(EW_C0_SinPN_PCo1_contmap, c("royalblue4", "slategray1", "white", "orange", "red4"))
plot(EW_C0_SinPN_PCo1_contmap, fsize=c(0.8), leg.txt="PCo1 values")

pdf("output/ParanasalAnalyses/EW_C0/Figures/EW_C0_paranasal_PCo1_contmap.pdf", width=17, height=9) 
plot(EW_C0_SinPN_PCo1_contmap)
dev.off()


### contmap PCo2
EW_C0_SinPN_PCo2n<-setNames(PCo2_SinPN, rownames(CrocSinPNData))
range(EW_C0_SinPN_PCo2n)
fit<-fastAnc(EW_C0_phylogeny_SinPN_Pruned_EqualMethod,EW_C0_SinPN_PCo2n,vars=TRUE,CI=TRUE)
fit
print(fit,printlen=10)
EW_C0_SinPN_PCo2_contmap<-contMap(EW_C0_phylogeny_SinPN_Pruned_EqualMethod,EW_C0_SinPN_PCo2n,plot=FALSE, res=300)
EW_C0_SinPN_PCo2_contmap<-setMap(EW_C0_SinPN_PCo2_contmap, c("royalblue4", "slategray1", "white", "orange", "red4"))
plot(EW_C0_SinPN_PCo2_contmap, fsize=c(0.8), leg.txt="PCo2 values")

pdf("output/ParanasalAnalyses/EW_C0/Figures/EW_C0_paranasal_PCo2_contmap.pdf", width=17, height=9) 
plot(EW_C0_SinPN_PCo2_contmap)
dev.off()


### contmap PCo3
EW_C0_SinPN_PCo3n<-setNames(PCo3_SinPN, rownames(CrocSinPNData))
range(EW_C0_SinPN_PCo3n)
fit<-fastAnc(EW_C0_phylogeny_SinPN_Pruned_EqualMethod,EW_C0_SinPN_PCo3n,vars=TRUE,CI=TRUE)
fit
print(fit,printlen=10)
EW_C0_SinPN_PCo3_contmap<-contMap(EW_C0_phylogeny_SinPN_Pruned_EqualMethod,EW_C0_SinPN_PCo3n,plot=FALSE, res=300)
EW_C0_SinPN_PCo3_contmap<-setMap(EW_C0_SinPN_PCo3_contmap, c("royalblue4", "slategray1", "white", "orange", "red4"))
plot(EW_C0_SinPN_PCo3_contmap, fsize=c(0.8), leg.txt="PCo3 values")

pdf("output/ParanasalAnalyses/EW_C0/Figures/EW_C0_paranasal_PCo3_contmap.pdf", width=17, height=9) 
plot(EW_C0_SinPN_PCo3_contmap)
dev.off()


# ----------------------------------------- #
# 7.1b Equal weights C0 analysis phenograms #
# ----------------------------------------- #

### phenogram PCo1 
EW_C0_SinPN_PCo1_phenogram<-phenogram(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, EW_C0_SinPN_PCo1n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0))

pdf("output/ParanasalAnalyses/EW_C0/Figures/EW_C0_paranasal_PCo1_phenogram.pdf", width=17, height=9) 
plot(phenogram(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, EW_C0_SinPN_PCo1n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0)))
dev.off()

ggsave("output/ParanasalAnalyses/EW_C0/Figures/EW_C0_paranasal_PCo1_phenogram.tiff",
plot(EW_C0_SinPN_PCo1_phenogram),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)


### phenogram PCo2
EW_C0_SinPN_PCo2_phenogram<-phenogram(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, EW_C0_SinPN_PCo2n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0))

pdf("output/ParanasalAnalyses/EW_C0/Figures/EW_C0_paranasal_PCo2_phenogram.pdf", width=17, height=9) 
plot(phenogram(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, EW_C0_SinPN_PCo2n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0)))
dev.off()

ggsave("output/ParanasalAnalyses/EW_C0/Figures/EW_C0_paranasal_PCo2_phenogram.tiff",
plot(EW_C0_SinPN_PCo2_phenogram),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)


### phenogram PCo3
EW_C0_SinPN_PCo3_phenogram<-phenogram(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, EW_C0_SinPN_PCo3n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0))

pdf("output/ParanasalAnalyses/EW_C0/Figures/EW_C0_paranasal_PCo3_phenogram.pdf", width=17, height=9) 
plot(phenogram(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, EW_C0_SinPN_PCo3n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0)))
dev.off()

ggsave("output/ParanasalAnalyses/EW_C0/Figures/EW_C0_paranasal_PCo3_phenogram.tiff",
plot(EW_C0_SinPN_PCo3_phenogram),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)


# --------------------------------------- #
# 7.2a Equal weights C1 analysis contmaps #
# --------------------------------------- #

### contmap PCo1
EW_C1_SinPN_PCo1n<-setNames(PCo1_SinPN, rownames(CrocSinPNData))
range(EW_C1_SinPN_PCo1n)
fit<-fastAnc(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, EW_C1_SinPN_PCo1n,vars=TRUE,CI=TRUE)
fit
print(fit,printlen=10)
EW_C1_SinPN_PCo1_contmap<-contMap(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, EW_C1_SinPN_PCo1n,plot=FALSE, res=300)
EW_C1_SinPN_PCo1_contmap<-setMap(EW_C1_SinPN_PCo1_contmap, c("royalblue4", "slategray1", "white", "orange", "red4"))
plot(EW_C1_SinPN_PCo1_contmap, fsize=c(0.8), leg.txt="PCo1 values")

pdf("output/ParanasalAnalyses/EW_C1/Figures/EW_C1_paranasal_PCo1_contmap.pdf", width=17, height=9) 
plot(EW_C1_SinPN_PCo1_contmap)
dev.off()


### contmap PCo2
EW_C1_SinPN_PCo2n<-setNames(PCo2_SinPN, rownames(CrocSinPNData))
range(EW_C1_SinPN_PCo2n)
fit<-fastAnc(EW_C1_phylogeny_SinPN_Pruned_EqualMethod,EW_C1_SinPN_PCo2n,vars=TRUE,CI=TRUE)
fit
print(fit,printlen=10)
EW_C1_SinPN_PCo2_contmap<-contMap(EW_C1_phylogeny_SinPN_Pruned_EqualMethod,EW_C1_SinPN_PCo2n,plot=FALSE, res=300)
EW_C1_SinPN_PCo2_contmap<-setMap(EW_C1_SinPN_PCo2_contmap, c("royalblue4", "slategray1", "white", "orange", "red4"))
plot(EW_C1_SinPN_PCo2_contmap, fsize=c(0.8), leg.txt="PCo2 values")

pdf("output/ParanasalAnalyses/EW_C1/Figures/EW_C1_paranasal_PCo2_contmap.pdf", width=17, height=9) 
plot(EW_C1_SinPN_PCo2_contmap)
dev.off()


### contmap PCo3
EW_C1_SinPN_PCo3n<-setNames(PCo3_SinPN, rownames(CrocSinPNData))
range(EW_C1_SinPN_PCo3n)
fit<-fastAnc(EW_C1_phylogeny_SinPN_Pruned_EqualMethod,EW_C1_SinPN_PCo3n,vars=TRUE,CI=TRUE)
fit
print(fit,printlen=10)
EW_C1_SinPN_PCo3_contmap<-contMap(EW_C1_phylogeny_SinPN_Pruned_EqualMethod,EW_C1_SinPN_PCo3n,plot=FALSE, res=300)
EW_C1_SinPN_PCo3_contmap<-setMap(EW_C1_SinPN_PCo3_contmap, c("royalblue4", "slategray1", "white", "orange", "red4"))
plot(EW_C1_SinPN_PCo3_contmap, fsize=c(0.8), leg.txt="PCo3 values")

pdf("output/ParanasalAnalyses/EW_C1/Figures/EW_C1_paranasal_PCo3_contmap.pdf", width=17, height=9) 
plot(EW_C1_SinPN_PCo3_contmap)
dev.off()


# ----------------------------------------- #
# 7.2b Equal weights C1 analysis phenograms #
# ----------------------------------------- #

### phenogram PCo1 
EW_C1_SinPN_PCo1_phenogram<-phenogram(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, EW_C1_SinPN_PCo1n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0))

pdf("output/ParanasalAnalyses/EW_C1/Figures/EW_C1_paranasal_PCo1_phenogram.pdf", width=17, height=9) 
plot(phenogram(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, EW_C1_SinPN_PCo1n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0)))
dev.off()

ggsave("output/ParanasalAnalyses/EW_C1/Figures/EW_C1_paranasal_PCo1_phenogram.tiff",
plot(EW_C1_SinPN_PCo1_phenogram),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)


### phenogram PCo2
EW_C1_SinPN_PCo2_phenogram<-phenogram(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, EW_C1_SinPN_PCo2n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0))

pdf("output/ParanasalAnalyses/EW_C1/Figures/EW_C1_paranasal_PCo2_phenogram.pdf", width=17, height=9) 
plot(phenogram(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, EW_C1_SinPN_PCo2n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0)))
dev.off()

ggsave("output/ParanasalAnalyses/EW_C1/Figures/EW_C1_paranasal_PCo2_phenogram.tiff",
plot(EW_C1_SinPN_PCo2_phenogram),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)


### phenogram PCo3
EW_C1_SinPN_PCo3_phenogram<-phenogram(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, EW_C1_SinPN_PCo3n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0))

pdf("output/ParanasalAnalyses/EW_C1/Figures/EW_C1_paranasal_PCo3_phenogram.pdf", width=17, height=9) 
plot(phenogram(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, EW_C1_SinPN_PCo3n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0)))
dev.off()

ggsave("output/ParanasalAnalyses/EW_C1/Figures/EW_C1_paranasal_PCo3_phenogram.tiff",
plot(EW_C1_SinPN_PCo3_phenogram),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)


# --------------------------------------- #
# 7.3a Equal weights C3 analysis contmaps #
# --------------------------------------- #

### contmap PCo1
EW_C3_SinPN_PCo1n<-setNames(PCo1_SinPN, rownames(CrocSinPNData))
range(EW_C3_SinPN_PCo1n)
fit<-fastAnc(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, EW_C3_SinPN_PCo1n,vars=TRUE,CI=TRUE)
fit
print(fit,printlen=10)
EW_C3_SinPN_PCo1_contmap<-contMap(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, EW_C3_SinPN_PCo1n,plot=FALSE, res=300)
EW_C3_SinPN_PCo1_contmap<-setMap(EW_C3_SinPN_PCo1_contmap, c("royalblue4", "slategray1", "white", "orange", "red4"))
plot(EW_C3_SinPN_PCo1_contmap, fsize=c(0.8), leg.txt="PCo1 values")

pdf("output/ParanasalAnalyses/EW_C3/Figures/EW_C3_paranasal_PCo1_contmap.pdf", width=17, height=9) 
plot(EW_C3_SinPN_PCo1_contmap)
dev.off()


### contmap PCo2
EW_C3_SinPN_PCo2n<-setNames(PCo2_SinPN, rownames(CrocSinPNData))
range(EW_C3_SinPN_PCo2n)
fit<-fastAnc(EW_C3_phylogeny_SinPN_Pruned_EqualMethod,EW_C3_SinPN_PCo2n,vars=TRUE,CI=TRUE)
fit
print(fit,printlen=10)
EW_C3_SinPN_PCo2_contmap<-contMap(EW_C3_phylogeny_SinPN_Pruned_EqualMethod,EW_C3_SinPN_PCo2n,plot=FALSE, res=300)
EW_C3_SinPN_PCo2_contmap<-setMap(EW_C3_SinPN_PCo2_contmap, c("royalblue4", "slategray1", "white", "orange", "red4"))
plot(EW_C3_SinPN_PCo2_contmap, fsize=c(0.8), leg.txt="PCo2 values")

pdf("output/ParanasalAnalyses/EW_C3/Figures/EW_C3_paranasal_PCo2_contmap.pdf", width=17, height=9) 
plot(EW_C3_SinPN_PCo2_contmap)
dev.off()


### contmap PCo3
EW_C3_SinPN_PCo3n<-setNames(PCo3_SinPN, rownames(CrocSinPNData))
range(EW_C3_SinPN_PCo3n)
fit<-fastAnc(EW_C3_phylogeny_SinPN_Pruned_EqualMethod,EW_C3_SinPN_PCo3n,vars=TRUE,CI=TRUE)
fit
print(fit,printlen=10)
EW_C3_SinPN_PCo3_contmap<-contMap(EW_C3_phylogeny_SinPN_Pruned_EqualMethod,EW_C3_SinPN_PCo3n,plot=FALSE, res=300)
EW_C3_SinPN_PCo3_contmap<-setMap(EW_C3_SinPN_PCo3_contmap, c("royalblue4", "slategray1", "white", "orange", "red4"))
plot(EW_C3_SinPN_PCo3_contmap, fsize=c(0.8), leg.txt="PCo3 values")

pdf("output/ParanasalAnalyses/EW_C3/Figures/EW_C3_paranasal_PCo3_contmap.pdf", width=17, height=9) 
plot(EW_C3_SinPN_PCo3_contmap)
dev.off()


# ----------------------------------------- #
# 7.3b Equal weights C3 analysis phenograms #
# ----------------------------------------- #

### phenogram PCo1 
EW_C3_SinPN_PCo1_phenogram<-phenogram(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, EW_C3_SinPN_PCo1n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0))

pdf("output/ParanasalAnalyses/EW_C3/Figures/EW_C3_paranasal_PCo1_phenogram.pdf", width=17, height=9) 
plot(phenogram(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, EW_C3_SinPN_PCo1n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0)))
dev.off()

ggsave("output/ParanasalAnalyses/EW_C3/Figures/EW_C3_paranasal_PCo1_phenogram.tiff",
plot(EW_C3_SinPN_PCo1_phenogram),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)


### phenogram PCo2
EW_C3_SinPN_PCo2_phenogram<-phenogram(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, EW_C3_SinPN_PCo2n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0))

pdf("output/ParanasalAnalyses/EW_C3/Figures/EW_C3_paranasal_PCo2_phenogram.pdf", width=17, height=9) 
plot(phenogram(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, EW_C3_SinPN_PCo2n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0)))
dev.off()

ggsave("output/ParanasalAnalyses/EW_C3/Figures/EW_C3_paranasal_PCo2_phenogram.tiff",
plot(EW_C3_SinPN_PCo2_phenogram),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)


### phenogram PCo3
EW_C3_SinPN_PCo3_phenogram<-phenogram(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, EW_C3_SinPN_PCo3n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0))

pdf("output/ParanasalAnalyses/EW_C3/Figures/EW_C3_paranasal_PCo3_phenogram.pdf", width=17, height=9) 
plot(phenogram(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, EW_C3_SinPN_PCo3n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0)))
dev.off()

ggsave("output/ParanasalAnalyses/EW_C3/Figures/EW_C3_paranasal_PCo3_phenogram.tiff",
plot(EW_C3_SinPN_PCo3_phenogram),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)


# --------------------------------------- #
# 7.4a Equal weights C4 analysis contmaps #
# --------------------------------------- #

### contmap PCo1
EW_C4_SinPN_PCo1n<-setNames(PCo1_SinPN, rownames(CrocSinPNData))
range(EW_C4_SinPN_PCo1n)
fit<-fastAnc(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, EW_C4_SinPN_PCo1n,vars=TRUE,CI=TRUE)
fit
print(fit,printlen=10)
EW_C4_SinPN_PCo1_contmap<-contMap(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, EW_C4_SinPN_PCo1n,plot=FALSE, res=300)
EW_C4_SinPN_PCo1_contmap<-setMap(EW_C4_SinPN_PCo1_contmap, c("royalblue4", "slategray1", "white", "orange", "red4"))
plot(EW_C4_SinPN_PCo1_contmap, fsize=c(0.8), leg.txt="PCo1 values")

pdf("output/ParanasalAnalyses/EW_C4/Figures/EW_C4_paranasal_PCo1_contmap.pdf", width=17, height=9) 
plot(EW_C4_SinPN_PCo1_contmap)
dev.off()


### contmap PCo2
EW_C4_SinPN_PCo2n<-setNames(PCo2_SinPN, rownames(CrocSinPNData))
range(EW_C4_SinPN_PCo2n)
fit<-fastAnc(EW_C4_phylogeny_SinPN_Pruned_EqualMethod,EW_C4_SinPN_PCo2n,vars=TRUE,CI=TRUE)
fit
print(fit,printlen=10)
EW_C4_SinPN_PCo2_contmap<-contMap(EW_C4_phylogeny_SinPN_Pruned_EqualMethod,EW_C4_SinPN_PCo2n,plot=FALSE, res=300)
EW_C4_SinPN_PCo2_contmap<-setMap(EW_C4_SinPN_PCo2_contmap, c("royalblue4", "slategray1", "white", "orange", "red4"))
plot(EW_C4_SinPN_PCo2_contmap, fsize=c(0.8), leg.txt="PCo2 values")

pdf("output/ParanasalAnalyses/EW_C4/Figures/EW_C4_paranasal_PCo2_contmap.pdf", width=17, height=9) 
plot(EW_C4_SinPN_PCo2_contmap)
dev.off()


### contmap PCo3
EW_C4_SinPN_PCo3n<-setNames(PCo3_SinPN, rownames(CrocSinPNData))
range(EW_C4_SinPN_PCo3n)
fit<-fastAnc(EW_C4_phylogeny_SinPN_Pruned_EqualMethod,EW_C4_SinPN_PCo3n,vars=TRUE,CI=TRUE)
fit
print(fit,printlen=10)
EW_C4_SinPN_PCo3_contmap<-contMap(EW_C4_phylogeny_SinPN_Pruned_EqualMethod,EW_C4_SinPN_PCo3n,plot=FALSE, res=300)
EW_C4_SinPN_PCo3_contmap<-setMap(EW_C4_SinPN_PCo3_contmap, c("royalblue4", "slategray1", "white", "orange", "red4"))
plot(EW_C4_SinPN_PCo3_contmap, fsize=c(0.8), leg.txt="PCo3 values")

pdf("output/ParanasalAnalyses/EW_C4/Figures/EW_C4_paranasal_PCo3_contmap.pdf", width=17, height=9) 
plot(EW_C4_SinPN_PCo3_contmap)
dev.off()


# ----------------------------------------- #
# 7.4b Equal weights C4 analysis phenograms #
# ----------------------------------------- #

### phenogram PCo1 
EW_C4_SinPN_PCo1_phenogram<-phenogram(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, EW_C4_SinPN_PCo1n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0))

pdf("output/ParanasalAnalyses/EW_C4/Figures/EW_C4_paranasal_PCo1_phenogram.pdf", width=17, height=9) 
plot(phenogram(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, EW_C4_SinPN_PCo1n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0)))
dev.off()

ggsave("output/ParanasalAnalyses/EW_C4/Figures/EW_C4_paranasal_PCo1_phenogram.tiff",
plot(EW_C4_SinPN_PCo1_phenogram),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)


### phenogram PCo2
EW_C4_SinPN_PCo2_phenogram<-phenogram(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, EW_C4_SinPN_PCo2n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0))

pdf("output/ParanasalAnalyses/EW_C4/Figures/EW_C4_paranasal_PCo2_phenogram.pdf", width=17, height=9) 
plot(phenogram(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, EW_C4_SinPN_PCo2n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0)))
dev.off()

ggsave("output/ParanasalAnalyses/EW_C4/Figures/EW_C4_paranasal_PCo2_phenogram.tiff",
plot(EW_C4_SinPN_PCo2_phenogram),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)


### phenogram PCo3
EW_C4_SinPN_PCo3_phenogram<-phenogram(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, EW_C4_SinPN_PCo3n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0))

pdf("output/ParanasalAnalyses/EW_C4/Figures/EW_C4_paranasal_PCo3_phenogram.pdf", width=17, height=9) 
plot(phenogram(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, EW_C4_SinPN_PCo3n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0)))
dev.off()

ggsave("output/ParanasalAnalyses/EW_C4/Figures/EW_C4_paranasal_PCo3_phenogram.tiff",
plot(EW_C4_SinPN_PCo3_phenogram),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)


# --------------------------------------------------------- #
# 7.5a Extended implied weights (k=15) C0 analysis contmaps #
# --------------------------------------------------------- #

### contmap PCo1
EIWk15_C0_SinPN_PCo1n<-setNames(PCo1_SinPN, rownames(CrocSinPNData))
range(EIWk15_C0_SinPN_PCo1n)
fit<-fastAnc(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C0_SinPN_PCo1n,vars=TRUE,CI=TRUE)
fit
print(fit,printlen=10)
EIWk15_C0_SinPN_PCo1_contmap<-contMap(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C0_SinPN_PCo1n,plot=FALSE, res=300)
EIWk15_C0_SinPN_PCo1_contmap<-setMap(EIWk15_C0_SinPN_PCo1_contmap, c("royalblue4", "slategray1", "white", "orange", "red4"))
plot(EIWk15_C0_SinPN_PCo1_contmap, fsize=c(0.8), leg.txt="PCo1 values")

pdf("output/ParanasalAnalyses/EIWk15_C0/Figures/EIWk15_C0_paranasal_PCo1_contmap.pdf", width=17, height=9) 
plot(EIWk15_C0_SinPN_PCo1_contmap)
dev.off()


### contmap PCo2
EIWk15_C0_SinPN_PCo2n<-setNames(PCo2_SinPN, rownames(CrocSinPNData))
range(EIWk15_C0_SinPN_PCo2n)
fit<-fastAnc(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod,EIWk15_C0_SinPN_PCo2n,vars=TRUE,CI=TRUE)
fit
print(fit,printlen=10)
EIWk15_C0_SinPN_PCo2_contmap<-contMap(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod,EIWk15_C0_SinPN_PCo2n,plot=FALSE, res=300)
EIWk15_C0_SinPN_PCo2_contmap<-setMap(EIWk15_C0_SinPN_PCo2_contmap, c("royalblue4", "slategray1", "white", "orange", "red4"))
plot(EIWk15_C0_SinPN_PCo2_contmap, fsize=c(0.8), leg.txt="PCo2 values")

pdf("output/ParanasalAnalyses/EIWk15_C0/Figures/EIWk15_C0_paranasal_PCo2_contmap.pdf", width=17, height=9) 
plot(EIWk15_C0_SinPN_PCo2_contmap)
dev.off()


### contmap PCo3
EIWk15_C0_SinPN_PCo3n<-setNames(PCo3_SinPN, rownames(CrocSinPNData))
range(EIWk15_C0_SinPN_PCo3n)
fit<-fastAnc(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod,EIWk15_C0_SinPN_PCo3n,vars=TRUE,CI=TRUE)
fit
print(fit,printlen=10)
EIWk15_C0_SinPN_PCo3_contmap<-contMap(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod,EIWk15_C0_SinPN_PCo3n,plot=FALSE, res=300)
EIWk15_C0_SinPN_PCo3_contmap<-setMap(EIWk15_C0_SinPN_PCo3_contmap, c("royalblue4", "slategray1", "white", "orange", "red4"))
plot(EIWk15_C0_SinPN_PCo3_contmap, fsize=c(0.8), leg.txt="PCo3 values")

pdf("output/ParanasalAnalyses/EIWk15_C0/Figures/EIWk15_C0_paranasal_PCo3_contmap.pdf", width=17, height=9) 
plot(EIWk15_C0_SinPN_PCo3_contmap)
dev.off()


# ----------------------------------------------------------- #
# 7.5b Extended implied weights (k=15) C0 analysis phenograms #
# ----------------------------------------------------------- #

### phenogram PCo1 
EIWk15_C0_SinPN_PCo1_phenogram<-phenogram(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C0_SinPN_PCo1n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0))

pdf("output/ParanasalAnalyses/EIWk15_C0/Figures/EIWk15_C0_paranasal_PCo1_phenogram.pdf", width=17, height=9) 
plot(phenogram(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C0_SinPN_PCo1n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0)))
dev.off()

ggsave("output/ParanasalAnalyses/EIWk15_C0/Figures/EIWk15_C0_paranasal_PCo1_phenogram.tiff",
plot(EIWk15_C0_SinPN_PCo1_phenogram),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)


### phenogram PCo2
EIWk15_C0_SinPN_PCo2_phenogram<-phenogram(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C0_SinPN_PCo2n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0))

pdf("output/ParanasalAnalyses/EIWk15_C0/Figures/EIWk15_C0_paranasal_PCo2_phenogram.pdf", width=17, height=9) 
plot(phenogram(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C0_SinPN_PCo2n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0)))
dev.off()

ggsave("output/ParanasalAnalyses/EIWk15_C0/Figures/EIWk15_C0_paranasal_PCo2_phenogram.tiff",
plot(EIWk15_C0_SinPN_PCo2_phenogram),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)


### phenogram PCo3
EIWk15_C0_SinPN_PCo3_phenogram<-phenogram(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C0_SinPN_PCo3n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0))

pdf("output/ParanasalAnalyses/EIWk15_C0/Figures/EIWk15_C0_paranasal_PCo3_phenogram.pdf", width=17, height=9) 
plot(phenogram(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C0_SinPN_PCo3n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0)))
dev.off()

ggsave("output/ParanasalAnalyses/EIWk15_C0/Figures/EIWk15_C0_paranasal_PCo3_phenogram.tiff",
plot(EIWk15_C0_SinPN_PCo3_phenogram),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)


# --------------------------------------------------------- #
# 7.6a Extended implied weights (k=15) C1 analysis contmaps #
# --------------------------------------------------------- #

### contmap PCo1
EIWk15_C1_SinPN_PCo1n<-setNames(PCo1_SinPN, rownames(CrocSinPNData))
range(EIWk15_C1_SinPN_PCo1n)
fit<-fastAnc(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C1_SinPN_PCo1n,vars=TRUE,CI=TRUE)
fit
print(fit,printlen=10)
EIWk15_C1_SinPN_PCo1_contmap<-contMap(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C1_SinPN_PCo1n,plot=FALSE, res=300)
EIWk15_C1_SinPN_PCo1_contmap<-setMap(EIWk15_C1_SinPN_PCo1_contmap, c("royalblue4", "slategray1", "white", "orange", "red4"))
plot(EIWk15_C1_SinPN_PCo1_contmap, fsize=c(0.8), leg.txt="PCo1 values")

pdf("output/ParanasalAnalyses/EIWk15_C1/Figures/EIWk15_C1_paranasal_PCo1_contmap.pdf", width=17, height=9) 
plot(EIWk15_C1_SinPN_PCo1_contmap)
dev.off()


### contmap PCo2
EIWk15_C1_SinPN_PCo2n<-setNames(PCo2_SinPN, rownames(CrocSinPNData))
range(EIWk15_C1_SinPN_PCo2n)
fit<-fastAnc(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod,EIWk15_C1_SinPN_PCo2n,vars=TRUE,CI=TRUE)
fit
print(fit,printlen=10)
EIWk15_C1_SinPN_PCo2_contmap<-contMap(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod,EIWk15_C1_SinPN_PCo2n,plot=FALSE, res=300)
EIWk15_C1_SinPN_PCo2_contmap<-setMap(EIWk15_C1_SinPN_PCo2_contmap, c("royalblue4", "slategray1", "white", "orange", "red4"))
plot(EIWk15_C1_SinPN_PCo2_contmap, fsize=c(0.8), leg.txt="PCo2 values")

pdf("output/ParanasalAnalyses/EIWk15_C1/Figures/EIWk15_C1_paranasal_PCo2_contmap.pdf", width=17, height=9) 
plot(EIWk15_C1_SinPN_PCo2_contmap)
dev.off()


### contmap PCo3
EIWk15_C1_SinPN_PCo3n<-setNames(PCo3_SinPN, rownames(CrocSinPNData))
range(EIWk15_C1_SinPN_PCo3n)
fit<-fastAnc(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod,EIWk15_C1_SinPN_PCo3n,vars=TRUE,CI=TRUE)
fit
print(fit,printlen=10)
EIWk15_C1_SinPN_PCo3_contmap<-contMap(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod,EIWk15_C1_SinPN_PCo3n,plot=FALSE, res=300)
EIWk15_C1_SinPN_PCo3_contmap<-setMap(EIWk15_C1_SinPN_PCo3_contmap, c("royalblue4", "slategray1", "white", "orange", "red4"))
plot(EIWk15_C1_SinPN_PCo3_contmap, fsize=c(0.8), leg.txt="PCo3 values")

pdf("output/ParanasalAnalyses/EIWk15_C1/Figures/EIWk15_C1_paranasal_PCo3_contmap.pdf", width=17, height=9) 
plot(EIWk15_C1_SinPN_PCo3_contmap)
dev.off()


# ----------------------------------------------------------- #
# 7.6b Extended implied weights (k=15) C1 analysis phenograms #
# ----------------------------------------------------------- #

### phenogram PCo1 
EIWk15_C1_SinPN_PCo1_phenogram<-phenogram(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C1_SinPN_PCo1n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0))

pdf("output/ParanasalAnalyses/EIWk15_C1/Figures/EIWk15_C1_paranasal_PCo1_phenogram.pdf", width=17, height=9) 
plot(phenogram(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C1_SinPN_PCo1n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0)))
dev.off()

ggsave("output/ParanasalAnalyses/EIWk15_C1/Figures/EIWk15_C1_paranasal_PCo1_phenogram.tiff",
plot(EIWk15_C1_SinPN_PCo1_phenogram),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)


### phenogram PCo2
EIWk15_C1_SinPN_PCo2_phenogram<-phenogram(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C1_SinPN_PCo2n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0))

pdf("output/ParanasalAnalyses/EIWk15_C1/Figures/EIWk15_C1_paranasal_PCo2_phenogram.pdf", width=17, height=9) 
plot(phenogram(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C1_SinPN_PCo2n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0)))
dev.off()

ggsave("output/ParanasalAnalyses/EIWk15_C1/Figures/EIWk15_C1_paranasal_PCo2_phenogram.tiff",
plot(EIWk15_C1_SinPN_PCo2_phenogram),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)


### phenogram PCo3
EIWk15_C1_SinPN_PCo3_phenogram<-phenogram(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C1_SinPN_PCo3n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0))

pdf("output/ParanasalAnalyses/EIWk15_C1/Figures/EIWk15_C1_paranasal_PCo3_phenogram.pdf", width=17, height=9) 
plot(phenogram(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C1_SinPN_PCo3n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0)))
dev.off()

ggsave("output/ParanasalAnalyses/EIWk15_C1/Figures/EIWk15_C1_paranasal_PCo3_phenogram.tiff",
plot(EIWk15_C1_SinPN_PCo3_phenogram),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)


# --------------------------------------------------------- #
# 7.7a Extended implied weights (k=15) C3 analysis contmaps #
# --------------------------------------------------------- #

### contmap PCo1
EIWk15_C3_SinPN_PCo1n<-setNames(PCo1_SinPN, rownames(CrocSinPNData))
range(EIWk15_C3_SinPN_PCo1n)
fit<-fastAnc(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C3_SinPN_PCo1n,vars=TRUE,CI=TRUE)
fit
print(fit,printlen=10)
EIWk15_C3_SinPN_PCo1_contmap<-contMap(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C3_SinPN_PCo1n,plot=FALSE, res=300)
EIWk15_C3_SinPN_PCo1_contmap<-setMap(EIWk15_C3_SinPN_PCo1_contmap, c("royalblue4", "slategray1", "white", "orange", "red4"))
plot(EIWk15_C3_SinPN_PCo1_contmap, fsize=c(0.8), leg.txt="PCo1 values")

pdf("output/ParanasalAnalyses/EIWk15_C3/Figures/EIWk15_C3_paranasal_PCo1_contmap.pdf", width=17, height=9) 
plot(EIWk15_C3_SinPN_PCo1_contmap)
dev.off()


### contmap PCo2
EIWk15_C3_SinPN_PCo2n<-setNames(PCo2_SinPN, rownames(CrocSinPNData))
range(EIWk15_C3_SinPN_PCo2n)
fit<-fastAnc(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod,EIWk15_C3_SinPN_PCo2n,vars=TRUE,CI=TRUE)
fit
print(fit,printlen=10)
EIWk15_C3_SinPN_PCo2_contmap<-contMap(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod,EIWk15_C3_SinPN_PCo2n,plot=FALSE, res=300)
EIWk15_C3_SinPN_PCo2_contmap<-setMap(EIWk15_C3_SinPN_PCo2_contmap, c("royalblue4", "slategray1", "white", "orange", "red4"))
plot(EIWk15_C3_SinPN_PCo2_contmap, fsize=c(0.8), leg.txt="PCo2 values")

pdf("output/ParanasalAnalyses/EIWk15_C3/Figures/EIWk15_C3_paranasal_PCo2_contmap.pdf", width=17, height=9) 
plot(EIWk15_C3_SinPN_PCo2_contmap)
dev.off()


### contmap PCo3
EIWk15_C3_SinPN_PCo3n<-setNames(PCo3_SinPN, rownames(CrocSinPNData))
range(EIWk15_C3_SinPN_PCo3n)
fit<-fastAnc(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod,EIWk15_C3_SinPN_PCo3n,vars=TRUE,CI=TRUE)
fit
print(fit,printlen=10)
EIWk15_C3_SinPN_PCo3_contmap<-contMap(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod,EIWk15_C3_SinPN_PCo3n,plot=FALSE, res=300)
EIWk15_C3_SinPN_PCo3_contmap<-setMap(EIWk15_C3_SinPN_PCo3_contmap, c("royalblue4", "slategray1", "white", "orange", "red4"))
plot(EIWk15_C3_SinPN_PCo3_contmap, fsize=c(0.8), leg.txt="PCo3 values")

pdf("output/ParanasalAnalyses/EIWk15_C3/Figures/EIWk15_C3_paranasal_PCo3_contmap.pdf", width=17, height=9) 
plot(EIWk15_C3_SinPN_PCo3_contmap)
dev.off()


# ----------------------------------------------------------- #
# 7.7b Extended implied weights (k=15) C3 analysis phenograms #
# ----------------------------------------------------------- #

### phenogram PCo1 
EIWk15_C3_SinPN_PCo1_phenogram<-phenogram(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C3_SinPN_PCo1n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0))

pdf("output/ParanasalAnalyses/EIWk15_C3/Figures/EIWk15_C3_paranasal_PCo1_phenogram.pdf", width=17, height=9) 
plot(phenogram(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C3_SinPN_PCo1n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0)))
dev.off()

ggsave("output/ParanasalAnalyses/EIWk15_C3/Figures/EIWk15_C3_paranasal_PCo1_phenogram.tiff",
plot(EIWk15_C3_SinPN_PCo1_phenogram),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)


### phenogram PCo2
EIWk15_C3_SinPN_PCo2_phenogram<-phenogram(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C3_SinPN_PCo2n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0))

pdf("output/ParanasalAnalyses/EIWk15_C3/Figures/EIWk15_C3_paranasal_PCo2_phenogram.pdf", width=17, height=9) 
plot(phenogram(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C3_SinPN_PCo2n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0)))
dev.off()

ggsave("output/ParanasalAnalyses/EIWk15_C3/Figures/EIWk15_C3_paranasal_PCo2_phenogram.tiff",
plot(EIWk15_C3_SinPN_PCo2_phenogram),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)


### phenogram PCo3
EIWk15_C3_SinPN_PCo3_phenogram<-phenogram(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C3_SinPN_PCo3n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0))

pdf("output/ParanasalAnalyses/EIWk15_C3/Figures/EIWk15_C3_paranasal_PCo3_phenogram.pdf", width=17, height=9) 
plot(phenogram(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C3_SinPN_PCo3n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0)))
dev.off()

ggsave("output/ParanasalAnalyses/EIWk15_C3/Figures/EIWk15_C3_paranasal_PCo3_phenogram.tiff",
plot(EIWk15_C3_SinPN_PCo3_phenogram),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)


# --------------------------------------------------------- #
# 7.8a Extended implied weights (k=15) C4 analysis contmaps #
# --------------------------------------------------------- #

### contmap PCo1
EIWk15_C4_SinPN_PCo1n<-setNames(PCo1_SinPN, rownames(CrocSinPNData))
range(EIWk15_C4_SinPN_PCo1n)
fit<-fastAnc(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C4_SinPN_PCo1n,vars=TRUE,CI=TRUE)
fit
print(fit,printlen=10)
EIWk15_C4_SinPN_PCo1_contmap<-contMap(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C4_SinPN_PCo1n,plot=FALSE, res=300)
EIWk15_C4_SinPN_PCo1_contmap<-setMap(EIWk15_C4_SinPN_PCo1_contmap, c("royalblue4", "slategray1", "white", "orange", "red4"))
plot(EIWk15_C4_SinPN_PCo1_contmap, fsize=c(0.8), leg.txt="PCo1 values")

pdf("output/ParanasalAnalyses/EIWk15_C4/Figures/EIWk15_C4_paranasal_PCo1_contmap.pdf", width=17, height=9) 
plot(EIWk15_C4_SinPN_PCo1_contmap)
dev.off()


### contmap PCo2
EIWk15_C4_SinPN_PCo2n<-setNames(PCo2_SinPN, rownames(CrocSinPNData))
range(EIWk15_C4_SinPN_PCo2n)
fit<-fastAnc(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod,EIWk15_C4_SinPN_PCo2n,vars=TRUE,CI=TRUE)
fit
print(fit,printlen=10)
EIWk15_C4_SinPN_PCo2_contmap<-contMap(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod,EIWk15_C4_SinPN_PCo2n,plot=FALSE, res=300)
EIWk15_C4_SinPN_PCo2_contmap<-setMap(EIWk15_C4_SinPN_PCo2_contmap, c("royalblue4", "slategray1", "white", "orange", "red4"))
plot(EIWk15_C4_SinPN_PCo2_contmap, fsize=c(0.8), leg.txt="PCo2 values")

pdf("output/ParanasalAnalyses/EIWk15_C4/Figures/EIWk15_C4_paranasal_PCo2_contmap.pdf", width=17, height=9) 
plot(EIWk15_C4_SinPN_PCo2_contmap)
dev.off()


### contmap PCo3
EIWk15_C4_SinPN_PCo3n<-setNames(PCo3_SinPN, rownames(CrocSinPNData))
range(EIWk15_C4_SinPN_PCo3n)
fit<-fastAnc(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod,EIWk15_C4_SinPN_PCo3n,vars=TRUE,CI=TRUE)
fit
print(fit,printlen=10)
EIWk15_C4_SinPN_PCo3_contmap<-contMap(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod,EIWk15_C4_SinPN_PCo3n,plot=FALSE, res=300)
EIWk15_C4_SinPN_PCo3_contmap<-setMap(EIWk15_C4_SinPN_PCo3_contmap, c("royalblue4", "slategray1", "white", "orange", "red4"))
plot(EIWk15_C4_SinPN_PCo3_contmap, fsize=c(0.8), leg.txt="PCo3 values")

pdf("output/ParanasalAnalyses/EIWk15_C4/Figures/EIWk15_C4_paranasal_PCo3_contmap.pdf", width=17, height=9) 
plot(EIWk15_C4_SinPN_PCo3_contmap)
dev.off()


# ----------------------------------------------------------- #
# 7.8b Extended implied weights (k=15) C4 analysis phenograms #
# ----------------------------------------------------------- #

### phenogram PCo1 
EIWk15_C4_SinPN_PCo1_phenogram<-phenogram(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C4_SinPN_PCo1n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0))

pdf("output/ParanasalAnalyses/EIWk15_C4/Figures/EIWk15_C4_paranasal_PCo1_phenogram.pdf", width=17, height=9) 
plot(phenogram(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C4_SinPN_PCo1n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0)))
dev.off()

ggsave("output/ParanasalAnalyses/EIWk15_C4/Figures/EIWk15_C4_paranasal_PCo1_phenogram.tiff",
plot(EIWk15_C4_SinPN_PCo1_phenogram),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)


### phenogram PCo2
EIWk15_C4_SinPN_PCo2_phenogram<-phenogram(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C4_SinPN_PCo2n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0))

pdf("output/ParanasalAnalyses/EIWk15_C4/Figures/EIWk15_C4_paranasal_PCo2_phenogram.pdf", width=17, height=9) 
plot(phenogram(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C4_SinPN_PCo2n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0)))
dev.off()

ggsave("output/ParanasalAnalyses/EIWk15_C4/Figures/EIWk15_C4_paranasal_PCo2_phenogram.tiff",
plot(EIWk15_C4_SinPN_PCo2_phenogram),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)


### phenogram PCo3
EIWk15_C4_SinPN_PCo3_phenogram<-phenogram(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C4_SinPN_PCo3n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0))

pdf("output/ParanasalAnalyses/EIWk15_C4/Figures/EIWk15_C4_paranasal_PCo3_phenogram.pdf", width=17, height=9) 
plot(phenogram(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, EIWk15_C4_SinPN_PCo3n,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0)))
dev.off()

ggsave("output/ParanasalAnalyses/EIWk15_C4/Figures/EIWk15_C4_paranasal_PCo3_phenogram.tiff",
plot(EIWk15_C4_SinPN_PCo3_phenogram),
 device = "tiff", dpi = 600, units = "cm", height = 17, width = 33)

 beep(1)


###################################
#                                 #
#  Multivariate statistical tests #
#                                 #
###################################


# ===================================== #
# 8. Phylogenetic Independent Contrasts #
# ===================================== #

### Required libraries:
 library(beepr)
 library(TeachingDemos)

# ----------------------------------- #
# 8.1 Equal weights C0 analysis: PICs #
# ----------------------------------- #

### writes to file:
txtStart("output/ParanasalAnalyses/EW_C0/Files/EW_C0_paranasal_PIC.txt")

 EW_C0_SinPN_PIC1<-pic(PCo1_SinPN, EW_C0_phylogeny_SinPN_Pruned_EqualMethod)
 EW_C0_SinPN_PIC1
 summary(EW_C0_SinPN_PIC1)
 plot(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, main= "PCo1", cex=0.5)
   nodelabels (round(EW_C0_SinPN_PIC1,5),adj=c(0,-0.5),frame="n")

 EW_C0_SinPN_PIC2<-pic(PCo2_SinPN, EW_C0_phylogeny_SinPN_Pruned_EqualMethod)
 EW_C0_SinPN_PIC2
 summary(EW_C0_SinPN_PIC2)
 plot(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, main= "PCo2", cex=0.5)
   nodelabels (round(EW_C0_SinPN_PIC2,5),adj=c(0,-0.5),frame="n")

 EW_C0_SinPN_PIC3<-pic(PCo3_SinPN, EW_C0_phylogeny_SinPN_Pruned_EqualMethod)
 EW_C0_SinPN_PIC3
 summary(EW_C0_SinPN_PIC3)
 plot(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, main= "PCo3", cex=0.5)
   nodelabels (round(EW_C0_SinPN_PIC3,5),adj=c(0,-0.5),frame="n")

### Schwab et al. 2020 habitat groupings:
 Habitat_SinSchwab<-as.factor(Habitat_SinSchwab)
 EW_C0_SinPN_PIC_SchwabHab<-pic(Habitat_SinSchwab, EW_C0_phylogeny_SinPN_Pruned_EqualMethod)
 EW_C0_SinPN_PIC_SchwabHab
 plot(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, main= "Habitat", cex=0.5)
   nodelabels (round(EW_C0_SinPN_PIC_SchwabHab,5),adj=c(0,-0.5),frame="n")

 EW_C0_SinPN_fit_PIC_SchwabHab<-lm(EW_C0_SinPN_PIC_SchwabHab ~ EW_C0_SinPN_PIC1+0)
 EW_C0_SinPN_fit_PIC_SchwabHab
 summary(EW_C0_SinPN_fit_PIC_SchwabHab)

 plot(EW_C0_SinPN_PIC1, EW_C0_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
 abline(EW_C0_SinPN_fit_PIC_SchwabHab,lwd=2,lty="dashed",col="red")

### Wilberg et al. 2019 habitat groupings:
 Habitat_SinWilberg <-as.factor(Habitat_SinWilberg)
 EW_C0_SinPN_PIC_WilbergHab<-pic(Habitat_SinWilberg, EW_C0_phylogeny_SinPN_Pruned_EqualMethod)
 EW_C0_SinPN_PIC_WilbergHab
 plot(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, main= "Habitat", cex=0.5)
   nodelabels (round(EW_C0_SinPN_PIC_WilbergHab,5),adj=c(0,-0.5),frame="n")

 EW_C0_SinPN_fit_PIC_WilbergHab<-lm(EW_C0_SinPN_PIC_WilbergHab ~ EW_C0_SinPN_PIC1 +0)
 EW_C0_SinPN_fit_PIC_WilbergHab
 summary(EW_C0_SinPN_fit_PIC_WilbergHab)

 plot(EW_C0_SinPN_PIC1, EW_C0_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
 abline(EW_C0_SinPN_fit_PIC_WilbergHab,lwd=2,lty="dashed",col="red")

txtStop()


### Generates PDFs:
pdf("output/ParanasalAnalyses/EW_C0/Figures/EW_C0_paranasal_PIC_SchwabHab_regression.pdf", width=17, height=9) 
plot(EW_C0_SinPN_PIC1, EW_C0_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
abline(EW_C0_SinPN_fit_PIC_SchwabHab,lwd=2,lty="dashed",col="red")
dev.off()

pdf("output/ParanasalAnalyses/EW_C0/Figures/EW_C0_paranasal_PIC_WilbergHab_regression.pdf", width=17, height=9) 
plot(EW_C0_SinPN_PIC1, EW_C0_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
abline(EW_C0_SinPN_fit_PIC_WilbergHab,lwd=2,lty="dashed",col="red")
dev.off()


# ----------------------------------- #
# 8.2 Equal weights C1 analysis: PICs #
# ----------------------------------- #

### Writes to file:
txtStart("output/ParanasalAnalyses/EW_C1/Files/EW_C1_paranasal_PIC.txt")

 EW_C1_SinPN_PIC1<-pic(PCo1_SinPN, EW_C1_phylogeny_SinPN_Pruned_EqualMethod)
 EW_C1_SinPN_PIC1
 summary(EW_C1_SinPN_PIC1)
 plot(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, main= "PCo1", cex=0.5)
   nodelabels (round(EW_C1_SinPN_PIC1,5),adj=c(0,-0.5),frame="n")

 EW_C1_SinPN_PIC2<-pic(PCo2_SinPN, EW_C1_phylogeny_SinPN_Pruned_EqualMethod)
 EW_C1_SinPN_PIC2
 summary(EW_C1_SinPN_PIC2)
 plot(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, main= "PCo2", cex=0.5)
   nodelabels (round(EW_C1_SinPN_PIC2,5),adj=c(0,-0.5),frame="n")

 EW_C1_SinPN_PIC3<-pic(PCo3_SinPN, EW_C1_phylogeny_SinPN_Pruned_EqualMethod)
 EW_C1_SinPN_PIC3
 summary(EW_C1_SinPN_PIC3)
 plot(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, main= "PCo3", cex=0.5)
   nodelabels (round(EW_C1_SinPN_PIC3,5),adj=c(0,-0.5),frame="n")

### Schwab et al. 2020 habitat groupings:
 Habitat_SinSchwab<-as.factor(Habitat_SinSchwab)
 EW_C1_SinPN_PIC_SchwabHab<-pic(Habitat_SinSchwab, EW_C1_phylogeny_SinPN_Pruned_EqualMethod)
 EW_C1_SinPN_PIC_SchwabHab
 plot(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, main= "Habitat", cex=0.5)
   nodelabels (round(EW_C1_SinPN_PIC_SchwabHab,5),adj=c(0,-0.5),frame="n")

 EW_C1_SinPN_fit_PIC_SchwabHab<-lm(EW_C1_SinPN_PIC_SchwabHab ~ EW_C1_SinPN_PIC1+0)
 EW_C1_SinPN_fit_PIC_SchwabHab
 summary(EW_C1_SinPN_fit_PIC_SchwabHab)

 plot(EW_C1_SinPN_PIC1, EW_C1_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
 abline(EW_C1_SinPN_fit_PIC_SchwabHab,lwd=2,lty="dashed",col="red")

### Wilberg et al. 2019 habitat groupings:
 Habitat_SinWilberg <-as.factor(Habitat_SinWilberg)
 EW_C1_SinPN_PIC_WilbergHab<-pic(Habitat_SinWilberg, EW_C1_phylogeny_SinPN_Pruned_EqualMethod)
 EW_C1_SinPN_PIC_WilbergHab
 plot(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, main= "Habitat", cex=0.5)
   nodelabels (round(EW_C1_SinPN_PIC_WilbergHab,5),adj=c(0,-0.5),frame="n")

 EW_C1_SinPN_fit_PIC_WilbergHab<-lm(EW_C1_SinPN_PIC_WilbergHab ~ EW_C1_SinPN_PIC1 +0)
 EW_C1_SinPN_fit_PIC_WilbergHab
 summary(EW_C1_SinPN_fit_PIC_WilbergHab)

 plot(EW_C1_SinPN_PIC1, EW_C1_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
 abline(EW_C1_SinPN_fit_PIC_WilbergHab,lwd=2,lty="dashed",col="red")

txtStop()


### Generates PDFs:
pdf("output/ParanasalAnalyses/EW_C1/Figures/EW_C1_paranasal_PIC_SchwabHab_regression.pdf", width=17, height=9) 
plot(EW_C1_SinPN_PIC1, EW_C1_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
abline(EW_C1_SinPN_fit_PIC_SchwabHab,lwd=2,lty="dashed",col="red")
dev.off()

pdf("output/ParanasalAnalyses/EW_C1/Figures/EW_C1_paranasal_PIC_WilbergHab_regression.pdf", width=17, height=9) 
plot(EW_C1_SinPN_PIC1, EW_C1_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
abline(EW_C1_SinPN_fit_PIC_WilbergHab,lwd=2,lty="dashed",col="red")
dev.off()


# ----------------------------------- #
# 8.3 Equal weights C3 analysis: PICs #
# ----------------------------------- #

### Writes to file:
txtStart("output/ParanasalAnalyses/EW_C3/Files/EW_C3_paranasal_PIC.txt")

 EW_C3_SinPN_PIC1<-pic(PCo1_SinPN, EW_C3_phylogeny_SinPN_Pruned_EqualMethod)
 EW_C3_SinPN_PIC1
 summary(EW_C3_SinPN_PIC1)
 plot(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, main= "PCo1", cex=0.5)
   nodelabels (round(EW_C3_SinPN_PIC1,5),adj=c(0,-0.5),frame="n")

 EW_C3_SinPN_PIC2<-pic(PCo2_SinPN, EW_C3_phylogeny_SinPN_Pruned_EqualMethod)
 EW_C3_SinPN_PIC2
 summary(EW_C3_SinPN_PIC2)
 plot(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, main= "PCo2", cex=0.5)
   nodelabels (round(EW_C3_SinPN_PIC2,5),adj=c(0,-0.5),frame="n")

 EW_C3_SinPN_PIC3<-pic(PCo3_SinPN, EW_C3_phylogeny_SinPN_Pruned_EqualMethod)
 EW_C3_SinPN_PIC3
 summary(EW_C3_SinPN_PIC3)
 plot(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, main= "PCo3", cex=0.5)
   nodelabels (round(EW_C3_SinPN_PIC3,5),adj=c(0,-0.5),frame="n")

### Schwab et al. 2020 habitat groupings:
 Habitat_SinSchwab<-as.factor(Habitat_SinSchwab)
 EW_C3_SinPN_PIC_SchwabHab<-pic(Habitat_SinSchwab, EW_C3_phylogeny_SinPN_Pruned_EqualMethod)
 EW_C3_SinPN_PIC_SchwabHab
 plot(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, main= "Habitat", cex=0.5)
   nodelabels (round(EW_C3_SinPN_PIC_SchwabHab,5),adj=c(0,-0.5),frame="n")

 EW_C3_SinPN_fit_PIC_SchwabHab<-lm(EW_C3_SinPN_PIC_SchwabHab ~ EW_C3_SinPN_PIC1+0)
 EW_C3_SinPN_fit_PIC_SchwabHab
 summary(EW_C3_SinPN_fit_PIC_SchwabHab)

 plot(EW_C3_SinPN_PIC1, EW_C3_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
 abline(EW_C3_SinPN_fit_PIC_SchwabHab,lwd=2,lty="dashed",col="red")

### Wilberg et al. 2019 habitat groupings:
 Habitat_SinWilberg <-as.factor(Habitat_SinWilberg)
 EW_C3_SinPN_PIC_WilbergHab<-pic(Habitat_SinWilberg, EW_C3_phylogeny_SinPN_Pruned_EqualMethod)
 EW_C3_SinPN_PIC_WilbergHab
 plot(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, main= "Habitat", cex=0.5)
   nodelabels (round(EW_C3_SinPN_PIC_WilbergHab,5),adj=c(0,-0.5),frame="n")

 EW_C3_SinPN_fit_PIC_WilbergHab<-lm(EW_C3_SinPN_PIC_WilbergHab ~ EW_C3_SinPN_PIC1 +0)
 EW_C3_SinPN_fit_PIC_WilbergHab
 summary(EW_C3_SinPN_fit_PIC_WilbergHab)

 plot(EW_C3_SinPN_PIC1, EW_C3_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
 abline(EW_C3_SinPN_fit_PIC_WilbergHab,lwd=2,lty="dashed",col="red")

txtStop()


### Generates PDFs:
pdf("output/ParanasalAnalyses/EW_C3/Figures/EW_C3_paranasal_PIC_SchwabHab_regression.pdf", width=17, height=9) 
plot(EW_C3_SinPN_PIC1, EW_C3_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
abline(EW_C3_SinPN_fit_PIC_SchwabHab,lwd=2,lty="dashed",col="red")
dev.off()

pdf("output/ParanasalAnalyses/EW_C3/Figures/EW_C3_paranasal_PIC_WilbergHab_regression.pdf", width=17, height=9) 
plot(EW_C3_SinPN_PIC1, EW_C3_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
abline(EW_C3_SinPN_fit_PIC_WilbergHab,lwd=2,lty="dashed",col="red")
dev.off()


# ----------------------------------- #
# 8.4 Equal weights C4 analysis: PICs #
# ----------------------------------- #

### Writes to file:
txtStart("output/ParanasalAnalyses/EW_C4/Files/EW_C4_paranasal_PIC.txt")

 EW_C4_SinPN_PIC1<-pic(PCo1_SinPN, EW_C4_phylogeny_SinPN_Pruned_EqualMethod)
 EW_C4_SinPN_PIC1
 summary(EW_C4_SinPN_PIC1)
 plot(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, main= "PCo1", cex=0.5)
   nodelabels (round(EW_C4_SinPN_PIC1,5),adj=c(0,-0.5),frame="n")

 EW_C4_SinPN_PIC2<-pic(PCo2_SinPN, EW_C4_phylogeny_SinPN_Pruned_EqualMethod)
 EW_C4_SinPN_PIC2
 summary(EW_C4_SinPN_PIC2)
 plot(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, main= "PCo2", cex=0.5)
   nodelabels (round(EW_C4_SinPN_PIC2,5),adj=c(0,-0.5),frame="n")

 EW_C4_SinPN_PIC3<-pic(PCo3_SinPN, EW_C4_phylogeny_SinPN_Pruned_EqualMethod)
 EW_C4_SinPN_PIC3
 summary(EW_C4_SinPN_PIC3)
 plot(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, main= "PCo3", cex=0.5)
   nodelabels (round(EW_C4_SinPN_PIC3,5),adj=c(0,-0.5),frame="n")

### Schwab et al. 2020 habitat groupings:
 Habitat_SinSchwab<-as.factor(Habitat_SinSchwab)
 EW_C4_SinPN_PIC_SchwabHab<-pic(Habitat_SinSchwab, EW_C4_phylogeny_SinPN_Pruned_EqualMethod)
 EW_C4_SinPN_PIC_SchwabHab
 plot(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, main= "Habitat", cex=0.5)
   nodelabels (round(EW_C4_SinPN_PIC_SchwabHab,5),adj=c(0,-0.5),frame="n")

 EW_C4_SinPN_fit_PIC_SchwabHab<-lm(EW_C4_SinPN_PIC_SchwabHab ~ EW_C4_SinPN_PIC1+0)
 EW_C4_SinPN_fit_PIC_SchwabHab
 summary(EW_C4_SinPN_fit_PIC_SchwabHab)

 plot(EW_C4_SinPN_PIC1, EW_C4_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
 abline(EW_C4_SinPN_fit_PIC_SchwabHab,lwd=2,lty="dashed",col="red")

### Wilberg et al. 2019 habitat groupings:
 Habitat_SinWilberg <-as.factor(Habitat_SinWilberg)
 EW_C4_SinPN_PIC_WilbergHab<-pic(Habitat_SinWilberg, EW_C4_phylogeny_SinPN_Pruned_EqualMethod)
 EW_C4_SinPN_PIC_WilbergHab
 plot(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, main= "Habitat", cex=0.5)
   nodelabels (round(EW_C4_SinPN_PIC_WilbergHab,5),adj=c(0,-0.5),frame="n")

 EW_C4_SinPN_fit_PIC_WilbergHab<-lm(EW_C4_SinPN_PIC_WilbergHab ~ EW_C4_SinPN_PIC1 +0)
 EW_C4_SinPN_fit_PIC_WilbergHab
 summary(EW_C4_SinPN_fit_PIC_WilbergHab)

 plot(EW_C4_SinPN_PIC1, EW_C4_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
 abline(EW_C4_SinPN_fit_PIC_WilbergHab,lwd=2,lty="dashed",col="red")

txtStop()


### Generates PDFs:
pdf("output/ParanasalAnalyses/EW_C4/Figures/EW_C4_paranasal_PIC_SchwabHab_regression.pdf", width=17, height=9) 
plot(EW_C4_SinPN_PIC1, EW_C4_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
abline(EW_C4_SinPN_fit_PIC_SchwabHab,lwd=2,lty="dashed",col="red")
dev.off()

pdf("output/ParanasalAnalyses/EW_C4/Figures/EW_C4_paranasal_PIC_WilbergHab_regression.pdf", width=17, height=9) 
plot(EW_C4_SinPN_PIC1, EW_C4_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
abline(EW_C4_SinPN_fit_PIC_WilbergHab,lwd=2,lty="dashed",col="red")
dev.off()


# ----------------------------------------------------- #
# 8.5 Extended implied weights (k=15) C0 analysis: PICs #
# ----------------------------------------------------- #

### Writes to file:
txtStart("output/ParanasalAnalyses/EIWk15_C0/Files/EIWk15_C0_paranasal_PIC.txt")

 EIWk15_C0_SinPN_PIC1<-pic(PCo1_SinPN, EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod)
 EIWk15_C0_SinPN_PIC1
 summary(EIWk15_C0_SinPN_PIC1)
 plot(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, main= "PCo1", cex=0.5)
   nodelabels (round(EIWk15_C0_SinPN_PIC1,5),adj=c(0,-0.5),frame="n")

 EIWk15_C0_SinPN_PIC2<-pic(PCo2_SinPN, EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod)
 EIWk15_C0_SinPN_PIC2
 summary(EIWk15_C0_SinPN_PIC2)
 plot(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, main= "PCo2", cex=0.5)
   nodelabels (round(EIWk15_C0_SinPN_PIC2,5),adj=c(0,-0.5),frame="n")

 EIWk15_C0_SinPN_PIC3<-pic(PCo3_SinPN, EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod)
 EIWk15_C0_SinPN_PIC3
 summary(EIWk15_C0_SinPN_PIC3)
 plot(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, main= "PCo3", cex=0.5)
   nodelabels (round(EIWk15_C0_SinPN_PIC3,5),adj=c(0,-0.5),frame="n")

### Schwab et al. 2020 habitat groupings:
 Habitat_SinSchwab<-as.factor(Habitat_SinSchwab)
 EIWk15_C0_SinPN_PIC_SchwabHab<-pic(Habitat_SinSchwab, EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod)
 EIWk15_C0_SinPN_PIC_SchwabHab
 plot(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, main= "Habitat", cex=0.5)
   nodelabels (round(EIWk15_C0_SinPN_PIC_SchwabHab,5),adj=c(0,-0.5),frame="n")

 EIWk15_C0_SinPN_fit_PIC_SchwabHab<-lm(EIWk15_C0_SinPN_PIC_SchwabHab ~ EIWk15_C0_SinPN_PIC1+0)
 EIWk15_C0_SinPN_fit_PIC_SchwabHab
 summary(EIWk15_C0_SinPN_fit_PIC_SchwabHab)

 plot(EIWk15_C0_SinPN_PIC1, EIWk15_C0_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
 abline(EIWk15_C0_SinPN_fit_PIC_SchwabHab,lwd=2,lty="dashed",col="red")

### Wilberg et al. 2019 habitat groupings:
 Habitat_SinWilberg <-as.factor(Habitat_SinWilberg)
 EIWk15_C0_SinPN_PIC_WilbergHab<-pic(Habitat_SinWilberg, EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod)
 EIWk15_C0_SinPN_PIC_WilbergHab
 plot(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, main= "Habitat", cex=0.5)
   nodelabels (round(EIWk15_C0_SinPN_PIC_WilbergHab,5),adj=c(0,-0.5),frame="n")

 EIWk15_C0_SinPN_fit_PIC_WilbergHab<-lm(EIWk15_C0_SinPN_PIC_WilbergHab ~ EIWk15_C0_SinPN_PIC1 +0)
 EIWk15_C0_SinPN_fit_PIC_WilbergHab
 summary(EIWk15_C0_SinPN_fit_PIC_WilbergHab)

 plot(EIWk15_C0_SinPN_PIC1, EIWk15_C0_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
 abline(EIWk15_C0_SinPN_fit_PIC_WilbergHab,lwd=2,lty="dashed",col="red")

txtStop()


### Generates PDFs:
pdf("output/ParanasalAnalyses/EIWk15_C0/Figures/EIWk15_C0_paranasal_PIC_SchwabHab_regression.pdf", width=17, height=9) 
plot(EIWk15_C0_SinPN_PIC1, EIWk15_C0_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
abline(EIWk15_C0_SinPN_fit_PIC_SchwabHab,lwd=2,lty="dashed",col="red")
dev.off()

pdf("output/ParanasalAnalyses/EIWk15_C0/Figures/EIWk15_C0_paranasal_PIC_WilbergHab_regression.pdf", width=17, height=9) 
plot(EIWk15_C0_SinPN_PIC1, EIWk15_C0_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
abline(EIWk15_C0_SinPN_fit_PIC_WilbergHab,lwd=2,lty="dashed",col="red")
dev.off()


# ----------------------------------------------------- #
# 8.6 Extended implied weights (k=15) C1 analysis: PICs #
# ----------------------------------------------------- #

### Writes to file:
txtStart("output/ParanasalAnalyses/EIWk15_C1/Files/EIWk15_C1_paranasal_PIC.txt")

 EIWk15_C1_SinPN_PIC1<-pic(PCo1_SinPN, EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod)
 EIWk15_C1_SinPN_PIC1
 summary(EIWk15_C1_SinPN_PIC1)
 plot(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, main= "PCo1", cex=0.5)
   nodelabels (round(EIWk15_C1_SinPN_PIC1,5),adj=c(0,-0.5),frame="n")

 EIWk15_C1_SinPN_PIC2<-pic(PCo2_SinPN, EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod)
 EIWk15_C1_SinPN_PIC2
 summary(EIWk15_C1_SinPN_PIC2)
 plot(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, main= "PCo2", cex=0.5)
   nodelabels (round(EIWk15_C1_SinPN_PIC2,5),adj=c(0,-0.5),frame="n")

 EIWk15_C1_SinPN_PIC3<-pic(PCo3_SinPN, EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod)
 EIWk15_C1_SinPN_PIC3
 summary(EIWk15_C1_SinPN_PIC3)
 plot(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, main= "PCo3", cex=0.5)
   nodelabels (round(EIWk15_C1_SinPN_PIC3,5),adj=c(0,-0.5),frame="n")

### Schwab et al. 2020 habitat groupings:
 Habitat_SinSchwab<-as.factor(Habitat_SinSchwab)
 EIWk15_C1_SinPN_PIC_SchwabHab<-pic(Habitat_SinSchwab, EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod)
 EIWk15_C1_SinPN_PIC_SchwabHab
 plot(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, main= "Habitat", cex=0.5)
   nodelabels (round(EIWk15_C1_SinPN_PIC_SchwabHab,5),adj=c(0,-0.5),frame="n")

 EIWk15_C1_SinPN_fit_PIC_SchwabHab<-lm(EIWk15_C1_SinPN_PIC_SchwabHab ~ EIWk15_C1_SinPN_PIC1+0)
 EIWk15_C1_SinPN_fit_PIC_SchwabHab
 summary(EIWk15_C1_SinPN_fit_PIC_SchwabHab)

 plot(EIWk15_C1_SinPN_PIC1, EIWk15_C1_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
 abline(EIWk15_C1_SinPN_fit_PIC_SchwabHab,lwd=2,lty="dashed",col="red")

### Wilberg et al. 2019 habitat groupings:
 Habitat_SinWilberg <-as.factor(Habitat_SinWilberg)
 EIWk15_C1_SinPN_PIC_WilbergHab<-pic(Habitat_SinWilberg, EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod)
 EIWk15_C1_SinPN_PIC_WilbergHab
 plot(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, main= "Habitat", cex=0.5)
   nodelabels (round(EIWk15_C1_SinPN_PIC_WilbergHab,5),adj=c(0,-0.5),frame="n")

 EIWk15_C1_SinPN_fit_PIC_WilbergHab<-lm(EIWk15_C1_SinPN_PIC_WilbergHab ~ EIWk15_C1_SinPN_PIC1 +0)
 EIWk15_C1_SinPN_fit_PIC_WilbergHab
 summary(EIWk15_C1_SinPN_fit_PIC_WilbergHab)

 plot(EIWk15_C1_SinPN_PIC1, EIWk15_C1_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
 abline(EIWk15_C1_SinPN_fit_PIC_WilbergHab,lwd=2,lty="dashed",col="red")

txtStop()


### Generates PDFs:
pdf("output/ParanasalAnalyses/EIWk15_C1/Figures/EIWk15_C1_paranasal_PIC_SchwabHab_regression.pdf", width=17, height=9) 
plot(EIWk15_C1_SinPN_PIC1, EIWk15_C1_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
abline(EIWk15_C1_SinPN_fit_PIC_SchwabHab,lwd=2,lty="dashed",col="red")
dev.off()

pdf("output/ParanasalAnalyses/EIWk15_C1/Figures/EIWk15_C1_paranasal_PIC_WilbergHab_regression.pdf", width=17, height=9) 
plot(EIWk15_C1_SinPN_PIC1, EIWk15_C1_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
abline(EIWk15_C1_SinPN_fit_PIC_WilbergHab,lwd=2,lty="dashed",col="red")
dev.off()


# ----------------------------------------------------- #
# 8.7 Extended implied weights (k=15) C3 analysis: PICs #
# ----------------------------------------------------- #

### Writes to file:
txtStart("output/ParanasalAnalyses/EIWk15_C3/Files/EIWk15_C3_paranasal_PIC.txt")

 EIWk15_C3_SinPN_PIC1<-pic(PCo1_SinPN, EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod)
 EIWk15_C3_SinPN_PIC1
 summary(EIWk15_C3_SinPN_PIC1)
 plot(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, main= "PCo1", cex=0.5)
   nodelabels (round(EIWk15_C3_SinPN_PIC1,5),adj=c(0,-0.5),frame="n")

 EIWk15_C3_SinPN_PIC2<-pic(PCo2_SinPN, EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod)
 EIWk15_C3_SinPN_PIC2
 summary(EIWk15_C3_SinPN_PIC2)
 plot(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, main= "PCo2", cex=0.5)
   nodelabels (round(EIWk15_C3_SinPN_PIC2,5),adj=c(0,-0.5),frame="n")

 EIWk15_C3_SinPN_PIC3<-pic(PCo3_SinPN, EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod)
 EIWk15_C3_SinPN_PIC3
 summary(EIWk15_C3_SinPN_PIC3)
 plot(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, main= "PCo3", cex=0.5)
   nodelabels (round(EIWk15_C3_SinPN_PIC3,5),adj=c(0,-0.5),frame="n")

### Schwab et al. 2020 habitat groupings:
 Habitat_SinSchwab<-as.factor(Habitat_SinSchwab)
 EIWk15_C3_SinPN_PIC_SchwabHab<-pic(Habitat_SinSchwab, EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod)
 EIWk15_C3_SinPN_PIC_SchwabHab
 plot(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, main= "Habitat", cex=0.5)
   nodelabels (round(EIWk15_C3_SinPN_PIC_SchwabHab,5),adj=c(0,-0.5),frame="n")

 EIWk15_C3_SinPN_fit_PIC_SchwabHab<-lm(EIWk15_C3_SinPN_PIC_SchwabHab ~ EIWk15_C3_SinPN_PIC1+0)
 EIWk15_C3_SinPN_fit_PIC_SchwabHab
 summary(EIWk15_C3_SinPN_fit_PIC_SchwabHab)

 plot(EIWk15_C3_SinPN_PIC1, EIWk15_C3_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
 abline(EIWk15_C3_SinPN_fit_PIC_SchwabHab,lwd=2,lty="dashed",col="red")

### Wilberg et al. 2019 habitat groupings:
 Habitat_SinWilberg <-as.factor(Habitat_SinWilberg)
 EIWk15_C3_SinPN_PIC_WilbergHab<-pic(Habitat_SinWilberg, EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod)
 EIWk15_C3_SinPN_PIC_WilbergHab
 plot(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, main= "Habitat", cex=0.5)
   nodelabels (round(EIWk15_C3_SinPN_PIC_WilbergHab,5),adj=c(0,-0.5),frame="n")

 EIWk15_C3_SinPN_fit_PIC_WilbergHab<-lm(EIWk15_C3_SinPN_PIC_WilbergHab ~ EIWk15_C3_SinPN_PIC1 +0)
 EIWk15_C3_SinPN_fit_PIC_WilbergHab
 summary(EIWk15_C3_SinPN_fit_PIC_WilbergHab)

 plot(EIWk15_C3_SinPN_PIC1, EIWk15_C3_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
 abline(EIWk15_C3_SinPN_fit_PIC_WilbergHab,lwd=2,lty="dashed",col="red")

txtStop()


### Generates PDFs:
pdf("output/ParanasalAnalyses/EIWk15_C3/Figures/EIWk15_C3_paranasal_PIC_SchwabHab_regression.pdf", width=17, height=9) 
plot(EIWk15_C3_SinPN_PIC1, EIWk15_C3_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
abline(EIWk15_C3_SinPN_fit_PIC_SchwabHab,lwd=2,lty="dashed",col="red")
dev.off()

pdf("output/ParanasalAnalyses/EIWk15_C3/Figures/EIWk15_C3_paranasal_PIC_WilbergHab_regression.pdf", width=17, height=9) 
plot(EIWk15_C3_SinPN_PIC1, EIWk15_C3_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
abline(EIWk15_C3_SinPN_fit_PIC_WilbergHab,lwd=2,lty="dashed",col="red")
dev.off()


# ----------------------------------------------------- #
# 8.8 Extended implied weights (k=15) C4 analysis: PICs #
# ----------------------------------------------------- #

### Writes to file:
txtStart("output/ParanasalAnalyses/EIWk15_C4/Files/EIWk15_C4_paranasal_PIC.txt")

 EIWk15_C4_SinPN_PIC1<-pic(PCo1_SinPN, EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod)
 EIWk15_C4_SinPN_PIC1
 summary(EIWk15_C4_SinPN_PIC1)
 plot(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, main= "PCo1", cex=0.5)
   nodelabels (round(EIWk15_C4_SinPN_PIC1,5),adj=c(0,-0.5),frame="n")

 EIWk15_C4_SinPN_PIC2<-pic(PCo2_SinPN, EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod)
 EIWk15_C4_SinPN_PIC2
 summary(EIWk15_C4_SinPN_PIC2)
 plot(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, main= "PCo2", cex=0.5)
   nodelabels (round(EIWk15_C4_SinPN_PIC2,5),adj=c(0,-0.5),frame="n")

 EIWk15_C4_SinPN_PIC3<-pic(PCo3_SinPN, EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod)
 EIWk15_C4_SinPN_PIC3
 summary(EIWk15_C4_SinPN_PIC3)
 plot(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, main= "PCo3", cex=0.5)
   nodelabels (round(EIWk15_C4_SinPN_PIC3,5),adj=c(0,-0.5),frame="n")

### Schwab et al. 2020 habitat groupings:
 Habitat_SinSchwab<-as.factor(Habitat_SinSchwab)
 EIWk15_C4_SinPN_PIC_SchwabHab<-pic(Habitat_SinSchwab, EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod)
 EIWk15_C4_SinPN_PIC_SchwabHab
 plot(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, main= "Habitat", cex=0.5)
   nodelabels (round(EIWk15_C4_SinPN_PIC_SchwabHab,5),adj=c(0,-0.5),frame="n")

 EIWk15_C4_SinPN_fit_PIC_SchwabHab<-lm(EIWk15_C4_SinPN_PIC_SchwabHab ~ EIWk15_C4_SinPN_PIC1+0)
 EIWk15_C4_SinPN_fit_PIC_SchwabHab
 summary(EIWk15_C4_SinPN_fit_PIC_SchwabHab)

 plot(EIWk15_C4_SinPN_PIC1, EIWk15_C4_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
 abline(EIWk15_C4_SinPN_fit_PIC_SchwabHab,lwd=2,lty="dashed",col="red")

### Wilberg et al. 2019 habitat groupings:
 Habitat_SinWilberg <-as.factor(Habitat_SinWilberg)
 EIWk15_C4_SinPN_PIC_WilbergHab<-pic(Habitat_SinWilberg, EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod)
 EIWk15_C4_SinPN_PIC_WilbergHab
 plot(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, main= "Habitat", cex=0.5)
   nodelabels (round(EIWk15_C4_SinPN_PIC_WilbergHab,5),adj=c(0,-0.5),frame="n")

 EIWk15_C4_SinPN_fit_PIC_WilbergHab<-lm(EIWk15_C4_SinPN_PIC_WilbergHab ~ EIWk15_C4_SinPN_PIC1 +0)
 EIWk15_C4_SinPN_fit_PIC_WilbergHab
 summary(EIWk15_C4_SinPN_fit_PIC_WilbergHab)

 plot(EIWk15_C4_SinPN_PIC1, EIWk15_C4_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
 abline(EIWk15_C4_SinPN_fit_PIC_WilbergHab,lwd=2,lty="dashed",col="red")

txtStop()


### Generates PDFs:
pdf("output/ParanasalAnalyses/EIWk15_C4/Figures/EIWk15_C4_paranasal_PIC_SchwabHab_regression.pdf", width=17, height=9) 
plot(EIWk15_C4_SinPN_PIC1, EIWk15_C4_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
abline(EIWk15_C4_SinPN_fit_PIC_SchwabHab,lwd=2,lty="dashed",col="red")
dev.off()

pdf("output/ParanasalAnalyses/EIWk15_C4/Figures/EIWk15_C4_paranasal_PIC_WilbergHab_regression.pdf", width=17, height=9) 
plot(EIWk15_C4_SinPN_PIC1, EIWk15_C4_SinPN_PIC2,bg="grey",cex=1.4,pch=21)
abline(EIWk15_C4_SinPN_fit_PIC_WilbergHab,lwd=2,lty="dashed",col="red")
dev.off()

 beep(1)


# =============== #
# 9. Blomberg's K #
# =============== #

### Required libraries:
 library(beepr)
 library(TeachingDemos)

# ----------------------------- #
# 9.1 Equal weights C0 analysis #
# ----------------------------- #

### Writes to file:
txtStart("output/ParanasalAnalyses/EW_C0/Files/EW_C0_paranasal_BlombergK.txt")

 EW_C0_Kvalue_PC1 <- phylosig(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, PCo1_SinPN, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EW_C0_Kvalue_PC1
 EW_C0_Kvalue_PC2 <- phylosig(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, PCo2_SinPN, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EW_C0_Kvalue_PC2
 EW_C0_Kvalue_PC3 <- phylosig(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, PCo3_SinPN, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EW_C0_Kvalue_PC3

txtStop()


# ----------------------------- #
# 9.2 Equal weights C1 analysis #
# ----------------------------- #

### Writes to file:
txtStart("output/ParanasalAnalyses/EW_C1/Files/EW_C1_paranasal_BlombergK.txt")

 EW_C1_Kvalue_PC1 <- phylosig(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, PCo1_SinPN, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EW_C1_Kvalue_PC1
 EW_C1_Kvalue_PC2 <- phylosig(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, PCo2_SinPN, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EW_C1_Kvalue_PC2
 EW_C1_Kvalue_PC3 <- phylosig(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, PCo3_SinPN, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EW_C1_Kvalue_PC3

txtStop()


# ----------------------------- #
# 9.3 Equal weights C3 analysis #
# ----------------------------- #

### Writes to file:
txtStart("output/ParanasalAnalyses/EW_C3/Files/EW_C3_paranasal_BlombergK.txt")

 EW_C3_Kvalue_PC1 <- phylosig(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, PCo1_SinPN, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EW_C3_Kvalue_PC1
 EW_C3_Kvalue_PC2 <- phylosig(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, PCo2_SinPN, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EW_C3_Kvalue_PC2
 EW_C3_Kvalue_PC3 <- phylosig(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, PCo3_SinPN, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EW_C3_Kvalue_PC3

txtStop()


# ----------------------------- #
# 9.4 Equal weights C4 analysis #
# ----------------------------- #

### Writes to file:
txtStart("output/ParanasalAnalyses/EW_C4/Files/EW_C4_paranasal_BlombergK_EQ.txt")

 EW_C4_Kvalue_PC1 <- phylosig(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, PCo1_SinPN, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EW_C4_Kvalue_PC1
 EW_C4_Kvalue_PC2 <- phylosig(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, PCo2_SinPN, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EW_C4_Kvalue_PC2
 EW_C4_Kvalue_PC3 <- phylosig(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, PCo3_SinPN, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EW_C4_Kvalue_PC3

txtStop()


# ------------------------------------------------- #
# 9.5 Extended implied weighting (k=15) C0 analysis #
# ------------------------------------------------- #

### Writes to file:
txtStart("output/ParanasalAnalyses/EIWk15_C0/Files/EIWk15_C0_paranasal_BlombergK.txt")

 EIWk15_C0_Kvalue_PC1 <- phylosig(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, PCo1_SinPN, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EIWk15_C0_Kvalue_PC1
 EIWk15_C0_Kvalue_PC2 <- phylosig(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, PCo2_SinPN, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EIWk15_C0_Kvalue_PC2
 EIWk15_C0_Kvalue_PC3 <- phylosig(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, PCo3_SinPN, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EIWk15_C0_Kvalue_PC3

txtStop()


# ------------------------------------------------- #
# 9.6 Extended implied weighting (k=15) C1 analysis #
# ------------------------------------------------- #

### Writes to file:
txtStart("output/ParanasalAnalyses/EIWk15_C1/Files/EIWk15_C1_paranasal_BlombergK.txt")

 EIWk15_C1_Kvalue_PC1 <- phylosig(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, PCo1_SinPN, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EIWk15_C1_Kvalue_PC1
 EIWk15_C1_Kvalue_PC2 <- phylosig(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, PCo2_SinPN, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EIWk15_C1_Kvalue_PC2
 EIWk15_C1_Kvalue_PC3 <- phylosig(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, PCo3_SinPN, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EIWk15_C1_Kvalue_PC3

txtStop()


# ------------------------------------------------- #
# 9.7 Extended implied weighting (k=15) C3 analysis #
# ------------------------------------------------- #

### Writes to file:
txtStart("output/ParanasalAnalyses/EIWk15_C3/Files/EIWk15_C3_paranasal_BlombergK.txt")

 EIWk15_C3_Kvalue_PC1 <- phylosig(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, PCo1_SinPN, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EIWk15_C3_Kvalue_PC1
 EIWk15_C3_Kvalue_PC2 <- phylosig(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, PCo2_SinPN, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EIWk15_C3_Kvalue_PC2
 EIWk15_C3_Kvalue_PC3 <- phylosig(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, PCo3_SinPN, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EIWk15_C3_Kvalue_PC3

txtStop()


# ------------------------------------------------- #
# 9.8 Extended implied weighting (k=15) C4 analysis #
# ------------------------------------------------- #

### Writes to file:
txtStart("output/ParanasalAnalyses/EIWk15_C4/Files/EIWk15_C4_paranasal_BlombergK.txt")

 EIWk15_C4_Kvalue_PC1 <- phylosig(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, PCo1_SinPN, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EIWk15_C4_Kvalue_PC1
 EIWk15_C4_Kvalue_PC2 <- phylosig(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, PCo2_SinPN, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EIWk15_C4_Kvalue_PC2
 EIWk15_C4_Kvalue_PC3 <- phylosig(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, PCo3_SinPN, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EIWk15_C4_Kvalue_PC3

txtStop()

 beep(1)


# ================== #
# 10. Pagel's Lambda #
# ================== #

### Required libraries:
 library(beepr)
 library(TeachingDemos)

# ------------------------------ #
# 10.1 Equal weights C0 analysis #
# ------------------------------ #

### Writes to file:
txtStart("output/ParanasalAnalyses/EW_C0/Files/EW_C0_paranasal_PagelsL.txt")

 EW_C0_lambda_PC1 <- phylosig(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, PCo1_SinPN, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EW_C0_lambda_PC1
 EW_C0_lambda_PC2 <- phylosig(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, PCo2_SinPN, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EW_C0_lambda_PC2
 EW_C0_lambda_PC3 <- phylosig(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, PCo3_SinPN, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EW_C0_lambda_PC3

txtStop()


# ------------------------------ #
# 10.2 Equal weights C1 analysis #
# ------------------------------ #

### Writes to file:
txtStart("output/ParanasalAnalyses/EW_C1/Files/EW_C1_paranasal_PagelsL.txt")

 EW_C1_lambda_PC1 <- phylosig(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, PCo1_SinPN, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EW_C1_lambda_PC1
 EW_C1_lambda_PC2 <- phylosig(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, PCo2_SinPN, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EW_C1_lambda_PC2
 EW_C1_lambda_PC3 <- phylosig(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, PCo3_SinPN, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EW_C1_lambda_PC3

txtStop()


# ------------------------------ #
# 10.3 Equal weights C3 analysis #
# ------------------------------ #

### Writes to file:
txtStart("output/ParanasalAnalyses/EW_C3/Files/EW_C3_paranasal_PagelsL.txt")

 EW_C3_lambda_PC1 <- phylosig(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, PCo1_SinPN, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EW_C3_lambda_PC1
 EW_C3_lambda_PC2 <- phylosig(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, PCo2_SinPN, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EW_C3_lambda_PC2
 EW_C3_lambda_PC3 <- phylosig(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, PCo3_SinPN, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EW_C3_lambda_PC3

txtStop()


# ------------------------------ #
# 10.4 Equal weights C4 analysis #
# ------------------------------ #

### Writes to file:
txtStart("output/ParanasalAnalyses/EW_C4/Files/EW_C4_paranasal_PagelsL.txt")

 EW_C4_lambda_PC1 <- phylosig(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, PCo1_SinPN, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EW_C4_lambda_PC1
 EW_C4_lambda_PC2 <- phylosig(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, PCo2_SinPN, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EW_C4_lambda_PC2
 EW_C4_lambda_PC3 <- phylosig(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, PCo3_SinPN, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EW_C4_lambda_PC3

txtStop()


# -------------------------------------------------- #
# 10.5 Extended implied weighting (k=15) C0 analysis #
# -------------------------------------------------- #

### Writes to file:
txtStart("output/ParanasalAnalyses/EIWk15_C0/Files/EIWk15_C0_paranasal_PagelsL.txt")

 EIWk15_C0_lambda_PC1 <- phylosig(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, PCo1_SinPN, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EIWk15_C0_lambda_PC1
 EIWk15_C0_lambda_PC2 <- phylosig(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, PCo2_SinPN, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EIWk15_C0_lambda_PC2
 EIWk15_C0_lambda_PC3 <- phylosig(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, PCo3_SinPN, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EIWk15_C0_lambda_PC3

txtStop()


# -------------------------------------------------- #
# 10.6 Extended implied weighting (k=15) C1 analysis #
# -------------------------------------------------- #

### Writes to file:
txtStart("output/ParanasalAnalyses/EIWk15_C1/Files/EIWk15_C1_paranasal_PagelsL.txt")

 EIWk15_C1_lambda_PC1 <- phylosig(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, PCo1_SinPN, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EIWk15_C1_lambda_PC1
 EIWk15_C1_lambda_PC2 <- phylosig(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, PCo2_SinPN, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EIWk15_C1_lambda_PC2
 EIWk15_C1_lambda_PC3 <- phylosig(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, PCo3_SinPN, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EIWk15_C1_lambda_PC3

txtStop()


# -------------------------------------------------- #
# 10.7 Extended implied weighting (k=15) C3 analysis #
# -------------------------------------------------- #

### Writes to file:
txtStart("output/ParanasalAnalyses/EIWk15_C3/Files/EIWk15_C3_paranasal_PagelsL.txt")

 EIWk15_C3_lambda_PC1 <- phylosig(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, PCo1_SinPN, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EIWk15_C3_lambda_PC1
 EIWk15_C3_lambda_PC2 <- phylosig(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, PCo2_SinPN, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EIWk15_C3_lambda_PC2
 EIWk15_C3_lambda_PC3 <- phylosig(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, PCo3_SinPN, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EIWk15_C3_lambda_PC3

txtStop()


# -------------------------------------------------- #
# 10.8 Extended implied weighting (k=15) C4 analysis #
# -------------------------------------------------- #

### Writes to file:
txtStart("output/ParanasalAnalyses/EIWk15_C4/Files/EIWk15_C4_paranasal_PagelsL.txt")

 EIWk15_C4_lambda_PC1 <- phylosig(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, PCo1_SinPN, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EIWk15_C4_lambda_PC1
 EIWk15_C4_lambda_PC2 <- phylosig(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, PCo2_SinPN, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EIWk15_C4_lambda_PC2
 EIWk15_C4_lambda_PC3 <- phylosig(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, PCo3_SinPN, method="lambda", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
 EIWk15_C4_lambda_PC3

txtStop()

 beep(1)


# ============= #
# 11. PERMANOVA #
# ============= #

### Required libraries:
 library(beepr)
 library(TeachingDemos)

### Writes to file:
txtStart("output/ParanasalAnalyses/NonPhylo/Files/Paranasal_PERMANOVA.txt")

### PCo scores against Schwab et al. 2020 habit groupings:
pairwise.adonis <- function(x, factors, sim.method = 'euclidean', p.adjust.m ='fdr')
{
  library(vegan)
  co = combn(unique(factors),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(co[1,elem],co[2,elem]),] ~ factors[factors %in% c(co[1,elem],co[2,elem])] , method =sim.method);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}


PERMANOVA.results <- pairwise.adonis(PCos_SinPN, factors=as.vector(Habitat_SinSchwab), sim.method = 'euclidean', p.adjust.m ='fdr')
PERMANOVA.results


### PCo scores against Wilberg et al. 2019 habitat groupings:
pairwise.adonis <- function(x, factors, sim.method = 'euclidean', p.adjust.m ='fdr')
{
  library(vegan)
  co = combn(unique(factors),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(co[1,elem],co[2,elem]),] ~ factors[factors %in% c(co[1,elem],co[2,elem])] , method =sim.method);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}


PERMANOVA.results <- pairwise.adonis(PCos_SinPN, factors=as.vector(Habitat_SinWilberg), sim.method = 'euclidean', p.adjust.m ='fdr')
PERMANOVA.results


### PCo scores against 'simple habitat' groupings:
pairwise.adonis <- function(x, factors, sim.method = 'euclidean', p.adjust.m ='fdr')
{
  library(vegan)
  co = combn(unique(factors),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(co[1,elem],co[2,elem]),] ~ factors[factors %in% c(co[1,elem],co[2,elem])] , method =sim.method);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}


PERMANOVA.results <- pairwise.adonis(PCos_SinPN, factors=as.vector(Habitat_SinSimple), sim.method = 'euclidean', p.adjust.m ='fdr')
PERMANOVA.results


### PCo scores against Thalattosuchia clade grouping:
pairwise.adonis <- function(x, factors, sim.method = 'euclidean', p.adjust.m ='fdr')
{
  library(vegan)
  co = combn(unique(factors),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(co[1,elem],co[2,elem]),] ~ factors[factors %in% c(co[1,elem],co[2,elem])] , method =sim.method);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}


PERMANOVA.results <- pairwise.adonis(PCos_SinPN, factors=as.vector(Thalattosuchian_Sin), sim.method = 'euclidean', p.adjust.m ='fdr')
PERMANOVA.results


### PCo scores against Metriorhynchidae clade grouping:
pairwise.adonis <- function(x, factors, sim.method = 'euclidean', p.adjust.m ='fdr')
{
  library(vegan)
  co = combn(unique(factors),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(co[1,elem],co[2,elem]),] ~ factors[factors %in% c(co[1,elem],co[2,elem])] , method =sim.method);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}


PERMANOVA.results <- pairwise.adonis(PCos_SinPN, factors=as.vector(Metriorhynchid_Sin), sim.method = 'euclidean', p.adjust.m ='fdr')
PERMANOVA.results

txtStop()

 beep(1)


# ====================================== #
# 12. Phylogenetic General Linear Models #
# ====================================== #

### Required libraries:
 library(ape)
 library(beepr)
 library(nlme)
 library(rr2)


# ---------------------------- #
# 12.1 GLMS on the PCoA scores #
# ---------------------------- #

### Writes to file:
txtStart("output/ParanasalAnalyses/NonPhylo/Files/Paranasal_OLS.txt")

### PCo1 OLS
# PCo1 v terrestrial:
 PCo1_SinPN_OLS<-gls(PCo1_SinPN~Terrestrial_Sin, method= "ML", data= CrocSinPNData)
  summary(PCo1_SinPN_OLS)
  plot(PCo1_SinPN_OLS)
  AIC(PCo1_SinPN_OLS)
  R2_lik(PCo1_SinPN_OLS)

# PCo1 v semiaquatic:
 PCo1_SinPN_OLS<-gls(PCo1_SinPN~Semiaquatic_Sin, method= "ML", data= CrocSinPNData)
  summary(PCo1_SinPN_OLS)
  plot(PCo1_SinPN_OLS)
  AIC(PCo1_SinPN_OLS)
  R2_lik(PCo1_SinPN_OLS)

# PCo1 v pelagic:
 PCo1_SinPN_OLS<-gls(PCo1_SinPN~Pelagic_Sin, method= "ML", data= CrocSinPNData)
  summary(PCo1_SinPN_OLS)
  plot(PCo1_SinPN_OLS)
  AIC(PCo1_SinPN_OLS)
  R2_lik(PCo1_SinPN_OLS)

# PCo1 v marine:
 PCo1_SinPN_OLS<-gls(PCo1_SinPN~Marine_Sin, method= "ML", data= CrocSinPNData)
  summary(PCo1_SinPN_OLS)
  plot(PCo1_SinPN_OLS)
  AIC(PCo1_SinPN_OLS)
  R2_lik(PCo1_SinPN_OLS)

# PCo1 v freshwater:
 PCo1_SinPN_OLS<-gls(PCo1_SinPN~Freshwater_Sin, method= "ML", data= CrocSinPNData)
  summary(PCo1_SinPN_OLS)
  plot(PCo1_SinPN_OLS)
  AIC(PCo1_SinPN_OLS)
  R2_lik(PCo1_SinPN_OLS)

# PCo1 v aquatic:
 PCo1_SinPN_OLS<-gls(PCo1_SinPN~Aquatic_Sin, method= "ML", data= CrocSinPNData)
  summary(PCo1_SinPN_OLS)
  plot(PCo1_SinPN_OLS)
  AIC(PCo1_SinPN_OLS)
  R2_lik(PCo1_SinPN_OLS)

# PCo1 v Thalattosuchia:
 PCo1_SinPN_OLS<-gls(PCo1_SinPN~Thalattosuchian_Sin, method= "ML", data= CrocSinPNData)
  summary(PCo1_SinPN_OLS)
  plot(PCo1_SinPN_OLS)
  AIC(PCo1_SinPN_OLS)
  R2_lik(PCo1_SinPN_OLS)

# PCo1 v Metriorhynchidae:
 PCo1_SinPN_OLS<-gls(PCo1_SinPN~Metriorhynchid_Sin, method= "ML", data= CrocSinPNData)
  summary(PCo1_SinPN_OLS)
  plot(PCo1_SinPN_OLS)
  AIC(PCo1_SinPN_OLS)
  R2_lik(PCo1_SinPN_OLS)


### PCo2 OLS
# PCo2 v terrestrial:
 PCo2_SinPN_OLS<-gls(PCo2_SinPN~Terrestrial_Sin, method= "ML", data= CrocSinPNData)
  summary(PCo2_SinPN_OLS)
  plot(PCo2_SinPN_OLS)
  AIC(PCo2_SinPN_OLS)
  R2_lik(PCo2_SinPN_OLS)

# PCo2 v semiaquatic:
 PCo2_SinPN_OLS<-gls(PCo2_SinPN~Semiaquatic_Sin, method= "ML", data= CrocSinPNData)
  summary(PCo2_SinPN_OLS)
  plot(PCo2_SinPN_OLS)
  AIC(PCo2_SinPN_OLS)
  R2_lik(PCo2_SinPN_OLS)

# PCo2 v pelagic:
 PCo2_SinPN_OLS<-gls(PCo2_SinPN~Pelagic_Sin, method= "ML", data= CrocSinPNData)
  summary(PCo2_SinPN_OLS)
  plot(PCo2_SinPN_OLS)
  AIC(PCo2_SinPN_OLS)
  R2_lik(PCo2_SinPN_OLS)

# PCo2 v marine:
 PCo2_SinPN_OLS<-gls(PCo2_SinPN~Marine_Sin, method= "ML", data= CrocSinPNData)
  summary(PCo2_SinPN_OLS)
  plot(PCo2_SinPN_OLS)
  AIC(PCo2_SinPN_OLS)
  R2_lik(PCo2_SinPN_OLS)

# PCo2 v freshwater:
 PCo2_SinPN_OLS<-gls(PCo2_SinPN~Freshwater_Sin, method= "ML", data= CrocSinPNData)
  summary(PCo2_SinPN_OLS)
  plot(PCo2_SinPN_OLS)
  AIC(PCo2_SinPN_OLS)
  R2_lik(PCo2_SinPN_OLS)

# PCo2 v aquatic:
 PCo2_SinPN_OLS<-gls(PCo2_SinPN~Aquatic_Sin, method= "ML", data= CrocSinPNData)
  summary(PCo2_SinPN_OLS)
  plot(PCo2_SinPN_OLS)
  AIC(PCo2_SinPN_OLS)
  R2_lik(PCo2_SinPN_OLS)

# PCo2 v Thalattosuchia:
 PCo2_SinPN_OLS<-gls(PCo2_SinPN~Thalattosuchian_Sin, method= "ML", data= CrocSinPNData)
  summary(PCo2_SinPN_OLS)
  plot(PCo2_SinPN_OLS)
  AIC(PCo2_SinPN_OLS)
  R2_lik(PCo2_SinPN_OLS)

# PCo2 v Metriorhynchidae:
 PCo2_SinPN_OLS<-gls(PCo2_SinPN~Metriorhynchid_Sin, method= "ML", data= CrocSinPNData)
  summary(PCo2_SinPN_OLS)
  plot(PCo2_SinPN_OLS)
  AIC(PCo2_SinPN_OLS)
  R2_lik(PCo2_SinPN_OLS)


### PCo3 OLS
# PCo3 v terrestrial:
 PCo3_SinPN_OLS<-gls(PCo3_SinPN~Terrestrial_Sin, method= "ML", data= CrocSinPNData)
  summary(PCo3_SinPN_OLS)
  plot(PCo3_SinPN_OLS)
  AIC(PCo3_SinPN_OLS)
  R2_lik(PCo3_SinPN_OLS)

# PCo3 v semiaquatic:
 PCo3_SinPN_OLS<-gls(PCo3_SinPN~Semiaquatic_Sin, method= "ML", data= CrocSinPNData)
  summary(PCo3_SinPN_OLS)
  plot(PCo3_SinPN_OLS)
  AIC(PCo3_SinPN_OLS)
  R2_lik(PCo3_SinPN_OLS)

# PCo3 v pelagic:
 PCo3_SinPN_OLS<-gls(PCo3_SinPN~Pelagic_Sin, method= "ML", data= CrocSinPNData)
  summary(PCo3_SinPN_OLS)
  plot(PCo3_SinPN_OLS)
  AIC(PCo3_SinPN_OLS)
  R2_lik(PCo3_SinPN_OLS)

# PCo3 v marine:
 PCo3_SinPN_OLS<-gls(PCo3_SinPN~Marine_Sin, method= "ML", data= CrocSinPNData)
  summary(PCo3_SinPN_OLS)
  plot(PCo3_SinPN_OLS)
  AIC(PCo3_SinPN_OLS)
  R2_lik(PCo3_SinPN_OLS)

# PCo3 v freshwater:
 PCo3_SinPN_OLS<-gls(PCo3_SinPN~Freshwater_Sin, method= "ML", data= CrocSinPNData)
  summary(PCo3_SinPN_OLS)
  plot(PCo3_SinPN_OLS)
  AIC(PCo3_SinPN_OLS)
  R2_lik(PCo3_SinPN_OLS)

# PCo3 v aquatic:
 PCo3_SinPN_OLS<-gls(PCo3_SinPN~Aquatic_Sin, method= "ML", data= CrocSinPNData)
  summary(PCo3_SinPN_OLS)
  plot(PCo3_SinPN_OLS)
  AIC(PCo3_SinPN_OLS)
  R2_lik(PCo3_SinPN_OLS)

# PCo3 v Thalattosuchia:
 PCo3_SinPN_OLS<-gls(PCo3_SinPN~Thalattosuchian_Sin, method= "ML", data= CrocSinPNData)
  summary(PCo3_SinPN_OLS)
  plot(PCo3_SinPN_OLS)
  AIC(PCo3_SinPN_OLS)
  R2_lik(PCo3_SinPN_OLS)

# PCo3 v Metriorhynchidae:
 PCo3_SinPN_OLS<-gls(PCo3_SinPN~Metriorhynchid_Sin, method= "ML", data= CrocSinPNData)
  summary(PCo3_SinPN_OLS)
  plot(PCo3_SinPN_OLS)
  AIC(PCo3_SinPN_OLS)
  R2_lik(PCo3_SinPN_OLS)

txtStop()


# ------------------------------------ #
# 12.2 Equal weights C0 paranasal GLMs #
# ------------------------------------ #

### Writes to file:
txtStart("output/ParanasalAnalyses/EW_C0/Files/EW_C0_paranasal_pGLS.txt")


EW_C0_SinPN_bm<-corBrownian(1,EW_C0_phylogeny_SinPN_Pruned_EqualMethod, form= ~Species_Sin)


### PCo1 pGLS
# PCo1 v terrestrial:
 EW_C0_PCo1_SinPN_pGLS_Terrestrial<-gls(PCo1_SinPN~ Terrestrial_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C0_SinPN_bm)
  summary(EW_C0_PCo1_SinPN_pGLS_Terrestrial)
  plot(EW_C0_PCo1_SinPN_pGLS_Terrestrial)
  AIC(EW_C0_PCo1_SinPN_pGLS_Terrestrial)
  R2_lik(EW_C0_PCo1_SinPN_pGLS_Terrestrial)

# PCo1 v semiaquatic:
 EW_C0_PCo1_SinPN_pGLS_semiaquatic<-gls(PCo1_SinPN~ Semiaquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C0_SinPN_bm)
  summary(EW_C0_PCo1_SinPN_pGLS_semiaquatic)
  plot(EW_C0_PCo1_SinPN_pGLS_semiaquatic)
  AIC(EW_C0_PCo1_SinPN_pGLS_semiaquatic)
  R2_lik(EW_C0_PCo1_SinPN_pGLS_semiaquatic)

# PCo1 v pelagic:
 EW_C0_PCo1_SinPN_pGLS_pelagic<-gls(PCo1_SinPN~ Pelagic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C0_SinPN_bm)
  summary(EW_C0_PCo1_SinPN_pGLS_pelagic)
  plot(EW_C0_PCo1_SinPN_pGLS_pelagic)
  AIC(EW_C0_PCo1_SinPN_pGLS_pelagic)
  R2_lik(EW_C0_PCo1_SinPN_pGLS_pelagic)

# PCo1 v marine:
 EW_C0_PCo1_SinPN_pGLS_marine<-gls(PCo1_SinPN~ Marine_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C0_SinPN_bm)
  summary(EW_C0_PCo1_SinPN_pGLS_marine)
  plot(EW_C0_PCo1_SinPN_pGLS_marine)
  AIC(EW_C0_PCo1_SinPN_pGLS_marine)
  R2_lik(EW_C0_PCo1_SinPN_pGLS_marine)

# PCo1 v freshwater:
 EW_C0_PCo1_SinPN_pGLS_freshwater<-gls(PCo1_SinPN~ Freshwater_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C0_SinPN_bm)
  summary(EW_C0_PCo1_SinPN_pGLS_freshwater)
  plot(EW_C0_PCo1_SinPN_pGLS_freshwater)
  AIC(EW_C0_PCo1_SinPN_pGLS_freshwater)
  R2_lik(EW_C0_PCo1_SinPN_pGLS_freshwater)

# PCo1 v aquatic:
 EW_C0_PCo1_SinPN_pGLS_aquatic<-gls(PCo1_SinPN~ Aquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C0_SinPN_bm)
  summary(EW_C0_PCo1_SinPN_pGLS_aquatic)
  plot(EW_C0_PCo1_SinPN_pGLS_aquatic)
  AIC(EW_C0_PCo1_SinPN_pGLS_aquatic)
  R2_lik(EW_C0_PCo1_SinPN_pGLS_aquatic)

# PCo1 v Thalattosuchia:
 EW_C0_PCo1_SinPN_pGLS_thalattosuchia<-gls(PCo1_SinPN~ Thalattosuchian_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C0_SinPN_bm)
  summary(EW_C0_PCo1_SinPN_pGLS_thalattosuchia)
  plot(EW_C0_PCo1_SinPN_pGLS_thalattosuchia)
  AIC(EW_C0_PCo1_SinPN_pGLS_thalattosuchia)
  R2_lik(EW_C0_PCo1_SinPN_pGLS_thalattosuchia)

# PCo1 v Metriorhynchidae:
 EW_C0_PCo1_SinPN_pGLS_metriorhynchidae<-gls(PCo1_SinPN~ Metriorhynchid_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C0_SinPN_bm)
  summary(EW_C0_PCo1_SinPN_pGLS_metriorhynchidae)
  plot(EW_C0_PCo1_SinPN_pGLS_metriorhynchidae)
  AIC(EW_C0_PCo1_SinPN_pGLS_metriorhynchidae)
  R2_lik(EW_C0_PCo1_SinPN_pGLS_metriorhynchidae)


### PCo2 pGLS
# PCo2 v terrestrial:
 EW_C0_PCo2_SinPN_pGLS_terrestrial<-gls(PCo2_SinPN~ Terrestrial_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C0_SinPN_bm)
  summary(EW_C0_PCo2_SinPN_pGLS_terrestrial)
  plot(EW_C0_PCo2_SinPN_pGLS_terrestrial)
  AIC(EW_C0_PCo2_SinPN_pGLS_terrestrial)
  R2_lik(EW_C0_PCo2_SinPN_pGLS_terrestrial) 

# PCo2 v semiaquatic:
 EW_C0_PCo2_SinPN_pGLS_semiaquatic<-gls(PCo2_SinPN~ Semiaquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C0_SinPN_bm)
  summary(EW_C0_PCo2_SinPN_pGLS_semiaquatic)
  plot(EW_C0_PCo2_SinPN_pGLS_semiaquatic)
  AIC(EW_C0_PCo2_SinPN_pGLS_semiaquatic)
  R2_lik(EW_C0_PCo2_SinPN_pGLS_semiaquatic) 

# PCo2 v pelagic:
 EW_C0_PCo2_SinPN_pGLS_pelagic<-gls(PCo2_SinPN~ Pelagic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C0_SinPN_bm)
  summary(EW_C0_PCo2_SinPN_pGLS_pelagic)
  plot(EW_C0_PCo2_SinPN_pGLS_pelagic)
  AIC(EW_C0_PCo2_SinPN_pGLS_pelagic)
  R2_lik(EW_C0_PCo2_SinPN_pGLS_pelagic) 

# PCo2 v marine:
 EW_C0_PCo2_SinPN_pGLS_marine<-gls(PCo2_SinPN~ Marine_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C0_SinPN_bm)
  summary(EW_C0_PCo2_SinPN_pGLS_marine)
  plot(EW_C0_PCo2_SinPN_pGLS_marine)
  AIC(EW_C0_PCo2_SinPN_pGLS_marine)
  R2_lik(EW_C0_PCo2_SinPN_pGLS_marine) 

# PCo2 v freshwater:
 EW_C0_PCo2_SinPN_pGLS_freshwater<-gls(PCo2_SinPN~ Freshwater_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C0_SinPN_bm)
  summary(EW_C0_PCo2_SinPN_pGLS_freshwater)
  plot(EW_C0_PCo2_SinPN_pGLS_freshwater)
  AIC(EW_C0_PCo2_SinPN_pGLS_freshwater)
  R2_lik(EW_C0_PCo2_SinPN_pGLS_freshwater) 

# PCo2 v aquatic:
 EW_C0_PCo2_SinPN_pGLS_aquatic<-gls(PCo2_SinPN~ Aquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C0_SinPN_bm)
  summary(EW_C0_PCo2_SinPN_pGLS_aquatic)
  plot(EW_C0_PCo2_SinPN_pGLS_aquatic)
  AIC(EW_C0_PCo2_SinPN_pGLS_aquatic)
  R2_lik(EW_C0_PCo2_SinPN_pGLS_aquatic) 

# PCo2 v Thalattosuchia:
 EW_C0_PCo2_SinPN_pGLS_thalattosuchia<-gls(PCo2_SinPN~ Thalattosuchian_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C0_SinPN_bm)
  summary(EW_C0_PCo2_SinPN_pGLS_thalattosuchia)
  plot(EW_C0_PCo2_SinPN_pGLS_thalattosuchia)
  AIC(EW_C0_PCo2_SinPN_pGLS_thalattosuchia)
  R2_lik(EW_C0_PCo2_SinPN_pGLS_thalattosuchia) 

# PCo2 v Metriorhynchidae:
 EW_C0_PCo2_SinPN_pGLS_metriorhynchidae<-gls(PCo2_SinPN~ Metriorhynchid_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C0_SinPN_bm)
  summary(EW_C0_PCo2_SinPN_pGLS_metriorhynchidae)
  plot(EW_C0_PCo2_SinPN_pGLS_metriorhynchidae)
  AIC(EW_C0_PCo2_SinPN_pGLS_metriorhynchidae)
  R2_lik(EW_C0_PCo2_SinPN_pGLS_metriorhynchidae) 


### PCo3 pGLS
# PCo3 v terrestrial:
 EW_C0_PCo3_SinPN_pGLS_terrestrial<-gls(PCo3_SinPN~ Terrestrial_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C0_SinPN_bm)
  summary(EW_C0_PCo3_SinPN_pGLS_terrestrial)
  plot(EW_C0_PCo3_SinPN_pGLS_terrestrial)
  AIC(EW_C0_PCo3_SinPN_pGLS_terrestrial)
  R2_lik(EW_C0_PCo3_SinPN_pGLS_terrestrial)

# PCo3 v semiaquatic:
 EW_C0_PCo3_SinPN_pGLS_semiaquatic<-gls(PCo3_SinPN~ Semiaquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C0_SinPN_bm)
  summary(EW_C0_PCo3_SinPN_pGLS_semiaquatic)
  plot(EW_C0_PCo3_SinPN_pGLS_semiaquatic)
  AIC(EW_C0_PCo3_SinPN_pGLS_semiaquatic)
  R2_lik(EW_C0_PCo3_SinPN_pGLS_semiaquatic)

# PCo3 v pelagic:
 EW_C0_PCo3_SinPN_pGLS_pelagic<-gls(PCo3_SinPN~ Pelagic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C0_SinPN_bm)
  summary(EW_C0_PCo3_SinPN_pGLS_pelagic)
  plot(EW_C0_PCo3_SinPN_pGLS_pelagic)
  AIC(EW_C0_PCo3_SinPN_pGLS_pelagic)
  R2_lik(EW_C0_PCo3_SinPN_pGLS_pelagic)

# PCo3 v marine:
 EW_C0_PCo3_SinPN_pGLS_marine<-gls(PCo3_SinPN~ Marine_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C0_SinPN_bm)
  summary(EW_C0_PCo3_SinPN_pGLS_marine)
  plot(EW_C0_PCo3_SinPN_pGLS_marine)
  AIC(EW_C0_PCo3_SinPN_pGLS_marine)
  R2_lik(EW_C0_PCo3_SinPN_pGLS_marine)

# PCo3 v freshwater:
 EW_C0_PCo3_SinPN_pGLS_freshwater<-gls(PCo3_SinPN~ Freshwater_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C0_SinPN_bm)
  summary(EW_C0_PCo3_SinPN_pGLS_freshwater)
  plot(EW_C0_PCo3_SinPN_pGLS_freshwater)
  AIC(EW_C0_PCo3_SinPN_pGLS_freshwater)
  R2_lik(EW_C0_PCo3_SinPN_pGLS_freshwater)

# PCo3 v aquatic:
 EW_C0_PCo3_SinPN_pGLS_aquatic<-gls(PCo3_SinPN~ Aquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C0_SinPN_bm)
  summary(EW_C0_PCo3_SinPN_pGLS_aquatic)
  plot(EW_C0_PCo3_SinPN_pGLS_aquatic)
  AIC(EW_C0_PCo3_SinPN_pGLS_aquatic)
  R2_lik(EW_C0_PCo3_SinPN_pGLS_aquatic)

# PCo3 v Thalattosuchia:
 EW_C0_PCo3_SinPN_pGLS_thalattosuchia<-gls(PCo3_SinPN~ Thalattosuchian_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C0_SinPN_bm)
  summary(EW_C0_PCo3_SinPN_pGLS_thalattosuchia)
  plot(EW_C0_PCo3_SinPN_pGLS_thalattosuchia)
  AIC(EW_C0_PCo3_SinPN_pGLS_thalattosuchia)
  R2_lik(EW_C0_PCo3_SinPN_pGLS_thalattosuchia)

# PCo3 v Metriorhynchidae:
 EW_C0_PCo3_SinPN_pGLS_metriorhynchidae<-gls(PCo3_SinPN~ Metriorhynchid_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C0_SinPN_bm)
  summary(EW_C0_PCo3_SinPN_pGLS_metriorhynchidae)
  plot(EW_C0_PCo3_SinPN_pGLS_metriorhynchidae)
  AIC(EW_C0_PCo3_SinPN_pGLS_metriorhynchidae)
  R2_lik(EW_C0_PCo3_SinPN_pGLS_metriorhynchidae)

txtStop()


# ------------------------------------ #
# 12.3 Equal weights C1 paranasal GLMs #
# ------------------------------------ #

### Writes to file:
txtStart("output/ParanasalAnalyses/EW_C1/Files/EW_C1_paranasal_pGLS.txt")


EW_C1_SinPN_bm<-corBrownian(1,EW_C1_phylogeny_SinPN_Pruned_EqualMethod, form= ~Species_Sin)


### PCo1 pGLS
# PCo1 v terrestrial:
 EW_C1_PCo1_SinPN_pGLS_Terrestrial<-gls(PCo1_SinPN~ Terrestrial_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C1_SinPN_bm)
  summary(EW_C1_PCo1_SinPN_pGLS_Terrestrial)
  plot(EW_C1_PCo1_SinPN_pGLS_Terrestrial)
  AIC(EW_C1_PCo1_SinPN_pGLS_Terrestrial)
  R2_lik(EW_C1_PCo1_SinPN_pGLS_Terrestrial)

# PCo1 v semiaquatic:
 EW_C1_PCo1_SinPN_pGLS_semiaquatic<-gls(PCo1_SinPN~ Semiaquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C1_SinPN_bm)
  summary(EW_C1_PCo1_SinPN_pGLS_semiaquatic)
  plot(EW_C1_PCo1_SinPN_pGLS_semiaquatic)
  AIC(EW_C1_PCo1_SinPN_pGLS_semiaquatic)
  R2_lik(EW_C1_PCo1_SinPN_pGLS_semiaquatic)

# PCo1 v pelagic:
 EW_C1_PCo1_SinPN_pGLS_pelagic<-gls(PCo1_SinPN~ Pelagic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C1_SinPN_bm)
  summary(EW_C1_PCo1_SinPN_pGLS_pelagic)
  plot(EW_C1_PCo1_SinPN_pGLS_pelagic)
  AIC(EW_C1_PCo1_SinPN_pGLS_pelagic)
  R2_lik(EW_C1_PCo1_SinPN_pGLS_pelagic)

# PCo1 v marine:
 EW_C1_PCo1_SinPN_pGLS_marine<-gls(PCo1_SinPN~ Marine_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C1_SinPN_bm)
  summary(EW_C1_PCo1_SinPN_pGLS_marine)
  plot(EW_C1_PCo1_SinPN_pGLS_marine)
  AIC(EW_C1_PCo1_SinPN_pGLS_marine)
  R2_lik(EW_C1_PCo1_SinPN_pGLS_marine)

# PCo1 v freshwater:
 EW_C1_PCo1_SinPN_pGLS_freshwater<-gls(PCo1_SinPN~ Freshwater_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C1_SinPN_bm)
  summary(EW_C1_PCo1_SinPN_pGLS_freshwater)
  plot(EW_C1_PCo1_SinPN_pGLS_freshwater)
  AIC(EW_C1_PCo1_SinPN_pGLS_freshwater)
  R2_lik(EW_C1_PCo1_SinPN_pGLS_freshwater)

# PCo1 v aquatic:
 EW_C1_PCo1_SinPN_pGLS_aquatic<-gls(PCo1_SinPN~ Aquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C1_SinPN_bm)
  summary(EW_C1_PCo1_SinPN_pGLS_aquatic)
  plot(EW_C1_PCo1_SinPN_pGLS_aquatic)
  AIC(EW_C1_PCo1_SinPN_pGLS_aquatic)
  R2_lik(EW_C1_PCo1_SinPN_pGLS_aquatic)

# PCo1 v Thalattosuchia:
 EW_C1_PCo1_SinPN_pGLS_thalattosuchia<-gls(PCo1_SinPN~ Thalattosuchian_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C1_SinPN_bm)
  summary(EW_C1_PCo1_SinPN_pGLS_thalattosuchia)
  plot(EW_C1_PCo1_SinPN_pGLS_thalattosuchia)
  AIC(EW_C1_PCo1_SinPN_pGLS_thalattosuchia)
  R2_lik(EW_C1_PCo1_SinPN_pGLS_thalattosuchia)

# PCo1 v Metriorhynchidae:
 EW_C1_PCo1_SinPN_pGLS_metriorhynchidae<-gls(PCo1_SinPN~ Metriorhynchid_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C1_SinPN_bm)
  summary(EW_C1_PCo1_SinPN_pGLS_metriorhynchidae)
  plot(EW_C1_PCo1_SinPN_pGLS_metriorhynchidae)
  AIC(EW_C1_PCo1_SinPN_pGLS_metriorhynchidae)
  R2_lik(EW_C1_PCo1_SinPN_pGLS_metriorhynchidae)


### PCo2 pGLS
# PCo2 v terrestrial:
 EW_C1_PCo2_SinPN_pGLS_terrestrial<-gls(PCo2_SinPN~ Terrestrial_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C1_SinPN_bm)
  summary(EW_C1_PCo2_SinPN_pGLS_terrestrial)
  plot(EW_C1_PCo2_SinPN_pGLS_terrestrial)
  AIC(EW_C1_PCo2_SinPN_pGLS_terrestrial)
  R2_lik(EW_C1_PCo2_SinPN_pGLS_terrestrial) 

# PCo2 v semiaquatic:
 EW_C1_PCo2_SinPN_pGLS_semiaquatic<-gls(PCo2_SinPN~ Semiaquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C1_SinPN_bm)
  summary(EW_C1_PCo2_SinPN_pGLS_semiaquatic)
  plot(EW_C1_PCo2_SinPN_pGLS_semiaquatic)
  AIC(EW_C1_PCo2_SinPN_pGLS_semiaquatic)
  R2_lik(EW_C1_PCo2_SinPN_pGLS_semiaquatic) 

# PCo2 v pelagic:
 EW_C1_PCo2_SinPN_pGLS_pelagic<-gls(PCo2_SinPN~ Pelagic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C1_SinPN_bm)
  summary(EW_C1_PCo2_SinPN_pGLS_pelagic)
  plot(EW_C1_PCo2_SinPN_pGLS_pelagic)
  AIC(EW_C1_PCo2_SinPN_pGLS_pelagic)
  R2_lik(EW_C1_PCo2_SinPN_pGLS_pelagic) 

# PCo2 v marine:
 EW_C1_PCo2_SinPN_pGLS_marine<-gls(PCo2_SinPN~ Marine_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C1_SinPN_bm)
  summary(EW_C1_PCo2_SinPN_pGLS_marine)
  plot(EW_C1_PCo2_SinPN_pGLS_marine)
  AIC(EW_C1_PCo2_SinPN_pGLS_marine)
  R2_lik(EW_C1_PCo2_SinPN_pGLS_marine) 

# PCo2 v freshwater:
 EW_C1_PCo2_SinPN_pGLS_freshwater<-gls(PCo2_SinPN~ Freshwater_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C1_SinPN_bm)
  summary(EW_C1_PCo2_SinPN_pGLS_freshwater)
  plot(EW_C1_PCo2_SinPN_pGLS_freshwater)
  AIC(EW_C1_PCo2_SinPN_pGLS_freshwater)
  R2_lik(EW_C1_PCo2_SinPN_pGLS_freshwater) 

# PCo2 v aquatic:
 EW_C1_PCo2_SinPN_pGLS_aquatic<-gls(PCo2_SinPN~ Aquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C1_SinPN_bm)
  summary(EW_C1_PCo2_SinPN_pGLS_aquatic)
  plot(EW_C1_PCo2_SinPN_pGLS_aquatic)
  AIC(EW_C1_PCo2_SinPN_pGLS_aquatic)
  R2_lik(EW_C1_PCo2_SinPN_pGLS_aquatic) 

# PCo2 v Thalattosuchia:
 EW_C1_PCo2_SinPN_pGLS_thalattosuchia<-gls(PCo2_SinPN~ Thalattosuchian_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C1_SinPN_bm)
  summary(EW_C1_PCo2_SinPN_pGLS_thalattosuchia)
  plot(EW_C1_PCo2_SinPN_pGLS_thalattosuchia)
  AIC(EW_C1_PCo2_SinPN_pGLS_thalattosuchia)
  R2_lik(EW_C1_PCo2_SinPN_pGLS_thalattosuchia) 

# PCo2 v Metriorhynchidae:
 EW_C1_PCo2_SinPN_pGLS_metriorhynchidae<-gls(PCo2_SinPN~ Metriorhynchid_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C1_SinPN_bm)
  summary(EW_C1_PCo2_SinPN_pGLS_metriorhynchidae)
  plot(EW_C1_PCo2_SinPN_pGLS_metriorhynchidae)
  AIC(EW_C1_PCo2_SinPN_pGLS_metriorhynchidae)
  R2_lik(EW_C1_PCo2_SinPN_pGLS_metriorhynchidae) 


### PCo3 pGLS
# PCo3 v terrestrial:
 EW_C1_PCo3_SinPN_pGLS_terrestrial<-gls(PCo3_SinPN~ Terrestrial_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C1_SinPN_bm)
  summary(EW_C1_PCo3_SinPN_pGLS_terrestrial)
  plot(EW_C1_PCo3_SinPN_pGLS_terrestrial)
  AIC(EW_C1_PCo3_SinPN_pGLS_terrestrial)
  R2_lik(EW_C1_PCo3_SinPN_pGLS_terrestrial)

# PCo3 v semiaquatic:
 EW_C1_PCo3_SinPN_pGLS_semiaquatic<-gls(PCo3_SinPN~ Semiaquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C1_SinPN_bm)
  summary(EW_C1_PCo3_SinPN_pGLS_semiaquatic)
  plot(EW_C1_PCo3_SinPN_pGLS_semiaquatic)
  AIC(EW_C1_PCo3_SinPN_pGLS_semiaquatic)
  R2_lik(EW_C1_PCo3_SinPN_pGLS_semiaquatic)

# PCo3 v pelagic:
 EW_C1_PCo3_SinPN_pGLS_pelagic<-gls(PCo3_SinPN~ Pelagic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C1_SinPN_bm)
  summary(EW_C1_PCo3_SinPN_pGLS_pelagic)
  plot(EW_C1_PCo3_SinPN_pGLS_pelagic)
  AIC(EW_C1_PCo3_SinPN_pGLS_pelagic)
  R2_lik(EW_C1_PCo3_SinPN_pGLS_pelagic)

# PCo3 v marine:
 EW_C1_PCo3_SinPN_pGLS_marine<-gls(PCo3_SinPN~ Marine_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C1_SinPN_bm)
  summary(EW_C1_PCo3_SinPN_pGLS_marine)
  plot(EW_C1_PCo3_SinPN_pGLS_marine)
  AIC(EW_C1_PCo3_SinPN_pGLS_marine)
  R2_lik(EW_C1_PCo3_SinPN_pGLS_marine)

# PCo3 v freshwater:
 EW_C1_PCo3_SinPN_pGLS_freshwater<-gls(PCo3_SinPN~ Freshwater_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C1_SinPN_bm)
  summary(EW_C1_PCo3_SinPN_pGLS_freshwater)
  plot(EW_C1_PCo3_SinPN_pGLS_freshwater)
  AIC(EW_C1_PCo3_SinPN_pGLS_freshwater)
  R2_lik(EW_C1_PCo3_SinPN_pGLS_freshwater)

# PCo3 v aquatic:
 EW_C1_PCo3_SinPN_pGLS_aquatic<-gls(PCo3_SinPN~ Aquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C1_SinPN_bm)
  summary(EW_C1_PCo3_SinPN_pGLS_aquatic)
  plot(EW_C1_PCo3_SinPN_pGLS_aquatic)
  AIC(EW_C1_PCo3_SinPN_pGLS_aquatic)
  R2_lik(EW_C1_PCo3_SinPN_pGLS_aquatic)

# PCo3 v Thalattosuchia:
 EW_C1_PCo3_SinPN_pGLS_thalattosuchia<-gls(PCo3_SinPN~ Thalattosuchian_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C1_SinPN_bm)
  summary(EW_C1_PCo3_SinPN_pGLS_thalattosuchia)
  plot(EW_C1_PCo3_SinPN_pGLS_thalattosuchia)
  AIC(EW_C1_PCo3_SinPN_pGLS_thalattosuchia)
  R2_lik(EW_C1_PCo3_SinPN_pGLS_thalattosuchia)

# PCo3 v Metriorhynchidae:
 EW_C1_PCo3_SinPN_pGLS_metriorhynchidae<-gls(PCo3_SinPN~ Metriorhynchid_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C1_SinPN_bm)
  summary(EW_C1_PCo3_SinPN_pGLS_metriorhynchidae)
  plot(EW_C1_PCo3_SinPN_pGLS_metriorhynchidae)
  AIC(EW_C1_PCo3_SinPN_pGLS_metriorhynchidae)
  R2_lik(EW_C1_PCo3_SinPN_pGLS_metriorhynchidae)

txtStop()


# ------------------------------------ #
# 12.4 Equal weights C3 paranasal GLMs #
# ------------------------------------ #

### Writes to file:
txtStart("output/ParanasalAnalyses/EW_C3/Files/EW_C3_paranasal_pGLS.txt")


EW_C3_SinPN_bm<-corBrownian(1,EW_C3_phylogeny_SinPN_Pruned_EqualMethod, form= ~Species_Sin)


### PCo1 pGLS
# PCo1 v terrestrial:
 EW_C3_PCo1_SinPN_pGLS_Terrestrial<-gls(PCo1_SinPN~ Terrestrial_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C3_SinPN_bm)
  summary(EW_C3_PCo1_SinPN_pGLS_Terrestrial)
  plot(EW_C3_PCo1_SinPN_pGLS_Terrestrial)
  AIC(EW_C3_PCo1_SinPN_pGLS_Terrestrial)
  R2_lik(EW_C3_PCo1_SinPN_pGLS_Terrestrial)

# PCo1 v semiaquatic:
 EW_C3_PCo1_SinPN_pGLS_semiaquatic<-gls(PCo1_SinPN~ Semiaquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C3_SinPN_bm)
  summary(EW_C3_PCo1_SinPN_pGLS_semiaquatic)
  plot(EW_C3_PCo1_SinPN_pGLS_semiaquatic)
  AIC(EW_C3_PCo1_SinPN_pGLS_semiaquatic)
  R2_lik(EW_C3_PCo1_SinPN_pGLS_semiaquatic)

# PCo1 v pelagic:
 EW_C3_PCo1_SinPN_pGLS_pelagic<-gls(PCo1_SinPN~ Pelagic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C3_SinPN_bm)
  summary(EW_C3_PCo1_SinPN_pGLS_pelagic)
  plot(EW_C3_PCo1_SinPN_pGLS_pelagic)
  AIC(EW_C3_PCo1_SinPN_pGLS_pelagic)
  R2_lik(EW_C3_PCo1_SinPN_pGLS_pelagic)

# PCo1 v marine:
 EW_C3_PCo1_SinPN_pGLS_marine<-gls(PCo1_SinPN~ Marine_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C3_SinPN_bm)
  summary(EW_C3_PCo1_SinPN_pGLS_marine)
  plot(EW_C3_PCo1_SinPN_pGLS_marine)
  AIC(EW_C3_PCo1_SinPN_pGLS_marine)
  R2_lik(EW_C3_PCo1_SinPN_pGLS_marine)

# PCo1 v freshwater:
 EW_C3_PCo1_SinPN_pGLS_freshwater<-gls(PCo1_SinPN~ Freshwater_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C3_SinPN_bm)
  summary(EW_C3_PCo1_SinPN_pGLS_freshwater)
  plot(EW_C3_PCo1_SinPN_pGLS_freshwater)
  AIC(EW_C3_PCo1_SinPN_pGLS_freshwater)
  R2_lik(EW_C3_PCo1_SinPN_pGLS_freshwater)

# PCo1 v aquatic:
 EW_C3_PCo1_SinPN_pGLS_aquatic<-gls(PCo1_SinPN~ Aquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C3_SinPN_bm)
  summary(EW_C3_PCo1_SinPN_pGLS_aquatic)
  plot(EW_C3_PCo1_SinPN_pGLS_aquatic)
  AIC(EW_C3_PCo1_SinPN_pGLS_aquatic)
  R2_lik(EW_C3_PCo1_SinPN_pGLS_aquatic)

# PCo1 v Thalattosuchia:
 EW_C3_PCo1_SinPN_pGLS_thalattosuchia<-gls(PCo1_SinPN~ Thalattosuchian_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C3_SinPN_bm)
  summary(EW_C3_PCo1_SinPN_pGLS_thalattosuchia)
  plot(EW_C3_PCo1_SinPN_pGLS_thalattosuchia)
  AIC(EW_C3_PCo1_SinPN_pGLS_thalattosuchia)
  R2_lik(EW_C3_PCo1_SinPN_pGLS_thalattosuchia)

# PCo1 v Metriorhynchidae:
 EW_C3_PCo1_SinPN_pGLS_metriorhynchidae<-gls(PCo1_SinPN~ Metriorhynchid_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C3_SinPN_bm)
  summary(EW_C3_PCo1_SinPN_pGLS_metriorhynchidae)
  plot(EW_C3_PCo1_SinPN_pGLS_metriorhynchidae)
  AIC(EW_C3_PCo1_SinPN_pGLS_metriorhynchidae)
  R2_lik(EW_C3_PCo1_SinPN_pGLS_metriorhynchidae)


### PCo2 pGLS
# PCo2 v terrestrial:
 EW_C3_PCo2_SinPN_pGLS_terrestrial<-gls(PCo2_SinPN~ Terrestrial_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C3_SinPN_bm)
  summary(EW_C3_PCo2_SinPN_pGLS_terrestrial)
  plot(EW_C3_PCo2_SinPN_pGLS_terrestrial)
  AIC(EW_C3_PCo2_SinPN_pGLS_terrestrial)
  R2_lik(EW_C3_PCo2_SinPN_pGLS_terrestrial) 

# PCo2 v semiaquatic:
 EW_C3_PCo2_SinPN_pGLS_semiaquatic<-gls(PCo2_SinPN~ Semiaquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C3_SinPN_bm)
  summary(EW_C3_PCo2_SinPN_pGLS_semiaquatic)
  plot(EW_C3_PCo2_SinPN_pGLS_semiaquatic)
  AIC(EW_C3_PCo2_SinPN_pGLS_semiaquatic)
  R2_lik(EW_C3_PCo2_SinPN_pGLS_semiaquatic) 

# PCo2 v pelagic:
 EW_C3_PCo2_SinPN_pGLS_pelagic<-gls(PCo2_SinPN~ Pelagic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C3_SinPN_bm)
  summary(EW_C3_PCo2_SinPN_pGLS_pelagic)
  plot(EW_C3_PCo2_SinPN_pGLS_pelagic)
  AIC(EW_C3_PCo2_SinPN_pGLS_pelagic)
  R2_lik(EW_C3_PCo2_SinPN_pGLS_pelagic) 

# PCo2 v marine:
 EW_C3_PCo2_SinPN_pGLS_marine<-gls(PCo2_SinPN~ Marine_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C3_SinPN_bm)
  summary(EW_C3_PCo2_SinPN_pGLS_marine)
  plot(EW_C3_PCo2_SinPN_pGLS_marine)
  AIC(EW_C3_PCo2_SinPN_pGLS_marine)
  R2_lik(EW_C3_PCo2_SinPN_pGLS_marine) 

# PCo2 v freshwater:
 EW_C3_PCo2_SinPN_pGLS_freshwater<-gls(PCo2_SinPN~ Freshwater_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C3_SinPN_bm)
  summary(EW_C3_PCo2_SinPN_pGLS_freshwater)
  plot(EW_C3_PCo2_SinPN_pGLS_freshwater)
  AIC(EW_C3_PCo2_SinPN_pGLS_freshwater)
  R2_lik(EW_C3_PCo2_SinPN_pGLS_freshwater) 

# PCo2 v aquatic:
 EW_C3_PCo2_SinPN_pGLS_aquatic<-gls(PCo2_SinPN~ Aquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C3_SinPN_bm)
  summary(EW_C3_PCo2_SinPN_pGLS_aquatic)
  plot(EW_C3_PCo2_SinPN_pGLS_aquatic)
  AIC(EW_C3_PCo2_SinPN_pGLS_aquatic)
  R2_lik(EW_C3_PCo2_SinPN_pGLS_aquatic) 

# PCo2 v Thalattosuchia:
 EW_C3_PCo2_SinPN_pGLS_thalattosuchia<-gls(PCo2_SinPN~ Thalattosuchian_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C3_SinPN_bm)
  summary(EW_C3_PCo2_SinPN_pGLS_thalattosuchia)
  plot(EW_C3_PCo2_SinPN_pGLS_thalattosuchia)
  AIC(EW_C3_PCo2_SinPN_pGLS_thalattosuchia)
  R2_lik(EW_C3_PCo2_SinPN_pGLS_thalattosuchia) 

# PCo2 v Metriorhynchidae:
 EW_C3_PCo2_SinPN_pGLS_metriorhynchidae<-gls(PCo2_SinPN~ Metriorhynchid_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C3_SinPN_bm)
  summary(EW_C3_PCo2_SinPN_pGLS_metriorhynchidae)
  plot(EW_C3_PCo2_SinPN_pGLS_metriorhynchidae)
  AIC(EW_C3_PCo2_SinPN_pGLS_metriorhynchidae)
  R2_lik(EW_C3_PCo2_SinPN_pGLS_metriorhynchidae) 


### PCo3 pGLS
# PCo3 v terrestrial:
 EW_C3_PCo3_SinPN_pGLS_terrestrial<-gls(PCo3_SinPN~ Terrestrial_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C3_SinPN_bm)
  summary(EW_C3_PCo3_SinPN_pGLS_terrestrial)
  plot(EW_C3_PCo3_SinPN_pGLS_terrestrial)
  AIC(EW_C3_PCo3_SinPN_pGLS_terrestrial)
  R2_lik(EW_C3_PCo3_SinPN_pGLS_terrestrial)

# PCo3 v semiaquatic:
 EW_C3_PCo3_SinPN_pGLS_semiaquatic<-gls(PCo3_SinPN~ Semiaquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C3_SinPN_bm)
  summary(EW_C3_PCo3_SinPN_pGLS_semiaquatic)
  plot(EW_C3_PCo3_SinPN_pGLS_semiaquatic)
  AIC(EW_C3_PCo3_SinPN_pGLS_semiaquatic)
  R2_lik(EW_C3_PCo3_SinPN_pGLS_semiaquatic)

# PCo3 v pelagic:
 EW_C3_PCo3_SinPN_pGLS_pelagic<-gls(PCo3_SinPN~ Pelagic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C3_SinPN_bm)
  summary(EW_C3_PCo3_SinPN_pGLS_pelagic)
  plot(EW_C3_PCo3_SinPN_pGLS_pelagic)
  AIC(EW_C3_PCo3_SinPN_pGLS_pelagic)
  R2_lik(EW_C3_PCo3_SinPN_pGLS_pelagic)

# PCo3 v marine:
 EW_C3_PCo3_SinPN_pGLS_marine<-gls(PCo3_SinPN~ Marine_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C3_SinPN_bm)
  summary(EW_C3_PCo3_SinPN_pGLS_marine)
  plot(EW_C3_PCo3_SinPN_pGLS_marine)
  AIC(EW_C3_PCo3_SinPN_pGLS_marine)
  R2_lik(EW_C3_PCo3_SinPN_pGLS_marine)

# PCo3 v freshwater:
 EW_C3_PCo3_SinPN_pGLS_freshwater<-gls(PCo3_SinPN~ Freshwater_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C3_SinPN_bm)
  summary(EW_C3_PCo3_SinPN_pGLS_freshwater)
  plot(EW_C3_PCo3_SinPN_pGLS_freshwater)
  AIC(EW_C3_PCo3_SinPN_pGLS_freshwater)
  R2_lik(EW_C3_PCo3_SinPN_pGLS_freshwater)

# PCo3 v aquatic:
 EW_C3_PCo3_SinPN_pGLS_aquatic<-gls(PCo3_SinPN~ Aquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C3_SinPN_bm)
  summary(EW_C3_PCo3_SinPN_pGLS_aquatic)
  plot(EW_C3_PCo3_SinPN_pGLS_aquatic)
  AIC(EW_C3_PCo3_SinPN_pGLS_aquatic)
  R2_lik(EW_C3_PCo3_SinPN_pGLS_aquatic)

# PCo3 v Thalattosuchia:
 EW_C3_PCo3_SinPN_pGLS_thalattosuchia<-gls(PCo3_SinPN~ Thalattosuchian_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C3_SinPN_bm)
  summary(EW_C3_PCo3_SinPN_pGLS_thalattosuchia)
  plot(EW_C3_PCo3_SinPN_pGLS_thalattosuchia)
  AIC(EW_C3_PCo3_SinPN_pGLS_thalattosuchia)
  R2_lik(EW_C3_PCo3_SinPN_pGLS_thalattosuchia)

# PCo3 v Metriorhynchidae:
 EW_C3_PCo3_SinPN_pGLS_metriorhynchidae<-gls(PCo3_SinPN~ Metriorhynchid_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C3_SinPN_bm)
  summary(EW_C3_PCo3_SinPN_pGLS_metriorhynchidae)
  plot(EW_C3_PCo3_SinPN_pGLS_metriorhynchidae)
  AIC(EW_C3_PCo3_SinPN_pGLS_metriorhynchidae)
  R2_lik(EW_C3_PCo3_SinPN_pGLS_metriorhynchidae)

txtStop()


# ------------------------------------ #
# 12.5 Equal weights C4 paranasal GLMs #
# ------------------------------------ #

### Writes to file:
txtStart("output/ParanasalAnalyses/EW_C4/Files/EW_C4_paranasal_pGLS.txt")


EW_C4_SinPN_bm<-corBrownian(1,EW_C4_phylogeny_SinPN_Pruned_EqualMethod, form= ~Species_Sin)


### PCo1 pGLS
# PCo1 v terrestrial:
 EW_C4_PCo1_SinPN_pGLS_Terrestrial<-gls(PCo1_SinPN~ Terrestrial_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C4_SinPN_bm)
  summary(EW_C4_PCo1_SinPN_pGLS_Terrestrial)
  plot(EW_C4_PCo1_SinPN_pGLS_Terrestrial)
  AIC(EW_C4_PCo1_SinPN_pGLS_Terrestrial)
  R2_lik(EW_C4_PCo1_SinPN_pGLS_Terrestrial)

# PCo1 v semiaquatic:
 EW_C4_PCo1_SinPN_pGLS_semiaquatic<-gls(PCo1_SinPN~ Semiaquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C4_SinPN_bm)
  summary(EW_C4_PCo1_SinPN_pGLS_semiaquatic)
  plot(EW_C4_PCo1_SinPN_pGLS_semiaquatic)
  AIC(EW_C4_PCo1_SinPN_pGLS_semiaquatic)
  R2_lik(EW_C4_PCo1_SinPN_pGLS_semiaquatic)

# PCo1 v pelagic:
 EW_C4_PCo1_SinPN_pGLS_pelagic<-gls(PCo1_SinPN~ Pelagic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C4_SinPN_bm)
  summary(EW_C4_PCo1_SinPN_pGLS_pelagic)
  plot(EW_C4_PCo1_SinPN_pGLS_pelagic)
  AIC(EW_C4_PCo1_SinPN_pGLS_pelagic)
  R2_lik(EW_C4_PCo1_SinPN_pGLS_pelagic)

# PCo1 v marine:
 EW_C4_PCo1_SinPN_pGLS_marine<-gls(PCo1_SinPN~ Marine_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C4_SinPN_bm)
  summary(EW_C4_PCo1_SinPN_pGLS_marine)
  plot(EW_C4_PCo1_SinPN_pGLS_marine)
  AIC(EW_C4_PCo1_SinPN_pGLS_marine)
  R2_lik(EW_C4_PCo1_SinPN_pGLS_marine)

# PCo1 v freshwater:
 EW_C4_PCo1_SinPN_pGLS_freshwater<-gls(PCo1_SinPN~ Freshwater_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C4_SinPN_bm)
  summary(EW_C4_PCo1_SinPN_pGLS_freshwater)
  plot(EW_C4_PCo1_SinPN_pGLS_freshwater)
  AIC(EW_C4_PCo1_SinPN_pGLS_freshwater)
  R2_lik(EW_C4_PCo1_SinPN_pGLS_freshwater)

# PCo1 v aquatic:
 EW_C4_PCo1_SinPN_pGLS_aquatic<-gls(PCo1_SinPN~ Aquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C4_SinPN_bm)
  summary(EW_C4_PCo1_SinPN_pGLS_aquatic)
  plot(EW_C4_PCo1_SinPN_pGLS_aquatic)
  AIC(EW_C4_PCo1_SinPN_pGLS_aquatic)
  R2_lik(EW_C4_PCo1_SinPN_pGLS_aquatic)

# PCo1 v Thalattosuchia:
 EW_C4_PCo1_SinPN_pGLS_thalattosuchia<-gls(PCo1_SinPN~ Thalattosuchian_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C4_SinPN_bm)
  summary(EW_C4_PCo1_SinPN_pGLS_thalattosuchia)
  plot(EW_C4_PCo1_SinPN_pGLS_thalattosuchia)
  AIC(EW_C4_PCo1_SinPN_pGLS_thalattosuchia)
  R2_lik(EW_C4_PCo1_SinPN_pGLS_thalattosuchia)

# PCo1 v Metriorhynchidae:
 EW_C4_PCo1_SinPN_pGLS_metriorhynchidae<-gls(PCo1_SinPN~ Metriorhynchid_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C4_SinPN_bm)
  summary(EW_C4_PCo1_SinPN_pGLS_metriorhynchidae)
  plot(EW_C4_PCo1_SinPN_pGLS_metriorhynchidae)
  AIC(EW_C4_PCo1_SinPN_pGLS_metriorhynchidae)
  R2_lik(EW_C4_PCo1_SinPN_pGLS_metriorhynchidae)


### PCo2 pGLS
# PCo2 v terrestrial:
 EW_C4_PCo2_SinPN_pGLS_terrestrial<-gls(PCo2_SinPN~ Terrestrial_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C4_SinPN_bm)
  summary(EW_C4_PCo2_SinPN_pGLS_terrestrial)
  plot(EW_C4_PCo2_SinPN_pGLS_terrestrial)
  AIC(EW_C4_PCo2_SinPN_pGLS_terrestrial)
  R2_lik(EW_C4_PCo2_SinPN_pGLS_terrestrial) 

# PCo2 v semiaquatic:
 EW_C4_PCo2_SinPN_pGLS_semiaquatic<-gls(PCo2_SinPN~ Semiaquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C4_SinPN_bm)
  summary(EW_C4_PCo2_SinPN_pGLS_semiaquatic)
  plot(EW_C4_PCo2_SinPN_pGLS_semiaquatic)
  AIC(EW_C4_PCo2_SinPN_pGLS_semiaquatic)
  R2_lik(EW_C4_PCo2_SinPN_pGLS_semiaquatic) 

# PCo2 v pelagic:
 EW_C4_PCo2_SinPN_pGLS_pelagic<-gls(PCo2_SinPN~ Pelagic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C4_SinPN_bm)
  summary(EW_C4_PCo2_SinPN_pGLS_pelagic)
  plot(EW_C4_PCo2_SinPN_pGLS_pelagic)
  AIC(EW_C4_PCo2_SinPN_pGLS_pelagic)
  R2_lik(EW_C4_PCo2_SinPN_pGLS_pelagic) 

# PCo2 v marine:
 EW_C4_PCo2_SinPN_pGLS_marine<-gls(PCo2_SinPN~ Marine_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C4_SinPN_bm)
  summary(EW_C4_PCo2_SinPN_pGLS_marine)
  plot(EW_C4_PCo2_SinPN_pGLS_marine)
  AIC(EW_C4_PCo2_SinPN_pGLS_marine)
  R2_lik(EW_C4_PCo2_SinPN_pGLS_marine) 

# PCo2 v freshwater:
 EW_C4_PCo2_SinPN_pGLS_freshwater<-gls(PCo2_SinPN~ Freshwater_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C4_SinPN_bm)
  summary(EW_C4_PCo2_SinPN_pGLS_freshwater)
  plot(EW_C4_PCo2_SinPN_pGLS_freshwater)
  AIC(EW_C4_PCo2_SinPN_pGLS_freshwater)
  R2_lik(EW_C4_PCo2_SinPN_pGLS_freshwater) 

# PCo2 v aquatic:
 EW_C4_PCo2_SinPN_pGLS_aquatic<-gls(PCo2_SinPN~ Aquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C4_SinPN_bm)
  summary(EW_C4_PCo2_SinPN_pGLS_aquatic)
  plot(EW_C4_PCo2_SinPN_pGLS_aquatic)
  AIC(EW_C4_PCo2_SinPN_pGLS_aquatic)
  R2_lik(EW_C4_PCo2_SinPN_pGLS_aquatic) 

# PCo2 v Thalattosuchia:
 EW_C4_PCo2_SinPN_pGLS_thalattosuchia<-gls(PCo2_SinPN~ Thalattosuchian_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C4_SinPN_bm)
  summary(EW_C4_PCo2_SinPN_pGLS_thalattosuchia)
  plot(EW_C4_PCo2_SinPN_pGLS_thalattosuchia)
  AIC(EW_C4_PCo2_SinPN_pGLS_thalattosuchia)
  R2_lik(EW_C4_PCo2_SinPN_pGLS_thalattosuchia) 

# PCo2 v Metriorhynchidae:
 EW_C4_PCo2_SinPN_pGLS_metriorhynchidae<-gls(PCo2_SinPN~ Metriorhynchid_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C4_SinPN_bm)
  summary(EW_C4_PCo2_SinPN_pGLS_metriorhynchidae)
  plot(EW_C4_PCo2_SinPN_pGLS_metriorhynchidae)
  AIC(EW_C4_PCo2_SinPN_pGLS_metriorhynchidae)
  R2_lik(EW_C4_PCo2_SinPN_pGLS_metriorhynchidae) 


### PCo3 pGLS
# PCo3 v terrestrial:
 EW_C4_PCo3_SinPN_pGLS_terrestrial<-gls(PCo3_SinPN~ Terrestrial_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C4_SinPN_bm)
  summary(EW_C4_PCo3_SinPN_pGLS_terrestrial)
  plot(EW_C4_PCo3_SinPN_pGLS_terrestrial)
  AIC(EW_C4_PCo3_SinPN_pGLS_terrestrial)
  R2_lik(EW_C4_PCo3_SinPN_pGLS_terrestrial)

# PCo3 v semiaquatic:
 EW_C4_PCo3_SinPN_pGLS_semiaquatic<-gls(PCo3_SinPN~ Semiaquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C4_SinPN_bm)
  summary(EW_C4_PCo3_SinPN_pGLS_semiaquatic)
  plot(EW_C4_PCo3_SinPN_pGLS_semiaquatic)
  AIC(EW_C4_PCo3_SinPN_pGLS_semiaquatic)
  R2_lik(EW_C4_PCo3_SinPN_pGLS_semiaquatic)

# PCo3 v pelagic:
 EW_C4_PCo3_SinPN_pGLS_pelagic<-gls(PCo3_SinPN~ Pelagic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C4_SinPN_bm)
  summary(EW_C4_PCo3_SinPN_pGLS_pelagic)
  plot(EW_C4_PCo3_SinPN_pGLS_pelagic)
  AIC(EW_C4_PCo3_SinPN_pGLS_pelagic)
  R2_lik(EW_C4_PCo3_SinPN_pGLS_pelagic)

# PCo3 v marine:
 EW_C4_PCo3_SinPN_pGLS_marine<-gls(PCo3_SinPN~ Marine_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C4_SinPN_bm)
  summary(EW_C4_PCo3_SinPN_pGLS_marine)
  plot(EW_C4_PCo3_SinPN_pGLS_marine)
  AIC(EW_C4_PCo3_SinPN_pGLS_marine)
  R2_lik(EW_C4_PCo3_SinPN_pGLS_marine)

# PCo3 v freshwater:
 EW_C4_PCo3_SinPN_pGLS_freshwater<-gls(PCo3_SinPN~ Freshwater_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C4_SinPN_bm)
  summary(EW_C4_PCo3_SinPN_pGLS_freshwater)
  plot(EW_C4_PCo3_SinPN_pGLS_freshwater)
  AIC(EW_C4_PCo3_SinPN_pGLS_freshwater)
  R2_lik(EW_C4_PCo3_SinPN_pGLS_freshwater)

# PCo3 v aquatic:
 EW_C4_PCo3_SinPN_pGLS_aquatic<-gls(PCo3_SinPN~ Aquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C4_SinPN_bm)
  summary(EW_C4_PCo3_SinPN_pGLS_aquatic)
  plot(EW_C4_PCo3_SinPN_pGLS_aquatic)
  AIC(EW_C4_PCo3_SinPN_pGLS_aquatic)
  R2_lik(EW_C4_PCo3_SinPN_pGLS_aquatic)

# PCo3 v Thalattosuchia:
 EW_C4_PCo3_SinPN_pGLS_thalattosuchia<-gls(PCo3_SinPN~ Thalattosuchian_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C4_SinPN_bm)
  summary(EW_C4_PCo3_SinPN_pGLS_thalattosuchia)
  plot(EW_C4_PCo3_SinPN_pGLS_thalattosuchia)
  AIC(EW_C4_PCo3_SinPN_pGLS_thalattosuchia)
  R2_lik(EW_C4_PCo3_SinPN_pGLS_thalattosuchia)

# PCo3 v Metriorhynchidae:
 EW_C4_PCo3_SinPN_pGLS_metriorhynchidae<-gls(PCo3_SinPN~ Metriorhynchid_Sin, method= "ML", data= CrocSinPNData, correlation = EW_C4_SinPN_bm)
  summary(EW_C4_PCo3_SinPN_pGLS_metriorhynchidae)
  plot(EW_C4_PCo3_SinPN_pGLS_metriorhynchidae)
  AIC(EW_C4_PCo3_SinPN_pGLS_metriorhynchidae)
  R2_lik(EW_C4_PCo3_SinPN_pGLS_metriorhynchidae)

txtStop()


# --------------------------------------------------------- #
# 12.6 Extended implied weights (k=15) C0 paranasal GLMs    #
# --------------------------------------------------------- #

### Writes to file:
txtStart("output/ParanasalAnalyses/EIWk15_C0/Files/EIWk15_C0_paranasal_pGLS.txt")


EIWk15_C0_SinPN_bm<-corBrownian(1,EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, form= ~Species_Sin)


### PCo1 pGLS
# PCo1 v terrestrial:
 EIWk15_C0_PCo1_SinPN_pGLS_Terrestrial<-gls(PCo1_SinPN~ Terrestrial_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C0_SinPN_bm)
  summary(EIWk15_C0_PCo1_SinPN_pGLS_Terrestrial)
  plot(EIWk15_C0_PCo1_SinPN_pGLS_Terrestrial)
  AIC(EIWk15_C0_PCo1_SinPN_pGLS_Terrestrial)
  R2_lik(EIWk15_C0_PCo1_SinPN_pGLS_Terrestrial)

# PCo1 v semiaquatic:
 EIWk15_C0_PCo1_SinPN_pGLS_semiaquatic<-gls(PCo1_SinPN~ Semiaquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C0_SinPN_bm)
  summary(EIWk15_C0_PCo1_SinPN_pGLS_semiaquatic)
  plot(EIWk15_C0_PCo1_SinPN_pGLS_semiaquatic)
  AIC(EIWk15_C0_PCo1_SinPN_pGLS_semiaquatic)
  R2_lik(EIWk15_C0_PCo1_SinPN_pGLS_semiaquatic)

# PCo1 v pelagic:
 EIWk15_C0_PCo1_SinPN_pGLS_pelagic<-gls(PCo1_SinPN~ Pelagic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C0_SinPN_bm)
  summary(EIWk15_C0_PCo1_SinPN_pGLS_pelagic)
  plot(EIWk15_C0_PCo1_SinPN_pGLS_pelagic)
  AIC(EIWk15_C0_PCo1_SinPN_pGLS_pelagic)
  R2_lik(EIWk15_C0_PCo1_SinPN_pGLS_pelagic)

# PCo1 v marine:
 EIWk15_C0_PCo1_SinPN_pGLS_marine<-gls(PCo1_SinPN~ Marine_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C0_SinPN_bm)
  summary(EIWk15_C0_PCo1_SinPN_pGLS_marine)
  plot(EIWk15_C0_PCo1_SinPN_pGLS_marine)
  AIC(EIWk15_C0_PCo1_SinPN_pGLS_marine)
  R2_lik(EIWk15_C0_PCo1_SinPN_pGLS_marine)

# PCo1 v freshwater:
 EIWk15_C0_PCo1_SinPN_pGLS_freshwater<-gls(PCo1_SinPN~ Freshwater_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C0_SinPN_bm)
  summary(EIWk15_C0_PCo1_SinPN_pGLS_freshwater)
  plot(EIWk15_C0_PCo1_SinPN_pGLS_freshwater)
  AIC(EIWk15_C0_PCo1_SinPN_pGLS_freshwater)
  R2_lik(EIWk15_C0_PCo1_SinPN_pGLS_freshwater)

# PCo1 v aquatic:
 EIWk15_C0_PCo1_SinPN_pGLS_aquatic<-gls(PCo1_SinPN~ Aquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C0_SinPN_bm)
  summary(EIWk15_C0_PCo1_SinPN_pGLS_aquatic)
  plot(EIWk15_C0_PCo1_SinPN_pGLS_aquatic)
  AIC(EIWk15_C0_PCo1_SinPN_pGLS_aquatic)
  R2_lik(EIWk15_C0_PCo1_SinPN_pGLS_aquatic)

# PCo1 v Thalattosuchia:
 EIWk15_C0_PCo1_SinPN_pGLS_thalattosuchia<-gls(PCo1_SinPN~ Thalattosuchian_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C0_SinPN_bm)
  summary(EIWk15_C0_PCo1_SinPN_pGLS_thalattosuchia)
  plot(EIWk15_C0_PCo1_SinPN_pGLS_thalattosuchia)
  AIC(EIWk15_C0_PCo1_SinPN_pGLS_thalattosuchia)
  R2_lik(EIWk15_C0_PCo1_SinPN_pGLS_thalattosuchia)

# PCo1 v Metriorhynchidae:
 EIWk15_C0_PCo1_SinPN_pGLS_metriorhynchidae<-gls(PCo1_SinPN~ Metriorhynchid_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C0_SinPN_bm)
  summary(EIWk15_C0_PCo1_SinPN_pGLS_metriorhynchidae)
  plot(EIWk15_C0_PCo1_SinPN_pGLS_metriorhynchidae)
  AIC(EIWk15_C0_PCo1_SinPN_pGLS_metriorhynchidae)
  R2_lik(EIWk15_C0_PCo1_SinPN_pGLS_metriorhynchidae)


### PCo2 pGLS
# PCo2 v terrestrial:
 EIWk15_C0_PCo2_SinPN_pGLS_terrestrial<-gls(PCo2_SinPN~ Terrestrial_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C0_SinPN_bm)
  summary(EIWk15_C0_PCo2_SinPN_pGLS_terrestrial)
  plot(EIWk15_C0_PCo2_SinPN_pGLS_terrestrial)
  AIC(EIWk15_C0_PCo2_SinPN_pGLS_terrestrial)
  R2_lik(EIWk15_C0_PCo2_SinPN_pGLS_terrestrial) 

# PCo2 v semiaquatic:
 EIWk15_C0_PCo2_SinPN_pGLS_semiaquatic<-gls(PCo2_SinPN~ Semiaquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C0_SinPN_bm)
  summary(EIWk15_C0_PCo2_SinPN_pGLS_semiaquatic)
  plot(EIWk15_C0_PCo2_SinPN_pGLS_semiaquatic)
  AIC(EIWk15_C0_PCo2_SinPN_pGLS_semiaquatic)
  R2_lik(EIWk15_C0_PCo2_SinPN_pGLS_semiaquatic) 

# PCo2 v pelagic:
 EIWk15_C0_PCo2_SinPN_pGLS_pelagic<-gls(PCo2_SinPN~ Pelagic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C0_SinPN_bm)
  summary(EIWk15_C0_PCo2_SinPN_pGLS_pelagic)
  plot(EIWk15_C0_PCo2_SinPN_pGLS_pelagic)
  AIC(EIWk15_C0_PCo2_SinPN_pGLS_pelagic)
  R2_lik(EIWk15_C0_PCo2_SinPN_pGLS_pelagic) 

# PCo2 v marine:
 EIWk15_C0_PCo2_SinPN_pGLS_marine<-gls(PCo2_SinPN~ Marine_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C0_SinPN_bm)
  summary(EIWk15_C0_PCo2_SinPN_pGLS_marine)
  plot(EIWk15_C0_PCo2_SinPN_pGLS_marine)
  AIC(EIWk15_C0_PCo2_SinPN_pGLS_marine)
  R2_lik(EIWk15_C0_PCo2_SinPN_pGLS_marine) 

# PCo2 v freshwater:
 EIWk15_C0_PCo2_SinPN_pGLS_freshwater<-gls(PCo2_SinPN~ Freshwater_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C0_SinPN_bm)
  summary(EIWk15_C0_PCo2_SinPN_pGLS_freshwater)
  plot(EIWk15_C0_PCo2_SinPN_pGLS_freshwater)
  AIC(EIWk15_C0_PCo2_SinPN_pGLS_freshwater)
  R2_lik(EIWk15_C0_PCo2_SinPN_pGLS_freshwater) 

# PCo2 v aquatic:
 EIWk15_C0_PCo2_SinPN_pGLS_aquatic<-gls(PCo2_SinPN~ Aquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C0_SinPN_bm)
  summary(EIWk15_C0_PCo2_SinPN_pGLS_aquatic)
  plot(EIWk15_C0_PCo2_SinPN_pGLS_aquatic)
  AIC(EIWk15_C0_PCo2_SinPN_pGLS_aquatic)
  R2_lik(EIWk15_C0_PCo2_SinPN_pGLS_aquatic) 

# PCo2 v Thalattosuchia:
 EIWk15_C0_PCo2_SinPN_pGLS_thalattosuchia<-gls(PCo2_SinPN~ Thalattosuchian_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C0_SinPN_bm)
  summary(EIWk15_C0_PCo2_SinPN_pGLS_thalattosuchia)
  plot(EIWk15_C0_PCo2_SinPN_pGLS_thalattosuchia)
  AIC(EIWk15_C0_PCo2_SinPN_pGLS_thalattosuchia)
  R2_lik(EIWk15_C0_PCo2_SinPN_pGLS_thalattosuchia) 

# PCo2 v Metriorhynchidae:
 EIWk15_C0_PCo2_SinPN_pGLS_metriorhynchidae<-gls(PCo2_SinPN~ Metriorhynchid_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C0_SinPN_bm)
  summary(EIWk15_C0_PCo2_SinPN_pGLS_metriorhynchidae)
  plot(EIWk15_C0_PCo2_SinPN_pGLS_metriorhynchidae)
  AIC(EIWk15_C0_PCo2_SinPN_pGLS_metriorhynchidae)
  R2_lik(EIWk15_C0_PCo2_SinPN_pGLS_metriorhynchidae) 


### PCo3 pGLS
# PCo3 v terrestrial:
 EIWk15_C0_PCo3_SinPN_pGLS_terrestrial<-gls(PCo3_SinPN~ Terrestrial_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C0_SinPN_bm)
  summary(EIWk15_C0_PCo3_SinPN_pGLS_terrestrial)
  plot(EIWk15_C0_PCo3_SinPN_pGLS_terrestrial)
  AIC(EIWk15_C0_PCo3_SinPN_pGLS_terrestrial)
  R2_lik(EIWk15_C0_PCo3_SinPN_pGLS_terrestrial)

# PCo3 v semiaquatic:
 EIWk15_C0_PCo3_SinPN_pGLS_semiaquatic<-gls(PCo3_SinPN~ Semiaquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C0_SinPN_bm)
  summary(EIWk15_C0_PCo3_SinPN_pGLS_semiaquatic)
  plot(EIWk15_C0_PCo3_SinPN_pGLS_semiaquatic)
  AIC(EIWk15_C0_PCo3_SinPN_pGLS_semiaquatic)
  R2_lik(EIWk15_C0_PCo3_SinPN_pGLS_semiaquatic)

# PCo3 v pelagic:
 EIWk15_C0_PCo3_SinPN_pGLS_pelagic<-gls(PCo3_SinPN~ Pelagic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C0_SinPN_bm)
  summary(EIWk15_C0_PCo3_SinPN_pGLS_pelagic)
  plot(EIWk15_C0_PCo3_SinPN_pGLS_pelagic)
  AIC(EIWk15_C0_PCo3_SinPN_pGLS_pelagic)
  R2_lik(EIWk15_C0_PCo3_SinPN_pGLS_pelagic)

# PCo3 v marine:
 EIWk15_C0_PCo3_SinPN_pGLS_marine<-gls(PCo3_SinPN~ Marine_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C0_SinPN_bm)
  summary(EIWk15_C0_PCo3_SinPN_pGLS_marine)
  plot(EIWk15_C0_PCo3_SinPN_pGLS_marine)
  AIC(EIWk15_C0_PCo3_SinPN_pGLS_marine)
  R2_lik(EIWk15_C0_PCo3_SinPN_pGLS_marine)

# PCo3 v freshwater:
 EIWk15_C0_PCo3_SinPN_pGLS_freshwater<-gls(PCo3_SinPN~ Freshwater_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C0_SinPN_bm)
  summary(EIWk15_C0_PCo3_SinPN_pGLS_freshwater)
  plot(EIWk15_C0_PCo3_SinPN_pGLS_freshwater)
  AIC(EIWk15_C0_PCo3_SinPN_pGLS_freshwater)
  R2_lik(EIWk15_C0_PCo3_SinPN_pGLS_freshwater)

# PCo3 v aquatic:
 EIWk15_C0_PCo3_SinPN_pGLS_aquatic<-gls(PCo3_SinPN~ Aquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C0_SinPN_bm)
  summary(EIWk15_C0_PCo3_SinPN_pGLS_aquatic)
  plot(EIWk15_C0_PCo3_SinPN_pGLS_aquatic)
  AIC(EIWk15_C0_PCo3_SinPN_pGLS_aquatic)
  R2_lik(EIWk15_C0_PCo3_SinPN_pGLS_aquatic)

# PCo3 v Thalattosuchia:
 EIWk15_C0_PCo3_SinPN_pGLS_thalattosuchia<-gls(PCo3_SinPN~ Thalattosuchian_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C0_SinPN_bm)
  summary(EIWk15_C0_PCo3_SinPN_pGLS_thalattosuchia)
  plot(EIWk15_C0_PCo3_SinPN_pGLS_thalattosuchia)
  AIC(EIWk15_C0_PCo3_SinPN_pGLS_thalattosuchia)
  R2_lik(EIWk15_C0_PCo3_SinPN_pGLS_thalattosuchia)

# PCo3 v Metriorhynchidae:
 EIWk15_C0_PCo3_SinPN_pGLS_metriorhynchidae<-gls(PCo3_SinPN~ Metriorhynchid_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C0_SinPN_bm)
  summary(EIWk15_C0_PCo3_SinPN_pGLS_metriorhynchidae)
  plot(EIWk15_C0_PCo3_SinPN_pGLS_metriorhynchidae)
  AIC(EIWk15_C0_PCo3_SinPN_pGLS_metriorhynchidae)
  R2_lik(EIWk15_C0_PCo3_SinPN_pGLS_metriorhynchidae)

txtStop()


# ------------------------------------------------------ #
# 12.7 Extended implied weights (k=15) C1 paranasal GLMs #
# ------------------------------------------------------ #

### Writes to file:
txtStart("output/ParanasalAnalyses/EIWk15_C1/Files/EIWk15_C1_paranasal_pGLS.txt")


EIWk15_C1_SinPN_bm<-corBrownian(1,EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, form= ~Species_Sin)


### PCo1 pGLS
# PCo1 v terrestrial:
 EIWk15_C1_PCo1_SinPN_pGLS_Terrestrial<-gls(PCo1_SinPN~ Terrestrial_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C1_SinPN_bm)
  summary(EIWk15_C1_PCo1_SinPN_pGLS_Terrestrial)
  plot(EIWk15_C1_PCo1_SinPN_pGLS_Terrestrial)
  AIC(EIWk15_C1_PCo1_SinPN_pGLS_Terrestrial)
  R2_lik(EIWk15_C1_PCo1_SinPN_pGLS_Terrestrial)

# PCo1 v semiaquatic:
 EIWk15_C1_PCo1_SinPN_pGLS_semiaquatic<-gls(PCo1_SinPN~ Semiaquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C1_SinPN_bm)
  summary(EIWk15_C1_PCo1_SinPN_pGLS_semiaquatic)
  plot(EIWk15_C1_PCo1_SinPN_pGLS_semiaquatic)
  AIC(EIWk15_C1_PCo1_SinPN_pGLS_semiaquatic)
  R2_lik(EIWk15_C1_PCo1_SinPN_pGLS_semiaquatic)

# PCo1 v pelagic:
 EIWk15_C1_PCo1_SinPN_pGLS_pelagic<-gls(PCo1_SinPN~ Pelagic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C1_SinPN_bm)
  summary(EIWk15_C1_PCo1_SinPN_pGLS_pelagic)
  plot(EIWk15_C1_PCo1_SinPN_pGLS_pelagic)
  AIC(EIWk15_C1_PCo1_SinPN_pGLS_pelagic)
  R2_lik(EIWk15_C1_PCo1_SinPN_pGLS_pelagic)

# PCo1 v marine:
 EIWk15_C1_PCo1_SinPN_pGLS_marine<-gls(PCo1_SinPN~ Marine_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C1_SinPN_bm)
  summary(EIWk15_C1_PCo1_SinPN_pGLS_marine)
  plot(EIWk15_C1_PCo1_SinPN_pGLS_marine)
  AIC(EIWk15_C1_PCo1_SinPN_pGLS_marine)
  R2_lik(EIWk15_C1_PCo1_SinPN_pGLS_marine)

# PCo1 v freshwater:
 EIWk15_C1_PCo1_SinPN_pGLS_freshwater<-gls(PCo1_SinPN~ Freshwater_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C1_SinPN_bm)
  summary(EIWk15_C1_PCo1_SinPN_pGLS_freshwater)
  plot(EIWk15_C1_PCo1_SinPN_pGLS_freshwater)
  AIC(EIWk15_C1_PCo1_SinPN_pGLS_freshwater)
  R2_lik(EIWk15_C1_PCo1_SinPN_pGLS_freshwater)

# PCo1 v aquatic:
 EIWk15_C1_PCo1_SinPN_pGLS_aquatic<-gls(PCo1_SinPN~ Aquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C1_SinPN_bm)
  summary(EIWk15_C1_PCo1_SinPN_pGLS_aquatic)
  plot(EIWk15_C1_PCo1_SinPN_pGLS_aquatic)
  AIC(EIWk15_C1_PCo1_SinPN_pGLS_aquatic)
  R2_lik(EIWk15_C1_PCo1_SinPN_pGLS_aquatic)

# PCo1 v Thalattosuchia:
 EIWk15_C1_PCo1_SinPN_pGLS_thalattosuchia<-gls(PCo1_SinPN~ Thalattosuchian_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C1_SinPN_bm)
  summary(EIWk15_C1_PCo1_SinPN_pGLS_thalattosuchia)
  plot(EIWk15_C1_PCo1_SinPN_pGLS_thalattosuchia)
  AIC(EIWk15_C1_PCo1_SinPN_pGLS_thalattosuchia)
  R2_lik(EIWk15_C1_PCo1_SinPN_pGLS_thalattosuchia)

# PCo1 v Metriorhynchidae:
 EIWk15_C1_PCo1_SinPN_pGLS_metriorhynchidae<-gls(PCo1_SinPN~ Metriorhynchid_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C1_SinPN_bm)
  summary(EIWk15_C1_PCo1_SinPN_pGLS_metriorhynchidae)
  plot(EIWk15_C1_PCo1_SinPN_pGLS_metriorhynchidae)
  AIC(EIWk15_C1_PCo1_SinPN_pGLS_metriorhynchidae)
  R2_lik(EIWk15_C1_PCo1_SinPN_pGLS_metriorhynchidae)


### PCo2 pGLS
# PCo2 v terrestrial:
 EIWk15_C1_PCo2_SinPN_pGLS_terrestrial<-gls(PCo2_SinPN~ Terrestrial_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C1_SinPN_bm)
  summary(EIWk15_C1_PCo2_SinPN_pGLS_terrestrial)
  plot(EIWk15_C1_PCo2_SinPN_pGLS_terrestrial)
  AIC(EIWk15_C1_PCo2_SinPN_pGLS_terrestrial)
  R2_lik(EIWk15_C1_PCo2_SinPN_pGLS_terrestrial) 

# PCo2 v semiaquatic:
 EIWk15_C1_PCo2_SinPN_pGLS_semiaquatic<-gls(PCo2_SinPN~ Semiaquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C1_SinPN_bm)
  summary(EIWk15_C1_PCo2_SinPN_pGLS_semiaquatic)
  plot(EIWk15_C1_PCo2_SinPN_pGLS_semiaquatic)
  AIC(EIWk15_C1_PCo2_SinPN_pGLS_semiaquatic)
  R2_lik(EIWk15_C1_PCo2_SinPN_pGLS_semiaquatic) 

# PCo2 v pelagic:
 EIWk15_C1_PCo2_SinPN_pGLS_pelagic<-gls(PCo2_SinPN~ Pelagic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C1_SinPN_bm)
  summary(EIWk15_C1_PCo2_SinPN_pGLS_pelagic)
  plot(EIWk15_C1_PCo2_SinPN_pGLS_pelagic)
  AIC(EIWk15_C1_PCo2_SinPN_pGLS_pelagic)
  R2_lik(EIWk15_C1_PCo2_SinPN_pGLS_pelagic) 

# PCo2 v marine:
 EIWk15_C1_PCo2_SinPN_pGLS_marine<-gls(PCo2_SinPN~ Marine_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C1_SinPN_bm)
  summary(EIWk15_C1_PCo2_SinPN_pGLS_marine)
  plot(EIWk15_C1_PCo2_SinPN_pGLS_marine)
  AIC(EIWk15_C1_PCo2_SinPN_pGLS_marine)
  R2_lik(EIWk15_C1_PCo2_SinPN_pGLS_marine) 

# PCo2 v freshwater:
 EIWk15_C1_PCo2_SinPN_pGLS_freshwater<-gls(PCo2_SinPN~ Freshwater_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C1_SinPN_bm)
  summary(EIWk15_C1_PCo2_SinPN_pGLS_freshwater)
  plot(EIWk15_C1_PCo2_SinPN_pGLS_freshwater)
  AIC(EIWk15_C1_PCo2_SinPN_pGLS_freshwater)
  R2_lik(EIWk15_C1_PCo2_SinPN_pGLS_freshwater) 

# PCo2 v aquatic:
 EIWk15_C1_PCo2_SinPN_pGLS_aquatic<-gls(PCo2_SinPN~ Aquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C1_SinPN_bm)
  summary(EIWk15_C1_PCo2_SinPN_pGLS_aquatic)
  plot(EIWk15_C1_PCo2_SinPN_pGLS_aquatic)
  AIC(EIWk15_C1_PCo2_SinPN_pGLS_aquatic)
  R2_lik(EIWk15_C1_PCo2_SinPN_pGLS_aquatic) 

# PCo2 v Thalattosuchia:
 EIWk15_C1_PCo2_SinPN_pGLS_thalattosuchia<-gls(PCo2_SinPN~ Thalattosuchian_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C1_SinPN_bm)
  summary(EIWk15_C1_PCo2_SinPN_pGLS_thalattosuchia)
  plot(EIWk15_C1_PCo2_SinPN_pGLS_thalattosuchia)
  AIC(EIWk15_C1_PCo2_SinPN_pGLS_thalattosuchia)
  R2_lik(EIWk15_C1_PCo2_SinPN_pGLS_thalattosuchia) 

# PCo2 v Metriorhynchidae:
 EIWk15_C1_PCo2_SinPN_pGLS_metriorhynchidae<-gls(PCo2_SinPN~ Metriorhynchid_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C1_SinPN_bm)
  summary(EIWk15_C1_PCo2_SinPN_pGLS_metriorhynchidae)
  plot(EIWk15_C1_PCo2_SinPN_pGLS_metriorhynchidae)
  AIC(EIWk15_C1_PCo2_SinPN_pGLS_metriorhynchidae)
  R2_lik(EIWk15_C1_PCo2_SinPN_pGLS_metriorhynchidae) 


### PCo3 pGLS
# PCo3 v terrestrial:
 EIWk15_C1_PCo3_SinPN_pGLS_terrestrial<-gls(PCo3_SinPN~ Terrestrial_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C1_SinPN_bm)
  summary(EIWk15_C1_PCo3_SinPN_pGLS_terrestrial)
  plot(EIWk15_C1_PCo3_SinPN_pGLS_terrestrial)
  AIC(EIWk15_C1_PCo3_SinPN_pGLS_terrestrial)
  R2_lik(EIWk15_C1_PCo3_SinPN_pGLS_terrestrial)

# PCo3 v semiaquatic:
 EIWk15_C1_PCo3_SinPN_pGLS_semiaquatic<-gls(PCo3_SinPN~ Semiaquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C1_SinPN_bm)
  summary(EIWk15_C1_PCo3_SinPN_pGLS_semiaquatic)
  plot(EIWk15_C1_PCo3_SinPN_pGLS_semiaquatic)
  AIC(EIWk15_C1_PCo3_SinPN_pGLS_semiaquatic)
  R2_lik(EIWk15_C1_PCo3_SinPN_pGLS_semiaquatic)

# PCo3 v pelagic:
 EIWk15_C1_PCo3_SinPN_pGLS_pelagic<-gls(PCo3_SinPN~ Pelagic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C1_SinPN_bm)
  summary(EIWk15_C1_PCo3_SinPN_pGLS_pelagic)
  plot(EIWk15_C1_PCo3_SinPN_pGLS_pelagic)
  AIC(EIWk15_C1_PCo3_SinPN_pGLS_pelagic)
  R2_lik(EIWk15_C1_PCo3_SinPN_pGLS_pelagic)

# PCo3 v marine:
 EIWk15_C1_PCo3_SinPN_pGLS_marine<-gls(PCo3_SinPN~ Marine_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C1_SinPN_bm)
  summary(EIWk15_C1_PCo3_SinPN_pGLS_marine)
  plot(EIWk15_C1_PCo3_SinPN_pGLS_marine)
  AIC(EIWk15_C1_PCo3_SinPN_pGLS_marine)
  R2_lik(EIWk15_C1_PCo3_SinPN_pGLS_marine)

# PCo3 v freshwater:
 EIWk15_C1_PCo3_SinPN_pGLS_freshwater<-gls(PCo3_SinPN~ Freshwater_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C1_SinPN_bm)
  summary(EIWk15_C1_PCo3_SinPN_pGLS_freshwater)
  plot(EIWk15_C1_PCo3_SinPN_pGLS_freshwater)
  AIC(EIWk15_C1_PCo3_SinPN_pGLS_freshwater)
  R2_lik(EIWk15_C1_PCo3_SinPN_pGLS_freshwater)

# PCo3 v aquatic:
 EIWk15_C1_PCo3_SinPN_pGLS_aquatic<-gls(PCo3_SinPN~ Aquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C1_SinPN_bm)
  summary(EIWk15_C1_PCo3_SinPN_pGLS_aquatic)
  plot(EIWk15_C1_PCo3_SinPN_pGLS_aquatic)
  AIC(EIWk15_C1_PCo3_SinPN_pGLS_aquatic)
  R2_lik(EIWk15_C1_PCo3_SinPN_pGLS_aquatic)

# PCo3 v Thalattosuchia:
 EIWk15_C1_PCo3_SinPN_pGLS_thalattosuchia<-gls(PCo3_SinPN~ Thalattosuchian_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C1_SinPN_bm)
  summary(EIWk15_C1_PCo3_SinPN_pGLS_thalattosuchia)
  plot(EIWk15_C1_PCo3_SinPN_pGLS_thalattosuchia)
  AIC(EIWk15_C1_PCo3_SinPN_pGLS_thalattosuchia)
  R2_lik(EIWk15_C1_PCo3_SinPN_pGLS_thalattosuchia)

# PCo3 v Metriorhynchidae:
 EIWk15_C1_PCo3_SinPN_pGLS_metriorhynchidae<-gls(PCo3_SinPN~ Metriorhynchid_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C1_SinPN_bm)
  summary(EIWk15_C1_PCo3_SinPN_pGLS_metriorhynchidae)
  plot(EIWk15_C1_PCo3_SinPN_pGLS_metriorhynchidae)
  AIC(EIWk15_C1_PCo3_SinPN_pGLS_metriorhynchidae)
  R2_lik(EIWk15_C1_PCo3_SinPN_pGLS_metriorhynchidae)

txtStop()


# ------------------------------------------------------ #
# 12.8 Extended implied weights (k=15) C3 paranasal GLMs #
# ------------------------------------------------------ #

### Writes to file:
txtStart("output/ParanasalAnalyses/EIWk15_C3/Files/EIWk15_C3_paranasal_pGLS.txt")


EIWk15_C3_SinPN_bm<-corBrownian(1,EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, form= ~Species_Sin)


### PCo1 pGLS
# PCo1 v terrestrial:
 EIWk15_C3_PCo1_SinPN_pGLS_Terrestrial<-gls(PCo1_SinPN~ Terrestrial_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C3_SinPN_bm)
  summary(EIWk15_C3_PCo1_SinPN_pGLS_Terrestrial)
  plot(EIWk15_C3_PCo1_SinPN_pGLS_Terrestrial)
  AIC(EIWk15_C3_PCo1_SinPN_pGLS_Terrestrial)
  R2_lik(EIWk15_C3_PCo1_SinPN_pGLS_Terrestrial)

# PCo1 v semiaquatic:
 EIWk15_C3_PCo1_SinPN_pGLS_semiaquatic<-gls(PCo1_SinPN~ Semiaquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C3_SinPN_bm)
  summary(EIWk15_C3_PCo1_SinPN_pGLS_semiaquatic)
  plot(EIWk15_C3_PCo1_SinPN_pGLS_semiaquatic)
  AIC(EIWk15_C3_PCo1_SinPN_pGLS_semiaquatic)
  R2_lik(EIWk15_C3_PCo1_SinPN_pGLS_semiaquatic)

# PCo1 v pelagic:
 EIWk15_C3_PCo1_SinPN_pGLS_pelagic<-gls(PCo1_SinPN~ Pelagic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C3_SinPN_bm)
  summary(EIWk15_C3_PCo1_SinPN_pGLS_pelagic)
  plot(EIWk15_C3_PCo1_SinPN_pGLS_pelagic)
  AIC(EIWk15_C3_PCo1_SinPN_pGLS_pelagic)
  R2_lik(EIWk15_C3_PCo1_SinPN_pGLS_pelagic)

# PCo1 v marine:
 EIWk15_C3_PCo1_SinPN_pGLS_marine<-gls(PCo1_SinPN~ Marine_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C3_SinPN_bm)
  summary(EIWk15_C3_PCo1_SinPN_pGLS_marine)
  plot(EIWk15_C3_PCo1_SinPN_pGLS_marine)
  AIC(EIWk15_C3_PCo1_SinPN_pGLS_marine)
  R2_lik(EIWk15_C3_PCo1_SinPN_pGLS_marine)

# PCo1 v freshwater:
 EIWk15_C3_PCo1_SinPN_pGLS_freshwater<-gls(PCo1_SinPN~ Freshwater_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C3_SinPN_bm)
  summary(EIWk15_C3_PCo1_SinPN_pGLS_freshwater)
  plot(EIWk15_C3_PCo1_SinPN_pGLS_freshwater)
  AIC(EIWk15_C3_PCo1_SinPN_pGLS_freshwater)
  R2_lik(EIWk15_C3_PCo1_SinPN_pGLS_freshwater)

# PCo1 v aquatic:
 EIWk15_C3_PCo1_SinPN_pGLS_aquatic<-gls(PCo1_SinPN~ Aquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C3_SinPN_bm)
  summary(EIWk15_C3_PCo1_SinPN_pGLS_aquatic)
  plot(EIWk15_C3_PCo1_SinPN_pGLS_aquatic)
  AIC(EIWk15_C3_PCo1_SinPN_pGLS_aquatic)
  R2_lik(EIWk15_C3_PCo1_SinPN_pGLS_aquatic)

# PCo1 v Thalattosuchia:
 EIWk15_C3_PCo1_SinPN_pGLS_thalattosuchia<-gls(PCo1_SinPN~ Thalattosuchian_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C3_SinPN_bm)
  summary(EIWk15_C3_PCo1_SinPN_pGLS_thalattosuchia)
  plot(EIWk15_C3_PCo1_SinPN_pGLS_thalattosuchia)
  AIC(EIWk15_C3_PCo1_SinPN_pGLS_thalattosuchia)
  R2_lik(EIWk15_C3_PCo1_SinPN_pGLS_thalattosuchia)

# PCo1 v Metriorhynchidae:
 EIWk15_C3_PCo1_SinPN_pGLS_metriorhynchidae<-gls(PCo1_SinPN~ Metriorhynchid_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C3_SinPN_bm)
  summary(EIWk15_C3_PCo1_SinPN_pGLS_metriorhynchidae)
  plot(EIWk15_C3_PCo1_SinPN_pGLS_metriorhynchidae)
  AIC(EIWk15_C3_PCo1_SinPN_pGLS_metriorhynchidae)
  R2_lik(EIWk15_C3_PCo1_SinPN_pGLS_metriorhynchidae)


### PCo2 pGLS
# PCo2 v terrestrial:
 EIWk15_C3_PCo2_SinPN_pGLS_terrestrial<-gls(PCo2_SinPN~ Terrestrial_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C3_SinPN_bm)
  summary(EIWk15_C3_PCo2_SinPN_pGLS_terrestrial)
  plot(EIWk15_C3_PCo2_SinPN_pGLS_terrestrial)
  AIC(EIWk15_C3_PCo2_SinPN_pGLS_terrestrial)
  R2_lik(EIWk15_C3_PCo2_SinPN_pGLS_terrestrial) 

# PCo2 v semiaquatic:
 EIWk15_C3_PCo2_SinPN_pGLS_semiaquatic<-gls(PCo2_SinPN~ Semiaquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C3_SinPN_bm)
  summary(EIWk15_C3_PCo2_SinPN_pGLS_semiaquatic)
  plot(EIWk15_C3_PCo2_SinPN_pGLS_semiaquatic)
  AIC(EIWk15_C3_PCo2_SinPN_pGLS_semiaquatic)
  R2_lik(EIWk15_C3_PCo2_SinPN_pGLS_semiaquatic) 

# PCo2 v pelagic:
 EIWk15_C3_PCo2_SinPN_pGLS_pelagic<-gls(PCo2_SinPN~ Pelagic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C3_SinPN_bm)
  summary(EIWk15_C3_PCo2_SinPN_pGLS_pelagic)
  plot(EIWk15_C3_PCo2_SinPN_pGLS_pelagic)
  AIC(EIWk15_C3_PCo2_SinPN_pGLS_pelagic)
  R2_lik(EIWk15_C3_PCo2_SinPN_pGLS_pelagic) 

# PCo2 v marine:
 EIWk15_C3_PCo2_SinPN_pGLS_marine<-gls(PCo2_SinPN~ Marine_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C3_SinPN_bm)
  summary(EIWk15_C3_PCo2_SinPN_pGLS_marine)
  plot(EIWk15_C3_PCo2_SinPN_pGLS_marine)
  AIC(EIWk15_C3_PCo2_SinPN_pGLS_marine)
  R2_lik(EIWk15_C3_PCo2_SinPN_pGLS_marine) 

# PCo2 v freshwater:
 EIWk15_C3_PCo2_SinPN_pGLS_freshwater<-gls(PCo2_SinPN~ Freshwater_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C3_SinPN_bm)
  summary(EIWk15_C3_PCo2_SinPN_pGLS_freshwater)
  plot(EIWk15_C3_PCo2_SinPN_pGLS_freshwater)
  AIC(EIWk15_C3_PCo2_SinPN_pGLS_freshwater)
  R2_lik(EIWk15_C3_PCo2_SinPN_pGLS_freshwater) 

# PCo2 v aquatic:
 EIWk15_C3_PCo2_SinPN_pGLS_aquatic<-gls(PCo2_SinPN~ Aquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C3_SinPN_bm)
  summary(EIWk15_C3_PCo2_SinPN_pGLS_aquatic)
  plot(EIWk15_C3_PCo2_SinPN_pGLS_aquatic)
  AIC(EIWk15_C3_PCo2_SinPN_pGLS_aquatic)
  R2_lik(EIWk15_C3_PCo2_SinPN_pGLS_aquatic) 

# PCo2 v Thalattosuchia:
 EIWk15_C3_PCo2_SinPN_pGLS_thalattosuchia<-gls(PCo2_SinPN~ Thalattosuchian_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C3_SinPN_bm)
  summary(EIWk15_C3_PCo2_SinPN_pGLS_thalattosuchia)
  plot(EIWk15_C3_PCo2_SinPN_pGLS_thalattosuchia)
  AIC(EIWk15_C3_PCo2_SinPN_pGLS_thalattosuchia)
  R2_lik(EIWk15_C3_PCo2_SinPN_pGLS_thalattosuchia) 

# PCo2 v Metriorhynchidae:
 EIWk15_C3_PCo2_SinPN_pGLS_metriorhynchidae<-gls(PCo2_SinPN~ Metriorhynchid_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C3_SinPN_bm)
  summary(EIWk15_C3_PCo2_SinPN_pGLS_metriorhynchidae)
  plot(EIWk15_C3_PCo2_SinPN_pGLS_metriorhynchidae)
  AIC(EIWk15_C3_PCo2_SinPN_pGLS_metriorhynchidae)
  R2_lik(EIWk15_C3_PCo2_SinPN_pGLS_metriorhynchidae) 


### PCo3 pGLS
# PCo3 v terrestrial:
 EIWk15_C3_PCo3_SinPN_pGLS_terrestrial<-gls(PCo3_SinPN~ Terrestrial_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C3_SinPN_bm)
  summary(EIWk15_C3_PCo3_SinPN_pGLS_terrestrial)
  plot(EIWk15_C3_PCo3_SinPN_pGLS_terrestrial)
  AIC(EIWk15_C3_PCo3_SinPN_pGLS_terrestrial)
  R2_lik(EIWk15_C3_PCo3_SinPN_pGLS_terrestrial)

# PCo3 v semiaquatic:
 EIWk15_C3_PCo3_SinPN_pGLS_semiaquatic<-gls(PCo3_SinPN~ Semiaquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C3_SinPN_bm)
  summary(EIWk15_C3_PCo3_SinPN_pGLS_semiaquatic)
  plot(EIWk15_C3_PCo3_SinPN_pGLS_semiaquatic)
  AIC(EIWk15_C3_PCo3_SinPN_pGLS_semiaquatic)
  R2_lik(EIWk15_C3_PCo3_SinPN_pGLS_semiaquatic)

# PCo3 v pelagic:
 EIWk15_C3_PCo3_SinPN_pGLS_pelagic<-gls(PCo3_SinPN~ Pelagic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C3_SinPN_bm)
  summary(EIWk15_C3_PCo3_SinPN_pGLS_pelagic)
  plot(EIWk15_C3_PCo3_SinPN_pGLS_pelagic)
  AIC(EIWk15_C3_PCo3_SinPN_pGLS_pelagic)
  R2_lik(EIWk15_C3_PCo3_SinPN_pGLS_pelagic)

# PCo3 v marine:
 EIWk15_C3_PCo3_SinPN_pGLS_marine<-gls(PCo3_SinPN~ Marine_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C3_SinPN_bm)
  summary(EIWk15_C3_PCo3_SinPN_pGLS_marine)
  plot(EIWk15_C3_PCo3_SinPN_pGLS_marine)
  AIC(EIWk15_C3_PCo3_SinPN_pGLS_marine)
  R2_lik(EIWk15_C3_PCo3_SinPN_pGLS_marine)

# PCo3 v freshwater:
 EIWk15_C3_PCo3_SinPN_pGLS_freshwater<-gls(PCo3_SinPN~ Freshwater_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C3_SinPN_bm)
  summary(EIWk15_C3_PCo3_SinPN_pGLS_freshwater)
  plot(EIWk15_C3_PCo3_SinPN_pGLS_freshwater)
  AIC(EIWk15_C3_PCo3_SinPN_pGLS_freshwater)
  R2_lik(EIWk15_C3_PCo3_SinPN_pGLS_freshwater)

# PCo3 v aquatic:
 EIWk15_C3_PCo3_SinPN_pGLS_aquatic<-gls(PCo3_SinPN~ Aquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C3_SinPN_bm)
  summary(EIWk15_C3_PCo3_SinPN_pGLS_aquatic)
  plot(EIWk15_C3_PCo3_SinPN_pGLS_aquatic)
  AIC(EIWk15_C3_PCo3_SinPN_pGLS_aquatic)
  R2_lik(EIWk15_C3_PCo3_SinPN_pGLS_aquatic)

# PCo3 v Thalattosuchia:
 EIWk15_C3_PCo3_SinPN_pGLS_thalattosuchia<-gls(PCo3_SinPN~ Thalattosuchian_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C3_SinPN_bm)
  summary(EIWk15_C3_PCo3_SinPN_pGLS_thalattosuchia)
  plot(EIWk15_C3_PCo3_SinPN_pGLS_thalattosuchia)
  AIC(EIWk15_C3_PCo3_SinPN_pGLS_thalattosuchia)
  R2_lik(EIWk15_C3_PCo3_SinPN_pGLS_thalattosuchia)

# PCo3 v Metriorhynchidae:
 EIWk15_C3_PCo3_SinPN_pGLS_metriorhynchidae<-gls(PCo3_SinPN~ Metriorhynchid_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C3_SinPN_bm)
  summary(EIWk15_C3_PCo3_SinPN_pGLS_metriorhynchidae)
  plot(EIWk15_C3_PCo3_SinPN_pGLS_metriorhynchidae)
  AIC(EIWk15_C3_PCo3_SinPN_pGLS_metriorhynchidae)
  R2_lik(EIWk15_C3_PCo3_SinPN_pGLS_metriorhynchidae)

txtStop()


# ------------------------------------------------------ #
# 12.9 Extended implied weights (k=15) C4 paranasal GLMs #
# ------------------------------------------------------ #

### Writes to file:
txtStart("output/ParanasalAnalyses/EIWk15_C4/Files/EIWk15_C4_paranasal_pGLS.txt")


EIWk15_C4_SinPN_bm<-corBrownian(1,EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, form= ~Species_Sin)


### PCo1 pGLS
# PCo1 v terrestrial:
 EIWk15_C4_PCo1_SinPN_pGLS_Terrestrial<-gls(PCo1_SinPN~ Terrestrial_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C4_SinPN_bm)
  summary(EIWk15_C4_PCo1_SinPN_pGLS_Terrestrial)
  plot(EIWk15_C4_PCo1_SinPN_pGLS_Terrestrial)
  AIC(EIWk15_C4_PCo1_SinPN_pGLS_Terrestrial)
  R2_lik(EIWk15_C4_PCo1_SinPN_pGLS_Terrestrial)

# PCo1 v semiaquatic:
 EIWk15_C4_PCo1_SinPN_pGLS_semiaquatic<-gls(PCo1_SinPN~ Semiaquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C4_SinPN_bm)
  summary(EIWk15_C4_PCo1_SinPN_pGLS_semiaquatic)
  plot(EIWk15_C4_PCo1_SinPN_pGLS_semiaquatic)
  AIC(EIWk15_C4_PCo1_SinPN_pGLS_semiaquatic)
  R2_lik(EIWk15_C4_PCo1_SinPN_pGLS_semiaquatic)

# PCo1 v pelagic:
 EIWk15_C4_PCo1_SinPN_pGLS_pelagic<-gls(PCo1_SinPN~ Pelagic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C4_SinPN_bm)
  summary(EIWk15_C4_PCo1_SinPN_pGLS_pelagic)
  plot(EIWk15_C4_PCo1_SinPN_pGLS_pelagic)
  AIC(EIWk15_C4_PCo1_SinPN_pGLS_pelagic)
  R2_lik(EIWk15_C4_PCo1_SinPN_pGLS_pelagic)

# PCo1 v marine:
 EIWk15_C4_PCo1_SinPN_pGLS_marine<-gls(PCo1_SinPN~ Marine_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C4_SinPN_bm)
  summary(EIWk15_C4_PCo1_SinPN_pGLS_marine)
  plot(EIWk15_C4_PCo1_SinPN_pGLS_marine)
  AIC(EIWk15_C4_PCo1_SinPN_pGLS_marine)
  R2_lik(EIWk15_C4_PCo1_SinPN_pGLS_marine)

# PCo1 v freshwater:
 EIWk15_C4_PCo1_SinPN_pGLS_freshwater<-gls(PCo1_SinPN~ Freshwater_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C4_SinPN_bm)
  summary(EIWk15_C4_PCo1_SinPN_pGLS_freshwater)
  plot(EIWk15_C4_PCo1_SinPN_pGLS_freshwater)
  AIC(EIWk15_C4_PCo1_SinPN_pGLS_freshwater)
  R2_lik(EIWk15_C4_PCo1_SinPN_pGLS_freshwater)

# PCo1 v aquatic:
 EIWk15_C4_PCo1_SinPN_pGLS_aquatic<-gls(PCo1_SinPN~ Aquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C4_SinPN_bm)
  summary(EIWk15_C4_PCo1_SinPN_pGLS_aquatic)
  plot(EIWk15_C4_PCo1_SinPN_pGLS_aquatic)
  AIC(EIWk15_C4_PCo1_SinPN_pGLS_aquatic)
  R2_lik(EIWk15_C4_PCo1_SinPN_pGLS_aquatic)

# PCo1 v Thalattosuchia:
 EIWk15_C4_PCo1_SinPN_pGLS_thalattosuchia<-gls(PCo1_SinPN~ Thalattosuchian_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C4_SinPN_bm)
  summary(EIWk15_C4_PCo1_SinPN_pGLS_thalattosuchia)
  plot(EIWk15_C4_PCo1_SinPN_pGLS_thalattosuchia)
  AIC(EIWk15_C4_PCo1_SinPN_pGLS_thalattosuchia)
  R2_lik(EIWk15_C4_PCo1_SinPN_pGLS_thalattosuchia)

# PCo1 v Metriorhynchidae:
 EIWk15_C4_PCo1_SinPN_pGLS_metriorhynchidae<-gls(PCo1_SinPN~ Metriorhynchid_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C4_SinPN_bm)
  summary(EIWk15_C4_PCo1_SinPN_pGLS_metriorhynchidae)
  plot(EIWk15_C4_PCo1_SinPN_pGLS_metriorhynchidae)
  AIC(EIWk15_C4_PCo1_SinPN_pGLS_metriorhynchidae)
  R2_lik(EIWk15_C4_PCo1_SinPN_pGLS_metriorhynchidae)


### PCo2 pGLS
# PCo2 v terrestrial:
 EIWk15_C4_PCo2_SinPN_pGLS_terrestrial<-gls(PCo2_SinPN~ Terrestrial_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C4_SinPN_bm)
  summary(EIWk15_C4_PCo2_SinPN_pGLS_terrestrial)
  plot(EIWk15_C4_PCo2_SinPN_pGLS_terrestrial)
  AIC(EIWk15_C4_PCo2_SinPN_pGLS_terrestrial)
  R2_lik(EIWk15_C4_PCo2_SinPN_pGLS_terrestrial) 

# PCo2 v semiaquatic:
 EIWk15_C4_PCo2_SinPN_pGLS_semiaquatic<-gls(PCo2_SinPN~ Semiaquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C4_SinPN_bm)
  summary(EIWk15_C4_PCo2_SinPN_pGLS_semiaquatic)
  plot(EIWk15_C4_PCo2_SinPN_pGLS_semiaquatic)
  AIC(EIWk15_C4_PCo2_SinPN_pGLS_semiaquatic)
  R2_lik(EIWk15_C4_PCo2_SinPN_pGLS_semiaquatic) 

# PCo2 v pelagic:
 EIWk15_C4_PCo2_SinPN_pGLS_pelagic<-gls(PCo2_SinPN~ Pelagic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C4_SinPN_bm)
  summary(EIWk15_C4_PCo2_SinPN_pGLS_pelagic)
  plot(EIWk15_C4_PCo2_SinPN_pGLS_pelagic)
  AIC(EIWk15_C4_PCo2_SinPN_pGLS_pelagic)
  R2_lik(EIWk15_C4_PCo2_SinPN_pGLS_pelagic) 

# PCo2 v marine:
 EIWk15_C4_PCo2_SinPN_pGLS_marine<-gls(PCo2_SinPN~ Marine_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C4_SinPN_bm)
  summary(EIWk15_C4_PCo2_SinPN_pGLS_marine)
  plot(EIWk15_C4_PCo2_SinPN_pGLS_marine)
  AIC(EIWk15_C4_PCo2_SinPN_pGLS_marine)
  R2_lik(EIWk15_C4_PCo2_SinPN_pGLS_marine) 

# PCo2 v freshwater:
 EIWk15_C4_PCo2_SinPN_pGLS_freshwater<-gls(PCo2_SinPN~ Freshwater_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C4_SinPN_bm)
  summary(EIWk15_C4_PCo2_SinPN_pGLS_freshwater)
  plot(EIWk15_C4_PCo2_SinPN_pGLS_freshwater)
  AIC(EIWk15_C4_PCo2_SinPN_pGLS_freshwater)
  R2_lik(EIWk15_C4_PCo2_SinPN_pGLS_freshwater) 

# PCo2 v aquatic:
 EIWk15_C4_PCo2_SinPN_pGLS_aquatic<-gls(PCo2_SinPN~ Aquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C4_SinPN_bm)
  summary(EIWk15_C4_PCo2_SinPN_pGLS_aquatic)
  plot(EIWk15_C4_PCo2_SinPN_pGLS_aquatic)
  AIC(EIWk15_C4_PCo2_SinPN_pGLS_aquatic)
  R2_lik(EIWk15_C4_PCo2_SinPN_pGLS_aquatic) 

# PCo2 v Thalattosuchia:
 EIWk15_C4_PCo2_SinPN_pGLS_thalattosuchia<-gls(PCo2_SinPN~ Thalattosuchian_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C4_SinPN_bm)
  summary(EIWk15_C4_PCo2_SinPN_pGLS_thalattosuchia)
  plot(EIWk15_C4_PCo2_SinPN_pGLS_thalattosuchia)
  AIC(EIWk15_C4_PCo2_SinPN_pGLS_thalattosuchia)
  R2_lik(EIWk15_C4_PCo2_SinPN_pGLS_thalattosuchia) 

# PCo2 v Metriorhynchidae:
 EIWk15_C4_PCo2_SinPN_pGLS_metriorhynchidae<-gls(PCo2_SinPN~ Metriorhynchid_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C4_SinPN_bm)
  summary(EIWk15_C4_PCo2_SinPN_pGLS_metriorhynchidae)
  plot(EIWk15_C4_PCo2_SinPN_pGLS_metriorhynchidae)
  AIC(EIWk15_C4_PCo2_SinPN_pGLS_metriorhynchidae)
  R2_lik(EIWk15_C4_PCo2_SinPN_pGLS_metriorhynchidae) 


### PCo3 pGLS
# PCo3 v terrestrial:
 EIWk15_C4_PCo3_SinPN_pGLS_terrestrial<-gls(PCo3_SinPN~ Terrestrial_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C4_SinPN_bm)
  summary(EIWk15_C4_PCo3_SinPN_pGLS_terrestrial)
  plot(EIWk15_C4_PCo3_SinPN_pGLS_terrestrial)
  AIC(EIWk15_C4_PCo3_SinPN_pGLS_terrestrial)
  R2_lik(EIWk15_C4_PCo3_SinPN_pGLS_terrestrial)

# PCo3 v semiaquatic:
 EIWk15_C4_PCo3_SinPN_pGLS_semiaquatic<-gls(PCo3_SinPN~ Semiaquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C4_SinPN_bm)
  summary(EIWk15_C4_PCo3_SinPN_pGLS_semiaquatic)
  plot(EIWk15_C4_PCo3_SinPN_pGLS_semiaquatic)
  AIC(EIWk15_C4_PCo3_SinPN_pGLS_semiaquatic)
  R2_lik(EIWk15_C4_PCo3_SinPN_pGLS_semiaquatic)

# PCo3 v pelagic:
 EIWk15_C4_PCo3_SinPN_pGLS_pelagic<-gls(PCo3_SinPN~ Pelagic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C4_SinPN_bm)
  summary(EIWk15_C4_PCo3_SinPN_pGLS_pelagic)
  plot(EIWk15_C4_PCo3_SinPN_pGLS_pelagic)
  AIC(EIWk15_C4_PCo3_SinPN_pGLS_pelagic)
  R2_lik(EIWk15_C4_PCo3_SinPN_pGLS_pelagic)

# PCo3 v marine:
 EIWk15_C4_PCo3_SinPN_pGLS_marine<-gls(PCo3_SinPN~ Marine_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C4_SinPN_bm)
  summary(EIWk15_C4_PCo3_SinPN_pGLS_marine)
  plot(EIWk15_C4_PCo3_SinPN_pGLS_marine)
  AIC(EIWk15_C4_PCo3_SinPN_pGLS_marine)
  R2_lik(EIWk15_C4_PCo3_SinPN_pGLS_marine)

# PCo3 v freshwater:
 EIWk15_C4_PCo3_SinPN_pGLS_freshwater<-gls(PCo3_SinPN~ Freshwater_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C4_SinPN_bm)
  summary(EIWk15_C4_PCo3_SinPN_pGLS_freshwater)
  plot(EIWk15_C4_PCo3_SinPN_pGLS_freshwater)
  AIC(EIWk15_C4_PCo3_SinPN_pGLS_freshwater)
  R2_lik(EIWk15_C4_PCo3_SinPN_pGLS_freshwater)

# PCo3 v aquatic:
 EIWk15_C4_PCo3_SinPN_pGLS_aquatic<-gls(PCo3_SinPN~ Aquatic_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C4_SinPN_bm)
  summary(EIWk15_C4_PCo3_SinPN_pGLS_aquatic)
  plot(EIWk15_C4_PCo3_SinPN_pGLS_aquatic)
  AIC(EIWk15_C4_PCo3_SinPN_pGLS_aquatic)
  R2_lik(EIWk15_C4_PCo3_SinPN_pGLS_aquatic)

# PCo3 v Thalattosuchia:
 EIWk15_C4_PCo3_SinPN_pGLS_thalattosuchia<-gls(PCo3_SinPN~ Thalattosuchian_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C4_SinPN_bm)
  summary(EIWk15_C4_PCo3_SinPN_pGLS_thalattosuchia)
  plot(EIWk15_C4_PCo3_SinPN_pGLS_thalattosuchia)
  AIC(EIWk15_C4_PCo3_SinPN_pGLS_thalattosuchia)
  R2_lik(EIWk15_C4_PCo3_SinPN_pGLS_thalattosuchia)

# PCo3 v Metriorhynchidae:
 EIWk15_C4_PCo3_SinPN_pGLS_metriorhynchidae<-gls(PCo3_SinPN~ Metriorhynchid_Sin, method= "ML", data= CrocSinPNData, correlation = EIWk15_C4_SinPN_bm)
  summary(EIWk15_C4_PCo3_SinPN_pGLS_metriorhynchidae)
  plot(EIWk15_C4_PCo3_SinPN_pGLS_metriorhynchidae)
  AIC(EIWk15_C4_PCo3_SinPN_pGLS_metriorhynchidae)
  R2_lik(EIWk15_C4_PCo3_SinPN_pGLS_metriorhynchidae)

txtStop()

 beep(1)


# ========================== #
# 13. Model fitting analyses #
# ========================== #

### Required libraries:
 library(beepr)
 library(geiger)
 library(TeachingDemos)

clusterCall(cluster, function () {
 library(geiger)
})


# ------------------------------ #
# 13.1 Equal weights C0 analysis #
# ------------------------------ #

### Writes to file:
txtStart("output/ParanasalAnalyses/EW_C0/Files/EW_C0_paranasal_fittingmodel.txt")

 EW_C0_paranasal_fitBM     <- fitContinuous(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="BM")
 EW_C0_paranasal_fitOU     <- fitContinuous(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="OU")
 EW_C0_paranasal_fitEB     <- fitContinuous(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="EB")
 EW_C0_paranasal_fittrend  <- fitContinuous(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="trend")
 EW_C0_paranasal_fitlambda <- fitContinuous(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="lambda")
 EW_C0_paranasal_fitkappa  <- fitContinuous(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="kappa")
 EW_C0_paranasal_fitdelta  <- fitContinuous(EW_C0_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="delta")

# Axis 1:
 EW_C0_paranasal_aic_vals1<-setNames(c(EW_C0_paranasal_fitBM$Axis.1$opt$aicc,EW_C0_paranasal_fitOU$Axis.1$opt$aicc,EW_C0_paranasal_fitEB$Axis.1$opt$aicc,EW_C0_paranasal_fittrend$Axis.1$opt$aicc,EW_C0_paranasal_fitlambda$Axis.1$opt$aicc,EW_C0_paranasal_fitkappa$Axis.1$opt$aicc,EW_C0_paranasal_fitdelta$Axis.1$opt$aicc),
                   c("BM","OU","EB","trend","lambda","kappa","delta"))
 EW_C0_paranasal_aic_vals1

# Axis 2:
 EW_C0_paranasal_aic_vals2<-setNames(c(EW_C0_paranasal_fitBM$Axis.2$opt$aicc,EW_C0_paranasal_fitOU$Axis.2$opt$aicc,EW_C0_paranasal_fitEB$Axis.2$opt$aicc,EW_C0_paranasal_fittrend$Axis.2$opt$aicc,EW_C0_paranasal_fitlambda$Axis.2$opt$aicc,EW_C0_paranasal_fitkappa$Axis.2$opt$aicc,EW_C0_paranasal_fitdelta$Axis.2$opt$aicc),
                   c("BM","OU","EB","trend","lambda","kappa","delta"))
 EW_C0_paranasal_aic_vals2

# Axis 3:
 EW_C0_paranasal_aic_vals3<-setNames(c(EW_C0_paranasal_fitBM$Axis.3$opt$aicc,EW_C0_paranasal_fitOU$Axis.3$opt$aicc,EW_C0_paranasal_fitEB$Axis.3$opt$aicc,EW_C0_paranasal_fittrend$Axis.3$opt$aicc,EW_C0_paranasal_fitlambda$Axis.3$opt$aicc,EW_C0_paranasal_fitkappa$Axis.3$opt$aicc,EW_C0_paranasal_fitdelta$Axis.3$opt$aicc),
                   c("BM","OU","EB","trend","lambda","kappa","delta"))
 EW_C0_paranasal_aic_vals3

txtStop()


# ------------------------------ #
# 13.2 Equal weights C1 analysis #
# ------------------------------ #

### Writes to file:
txtStart("output/ParanasalAnalyses/EW_C1/Files/EW_C1_paranasal_fittingmodel.txt")

 EW_C1_paranasal_fitBM     <- fitContinuous(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="BM")
 EW_C1_paranasal_fitOU     <- fitContinuous(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="OU")
 EW_C1_paranasal_fitEB     <- fitContinuous(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="EB")
 EW_C1_paranasal_fittrend  <- fitContinuous(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="trend")
 EW_C1_paranasal_fitlambda <- fitContinuous(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="lambda")
 EW_C1_paranasal_fitkappa  <- fitContinuous(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="kappa")
 EW_C1_paranasal_fitdelta  <- fitContinuous(EW_C1_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="delta")

# Axis 1:
 EW_C1_paranasal_aic_vals1<-setNames(c(EW_C1_paranasal_fitBM$Axis.1$opt$aicc,EW_C1_paranasal_fitOU$Axis.1$opt$aicc,EW_C1_paranasal_fitEB$Axis.1$opt$aicc,EW_C1_paranasal_fittrend$Axis.1$opt$aicc,EW_C1_paranasal_fitlambda$Axis.1$opt$aicc,EW_C1_paranasal_fitkappa$Axis.1$opt$aicc,EW_C1_paranasal_fitdelta$Axis.1$opt$aicc),
                   c("BM","OU","EB","trend","lambda","kappa","delta"))
 EW_C1_paranasal_aic_vals1

# Axis 2:
 EW_C1_paranasal_aic_vals2<-setNames(c(EW_C1_paranasal_fitBM$Axis.2$opt$aicc,EW_C1_paranasal_fitOU$Axis.2$opt$aicc,EW_C1_paranasal_fitEB$Axis.2$opt$aicc,EW_C1_paranasal_fittrend$Axis.2$opt$aicc,EW_C1_paranasal_fitlambda$Axis.2$opt$aicc,EW_C1_paranasal_fitkappa$Axis.2$opt$aicc,EW_C1_paranasal_fitdelta$Axis.2$opt$aicc),
                   c("BM","OU","EB","trend","lambda","kappa","delta"))
 EW_C1_paranasal_aic_vals2

# Axis 3:
 EW_C1_paranasal_aic_vals3<-setNames(c(EW_C1_paranasal_fitBM$Axis.3$opt$aicc,EW_C1_paranasal_fitOU$Axis.3$opt$aicc,EW_C1_paranasal_fitEB$Axis.3$opt$aicc,EW_C1_paranasal_fittrend$Axis.3$opt$aicc,EW_C1_paranasal_fitlambda$Axis.3$opt$aicc,EW_C1_paranasal_fitkappa$Axis.3$opt$aicc,EW_C1_paranasal_fitdelta$Axis.3$opt$aicc),
                   c("BM","OU","EB","trend","lambda","kappa","delta"))
 EW_C1_paranasal_aic_vals3

txtStop()


# ------------------------------ #
# 13.3 Equal weights C3 analysis #
# ------------------------------ #

### Writes to file:
txtStart("output/ParanasalAnalyses/EW_C3/Files/EW_C3_paranasal_fittingmodel.txt")

 EW_C3_paranasal_fitBM     <- fitContinuous(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="BM")
 EW_C3_paranasal_fitOU     <- fitContinuous(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="OU")
 EW_C3_paranasal_fitEB     <- fitContinuous(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="EB")
 EW_C3_paranasal_fittrend  <- fitContinuous(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="trend")
 EW_C3_paranasal_fitlambda <- fitContinuous(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="lambda")
 EW_C3_paranasal_fitkappa  <- fitContinuous(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="kappa")
 EW_C3_paranasal_fitdelta  <- fitContinuous(EW_C3_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="delta")

# Axis 1:
 EW_C3_paranasal_aic_vals1<-setNames(c(EW_C3_paranasal_fitBM$Axis.1$opt$aicc,EW_C3_paranasal_fitOU$Axis.1$opt$aicc,EW_C3_paranasal_fitEB$Axis.1$opt$aicc,EW_C3_paranasal_fittrend$Axis.1$opt$aicc,EW_C3_paranasal_fitlambda$Axis.1$opt$aicc,EW_C3_paranasal_fitkappa$Axis.1$opt$aicc,EW_C3_paranasal_fitdelta$Axis.1$opt$aicc),
                   c("BM","OU","EB","trend","lambda","kappa","delta"))
 EW_C3_paranasal_aic_vals1

# Axis 2:
 EW_C3_paranasal_aic_vals2<-setNames(c(EW_C3_paranasal_fitBM$Axis.2$opt$aicc,EW_C3_paranasal_fitOU$Axis.2$opt$aicc,EW_C3_paranasal_fitEB$Axis.2$opt$aicc,EW_C3_paranasal_fittrend$Axis.2$opt$aicc,EW_C3_paranasal_fitlambda$Axis.2$opt$aicc,EW_C3_paranasal_fitkappa$Axis.2$opt$aicc,EW_C3_paranasal_fitdelta$Axis.2$opt$aicc),
                   c("BM","OU","EB","trend","lambda","kappa","delta"))
 EW_C3_paranasal_aic_vals2

# Axis 3:
 EW_C3_paranasal_aic_vals3<-setNames(c(EW_C3_paranasal_fitBM$Axis.3$opt$aicc,EW_C3_paranasal_fitOU$Axis.3$opt$aicc,EW_C3_paranasal_fitEB$Axis.3$opt$aicc,EW_C3_paranasal_fittrend$Axis.3$opt$aicc,EW_C3_paranasal_fitlambda$Axis.3$opt$aicc,EW_C3_paranasal_fitkappa$Axis.3$opt$aicc,EW_C3_paranasal_fitdelta$Axis.3$opt$aicc),
                   c("BM","OU","EB","trend","lambda","kappa","delta"))
 EW_C3_paranasal_aic_vals3

txtStop()


# ------------------------------ #
# 13.4 Equal weights C4 analysis #
# ------------------------------ #

### Writes to file:
txtStart("output/ParanasalAnalyses/EW_C4/Files/EW_C4_paranasal_fittingmodel.txt")

 EW_C4_paranasal_fitBM     <- fitContinuous(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="BM")
 EW_C4_paranasal_fitOU     <- fitContinuous(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="OU")
 EW_C4_paranasal_fitEB     <- fitContinuous(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="EB")
 EW_C4_paranasal_fittrend  <- fitContinuous(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="trend")
 EW_C4_paranasal_fitlambda <- fitContinuous(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="lambda")
 EW_C4_paranasal_fitkappa  <- fitContinuous(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="kappa")
 EW_C4_paranasal_fitdelta  <- fitContinuous(EW_C4_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="delta")

# Axis 1:
 EW_C4_paranasal_aic_vals1<-setNames(c(EW_C4_paranasal_fitBM$Axis.1$opt$aicc,EW_C4_paranasal_fitOU$Axis.1$opt$aicc,EW_C4_paranasal_fitEB$Axis.1$opt$aicc,EW_C4_paranasal_fittrend$Axis.1$opt$aicc,EW_C4_paranasal_fitlambda$Axis.1$opt$aicc,EW_C4_paranasal_fitkappa$Axis.1$opt$aicc,EW_C4_paranasal_fitdelta$Axis.1$opt$aicc),
                   c("BM","OU","EB","trend","lambda","kappa","delta"))
 EW_C4_paranasal_aic_vals1

# Axis 2:
 EW_C4_paranasal_aic_vals2<-setNames(c(EW_C4_paranasal_fitBM$Axis.2$opt$aicc,EW_C4_paranasal_fitOU$Axis.2$opt$aicc,EW_C4_paranasal_fitEB$Axis.2$opt$aicc,EW_C4_paranasal_fittrend$Axis.2$opt$aicc,EW_C4_paranasal_fitlambda$Axis.2$opt$aicc,EW_C4_paranasal_fitkappa$Axis.2$opt$aicc,EW_C4_paranasal_fitdelta$Axis.2$opt$aicc),
                   c("BM","OU","EB","trend","lambda","kappa","delta"))
 EW_C4_paranasal_aic_vals2

# Axis 3:
 EW_C4_paranasal_aic_vals3<-setNames(c(EW_C4_paranasal_fitBM$Axis.3$opt$aicc,EW_C4_paranasal_fitOU$Axis.3$opt$aicc,EW_C4_paranasal_fitEB$Axis.3$opt$aicc,EW_C4_paranasal_fittrend$Axis.3$opt$aicc,EW_C4_paranasal_fitlambda$Axis.3$opt$aicc,EW_C4_paranasal_fitkappa$Axis.3$opt$aicc,EW_C4_paranasal_fitdelta$Axis.3$opt$aicc),
                   c("BM","OU","EB","trend","lambda","kappa","delta"))
 EW_C4_paranasal_aic_vals3

txtStop()


# ------------------------------------------------ #
# 13.5 Extended implied weights (k=15) C0 analysis #
# ------------------------------------------------ #

### Writes to file:
txtStart("output/ParanasalAnalyses/EIWk15_C0/Files/EIWk15_C0_paranasal_fittingmodel.txt")

 EIWk15_C0_paranasal_fitBM     <- fitContinuous(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="BM")
 EIWk15_C0_paranasal_fitOU     <- fitContinuous(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="OU")
 EIWk15_C0_paranasal_fitEB     <- fitContinuous(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="EB")
 EIWk15_C0_paranasal_fittrend  <- fitContinuous(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="trend")
 EIWk15_C0_paranasal_fitlambda <- fitContinuous(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="lambda")
 EIWk15_C0_paranasal_fitkappa  <- fitContinuous(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="kappa")
 EIWk15_C0_paranasal_fitdelta  <- fitContinuous(EIWk15_C0_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="delta")

# Axis 1:
 EIWk15_C0_paranasal_aic_vals1<-setNames(c(EIWk15_C0_paranasal_fitBM$Axis.1$opt$aicc,EIWk15_C0_paranasal_fitOU$Axis.1$opt$aicc,EIWk15_C0_paranasal_fitEB$Axis.1$opt$aicc,EIWk15_C0_paranasal_fittrend$Axis.1$opt$aicc,EIWk15_C0_paranasal_fitlambda$Axis.1$opt$aicc,EIWk15_C0_paranasal_fitkappa$Axis.1$opt$aicc,EIWk15_C0_paranasal_fitdelta$Axis.1$opt$aicc),
                   c("BM","OU","EB","trend","lambda","kappa","delta"))
 EIWk15_C0_paranasal_aic_vals1

# Axis 2:
 EIWk15_C0_paranasal_aic_vals2<-setNames(c(EIWk15_C0_paranasal_fitBM$Axis.2$opt$aicc,EIWk15_C0_paranasal_fitOU$Axis.2$opt$aicc,EIWk15_C0_paranasal_fitEB$Axis.2$opt$aicc,EIWk15_C0_paranasal_fittrend$Axis.2$opt$aicc,EIWk15_C0_paranasal_fitlambda$Axis.2$opt$aicc,EIWk15_C0_paranasal_fitkappa$Axis.2$opt$aicc,EIWk15_C0_paranasal_fitdelta$Axis.2$opt$aicc),
                   c("BM","OU","EB","trend","lambda","kappa","delta"))
 EIWk15_C0_paranasal_aic_vals2

# Axis 3:
 EIWk15_C0_paranasal_aic_vals3<-setNames(c(EIWk15_C0_paranasal_fitBM$Axis.3$opt$aicc,EIWk15_C0_paranasal_fitOU$Axis.3$opt$aicc,EIWk15_C0_paranasal_fitEB$Axis.3$opt$aicc,EIWk15_C0_paranasal_fittrend$Axis.3$opt$aicc,EIWk15_C0_paranasal_fitlambda$Axis.3$opt$aicc,EIWk15_C0_paranasal_fitkappa$Axis.3$opt$aicc,EIWk15_C0_paranasal_fitdelta$Axis.3$opt$aicc),
                   c("BM","OU","EB","trend","lambda","kappa","delta"))
 EIWk15_C0_paranasal_aic_vals3

txtStop()


# ------------------------------------------------ #
# 13.6 Extended implied weights (k=15) C1 analysis #
# ------------------------------------------------ #

### Writes to file:
txtStart("output/ParanasalAnalyses/EIWk15_C1/Files/EIWk15_C1_paranasal_fittingmodel.txt")

 EIWk15_C1_paranasal_fitBM     <- fitContinuous(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="BM")
 EIWk15_C1_paranasal_fitOU     <- fitContinuous(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="OU")
 EIWk15_C1_paranasal_fitEB     <- fitContinuous(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="EB")
 EIWk15_C1_paranasal_fittrend  <- fitContinuous(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="trend")
 EIWk15_C1_paranasal_fitlambda <- fitContinuous(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="lambda")
 EIWk15_C1_paranasal_fitkappa  <- fitContinuous(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="kappa")
 EIWk15_C1_paranasal_fitdelta  <- fitContinuous(EIWk15_C1_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="delta")

# Axis 1:
 EIWk15_C1_paranasal_aic_vals1<-setNames(c(EIWk15_C1_paranasal_fitBM$Axis.1$opt$aicc,EIWk15_C1_paranasal_fitOU$Axis.1$opt$aicc,EIWk15_C1_paranasal_fitEB$Axis.1$opt$aicc,EIWk15_C1_paranasal_fittrend$Axis.1$opt$aicc,EIWk15_C1_paranasal_fitlambda$Axis.1$opt$aicc,EIWk15_C1_paranasal_fitkappa$Axis.1$opt$aicc,EIWk15_C1_paranasal_fitdelta$Axis.1$opt$aicc),
                   c("BM","OU","EB","trend","lambda","kappa","delta"))
 EIWk15_C1_paranasal_aic_vals1

# Axis 2:
 EIWk15_C1_paranasal_aic_vals2<-setNames(c(EIWk15_C1_paranasal_fitBM$Axis.2$opt$aicc,EIWk15_C1_paranasal_fitOU$Axis.2$opt$aicc,EIWk15_C1_paranasal_fitEB$Axis.2$opt$aicc,EIWk15_C1_paranasal_fittrend$Axis.2$opt$aicc,EIWk15_C1_paranasal_fitlambda$Axis.2$opt$aicc,EIWk15_C1_paranasal_fitkappa$Axis.2$opt$aicc,EIWk15_C1_paranasal_fitdelta$Axis.2$opt$aicc),
                   c("BM","OU","EB","trend","lambda","kappa","delta"))
 EIWk15_C1_paranasal_aic_vals2

# Axis 3:
 EIWk15_C1_paranasal_aic_vals3<-setNames(c(EIWk15_C1_paranasal_fitBM$Axis.3$opt$aicc,EIWk15_C1_paranasal_fitOU$Axis.3$opt$aicc,EIWk15_C1_paranasal_fitEB$Axis.3$opt$aicc,EIWk15_C1_paranasal_fittrend$Axis.3$opt$aicc,EIWk15_C1_paranasal_fitlambda$Axis.3$opt$aicc,EIWk15_C1_paranasal_fitkappa$Axis.3$opt$aicc,EIWk15_C1_paranasal_fitdelta$Axis.3$opt$aicc),
                   c("BM","OU","EB","trend","lambda","kappa","delta"))
 EIWk15_C1_paranasal_aic_vals3

txtStop()


# ------------------------------------------------ #
# 13.7 Extended implied weights (k=15) C3 analysis #
# ------------------------------------------------ #

### Writes to file:
txtStart("output/ParanasalAnalyses/EIWk15_C3/Files/EIWk15_C3_paranasal_fittingmodel.txt")

 EIWk15_C3_paranasal_fitBM     <- fitContinuous(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="BM")
 EIWk15_C3_paranasal_fitOU     <- fitContinuous(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="OU")
 EIWk15_C3_paranasal_fitEB     <- fitContinuous(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="EB")
 EIWk15_C3_paranasal_fittrend  <- fitContinuous(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="trend")
 EIWk15_C3_paranasal_fitlambda <- fitContinuous(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="lambda")
 EIWk15_C3_paranasal_fitkappa  <- fitContinuous(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="kappa")
 EIWk15_C3_paranasal_fitdelta  <- fitContinuous(EIWk15_C3_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="delta")

# Axis 1:
 EIWk15_C3_paranasal_aic_vals1<-setNames(c(EIWk15_C3_paranasal_fitBM$Axis.1$opt$aicc,EIWk15_C3_paranasal_fitOU$Axis.1$opt$aicc,EIWk15_C3_paranasal_fitEB$Axis.1$opt$aicc,EIWk15_C3_paranasal_fittrend$Axis.1$opt$aicc,EIWk15_C3_paranasal_fitlambda$Axis.1$opt$aicc,EIWk15_C3_paranasal_fitkappa$Axis.1$opt$aicc,EIWk15_C3_paranasal_fitdelta$Axis.1$opt$aicc),
                   c("BM","OU","EB","trend","lambda","kappa","delta"))
 EIWk15_C3_paranasal_aic_vals1

# Axis 2:
 EIWk15_C3_paranasal_aic_vals2<-setNames(c(EIWk15_C3_paranasal_fitBM$Axis.2$opt$aicc,EIWk15_C3_paranasal_fitOU$Axis.2$opt$aicc,EIWk15_C3_paranasal_fitEB$Axis.2$opt$aicc,EIWk15_C3_paranasal_fittrend$Axis.2$opt$aicc,EIWk15_C3_paranasal_fitlambda$Axis.2$opt$aicc,EIWk15_C3_paranasal_fitkappa$Axis.2$opt$aicc,EIWk15_C3_paranasal_fitdelta$Axis.2$opt$aicc),
                   c("BM","OU","EB","trend","lambda","kappa","delta"))
 EIWk15_C3_paranasal_aic_vals2

# Axis 3:
 EIWk15_C3_paranasal_aic_vals3<-setNames(c(EIWk15_C3_paranasal_fitBM$Axis.3$opt$aicc,EIWk15_C3_paranasal_fitOU$Axis.3$opt$aicc,EIWk15_C3_paranasal_fitEB$Axis.3$opt$aicc,EIWk15_C3_paranasal_fittrend$Axis.3$opt$aicc,EIWk15_C3_paranasal_fitlambda$Axis.3$opt$aicc,EIWk15_C3_paranasal_fitkappa$Axis.3$opt$aicc,EIWk15_C3_paranasal_fitdelta$Axis.3$opt$aicc),
                   c("BM","OU","EB","trend","lambda","kappa","delta"))
 EIWk15_C3_paranasal_aic_vals3

txtStop()


# ------------------------------------------------ #
# 13.8 Extended implied weights (k=15) C4 analysis #
# ------------------------------------------------ #

### Writes to file:
txtStart("output/ParanasalAnalyses/EIWk15_C4/Files/EIWk15_C4_paranasal_fittingmodel.txt")

 EIWk15_C4_paranasal_fitBM     <- fitContinuous(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="BM")
 EIWk15_C4_paranasal_fitOU     <- fitContinuous(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="OU")
 EIWk15_C4_paranasal_fitEB     <- fitContinuous(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="EB")
 EIWk15_C4_paranasal_fittrend  <- fitContinuous(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="trend")
 EIWk15_C4_paranasal_fitlambda <- fitContinuous(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="lambda")
 EIWk15_C4_paranasal_fitkappa  <- fitContinuous(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="kappa")
 EIWk15_C4_paranasal_fitdelta  <- fitContinuous(EIWk15_C4_phylogeny_SinPN_Pruned_EqualMethod, PCos_SinPN, model="delta")

# Axis 1:
 EIWk15_C4_paranasal_aic_vals1<-setNames(c(EIWk15_C4_paranasal_fitBM$Axis.1$opt$aicc,EIWk15_C4_paranasal_fitOU$Axis.1$opt$aicc,EIWk15_C4_paranasal_fitEB$Axis.1$opt$aicc,EIWk15_C4_paranasal_fittrend$Axis.1$opt$aicc,EIWk15_C4_paranasal_fitlambda$Axis.1$opt$aicc,EIWk15_C4_paranasal_fitkappa$Axis.1$opt$aicc,EIWk15_C4_paranasal_fitdelta$Axis.1$opt$aicc),
                   c("BM","OU","EB","trend","lambda","kappa","delta"))
 EIWk15_C4_paranasal_aic_vals1

# Axis 2:
 EIWk15_C4_paranasal_aic_vals2<-setNames(c(EIWk15_C4_paranasal_fitBM$Axis.2$opt$aicc,EIWk15_C4_paranasal_fitOU$Axis.2$opt$aicc,EIWk15_C4_paranasal_fitEB$Axis.2$opt$aicc,EIWk15_C4_paranasal_fittrend$Axis.2$opt$aicc,EIWk15_C4_paranasal_fitlambda$Axis.2$opt$aicc,EIWk15_C4_paranasal_fitkappa$Axis.2$opt$aicc,EIWk15_C4_paranasal_fitdelta$Axis.2$opt$aicc),
                   c("BM","OU","EB","trend","lambda","kappa","delta"))
 EIWk15_C4_paranasal_aic_vals2

# Axis 3:
 EIWk15_C4_paranasal_aic_vals3<-setNames(c(EIWk15_C4_paranasal_fitBM$Axis.3$opt$aicc,EIWk15_C4_paranasal_fitOU$Axis.3$opt$aicc,EIWk15_C4_paranasal_fitEB$Axis.3$opt$aicc,EIWk15_C4_paranasal_fittrend$Axis.3$opt$aicc,EIWk15_C4_paranasal_fitlambda$Axis.3$opt$aicc,EIWk15_C4_paranasal_fitkappa$Axis.3$opt$aicc,EIWk15_C4_paranasal_fitdelta$Axis.3$opt$aicc),
                   c("BM","OU","EB","trend","lambda","kappa","delta"))
 EIWk15_C4_paranasal_aic_vals3

txtStop()

 beep(1)


# ======================== #
# 0.1 Parallelisation ends #
# ======================== #

# stop cluster
stopImplicitCluster()


# ========= #
# 0.2 Done! #
# ========= #

### Noise when script has finished (useful if you are running this in the background):
 library(beepr)

 beep(8)
