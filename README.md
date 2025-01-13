# Santalaceae
Analyses of target sequence capture data for Santalaceae _sensu lato_

This repository contains the scripts and steps for analysing target sequence capture data in Anderson _et al._ (2025) "Evolutionary relationships in Santalales inferred using target capture with Angiosperms353, focusing on Australasian Santalaceae _sensu lato_"  

Scripts are located in the `scripts` folder (and subfolders), an R markdown file for tree plotting in the `rmd` folder, and Singularity recipes in `singularity`  
Some additional label files and the HybPiper targets file are in the `files` folder  

Raw sequencing data are available at the European Nucleotide Archive under projects PRJEB49212 (GAP stage 1), PRJEB78980 (GAP stage 2) and PRJEB79126 (samples sequenced outside Australia)  
Note that uploading data to ENA removes Illumina information from read headers, so the filtering step to remove optical duplicates needs to be turned off  

The full set of analysis steps can be followed in the `Sant_analyses.md` file, though many steps are specific to the NCI Gadi supercomputer where these analyses were run  
