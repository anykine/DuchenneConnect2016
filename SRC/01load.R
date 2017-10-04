#######################################
# DuchenneConnect November 2016 data
#
#
# following the example here: http://stackoverflow.com/questions/1429907/workflow-for-statistical-analysis-and-report-writing
# load.R
# clean.R
# func.R
# do.R

# analysis of DuchenneConnect Data from November 2016

# NOTES 11/2/16
#  Seems to be too few exon 44 skippable patients in data set. 
# 11/29/16
#   Compute your own Exon44 skip and do not rely on amenable to exon 44 skip

# load data sets, 

## 

library(ggplot2)
library(dplyr)
library(DT)
library(data.table)
library(dtplyr)
library(RColorBrewer)
library(broom)

setwd("/home/rwang/projects/DuchenneConnect2016")  # lab desktop
#setwd("C://Users//Richard//Box Sync/projects/DuchenneConnect/")
#setwd("/home/rwang/box.com/projects/DuchenneConnect/")

#load("SRC/DCdata2016.RData")

list.files("DATA/2016October/")
data_dir = "DATA/2016October/"


#
# ----------- skip -------------
skips= read.csv("DATA/exon skippable mutations1.csv")
skips=skips %>% mutate(startstop = paste(exon_start, exon_stop, sep=","))
#

#
#
# ----------- MUSCLE DATA -------------
#
muscle = fread(paste(data_dir, "MuscleFxn Module Oct 2016.csv", sep="/"))
setnames(muscle, make.names(colnames(muscle)))   # fix the spaces in names

#muscle = read.csv(paste(data_dir, "MuscleFxn Module Oct 2016.csv", sep="/"))


#
# ----------- Genetic DATA -------------
#

gen = fread(paste(data_dir, "GeneticTest_MutationInfo Oct 2016.csv", sep="/"))
setnames(gen, make.names(colnames(gen)))   # fix the spaces in names
#gen = read.csv(paste(data_dir, "GeneticTest_MutationInfo Oct 2016.csv", sep="/"))





#
# ----------- Steroid DATA -------------
#

ster = fread(paste(data_dir, "Corticosteroid Module Oct 2016.csv", sep="/"))
setnames(ster, make.names(colnames(ster)))   # fix the spaces in names
#ster = read.csv(paste(data_dir, "Corticosteroid Module Oct 2016.csv", sep="/"))


#
# ----------- cardio DATA -------------
#
cardio = fread(paste(data_dir, "Cardiac Module Oct 2016.csv", sep="/"))
setnames(cardio, make.names(colnames(cardio)))   # fix the spaces in names


# save.image("SRC/DCdata2016.RData")
