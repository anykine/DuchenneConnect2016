#######################################
# Project: DuchenneConnect2016
# Purpose: Analyze DuchenneConnect Data from October 2016 & October 2017
# Module: DuchenneConnect2016/SRC/rclean/104do.R
# Author: Richard T Wang
#

# load data sets


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

#load("SRC/rclean/DCdata201710.RData")

list.files("DATA/2016October/")
data_dir = "DATA/2016October/"

source("SRC/rclean/103functions.R")

#
# ----------- skip -------------
skips= read.csv("DATA/exon skippable mutations1.csv")
skips=skips %>% mutate(startstop = paste(exon_start, exon_stop, sep=","))


#
#
# ----------- MUSCLE DATA -------------
#
muscle = fread(paste(data_dir, "MuscleFxn Module Oct 2016.csv", sep="/"))
setnames(muscle, make.names.prefix(colnames(muscle),"Mus"))   # fix the spaces in names



#
# ----------- Genetic DATA -------------
#

gen = fread(paste(data_dir, "GeneticTest_MutationInfo Oct 2016.csv", sep="/"))
setnames(gen, make.names.prefix(colnames(gen), "Gen"))   # fix the spaces in names






#
# ----------- Steroid DATA -------------
#

ster = fread(paste(data_dir, "Corticosteroid Module Oct 2016.csv", sep="/"))
setnames(ster, make.names.prefix(colnames(ster), "Ste"))   # fix the spaces in names



#
# ----------- cardio DATA -------------
#
cardio = fread(paste(data_dir, "Cardiac Module Oct 2016.csv", sep="/"))
setnames(cardio, make.names.prefix(colnames(cardio), "Car"))   # fix the spaces in names

#save.image("SRC/rclean/DCdata201710.RData")
