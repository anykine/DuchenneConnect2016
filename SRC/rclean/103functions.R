#######################################
# Project: DuchenneConnect2016
# Purpose: Useful functions for analysis of DuchenneConnect Data 2016/2017
# Module: DuchenneConnect2016/SRC/rclean/103functions.R
# Author: Richard T Wang
#

# Append a prefix to column names
make.names.prefix <- function(columnNames, prefix){
  var1=make.names(columnNames)
  var2=paste0(prefix, var1[2:length(var1)])
  c(var1[1], var2)
}
  