# work with the develop proteome of Manduca sexta
filename <- 'C:\\Users\\ATPs\\OneDrive\\Lab\\works\\2017DmSPSPH\\20180166ManducaMassSpectDevelopMembrane\\20190209MsMemGroup_ori.xlsx'


library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)



mydata <- read_excel(filename)
