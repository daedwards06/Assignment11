setwd("~/SMUWork")
library(affy)

source("~/SMUWork/SimulationMain.R")
values = simulate.array("expaffy.ini","configaffy.ini")

#cel = ReadAffy("GSM269808.CEL",celfile.path = "./Exp_Raw/GSE10685_RAW")


#pm.index = unlist(indexProbes(cel,which = "pm"))
#mm.index = unlist(indexProbes(cel,which = "mm"))

#writeBin(pm.index,"PMIndex.bin")
#writeBin(mm.index,"MMIndex.bin")

#rm(pm.index)
#rm(mm.index)

#pm.index = readBin("PMIndex.bin",what="integer",n=10000000)
#mm.index = readBin("MMIndex.bin",what="integer",n=10000000)