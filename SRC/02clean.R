# Some data sets have multiple rows per individual (time course)
# so I use data.tables() to get last row to dedup

#
#
# ----------- MUSCLE DATA -------------
#
# Data cleanup
setkey(muscle, Patient.ID)
muscle2 = muscle[, .SD[c(.N)], by=Patient.ID]  # this gets the last row for each patient record 


#
#
# ----------- GENETIC DATA -------------
#
# get the last observation for each group
setkey(gen, Patient.ID)
gen2= gen[, .SD[c(.N)], by=Patient.ID]  # this gets the last row for each patient record 

# REQUIRED fix columns with same name, for merge later
# REQUIRED There are 3 columns with name "Other Value"
tmp.names = make.names(colnames(gen))
tmp.names[c(22,24,26)] = c("Laboratory.Other.Value", "Test.Method.Other.Value", "Category.Other.Value")
setnames(gen, tmp.names)
setnames(gen2, tmp.names)

#
#
# ----------- STEROID DATA -------------
#
setkey(ster, Patient.ID)
ster2 = ster[ , .SD[c(.N)], by=Patient.ID]

#
#
# ----------- CARDIO DATA -------------
#
setkey(cardio, Patient.ID)
cardio2 = cardio[ , .SD[c(.N)], by=Patient.ID]

#
# ----------- combine muscle2, ster2, gen2 into 'all' ------------
#
#check no dups in each table, should be 0rows
lapply(list(muscle2,ster2, gen2), checkUniquePatientsPerRow, myKey="Patient.ID")

# Age.m = muscle, Age.s = steroid, Age = Genetic Excel table
all = merge(merge(muscle2, ster2, by.x="Patient.ID", by.y="Patient.ID", suffixes=c(".m",".s")),
            gen2, by.x="Patient.ID", by.y="Patient.ID", suffixes=c(".ms", ".g"))

#
# ----------- combine 'all' with cardio2 ------------
# Don't use yet
#all = merge(all, cardio2, by.x="Patient.ID", by.y="Patient.ID", suffixes=c("", ".c"))

# save.image("SRC/DCdata2016.RData")

