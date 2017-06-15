# analysis of DuchenneConnect Data from November 2016

# NOTES 11/2/16
#  Seems to be too few exon 44 skippable patients in data set. 
# 11/29/16
#   Compute your own Exon44 skip and do not rely on amenable to exon 44 skip

# load data sets, get last row to dedup, merge, filter

## 

library(ggplot2)
library(dplyr)
library(DT)
library(data.table)
library(dtplyr)

setwd("C://Users//Richard//Box Sync/projects/DuchenneConnect/")
#setwd("/home/rwang/box.com/projects/DuchenneConnect/")
list.files("DATA/2016October/")
data_dir = "DATA/2016October/"


#
# ----------- skip -------------
skips= read.csv("DATA/exon skippable mutations1.csv")
skips=skips %>% mutate(startstop = paste(exon_start, exon_stop, sep=","))

#
#
# ----------- MUSCLE DATA -------------
#
muscle = fread(paste(data_dir, "MuscleFxn Module Oct 2016.csv", sep="/"))
setnames(muscle, make.names(colnames(muscle)))   # fix the spaces in names

#muscle = read.csv(paste(data_dir, "MuscleFxn Module Oct 2016.csv", sep="/"))

# basic stats
names(muscle)
table(muscle$Gender)
table(muscle$Diagnosis)
table(muscle$Ethnicity)                      
xtabs( ~Diagnosis, data = muscle) %>% t

# how many nonwalking boys?
# q1 = 3, DMD, male
muscle %>% filter(Q1 == 3, Diagnosis=="Duchenne", Gender=="Male") %>% dim
muscle %>% filter(Q1 == 3, Diagnosis=="Duchenne", Gender=="Male") %>% 
  xtabs(~Country, data = .)


# 235 folks have multiple profiles
muscle %>% group_by(Patient.ID) %>% summarize(count=n()) %>%
  filter(count>1) %>% arrange(desc(count))

# female carriers check
muscle %>% filter(Gender=="Female") %>% datatable
muscle %>% filter(Gender=="Female") %>% xtabs( ~Diagnosis, data=.)

# Data cleanup
setkey(muscle, Patient.ID)
muscle2 = muscle[, .SD[c(.N)], by=Patient.ID]  # this gets the last row for each patient record 

#
# ----------- Genetic DATA -------------
#

gen = fread(paste(data_dir, "GeneticTest_MutationInfo Oct 2016.csv", sep="/"))
setnames(gen, make.names(colnames(gen)))   # fix the spaces in names
#gen = read.csv(paste(data_dir, "GeneticTest_MutationInfo Oct 2016.csv", sep="/"))

MultipleRowsForPatient(gen, returnData = T)  #35 dups

# get the last observation for each group
setkey(gen, Patient.ID)
gen2= gen[, .SD[c(.N)], by=Patient.ID]  # this gets the last row for each patient record 

# REQUIRED fix columns with same name, for merge later
# REQUIRED There are 3 columns with name "Other Value"
tmp.names = make.names(colnames(gen))
tmp.names[c(22,24,26)] = c("Laboratory.Other.Value", "Test.Method.Other.Value", "Category.Other.Value")
setnames(gen, tmp.names)
setnames(gen2, tmp.names)


# ignore this step
#how many exon44 skippables only 14? Check with previous data set
# Field AmenableToSkipExon1 is NOT COMPLETED for all patients
gen2 %>% filter(Amendable.To.Skip.Exon.1==44)



#
# ----------- Steroid DATA -------------
#

ster = fread(paste(data_dir, "Corticosteroid Module Oct 2016.csv", sep="/"))
setnames(ster, make.names(colnames(ster)))   # fix the spaces in names
#ster = read.csv(paste(data_dir, "Corticosteroid Module Oct 2016.csv", sep="/"))

MultipleRowsForPatient(ster, returnData = T)  #142 dups

setkey(ster, Patient.ID)
ster2 = ster[ , .SD[c(.N)], by=Patient.ID]

#
# ---------- read in skips ------------
# 
skips = read.csv("DATA/exon skippable mutations1.csv")

#
# ----------- combine all ------------
#
#check no dups
lapply(list(muscle2,ster2, gen2), checkUniquePatientsPerRow, myKey="Patient.ID")

# Age.m = muscle, Age.s = steroid, Age = Genetic Excel table
all = merge(merge(muscle2, ster2, by.x="Patient.ID", by.y="Patient.ID", suffixes=c(".m",".s")),
            gen2, by.x="Patient.ID", by.y="Patient.ID", suffixes=c(".ms", ".g"))

# filters
all %>% group_by(Gender) %>% summarize(count=n())

# preliminary filter
# oecd,duchenne,male,validIntron,deletion,validSteroid,validWalking
all.0 = all %>% filter(
  Gender == "Male",
  Health.Reason.for.Registration.g == "Duchenne",
  Diagnosis == "Duchenne",
  Country %in% c("AUSTRALIA", "AUSTRIA", "BELGIUM", "CANADA", "CHILE", "CZECH REPUBLIC", "DENMARK", "ESTONIA",
                 "FINALND", "FRANCE", "GERMANY", "GREECE", "HUNGARY", "ICELAND", "IRELAND", "ISREAL", "ITALY",
                 "JAPAN", "KOREA", "LATVIA", "LUXEMBOURG","MEXICO","NETHERLANDS","NEW ZEALAND","NORWAY",
                 "POLAND", "PORTUGAL", "SLOVAK REPUBLIC", "SLOVENIA", "SPAIN", "SWEDEN", "SWITZERLAND", "TURKEY",
                 "UNITED KINGDOM", "UNITED STATES"),
  Category == "Deletion",
  Start.Exon.Intron != "",
  End.Exon.Intron != "",
  # steroid 1=pred,2=dflz,3=stopped,4=never
  Q1.s %in% c(1,2,3,4),
  #Walking 1=own,2=partial,3=WC,4=infant
  Q1.m %in% c(1,2,3)
) %>% 
  mutate(
    Start.Exon.Intron = gsub("^Ex", "", Start.Exon.Intron),
    End.Exon.Intron = gsub("^Ex", "", End.Exon.Intron),
    startstop = paste(Start.Exon.Intron, End.Exon.Intron, sep=",")) %>% 
  merge(skips, by.x="startstop", by.y="startstop", all.x=T) %>%
  mutate(
    walking_01 = ifelse(Q3.2.m > 5 & !is.na(Q3.2.m), 0, 1)
  ) %>% 
  mutate(
    time_to_wheelchair = ifelse(walking_01==0, Q3.2.m, Age.m)
  )


all.01 = all.0 %>% filter(
    # steroid 1=pred,2=dflz,3=stopped,4=never
    Q1.s %in% c(1,2),
    #Walking 1=own,2=partial,3=WC,4=infant
    Q1.m %in% c(1,2,3),
    Age.m < 25,
    Age.m > 5
  ) %>% 
  mutate(
    strat = ifelse(skip_to_render_inframe=="44" & Start.Exon.Intron==45 & End.Exon.Intron==45, "skip44", "noskip44")
    #strat = ifelse(skip_to_render_inframe=="44" , "skip44", "noskip44")
  )

all.ster.44 = all.0 %>% filter(
  # steroid 1=pred,2=dflz,3=stopped,4=never
  Q1.s %in% c(1,2),
  #Walking 1=own,2=partial,3=WC,4=infant
  Q1.m %in% c(1,2,3),
  #Age.m < 25,
  Age.m > 5
) %>% 
  mutate(
    strat = ifelse(skip_to_render_inframe=="44", "skip44", "noskip44"),
    s44 = ifelse(skip_to_render_inframe==44, 1, 0),
    s45 = ifelse(skip_to_render_inframe==45, 1, 0),
    s50 = ifelse(skip_to_render_inframe==50, 1, 0),
    s51 = ifelse(skip_to_render_inframe==51 , 1, 0),
    s53 = ifelse(skip_to_render_inframe==53, 1, 0)
  )
formula = Surv(time_to_wheelchair, walking_01==0) ~ strat
#formula = Surv(time_to_wheelchair, walking_01==0) ~ s44
plot(survfit(formula, data = all.ster.44, conf.type="log-log"), lty=1:5, col=1:5)
fit.diff = survdiff(formula, data = all.ster.44)
fit.diff


# ***44 more mild; limit age < 25 seems make difference
# 1. Looking at all 44 skippables, p-0.06 survival difference
# 2. If we just look at 44 skippables, E45del, p is 0.02
# & Start.Exon.Intron==45 & End.Exon.Intron==45

# ***51 more severe? look at start/stop exons
all.ster.44 %>% filter(skip_to_render_inframe==51) %>% xtabs(~Start.Exon.Intron+End.Exon.Intron,data=.)
all.ster.51 = all.ster.44 %>% filter(
  skip_to_render_inframe==51
) %>% mutate(
  reg4550 = ifelse(Start.Exon.Intron==45, 1, 0),
  reg4850 = ifelse(Start.Exon.Intron==48, 1, 0),
  reg4950 = ifelse(Start.Exon.Intron==49, 1, 0),  #this one is severe
  reg5050 = ifelse(Start.Exon.Intron==50, 1, 0),
  reg5252 = ifelse(Start.Exon.Intron==52, 1, 0)
)
formula = Surv(time_to_wheelchair, walking_01==0) ~ reg4950
plot(survfit(formula, data = all.ster.51, conf.type="log-log"), lty=1:5, col=1:5)
fit.diff = survdiff(formula, data = all.ster.51)
fit.diff

# ***keep this plot. boxplot by startstop
all.ster.51 %>% filter(walking_01==0) %>% ggplot(aes(x=startstop,y=time_to_wheelchair))+
  geom_boxplot()

# scratch
all.0 %>% filter(skip_to_render_inframe==44) %>% select( Start.Exon.Intron, End.Exon.Intron) %>% 
  xtabs(~Start.Exon.Intron+ End.Exon.Intron, data = .)

all.1 = all.0 %>% select(Q1.m, Q3.2.m, Age.m)
all.1 %>% filter(!is.na(Q3.2.m)) %>% xtabs( ~ Q1.m + Age.m, data = .)
all.2 = all.1 %>% mutate(walking_01=ifelse(Q1.m==3 & Q3.2.m > 5, 0, 1))


validAge = function(age){
  if age > 4
}

# 44 skippable, steroid+, age stopped walking mean = 12.44, SD=4.5, N=20
all %>% mutate(
  Start.Exon.Intron = gsub("^Ex", "", Start.Exon.Intron),
  End.Exon.Intron = gsub("^Ex", "", End.Exon.Intron),
  startstop = paste(Start.Exon.Intron, End.Exon.Intron, sep=",")) %>% 
merge(skips, by.x="startstop", by.y="startstop", all.x=T) %>%
filter(
  Diagnosis=="Duchenne",
  Gender == "Male",
  Category == "Deletion",
  Start.Exon.Intron != "",
  End.Exon.Intron != "",
  # steroid 1=pred,2=dflz,3=stopped,4=never
  Q1.y %in% c(1,2),
  #Walking 1=own,2=partial,3=WC,4=infant
  Q1.x %in% c(1,2,3,4),
  Q1.x==3,
  #ex44 skippable
  skip_to_render_inframe==44
) %>% 
select(
    Q3.2.x
  ) %>% na.omit() %>% 
summarize(
  avg = mean(Q3.2.x),
  stdev = sd(Q3.2.x)
)

# # 44 NOT skippable, steroid+, age stopped walking mean = 11.33, sd=2.65, N=139
all %>% mutate(
  Start.Exon.Intron = gsub("^Ex", "", Start.Exon.Intron),
  End.Exon.Intron = gsub("^Ex", "", End.Exon.Intron),
  startstop = paste(Start.Exon.Intron, End.Exon.Intron, sep=",")) %>% 
  merge(skips, by.x="startstop", by.y="startstop", all.x=T) %>%
  filter(
    Diagnosis=="Duchenne",
    Gender == "Male",
    Category == "Deletion",
    Start.Exon.Intron != "",
    End.Exon.Intron != "",
    # steroid 1=pred,2=dflz,3=stopped,4=never
    Q1.y %in% c(1,2),
    #Walking 1=own,2=partial,3=WC,4=infant
    Q1.x %in% c(1,2,3,4),
    Q1.x==3,
    #ex44 skippable
    skip_to_render_inframe!=44
  ) %>% 
select(
    Q3.2.x
  ) %>% na.omit() %>% 
summarize(
  avg = mean(Q3.2.x), 
  stdev=sd(Q3.2.x)
)

#
# ----- filter data set -----
#
all %>% select(Country)

# DO the survival analysis
# is WCage survival different between 44skippable and non44skippable?
#   create time_to_wheelchair = if walking then,curr age; if not walking, WCage
#          
all44 = all %>% mutate(
  Start.Exon.Intron = gsub("^Ex", "", Start.Exon.Intron),
  End.Exon.Intron = gsub("^Ex", "", End.Exon.Intron),
  startstop = paste(Start.Exon.Intron, End.Exon.Intron, sep=",")) %>% 
  merge(skips, by.x="startstop", by.y="startstop", all.x=T) %>%
  filter(
    Diagnosis=="Duchenne",
    Gender == "Male",
    Category == "Deletion",
    Start.Exon.Intron != "",
    End.Exon.Intron != "",
    # steroid 1=pred,2=dflz,3=stopped,4=never
    Q1.y %in% c(1,2),
    #Walking 1=own,2=partial,3=WC,4=infant
    Q1.x %in% c(1,2,3,4),
    Q1.x==3,
    #ex44 skippable
    skip_to_render_inframe==44
  ) %>% 
  mutate(
    walking_01 = ifelse(Q1.x==3,0,1)
  ) %>%
  mutate(
    time_to_wheelchair = ifelse(walking_01==0, Q3.2.x, Age.x),
    strat=ifelse(skip_to_render_inframe==44 , as.character(skip_to_render_inframe), "not44")
  )


names(all)[grep("Q3.2", names(all))]

# ----- functions ------


MultipleRowsForPatient = function(data, returnData = F) {
  tab = data %>% group_by(Patient.ID) %>% summarize(count=n()) %>%
    filter(count > 1)
  if (returnData == TRUE) {
    tab
  } else {
    dim(tab)
  }
}

checkUniquePatientsPerRow = function(data, myKey="Patient.ID") {
  data %>% group_by_(myKey) %>% summarize(count=n()) %>% filter(count>1) %>% as.data.frame
}

## ---- EDA ----
gen %>% group_by(Amendable.To.Skip.Exon.1) %>% summarize(count=n()) %>%
                                                           as.data.frame
new1 = gen %>% filter(Amendable.To.Skip.Exon.1==44) %>% select(Patient.ID)


oldgen = read.csv("C:\\Users\\Richard\\Documents\\nelsonMiceliLabs\\projDuchenneConnect\\2013data\\coded 2013-08-13\\StanNelsonDataExport_5-7-2013_coded_08132013.csv", header=T)
names(oldgen) = sapply(names(oldgen), tolower)
skips = read.csv("C:\\Users\\Richard\\Documents\\nelsonMiceliLabs\\projDuchenneConnect\\skips\\exon skippable mutations1.csv")
skips$startstop = paste(skips[,1], skips[,2], sep=",")

x.cols = oldgen %>% select(patientid, time_to_wheelchair, corticosteroids_01, walking_01, category, start.exon.intron, end.exon.intron, km) %>%
  filter(km==1, category=="Deletion") %>% 
  mutate(start.exon.intron = gsub("^Ex", "", start.exon.intron),
         end.exon.intron = gsub("^Ex", "", end.exon.intron),
         startstop = paste(start.exon.intron, end.exon.intron, sep=",")) %>% 
  merge(skips, by.x="startstop", by.y="startstop", all.x=T)

old1 = x.cols %>% filter(skip_to_render_inframe==44) %>% select(patientid)


## ---- steroid effect -----
# let's check the steroid effect to see if our data is replicable. YES REPLICATES
# 1. steroid versus no steroid
# 2. pred versus deflaz
# 3. pred/dflz versus previous versus never
all.ster = all.0 %>% filter(
  Age.m < 25
) %>% 
  mutate(
    strat = ifelse(Q1.s == 1, "pred",  ifelse(Q1.s==2, "dflz", ifelse(Q1.s==3,"prev","never")))
  ) 

#strat=ifelse(skip_to_render_inframe==44 , as.character(skip_to_render_inframe), "not44")
formula = Surv(time_to_wheelchair, walking_01==0) ~ strat
plot(survfit(formula, data = all.ster, conf.type="log-log"), lty=1:5, col=1:5)
fit.diff = survdiff(formula, data = all.ster)
fit.diff
legend("bottomleft",legend=c("dflz", "never","pred","prev"),fill=1:5)


survivalVis(cobj, data=veteran, plot.title="Veteran Survival Data", group="trt", 
            group.names=c("Treatment", "No Treatment"), line.col=c("#E495A5","#39BEB1"))