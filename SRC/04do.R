library(survival)
#
#
# ----------- MUSCLE DATA -------------
#
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

#
#
# ----------- GENETIC DATA -------------
#
# 

MultipleRowsForPatient(gen, returnData = T)  #35 dups

# ignore this step
#how many exon44 skippables only 14? Check with previous data set
# Field AmenableToSkipExon1 is NOT COMPLETED for all patients
gen2 %>% filter(Amendable.To.Skip.Exon.1==44)



#
# ----------- Steroid DATA -------------
#

MultipleRowsForPatient(ster, returnData = T)  #142 dups





#
# ----------- combine all ------------
#
#check no dups
lapply(list(muscle2,ster2, gen2), checkUniquePatientsPerRow, myKey="Patient.ID")

lapply(list(muscle2,ster2, gen2), MultipleRowsForPatient)


#
# ----------- Basic question ------------
#

# sex breakdown
all %>% group_by(Gender) %>% summarize(count=n())


#
# ----------- Basic filter ------------
#

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


# ------------ what is the mutation class -------------
all.0 %>% group_by(Category) %>% summarize(cnt=n())
all.0 %>% group_by(Frame) %>% summarize(cnt=n())
# deletion mutation spectrum where N>5
# deletions_n5.png
all.0 %>% group_by(startstop) %>% summarize(cnt=n()) %>% 
  filter(cnt>5) %>% as.data.frame %>%
  ggplot(aes(x=factor(startstop), y=cnt)) + geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle=90)) + ggtitle("mutations with N > 5")
# barchar: exon skip required to put in frame
# skip_reqd_for_inframe_bar.png
all.0 %>% group_by(skip_to_render_inframe) %>% summarize(cnt=n()) %>% 
  arrange(skip_to_render_inframe) %>% as.data.frame %>%
  ggplot(aes(x=factor(skip_to_render_inframe), y=cnt)) + geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle=90)) + ggtitle("count of skip reqd to put inframe ")




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

# ------------ DATASET all steroid valid_walking strats  --------

# IFELSE ignores NA values in skip_to_render_inframe, so it effectively compares
# for example, 44 skippables ('skip44') to 8/45/50/51/52/53/55 BUT it does not set NA values to 
# 'noskip44'. 
# SOLVED the NA problem (eg s8na). flags wigh "na" on the end convert NA's to zeros
all.ster.44 = all.0 %>% filter(
  # steroid 1=pred,2=dflz,3=stopped,4=never
  Q1.s %in% c(1,2),
  #Walking 1=own,2=partial,3=WC,4=infant
  Q1.m %in% c(1,2,3),
  Age.m < 25,
  Age.m > 5
) %>% 
  mutate(
    #strat = ifelse(skip_to_render_inframe=="44", "skip44", "noskip44"),
    strat = ifelse(is.na(skip_to_render_inframe) | skip_to_render_inframe!=44, "noskip44", "skip44"),
    s8 = ifelse(skip_to_render_inframe==8, 1, 0),
    s8na = ifelse(is.na(skip_to_render_inframe) | skip_to_render_inframe!=8, 0, 1),
    s44 = ifelse(skip_to_render_inframe==44, 1, 0),
    s44na = ifelse(is.na(skip_to_render_inframe) | skip_to_render_inframe!=44, 0, 1),
    s45 = ifelse(skip_to_render_inframe==45, 1, 0),
    s45na = ifelse(is.na(skip_to_render_inframe) | skip_to_render_inframe!=45, 0, 1),
    s50 = ifelse(skip_to_render_inframe==50, 1, 0),
    s50na = ifelse(is.na(skip_to_render_inframe) | skip_to_render_inframe!=50, 0, 1),
    s51 = ifelse(skip_to_render_inframe==51 , 1, 0),
    s51na = ifelse(is.na(skip_to_render_inframe) | skip_to_render_inframe!=51, 0, 1),
    s52 = ifelse(skip_to_render_inframe==52, 1, 0),
    s52na = ifelse(is.na(skip_to_render_inframe) | skip_to_render_inframe!=52, 0, 1),
    s53 = ifelse(skip_to_render_inframe==53, 1, 0),
    s53na = ifelse(is.na(skip_to_render_inframe) | skip_to_render_inframe!=53, 0, 1),
    s55 = ifelse(skip_to_render_inframe==55, 1, 0),
    s55na = ifelse(is.na(skip_to_render_inframe) | skip_to_render_inframe!=55, 0, 1)
  )
makeSurvPlot("s8", all.ster.44, "8 skippable")       # Significant
makeSurvPlot("s8na", all.ster.44, "8 skippable NA")  # Significant
makeSurvPlot("s44", all.ster.44, "44 skippable")
makeSurvPlot("s44na", all.ster.44, "44 skippable NA")  # Significant
makeSurvPlot("s45", all.ster.44, "45 skippable")
makeSurvPlot("s45na", all.ster.44, "45 skippable NA")
makeSurvPlot("s50", all.ster.44, "50 skippable")
makeSurvPlot("s50na", all.ster.44, "50 skippable NA")
makeSurvPlot("s51", all.ster.44, "51 skippable")       # Significant
makeSurvPlot("s51na", all.ster.44, "51 skippable NA")  # Significant
makeSurvPlot("s52", all.ster.44, "52 skippable")
makeSurvPlot("s52na", all.ster.44, "52 skippable NA")
makeSurvPlot("s53", all.ster.44, "53 skippable")
makeSurvPlot("s53na", all.ster.44, "53 skippable NA")
makeSurvPlot("s55", all.ster.44, "55 skippable")
makeSurvPlot("s55na", all.ster.44, "55 skippable NA")


# EDA
all.ster.44 %>% filter(s44na==1) %>% select(skip_to_render_inframe, Category, Start.Exon.Intron, 
                                                   End.Exon.Intron, walking_01, time_to_wheelchair) %>%
  datatable()

all.ster.44 %>% select(skip_to_render_inframe) %>% table(useNA="ifany")
all.ster.44 %>% select(s55na) %>% table(useNA="ifany")

#######  todo: Exon 44, E45 del ######

# ---------- E44 mild in depth ------------
# ***44 more mild; limit age < 25 seems make difference
# 1. Looking at all 44 skippables, p-0.06 survival difference
# 2. If we just look at 44 skippables, E45del, p is 0.02
# & Start.Exon.Intron==45 & End.Exon.Intron==45

# list all start/stop for E44 skippables
all.0 %>% filter(skip_to_render_inframe==44) %>% select( Start.Exon.Intron, End.Exon.Intron) %>% 
  xtabs(~Start.Exon.Intron+ End.Exon.Intron, data = .)

# this compares E45del LOA only WITHIN all 44 skippables
all.ster.44.detail = all.ster.44 %>% filter(
  skip_to_render_inframe==44
) %>% mutate(
  reg4545 = ifelse(Start.Exon.Intron==45 & End.Exon.Intron==45, 1, 0),
  reg4554 = ifelse(Start.Exon.Intron==45 & End.Exon.Intron==54, 1, 0)
)
makeSurvPlot("reg4545", all.ster.44.detail, "44 skippable E45del")   # only compared E45del to nonE45del that are 44 skippable
makeSurvPlot("reg4554", all.ster.44.detail, "44 skippable E45del")   # only compared E45_54del to nonE45del that are 44 skippable

# Compares E44 skippable (E45del) to ALL other deletions
# Significant
all.0 %>% filter(
  # steroid 1=pred,2=dflz,3=stopped,4=never
  Q1.s %in% c(1,2),
  #Walking 1=own,2=partial,3=WC,4=infant
  Q1.m %in% c(1,2,3),
  Age.m < 25,
  Age.m > 5
) %>% 
  mutate(
    strat = ifelse(skip_to_render_inframe=="44" & Start.Exon.Intron==45 & End.Exon.Intron==45, "skip44", "noskip44"),
    reg4545 = ifelse(Start.Exon.Intron==45 & End.Exon.Intron==45, 1, 0),
    reg4554 = ifelse(Start.Exon.Intron==45 & End.Exon.Intron==54, 1, 0)
    #strat = ifelse(skip_to_render_inframe=="44" , "skip44", "noskip44")
) %>% makeSurvPlot("strat", ., "44 skippable (E45del) versus all other deletions")

# ---------- E51 severe in depth --------
# ***51 more severe? look at start/stop exons
# list all start/stop for E51 skippables
all.ster.44 %>% filter(skip_to_render_inframe==51) %>% xtabs(~Start.Exon.Intron+End.Exon.Intron,data=.)

# Compares subgroups of E51 skippables to all E51 skippables
all.ster.51 = all.ster.44 %>% filter(
  skip_to_render_inframe==51
) %>% mutate(
  reg4550 = ifelse(Start.Exon.Intron==45, 1, 0),
  reg4850 = ifelse(Start.Exon.Intron==48, 1, 0),
  reg4950 = ifelse(Start.Exon.Intron==49, 1, 0),  #this one is severe
  reg5050 = ifelse(Start.Exon.Intron==50, 1, 0),
  reg5252 = ifelse(Start.Exon.Intron==52, 1, 0)
)

makeSurvPlot("reg4550", all.ster.51, "51 skippable E45_50del") 
makeSurvPlot("reg4850", all.ster.51, "51 skippable E48_50del") 
makeSurvPlot("reg4950", all.ster.51, "51 skippable E49_50del") # Significant
makeSurvPlot("reg5050", all.ster.51, "51 skippable E50del") 
makeSurvPlot("reg5252", all.ster.51, "51 skippable E52del") 


# Compares E51 skippable (E49_50del) to all other deletions
# SIGNIFICANT
all.0 %>% filter(
  # steroid 1=pred,2=dflz,3=stopped,4=never
  Q1.s %in% c(1,2),
  #Walking 1=own,2=partial,3=WC,4=infant
  Q1.m %in% c(1,2,3),
  Age.m < 25,
  Age.m > 5
) %>% 
  mutate(
    strat = ifelse(skip_to_render_inframe=="51" & Start.Exon.Intron==49 & End.Exon.Intron==50, "skip44", "noskip44")
    #strat = ifelse(skip_to_render_inframe=="44" , "skip44", "noskip44")
  ) %>% makeSurvPlot("strat", ., "51 skippable E49_50del versus all other deletions")


# ***keep this plot. boxplot by startstop, WC age for E51 skippables
all.ster.51 %>% filter(walking_01==0) %>% ggplot(aes(x=startstop,y=time_to_wheelchair))+
  geom_boxplot()



# delete?
all.1 = all.0 %>% select(Q1.m, Q3.2.m, Age.m)
all.1 %>% filter(!is.na(Q3.2.m)) %>% xtabs( ~ Q1.m + Age.m, data = .)
all.2 = all.1 %>% mutate(walking_01=ifelse(Q1.m==3 & Q3.2.m > 5, 0, 1))




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

makeSurvPlot("strat", all.ster, "steroid type effect")  # steroid_type_survival.png

# #strat=ifelse(skip_to_render_inframe==44 , as.character(skip_to_render_inframe), "not44")
# formula = Surv(time_to_wheelchair, walking_01==0) ~ strat
# plot(survfit(formula, data = all.ster, conf.type="log-log"), lty=1:5, col=1:5)
# fit.diff = survdiff(formula, data = all.ster)
# fit.diff
# legend("bottomleft",legend=c("dflz", "never","pred","prev"),fill=1:5)


# steroid_v_nosteroid.png
all.0 %>% filter(
  Age.m < 25
) %>% 
  mutate(
    strat = ifelse(Q1.s == 1 | Q1.s==2, "ster",  ifelse(Q1.s==3|Q1.s==4, "noster", NA))
) %>% makeSurvPlot("strat", . , "steroid v no steroid")



# healthvis function
survivalVis(cobj, data=veteran, plot.title="Veteran Survival Data", group="trt", 
            group.names=c("Treatment", "No Treatment"), line.col=c("#E495A5","#39BEB1"))

### ----- WCage by exon skippable mutation ----
# exon50/51 skippables have earlier WCage
# exon44 skippables have later WCage  
all.ster.any = all.0 %>% filter(
  walking_01 == 0,
  Age.m < 25
)

all.ster.any %>% ggplot(aes(x=factor(skip_to_render_inframe), y=time_to_wheelchair)) +
  geom_boxplot() + geom_jitter() + ggtitle("WCage by exon_skip_to_render_inframe")+
  coord_flip()

all.ster.any %>% xtabs(~ skip_to_render_inframe + time_to_wheelchair, data = .)
  
# exon51





# -------------- Cardio --------------
# examine effect of steroid on cardio
# examine effect of cardiac drugs on cardio
# examine WC age effects due to cardio drugs?
# is cardio predictive of WCage?

#how many ppl had decreased heart function (DHF) diagnosis 
all.0 %>% filter(Q3==1) %>% dim #147

# age steroid vs age heart problem strat by steroid status(1=pred,2=dflz)
all.0 %>% select(Q1.s, Q1.1, Q3.1) %>% ggplot(aes(x=Q1.1, y=Q3.1))+
  geom_point(aes(color=factor(Q1.s)), size=3) + geom_jitter() +
  xlab("Age started corticosteroids") + ylab("Age decreas heart function") +
  ylim(0,25)

# hist of age decres heart function strat by steroid type
all.0 %>% filter(Q1.s==1 | Q1.s==2) %>% select(Q1.s, Q1.1, Q3.1) %>% 
  ggplot(aes(Q3.1)) + geom_density(aes(color=as.factor(Q1.s))) +
  xlim(0,25) + xlab("age DHF") + ggtitle("age DHF strat steroid type")
all.0 %>% filter(Q1.s==1 | Q1.s==2, Q1.1 < Q3.1) %>% dim #41 ppl started ster before DHF
all.0 %>% filter(Q1.s==1 | Q1.s==2, Q1.1 > Q3.1) %>% dim #4 ppl started ster AFTER DHF
all.0 %>% filter(Q1.s==1 | Q1.s==2) %>% dim

# 7 ppl beta blockers (start date before DHF)
all.0 %>% filter(Q1a_start < Q3.1) %>% dim
# 17 ppl ACE
all.0 %>% filter(Q1b._start < Q3.1) %>% dim
# 4  ppl ARB
all.0 %>% filter(Q1c_start < Q3.1) %>% dim
# 0  ppl Ca block
all.0 %>% filter(Q1d_start < Q3.1) %>% dim
# 1 diuretics
all.0 %>% filter(Q1e_start < Q3.1) %>% dim
# 0 adrenergic
all.0 %>% filter(Q1f_start < Q3.1) %>% dim
# 2 heart failure
all.0 %>% filter(Q1g_start < Q3.1) %>% dim
# 0 antiarrhythmia 
all.0 %>% filter(Q1h_start < Q3.1) %>% dim
# 0 mineralcorticoid 
all.0 %>% filter(Q1i_start < Q3.1) %>% dim

# 17 ppl ACE
# compare this plot to age DHF strat by steroid type - looks the same

ace = all.0 %>% filter(Q1b._start < Q3.1)
ace %>% filter(Q1.s==1 | Q1.s==2) %>% select(Q1.s, Q1.1, Q3.1) %>% 
  ggplot(aes(Q3.1)) + geom_density(aes(color=as.factor(Q1.s))) +
  xlim(0,25) + xlab("age DHF") + ggtitle('DHF for ACE users N=17 (strat by steroid type)')

# are all ACE users also steroid users (Q1.s=1,2)?
all.0 %>% filter(!is.na(Q1b._start)) %>% select(Q1b._start, Q1.s) %>%
  xtabs(~Q1.s + Q1b._start, data = .)
