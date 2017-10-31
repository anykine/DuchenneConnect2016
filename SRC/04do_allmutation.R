#----basic dataset wtih Deletion, Duplication, or Nonsense----
# Kaplan Meier Analysis and Log Rank Test
# repeat 04do.R using more DMD mutation categories

# this is equivalent to all.0 on 04do.R but includes Dups and Nonsense
all.anymutation = all %>% dplyr::filter(
  Gender.m == "Male",
  Health.Reason.for.Registration.g == "Duchenne",
  Diagnosis == "Duchenne",
  Country %in% c("AUSTRALIA", "AUSTRIA", "BELGIUM", "CANADA", "CHILE", "CZECH REPUBLIC", "DENMARK", "ESTONIA",
                 "FINALND", "FRANCE", "GERMANY", "GREECE", "HUNGARY", "ICELAND", "IRELAND", "ISREAL", "ITALY",
                 "JAPAN", "KOREA", "LATVIA", "LUXEMBOURG","MEXICO","NETHERLANDS","NEW ZEALAND","NORWAY",
                 "POLAND", "PORTUGAL", "SLOVAK REPUBLIC", "SLOVENIA", "SPAIN", "SWEDEN", "SWITZERLAND", "TURKEY",
                 "UNITED KINGDOM", "UNITED STATES"),
  Category %in% c("Deletion", "Duplication", "Nonsense"),
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

# ----Analyze Deletion, Duplication and Nonsense) + Steroids----
# this replicates 04do.R (all.ster.44 survival plots)
# stratifications are tricky! For "NA" groups, use NOT(deletion AND skip_to_render)
all.anymutation.ster = all.anymutation %>% filter(
  # steroid 1=pred,2=dflz,3=stopped,4=never
  Q1.s %in% c(1,2),
  #Walking 1=own,2=partial,3=WC,4=infant
  Q1.m %in% c(1,2,3),
  Age.m < 25,
  Age.m > 5
) %>% 
  mutate(
    s8Del = ifelse(skip_to_render_inframe==8 & Category=="Deletion", 1, 0),
    s8naDel = ifelse(is.na(skip_to_render_inframe) | !(skip_to_render_inframe==8 & Category == "Deletion"), 0, 1),
    s44Del = ifelse(skip_to_render_inframe==44 & Category=="Deletion", 1, 0),
    s44naDel = ifelse(is.na(skip_to_render_inframe) | !(skip_to_render_inframe==44 & Category=="Deletion") , 0, 1),
    s45Del = ifelse(skip_to_render_inframe==45 & Category=="Deletion", 1, 0),
    s45naDel = ifelse(is.na(skip_to_render_inframe) | !(skip_to_render_inframe==45 & Category=="Deletion"), 0, 1),
    s50Del = ifelse(skip_to_render_inframe==50 & Category=="Deletion", 1, 0),
    s50naDel = ifelse(is.na(skip_to_render_inframe) | !(skip_to_render_inframe==50 & Category=="Deletion"), 0, 1),
    s51Del = ifelse(skip_to_render_inframe==51  & Category=="Deletion", 1, 0),
    s51naDel = ifelse(is.na(skip_to_render_inframe) | !(skip_to_render_inframe==51 & Category=="Deletion"), 0, 1),
    s52Del = ifelse(skip_to_render_inframe==52 & Category=="Deletion", 1, 0),
    s52naDel = ifelse(is.na(skip_to_render_inframe) | !(skip_to_render_inframe==52 & Category=="Deletion"), 0, 1),
    s53Del = ifelse(skip_to_render_inframe==53 & Category=="Deletion", 1, 0),
    s53naDel = ifelse(is.na(skip_to_render_inframe) | !(skip_to_render_inframe==53 & Category=="Deletion"), 0, 1),
    s55Del = ifelse(skip_to_render_inframe==55 & Category=="Deletion", 1, 0),
    s55naDel = ifelse(is.na(skip_to_render_inframe) | !(skip_to_render_inframe==55 & Category=="Deletion"), 0, 1),
    Dup = ifelse(Category=="Duplication", 1, 0),
    Nonsense = ifelse(Category=="Nonsense", 1, 0)
  )
# LOA survival for mutatin type) + STEROID
makeSurvPlot("s8Del", all.anymutation.ster, "8 skippable")       # Significant
makeSurvPlot("s8naDel", all.anymutation.ster, "8 skippable NA")  # Significant
makeSurvPlot("s44Del", all.anymutation.ster, "44 skippable")
makeSurvPlot("s44naDel", all.anymutation.ster, "44 skippable NA")  # Significant
makeSurvPlot("s45Del", all.anymutation.ster, "45 skippable")
makeSurvPlot("s45naDel", all.anymutation.ster, "45 skippable NA")
makeSurvPlot("s50Del", all.anymutation.ster, "50 skippable")
makeSurvPlot("s50naDel", all.anymutation.ster, "50 skippable NA")
makeSurvPlot("s51Del", all.anymutation.ster, "51 skippable")       # Significant
makeSurvPlot("s51naDel", all.anymutation.ster, "51 skippable NA")  # Significant but compard to 04do.R where nonsense/dup NOT included, this is less significant
makeSurvPlot("s52Del", all.anymutation.ster, "52 skippable")
makeSurvPlot("s52naDel", all.anymutation.ster, "52 skippable NA")
makeSurvPlot("s53Del", all.anymutation.ster, "53 skippable")
makeSurvPlot("s53naDel", all.anymutation.ster, "53 skippable NA")
makeSurvPlot("s55Del", all.anymutation.ster, "55 skippable")
makeSurvPlot("s55naDel", all.anymutation.ster, "55 skippable NA")
makeSurvPlot("Dup", all.anymutation.ster, "Dupications")
makeSurvPlot("Nonsense", all.anymutation.ster, "Nonsense")

#rought counts ... numbers are small
all.anymutation.ster %>% filter(Category=="Deletion") %>% select(skip_to_render_inframe) %>% table()

# E44 skippable: Compare E45del  or E45-54del Deletions) + Steroids to everyone else
# E45dels are milder! Significant! P=0.029
all.anymutation.ster.44.detail = all.anymutation.ster %>% mutate(
  reg4545 = ifelse(Start.Exon.Intron==45 & End.Exon.Intron==45 & skip_to_render_inframe==44 & Category=="Deletion", 1, 0),
  reg4554 = ifelse(Start.Exon.Intron==45 & End.Exon.Intron==54 & skip_to_render_inframe==44 & Category=="Deletion", 1, 0)
)
makeSurvPlot("reg4545", all.anymutation.ster.44.detail, "44 skippable E45del")   # significant
makeSurvPlot("reg4554", all.anymutation.ster.44.detail, "44 skippable E45del")   

# E44 skippable: split into E45del, not E45del, and all other groups
# for 3-way survival plot 
# Use for paper.
all.anymutation.ster %>% filter(skip_to_render_inframe==44 & Category=="Deletion") %>%  
  xtabs(~startstop, data = .)
all.anymutation.ster.44.subgroups = all.anymutation.ster %>% mutate(
  s44_subgroups = ifelse(s44naDel==1 & Start.Exon.Intron==45 & End.Exon.Intron==45, "44skippable E45del",
                         ifelse(s44naDel==1 & Start.Exon.Intron!=45 & End.Exon.Intron!=45, "not44skippable E45del", "all others"))
) #%>% makeSurvPlot("s44_subgroups", data = . , "44 skippable by E45 subgroups")   


# E51 skippable: Compare subgroups of E51 skippables) + Steroids to everyone else
# E49_50dels are more severe! Significant P=0.00791
all.anymutation.ster.51.detail = all.anymutation.ster %>% mutate(
  reg4550 = ifelse(Start.Exon.Intron==45 & skip_to_render_inframe==51 & Category=="Deletion" , 1, 0),
  reg4850 = ifelse(Start.Exon.Intron==48 & skip_to_render_inframe==51 & Category=="Deletion", 1, 0),
  reg4950 = ifelse(Start.Exon.Intron==49 & skip_to_render_inframe==51 & Category=="Deletion", 1, 0),  #this one is severe
  reg5050 = ifelse(Start.Exon.Intron==50 & skip_to_render_inframe==51 & Category=="Deletion", 1, 0),
  reg5252 = ifelse(Start.Exon.Intron==52 & skip_to_render_inframe==51 & Category=="Deletion", 1, 0)
)

makeSurvPlot("reg4550", all.anymutation.ster.51.detail, "51 skippable E45_50del") 
makeSurvPlot("reg4850", all.anymutation.ster.51.detail, "51 skippable E48_50del") 
makeSurvPlot("reg4950", all.anymutation.ster.51.detail, "51 skippable E49_50del") # Significant
makeSurvPlot("reg5050", all.anymutation.ster.51.detail, "51 skippable E50del") 
makeSurvPlot("reg5252", all.anymutation.ster.51.detail, "51 skippable E52del") 

# E51 skippable by subgroups
all.anymutation.ster %>% filter(skip_to_render_inframe==51 & Category=="Deletion") %>%  
  xtabs(~startstop, data = .)
all.anymutation.ster.51.subgroups = all.anymutation.ster %>% mutate(
  s51_subgroups = ifelse(s51naDel==1 & Start.Exon.Intron == 45, "45-50del",
                         ifelse(s51naDel==1 & Start.Exon.Intron == 48 , "48-50del",
                                ifelse(s51naDel==1 & Start.Exon.Intron == 49 , "49-50del",
                                       ifelse(s51naDel==1 & Start.Exon.Intron == 50 , "50-50del",
                                              ifelse(s51naDel==1 & Start.Exon.Intron == 52, "52-52del", "all others"
                                              )))))
  )
  
  

#---- NO Steroid Use ----
# repeat above but using no steroid group
all.anymutation.noster = all.anymutation %>% filter(
  # steroid 1=pred,2=dflz,3=stopped,4=never
  Q1.s %in% c(3,4),
  #Walking 1=own,2=partial,3=WC,4=infant
  Q1.m %in% c(1,2,3),
  Age.m < 25,
  Age.m > 5
) %>% 
  mutate(
    s8Del = ifelse(skip_to_render_inframe==8 & Category=="Deletion", 1, 0),
    s8naDel = ifelse(is.na(skip_to_render_inframe) | !(skip_to_render_inframe==8 & Category == "Deletion"), 0, 1),
    s44Del = ifelse(skip_to_render_inframe==44 & Category=="Deletion", 1, 0),
    s44naDel = ifelse(is.na(skip_to_render_inframe) | !(skip_to_render_inframe==44 & Category=="Deletion") , 0, 1),
    s45Del = ifelse(skip_to_render_inframe==45 & Category=="Deletion", 1, 0),
    s45naDel = ifelse(is.na(skip_to_render_inframe) | !(skip_to_render_inframe==45 & Category=="Deletion"), 0, 1),
    s50Del = ifelse(skip_to_render_inframe==50 & Category=="Deletion", 1, 0),
    s50naDel = ifelse(is.na(skip_to_render_inframe) | !(skip_to_render_inframe==50 & Category=="Deletion"), 0, 1),
    s51Del = ifelse(skip_to_render_inframe==51  & Category=="Deletion", 1, 0),
    s51naDel = ifelse(is.na(skip_to_render_inframe) | !(skip_to_render_inframe==51 & Category=="Deletion"), 0, 1),
    s52Del = ifelse(skip_to_render_inframe==52 & Category=="Deletion", 1, 0),
    s52naDel = ifelse(is.na(skip_to_render_inframe) | !(skip_to_render_inframe==52 & Category=="Deletion"), 0, 1),
    s53Del = ifelse(skip_to_render_inframe==53 & Category=="Deletion", 1, 0),
    s53naDel = ifelse(is.na(skip_to_render_inframe) | !(skip_to_render_inframe==53 & Category=="Deletion"), 0, 1),
    s55Del = ifelse(skip_to_render_inframe==55 & Category=="Deletion", 1, 0),
    s55naDel = ifelse(is.na(skip_to_render_inframe) | !(skip_to_render_inframe==55 & Category=="Deletion"), 0, 1),
    Dup = ifelse(Category=="Duplication", 1, 0),
    Nonsense = ifelse(Category=="Nonsense", 1, 0)
  )
#rought counts ... numbers are small
all.anymutation.noster %>% filter(Category=="Deletion") %>% select(skip_to_render_inframe) %>% table()

# none of these are significant (small N)
makeSurvPlot("s8Del", all.anymutation.noster, "8 skippable")       
makeSurvPlot("s8naDel", all.anymutation.noster, "8 skippable NA")  
makeSurvPlot("s44Del", all.anymutation.noster, "44 skippable")
makeSurvPlot("s44naDel", all.anymutation.noster, "44 skippable NA") 
makeSurvPlot("s45Del", all.anymutation.noster, "45 skippable")
makeSurvPlot("s45naDel", all.anymutation.noster, "45 skippable NA")
makeSurvPlot("s50Del", all.anymutation.noster, "50 skippable")
makeSurvPlot("s50naDel", all.anymutation.noster, "50 skippable NA")
makeSurvPlot("s51Del", all.anymutation.noster, "51 skippable")       
makeSurvPlot("s51naDel", all.anymutation.noster, "51 skippable NA") 
makeSurvPlot("s52Del", all.anymutation.noster, "52 skippable")
makeSurvPlot("s52naDel", all.anymutation.noster, "52 skippable NA")
makeSurvPlot("s53Del", all.anymutation.noster, "53 skippable")
makeSurvPlot("s53naDel", all.anymutation.noster, "53 skippable NA")
makeSurvPlot("s55Del", all.anymutation.noster, "55 skippable")
makeSurvPlot("s55naDel", all.anymutation.noster, "55 skippable NA")
makeSurvPlot("Dup", all.anymutation.noster, "Dupications")
makeSurvPlot("Nonsense", all.anymutation.noster, "Nonsense")




# ---- FIGURES ----

# Fig 1 survival: any mutation, + steroid
png(file="RESULTS/2016/201604/fig1_surv1.png")
layout(matrix(1:6,3,2))
makeSurvPlot("s8naDel", all.anymutation.ster, "8 skippable NA")  # Significant
makeSurvPlot("s44naDel", all.anymutation.ster, "44 skippable NA")  # Significant
makeSurvPlot("s45naDel", all.anymutation.ster, "45 skippable NA")
makeSurvPlot("s50naDel", all.anymutation.ster, "50 skippable NA")
makeSurvPlot("s51naDel", all.anymutation.ster, "51 skippable NA")  # Significant but compard to 04do.R where nonsense/dup NOT included, this is less significant
makeSurvPlot("s52naDel", all.anymutation.ster, "52 skippable NA")
dev.off()
png(file="RESULTS/2016/201604/fig1_surv2.png")
layout(matrix(1:6,3,2))
makeSurvPlot("s53naDel", all.anymutation.ster, "53 skippable NA")
makeSurvPlot("s55naDel", all.anymutation.ster, "55 skippable NA")
makeSurvPlot("Dup", all.anymutation.ster, "Dupications")
makeSurvPlot("Nonsense", all.anymutation.ster, "Nonsense")
dev.off()

# Fig2 Survival Detail fo E45del and E49_50del
png(file="RESULTS/2016/201604/fig2_survDetail.png")
layout(matrix(1:2,2,1))
makeSurvPlot("reg4545", all.anymutation.ster.44.detail, "44 skippable E45del")   
makeSurvPlot("reg4950", all.anymutation.ster.51.detail, "51 skippable E49_50del") # Significant
dev.off()

# Fig 2 with all subgroups
png(file="RESULTS/2016/201604/fig2_survDetailsubgroup.png")
layout(matrix(1:2,2,1))
makeSurvPlot("s44_subgroups", all.anymutation.ster.44.subgroups , "44 skippable by E45 subgroups")  
makeSurvPlot("s51_subgroups", all.anymutation.ster.51.subgroups , "51 skippable by subgroups")
dev.off()