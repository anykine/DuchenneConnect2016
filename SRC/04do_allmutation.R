#----basic dataset wtih Deletion, Duplication, or Nonsense----
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



# E44 skippable: Compare E45del  or E45-54del Deletions) + Steroids to everyone else
# E45dels are milder! Significant! P=0.029
all.anymutation.ster.44.detail = all.anymutation.ster %>% mutate(
  reg4545 = ifelse(Start.Exon.Intron==45 & End.Exon.Intron==45 & skip_to_render_inframe==44 & Category=="Deletion", 1, 0),
  reg4554 = ifelse(Start.Exon.Intron==45 & End.Exon.Intron==54 & skip_to_render_inframe==44 & Category=="Deletion", 1, 0)
)
makeSurvPlot("reg4545", all.anymutation.ster.44.detail, "44 skippable E45del")   # significant
makeSurvPlot("reg4554", all.anymutation.ster.44.detail, "44 skippable E45del")   

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

# ---- Cox Regression ----
# according to Ake, we should create 1 variable for mutation type with multiple levels
# let's also relevel steroids so that no steroid group is the reference


#----Cox regression----
all.anymutation.cox = all.anymutation %>% mutate(
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
# collapse all mutation types into one factor "muttype"
all.anymutation.cox2 = all.anymutation.cox %>% mutate(
  muttype=ifelse(s8naDel==1, "s8",
                 ifelse(s44naDel==1, "s44",
                        ifelse(s45naDel==1, "s45", 
                               ifelse(s50naDel==1, "s50",
                                      ifelse(s51naDel==1, "s51",
                                             ifelse(s52naDel==1,"s52",
                                                    ifelse(s53naDel==1, "s53",
                                                           ifelse(s55naDel==1, "s55",
                                                                  ifelse(Dup==1, "dup",
                                                                         ifelse(Nonsense==1, "nonsense", "other"))))))))))
)
# make NO STEROID, the reference level
all.anymutation.cox2$Q1.s = factor(all.anymutation.cox2$Q1.s)
all.anymutation.cox2$Q1.s = relevel(all.anymutation.cox2$Q1.s, 4)
# collapsed mutation type to one variable
coxfit2 = coxph(Surv(all.anymutation.cox2$time_to_wheelchair, all.anymutation.cox2$walking_01==0) ~ 
                  as.factor(all.anymutation.cox2$Q1.s)+ as.factor(all.anymutation.cox2$muttype))
summary(coxfit2)                 
cox.zph(coxfit2) # pval small means cox assumptions violated...due to deflazacort?
plot(cox.zph(coxfit2))  #diagnostic, is line flat? No.

# add time variable age:mutation interaction
coxfit3 = coxph(Surv(all.anymutation.cox2$time_to_wheelchair, all.anymutation.cox2$walking_01==0) ~ 
                  as.factor(all.anymutation.cox2$Q1.s)+ as.factor(all.anymutation.cox2$muttype) + 
                  as.factor(all.anymutation.cox2$muttype):all.anymutation.cox2$Age)
summary(coxfit3)   
# stratify by setroid, (will calc different baselines for pred/dflz/none/prev)
coxfit4 = coxph(Surv(all.anymutation.cox2$time_to_wheelchair, all.anymutation.cox2$walking_01==0) ~ 
                  strata(as.factor(all.anymutation.cox2$Q1.s))+ as.factor(all.anymutation.cox2$muttype) ) 
                  #as.factor(all.anymutation.cox2$muttype):all.anymutation.cox2$Age)
summary(coxfit4)   

# stratify by setroid, (will calc different baselines for pred/dflz/none/prev)
coxfit5 = coxph(Surv(all.anymutation.cox2$time_to_wheelchair, all.anymutation.cox2$walking_01==0) ~ 
                  as.factor(all.anymutation.cox2$Q1.s)+ strata(as.factor(all.anymutation.cox2$muttype) ) )
#as.factor(all.anymutation.cox2$muttype):all.anymutation.cox2$Age)
summary(coxfit5)
plot(cox.zph(coxfit5))

# try to adjust w/ time variying cohort, strata(), other ideas?

# ---- ANOVA ----
# compare mean/variance of E44 and other skippables
# plot means (steroid effect NOT ACCOUNTED)
all.anymutation.cox2 %>% filter(Q1.s %in% c(1,2)) %>% ggplot(aes(x=muttype,y=time_to_wheelchair)) +
  geom_boxplot(fill = "grey80", colour = "blue") +
  scale_x_discrete() + xlab("Mutation group") +
  ylab("Age LOA") + ggtitle("age LOA by mutation (currently on steroids)")
#ggsave("RESULTS/2016/201604/boxplot_LOA_mutation_steroid1.png")


# 1-way anova, using lm()
anova.model = lm(time_to_wheelchair ~ muttype, data = all.anymutation.cox2)
summary(anova.model)
anova(anova.model)
confint(anova.model)
# 1-way anova, using aov()
anova.model2 = aov(time_to_wheelchair ~ as.factor(muttype), data = all.anymutation.cox2)
summary(anova.model2)
#pairwise.t.test(all.anymutation.cox2$muttype, all.anymutation.cox2$time_to_wheelchair,p.adj="none",na.rm=T)

# 2 factor anova
aov2 = aov(time_to_wheelchair ~ as.factor(muttype) + as.factor(Q1.s) + as.factor(muttype)*as.factor(Q1.s), data = all.anymutation.cox2)
summary(aov2 )
print(model.tables(aov2, "means"), digits=3) # print means of each grouup and number of samples per cell
res = TukeyHSD(aov2, c("as.factor(muttype)"), ordered=T)

# nice boxplot
png(file="RESULTS/2016/201604/boxplot_LOA_mutation_steroid.png",width=1024,height=800)
p=boxplot(time_to_wheelchair ~ as.factor(muttype) + as.factor(Q1.s), data = all.anymutation.cox2)
boxplot(time_to_wheelchair ~ as.factor(muttype) + as.factor(Q1.s), data = all.anymutation.cox2, axes=F, col=brewer.pal(11,"Set3"))
axis(1, at=1:44, labels=p$names, las=2)
axis(2, at=1:25, las=1)
title("Age LOA by mutation x steroid", ylab="age LOA")
#dev.off()

# 1way anova, ster+ only
# the only sig differences are between exon8 skippables and other groups
all.anymutation.cox2.ster = all.anymutation.cox2 %>% filter(Q1.s %in% c(1,2))
anova.model2ster = aov(time_to_wheelchair ~ as.factor(muttype), data = all.anymutation.cox2.ster)
summary(anova.model2ster)
pairwise.t.test(as.factor(all.anymutation.cox2.ster$muttype), all.anymutation.cox2.ster$time_to_wheelchair,p.adj="none")
TukeyHSD(anova.model2ster, "as.factor(muttype)")

# (old)Cox regression (old but should be equivalent to cox2 above)
# significant for steroids, E8 skippable, E44 skippable. Not significant for E51 skippables
coxfit = coxph(Surv(all.anymutation.cox$time_to_wheelchair, all.anymutation.cox$walking_01==0) ~ as.factor(all.anymutation.cox$Q1.s)+
                 as.factor(all.anymutation.cox$s8naDel) + 
                 as.factor(all.anymutation.cox$s44naDel) + 
                 as.factor(all.anymutation.cox$s45naDel) + 
                 as.factor(all.anymutation.cox$s50naDel) + 
                 as.factor(all.anymutation.cox$s51naDel) + 
                 as.factor(all.anymutation.cox$s52naDel) + 
                 as.factor(all.anymutation.cox$s53naDel) + 
                 as.factor(all.anymutation.cox$s55naDel) + 
                 as.factor(all.anymutation.cox$Dup) + 
                 as.factor(all.anymutation.cox$Nonsense ))
summary(coxfit)


tmp1=all.anymutation.cox %>% filter(
  # steroid 1=pred,2=dflz,3=stopped,4=never
  Q1.s %in% c(1,2),
  #Walking 1=own,2=partial,3=WC,4=infant
  Q1.m %in% c(1,2,3),
  Age.m < 25,
  Age.m > 5,
  skip_to_render_inframe %in% c(8,44,45,50,51,52,53,55),
  Category %in% c("Deletion", "Duplication", "Nonsesnse"
                  ))
tmp1coxfit = coxph(Surv(tmp1$time_to_wheelchair, tmp1$walking_01==0) ~ as.factor(tmp1$Q1.s)+
                 as.factor(tmp1$s8naDel) + 
                 as.factor(tmp1$s44naDel) + 
                 as.factor(tmp1$s45naDel) + 
                 as.factor(tmp1$s50naDel) + 
                 as.factor(tmp1$s51naDel) + 
                 as.factor(tmp1$s52naDel) + 
                 as.factor(tmp1$s53naDel) + 
                 as.factor(tmp1$s55naDel) )
                 #as.factor(tmp1$Dup) )
                 #as.factor(tmp1$Nonsense ))
                 
summary(tmp1coxfit)


# basic status for table1
all.anymutation.cox %>% dim
all.anymutation.cox %>% count(s8naDel, s44naDel,Nonsense) 


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