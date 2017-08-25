# ---- Cox Regression ----
# according to Ake, we should create 1 variable for mutation type with multiple levels
# let's also relevel steroids so that no steroid group is the reference

all.anymutation.coxNew = all.anymutation %>% filter(
  Category=="Deletion" | Category=="Duplication" | Category=="Nonsense"
) %>% dplyr::mutate( 
  mutationType = case_when(
    .$skip_to_render_inframe==8 ~ "s8",
    .$skip_to_render_inframe==44 ~ "s44",
    .$skip_to_render_inframe==45 ~ "s45",
    .$skip_to_render_inframe==50 ~ "s50",
    .$skip_to_render_inframe==51 ~ "s51",
    .$skip_to_render_inframe==52 ~ "s52",
    .$skip_to_render_inframe==53 ~ "s53",
    .$skip_to_render_inframe==55 ~ "s55")
)  %>% dplyr::mutate( mutationType = ifelse(Category=="Duplication", "Dup", mutationType),  #weird workaround
                      mutationType = ifelse(Category=="Nonsense", "Nonsense", mutationType)) 

# check
all.anymutation.coxNew %>% select(mutationType) %>% as.data.frame %>% table()
all.anymutation.coxNew %>% filter(mutationType=="s44") %>% select(mutationType, Category) %>% table()

# steroid status, mutation must be factors
all.anymutation.coxNew = all.anymutation.coxNew %>% mutate(steroidStatus = factor(all.anymutation.coxNew$Q1.s))
all.anymutation.coxNew = all.anymutation.coxNew %>% mutate(mutationType = factor(mutationType))
# make no steroids the reference level
all.anymutation.coxNew$steroidStatus = relevel(all.anymutation.coxNew$steroidStatus, 4)
# make nonsense the reference level
all.anymutation.coxNew$mutationType = relevel(all.anymutation.coxNew$mutationType, "Nonsense")

# collapsed mutation type to one variable
coxfitNew = coxph(Surv(all.anymutation.coxNew$time_to_wheelchair, all.anymutation.coxNew$walking_01==0) ~ 
                    as.factor(all.anymutation.coxNew$steroidStatus)+ as.factor(all.anymutation.coxNew$mutationType))
summary(coxfitNew)
survfit(coxfitNew)
#diagnostic
cox.zph(coxfitNew) # pval small means cox assumptions violated...due to deflazacort?
plot(cox.zph(coxfitNew))  #diagnostic, is line flat? 

# plot residuals. Not working, try viewing videos
#https://github.com/ryandata/Survival
#res_martingale<-residuals(coxfitNew, type="martingale")
#scatter.smooth(all.anymutation.coxNew$steroidStatus,res_martingale)
# Plot martingale residuals. See
# http://www.math.wustl.edu/~jmding/math434/R_model_diag.R

#----cox regression, split E44 skippable into s44E45del and s44notE45del-----


all.anymutation.coxNew2 = all.anymutation.coxNew %>% 
  mutate(mutationType2 = ifelse(
    Category=="Deletion" & skip_to_render_inframe==44 & startstop=="45,45", 
    "s44E45del",
      ifelse(
        Category=="Deletion" & skip_to_render_inframe==44 & startstop!="45,45",
        "s44notE45del",
        as.character(mutationType)
      )
    )
)
# lost our labels for some reason
#check
all.anymutation.coxNew2 %>% select(mutationType2) %>% table()
all.anymutation.coxNew2$mutationType2 %>% class
# must be factor
all.anymutation.coxNew2 = all.anymutation.coxNew2 %>% mutate(mutationType2 = factor(mutationType2))
# relevel mutationType2
all.anymutation.coxNew2$mutationType2 = relevel(all.anymutation.coxNew2$mutationType2, "Nonsense")

coxfitNew2 = coxph(Surv(all.anymutation.coxNew2$time_to_wheelchair, all.anymutation.coxNew2$walking_01==0) ~ 
                    as.factor(all.anymutation.coxNew2$steroidStatus)+ as.factor(all.anymutation.coxNew2$mutationType2))
summary(coxfitNew2)
survfit(coxfitNew2)

#----cox regression, split E51 skippable into E49_50del and everyone else -----


all.anymutation.coxNew3 = all.anymutation.coxNew2 %>% 
  mutate(mutationType3 = ifelse(
    Category=="Deletion" & skip_to_render_inframe==51 & startstop=="49,50", 
    "s51E4951del",
    ifelse(
      Category=="Deletion" & skip_to_render_inframe==51 & startstop!="49,50",
      "s51notE4951del",
      as.character(mutationType2)
    )
  )
)
all.anymutation.coxNew3 %>% select(mutationType3) %>% table()
all.anymutation.coxNew3$mutationType3 %>% class

# must be factor
all.anymutation.coxNew3 = all.anymutation.coxNew3 %>% mutate(mutationType3 = factor(mutationType3))
# relevel mutationType2
all.anymutation.coxNew3$mutationType3 = relevel(all.anymutation.coxNew3$mutationType3, "Nonsense")

coxfitNew3 = coxph(Surv(all.anymutation.coxNew3$time_to_wheelchair, all.anymutation.coxNew3$walking_01==0) ~ 
                     as.factor(all.anymutation.coxNew3$steroidStatus)+ as.factor(all.anymutation.coxNew3$mutationType3))
summary(coxfitNew3)
survfit(coxfitNew3)

#----cox regression, split E8 skippable into E3_7del and everyone else -----
# this one start with all.anymutation.coxNew2 (e44 split because e51 is not informative)
# this one CONFIRMS that steroids, E45del, and E3_7del are significant covariates 
# FINAL
all.anymutation.coxNew2 %>% filter(Category=="Deletion", skip_to_render_inframe==8) %>% select(startstop) %>% table()

all.anymutation.coxNew4 = all.anymutation.coxNew2 %>% 
  mutate(mutationType4 = ifelse(
    Category=="Deletion" & skip_to_render_inframe==8 & startstop=="3,7", 
    "s8E37del",
    ifelse(
      Category=="Deletion" & skip_to_render_inframe==8 & startstop!="3,7",
      "s8notE37del",
      as.character(mutationType2)
    )
  )
  )
all.anymutation.coxNew4 %>% select(mutationType4) %>% table()
all.anymutation.coxNew4$mutationType4 %>% class
# must be factor
all.anymutation.coxNew4 = all.anymutation.coxNew4 %>% mutate(mutationType4 = factor(mutationType4))
# relevel mutationType2
all.anymutation.coxNew4$mutationType4 = relevel(all.anymutation.coxNew4$mutationType4, "Nonsense")

coxfitNew4 = coxph(Surv(all.anymutation.coxNew4$time_to_wheelchair, all.anymutation.coxNew4$walking_01==0) ~ 
                     as.factor(all.anymutation.coxNew4$steroidStatus)+ as.factor(all.anymutation.coxNew4$mutationType4))
summary(coxfitNew4)
survfit(coxfitNew4)


#----cox regression on steroid only patients, what mutations are signif?----
all.anymutation.coxNewSter = all.anymutation.coxNew %>% filter(
  Q1.s %in% c(1,2)
)
coxfitNewSter = coxph(Surv(all.anymutation.coxNewSter$time_to_wheelchair, all.anymutation.coxNewSter$walking_01==0) ~ 
                        as.factor(all.anymutation.coxNewSter$mutationType))
summary(coxfitNewSter)
survfit(coxfitNewSter)


#----Cox regression (previous attemp but probably equivalent to the above all.anymutation.coxNew)----
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
# stratify by steroid, (will calc different baselines for pred/dflz/none/prev)
coxfit4 = coxph(Surv(all.anymutation.cox2$time_to_wheelchair, all.anymutation.cox2$walking_01==0) ~ 
                  strata(as.factor(all.anymutation.cox2$Q1.s))+ as.factor(all.anymutation.cox2$muttype) ) 
#as.factor(all.anymutation.cox2$muttype):all.anymutation.cox2$Age)
summary(coxfit4)   

# stratify by steroid, (will calc different baselines for pred/dflz/none/prev)
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