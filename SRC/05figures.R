library(gridExtra)

# figures & data
#
#


# how many people in theses modules?
dim(gen2)
dim(ster2)
dim(muscle2)


# When was the module survey most recently completed?
gen2 %>% select(Survey.Time) %>% unlist(use.names=F) %>% strsplit( split=" ") %>% 
  lapply(function(x) x[1]) %>% unlist %>% as.Date("%m/%d/%Y") %>%
  format( "%Y") %>% table() %>% as.data.frame.table() %>% ggplot(aes_string(x=".", y="Freq")) +
  geom_bar(stat="identity")

p1 = histSurveyTime(gen2) + ggtitle("genetic")
p2 = histSurveyTime(muscle2) + ggtitle("muscle")
p3 = histSurveyTime(ster2) + ggtitle("steroid")

grid.arrange(p1,p2,p3, ncol=2, nrow=2)
g = arrangeGrob(p1,p2,p3,nrow=2)
ggsave("../RESULTS/2016/date_completed_survey.png", g)


# loss of ambulation, N=366
# mean 10, median 10.85
all.0 %>% filter(walking_01==0) %>% select(Q3.2.m) %>% summary
all.0 %>% filter(walking_01==0) %>% select(Q3.2.m) %>% ggplot(aes(x=Q3.2.m)) + 
  geom_histogram(binwidth=1, aes(fill=..count..)) +
  geom_vline(xintercept=10, col="green") + geom_vline(xintercept = 10.8, col="red") +
  ggtitle("Age at loss of ambulation") + xlab("years")

# --- Survival PLOTS for E8/45/50/51/52/53/55 (NO DUPLICATION or NONSENSE) + Steroids ----
# Survival plots exon skippable 8, 44, 45,50
layout(matrix(c(1,2,3,4), 2, 2))
makeSurvPlot("s8na", all.ster.44, "8 skippable NA")  # Significant
makeSurvPlot("s44na", all.ster.44, "44 skippable NA")  # Significant
makeSurvPlot("s45na", all.ster.44, "45 skippable NA")
makeSurvPlot("s50na", all.ster.44, "50 skippable NA")
# Survival plots exon skippable 51,52,53,55
makeSurvPlot("s51na", all.ster.44, "51 skippable NA")  # Significant
makeSurvPlot("s52na", all.ster.44, "52 skippable NA")
makeSurvPlot("s53na", all.ster.44, "53 skippable NA")
makeSurvPlot("s55na", all.ster.44, "55 skippable NA")




# ----FINAL Survival PLOT  Deletion, Duplication and Nonsense) + Steroids----
# from 04do_allmutation.R
# FINAL plots for publication: 
# E8 skippable
# E44 skippable
# E45 skippable
# E50 skippable
# E51 skippable
# E52 skippable
# E53 skippable
# E55 skippable
# Duplications
# Nonsense


pdf(file="RESULTS/2016/201604/fig1_surv1.pdf")
layout(matrix(1:6,3,2))
makeSurvPlot("s8naDel", all.anymutation.ster, "8 skippable NA")  # Significant
makeSurvPlot("s44naDel", all.anymutation.ster, "44 skippable NA")  # Significant
makeSurvPlot("s45naDel", all.anymutation.ster, "45 skippable NA")
makeSurvPlot("s50naDel", all.anymutation.ster, "50 skippable NA")
makeSurvPlot("s51naDel", all.anymutation.ster, "51 skippable NA")  # Significant but compard to 04do.R where nonsense/dup NOT included, this is less significant
makeSurvPlot("s52naDel", all.anymutation.ster, "52 skippable NA")
dev.off()
pdf(file="RESULTS/2016/201604/fig1_surv2.pdf")
layout(matrix(1:6,3,2))
makeSurvPlot("s53naDel", all.anymutation.ster, "53 skippable NA")
makeSurvPlot("s55naDel", all.anymutation.ster, "55 skippable NA")
makeSurvPlot("Dup", all.anymutation.ster, "Duplications")
makeSurvPlot("Nonsense", all.anymutation.ster, "Nonsense")
dev.off()


pdf(file="RESULTS/2016/201604/fig2_survDetail.pdf")
layout(matrix(1:2,1,2))
# E44 skippable, E45 deletions only, more mild
# compare E45dels to ALL OTHERS
makeSurvPlot("reg4545", all.anymutation.ster.44.detail, "44 skippable E45del")   # significant
# E51 skippable, E49_50dels are more severe! Significant P=0.00791
# compares E49_50dels to ALL OTHERS
makeSurvPlot("reg4950", all.anymutation.ster.51.detail, "51 skippable E49_50del") # Significant

dev.off()


#----- alt format Figure 1: all plots on same plot ----
# reshape2 a new variable
skip.cols = grep("naDel|Dup|Nonsense", names(all.anymutation.ster))
z = all.anymutation.ster %>% mutate(
  survPlotCat = case_when(
    s8naDel == 1 ~ "s8naDel",
    s44naDel == 1 ~ "s44naDel",
    s45naDel == 1 ~ "s45naDel",
    s50naDel == 1 ~ "s50naDel",
    s51naDel == 1 ~ "s51naDel",
    s52naDel == 1 ~ "s52naDel",
    s53naDel == 1 ~ "s53naDel",
    s55naDel == 1 ~ "s55naDel",
    Dup == 1 ~ "Dup",
    Nonsense==1 ~ "Nonsense"
  )
)

pdf(file="RESULTS/2016/201604/fig1alt_allPlotsTogether.pdf")
makeSurvPlot("survPlotCat", z, "all survival plots", conf.int=F) 
dev.off()

#----- alt format Figure 1: only 8,44,51 separate and all other as average // using RMS to plot----
z = all.anymutation.ster %>% mutate(
  survPlotCat = case_when(
    s8naDel == 1 ~ "s8naDel",
    s44naDel == 1 ~ "s44naDel",
    s51naDel == 1 ~ "s51naDel",
    TRUE  ~ "others"
  )
)

pdf(file="RESULTS/2016/201604/fig1alt_3PlotsPlusOthers.pdf")
makeSurvPlot("survPlotCat", z, "8/44/51/others merge") 
dev.off()


km.by.cat <- survfit(Surv(time_to_wheelchair, walking_01==0) ~ survPlotCat, data = z, conf.type = "log-log")
plot(km.by.cat, lty=1:10, col=1:10, main='title',
     xaxt="n", yaxt="n", bty="l", xlab="Age (years)", ylab="Percent Survival", mark.time=T, conf.int=TRUE)


# fix rms::survplot() use npsurv()
# RMS gives ugly conf intervals
library(rms)
km.by.catNp = npsurv(formula=Surv(time_to_wheelchair, walking_01==0) ~ survPlotCat, data = z)
survplot(km.by.catNp, conf.type="log-log", conf="bands", col=1:10, conf.int=0.99)
survplot(km.by.catNp, conf.type="log-log", conf="bars", col=1:10, conf.int=0.99) 
survplotp(km.by.catNp, conf.type="log-log", conf="bands") #interactive

#----- alt format Figure 1: only 8,44,51 separate and all other as average // using survminer----
library(survminer)
require(survival)

# remove duplicates? argh. some deletions can be solved by skipping 2 different exons...ignoring for now
dups = all.anymutation$Patient.ID[ duplicated(all.anymutation$Patient.ID)]
zz = z %>% filter(!Patient.ID %in% dups)

# create the Survfit
fit <- survfit(Surv(time_to_wheelchair, walking_01==0) ~ survPlotCat, data = z, conf.type = "log-log")
ggsurvplot(fit, data = z,
           risk.table=F,
           pval=TRUE,
           conf.int=0.99,
           break.time.by=3,
           risk.table.y.text.col=T,
           risk.table.y.text = F,
           xlim=c(4,20),
           conf.int.style="ribbon")

ggsave("fig1_condensed.pdf") # edit in illustrator to remove CI bands around 8/44/51
#ggsave("../RESULTS/2016/201604/fig1_condensed.pdf") #can't write to directory for some reason
ggsurvplot_add_all(fit, data = z)

#even after removing duplicates, its still significant. too hard to explain.
fit2 = survfit(Surv(time_to_wheelchair, walking_01==0)~survPlotCat, data = zz, conf.type="log-log")
ggsurvplot(fit2, data = zz,
           risk.table=F,
           pval=TRUE,
           conf.int=0.99,
           break.time.by=3,
           risk.table.y.text.col=T,
           risk.table.y.text = T,
           xlim=c(0,20),
           conf.int.style="ribbon")
