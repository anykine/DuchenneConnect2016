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