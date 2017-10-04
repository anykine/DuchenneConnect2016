# exploratory questions and answers

# ---- what % of E44 skippable are E45del ----
# 50/78 = 64% of E44 skippable are E45del in steroid+
all.anymutation.ster %>% filter(skip_to_render_inframe==44) %>% dim
all.anymutation.ster %>% filter(skip_to_render_inframe==44, startstop=="45,45") %>% dim()

# ---- what % of E8 skippable are E-7del ----
# 50/78 = 64% of E44 skippable are E45del in steroid+
all.anymutation.ster %>% filter(skip_to_render_inframe==8) %>% dim()
all.anymutation.ster %>% filter(skip_to_render_inframe==8, startstop=="3,7") %>% dim()

# ---- what % of E51 skippable are E49-50del ----
# E49-50del is 22% of E51 skippables
all.anymutation.ster %>% filter(skip_to_render_inframe==51) %>% dim()
all.anymutation.ster %>% filter(skip_to_render_inframe==51, startstop=="49,50") %>% dim()

# EDA

# hist of WCage E45-50del, is it gaussian?
all.anymutation.ster %>% filter(startstop=="45,50") %>% 
  ggplot(aes(time_to_wheelchair)) + geom_histogram(binwidth=0.5)+
  ylab("WC age") + ggtitle("E45-50 dels")


# hist of WCage E45del, is it bimodal?
# TODO: this includes del/dup/etc. color by mutation type?
all.anymutation.ster %>% filter(startstop=="45,45", Category=="Deletion") %>% 
  ggplot(aes(time_to_wheelchair)) + geom_histogram(binwidth=0.5) +
  ylab("WC age") + ggtitle("E45 dels")

# try histogram by deletion category
all.anymutation.ster %>% filter(skip_to_render_inframe==44) %>% 
  ggplot(aes(time_to_wheelchair)) + geom_histogram(binwidth=0.5)+
  ylab("WC age") + ggtitle("44 skippables")

all.anymutation.ster %>% filter(skip_to_render_inframe==44) %>% 
  ggplot() + geom_density(aes(x=time_to_wheelchair, color=Category))+
 ggtitle("44 skippables")

# histogram by start,stop
all.anymutation.ster %>% filter(skip_to_render_inframe==44) %>% 
  ggplot() + geom_density(aes(x=time_to_wheelchair, color=startstop))+
   ggtitle("44 skippables")
