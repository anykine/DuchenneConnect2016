#######################################
# Project: DuchenneConnect2016
# Purpose: Analyze DuchenneConnect Data from October 2016 & October 2017
# Module: DuchenneConnect2016/SRC/rclean/104do.R
# Author: Richard T Wang
#

all.anymutation = all %>% dplyr::filter(
  MusGender == "Male",
  SteHealth.Reason.for.Registration == "Duchenne",
  MusDiagnosis == "Duchenne",
  MusCountry %in% c("AUSTRALIA", "AUSTRIA", "BELGIUM", "CANADA", "CHILE", "CZECH REPUBLIC", "DENMARK", "ESTONIA",
                 "FINALND", "FRANCE", "GERMANY", "GREECE", "HUNGARY", "ICELAND", "IRELAND", "ISREAL", "ITALY",
                 "JAPAN", "KOREA", "LATVIA", "LUXEMBOURG","MEXICO","NETHERLANDS","NEW ZEALAND","NORWAY",
                 "POLAND", "PORTUGAL", "SLOVAK REPUBLIC", "SLOVENIA", "SPAIN", "SWEDEN", "SWITZERLAND", "TURKEY",
                 "UNITED KINGDOM", "UNITED STATES"),
  GenCategory %in% c("Deletion", "Duplication", "Nonsense"),
  GenStart.Exon.Intron != "",
  GenEnd.Exon.Intron != "",
  # steroid 1=pred,2=dflz,3=stopped,4=never
  SteQ1 %in% c(1,2,3,4),
  #Walking 1=own,2=partial,3=WC,4=infant
  MusQ1 %in% c(1,2,3)
) %>% mutate(
  GenStart.Exon.Intron = gsub("^Ex", "", GenStart.Exon.Intron),
  GenEnd.Exon.Intron = gsub("^Ex", "", GenEnd.Exon.Intron),
  startstop = paste(GenStart.Exon.Intron, GenEnd.Exon.Intron, sep=",")) %>% 
  merge(skips, by.x="startstop", by.y="startstop", all.x=T) %>%
  mutate(
    walking_01 = ifelse(MusQ3.2 > 5 & !is.na(MusQ3.2), 0, 1)
  ) %>% 
  mutate(
    time_to_wheelchair = ifelse(walking_01==0, MusQ3.2, MusAge)
  )


all.anymutation.ster = all.anymutation %>% filter(
  # steroid 1=pred,2=dflz,3=stopped,4=never
  SteQ1 %in% c(1,2),
  #Walking 1=own,2=partial,3=WC,4=infant
  MusQ1 %in% c(1,2,3),
  MusAge < 25,
  MusAge > 5
) %>% 
  mutate(
    s8Del = ifelse(skip_to_render_inframe==8 & GenCategory=="Deletion", 1, 0),
    s8naDel = ifelse(is.na(skip_to_render_inframe) | !(skip_to_render_inframe==8 & GenCategory == "Deletion"), 0, 1),
    s44Del = ifelse(skip_to_render_inframe==44 & GenCategory=="Deletion", 1, 0),
    s44naDel = ifelse(is.na(skip_to_render_inframe) | !(skip_to_render_inframe==44 & GenCategory=="Deletion") , 0, 1),
    s45Del = ifelse(skip_to_render_inframe==45 & GenCategory=="Deletion", 1, 0),
    s45naDel = ifelse(is.na(skip_to_render_inframe) | !(skip_to_render_inframe==45 & GenCategory=="Deletion"), 0, 1),
    s50Del = ifelse(skip_to_render_inframe==50 & GenCategory=="Deletion", 1, 0),
    s50naDel = ifelse(is.na(skip_to_render_inframe) | !(skip_to_render_inframe==50 & GenCategory=="Deletion"), 0, 1),
    s51Del = ifelse(skip_to_render_inframe==51  & GenCategory=="Deletion", 1, 0),
    s51naDel = ifelse(is.na(skip_to_render_inframe) | !(skip_to_render_inframe==51 & GenCategory=="Deletion"), 0, 1),
    s52Del = ifelse(skip_to_render_inframe==52 & GenCategory=="Deletion", 1, 0),
    s52naDel = ifelse(is.na(skip_to_render_inframe) | !(skip_to_render_inframe==52 & GenCategory=="Deletion"), 0, 1),
    s53Del = ifelse(skip_to_render_inframe==53 & GenCategory=="Deletion", 1, 0),
    s53naDel = ifelse(is.na(skip_to_render_inframe) | !(skip_to_render_inframe==53 & GenCategory=="Deletion"), 0, 1),
    s55Del = ifelse(skip_to_render_inframe==55 & GenCategory=="Deletion", 1, 0),
    s55naDel = ifelse(is.na(skip_to_render_inframe) | !(skip_to_render_inframe==55 & GenCategory=="Deletion"), 0, 1),
    Dup = ifelse(GenCategory=="Duplication", 1, 0),
    Nonsense = ifelse(GenCategory=="Nonsense", 1, 0)
  )

makeSurvPlot("s8Del", all.anymutation.ster, "8 skippable")       # Significant
