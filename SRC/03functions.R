# ----- functions ------

# Return number of rows per person
MultipleRowsForPatient = function(data, returnData = F) {
  tab = data %>% group_by(Patient.ID) %>% summarize(count=n()) %>%
    filter(count > 1)
  if (returnData == TRUE) {
    tab
  } else {
    dim(tab)
  }
}

# Get rows with >1 row per patient
checkUniquePatientsPerRow = function(data, myKey="Patient.ID") {
  data %>% group_by_(myKey) %>% summarize(count=n()) %>% filter(count>1) %>% as.data.frame
}

# draw histogram of survey time for each module data set
# eg histSurveyTime(gen2)
histSurveyTime = function(dt, column="Survey.Time") {
  dt %>% select_(column) %>% unlist(use.names=F) %>% strsplit( split=" ") %>% 
    lapply(function(x) x[1]) %>% unlist %>% as.Date("%m/%d/%Y") %>%
    format( "%Y") %>% table() %>% as.data.frame.table() %>% ggplot(aes_string(x=".", y="Freq")) +
    geom_bar(stat="identity") + theme(axis.text.x = element_text(angle=90))
}


# plot survival curve
# makeSurvPlot("strat", all.ster, "Steroid effect")
# data is a data.table
makeSurvPlot = function(strat, data, title,...) {
  require(survival)
  formula = eval(
    parse(text = paste0("Surv(time_to_wheelchair, walking_01==0) ~", strat))
  )
  plot(survfit(formula, data = data, conf.type="log-log"), lty=1:10, col=1:10, main=title,
       xaxt="n", yaxt="n", bty="l", xlab="Age (years)", ylab="Percent Survival", mark.time=T,...)
  axis(1, las=1)
  axis(2, las=2)
  fit.diff = survdiff(formula, data = data)
  fit.diff
  sf = survfit(formula, data = data, conf.type="log-log")
  print(sf)
  
  pval <- signif( 1-pchisq(fit.diff$chisq, length(fit.diff$n)-1), digits=3)
  text(3, 0.6, paste0("p=",pval))
  legend("bottomleft", bty="n",
         legend=apply(as.data.frame(table(eval(data[[strat]]))), 1, paste, collapse=" n="),
         lty=1:10, col=1:10)
  print(fit.diff)
}

makeSurvPlotNoStrat = function(data, title) { #not working
  require(survival)
  formula = eval(
    parse(text = "Surv(time_to_wheelchair, walking_01==0)~1")
  )
  plot(survfit(formula, data = data, conf.type="log-log"), lty=1:10, col=1:10, main=title,
       xaxt="n", yaxt="n", bty="l", xlab="Age (years)", ylab="Percent Survival", mark.time=T)
  axis(1, las=1)
  axis(2, las=2)
  fit.diff = survdiff(formula, data = data)
  fit.diff
  sf = survfit(formula, data = data, conf.type="log-log")
  print(sf)
  
  pval <- signif( 1-pchisq(fit.diff$chisq, length(fit.diff$n)-1), digits=3)
  text(3, 0.6, paste0("p=",pval))
  # legend("bottomleft", bty="n",
  #        legend=apply(as.data.frame(table(eval(data[[strat]]))), 1, paste, collapse=" n="),
  #        lty=1:10, col=1:10)
  print(fit.diff)
}
