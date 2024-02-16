library("devtools")
#requires dependency Cmake
#devtools::install_github("NIEHS/ToxicR")
library("ToxicR")
library("dplyr")
library("PMCMRplus")
library(parallel)
library("pbmcapply")
#bmdoutput <- read.table("~/storage/biomarker_input_.txt",
#                        header=FALSE, sep = "\t", row.names = 1)
bmdoutput <- read.table("~/storage/ToxicR mod/data for BMDExpress log2 CPM Males.txt",
                        header=FALSE, sep = "\t", row.names = 1)
#
#bmdoutput <- bmdoutput[-3,-c(1:3,5:9,11:15,17:18)]
#
bmdoutput <- as.data.frame(t(bmdoutput))
bmdoutput <-bmdoutput[,-1]
bmdoutput[-1] <- lapply(bmdoutput[-1], function(x) {
  as.numeric(as.character(x))
})

dose_column <-as.numeric(bmdoutput[,1]) 
start.time <- Sys.time()
Func <- function(resp_bygene_column){ 
  trendtestres<- williamsTest(resp_bygene_column ~ dose_column, alternative = c("greater", "less"))
  trendestresdf <- data_frame("pass/fail value" = trendtestres$crit.value[,1] - trendtestres$statistic[,1]) 
  #only generate model when significant based on williams test 
  if (any(trendestresdf$`pass/fail value` <0)){ 
    model<-single_continuous_fit(dose_column,resp_bygene_column)
    #capture p-value output from ToxicR for goodness of fit
    modeloutput<- capture.output(summary(model))
    all_numbers_mean_adequate_pvalue <- lapply(modeloutput[9:10], function(line) {
      as.numeric(tail(regmatches(line, gregexpr("\\b\\d+\\.\\d+\\b", line))[[1]],n=1))
    })
    #only return model if goodness of fit is significant 
    if(!any(as.data.frame(all_numbers_mean_adequate_pvalue) >0.05)){ 
      return(as.numeric(model$bmd[1]))
    }
  }
}

result <- pbmclapply(bmdoutput[,2:ncol(bmdoutput)], function(col) {
  return(Func(col))
})

#structuring of results to remove NULLs and have rows for output
signficant_BMD_table <- t(as.data.frame(Filter(Negate(is.null), result)))
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken
write.table(signficant_BMD_table, file ="~/storage/toxicR output.txt", sep=" ", col.names = TRUE)
