library("devtools")
#requires dependency Cmake
#devtools::install_github("NIEHS/ToxicR")
library("ToxicR")
library("dplyr")
library("PMCMRplus")
library(parallel)
library("pbmcapply")
library(foreach)
library(doParallel)

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
  trendtestresgreater<- williamsTest(resp_bygene_column ~ dose_column, alternative = "greater")
  trendtestresdfgreater <- data_frame("pass/fail value" = trendtestresgreater$crit.value[,1] - trendtestresgreater$statistic[,1]) 
  
  trendtestresless<- williamsTest(resp_bygene_column ~ dose_column, alternative = "less")
  trendtestresdfless <- data_frame("pass/fail value" = trendtestresgreater$crit.value[,1] - trendtestresgreater$statistic[,1])

  trendtestresdf <- rbind(trendtestresdfgreater, trendtestresdfless)
  
    #only generate model when significant based on williams test 
  if (any(trendtestresdf$`pass/fail value` <0)){ 
    
    #cl <- makeCluster(detectCores())
    #registerDoParallel(cl)
    #clusterEvalQ(cl, {
    #  library("ToxicR")
    #})
    #clusterExport(cl, c("dose_column"))
    modeltypes = c("hill","exp-3","exp-5","power","polynomial")
    models <- foreach(model_type = modeltypes) %dopar% {
      single_continuous_fit(dose_column, resp_bygene_column, model_type = model_type)
    }
    #stopCluster(cl)
    
    
    
    #runmodels <- function(modeltype){ 
    #  single_continuous_fit(dose_column, resp_bygene_column, model_type = modeltype)
    #}
    #models <- pbmclapply(modeltypes,runmodels, mc.cores = detectCores())
    
    modeloptions <- sapply(models, function(model) model$maximum)
    index_of_best <- which.max(modeloptions)
    model <- sapply(models[[index_of_best]], function(model) model)
    #capture p-value output from ToxicR for goodness of fit
    modeloutput<- capture.output(summary(model))
    all_numbers_mean_adequate_pvalue <- lapply(modeloutput[9:10], function(line) {
      as.numeric(tail(regmatches(line, gregexpr("\\b\\d+\\.\\d+\\b", line))[[1]],n=1))
    })
    #only return model if goodness of fit is significant 
    if(!any(as.data.frame(all_numbers_mean_adequate_pvalue) >0.05)){ 
      return(list(as.numeric(model$bmd[1]),model$model))
    }
  }
}

result <-pbmclapply(bmdoutput[,2:ncol(bmdoutput)], Func, mc.cores = detectCores())

filtered_list <- Filter(Negate(is.null), result)
BMDs <- lapply(filtered_list, function(gene) as.numeric(gene[1]))
modeltypesresults <- lapply(filtered_list, function(gene) as.character(gene[2]))
#structuring of results to remove NULLs and have rows for output
significant_BMD_table <- data.frame("Gene" = names(filtered_list), 
                                    "BMD" = unlist(BMDs),
                                    "Modeltype" = unlist(modeltypesresults),row.names = NULL)

end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken
write.table(significant_BMD_table, file ="~/storage/toxicR output.txt", sep=" ", col.names = TRUE)
