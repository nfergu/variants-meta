model.description = "None"

names = c("IS1", "IS2", "IS3", "IS4")  

env.results <- new.env()

env.results$results <- data.frame(matrix(ncol = 0, nrow = length(names)))
env.results$results$NAME <- NA
env.results$results$ACCURACY <- NA

i <- 1

for (name in names) {
  
  vcfName <- paste("synthetic_", substr(name, nchar(name), nchar(name)), sep="")
  
  clean.data <- loadCleanData(base.dir, name)
  
  df.clean <- clean.data$df
  genomic.features.clean <- clean.data$genomic.features
  vcf.clean <- clean.data$vcf
  
  # Restart a clean run for an individual sample here
  
  env <- new.env()
  
  env$df <- df.clean
  env$genomic.features <- genomic.features.clean
  env$vcf <- vcf.clean
  env$vcf.set = TRUE
  
  prepareData(name, votes.threshold, votes.rate.threshold)
  
  #outputStats(base.dir, name, num.variants)
  
  # Filter out columns that will not be used for prediction (like the unranked submissions)
  filterColumns(TRUE)
  
  test <- env$df
  
  test[is.na(test)] <- 0
  
  model <- readRDS(file = paste(base.dir, "/model/latestboostingall.rds", sep = ""))
  
  print(paste("RESULTS FOR ", name, sep=""))
  
  predicted <- validateModel(test, model, paste("TEST DATA for ", name, sep=""), name, 
                             num.variants, i, threshold = 0.5)

  predicted$Sample <- vcfName
  
  predicted.out <- select(predicted, CHROM, POS, Sample, PREDICTED_IS_VARIANT_BINARY, PREDICTED_IS_VARIANT)
  predicted.out <- rename(predicted.out, Predicted = PREDICTED_IS_VARIANT_BINARY)
  predicted.out <- rename(predicted.out, Probability = PREDICTED_IS_VARIANT)
  
  # Add one to conform to VCF format
  predicted.out$POS = predicted.out$POS + 1
  
  write.table(predicted.out, file = paste(base.dir, "/", name, "/predictions.csv", sep=""), sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names = FALSE)

  print(colnames(test))
    
  i <- i+1
  
}

print(paste("Results for ", model.description, sep=""))
print(env.results$results)
print("MEDIAN:")
print(median(env.results$results$ACCURACY))