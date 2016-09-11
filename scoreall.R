model.description = "None"

names = c(
"IS1", 
"IS2", 
"IS3", 
"IS4",
"CPCG0100",
"CPCG0183",
"CPCG0184",
"CPCG0196",
"CPCG0235",
"PCSI0023",
"PCSI0044",
"PCSI0046",
"PCSI0048",
"PCSI0072")

header.file <- paste(base.dir, "/predictions/header", sep="")
file.remove(header.file)

for (name in names) {
  
  vcfName <- name
  if (startsWith(name, "IS")) {
    vcfName <- paste("synthetic_", substr(name, nchar(name), nchar(name)), sep="")
  }
  
  clean.data <- loadCleanData(base.dir, name, load.vcf = FALSE)
  
  df.clean <- clean.data$df
  genomic.features.clean <- clean.data$genomic.features
  
  # Restart a clean run for an individual sample here
  
  env <- new.env()
  
  env$df <- df.clean
  env$genomic.features <- genomic.features.clean
  env$vcf.set = FALSE
  
  print(paste("Processing data frame for ", name, " with ", length(colnames(env$df)), " columns and ", nrow(env$df), " rows", sep=""))
  
  prepareData(name, votes.threshold, votes.rate.threshold)
  
  # Filter out columns that will not be used for prediction (like the unranked submissions)
  filterColumns(includeKey = TRUE, includeDv = FALSE)
  
  test <- env$df
  
  test[is.na(test)] <- 0
  
  model <- readRDS(file = paste(base.dir, "/model/latestboostingall.rds", sep = ""))
  
  predicted <- doPredictions(test, model, threshold = 0.5)
  
  predicted$Sample <- vcfName
  
  predicted.out <- select(predicted, CHROM, POS, Sample, PREDICTED_IS_VARIANT_BINARY, PREDICTED_IS_VARIANT)
  predicted.out <- rename(predicted.out, Predicted = PREDICTED_IS_VARIANT_BINARY)
  predicted.out <- rename(predicted.out, Probability = PREDICTED_IS_VARIANT)
  
  # Add one to conform to VCF format
  predicted.out$POS = predicted.out$POS + 1
  
  predicted.out$CHROM <- sapply(predicted.out$CHROM, function(chr.name) {
    chr.name <- as.character(chr.name)
    if (startsWith(chr.name, "chr")) {
      return (substr(chr.name, nchar("chr")+1, nchar(chr.name)))
    }
    else {
      return (chr.name)
    }
  })
  
  write.table(predicted.out, file = paste(base.dir, "/predictions/predictions-", name, ".csv", sep=""), sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names = FALSE)
  
  callers.used <- env$callers.used
  callers.used.string <- paste(callers.used, collapse=",")
  callers.used.string <- paste("##", vcfName, "_Pipelines=", callers.used.string, sep="")
  write(callers.used.string,file=header.file,append=TRUE)
  
  print(colnames(test))
  
}

write("#CHROM\tPOS\tSample\tPredicted\tProbability",file=header.file,append=TRUE)
