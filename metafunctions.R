library(VariantAnnotation)
library(dplyr)
library(glmnet)
library(pROC)
library(data.table)
library(car)
library(caret)
library(Causata)
library(xgboost)
library(dummies)
library(gdata)

base.dir <- "/Users/nferguson/dream/variants-meta"

dvname <- 'IS_VARIANT'

num.variants <- 0
# These thresholds were calculated based on ALL
votes.threshold <- 38
votes.rate.threshold <- 0.27

loadCleanData <- function(base.dir, name, load.vcf = TRUE) {
  df.clean <- read.delim(paste(base.dir, "/", name, "/SNVCalls_", name, ".txt", sep=""), header = TRUE)
  genomic.features.clean <- read.delim(paste(base.dir, "/", name, "/GenomicFeatures_", name, ".txt", sep=""), header = TRUE)
  if (load.vcf) {
    vcf.clean <- readVcf(paste(base.dir, "/", name, "/truth.vcf", sep = ""), name)
    return(list("df" = df.clean, "genomic.features" = genomic.features.clean, "vcf" = vcf.clean))
  }
  else {
    return(list("df" = df.clean, "genomic.features" = genomic.features.clean))
  }
  
}

calculateVotesAndConsensus <- function(threshold, rate.threshold) {
  
  names <- colnames(env$df)[grep("X.*", colnames(env$df))]
  
  submissions <- env$df[,names]
  env$df$votes = rowSums(submissions, na.rm = TRUE)
  
  submission.count <- length(colnames(submissions))
  
  env$df <- mutate(env$df, naive.consensus = votes >= (submission.count / 2))
  env$df$naive.consensus <- as.numeric(env$df$naive.consensus)
  env$df <- mutate(env$df, optimal.consensus = votes >= threshold)
  env$df$optimal.consensus <- as.numeric(env$df$optimal.consensus)
  
  env$df$submission.count <- submission.count
  env$df$vote.rate <- env$df$votes / env$df$submission.count
  
  env$df <- mutate(env$df, naive.consensus.from.rate = vote.rate >= 0.5)
  env$df$naive.consensus.from.rate <- as.numeric(env$df$naive.consensus.from.rate)
  env$df <- mutate(env$df, optimal.consensus.from.rate = vote.rate >= rate.threshold)
  env$df$optimal.consensus.from.rate <- as.numeric(env$df$optimal.consensus.from.rate)
  
}

joinData <- function() {
  if (env$vcf.set) {
    rowranges <- as.data.frame(rowRanges(env$vcf))
    rowranges <- rename(rowranges, CHROM = seqnames)
    rowranges <- rename(rowranges, POS = start)
    # BED is zero-based, but VCF is one-based, so subtract one here 
    # reference: https://genome.ucsc.edu/FAQ/FAQformat.html#format1 and https://www.biostars.org/p/18357/
    rowranges$POS <- rowranges$POS - 1
    variants <- select(rowranges, CHROM, POS)
    # Create a new column called IS_VARIANT containing the ground truth. 
    # Everything in the VCF file has a value of 1 for this
    variants$IS_VARIANT = 1
    env$df <- left_join(env$df, variants, by = c("CHROM", "POS"))
    env$df$IS_VARIANT[is.na(env$df$IS_VARIANT)] <- 0
    not.present <- anti_join(rowranges, env$df, by = c("CHROM", "POS"))
    #print("NOT PRESENT IN DATA:")
    #print(not.present)
  }
  env$df <- left_join(env$df, env$genomic.features, by = c("CHROM", "POS"))
}

calculateAll <- function(df, FUN, truthCol = "IS_VARIANT") {
  result.vec <- sapply(colnames(df), function(col.name) {
    if (grepl("X.*", col.name)) {
      c(col.name, FUN(df[[col.name]], df[[truthCol]]))
    }
    else {
      NA
    }
  })
  result.vec <- result.vec[!is.na(result.vec)]
  result.df <- data.frame(matrix(unlist(result.vec), nrow=length(result.vec), byrow = TRUE))
  result.df$X2 = as.numeric(as.character(result.df$X2))
  result.df <- rename(result.df, COLUMN = X1)
  result.df <- rename(result.df, VALUE = X2)
  result.df <- filter(result.df, !is.na(VALUE))
  result.df
}

balancedAccuracy <- function(prediction, truth) {
  truth <- as.factor(truth)
  prediction <- as.factor(prediction)
  sensitivity = sensitivity(prediction, truth)
  specificity = specificity(prediction, truth)
  (sensitivity + specificity) / 2.0
}

balancedAccuracyNew <- function(prediction, truth, trurecs) {
  return (f1Score(prediction, truth, trurecs))
}

balancedAccuracyFromEvaluatorScript <- function(prediction, truth, trurecs) {
  sensitivity <- sensitivityNew(prediction, truth, trurecs)
  specificity <- specificityNew(prediction, truth)
  return ((sensitivity + specificity) / 2)
}

f1Score <- function(prediction, truth, trurecs) {
  sensitivity <- sensitivityNew(prediction, truth, trurecs)
  precision <- precisionValue(prediction, truth)
  return (2 * ((precision * sensitivity) / (precision + sensitivity))) 
}

precisionValue <- function(prediction, truth) {
  levels(prediction) <- c(0, 1) # Some eejit has predicted 1 for everything so it doesn't have the right factor levels
  evaluation.df <- data.frame(PREDICTION = prediction, TRUTH = truth)
  tpcount <- nrow(filter(evaluation.df, PREDICTION == 1 & PREDICTION == TRUTH))
  fpcount <- nrow(filter(evaluation.df, PREDICTION == 1 & PREDICTION != TRUTH))
  return (tpcount / (tpcount + fpcount))
}

sensitivityNew <- function(prediction, truth, trurecs) {
  levels(prediction) <- c(0, 1) # Some eejit has predicted 1 for everything so it doesn't have the right factor levels
  evaluation.df <- data.frame(PREDICTION = prediction, TRUTH = truth)
  tpcount <- nrow(filter(evaluation.df, PREDICTION == 1 & PREDICTION == TRUTH))
  return (tpcount / trurecs)
}

specificityNew <- function(prediction, truth) {
  levels(prediction) <- c(0, 1) # Some eejit has predicted 1 for everything so it doesn't have the right factor levels
  evaluation.df <- data.frame(PREDICTION = prediction, TRUTH = truth)
  subrecs <- nrow(filter(evaluation.df, PREDICTION == 1))
  fpcount <- nrow(filter(evaluation.df, PREDICTION == 1 & PREDICTION != TRUTH))
  return (1 - (fpcount / subrecs))
}

binaryPredictions <- function(rootDir) {
  df <- env$df
  binaryPredictorSummary <- BinaryPredictor(df, df[[dvname]], min.robustness = 0.1)
  return (print(binaryPredictorSummary, file = paste(rootDir, "/binarypredictors", "-", name, ".yaml", sep="")))
}

outputStats <- function(rootDir, name, num.variants, calculate.max.min = TRUE) {
  
  df <- env$df
  
  if (num.variants == 0) {
    num.variants = nrow(filter(df, IS_VARIANT == 1))
  }
  
  print(paste("NUM VARIANTS: ", num.variants, sep=""))
  
  max.accuracy <- 0
  best.threshold <- 0
  for (i in 1:100) {
    threshold.attempt <- i
    df.new <- mutate(df, IS_VARIANT_NEW = votes >= threshold.attempt)
    accuracy <- balancedAccuracyNew(df.new$IS_VARIANT_NEW, df.new$IS_VARIANT, num.variants)
    if (!is.na(accuracy) & accuracy > max.accuracy) {
      max.accuracy = accuracy
      best.threshold <- threshold.attempt
    }
  }
  print(paste("BEST ACCURACY FOR ", name, ": ", max.accuracy, sep=""))
  print(paste("BEST THRESHOLD FOR ", name, ": ", best.threshold, sep=""))
  
  max.accuracy <- 0
  best.threshold <- 0
  for (i in 1:100) {
    threshold.attempt <- i / 100
    df.new <- mutate(df, IS_VARIANT_NEW = vote.rate >= threshold.attempt)
    accuracy <- balancedAccuracyNew(df.new$IS_VARIANT_NEW, df.new$IS_VARIANT, num.variants)
    if (!is.na(accuracy) & accuracy > max.accuracy) {
      max.accuracy = accuracy
      best.threshold <- threshold.attempt
    }
  }
  print(paste("BEST ACCURACY FOR RATE FOR ", name, ": ", max.accuracy, sep=""))
  print(paste("BEST THRESHOLD FOR RATE FOR ", name, ": ", best.threshold, sep=""))
  
  max.accuracy <- 0
  best.threshold <- 0
  for (i in 1:100) {
    threshold.attempt <- i / 100
    df.new <- mutate(df, IS_VARIANT_NEW = vote.rate.from.ranked >= threshold.attempt)
    accuracy <- balancedAccuracyNew(df.new$IS_VARIANT_NEW, df.new$IS_VARIANT, num.variants)
    if (!is.na(accuracy) & accuracy > max.accuracy) {
      max.accuracy = accuracy
      best.threshold <- threshold.attempt
    }
  }
  print(paste("BEST ACCURACY FOR RATE FROM RANKED FOR ", name, ": ", max.accuracy, sep=""))
  print(paste("BEST THRESHOLD FOR RATE FROM RANKED FOR ", name, ": ", best.threshold, sep=""))
  
  #roc.votes <- roc(df$IS_VARIANT, df$votes)
  #print("ROC for votes:")
  #plot(roc.votes, print.thres = TRUE)
  #print(roc.votes)
  
  print(paste("NEW NAIVE CONSENSUS ACCURACY: ", balancedAccuracyNew(df$naive.consensus, df$IS_VARIANT, num.variants)))
  print(paste("NEW OPTIMAL CONSENSUS ACCURACY: ", balancedAccuracyNew(df$optimal.consensus, df$IS_VARIANT, num.variants)))
  
  print(paste("NEW NAIVE CONSENSUS RATE ACCURACY: ", balancedAccuracyNew(df$naive.consensus.from.rate, df$IS_VARIANT, num.variants)))
  print(paste("NEW OPTIMAL CONSENSUS RATE ACCURACY: ", balancedAccuracyNew(df$optimal.consensus.from.rate, df$IS_VARIANT, num.variants)))
  
  if (calculate.max.min) {
    
    pearsons.df <- calculateAll(df, cor)
    print(paste("MAX PEARSONS: ", max(pearsons.df$VALUE)))
    print(paste("MIN PEARSONS: ", min(pearsons.df$VALUE)))
    
    error.df <- calculateAll(df, function(x, y) {sum(x != y)})
    print(paste("MIN ERRORS: ", min(error.df$VALUE)))
    print(paste("MIN ERROR %: ", as.character((min(error.df$VALUE) / nrow(df)) * 100)))
    
    sensitivity.df <- calculateAll(df, function(x, y) {sensitivityNew(as.factor(x), as.factor(y), num.variants)})
    print(paste("NEW MAX SENSITIVITY: ", max(sensitivity.df$VALUE)))
    
    sensitivity.df <- calculateAll(df, function(x, y) {precisionValue(as.factor(x), as.factor(y))})
    print(paste("NEW MAX PRECISION: ", max(sensitivity.df$VALUE)))
    
    specificity.df <- calculateAll(df, function(x, y) {specificityNew(as.factor(x), as.factor(y))})
    print(paste("NEW MAX SPECIFICITY: ", max(specificity.df$VALUE)))
    
    accuracy.df <- calculateAll(df, function(x, y) {balancedAccuracyNew(x, y, num.variants)})
    print(paste("NEW MAX BALANCED ACCURACY: ", max(accuracy.df$VALUE)))

  }
    
}

createRankedBins <- function(df, rank.count, bin.size) {
  
  # Always round down here, so that we have complete bins. We may miss some data though.
  num.bins <- floor(rank.count / bin.size)
  bins.range <- 1:num.bins
  ranges <- as.data.frame(sapply(bins.range, function(a) {
    multiplier <- a-1
    start <- (multiplier * bin.size) + 1
    end <- (multiplier * bin.size) + bin.size
    return(start:end)
  }))
  i <- 1
  for (col in colnames(ranges)) {
    ranks <- select(df, num_range("RANK_", ranges[,col]))
    rank.sums <- rowSums(ranks)
    threshold <- bin.size / 2
    bin.consensus <- as.numeric(rank.sums >= threshold)
    bin.rate <- rank.sums / bin.size
    df[,paste("RANKBIN_", i, sep="")] <- bin.consensus
    df[,paste("RANKBINRATE_", i, sep="")] <- bin.rate
    i <- i+1
  }
  
  return(df)
  
}

linearize <- function() {
  df <- env$df
  df$ref.allele.percentage <- BinaryCut(df$ref.allele.percentage, df$IS_VARIANT)
  df$non.ref.allele.percentage <- BinaryCut(df$non.ref.allele.percentage, df$IS_VARIANT)
  #df$tumor.coverage.percentage <- BinaryCut(df$tumor.coverage.percentage, df$IS_VARIANT)
  #df$normal.coverage.percentage <- BinaryCut(df$normal.coverage.percentage, df$IS_VARIANT)
  env$df <- df
}

calculateNewFeatures <- function() {
  df <- env$df
  callers.to.use <- 5
  bin.size <- 2
  # Add the top N callers based on their accuracy against the consensus
  consensus.accuracy.df <- calculateAll(df, function(x, y) {balancedAccuracy(x, y)}, truthCol = "optimal.consensus.from.rate")
  consensus.accuracy.df <- arrange(consensus.accuracy.df, desc(VALUE))
  iterations <- min(callers.to.use, nrow(consensus.accuracy.df))
  callers.used <- as.character(vector(length=iterations))
  for(i in 1:iterations) {
    submission.id <- consensus.accuracy.df$COLUMN[i]
    column.name <- paste("RANK_", i, sep="")
    df[,column.name] <- as.numeric(df[,as.character(submission.id)])
    callers.used[i] <- as.character(submission.id)
  }
  names <- colnames(df)[grep("RANK.*", colnames(df))]
  submissions <- df[,names]
  df$votes.from.ranked = rowSums(submissions, na.rm = TRUE)
  df <- mutate(df, vote.rate.from.ranked = votes.from.ranked / length(names))
  df <- createRankedBins(df, length(names), bin.size)
  df <- mutate(df, ref.allele.percentage = (RefAllele / (RefAllele + NonRefAllele)) * 100)
  df <- mutate(df, non.ref.allele.percentage = (NonRefAllele / (RefAllele + NonRefAllele)) * 100)
  df <- mutate(df, tumor.coverage.percentage = (TumourCoverage / (TumourCoverage + NormalCoverage)) * 100)
  df <- mutate(df, normal.coverage.percentage = (NormalCoverage / (TumourCoverage + NormalCoverage)) * 100)
  env$df <- df
  env$callers.used <- callers.used
}

calculatePredictors <- function(df, dvname, base.dir) {
  binaryPredictorSummary <- BinaryPredictor(df, df[[dvname]], min.robustness = 0.1)
  df.summary <- print(binaryPredictorSummary, file = paste(base.dir, "/binarypredictors.yaml", sep=""))
}

filterColumns <- function(includeKey = FALSE, includeDv = TRUE) {
  df <- env$df
  names <- c()
  if (includeDv) {
    names <- c(names, dvname)  
  }
  names <- c(names, "vote.rate.from.ranked")
  #names <- c(names, "vote.rate")
  #names <- c(names, "votes", "vote.rate")
  #names <- c(names, "naive.consensus", "optimal.consensus", "naive.consensus.from.rate", "optimal.consensus.from.rate")
  if (includeKey) {
    names <- c(names, "CHROM", "POS")
  }
  names <- c(names, colnames(df)[grep("RANKBINRATE_.*", colnames(df))])
  names <- c(names, colnames(df)[grep("RANKBIN_.*", colnames(df))])
  names <- c(names, colnames(df)[grep("RANK_.*", colnames(df))])
  #names <- c(names, "TumourCoverage", "NormalCoverage")
  names <- c(names, "tumor.coverage.percentage")
  names <- c(names, "non.ref.allele.percentage")
  #names <- c(names, "RefAllele", "NonRefAllele")
  #names <- c(names, "ref.allele.percentage", "non.ref.allele.percentage", "tumor.coverage.percentage", "normal.coverage.percentage")
  names <- c(names, "BaseQual", "MapQual", "Trinucleotide", "GenomicLocation")
  #names <- c(names, "RefAllele", "NonRefAllele", "BaseQual", "TumourCoverage", "NormalCoverage",
             #"MapQual", "ReadPosition", "HomopolymerRate", "GCcontent", "Distance", "StrandBias",
             #"Trinucleotide", "GenomicLocation")
  #print(paste("Missing column: ", setdiff(names, colnames(df))))
  env$df <- df[,names]
}

splitData <- function(training.proportion) {
  df <- env$df
  df[is.na(df)] <- 0
  smp_size <- floor(training.proportion * nrow(df))
  set.seed(238476238)
  train_ind <- base::sample(seq_len(nrow(df)), size = smp_size)
  env$train <- df[train_ind, ]
  env$test <- df[-train_ind, ]
}

prepareData <- function(name, threshold, rate.threshold) {
  
  # and some consensus values based on these votes
  # This calculates additional columns on the data frame based on the number of votes from different callers
  calculateVotesAndConsensus(threshold, rate.threshold)
  
  # Join-in the genomic features and the ground truth to the main data frame 
  joinData()

  # Calculate new features based on the existing data
  calculateNewFeatures()
  
  # Convert non-linear numeric variables to factors
  #linearize()

  env$df$data.set.name <- name 
  
}

evalerror <- function(preds, dtrain) {
  truth <- getinfo(dtrain, "label")
  trurecs <- nrow(filter(as.data.frame(truth), truth == 1))
  preds.df <- mutate(as.data.frame(preds), pred_binary = preds >= 0.5)
  accuracy <- balancedAccuracyNew(preds.df$pred_binary, truth, trurecs)
  f1.score <- f1Score(preds.df$pred_binary, truth, trurecs)
  return(list(metric = "accuracy", value = accuracy))
}

to.training.matrix <- function(train) {
  train.dummies <- dummy.data.frame(train)
  return (as.matrix(train.dummies[, names(train.dummies) != dvname]))
}

trainXgb <- function(train, iterations) {
  yResponseVariable <- train[[dvname]]
  train.matrix <- to.training.matrix(train)
  xgboost(data       = train.matrix,
         label       = yResponseVariable,
         nrounds     = iterations,
         objective   = "binary:logistic",
         eval_metric = evalerror,
         nfold       = 5, 
         maximize    = TRUE)
}

trainXgbCv <- function(train, dummy) {
  yResponseVariable <- train[[dvname]]
  train.matrix <- to.training.matrix(train)
  xgb.cv(data       = train.matrix,
          label       = yResponseVariable,
          nrounds     = 1000000,
          objective   = "binary:logistic",
          eval_metric = evalerror,
          nfold       = 5, 
          maximize    = TRUE,
          early.stop.round = 100)
}

trainGlmnet <- function(train) {
  yResponseVariable <- train[[dvname]]
  train.matrix <- to.training.matrix(train)
  return(cv.glmnet(train.matrix, yResponseVariable, family="binomial"))
}

doPredictions <- function(test, model, threshold = 0.5) {
  test.new <- test
  test.new$CHROM <- NULL
  test.new$POS <- NULL
  test.matrix <- to.training.matrix(test.new)
  predicted.test <- predict(model, test.matrix)
  predicted.df = as.data.frame(predicted.test)
  colnames(predicted.df) <- "predicted.test"
  predicted.df <- rename(predicted.df, PREDICTED_IS_VARIANT = predicted.test)
  # The threshold here was chosen by experimentation
  predicted.df <- mutate(predicted.df, PREDICTED_IS_VARIANT_BINARY = predicted.df$PREDICTED_IS_VARIANT > threshold)
  predicted.df$PREDICTED_IS_VARIANT_BINARY = as.numeric(predicted.df$PREDICTED_IS_VARIANT_BINARY)
  predicted.df$CHROM = test$CHROM
  predicted.df$POS = test$POS
  return (predicted.df)
}

validateModel <- function(test, model, dataset.name, name, num.variants, result.index, threshold = 0.5) {

  results <- env.results$results
  
  results[result.index, "NAME"] <- name
  
  print(paste("THRESHOLD FOR ", dataset.name, ": ", as.character(threshold), sep =""))
  
  print(paste("COLUMN NAMES FOR ", dataset.name, ": ", sep =""))
  print(str(test))
  
  if (num.variants == 0) {
    num.variants = nrow(filter(test, IS_VARIANT == 1))
  }
  
  #print(paste("Correlation for ", dataset.name, sep = ""))
  #print(cor.test(test$IS_VARIANT, predicted.test))
  
  predicted.df = doPredictions(test, model, threshold) 
  
  predicted.df$IS_VARIANT = test$IS_VARIANT
  
  #print(paste("ROC for ", dataset.name, sep = ""))
  #roc.test <- roc(predicted.df$IS_VARIANT, predicted.df$PREDICTED_IS_VARIANT)
  #plot(roc.test, print.thres = TRUE)
  #print(roc.test)
  
  predicted.df <- mutate(predicted.df, NOT_EQUALS = IS_VARIANT != PREDICTED_IS_VARIANT_BINARY) 
  predicted.not.equals <- sum(as.numeric(predicted.df$NOT_EQUALS))
  print(paste("PREDICTION ERRORS for ", dataset.name, ": ", predicted.not.equals, sep = ""))
  predicted.not.equals.percent = (predicted.not.equals / nrow(predicted.df)) * 100
  print(paste("PREDICTION ERROR % for ", dataset.name, ": ", predicted.not.equals.percent))
  
  print(paste("PREDICTION SENSITIVITY for ", dataset.name, ": ", sensitivity(as.factor(predicted.df$PREDICTED_IS_VARIANT_BINARY), as.factor(predicted.df$IS_VARIANT)), sep=""))
  print(paste("PREDICTION SPECIFICITY for ", dataset.name, ":", specificity(as.factor(predicted.df$PREDICTED_IS_VARIANT_BINARY), as.factor(predicted.df$IS_VARIANT)), sep=""))
  print(paste("PREDICTION ACCURACY for ", dataset.name, ":", balancedAccuracy(predicted.df$PREDICTED_IS_VARIANT_BINARY, predicted.df$IS_VARIANT), sep = ""))
  
  trurecs <- num.variants
  
  sensitivity <- sensitivityNew(predicted.df$PREDICTED_IS_VARIANT_BINARY, predicted.df$IS_VARIANT, trurecs)
  precision <- precisionValue(predicted.df$PREDICTED_IS_VARIANT_BINARY, predicted.df$IS_VARIANT)
  specificity <- specificityNew(predicted.df$PREDICTED_IS_VARIANT_BINARY, predicted.df$IS_VARIANT)
  
  accuracy <- balancedAccuracyNew(predicted.df$PREDICTED_IS_VARIANT_BINARY, predicted.df$IS_VARIANT, trurecs)
  
  print(paste("NEW PREDICTION SENSITIVITY for ", dataset.name, ": ", sensitivity, sep=""))
  print(paste("NEW PREDICTION PRECISION for ", dataset.name, ": ", precision, sep=""))
  print(paste("NEW PREDICTION SPECIFICITY for ", dataset.name, ":", specificity, sep=""))
  print(paste("NEW PREDICTION ACCURACY for ", dataset.name, ":", accuracy, sep = ""))
  
  results[result.index, "ACCURACY"] <- accuracy
  
  max.accuracy <- 0
  best.threshold <- 0
  for (i in 1:100) {
    threshold.attempt <- i / 100
    trurecs <- num.variants
    predicted.df.new <- mutate(predicted.df, IS_VARIANT_NEW = PREDICTED_IS_VARIANT > threshold.attempt)
    accuracy <- f1Score(predicted.df.new$IS_VARIANT_NEW, predicted.df.new$IS_VARIANT, trurecs)
    if (!is.na(accuracy) & accuracy > max.accuracy) {
      max.accuracy = accuracy
      best.threshold <- threshold.attempt
    }
  }
  print(paste("BEST ACCURACY FOR ", dataset.name, ": ", max.accuracy, sep=""))
  print(paste("BEST THRESHOLD FOR ", dataset.name, ": ", best.threshold, sep=""))
  
  env.results$results <- results
  
  return (predicted.df)

}