name <- "MODEL"

names <- as.vector(c("IS1", "IS2", "IS3", "IS4")) 

env <- new.env()

data.frames <- sapply(names, function(name) {
  
  print(paste("Preparing data for", name))
  
  clean.data <- loadCleanData(base.dir, name)
  
  df.clean <- clean.data$df
  genomic.features.clean <- clean.data$genomic.features
  vcf.clean <- clean.data$vcf
  
  env <<- new.env()
  env$df <- df.clean
  env$genomic.features <- genomic.features.clean
  env$vcf <- vcf.clean
  env$vcf.set = TRUE
  
  prepareData(name, votes.threshold, votes.rate.threshold)
  
  return (env$df)

})

all.clean <- as.data.frame(rbindlist(data.frames, use.names=TRUE, fill=TRUE))

# Get rid of NAs
all.clean[is.na(all.clean)] <- 0

all.clean.retrain <- all.clean

# Restart a clean training run here

env$df <- all.clean.retrain

outputStats(base.dir, "ALL", num.variants, calculate.max.min = FALSE)

training.proportion <- 0.5

# Filter out columns that will not be used for prediction (like the unranked submissions)
filterColumns()

write.table(env$df, file = paste(base.dir, "/all_synthetic.csv", sep=""), sep=",", quote = TRUE, row.names = FALSE, col.names = TRUE)

#binary.preds <- binaryPredictions(base.dir)

splitData(training.proportion)
train = env$train
test = env$test

print(str(env$train))
model <- trainXgb(train, 45)
#model <- trainGlmnet(train)

env.results <- new.env()
env.results$results <- data.frame(matrix(ncol = 0, nrow = 2))
env.results$results$NAME <- NA
env.results$results$ACCURACY <- NA

predicted.df.train <- validateModel(train, model, "TRAINING DATA", name, num.variants, 1, threshold = 0.5)

if (nrow(test) > 0) {
  predicted.df.test <- validateModel(test, model, "TEST DATA", name, num.variants, 2, threshold = 0.5)
}

print(env.results$results)

saveRDS(model, file = paste(base.dir, "/model/latestboostingall.rds", sep = ""))