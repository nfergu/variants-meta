# Variants-Meta

Builds a model for the [ICGC-TCGA DREAM Somatic Mutation Calling Meta-pipeline Challenge](https://www.synapse.org/#!Synapse:syn4588939/wiki/) 
(SMC-DNA Meta), and uses this model to make predictions against the tumor data provided as part of the challenge. The model is built with 
[gradient-boosted trees](http://xgboost.readthedocs.io/en/latest/model.html) using the [XGBoost library] in R.

To build a new model and score the tumor data against it:

 * In `metafunctions.R` change the `base.dir` variable at the top to point to the directory with the tumor data (synthetic and real) in it.
 * Run the `metafunctions.R` script, which creates some shared functions.
 * Run the `retrain-multi.R` script. This build a new model against 50% of the synthetic tumor data. The model is saved to the base directory
   in `/model/latestboostingall.rds`.
 * Run the `validate.R` script. This loads the model saved in the previous step, and checks its accuracy against all of the synthetic tunor data.
 * Run the `scoreall.R` script. This loads the model saved previously, and uses it to make predictions based on all of the tumor data ](synthetic and real). These predictions are saved to the `predictions` directory within the base directory.

Note: the code needs to be tidied-up. Apologies!

# TODO

 * Tidy-up code
 * Better documentation

# License

Variants-Meta is licensed under the [Apache License, Version 2.0](LICENSE.md)



