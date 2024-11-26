# Spaethe_et_al_2024_beeDecisions
Analysis files for the publication "Bumblebees flexibly adjust learning and decision strategies to sensory information content in a multi-cue foraging task" , DOI follows soon. 

# Plotting of data and generation of additional files for statistical analysis
1) plotChoiceDistr.m plots the fraction of choices in the conflict, and second attribute tests, as well as of single attribute tests
2) plotNaiveChoices.m plots the fraction of choices in the naive tests, before training
3) plot_learningCurves.m plots the learning curves of the different conditions, and extracts choice data from individual training blocks (selectable) to be analysed separately

# Statistical analysis of choice data and training data
4) GLMM Mean Test Binomial learningBlocks.R performs a generalised linear model with random effects (glmer), with a binomial family and "logit" link on the choices for the rewarded and non-rewarded stimuli, included animal identity as random effects, and the stimulus type (i.e. distant and close colours) as a fixed effect
5) GLMM Mean Test Binomial patternVScol.R performs a generalised linear model with random effects (glmer), with a binomial family and "logit" link on the choices for the rewarded and non-rewarded stimuli, included animal identity as random effects, and the stimulus type (colour combination) as a fixed effect
6) GLMM Mean Test Naive patternVScol.R performs a generalised linear model with random effects (glmer), with a binomial family and "logit" link on the choices for the rewarded and non-rewarded stimuli, included animal identity as random effects, and the stimulus types (colour and pattern or shape) as a fixed effects

# Bayesian model of conflict choices based on single attribute learning as priors
7) bayesianModel.m reads in the choice distributions of single attribute tests, and uses these as priors to generate a Bayesian model, which predicts the fraction of choices for colour in the conflict tests. The priors are combined with a weighting factor, which is fitted to the observed choice data in the conflict tests.
