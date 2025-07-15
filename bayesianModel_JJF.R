# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2025 07 07
#     MODIFIED:	James Foster              DATE: 2025 07 14
#
#  DESCRIPTION: Modified from "bayesianModel.m" by Anna Stöckl.
#               Uses data from the single cue experiments to predict performance
#               in the conflict trials, fitting a weighting variable to describe
#               the relative weighting of each cue.
#               
#       INPUTS: (Data in script)
#               
#      OUTPUTS: Plots
#
#	   CHANGES: - 
#
#   REFERENCES: Berger, J.O. (1985).
#               Statistical Decision Theory and Bayesian Analysis (Springer New York)..
# 
#
#    EXAMPLES:  
#
# 
#TODO   ---------------------------------------------
#TODO   
#- LMER fit to individuals +
#- BRMS estimate priors and weights simultaneously +
#- Reorganise functions
#- Comment
#- BRMS joint model ?


# Load packages -----------------------------------------------------------
require(lme4)#for mixed effects modelling
require(glmmTMB)#for beta-binomial models
require(brms)#for (non-linear) Bayesian estimation

print(R.version.string)
citation("lme4")
citation("brms")
rstan::stan_version()
# Set up useful functions -------------------------------------------------


## Prior estimation functions --------------------------------------------


#set up function for calculating MLE binomial distribution
#N.B. when all proportions are calculated from the same number of observations, 
#this is equivalent calculating the mean probability.
BinomMLE = function(obs_p,#observed proportion
                    n,#total number of trials for each 
                    ...)
{
  #check that the original binomial data can be recovered
  if( any(as.logical(obs_p*n %% 1)) )
  {
    stop(
      "\n obs_p*n doesn't produce an integer
either obs_p or n are incorrect"
    ) 
  }
  mod = glm(cbind(obs_p*n, n - obs_p*n) ~ 1, #fit a binomial GLM
            family = binomial(link = 'logit'))
  lodds = coef(mod) #find the log odds
  prob = plogis(lodds) #find the probability
  ic = AIC(mod)
  prm = c(probability = prob, 
          AIC = ic)
  return(prm) #return the maximum likelihood probability
}

#set up function for calculating MLE beta-binomial distribution
#This is similar, but accounts for different variance structures around the mean.
BetaBinomMLE = function(obs_p,#observed proportion
                        n,#total number of trials for each 
                        ...)
{
  #check that the original binomial data can be recovered
  if( any(as.logical(obs_p*n %% 1)) )
  {
    stop(
      "\n obs_p*n doesn't produce an integer
either obs_p or n are incorrect"
    ) 
  }
  dt = data.frame(observed = obs_p*n,
                  n = n)
  mod = glmmTMB::glmmTMB(cbind(observed, n-observed) ~ 1, #fit a binomial GLM
                         dispformula = ~1, #fit a "dispersion" (1/variance) parameter
                         family = glmmTMB::betabinomial(link = 'logit'),
                         data = dt)
  prm = mod$fit$par #find the log odds and dispersion
  prm = c(prm,
          probability = plogis(prm[1]),
          AIC = AIC(mod))
  return(prm) #return the maximum likelihood probability
}

BinomLMER = function(obs_p,#observed proportion
                    n,#total number of trials for each 
                    ID = NULL,#individual identities (guessed by default)
                    ...)
{
  #check that the original binomial data can be recovered
  if( any(as.logical(obs_p*n %% 1)) )
  {
    stop(
      "\n obs_p*n doesn't produce an integer
either obs_p or n are incorrect"
    ) 
  }
  if(is.null(ID))
  {
    ID = 1:length(obs_p)#by default, all observations are different individuals
  }
  
  mod = glmer(cbind(obs_p*n, n - obs_p*n) ~ 1 + (1|ID), #fit a binomial GLMM
            family = binomial(link = 'logit'))
  lodds = fixef(mod) #find the log odds
  prob = plogis(lodds) #find the probability
  fe_sigma = coef(summary(mod))[2]#extract population variance
  re_sigma = VarCorr(mod)#extract individual variation
  ic = AIC(mod)
  prm = c(probability = prob, 
          sem = fe_sigma,
          ind_sd = re_sigma,
          AIC = ic)
  return(prm) #return the maximum likelihood probability
}

#set up function for calculating MLE beta-binomial distribution
#This is similar, but accounts for different variance structures around the mean.
#In this dataset, it is not possible to distinguish individual variance from 
#population dispersion
BetaBinomLMER = function(obs_p,#observed proportion
                        n,#total number of trials for each 
                        ID = NULL,#individual identities (guessed by default)
                        ...)
{
  #check that the original binomial data can be recovered
  if( any(as.logical(obs_p*n %% 1)) )
  {
    stop(
      "\n obs_p*n doesn't produce an integer
either obs_p or n are incorrect"
    ) 
  }
  
  if(is.null(ID))
  {
    ID = 1:length(obs_p)#by default, all observations are different individuals
  }
  
  dt = data.frame(observed = obs_p*n,
                  n = n,
                  ID = ID)
  mod = glmmTMB::glmmTMB(cbind(observed, n-observed) ~ 1 + (1|ID), #fit a binomial GLM
                         dispformula = ~1, #fit a "dispersion" (1/variance) parameter
                         family = glmmTMB::betabinomial(link = 'logit'),
                         data = dt)
  prm = mod$fit$par #find the log odds and dispersion
  fe_sigma = coef(summary(mod))$cond[2]#extract population variance
  re_sigma = VarCorr(mod)#extract individual variation
  prm = c(prm,
          probability = plogis(prm[1]),
          sem = fe_sigma,
          ind_sd = re_sigma,
          AIC = AIC(mod))
  return(prm) #return the maximum likelihood probability
}


## Set up the Bayesian choice function -------------------------------------
#this can be fitted via optimisation
EstPost = function(weights,
                   prior1,
                   prior2,
                   observation = 0.5)
{
  #construct the posterior probability equation using the 
  #priors, weights, and expected observation ratio
  numerator = ( (prior1^weights)*
                  (1-prior2)^(1-weights) )*
                  observation
  
  denominator =  ( (prior1^weights)*
                     (1-prior2)^(1-weights) )*
                  observation +
                ( (1-prior1)^weights )*
                ( (prior2)^(1-weights) )*
                observation     
  
  estPosterior = numerator/denominator
  return(estPosterior)
}

BayesianFun = function(weights,
                       prior1,
                       prior2,
                       posterior,
                       posterior_num = 10, #default 10 choices
                       observation = 0.5) #null hypothesis 50:50
{
  #estimate the posterior probability of choosing option "1"
  estPosterior = EstPost(weights,
                         prior1,
                         prior2,
                         observation)
  #calculate the log-likelihood of this estimate, given the observed posterior choices
  likelihood_estPosterior = sum( #add all log likelihoods (i.e. multiply probabilities)
    dbinom(x = posterior*posterior_num, #recover choice number
           size = posterior_num, #total choices
           prob = estPosterior, #true mean probability is our estimate
           log = TRUE) #log scaled for numerical stability
  )
  neg_ll = -likelihood_estPosterior #make negative for optimisation algorithm
  return(
    if(is.infinite(neg_ll) | is.na(neg_ll))
    {1e16}else
    {neg_ll}
  )
}

#Fitting via optimiser
Optim_Bfun = function(posterior,
                      prior1,
                      prior2,
                      wt = 0.5, #initial estimate of weighting (default 50:50)
                      method = "L-BFGS-B",
                      lower = 0,
                      upper = 1,
                      ... #passed to optim
)
{
  opt = optim(par = wt,
              fn = BayesianFun,
              prior1 = prior1,
              prior2 = prior2,
              posterior = posterior,
              method = method,
              lower = lower,
              upper = upper,
              ...
  )
  return( opt$par )
}

#Resample posterior for confidence intervals
Resample_Bfun = function(posterior,
                         wt = 0.5, #initial estimate of weighting (default 50:50)
                         prior1,
                         prior2,
                         nsamples = 1e4,
                         method = "L-BFGS-B",
                         lower = 0,
                         upper = 1,
                         ... #can't pass to optim this way?
)
{
  replicate(n = nsamples,
            expr = 
              {
                Optim_Bfun(posterior= sample(x = posterior,
                                             size = length(posterior),
                                             replace = TRUE),
                           prior1 = prior1,
                           prior2 = prior2,
                           wt = wt, #initial estimate of weighting (default 50:50)
                           method = method,
                           lower = lower,
                           upper = upper
                )
              }
  )
}


# Functions for BRMS ------------------------------------------------------

inv_logit = inv_logit_scaled


# fitting Bayesian decision model and estimate weights for combine --------


# Data --------------------------------------------------------------------

## Priors for each cue --------------------------------------------------
#animals' prior for the two cues
#blue-orange
prior_col1=c(1.0000,    0.9000,    1.0000,    0.9000,    1.0000,
             1.0000,    1.0000,    1.0000,    1.0000,    1.0000)
#blue-teal
prior_col2=c(0.8000,    0.6000,    0.6000,    0.9000,    0.9000,
             0.9000,    0.6000,    0.9000,    0.8000,    0.7000)
#prior of pattern 
prior_pat=c(0.8000,    0.5000,    0.7000,    0.8000,    0.9000,
            0.7000,    0.6000,    0.8000,    0.6000,    0.8000)
#prior of shape
prior_shap=c(0.8000,    0.9000,    0.8000,    0.6000,    0.8000,
            0.8000,    0.7000,    0.5000,    0.6000,    0.7000)

## posterior - observed choices of animals in cue conflict -----------------
#pattern
#blue orange pattern
posterior_col1_pat=c(1.0000,   1.0000,    1.0000,    1.0000,    0.9000,
                    0.9000,    1.0000,    1.0000,    1.0000,    1.0000,
                    1.0000,    1.0000,    1.0000,    1.0000,    1.0000,
                    1.0000,    0.9000,    0.9000,    0.9000,    1.0000)

#blue teal pattern
posterior_col2_pat=c(0.6000,    0.6000,    0.8000,    0.6000,    0.4000,
                     0.4000,    0.3000,    0.9000,    0.7000,    0.8000,
                     0.4000,    0.6000,    0.7000,    0.8000,    0.6000,
                     0.4000,    0.7000,    0.7000,    0.6000,    0.8000)

#shape
#blue orange pattern
posterior_col1_shap=c(1.0000,    1.0000,    1.0000,    1.0000,    1.0000,
                      1.0000,    0.9000,    1.0000,    1.0000,    1.0000,
                      0.8000,    1.0000,    1.0000,    0.9000,    1.0000)

#blue teal pattern
posterior_col2_shap=c(0.8000,    0.7000,    0.8000,    0.4000,    0.6000,
                      0.9000,    0.3000,    0.4000,    0.7000,    0.6000,
                      0.8000,    0.6000,    0.5000,    0.7000,    0.7000)

posterior_num = 10#data represents a set of n draws from a binomial distribution

## probability of each observation -----------------------------------------
observation = 0.5

### Estimate prior distribution ------------------------------------------
prior_num = 10#data represents a set of n draws from a binomial distribution

#blue-orange
p_col1 = BinomLMER(obs = prior_col1,
                  n = prior_num)
beta_col1 = BetaBinomMLE(obs = prior_col1,
                  n = prior_num)
#blue-teal
p_col2 = BinomLMER(obs = prior_col2,
                  n = prior_num)
beta_col2 = BetaBinomMLE(obs = prior_col2,
                  n = prior_num)
#pattern 
p_pat = BinomLMER(obs = prior_pat,
                  n = prior_num)
beta_pat = BetaBinomMLE(obs = prior_pat,
                  n = prior_num)
#shape
p_shap = BinomLMER(obs = prior_shap,
                  n = prior_num)
beta_shap = BetaBinomMLE(obs = prior_shap,
                  n = prior_num)
#Inspect and compare
cbind(mean = sapply(X = list(prior_col1, prior_col2, prior_pat, prior_shap),
             FUN = mean),
      binom = 
            rbind(p_col1, p_col2, p_pat, p_shap),
    beta.bino = 
rbind(beta_col1, beta_col2, beta_pat, beta_shap) )
#Beta binomial population model is a better fit than the binomial
#but the binomial with individual effects appears similar.


# Fit via optimiser -------------------------------------------------------

## Colour 1 vs pattern ---------------------------------------------------


#Find MLE weighting
wt_col1_pat = Optim_Bfun(
                    prior1 = p_col1[[1]],
                    prior2 = p_pat[[1]],
                    posterior = posterior_col1_pat
                    )
#Resample to get CI
wt_col1_pat_ci = quantile(x = 
                              Resample_Bfun(
                                            prior1 = p_col1[[1]],
                                            prior2 = p_pat[[1]],
                                            posterior = posterior_col1_pat
                                        ),
                         probs = c(0,1)+c(1,-1)*(0.05/2))
#Calculate estimates using MLE and CI for the weight
est_col1_pat = sapply(X = list(mean = wt_col1_pat,
                               lower = wt_col1_pat_ci[1],
                               upper = wt_col1_pat_ci[2]),
                      FUN = EstPost,
                      prior1 = prior_col1[1],
                      prior2 = prior_pat[1],
                      observation = 0.5
                      )


#Or simulate and plot the full 
wt_col1_pat_sim = Resample_Bfun(
    prior1 = p_col1[[1]],
    prior2 = p_pat[[1]],
    posterior = posterior_col1_pat
)
#in this case, posterior is always nearly 1.0
col1_pat_sim = EstPost(weights = wt_col1_pat_sim, 
                       prior1 = prior_col1[[1]],
                       prior2 = prior_pat[[1]],
                       observation = 0.5)

## Colour 2 vs pattern ---------------------------------------------------


#Find MLE weighting
wt_col2_pat = Optim_Bfun(
                    prior1 = p_col2[[1]],
                    prior2 = p_pat[[1]],
                    posterior = posterior_col2_pat
                    )
#Resample to get CI
wt_col2_pat_ci = quantile(x = 
                              Resample_Bfun(
                                            prior1 = p_col2[[1]],
                                            prior2 = p_pat[[1]],
                                            posterior = posterior_col2_pat
                                        ),
                         probs = c(0,1)+c(1,-1)*(0.05/2))
#Calculate estimates using MLE and CI for the weight
est_col2_pat = sapply(X = list(mean = wt_col2_pat,
                               lower = wt_col2_pat_ci[1],
                               upper = wt_col2_pat_ci[2]),
                      FUN = EstPost,
                      prior1 = prior_col2[[1]],
                      prior2 = prior_pat[[1]],
                      observation = 0.5
                      )


#Or simulate and plot the full 
wt_col2_pat_sim = Resample_Bfun(
    prior1 = p_col2[[1]],
    prior2 = p_pat[[1]],
    posterior = posterior_col2_pat
)


col2_pat_sim = EstPost(weights = wt_col2_pat_sim, 
                       prior1 = prior_col2[[1]],
                       prior2 = prior_pat[[1]],
                       observation = 0.5)


## Colour 1 vs shape ---------------------------------------------------


#Find MLE weighting
wt_col1_shap = Optim_Bfun(
                    prior1 = p_col1[[1]],
                    prior2 = p_shap[[1]],
                    posterior = posterior_col1_shap
                    )
#Resample to get CI
wt_col1_shap_ci = quantile(x = 
                              Resample_Bfun(
                                            prior1 = p_col1[[1]],
                                            prior2 = p_shap[[1]],
                                            posterior = posterior_col1_shap
                                        ),
                         probs = c(0,1)+c(1,-1)*(0.05/2))
#Calculate estimates using MLE and CI for the weight
est_col1_shap = sapply(X = list(mean = wt_col1_shap,
                               lower = wt_col1_shap_ci[1],
                               upper = wt_col1_shap_ci[2]),
                      FUN = EstPost,
                      prior1 = prior_col1[1],
                      prior2 = prior_shap[1],
                      observation = 0.5
                      )


#Or simulate and plot the full 
wt_col1_shap_sim = Resample_Bfun(
    prior1 = p_col1[[1]],
    prior2 = p_shap[[1]],
    posterior = posterior_col1_shap
)

col1_shap_sim = EstPost(weights = wt_col1_shap_sim, 
                       prior1 = prior_col1[[1]],
                       prior2 = prior_shap[[1]],
                       observation = 0.5)


## Colour 2 vs shape ---------------------------------------------------


#Find MLE weighting
wt_col2_shap = Optim_Bfun(
                    prior1 = p_col2[[1]],
                    prior2 = p_shap[[1]],
                    posterior = posterior_col2_shap
                    )
#Resample to get CI
wt_col2_shap_ci = quantile(x = 
                              Resample_Bfun(
                                            prior1 = p_col2[[1]],
                                            prior2 = p_shap[[1]],
                                            posterior = posterior_col2_shap
                                        ),
                         probs = c(0,1)+c(1,-1)*(0.05/2))
#Calculate estimates using MLE and CI for the weight
est_col2_shap = sapply(X = list(mean = wt_col2_shap,
                               lower = wt_col2_shap_ci[1],
                               upper = wt_col2_shap_ci[2]),
                      FUN = EstPost,
                      prior1 = prior_col2[[1]],
                      prior2 = prior_shap[[1]],
                      observation = 0.5
                      )


#Or simulate and plot the full 
wt_col2_shap_sim = Resample_Bfun(
    prior1 = p_col2[[1]],
    prior2 = p_shap[[1]],
    posterior = posterior_col2_shap
)

col2_shap_sim = EstPost(weights = wt_col2_shap_sim, 
                       prior1 = prior_col2[[1]],
                       prior2 = prior_shap[[1]],
                       observation = 0.5)


## Plot predictions ------------------------------------------------------

sim_data = data.frame( prop = c(col1_pat_sim,
                                col1_shap_sim,
                                col2_pat_sim,
                                col2_shap_sim),
                       condition = c(rep('col1_pat', length(col1_pat_sim)),
                                     rep('col1_shap', length(col1_shap_sim)),
                                     rep('col2_pat', length(col2_pat_sim)),
                                     rep('col2_shap', length(col2_shap_sim)))
)

post_data = data.frame( prop = c(posterior_col1_pat,
                                 posterior_col1_shap,
                                 posterior_col2_pat,
                                 posterior_col2_shap),
                        condition = c(rep('col1_pat', length(posterior_col1_pat)),
                                      rep('col1_shap', length(posterior_col1_shap)),
                                      rep('col2_pat', length(posterior_col2_pat)),
                                      rep('col2_shap', length(posterior_col2_shap)))
)

par(mar = c(3,3,0.5,0))
boxplot(prop ~ condition,
        ylim = c(0,1),
        data = post_data,
        pch = 19,
        col = c('darkgreen', 'darkgreen', 'green', 'green')
)

stripchart(prop ~ condition,
           vertical = TRUE,
           method = 'jitter',
          data = sim_data,
           jitter = 0.3,
           ylim = c(0,1),
           pch = 19,
           add = TRUE,
           col = adjustcolor('gray50', 0.1)
)
abline(h = c(0,1))
legend(x = 'bottomleft',
       legend = c('Observed',
                  'Simulated'),
       col = c('green', 'gray50'),
       pch = c(15, 19),
       cex = 0.7
)


# Fit via BRMS ------------------------------------------------------------


## Specify decision formula ----------------------------------------------


decision_formula = bf(formula = 
                        choice_count | trials(total_trials) ~ # response (binomial)
                         ( ( inv_logit(lPrior1)^inv_logit(lWt) ) * #p(H | C1)^w
                        ( 1-inv_logit(lPrior2))^(1-inv_logit(lWt) ) * #p(H | C2)^(1-w)
                                                   0.5)/ ( #p(O | H)
                         (inv_logit(lPrior1)^inv_logit(lWt)) * #p(H | C1)^w #1
                           (1-inv_logit(lPrior2))^(1-inv_logit(lWt)) * #p(H | C2)^(1-w) #1
                           0.5 + #p(H | O) #1
                         ((1-inv_logit(lPrior1))^inv_logit(lWt)) * #p(H | C1)^w #0
                           (inv_logit(lPrior2)^(1-inv_logit(lWt))) * #p(H | C2)^(1-w) #0
                           0.5 ), #p(H | O) #0
                      lWt ~ 1,  #single estimate per dataset
                      lPrior1 ~ 1,  #single estimate per dataset
                      lPrior2 ~ 1,  #single estimate per dataset
                      family = binomial(link = 'identity'),
                      nl = TRUE)


##Organise data for BRMS -------------------------------------------------
#organise data
dt_col1_pat = data.frame(choice_count = posterior_col1_pat*prior_num,
                         total_trials = prior_num)

dt_col2_pat = data.frame(choice_count = posterior_col2_pat*prior_num,
                         total_trials = prior_num)

dt_col1_shap = data.frame(choice_count = posterior_col1_shap*prior_num,
                         total_trials = prior_num)

dt_col2_shap = data.frame(choice_count = posterior_col2_shap*prior_num,
                         total_trials = prior_num)



## Set up priors based on test data --------------------------------------

                         
pr_col1_pat = get_prior(formula = decision_formula,
                         data = dt_col1_pat)
                         
pr_col2_pat = get_prior(formula = decision_formula,
                         data = dt_col2_pat)

pr_col1_shap = get_prior(formula = decision_formula,
                         data = dt_col1_shap)
                         
pr_col2_shap = get_prior(formula = decision_formula,
                         data = dt_col2_shap)
# prior class      coef group resp dpar   nlpar lb ub       source
# (flat)     b                           lPrior1            default
# (flat)     b Intercept                 lPrior1       (vectorized)
# (flat)     b                           lPrior2            default
# (flat)     b Intercept                 lPrior2       (vectorized)
# (flat)     b                               lWt            default
# (flat)     b Intercept                     lWt       (vectorized)


### Weighting priors -----------------------------------------------------


#add a broad prior for weightings
#general to all 
pr_col1_pat = within(pr_col1_pat, 
                      { prior[
                        class %in% 'b' & #just the fixed effects
                          nlpar %in% 'lWt'
                      ] = 'normal(0,4)' # most likely weightings in (0.02, 0.98), extremes (0.0003, 0.9997) 
                      }
)
pr_col2_pat = within(pr_col2_pat, 
                      { prior[
                        class %in% 'b' & #just the fixed effects
                          nlpar %in% 'lWt'
                      ] = 'normal(0,4)' # most likely weightings in (0.02, 0.98), extremes (0.0003, 0.9997)
                      }
)
pr_col1_shap = within(pr_col1_shap, 
                      { prior[
                        class %in% 'b' & #just the fixed effects
                          nlpar %in% 'lWt'
                      ] = 'normal(0,4)' # most likely weightings in (0.02, 0.98), extremes (0.0003, 0.9997)
                      }
)
pr_col2_shap = within(pr_col2_shap, 
                      { prior[
                        class %in% 'b' & #just the fixed effects
                          nlpar %in% 'lWt'
                      ] = 'normal(0,4)' # most likely weightings in (0.02, 0.98), extremes (0.0003, 0.9997)
                      }
)


### Cue 1 priors ---------------------------------------------------------


#add distribution from colour 1 test as a prior
pr_col1_pat = within(pr_col1_pat, 
                      { prior[
                        class %in% 'b' & #just the fixed effects
                          nlpar %in% 'lPrior1'
                      ] = 
                        with(p_col1,
                             {
                        paste0('normal(',
                                 qlogis(`probability.(Intercept)`), ',',
                               sum(sem, ind_sd.ID), ')') #the distribution of the prior
                      }
                      )
                      }
)
#add distribution from colour 2 test as a prior
pr_col2_pat = within(pr_col2_pat, 
                      { prior[
                        class %in% 'b' & #just the fixed effects
                          nlpar %in% 'lPrior1'
                      ] = 
                        with(p_col2,
                             {
                        paste0('normal(',
                                 qlogis(`probability.(Intercept)`), ',',
                               sum(sem, ind_sd.ID), ')') #the distribution of the prior
                      }
                      )
                      }
)
#add distribution from colour 1 test as a prior
pr_col1_shap = within(pr_col1_shap, 
                      { prior[
                        class %in% 'b' & #just the fixed effects
                          nlpar %in% 'lPrior1'
                      ] = 
                        with(p_col1,
                             {
                        paste0('normal(',
                                 qlogis(`probability.(Intercept)`), ',',
                               sum(sem, ind_sd.ID), ')') #the distribution of the prior
                      }
                      )
                      }
)
#add distribution from colour 2 test as a prior
pr_col2_shap = within(pr_col2_shap, 
                      { prior[
                        class %in% 'b' & #just the fixed effects
                          nlpar %in% 'lPrior1'
                      ] = 
                        with(p_col2,
                             {
                        paste0('normal(',
                                 qlogis(`probability.(Intercept)`), ',',
                               sum(sem, ind_sd.ID), ')') #the distribution of the prior
                      }
                      )
                      }
)


### Cue 2 priors ---------------------------------------------------------


#add distribution from pattern test as a prior
pr_col1_pat = within(pr_col1_pat, 
                      { prior[
                        class %in% 'b' & #just the fixed effects
                          nlpar %in% 'lPrior2'
                      ] = 
                        with(p_pat,
                             {
                        paste0('normal(',
                                 qlogis(`probability.(Intercept)`), ',',
                               max(sem, ind_sd.ID), ')') #the distribution of the prior
                      }
                      )
                      }
)
#add distribution from pattern test as a prior
pr_col2_pat = within(pr_col2_pat, 
                      { prior[
                        class %in% 'b' & #just the fixed effects
                          nlpar %in% 'lPrior2'
                      ] = 
                        with(p_pat,
                             {
                        paste0('normal(',
                                 qlogis(`probability.(Intercept)`), ',',
                               max(sem, ind_sd.ID), ')') #the distribution of the prior
                      }
                      )
                      }
)

#add distribution from shape test as a prior
pr_col1_shap = within(pr_col1_shap, 
                      { prior[
                        class %in% 'b' & #just the fixed effects
                          nlpar %in% 'lPrior2'
                      ] = 
                        with(p_shap,
                             {
                        paste0('normal(',
                                 qlogis(`probability.(Intercept)`), ',',
                               max(sem, ind_sd.ID), ')') #the distribution of the prior
                      }
                      )
                      }
)
#add distribution from pattern test as a prior
pr_col2_shap = within(pr_col2_shap, 
                      { prior[
                        class %in% 'b' & #just the fixed effects
                          nlpar %in% 'lPrior2'
                      ] = 
                        with(p_shap,
                             {
                        paste0('normal(',
                                 qlogis(`probability.(Intercept)`), ',',
                               max(sem, ind_sd.ID), ')') #the distribution of the prior
                      }
                      )
                      }
)


## Fit all models --------------------------------------------------------



system.time( # takes less than two minutes
  {
    fit_col1_pat = brm( formula = decision_formula,
                          data = dt_col1_pat,
                          family = binomial(link = "identity"), 
                          prior = pr_col1_pat,
                          control = list(adapt_delta = 0.9),
                          silent = 2, #don't print lots of iteration information
                          iter = 2000,
                          chains= 4,
                          cores = 4,
                          backend = 'cmdstanr')
    
    fit_col2_pat = brm( formula = decision_formula,
                          data = dt_col2_pat,
                          family = binomial(link = "identity"), 
                          prior = pr_col2_pat,
                          control = list(adapt_delta = 0.9),
                          silent = 2, #don't print lots of iteration information
                          iter = 2000,
                          chains= 4,
                          cores = 4,
                          backend = 'cmdstanr')
    fit_col1_shap = brm( formula = decision_formula,
                          data = dt_col1_shap,
                          family = binomial(link = "identity"), 
                          prior = pr_col1_shap,
                          control = list(adapt_delta = 0.9),
                          silent = 2, #don't print lots of iteration information
                          iter = 2000,
                          chains= 4,
                          cores = 4,
                          backend = 'cmdstanr')
    
    fit_col2_shap = brm( formula = decision_formula,
                          data = dt_col2_shap,
                          family = binomial(link = "identity"), 
                          prior = pr_col2_shap,
                          control = list(adapt_delta = 0.9),
                          silent = 2, #don't print lots of iteration information
                          iter = 2000,
                          chains= 4,
                          cores = 4,
                          backend = 'cmdstanr')
  }
)


## Extract model predictions ---------------------------------------------
mod_lst = list(
  col1_pat = fit_col1_pat,
  col1_shap = fit_col1_shap,
  col2_pat = fit_col2_pat,
  col2_shap = fit_col2_shap
)

sm_lst = lapply(mod_lst,
                summary,
                robust = TRUE)

QuantEffects = function(mod)
{
  sm = summary(object = mod,
               robust = TRUE)
  qe = data.frame(apply(X = t(sm$fixed[,c(1,3:4)]),
                        MARGIN = 2,
                        FUN = plogis)
  )
  return( qe )
}

qe_lst = lapply(mod_lst,
                QuantEffects)

draw_lst = lapply(mod_lst,
                  prepare_predictions, 
                  ndraws = NULL)#use all


Draws2Post = function(draw)
{
  with(draw$nlpars,
       {
         EstPost(weights = plogis(c(lWt$fe$b)), 
                 prior1 = plogis(c(lPrior1$fe$b)),
                 prior2 = plogis(c(lPrior2$fe$b)),
                 observation = 0.5)
       }
  )
}

ppo_lst = sapply(draw_lst,
                 Draws2Post
                  )
#in this case, posterior is always nearly 1.0


## Plot predictions ------------------------------------------------------

post_data = data.frame( prop = c(posterior_col1_pat,
                                  posterior_col1_shap,
                                  posterior_col2_pat,
                                  posterior_col2_shap),
                       condition = c(rep('col1_pat', length(posterior_col1_pat)),
                                     rep('col1_shap', length(posterior_col1_shap)),
                                     rep('col2_pat', length(posterior_col2_pat)),
                                     rep('col2_shap', length(posterior_col2_shap)))
                       )

par(mar = c(3,3,0.5,0))
boxplot(prop ~ condition,
        ylim = c(0,1),
        data = post_data,
        pch = 19,
        col = c('darkgreen', 'darkgreen', 'green', 'green')
        )

stripchart(x = data.frame(ppo_lst),
           vertical = TRUE,
           method = 'jitter',
           jitter = 0.3,
           ylim = c(0,1),
           pch = 19,
           add = TRUE,
           col = adjustcolor('gray50', 0.1)
)
abline(h = c(0,1))
legend(x = 'bottomleft',
       legend = c('Observed',
                  'Predicted'),
       col = c('green', 'gray50'),
       pch = c(15, 19),
       cex = 0.7
       )

CalcDensWt = function(draw)
{
  dd = density(plogis(draw$nlpars$lWt$fe$b),
          from = 0, to = 1)
}

dens_lst = lapply(draw_lst, 
                  CalcDensWt)

PlotDens = function(dens,
                    col_fill = 'darkgreen',
                    scaling = 'original',
                    normby = 4000,
                    ... # passed to plot
                    )
{
  plot(x = NULL,
       ylab = switch(EXPR = scaling,
                     original = 'estimate density (rel. units)',
                     auc_1 = 'probability density (rel. units)',
                     bymax = 'normalised density (rel. units)',
                     bynorm = 'normalised density (rel. units)',
                     'estimate density (rel. units)'),
       xlab ='weight (rel. units)',
       ...
  )
  xx = with(dens, {c(x, rev(x))})
  yy = with(dens, 
            {
              switch(EXPR = scaling,
              original = c(y, rep(0, length(y))),
              auc_1 = c(y/ #smoothing makes this less than 1.0, integrate to estimate smoothed area
                          sfsmisc::integrate.xy(x = x, fx = y), 
                        rep(0, length(y))),
              bymax = c(y/max(y, na.rm = TRUE), 
                        rep(0, length(y))),
              bynorm = c(y/normby, 
                        rep(0, length(y))),
              c(y, rep(0, length(y)))
              )
            }
            )
  polygon(x = xx,
          y = yy,
          col = col_fill,
          border = NA
          )
  abline(h = 0)
}

#original version
par(mfrow = c(2,2),
    mar = c(5,5,3,0))
mapply(FUN = PlotDens,
       dens = dens_lst,
       col_fill = c('darkgreen', 
                    'darkgreen', 
                    'green', 
                    'green'),
       main = c('Distant colours vs patterns',
                'Distant colours vs shapes',
                'Close colours vs patterns',
                'Close colours vs shapes'),
       xlim = list(c(0,1)),
       ylim = list(c(0,12))
       )

#normalised to AUC version
par(mfrow = c(2,2),
    mar = c(5,5,3,0))
mapply(FUN = PlotDens,
       dens = dens_lst,
       scaling = 'auc_1',
       col_fill = c('darkgreen', 
                    'darkgreen', 
                    'green', 
                    'green'),
       main = c('Distant colours vs patterns',
                'Distant colours vs shapes',
                'Close colours vs patterns',
                'Close colours vs shapes'),
       xlim = list(c(0,1)),
       ylim = list(c(0,17))
       )
#normalised to maximum density version
par(mfrow = c(2,2),
    mar = c(5,5,3,0))
mapply(FUN = PlotDens,
       dens = dens_lst,
       scaling = 'bynorm',
       normby = max(unlist(do.call(rbind, dens_lst)[,'y'])),
       col_fill = c('darkgreen', 
                    'darkgreen', 
                    'green', 
                    'green'),
       main = c('Distant colours vs patterns',
                'Distant colours vs shapes',
                'Close colours vs patterns',
                'Close colours vs shapes'),
       xlim = list(c(0,1)),
       ylim = list(c(0,1))
       )

# Plot as bars ------------------------------------------------------------
#generate quantiles
quant_ppo = apply(X = data.frame(ppo_lst),
                  MARGIN = 2,
                  FUN = quantile,
                  probs = c(0,0.25,0.5,0.75,1)+c(1,0,0,0,-1)*0.05/2 )

#generates the same spacing as R's default barplot function
BarSpacer = function(n, 
                     spa = 0.2,
                     wdt = 1.0)
{
  seq(from = spa+1-wdt/2,
      to = n*(1+spa)-wdt/2,
      length.out = n )
}

par(mfrow = c(1,1))
barplot(height = quant_ppo[3,],
        col = c('darkgreen', 
                     'darkgreen', 
                     'green', 
                     'green'),
        ylim = c(0,1),
        ylab = 'weights (rel. units)'
        )
abline(h = 0.5)
arrows(x0 = BarSpacer(4),
       y0 = quant_ppo[1,],
       y1 = quant_ppo[5,],
       code = 3,
       angle = 90,
       length = 0.2,
       lwd = 3,
       lend ='butt',
       col = 'gray'
       )
arrows(x0 = BarSpacer(4),
       y0 = quant_ppo[2,],
       y1 = quant_ppo[4,],
       code = 3,
       angle = 90,
       length = 0.2,
       lwd = 5,
       lend = 'butt'
       )
