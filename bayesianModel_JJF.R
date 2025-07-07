
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


### Estimate prior distribution ------------------------------------------
prior_num = 10#data represents a set of n draws from a binomial distribution

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

#blue-orange
p_col1 = BinomMLE(obs = prior_col1,
                  n = prior_num)
beta_col1 = BetaBinomMLE(obs = prior_col1,
                  n = prior_num)
#blue-teal
p_col2 = BinomMLE(obs = prior_col2,
                  n = prior_num)
beta_col2 = BetaBinomMLE(obs = prior_col2,
                  n = prior_num)
#pattern 
p_pat = BinomMLE(obs = prior_pat,
                  n = prior_num)
beta_pat = BetaBinomMLE(obs = prior_pat,
                  n = prior_num)
#shape
p_shap = BinomMLE(obs = prior_shap,
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
#N.B. Beta-binomial appears a better fit (lower AIC)
#betadisp >2 suggests underdispersion (low variance)
        # plot(x = 0:10,
        #      y = brms::dbeta_binomial(x = 0:10, size = 10, mu = 0.5, phi = 2.0),
        #      type = 'b',
        #      ylim = c(0, 1)*0.2,
        #      ylab = 'probability mass',
        #      xlab = 'choices')
        # lines(x = 0:10,
        #      y = brms::dbeta_binomial(x = 0:10, size = 10, mu = 0.5, phi = 1.0),
        #      lty = 3,
        #      pch = 1
        #      )
        # lines(x = 0:10,
        #      y = brms::dbeta_binomial(x = 0:10, size = 10, mu = 0.5, phi = 3.0),
        #      pch = 1,
        #      lty = 2) 
        # legend(x = 'top',
        #        legend = paste('phi =', c('1.0', '2.0', '3.0')),
        #        lty = c(3,1,2)
        #        )
#for binomial probability we would use:
        # dbinom(x = 10,
        #        size = 10,
        #        prob = p_col1[1])
#for beta binomial probability we would use
        # brms::dbeta_binomial(x = 10,
        #                      size =  10, 
        #                      mu = beta_col1[3],
        #                      phi = beta_col1[2])

# posterior - observed choices of animals in cue conflict -----------------
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

# probability of each observation -----------------------------------------
observation = 0.5


# Set up the Bayesian choice function -------------------------------------
#this can be fitted via optimisation
EstPost = function(weights,
                   prior1,
                   prior2,
                   observation = 0.5)
{
    #construct the posterior probability equation using the 
    #priors, weights, and expected observation ratio
    numerator = (prior1^weights)*(1-prior2)^(1-weights)*observation
    denominator =  (prior1^weights)*(1-prior2)^(1-weights)*observation+
        ((1-prior1)^weights)*((prior2)^(1-weights))*observation     
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


## Colour 1 vs pattern ---------------------------------------------------


#Find MLE weighting
wt_col1_pat = Optim_Bfun(
                    prior1 = p_col1[1],
                    prior2 = p_pat[1],
                    posterior = posterior_col1_pat
                    )
#Resample to get CI
wt_col1_pat_ci = quantile(x = 
                              Resample_Bfun(
                                            prior1 = p_col1[1],
                                            prior2 = p_pat[1],
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
    prior1 = p_col1[1],
    prior2 = p_pat[1],
    posterior = posterior_col1_pat
)
#in this case, posterior is always nearly 1.0
stripchart(x = EstPost(weights = wt_col1_pat_sim, 
                       prior1 = prior_col1[1],
                       prior2 = prior_pat[1],
                       observation = 0.5),
            vertical = TRUE,
            method = 'jitter',
            jitter = 0.2,
            ylim = c(0,1),
            pch = 19,
            col = gray(0, 0.005)
)

## Colour 2 vs pattern ---------------------------------------------------


#Find MLE weighting
wt_col2_pat = Optim_Bfun(
                    prior1 = p_col2[1],
                    prior2 = p_pat[1],
                    posterior = posterior_col2_pat
                    )
#Resample to get CI
wt_col2_pat_ci = quantile(x = 
                              Resample_Bfun(
                                            prior1 = p_col2[1],
                                            prior2 = p_pat[1],
                                            posterior = posterior_col2_pat
                                        ),
                         probs = c(0,1)+c(1,-1)*(0.05/2))
#Calculate estimates using MLE and CI for the weight
est_col2_pat = sapply(X = list(mean = wt_col2_pat,
                               lower = wt_col2_pat_ci[1],
                               upper = wt_col2_pat_ci[2]),
                      FUN = EstPost,
                      prior1 = prior_col2[1],
                      prior2 = prior_pat[1],
                      observation = 0.5
                      )


#Or simulate and plot the full 
wt_col2_pat_sim = Resample_Bfun(
    prior1 = p_col1[1],
    prior2 = p_pat[1],
    posterior = posterior_col1_pat
)

stripchart(x = EstPost(weights = wt_col2_pat_sim, 
                       prior1 = prior_col2[1],
                       prior2 = prior_pat[1],
                       observation = 0.5),
            vertical = TRUE,
            method = 'jitter',
            jitter = 0.2,
            ylim = c(0,1),
            pch = 19,
           # cex = 0.002,
            col = gray(0, 0.002)
)
# 
# 
# %% loop to calculate results for all 4 tested conditions
# reps=10000;%repetition of fits with different subsets
# all_params=nan(reps,4);
# all_estPosteriors=nan(100,4);
# 
# for u=1:4
# %select the right combination of inputs for all 4 conditions
# if u==1; prior_col=prior_col1;prior_sec=prior_pat;posterior=posterior_col1_pat;end
# if u==3; prior_col=prior_col2;prior_sec=prior_pat;posterior=posterior_col2_pat;end
# if u==2; prior_col=prior_col1;prior_sec=prior_shap;posterior=posterior_col1_shap;end
# if u==4; prior_col=prior_col2;prior_sec=prior_shap;posterior=posterior_col2_shap;end
# 
# prior_col(prior_col==1)=prior_col(prior_col==1)-0.00001;%since calculation does not work with 1 or 0 integer
# posterior(posterior==1)=posterior(posterior==1)-0.00001;%since calculation does not work with 1 or 0 integer
# 
# %% set up Bayesian model fit for conflict situation
# %this is for the option of choosing colour over the second feature
# prior1=prior_col;
# prior2=prior_sec;
# 
# lb=[0];ub=[1];
# start_param=0.5;%initialise at equal weights
# options = optimset('Algorithm','active-set','Display','off','MaxIter',10^6);
# 
# %generate prior and observation matrices, to model all possible
# %combinations
# prior1_m=repmat(prior1',size(prior2,1),1);
# prior2_m=repmat(prior2,1,size(prior1,1));
# 
# est_params=nan(reps,1);
# 
# for i=1:reps
#     subsample=randsample(numel(prior2_m),length(posterior));%initialise a different subsample for each run
# est_params(i)=fmincon(@(params) BayesianFun(params,prior1_m,prior2_m,posterior,observation,subsample),start_param,[],[],[],[],lb,ub,[],options);
# end
# 
# %use mean for posterior estimation
# mean_param=nanmean(est_params);
# 
# estPosteriors_col=prior1_m.^mean_param.*(1-prior2_m).^(1-mean_param).*observation...
#     ./(prior1_m.^mean_param.*(1-prior2_m).^(1-mean_param).*observation+(1-prior1_m).^mean_param.*(prior2_m).^(1-mean_param).*observation);
# 
# %collect data from all conditions
# all_params(:,u)=est_params;
# all_estPosteriors(:,u)=estPosteriors_col(:);
# 
# %% plot data
# 
# figure;hold on;
# subplot(1,2,1);hold on;
# data=[posterior;estPosteriors_col(:)];
# groups=[ones(size(posterior));2*ones(size(estPosteriors_col(:)))];
# b1=boxplot(data,[ones(size(posterior));2*ones(size(estPosteriors_col(:)))]);
# for u=1:length(unique(groups))
#     n=sum(groups==u);
# plot(0.05*randn(n,1)+u*ones(n,1),data(groups==u,1),'.','MarkerSize',12,'color','k');hold on
# end
# ylabel('posterior probabilities')
# ylim([0 1])
# plot([0,2.5],[0.5,0.5],'k--')
# 
# subplot(1,2,2);hold on
# mean_param=nanmean(est_params,1);
# bar([1,2],[mean_param,1-mean_param])
# % sem_est_params=nanstd(est_params,1,1)/sqrt(numel(est_params)-1);
# %95% credible interval of the estimate
# cred_interval=abs(mean_param-[quantile(est_params,.025) quantile(est_params,.975)]);
# errorbar([1,2],[mean_param,1-mean_param],[cred_interval(1),cred_interval(1)],[cred_interval(2),cred_interval(2)],'.k')
# ylim([0 1])
# ylabel('weights (rel. units)')
# plot([0,2.5],[0.5,0.5],'k--')
# 
# end
# 
# %% summary plot for all
# figure('position',[100 100 900 500]);hold on;
# subplot(1,3,[1 2]);hold on;
# data=[all_estPosteriors(:)];
# groups=[ones(size(all_estPosteriors,1),1);2*ones(size(all_estPosteriors,1),1);3*ones(size(all_estPosteriors,1),1);4*ones(size(all_estPosteriors,1),1)];
# b1=boxplot(all_estPosteriors);
# for u=1:length(unique(groups))
#     n=sum(groups==u);
# plot(0.05*randn(n,1)+u*ones(n,1),data(groups==u,1),'.','MarkerSize',12,'color','k');hold on
# end
# ylabel('posterior probabilities')
# ylim([0 1])
# plot([0,5],[0.5,0.5],'k--')
# xlim([0 5])
# 
# subplot(1,3,3);hold on
# mean_param=nanmean(all_params,1);
# bar([1,2,3,4],[mean_param])
# % sem_est_params=nanstd(est_params,1,1)/sqrt(numel(est_params)-1);
# %95% credible interval of the estimate / use IQR
# cred_interval_low=abs(mean_param-[quantile(all_params,.25,1)]);
# cred_interval_up=abs(mean_param-[quantile(all_params,.75,1)]);
# errorbar([1,2,3,4],mean_param,cred_interval_low,cred_interval_up,'.k')
# ylim([0 1])
# ylabel('weights (rel. units)')
# plot([0,5],[0.5,0.5],'k--')
# xlim([0 5])
# 
# %% in-script functions
# function sse=BayesianFun(weights,prior1,prior2,posterior,observation,subsample)
# 
# estPosterior=prior1.^weights.*(1-prior2).^(1-weights).*observation...
#     ./(prior1.^weights.*(1-prior2).^(1-weights).*observation+(1-prior1).^weights.*(prior2).^(1-weights).*observation);     
# 
# %if estPosterior is a matrix, select a subset for fitting:
# subset_estPosterior=estPosterior(subsample);
# 
# %fit sse on probabilities (not ideal, because variance does not scale
# %equally from 0 to 1)
# sse=nansum((subset_estPosterior-posterior).^2);
# 
# %generate max likelihood estimate by calculating log odds ratio instead of
# %probability, and taking sse from there
# %odds ratio is not treating entries similarly, but massively over-punishing
# %large entries (as not linear distance measured, but the ratio, which is
# %extremely high for values very close to 1)
# % logodds_estPosterior=log(subset_estPosterior./(1-subset_estPosterior));
# % logodds_posterior=log(posterior./(1-posterior));
# % sse=nansum((logodds_estPosterior-logodds_posterior).^2);
# 
# sse(isnan(sse))=10^6;
# sse(isinf(sse))=10^6;
# end

