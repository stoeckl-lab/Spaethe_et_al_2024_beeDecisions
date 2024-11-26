%fitting Bayesian decision model and estimate weights for combined cues
close all
clear all

%%inputs

%animals prior for the two cues
prior_col1=[1.0000    0.9000    1.0000    0.9000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000]';%blue-orange
prior_col2=[0.8000    0.6000    0.6000    0.9000    0.9000    0.9000    0.6000    0.9000    0.8000    0.7000]';%blue-teal

prior_pat=[0.8000    0.5000    0.7000    0.8000    0.9000    0.7000    0.6000    0.8000    0.6000    0.8000]';%prior of pattern 
prior_shap=[0.8000    0.9000    0.8000    0.6000    0.8000    0.8000    0.7000    0.5000    0.6000    0.7000]';%prior of shape

%posterior - observed choices of animals in cue conflict
%pattern
posterior_col1_pat=[1.0000   1.0000    1.0000    1.0000    0.9000    0.9000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000 1.0000    1.0000    1.0000    0.9000    0.9000    0.9000    1.0000]';%blue orange pattern
posterior_col2_pat=[0.6000    0.6000    0.8000    0.6000    0.4000    0.4000    0.3000    0.9000    0.7000    0.8000    0.4000    0.6000    0.7000 0.8000    0.6000    0.4000    0.7000    0.7000    0.6000    0.8000]';%blue teal pattern
%shape
posterior_col1_shap=[1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    0.9000    1.0000    1.0000    1.0000    0.8000    1.0000    1.0000 0.9000    1.0000]';%blue orange pattern
posterior_col2_shap=[0.8000    0.7000    0.8000    0.4000    0.6000    0.9000    0.3000    0.4000    0.7000    0.6000    0.8000    0.6000    0.5000 0.7000    0.7000]';%blue teal pattern

%probability of each observation
observation = 0.5;

%% loop to calculate results for all 4 tested conditions
reps=10000;%repetition of fits with different subsets
all_params=nan(reps,4);
all_estPosteriors=nan(100,4);

for u=1:4
%select the right combination of inputs for all 4 conditions
if u==1; prior_col=prior_col1;prior_sec=prior_pat;posterior=posterior_col1_pat;end
if u==3; prior_col=prior_col2;prior_sec=prior_pat;posterior=posterior_col2_pat;end
if u==2; prior_col=prior_col1;prior_sec=prior_shap;posterior=posterior_col1_shap;end
if u==4; prior_col=prior_col2;prior_sec=prior_shap;posterior=posterior_col2_shap;end

prior_col(prior_col==1)=prior_col(prior_col==1)-0.00001;%since calculation does not work with 1 or 0 integer
posterior(posterior==1)=posterior(posterior==1)-0.00001;%since calculation does not work with 1 or 0 integer

%% set up Bayesian model fit for conflict situation
%this is for the option of choosing colour over the second feature
prior1=prior_col;
prior2=prior_sec;

lb=[0];ub=[1];
start_param=0.5;%initialise at equal weights
options = optimset('Algorithm','active-set','Display','off','MaxIter',10^6);

%generate prior and observation matrices, to model all possible
%combinations
prior1_m=repmat(prior1',size(prior2,1),1);
prior2_m=repmat(prior2,1,size(prior1,1));

est_params=nan(reps,1);

for i=1:reps
    subsample=randsample(numel(prior2_m),length(posterior));%initialise a different subsample for each run
est_params(i)=fmincon(@(params) BayesianFun(params,prior1_m,prior2_m,posterior,observation,subsample),start_param,[],[],[],[],lb,ub,[],options);
end

%use mean for posterior estimation
mean_param=nanmean(est_params);

estPosteriors_col=prior1_m.^mean_param.*(1-prior2_m).^(1-mean_param).*observation...
    ./(prior1_m.^mean_param.*(1-prior2_m).^(1-mean_param).*observation+(1-prior1_m).^mean_param.*(prior2_m).^(1-mean_param).*observation);

%collect data from all conditions
all_params(:,u)=est_params;
all_estPosteriors(:,u)=estPosteriors_col(:);

%% plot data

figure;hold on;
subplot(1,2,1);hold on;
data=[posterior;estPosteriors_col(:)];
groups=[ones(size(posterior));2*ones(size(estPosteriors_col(:)))];
b1=boxplot(data,[ones(size(posterior));2*ones(size(estPosteriors_col(:)))]);
for u=1:length(unique(groups))
    n=sum(groups==u);
plot(0.05*randn(n,1)+u*ones(n,1),data(groups==u,1),'.','MarkerSize',12,'color','k');hold on
end
ylabel('posterior probabilities')
ylim([0 1])
plot([0,2.5],[0.5,0.5],'k--')

subplot(1,2,2);hold on
mean_param=nanmean(est_params,1);
bar([1,2],[mean_param,1-mean_param])
% sem_est_params=nanstd(est_params,1,1)/sqrt(numel(est_params)-1);
%95% credible interval of the estimate
cred_interval=abs(mean_param-[quantile(est_params,.025) quantile(est_params,.975)]);
errorbar([1,2],[mean_param,1-mean_param],[cred_interval(1),cred_interval(1)],[cred_interval(2),cred_interval(2)],'.k')
ylim([0 1])
ylabel('weights (rel. units)')
plot([0,2.5],[0.5,0.5],'k--')

end

%% summary plot for all
figure('position',[100 100 900 500]);hold on;
subplot(1,3,[1 2]);hold on;
data=[all_estPosteriors(:)];
groups=[ones(size(all_estPosteriors,1),1);2*ones(size(all_estPosteriors,1),1);3*ones(size(all_estPosteriors,1),1);4*ones(size(all_estPosteriors,1),1)];
b1=boxplot(all_estPosteriors);
for u=1:length(unique(groups))
    n=sum(groups==u);
plot(0.05*randn(n,1)+u*ones(n,1),data(groups==u,1),'.','MarkerSize',12,'color','k');hold on
end
ylabel('posterior probabilities')
ylim([0 1])
plot([0,5],[0.5,0.5],'k--')
xlim([0 5])

subplot(1,3,3);hold on
mean_param=nanmean(all_params,1);
bar([1,2,3,4],[mean_param])
% sem_est_params=nanstd(est_params,1,1)/sqrt(numel(est_params)-1);
%95% credible interval of the estimate / use IQR
cred_interval_low=abs(mean_param-[quantile(all_params,.25,1)]);
cred_interval_up=abs(mean_param-[quantile(all_params,.75,1)]);
errorbar([1,2,3,4],mean_param,cred_interval_low,cred_interval_up,'.k')
ylim([0 1])
ylabel('weights (rel. units)')
plot([0,5],[0.5,0.5],'k--')
xlim([0 5])

%% in-script functions
function sse=BayesianFun(weights,prior1,prior2,posterior,observation,subsample)

estPosterior=prior1.^weights.*(1-prior2).^(1-weights).*observation...
    ./(prior1.^weights.*(1-prior2).^(1-weights).*observation+(1-prior1).^weights.*(prior2).^(1-weights).*observation);     

%if estPosterior is a matrix, select a subset for fitting:
subset_estPosterior=estPosterior(subsample);

%fit sse on probabilities (not ideal, because variance does not scale
%equally from 0 to 1)
sse=nansum((subset_estPosterior-posterior).^2);

%generate max likelihood estimate by calculating log odds ratio instead of
%probability, and taking sse from there
%odds ratio is not treating entries similarly, but massively over-punishing
%large entries (as not linear distance measured, but the ratio, which is
%extremely high for values very close to 1)
% logodds_estPosterior=log(subset_estPosterior./(1-subset_estPosterior));
% logodds_posterior=log(posterior./(1-posterior));
% sse=nansum((logodds_estPosterior-logodds_posterior).^2);

sse(isnan(sse))=10^6;
sse(isinf(sse))=10^6;
end

