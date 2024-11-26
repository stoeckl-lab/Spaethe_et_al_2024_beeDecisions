rm(list = ls())
graphics.off()
# Details ---------------------------------------------------------------
#       AUTHOR:	James Foster              DATE: 2020 03 01
#     MODIFIED:	James Foster              DATE: 2021 03 01
#     MODIFIED: Anna Stoeckl               DATE: 2023 10 16
#
#  DESCRIPTION: Loads a text files, fits a mixed effects logistic model and
#               estimates the p-value for mean ≤ reference mean.
#               
#       INPUTS: A ".xlsx" table with columns for experiment type ("Test"),
#               individual ID ("Moth number"), proportion of correct choices ("%")
#               total number of choices ("total Choice number"), correct choices
#               ("choices cross"), incorrect choices ("choices circle").
#               User should specify h0 reference correct choice rate (line 50).
#               
#      OUTPUTS: Plots of confidence intervals (.pdf).
#
#	   CHANGES: 
#             
#             
#
#   REFERENCES: Bates D, Maechler M, Bolker B, Walker S (2015).
#               Fitting Linear Mixed-Effects Models Using lme4.
#               Journal of Statistical Software, 67(1), 1-48.
#               doi:10.18637/jss.v067.i01.
#
#    EXAMPLES:  Fill out user input (line 50), then press ctrl+shift+s to run
#
# 
#TODO   ---------------------------------------------

# Useful functions --------------------------------------------------------
# . Load package ----------------------------------------------------------
#needs installing before first use (in Rstudio, see automatic message)
require(lme4)#package for fitting GLMMs (generalised linear mixed-effects models)
require(readxl)
require(lmerTest)
require(emmeans)
require(ggplot2)
require(beeswarm)
require(multcomp)

 
# Input Variables ----------------------------------------------------------

#  .  User input -----------------------------------------------------------
#further below code checks if data contains threshold for choices itself and uses that if so
h0 = 0.5#reference choice rate 

h1 = 'greater'#alternative hypothesis, mean "greater" than reference or "less"
csv_sep = ','#Is the csv comma separated or semicolon separated? For tab sep, use "\t"


#Check the operating system and assign a logical flag (T or F)
sys_win <- Sys.info()[['sysname']] == 'Windows'
#On computers set up by JMU Würzburg, use user profile instead of home directory
if(sys_win){
  #get rid of all the backslashes
  ltp <- gsub('\\\\', '/', Sys.getenv('USERPROFILE'))#Why does windows have to make this so difficult
}else{#Root directory should be the "HOME" directory on a Mac (or Linux?)
  ltp <- Sys.getenv('HOME')#Life was easier on Mac
}

# . Select files ---------------------------------------------------------


# manually set path file
wd="D:/Nextcloud/Home/Behaviour/PatternDiscrimination/Bee_Pattern_vs_Colour/col_vs_pattern/"; add=""

filename <- 'naiveChoices'
# filename <- 'naiveChoices_shape'


# Read in file ------------------------------------------------------------
#read in xlsx file with multiple sheets
path_file <- paste0(wd, add,filename,'.xlsx')
mdata = read_excel(path_file,sheet = "All")

# get proportion of choice1 out of all
mdata$choice1R <- mdata$choice1/(mdata$choice1 + mdata$choice2)


temp1 <- strsplit(mdata$ColPat1, "")

mdataNew = within(mdata,
               {
                 col1 = sapply(X = ColPat1,
                               FUN = function(x)
                                 {
                                 strsplit(x, '')[[1]][1]
                               }
                               )
                 patshap = sapply(X = ColPat1,
                               FUN = function(x)
                                 {
                                 strsplit(x, '')[[1]][2]
                               }
                               )
                 col2 = sapply(X = ColPat2,
                                  FUN = function(x)
                                  {
                                    strsplit(x, '')[[1]][1]
                                  }
                 )
               }
               )


# switch conditions in data so that they are all unique, and not split into O-Y vs Y-O for example
# hard code this for now, as it affects only the orange and yellow

indO <- which(mdataNew$col1 == "O")
for (i in 1:length(indO)) {
  mdataNew$col1[indO[i]] <- "Y"
  mdataNew$col2[indO[i]] <- "O"
  mdataNew$choice1R[indO[i]] <- 1-mdataNew$choice1R[indO[i]]
  if (mdataNew$patshap[indO[i]] == "K") {
    mdataNew$patshap[indO[i]] <- "S"
  }else {
    mdataNew$patshap[indO[i]] <- "K"
  }
}

#make a new variable which combines the colour conditions
mdataNew = within(mdataNew,{colourCond=paste0(col1,col2)})

# reorder conditions to fit the order in matlaplot.spec.coherency
reorderCons = c("BY", "BO", "BT", "YO")
mdataNew$colourCond <- factor(mdataNew$colourCond, levels=reorderCons) # reorder factors in open, semi, closed

# add animal ID as factor
# mdataNew$AnimalID <- factor(mdataNew$AnimalID) 


# plot data

setEPS()
postscript(paste0(wd, filename,add,'_patShap.eps'))
par(mar = c(2, 12, 1, 9))  # Set margins (bottom, left, top, right)
with(mdataNew, boxplot(choice1R ~ patshap,ylim=c(0,1))) 
with(mdataNew, beeswarm(choice1R ~ patshap, add = TRUE,
                     col = 3,   # Color
                     pch = 19,  # Symbol
                     cex = .7,
                     corral = "random"))
abline(h=c(0,1,0.5),lty=c(1,1,2))
dev.off()

setEPS()
postscript(paste0(wd, filename,add,'_fullInteract.eps'))
#plot shape / pattern choices within each colour combination
with(mdataNew, boxplot(choice1R ~ patshap * colourCond,ylim=c(0,1))) 
with(mdataNew, beeswarm(choice1R ~ patshap * colourCond, add = TRUE,
                        col = 3,   # Color
                        pch = 19,  # Symbol
                        cex = .5,
                        corral = "random"))
abline(h=c(0,1,0.5),lty=c(1,1,2))
dev.off()

#plot colour conbdition only
setEPS()
postscript(paste0(wd, filename,add,'_col.eps'))
par(mar = c(2, 9, 1, 6))  # Set margins (bottom, left, top, right)
with(mdataNew, boxplot(choice1R ~ colourCond,ylim=c(0,1))) 
with(mdataNew, beeswarm(choice1R ~ colourCond, add = TRUE,
                        col = 3,   # Color
                        pch = 19,  # Symbol
                        cex = .5,
                        corral = "random"))
abline(h=c(0,1,0.5),lty=c(1,1,2))
dev.off()

# plot all interactions present in data

# this lets us drop all interactions that are not present
mdataNew$interaction_var <- interaction(mdataNew$patshap, mdataNew$col1, mdataNew$col2, drop = TRUE)
boxplot(choice1R ~ interaction_var, data = mdataNew, main = "Boxplot with Interactions",ylim=c(0,1))
abline(h=c(0,1,0.5),lty=c(1,1,2))

with(mdataNew, beeswarm(choice1R ~ interaction_var, add = TRUE,
                        col = 3,   # Color
                        pch = 19,  # Symbol
                        cex = .5))

# plot colour interactions present in data

# this lets us drop all interactions that are not present
mdataNew$interaction_var <- interaction(mdataNew$col1, mdataNew$col2, drop = TRUE)
boxplot(choice1R ~ interaction_var, data = mdataNew, main = "Boxplot with Interactions",ylim=c(0,1))
abline(h=c(0,1,0.5),lty=c(1,1,2))

with(mdataNew, beeswarm(choice1R ~ interaction_var, add = TRUE,
                        col = 3,   # Color
                        pch = 19,  # Symbol
                        cex = .5))

# #read in csv file 
# path_file <- paste0(wd, add,filename,'.csv')
# mdata = read.table(file = path_file,#read from user-selected file
#                     header = T,#read the file header to use for variable names
#                     sep = csv_sep#,#values are separated by the user-specified character
#                     #other parameters can be added here for troubleshooting
#                     )


# Model -------------------------------------------------------------------


#lmm with data and all interactions
# lmm.all <- lmer(choice1R ~ patshap * col1 * col2 + (1 + patshap * col1 * col2 | AnimalID), data = mdataNew)
# no random effects (as animals in all conditions are different) - and too many random effects for number of observations
# lmm.all <- lm(choice1R ~ patshap * col1 * col2, data = mdataNew) 

# need a binomial model because data between 0 and 1
# lmm.all <- glm(choice1R ~ patshap * col1 + patshap * col2, data = mdataNew,family = quasibinomial()) 

#make a model to test for patshap choice, given the colour conditions
lmm.all <- glm(choice1R ~ patshap * colourCond, data = mdataNew,family = quasibinomial()) 


# glmm.all <- glmer(choice1R ~ patshap * col1 * col2 + (1 + patshap * col1 * col2 | AnimalID), data = mdataNew,family = Gamma)

#check normality of residuals of selected model (only relevant for Gaussian models)
qqnorm(residuals(lmm.all))
qqline(residuals(lmm.all))

plot(density(residuals(lmm.all)), col = 'darkblue', main = 'Residuals', xlab = 'Residual size', ylab = 'Frequency Density', xlim = c(-10, 10))
# plot(density(summary(lmm.allbg)$residuals), col = 'darkblue', main = 'Residuals', xlab = 'Residual size', ylab = 'Frequency Density', xlim = c(-10, 10))
lines(rep(0,2), c(0,1), lwd = .25)


#also check that this model explains the data at all
lmm.null <- glm(choice1R ~ patshap, data = mdataNew,family = quasibinomial())#model with only random effects, how well is the data explained by just splitting it into individuals and ignoring all experimental manipulation?
AIC(lmm.all)
AIC(lmm.null)


# # for the sake of publications, you can test for change in "deviance"
print(" ")
print("Model Comparison with null model")
print(anova(lmm.all, lmm.null,test="Chisq"))
#this does a chi-squared test on whether deviance between the model and the data is significantly lower for lmm.all
#the degrees of freedom are the difference in number of factors (1)
# the statistic is a deviance and it gives a p-value that reviewers like to see
#Now you can do some post-hoc tests

#for completion put here again
summary(lmm.all)

# print results of post-hoc tests for relevant factors (include * bg Color if tests suggest it is a relevant factor)
resultEMS <- emmeans(lmm.all, ~ patshap * colourCond,type="response")
summary(resultEMS)

print(summary(emmeans(lmm.all, list(pairwise ~ patshap * colourCond))))
print(summary(emmeans(lmm.all, list(pairwise ~ colourCond))))
print(summary(emmeans(lmm.all, list(pairwise ~ patshap),type="response")))
print(summary(emmeans(lmm.all, ~ colourCond,type="response")))
colCond_tst = test(emmeans(lmm.all, ~ colourCond,type="response"))

# colCond_tst = within(colCond_tst, sig = ifelse(p.value < 0.05, '*', ''))

#test if subsets of data defined by condition are sig different from threshold (0.5)
emm <- emmeans(lmm.all, ~ colourCond)
test_results <- test(emm, type = "response", null = 0)
summaryChoiceFreqCond <- summary(test_results)
summary(test_results)
