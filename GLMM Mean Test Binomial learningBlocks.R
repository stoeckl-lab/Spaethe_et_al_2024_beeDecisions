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

# Input Variables ----------------------------------------------------------

#  .  User input -----------------------------------------------------------
#further below code checks if data contains threshold for choices itself and uses that if so
# h0 = 0.25#reference choice rate line sectors
h0 = 0.5#reference choice rate cross sectors

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


# # manually set path file
wd="D:/Nextcloud/Home/Behaviour/PatternDiscrimination/Bee_Pattern_vs_Colour/col_vs_pattern/"; add=""

filename <- 'allLastBlock' #to compare last block of learning data
# filename <- 'conflictChoicesAll_single' #to compare test data from single attributes
#filename <- 'naiveChoices_All_single'

# wd="D:/Nextcloud/Home/Behaviour/PatternDiscrimination/Bee_Pattern_vs_Colour/col_vs_pattern/blockTests/"; add=""
# 
# filename <- 'allBlock_pattern_6' #to compare either block of pattern data
# filename <- 'allBlock_shape_1'  #to compare either block of shape data
# filename <- 'allBlock_single_5' #to compare either block of pattern data


# Read in file ------------------------------------------------------------
#read in xlsx file with multiple sheets
path_file <- paste0(wd, add,filename,'.xlsx')
mdata = read_excel(path_file)

# #read in csv file 
# path_file <- paste0(wd, add,filename,'.csv')
# mdata = read.table(file = path_file,#read from user-selected file
#                     header = T,#read the file header to use for variable names
#                     sep = csv_sep#,#values are separated by the user-specified character
#                     #other parameters can be added here for troubleshooting
#                     )

# View(mdata)#show the user the data that was
#find threshold if present in data
mdata$threshold=0.5 #if threshold is not defined
if("threshold" %in% colnames(mdata))
{
  #h0 <- mdata[1,5];
  h0 <- mdata$threshold[1];
}


#TODO, make this more flexible
mdata$Moth.number <- as.factor(mdata$AnimalID)#column named animalID contains ID information
mdata$colourCond <- as.factor(mdata$group)

#subset data (either only orange-blue or teal-blue for last block comparions)
# mdata <- mdata[ which(mdata$colourCond=="1" | mdata$colourCond=="2" | mdata$colourCond=="3"), ]; choiceType="1" #for orange blue
mdata <- mdata[ which(mdata$colourCond=="4" | mdata$colourCond=="5" | mdata$colourCond=="6" | mdata$colourCond=="7" | mdata$colourCond=="8"), ]; choiceType="2" #for blue blue

#generate a new attribute, close or distant colour, which is used as the colour condition
#for comparisons of learning blocks
# mdata$colourCond <- factor(ifelse(mdata$group %in% c(1, 2), "distant","close" ))

#for final test of single attribute data
choiceType=" " #for 

#need to drop the other levels of the factor, because otherwise retains them for plotting
mdata$colourCond <- droplevels(mdata$colourCond)

#for naive choices of single test, reorder conditions
# reorderCons = c("BO", "TB", "RK", "SC")
# mdata$colourCond <- factor(mdata$colourCond, levels=reorderCons) # reorder factors in open, semi, closed

# Fit model ---------------------------------------------------------------

# 
mdata$choice1 <- mdata$correctChoices  #select what to treat as choice option 1
mdata$choice2 <- mdata$incorrectChoices #select what to treat as choice option 2


frm.1 = formula(#response variable is a pair of columns (correct, incorrect)
          cbind(choice1, choice2) ~ #responses are binomial choice counts
            1 + #the predictor is a single mean
                (1|Moth.number) # individuals have their own random-effects means
               )

frm.1bg = formula(#response variable is a pair of columns (correct, incorrect)
  cbind(choice1, choice2) ~ #responses are binomial choice counts
    1 + colourCond + #the predictor is a single mean, CHANGED IT TO COLOUR< WAS 1 BEFORE
    (1 | Moth.number) # individuals have their own random-effects means
)


mod.1 = glmer(formula = frm.1,#fit this formula to a _generalised_ linear model
               data = mdata,#using this data (i.e. header names match formula variables)
               family = binomial(link = 'logit')#mean and s.d. are fitted on logit scale https://en.wikipedia.org/wiki/Logit 
               )

mod.1bg = glmer(formula = frm.1bg,#fit this formula to a _generalised_ linear model
              data = mdata,#using this data (i.e. header names match formula variables)
              family = binomial(link = 'logit')#mean and s.d. are fitted on logit scale https://en.wikipedia.org/wiki/Logit
)



#N.B. this is the simplest possible model including random effects,
#so no model selection is necessary

#inspect model
summary(mod.1)
summary(mod.1bg)
# print("Model Comparison with and without background")

anovaAllBg <- anova(mod.1, mod.1bg)
print(anovaAllBg)
#AIC:   Akaike Information Criterion, lower = better
#Random effects; Std. Dev.: should be >0 or individual effects were not fitted
#Fixed effects; Estimate: logit(mean success rate), 1/(1+exp(-Estimate)) = mean success rate

#perform post-hoc contrasts for the colour condition
emm <- emmeans(mod.1bg, ~ colourCond)
pairwise <- pairs(emm, adjust = "sidak")
summaryCondPosthoc <- summary(pairwise)

#test if subsets of data defined by condition are sig different from threshold (0.5)
test_results <- test(emm, type = "response", null = 0, side = ">")
summaryChoiceFreqCond <- summary(test_results)
summary(test_results)
#there is a conversion problem between linear scale and response scale (logit),
# and alledgedly the null that is put in is interpreted on the response scale (in this case probability),
# but in the summary it is not, it is shown as the conversion in linear scale (so 0.5 turns to 0.622).
# 0 input corresponds to 0.5 null value in the summary result with response scale, so this makes me think
# it is input in linear scale = 0 turns into null in logit space = 0.5

#write results into txt file
summaryCondPosthoc_text <- capture.output(summaryCondPosthoc)
summaryChoiceFreqCond_text <- capture.output(summaryChoiceFreqCond)

file_path_summary <- paste0(wd,filename,choiceType,"_binModel.txt")
writeLines(c(summaryCondPosthoc_text, summaryChoiceFreqCond_text), file_path_summary)

#plot choice data
mdata$relChoice <- mdata$choice1 / (mdata$choice1+mdata$choice2)

ggplot(mdata, aes(x = factor(colourCond), y = relChoice)) +
  geom_boxplot() +
  labs(x = "Group", y = "Choice Frequency") +
  ggtitle("Relative Choice Frequencies by Group") +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), alpha = 0.5)

setEPS()
postscript(paste0(wd, filename,add,choiceType,'_comp.eps'))
par(mar = c(2, 4, 1, 1))  # Set margins (bottom, left, top, right)
with(mdata, boxplot(relChoice ~ colourCond,ylim=c(0,1)))
with(mdata, beeswarm(relChoice ~ colourCond, add = TRUE,
                        col = 3,   # Color
                        pch = 19,  # Symbol
                        cex = .7,
                        corral = "random"))
abline(h=c(0,1,0.5),lty=c(1,1,2))
dev.off()



#now how do I analyse all the individual conditions' results and whether the choices are different from 0.5??

#check if full model with colour fits better, then use that for further analysis
if (anovaAllBg[2,8]<0.05 && anovaAllBg[2,5]<anovaAllBg[1,5])
{mod.1 <- mod.1bg}

#  .  Compare between conditions	----------------------------------------------
glm.all <- glm(relChoice ~ colourCond, data = mdata,family = quasibinomial()) 
#check normality of residuals of selected model (only relevant for Gaussian models)
qqnorm(residuals(glm.all))
qqline(residuals(glm.all))
# do posthoc tests
print(summary(emmeans(glm.all, list(pairwise ~ colourCond))))
summaryGLM_text <- capture.output(summary(emmeans(glm.all, list(pairwise ~ colourCond))))

file_path_summary <- paste0(wd,filename,choiceType,"_glm.txt")
writeLines(summaryGLM_text, file_path_summary)
# 

#  .  Derive variables	----------------------------------------------

#transform fixed effects mean to odds, logit is equivalent to log(odds)
exp(fixef(mod.1))#X correct choices for every incorrect choice
#transform fixed effects mean to logistic probability (proportion correct)
mod.mu = fixef(mod.1)
message('mean cross choice rate\n',
        round(plogis(mod.mu), 3)*100,
        '%')
#extract fixed-effects standard deviation
mod.sigma = coef(summary(mod.1))[ , "Std. Error"]
#N.B. this can also be calculated from the mode, and for each individual including individual effects
    # modmat  =  model.matrix(terms(mod.1), mdata)#extract model matrix
    # fixedfx = modmat %*% fixef(mod.1) #calulate fixed effects
    # pvar = diag(modmat %*% tcrossprod(vcov(mod.1), modmat)) #use variance-covariance matrix to calculate model variance
    # tvar = pvar+as.numeric(VarCorr(mod.1))
    # p_upper = fixedfx + 1.96*sqrt(tvar) 
    # p_lower = fixedfx - 1.96*sqrt(tvar)
#extract fixed-effects confidence interval (explained below)
#confidence interval
mod.fixef.ci = qnorm(c(0.05/2,1-0.05/2),
                    mean = mod.mu,
                    sd = mod.sigma
                    )
message('two-tailed 95% confidence interval for mean cross choice rate\n',
        round(plogis(mod.fixef.ci[1]), 3)*100,
        ' - ',
        round(plogis(mod.fixef.ci[2]), 3)*100,
        '%')

# Calculate p value -------------------------------------------------------

#Logistic regression assumes a normal distribution on the logit scale
#The fitted model has the probability density function
plot(seq(-5,5,0.01),
     dnorm(seq(-5,5,0.01),
           mean = mod.mu,
           sd = mod.sigma),
     type = 'l',
     xlab = 'log(odds)',
     ylab = 'probability density'
     )
#For which the cumulative probability of the true mean
# being at a specific value, or lower, is
plot(seq(-5,5,0.01),
     pnorm(seq(-5,5,0.01),
         mean = mod.mu,
         sd = mod.sigma),
     type = 'l',
     xlab = 'log(odds)',
     ylab = 'probility of true mean ≤ value'
     )
abline(h = 0.05,
       col = 'red',
       lty = 2#dashed line
       )
#below the 5th percentile for this normal distribution (quantile function),
#the probability of a true mean at that value or lower is less than 0.05
abline(v = qnorm(0.05,
                 mean = mod.mu,
                 sd = mod.sigma),
       col = 'red',
       lty = 1#solid line
       )

#Convert reference value to logit space
h0.logit = qlogis(h0)#quantile function for the logistic distribution performs logit tranform
p.h0 = ifelse(h1 == 'greater',#if alternative hypothesis that true mean is greater than reference
              pnorm(h0.logit,#cumulative probability up to reference
                    mean = mod.mu,#distribution with model mean
                    sd = mod.sigma#and model standard deviation
                    ),
              1-pnorm(h0.logit,#cumulative probability down to reference (1-p)
                      mean = mod.mu,#distribution with model mean
                      sd = mod.sigma#and model standard deviation
                      )
              )
message('Probability mean correct choices are NOT ',
        h1,
        ' than ',
        h0,
        '\n p = ',
        round(p.h0,5)
        )
#calculate z-score
z_score = (qlogis(h0)-mod.mu)/mod.sigma
message('z-score = ',
        round(z_score,3)
        )
#save result
write.table(x = cbind(h0 = h0,#save relevant info in a table
                      h1 = h1,
                      mean_choice = plogis(mod.mu),
                      ci_02.5 = plogis(mod.fixef.ci[1]),
                      ci_97.5 = plogis(mod.fixef.ci[2]),
                      z_score = z_score,
                      p = p.h0,
                      mod1 = c("noBG"),
                      mod2 = c("BG"),
                      Dev1 = anovaAllBg[1,5],
                      Dev2 = anovaAllBg[2,5],
                      anovaChi = anovaAllBg[2,6],
                      anovap = anovaAllBg[2,8]
                      ),
            file = file.path(dirname(path_file),#same place as original with
                             paste(sub(pattern='.csv',#similar name
                                       replacement = '',
                                       x = basename(path_file)
                                       ),
                                  'GLMM mean test.csv'#ending in mean test
                                  )
                             ),
            row.names = FALSE,#rows do not need names
            sep = csv_sep #Use same separator as original csv file
  
)

# Plot mean and confidence intervals --------------------------------------

just_path = dirname(path_file)#folder
file_name = basename(path_file)#filename
path_to_newfile = file.path(just_path, paste0(file_name,'.pdf'))#original filename + '.pdf'
#open a pdf
pdf(file = path_to_newfile)#open file

#fixed effects model in probability space
#mean correct
barplot(height = plogis(mod.mu),
        width = 0.1,
        xlim = c(0,1),
        ylim = c(0,1),
        space = 0.5,
        main = paste0('p(mean',
                      ifelse(h1=='greater','≤','≥'),
                      'reference) = ',
                      round(p.h0,3)
                      )
        )
#reference correct choice rate
abline(h = h0,
       col = 'red'
       )
#Fixed effects 95% confidence interval for the mean
#N.B. This confidence interval is two-tailed, it extends to the 2.5th percentile
# whereas our test is one-tailed, its "confidence interval" extends to the 5th percentile
arrows(x0 = 0.1,
      x1 = 0.1,
      y0 = plogis(mod.fixef.ci[1]),
      y1 = plogis(mod.fixef.ci[2]),
      code = 3,
      angle = 90,
      length = 0.1
      )
#add mean back for appearances
points(x = 0.1,
       y = plogis(mod.mu),
       pch = 19#a small solid dot
       )

dev.off()#close file to save