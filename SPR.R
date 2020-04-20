#This assumes data collected in IbexFarm, https://github.com/addrummond/ibex. However, the analysis and plotting are more widely applicable.

#The bayesian analysis uses brms, a wrapper for STAN: https://github.com/paul-buerkner/brms


setwd("mywd")
library(ggplot2)
library(ggthemes)
library(magrittr)
library(dplyr)
library(plotrix)
library(reshape)
library(blmer)
library(rstan) 
rstan_options(auto_write=TRUE)
library(rstanarm)
options(mc.cores=parallel::detectCores())
library(stringr)


#read in IbexFarm output
allData_SPR <- read.csv('IbexResults.csv', header = 0, comment.char = "#")
#this uses the submission time as Participant ID bc it was run in a lab
colnames(allData_SPR) <- c("Participant", "MD5","Controller", "Group", "Element", "Condition", "Item", "Seq", "Word","RT","Newline","FullSent")

#remove any participants necessary
allData_SPR3 <- subset(allData_SPR3, allData_SPR3$Participant != 3)


#pull exp. items
#this assumes 4 conditions as in a 2x2 factorial design: NL, ND, EL, ED
doi_SPR <- subset(allData_SPR, allData_SPR$Condition == "NL" | allData_SPR$Condition == "ND" | allData_SPR$Condition == "EL" | allData_SPR$Condition == "ED")

#pull dashed sentences
#pulls only the self-paced reading sentences using the IbexFarm controller label
dash_SPR <- filter(doi_SPR, Controller == "DashedSentence")

#look at a summary of the items by experimental condition
dash_SPR %>% 
  group_by(Item, Condition) %>% 
  summarize(n=n()) -> items_SPR_summary


#Look at question performance, can use to exclude participants
#this assumes comprehension questions were used as well as SPR sentences
#pull questions using the IbexFarm question controller
questions_SPR <- filter(doi_SPR, Controller == "Question")
#label the question columns, which are not the same as the DashedSentence columns
colnames(questions_SPR) <- c("Participant", "MD5.Q","Controller.Q", "Group.Q", "Element.Q", "Condition", "Item", "FullQues", "Resp","Correct","RT.Q","Blank")
questions_SPR <- questions_SPR[c(1:11)]

#includes null responses (time-outs)
#look at overall by-participant performance
questions_SPR %>%
  group_by(Participant, Correct) %>%
  summarize(n=n()) -> questions_SPR_byPart

#remove null responses
questions_SPR_noNA <- subset(questions_SPR, questions_SPR$Correct != "NULL")
questions_SPR_noNA %>%
  group_by(Participant, Correct) %>%
  summarize(n=n()) -> questions_SPR_byPart_noNA

#use spread to reshape data to get average correct by participant
questions_SPR_byPart_noNA <- spread(questions_SPR_byPart_noNA, Correct, n)
questions_SPR_byPart_noNA$total <- questions_SPR_byPart_noNA$"0" + questions_SPR_byPart_noNA$"1"
questions_SPR_byPart_noNA$percCorr <- questions_SPR_byPart_noNA$"1" / questions_SPR_byPart_noNA$total

#combine both dataframes into one large dataframe called 'dataPlaus_SPR'
#this can be helpful for analysis purposes
#remember here the two dataframes need the same label for whatever the 'by' label is
dataPlaus_SPR <-merge(questions_SPR, dash_SPR, by=c("Participant","Item", "Condition"))

#end looking at participant question data



#label critical regions
#this is only if your critical regions are not the same Seq across items/conditions
#remember that SPR indexes Seq from 1, while Maze indexes from 0
dash_SPR$Region <- if_else(dash_SPR$Seq == 6 & (dash_SPR$Condition == "ND" | dash_SPR$Condition == "NL"), 2, dash_SPR$Region)
dash_SPR$Region <- if_else(dash_SPR$Seq == 7 & (dash_SPR$Condition == "ND" | dash_SPR$Condition == "NL"), 1, dash_SPR$Region)

dash_SPR$Region <- if_else(dash_SPR$Seq == 5 & (dash_SPR$Condition == "ED" | dash_SPR$Condition == "EL"), 2, dash_SPR$Region)
dash_SPR$Region <- if_else(dash_SPR$Seq == 6 & (dash_SPR$Condition == "ED" | dash_SPR$Condition == "EL"), 1, dash_SPR$Region)

#check regions
dash_SPR %>% 
  group_by(Region, Word) %>% 
  summarize(n=n()) -> regionCheck_SPR

#helpful to create a new df in case something goes wrong from here on out
usableDOI_SPR <- dash_SPR


#TRIM DATA
#remove sequence 0, since it's not indicative of actual reading time
usableDOI_SPR <- subset(usableDOI_SPR, usableDOI_SPR$Seq != 0)

usableDOI_SPR$RT <- as.numeric(as.character(usableDOI_SPR$RT))
is.numeric(usableDOI_SPR$RT)

#trim off the top .3% of the RTs
#this is an arbitrary cut-off in some sense
upperTrim_SPR <- quantile(x = usableDOI_SPR$RT, probs = .997)
upperTrim_SPR
usableDOI_SPR <- subset(usableDOI_SPR, usableDOI_SPR$RT < 2287.568)  

#trim fast end: anything below 100ms
#this is physically-bounded by possible reaction times and is the standard cut-off in single-word reaction time studies
usableDOI_SPR <- subset(usableDOI_SPR,usableDOI_SPR$RT > 99)


#make factor levels
#this uses contrast coding based on Condition labels, which I use to represent factor levels
usableDOI_SPR$Condition %>%
substr(1,1) %>% factor -> usableDOI_SPR$Ellipsis
contrasts(usableDOI_SPR$Ellipsis) <- contr.sum(2) * (1/2)
#defaults at 1 and -1, multiplying by .5 gives you contrasts .5 and -.5
levels(usableDOI_SPR$Ellipsis) <- c("Ellipsis", "NoEllipsis")
#sets your level labels

usableDOI_SPR$Condition %>%
substr(2,2) %>% factor -> usableDOI_SPR$Locality
contrasts(usableDOI_SPR$Locality)  <- contr.sum(2) * (-1/2)
#-.5 here will switch the sign of your label contrasts
#can always check your contrast labelling using: contrasts(usableDOI_SPR$Locality)
levels(usableDOI_SPR$Locality) <- c("Distant", "Local")


#find region means (using Region 2 here as an example)

region2_SPR <- subset(usableDOI_SPR, usableDOI_SPR$Region == 2)

#calculate means by item
region2_SPR %>% 
  group_by(Locality, Ellipsis, Item) %>% 
  summarize(m=mean(RT), med=median(RT), sd = sd(RT), sd2 = sd(RT)*2, n=n(), se=std.error(RT)) -> summary.by.items_R2_SPR

#calculate means for condition using by-item means
summary.by.items_R2_SPR %>% 
  group_by(Locality, Ellipsis) %>% 
  summarize(RT=mean(m), med=median(med), sd = sd(m), sd2 = sd(m)*2, n=n(), se=sd/sqrt(n)) -> summary.R2_SPR

summary.R2_SPR$Region <- "NP"

#plot means
png("myPlot.png", width=20, height=10, units="cm", res=600, bg="transparent")

summary.R2_SPR %>% ggplot(aes(y=RT)) -> SPR_R2

SPR_R2 + geom_pointrange(aes(x=Ellipsis, shape=Locality, col=Locality,ymin=RT-se, ymax=RT+se), 
                         position=position_dodge(width=0.4)) +
  theme_minimal(base_family="Times New Roman") + scale_shape_circlefill() + scale_color_manual(values=c("#09406E","#027495")) + 
  ylim(300,600) + ylab("Mean RT") + labs(title="Region 2") +
  labs(caption="SPR Experiment, Region 2. Error bars indicate standard error of the mean.") +
  theme(axis.text.x=element_text(size=6)) + theme(panel.grid.minor=element_blank())
dev.off()


#run Bayesian analysis using brms
#https://github.com/paul-buerkner/brms

#good to double-check your contrasts bc dplyr manipulations can mess them up
contrasts(region2_SPR$Ellipsis)
contrasts(region2_SPR$Locality)

#always run a fixed effects only model first
#then compare the results with the full model to ensure no craziness has occurred
# see 10.1016/j.jml.2012.11.001 for why using models with full-random effects structure matters
brmR2_SPR <- brm(RT ~ Locality*Ellipsis,  #fixed effects only
                 data = region2_SPR, 
                 family = shifted_lognormal(link=identity),
                 control = list(adapt_delta = 0.8, max_treedepth = 10),
                 warmup = 1000, iter = 5000, chains = 4) #good defaults, can change if warnings occur
summary(brmR2_SPR, waic = TRUE)

#note: if you increase adapt_delta to .99 and still get divergence warnings, the warnings will send you to STAN
#help pages telling you to re-paramaterize your model. brms already runs this parameter as a default, so your 
#problem is elsewhere. try using a different family to fit your model.
#gaussian is a good family for log-transformed RTs
#exgaussian and shifted_lognormal are good families for raw RTs (account for the gaussiandistribution+exponential tail of reaction times)

brmR2_SPR <- brm(RT ~ Locality*Ellipsis + #fixed effects
                        (0 + Locality*Ellipsis|Participant) + #random Participant slopes
                        (0 + Locality*Ellipsis|Item) +  #random Item slopes
                        (1|Participant) + #random Participant intercept
                        (1|Item),   #random Item intercept
                      data = region2_SPR, 
                      family = shifted_lognormal(link=identity),
                      control = list(adapt_delta = 0.8, max_treedepth = 10),
                      warmup = 1000, iter = 5000, chains = 4)
summary(brmR2_SPR, waic = TRUE)
plot(brmR2_SPR) #get all your pretty plots

pp_check(brmR2_SPR,nsamples = 300) #check your posterior probability to see if your model matches the observed data


#Plot the RT results
#combine all region means together into a single dataframe
region_RTs_SPR <- full_join(summary.RN2_SPR, summary.RN1_SPR)
#...
region_RTs_SPR <- full_join(region_RTs_SPR, summary.R0_SPR)


png("myRTmeans.png", width=31.75, height=17.46, units="cm", res=800, bg="transparent")

#plot RTs from each Region
#assumes you computed means for each region as given above for Region 2
#thanks to Steven Foley for sharing his version of this code
# Reorder factors so they don't appear alphabetically
region_RTs_SPR$Region %<>% factor(levels=c('MC-2', 'MC-1', 'ClauseEnd', 'Quant','NP', 'be', 'Adverb', 'Spill'))

region_RTs_SPR$Locality %<>% factor(levels=c('Distant','Local'))

# Make an all-caps legend for the Ellipsis factor
Ellipsis_names <- c('Ellipsis' = "ELLIPSIS", 'NoEllipsis' = "NO ELLIPSIS")

region_plot_NPE3_SPR <- region_RTs_SPR %>% 
  ggplot(aes(x=Region, y=RT,
             col=Ellipsis, group=Locality, alpha=Locality))
region_plot_NPE3_SPR + 
  # Facets will be labelled with what we defined
  facet_grid(Ellipsis ~ ., labeller = as_labeller(Ellipsis_names)) +
  # Fiddle with the point/line sizes to make them an appropriate size
  # This gives points at mean RTs and confidence intervals around them
  geom_pointrange(size=1, aes(ymin=RT-se,
                              ymax=RT+se)) +
  geom_line(lwd=1) + theme_minimal() + 
  scale_shape_circlefill() +
  ylab("Mean RTs (ms)") +
  theme(text=element_text(size=18, family = "Times New Roman"), 
        element_line(size = 2),
        panel.grid = element_line(size=0.5),
        axis.title.x = element_blank()) +
  # Define pretty colors
  scale_color_manual(values=c("#A50D69","#7A0D87","#500DA5")) +
  # Alpha is transparency
  scale_alpha_manual(values=c(1, 0.5),
                     name="Locality",
                     breaks=c("Distant","Local"),
                     labels=c("Distant","Local")) +
  guides(color=FALSE) +
  labs(title ="My RT Data")

dev.off()





