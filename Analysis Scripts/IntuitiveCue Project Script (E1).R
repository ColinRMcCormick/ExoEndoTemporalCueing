# ----------------------------------- #
#    Temporal Attention Experiment
#     Discrimination temporal orienting
#       Intuitive Cues
#  N=40 (goal)
# ----------------------------------- #

#Current WD#
setwd()
#    Load Librarys
# ----------------------------------- #

library(tidyverse)
library(plyr)
library(calibrate)
library(dplyr)
library(ez)
library(readr)
library(stringr)
library(tidyr)
library(tidyverse)
library(ggpubr)
library(lme4)

###
# Read & process the trigger data ----
###
trigger_files = list.files(
  pattern = 'trigger'
  , recursive = TRUE
  , full.names = TRUE
)

data_filesrigger_files = trigger_files[!str_detect(trigger_files,'test')]

trigger_data = map_df(
  .x = trigger_files
  , .f = read_tsv
  , col_types = cols(
    minute = col_integer()
    , age = col_character()
    , hour = col_character()
    , id = col_character()
    , month = col_integer()
    ,day = col_character()
    ,sex = col_character()
  )
)

trigger_data = trigger_data[trigger_data$block!='practice',]
trigger_data$block = as.numeric(trigger_data$block)

trigger_data$value[trigger_data$value<(-1)] = -1
trigger_data$value[trigger_data$value>1] = 1
trigger_data$value = (trigger_data$value+1)/2*100
trigger_data$time = trigger_data$time*1000

trigger_stats = ddply(
  .data = trigger_data
  , .variables = .(id,block,trial_num,signal)
  , .fun = function(x){
    no_max = (sum(diff(sign(diff(x$value[x$trigger=='right'])))==-2)==0) & (sum(diff(sign(diff(x$value[x$trigger=='left'])))==-2)==0)
    both_triggers = length(unique(x$trigger))>1
    double_trigger = (sum(diff(sign(diff(x$value[x$trigger=='right'])))==-2)>1)|(sum(diff(sign(diff(x$value[x$trigger=='left'])))==-2)>1)
    start_slope = NA
    start_r = NA
    end_slope = NA
    end_r = NA
    full = NA
    hold = NA
    auc = NA
    if((!no_max)&(!both_triggers)&(!double_trigger)){
      first_max_index = 1+which(diff(sign(diff(x$value)))==-2)[1]
      if(!is.na(first_max_index)){
        if(first_max_index>1){
          start = x[1:first_max_index,]
          start_fit = lm(value~time,data=start)
          start_slope = as.numeric(start_fit$coef[2])
          start_r = summary(start_fit)$r.squared
          end = x[(first_max_index+1):nrow(x),]
          end_fit = lm(value~time,data=end)
          end_slope = as.numeric(end_fit$coef[2])
          end_r = summary(end_fit)$r.squared
          full = x$time[nrow(x)]-x$time[1]
          hold = x$time[first_max_index+1]-x$time[first_max_index]
          auc = 0
          # for(i in 1:(nrow(x)-1)){
          # auc = auc+x$value[i]*(x$time[i+1]-x$time[i])+abs(x$value[i+1]-x$value[i])*(x$time[i+1]-x$time[i])/2
          # }
        }
      }
    }
    to_return = data.frame(
      no_max = no_max
      , both_triggers = both_triggers
      , double_trigger = double_trigger
      , start_slope = start_slope
      , start_r = start_r
      , end_slope = end_slope
      , end_r = end_r
      , full = full
      , hold = hold
      , auc = auc
    )
    return(to_return)
  }
  , .progress = 'time'
)
summary(trigger_stats)

start_r = trigger_stats$start_r
start_r = start_r[!is.na(start_r)]
start_r = start_r[start_r<1]
start_z = atanh(start_r)
start_lo = tanh(mean(start_z)-2*sd(start_z))
trigger_stats$bad_start = trigger_stats$start_r<start_lo

end_r = trigger_stats$end_r
end_r = end_r[!is.na(end_r)]
end_r = end_r[end_r<1]
end_z = atanh(end_r)
end_lo = tanh(mean(end_z)-2*sd(end_z))
trigger_stats$bad_end = trigger_stats$end_r<end_lo

summary(trigger_stats)

###
# Read in the regular data ----
###

data_files = list.files(
  pattern = 'data'
  , recursive = TRUE
  , full.names = TRUE
)

#read all data files to a single data frame
c = purrr::map_df(
  .x = data_files
  , .f = read_tsv
  , col_types = cols(
    minute = col_integer()
    , age = col_character()
    , id = col_character()
    , cue = col_character()
    , month =col_integer()
    , hour = col_character()
    , day = col_character()
    , sex = col_character()
  )
)

c = c[c$block!='practice',]
c$block = as.numeric(c$block)


c = dplyr::full_join(
  x = c
  , y = trigger_stats
)


#####
#   Exclusion Criteria
####
#before getting rid of error trials, store them here
c_error<-c

#applying criteria
#get rates of trials we'll toss
hist(c$rt)
###Anticipation Errors##
mean(c$pre_target_response)
#
c = subset(c, pre_target_response==FALSE)
###Misses###
mean(is.na(c$rt[!c$pre_target_response]))
#
c = c[!is.na(c$rt),]
###ITI###
mean(c$ITI_response,na.rm=T)
#
c = subset(c, ITI_response==FALSE)
###Both Triggers###
summary(c)
#
c = subset(c, both_triggers==FALSE)

#no NO_MAX remain
c = subset(c, no_max==FALSE)
c = subset(c, double_trigger!=TRUE)
c$rt<-c$rt *1000

#Reaction Time Cut Off
summary(c)
# Define the bin edges and labels
breaks <- c(-Inf, 50, 100, 150, 200, 210, 220, 230, 240, 250, 300, 350, 400, 450, 500,
            550, 600, 650, 700, 750, 800, 810, 820, 830, 840, 850, 900, 950, 1000, Inf)
labels <- c("50", "50-100", "100-150", "150-200", "200-210", "210-220", "220-230",
            "230-240", "240-250", "250-300", "300-350", "350-400", "400-450",
            "450-500", "500-550", "550-600", "600-650", "650-700", "700-750",
            "750-800", "800-810", "810-820", "820-830", "830-840", "840-850",
            "850-900", "900-950", "950-1000")

# Bin the data using cut()
c$rt_bin <- cut(c$rt, breaks = breaks, labels = labels, right = FALSE)

# Aggregate error by rt_bin
bin_er <- aggregate(error ~ rt_bin, data = c, mean)
bin_er1 <- aggregate(error ~ rt_bin, data = c, length)

# Merge results
a <- merge(bin_er, bin_er1, by = "rt_bin")
a
c<-subset(c, rt > 220)
###Too slow###
c<-subset(c, rt < 800)


#Exclusion Criteria
#---------------------------
#Making ErrorRate
c$error [c$error == TRUE] <- 1
c$error [c$error == FALSE] <- 0
c$error <- as.numeric(c$error)
c_error<-c
summary(c_error)
#All of these prior errors represent not true error performance. It's now all set
#up, so we'll cut the errors from only the C dataset
c <- subset (c, error!=TRUE)

##################3
####Checking usable trials per participant#######
##################3

summary(c)
chtr<-aggregate(rt~id, data= c_error, length )
chpar<- transform(chtr, completed_trials =rt/480)
chpar
c<-subset(c, id!='01309') # 60
c<-subset(c, id!='47523') # 56
c<-subset(c, id!='1000') # 47
c<-subset(c, id!='86457') # 36
c<-subset(c, id!='1007') # 66
c<-subset(c, id!='2146545') # 69
c<-subset(c, id!='12496') # 60
c<-subset(c, id!='91826') #56

summary(c_error)

c_error<-subset(c_error, id!='47523') #
c_error<-subset(c_error, id!='1000') #
c_error<-subset(c_error, id!='86457') #
c_error<-subset(c_error, id!='1007') #
c_error<-subset(c_error, id!='2146545') #
c_error<-subset(c_error, id!='01309') #
c_error<-subset(c_error, id!='12496') #
c_error<-subset(c_error, id!='91826')

####
#Separating short and long SOAs to allow for separate analyses#
####
ShSOAc<-subset(c, soa!=1.6)
LoSOAc<-subset(c, soa!=.4)
ShSOAc_er<-subset(c_error, soa!=1.6)
LoSOAc_er<-subset(c_error, soa!=.4)

#inverse reaction time
ShSOAc$rt_inv<-NA
ShSOAc$rt_inv<-1/ShSOAc$rt
LoSOAc$rt_inv<-NA
LoSOAc$rt_inv<-1/LoSOAc$rt


####################
#Reaction Time# (Short)
unrestrictM_target <- lmer( rt_inv ~target + (1|id), data = ShSOAc, REML = F ) #influence of target
unrestrictM_cue <- lmer( rt_inv ~cue + (1|id), data = ShSOAc, REML = F ) #influence of cue
unrestrictM_signal <- lmer( rt_inv ~ signal + (1|id), data = ShSOAc, REML = F ) #influence of signal
unrestrictM_int <- lmer( rt_inv ~ signal+ cue+ signal * cue + (1|id), data = ShSOAc, REML = F ) #full
restrictM <- lmer( rt_inv ~ (1|id), data = ShSOAc, REML = F ) #only participant error (for comparison to other ME models)
restrictM_interaction <- lmer( rt_inv ~ (1|id) +cue + signal, data = ShSOAc, REML = F )#everything but interaction (for comparison to interaction)
model_interact<-lmer( rt_inv ~ (1|id) +cue * signal, data = ShSOAc, REML = F )
###Doing AIC comparisons
(AIC(restrictM)-AIC(unrestrictM_target))*log2(exp(1))
(AIC(restrictM)-AIC(unrestrictM_cue))*log2(exp(1))
(AIC(restrictM)-AIC(unrestrictM_signal))*log2(exp(1))
(AIC(restrictM_interaction)-AIC(unrestrictM_int))*log2(exp(1))

#Calculating confidence intervals for main effects and interactions (ask RJ why)
#confint(restrictM_interaction, method = 'boot', .progress="txt", nsim=500)
#confint(model_interact, method = 'boot', .progress="txt", nsim=500)

#Error Rate (short)#
unrestrictM_target <- glmer( error ~target + (1|id), data = ShSOAc_er, REML = F, family = 'binomial' ) #influence of cue
unrestrictM_cue <- glmer( error ~cue + (1|id), data = ShSOAc_er, REML = F, family = 'binomial' ) #influence of cue
unrestrictM_signal <- glmer( error ~ signal + (1|id), data = ShSOAc_er, REML = F, family = 'binomial' ) #influence of signal
unrestrictM_int <- glmer( error ~ signal+ cue+ signal * cue + (1|id), data = ShSOAc_er, REML = F, family = 'binomial' ) #full
restrictM <- glmer( error ~ (1|id), data = ShSOAc_er, REML = F, family = 'binomial' ) #only participant error (for comparison to other ME models)
restrictM_interaction <- glmer( error ~ (1|id) +cue + signal, data = ShSOAc_er, REML = F, family = 'binomial' )#everything but interaction (for comparison to interaction)
model_interact<-glmer( error ~ (1|id) +cue * signal, data = ShSOAc_er, REML = F, family = 'binomial' )
(AIC(restrictM)-AIC(unrestrictM_target))*log2(exp(1))
(AIC(restrictM)-AIC(unrestrictM_cue))*log2(exp(1))
(AIC(restrictM)-AIC(unrestrictM_signal))*log2(exp(1))
(AIC(restrictM_interaction)-AIC(unrestrictM_int))*log2(exp(1))
summary(c_error)
#confint(restrictM_interaction, method = 'boot', .progress="txt", nsim=500)
#confint(model_interact, method = 'boot', .progress="txt", nsim=500)


#Long SOA RT
library(lme4)
unrestrictM_cue <- lmer( rt_inv ~cue + (1|id), data = LoSOAc, REML = F ) #influence of cue
unrestrictM_signal <- lmer( rt_inv ~ signal + (1|id), data = LoSOAc, REML = F ) #influence of signal
unrestrictM_int <- lmer( rt_inv ~ signal+ cue+ signal * cue + (1|id), data = LoSOAc, REML = F ) #full
restrictM <- lmer( rt_inv ~ (1|id), data = LoSOAc, REML = F ) #only participant error (for comparison to other ME models)
restrictM_interaction <- lmer( rt_inv ~ (1|id) +cue + signal, data = LoSOAc, REML = F )#everything but interaction (for comparison to interaction)
model_interact<-lmer( rt_inv ~ (1|id) +cue * signal, data = LoSOAc, REML = F )
###Doing AIC comparisons
(AIC(restrictM)-AIC(unrestrictM_cue))*log2(exp(1))
(AIC(restrictM)-AIC(unrestrictM_signal))*log2(exp(1))
(AIC(restrictM_interaction)-AIC(unrestrictM_int))*log2(exp(1))

#Long SOA Error Rate#
LoSOAc_er$error<- as.numeric(LoSOAc_er$error)
unrestrictM_cue <- glmer( error ~cue + (1|id), data = LoSOAc_er, REML = F, family = 'binomial' ) #influence of cue
unrestrictM_signal <- glmer( error ~ signal + (1|id), data = LoSOAc_er, REML = F, family = 'binomial' ) #influence of signal
unrestrictM_int <- glmer( error ~ signal+ cue+ signal * cue + (1|id), data = LoSOAc_er, REML = F, family = 'binomial' ) #full
restrictM <- glmer( error ~ (1|id), data = LoSOAc_er, REML = F, family = 'binomial' ) #only participant error (for comparison to other ME models)
restrictM_interaction <- glmer( error ~ (1|id) +cue + signal, data = LoSOAc_er, REML = F, family = 'binomial' )#everything but interaction (for comparison to interaction)
model_interact<-glmer( error ~ (1|id) +cue * signal, data = LoSOAc_er, REML = F, family = 'binomial' )
(AIC(restrictM)-AIC(unrestrictM_cue))*log2(exp(1))
(AIC(restrictM)-AIC(unrestrictM_signal))*log2(exp(1))
(AIC(restrictM_interaction)-AIC(unrestrictM_int))*log2(exp(1))


#######################################
#MEANS of REACTION TIME AND ERROR RATE#
#######################################
# Use participant means and not overall means #
ShSOAc_er <- ShSOAc_er %>% dplyr::filter(!is.na(error))
LoSOAc_er <- LoSOAc_er %>% dplyr::filter(!is.na(error))

#Short SOA reaction time/ER#
ezStats(
  data = ShSOAc
  , dv = .(rt)
  , wid = .(id)
  , within = .(signal, cue)

)
ezStats(
  data = ShSOAc_er
  , wid = id
  , dv = .(error)
  , within = .(signal, cue)
)
#Long SOA reaction time/ER
ezStats(
  data = LoSOAc
  , dv = .(rt)
  , wid = .(id)
  , within = .(signal, cue)

)
ezStats(
  data = LoSOAc_er
  , wid = id
  , dv = .(error)
  , within = .(signal, cue)
)

summary(LoSOAc_er)

#Target check
#Short#
#RT
ezStats(
  data = ShSOAc
  , wid = id
  , dv = .(rt)
  , within = .(target, cue)
)
ezStats(
  data = ShSOAc
  , wid = id
  , dv = .(rt)
  , within = .(target, signal)
)

#ER
ezStats(
  data = ShSOAc_er
  , wid = id
  , dv = .(error)
  , within = .(target, cue)
)
ezStats(
  data = ShSOAc_er
  , wid = id
  , dv = .(error)
  , within = .(target, signal)
)

#######################
#########Plots#########
########################
cbPalette <- c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#SHORT SOA RT
gr_a = ezPlot(
  data = ShSOAc
  , dv = .(rt)
  , wid = .(id)
  , within = .(signal, cue)
  ,do_lines = TRUE
  ,split= .(cue)
  , x = .(signal)
  , x_lab = 'Signal Type'
  , y_lab = 'Reaction Time (ms)'
  ,split_lab = 'Cues'

)
gr_a+ theme_pubr(base_size = 14, base_family = "", border = FALSE,
                 margin = TRUE, legend = c("top", "bottom", "left", "right", "none"),
                 x.text.angle = 0)+
  scale_colour_manual(values=cbPalette)

############################################
#ER plot for SHORT SOA
gr_a = ezPlot(
  data = ShSOAc_er
  , dv = .(error)
  , wid = .(id)
  , within = .(signal, cue)
  ,do_lines = TRUE
  ,split= .(cue)
  , x = .(signal)
  , x_lab = 'Signal Type'
  , y_lab = 'Reaction Time (ms)'
  ,split_lab = 'Cues'

)
gr_a+ theme_pubr(base_size = 14, base_family = "", border = FALSE,
                 margin = TRUE, legend = c("top", "bottom", "left", "right", "none"),
                 x.text.angle = 0)+
  scale_colour_manual(values=cbPalette)

#LONG SOA ERROR RATE
gr_a = ezPlot(
  data = LoSOAc_er
  , dv = .(error)
  , wid = .(id)
  , within = .(signal, cue)
  ,do_lines = TRUE
  ,split= .(cue)
  , x = .(signal)
  , x_lab = 'Signal Type'
  , y_lab = 'Reaction Time (ms)'
  ,split_lab = 'Cues'

)
gr_a+ theme_pubr(base_size = 14, base_family = "", border = FALSE,
                 margin = TRUE, legend = c("top", "bottom", "left", "right", "none"),
                 x.text.angle = 0)+
  scale_colour_manual(values=cbPalette)

############################################
#RT plot for Long SOA RT
gr_a = ezPlot(
  data = LoSOAc
  , dv = .(rt)
  , wid = .(id)
  , within = .(signal, cue)
  ,do_lines = TRUE
  ,split= .(cue)
  , x = .(signal)
  , x_lab = 'Signal Type'
  , y_lab = 'Reaction Time (ms)'
  ,split_lab = 'Cues'

)
gr_a+ theme_pubr(base_size = 14, base_family = "", border = FALSE,
                 margin = TRUE, legend = c("top", "bottom", "left", "right", "none"),
                 x.text.angle = 0)+
  scale_colour_manual(values=cbPalette)



###
#########BONUS CODE Analyses##############
###

#OTHER Error Rate#
unrestrictM_target <- glmer( error_t ~target + (1|id), data = ShSOAc_er, REML = F, family = 'binomial' ) #influence of cue
unrestrictM_cue <- glmer( error_t ~cue + (1|id), data = ShSOAc_er, REML = F, family = 'binomial' ) #influence of cue
unrestrictM_signal <- glmer( error_t ~ signal + (1|id), data = ShSOAc_er, REML = F, family = 'binomial' ) #influence of signal
unrestrictM_int <- glmer( error_t ~ signal+ cue+ signal * cue + (1|id), data = ShSOAc_er, REML = F, family = 'binomial' ) #full
restrictM <- glmer( error_t ~ (1|id), data = ShSOAc_er, REML = F, family = 'binomial' ) #only participant error (for comparison to other ME models)
restrictM_interaction <- glmer( error_t ~ (1|id) +cue + signal, data = ShSOAc_er, REML = F, family = 'binomial' )#everything but interaction (for comparison to interaction)
model_interact<-glmer( error_t ~ (1|id) +cue * signal, data = ShSOAc_er, REML = F, family = 'binomial' )
(AIC(restrictM)-AIC(unrestrictM_target))*log2(exp(1))
(AIC(restrictM)-AIC(unrestrictM_cue))*log2(exp(1))
(AIC(restrictM)-AIC(unrestrictM_signal))*log2(exp(1))
(AIC(restrictM_interaction)-AIC(unrestrictM_int))*log2(exp(1))
summary(c_error

