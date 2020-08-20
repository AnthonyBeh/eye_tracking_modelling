# yet another blooming thing to learn.
#
# R version of psychometric fitting (probably long term good idea for plotting this kind of thing)


# installing package that reads HDF5 --------------------------------------

install.packages('BiocManager')
BiocManager::install("rhdf5")

# installing package that deals with eye tracking data --------------------

install.packages('devtools')
remotes::install_github("tmalsburg/saccades/saccades", dependencies=TRUE)


install.packages('imager')
install.packages("dplyr")
install.packages("tidyverse")
install.packages("quickpsy")
install.packages('janitor')
install.packages('ggplot2')
install.packages('rstatix')
install.packages('hdf5r')
library(purrr)  # functional programming
library(imager)
library("devtools")
library(dplyr)
library(tidyverse)
library(quickpsy)
library(janitor)
library(ggplot2)
library(rhdf5)
library(saccades)
library(rstatix)
library(hdf5r)

# Analyzing single read.csv files ------------------------------------------------

# # save all read.csv files and read them into all_data
# files <- gsub("\\.csv$","",list.files(path = "./000", pattern = "*read.csv"))  # This takes the file names without the extension (?gsub)   
# all_data = quickreadfiles(path = './000',
#                           session = files, extension ='csv')
# 
# # Clean up names and add a few more variables
# all_data <- all_data %>%
#   janitor::clean_names() %>% 
#   mutate(inverse_reading_speed = 1/reading_speed,
#          blurring = 1 - cutoff_value) %>% 
#   mutate(percent_readingSpeed = inverse_reading_speed/max(inverse_reading_speed))
# 

# Now do it for all subjects ----------------------------------------------
setwd('E://thumbdrive/ageMatchControls/')

# Get all subject folders
subDir <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
pullRead <- function(fname){
  files <- gsub("\\.csv$","",list.files(path = fname, pattern = "*read.csv"))
  all_data = quickreadfiles(path = fname,
                            session = files, extension ='csv')
  all_data <- all_data %>% 
    janitor::clean_names() %>% 
    # filter(cutoff_value <1) %>%
    mutate(iSub = as.numeric(str_remove(fname,'./'))) %>% 
    mutate(inverse_reading_speed = 1/reading_speed,
           blurring = 1 - cutoff_value) %>%

    # get the baseline as a new column
    mutate(baseline=case_when( cutoff_value == 1.00 ~ inverse_reading_speed,
                               TRUE ~ 0)) %>% 
    # change all values to the baseline
    mutate(baseline = max(baseline)) %>% 
    
    # get percent change as a function of baseline
    mutate(percent_readingSpeed = inverse_reading_speed/baseline) %>% 
    rename(iTrial = cutoff_value) %>% 
    rename(iVfl = visual_field_loss_side)
  
}
readData <- map_df(subDir, pullRead)



## Now the current dataset
# Read specific values
mono_eyeSamples ="/data_collection/events/eyetracker/MonocularEyeSampleEvent"
message_events = "/data_collection/events/experiment/MessageEvent"

# get all hdf5 files
filenames <- list.files(path = , pattern = "*.hdf5",full.names = TRUE, recursive = TRUE)

# Pulls event messages
pullEvents <- function(fname){
  h5read(file = fname, name = message_events) %>% 
    # filter for the latest session, as sometimes it saves previous settings
    filter(session_id == max(session_id))
}

# Pulls eye tracking traces
pullTraces <- function(fname){
  h5read(file = fname, name = mono_eyeSamples) %>% 
    # filter for the latest session, as sometimes it saves previous settings
    filter(session_id == max(session_id)) %>% 
    select(time,gaze_x,gaze_y,pupil_measure1)
}

# Finds stimulus start and end time using pullEvents
findStimTime <- function(fname){
  events = pullEvents(fname)
  
  # Find when the stimulus starts and ends
  stimStart = which(grepl("Stimulus Start", events$text))
  stimStart = events$time[tail(stimStart,n=1)+1]
  
  cutOff = which(grepl("Cutoff level", events$text))
  cutOff = events$text[cutOff+1]
  cutOff= cutOff[1] #hack because sometimes experiment crashed and reinput params, causing multiple versions\
  cutOff = as.numeric(cutOff)
  checkVFLSide = which(grepl("Heminopia Side", events$text))
  checkVFLSide = as.numeric(events$text[checkVFLSide+1])
  if (checkVFLSide ==1){
    vflSide = "R"
  } else{
    vflSide = "L"
  }
  
  stimEnd = which(grepl("Stimulus End", events$text))
  stimEnd = events$time[stimEnd]
  
  df <- data.frame(cutOff,stimStart, stimEnd,vflSide)
  colnames(df) <- c("cutoff","stimulus_start","stimulus_end","vflSide")
  return(df)
}

# Filter traces from pullTraces using events from findStimTime
filterTraces <- function(fname){
  
  # Get subject ID
  subID = (gsub(".hdf5Files*(.*)","",fname))
  subID = as.numeric(sub("./","",subID))
  
  dfFilter=findStimTime(fname)
  dfData = pullTraces(fname)
  dfData = subset(dfData, dfData$time > dfFilter$stimulus_start[1] &
                    dfData$time < dfFilter$stimulus_end[1] )
  dfData = dfData %>%
    mutate(cutoff = dfFilter$cutoff) %>%
    mutate(vflSide = dfFilter$vflSide) %>% 
    mutate(subject = subID) %>% 
    #mutate(time = time-min(time)) %>%
    
    # filter out samples that are delayed in time (comment out to check how many samples have been filtered out)
    mutate(diff = c(0.1,diff(time))) %>% 
    filter(diff>0) %>% 

    rename(x = gaze_x, y = gaze_y, trial = cutoff) %>%
    select(subject, vflSide, time, x, y, trial,diff) 
}

# Get all eye traces for each cut off
all_eye <- map_df(filenames, filterTraces)

library(signal)
## Helper function for adaptive algorithm
# this function runs an interative process to find the estimated velocity threshold based on the data
vec_Estimate <- function(x) {
  prevThres <-  0
  newThres <- 100 # give a number higher than the current threshold
  while (abs(newThres - prevThres)>1){
    prevThres <- newThres # assign the new threshold onto the previous threshold
    # get the mean and standard deviation and calculate the new threshold
    filtData <- x[x<prevThres]
    meanVec <- mean(filtData)
    stdVec <- sd(filtData)
    newThres <- meanVec + 6*stdVec
    # print(paste("the current threshold is ",newThres))
    # print(paste("Difference from previous threshold ", abs(newThres-prevThres)))
  }
  # print(paste("The estimated threshold is ", newThres))
  est_params <- c(newThres, meanVec, stdVec)
  return(est_params)
}


## Detection paramters
# Eye tracker sample rate
sampleRate <- 500

# Saccade parameters
sac_vec_thres <- 100 # Initial saccade velocity threshold (degrees)
sac_min_dur <- 0.01*sampleRate # Filter out saccades that are too small (10msec)
# Filter parameter (Satvitzky Golay) ------------------------------------
filter_order = 2
filter_window = 11
# Fixation paramters
fix_min_dur <- 0.04*sampleRate # Minimum fixation duration
# Glissade parameters (refer Nystrom 2010)
alpha <-  0.7 # global noise weighting
beta <-  0.3 # local noise weighting


eye_events <- data.frame(matrix(ncol = 15, nrow = 0))
colnames(eye_events) <- c("subject","vflSide","trial", 
                          "num_fix", "mean_fixDur","sd_fixDur",
                          "num_sac", "mean_sacDur","sd_sacDur",
                          "num_highGlis", "mean_highGlisDur","sd_highGlisDur",
                          "num_lowGlis","mean_lowGlisDur","sd_lowGlisDur")

count= 0
for (iSub in unique(all_eye$subject)){
  for (iVfl in unique(all_eye$vflSide)){
    for (iTrial in unique(all_eye$trial)){
      if (iTrial == 1){
        tempDf <- all_eye %>% 
          dplyr::filter(subject == iSub) %>% 
          dplyr::filter(trial == iTrial) 
      } else {
        tempDf <- all_eye %>% 
          dplyr::filter(subject == iSub) %>% 
          dplyr::filter(vflSide == iVfl) %>%
          dplyr::filter(trial == iTrial)
      }
      count = count +1
      print(paste("run no:", count))
      tempDf <- tempDf %>% 
        dplyr::filter(x !=0 | y !=0) %>% # remove all samples where x and y is 0, maybe look to see if there is also pupil data to make this better
        mutate(time_sec = time-tempDf$time[1]) %>% 
        mutate(filt_x = sgolayfilt(x, p=filter_order , n =filter_window)) %>% 
        mutate(filt_y = sgolayfilt(y, p=filter_order , n =filter_window)) %>% 
        mutate(filt_D = sqrt(filt_x^2+filt_y^2)) %>% 
        mutate(filt_vec = c(NA,abs((diff(filt_D)/0.002))*0.01972)) %>% 
        mutate(filt_acc = c(NA,(diff(filt_vec)/0.002))) %>% 
        dplyr::filter(filt_vec < 1000) %>%  # filter for velocity above physiological limitation
        dplyr::filter(filt_acc < 10000) %>%  # filter for acceleration above physiological limitation
        select(subject, vflSide, trial, time_sec, x, y, filt_x, filt_y, filt_D, filt_vec, filt_acc)
      
      # find the adapted velocity threshold
      # var[1] = velocity threshold, var[2] =  mean, var[3] = std  
      est_params <- vec_Estimate(tempDf$filt_vec) 
      
      # Finding saccade onsets/offsets
      # find boolean vector for samples above threshold
      upThres <- tempDf$filt_vec > est_params[1]
      diffOnOff <- c(1,diff(upThres))
      thresOn <- which(diffOnOff %in% 1)
      thresOff <- which(diffOnOff %in% -1)-1
      if (length(thresOn) > length(thresOff)){
        thresOn <- thresOn[2:length(thresOn)] # remove first dummy number
      }
      # sometimes the experiment cuts off before offset
      if (length(thresOn) > length(thresOff)){
        thresOn <- thresOn[1:length(thresOn)-1]
      }
      
      # get rid of intervals that are too short and therefore unlikely to be fixations
      on_off_interval <- thresOn[2:length(thresOn)]-thresOff[1:length(thresOff)-1] # shift "on" forward, shift "off" backward, can work out mathematically
      longEnough <- c(TRUE, (on_off_interval >=fix_min_dur )) # added dummy variable in front to compensate for shift
      thresOff <- thresOff[longEnough]
      thresOn <- thresOn[longEnough]
      
      # now let's find the saccade onset and offset
      sac_onset_thres <- est_params[2] + 3*est_params[3] # mean + 3*std
      
      onsets <- as.numeric() # create empty data frame to store onsets
      offsets <- as.numeric() # create empty data frame to store offsets
      high_glissOn <- as.numeric()
      high_glissOff <- as.numeric()
      low_glissOn <- as.numeric()
      low_glissOff <- as.numeric()
      
      # loop for number of windows
      for (iRange in 1:length(thresOn)){
        # set iOn, iOff and local noise window
        iOn <- thresOn[iRange]
        iOff <- thresOff[iRange]
        # sometimes the saccades starts at the start of the run, so there is no way to get local noise
        if (iOn<fix_min_dur){
          next
        } else{
          local_offset_window <- c((iOn-fix_min_dur):iOn)
        }
        
        local_offset_window <- tempDf$filt_vec[local_offset_window]
        # onsets ------------------------------------------------------------------
        currOn <-  tempDf$filt_vec[iOn]
        nextOn <-  tempDf$filt_vec[iOn-1]
        # looking for the first sample below saccade onset threshold and (Xn - Xn+1 >=0)
        while (currOn > sac_onset_thres & currOn - nextOn >= 0) {
          iOn <-  iOn - 1
          currOn <- nextOn
          nextOn <-  tempDf$filt_vec[iOn-1]
        }
        # offsets -----------------------------------------------------------------
        currOff <-  tempDf$filt_vec[iOff]
        nextOff <-  tempDf$filt_vec[iOff+1]
        # weighted global noise + weighted local noise before saccade
        localNoise <- mean(local_offset_window) + 3 *sd(local_offset_window)
        sac_offset_thres <-  alpha * sac_onset_thres + beta * localNoise
        # looking for the first sample below saccade offset threshold and (Xn - Xn+1 >=0)
        while (currOff > sac_offset_thres & currOff -nextOff >= 0) {
          iOff <- iOff +1
          currOff <-  nextOff 
          nextOff <-  tempDf$filt_vec[iOff+1]
        }
        
        # Filter out saccades that are too small (10msec)
        if (iOff - iOn>=sac_min_dur ){
          # append list
          onsets <- c(onsets, iOn) 
          offsets <- c(offsets, iOff)
          
          # look for glissades too
          # window of minimum fixation duration after saccadic offset
          # first let's add a break clause if this happens towards the end of the run
          if (iOff +fix_min_dur-1 > length(tempDf$filt_vec)){
            next
          }
          glissadeWindow <- tempDf$filt_vec[iOff:(iOff+fix_min_dur-1)] # -1 to get the right number
          # check for high velocity glissades then low velocity glissades
          if (any(glissadeWindow > est_params[1])){
            # onset of glissade == offset of preceeding saccade
            high_glissOn <- c(high_glissOn, iOff)
            diffOnOff <- c(1,diff(tempDf$filt_vec[iOff:(iOff+100)] > est_params[1]))
            findOff <- (which(diffOnOff %in% -1)-1)[1] # take the first offset
            # if (is.na(findOff)){
            #   # break if it's NA, means that its not a glissade
            #   next #over here
            # }
            # find the next local velocity minimum
            while (tempDf$filt_vec[iOff+findOff-1]-tempDf$filt_vec[iOff+findOff]>=0){
              findOff <- findOff + 1
            }
            high_glissOff <- c(high_glissOff, iOff+findOff-1) # -1 to get the right number
          } else if (any(glissadeWindow >sac_offset_thres)){
            # onset of glissade == offset of preceeding saccade
            low_glissOn <- c(low_glissOn, iOff)
            diffOnOff <- c(1,diff(tempDf$filt_vec[iOff:(iOff+100)] > sac_offset_thres))
            findOff <- (which(diffOnOff %in% -1)-1)[1] # take the first offset
            # if (is.na(findOff)){
            #   # break if it's NA, means that its not a glissade
            #   next #over here
            # }
            # find the next local velocity minimum
            while (tempDf$filt_vec[iOff+findOff-1]-tempDf$filt_vec[iOff+findOff]>=0){
              findOff <- findOff + 1
            }
            low_glissOff <- c(low_glissOff, iOff+findOff-1) # -1 to get the right number
          }
        }
      }
      
      # fixation
      fixations <-as.numeric()
      fixCount <-  0
      
      for (iFix in 1:length(tempDf$filt_vec)){
        fixCount <- fixCount + 1
        # if fixCount matches a number in onsets
        if (fixCount %in% onsets){
          # find index to match to offset 
          idx <- match(fixCount, onsets)
          # check to see there is glissade, if not use saccade offset
          if (offsets[idx] %in% high_glissOn){
            idx <- match(offsets[idx], high_glissOn)
            fixCount <- high_glissOff[idx]
          } else if (offsets[idx] %in% low_glissOn) {
            idx <- match(offsets[idx], low_glissOn)
            fixCount <- low_glissOff[idx]
          }else{
            fixCount <- offsets[idx]
          }
        } else{
          
          fixations <- c(fixations, fixCount)
        }
      }
      
      # filter out fixations that are too short, amazing function below I found online, using a mix of
      # diff, not 1, and cumsum to segregate bins
      fixations <- split(fixations, cumsum(c(1, diff(fixations) != 1)))
      fixations <- fixations[lengths(fixations) >= fix_min_dur]
      num_fix <- length(fixations)# number of fixations
      mean_fixDur <- mean(lengths(fixations))/sampleRate # mean fixation duration
      sd_fixDur <- sd(lengths(fixations))/sampleRate # standard deviation fixation duration
      
      # save saccades 
      sacOn <- tempDf$time_sec[onsets]
      sacOff <- tempDf$time_sec[offsets]
      num_sac <- length(sacOn) # number of saccades
      mean_sacDur <- mean(sacOff-sacOn) # mean saccade duration
      sd_sacDur <-sd(sacOff-sacOn) # standard deviation saccade duration
      
      # save glissades
      high_glissOn <- tempDf$time_sec[high_glissOn]
      high_glissOff <- tempDf$time_sec[high_glissOff]
      # remove incomplete glissades, offset is not recorded*
      if (length(high_glissOff) < length(high_glissOn)){
        high_glissOn <- high_glissOn[1:(length(high_glissOn)-1)]
      }
      num_highGlis <- length(high_glissOn) # number of glissades
      mean_highGlisDur <- mean(high_glissOff-high_glissOn) # mean glissade duration
      sd_highGlisDur <- sd(high_glissOff-high_glissOn) # standard deviation glissade duration
      
      low_glissOn <- tempDf$time_sec[low_glissOn]
      low_glissOff <- tempDf$time_sec[low_glissOff]
      # remove incomplete glissades, offset is not recorded*
      if (length(low_glissOff) < length(low_glissOn)){
        low_glissOn <- low_glissOn[1:(length(low_glissOn)-1)] 
      }
      num_lowGlis <- length(low_glissOn) # number of glissades
      mean_lowGlisDur <- mean(low_glissOff-low_glissOn) # mean glissade duration
      sd_lowGlisDur <- sd(low_glissOff-low_glissOn) # standard deviation glissade duration
      
      df <- data.frame(iSub, iVfl, iTrial,
                       num_fix, mean_fixDur,sd_fixDur,
                       num_sac, mean_sacDur,sd_sacDur,
                       num_highGlis, mean_highGlisDur,sd_highGlisDur,
                       num_lowGlis,mean_lowGlisDur,sd_lowGlisDur)
      
      eye_events <- rbind(eye_events, df)
      
      
      # return(eye_events)
      
    }
  }
}


finalData <- merge(eye_events, readData, by=c('iSub','iTrial','iVfl'))
finalData <- finalData %>% 
  select(iSub, iVfl, iTrial, reading_speed, num_fix, mean_fixDur, sd_fixDur,
         num_sac, mean_sacDur, sd_sacDur, num_highGlis, mean_highGlisDur, sd_highGlisDur,
         num_lowGlis, mean_lowGlisDur, sd_lowGlisDur) %>% 
  mutate(iVfl = case_when(iTrial == 1 ~ 'R', # if trial is 1 (baseline), set the vfl to 'R', otherwise, leave it (TRUE)
                          TRUE ~ as.character(iVfl)))




# reading speed -----------------------------------------------------------

reading_behaviour <- finalData %>% 
  mutate(baseline = case_when(iTrial ==  1 ~ reading_speed,
                              TRUE ~ 0)) %>% 
  group_by(iSub) %>% 
  mutate(baseline = max(baseline)) %>% 
  mutate(percent_readSpeed = baseline/reading_speed) %>% 
  select(iTrial, iVfl, percent_readSpeed)

# set up plot variables
plot <- reading_behaviour %>%
  group_by(iTrial,iVfl) %>%
  dplyr::filter(iTrial <1) %>% 
  dplyr::summarize(count=n(),mean = mean(percent_readSpeed, na.rm=TRUE),
                   sd = sd(percent_readSpeed, na.rm=TRUE)) 
  
                   
ggplot(plot, aes(x = factor(iTrial), y = mean, fill = iVfl))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9))+
  geom_hline(yintercept = 1, linetype='dashed')+
  theme_minimal()+
  theme(text=element_text(size =12))+
  labs(title = "Reading speed loss across different cut-off frequencies",
         x = 'relative cut-off frequency', y = 'average percent reading speed (%)')+
  scale_fill_discrete(name = "Simulated VFL", 
                       labels=c("Left", "Right"))
  
# running aov
reading_behaviour <- reading_behaviour %>% 
  dplyr::filter(iTrial <1)
# readSpeed.aov <- aov(percent_readSpeed ~ factor(iTrial) * iVfl, data = reading_behaviour)
# summary(readSpeed.aov)
# model.tables(readSpeed.aov, type="means", se = TRUE)
# TukeyHSD(readSpeed.aov, "factor(iTrial)")




linearRead <- lm(scale(percent_readSpeed) ~ scale(iTrial) + iVfl + iTrial:iVfl, data=reading_behaviour) 
summary(linearRead)



# Oculometric measures ---------------------------------------------------


# number of fixations -----------------------------------------------------
analyzeDf <- finalData %>% 
  mutate(baseline = case_when(iTrial ==  1 ~ num_fix,
                              TRUE ~ as.integer(0))) %>% 
  group_by(iSub) %>% 
  mutate(baseline = max(baseline)) %>% 
  select(iTrial, iVfl, num_fix,baseline)
# set up plot variables
plot <- analyzeDf %>%
  group_by(iTrial,iVfl) %>%
  dplyr::filter(iTrial <1) %>% 
  dplyr::summarize(count=n(),mean = mean(num_fix, na.rm=TRUE),
                   sd = sd(num_fix, na.rm=TRUE)) 
ggplot(plot, aes(x = factor(iTrial), y = mean, fill = iVfl))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9))+
  geom_hline(yintercept = analyzeDf$baseline[1], linetype='dashed')+
  theme_minimal()+
  theme(text=element_text(size =12))+
  labs(title = "Average number of fixations \n across different cut-off frequencies",
       x = 'relative cut-off frequency', y = 'average number of fixations')+
  scale_fill_discrete(name = "Simulated VFL", 
                      labels=c("Left", "Right"))


analyzeDf <- analyzeDf %>% 
  dplyr::filter(iTrial <1)
linearRead <- lm(num_fix ~ iTrial + iVfl + iTrial:iVfl, data=analyzeDf) 
summary(linearRead)



# mean fixation duration -------------------------------------------------------

analyzeDf <- finalData %>% 
  mutate(baseline = case_when(iTrial ==  1 ~ mean_fixDur,
                              TRUE ~ 0)) %>% 
  group_by(iSub) %>% 
  mutate(baseline = max(baseline)) %>% 
  select(iTrial, iVfl, mean_fixDur,baseline)
# set up plot variables
plot <- analyzeDf %>%
  group_by(iTrial,iVfl) %>%
  dplyr::filter(iTrial <1) %>% 
  dplyr::summarize(count=n(),mean = mean(mean_fixDur, na.rm=TRUE),
                   sd = sd(mean_fixDur, na.rm=TRUE)) 
ggplot(plot, aes(x = factor(iTrial), y = mean, fill = iVfl))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9))+
  geom_hline(yintercept = analyzeDf$baseline[1], linetype='dashed')+
  theme_minimal()+
  theme(text=element_text(size =12))+
  labs(title = "Average mean fixation duration \n across different cut-off frequencies",
       x = 'relative cut-off frequency', y = 'average mean fixation duration (s)')+
  scale_fill_discrete(name = "Simulated VFL", 
                      labels=c("Left", "Right"))


analyzeDf <- analyzeDf %>% 
  dplyr::filter(iTrial <1)
linearRead <- lm(mean_fixDur ~ iTrial + iVfl + iTrial:iVfl, data=analyzeDf) 
summary(linearRead)


# number of saccades  -------------------------------------------------------

analyzeDf <- finalData %>% 
  mutate(baseline = case_when(iTrial ==  1 ~ num_sac,
                              TRUE ~ as.integer(0))) %>% 
  group_by(iSub) %>% 
  mutate(baseline = max(baseline)) %>% 
  select(iTrial, iVfl, num_sac,baseline)
# set up plot variables
plot <- analyzeDf %>%
  group_by(iTrial,iVfl) %>%
  dplyr::filter(iTrial <1) %>% 
  dplyr::summarize(count=n(),mean = mean(num_sac, na.rm=TRUE),
                   sd = sd(num_sac, na.rm=TRUE)) 
ggplot(plot, aes(x = factor(iTrial), y = mean, fill = iVfl))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9))+
  geom_hline(yintercept = analyzeDf$baseline[1], linetype='dashed')+
  theme_minimal()+
  theme(text=element_text(size =12))+
  labs(title = "Average number of saccades \n across different cut-off frequencies",
       x = 'relative cut-off frequency', y = 'average number of saccades')+
  scale_fill_discrete(name = "Simulated VFL", 
                      labels=c("Left", "Right"))


analyzeDf <- analyzeDf %>% 
  dplyr::filter(iTrial <1)
linearRead <- lm(num_sac ~ iTrial + iVfl + iTrial:iVfl, data=analyzeDf) 
summary(linearRead)

# mean saccade duration -------------------------------------------------------

analyzeDf <- finalData %>% 
  mutate(baseline = case_when(iTrial ==  1 ~ mean_sacDur,
                              TRUE ~ 0)) %>% 
  group_by(iSub) %>% 
  mutate(baseline = max(baseline)) %>% 
  select(iTrial, iVfl, mean_sacDur,baseline)
# set up plot variables
plot <- analyzeDf %>%
  group_by(iTrial,iVfl) %>%
  dplyr::filter(iTrial <1) %>% 
  dplyr::summarize(count=n(),mean = mean(mean_sacDur, na.rm=TRUE),
                   sd = sd(mean_sacDur, na.rm=TRUE)) 
ggplot(plot, aes(x = factor(iTrial), y = mean, fill = iVfl))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9))+
  geom_hline(yintercept = analyzeDf$baseline[1], linetype='dashed')+
  theme_minimal()+
  theme(text=element_text(size =12))+
  labs(title = "Average mean saccade duration \n across different cut-off frequencies",
       x = 'relative cut-off frequency', y = 'average mean saccade duration (s)')+
  scale_fill_discrete(name = "Simulated VFL", 
                      labels=c("Left", "Right"))


analyzeDf <- analyzeDf %>% 
  dplyr::filter(iTrial <1)
linearRead <- lm(mean_sacDur ~ iTrial + iVfl + iTrial:iVfl, data=analyzeDf) 
summary(linearRead)
# number of high velocity glissades  -------------------------------------------------------

analyzeDf <- finalData %>% 
  mutate(baseline = case_when(iTrial ==  1 ~ num_highGlis,
                              TRUE ~ as.integer(0))) %>% 
  group_by(iSub) %>% 
  mutate(baseline = max(baseline)) %>% 
  select(iTrial, iVfl, num_highGlis,baseline)
# set up plot variables
plot <- analyzeDf %>%
  group_by(iTrial,iVfl) %>%
  dplyr::filter(iTrial <1) %>% 
  dplyr::summarize(count=n(),mean = mean(num_highGlis, na.rm=TRUE),
                   sd = sd(num_highGlis, na.rm=TRUE)) 
ggplot(plot, aes(x = factor(iTrial), y = mean, fill = iVfl))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9))+
  geom_hline(yintercept = analyzeDf$baseline[1], linetype='dashed')+
  theme_minimal()+
  theme(text=element_text(size =12))+
  labs(title = "Average number of high velocity glissades \n across different cut-off frequencies",
       x = 'relative cut-off frequency', y = 'average number of high velocity glissades')+
  scale_fill_discrete(name = "Simulated VFL", 
                      labels=c("Left", "Right"))


analyzeDf <- analyzeDf %>% 
  dplyr::filter(iTrial <1)
linearRead <- lm(num_highGlis ~ iTrial + iVfl + iTrial:iVfl, data=analyzeDf) 
summary(linearRead)

# mean high velocity glissade duration -------------------------------------------------------

analyzeDf <- finalData %>% 
  mutate(baseline = case_when(iTrial ==  1 ~ mean_highGlisDur,
                              TRUE ~ 0)) %>% 
  group_by(iSub) %>% 
  mutate(baseline = max(baseline)) %>% 
  select(iTrial, iVfl, mean_highGlisDur,baseline)
# set up plot variables
plot <- analyzeDf %>%
  group_by(iTrial,iVfl) %>%
  dplyr::filter(iTrial <1) %>% 
  dplyr::summarize(count=n(),mean = mean(mean_highGlisDur, na.rm=TRUE),
                   sd = sd(mean_highGlisDur, na.rm=TRUE)) 
ggplot(plot, aes(x = factor(iTrial), y = mean, fill = iVfl))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9))+
  geom_hline(yintercept = analyzeDf$baseline[1], linetype='dashed')+
  theme_minimal()+
  theme(text=element_text(size =12))+
  labs(title = "Average mean high velocity glissade duration \n across different cut-off frequencies",
       x = 'relative cut-off frequency', y = 'average mean high velocity glissade duration (s)')+
  scale_fill_discrete(name = "Simulated VFL", 
                      labels=c("Left", "Right"))


analyzeDf <- analyzeDf %>% 
  dplyr::filter(iTrial <1)
linearRead <- lm(mean_highGlisDur ~ iTrial + iVfl + iTrial:iVfl, data=analyzeDf) 
summary(linearRead)
# number of low velocity glissades  -------------------------------------------------------

analyzeDf <- finalData %>% 
  mutate(baseline = case_when(iTrial ==  1 ~ num_lowGlis,
                              TRUE ~ as.integer(0))) %>% 
  group_by(iSub) %>% 
  mutate(baseline = max(baseline)) %>% 
  select(iTrial, iVfl, num_lowGlis,baseline)
# set up plot variables
plot <- analyzeDf %>%
  group_by(iTrial,iVfl) %>%
  dplyr::filter(iTrial <1) %>% 
  dplyr::summarize(count=n(),mean = mean(num_lowGlis, na.rm=TRUE),
                   sd = sd(num_lowGlis, na.rm=TRUE)) 
ggplot(plot, aes(x = factor(iTrial), y = mean, fill = iVfl))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9))+
  geom_hline(yintercept = analyzeDf$baseline[1], linetype='dashed')+
  theme_minimal()+
  theme(text=element_text(size =12))+
  labs(title = "Average number of low velocity glissades \n across different cut-off frequencies",
       x = 'relative cut-off frequency', y = 'average number of low velocity glissades')+
  scale_fill_discrete(name = "Simulated VFL", 
                      labels=c("Left", "Right"))


analyzeDf <- analyzeDf %>% 
  dplyr::filter(iTrial <1)
linearRead <- lm(num_lowGlis ~ iTrial + iVfl + iTrial:iVfl, data=analyzeDf) 
summary(linearRead)

# mean low velocity glissade duration -------------------------------------------------------

analyzeDf <- finalData %>% 
  mutate(baseline = case_when(iTrial ==  1 ~ mean_lowGlisDur,
                              TRUE ~ 0)) %>% 
  group_by(iSub) %>% 
  mutate(baseline = max(baseline)) %>% 
  select(iTrial, iVfl, mean_lowGlisDur,baseline)
# set up plot variables
plot <- analyzeDf %>%
  group_by(iTrial,iVfl) %>%
  dplyr::filter(iTrial <1) %>% 
  dplyr::summarize(count=n(),mean = mean(mean_lowGlisDur, na.rm=TRUE),
                   sd = sd(mean_lowGlisDur, na.rm=TRUE)) 
ggplot(plot, aes(x = factor(iTrial), y = mean, fill = iVfl))+
  geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9))+
  geom_hline(yintercept = analyzeDf$baseline[1], linetype='dashed')+
  theme_minimal()+
  theme(text=element_text(size =12))+
  labs(title = "Average mean low velocity glissade duration \n across different cut-off frequencies",
       x = 'relative cut-off frequency', y = 'average mean low velocity glissade duration (s)')+
  scale_fill_discrete(name = "Simulated VFL", 
                      labels=c("Left", "Right"))


analyzeDf <- analyzeDf %>% 
  dplyr::filter(iTrial <1)
linearRead <- lm(mean_lowGlisDur ~ iTrial + iVfl + iTrial:iVfl, data=analyzeDf) 
summary(linearRead)





#
# install.packages('lsmeans')
# library(lsmeans)
# # construct linear regression model
# m.interaction <- lm(iTrial ~ percent_readSpeed*iVfl, data = reading_behaviour)
# # run anova
# anova(m.interaction)
# # check coefficients
# m.interaction$coefficients
# # obtain slopes
# m.lst <- lstrends(m.interaction, "iVfl", var="percent_readSpeed")
# m.lst
# # compare slopes
# pairs(m.lst)
# 
# 
# install.packages("psych")
# library(psych)
# library(data.table)
# readTable <- as.data.table(reading_behaviour)
# # Calculate Pearson's R
# m.correlations <- readTable[, cor(iTrial, percent_readSpeed), by = iVfl]
# m.correlations
# # Compare R values with Fisher's R to Z
# psych::paired.r(m.correlations[iVfl=="R", V1], m.correlations[iVfl=="L", V1], 
#          n = readTable[iVfl %in% c("R", "L"), .N])
# 
# install.packages("mnormt")
# 
# writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")






# Visualize data ----------------------------------------------------------

plot(tempDf$time_sec, tempDf$filt_vec, type = 'l',main = "Event detection",
     cex.main = 1.5,
     xlab = "Time (s)",
     ylab = "Angular velocity (\u00B0/sec)",
     cex.lab =1.2,
     xlim= c(0,1),
     cex = 4)
abline(h =est_params[1], col = 'red', lty = 2)
abline(v = sacOn, col = 'red')
abline(v = sacOff, col = 'blue')
abline(v = high_glissOff, col = 'blue', lty = 2)
abline(v = low_glissOff, col = 'blue', lty = 2)

legend(0.6,300, legend = c("saccadic onset", "saccadic offset", "glissade offset","velocity threshold"), col = c("red","blue","blue","red"),
       lty = c(1,1,2,2), cex = 0.8)


#ggplot(fixStats, aes(x=trial, y =totalMeanFix, colour = vflSide))+
  #   geom_bar(stat="identity",position=position_dodge(),fill="white")+
  #   geom_errorbar(aes(ymin=totalMeanFix -sdMeanFix ,
  #                     ymax=totalMeanFix + sdMeanFix ),
  #                 width=.01, size =1, position=position_dodge(.045)) +
  #   theme_minimal()+
  #   labs(x= "cutoff values", y = "average mean fixations (s)", colour = 'visual field loss side')
# Time to run some analysis -----------------------------------------------




# # Plot out full eye traces
# ggplot(single_eye, aes(x=x, y= y))+
#   geom_point(size=0.005)+
#   geom_line(size = 0.005)+
#   coord_cartesian(xlim =c(-512,512), ylim = c(-384, 384))

# # Loop through to look for fixations
# fixData = data.frame()
# for (iSub in 1:length(unique(all_eye$subject))){
#   for(iTrial in 1:length(unique(all_eye$trial))){
#     for(iVFL in 1:length(unique(all_eye$vflSide))){
#       
#       subject = unique(all_eye$subject)[iSub]
#       trial = unique(all_eye$trial)[iTrial]
#       vfl = unique(all_eye$vflSide)[iVFL]
#       df <- getFixStats(all_eye, subject, trial, vfl)
#       fixData <- rbind(fixData, df)
#       print(fixData)
# 
#     }
#   }
# }
# 
# # Debug code to look into fixation dataframe
# debugFix <-  function(df,sub,cutoff,vfl){
#   # Filter data
#   if (cutoff==1.00){
#     single_eye <- subset(df, subject == sub & trial == cutoff )
#   } else {
#     single_eye <- subset(df, subject == sub & trial == cutoff & vflSide == vfl )
#   }
#   fixations <-detect.fixations(single_eye)
# }
# single_eye <- subset(all_eye, subject == 11 & trial == 0.15 & vflSide == 'R' )
# fixations <- subset(detect.fixations(single_eye), event =='fixation')
# debug <- debugFix(all_eye, 3, 0.05, 'L')
# singleEyeRaw <- list(single_eye,fixations)
# nms <- c("raw_data","fixation_data")
# names(singleEyeRaw) <- nms
# save(singleEyeRaw,file='sub11_R_text1.Rdata')
# 
# # Visualise fixation data
# filterFix <- fixData %>% 
#   filter(trial<1)
# baseline <- fixData %>% 
#   filter(trial==1)
# 
# # mean fixation duration
# ggplot(filterFix, aes(x = trial, y=mean_fixation, colour = vflSide))+
#   geom_bar(stat="identity",position=position_dodge(),fill="white")+
#   geom_hline(data = baseline, aes(yintercept=mean_fixation),linetype='dashed')+
#   theme_minimal()+
#   facet_wrap(~subject)+
#   labs(x= "cutoff values", y = "mean fixation duration (s)", colour = 'visual field loss side') 
# 
# # number of fixations
# ggplot(filterFix, aes(x = trial, y=number_fixation, colour = vflSide))+
#   geom_bar(stat="identity",position=position_dodge(),fill="white")+
#   geom_hline(data = baseline, aes(yintercept=number_fixation),linetype='dashed')+
#   theme_minimal()+
#   facet_wrap(~subject)+
#   labs(x= "cutoff values", y = "fixation duration (s)", colour = 'visual field loss side') 
# 
# 
# # Group by the two independent variables
# fixStats <- filterFix %>% 
#   group_by(trial, vflSide) %>% 
#   summarise(totalMeanFix = mean(mean_fixation), sdMeanFix = sd(mean_fixation))
# 
# # Visualise data
# ggplot(fixStats, aes(x=trial, y =totalMeanFix, colour = vflSide))+
#   geom_bar(stat="identity",position=position_dodge(),fill="white")+
#   geom_errorbar(aes(ymin=totalMeanFix -sdMeanFix ,
#                     ymax=totalMeanFix + sdMeanFix ),
#                 width=.01, size =1, position=position_dodge(.045)) +
#   theme_minimal()+
#   labs(x= "cutoff values", y = "average mean fixations (s)", colour = 'visual field loss side')
# 
# 
# 
# # Merging datasets together (tidy data) -----------------------------------
# 
# test <- merge(readData, fixData, by = c('subject','trial','vflSide'))
# test <- test %>% 
#   select(subject,trial,vflSide,text_used,reading_speed,inverse_reading_speed,
#          baseline,percent_readingSpeed, number_fixation,mean_fixation,std_fixation)
# 
# UG_data <- test %>% 
#   select(subject,trial,vflSide,text_used,reading_speed,percent_readingSpeed, number_fixation,mean_fixation,std_fixation)
# 
# write.csv(UG_data, "UG.csv")
# # Maybe do a full check for when samples go back in time? (figure it out tomorrow)
# check <- getFixStats(all_eye,0,0.2,'L')
# 
# check2 <- pullTraces(filenames[1])
# check3 <- filterTraces(filenames[1])
# check4 <-  h5read(file = filenames[10], name = mono_eyeSamples)
# check5 <- pullEvents(filenames[1])
# check6 <- findStimTime(filenames[1])
# check7 <- all_eye %>% 
#   filter(subject == 0, trial == 0.2, vflSide == 'L')
# test <- data.frame(diff(as.matrix(check$time)))
# which(test<=0)
# test2 <- data.frame(test[which(test>=0),1])
# test3 <- check %>% 
#   # filter(diff(time)>0)
#   mutate(diff2 = c(0.1,diff(time))) %>% 
#   filter(diff2<0)
# 
# 
# 
# fixations <- detect.fixations(check)
# plot(check$time)
# # Draw some figures to illustrate fixations
# ggplot(fixations, aes(x=x, y= y))+
#   geom_point(size=1)+
#   geom_line(size = 1)+
#   coord_cartesian(xlim =c(-512,512), ylim = c(-384, 384))
# 
# par(mfrow=c(1,2))
# 
# plot(single_eye$time, single_eye$x, cex = .05, col = 'red', xlab = 'time (s)', ylab = 'gaze x position')
# lines(fixations$start,fixations$x, lty = 2,lwd = 3)
# 
# plot(single_eye$time, single_eye$y,cex = .05, col = 'red', xlab = 'time(s)', ylab = 'gaze y position')
# lines(fixations$start,fixations$y, lty = 2, lwd = 3)
# 
# 
# 
# 
# 
# 
# diagnostic.plot(sub_eye$`0`, fixations)
# ggplot(single_eye, aes(x, y)) +
#   geom_point(size=0.2) +
#   geom_line(size=0.2)+
#   coord_fixed() +
#   facet_wrap(~trial)


















# lets work with the first 1000 samples

p1 <- ggplot(tempDf,aes(x,y))+
  geom_point(size = 1)+
  # geom_line(lty= 1,alpha= 0.7,colour ='red')+
  coord_cartesian(xlim = c(-512,512),ylim=c(-384,384))+
  theme_minimal()

p2 <- ggplot(tempDf[first(fixTest):length(fixTest),],aes(x[fixTest]+512,y[fixTest]+384))+
  geom_point(size = 1)+
  # geom_line(lty= 1,alpha= 0.7,colour ='red')+
  coord_cartesian(xlim = c(0,1024),ylim=c(0,768))+
  theme_minimal()


multiplot(p1, p2, cols = 2)

# process image and read text!
install.packages("tesseract")
install.packages("magick")
library(tesseract)
library(magick)
img <- image_read("E:\\thumbdrive\\gazeFilter\\textImg\\4.png")
cat(image_ocr(img))



# creating AOI in reading text
install.packages("imager")
library(imager)
img <- load.image("E:\\thumbdrive\\gazeFilter\\textImg\\4.png")
img <- grayscale(img) 
img <- threshold(img,"30%") # threshold to get right values
img_px <- as.pixset(img)
plot(img_px)
highlight(img_px)
# grow some pixels
img_blob <- grow(img_px,50) %>% 
  shrink(3) %>% 
  grow(4) %>% 
  shrink(4) %>% 
  plot(main = "grow then shrink")

  # get bounding box
box <- plot(img_blob) %>% highlight

test <- mask(x = box, mask = denpasar)


install.packages("BiocManager")
BiocManager::install("EBImage")
library("EBImage")

img = readImage("E:\\thumbdrive\\gazeFilter\\textImg\\4.png")
display(img, method="browser")
plot(img)
test = bwlabel(img)
table(test)
max(test)
display( normalize(test) )


length(tempDf$x[fixTest])
length(tempDf[first(fixTest):last(fixTest),])
fixations[2]
fixTest <- unlist(fixations,use.names=FALSE)
fixTest

fixations <- split(fixations, cumsum(c(1, diff(fixations) != 1)))
fixations <- fixations[lengths(fixations) >= fix_min_dur]
num_fix <- length(fixations)# number of fixations
mean_fixDur <- mean(lengths(fixations))/sampleRate # mean fixation duration
sd_fixDur <- sd(lengths(fixations))/sampleRate # standard deviation fixation duration

# save saccades 
sacOn <- tempDf$time_sec[onsets]
sacOff <- tempDf$time_sec[offsets]
num_sac <- length(sacOn) # number of saccades
mean_sacDur <- mean(sacOff-sacOn) # mean saccade duration
sd_sacDur <-sd(sacOff-sacOn) # standard deviation saccade duration

# save glissades
high_glissOn <- tempDf$time_sec[high_glissOn]
high_glissOff <- tempDf$time_sec[high_glissOff]











# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}




