## This script processes data, performs statistical tests, and generates figures
## for the revised submission of:
##    Title: Motor learning without movement
##    Authors: Kim, Forrence, McDougle

### USER TODOS: ----------------------------------------------------------------
homedir <- 'C:/Users/kimol/Documents/GitHub/LearningFromThePathNotTaken' ## TODO: Set home directory -- change this to the directory you
                                                                        ## saved the data in
savedir <- paste(homedir,'/output',sep='')
setwd(homedir) 
excludeBasedOnInstrRecall <- 1                                       ## TODO: Set excludeBasedOnInstrRecall to 1 to conduct data analyses
                                                                          ## without data from participants who failed to completely recall
                                                                          ## the task instructions after online experiments. Set to any other
                                                                          ## value to include data from all participants in analyses.


# import required packages -----------------------------------------------------
library(magrittr)
library(rstatix)
library(effectsize)
library(lmerTest)
library(ggplot2)
library(ggeffects)
library(MuMIn)
library(emmeans)
library(effsize)
library(ggpubr)
library(r2glmm)



# set up some custom functions -------------------------------------------------
#------------------ for data organization & repeated statistical tests ---------

splitRotnMvtMeans <- function(data) {
  stop_cw = tapply(data$Single_triplet_learning[
    data$Rotation == 'cw' & data$Movement == 'stop'],
    data$Subject_ID[data$Rotation == 'cw' & data$Movement == 'stop'],
    mean, na.rm = T)
  stop_ccw = tapply(data$Single_triplet_learning[
    data$Rotation == 'ccw' & data$Movement == 'stop'],
    data$Subject_ID[data$Rotation == 'ccw' & data$Movement == 'stop'],
    mean, na.rm = T)
  stop_0 = tapply(data$Single_triplet_learning[
    data$Rotation == '0' & data$Movement == 'stop'],
    data$Subject_ID[data$Rotation == '0' & data$Movement == 'stop'],
    mean, na.rm = T)
  go_0 = tapply(data$Single_triplet_learning[
    data$Rotation == '0' & data$Movement == 'go'],
    data$Subject_ID[data$Rotation == '0' & data$Movement == 'go'],
    mean, na.rm = T)
  go_cw = tapply(data$Single_triplet_learning[
    data$Rotation == 'cw' & data$Movement == 'go'],
    data$Subject_ID[data$Rotation == 'cw' & data$Movement == 'go'],
    mean, na.rm = T)
  go_ccw = tapply(data$Single_triplet_learning[
    data$Rotation == 'ccw' & data$Movement == 'go'],
    data$Subject_ID[data$Rotation == 'ccw' & data$Movement == 'go'],
    mean, na.rm = T)
  output = data.frame(go_cw, go_ccw, go_0, stop_cw, stop_ccw, stop_0)
  return(output)
}

meanDfToLongForm <- function(meandf) {
  meanvals = c()
  mvt = c()
  rotn = c()
  subj = c()
  for (m in 1:length(meandf$stop_cw)) {
    subj = c(subj, m)
    mvt = c(mvt, 'go')
    rotn = c(rotn, 'cw')
    meanvals = c(meanvals, meandf$go_cw[m])
    
    subj = c(subj, m)
    mvt = c(mvt, 'go')
    rotn = c(rotn, 'ccw')
    meanvals = c(meanvals, meandf$go_ccw[m])
    
    subj = c(subj, m)
    mvt = c(mvt, 'go')
    rotn = c(rotn, '0')
    meanvals = c(meanvals, meandf$go_0[m])
    
    subj = c(subj, m)
    mvt = c(mvt, 'stop')
    rotn = c(rotn, 'cw')
    meanvals = c(meanvals, meandf$stop_cw[m])
    
    subj = c(subj, m)
    mvt = c(mvt, 'stop')
    rotn = c(rotn, 'ccw')
    meanvals = c(meanvals, meandf$stop_ccw[m])
    
    subj = c(subj, m)
    mvt = c(mvt, 'stop')
    rotn = c(rotn, '0')
    meanvals = c(meanvals, meandf$stop_0[m])
  }
  longfmt = data.frame(subj, mvt, rotn, meanvals)
  longfmt$mvt = factor(longfmt$mvt, levels = c('go', 'stop'))
  longfmt$rotn = factor(longfmt$rotn, levels = c('cw', '0', 'ccw'))
  longfmt$subj = factor(longfmt$subj)
  return(longfmt)
}

organizeFactors <- function(data) {
  rotation = c()
  movement = c()
  for (i in 1:length(data$tripType)) {
    if (data$tripType[i] <= 3) {
      movement = c(movement, 'go')
    } else {
      movement = c(movement, 'stop')
    }
    if (data$tripType[i] == 1 | data$tripType[i] == 4) {
      rotation = c(rotation, 'ccw')
    } else if (data$tripType[i] == 2 | data$tripType[i] == 5){
      rotation = c(rotation, '0')
    } else if (data$tripType[i] == 3 | data$tripType[i] == 6) {
      rotation = c(rotation, 'cw')
    }
  }
  data$Rotation = factor(rotation, levels = c('cw', '0', 'ccw'))
  data$Movement = factor(movement, levels = c('go', 'stop'))
  data$Subject_ID = factor(data$Subject_ID)
  return(data)
}

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sem = sd(x[[col]], na.rm=TRUE)/sqrt(length(x[[col]])))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

makeMedianXTmDf <- function(data) {
  medians_go_STL_cw <- tapply(data$Single_triplet_learning[data$Movement == 'go'
                                                           & data$Rotation == 'cw'],
                              data$Subject_ID[data$Movement == 'go'
                                              & data$Rotation == 'cw'],
                              median, na.rm = T)
  medians_go_STL_ccw <- tapply(data$Single_triplet_learning[data$Movement == 'go'
                                                            & data$Rotation == 'ccw'],
                               data$Subject_ID[data$Movement == 'go'
                                               & data$Rotation == 'ccw'],
                               median, na.rm = T)
  medians_go_RET_cw <- tapply(data$Retention_to_next_triplet[data$Movement == 'go'
                                                             & data$Rotation == 'cw'],
                              data$Subject_ID[data$Movement == 'go'
                                              & data$Rotation == 'cw'],
                              median, na.rm = T)
  medians_go_RET_ccw <- tapply(data$Retention_to_next_triplet[data$Movement == 'go'
                                                              & data$Rotation == 'ccw'],
                               data$Subject_ID[data$Movement == 'go'
                                               & data$Rotation == 'ccw'],
                               median, na.rm = T)
  
  medians_stop_STL_cw <- tapply(data$Single_triplet_learning[data$Movement == 'stop'
                                                             & data$Rotation == 'cw'],
                                data$Subject_ID[data$Movement == 'stop'
                                                & data$Rotation == 'cw'],
                                median, na.rm = T)
  medians_stop_STL_ccw <- tapply(data$Single_triplet_learning[data$Movement == 'stop'
                                                              & data$Rotation == 'ccw'],
                                 data$Subject_ID[data$Movement == 'stop'
                                                 & data$Rotation == 'ccw'],
                                 median, na.rm = T)
  medians_stop_RET_cw <- tapply(data$Retention_to_next_triplet[data$Movement == 'stop'
                                                               & data$Rotation == 'cw'],
                                data$Subject_ID[data$Movement == 'stop'
                                                & data$Rotation == 'cw'],
                                median, na.rm = T)
  medians_stop_RET_ccw <- tapply(data$Retention_to_next_triplet[data$Movement == 'stop'
                                                                & data$Rotation == 'ccw'],
                                 data$Subject_ID[data$Movement == 'stop'
                                                 & data$Rotation == 'ccw'],
                                 median, na.rm = T)
  
  bl = c()
  for (i in 1:length(medians_go_STL_cw)) bl <- c(bl,0)
  perf = c()
  trialsOut = c()
  subj = c()
  mvt = c()
  for (i in 1:length(medians_go_STL_cw)) {
    perf = c(perf, bl[i])
    trialsOut = c(trialsOut, -1)
    subj = c(subj, names(medians_go_STL_cw[i]))
    mvt = c(mvt, 'go')
    
    perf = c(perf, (medians_go_STL_cw[i] + medians_go_STL_ccw[i]*-1)/2)
    trialsOut = c(trialsOut, 1)
    subj = c(subj, names(medians_go_STL_cw[i]))
    mvt = c(mvt, 'go')
    
    perf = c(perf, (medians_go_RET_cw[i] + medians_go_RET_ccw[i]*-1)/2)
    trialsOut = c(trialsOut, 2)
    subj = c(subj, names(medians_go_RET_cw[i]))
    mvt = c(mvt, 'go')
    
    perf = c(perf, bl[i])
    trialsOut = c(trialsOut, -1)
    subj = c(subj, names(medians_stop_STL_cw[i]))
    mvt = c(mvt, 'stop')
    
    perf = c(perf, (medians_stop_STL_cw[i] + medians_stop_STL_ccw[i]*-1)/2)
    trialsOut = c(trialsOut, 1)
    subj = c(subj, names(medians_stop_STL_cw[i]))
    mvt = c(mvt, 'stop')
    
    perf = c(perf, (medians_stop_RET_cw[i] + medians_stop_RET_ccw[i]*-1)/2)
    trialsOut = c(trialsOut, 2)
    subj = c(subj, names(medians_stop_RET_cw[i]))
    mvt = c(mvt, 'stop')
  }
  xtmdf = data.frame(subj, mvt, trialsOut, perf)
}

getSubjectwiseRetRatios <- function(data) {
  ret = data[data$retRatio>-Inf & data$retRatio<Inf,]
  medians_go_cw = tapply(ret$retRatio[ret$Movement == 'go' & ret$Rotation == 'cw'],
                         ret$Subject_ID[ret$Movement == 'go' & ret$Rotation == 'cw'],
                         median, na.rm = T)
  medians_go_ccw = tapply(ret$retRatio[ret$Movement == 'go' & ret$Rotation == 'ccw'],
                          ret$Subject_ID[ret$Movement == 'go' & ret$Rotation == 'ccw'],
                          median, na.rm = T)
  medians_stop_cw = tapply(ret$retRatio[ret$Movement == 'stop' & ret$Rotation == 'cw'],
                           ret$Subject_ID[ret$Movement == 'stop' & ret$Rotation == 'cw'],
                           median, na.rm = T)
  medians_stop_ccw = tapply(ret$retRatio[ret$Movement == 'stop' & ret$Rotation == 'ccw'],
                            ret$Subject_ID[ret$Movement == 'stop' & ret$Rotation == 'ccw'],
                            median, na.rm = T)
  medvals = c()
  mvt = c()
  subj = c()
  for (m in 1:length(medians_go_cw)) {
    medvals = c(medvals, (medians_go_cw[m] + medians_go_ccw[m])/2)
    mvt= c(mvt, 'go')
    subj = c(subj, names(medians_go_cw[m]))
    medvals = c(medvals, (medians_stop_cw[m] + medians_stop_ccw[m])/2)
    mvt= c(mvt, 'stop')
    subj = c(subj, names(medians_stop_cw[m]))
  }
  retratiodf = data.frame(subj, mvt, medvals)
  retratiodf$mvt = factor(retratiodf$mvt, c("go", "stop"))
  return(retratiodf)
}

pairedTAndCohensD <- function(dataA, dataB) {
  testout = t.test(dataA, dataB, paired = TRUE)
  t = testout$statistic
  df = testout$parameter
  p = testout$p.value
  
  testout = cohen.d(dataA, dataB, paired = TRUE, na.rm = TRUE)
  d = testout$estimate
  return(c(t, df, p, d))
}

oneSampTAndCohensD <- function(dataA, mu) {
  testout = t.test(dataA, m = mu)
  t = testout$statistic
  df = testout$parameter
  p = testout$p.value
  
  d = (mu - mean(dataA[is.na(dataA)==FALSE]))/sd(dataA[is.na(dataA)==FALSE])
  return(c(t, df, p, d))
}

retentionRatioTTests <- function(data) {
  govals = data[data$mvt=='go',]
  stopvals = data[data$mvt=='stop',]
  groupA = c()
  groupB = c()
  test = c()
  t = c()
  df = c()
  p = c()
  d = c()
  for (f in 1:3){
    if (f == 1){ 
      if (shapiro.test(govals$medvals)$p.value > 0.05) {
        uset = 1
        statlist = oneSampTAndCohensD(govals$medvals, 0)
      } else {
        uset = 0
        statlist = rstatix::wilcox_test(formula = medvals ~ 1, data = govals)
        r = wilcox_effsize(formula = medvals ~ 1, data = govals)$effsize
      }
      groupA = c(groupA, 'go')
      groupB = c(groupB, '1 samp v zero')
    } else if (f == 2) {
      if (shapiro.test(stopvals$medvals)$p.value > 0.05) {
        uset = 1
        statlist = oneSampTAndCohensD(stopvals$medvals, 0)
      } else {
        uset = 0
        statlist = rstatix::wilcox_test(formula = medvals ~ 1, data = stopvals)
        r = wilcox_effsize(formula = medvals ~ 1, data = stopvals)$effsize
      }
      groupA = c(groupA, 'stop')
      groupB = c(groupB, '1 samp v zero')
    } else if (f == 3) {
      if (shapiro.test(govals$medvals)$p.value > 0.05) {
        uset = 1
        statlist = pairedTAndCohensD(govals$medvals, stopvals$medvals)
      } else {
        uset = 0
        statlist = rstatix::wilcox_test(formula = medvals ~ mvt,
                                        data = data, paired=TRUE)
        r = wilcox_effsize(formula = medvals ~ mvt, data = data,
                           paired = TRUE)$effsize
      }
      groupA = c(groupA, 'go')
      groupB = c(groupB, 'stop')
    }
    if (uset==1) {
      test = c(test, 't')
      t = c(t, statlist[1])
      df = c(df, statlist[2])
      p = c(p, statlist[3])
      d = c(d, statlist[4])
    } else {
      test = c(test, 'wilcox')
      t = c(t, statlist$statistic)
      df = c(df, NA)
      p = c(p, statlist$p)
      d = c(d, r)
    }
  }
  padj = p.adjust(p, method='fdr', n=3)
  retRatioTTests = data.frame(test, groupA, groupB, t, df, p, padj, d)
  return(retRatioTTests)
}

establishFStyleStatsOutput <- function(aovTypeObj, expNum) {
  Experiment <- c()
  Test <- c()
  Effect <- c()
  DFn <- c()
  DFd <- c()
  Fval <- c()
  pval <- c()
  effsize <- c()
  
  for (r in 1:length(aovTypeObj$Effect)) {
    Experiment = c(Experiment, expNum)
    Test = c(Test, "RM ANOVA")
    Effect = c(Effect, aovTypeObj$Effect[r])
    DFn = c(DFn, aovTypeObj$DFn[r])
    DFd = c(DFd, aovTypeObj$DFd[r])
    Fval = c(Fval, aovTypeObj$F[r])
    pval = c(pval, aovTypeObj$p[r])
    effsize = c(effsize, aovTypeObj$ges[r])
  }
  
  FStyleStatsDf <- data.frame(Experiment, Test, Effect, DFn, DFd, Fval,
                       pval, effsize)
  
  return(FStyleStatsDf)
}

addLMERToFStyleStatsOutput <- function(aovTypeObj, r2Info, expNum, FStyleStatsDf) {
  Experiment <- c()
  Test <- c()
  Effect <- c()
  DFn <- c()
  DFd <- c()
  Fval <- c()
  pval <- c()
  effsize <- c()
  
  rows = row.names(aovTypeObj)
  for (r in 1:length(rows)) {
    Experiment = c(Experiment, expNum)
    Test = c(Test, "LMM")
    Effect = c(Effect, rows[r])
    DFn = c(DFn, aovTypeObj$NumDF[r])
    DFd = c(DFd, aovTypeObj$DenDF[r])
    Fval = c(Fval, aovTypeObj$F[r])
    pval = c(pval, aovTypeObj$Pr[r])
    effsize = c(effsize, r2Info$Rsq[r2Info$Effect==rows[r]])
  }
  
  tempdf <- data.frame(Experiment, Test, Effect, DFn, DFd, Fval,
                       pval, effsize)
  FStyleStatsDf <- rbind(FStyleStatsDf, tempdf)
  
  return(FStyleStatsDf)
}

establishPostHocStatsOutput <- function(pairwiseOutput, expNum) {
  Experiment <- c()
  Test <- c()
  groupA <- c()
  groupB <- c()
  Statistic <- c()
  df <- c()
  padjval <- c()
  effsize <- c()
  
  for (r in 1:length(pairwiseOutput$test)) {
    Experiment <- c(Experiment, expNum)
    Test <- c(Test, pairwiseOutput$test[r])
    groupA <- c(groupA, pairwiseOutput$groupA[r])
    groupB <- c(groupB, pairwiseOutput$groupB[r])
    Statistic <- c(Statistic, pairwiseOutput$t[r])
    df <- c(df, pairwiseOutput$df[r])
    padjval <- c(padjval, pairwiseOutput$padj[r])
    effsize <- c(effsize, pairwiseOutput$d[r])
  }
  
  postHocStatsOutput <- data.frame(Experiment, Test, groupA, groupB, Statistic,
                                   df, padjval, effsize)
    
  return(postHocStatsOutput)
}

addEMMContrastsToPostHocStatsOutput <- function(EMMContrasts, contrastEffsize,
                                                expNum, postHocStatsOutput) {
  Experiment <- c()
  Test <- c()
  groupA <- c()
  groupB <- c()
  Statistic <- c()
  df <- c()
  padjval <- c()
  effsize <- c()
  
  for (r in 1:length(EMMContrasts$contrast)) {
    tempstrings = strsplit(EMMContrasts$contrast[r], ' - ')
    thisGroupA = tempstrings[[1]][1]
    thisGroupB = tempstrings[[1]][2]
    
    Experiment <- c(Experiment, expNum)
    Test <- c(Test, 'EMM Contrast (t)')
    groupA <- c(groupA, thisGroupA)
    groupB <- c(groupB, thisGroupB)
    Statistic <- c(Statistic, EMMContrasts$t.ratio[r])
    df <- c(df, EMMContrasts$df[r])
    padjval <- c(padjval, EMMContrasts$p.value[r])
    
    found = FALSE
    idx = 0
    while (idx < length(contrastEffsize$contrast) && found == FALSE) {
      idx = idx + 1
      tempstrings = strsplit(contrastEffsize$contrast[idx], ' - ')
      checkA = tempstrings[[1]][1]
      checkB = tempstrings[[1]][2]
      if ((checkA == thisGroupA && checkB == thisGroupB) ||
          (checkA == thisGroupB && checkB == thisGroupA))  {
        found = TRUE
      }
    }
    
    if (found) {
      effsize <- c(effsize, contrastEffsize$effect.size[idx])
    } else {
      effsize <- c(effsize, NA)
    }
  }
  
  tempdf <- data.frame(Experiment, Test, groupA, groupB, Statistic,
                                   df, padjval, effsize)
  postHocStatsOutput <- rbind(postHocStatsOutput, tempdf)
  
  return(postHocStatsOutput)
}

logRetentionCompStatsOutput <- function(pairwiseOutput, expNum, retCompStatsOutput) {
  Experiment <- c()
  Test <- c()
  groupA <- c()
  groupB <- c()
  Statistic <- c()
  df <- c()
  padjval <- c()
  effsize <- c()
  
  for (r in 1:length(pairwiseOutput$test)) {
    Experiment <- c(Experiment, expNum)
    Test <- c(Test, pairwiseOutput$test[r])
    groupA <- c(groupA, pairwiseOutput$groupA[r])
    groupB <- c(groupB, pairwiseOutput$groupB[r])
    Statistic <- c(Statistic, pairwiseOutput$t[r])
    df <- c(df, pairwiseOutput$df[r])
    padjval <- c(padjval, pairwiseOutput$padj[r])
    effsize <- c(effsize, pairwiseOutput$d[r])
  }
  
  if (length(retCompStatsOutput) == 0) {
    retCompStatsOutput <- data.frame(Experiment, Test, groupA, groupB, Statistic,
                                     df, padjval, effsize)
  } else {
    tempdf <- data.frame(Experiment, Test, groupA, groupB, Statistic,
                                     df, padjval, effsize)
    retCompStatsOutput <- rbind(retCompStatsOutput, tempdf)
  }
  
  return(retCompStatsOutput)
}

### LOAD DATA ------------------------------------------------------------------
E1_anovadata <-  read.csv('STL_IRLRotn_ANOVAForm.csv')  # Experiment 1 data
E1_longformdata <-  read.csv('STL_IRLRotn_longForm.csv')

tripletdata <- read.csv('tripletdata_onlineStudies.csv') # Experiments 2-5 data
if (excludeBasedOnInstrRecall == 1) {
  E2Data = tripletdata[tripletdata$Experiment=='Rotation_Online'
                               & tripletdata$includeSubj == 1
                               & tripletdata$includeTrip == 1,]
  E3Data = tripletdata[tripletdata$Experiment=='Clamp_Online'
                          & tripletdata$includeSubj == 1
                          & tripletdata$includeTrip == 1,]
  E4Data = tripletdata[tripletdata$Experiment=='Rotation_Online_0go'
                           & tripletdata$includeSubj == 1
                           & tripletdata$includeTrip == 1,]
  E5Data = tripletdata[tripletdata$Experiment=='Clamp_Online_0go'
                             & tripletdata$includeSubj == 1
                             & tripletdata$includeTrip == 1,]
} else {
  E2Data = tripletdata[tripletdata$Experiment=='Rotation_Online'
                                & tripletdata$includeTrip == 1,]
  E3Data = tripletdata[tripletdata$Experiment=='Clamp_Online'
                                & tripletdata$includeTrip == 1,]
  E4Data = tripletdata[tripletdata$Experiment=='Rotation_Online_0go'
                                & tripletdata$includeTrip == 1,]
  E5Data = tripletdata[tripletdata$Experiment=='Clamp_Online_0go'
                                & tripletdata$includeTrip == 1,]
}



#### ----------------------- EXPERIMENT 1 ----------------------------------####
E1_anovadata$PID <- as.factor(E1_anovadata$Subject)                    # prepare data for ANOVA
E1_anovadata$MovementCondition[
  E1_anovadata$MovementCondition=='nogo'] = 'stop'                           # recoding data tags for legibility
E1_anovadata$Movement <- factor(E1_anovadata$MovementCondition,
                                       c('go','stop'))
E1_anovadata$Rotation[E1_anovadata$Rotation==15] = 'ccw'
E1_anovadata$Rotation[E1_anovadata$Rotation==-15] = 'cw'
E1_anovadata$Rotation <- factor(E1_anovadata$Rotation, c('cw', 'ccw'))

aov_E1 = anova_test(
  data = E1_anovadata,
  formula = STL ~ Rotation * Movement + Error(PID/(Rotation*Movement)),
  dv = STL, wid = PID, within = c(Rotation, Movement)
)                                                                                  # Run 2-way within-within ANOVA
get_anova_table(aov_E1)
anovasummarydf = data_summary(E1_anovadata, varname='STL',
                              groupnames=c('Movement', 'Rotation'))
outputdf_FStyleStats <- establishFStyleStatsOutput(aov_E1, 1) # log anova style stats for output


groupA = c()
groupB = c()
test = c()
t = c()
df = c()
p = c()
d = c()
for (i in 1:4) # compute post hoc comparisons
{
  if (i == 1) {
    rotA = 'ccw'
    movA = 'go'
    rotB = 'ccw'
    movB = 'stop'
  } else if (i == 2) {
    rotA = 'cw'
    movA = 'go'
    rotB = 'cw'
    movB = 'stop'
  } else if (i == 3) {
    rotA = 'ccw'
    movA = 'go'
    rotB = 'cw'
    movB = 'go'
  } else if (i == 4) {
    rotA = 'ccw'
    movA = 'stop'
    rotB = 'cw'
    movB = 'stop'
  }
  
  groupA = c(groupA, paste(movA, rotA))
  groupB = c(groupB, paste(movB, rotB))
  adata = E1_anovadata[
    (E1_anovadata$Rotation==rotA)
    & (E1_anovadata$Movement==movA),]$STL
  bdata = E1_anovadata[(E1_anovadata$Rotation==rotB) & 
                               (E1_anovadata$Movement==movB),]$STL
  bothdata = E1_anovadata[((E1_anovadata$Rotation==rotA)
                                 & (E1_anovadata$Movement==movA)) |
                                  ((E1_anovadata$Rotation==rotB) & 
                                  (E1_anovadata$Movement==movB)),]
  if (shapiro.test(adata -bdata)$p.value > 0.05) {
    statlist = pairedTAndCohensD(adata, bdata)
    test = c(test, 't')
    t = c(t, statlist[1])
    df = c(df, statlist[2])
    p = c(p, statlist[3])
    d = c(d, statlist[4])
  } else {
    if (movA == movB) {
      # groups differ by rotation condition
      statlist = rstatix::wilcox_test(formula = STL ~ Rotation,
                                      data = bothdata, paired = TRUE)
      r = wilcox_effsize(formula = STL ~ Rotation, data = bothdata,
                         paired = TRUE)$effsize
    } else if (rotA == rotB) {
      # groups differ by movement condition
      statlist = rstatix::wilcox_test(formula = STL ~ Movement,
                                      data = bothdata, paired = TRUE)
      r = wilcox_effsize(formula = STL ~ Movement, data = bothdata,
                         paired = TRUE)$effsize
    } else {
      print('PROBLEM HERE')
    }
    test = c(test, 'wilcox')
    t = c(t, statlist$statistic)
    df = c(df, NA)
    p = c(p, statlist$p)
    d = c(d, r)
  }
}
padj = p.adjust(p, method='fdr', n=4)
postHocs_E1 = data.frame(test, groupA, groupB,
                                    t, df, p, padj, d)
postHocs_E1
outputdf_postHocTests <- establishPostHocStatsOutput(postHocs_E1, 1)  # log post-hoc tests for output later

# now get performance x time
E1_longformdata$Subject_ID = as.factor(E1_longformdata$Subject)
E1_longformdata$Movement = factor(E1_longformdata$MovementCondition,
                                        c('go', 'stop'))
E1_longformdata$Rotation = factor(E1_longformdata$Rotation,
                                        c('cw', 'ccw'))
E1_longformdata$Single_triplet_learning = E1_longformdata$STL
E1_longformdata$Retention_to_next_triplet = E1_longformdata$RetainedRaw
E1_longformdata$retRatio = E1_longformdata$Retained

xtmdf_E1 = makeMedianXTmDf(E1_longformdata)
df2_E1 = data_summary(xtmdf_E1, varname='perf', groupnames=c('mvt', 'trialsOut'))


## now check retention values
retratiodf_E1 = getSubjectwiseRetRatios(E1_longformdata)
retratiodf2_E1 = data_summary(retratiodf_E1, varname='medvals', groupnames=c('mvt'))
retRatioTtests_E1 = retentionRatioTTests(retratiodf_E1)
retRatioTtests_E1
outputdf_retentionComparisons <- logRetentionCompStatsOutput(retRatioTtests_E1,
                                                             1, c()) # log one-off retention comparisons for output


# set up plots for Figure 2d and e
fig2_d_left <- ggplot(anovasummarydf[anovasummarydf$Movement == 'go',],
                       aes(x=Rotation, y=STL, group=Movement, color=Movement)) +
  geom_hline(yintercept=0, linetype='dashed', color='black') +
  geom_pointrange(aes(ymin=STL-sem, ymax=STL+sem))+
  theme(legend.position = c(0.8, 0.8)) +
  ggtitle("d. E1 Movement") +
  ylim(-6,6)+
  ylab('Single-Trial Learning (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

fig2_d_right <- ggplot(anovasummarydf[anovasummarydf$Movement == 'stop',],
                       aes(x=Rotation, y=STL, group=Movement, color=Movement)) + 
  geom_hline(yintercept=0, linetype='dashed', color='black') +
  geom_pointrange(aes(ymin=STL-sem, ymax=STL+sem))+
  theme(legend.position = c(0.8, 0.8)) +
  ggtitle("No-Movement") +
  ylim(-6,6)+
  ylab('Single-Trial Learning (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

fig2_e_left <- ggplot(df2_E1, aes(x=trialsOut, y=perf, group=mvt, color=mvt)) + 
  geom_line()+
  geom_pointrange(aes(ymin=perf-sem, ymax=perf+sem))+
  theme(legend.position = c(0.25, 0.8)) +
  ggtitle("e. E1 Retention") +
  ylim(0,4)+
  ylab('delta Hand Angle (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

fig2_e_right <- ggplot(retratiodf2_E1, aes(x=mvt, y=medvals, color=mvt)) +
  geom_pointrange(aes(ymin=medvals-sem, ymax=medvals+sem))+
  scale_y_continuous(breaks=seq(0,1,by=0.25))+
  expand_limits(y = c(0, 1)) +
  theme(legend.position = c(0.8, 0.8)) +
  ggtitle("Ratio") +
  ylab('Remembered STL (ratio)') +
  xlab('Movement') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))



#### ------------------------- EXPERIMENT 2 ------------------------------- ####
E2Data = organizeFactors(E2Data)

meandf_E2 = splitRotnMvtMeans(E2Data)
longfmt_E2 = meanDfToLongForm(meandf_E2)

E2_crossedRE = lmerTest::lmer(
  Single_triplet_learning ~ 1 + Rotation * Movement
  + (1 | Subject_ID),
  data = E2Data
)
E2.stl.aov<-anova(E2_crossedRE, type=3, ddf="Kenward-Roger")
E2.stl.aov
E2.r2 <- r2beta(E2_crossedRE, method = 'kr')
E2.r2
outputdf_FStyleStats <- addLMERToFStyleStatsOutput(E2.stl.aov, E2.r2, 2,
                                                   outputdf_FStyleStats) # log anova style stats for output

# post hoc tests
E2.emm.s <- emmeans(E2_crossedRE, c("Movement", "Rotation"))

# specify emmeans to test
gocw_lookup = c(1, 0, 0, 0, 0, 0) # lookup vector for first emmean in grid
stopcw_lookup = c(0, 1, 0, 0, 0, 0)
go0_lookup = c(0, 0, 1, 0, 0, 0)
stop0_lookup = c(0, 0, 0, 1, 0, 0)
goccw_lookup = c(0, 0, 0, 0, 1, 0)
stopccw_lookup = c(0, 0, 0, 0, 0, 1)
contrastOut <- summary(contrast(E2.emm.s, method = list("go cw - stop cw" = gocw_lookup - stopcw_lookup,
                                                              "go cw - go 0" = gocw_lookup - go0_lookup,
                                                              "go cw - go ccw" = gocw_lookup - goccw_lookup,
                                                              "go 0 - go ccw" = go0_lookup - goccw_lookup,
                                                              "go 0 - stop 0" = go0_lookup - stop0_lookup,
                                                              "go ccw - stop ccw" = goccw_lookup - stopccw_lookup,
                                                              "stop cw - stop 0" = stopcw_lookup - stop0_lookup,
                                                              "stop cw - stop ccw" = stopcw_lookup - stopccw_lookup,
                                                              "stop 0 - stop ccw" = stop0_lookup - stopccw_lookup),
                                adjust = "fdr"))
contrastOut
effectSizeOut <- summary(eff_size(E2.emm.s,
                                  sigma = sigma(E2_crossedRE),
                                  edf = df.residual(E2_crossedRE))) # returns cohen's d: https://cran.r-project.org/web/packages/emmeans/vignettes/
effectSizeOut
outputdf_postHocTests <- addEMMContrastsToPostHocStatsOutput(contrastOut, effectSizeOut,
                                                2, outputdf_postHocTests) # log post hoc stats for output

## collapse STL and remembered STL across rotation directions
xtmdf_E2 = makeMedianXTmDf(E2Data)
df2_E2 = data_summary(xtmdf_E2, varname='perf', groupnames=c('mvt', 'trialsOut'))

## now check retention values
retratiodf_E2 = getSubjectwiseRetRatios(E2Data)
retratiodf2_E2 = data_summary(retratiodf_E2, varname='medvals', groupnames=c('mvt'))
retRatioTtests_E2 = retentionRatioTTests(retratiodf_E2)
retRatioTtests_E2
outputdf_retentionComparisons <- logRetentionCompStatsOutput(retRatioTtests_E2,
                                                             2, outputdf_retentionComparisons) # log one-off retention comparisons for output

# Set up plots for Figure 2g and h, also plot data for SI Appendix Supplemental Figure 2a
SupFig2_a <- ggplot(longfmt_E2, aes(x=rotn, y=meanvals, colour=mvt)) +
  geom_hline(yintercept=0, linetype='dashed', color='black') +
  geom_boxplot(outlier.shape = NA, notch=TRUE, width=0.5) +
  scale_y_continuous(breaks=seq(-10,10,by=5))+
  expand_limits(y = c(-10, 10)) +
  theme(legend.position = c(0.9, 0.9)) +
  ggtitle("a. E2 STL") +
  ylab('STL (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

Fig2_g_left <- ggplot(as.data.frame(E2.emm.s)[as.data.frame(E2.emm.s)$Movement == 'go',],
                           aes(x = Rotation, y = emmean, color = Movement)) +
  geom_hline(yintercept=0, linetype='dashed', color='black') +
  geom_pointrange(aes(ymin = emmean - (emmean-lower.CL), ymax = emmean + (upper.CL-emmean))) +
  ylim(-6, 6) +
  theme(legend.position = c(0.8, 0.8)) +
  ggtitle("g. E2 Movement") +
  ylab('STL (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

Fig2_g_right <- ggplot(as.data.frame(E2.emm.s)[as.data.frame(E2.emm.s)$Movement == 'stop',],
                             aes(x = Rotation, y = emmean, color = Movement)) +
  geom_hline(yintercept=0, linetype='dashed', color='black') +
  geom_pointrange(aes(ymin = emmean - (emmean-lower.CL), ymax = emmean + (upper.CL-emmean))) +
  ylim(-6, 6) +
  theme(legend.position = c(0.8, 0.8)) +
  ggtitle("No-Movement") +
  ylab('STL (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

Fig2_h_left <- ggplot(df2_E2, aes(x=trialsOut, y=perf, group=mvt, color=mvt)) + 
  geom_line()+
  geom_pointrange(aes(ymin=perf-sem, ymax=perf+sem))+
  theme(legend.position = c(0.25, 0.8)) +
  ggtitle("h. E2 Retention") +
  ylim(0,4)+
  ylab('delta Hand Angle (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

Fig2_h_right <- ggplot(retratiodf2_E2, aes(x=mvt, y=medvals, color=mvt)) +
  geom_pointrange(aes(ymin=medvals-sem, ymax=medvals+sem))+
  scale_y_continuous(breaks=seq(0,1,by=0.25))+
  expand_limits(y = c(0, 1)) +
  theme(legend.position = c(0.8, 0.8)) +
  ggtitle("Ratio") +
  ylab('Remembered STL/STL') +
  xlab('Movement') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))


#### ------------------------- EXPERIMENT 3  ------------------------------ ####
E3Data = organizeFactors(E3Data)

meandf_E3 = splitRotnMvtMeans(E3Data)
longfmt_E3 = meanDfToLongForm(meandf_E3)

E3_crossedRE = lmerTest::lmer(
  Single_triplet_learning ~ 1 + Rotation * Movement
  + (1 | Subject_ID),
  data = E3Data
)
E3.stl.aov<-anova(E3_crossedRE, type=3, ddf="Kenward-Roger")
E3.stl.aov
E3.r2 = r2beta(E3_crossedRE, method = 'kr')
E3.r2
outputdf_FStyleStats <- addLMERToFStyleStatsOutput(E3.stl.aov, E3.r2, 3,
                                                   outputdf_FStyleStats) # log anova style stats for output

# post hoc tests
E3.emm.s <- emmeans(E3_crossedRE, c("Movement", "Rotation"))
contrastOut <- summary(contrast(E3.emm.s, method = list("go cw - stop cw" = gocw_lookup - stopcw_lookup,
                                                              "go cw - go 0" = gocw_lookup - go0_lookup,
                                                              "go cw - go ccw" = gocw_lookup - goccw_lookup,
                                                              "go 0 - go ccw" = go0_lookup - goccw_lookup,
                                                              "go 0 - stop 0" = go0_lookup - stop0_lookup,
                                                              "go ccw - stop ccw" = goccw_lookup - stopccw_lookup,
                                                              "stop cw - stop 0" = stopcw_lookup - stop0_lookup,
                                                              "stop cw - stop ccw" = stopcw_lookup - stopccw_lookup,
                                                              "stop 0 - stop ccw" = stop0_lookup - stopccw_lookup),
                                adjust = "fdr"))
contrastOut
effectSizeOut <- summary(eff_size(E3.emm.s, sigma = sigma(E3_crossedRE),
                                  edf = df.residual(E3_crossedRE)))
outputdf_postHocTests <- addEMMContrastsToPostHocStatsOutput(contrastOut, effectSizeOut,
                                                             3, outputdf_postHocTests) # log post hoc stats for output

## collapse STL and remembered STL across rotation directions
xtmdf_E3 = makeMedianXTmDf(E3Data)
df2_E3 = data_summary(xtmdf_E3, varname='perf', groupnames=c('mvt', 'trialsOut'))

## now check retention values
retratiodf_E3 = getSubjectwiseRetRatios(E3Data)
retratiodf2_E3 = data_summary(retratiodf_E3, varname='medvals', groupnames=c('mvt'))
retRatioTtests_E3 = retentionRatioTTests(retratiodf_E3)
retRatioTtests_E3
outputdf_retentionComparisons <- logRetentionCompStatsOutput(retRatioTtests_E3,
                                                             3, outputdf_retentionComparisons) # log retention comparisons for output later

# PLOTTING - Fig2j and k, as well as SI Appendix Supplemental Fig 2 b
SupFig2_b <- ggplot(longfmt_E3, aes(x=rotn, y=meanvals, colour=mvt)) +
  geom_hline(yintercept=0, linetype='dashed', color='black') +
  geom_boxplot(outlier.shape = NA, notch=TRUE, width=0.5) +
  scale_y_continuous(breaks=seq(-10,10,by=5), limits=c(-10,10))+
  expand_limits(y = c(-10, 10)) +
  theme(legend.position = c(0.9, 0.9)) +
  ggtitle("b. E3 STL") +
  ylab('STL (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

Fig2_j_left <- ggplot(as.data.frame(E3.emm.s)[as.data.frame(E3.emm.s)$Movement == 'go',],
                           aes(x = Rotation, y = emmean, color = Movement)) +
  geom_hline(yintercept=0, linetype='dashed', color='black') +
  geom_pointrange(aes(ymin = emmean - (emmean-lower.CL), ymax = emmean + (upper.CL-emmean))) +
  ylim(-6, 6) +
  theme(legend.position = c(0.8, 0.8)) +
  ggtitle("j. E3 Movement") +
  ylab('STL (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

Fig2_j_right <- ggplot(as.data.frame(E3.emm.s)[as.data.frame(E3.emm.s)$Movement == 'stop',],
                             aes(x = Rotation, y = emmean, color = Movement)) +
  geom_hline(yintercept=0, linetype='dashed', color='black') +
  geom_pointrange(aes(ymin = emmean - (emmean-lower.CL), ymax = emmean + (upper.CL-emmean))) +
  ylim(-6, 6) +
  theme(legend.position = c(0.8, 0.8)) +
  ggtitle("No-Movement") +
  ylab('STL (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

Fig2_k_left <- ggplot(df2_E3, aes(x=trialsOut, y=perf, group=mvt, color=mvt)) + 
  geom_line()+
  geom_pointrange(aes(ymin=perf-sem, ymax=perf+sem))+
  theme(legend.position = c(0.25, 0.8)) +
  ggtitle("k. E3 Retention") +
  ylim(0,4)+
  ylab('delta Hand Angle (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

Fig2_k_right <- ggplot(retratiodf2_E3, aes(x=mvt, y=medvals, color=mvt)) +
  geom_pointrange(aes(ymin=medvals-sem, ymax=medvals+sem))+
  scale_y_continuous(breaks=seq(0,1,by=0.25))+
  expand_limits(y = c(0, 1)) +
  theme(legend.position = c(0.8, 0.8)) +
  ggtitle("Ratio") +
  ylab('Remembered STL/STL') +
  xlab('Movement') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

### ---------------------------- EXPERIMENT 4 ------------------------------ ###
E4Data = organizeFactors(E4Data)

meandf_E4 = splitRotnMvtMeans(E4Data)
longfmt_E4 = meanDfToLongForm(meandf_E4)

E4_model = lmerTest::lmer(
  Single_triplet_learning ~ 1 + Rotation
  + (1 | Subject_ID),
  data = E4Data[E4Data$Movement == 'stop',]
)
E4.stl.aov<-anova(E4_model, type=3, ddf="Kenward-Roger")
E4.stl.aov
E4.r2 = r2beta(E4_model, method = 'kr')
E4.r2
outputdf_FStyleStats <- addLMERToFStyleStatsOutput(E4.stl.aov, E4.r2, 4,
                                                   outputdf_FStyleStats) # log anova style stats for output

# post hoc tests
E4.emm.s <- emmeans(E4_model, c("Rotation"))
zero_lookup = c(0, 1, 0) # lookup vector for first emmean in grid
cw_lookup = c(1, 0)
ccw_lookup = c(0, 1)
contrastOut <- summary(contrast(E4.emm.s, method = list("cw - ccw" = cw_lookup - ccw_lookup),
                                adjust = "fdr"))
effectSizeOut <- summary(eff_size(E4.emm.s, sigma = sigma(E4_model),
                                  edf = df.residual(E4_model)))
contrastOut
effectSizeOut
outputdf_postHocTests <- addEMMContrastsToPostHocStatsOutput(contrastOut, effectSizeOut,
                                                             4, outputdf_postHocTests) # log post hoc stats for output

## collapse STL and remembered STL across rotation directions
xtmdf_E4 = makeMedianXTmDf(E4Data)
df2_E4 = data_summary(xtmdf_E4, varname='perf', groupnames=c('mvt', 'trialsOut'))

## now check retention values
retratiodf_E4 = getSubjectwiseRetRatios(E4Data)
retratiodf2_E4 = data_summary(retratiodf_E4, varname='medvals', groupnames=c('mvt'))
Experiment <- 4
groupA <- 'stop'
groupB <- '1 samp v zero'
shapiro.test(retratiodf_E4$medvals[retratiodf_E4$mvt=='stop'])
Test = 't'
e4retcomp = t.test(retratiodf_E4$medvals[retratiodf_E4$mvt=='stop'], mu=0)
Statistic <- e4retcomp$statistic
df <- e4retcomp$parameter
padjval <- e4retcomp$p.value # no other comparisons in this test family, so no need to adjust p value here
effsize = (0 - mean(retratiodf_E4$medvals[retratiodf_E4$mvt=='stop' & 
                                          is.na(
                                            retratiodf_E4$medvals)==FALSE])
     )/sd(retratiodf_E4$medvals[retratiodf_E4$mvt=='stop' & 
                                      is.na(
                                        retratiodf_E4$medvals)==FALSE])
tempdf <- data.frame(Experiment, Test, groupA, groupB, Statistic, df, padjval, effsize)
outputdf_retentionComparisons <- rbind(outputdf_retentionComparisons, tempdf)  # log retention test results for output

# PLOTTING - FIGURE 3, TOP PANEL
fig3_c_left <- ggplot(longfmt_E4, aes(x=rotn, y=meanvals, colour=mvt)) +
  geom_hline(yintercept=0, linetype='dashed', color='black') +
  geom_boxplot(outlier.shape = NA, notch=TRUE, width=0.5) +
  scale_y_continuous(breaks=seq(-6,6,by=4), limits = c(-6,6))+
  expand_limits(y = c(-6, 6)) +
  theme(legend.position = c(0.9, 0.9)) +
  ggtitle("c. E4 Rotation") +
  ylab('STL (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

fig3_d_left <- ggplot(as.data.frame(E4.emm.s),
                              aes(x = Rotation, y = emmean)) +
  geom_hline(yintercept=0, linetype='dashed', color='black') +
  geom_pointrange(aes(ymin = emmean - (emmean-lower.CL), ymax = emmean + (upper.CL-emmean))) +
  ylim(-4, 4) +
  theme(legend.position = c(0.8, 0.8)) +
  ggtitle("d. E4 Rotation") +
  ylab('STL (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

fig3_e_left <- ggplot(df2_E4[df2_E4$mvt=='stop',], aes(x=trialsOut, y=perf)) + 
  geom_line()+
  geom_pointrange(aes(ymin=perf-sem, ymax=perf+sem))+
  theme(legend.position = c(0.25, 0.8)) +
  ggtitle("e. E4 Rotation") +
  ylim(-0.1,4)+
  ylab('delta Hand Angle (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

fig3_retRatioPlot <- ggplot(retratiodf2_E4[retratiodf2_E4$mvt=='stop',], aes(x=mvt, y=medvals)) +
  geom_pointrange(aes(ymin=medvals-sem, ymax=medvals+sem))+
  scale_y_continuous(breaks=seq(0,1,by=0.25))+
  expand_limits(y = c(0, 1)) +
  theme(legend.position = c(0.8, 0.8)) +
  ylab('Remembered STL/STL') +
  xlab('Movement') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))


### ----------------------------- EXPERIMENT 5 ----------------------------- ###
E5Data = organizeFactors(E5Data)

meandf_E5 = splitRotnMvtMeans(E5Data)
longfmt_E5 = meanDfToLongForm(meandf_E5)

E5_model = lmerTest::lmer(
  Single_triplet_learning ~ 1 + Rotation
  + (1 | Subject_ID),
  data = E5Data[E5Data$Movement == 'stop',]
)
E5.stl.aov<-anova(E5_model, type=3, ddf="Kenward-Roger")
E5.stl.aov
E5.r2 = r2beta(E5_model, method = 'kr')
E5.r2
outputdf_FStyleStats <- addLMERToFStyleStatsOutput(E5.stl.aov, E5.r2, 5,
                                                   outputdf_FStyleStats) # log anova style stats for output

# post hoc tests
E5.emm.s <- emmeans(E5_model, c("Rotation"))
zero_lookup = c(0, 1, 0) # lookup vector for first emmean in grid
cw_lookup = c(1, 0, 0)
ccw_lookup = c(0, 0, 1)
contrastOut <- summary(contrast(E5.emm.s, method = list("cw - ccw" = cw_lookup - ccw_lookup,
                                                              "cw - 0" = cw_lookup - zero_lookup,
                                                              "0 - ccw" = zero_lookup - ccw_lookup),
                                adjust = "fdr"))
effectSizeOut <- summary(eff_size(E5.emm.s, sigma = sigma(E5_model),
                                  edf = df.residual(E5_model)))
contrastOut
effectSizeOut
outputdf_postHocTests <- addEMMContrastsToPostHocStatsOutput(contrastOut, effectSizeOut,
                                                             5, outputdf_postHocTests) # log post hoc stats for output

## collapse STL and remembered STL across rotation directions
xtmdf_E5 = makeMedianXTmDf(E5Data)
df2_E5 = data_summary(xtmdf_E5, varname='perf', groupnames=c('mvt', 'trialsOut'))

## now check retention values
retratiodf_E5 = getSubjectwiseRetRatios(E5Data)
retratiodf2_E5 = data_summary(retratiodf_E5, varname='medvals', groupnames=c('mvt'))
shapiro.test(retratiodf_E5$medvals[retratiodf_E5$mvt=='stop'])
e5retcomp = rstatix::wilcox_test(formula = medvals ~ 1,
                    data = retratiodf_E5[retratiodf_E5$mvt=='stop',])
e5reteffsize = rstatix::wilcox_effsize(formula = medvals ~ 1,
                     data = retratiodf_E5[retratiodf_E5$mvt=='stop',])
quantile(retratiodf_E5$medvals[retratiodf_E5$mvt=='stop'],
         c(0.25, 0.5, 0.75), na.rm = TRUE)

Experiment <- 5
groupA <- 'stop'
groupB <- '1 samp v zero'
Test = 'wilcox'
Statistic <- e5retcomp$statistic
df <- NA
padjval <- e5retcomp$p # no other comparisons in this test family, so no need to adjust p value here
effsize = e5reteffsize$effsize
tempdf <- data.frame(Experiment, Test, groupA, groupB, Statistic, df, padjval, effsize)
outputdf_retentionComparisons <- rbind(outputdf_retentionComparisons, tempdf) # log retention test results for output



# PLOTTING - FIGURE 3, BOTTOM ROW
fig3_c_right <- ggplot(longfmt_E5, aes(x=rotn, y=meanvals)) +
  geom_hline(yintercept=0, linetype='dashed', color='black') +
  geom_boxplot(outlier.shape = NA, notch=TRUE, width=0.5) +
  scale_y_continuous(breaks=seq(-8,8,by=4))+
  expand_limits(y = c(-8, 8)) +
  theme(legend.position = c(0.9, 0.9)) +
  ggtitle("E5 Error-Clamp") +
  ylab('STL (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

fig3_d_right <- ggplot(as.data.frame(E5.emm.s),
                                    aes(x = Rotation, y = emmean)) +
  geom_hline(yintercept=0, linetype='dashed', color='black') +
  geom_pointrange(aes(ymin = emmean - (emmean-lower.CL), ymax = emmean + (upper.CL-emmean))) +
  ylim(-4, 4) +
  theme(legend.position = c(0.8, 0.8)) +
  ggtitle("E5 Error-Clamp") +
  ylab('STL (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

fig3_e_right <- ggplot(df2_E5[df2_E5$mvt=='stop',], aes(x=trialsOut, y=perf)) + 
  geom_line()+
  geom_pointrange(aes(ymin=perf-sem, ymax=perf+sem))+
  theme(legend.position = c(0.25, 0.8)) +
  ylim(-0.1,4)+
  ggtitle("E5 Error-Clamp") +
  ylab('delta Hand Angle (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

fig3_retRatioPlot_E5 <- ggplot(retratiodf2_E5[retratiodf2_E5$mvt=='stop',], aes(x=mvt, y=medvals)) +
  geom_pointrange(aes(ymin=medvals-sem, ymax=medvals+sem))+
  scale_y_continuous(breaks=seq(0,1,by=0.25))+
  expand_limits(y = c(0, 1)) +
  theme(legend.position = c(0.8, 0.8)) +
  ylab('Remembered STL/STL') +
  xlab('Movement') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

### -------------------------- SAVE STATS OUTPUT --------------------------- ###
setwd(savedir)
if (excludeBasedOnInstrRecall == 1) {
  fstylefilename = "FStyleStatsOutput_onlyIncludeIfRecalledInstr.csv"
  posthocfilename = "postHocTestOutput_onlyIncludeIfRecalledInstr.csv"
  retentionfilename = "retentionRatioTestOutput_onlyIncludeIfRecalledInstr.csv"
  figpaste = "_onlyIncludeIfRecalledInstr.pdf"
} else {
  fstylefilename = "FStyleStatsOutput_allParticipantsIncluded.csv"
  posthocfilename = "postHocTestOutput_allParticipantsIncluded.csv"
  retentionfilename = "retentionRatioTestOutputl_alParticipantsIncluded.csv"
  figpaste = "_alParticipantsIncluded.pdf"
}
write.csv(outputdf_FStyleStats, fstylefilename, row.names = FALSE)
write.csv(outputdf_postHocTests, posthocfilename, row.names = FALSE)
write.csv(outputdf_retentionComparisons, retentionfilename, row.names = FALSE)


### ----------------- PLOT FIG 2 AND SAVE -------------------- ###
fig2 <- ggarrange(fig2_d_left, fig2_d_right, fig2_e_left, fig2_e_right,
                  Fig2_g_left, Fig2_g_right, Fig2_h_left, Fig2_h_right,
                  Fig2_j_left, Fig2_j_right, Fig2_k_left, Fig2_k_right,
                  ncol = 4, nrow = 3, widths = c(0.5, 0.5, 0.4, 0.3,
                                                 0.5, 0.5, 0.4, 0.3,
                                                 0.5, 0.5, 0.4, 0.3))
ggsave(paste("fig2",figpaste,sep=""),fig2, width = 9, height = 12, units = "in", device = 'pdf')
fig2


### ------------- PLOT SI Fig S2 AND SAVE --------------------- ###
supfig2 <- ggarrange(SupFig2_a, SupFig2_b,
                  ncol = 2, nrow = 1, widths = c(0.5, 0.5))
ggsave(paste("supfig2",figpaste,sep=""),supfig2, width = 6, height = 3, units = "in", device = 'pdf')
supfig2


### --------------- PLOT FIG 3 AND SAVE -------------------- ###
fig3_bottomRow <- ggarrange(fig3_c_left, fig3_c_right, fig3_d_left, fig3_d_right,
                            fig3_e_left, fig3_e_right,
                            ncol = 6, nrow = 1, widths = c(0.3, 0.3, 0.2, 0.2, 0.2, 0.2))
ggsave(paste("fig3",figpaste,sep=""),fig3_bottomRow, width = 14, height = 4, units = "in", device = 'pdf')
fig3_bottomRow
