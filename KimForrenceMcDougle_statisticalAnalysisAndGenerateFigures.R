## This script processes data, performs statistical tests, and generates figures
## for the initial submission of:
##    Title: Learning from the path not taken: Sensory prediction error computation and implicit motor adaptation proceed without movement
##    Authors: Kim, Forrence, McDougle
##
## TODO: make sure you set the working directory to the folder where you
## downloaded this file and the accompanying data


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
#------------------------------ these are mainly for data organization ---------

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

getRightwayMeans <- function(data) {
  rightway_cw_go = data[data$Rotation=='cw' & data$Single_triplet_learning>0
                        & data$Movement == 'go',]
  rightway_ccw_go = data[data$Rotation=='ccw' & data$Single_triplet_learning<0
                         & data$Movement == 'go',]
  rightway_cw_stop = data[data$Rotation=='cw' & data$Single_triplet_learning>0
                          & data$Movement == 'stop',]
  rightway_ccw_stop = data[data$Rotation=='ccw' & data$Single_triplet_learning<0
                           & data$Movement == 'stop',]
  means_cw_go = tapply(rightway_cw_go$Single_triplet_learning, rightway_cw_go$Subject_ID,
                       median, na.rm = T)
  means_ccw_go = tapply(rightway_ccw_go$Single_triplet_learning, rightway_ccw_go$Subject_ID,
                        median, na.rm = T)
  means_cw_stop = tapply(rightway_cw_stop$Single_triplet_learning, rightway_cw_stop$Subject_ID,
                         median, na.rm = T)
  means_ccw_stop = tapply(rightway_ccw_stop$Single_triplet_learning, rightway_ccw_stop$Subject_ID,
                          median, na.rm = T)
  means_go = c()
  means_stop = c()
  for (i in 1:length(means_cw_go)) {
    subject = names(means_cw_go)[i]
    if (is.na(means_cw_go[subject])) {
      means_go = c(means_go, means_ccw_go[subject]*-1)
    } else if (is.na(means_ccw_go[subject])) {
      means_go = c(means_go, means_cw_go[subject])
    } else {
      means_go = c(means_go, (means_cw_go[subject] + means_ccw_go[subject]*-1)/2)
    }
    
    if (is.na(means_cw_stop[subject])) {
      means_stop = c(means_stop, means_ccw_stop[subject]*-1)
    } else if (is.na(means_ccw_stop[subject])) {
      means_stop = c(means_stop, means_cw_stop[subject])
    } else {
      means_stop = c(means_stop, (means_cw_stop[subject] + means_ccw_stop[subject]*-1)/2)
    }
  }
  rightwaySTLdf = data.frame(means_go, means_stop)
  return(rightwaySTLdf)
}

getWrongwayMeans <- function(data) {
  wrongway_cw_go = data[data$Rotation=='cw' & data$Single_triplet_learning<0
                        & data$Movement == 'go',]
  wrongway_ccw_go = data[data$Rotation=='ccw' & data$Single_triplet_learning>0
                         & data$Movement == 'go',]
  wrongway_cw_stop = data[data$Rotation=='cw' & data$Single_triplet_learning<0
                          & data$Movement == 'stop',]
  wrongway_ccw_stop = data[data$Rotation=='ccw' & data$Single_triplet_learning>0
                           & data$Movement == 'stop',]
  means_cw_go = tapply(wrongway_cw_go$Single_triplet_learning, wrongway_cw_go$Subject_ID,
                       median, na.rm = T)
  means_ccw_go = tapply(wrongway_ccw_go$Single_triplet_learning, wrongway_ccw_go$Subject_ID,
                        median, na.rm = T)
  means_cw_stop = tapply(wrongway_cw_stop$Single_triplet_learning, wrongway_cw_stop$Subject_ID,
                         median, na.rm = T)
  means_ccw_stop = tapply(wrongway_ccw_stop$Single_triplet_learning, wrongway_ccw_stop$Subject_ID,
                          median, na.rm = T)
  means_go = c()
  means_stop = c()
  for (i in 1:length(means_cw_go)) {
    subject = names(means_cw_go)[i]
    if (is.na(means_cw_go[subject])) {
      means_go = c(means_go, means_ccw_go[subject]*-1)
    } else if (is.na(means_ccw_go[subject])) {
      means_go = c(means_go, means_cw_go[subject])
    } else {
      means_go = c(means_go, (means_cw_go[subject] + means_ccw_go[subject]*-1)/2)
    }
    
    if (is.na(means_cw_stop[subject])) {
      means_stop = c(means_stop, means_ccw_stop[subject]*-1)
    } else if (is.na(means_ccw_stop[subject])) {
      means_stop = c(means_stop, means_cw_stop[subject])
    } else {
      means_stop = c(means_stop, (means_cw_stop[subject] + means_ccw_stop[subject]*-1)/2)
    }
  }
  wrongwaySTLdf = data.frame(means_go, means_stop)
  return(wrongwaySTLdf)
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

# navigate to directory for the rotation experiment and load files--------------
# load in data for Figure 2
setwd('C:/Users/kimol/Documents/GitHub/LearningFromThePathNotTaken') ## TODO: set working directory
anovadata_rotation <-  read.csv('STL_IRLRotn_ANOVAForm.csv')   
longformdata_rotation <-  read.csv('STL_IRLRotn_longForm.csv')  

# load in data from online experiments
tripletdata <- read.csv('tripletdata_onlineStudies.csv')



#### -------------------- R O T A T I O N ----------------------------------####
anovadata_rotation$PID <- as.factor(anovadata_rotation$Subject)                    # prepare data for ANOVA
anovadata_rotation$MovementCondition[
  anovadata_rotation$MovementCondition=='nogo'] = 'stop' # recoding data tag for legibility
anovadata_rotation$Movement <- factor(anovadata_rotation$MovementCondition,
                                       c('go','stop'))
anovadata_rotation$Rotation[anovadata_rotation$Rotation==15] = 'ccw'
anovadata_rotation$Rotation[anovadata_rotation$Rotation==-15] = 'cw'
anovadata_rotation$Rotation <- factor(anovadata_rotation$Rotation, c('cw', 'ccw'))

aov_Rotation = anova_test(
  data = anovadata_rotation,
  formula = STL ~ Rotation * Movement + Error(PID/(Rotation*Movement)),
  dv = STL, wid = PID, within = c(Rotation, Movement)
)                                                                                    # Run 2-way within-within ANOVA
get_anova_table(aov_Rotation)
anovasummarydf = data_summary(anovadata_rotation, varname='STL',
                              groupnames=c('Movement', 'Rotation'))


groupA = c()
groupB = c()
test = c()
t = c()
df = c()
p = c()
d = c()
for (i in 1:4) # compute post hocs
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
  adata = anovadata_rotation[
    (anovadata_rotation$Rotation==rotA)
    & (anovadata_rotation$Movement==movA),]$STL
  bdata = anovadata_rotation[(anovadata_rotation$Rotation==rotB) & 
                               (anovadata_rotation$Movement==movB),]$STL
  bothdata = anovadata_rotation[((anovadata_rotation$Rotation==rotA)
                                 & (anovadata_rotation$Movement==movA)) |
                                  ((anovadata_rotation$Rotation==rotB) & 
                                  (anovadata_rotation$Movement==movB)),]
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
postHocs_RotationInLab = data.frame(test, groupA, groupB,
                                    t, df, p, padj, d)
postHocs_RotationInLab


# now get performance x time
longformdata_rotation$Subject_ID = as.factor(longformdata_rotation$Subject)
longformdata_rotation$Movement = factor(longformdata_rotation$MovementCondition,
                                        c('go', 'stop'))
longformdata_rotation$Rotation = factor(longformdata_rotation$Rotation,
                                        c('cw', 'ccw'))
longformdata_rotation$Single_triplet_learning = longformdata_rotation$STL
longformdata_rotation$Retention_to_next_triplet = longformdata_rotation$RetainedRaw
longformdata_rotation$retRatio = longformdata_rotation$Retained

xtmdf_rotation = makeMedianXTmDf(longformdata_rotation)
df2_rotation = data_summary(xtmdf_rotation, varname='perf', groupnames=c('mvt', 'trialsOut'))


## now check retention values
retratiodf_rotation = getSubjectwiseRetRatios(longformdata_rotation)
retratiodf2_rotation = data_summary(retratiodf_rotation, varname='medvals', groupnames=c('mvt'))
retRatioTtests_RotationInLab = retentionRatioTTests(retratiodf_rotation)
retRatioTtests_RotationInLab


## now split values by whether or not adaptation was in the correct direction
rightwaySTLdf_rotation = getRightwayMeans(longformdata_rotation)
cor.test(rightwaySTLdf_rotation$means_go,
         rightwaySTLdf_rotation$means_stop, method = "pearson")

wrongwaySTLdf_rotation = getWrongwayMeans(longformdata_rotation)
cor.test(wrongwaySTLdf_rotation$means_go,
         wrongwaySTLdf_rotation$means_stop, method = "pearson")



# set up plots for Figure 2
fig2_ANOVAplot_go <- ggplot(anovasummarydf[anovasummarydf$Movement == 'go',],
                       aes(x=Rotation, y=STL, group=Movement, color=Movement)) +
  geom_hline(yintercept=0, linetype='dashed', color='black') +
  geom_pointrange(aes(ymin=STL-sem, ymax=STL+sem))+
  theme(legend.position = c(0.8, 0.8)) +
  ylim(-6,6)+
  ylab('Single-Trial Learning (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

fig2_ANOVAplot_stop <- ggplot(anovasummarydf[anovasummarydf$Movement == 'stop',],
                       aes(x=Rotation, y=STL, group=Movement, color=Movement)) + 
  geom_hline(yintercept=0, linetype='dashed', color='black') +
  geom_pointrange(aes(ymin=STL-sem, ymax=STL+sem))+
  theme(legend.position = c(0.8, 0.8)) +
  ylim(-6,6)+
  ylab('Single-Trial Learning (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

fig2_perfXTrials <- ggplot(df2_rotation, aes(x=trialsOut, y=perf, group=mvt, color=mvt)) + 
  geom_line()+
  geom_pointrange(aes(ymin=perf-sem, ymax=perf+sem))+
  theme(legend.position = c(0.25, 0.8)) +
  ylim(0,4)+
  ylab('delta Hand Angle (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

fig2_retRatioPlot <- ggplot(retratiodf2_rotation, aes(x=mvt, y=medvals, color=mvt)) +
  geom_pointrange(aes(ymin=medvals-sem, ymax=medvals+sem))+
  scale_y_continuous(breaks=seq(0,1,by=0.25))+
  expand_limits(y = c(0, 1)) +
  theme(legend.position = c(0.8, 0.8)) +
  ylab('Remembered STL (ratio)') +
  xlab('Movement') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

fig2_rreg<-lm(means_stop ~ means_go, data = rightwaySTLdf_rotation)
fig2_rightWayCorrPlot <- ggplot(rightwaySTLdf_rotation, aes(x = means_go, y = means_stop)) +
  geom_point() +
  geom_abline(intercept = fig2_rreg$coefficients[1],
              slope = fig2_rreg$coefficients[2], 
              linetype="dashed") +
  scale_y_continuous(breaks=seq(0,8,by=2))+
  scale_x_continuous(breaks=seq(0,8,by=2))+
  expand_limits(y = c(0, 8), x = c(0, 8)) +
  ylab('Right-Way STL on Stop Trials') +
  xlab('Right-Way STL on Go Trials') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  coord_fixed()

fig2_wreg<-lm(means_stop ~ means_go, data = wrongwaySTLdf_rotation)
fig2_wrongWayCorrPlot <- ggplot(wrongwaySTLdf_rotation, aes(x = means_go, y = means_stop)) +
  geom_point() +
  geom_abline(intercept = fig2_wreg$coefficients[1],
              slope = fig2_wreg$coefficients[2], 
              linetype="dashed") +
  scale_y_continuous(breaks=seq(-6,0,by=2))+
  scale_x_continuous(breaks=seq(-6,0,by=2))+
  expand_limits(y = c(-6, 0), x = c(-6, 0)) +
  ylab('Wrong-Way STL on Stop Trials') +
  xlab('Wrong-Way STL on Go Trials') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  coord_fixed()

blankPlot <- ggplot()+geom_blank(aes(1,1))


# first row of panels
fig2_upperRow <- ggarrange(blankPlot, fig2_ANOVAplot_go, fig2_ANOVAplot_stop, fig2_perfXTrials,
          fig2_retRatioPlot,
          ncol = 5, nrow = 1, widths = c(0.5, 0.5, 0.5, 0.5, 0.3))
ggsave("fig2_upperRow.pdf",fig2_upperRow)
fig2_upperRow

# second row of panels
fig2_lowerRow <- ggarrange(blankPlot, fig2_rightWayCorrPlot, fig2_wrongWayCorrPlot,
          ncol = 3, nrow = 1, widths = c(0.4, 0.6, 0.6))
ggsave("fig2_lowerRow.pdf",fig2_lowerRow)
fig2_lowerRow


#### -------------------- R O T A T I O N ----------------------------------####
rotnOnlineData = tripletdata[tripletdata$Experiment=='Rotation_Online'
                       & tripletdata$includeSubj == 1
                       & tripletdata$includeTrip == 1,]
rotnOnlineData = organizeFactors(rotnOnlineData)

meandf_rotnOnline = splitRotnMvtMeans(rotnOnlineData)
longfmt_rotnOnline = meanDfToLongForm(meandf_rotnOnline)

rotnOnline_crossedRE = lmerTest::lmer(
  Single_triplet_learning ~ 1 + Rotation * Movement
  + (1 | Subject_ID),
  data = rotnOnlineData
)
rotn.stl.aov<-anova(rotnOnline_crossedRE, type=3, ddf="Kenward-Roger")
rotn.stl.aov
r2beta(rotnOnline_crossedRE, method = 'kr')

# post hoc tests
rotnOnline.emm.s <- emmeans(rotnOnline_crossedRE, c("Movement", "Rotation"))

# specify emmeans to test
gocw_lookup = c(1, 0, 0, 0, 0, 0) # lookup vector for first emmean in grid
stopcw_lookup = c(0, 1, 0, 0, 0, 0)
go0_lookup = c(0, 0, 1, 0, 0, 0)
stop0_lookup = c(0, 0, 0, 1, 0, 0)
goccw_lookup = c(0, 0, 0, 0, 1, 0)
stopccw_lookup = c(0, 0, 0, 0, 0, 1)
contrastOut <- summary(contrast(rotnOnline.emm.s, method = list("go cw - stop cw" = gocw_lookup - stopcw_lookup,
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
effectSizeOut <- summary(eff_size(rotnOnline.emm.s,
                                  sigma = sigma(rotnOnline_crossedRE),
                                  edf = df.residual(rotnOnline_crossedRE))) # returns cohen's d: https://cran.r-project.org/web/packages/emmeans/vignettes/
effectSizeOut


## collapse STL and remembered STL across rotation directions
xtmdf_rotnOnline = makeMedianXTmDf(rotnOnlineData)
df2_rotnOnline = data_summary(xtmdf_rotnOnline, varname='perf', groupnames=c('mvt', 'trialsOut'))

## now check retention values
retratiodf_rotnOnline = getSubjectwiseRetRatios(rotnOnlineData)
retratiodf2_rotnOnline = data_summary(retratiodf_rotnOnline, varname='medvals', groupnames=c('mvt'))
retRatioTtests_RotnOnline = retentionRatioTTests(retratiodf_rotnOnline)
retRatioTtests_RotnOnline


## now split values by whether or not adaptation was in the correct direction
rightwaySTLdf_rotnOnline = getRightwayMeans(rotnOnlineData)
cor.test(rightwaySTLdf_rotnOnline$means_go,
         rightwaySTLdf_rotnOnline$means_stop, method = "pearson")

wrongwaySTLdf_rotnOnline = getWrongwayMeans(rotnOnlineData)
cor.test(wrongwaySTLdf_rotnOnline$means_go,
         wrongwaySTLdf_rotnOnline$means_stop, method = "pearson")


# PLOTTING - SUPPLEMENTARY FIGURE 2
SFig2_boxplot_rotn <- ggplot(longfmt_rotnOnline, aes(x=rotn, y=meanvals, colour=mvt)) +
  geom_hline(yintercept=0, linetype='dashed', color='black') +
  geom_boxplot(outlier.shape = NA, notch=TRUE, width=0.5) +
  scale_y_continuous(breaks=seq(-10,10,by=5))+
  expand_limits(y = c(-10, 10)) +
  theme(legend.position = c(0.9, 0.9)) +
  ylab('STL (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

SFig2_emmsplot_rotn_go <- ggplot(as.data.frame(rotnOnline.emm.s)[as.data.frame(rotnOnline.emm.s)$Movement == 'go',],
                           aes(x = Rotation, y = emmean, color = Movement)) +
  geom_hline(yintercept=0, linetype='dashed', color='black') +
  geom_pointrange(aes(ymin = emmean - (emmean-lower.CL), ymax = emmean + (upper.CL-emmean))) +
  ylim(-6, 6) +
  theme(legend.position = c(0.8, 0.8)) +
  ylab('STL (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

SFig2_emmsplot_rotn_stop <- ggplot(as.data.frame(rotnOnline.emm.s)[as.data.frame(rotnOnline.emm.s)$Movement == 'stop',],
                             aes(x = Rotation, y = emmean, color = Movement)) +
  geom_hline(yintercept=0, linetype='dashed', color='black') +
  geom_pointrange(aes(ymin = emmean - (emmean-lower.CL), ymax = emmean + (upper.CL-emmean))) +
  ylim(-6, 6) +
  theme(legend.position = c(0.8, 0.8)) +
  ylab('STL (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

SFig2_perfXTrials <- ggplot(df2_rotnOnline, aes(x=trialsOut, y=perf, group=mvt, color=mvt)) + 
  geom_line()+
  geom_pointrange(aes(ymin=perf-sem, ymax=perf+sem))+
  theme(legend.position = c(0.25, 0.8)) +
  ylim(0,4)+
  ylab('delta Hand Angle (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

SFig2_retRatioPlot <- ggplot(retratiodf2_rotnOnline, aes(x=mvt, y=medvals, color=mvt)) +
  geom_pointrange(aes(ymin=medvals-sem, ymax=medvals+sem))+
  scale_y_continuous(breaks=seq(0,1,by=0.25))+
  expand_limits(y = c(0, 1)) +
  theme(legend.position = c(0.8, 0.8)) +
  ylab('Remembered STL/STL') +
  xlab('Movement') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

rreg<-lm(means_stop ~ means_go, data = rightwaySTLdf_rotnOnline)
SFig2_rightWayCorrPlot <- ggplot(rightwaySTLdf_rotnOnline, aes(x = means_go, y = means_stop)) +
  geom_point() +
  geom_abline(intercept = rreg$coefficients[1],
              slope = rreg$coefficients[2], 
              linetype="dashed") +
  scale_y_continuous(breaks=seq(2,12,by=2.5))+
  scale_x_continuous(breaks=seq(2,12,by=2.5))+
  expand_limits(y = c(2, 12), x = c(2, 12)) +
  ylab('Right-Way STL on Stop Trials') +
  xlab('Right-Way STL on Go Trials') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  coord_fixed()

wreg<-lm(means_stop ~ means_go, data = wrongwaySTLdf_rotnOnline)
SFig2_wrongWayCorrPlot <- ggplot(wrongwaySTLdf_rotnOnline, aes(x = means_go, y = means_stop)) +
  geom_point() +
  geom_abline(intercept = wreg$coefficients[1],
              slope = wreg$coefficients[2], 
              linetype="dashed") +
  scale_y_continuous(breaks=seq(-12,0,by=3))+
  scale_x_continuous(breaks=seq(-12,0,by=3))+
  expand_limits(y = c(-12, 0), x = c(-12, 0)) +
  ylab('Wrong-Way STL on Stop Trials') +
  xlab('Wrong-Way STL on Go Trials') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  coord_fixed()

# first row of panels
sfig2_upperRow <- ggarrange(blankPlot, SFig2_boxplot_rotn, SFig2_emmsplot_rotn_go, SFig2_emmsplot_rotn_stop,
          ncol = 4, nrow = 1, widths = c(0.6, 0.8, 0.5, 0.5))
ggsave("sfig2_upperRow.pdf",sfig2_upperRow)
sfig2_upperRow
# second row of panels
sfig2_lowerRow <- ggarrange(SFig2_perfXTrials, SFig2_retRatioPlot, SFig2_rightWayCorrPlot, SFig2_wrongWayCorrPlot,
          ncol = 4, nrow = 1, widths = c(0.5, 0.3, 0.6, 0.6))
ggsave("sfig2_lowerRow.pdf",sfig2_lowerRow)
sfig2_lowerRow


#### ------------------------- C L A M P  ----------------------------------####
clampdata = tripletdata[tripletdata$Experiment=='Clamp_Online'
                       & tripletdata$includeSubj == 1
                       & tripletdata$includeTrip == 1,]
clampdata = organizeFactors(clampdata)

meandf_clampOnline = splitRotnMvtMeans(clampdata)
longfmt_clampOnline = meanDfToLongForm(meandf_clampOnline)

clamp_crossedRE = lmerTest::lmer(
  Single_triplet_learning ~ 1 + Rotation * Movement
  + (1 | Subject_ID),
  data = clampdata
)
clamp.stl.aov<-anova(clamp_crossedRE, type=3, ddf="Kenward-Roger")
clamp.stl.aov
r2beta(clamp_crossedRE, method = 'kr')

# post hoc tests
clamp.emm.s <- emmeans(clamp_crossedRE, c("Movement", "Rotation"))
contrastOut <- summary(contrast(clamp.emm.s, method = list("go cw - stop cw" = gocw_lookup - stopcw_lookup,
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
effectSizeOut <- summary(eff_size(clamp.emm.s, sigma = sigma(clamp_crossedRE),
                                  edf = df.residual(clamp_crossedRE)))

## collapse STL and remembered STL across rotation directions
xtmdf_clamp = makeMedianXTmDf(clampdata)
df2_clamp = data_summary(xtmdf_clamp, varname='perf', groupnames=c('mvt', 'trialsOut'))

## now check retention values
retratiodf_clamp = getSubjectwiseRetRatios(clampdata)
retratiodf2_clamp = data_summary(retratiodf_clamp, varname='medvals', groupnames=c('mvt'))
retRatioTtests_Clamp = retentionRatioTTests(retratiodf_clamp)
retRatioTtests_Clamp

## now split values by whether or not adaptation was in the correct direction
rightwaySTLdf_clamp = getRightwayMeans(clampdata)
cor.test(rightwaySTLdf_clamp$means_go, rightwaySTLdf_clamp$means_stop, method = "pearson")

wrongwaySTLdf_clamp = getWrongwayMeans(clampdata)
cor.test(wrongwaySTLdf_clamp$means_go, wrongwaySTLdf_clamp$means_stop, method = "pearson")


# PLOTTING - SUPPLEMENTARY FIGURE 4
sfig4_boxplot_clamp <- ggplot(longfmt_clampOnline, aes(x=rotn, y=meanvals, colour=mvt)) +
  geom_hline(yintercept=0, linetype='dashed', color='black') +
  geom_boxplot(outlier.shape = NA, notch=TRUE, width=0.5) +
  scale_y_continuous(breaks=seq(-10,10,by=5), limits=c(-10,10))+
  expand_limits(y = c(-10, 10)) +
  theme(legend.position = c(0.9, 0.9)) +
  ylab('STL (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

sfig4_emmsplot_clamp_go <- ggplot(as.data.frame(clamp.emm.s)[as.data.frame(clamp.emm.s)$Movement == 'go',],
                           aes(x = Rotation, y = emmean, color = Movement)) +
  geom_hline(yintercept=0, linetype='dashed', color='black') +
  geom_pointrange(aes(ymin = emmean - (emmean-lower.CL), ymax = emmean + (upper.CL-emmean))) +
  ylim(-6, 6) +
  theme(legend.position = c(0.8, 0.8)) +
  ylab('STL (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
sfig4_emmsplot_clamp_stop <- ggplot(as.data.frame(clamp.emm.s)[as.data.frame(clamp.emm.s)$Movement == 'stop',],
                             aes(x = Rotation, y = emmean, color = Movement)) +
  geom_hline(yintercept=0, linetype='dashed', color='black') +
  geom_pointrange(aes(ymin = emmean - (emmean-lower.CL), ymax = emmean + (upper.CL-emmean))) +
  ylim(-6, 6) +
  theme(legend.position = c(0.8, 0.8)) +
  ylab('STL (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

sfig4_perfXTrials <- ggplot(df2_clamp, aes(x=trialsOut, y=perf, group=mvt, color=mvt)) + 
  geom_line()+
  geom_pointrange(aes(ymin=perf-sem, ymax=perf+sem))+
  theme(legend.position = c(0.25, 0.8)) +
  ylim(0,4)+
  ylab('delta Hand Angle (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

sfig4_retRatioPlot <- ggplot(retratiodf2_clamp, aes(x=mvt, y=medvals, color=mvt)) +
  geom_pointrange(aes(ymin=medvals-sem, ymax=medvals+sem))+
  scale_y_continuous(breaks=seq(0,1,by=0.25))+
  expand_limits(y = c(0, 1)) +
  theme(legend.position = c(0.8, 0.8)) +
  ylab('Remembered STL/STL') +
  xlab('Movement') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

rreg<-lm(means_stop ~ means_go, data = rightwaySTLdf_clamp)
sfig4_rightWayCorrPlot <- ggplot(rightwaySTLdf_clamp, aes(x = means_go, y = means_stop)) +
  geom_point() +
  geom_abline(intercept = rreg$coefficients[1],
              slope = rreg$coefficients[2], 
              linetype="dashed") +
  scale_y_continuous(breaks=seq(0,15,by=5))+
  scale_x_continuous(breaks=seq(0,15,by=5))+
  expand_limits(y = c(0, 15), x = c(0, 15)) +
  ylab('Right-Way STL on Stop Trials') +
  xlab('Right-Way STL on Go Trials') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  coord_fixed()

wreg<-lm(means_stop ~ means_go, data = wrongwaySTLdf_clamp)
sfig4_wrongWayCorrPlot <- ggplot(wrongwaySTLdf_clamp, aes(x = means_go, y = means_stop)) +
  geom_point() +
  geom_abline(intercept = wreg$coefficients[1],
              slope = wreg$coefficients[2], 
              linetype="dashed") +
  scale_y_continuous(breaks=seq(-20,0,by=5))+
  scale_x_continuous(breaks=seq(-20,0,by=5))+
  expand_limits(y = c(-20, 0), x = c(-20, 0)) +
  ylab('Wrong-Way STL on Stop Trials') +
  xlab('Wrong-Way STL on Go Trials') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  coord_fixed()

# first row of panels
sfig4_upperRow<-ggarrange(blankPlot, sfig4_boxplot_clamp,
                          sfig4_emmsplot_clamp_go, sfig4_emmsplot_clamp_stop,
          ncol = 4, nrow = 1, widths = c(0.6, 0.8, 0.5, 0.5))
ggsave("sfig4_upperRow.pdf",sfig4_upperRow)
sfig4_upperRow
# second row of panels
sfig4_lowerRow <- ggarrange(sfig4_perfXTrials, sfig4_retRatioPlot,
                            sfig4_rightWayCorrPlot, sfig4_wrongWayCorrPlot,
          ncol = 4, nrow = 1, widths = c(0.5, 0.3, 0.6, 0.6))
ggsave("sfig4_lowerRow.pdf",sfig4_lowerRow)
sfig4_lowerRow

# ---------------------------------- ROTATION 0 GO ---------------------------#
rot0godata = tripletdata[tripletdata$Experiment=='Rotation_Online_0go'
                        & tripletdata$includeSubj == 1
                        & tripletdata$includeTrip == 1,]
rot0godata = organizeFactors(rot0godata)

meandf_rot0goOnline = splitRotnMvtMeans(rot0godata)
longfmt_rot0goOnline = meanDfToLongForm(meandf_rot0goOnline)

rot0go_model = lmerTest::lmer(
  Single_triplet_learning ~ 1 + Rotation
  + (1 | Subject_ID),
  data = rot0godata[rot0godata$Movement == 'stop',]
)
rot0go.stl.aov<-anova(rot0go_model, type=3, ddf="Kenward-Roger")
rot0go.stl.aov
r2beta(rot0go_model, method = 'kr')

# post hoc tests
rot0go.emm.s <- emmeans(rot0go_model, c("Rotation"))
zero_lookup = c(0, 1, 0) # lookup vector for first emmean in grid
cw_lookup = c(1, 0)
ccw_lookup = c(0, 1)
contrastOut <- summary(contrast(rot0go.emm.s, method = list("cw - ccw" = cw_lookup - ccw_lookup),
                                adjust = "fdr"))
effectSizeOut <- summary(eff_size(rot0go.emm.s, sigma = sigma(rot0go_model),
                                  edf = df.residual(rot0go_model)))
contrastOut
effectSizeOut

## collapse STL and remembered STL across rotation directions
xtmdf_rot0go = makeMedianXTmDf(rot0godata)
df2_rot0go = data_summary(xtmdf_rot0go, varname='perf', groupnames=c('mvt', 'trialsOut'))

## now check retention values
retratiodf_rot0go = getSubjectwiseRetRatios(rot0godata)
retratiodf2_rot0go = data_summary(retratiodf_rot0go, varname='medvals', groupnames=c('mvt'))
shapiro.test(retratiodf_rot0go$medvals[retratiodf_rot0go$mvt=='stop'])
t.test(retratiodf_rot0go$medvals[retratiodf_rot0go$mvt=='stop'], mu=0)
d = (0 - mean(retratiodf_rot0go$medvals[retratiodf_rot0go$mvt=='stop' & 
                                          is.na(
                                            retratiodf_rot0go$medvals)==FALSE])
     )/sd(retratiodf_rot0go$medvals[retratiodf_rot0go$mvt=='stop' & 
                                      is.na(
                                        retratiodf_rot0go$medvals)==FALSE])

# PLOTTING - FIGURE 3, TOP PANEL
fig3_boxplot_rot0go <- ggplot(longfmt_rot0goOnline, aes(x=rotn, y=meanvals, colour=mvt)) +
  geom_hline(yintercept=0, linetype='dashed', color='black') +
  geom_boxplot(outlier.shape = NA, notch=TRUE, width=0.5) +
  scale_y_continuous(breaks=seq(-6,6,by=4), limits = c(-6,6))+
  expand_limits(y = c(-6, 6)) +
  theme(legend.position = c(0.9, 0.9)) +
  ylab('STL (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

fig3_emmsplot_rot0go_stop <- ggplot(as.data.frame(rot0go.emm.s),
                              aes(x = Rotation, y = emmean)) +
  geom_hline(yintercept=0, linetype='dashed', color='black') +
  geom_pointrange(aes(ymin = emmean - (emmean-lower.CL), ymax = emmean + (upper.CL-emmean))) +
  ylim(-4, 4) +
  theme(legend.position = c(0.8, 0.8)) +
  ylab('STL (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

fig3_perfXTrials <- ggplot(df2_rot0go[df2_rot0go$mvt=='stop',], aes(x=trialsOut, y=perf)) + 
  geom_line()+
  geom_pointrange(aes(ymin=perf-sem, ymax=perf+sem))+
  theme(legend.position = c(0.25, 0.8)) +
  ylim(-0.1,4)+
  ylab('delta Hand Angle (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

fig3_retRatioPlot <- ggplot(retratiodf2_rot0go[retratiodf2_rot0go$mvt=='stop',], aes(x=mvt, y=medvals)) +
  geom_pointrange(aes(ymin=medvals-sem, ymax=medvals+sem))+
  scale_y_continuous(breaks=seq(0,1,by=0.25))+
  expand_limits(y = c(0, 1)) +
  theme(legend.position = c(0.8, 0.8)) +
  ylab('Remembered STL/STL') +
  xlab('Movement') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

# first row of panels
fig3_upperRow <- ggarrange(blankPlot, fig3_boxplot_rot0go, fig3_emmsplot_rot0go_stop,
          fig3_perfXTrials, fig3_retRatioPlot,
          ncol = 5, nrow = 1, widths = c(0.6, 0.5, 0.5, 0.3, 0.2))
ggsave("fig3_upperRow.pdf",fig3_upperRow)
fig3_upperRow


# ---------------------------------- CLAMP 0 GO ---------------------------#
clamp0godata = tripletdata[tripletdata$Experiment=='Clamp_Online_0go'
                         & tripletdata$includeSubj == 1
                         & tripletdata$includeTrip == 1,]
clamp0godata = organizeFactors(clamp0godata)

meandf_clamp0goOnline = splitRotnMvtMeans(clamp0godata)
longfmt_clamp0goOnline = meanDfToLongForm(meandf_clamp0goOnline)

clamp0go_model = lmerTest::lmer(
  Single_triplet_learning ~ 1 + Rotation
  + (1 | Subject_ID),
  data = clamp0godata[clamp0godata$Movement == 'stop',]
)
clamp0go.stl.aov<-anova(clamp0go_model, type=3, ddf="Kenward-Roger")
clamp0go.stl.aov
r2beta(clamp0go_model, method = 'kr')

# post hoc tests
clamp0go.emm.s <- emmeans(clamp0go_model, c("Rotation"))
zero_lookup = c(0, 1, 0) # lookup vector for first emmean in grid
cw_lookup = c(1, 0, 0)
ccw_lookup = c(0, 0, 1)
contrastOut <- summary(contrast(clamp0go.emm.s, method = list("cw - ccw" = cw_lookup - ccw_lookup,
                                                              "cw - 0" = cw_lookup - zero_lookup,
                                                              "0 - ccw" = zero_lookup - ccw_lookup),
                                adjust = "fdr"))
effectSizeOut <- summary(eff_size(clamp0go.emm.s, sigma = sigma(clamp0go_model),
                                  edf = df.residual(clamp0go_model)))
contrastOut
effectSizeOut

## collapse STL and remembered STL across rotation directions
xtmdf_clamp0go = makeMedianXTmDf(clamp0godata)
df2_clamp0go = data_summary(xtmdf_clamp0go, varname='perf', groupnames=c('mvt', 'trialsOut'))

## now check retention values
retratiodf_clamp0go = getSubjectwiseRetRatios(clamp0godata)
retratiodf2_clamp0go = data_summary(retratiodf_clamp0go, varname='medvals', groupnames=c('mvt'))
shapiro.test(retratiodf_clamp0go$medvals[retratiodf_clamp0go$mvt=='stop'])
rstatix::wilcox_test(formula = medvals ~ 1,
                    data = retratiodf_clamp0go[retratiodf_clamp0go$mvt=='stop',])
rstatix::wilcox_effsize(formula = medvals ~ 1,
                     data = retratiodf_clamp0go[retratiodf_clamp0go$mvt=='stop',])
quantile(retratiodf_clamp0go$medvals[retratiodf_clamp0go$mvt=='stop'],
         c(0.25, 0.5, 0.75), na.rm = TRUE)



# PLOTTING - FIGURE 3, BOTTOM ROW
fig3_boxplot_clamp0go <- ggplot(longfmt_clamp0goOnline, aes(x=rotn, y=meanvals)) +
  geom_hline(yintercept=0, linetype='dashed', color='black') +
  geom_boxplot(outlier.shape = NA, notch=TRUE, width=0.5) +
  scale_y_continuous(breaks=seq(-8,8,by=4))+
  expand_limits(y = c(-8, 8)) +
  theme(legend.position = c(0.9, 0.9)) +
  ylab('STL (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

fig3_emmsplot_clamp0go_stop <- ggplot(as.data.frame(clamp0go.emm.s),
                                    aes(x = Rotation, y = emmean)) +
  geom_hline(yintercept=0, linetype='dashed', color='black') +
  geom_pointrange(aes(ymin = emmean - (emmean-lower.CL), ymax = emmean + (upper.CL-emmean))) +
  ylim(-4, 4) +
  theme(legend.position = c(0.8, 0.8)) +
  ylab('STL (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

fig3_perfXTrials_clamp0go <- ggplot(df2_clamp0go[df2_clamp0go$mvt=='stop',], aes(x=trialsOut, y=perf)) + 
  geom_line()+
  geom_pointrange(aes(ymin=perf-sem, ymax=perf+sem))+
  theme(legend.position = c(0.25, 0.8)) +
  ylim(-0.1,4)+
  ylab('delta Hand Angle (degrees)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

fig3_retRatioPlot_clamp0go <- ggplot(retratiodf2_clamp0go[retratiodf2_clamp0go$mvt=='stop',], aes(x=mvt, y=medvals)) +
  geom_pointrange(aes(ymin=medvals-sem, ymax=medvals+sem))+
  scale_y_continuous(breaks=seq(0,1,by=0.25))+
  expand_limits(y = c(0, 1)) +
  theme(legend.position = c(0.8, 0.8)) +
  ylab('Remembered STL/STL') +
  xlab('Movement') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

# first row of panels
fig3_bottomRow <- ggarrange(blankPlot, fig3_boxplot_clamp0go, fig3_emmsplot_clamp0go_stop,
          fig3_perfXTrials_clamp0go, fig3_retRatioPlot_clamp0go,
          ncol = 5, nrow = 1, widths = c(0.6, 0.5, 0.5, 0.3, 0.2))
ggsave("fig3_bottomRow.pdf",fig3_bottomRow)
fig3_bottomRow