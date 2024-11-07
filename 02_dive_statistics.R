# Script 002: Dive statistics
# - Calculate simple statistics for parameters in tables
# - Make additional data plots
# - Calculate correlation coefficients etc

############ User defined ############
ptts <-c(134663,134664,134666,134667,134668,134669,134670,161587,161588,161590,161591,161592,161593)
ptts_w_data = ptts[-3] # One ptt does not have any dive data (134666)

############ Table 1: deployment overview ############
load("03_data_tables.Rdata")
load("tracks.Rdata")

### Duration: Time interval from tag-on to last BL data point
totaldur = 0 
for(i in 1:length(ptts)){ 
  ttagon = tracks$date_time[which(tracks$ptt==ptts[i])[1]]
  tlast = max(btab$End[btab$Ptt==ptts[i]])
  dtdays = floor(difftime(tlast,ttagon,units="days"))
  dthours = round((difftime(tlast,ttagon,units="hour") - dtdays)/3600)
  print(paste(dtdays,'d ',dthours,'h'))
  ttmp = difftime(tlast,ttagon,units="days")
  if(!is.infinite(ttmp)){
    totaldur = totaldur + difftime(tlast,ttagon,units="days") 
  }
}
print(totaldur) # total duration in days

############ Table 3: Dives to oceanographic zones ############

# Remove dives that are too shallow or during post-tagging and sonar periods
btab = btab[!btab$ishallow & !btab$itagging & !btab$isonar,]

tab3 = matrix(NA, 12, 10)
for(i in 1:length(ptts_w_data)){
  n.dives = sum(btab$What=="Dive" & btab$Ptt==ptts_w_data[i])
  n.sdives = sum(btab$What=="Dive" & btab$Ptt==ptts_w_data[i] & btab$cluster==1)
  n.ddives = sum(btab$What=="Dive" & btab$Ptt==ptts_w_data[i] & btab$cluster==2)
  tab3[i,1] = n.sdives
  tab3[i,2] = n.ddives  
  ind = btab$What=="Dive" & btab$Zone=="Epipelagic" & !btab$ibentic & btab$Ptt==ptts_w_data[i]
  tab3[i,3] = sum(ind & btab$cluster==1) / n.sdives
  tab3[i,4] = sum(ind & btab$cluster==2) / n.ddives
  ind = btab$What=="Dive" & btab$Zone=="Mesopelagic" & !btab$ibentic & btab$Ptt==ptts_w_data[i]
  tab3[i,5] = sum(ind & btab$cluster==1) / n.sdives
  tab3[i,6] = sum(ind & btab$cluster==2) / n.ddives
  ind = btab$What=="Dive" & btab$Zone=="Bathypelagic" & !btab$ibentic & btab$Ptt==ptts_w_data[i]
  tab3[i,7] = sum(ind & btab$cluster==1) / n.sdives
  tab3[i,8] = sum(ind & btab$cluster==2) / n.ddives
  ind = btab$What=="Dive" & btab$ibentic & btab$Ptt==ptts_w_data[i]
  tab3[i,9] = sum(ind & btab$cluster==1) / n.sdives
  tab3[i,10] = sum(ind & btab$cluster==2) / n.ddives
}
tab3[,3:10] = round(tab3[,3:10]*100)
tab3

# Mean and SD of max depth of benthic dives
mean(btab$DepthMean[btab$What=="Dive" & btab$ibentic]) # 1115 m
sd(btab$DepthMean[btab$What=="Dive" & btab$ibentic]) # 485 m


############ Longest and deepest dive ############
durmax = max(btab$DurationMean[btab$What=="Dive"], na.rm=T)/60 # 97.8333 min
btab$Start[which(btab$DurationMean==durmax*60)] # "2016-07-06 06:08:00 UTC"
btab$lat.pred[which(btab$DurationMean==durmax*60)] # 64.17344
btab$lon.pred[which(btab$DurationMean==durmax*60)] # -9.085513
btab$IDDI[which(btab$DurationMean==durmax*60)] # 88 min
depthmax = max(btab$DepthMean[btab$What=="Dive"], na.rm=T) # 2288 m
btab$Start[which(btab$DepthMean==depthmax)] # "2016-06-30 15:21:44 UTC" "2016-06-30 18:34:58 UTC"
btab$lat.pred[which(btab$DepthMean==depthmax)[1]] # 70.9433
btab$lon.pred[which(btab$DepthMean==depthmax)[1]] # -6.6149
btab$lat.pred[which(btab$DepthMean==depthmax)[2]] # 70.9776
btab$lon.pred[which(btab$DepthMean==depthmax)[2]] # -6.4692

# More info deepest dives
load("01_data_tables.Rdata")
btab$DepthMin[which(btab$DepthMean==depthmax)] # 2256 m
btab$DepthMax[which(btab$DepthMean==depthmax)] # 2319 m
btab$DurationMean[which(btab$DepthMean==depthmax)]/60 # 78 min, 73 min
difftime((as.POSIXct("2016-06-30 15:21:44 UTC") + (78*60)),"2016-06-30 18:34:58 UTC") # "IDDI"

# More info longest dive
btab$DurationMin[which(btab$DurationMean==durmax*60)]/60 # 97.82 min
btab$DurationMax[which(btab$DurationMean==durmax*60)]/60 # 97.85 min
btab$DepthMean[which(btab$DurationMean==durmax*60)] # 1231.5 m
btab$DepthMin[which(btab$DurationMean==durmax*60)] # 1216 m
btab$DepthMax[which(btab$DurationMean==durmax*60)] # 1247 m
btab$Bathy[which(btab$DurationMean==durmax*60)] # 1114 m

############ Table 2: Descriptive statistics of dive variables ############
load("03_data_tables.Rdata")

# Remove dives during tagging and sonar periods
btab = btab[!btab$itagging & !btab$isonar,]

# Replace dives that are too shallow by surface periods
btab$What[btab$ishallow]="Surface"

# Add surface duration to previous surface duration if needed
for(i in which(btab$ishallow)){
    if(btab$What[i-1]=="Surface" & btab$Ptt[i-1]==btab$Ptt[i] & !is.na(btab$DurationMean[i])){
      btab$DurationMean[i-1] = btab$DurationMean[i-1] + btab$DurationMean[i]
      btab$DurationMean[i] = NA
    }
}

# Indices of the 2 subsets 2014-2015 and 2016
iold = btab$Datetime<"2016-01-01"
inew = btab$Datetime>"2016-01-01"

### Number of dives
n = matrix(NA, 3, 2)
n[1,1] = nrow(btab[iold & btab$What=="Dive" & btab$cluster==1,]) # n=1567 dives
n[1,2] = nrow(btab[iold & btab$What=="Dive" & btab$cluster==2,]) # n=648 dives
n[2,1] = nrow(btab[inew & btab$What=="Dive" & btab$cluster==1,]) # n=4345 dives
n[2,2] = nrow(btab[inew & btab$What=="Dive" & btab$cluster==2,]) # n=1812 dives
n[3,1] = sum(n[1:2,1]) # n=5912 dives
n[3,2] = sum(n[1:2,2]) # n=2460 dives
n

### Dive rate 2016: number of shallow and deep dives / data hours
diverate = n[2,] / sum(btab$DurationMean[inew],na.rm=T)/3600

### Maximum dive depth (mean,sd,min,max)
tab2 = matrix(NA, 6, 4)
i1 = btab$DepthMean[iold & btab$What=="Dive" & btab$cluster==1]
i2 = btab$DepthMean[iold & btab$What=="Dive" & btab$cluster==2]
tab2[1, ] = c(mean(i1),sd(i1),mean(i2),sd(i2))
tab2[2, ] = c(min(i1),max(i1),min(i2),max(i2))
i1 = btab$DepthMean[inew & btab$What=="Dive" & btab$cluster==1]
i2 = btab$DepthMean[inew & btab$What=="Dive" & btab$cluster==2]
tab2[3, ] = c(mean(i1),sd(i1),mean(i2),sd(i2))
tab2[4, ] = c(min(i1),max(i1),min(i2),max(i2))
i1 = btab$DepthMean[btab$What=="Dive" & btab$cluster==1]
i2 = btab$DepthMean[btab$What=="Dive" & btab$cluster==2]
tab2[5, ] = c(mean(i1),sd(i1),mean(i2),sd(i2))
tab2[6, ] = c(min(i1),max(i1),min(i2),max(i2))
round(tab2)

# Mature male truncated due to tag's limited depth range
ind = btab$What=="Dive" & btab$Ptt=="134668" & btab$cluster==2 
length(btab$DepthMean[ind])  # n=423 long/deep dives
table(btab$DepthMean[ind])   # n=21 (5%) of those to maximum depth of 1847.75 m
iben = btab$What=="Dive" & btab$Ptt=="134668" & btab$ibentic
length(btab$DepthMean[iben]) # n=32 of (long/deep) dives are benthic dives
table(btab$DepthMean[iben])  # n=5 (1%) of those to maximum depth of 1847.75 m 

### Dive duration (mean,sd,min,max) in 2016
tab2 = matrix(NA, 2, 4)
i1 = btab$DurationMean[inew & btab$What=="Dive" & btab$cluster==1]/60
i2 = btab$DurationMean[inew & btab$What=="Dive" & btab$cluster==2]/60
tab2[1, ] = c(mean(i1),sd(i1),mean(i2),sd(i2))
tab2[2, ] = c(min(i1),max(i1),min(i2),max(i2))
round(tab2)

### Surface duration (mean,sd,min,max) in 2016
tab2 = matrix(NA, 2, 2)
i1 = btab$DurationMean[inew & btab$What=="Surface"]/60
tab2[1, ] = c(mean(i1,na.rm=T),sd(i1,na.rm=T))
tab2[2, ] = c(min(i1,na.rm=T),max(i1,na.rm=T))
round(tab2*10)/10
sum(inew & btab$What=="Surface" & !is.na(btab$DurationMean)) # n in 2016

### IDDI (mean,sd,min,max)
tab2 = matrix(NA, 6, 2)
i1 = btab$IDDI[iold & !is.na(btab$IDDI)]
tab2[1, ] = c(mean(i1,na.rm=T),sd(i1,na.rm=T))
#tab2[2, ] = c(min(i1,na.rm=T),max(i1,na.rm=T))
i1 = btab$IDDI[inew & !is.na(btab$IDDI)]
tab2[3, ] = c(mean(i1,na.rm=T),sd(i1,na.rm=T))
tab2[4, ] = c(min(i1,na.rm=T),max(i1,na.rm=T))
i1 = btab$IDDI[!is.na(btab$IDDI)]
tab2[5, ] = c(mean(i1,na.rm=T),sd(i1,na.rm=T))
tab2[6, ] = c(min(i1,na.rm=T),max(i1,na.rm=T))
round(tab2)

### % of time in each behavior
# Load hours in short/shallow, long/deep and surface, data gaps
load("05_time_spent_in_clusters.Rdata")

# Add hours within clusters across tags
dur.surf = 0
dur.clus1 = 0
dur.clus2 = 0
for(i in c(7:12)){ # 2016 only
  dur.surf = dur.surf + sum(dur_list[[i]][,1])
  dur.clus1 = dur.clus1 + sum(dur_list[[i]][,2])
  dur.clus2 = dur.clus2 + sum(dur_list[[i]][,3]) 
}
dur.tot = sum(dur.surf,dur.clus1,dur.clus2) # 2514 h (close to 2513 h used for dive rates)
round(c(dur.surf,dur.clus1,dur.clus2)/dur.tot*100)/100

############ Transition probabilities within and between dive types ############
load("03_data_tables.Rdata")

# Remove dives during tagging and sonar periods
btab = btab[!btab$itagging & !btab$isonar,]

# Replace dives that are too shallow by surface periods
btab$What[btab$ishallow]="Surface"

# Add surface duration to previous surface duration
for(i in which(btab$ishallow)){
  if(btab$What[i-1]=="Surface" & btab$Ptt[i-1]==btab$Ptt[i] & !is.na(btab$DurationMean[i])){
    btab$DurationMean[i-1] = btab$DurationMean[i-1] + btab$DurationMean[i]
    btab$DurationMean[i] = NA
  }
}
# Split off surface periods
btab.surf = btab[btab$What=="Surface",]
btab = btab[btab$What=="Dive",]

# Probability of following deep dive
n.deepshallow = 0
n.deepdeep = 0 
for(i in 1:(nrow(btab)-1)){
  if(btab$Ptt[i]==btab$Ptt[i+1] & !btab$igap[i]){ # row is not new tag or data gap
    if(diff(btab$cluster[i:(i+1)])==-1){ # if dive type changed from 2 to 1
      n.deepshallow = n.deepshallow + 1
    }
    if(sum(btab$cluster[i:(i+1)])==4){ # if dive type stayed 2
      n.deepdeep = n.deepdeep + 1
    }
  }
}
n.deepshallow/(n.deepdeep+n.deepshallow) # P of deep followed by shallow

# Probability of following shallow dive
n.shallowdeep = 0
n.shallowshallow = 0 
for(i in 1:(nrow(btab)-1)){
  if(btab$Ptt[i]==btab$Ptt[i+1] & !btab$igap[i]){ # row is not new tag or data gap
    if(diff(btab$cluster[i:(i+1)])==1){ # if dive type changed from 1 to 2
      n.shallowdeep = n.shallowdeep + 1
    }
    if(sum(btab$cluster[i:(i+1)])==2){ # if dive type stayed 2
      n.shallowshallow = n.shallowshallow + 1
    }
  }
}
n.shallowshallow/(n.shallowshallow+n.shallowdeep) # P of shallow followed by shallow

############ Correlation tests ############
# Uses dives-only btab and btab.surf from above

# Max dive depth vs duration for 2016
cor.test(btab$DepthMean, btab$DurationMean, method = "spearman", exact=FALSE) 

# Max dive depth vs water depth of non-benthic deep dives
cor.test(btab$DepthMean[btab$cluster==2 & !btab$ibentic], btab$Bathy[btab$cluster==2 & !btab$ibentic], method = "spearman", exact=FALSE) 

# Plot of Max depth vs seafloor depth
ggplot(aes(Bathy, DepthMean, groupColour = T, groupFill = T), data=btab[btab$cluster==2 & !btab$ibentic,]) +
  geom_point(colour = "black", fill = "grey", size = 4, alpha = 0.1, shape = 21) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  xlab("Seafloor depth (m)") +  
  ylab("Maximum depth of long/deep dives (m)") +
  stat_smooth(method = "lm", formula = y ~ x, geom = "smooth", colour = "#111") 


# Assign IDDI to last dive instead of first dive
btab$IDDI2 <- NA
ideep = which(btab$cluster==2)
for(i in 1:nrow(btab)){
  if(!is.na(btab$IDDI[i])){
    btab$IDDI2[ideep[min(which(ideep>i))]] = btab$IDDI[i]
  }
}

# Max dive depth vs IDDI
i2016 = btab$Datetime>"2016-01-01"
cor.test(btab$DepthMean, btab$IDDI, method = "spearman", exact=FALSE)
cor.test(btab$DepthMean[!i2016], btab$IDDI[!i2016], method = "spearman", exact=FALSE)
cor.test(btab$DepthMean[i2016], btab$IDDI[i2016], method = "spearman", exact=FALSE)
cor.test(btab$DepthMean, btab$IDDI2, method = "spearman", exact=FALSE)

# Dive duration vs IDDI (2016 only)
r12 = cor.test(btab$DurationMean, btab$IDDI, method = "spearman", exact=FALSE)
cor.test(btab$DurationMean2[!i2016], btab$IDDI[!i2016], method = "spearman", exact=FALSE) # imputed durations not correlated
r13 = cor.test(btab$DurationMean, btab$IDDI2, method = "spearman", exact=FALSE)
r23 = cor.test(btab$IDDI[i2016], btab$IDDI2[i2016], method = "spearman", exact=FALSE)
diffcor::diffcor.dep(r12$estimate, r13$estimate, r23$estimate, sum(!is.na(btab$DurationMean) & !is.na(btab$IDDI)), digit = 8)

# Dive duration vs IDDI (2016 only) - excluding shorter ones
i30 = btab$IDDI>=30
j30 = btab$IDDI2>=30
r12 = cor.test(btab$DurationMean[i30], btab$IDDI[i30], method = "spearman", exact=FALSE)
r13 = cor.test(btab$DurationMean[i30], btab$IDDI2[j30], method = "spearman", exact=FALSE)
r23 = cor.test(btab$IDDI[i2016 & i30], btab$IDDI2[i2016 & j30], method = "spearman", exact=FALSE)
diffcor::diffcor.dep(r12$estimate, r13$estimate, r23$estimate, sum(!is.na(btab$DurationMean) & !is.na(btab$IDDI)), digit = 8)


### Plots of Dive duration vs IDDI before and after dive
a <- btab[!is.na(btab$IDDI2),] %>% 
  ggplot(aes(IDDI2, DurationMean/60, groupColour = T, groupFill = T)) +
  geom_point(colour = "black", fill = "grey", size = 4, alpha = 0.1, shape = 21) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_x_continuous(trans='log10') +
  xlab("IDDI before dive (min)") +
  ylab("Dive duration (min)") +
  stat_smooth(method = "lm", formula = y ~ x, geom = "smooth", colour = "#111") 
b <- btab[!is.na(btab$IDDI),] %>% 
  ggplot(aes(IDDI, DurationMean/60, groupColour = T, groupFill = T)) +
  geom_point(colour = "black", fill = "grey", size = 4, alpha = 0.1, shape = 21) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_x_continuous(trans='log10') +
  xlab("IDDI after dive (min)") +
  ylab("Dive duration (min)") +
  stat_smooth(method = "lm", formula = y ~ x, geom = "smooth", colour = "#111") 
g <- ggpubr::ggarrange(a, b, ncol = 2)
print(g)



