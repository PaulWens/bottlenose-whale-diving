# Code for conducting the descriptive analysis of dive behaviour described in
# the manuscript "Long-term depth records of satellite-tagged northern 
# bottlenose whales reveal extraordinary dive capabilities" by Neubarth et al.

############ User defined ############
ptts <-c(134663,134664,134666,134667,134668,134669,134670,161587,161588,161590,161591,161592,161593)
ptts_w_data = ptts[-3] # One ptt does not have any dive data (134666)
ptts2016 <-c(161587,161588,161590,161592,161593,161591)
sstart2015 <- "2015-06-29 02:48:00" # Start of sonar exposure experiment 2015
sstart2016 <- "2016-06-18 12:16:00" # Start of sonar exposure experiment 2016


############ Load R packages ############
library(dplyr)
library(ggplot2)


############ Load data tables ############
# "btab" has the relevant fields in the original Behavior Log, vertical velocity, 
# data gaps, latitude and longitude, bathymetric depth, oceanographic zone.
# "stab" includes the Timeseries depth data with similar additional variables.
load("data_tables.Rdata")
load("tracks.Rdata")

############ Tracking durations in Table 1 ############
# Time interval from tag-on to last BL data point. Rerun after removal of 
# potential migration data to get tracking durations for IDs 134670 and 161593.
totaldur = 0 
for(i in 1:length(ptts)){ 
  ttagon = tracks$date_time[which(tracks$ptt==ptts[i])[1]]
  tlast = max(btab$End[btab$Ptt==ptts[i]])
  dtdays = floor(difftime(tlast,ttagon,units="days"))
  dthours = round((difftime(tlast,ttagon,units="hour") - dtdays)/3600)
  print(cbind(ptts[i],paste(dtdays,'d ',dthours,'h')))
  ttmp = difftime(tlast,ttagon,units="days")
  if(!is.infinite(ttmp)){
    totaldur = totaldur + difftime(tlast,ttagon,units="days") 
  }
}
print(totaldur) # total duration in days


############ Data filters and IDDI calculation ############
# Remove data for two whales above and south of the Iceland-Faroe Ridge. Whales
# likely switched to a migratory state.
ind = (btab$Ptt==134670 & btab$Datetime>"2015-07-05 02:28:50") | (btab$Ptt==161593 & btab$Datetime>"2016-07-06 12:39:00")
btab = btab[!ind,]

# Remove dives during tagging and sonar periods
btab = btab[!btab$itagging & !btab$isonar,]

# Relabel dives that are too shallow as surface periods
btab$What[btab$ishallow]="Surface"
btab$Zone[btab$ishallow]=""

# Add surface duration to previous surface duration if needed
for(i in which(btab$ishallow)){
  if(btab$What[i-1]=="Surface" & btab$Ptt[i-1]==btab$Ptt[i] & !is.na(btab$DurationMean[i])){
    btab$DurationMean[i-1] = btab$DurationMean[i-1] + btab$DurationMean[i]
    btab$DurationMean[i] = NA
  }
}

# Exclude erroneous dive durations from 2014-2015
ind = btab$Datetime < "2016-01-01 00:00:00"
btab$DurationMean[ind] = NA

### Calculate IDDI: Time elapsed between the end of a deep dive and start of the next one
btab$IDDI<-NA
ind = which((btab$Zone == "Mesopelagic" & !is.na(btab$Zone)) | (btab$Zone == "Bathypelagic" & !is.na(btab$Zone))) # indices of deep dives
iddi = -(as.numeric(btab$End[ind[1:(length(ind)-1)]] - btab$Start[ind[2:(length(ind))]]))/60
btab$IDDI[(btab$Zone == "Mesopelagic" & !is.na(btab$Zone)) | (btab$Zone == "Bathypelagic" & !is.na(btab$Zone))] = c(iddi,NA)

# Exclude IDDI if interval contains data gap or is last IDDI of deployment
for(i in (1:(length(ind)-1))){
  if(sum(btab$igap[(ind[i]+1):ind[i+1]])>0){
    btab$IDDI[ind[i]] = NA 
  }
  if(length(unique(btab$Ptt[(ind[i]+1):ind[i+1]]))>1){
    btab$IDDI[ind[i]] = NA
  }
}

# Exclude IDDI if interval overlapped with start of exposure
for(j in 1:length(ptts)){
  if(ptts[j]=="134668" | ptts[j]=="134670"){
    ilast_before_sonar = max(which(btab$Start[ind]<sstart2015 & btab$Ptt[ind]==ptts[j])) 
    btab$IDDI[ind[ilast_before_sonar]] = NA
  }
  if(ptts[j]=="161587" | ptts[j]=="161588" | ptts[j]=="161590"){
    ilast_before_sonar = max(which(btab$Start[ind]<sstart2016 & btab$Ptt[ind]==ptts[j])) 
    btab$IDDI[ind[ilast_before_sonar]] = NA  
  }
}


############ Summary statistics ############
# Note that all tables are ordered by PTT/tag ID but chronologically in the ms

### Table 2. No. of dives, and percentage benthic, and mean, sd, max dive depth
tab2 = data.frame(matrix(NA, 12, 9))
rownames(tab2) <- ptts_w_data
for(i in 1:length(ptts_w_data)){
  n.dives = sum(btab$What=="Dive" & btab$Ptt==ptts_w_data[i])
  n.edives = sum(btab$What=="Dive" & btab$Ptt==ptts_w_data[i] & btab$Zone=="Epipelagic")
  n.mdives = sum(btab$What=="Dive" & btab$Ptt==ptts_w_data[i] & btab$Zone=="Mesopelagic")
  n.bdives = sum(btab$What=="Dive" & btab$Ptt==ptts_w_data[i] & btab$Zone=="Bathypelagic")
  tab2[i,1] = n.edives
  tab2[i,2] = n.mdives 
  tab2[i,3] = n.bdives
  ind = btab$What=="Dive" & btab$ibentic & btab$Ptt==ptts_w_data[i]
  tab2[i,4] = sum(ind & btab$Zone=="Epipelagic") / n.edives
  tab2[i,5] = sum(ind & btab$Zone=="Mesopelagic") / n.mdives
  tab2[i,6] = sum(ind & btab$Zone=="Bathypelagic") / n.bdives
  tab2[i,7] = mean(btab$DepthMean[btab$What=="Dive" & btab$Ptt==ptts_w_data[i]])
  tab2[i,8] = sd(btab$DepthMean[btab$What=="Dive" & btab$Ptt==ptts_w_data[i]])
  tab2[i,9] = max(btab$DepthMean[btab$What=="Dive" & btab$Ptt==ptts_w_data[i]])
}
tab2[,4:6] = round(tab2[,4:6]*100)
tab2[,7:9] = round(tab2[,7:9])
tab2

### Max dive depth and duration per dive category
mean(btab$DepthMean[btab$What=="Dive" & btab$Zone == "Epipelagic"],na.rm=T) # 123 m
sd(btab$DepthMean[btab$What=="Dive" & btab$Zone == "Epipelagic"],na.rm=T) 
mean(btab$DurationMean[btab$What=="Dive" & btab$Zone == "Epipelagic"],na.rm=T)/60 # 11.45 min
sd(btab$DurationMean[btab$What=="Dive" & btab$Zone == "Epipelagic"],na.rm=T)/60
mean(btab$DepthMean[btab$What=="Dive" & btab$Zone == "Mesopelagic"],na.rm=T) # 441 m
sd(btab$DepthMean[btab$What=="Dive" & btab$Zone == "Mesopelagic"],na.rm=T)
mean(btab$DurationMean[btab$What=="Dive" & btab$Zone == "Mesopelagic"],na.rm=T)/60 # 24.01 min
sd(btab$DurationMean[btab$What=="Dive" & btab$Zone == "Mesopelagic"],na.rm=T)/60 
mean(btab$DepthMean[btab$What=="Dive" & btab$Zone == "Bathypelagic"],na.rm=T) # 1487 m
sd(btab$DepthMean[btab$What=="Dive" & btab$Zone == "Bathypelagic"],na.rm=T)
mean(btab$DurationMean[btab$What=="Dive" & btab$Zone == "Bathypelagic"],na.rm=T)/60 # 55.50 min
sd(btab$DurationMean[btab$What=="Dive" & btab$Zone == "Bathypelagic"],na.rm=T)/60 

### Mean and SD of max depth of benthic dives
mean(btab$DepthMean[btab$What=="Dive" & btab$ibentic]) # 1188 m
sd(btab$DepthMean[btab$What=="Dive" & btab$ibentic]) # 471 m 

### Number of truncated dives due to tag's limited depth range
ind = btab$What=="Dive" & btab$Ptt=="134668" & btab$Zone=="Bathypelagic" 
length(btab$DepthMean[ind])  # n=149 bathypelagic dives by mature male
table(btab$DepthMean[ind])   # n=21 of those to maximum depth of 1847.75 m

### Dive shapes
i2016 = btab$Datetime>"2016-01-01" # Indices for 2016 only

# Percentage across all dives
sum(btab$What=="Dive" & i2016 & btab$Shape == "V") / sum(btab$What=="Dive" & i2016) * 100
sum(btab$What=="Dive" & i2016 & btab$Shape == "U") / sum(btab$What=="Dive" & i2016) * 100
sum(btab$What=="Dive" & i2016 & btab$Shape == "Square") / sum(btab$What=="Dive" & i2016) * 100

# Percentages within dive categories
sum(btab$What=="Dive" & i2016 & btab$Zone=="Epipelagic" & btab$Shape == "Square") / sum(btab$What=="Dive" & i2016 & btab$Zone=="Epipelagic") * 100
sum(btab$What=="Dive" & i2016 & btab$Zone=="Mesopelagic" & btab$Shape == "Square") / sum(btab$What=="Dive" & i2016 & btab$Zone=="Mesopelagic") * 100
sum(btab$What=="Dive" & i2016 & btab$Zone=="Bathypelagic" & btab$Shape == "Square") / sum(btab$What=="Dive" & i2016 & btab$Zone=="Bathypelagic") * 100
sum(btab$What=="Dive" & i2016 & btab$Zone=="Epipelagic" & btab$Shape == "U") / sum(btab$What=="Dive" & i2016 & btab$Zone=="Epipelagic") * 100
sum(btab$What=="Dive" & i2016 & btab$Zone=="Mesopelagic" & btab$Shape == "U") / sum(btab$What=="Dive" & i2016 & btab$Zone=="Mesopelagic") * 100
sum(btab$What=="Dive" & i2016 & btab$Zone=="Bathypelagic" & btab$Shape == "U") / sum(btab$What=="Dive" & i2016 & btab$Zone=="Bathypelagic") * 100
sum(btab$What=="Dive" & i2016 & btab$Zone=="Epipelagic" & btab$Shape == "V") / sum(btab$What=="Dive" & i2016 & btab$Zone=="Epipelagic")  * 100
sum(btab$What=="Dive" & i2016 & btab$Zone=="Mesopelagic" & btab$Shape == "V") / sum(btab$What=="Dive" & i2016 & btab$Zone=="Mesopelagic") * 100
sum(btab$What=="Dive" & i2016 & btab$Zone=="Bathypelagic" & btab$Shape == "V") / sum(btab$What=="Dive" & i2016 & btab$Zone=="Bathypelagic") * 100

### Number of dives per category in 2014-2015 vs. 2016
n = matrix(NA, 3, 3)
n[1,1] = nrow(btab[!i2016 & btab$What=="Dive" & btab$Zone=="Epipelagic",]) # 937
n[1,2] = nrow(btab[!i2016 & btab$What=="Dive" & btab$Zone=="Mesopelagic",]) # 617
n[1,3] = nrow(btab[!i2016 & btab$What=="Dive" & btab$Zone=="Bathypelagic",]) # 209
n[2,1] = nrow(btab[i2016 & btab$What=="Dive" & btab$Zone=="Epipelagic",]) # 2692
n[2,2] = nrow(btab[i2016 & btab$What=="Dive" & btab$Zone=="Mesopelagic",]) # 2760
n[2,3] = nrow(btab[i2016 & btab$What=="Dive" & btab$Zone=="Bathypelagic",]) # 316
n[3,1] = sum(n[1:2,1]) # 3629
n[3,2] = sum(n[1:2,2]) # 3377
n[3,3] = sum(n[1:2,3]) # 525
n

#### Table 3: overall and individual dive rates, dive and surface durations, IDDIs 

### Dive rate in 2016: number dives / data hours
n[2,] / (sum(btab$DurationMean[i2016],na.rm=T)/3600) # 1.15, 1.18, 0.14

tab3 = matrix(NA, 6, 3)
for(i in 1:length(ptts2016)){
  n.edives = sum(btab$What=="Dive" & btab$Ptt==ptts2016[i] & btab$Zone=="Epipelagic")
  n.mdives = sum(btab$What=="Dive" & btab$Ptt==ptts2016[i] & btab$Zone=="Mesopelagic")
  n.bdives = sum(btab$What=="Dive" & btab$Ptt==ptts2016[i] & btab$Zone=="Bathypelagic")
  tab3[i,1] = n.edives / (sum(btab$DurationMean[btab$Ptt==ptts2016[i]],na.rm=T)/3600) 
  tab3[i,2] = n.mdives / (sum(btab$DurationMean[btab$Ptt==ptts2016[i]],na.rm=T)/3600) 
  tab3[i,3] = n.bdives / (sum(btab$DurationMean[btab$Ptt==ptts2016[i]],na.rm=T)/3600) 
}
round(tab3*100)/100

### Dive duration (mean,sd;min,max) in 2016
tab3 = matrix(NA, 2, 6)
i1 = btab$DurationMean[i2016 & btab$What=="Dive" & btab$Zone=="Epipelagic"]/60
i2 = btab$DurationMean[i2016 & btab$What=="Dive" & btab$Zone=="Mesopelagic"]/60
i3 = btab$DurationMean[i2016 & btab$What=="Dive" & btab$Zone=="Bathypelagic"]/60
tab3[1, ] = c(mean(i1),sd(i1),mean(i2),sd(i2), mean(i3),sd(i3))
tab3[2, ] = c(min(i1),max(i1),min(i2),max(i2),min(i3),max(i3))
round(tab3)

tab3 = matrix(NA, 12, 6)
for(i in 1:length(ptts2016)){
  i1 = btab$DurationMean[i2016 & btab$What=="Dive" & btab$Ptt==ptts2016[i] & btab$Zone=="Epipelagic"]/60
  i2 = btab$DurationMean[i2016 & btab$What=="Dive" & btab$Ptt==ptts2016[i] & btab$Zone=="Mesopelagic"]/60
  i3 = btab$DurationMean[i2016 & btab$What=="Dive" & btab$Ptt==ptts2016[i] & btab$Zone=="Bathypelagic"]/60
  tab3[(i*2)-1, ] = c(mean(i1),sd(i1),mean(i2),sd(i2), mean(i3),sd(i3))
  tab3[(i*2), ] = c(min(i1),max(i1),min(i2),max(i2),min(i3),max(i3))
}
round(tab3)

### Surface duration (mean,sd;min,max) in 2016
tab3 = matrix(NA, 2, 2)
i1 = btab$DurationMean[i2016 & btab$What=="Surface"]/60
tab3[1, ] = c(mean(i1,na.rm=T),sd(i1,na.rm=T))
tab3[2, ] = c(min(i1,na.rm=T),max(i1,na.rm=T))
round(tab3*10)/10
sum(i2016 & btab$What=="Surface" & !is.na(btab$DurationMean)) # n in 2016

tab3 = matrix(NA, 12, 2)
for(i in 1:length(ptts2016)){
  i1 = btab$DurationMean[i2016 & btab$What=="Surface" & btab$Ptt==ptts2016[i] & !is.na(btab$DurationMean)]/60
  tab3[(i*2)-1, ] = c(mean(i1),sd(i1))
  tab3[(i*2), ] = c(min(i1),max(i1))
}
round(tab3*10)/10

### IDDI in 2016 (mean,sd; min,max)
tab3 = matrix(NA, 2, 2)
i1 = btab$IDDI[i2016 & !is.na(btab$IDDI)]
tab3[1, ] = c(mean(i1,na.rm=T),sd(i1,na.rm=T))
tab3[2, ] = c(min(i1,na.rm=T),max(i1,na.rm=T))
round(tab3)

tab3 = matrix(NA, 12, 2)
for(i in 1:length(ptts2016)){
  i1 = btab$IDDI[i2016 & btab$Ptt==ptts2016[i] & !is.na(btab$IDDI)]
  tab3[(i*2)-1, ] = c(mean(i1),sd(i1))
  tab3[(i*2), ] = c(min(i1),max(i1))
}
round(tab3)

### % of time in each dive type (2016 only)
dur.tot = (sum(btab$DurationMean[i2016],na.rm=T)/3600) # 2339 hours (dives + surfs)
(sum(btab$DurationMean[i2016 & btab$What=="Surface" & !is.na(btab$DurationMean)])/3600) / dur.tot * 100 # 18%
(sum(btab$DurationMean[i2016 & btab$Zone=="Epipelagic" & !is.na(btab$DurationMean)])/3600) / dur.tot * 100 # 22%
(sum(btab$DurationMean[i2016 & btab$Zone=="Mesopelagic" & !is.na(btab$DurationMean)])/3600) / dur.tot * 100 # 47%
(sum(btab$DurationMean[i2016 & btab$Zone=="Bathypelagic" & !is.na(btab$DurationMean)])/3600) / dur.tot * 100 # 12%


############ Info on longest and deepest dive ############
# Longest dive
durmax = max(btab$DurationMean[btab$What=="Dive"], na.rm=T)/60 # 97.8333 min
btab$Start[which(btab$DurationMean==durmax*60)] # "2016-07-06 06:08:00 UTC"
btab$lat.pred[which(btab$DurationMean==durmax*60)] # 64.17344
btab$lon.pred[which(btab$DurationMean==durmax*60)] # -9.085513
btab$IDDI[which(btab$DurationMean==durmax*60)] # 88 min
btab$DepthMean[which(btab$DurationMean==durmax*60)] # 1232 m
btab$Bathy[which(btab$DurationMean==durmax*60)] # 1114 m bathymetric depth

# Deepest dives
depthmax = max(btab$DepthMean[btab$What=="Dive"], na.rm=T) # 2288 m
btab$Start[which(btab$DepthMean==depthmax)] # "2016-06-30 15:21:44 UTC" "2016-06-30 18:34:58 UTC"
btab$lat.pred[which(btab$DepthMean==depthmax)[1]] # 70.9433
btab$lon.pred[which(btab$DepthMean==depthmax)[1]] # -6.6149
btab$lat.pred[which(btab$DepthMean==depthmax)[2]] # 70.9776
btab$lon.pred[which(btab$DepthMean==depthmax)[2]] # -6.4692
btab$DurationMean[which(btab$DepthMean==depthmax)]/60 # 78 min, 73 min
difftime((as.POSIXct("2016-06-30 15:21:44 UTC") + (78*60)),"2016-06-30 18:34:58 UTC",units="mins") # dt = 115 mins


############ Correlation tests ############
btab.surf = btab[btab$What=="Surface",]
btab = btab[btab$What=="Dive",] # dives-only btab now!

# Max dive depth vs duration for 2016
cor.test(btab$DepthMean, btab$DurationMean, method = "spearman", exact=FALSE) # r=+0.86

# Max dive depth vs water depth of non-benthic deep dives # r=+0.43
cor.test(btab$DepthMean[btab$Zone=="Bathypelagic" & !btab$ibentic], btab$Bathy[btab$Zone=="Bathypelagic" & !btab$ibentic], method = "spearman", exact=FALSE) 

# Assign IDDI to dive at the end instead of the start
btab$IDDI2 <- NA
ideep = which(btab$Zone=="Mesopelagic" | btab$Zone=="Bathypelagic")
for(i in 1:nrow(btab)){
  if(!is.na(btab$IDDI[i])){
    btab$IDDI2[ideep[min(which(ideep>i))]] = btab$IDDI[i]
  }
}

# Compare IDDI with dive duration of dive before and after
btab0 = btab[!is.na(btab$DurationMean) & (btab$DurationMean>(20*60)),]

# Or only top 5% of deep dive durations following Quick et al. 2020
q95 = quantile(btab$DurationMean[(btab$Zone=="Mesopelagic" | btab$Zone=="Bathypelagic") & !is.na(btab$DurationMean)],0.95)
btab0 = btab[!is.na(btab$DurationMean) & (btab$DurationMean>=q95),]

# Dive duration vs IDDI (2016 only)
r12 = cor.test(btab0$DurationMean, btab0$IDDI, method = "spearman", exact=FALSE)
r13 = cor.test(btab0$DurationMean, btab0$IDDI2, method = "spearman", exact=FALSE)
r23 = cor.test(btab0$IDDI, btab0$IDDI2, method = "spearman", exact=FALSE)
diffcor::diffcor.dep(r12$estimate, r13$estimate, r23$estimate, sum(!is.na(btab0$IDDI) & !is.na(btab0$IDDI2)), digit = 4)

### Plots of Dive duration (top 5%) vs IDDI
a <- btab0[!is.na(btab$IDDI2),] %>% 
  ggplot(aes(IDDI2, DurationMean/60, groupColour = T, groupFill = T)) +
  geom_point(colour = "black", fill = "grey", size = 4, alpha = 0.1, shape = 21) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_x_continuous(trans='log10') +
  xlab("IDDI before dive (min)") +
  ylab("Dive duration (min)") +
  stat_smooth(method = "lm", formula = y ~ x, geom = "smooth", colour = "#111") 
b <- btab0[!is.na(btab$IDDI),] %>% 
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

# Percentages of consecutive epi- and mesopelagic dives
sum(diff(which(btab$Zone=="Epipelagic"))==1) / length(which(btab$Zone=="Epipelagic")) * 100
sum(diff(which(btab$Zone=="Mesopelagic"))==1) / length(which(btab$Zone=="Mesopelagic")) * 100

