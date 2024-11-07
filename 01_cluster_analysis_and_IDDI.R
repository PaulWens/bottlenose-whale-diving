# Script 001: Cluster analysis and IDDI
# - K-means cluster analysis on max depth and duration from 2016
# - Imputation clustering for 2014-2015 duration data
# - IDDI calculation

############ User defined ############
ptts <-c(134663,134664,134666,134667,134668,134669,134670,161587,161588,161590,161591,161592,161593)
sstart2015 <- "2015-06-29 02:48:00" # Start of sonar exposure experiment 2015
sstart2016 <- "2016-06-18 12:16:00" # Start of sonar exposure experiment 2016

############ Load R packages ############
library(dplyr)
library(ClustImpute)
library(factoextra)
library(ggplot2)

############ Load data tables ############
# "btab" has the relevant fields in the original Behavior Log, vertical velocity, 
# data gaps, latitude and longitude, bathymetric depth, oceanographic zone.
# "stab" includes the Timeseries depth data with similar additional variables.
load("02_data_tables.Rdata")

# Exclude sonar periods, tagging periods, very shallow dives, and surface periods
btab = btab[!btab$isonar & !btab$itagging & !btab$ishallow & btab$What == "Dive", ]

# Set to NA the dive durations from 2014-2015 that could be wrong
ind = btab$Datetime < "2016-01-01 00:00:00"
btab$DurationMean[ind] = NA

############ K-means cluster analysis on 2016 data ############
# Cluster variables: log(max depth) and sqrt(duration)
btab.clus <- btab[,c("DepthMean","DurationMean")]
btab.clus$DepthMean = log(btab.clus$DepthMean)
btab.clus$DurationMean = sqrt(btab.clus$DurationMean)

# Standardize these variables 
mu = rep(NA,length(btab.clus))
sds = rep(NA,length(btab.clus))
for(i in 1:length(btab.clus)){
  mu[i] = mean(btab.clus[,i],na.rm=T)
  sds[i] = sd(btab.clus[,i],na.rm=T)
  btab.clus[,i] = (btab.clus[,i] - mu[i]) / sds[i]
}

# Remove observations from 2014-2015 for which Duration can be wrong 
ind = !is.na(btab.clus$DurationMean)
clus = btab.clus[ind,] # transformed data for clustering 
clus.org = btab[ind,]  # non-transformed data for plotting

# Determine optimal K: Elbow method
fviz_nbclust(clus, kmeans, method = "wss", k.max = 8) 

# Determine optimal K: Silhouette method
a <- fviz_nbclust(clus, kmeans, method = "silhouette") +
  theme_minimal() +
  ggtitle("")

# K-means clustering with K=2 clusters
km <- kmeans(clus, centers = 2, nstart = 25)

# Check silhouette widths to validate clusters
# -All clusters should have silhouette score greater than the average
# -Cluster size shouldn't fluctuate widely
sil <- cluster::silhouette(km$cluster, dist(clus)) 
b <- fviz_silhouette(sil, palette = c("#56B4E9", "darkred")) + 
  theme_minimal() +
  coord_flip() +  
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_x_discrete(breaks = "") +
  ggtitle("") +
  labs(y = "Sihouette width", x = "")

# Plot Silhouette results
g <- ggpubr::ggarrange(a, b, widths = c(0.7,1), ncol = 2)
print(g)

# Max depth vs duration - original, non-transformed scale
margplot <- clus.org %>% 
  as_tibble() %>%
  mutate(cluster = km$cluster) %>%
  ggplot(aes(DurationMean/60, DepthMean, color = factor(cluster))) +
  geom_point(size = 2, alpha = 0.2) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) + 
  scale_y_reverse(breaks = seq(0, 2250, by = 250)) + # reverse y-axis
  labs(x = "Dive duration (min)", y = ("Maximum dive depth (m)"), color = "Cluster") +
  scale_color_manual(values = c("#56B4E9", "darkred","#999999"))
ggExtra::ggMarginal(margplot + theme(legend.position = "none"),
                    groupColour = T, groupFill = T)


############ K-means cluster analysis on 2014-2016 data ############
## Imputation clustering to impute missing dive durations for 2014-2015 

# Determine optimal K: Elbow method
set.seed(321)
ClustImpute2 <- function(dataFrame, nr_cluster, nr_iter=10, c_steps=1, wf=default_wf, n_end=10, seed_nr=150519) {
  return(ClustImpute(dataFrame,nr_cluster, nr_iter, c_steps, wf, n_end, seed_nr))
}
res_list <- lapply(X=1:8, FUN=ClustImpute2, dataFrame=btab.clus, c_steps=50)
tmp <- var_reduction(res_list[[1]])
var_red_total <- tmp$Variance_reduction 
var_by_clust <- tmp$Variance_by_cluster
for (k in 2:8) {
  tmp <- var_reduction(res_list[[k]])
  var_by_clust <- rbind(var_by_clust,tmp$Variance_by_cluster)
  var_red_total <- rbind(var_red_total,tmp$Variance_reduction)
}
var_by_clust$nr_clusters <- 1:8

# Plot ratio of the sum of all within cluster variances by the overall variance
plot(1:8, 1-var_red_total, type="b", ylab = "Variance ratio", xlab="No. of clusters")

# Keep results for K=2
res <- res_list[[2]]

# Create new data frames with imputed durations
btab.clus2 = res$complete_data # same as btab.clus but with imputed values
btab.clus3 = btab.clus2 # same as btab.clus2 but back-transformed:
btab.clus3$DepthMean = exp((btab.clus2$DepthMean * sds[1]) + mu[1])
btab.clus3$DurationMean = ((btab.clus2$DurationMean * sds[2]) + mu[2])^2

### Save cluster assignment (dive type) and imputed durations
# If needed, switch cluster ids so that 1=short/shallow and 2=long/deep (CHECK AFTER RERUNNING)
load("02_data_tables.Rdata")
btab$DurationMean2 = NA
btab$cluster = NA
km <- kmeans(btab.clus2, centers = 2, nstart = 25)
btab$DurationMean2[!btab$isonar & !btab$itagging & !btab$ishallow & btab$What == "Dive"] = btab.clus3$DurationMean
btab$cluster[!btab$isonar & !btab$itagging & !btab$ishallow & btab$What == "Dive"] = km$cluster
btab$cluster = btab$cluster - 1 # switch cluster id (2 -> 1)
btab$cluster[btab$cluster == 0] = 2 # switch cluster id (1 -> 2)
#save(btab, stab, file="03_data_tables.Rdata")

### Calculate mean MaxDepth and Duration per dive type in Table 2
btab0 = btab[!btab$isonar & !btab$ishallow & !btab$itagging, ]
mean(btab0$DepthMean[btab0$cluster==1],na.rm=T)
mean(btab0$DepthMean[btab0$cluster==2],na.rm=T)
mean(btab$DurationMean[btab$cluster==1],na.rm=T)/60
mean(btab$DurationMean[btab$cluster==2],na.rm=T)/60
sum(!is.na(btab$DurationMean[btab$cluster==1]))
sum(!is.na(btab$DurationMean[btab$cluster==2]))


############ Inter deep dive interval ############
# IDDI = Time elapsed between the end of a long dive and start of the next long dive
btab$IDDI<-NA
ind = which(btab$cluster == 2 & !is.na(btab$cluster)) # indices of deep dives
iddi = -(as.numeric(btab$End[ind[1:(length(ind)-1)]] - btab$Start[ind[2:(length(ind))]]))/60
btab$IDDI[btab$cluster == 2 & !is.na(btab$cluster)] = c(iddi,NA)

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
#save(btab, stab, file="03_data_tables.Rdata")
