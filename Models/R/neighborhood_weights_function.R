## Script by Andrew Crosby

library(sf)
library(spdep)
getwd()

covariates_path <- "C:/Users/ilebe/Documents/!Masters!/Analysis/0. Data/Processed/covariates"
obs_path <- "C:/Users/ilebe/Documents/!Masters!/Analysis/0. Data/Processed/All Processed Final Location"

pts <- st_read(file.path(covariates_path,"coordinatesonly.shp")) # An sf object of all the point locations that lines up with the rest of the data (e.g., the detections and the envirionmental variables)

site_det <- read.csv(file.path(obs_path,"OCCU.csv"))# A site x spp detection matrix 
dim(site_det)
head(site_det)

# Identify neighbors by Euclidean distance of 1000m
nb<-dnearneigh(pts, d1=0, d2=1000) # creates a list, where each list item is a site, and the values refer to the row numbers of the ponits within 1000m 

# Create the 1000-meter neighborhood matrix
nb_wb<-nb2mat(nb, style="B", zero.policy=TRUE)

  # Create the 1000-meter neighborhood weights matrix
n <- 6
R <- dim(nb_wb)[1]

num_neigh<-rep(NA, R)
for(i in 1:R){
  num_neigh[i]<-sum(nb_wb[i,])
}
 length(which(num_neigh==0))
 
 sum_neighbors<-sum(num_neigh)
 
 print(names(site_det)[-1])
 
site_spec_weights<-matrix(0, nrow = R, ncol =  6)
# dimnames(site_det)[2][[1]][-1]
head(site_spec_weights)
colnames(site_spec_weights) <- dimnames(site_det)[2][[1]][-1]
head(site_spec_weights)

for(j in 1:R){
  for(i in 1:n){
    if(num_neigh[j]==0){
      site_spec_weights[j, i]<-0
    }else{
      d<-0
      e<-num_neigh[j]
      for(k in 1:R){
        if(nb_wb[j, k]==0){
          temp<-0
        }else{
          temp<-site_det[k, i+ 1]
        }
        d<-d+temp
      }
      site_spec_weights[j, i] <- d/e
    }
  }
}
 
weights_1000 <- site_spec_weights
rownames(weights_1000) <- pts$location
head(weights_1000)

# Convert the matrix to a data frame
weights_1000_df <- as.data.frame(weights_1000)

# Add a 'location' column at the beginning
weights_1000_df <- data.frame(location = rownames(weights_1000_df), weights_1000_df)

write.csv(weights_1000, file.path(covariates_path, "weights_1000.csv"), row.names=TRUE)



