library(tidybayes)
library(tidyverse)
library(ncdf4)

idata_path <-"C:/Users/ilebe/Documents/!Masters!/Analysis/AbundanceConditionalOccupancy/Model/Inference Data"

nc_data <- nc_open(file.path(idata_path, "ZIP_linear_occupancy.nc"))

class(nc_data)

M_pred <- ncvar_get(nc_data, "M_pred")

data_frame <- as.data.frame(matrix(ncvar_get(nc_data, varid="variable_name"), ncol=variable_dimension))

library(ncdf4)

# Read in the posterior_predictive data
post_predict_nc <- nc_open(file.path(idata_path,"ZIP_linear_occupancy_post_predict.nc"))
M_pred <- ncvar_get(post_predict_nc, "M_pred")

# Read in the posterior data
posterior_nc <- nc_open(file.path(idata_path, "ZIP_linear_occupancy_posterior.nc"))
alpha_theta <- ncvar_get(posterior_nc, "alpha_theta")
beta_longitude <- ncvar_get(posterior_nc, "beta_longitude")

longitudes <- seq(from = -119.82829, to = -110.67583, length.out = length(M_pred))
theta <- 1 / (1 + exp(- (mean(alpha_theta) + mean(beta_longitude) * longitudes)))

# Assuming you've already loaded the necessary data from the NetCDF files...
library(tidybayes)

theta_samples <- posterior_nc$theta_rep
post_predict_nc
theta_samples
prediction_data <- gather_rvars(theta_samples) 
