############    Overview & Install packages     ######
#' Goal: Create Figure illustrating biome-specific C accumulation curves of ecosystem regeneration
#' with and without adjustments for the net CH4-N2O effect (in CO2 equivalent). 
#' 
#' Analysis steps:
#'   a) Obtain parameters for biome-specific C accumulation curves from Nathaniel Robinson et al. (in prep)
#'   b) Plot biome-specific C accumulation curves 
#'   c) Compute the integrals of the net CH4-N2O effect in CO2e for each biome (from t=0 to 100 years). 
#'                 Starting units: [Mg CO2e ha-1 year-1 ]
#'                 Final units: [Mg CO2e ha-1 ]
#'  d) Final plot: add the net CH4-N2O effect (in CO2-equivalent) to the biome-specific C accumulation curves (CO2-only).  
#'  The net CH4-N2O effect will be a constant number since the Age coefficient was not significant in model results. 
#'  There will be two curves per biome, the original CO2-only curve and the added net CO2-CH4-N2O effect . 
#'  
#'  
#' Robinson et al (in prep) Modified Chapman-Richards AGC model: 
#' Above ground C = a * pow(1-(b*exp(-k*t)), 1/(1-m))
#' m=2/3
#' a, b, k computed with spatially-explicit RF models by Robinson et al (in prep)

#'  
#'  Net CH4-N2O effect based on meta-analysis results from Cooley et al (in prep) 
#'  
# Biome	Net CH4-N2O effect (Mg C ha-1 year-1)	+/- Standard Error

library(raster)
library(data.table)
library(ggplot2); library(gridExtra)
library(dplyr)
library(tidyverse)
library(gridExtra)
library(reshape2) 
library(RColorBrewer)
library(FNN)
library(sjPlot)

dir<-"/Users/savannah/Documents/Ch4/Data/Robinson_et_al_AGC_model_100yrs/"

############    Functions for the scientific model   ##########
# Function for the scientific model based on Robinson et al (in prep). 

#' Model 4
yhat_AGB_m4 = function(A, age){
  #mu= b*age/ (d+age)
  return(mu)
}



# Other model options that I explored
#' Model 1
yhat_AGB_m1 = function(c, z, age){
  mu= exp(c)*age^z
} 

#' Model 2
yhat_AGB_m2 = function(b, d, age){
  mu= b*age/ (d+age)
  return(mu)
}

# Model 3: Negative exponential
yhat_AGB_m3 = function(k, phi, age){
  mu= k*(1-exp(-phi*age))
  return(mu)
}

############    Mock up Fig: Michaelis-Mentin AGB model with CH4-N2O results      ######

#' Model 2: Michaelis-Mentin AGB model 
yhat_AGB_m2 = function(b, d, age){
  mu= b*age/ (d+age)
  return(mu)
}

# Simulate data for the independent variable
x= runif(n=10000, min = 1, max = 100)# c(1:35) n=nrow(dt_all)
Biome <- c("Temperate conifer forest","Temperate broadleaf forest",
           "Subtropical/tropical savanna","Subtropical/tropical forest",
           "Boreal forest", "GLOBAL mean") 

# AGbiomass/2 = AG Carbon. Aboveground Carbon *44/12 = Aboveground CO2 mass
b_params <- c(85, 90, 95, 180,  50,  83)*44/12 # b_params are ~half of ch2 Amazon b param bc AGB is ~50% Carbon. To go from C mass to CO2 mass, multiply by 44/12
d_params <- c(7.2, 7.3, 8.5, 8.4, 7.0, 6.2)# d_params+4.5 
net_CH4_N2O<- c(1.07, 0.06, -0.01, 1.64, 2.77, 1.106) # Back-transformed B0 intercepts in [Mg CO2e ha-1 yr-1] from Cooley et al (in prep) - Biome	Net CH4-N2O effect in [Mg C ha-1 year-1]
# vJan4: net_CH4_N2O<- c(1.081, -0.069, 0.23, 1.313, 4.572/2, 1.4254) # Back-transformed B0 intercepts in [Mg CO2e ha-1 yr-1] from Cooley et al (in prep) - Biome	Net CH4-N2O effect in [Mg C ha-1 year-1]
age_coeff<- c(0.03, 0, 0.023, 0, 0, 0.009) # Back-transformed in [Mg CO2e ha-1 yr-1 yr-1 * Age (Year)]. Only consider significant coefficients for age 
# vJan4: age_coeff<- c(0.026,0.028,0.022,0.021,0.023, 0.024) # Include all coefficients for age, regardless of statistical significance
ag_CH4_N2O<- c(1.12, 2.43, 1.45, 3.73, 0.409, 1.8278) # mean agriculture CH4-N2O effect in [Mg CO2e ha-1 yr-1]

# Decided not to do this: Adjust b and d to make sense w starting y at 0 [will wait ]
# Decided not to do this for now: add error estimates to this figure

# Integrate x=0 to 100: net_CH4_N2O = B0 + B1*x. So net_CH4_N2O = B0*x + (B1*x^2)/2
# Do I want to weight the global mean by the area each biome occupies on earth?

par(mfrow = c(3, 2))

# Plot the scientific Model2 predictions across different biomes  
for(i in 1:length(Biome)){ 
  y_orig = yhat_AGB_m2(b=b_params[i], d=d_params[i], age = x)  
  y_CH4_N2O_int = y_orig - (age_coeff[i]*x^2)/2 - net_CH4_N2O[i]*x 
  y_CH4_N2O_ag = -ag_CH4_N2O[i]*x
  y0 = 0*x
  
  plot(x,y_orig, col="forestgreen", cex=0.6, pch=15, cex.axis=1.8, cex.lab=0.1, 
       xlab= NULL, #xlab= "Forest Age (Years)", 
       ylab= NULL, #ylab= expression(atop("Radiative cooling effect", "(CO"[2]*"e" ~"ha"^-1*")")), #"Radiative effect (Mg CO2e ha-1)", 
       main=Biome[i], cex.main=2, 
       ylim = c(-ag_CH4_N2O[i]*100, max(y_orig)+10) ) #  main="Scientific model predictions",
  
  points(x,y_CH4_N2O_int, col="green2", xlab= NULL,ylab= NULL,cex.axis=1.8, cex.lab=0.1, cex=0.5, pch=1) # pch = 0:25)
  points(x,y_CH4_N2O_ag, col="purple",  xlab= NULL,ylab= NULL,cex.axis=1.8, cex.lab=0.1, cex=0.5, pch=1)  
  points(x,y0, col="grey",  xlab= NULL,ylab= NULL,cex.axis=1.8, cex.lab=0.1, cex=0.05, pch=1)  
  }



############    Assemble data: Robinson et al AGC model with CH4-N2O results      ######

#' Aggregate AGC mean (y), SD, min, max by biome type in new df

# Import global biome data set (from 2017 resolve ecoregions) as data table https://resourcewatch.org/data/explore/bio042-Ecoregion-by-Biome
r<-raster("/Users/savannah/Documents/Ch3/Data/Biomes_ecoregions/Biomes_11km_res.tif")
plot(r, col=bpy.colors(14)) #crs(r)<-"epsg:3857"; crs(r)<-"epsg:4326"; plot(r)

# reclassify the values into biome groups 
# c(x1,x2,x3) all values > x1 and <= x2 become x3, etc.
m <- c(4, 5, 1001, # 'Temperate conifer forest'
       3, 4, 1002, # 'Temperate broadleaf forest'
       11, 12, 1003, 6, 7, 1003, # 'Subtropical/tropical savanna'
       0, 3, 1004,  # 'Subtropical/tropical forest'
       5, 6, 1005, # 'Boreal forest'
       7, 11, NA, # Other biomes to NA
       12, 14, NA) # Other biomes to NA

rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc <- reclassify(r, rclmat); plot(rc)

# Read in the raster files as a stack
ls<-list.files(path = dir, pattern = "*.tif", full.names=T)
r_stack <- stack(ls)

# Get the distinct categories from the 'rc' raster
rc <-resample(rc, r_stack$AGC.1)
categories <- unique(rc[])

# Initialize a list to store the results
mean_results <- list()
sd_results <- list()

# Loop through each raster layer in the stack
for (i in 1:nlayers(r_stack)) {
  print(i)
  # Extract the current raster layer
  r_layer <- r_stack[[i]]
  
  # Perform zonal statistics
  mean_zonal_stats <- zonal(r_layer, rc, fun = "mean", na.rm = TRUE)
  sd_zonal_stats <- zonal(r_layer, rc, fun = "sd", na.rm = TRUE)
  
  # Store the results in the list
  mean_results[[i]] <- mean_zonal_stats[, 2, drop = FALSE]  # Extract only the mean values
  sd_results[[i]] <- sd_zonal_stats[, 2, drop = FALSE]  # Extract only the mean values
  
}

# Convert the list of results to a data frame
mean_df <- do.call(cbind, mean_results); 
mean_df <- data.frame( cbind(mean_df[1,1:100],mean_df[2,1:100],
                             mean_df[3,1:100],mean_df[4,1:100],
                             mean_df[5,1:100]) ) #mean_df <- as.data.frame(mean_df)
# Set the column names
colnames(mean_df) <- c("Temperate conifer forest","Temperate broadleaf forest",
                       "Subtropical/tropical savanna","Subtropical/tropical forest",
                       "Boreal forest") 
mean_df$Global.mean <- rowMeans(mean_df)

sd_df <- do.call(cbind, sd_results)
sd_df <- data.frame( cbind(sd_df[1,1:100],sd_df[2,1:100],
                           sd_df[3,1:100],sd_df[4,1:100],
                           sd_df[5,1:100]) ) #sd_df <- as.data.frame(sd_df)
# Set the column names
colnames(sd_df) <- c("Temperate conifer forest","Temperate broadleaf forest",
                     "Subtropical/tropical savanna","Subtropical/tropical forest",
                     "Boreal forest") 
sd_df$Global.mean <- rowMeans(sd_df)

write.csv(mean_df, paste(dir, "mean_df_AGC.yrs1_100.csv"))
write.csv(sd_df, paste(dir, "sd_df_AGC.yrs1_100.csv"))


############    Plot Fig4: AGC accumulation & net CH4-N2O effect including uncertainties      ######

mean_df<- read.csv(paste(dir, "mean_df_AGC.yrs1_100.csv"))
sd_df<- read.csv(paste(dir, "sd_df_AGC.yrs1_100.csv"))

# The atomic weight of a carbon atom is 12 and the atomic weight of oxygen is 16, 
# so the total atomic weight of CO2 is 44 (12 + (16 * 2) = 44). This means that a quantity of CO2 can be expressed in terms of the amount of carbon it contains 
# by multiplying the amount of CO2 by 0.27 (12/44).
mean_df <- mean_df*44/12; mean_df$Age<-seq(1,100,1); mean_df$X<-NULL
sd_df <- sd_df*44/12; sd_df$Age<-seq(1,100,1); sd_df$X<-NULL

# Remove outlier obs (artifact of Robinson et al code? seems to happen every 10 yrs or so)
par(mfrow = c(1, 1)); plot(mean_df$Age,mean_df$Temperate.broadleaf.forest )
plot(mean_df$Age,mean_df$Temperate.conifer.forest )
plot(sd_df$Age,sd_df$Subtropical.tropical.forest )

# Upon taking a closer look, I see the outliers happen every 11 years starting at 13 years: 
# 13, 24, 35, 46, ... So I will simply remove these
outliers<-seq(13, 100, 11)
# Option1: Remove rows based on index ordering
mean_df <- mean_df[-outliers, ]
sd_df <- sd_df[-outliers, ]
# Option2: Change the ages to be between 2-3 years since the og values are between those of 2-3 years
# new_Ages<-seq(2.1, 2.8, 0.1)
# mean_df[outliers, ]$Age <- new_Ages
# sd_df[outliers, ]$Age <- new_Ages


# Simulate data for the independent variable
x= runif(n=1000, min = 1, max = 100)# c(1:35) n=nrow(dt_all)

# Create a new data frame with the new Age values
new_df <- data.frame(Age = x)
sd_new_df <- data.frame(Age = x)

# Interpolate values for each column based on the new Age values
for (col in names(mean_df)[-7]) {
  
  # Use the approx() function to interpolate the y values for each biome (column) 
  new_df[[col]] <- approx(mean_df$Age, mean_df[[col]], xout = x, method = "linear")$y
  sd_new_df[[col]] <- approx(sd_df$Age, sd_df[[col]], xout = x, method = "linear")$y
}

head(new_df)

######

# I followed this specific biome ordering for all variables below
Biome <- c("Temperate conifer forest","Temperate broadleaf forest",
           "Subtropical/tropical savanna","Subtropical/tropical forest",
           "Boreal forest", "GLOBAL mean") 

# Results from Table 1
# vJan4: net_CH4_N2O<- c(1.081, -0.069, 0.23, 1.313, 4.572/2, 1.4254) # Back-transformed B0 intercepts in [Mg CO2e ha-1 yr-1] from Cooley et al (in prep) - Biome	Net CH4-N2O effect in [Mg C ha-1 year-1]
# vMar2024: net_CH4_N2O<- c(0.65, -0.79, 0.24, 1.64, 2.44, 0.83) # Back-transformed B0 intercepts in [Mg CO2e ha-1 yr-1]  
# vNov2024: net_CH4_N2O<- c(0.633104045, 1.931022676, 0.230676892, 0.984369964, 1.446009465,1.045037) # Back-transformed B0 intercepts in [Mg CO2e ha-1 yr-1]  
net_CH4_N2O<- c(0.913609648, 1.655289318, 0.207768513, 0.657922094, 1.849368525,1.056792) # Back-transformed B0 intercepts in [Mg CO2e ha-1 yr-1]  

# Age B1 coefficients added together Back-transformed CH4.B1+N2O.B1 in [Mg CO2e ha-1 yr-1 yr-1 * Age (Year)] ie. CH4model_B1 + N2Omodel_B1 for each biome 

# vJan4: age_coeff<- c(0.026,0.028,0.022,0.021,0.023, 0.024) # Include all coefficients for age, regardless of statistical significance. # Considered: ?Only consider significant coefficients for age - Back-transformed B1 in [Mg CO2e ha-1 yr-1 yr-1 * Age (Year)].
# vMar2024: age_coeff<- c(0.031,0.030,0.022,0.02122,0.02288, 0.02542 )  
# vNov2024: age_coeff<- c(0.03085,0.03217,0.02213,0.02151,0.02527,0.026386) 
age_coeff<- c(0.030948332,0.032474767,0.022125239,0.021309671,0.025149067, 0.02640142) # Current approach: age_coeff_net_CH4_N2O = CH4model_B1_backtransformed + N2Omodel_B1_backtransformed 
# Equivalent to: age_coeff<- c(sum(0.012770294,0.018178038),etc, etc, etc, sum(0.012693903, 0.012455165), mean(age_coeff[1:5]) )

# v1 2023: ag_CH4_N2O<- c(1.12, 2.43, 1.45, 3.73, 0.409, 1.8278) # v1 mean agriculture CH4-N2O effect in [Mg CO2e ha-1 yr-1]
# vMar2024: ag_CH4_N2O<- c(3.48, 3.53, 3.01, 6.72, 1.87,3.722) # v2 w pasture adjusted mean agriculture CH4-N2O effect in [Mg CO2e ha-1 yr-1]
ag_CH4_N2O<- c(3.05, 2.53, 2.06, 1.55, 1.87,2.212) # Nov2024 v3: no rice CH4 & pasture adjusted mean agriculture CH4-N2O effect in [Mg CO2e ha-1 yr-1]
# vMar2024: ag_CH4_N2O_SD<- c(23.4, 21.8, 12.1, 224, 2.12 ) # Standard deviation of agriculture CH4-N2O effect in [Mg CO2e ha-1 yr-1]
ag_CH4_N2O_SD<- c(13.2, 23.1, 8.76, 14.5, 1.38,12.188) # Nov2024v3: Standard deviation of agriculture CH4-N2O effect in [Mg CO2e ha-1 yr-1]

# I quantified the propagation of uncertainty for the net CH4-N2O effect of 
# regeneration through time using an analytical approach: 

#' σ_f= (〖(t σ_B0_br) ) )〗^2 +〖(t^2/2)(σ_B1_br)〗^2 + t^3/2 Cov(B0_br,B1_br) ) ^ 0.5
  
#' σ_B0_br = (CH4_B0_SE^2 + N2O_B0_SE^2)^0.5
# vNov2024: net_B0_SE <- c(0.38, 1.59, 0.07, 0, 0, 0.4084837) # Back-transformed SE of B0 intercepts in [Mg CO2e ha-1 yr-1]   
net_B0_SE <- c(0.249477092, 0.691284465, 0.063308363, 0.317685313, 1.073155474, 0.4789821) # Back-transformed SE of B0 intercepts in [Mg CO2e ha-1 yr-1]   

#' σ_B1_br = net_B1_SE <- (CH4_B1_SE^2 + N2O_B1_SE^2)^0.5
net_B1_SE<- c(0.001138995, 0.000508702, 4.87086E-05, 0.000785088, 0.000310482, 0.0005583951) # 

#' Cov(B0_br,B1_br) = (Cov_CH4^2 + Cov_N2O^2)^0.5 
# I computed vcov(m3) for each biome-specific model and got the following covariance numbers for B0 and B1 matrix entries:
Cov_CH4 = c(2.955053e-07, 7.089406e-06, 2.021907e-08, 2.069347e-08, 5.848122e-06, 2.654789e-06)
Cov_N2O = c(-0.0003960923, 0.0004004415, 2.061554e-05, 0.002300802, 0.001104355, 0.0006860243)
net_Cov_B0_B1 = (Cov_CH4^2 + Cov_N2O^2)^0.5 

######
# Integrate x=0 to 100: net_CH4_N2O = B0 + B1*x. So integrated net_CH4_N2O = B0*x + (B1*x^2)/2

# Initialize an empty list to store the plots
plot_list <- list() # par(mfrow = c(3, 2))

# Plot the model predictions across biomes
for(i in 1:length(Biome)){
  y_orig = new_df[,i+1]
  sd_y_orig = sd_new_df[,i+1]
  y_CH4_N2O_int = y_orig - ((age_coeff[i]*x^2)/2 + net_CH4_N2O[i]*x )
  SE_f = ( (x*net_B0_SE[i])^2 + ((x^2/2)*net_B1_SE[i])^2 + (x^3/2)*net_Cov_B0_B1[i] ) ^ 0.5
  y_CH4_N2O_int_SEplus = y_CH4_N2O_int + (sd_y_orig^2 + SE_f^2)^0.5
  y_CH4_N2O_int_SEminus = y_CH4_N2O_int - (sd_y_orig^2 + SE_f^2)^0.5
  
  y_CH4_N2O_ag = -ag_CH4_N2O[i]*x
  y_CH4_N2O_ag_SDplus = y_CH4_N2O_ag + ag_CH4_N2O_SD[i]
  y_CH4_N2O_ag_SDminus = y_CH4_N2O_ag - ag_CH4_N2O_SD[i]
  
  y0 = 0*x
  
  # Create a data frame for plotting
  df <- data.frame(x = x,
                   y_orig = y_orig,
                   y_orig_upper = y_orig + sd_y_orig,
                   y_orig_lower = y_orig - sd_y_orig,
                   y_CH4_N2O_int = y_CH4_N2O_int,
                   y_CH4_N2O_int_SEplus = y_CH4_N2O_int_SEplus,
                   y_CH4_N2O_int_SEminus = y_CH4_N2O_int_SEminus,
                   y_CH4_N2O_ag = y_CH4_N2O_ag,
                   y_CH4_N2O_ag_SDplus = y_CH4_N2O_ag_SDplus,
                   y_CH4_N2O_ag_SDminus = y_CH4_N2O_ag_SDminus,
                   y0 = y0)
  
  p <- ggplot(df, aes(x = x)) +
    geom_ribbon(aes(ymin = y_orig_lower, ymax = y_orig_upper), fill = "forestgreen", alpha = 0.3) +
    geom_point(aes(y = y_orig), col = "forestgreen", size = 1.5, shape = 15) +
    geom_ribbon(aes(ymin = y_CH4_N2O_ag_SDminus, ymax = y_CH4_N2O_ag_SDplus), fill = "orange", alpha = 0.3) +
    geom_point(aes(y = y_CH4_N2O_ag), col = "orange", size = 1.5, shape = 1) +
    geom_point(aes(y = y0), col = "grey", size = 0.8, shape = 1) +
    geom_ribbon(aes(ymin = y_CH4_N2O_int_SEminus, ymax = y_CH4_N2O_int_SEplus), fill = "green2", alpha = 0.3) +
    geom_point(aes(y = y_CH4_N2O_int), col = "green2", size = 1.5, shape = 15) +
    
    labs(title = Biome[i], x = NULL, y = NULL) + # x= "Ecosystem Age"
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"), 
          axis.text = element_text(size = 12) ) + #,
          #axis.ticks = element_blank(),
          #panel.grid = element_blank()) + 
          coord_cartesian(ylim = c(-ag_CH4_N2O[i]*100, max(y_orig+sd_y_orig)+10))
  
  # Store the plot in the list
  plot_list[[i]] <- p
}

# Arrange the plots in a grid
grid.arrange(grobs = plot_list, nrow = 3, ncol = 2)

# 70/236 = 0.2/x --> x=0.2*236/70 = 0.6743


