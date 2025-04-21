#####################################################
# Notes  
#####################################################
'''
Goal: Evaluate statistical models predicting CH4 and N2O fluxes globally.
Steps:
a) Identify best-performing biome-specific glmers based on RF variable importance output
b) Create figures, including:

The code is organized as follows:

1.	Helper Functions:
o	extract_model_diagnostics(): Computes model diagnostics (normality, heteroscedasticity)
o	get_covariate_names(): Extracts the names of covariates from model formula
o	back_transform_ch4(), back_transform_n2o(): Transform coefficients back to original scale
o	create_model_summary_table(): Main function to process models by gas type

2.	Model Processing:
o	Refit all models (m0-m3) for each biome-gas combination
o	Extract all relevant statistics and diagnostics
o	Format results for presentation

3.	Output Generation:
o	Complete, detailed CSV files
o	Figures 2 and 3
o	README documentation
'''


# Determining the best model by AIC hardcoded as follows:
# 
# For CH4:
# - Boreal forest: m0
# - Subtropical/tropical forest: m2
# - Subtropical/tropical savanna: m3
# - Temperate broadleaf forest: m0
# - Temperate conifer forest: m1
#
# For N2O:
# - Boreal forest: m0
# - Subtropical/tropical forest: m0
# - Subtropical/tropical savanna: m3
# - Temperate broadleaf forest: m0
# - Temperate conifer forest: m3
#
# These models were selected based on having the lowest AIC score among models m0-m3
# for each GHG-biome combination.
#####################################################
# Prep: install packages, define root dir 
#####################################################
library(data.table)
library(dplyr)
library(ggplot2)
library(terra)
library(raster)
library(lme4)
library(sjPlot)
set_theme(base = theme_classic(), #axis.tickslen = 0, # hides tick marks
          axis.title.size = 1.8, axis.textsize = 1.5, #legend.size = .7,  legend.title.size = .8,
          geom.label.size = 4.5, title.size = 2.4)

root="/Users/savannah/Documents/Ch4/Data/Data_extraction/"
figsDir="/Users/savannah/Documents/Ch4/Figs_and_tables/"


#####     v3 post-GCB: Ag_net_CH4_N2O exclude rice cultivation areas for CH4     ####   

# Read rice cultivation data
rice<- raster("/Users/savannah/Documents/Ch4/Data/Ag_opportunity_cost/CROPGRIDSv1.08_rice.nc",
              varname = "croparea") # varname can set to one of: harvarea, croparea, qual, set

dt_v2<- fread("/Users/savannah/Documents/Ch4/Data/Ag_opportunity_cost/__dt_Ag_CH4_N2O_xy_pastureV2.csv")

# Read rice cultivation data
rice <- raster("/Users/savannah/Documents/Ch4/Data/Ag_opportunity_cost/CROPGRIDSv1.08_rice.nc",
               varname = "croparea")

# Extract rice area values for each point using raster::extract
dt_v2[!is.na(Lon) & !is.na(Lat), rice_area := raster::extract(rice, cbind(Lon, Lat))]

# Create new columns
dt_v2[, `:=`(
  # Set CH4 to 0 where rice is present (rice_area > 0)
  Ag_CH4_CO2e_noRice = fifelse(!is.na(rice_area) & rice_area > 0, 0, Ag_CH4_CO2e),
  
  # Calculate new total (CH4 excluding rice areas + original N2O)
  Ag_noRice_CH4_N2O = NA_real_
)]

# Fill in total where we have both components
dt_v2[!is.na(Ag_CH4_CO2e_noRice) & !is.na(Ag_N2O_CO2e), 
      Ag_noRice_CH4_N2O := Ag_CH4_CO2e_noRice + Ag_N2O_CO2e]

# Print summary of changes
print("Summary of changes:")
print(dt_v2[, .(
  Original_CH4_mean = mean(Ag_CH4_CO2e, na.rm=TRUE),
  NoRice_CH4_mean = mean(Ag_CH4_CO2e_noRice, na.rm=TRUE),
  Original_total_mean = mean(Ag_net_CH4_N2O, na.rm=TRUE),
  NoRice_total_mean = mean(Ag_noRice_CH4_N2O, na.rm=TRUE)
)])

# Save updated data
fwrite(dt_v2, "/Users/savannah/Documents/Ch4/Data/Ag_opportunity_cost/__dt_Ag_CH4_N2O_xy_noRice.csv")


dt_v2<-fread("/Users/savannah/Documents/Ch4/Data/Ag_opportunity_cost/__dt_Ag_CH4_N2O_xy_noRice.csv")
#dt$Ag_net_CH4_N2O[dt$Ag_net_CH4_N2O < 0.1] <- NA # adjust data to only compute mean for values > 0.3 Mg ha-1 yr-1

dt_v2 %>% group_by(Biome) %>%
  summarise(Ag_CH4 = mean(Ag_CH4_CO2e_noRice, na.rm=T),
            Ag_N2O = mean(Ag_N2O_CO2e, na.rm=T),                
            Ag_CH4_N2O = mean(Ag_noRice_CH4_N2O, na.rm=T))

# Test: Only keep Croplands: at least 60% of area is cultivated cropland. 
dt_crops<-subset(dt_v2, LC_num ==12) # Rice CH4 already excluded (set CH4 fluxes for rice cultivation areas to 0 when generating __dt_Ag_CH4_N2O_xy_noRice.csv)

Ag_crops<- dt_crops %>% group_by(Biome) %>%
  summarise(Ag_CH4 = mean(Ag_CH4_CO2e_noRice, na.rm=T),
            Ag_N2O = mean(Ag_N2O_CO2e, na.rm=T),                
            Ag_CH4_N2O = mean(Ag_noRice_CH4_N2O, na.rm=T))
dt<-dt_v2
# Remove biomes not in scope of our study
dt<-subset(dt, Biome_num != 8); dt<-subset(dt, Biome_num != 9); dt<-subset(dt, Biome_num != 10); dt<-subset(dt, Biome_num != 11)
dt<-subset(dt, Biome_num != 13); dt<-subset(dt, Biome_num != 14); 
dt<- subset(dt, Biome!= "") # Remove NA Biome row  
Ag_s<- dt %>% group_by(Biome) %>%
  summarise(Ag_CH4 = mean(Ag_CH4_CO2e_noRice, na.rm=T),
            Ag_N2O = mean(Ag_N2O_CO2e, na.rm=T),
            Ag_CH4_N2O = mean(Ag_noRice_CH4_N2O, na.rm=T)) # , #


# Adjust to match Eco dt: Biome, B0, B0.SE, GHG
tmp<- reshape2::melt(Ag_s, id.vars = "Biome" )
colnames(tmp)[2:3]<- c("GHG", "B0") 

# This approach results in huge SD values bc of sumNovizing across space
# Ag_s<- dt %>% group_by(Biome) %>%
#   summarise(Ag_CH4 = sd(Ag_CH4_CO2e, na.rm=T),
#             Ag_N2O = sd(Ag_N2O_CO2e, na.rm=T),
#             Ag_CH4_N2O = sd(Ag_net_CH4_N2O, na.rm=T))

# # This approach results in smaller SD values bc of sumNovizing across time (5 yrs)
Ag_s<- dt %>% group_by(Biome) %>%
  summarise(Ag_CH4 = mean(Ag_CH4_CO2e.sd, na.rm=T),
            Ag_N2O = mean(Ag_N2O_CO2e.sd, na.rm=T),
            Ag_CH4_N2O = mean(sd_Ag_net_CH4_N2O, na.rm=T)) # , #

t<- cbind(tmp, reshape2::melt(Ag_s, id.vars = "Biome" ) )
t<- subset(t, select=c("Biome","B0", "value", "GHG"  )); colnames(t)[3]<-"B0.SE" # Actually SD. Using "B0.SE" to match other df s

# Add Global mean rows
g<- t %>%
  group_by(GHG) %>%
  summarise(Ag_B0.mean = mean(B0, na.rm=T),
            Ag_B0.SE.mean = mean(B0.SE, na.rm=T))
g$Biome<- "Global mean"; names(g)[2:3]<-c("B0", "B0.SE")
g<-subset(g, select=c("Biome","B0", "B0.SE", "GHG"))
t<-rbind(t,g)
fwrite(t,"/Users/savannah/Documents/Ch4/Data/Ag_opportunity_cost/__summ_Ag_CH4_N2O_xy_noRice.csv")


#####################################################
# Helper Functions
#####################################################
# Need to convert to all positive values for gamma regression model - used dt_Mar2024 version to be consistent with past results. 
#dt$CH4.CO2e_20yr_conv_pos<- dt$CH4.CO2e_20yr + abs(min(dt$CH4.CO2e_20yr, na.rm = T)) +1.0001 
#dt$N2O.CO2e_conv_pos<- dt$N2O.CO2e + abs(min(dt$N2O.CO2e, na.rm = T)) +1.0001

# Helper functions with explicit error propagation for back-transformation
back_transform_ch4 <- function(value) {
  # For CH4: exp(x) - (78.67489 + 1.0001)
  return(exp(value) - 79.67499)
}
back_transform_ch4_se <- function(estimate, se) {
  # Using delta method for error propagation
  # If Y = exp(X) - c, then Var(Y) = exp(2*mu_X) * Var(X)
  mean_transformed <- exp(estimate)
  se_transformed <- mean_transformed * se
  return(se_transformed)
}

back_transform_n2o <- function(value) {
  # For N2O: exp(x) - (1.358029 + 1.0001)
  return(exp(value) - 2.358129)
}
back_transform_n2o_se <- function(estimate, se) {
  # Similar error propagation for N2O
  mean_transformed <- exp(estimate)
  se_transformed <- mean_transformed * se
  return(se_transformed)
}

#####################################################
# Run gamma regression for each ghg-biome pair
#####################################################
# vNov 2024: updated data set to exclude peatland ecosystems (142 studies to 128 studies) 
# Also corrected Citation #19: http://dx.doi.org/10.1023/A:1015856617553. Old extracted data included only CH4 oxidation from table NOT full fluxes from figure
#dt<-fread(paste(root,"__Agg_Nov2024_drainedPeat.csv", sep="")) 
dt<-fread(paste(root,"__Agg_Nov2024_noPeatland.csv", sep="")) 

# Re-do Boreal flux obs. vNov2024: Going to assume growing season is half the year not 1/3. vMar2024 already multiplied by 0.33 in 2_prep_stat_mod.R. 
dt$CH4.C.converted[dt$Biome == "Boreal forest"] <- dt$CH4.C.converted[dt$Biome == "Boreal forest"]*0.5/0.33
dt$N2O.N.converted[dt$Biome == "Boreal forest"] <- dt$N2O.N.converted[dt$Biome == "Boreal forest"]*0.5/0.33

dt$CH4.CO2e_20yr<- dt$CH4.C.converted*14/12*85*1e4*1e-9 # mg CH4-C m-2 y-1 * 85 CO2e/CH4 on 20yr time horizon * 1e4 m2 ha-1 * 1e-9 Mg mg-1
dt$N2O.CO2e_20yr<- dt$N2O.N.converted*44/28*273*1e4*1e-9 # mg N2O-N m-2 y-1 * 273 CO2e/N2O on 20yr time horizon * 1e4 m2 ha-1* 1e-9 Mg mg-1

# Need to convert to all positive values for gamma regression model
dt$CH4.CO2e_20yr_conv_pos<- dt$CH4.CO2e_20yr + abs(-78.67489) +1.0001 # 1.0001 instead of 0.0001 for family = Gamma(link = "log"), default is Gamma(link = "inverse") but this did not work for my data
dt$N2O.CO2e_20yr_conv_pos<- dt$N2O.CO2e_20yr + abs(-1.358029) +1.0001

dt$Biome<-as.factor(dt$Biome); levels(dt$Biome);  
dt$Previous.LC<-as.factor(dt$Previous.LC); levels(dt$Previous.LC)

# Remove citation 174 bc suspiciously low CH4, DOI: https://doi.org/10.1016/j.ejsobi.2017.08.004
dt <-subset(dt, Citation != 174); dt$Citation<-as.factor(dt$Citation); length(levels(dt$Citation)); # 114 unique studies (no peatland included)
#summary(dt$Age); dt %>% plot_frq(Biome); dt %>% plot_frq(Chronosequ); 

dt %>%
  ggplot(aes(x = CH4.CO2e_20yr_conv_pos, fill = Biome)) +
  geom_histogram(binwidth = 0.1, position = "stack", alpha = 0.7) +  # Adjust binwidth as needed
  scale_fill_viridis_d() +  # You can use other color palettes
  labs(x = expression("CH"[4]*"-C soil flux (mg" ~ m^-2 ~ yr^-1*")") ,
       y = "Count") +
  theme_minimal()

dt %>%
  group_by(Biome) %>%
  summarise(freq = sum(CH4.CO2e_20yr_conv_pos, na.rm = TRUE)) %>%
  ggplot(aes(x = Biome, y = freq, fill = Biome)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d() +  # Or any other color scale you prefer
  labs(x = "Biome", y = "Frequency", 
       title = expression("CH"[4]*"-C soil flux (mg" ~ m^-2 ~ yr^-1*")") )
       #title = expression("N"[2]*"O-N soil flux (mg" ~ m^-2 ~ yr^-1*")") ) #+theme_minimal()

# Create a list to store standardization parameters for all biomes
standardization_params <- list()

#####     v4: run models & generate Table S3 (GWP 20 yrs)     ####   

# Modified code to create Table S3

# Function to extract model diagnostics for normality and heteroscedasticity
extract_model_diagnostics <- function(model) {
  # Get residuals
  resids <- residuals(model)
  
  # Shapiro-Wilk test for normality
  sw_test <- tryCatch(
    shapiro.test(resids),
    error = function(e) list(p.value = NA)
  )
  
  # Breusch-Pagan test for heteroscedasticity using lmtest package
  # If not already installed, use: install.packages("lmtest")
  bp_test <- tryCatch({
    require(lmtest)
    # Create a dataset with residuals and fitted values
    df_resid <- data.frame(resid = resids^2, fitted = fitted(model))
    # Run linear model on squared residuals
    lm_resid <- lm(resid ~ fitted, data = df_resid)
    # Get p-value from F-test
    anova(lm_resid)$"Pr(>F)"[1]
  }, error = function(e) NA)
  
  # Calculate AIC
  aic_val <- AIC(model)
  
  return(list(
    shapiro_p = sw_test$p.value,
    bp_p = bp_test,
    aic = aic_val
  ))
}

# Function to identify covariate names from model formula
get_covariate_names <- function(model) {
  # Extract formula terms
  terms <- attr(terms(model), "term.labels")
  # Remove random effects term if present
  terms <- terms[!grepl("\\|", terms)]
  # Return only the covariates (excluding Age which is always first)
  if(length(terms) > 1) {
    return(terms[2:length(terms)])
  } else {
    return(character(0))
  }
}

# UPDATED version
create_model_summary_table <- function(gas_type) {
  if(gas_type == "CH4") {
    biomes <- levels(dt$Biome)
    all_results <- data.frame()
    
    for(biome in biomes) {
      print(paste("Processing", biome, "for CH4"))
      
      # Subset data
      sub_dt <- subset(dt, select = c(precip_seasn, perc_clay, aspect, slope, perc_sand,
                                      perc_silt, temp_seasn, Age, elevation, MAT, MAP))
      sub_dt <- as.data.frame(sapply(sub_dt, as.numeric))
      sub_dt <- as.data.frame(scale(sub_dt, center = TRUE, scale = TRUE))
      
      # Store original Age values before scaling for SD calculation
      age_orig <- dt$Age[dt$Biome == biome]
      age_sd <- sd(age_orig, na.rm = TRUE)
      print(paste("Age SD for", biome, ":", age_sd))
      
      # Use CH4.CO2e_20yr_conv_pos or CH4.CO2e_100yr_conv_pos based on what's being analyzed
      # if("CH4.CO2e_100yr_conv_pos" %in% names(dt) ) { # && any(grepl("100yr", deparse(sys.call(-1))))) {
      #   sub_dt$CH4.CO2e_conv_pos <- dt$CH4.CO2e_100yr_conv_pos
      #   print("Using 100yr GWP for CH4")
      # } else {
      #   sub_dt$CH4.CO2e_conv_pos <- dt$CH4.CO2e_20yr_conv_pos
      #   print("Using 20yr GWP for CH4")
      # }
      sub_dt$CH4.CO2e_conv_pos <- dt$CH4.CO2e_20yr_conv_pos
      print("Using 20yr GWP for CH4")
      # sub_dt$CH4.CO2e_conv_pos <- dt$CH4.CO2e_100yr_conv_pos
      # print("Using 100yr GWP for CH4")
      
      sub_dt$Citation <- dt$Citation
      sub_dt$Citation <- as.factor(sub_dt$Citation)
      sub_dt$Biome <- dt$Biome
      sub_dt$Age_orig <- dt$Age
      sub_dt <- subset(sub_dt, Biome == biome)
      
      # Refit models based on biome type (using the same formulas as in your original code)
      # Note: m0 models are removed
      if(biome == "Boreal forest") {
        m1 <- glmer(CH4.CO2e_conv_pos ~ Age + (1|Citation) + slope, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1000)))
        m2 <- glmer(CH4.CO2e_conv_pos ~ Age + (1|Citation) + slope + MAT, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
        m3 <- glmer(CH4.CO2e_conv_pos ~ Age + (1|Citation) + slope + MAT + perc_sand, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
      } else if(biome == "Subtropical/tropical forest") {
        m1 <- glmer(CH4.CO2e_conv_pos ~ Age + (1|Citation) + slope, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1000)))
        m2 <- glmer(CH4.CO2e_conv_pos ~ Age + (1|Citation) + slope + perc_sand, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1000)))
        m3 <- glmer(CH4.CO2e_conv_pos ~ Age + (1|Citation) + slope + perc_sand + aspect, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
      } else if(biome == "Subtropical/tropical savanna") {
        m1 <- glmer(CH4.CO2e_conv_pos ~ Age + (1|Citation) + temp_seasn, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1000)))
        m2 <- glmer(CH4.CO2e_conv_pos ~ Age + (1|Citation) + temp_seasn + elevation, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
        m3 <- glmer(CH4.CO2e_conv_pos ~ Age + (1|Citation) + temp_seasn + elevation + precip_seasn, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
      } else if(biome == "Temperate broadleaf forest") {
        m1 <- glmer(CH4.CO2e_conv_pos ~ Age + (1|Citation) + temp_seasn, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
        m2 <- glmer(CH4.CO2e_conv_pos ~ Age + (1|Citation) + temp_seasn + MAT, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
        m3 <- glmer(CH4.CO2e_conv_pos ~ Age + (1|Citation) + temp_seasn + MAT + elevation, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
      } else { # Temperate conifer forest
        m1 <- glmer(CH4.CO2e_conv_pos ~ Age + (1|Citation) + perc_clay, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
        m2 <- glmer(CH4.CO2e_conv_pos ~ Age + (1|Citation) + perc_clay + elevation, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
        m3 <- glmer(CH4.CO2e_conv_pos ~ Age + (1|Citation) + perc_clay + elevation + perc_sand, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
      }
      
      # Process each model (only m1, m2, m3)
      for(model_num in 1:3) {
        model <- get(paste0("m", model_num))
        
        # Get model summary
        model_summary <- summary(model)
        
        # Get model diagnostics
        diag <- extract_model_diagnostics(model)
        
        # Get covariate names
        cov_names <- get_covariate_names(model)
        cov_names_full <- c("Age", cov_names)
        
        # Extract coefficients
        coefs <- fixef(model)
        coef_names <- names(coefs)
        
        # Extract standard errors
        se <- sqrt(diag(vcov(model)))
        
        # Create result row for this model
        result <- data.frame(
          GHG = "CH4",
          Biome = biome,
          Model = paste0("m", model_num),
          N_obs = nobs(model),
          Num_Cit = model_summary$ngrps,
          AIC = diag$aic,
          LogLik = logLik(model)[1],
          Shapiro_p = diag$shapiro_p,
          BP_p = diag$bp_p
        )
        
        # Add intercept
        result$Intercept <- coefs["(Intercept)"]
        result$SE_Intercept <- se["(Intercept)"]
        
        print(result)
        #Done later in Fig2 prep function: 
        result$Intercept_transformed <- back_transform_ch4(result$Intercept) 
        result$SE_Intercept_transformed <- back_transform_ch4_se(result$Intercept, result$SE_Intercept) 
        
        
        # Add Age coefficient (back-transformed)
        age_coef <- coefs["Age"]
        age_se <- se["Age"]
        
        # Back transform Age coefficient for interpretability
        # For gamma regression with log link, if the original data is on a log scale:
        # 1. The coefficient represents the expected change in log(y) for a one-unit change in x
        # 2. To get the multiplicative effect: exp(coef)
        # 3. For percentage change: (exp(coef) - 1) * 100
        # 4. Since Age was scaled, we need to divide by the SD to get the effect per original unit
        
        # Compute backtransformed coefficients directly
        if(!is.na(age_coef) && !is.na(age_sd) && age_sd > 0) {
          # For CH4, implement the backtransformation based on gamma regression with log link
          # and scaled predictor (Age)
          age_coef_orig <- (exp(age_coef)) / age_sd # (exp(age_coef) - 1) / age_sd: Change in the response per unit of original Age (as a proportion)
          # while exp(age_coef) / age_sd: Multiplicative factor in the response per unit of original Age
          
          age_se_orig <- exp(age_se) / age_sd
        } else {
          age_coef_orig <- NA
          age_se_orig <- NA
        }
        
        result$Age_coeff <- age_coef
        result$SE_Age_coeff <- age_se
        result$Age_coeff_transformed <- age_coef_orig
        result$SE_Age_coeff_transformed <- age_se_orig
        result$p_Age <- model_summary$coefficients["Age", "Pr(>|z|)"]
        
        # Add other coefficients if they exist
        if(length(cov_names) >= 1) {
          result$B2_name <- cov_names[1]
          result$B2_coeff <- if(cov_names[1] %in% coef_names) coefs[cov_names[1]] else NA
          result$SE_B2_coeff <- if(cov_names[1] %in% coef_names) se[cov_names[1]] else NA
          result$p_B2 <- if(cov_names[1] %in% coef_names) model_summary$coefficients[cov_names[1], "Pr(>|z|)"] else NA
        } else {
          result$B2_name <- NA
          result$B2_coeff <- NA
          result$SE_B2_coeff <- NA
          result$p_B2 <- NA
        }
        
        if(length(cov_names) >= 2) {
          result$B3_name <- cov_names[2]
          result$B3_coeff <- if(cov_names[2] %in% coef_names) coefs[cov_names[2]] else NA
          result$SE_B3_coeff <- if(cov_names[2] %in% coef_names) se[cov_names[2]] else NA
          result$p_B3 <- if(cov_names[2] %in% coef_names) model_summary$coefficients[cov_names[2], "Pr(>|z|)"] else NA
        } else {
          result$B3_name <- NA
          result$B3_coeff <- NA
          result$SE_B3_coeff <- NA
          result$p_B3 <- NA
        }
        
        if(length(cov_names) >= 3) {
          result$B4_name <- cov_names[3]
          result$B4_coeff <- if(cov_names[3] %in% coef_names) coefs[cov_names[3]] else NA
          result$SE_B4_coeff <- if(cov_names[3] %in% coef_names) se[cov_names[3]] else NA
          result$p_B4 <- if(cov_names[3] %in% coef_names) model_summary$coefficients[cov_names[3], "Pr(>|z|)"] else NA
        } else {
          result$B4_name <- NA
          result$B4_coeff <- NA
          result$SE_B4_coeff <- NA
          result$p_B4 <- NA
        }
        
        all_results <- rbind(all_results, result)
      }
    }
    
    return(all_results)
  } else if(gas_type == "N2O") {
    biomes <- levels(dt$Biome)
    all_results <- data.frame()
    
    for(biome in biomes) {
      print(paste("Processing", biome, "for N2O"))
      
      # Subset data
      sub_dt <- subset(dt, select = c(MAT, perc_sand, slope, precip_seasn, aspect, Age,
                                      soil_porosity, temp_seasn, perc_silt, MAP, elevation, perc_clay))
      sub_dt <- as.data.frame(sapply(sub_dt, as.numeric))
      sub_dt <- as.data.frame(scale(sub_dt, center = TRUE, scale = TRUE))
      
      # Store original Age values before scaling for SD calculation
      age_orig <- dt$Age[dt$Biome == biome]
      age_sd <- sd(age_orig, na.rm = TRUE)
      print(paste("Age SD for", biome, ":", age_sd))
      
      # Use the correct variable name based on whether we're analyzing 20yr or 100yr GWP
      # if("N2O.CO2e_100yr_conv_pos" %in% names(dt)  ) { # && any(grepl("100yr", deparse(sys.call(-1))))) {
      #   sub_dt$N2O.CO2e_conv_pos <- dt$N2O.CO2e_100yr_conv_pos
      #   print("Using 100yr GWP for N2O")
      # } else {
      #   sub_dt$N2O.CO2e_conv_pos <- dt$N2O.CO2e_20yr_conv_pos
      #   print("Using 20yr GWP for N2O")
      # }
      sub_dt$N2O.CO2e_conv_pos <- dt$N2O.CO2e_20yr_conv_pos
      print("Using 20yr GWP for N2O")
      # sub_dt$N2O.CO2e_conv_pos <- dt$N2O.CO2e_100yr_conv_pos
      # print("Using 100yr GWP for N2O")
      
      sub_dt$Citation <- dt$Citation
      sub_dt$Citation <- as.factor(sub_dt$Citation)
      sub_dt$Biome <- dt$Biome
      sub_dt$Age_orig <- dt$Age
      sub_dt <- subset(sub_dt, Biome == biome)
      
      # Refit models based on biome type (removing m0 models)
      if(biome == "Boreal forest") {
        m1 <- glmer(N2O.CO2e_conv_pos ~ Age + (1|Citation) + aspect, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
        m2 <- glmer(N2O.CO2e_conv_pos ~ Age + (1|Citation) + aspect + perc_sand, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
        m3 <- glmer(N2O.CO2e_conv_pos ~ Age + (1|Citation) + aspect + perc_sand + MAP, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
      } else if(biome == "Subtropical/tropical forest") {
        m1 <- glmer(N2O.CO2e_conv_pos ~ Age + (1|Citation) + slope, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1000)))
        m2 <- glmer(N2O.CO2e_conv_pos ~ Age + (1|Citation) + slope + precip_seasn, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1000)))
        m3 <- glmer(N2O.CO2e_conv_pos ~ Age + (1|Citation) + slope + precip_seasn + perc_silt, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
      } else if(biome == "Subtropical/tropical savanna") {
        m1 <- glmer(N2O.CO2e_conv_pos ~ Age + (1|Citation) + temp_seasn, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
        m2 <- glmer(N2O.CO2e_conv_pos ~ Age + (1|Citation) + temp_seasn + precip_seasn, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
        m3 <- glmer(N2O.CO2e_conv_pos ~ Age + (1|Citation) + temp_seasn + precip_seasn + elevation, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
      } else if(biome == "Temperate broadleaf forest") {
        m1 <- glmer(N2O.CO2e_conv_pos ~ Age + (1|Citation) + MAT, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1000)))
        m2 <- glmer(N2O.CO2e_conv_pos ~ Age + (1|Citation) + MAT + aspect, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1000)))
        m3 <- glmer(N2O.CO2e_conv_pos ~ Age + (1|Citation) + MAT + aspect + perc_sand, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
      } else { # Temperate conifer forest
        m1 <- glmer(N2O.CO2e_conv_pos ~ Age + (1|Citation) + elevation, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
        m2 <- glmer(N2O.CO2e_conv_pos ~ Age + (1|Citation) + elevation + perc_clay, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
        m3 <- glmer(N2O.CO2e_conv_pos ~ Age + (1|Citation) + elevation + perc_clay + soil_porosity, data = sub_dt,
                    family = Gamma(link = "log"), na.action = na.exclude, nAGQ = 0,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
      }
      
      # Process each model (only m1, m2, m3)
      for(model_num in 1:3) {
        model <- get(paste0("m", model_num))
        
        # Get model summary
        model_summary <- summary(model)
        
        # Get model diagnostics
        diag <- extract_model_diagnostics(model)
        
        # Get covariate names
        cov_names <- get_covariate_names(model)
        cov_names_full <- c("Age", cov_names)
        
        # Extract coefficients
        coefs <- fixef(model)
        coef_names <- names(coefs)
        
        # Extract standard errors
        se <- sqrt(diag(vcov(model)))
        
        # Create result row for this model
        result <- data.frame(
          GHG = "N2O",
          Biome = biome,
          Model = paste0("m", model_num),
          N_obs = nobs(model),
          Num_Cit = model_summary$ngrps,
          AIC = diag$aic,
          LogLik = logLik(model)[1],
          Shapiro_p = diag$shapiro_p,
          BP_p = diag$bp_p
        )
        
        # Add intercept
        result$Intercept <- coefs["(Intercept)"]
        result$SE_Intercept <- se["(Intercept)"]
        #Done later in Fig2 prep function: 
        result$Intercept_transformed <- back_transform_n2o(result$Intercept) 
        result$SE_Intercept_transformed <- back_transform_n2o_se(result$Intercept, result$SE_Intercept) 
        
        # Add Age coefficient (back-transformed)
        age_coef <- coefs["Age"]
        age_se <- se["Age"]
        
        # Back transform Age coefficient for interpretability using the Excel formula approach:
        # Age_backtr = (EXP(age_coef))/1SD_Age_biome
        if(!is.na(age_coef) && !is.na(age_sd) && age_sd > 0) {
          # For N2O, implement the backtransformation based on gamma regression with log link
          # and scaled predictor (Age) 
          age_coef_orig <- (exp(age_coef)) / age_sd # (exp(age_coef) - 1) / age_sd: Change in the response per unit of original Age (as a proportion)
          # while exp(age_coef) / age_sd: Multiplicative factor in the response per unit of original Age
          age_se_orig <- exp(age_se) / age_sd
        } else {
          age_coef_orig <- NA
          age_se_orig <- NA
        }
        
        result$Age_coeff <- age_coef
        result$SE_Age_coeff <- age_se
        result$Age_coeff_transformed <- age_coef_orig  # Use calculated value instead of NA
        result$SE_Age_coeff_transformed <- age_se_orig  # Use calculated value instead of NA
        result$p_Age <- model_summary$coefficients["Age", "Pr(>|z|)"]
        
        # Add other coefficients if they exist
        if(length(cov_names) >= 1) {
          result$B2_name <- cov_names[1]
          result$B2_coeff <- if(cov_names[1] %in% coef_names) coefs[cov_names[1]] else NA
          result$SE_B2_coeff <- if(cov_names[1] %in% coef_names) se[cov_names[1]] else NA
          result$p_B2 <- if(cov_names[1] %in% coef_names) model_summary$coefficients[cov_names[1], "Pr(>|z|)"] else NA
        } else {
          result$B2_name <- NA
          result$B2_coeff <- NA
          result$SE_B2_coeff <- NA
          result$p_B2 <- NA
        }
        
        if(length(cov_names) >= 2) {
          result$B3_name <- cov_names[2]
          result$B3_coeff <- if(cov_names[2] %in% coef_names) coefs[cov_names[2]] else NA
          result$SE_B3_coeff <- if(cov_names[2] %in% coef_names) se[cov_names[2]] else NA
          result$p_B3 <- if(cov_names[2] %in% coef_names) model_summary$coefficients[cov_names[2], "Pr(>|z|)"] else NA
        } else {
          result$B3_name <- NA
          result$B3_coeff <- NA
          result$SE_B3_coeff <- NA
          result$p_B3 <- NA
        }
        
        if(length(cov_names) >= 3) {
          result$B4_name <- cov_names[3]
          result$B4_coeff <- if(cov_names[3] %in% coef_names) coefs[cov_names[3]] else NA
          result$SE_B4_coeff <- if(cov_names[3] %in% coef_names) se[cov_names[3]] else NA
          result$p_B4 <- if(cov_names[3] %in% coef_names) model_summary$coefficients[cov_names[3], "Pr(>|z|)"] else NA
        } else {
          result$B4_name <- NA
          result$B4_coeff <- NA
          result$SE_B4_coeff <- NA
          result$p_B4 <- NA
        }
        
        all_results <- rbind(all_results, result)
      }
    }
    return(all_results)
  }
}

# Create the summary tables
ch4_results <- create_model_summary_table("CH4")
n2o_results <- create_model_summary_table("N2O")

# Combine results
table_s3 <- rbind(ch4_results, n2o_results)

# Write the complete table to a CSV file
write.csv(table_s3, paste(root, "Table_S3_GHG_reg_models_20yr_noPeat_vMarch2025.csv", sep = ""), row.names = FALSE)

#####     v4: Summarize Table S3 + README     ####   
# 1. Update the table_s3_formatted to only include the models used for figures 2 and 3
# These are the best performing models with lowest AIC among m1, m2, and m3

# Make sure dplyr is loaded
library(dplyr)

# Define which models to use for each ghg-biome combination
best_m1_m3_models <- data.frame(
  GHG = c("CH4", "CH4", "CH4", "CH4", "CH4", 
          "N2O", "N2O", "N2O", "N2O", "N2O"),
  Biome = c("Boreal forest", "Subtropical/tropical forest", "Subtropical/tropical savanna", 
            "Temperate broadleaf forest", "Temperate conifer forest",
            "Boreal forest", "Subtropical/tropical forest", "Subtropical/tropical savanna", 
            "Temperate broadleaf forest", "Temperate conifer forest"),
  Model = c("m1", "m2", "m3", "m1", "m1",
            "m1", "m1", "m3", "m1", "m3")
)

# Replace the existing table_s3_formatted creation with this:
# Use base R instead of dplyr for safer compatibility
table_s3_formatted <- merge(table_s3, best_m1_m3_models, by = c("GHG", "Biome", "Model"))
table_s3_formatted <- table_s3_formatted[, c("GHG", "Biome", "Model", "N_obs", "Num_Cit", "AIC", 
                                             "LogLik", "Shapiro_p", "BP_p", "Age_coeff_transformed", 
                                             "SE_Age_coeff_transformed", "p_Age", "B2_name", "B2_coeff", 
                                             "SE_B2_coeff", "p_B2", "B3_name", "B3_coeff", "SE_B3_coeff", 
                                             "p_B3", "B4_name", "B4_coeff", "SE_B4_coeff", "p_B4")]

# Write the formatted table to a CSV file
write.csv(table_s3_formatted, paste(figsDir, "Table_S3_GHG_reg_models_20yr_manuscript.csv", sep = ""), row.names = FALSE)



# Create a comprehensive README for Table S3
table_s3_readme <- "
Table S3. Results of the mixed-effect gamma regression models showing the biome-specific relationships 
between soil greenhouse gas fluxes of (A) CH4 and (B) N2O and covariates for Age (B_1 coefficient) 
as well as environmental covariates (B_2, B_3 and B_4 coefficients). 

The table includes:
- Model complexity comparison (m0: Age only, m1: Age + 1 covariate, m2: Age + 2 covariates, m3: Age + 3 covariates)
- Covariates represented by coefficients B_2, B_3 and B_4 vary in each model based on biome-specific variable importance
- All coefficient estimates were calculated with standardized variables
- Age coefficient has been back-transformed for direct interpretation as change in GHG flux per year
- Model evaluation metrics: Number of plot-year observations (N_obs), number of studies (Num_Cit)
- Model fit metrics: Akaike's Information Criteria (AIC) and log likelihood (LogLik)
- Model validation tests:
  - Normality of residuals (Shapiro-Wilk test, pass if p > 0.05) 
  - Homoscedasticity (Breusch-Pagan test, pass if p > 0.05)

Significance codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
"

writeLines(table_s3_readme, paste(figsDir, "Table_S3_README.txt", sep = ""))

#####     v4: GWP 100 yrs for CH4     ####   

# New version of Table S3 needs to based on results with a 100 year GWP 
dt$CH4.CO2e_100yr<- dt$CH4.C.converted*14/12*30*1e4*1e-9 # mg CH4-C m-2 y-1 * 30 CO2e/CH4 on 100 year time horizon * 1e4 m2 ha-1 * 1e-9 Mg mg-1
dt$N2O.CO2e_100yr<- dt$N2O.N.converted*44/28*298*1e4*1e-9 # mg N2O-N m-2 y-1 * 298 CO2e/N2O on a 20 and 100 year time horizon (298 multiplier is the same as 20yr) * 1e4 m2 ha-1* 1e-9 Mg mg-1

# Need to convert to all positive values for gamma regression model
# No need to update these values bc absolute value of new min 100 yr GWP are smaller than abs values of 20 yr GWP
# > min(dt$N2O.CO2e_100yr, na.rm = T) [1] -1.358029
# > min(dt$CH4.CO2e_100yr, na.rm = T) [1] -2.478
dt$CH4.CO2e_100yr_conv_pos<- dt$CH4.CO2e_100yr + abs(-78.67489) +1.0001 
dt$N2O.CO2e_100yr_conv_pos<- dt$N2O.CO2e_100yr + abs(-1.358029) +1.0001 

# Reproduce the regression analysis with GWP 100 yrs instead of 20 yrs
# Function to run regression models with GWP 100 yrs

# Modified run_gwp100_analysis function to use the updated create_model_summary_table
run_gwp100_analysis <- function() {
  # Make sure the 100yr variables exist
  if(!"CH4.CO2e_100yr_conv_pos" %in% names(dt)) {
    dt$N2O.CO2e_100yr_conv_pos <- dt$N2O.CO2e_100yr + abs(-1.358029) + 1.0001
    print("Created N2O.CO2e_100yr_conv_pos variable")
  }
  
  # Run the analyses using the updated function
  print("Running CH4 analysis with 100yr GWP...")
  ch4_results_100yr <- create_model_summary_table("CH4")
  
  print("Running N2O analysis with 100yr GWP...")
  n2o_results_100yr <- create_model_summary_table("N2O")
  
  # Combine results
  table_s3_100yr <- rbind(ch4_results_100yr, n2o_results_100yr)
  
  # Write the complete table to a CSV file
  write.csv(table_s3_100yr, paste(root, "Table_S3_GHG_reg_models_100yr_noPeat_vMarch2025.csv", sep = ""), row.names = FALSE)
  
  # Create the formatted summary table for the best models (m1-m3) using base R
  table_s3_100yr_formatted <- merge(table_s3_100yr, best_m1_m3_models, by = c("GHG", "Biome", "Model"))
  table_s3_100yr_formatted <- table_s3_100yr_formatted[, c("GHG", "Biome", "Model", "N_obs", "Num_Cit", "AIC",
                                                           "LogLik", "Shapiro_p", "BP_p", 
                                                           "Intercept_transformed","SE_Intercept_transformed","Age_coeff_transformed",
                                                           "SE_Age_coeff_transformed", "p_Age", "B2_name", "B2_coeff",
                                                           "SE_B2_coeff", "p_B2", "B3_name", "B3_coeff", "SE_B3_coeff",
                                                           "p_B3", "B4_name", "B4_coeff", "SE_B4_coeff", "p_B4")]
  
  # Format p-values with stars for significance using base R
  table_s3_100yr_formatted$Age_sig <- ""
  table_s3_100yr_formatted$Age_sig[table_s3_100yr_formatted$p_Age < 0.1] <- "."
  table_s3_100yr_formatted$Age_sig[table_s3_100yr_formatted$p_Age < 0.05] <- "*"
  table_s3_100yr_formatted$Age_sig[table_s3_100yr_formatted$p_Age < 0.01] <- "**"
  table_s3_100yr_formatted$Age_sig[table_s3_100yr_formatted$p_Age < 0.001] <- "***"
  
  table_s3_100yr_formatted$B2_sig <- ""
  table_s3_100yr_formatted$B2_sig[table_s3_100yr_formatted$p_B2 < 0.1] <- "."
  table_s3_100yr_formatted$B2_sig[table_s3_100yr_formatted$p_B2 < 0.05] <- "*"
  table_s3_100yr_formatted$B2_sig[table_s3_100yr_formatted$p_B2 < 0.01] <- "**"
  table_s3_100yr_formatted$B2_sig[table_s3_100yr_formatted$p_B2 < 0.001] <- "***"
  
  table_s3_100yr_formatted$B3_sig <- ""
  table_s3_100yr_formatted$B3_sig[table_s3_100yr_formatted$p_B3 < 0.1] <- "."
  table_s3_100yr_formatted$B3_sig[table_s3_100yr_formatted$p_B3 < 0.05] <- "*"
  table_s3_100yr_formatted$B3_sig[table_s3_100yr_formatted$p_B3 < 0.01] <- "**"
  table_s3_100yr_formatted$B3_sig[table_s3_100yr_formatted$p_B3 < 0.001] <- "***"
  
  table_s3_100yr_formatted$B4_sig <- ""
  table_s3_100yr_formatted$B4_sig[table_s3_100yr_formatted$p_B4 < 0.1] <- "."
  table_s3_100yr_formatted$B4_sig[table_s3_100yr_formatted$p_B4 < 0.05] <- "*"
  table_s3_100yr_formatted$B4_sig[table_s3_100yr_formatted$p_B4 < 0.01] <- "**"
  table_s3_100yr_formatted$B4_sig[table_s3_100yr_formatted$p_B4 < 0.001] <- "***"
  
  table_s3_100yr_formatted$Norm_test <- ifelse(table_s3_100yr_formatted$Shapiro_p > 0.05, "Pass", "Fail")
  table_s3_100yr_formatted$Hetero_test <- ifelse(table_s3_100yr_formatted$BP_p > 0.05, "Pass", "Fail")
  
  # Write the formatted 100yr table to a CSV file
  write.csv(table_s3_100yr_formatted, paste(figsDir, "Table_S3_GHG_reg_models_100yr_manuscript.csv", sep = ""), row.names = FALSE)
  
  return(table_s3_100yr)
  
}

# Run the GWP 100yr analysis
table_s3_100yr <- run_gwp100_analysis()

table_s3_formatted <- merge(table_s3_100yr, best_m1_m3_models, by = c("GHG", "Biome", "Model"))
table_s3_formatted <- table_s3_formatted[, c("GHG", "Biome", "Model", "N_obs", "Num_Cit", "AIC", 
                                             "LogLik", "Shapiro_p", "BP_p", "Intercept_transformed","SE_Intercept_transformed",
                                             "Age_coeff_transformed", 
                                             "SE_Age_coeff_transformed", "p_Age", "B2_name", "B2_coeff", 
                                             "SE_B2_coeff", "p_B2", "B3_name", "B3_coeff", "SE_B3_coeff", 
                                             "p_B3", "B4_name", "B4_coeff", "SE_B4_coeff", "p_B4")]

# Write the formatted table to a CSV file
write.csv(table_s3_formatted, paste(figsDir, "Table_S3_GHG_reg_models_100yr_manuscript.csv", sep = ""), row.names = FALSE)



###################################################### 
# Functions for Fig 2: Prepare GHG data, plot results
######################################################

prepare_ghg_data <- function(model_data, ag_file) {
  # Hardcoded best models based on AIC values
  ch4_best_models <- data.table(
    Biome = c("Boreal forest", "Subtropical/tropical forest", "Subtropical/tropical savanna", 
              "Temperate broadleaf forest", "Temperate conifer forest"),
    Model = c("m1", "m2", "m3", "m1", "m1") # if considering m0: ("m0", "m2", "m3", "m0", "m1")
  )
  
  n2o_best_models <- data.table(
    Biome = c("Boreal forest", "Subtropical/tropical forest", "Subtropical/tropical savanna", 
              "Temperate broadleaf forest", "Temperate conifer forest"),
    Model = c("m1", "m1", "m3", "m1", "m3") # if considering m0: ("m0", "m0", "m3", "m0", "m3")
  )
  
  # Extract the selected model for each biome
  results_CH4 <- model_data[model_data$GHG == "CH4", ]
  results_N2O <- model_data[model_data$GHG == "N2O", ]
  
  # Filter to only include best models
  results_CH4 <- merge(results_CH4, ch4_best_models, by = c("Biome", "Model"))
  results_N2O <- merge(results_N2O, n2o_best_models, by = c("Biome", "Model"))
  
  # Create GHG_model_comparison.csv with only the best models
  model_comparison <- rbind(
    results_CH4[, .(GHG, Biome, Model, AIC, LogLik, Shapiro_p, BP_p)],
    results_N2O[, .(GHG, Biome, Model, AIC, LogLik, Shapiro_p, BP_p)]
  )
  
  # Write the best model selection to CSV
  fwrite(model_comparison, "GHG_model_comparison.csv")
  
  # Load agricultural data
  results_AG <- fread(ag_file)
  
  # Back transform ecosystem CH4 results with proper SE propagation
  results_CH4[, `:=`(
    B0_backtr = back_transform_ch4(Intercept),
    B0_SE_backtr = back_transform_ch4_se(Intercept, SE_Intercept),
    GHG = "Eco_CH4"
  )]
  
  # Back transform ecosystem N2O results with proper SE propagation
  results_N2O[, `:=`(
    B0_backtr = back_transform_n2o(Intercept),
    B0_SE_backtr = back_transform_n2o_se(Intercept, SE_Intercept),
    GHG = "Eco_N2O"
  )]
  
  # Prepare agricultural data - rename columns to match expected format
  setnames(results_AG,
           old = c("B0", "B0.SE"),
           new = c("B0_backtr", "B0_SE_backtr"))
  
  # Remove any X column if it exists
  if("X" %in% names(results_AG)) {
    results_AG[, X := NULL]
  }
  
  # Combine results
  results_CH4 <- results_CH4[, .(Biome, B0_backtr, B0_SE_backtr, GHG, Num_Cit, Model)]
  results_N2O <- results_N2O[, .(Biome, B0_backtr, B0_SE_backtr, GHG, Num_Cit, Model)]
  results_AG <- results_AG[, .(Biome, B0_backtr, B0_SE_backtr, GHG)]
  
  # Remove any existing "Global mean" rows from agricultural data
  results_AG <- results_AG[Biome != "Global mean"]
  
  ghg_data <- rbindlist(list(results_CH4, results_N2O, results_AG), fill = TRUE)
  
  # Add model info to results for reporting
  cat("Models selected for CH4:\n")
  for(i in 1:nrow(results_CH4)) {
    cat(sprintf("%s: %s\n", results_CH4$Biome[i], results_CH4$Model[i]))
  }
  
  cat("\nModels selected for N2O:\n")
  for(i in 1:nrow(results_N2O)) {
    cat(sprintf("%s: %s\n", results_N2O$Biome[i], results_N2O$Model[i]))
  }
  
  # Calculate net effects for ecosystems
  net_effect_eco <- ghg_data[GHG %in% c("Eco_CH4", "Eco_N2O"),
                             .(B0_backtr = sum(B0_backtr),
                               B0_SE_backtr = sqrt(sum(B0_SE_backtr^2, na.rm = TRUE)),
                               Num_Cit = first(Num_Cit),
                               GHG = "Eco_CH4-N2O"),
                             by = .(Biome)]
  
  ghg_data <- rbindlist(list(ghg_data, net_effect_eco), fill = TRUE)
  
  # Remove any existing global means
  ghg_data <- ghg_data[Biome != "Global mean"]
  
  # Calculate global means for each GHG category
  global_means <- ghg_data[!is.na(Biome),
                           .(B0_backtr = mean(B0_backtr, na.rm = TRUE),
                             B0_SE_backtr = mean(B0_SE_backtr, na.rm = TRUE),
                             Num_Cit = mean(Num_Cit, na.rm = TRUE)),
                           by = .(GHG)]
  global_means[, Biome := "Global mean"]
  
  # Combine with main dataset
  ghg_data <- rbindlist(list(ghg_data, global_means), fill = TRUE)
  
  # Convert to data.frame for ggplot
  ghg_data <- as.data.frame(ghg_data)
  
  # Set factor levels
  ghg_data$GHG <- factor(ghg_data$GHG,
                         levels = c("Eco_CH4", "Ag_CH4",
                                    "Eco_N2O", "Ag_N2O",
                                    "Eco_CH4-N2O", "Ag_CH4_N2O"))
  
  ghg_data$Biome <- factor(ghg_data$Biome,
                           levels = c("Temperate conifer forest",
                                      "Temperate broadleaf forest",
                                      "Subtropical/tropical savanna",
                                      "Subtropical/tropical forest",
                                      "Boreal forest",
                                      "Global mean"))
  
  return(ghg_data)
}

print_global_means <- function(ghg_data) {
  global_data <- subset(ghg_data, Biome == "Global mean")
  print("Global Means Summary:")
  print(global_data[order(global_data$GHG),
                    c("GHG", "B0_backtr", "B0_SE_backtr")])
} 

plot_ghg_comparison <- function(ghg_data, y_limit = NULL) {
  # Calculate reasonable y-axis limits if not provided
  if (is.null(y_limit)) {
    # Get the next largest value after excluding CH4 outliers
    non_ch4_max <- max(ghg_data$B0_backtr[!grepl("CH4", ghg_data$GHG)] +
                         ghg_data$B0_SE_backtr[!grepl("CH4", ghg_data$GHG)],
                       na.rm = TRUE)
    y_limit <- ceiling(non_ch4_max)
  }
  
  # Create main plot
  p <- ggplot(ghg_data, aes(x = GHG, y = B0_backtr, fill = GHG)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = B0_backtr - B0_SE_backtr,
                      ymax = B0_backtr + B0_SE_backtr),
                  width = 0.25, position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = c("purple4", "purple",  # CH4 colors
                                 "blue", "lightblue",        # N2O colors
                                 "gold2", "yellow")) +    # Net effect colors
    facet_wrap(~Biome, nrow = 3) +
    labs(y = expression("CO"[2]*"e soil flux (Mg" ~ ha^-1 ~ yr^-1*")")) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size=16),
          axis.text.y = element_text(size=16),
          axis.title = element_text(size = 18, face = "bold" ),
          strip.text.x = element_text(size = 18, colour = "black", face = "bold" ))
  
  # Add y-axis limits while preserving data integrity
  if (!is.null(y_limit)) {
    p <- p + 
      coord_cartesian(ylim = c(-0.8,3)) + #c(-y_limit/3.8, y_limit-0.8)) + # adjusted to minimize white space on plot
      # Add arrow annotation for out-of-bounds error bars
      geom_segment(data = subset(ghg_data, (B0_backtr + B0_SE_backtr) > y_limit),
                   aes(x = as.numeric(GHG),
                       xend = as.numeric(GHG),
                       y = y_limit * 0.95,
                       yend = y_limit),
                   arrow = arrow(length = unit(0.2, "cm")),
                   color = "black")
  }
  
  return(p)
}

###################################################### 
# Functions for Figure 3: Global predictions map
######################################################

# Updated function to create global predictions using best models
create_global_predictions <- function(df, model_data) {
  # Hardcoded best models based on AIC values (excluding m0)
  ch4_best_models <- data.table(
    Biome = c("Boreal forest", "Subtropical/tropical forest", "Subtropical/tropical savanna",
              "Temperate broadleaf forest", "Temperate conifer forest"),
    Model = c("m1", "m2", "m3", "m1", "m1")
  )
  
  n2o_best_models <- data.table(
    Biome = c("Boreal forest", "Subtropical/tropical forest", "Subtropical/tropical savanna",
              "Temperate broadleaf forest", "Temperate conifer forest"),
    Model = c("m1", "m1", "m3", "m1", "m3")
  )
  
  # Extract the coefficients for each best model
  results_CH4 <- model_data[model_data$GHG == "CH4", ]
  results_N2O <- model_data[model_data$GHG == "N2O", ]
  
  # Filter to only include best models
  results_CH4 <- merge(results_CH4, ch4_best_models, by = c("Biome", "Model"))
  results_N2O <- merge(results_N2O, n2o_best_models, by = c("Biome", "Model"))
  
  # Create predictions for CH4 based on the best models
  df$CH4.CO2e.pred <- NA_real_
  
  # Apply each biome's model to its respective data
  for (i in 1:nrow(ch4_best_models)) {
    biome <- ch4_best_models[i, Biome]
    model_type <- ch4_best_models[i, Model]
    model_coefs <- results_CH4[Biome == biome]
    
    # Filter data for this biome
    biome_rows <- which(df$Biome == biome)
    
    if (length(biome_rows) > 0) {
      # Always add intercept regardless of model complexity
      df$CH4.CO2e.pred[biome_rows] <- model_coefs$Intercept
      
      # Apply model based on its complexity
      if (model_type == "m1") {
        # Model 1: Intercept + 1 covariate
        if (biome == "Temperate conifer forest") {
          df$CH4.CO2e.pred[biome_rows] <- df$CH4.CO2e.pred[biome_rows] +
            model_coefs$B2_coeff * df$S_perc_clay[biome_rows]
        } else if (biome == "Boreal forest") {
          df$CH4.CO2e.pred[biome_rows] <- df$CH4.CO2e.pred[biome_rows] +
            model_coefs$B2_coeff * df$S_slope[biome_rows]
        } else if (biome == "Temperate broadleaf forest") {
          df$CH4.CO2e.pred[biome_rows] <- df$CH4.CO2e.pred[biome_rows] +
            model_coefs$B2_coeff * df$S_temp_seasn[biome_rows]
        }
      } else if (model_type == "m2") {
        # Model 2: Intercept + 2 covariates
        if (biome == "Subtropical/tropical forest") {
          df$CH4.CO2e.pred[biome_rows] <- df$CH4.CO2e.pred[biome_rows] +
            model_coefs$B2_coeff * df$S_slope[biome_rows] +
            model_coefs$B3_coeff * df$S_perc_sand[biome_rows]
        }
      } else if (model_type == "m3") {
        # Model 3: Intercept + 3 covariates
        if (biome == "Subtropical/tropical savanna") {
          df$CH4.CO2e.pred[biome_rows] <- df$CH4.CO2e.pred[biome_rows] +
            model_coefs$B2_coeff * df$S_temp_seasn[biome_rows] +
            model_coefs$B3_coeff * df$S_elevation[biome_rows] +
            model_coefs$B4_coeff * df$S_precip_seasn[biome_rows]
        }
      }
    }
  }
  
  # Create predictions for N2O based on the best models
  df$N2O.CO2e.pred <- NA_real_
  
  # Apply each biome's model to its respective data
  for (i in 1:nrow(n2o_best_models)) {
    biome <- n2o_best_models[i, Biome]
    model_type <- n2o_best_models[i, Model]
    model_coefs <- results_N2O[Biome == biome]
    
    # Filter data for this biome
    biome_rows <- which(df$Biome == biome)
    
    if (length(biome_rows) > 0) {
      # Always add intercept regardless of model complexity
      df$N2O.CO2e.pred[biome_rows] <- model_coefs$Intercept
      
      # Apply model based on its complexity
      if (model_type == "m1") {
        # Model 1: Intercept + 1 covariate
        if (biome == "Temperate broadleaf forest") {
          df$N2O.CO2e.pred[biome_rows] <- df$N2O.CO2e.pred[biome_rows] +
            model_coefs$B2_coeff * df$S_MAT[biome_rows]
        } else if (biome == "Boreal forest") {
          df$N2O.CO2e.pred[biome_rows] <- df$N2O.CO2e.pred[biome_rows] +
            model_coefs$B2_coeff * df$S_aspect[biome_rows]
        } else if (biome == "Subtropical/tropical forest") {
          df$N2O.CO2e.pred[biome_rows] <- df$N2O.CO2e.pred[biome_rows] +
            model_coefs$B2_coeff * df$S_slope[biome_rows]
        }
      } else if (model_type == "m3") {
        # Model 3: Intercept + 3 covariates
        if (biome == "Subtropical/tropical savanna") {
          df$N2O.CO2e.pred[biome_rows] <- df$N2O.CO2e.pred[biome_rows] +
            model_coefs$B2_coeff * df$S_temp_seasn[biome_rows] +
            model_coefs$B3_coeff * df$S_precip_seasn[biome_rows] +
            model_coefs$B4_coeff * df$S_elevation[biome_rows]
        } else if (biome == "Temperate conifer forest") {
          df$N2O.CO2e.pred[biome_rows] <- df$N2O.CO2e.pred[biome_rows] +
            model_coefs$B2_coeff * df$S_elevation[biome_rows] +
            model_coefs$B3_coeff * df$S_perc_clay[biome_rows] +
            model_coefs$B4_coeff * df$S_soil_porosity[biome_rows]
        }
      }
    }
  }
  
  # Back-transform results
  df$CH4.CO2e.BackTr <- exp(df$CH4.CO2e.pred) - 79.67499 # Adjust for log link and scaling to positive values
  df$N2O.CO2e.BackTr <- exp(df$N2O.CO2e.pred) - 2.358129 # Adjust for log link and scaling to positive values
  
  # Calculate net CH4-N2O effect
  df$net_effect <- df$CH4.CO2e.BackTr + df$N2O.CO2e.BackTr
  
  # Add diagnostic columns to check for missing predictions
  df$CH4_has_pred <- !is.na(df$CH4.CO2e.pred)
  df$N2O_has_pred <- !is.na(df$N2O.CO2e.pred)
  
  # Summarize prediction coverage by biome
  cat("CH4 prediction coverage by biome:\n")
  ch4_summary <- aggregate(CH4_has_pred ~ Biome, data=df, FUN=mean)
  print(ch4_summary)
  
  cat("\nN2O prediction coverage by biome:\n")
  n2o_summary <- aggregate(N2O_has_pred ~ Biome, data=df, FUN=mean)
  print(n2o_summary)
  
  return(df)
}

# You'll also need to update the plot_global_predictions function for better diagnostics:
plot_global_predictions <- function(pred_data) {
  # Report the number of NA vs non-NA values
  cat("CH4 predictions: ", sum(!is.na(pred_data$CH4.CO2e.BackTr)), " non-NA values, ", 
      sum(is.na(pred_data$CH4.CO2e.BackTr)), " NA values\n", sep="")
  
  cat("N2O predictions: ", sum(!is.na(pred_data$N2O.CO2e.BackTr)), " non-NA values, ", 
      sum(is.na(pred_data$N2O.CO2e.BackTr)), " NA values\n", sep="")
  
  cat("Net effect predictions: ", sum(!is.na(pred_data$net_effect)), " non-NA values, ", 
      sum(is.na(pred_data$net_effect)), " NA values\n", sep="")
  
  # Base color scheme
  cooling_colors <- c("darkgreen", "green4", "palegreen")
  warming_colors <- c("khaki1", "gold", "orange", "darkred")
  
  # Plot CH4
  p1 <- ggplot(pred_data) +
    # Add a light grey background layer for land masses
    geom_tile(aes(x = Lon, y = Lat), fill = "grey88") +
    # Add the predictions layer on top
    geom_tile(aes(x = Lon, y = Lat, fill = CH4.CO2e.BackTr)) +
    scale_fill_gradientn(
      colors = c(cooling_colors, "white", warming_colors),
      values = scales::rescale(c(-2, -1, -0.5, 0, 1, 2, 4)),
      limits = c(-2, 4),
      na.value = "grey88",  # This controls the color of NA values
      name = expression("CH"[4] ~ "CO"[2]*"e (Mg" ~ ha^-1 ~ yr^-1*")")
    ) +
    theme_minimal() +
    coord_equal() +
    ggtitle("A")
  
  # Plot N2O
  p2 <- ggplot(pred_data) +
    geom_tile(aes(x = Lon, y = Lat), fill = "grey88") +
    geom_tile(aes(x = Lon, y = Lat, fill = N2O.CO2e.BackTr)) +
    scale_fill_gradientn(
      colors = c("white", warming_colors),
      values = scales::rescale(c( 0, 1, 2, 3)),
      limits = c(0, 3),
      na.value = "grey88",
      name = expression("N"[2]*"O CO"[2]*"e (Mg" ~ ha^-1 ~ yr^-1*")")
    ) +
    theme_minimal() +
    coord_equal() +
    ggtitle("B")
  
  # Plot net effect
  p3 <- ggplot(pred_data) +
    geom_tile(aes(x = Lon, y = Lat), fill = "grey88") +
    geom_tile(aes(x = Lon, y = Lat, fill = net_effect)) +
    scale_fill_gradientn(
      colors = c(cooling_colors, "white", warming_colors),
      values = scales::rescale(c(-2, -1, 0, 1, 3, 6)),  #  -0.5,
      limits = c(-2, 6),
      na.value = "grey88",
      name = expression("Net CO"[2]*"e (Mg" ~ ha^-1 ~ yr^-1*")")
    ) +
    theme_minimal() +
    coord_equal() +
    ggtitle("C")
  
  return(list(ch4 = p1, n2o = p2, net = p3))
}

create_summary_tables <- function(ghg_data) {
  # Ensure we're working with a data.table
  if (!is.data.table(ghg_data)) {
    ghg_data <- as.data.table(ghg_data)
  }
  
  # Table 1A: CH4 and N2O model results
  table_1a <- ghg_data[Biome != "Global mean" &
                         GHG %in% c("Eco_CH4", "Eco_N2O")]
  
  # Create separate summaries for CH4 and N2O
  ch4_summary <- table_1a[GHG == "Eco_CH4",
                          .(Biome,
                            CH4_mean = B0_backtr,
                            CH4_SE = B0_SE_backtr,
                            CH4_N = Num_Cit)]
  
  n2o_summary <- table_1a[GHG == "Eco_N2O",
                          .(Biome,
                            N2O_mean = B0_backtr,
                            N2O_SE = B0_SE_backtr,
                            N2O_N = Num_Cit)]
  
  # Merge CH4 and N2O summaries
  table_1a_final <- merge(ch4_summary, n2o_summary, by = "Biome")
  setnames(table_1a_final,
           old = c("CH4_mean", "CH4_SE", "CH4_N",
                   "N2O_mean", "N2O_SE", "N2O_N"),
           new = c("Mean CH4", "CH4 SE", "N obs CH4",
                   "Mean N2O", "N2O SE", "N obs N2O"))
  
  # Table 1B: Net effects
  table_1b <- ghg_data[Biome != "Global mean" &
                         GHG == "Eco_CH4-N2O",
                       .(Biome,
                         "Net CH4-N2O effect" = B0_backtr,
                         "Standard Error" = B0_SE_backtr)]
  
  return(list(table_1a = table_1a_final, table_1b = table_1b))
}

###################################################### 
# Global Effect Calculation
######################################################
calculate_global_effect <- function(pred_data, root) {
  # Create a raster from the points with specified resolution
  r <- rast(
    xmin = min(pred_data$Lon, na.rm = TRUE),
    xmax = max(pred_data$Lon, na.rm = TRUE),
    ymin = min(pred_data$Lat, na.rm = TRUE),
    ymax = max(pred_data$Lat, na.rm = TRUE),
    resolution = 0.5,  # 0.5 degree resolution
    crs = "EPSG:4326"
  )
  
  # Rasterize the points
  net_effect_raster <- rasterize(
    x = cbind(pred_data$Lon, pred_data$Lat),
    y = r,
    values = pred_data$net_effect,
    fun = mean  # Use mean if multiple points fall in same cell
  )
  
  # Calculate pixel area in hectares
  area_raster <- terra::cellSize(net_effect_raster, unit = "ha")
  
  # Mask the area raster to only non-NA values in net_effect
  valid_areas <- area_raster * !is.na(net_effect_raster)
  
  # Sum up the total area
  total_area_ha <- sum(values(valid_areas), na.rm = TRUE)
  
  # Calculate Total global net CH4-N2O effect
  MgCO2e_yr <- valid_areas * net_effect_raster
  total_MgCO2e_yr <- sum(values(MgCO2e_yr), na.rm = TRUE)
  
  # Report results
  cat("Total area of non-NA pixels:", format(total_area_ha, big.mark=","), "ha\n")
  cat("Total global net CH4-N2O effect:", format(total_MgCO2e_yr, big.mark=","), "Mg CO2e / year\n")
  
  # Optional: Remove outliers and recalculate
  # Create reclassification matrix for terra
  rcl <- matrix(c(400, Inf, NA,  # Values >= 400 become NA
                  -Inf, -400, NA),  # Values <= -400 become NA
                ncol=3, byrow=TRUE)
  
  # Reclassify
  net_effect_r_NoOut <- classify(net_effect_raster, rcl)
  MgCO2e_yr_NoOut <- valid_areas * net_effect_r_NoOut
  total_MgCO2e_yr_NoOut <- sum(values(MgCO2e_yr_NoOut), na.rm = TRUE)
  
  cat("Total global net CH4-N2O effect (excluding outliers):", 
      format(total_MgCO2e_yr_NoOut, big.mark=","), "Mg CO2e / year\n")
  
  # Save results to file
  results <- data.frame(
    total_area_ha = total_area_ha,
    total_MgCO2e_yr = total_MgCO2e_yr,
    total_MgCO2e_yr_NoOut = total_MgCO2e_yr_NoOut
  )
  
  write.csv(results, file.path(root, "global_effect_results.csv"), row.names = FALSE)
  
  return(results)
}

###################################################### 
# Main Execution Function
######################################################

main_analysis <- function(root, figsDir, model_data_file = "Table_S3_GHG_reg_models.txt") {
  # Create directories if they don't exist
  if (!dir.exists(figsDir)) {
    dir.create(figsDir, recursive = TRUE)
    cat("Created output directory:", figsDir, "\n")
  }
  
  # Read in the model data from the provided file
  cat("Reading model data from:", model_data_file, "\n")
  model_data <- fread(model_data_file)
  
  # Prepare the GHG data for Figure 2
  cat("Preparing GHG data for Figure 2...\n")
  ag_file <- file.path(root, "__summ_Ag_CH4_N2O_xy_noRice.csv")
  if (!file.exists(ag_file)) {
    stop("Agricultural data file not found: ", ag_file)
  }
  
  ghg_data <- prepare_ghg_data(
    model_data,
    ag_file = ag_file
  )
  
  # Print global means
  cat("\nGlobal means summary:\n")
  print_global_means(ghg_data)
  
  # Create Figure 2
  cat("Creating Figure 2...\n")
  fig2 <- plot_ghg_comparison(ghg_data) #, y_limit = 4)
  
  # Read in the global predictions data
  preds_file <- file.path(root, "global_preds_vMar2025.csv")
  if (!file.exists(preds_file)) {
    stop("Global predictions file not found: ", preds_file)
  }
  
  cat("Reading global predictions data...\n")
  pred_data <- fread(preds_file)
  
  # Update predictions using the best models
  cat("Creating updated predictions with best AIC models...\n")
  updated_pred_data <- create_global_predictions(pred_data, model_data)
  
  # Write the updated predictions to file
  updated_preds_file <- file.path(root, "global_preds_updated.csv")
  cat("Writing updated predictions to:", updated_preds_file, "\n")
  fwrite(updated_pred_data, updated_preds_file)
  
  # Create Figure 3
  cat("Creating Figure 3...\n")
  plots <- plot_global_predictions(updated_pred_data)
  
  # Arrange plots vertically
  fig3 <- gridExtra::grid.arrange(plots$ch4, plots$n2o, plots$net, 
                                  ncol = 1, heights = c(1, 1, 1))
  
  # Save the figures
  fig2_file <- file.path(figsDir, "figure2.pdf")
  fig3_file <- file.path(figsDir, "figure3.pdf")
  
  cat("Saving Figure 2 to:", fig2_file, "\n")
  ggsave(fig2_file, fig2, width = 12, height = 10)
  
  cat("Saving Figure 3 to:", fig3_file, "\n")
  ggsave(fig3_file, fig3, width = 10, height = 15)
  
  # Create summary tables
  cat("Creating summary tables...\n")
  tables <- create_summary_tables(ghg_data)
  
  table1a_file <- file.path(figsDir, "table_1a.csv")
  table1b_file <- file.path(figsDir, "table_1b.csv")
  
  cat("Saving summary tables to:", table1a_file, "and", table1b_file, "\n")
  write.csv(tables$table_1a, table1a_file, row.names = FALSE)
  write.csv(tables$table_1b, table1b_file, row.names = FALSE)
  
  # Calculate global CH4-N2O effect per year
  cat("Calculating global CH4-N2O effects...\n")
  calculate_global_effect(updated_pred_data, root)
  
  cat("Analysis complete!\n")
  
  # Return the results for inspection
  return(list(
    ghg_data = ghg_data,
    fig2 = fig2,
    fig3 = fig3,
    pred_data = updated_pred_data,
    tables = tables
  ))
}

# Ensure I have these required files in root directory:
#    - Table_S3_GHG_reg_models.txt (model results data)
#    - __summ_Ag_CH4_N2O_xy_noRice.csv (agricultural data)
#    - global_preds_vMar2025.csv # used to be: global_preds_vNov24.csv (global predictions data)

setwd(root)
# Run the main analysis function
results <- main_analysis(root, figsDir, "Table_S3_GHG_reg_models_20yr_noPeat_vMarch2025.csv")

# Examine the Results

# The analysis will produce the following files:
#    - GHG_model_comparison.csv (best model for each GHG-biome pair)
#    - global_preds_updated.csv (updated predictions using best models)
#    - figure2.pdf (GHG boxplots with agriculture data)
#    - figure3.pdf (global predictions maps)
#    - table_1a.csv and table_1b.csv (summary tables)
#    - global_effect_results.csv (global CH4-N2O effects)


# ----------------------------------------------
# Heterogeneity and Sampling bias correction
# ----------------------------------------------

source("/Users/savannah/Documents/GitHub/ch4.Meta_analysis/6_ext_heterogeneity-functions.R") #("6_ext_heterogeneity-functions.R")
#source("/Users/savannah/Documents/GitHub/ch4.Meta_analysis/6_ext_sampling-bias-functions.R")
source("/Users/savannah/Documents/GitHub/ch4.Meta_analysis/6_ext_visualization-results.R")
source("/Users/savannah/Documents/GitHub/ch4.Meta_analysis/6_ext_implementation.R") # Add the enhanced implementation to my workflow

results_combined <- main_analysis_enhanced(root, figsDir, "Table_S3_GHG_reg_models_20yr_noPeat_vMarch2025.csv")
# Updated workflow for visualizing the enhanced analysis results

# 1. First, visualize heterogeneity test results
het_test_plot <- plot_heterogeneity_test_results(
  results_combined$enhanced$model_results,
  file.path(enhanced_results_dir, "heterogeneity_tests.pdf")
)

# 2. Visualize the impact of robust modeling on coefficients
het_impact_plot <- plot_heterogeneity_impact(
  results_combined$enhanced$model_results,
  file.path(enhanced_results_dir, "heterogeneity_impact.pdf")
)

# 3. Create a summary table of heterogeneity results
het_summary <- create_heterogeneity_table(
  results_combined$enhanced$model_results,
  file.path(enhanced_results_dir, "heterogeneity_summary.csv")
)

# 4. Create a comprehensive visualization comparing original vs robust results
comp_viz <- create_comprehensive_visualization(
  results_combined$enhanced,
  file.path(enhanced_results_dir, "comprehensive_comparison.pdf")
)

# 5. Create a summary table of the best models
# For each GHG-biome combination
best_models <- data.frame(
  GHG = c("CH4", "CH4", "CH4", "CH4", "CH4",
          "N2O", "N2O", "N2O", "N2O", "N2O"),
  Biome = c("Boreal forest", "Subtropical/tropical forest", "Subtropical/tropical savanna",
            "Temperate broadleaf forest", "Temperate conifer forest",
            "Boreal forest", "Subtropical/tropical forest", "Subtropical/tropical savanna",
            "Temperate broadleaf forest", "Temperate conifer forest"),
  Model = c("m1", "m2", "m3", "m1", "m1",
            "m1", "m1", "m3", "m1", "m3")
)

# Filter the model results to get only the best models
model_data <- results_combined$enhanced$model_results

# Create a comparison table where original and enhanced are the same 
# (since we're just comparing robust vs non-robust approaches)
model_comparison <- create_model_comparison_table(
  original_results = model_data, 
  enhanced_results = model_data,
  output_file = file.path(enhanced_results_dir, "model_comparison.csv")
)

# Print a summary of what was generated
cat("Enhanced analysis visualization complete.\n")
cat("Results saved to:", enhanced_results_dir, "\n")



###################################################### 
# ROBUST - Enhanced Analysis Code
######################################################

# This contains the modified functions needed to replace the original models with robust models where appropriate (based on heteroscedasticity detection). 

# Modified prepare_ghg_data function to use robust models where appropriate
prepare_ghg_data_robust <- function(model_data, enhanced_model_data, ag_file) {
  # Hardcoded best models based on AIC values
  ch4_best_models <- data.table(
    Biome = c("Boreal forest", "Subtropical/tropical forest", "Subtropical/tropical savanna",
              "Temperate broadleaf forest", "Temperate conifer forest"),
    Model = c("m1", "m2", "m3", "m1", "m1") # if considering m0: ("m0", "m2", "m3", "m0", "m1")
  )
  
  n2o_best_models <- data.table(
    Biome = c("Boreal forest", "Subtropical/tropical forest", "Subtropical/tropical savanna",
              "Temperate broadleaf forest", "Temperate conifer forest"),
    Model = c("m1", "m1", "m3", "m1", "m3") # if considering m0: ("m0", "m0", "m3", "m0", "m3")
  )
  
  # Extract the selected model for each biome
  results_CH4 <- model_data[model_data$GHG == "CH4", ]
  results_N2O <- model_data[model_data$GHG == "N2O", ]
  
  # Filter to only include best models
  results_CH4 <- merge(results_CH4, ch4_best_models, by = c("Biome", "Model"))
  results_N2O <- merge(results_N2O, n2o_best_models, by = c("Biome", "Model"))
  
  # Now get the enhanced model data for the same models
  enhanced_CH4 <- enhanced_model_data[enhanced_model_data$GHG == "CH4", ]
  enhanced_N2O <- enhanced_model_data[enhanced_model_data$GHG == "N2O", ]
  
  # Filter to only include best models
  enhanced_CH4 <- merge(enhanced_CH4, ch4_best_models, by = c("Biome", "Model"))
  enhanced_N2O <- merge(enhanced_N2O, n2o_best_models, by = c("Biome", "Model"))
  
  # Create GHG_model_comparison.csv with only the best models
  model_comparison <- rbind(
    results_CH4[, .(GHG, Biome, Model, AIC, LogLik, Shapiro_p, BP_p)],
    results_N2O[, .(GHG, Biome, Model, AIC, LogLik, Shapiro_p, BP_p)]
  )
  
  # Write the best model selection to CSV
  fwrite(model_comparison, file.path(enhanced_results_dir, "GHG_model_comparison.csv"))
  
  # Load agricultural data
  results_AG <- fread(ag_file)
  
  # For each model, check if we should use the robust version instead of the original
  # We'll modify the CH4 and N2O results accordingly
  for (i in 1:nrow(results_CH4)) {
    biome <- results_CH4$Biome[i]
    model <- results_CH4$Model[i]
    
    # Find the corresponding enhanced model
    enhanced_idx <- which(enhanced_CH4$Biome == biome & enhanced_CH4$Model == model)
    
    if (length(enhanced_idx) > 0 && !is.na(enhanced_CH4$Used_Robust[enhanced_idx]) && enhanced_CH4$Used_Robust[enhanced_idx]) {
      # Replace with robust values
      cat("Using robust model for CH4", biome, model, "\n")
      
      # Check if the robust model values are available
      if (!is.na(enhanced_CH4$Robust_Intercept[enhanced_idx])) {
        results_CH4$Intercept[i] <- enhanced_CH4$Robust_Intercept[enhanced_idx]
        results_CH4$SE_Intercept[i] <- enhanced_CH4$Robust_SE_Intercept[enhanced_idx]
      }
      
      if (!is.na(enhanced_CH4$Robust_Age_coeff[enhanced_idx])) {
        results_CH4$Age_coeff[i] <- enhanced_CH4$Robust_Age_coeff[enhanced_idx]
        results_CH4$SE_Age_coeff[i] <- enhanced_CH4$Robust_SE_Age_coeff[enhanced_idx]
      }
    }
  }
  
  for (i in 1:nrow(results_N2O)) {
    biome <- results_N2O$Biome[i]
    model <- results_N2O$Model[i]
    
    # Find the corresponding enhanced model
    enhanced_idx <- which(enhanced_N2O$Biome == biome & enhanced_N2O$Model == model)
    
    if (length(enhanced_idx) > 0 && !is.na(enhanced_N2O$Used_Robust[enhanced_idx]) && enhanced_N2O$Used_Robust[enhanced_idx]) {
      # Replace with robust values
      cat("Using robust model for N2O", biome, model, "\n")
      
      # Check if the robust model values are available
      if (!is.na(enhanced_N2O$Robust_Intercept[enhanced_idx])) {
        results_N2O$Intercept[i] <- enhanced_N2O$Robust_Intercept[enhanced_idx]
        results_N2O$SE_Intercept[i] <- enhanced_N2O$Robust_SE_Intercept[enhanced_idx]
      }
      
      if (!is.na(enhanced_N2O$Robust_Age_coeff[enhanced_idx])) {
        results_N2O$Age_coeff[i] <- enhanced_N2O$Robust_Age_coeff[enhanced_idx]
        results_N2O$SE_Age_coeff[i] <- enhanced_N2O$Robust_SE_Age_coeff[enhanced_idx]
      }
    }
  }
  
  # Back transform ecosystem CH4 results with proper SE propagation
  results_CH4[, `:=`(
    B0_backtr = back_transform_ch4(Intercept),
    B0_SE_backtr = back_transform_ch4_se(Intercept, SE_Intercept),
    GHG = "Eco_CH4"
  )]
  
  # Back transform ecosystem N2O results with proper SE propagation
  results_N2O[, `:=`(
    B0_backtr = back_transform_n2o(Intercept),
    B0_SE_backtr = back_transform_n2o_se(Intercept, SE_Intercept),
    GHG = "Eco_N2O"
  )]
  
  # Prepare agricultural data - rename columns to match expected format
  setnames(results_AG,
           old = c("B0", "B0.SE"),
           new = c("B0_backtr", "B0_SE_backtr"))
  
  # Remove any X column if it exists
  if("X" %in% names(results_AG)) {
    results_AG[, X := NULL]
  }
  
  # Combine results
  results_CH4 <- results_CH4[, .(Biome, B0_backtr, B0_SE_backtr, GHG, Num_Cit, Model)]
  results_N2O <- results_N2O[, .(Biome, B0_backtr, B0_SE_backtr, GHG, Num_Cit, Model)]
  results_AG <- results_AG[, .(Biome, B0_backtr, B0_SE_backtr, GHG)]
  
  # Remove any existing "Global mean" rows from agricultural data
  results_AG <- results_AG[Biome != "Global mean"]
  
  ghg_data <- rbindlist(list(results_CH4, results_N2O, results_AG), fill = TRUE)
  
  # Create summary of which robust models were used
  cat("\nModels replaced with robust versions:\n")
  cat("CH4 models:", sum(enhanced_CH4$Used_Robust, na.rm = TRUE), "of", nrow(enhanced_CH4), "\n")
  cat("N2O models:", sum(enhanced_N2O$Used_Robust, na.rm = TRUE), "of", nrow(enhanced_N2O), "\n")
  
  # Add model info to results for reporting
  cat("\nModels selected for CH4:\n")
  for(i in 1:nrow(results_CH4)) {
    robust_flag <- ""
    enhanced_idx <- which(enhanced_CH4$Biome == results_CH4$Biome[i] & enhanced_CH4$Model == results_CH4$Model[i])
    if (length(enhanced_idx) > 0 && !is.na(enhanced_CH4$Used_Robust[enhanced_idx]) && enhanced_CH4$Used_Robust[enhanced_idx]) {
      robust_flag <- " (robust)"
    }
    cat(sprintf("%s: %s%s\n", results_CH4$Biome[i], results_CH4$Model[i], robust_flag))
  }
  
  cat("\nModels selected for N2O:\n")
  for(i in 1:nrow(results_N2O)) {
    robust_flag <- ""
    enhanced_idx <- which(enhanced_N2O$Biome == results_N2O$Biome[i] & enhanced_N2O$Model == results_N2O$Model[i])
    if (length(enhanced_idx) > 0 && !is.na(enhanced_N2O$Used_Robust[enhanced_idx]) && enhanced_N2O$Used_Robust[enhanced_idx]) {
      robust_flag <- " (robust)"
    }
    cat(sprintf("%s: %s%s\n", results_N2O$Biome[i], results_N2O$Model[i], robust_flag))
  }
  
  # Calculate net effects for ecosystems
  net_effect_eco <- ghg_data[GHG %in% c("Eco_CH4", "Eco_N2O"),
                             .(B0_backtr = sum(B0_backtr),
                               B0_SE_backtr = sqrt(sum(B0_SE_backtr^2, na.rm = TRUE)),
                               Num_Cit = first(Num_Cit),
                               GHG = "Eco_CH4-N2O"),
                             by = .(Biome)]
  
  ghg_data <- rbindlist(list(ghg_data, net_effect_eco), fill = TRUE)
  
  # Remove any existing global means
  ghg_data <- ghg_data[Biome != "Global mean"]
  
  # Calculate global means for each GHG category
  global_means <- ghg_data[!is.na(Biome),
                           .(B0_backtr = mean(B0_backtr, na.rm = TRUE),
                             B0_SE_backtr = mean(B0_SE_backtr, na.rm = TRUE),
                             Num_Cit = mean(Num_Cit, na.rm = TRUE)),
                           by = .(GHG)]
  global_means[, Biome := "Global mean"]
  
  # Combine with main dataset
  ghg_data <- rbindlist(list(ghg_data, global_means), fill = TRUE)
  
  # Convert to data.frame for ggplot
  ghg_data <- as.data.frame(ghg_data)
  
  # Set factor levels
  ghg_data$GHG <- factor(ghg_data$GHG,
                         levels = c("Eco_CH4", "Ag_CH4",
                                    "Eco_N2O", "Ag_N2O",
                                    "Eco_CH4-N2O", "Ag_CH4_N2O"))
  
  ghg_data$Biome <- factor(ghg_data$Biome,
                           levels = c("Temperate conifer forest",
                                      "Temperate broadleaf forest",
                                      "Subtropical/tropical savanna",
                                      "Subtropical/tropical forest",
                                      "Boreal forest",
                                      "Global mean"))
  
  return(ghg_data)
}

# Updated function to create global predictions using robust versions of best models where appropriate
create_global_predictions_robust <- function(df, model_data, enhanced_model_data) {
  # Hardcoded best models based on AIC values (excluding m0)
  ch4_best_models <- data.table(
    Biome = c("Boreal forest", "Subtropical/tropical forest", "Subtropical/tropical savanna",
              "Temperate broadleaf forest", "Temperate conifer forest"),
    Model = c("m1", "m2", "m3", "m1", "m1")
  )
  
  n2o_best_models <- data.table(
    Biome = c("Boreal forest", "Subtropical/tropical forest", "Subtropical/tropical savanna",
              "Temperate broadleaf forest", "Temperate conifer forest"),
    Model = c("m1", "m1", "m3", "m1", "m3")
  )
  
  # Extract the coefficients for each best model
  results_CH4 <- model_data[model_data$GHG == "CH4", ]
  results_N2O <- model_data[model_data$GHG == "N2O", ]
  
  # Filter to only include best models
  results_CH4 <- merge(results_CH4, ch4_best_models, by = c("Biome", "Model"))
  results_N2O <- merge(results_N2O, n2o_best_models, by = c("Biome", "Model"))
  
  # Get enhanced model results
  enhanced_CH4 <- enhanced_model_data[enhanced_model_data$GHG == "CH4", ]
  enhanced_N2O <- enhanced_model_data[enhanced_model_data$GHG == "N2O", ]
  
  # Filter to only include best models
  enhanced_CH4 <- merge(enhanced_CH4, ch4_best_models, by = c("Biome", "Model"))
  enhanced_N2O <- merge(enhanced_N2O, n2o_best_models, by = c("Biome", "Model"))
  
  # For each model, check if we should use the robust version instead of the original
  # We'll modify the CH4 and N2O results accordingly
  cat("Checking which models should use robust versions for global predictions:\n")
  for (i in 1:nrow(results_CH4)) {
    biome <- results_CH4$Biome[i]
    model <- results_CH4$Model[i]
    
    # Find the corresponding enhanced model
    enhanced_idx <- which(enhanced_CH4$Biome == biome & enhanced_CH4$Model == model)
    
    if (length(enhanced_idx) > 0 && !is.na(enhanced_CH4$Used_Robust[enhanced_idx]) && enhanced_CH4$Used_Robust[enhanced_idx]) {
      # Replace with robust values
      cat("Using robust model for CH4", biome, model, "\n")
      
      # Update intercept with robust version
      if (!is.na(enhanced_CH4$Robust_Intercept[enhanced_idx])) {
        results_CH4$Intercept[i] <- enhanced_CH4$Robust_Intercept[enhanced_idx]
      }
      
      # The dataset doesn't have separate robust versions of B2_coeff, B3_coeff, B4_coeff
      # So we'll only update the Age coefficient
      # For improved accuracy, we may need to rerun the models to get proper robust coefficients for the covariates
    }
  }
  
  for (i in 1:nrow(results_N2O)) {
    biome <- results_N2O$Biome[i]
    model <- results_N2O$Model[i]
    
    # Find the corresponding enhanced model
    enhanced_idx <- which(enhanced_N2O$Biome == biome & enhanced_N2O$Model == model)
    
    if (length(enhanced_idx) > 0 && !is.na(enhanced_N2O$Used_Robust[enhanced_idx]) && enhanced_N2O$Used_Robust[enhanced_idx]) {
      # Replace with robust values
      cat("Using robust model for N2O", biome, model, "\n")
      
      # Update intercept with robust version
      if (!is.na(enhanced_N2O$Robust_Intercept[enhanced_idx])) {
        results_N2O$Intercept[i] <- enhanced_N2O$Robust_Intercept[enhanced_idx]
      }
      
      # The dataset doesn't have separate robust versions of B2_coeff, B3_coeff, B4_coeff
      # So we'll only update the Age coefficient
      # For improved accuracy, we may need to rerun the models to get proper robust coefficients for the covariates
    }
  }
  
  # Create predictions for CH4 based on the best models (potentially with robust coefficients)
  df$CH4.CO2e.pred <- NA_real_
  
  # Apply each biome's model to its respective data
  for (i in 1:nrow(ch4_best_models)) {
    biome <- ch4_best_models[i, Biome]
    model_type <- ch4_best_models[i, Model]
    model_coefs <- results_CH4[results_CH4$Biome == biome, ]
    
    # Filter data for this biome
    biome_rows <- which(df$Biome == biome)
    
    if (length(biome_rows) > 0 && nrow(model_coefs) > 0) {
      # Apply model based on its complexity
      if (model_type == "m1") {
        # Model 1: Intercept + 1 covariate
        if (biome == "Temperate conifer forest") {
          df$CH4.CO2e.pred[biome_rows] <- model_coefs$Intercept +
            model_coefs$B2_coeff * df$S_perc_clay[biome_rows]
        } else if (biome == "Boreal forest") {
          df$CH4.CO2e.pred[biome_rows] <- model_coefs$Intercept +
            model_coefs$B2_coeff * df$S_slope[biome_rows]
        } else if (biome == "Temperate broadleaf forest") {
          df$CH4.CO2e.pred[biome_rows] <- model_coefs$Intercept +
            model_coefs$B2_coeff * df$S_temp_seasn[biome_rows]
        }
      } else if (model_type == "m2") {
        # Model 2: Intercept + 2 covariates
        if (biome == "Subtropical/tropical forest") {
          df$CH4.CO2e.pred[biome_rows] <- model_coefs$Intercept +
            model_coefs$B2_coeff * df$S_slope[biome_rows] +
            model_coefs$B3_coeff * df$S_perc_sand[biome_rows]
        }
      } else if (model_type == "m3") {
        # Model 3: Intercept + 3 covariates
        if (biome == "Subtropical/tropical savanna") {
          df$CH4.CO2e.pred[biome_rows] <- model_coefs$Intercept +
            model_coefs$B2_coeff * df$S_temp_seasn[biome_rows] +
            model_coefs$B3_coeff * df$S_elevation[biome_rows] +
            model_coefs$B4_coeff * df$S_precip_seasn[biome_rows]
        }
      }
    }
  }
  
  # Create predictions for N2O based on the best models (potentially with robust coefficients)
  df$N2O.CO2e.pred <- NA_real_
  
  # Apply each biome's model to its respective data
  for (i in 1:nrow(n2o_best_models)) {
    biome <- n2o_best_models[i, Biome]
    model_type <- n2o_best_models[i, Model]
    model_coefs <- results_N2O[results_N2O$Biome == biome, ]
    
    # Filter data for this biome
    biome_rows <- which(df$Biome == biome)
    
    if (length(biome_rows) > 0 && nrow(model_coefs) > 0) {
      # Apply model based on its complexity
      if (model_type == "m1") {
        # Model 1: Intercept + 1 covariate
        if (biome == "Temperate broadleaf forest") {
          df$N2O.CO2e.pred[biome_rows] <- model_coefs$Intercept +
            model_coefs$B2_coeff * df$S_MAT[biome_rows]
        } else if (biome == "Boreal forest") {
          df$N2O.CO2e.pred[biome_rows] <- model_coefs$Intercept +
            model_coefs$B2_coeff * df$S_aspect[biome_rows]
        } else if (biome == "Subtropical/tropical forest") {
          df$N2O.CO2e.pred[biome_rows] <- model_coefs$Intercept +
            model_coefs$B2_coeff * df$S_slope[biome_rows]
        }
      } else if (model_type == "m3") {
        # Model 3: Intercept + 3 covariates
        if (biome == "Subtropical/tropical savanna") {
          df$N2O.CO2e.pred[biome_rows] <- model_coefs$Intercept +
            model_coefs$B2_coeff * df$S_temp_seasn[biome_rows] +
            model_coefs$B3_coeff * df$S_precip_seasn[biome_rows] +
            model_coefs$B4_coeff * df$S_elevation[biome_rows]
        } else if (biome == "Temperate conifer forest") {
          df$N2O.CO2e.pred[biome_rows] <- model_coefs$Intercept +
            model_coefs$B2_coeff * df$S_elevation[biome_rows] +
            model_coefs$B3_coeff * df$S_perc_clay[biome_rows] +
            model_coefs$B4_coeff * df$S_soil_porosity[biome_rows]
        }
      }
    }
  }
  
  # Back-transform results
  df$CH4.CO2e.BackTr <- exp(df$CH4.CO2e.pred) - 79.67499 # Adjust for log link and scaling to positive values
  df$N2O.CO2e.BackTr <- exp(df$N2O.CO2e.pred) - 2.358129 # Adjust for log link and scaling to positive values
  
  # Calculate net CH4-N2O effect
  df$net_effect <- df$CH4.CO2e.BackTr + df$N2O.CO2e.BackTr
  
  return(df)
}

# Modified calculate_global_effect function to save results to enhanced_results_dir
calculate_global_effect_robust <- function(pred_data, root, enhanced_results_dir) {
  # Create a raster from the points with specified resolution
  r <- rast(
    xmin = min(pred_data$Lon, na.rm = TRUE),
    xmax = max(pred_data$Lon, na.rm = TRUE),
    ymin = min(pred_data$Lat, na.rm = TRUE),
    ymax = max(pred_data$Lat, na.rm = TRUE),
    resolution = 0.5,  # 0.5 degree resolution
    crs = "EPSG:4326"
  )
  
  # Rasterize the points
  net_effect_raster <- rasterize(
    x = cbind(pred_data$Lon, pred_data$Lat),
    y = r,
    values = pred_data$net_effect,
    fun = mean  # Use mean if multiple points fall in same cell
  )
  
  # Calculate pixel area in hectares
  area_raster <- terra::cellSize(net_effect_raster, unit = "ha")
  
  # Mask the area raster to only non-NA values in net_effect
  valid_areas <- area_raster * !is.na(net_effect_raster)
  
  # Sum up the total area
  total_area_ha <- sum(values(valid_areas), na.rm = TRUE)
  
  # Calculate Total global net CH4-N2O effect
  MgCO2e_yr <- valid_areas * net_effect_raster
  total_MgCO2e_yr <- sum(values(MgCO2e_yr), na.rm = TRUE)
  
  # Report results
  cat("Total area of non-NA pixels:", format(total_area_ha, big.mark=","), "ha\n")
  cat("Total global net CH4-N2O effect:", format(total_MgCO2e_yr, big.mark=","), "Mg CO2e / year\n")
  
  # Optional: Remove outliers and recalculate
  # Create reclassification matrix for terra
  rcl <- matrix(c(400, Inf, NA,  # Values >= 400 become NA
                  -Inf, -400, NA),  # Values <= -400 become NA
                ncol=3, byrow=TRUE)
  
  # Reclassify
  net_effect_r_NoOut <- classify(net_effect_raster, rcl)
  MgCO2e_yr_NoOut <- valid_areas * net_effect_r_NoOut
  total_MgCO2e_yr_NoOut <- sum(values(MgCO2e_yr_NoOut), na.rm = TRUE)
  
  cat("Total global net CH4-N2O effect (excluding outliers):", 
      format(total_MgCO2e_yr_NoOut, big.mark=","), "Mg CO2e / year\n")
  
  # Save results to file
  results <- data.frame(
    total_area_ha = total_area_ha,
    total_MgCO2e_yr = total_MgCO2e_yr,
    total_MgCO2e_yr_NoOut = total_MgCO2e_yr_NoOut
  )
  
  # Save to enhanced results directory
  write.csv(results, file.path(enhanced_results_dir, "global_effect_results_robust.csv"), row.names = FALSE)
  
  return(results)
}

# Modified create_summary_tables function to use robust model results
create_summary_tables_robust <- function(ghg_data) {
  # Ensure we're working with a data.table
  if (!is.data.table(ghg_data)) {
    ghg_data <- as.data.table(ghg_data)
  }
  
  # Table 1A: CH4 and N2O model results
  table_1a <- ghg_data[Biome != "Global mean" &
                         GHG %in% c("Eco_CH4", "Eco_N2O")]
  
  # Create separate summaries for CH4 and N2O
  ch4_summary <- table_1a[GHG == "Eco_CH4",
                          .(Biome,
                            CH4_mean = B0_backtr,
                            CH4_SE = B0_SE_backtr,
                            CH4_N = Num_Cit)]
  
  n2o_summary <- table_1a[GHG == "Eco_N2O",
                          .(Biome,
                            N2O_mean = B0_backtr,
                            N2O_SE = B0_SE_backtr,
                            N2O_N = Num_Cit)]
  
  # Merge CH4 and N2O summaries
  table_1a_final <- merge(ch4_summary, n2o_summary, by = "Biome")
  setnames(table_1a_final,
           old = c("CH4_mean", "CH4_SE", "CH4_N",
                   "N2O_mean", "N2O_SE", "N2O_N"),
           new = c("Mean CH4", "CH4 SE", "N obs CH4",
                   "Mean N2O", "N2O SE", "N obs N2O"))
  
  # Table 1B: Net effects
  table_1b <- ghg_data[Biome != "Global mean" &
                         GHG == "Eco_CH4-N2O",
                       .(Biome,
                         "Net CH4-N2O effect" = B0_backtr,
                         "Standard Error" = B0_SE_backtr)]
  
  return(list(table_1a = table_1a_final, table_1b = table_1b))
}


###################################################### 
# ROBUST - Figs generation 
######################################################

source("/Users/savannah/Documents/GitHub/ch4.Meta_analysis/6_ext_heterogeneity-functions.R") #("6_ext_heterogeneity-functions.R")
#source("/Users/savannah/Documents/GitHub/ch4.Meta_analysis/6_ext_sampling-bias-functions.R")
source("/Users/savannah/Documents/GitHub/ch4.Meta_analysis/6_ext_visualization-results.R")
source("/Users/savannah/Documents/GitHub/ch4.Meta_analysis/6_ext_implementation.R") # Add the enhanced implementation to my workflow

# Create a directory for the enhanced results
enhanced_results_dir <- file.path(root, "enhanced_results")
if (!dir.exists(enhanced_results_dir)) {
  dir.create(enhanced_results_dir, recursive = TRUE)
}

# Path to the enhanced model file with robust model results
enhanced_model_file <- file.path(root, "Table_S3_Enhanced_GHG_reg_models_20yr_noPeat.csv")

# Check if the file exists
if (!file.exists(enhanced_model_file)) {
  stop("Enhanced model file not found: ", enhanced_model_file)
} else {
  cat("Found enhanced model file:", enhanced_model_file, "\n")
}

# Function to run the main analysis with enhanced/robust models

# Function to run the main analysis with enhanced/robust models
main_analysis_with_robust_models <- function(root, figsDir, enhanced_results_dir, model_data_file, enhanced_model_file) {
  # Create enhanced results directory if it doesn't exist
  if (!dir.exists(enhanced_results_dir)) {
    dir.create(enhanced_results_dir, recursive = TRUE)
    cat("Created enhanced results directory:", enhanced_results_dir, "\n")
  }
  
  # Read in the original model data
  cat("Reading original model data from:", model_data_file, "\n")
  model_data <- fread(model_data_file)
  
  # Read in the enhanced model data (with robust models)
  cat("Reading enhanced model data from:", enhanced_model_file, "\n")
  enhanced_model_data <- fread(enhanced_model_file)
  
  # Prepare the GHG data for Figure 2 using robust models where appropriate
  cat("Preparing GHG data for Figure 2 with robust models...\n")
  ag_file <- file.path(root, "__summ_Ag_CH4_N2O_xy_noRice.csv")
  if (!file.exists(ag_file)) {
    stop("Agricultural data file not found: ", ag_file)
  }
  
  ghg_data_robust <- prepare_ghg_data_robust(
    model_data,
    enhanced_model_data,
    ag_file = ag_file
  )
  
  # Print global means
  cat("\nGlobal means summary with robust models:\n")
  print_global_means(ghg_data_robust)
  
  # Create Figure 2 with robust models
  cat("Creating Figure 2 with robust models...\n")
  fig2_robust <- plot_ghg_comparison(ghg_data_robust, y_limit = 4)
  
  # Read in the global predictions data
  preds_file <- file.path(root, "global_preds_vMar2025.csv")
  if (!file.exists(preds_file)) {
    stop("Global predictions file not found: ", preds_file)
  }
  
  cat("Reading global predictions data...\n")
  pred_data <- fread(preds_file)
  
  # Update predictions using robust models where appropriate
  cat("Creating updated predictions with robust models...\n")
  updated_pred_data_robust <- create_global_predictions_robust(pred_data, model_data, enhanced_model_data)
  
  # Write the updated predictions to file
  updated_preds_file <- file.path(enhanced_results_dir, "global_preds_updated_robust.csv")
  cat("Writing updated predictions with robust models to:", updated_preds_file, "\n")
  fwrite(updated_pred_data_robust, updated_preds_file)
  
  # Create Figure 3 with robust models
  cat("Creating Figure 3 with robust models...\n")
  plots_robust <- plot_global_predictions(updated_pred_data_robust)
  
  # Arrange plots vertically
  fig3_robust <- gridExtra::grid.arrange(plots_robust$ch4, plots_robust$n2o, plots_robust$net,
                                         ncol = 1, heights = c(1, 1, 1))
  
  # Save the figures to enhanced results directory
  fig2_file_robust <- file.path(enhanced_results_dir, "figure2_robust.pdf")
  fig3_file_robust <- file.path(enhanced_results_dir, "figure3_robust.pdf")
  
  cat("Saving Figure 2 with robust models to:", fig2_file_robust, "\n")
  ggsave(fig2_file_robust, fig2_robust, width = 10, height = 12)
  
  cat("Saving Figure 3 with robust models to:", fig3_file_robust, "\n")
  ggsave(fig3_file_robust, fig3_robust, width = 10, height = 15)
  
  # Create summary tables using robust models
  cat("Creating summary tables with robust models...\n")
  tables_robust <- create_summary_tables_robust(ghg_data_robust)
  
  table1a_file_robust <- file.path(enhanced_results_dir, "table_1a_robust.csv")
  table1b_file_robust <- file.path(enhanced_results_dir, "table_1b_robust.csv")
  
  cat("Saving summary tables with robust models to:", table1a_file_robust, "and", table1b_file_robust, "\n")
  write.csv(tables_robust$table_1a, table1a_file_robust, row.names = FALSE)
  write.csv(tables_robust$table_1b, table1b_file_robust, row.names = FALSE)
  
  # Calculate global CH4-N2O effect per year using robust models
  cat("Calculating global CH4-N2O effects with robust models...\n")
  global_effects_robust <- calculate_global_effect_robust(updated_pred_data_robust, root, enhanced_results_dir)
  
  # Create a comparison of original vs robust model results
  cat("Creating comparison of original vs robust results...\n")
  
  # Compare global effect results
  # Original global effect results should be stored in global_effect_results.csv
  original_global_effect_file <- file.path(root, "global_effect_results.csv")
  if (file.exists(original_global_effect_file)) {
    original_global_effect <- fread(original_global_effect_file)
    
    # Calculate percent difference
    original_total <- original_global_effect$total_MgCO2e_yr[1]
    robust_total <- global_effects_robust$total_MgCO2e_yr
    
    pct_diff <- 100 * (robust_total - original_total) / abs(original_total)
    
    cat("Global Effect Comparison:\n")
    cat("Original:", format(original_total, big.mark=","), "Mg CO2e/yr\n")
    cat("With Robust Models:", format(robust_total, big.mark=","), "Mg CO2e/yr\n")
    cat("Percent Difference:", sprintf("%.2f%%", pct_diff), "\n")
    
    # Save comparison to file
    global_effect_comparison <- data.frame(
      Description = c("Original", "With Robust Models", "Absolute Difference", "Percent Difference"),
      Value = c(
        original_total,
        robust_total,
        robust_total - original_total,
        pct_diff
      )
    )
    
    write.csv(global_effect_comparison, 
              file.path(enhanced_results_dir, "global_effect_comparison.csv"),
              row.names = FALSE)
  }
  
  cat("Analysis with robust models complete!\n")
  
  # Return the results for inspection
  return(list(
    ghg_data_robust = ghg_data_robust,
    fig2_robust = fig2_robust,
    fig3_robust = fig3_robust,
    pred_data_robust = updated_pred_data_robust,
    tables_robust = tables_robust,
    global_effects_robust = global_effects_robust
  ))
}

# Run the enhanced analysis with robust models
cat("Running enhanced analysis with robust models...\n")
results_robust <- main_analysis_with_robust_models(
  root = root,
  figsDir = figsDir,
  enhanced_results_dir = enhanced_results_dir,
  model_data_file = file.path(root, "Table_S3_GHG_reg_models_20yr_noPeat_vMarch2025.csv"),
  enhanced_model_file = enhanced_model_file
)

# Create a robust models summary file
cat("Creating summary of robust model usage...\n")
robust_summary <- data.frame(
  GHG = c("CH4", "N2O"),
  Total_Models = c(0, 0),
  Robust_Models_Used = c(0, 0),
  Percent = c(0, 0)
)

# Count robust models used
enhanced_model_data <- fread(enhanced_model_file)
ch4_models <- enhanced_model_data[enhanced_model_data$GHG == "CH4", ]
n2o_models <- enhanced_model_data[enhanced_model_data$GHG == "N2O", ]

robust_summary$Total_Models[1] <- nrow(ch4_models)
robust_summary$Robust_Models_Used[1] <- sum(ch4_models$Used_Robust == TRUE, na.rm = TRUE)
robust_summary$Percent[1] <- 100 * robust_summary$Robust_Models_Used[1] / robust_summary$Total_Models[1]

robust_summary$Total_Models[2] <- nrow(n2o_models)
robust_summary$Robust_Models_Used[2] <- sum(n2o_models$Used_Robust == TRUE, na.rm = TRUE)
robust_summary$Percent[2] <- 100 * robust_summary$Robust_Models_Used[2] / robust_summary$Total_Models[2]

# Save robust summary
write.csv(robust_summary, 
          file.path(enhanced_results_dir, "robust_models_summary.csv"),
          row.names = FALSE)

cat("Enhanced analysis complete. Results saved to:", enhanced_results_dir, "\n")

# Print a comparison of key results
cat("\n====== COMPARISON OF ORIGINAL VS ROBUST MODELS ======\n")

# Compare global means
original_global_means <- ghg_data[ghg_data$Biome == "Global mean", ]
robust_global_means <- results_robust$ghg_data_robust[results_robust$ghg_data_robust$Biome == "Global mean", ]

cat("\nGlobal Mean Comparison:\n")
for (ghg_type in c("Eco_CH4", "Eco_N2O", "Eco_CH4-N2O")) {
  orig_val <- original_global_means$B0_backtr[original_global_means$GHG == ghg_type]
  robust_val <- robust_global_means$B0_backtr[robust_global_means$GHG == ghg_type]
  
  if (length(orig_val) > 0 && length(robust_val) > 0) {
    pct_diff <- 100 * (robust_val - orig_val) / abs(orig_val)
    cat(sprintf("%s: Original = %.4f, Robust = %.4f (%.2f%% difference)\n",
                ghg_type, orig_val, robust_val, pct_diff))
  }
}

cat("\nEnhanced analysis with robust models has been integrated successfully.\n")