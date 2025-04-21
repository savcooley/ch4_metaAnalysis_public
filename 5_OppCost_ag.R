#######   Overview & Install packages ######

#' Goal: CH4-N2O opportunity cost/benefit calculation: estimate biome-specific mean +/- SE CH4 and N2O fluxes for agricultural land

#' Steps:
#'   a) Retrieve data from Huddell et al (2020) and EDGARv7 (Crippa et al 2022)
#'   b) Data preparation & preliminary analysis 
#'   c) Estimate CH4 and N2O fluxes in agricultural land (from EDGARv8) 
#'   d) Create figures, including global maps & box plots of results (maybe will do box plots in another script)

library(data.table); library(sf); library(raster); library(terra); 
library(FNN); library(ggplot2); library(dplyr); 
library(sjPlot); library(tidyverse); library(ncdf4); library(colorspace); library(ggplot2);library(RColorBrewer)
root="/Users/savannah/Documents/Ch4/Data/Ag_opportunity_cost/"
figsDir="/Users/savannah/Documents/Ch4/Figs_and_tables/"

##    CH4: Total global Agriculture emissions w/ EDGARv8     #######  
# EDGARv8: Yearly Emissions gridmaps in ton substance / 0.1degree x 0.1degree / year 
# for the .txt files with longitude and latitude coordinates referring to the low-left corner of each grid-cell.
r_stack<-raster::stack() 
for(year in 2017:2022){ # v8.0_FT2022_GHG_GHG_CH4_2017_ENF.txt
  print(year)
  
  AWB_dt<-read.delim(sep = ";",skip=2,paste(root,"EDGARv8_CH4 - IPCC 4F - Agricultural waste burning/v8.0_FT2022_GHG_CH4_",year,"_AWB.txt", sep=""))
  AGS_dt<-read.delim(sep = ";",skip=2,paste(root,"EDGARv8_CH4 - IPCC 4C+4D - Agricultural soils/v8.0_FT2022_GHG_CH4_",year,"_AGS.txt", sep=""))
  ENF_dt<-read.delim(sep = ";",skip=2,paste(root,"EDGARv8_CH4 - IPCC 4A - Enteric fermentation/v8.0_FT2022_GHG_CH4_",year,"_ENF.txt", sep=""))
  MNM_dt<-read.delim(sep = ";",skip=2,paste(root,"EDGARv8_CH4 - IPCC 4B - Manure management/v8.0_FT2022_GHG_CH4_",year,"_MNM.txt", sep=""))
  
  dt<-merge(AWB_dt, AGS_dt, by=c("lat", "lon"), all=T)
  dt<-merge(dt, MNM_dt, by=c("lat", "lon"), all=T); names(dt)<-c("lat", "lon", "CH4_tons_AWB", "CH4_tons_AGS", "CH4_tons_NMN")
  dt<-merge(dt, ENF_dt, by=c("lat", "lon"), all=T); names(dt)<-c("lat", "lon", "CH4_tons_AWB", "CH4_tons_AGS", "CH4_tons_NMN", "CH4_tons_ENF")
  dt<-dt %>%
    mutate(sum_Ag_CH4 = rowSums(across( c(CH4_tons_AWB, CH4_tons_AGS, CH4_tons_NMN, CH4_tons_ENF)), na.rm=TRUE))
  
  rCH4_sum<-rasterFromXYZ(dt[, c("lon", "lat", "sum_Ag_CH4")]) # Convert table to raster
  raster::crs(rCH4_sum) <-"EPSG:4326" # deprecated PROJ4 string equivalent: "+proj=longlat +datum=WGS84 +no_defs", see https://gis.stackexchange.com/questions/374508/spproj4string-equivalent-in-the-proj6-framework
  
  r_stack<-raster::stack(r_stack, rCH4_sum)
}

# Calculate mean and SD for all years
rCH4_mean<- mean(r_stack, na.rm=T) 
#rCH4_sd<- sd(r_stack, na.rm = T) # didn't work for SD 
rCH4_sd <- calc(r_stack, fun = sd) # slow but works

plot(rCH4_mean, main="2017-2022 Mean Ag. CH4  (Mg CH4)") # ton substance / 0.1degree x 0.1degree / year 
writeRaster(rCH4_mean, paste(root,"_CH4_mean_ag_EDGARv8_2017-2022.tif", sep=""))
plot(rCH4_sd, main="2017-2022 SD Ag. CH4  (Mg CH4)")
writeRaster(rCH4_sd, paste(root,"_CH4_sd_ag_EDGARv8_2017-2022.tif", sep=""))

##    N2O: Total global Agriculture emissions w/ EDGARv8     #######  
# EDGARv8: Yearly Emissions gridmaps in ton substance / 0.1degree x 0.1degree / year for the .txt files with longitude and latitude coordinates referring to the low-left corner of each grid-cell.
r_stack<-raster::stack() 
for(year in 2017:2022){ 
  print(year)
  
  AWB_dt<-read.delim(sep = ";",skip=2,paste(root,"EDGARv8_N2O - IPCC 4F - Agricultural waste burning/v8.0_FT2022_GHG_N2O_",year,"_AWB.txt", sep=""))
  AGS_dt<-read.delim(sep = ";",skip=2,paste(root,"EDGARv8_N2O - IPCC 4C+4D - Agricultural soils/v8.0_FT2022_GHG_N2O_",year,"_AGS.txt", sep=""))
  INA_dt<-read.delim(sep = ";",skip=2,paste(root,"EDGARv8_N2O - IPCC 4D3 - Indirect N2O emissions from agriculture/v8.0_FT2022_GHG_N2O_",year,"_N2O.txt", sep=""))
  MNM_dt<-read.delim(sep = ";",skip=2,paste(root,"EDGARv8_N2O - IPCC 4B - Manure management/v8.0_FT2022_GHG_N2O_",year,"_MNM.txt", sep=""))
  
  dt<-merge(AWB_dt, AGS_dt, by=c("lat", "lon"), all=T)
  dt<-merge(dt, MNM_dt, by=c("lat", "lon"), all=T); names(dt)<-c("lat", "lon", "N2O_tons_AWB", "N2O_tons_AGS", "N2O_tons_NMN")
  dt<-merge(dt, INA_dt, by=c("lat", "lon"), all=T); names(dt)<-c("lat", "lon", "N2O_tons_AWB", "N2O_tons_AGS", "N2O_tons_NMN", "N2O_tons_INA")
  dt<-dt %>%
    mutate(sum_Ag_N2O = rowSums(across( c(N2O_tons_AWB, N2O_tons_AGS, N2O_tons_NMN, N2O_tons_INA)), na.rm=TRUE))
  
  rN2O_sum<-rasterFromXYZ(dt[, c("lon", "lat", "sum_Ag_N2O")]) # Convert table to raster
  raster::crs(rN2O_sum) <-"EPSG:4326" # deprecated PROJ4 string equivalent: "+proj=longlat +datum=WGS84 +no_defs", see https://gis.stackexchange.com/questions/374508/spproj4string-equivalent-in-the-proj6-framework
  
  r_stack<-raster::stack(r_stack, rN2O_sum)
}

# Test: "Global human-induced emissions, which are dominated by nitrogen additions to croplands

rN2O_tons_AGS<-rasterFromXYZ(dt[, c("lon", "lat", "N2O_tons_AGS")]) # Convert table to raster
plot(rN2O_tons_AGS, main="2022 Ag. soils N2O emissions  (Mg N2O)", xlim=400)

rN2O_tons_NMN<-rasterFromXYZ(dt[, c("lon", "lat", "N2O_tons_NMN")]) # Convert table to raster
plot(rN2O_tons_NMN, main="2022 Manure management N2O emissions  (Mg N2O)")

rN2O_tons_INA<-rasterFromXYZ(dt[, c("lon", "lat", "N2O_tons_INA")]) # Convert table to raster
plot(rN2O_tons_NMN, main="2022 Indirect Ag. N2O emissions  (Mg N2O)")



# Calculate mean and SD for all years
rN2O_mean<- mean(r_stack, na.rm=T) 
rN2O_sd <- calc(r_stack, fun = sd) # slow, but works

plot(rN2O_mean, main="2017-2022 Mean Ag. N2O  (Mg N2O)")
writeRaster(rN2O_mean, paste(root,"_N2O_mean_ag_EDGARv8_2017-2022.tif", sep=""))
plot(rN2O_sd, main="2017-2022 SD Ag. N2O  (Mg N2O)")
writeRaster(rN2O_sd, paste(root,"_N2O_sd_ag_EDGARv8_2017-2022.tif", sep=""))


##    Global map prep: CH4 and N2O Ag based on EDGAR only     #######  
#####     CH4 map     ####   
# EDGARv8: Yearly Emissions gridmaps in ton substance / 0.1degree x 0.1degree / year. 1 metric tonne = 1 Megagram 
rCH4_mean<-raster("/Users/savannah/Documents/Ch4/Data/Ag_opportunity_cost/_CH4_mean_ag_EDGARv8_2017-2022.tif")
Ag_CH4_CO2e<- rCH4_mean*86 # methane has a global warming potential of 25x CO2 over 100 years (GWP100 = 25) but 86 over 20 years (GWP20 = 86)
r <- rast(Ag_CH4_CO2e); crs(r)<-"epsg:4326"
a<- terra::cellSize(r, mask=T, unit="ha") # Calculate area of each pixel in hectares
minmax(a); plot(a, main="Area of each pixel (in ha) ")
a<- raster(a)#; raster::cellStats(a, sta="sum", na.rm=TRUE) # Sum area
Ag_CH4_CO2e<- Ag_CH4_CO2e/a # Convert from Mg 0.1x0.1 degree-1 yr-1 to Mg ha-1 yr-1

# Set rice cultivation areas to NA

tmp<- data.frame(rasterToPoints(Ag_CH4_CO2e, xy=T, na.rm=T)) #; dt<-subset(dt, layer<20); hist(dt$layer) 

## SD map
rCH4_sd<-raster("/Users/savannah/Documents/Ch4/Data/Ag_opportunity_cost/_CH4_sd_ag_EDGARv8_2017-2022.tif")
Ag_CH4_CO2e_sd<- rCH4_sd*86 # methane has a global warming potential of 25x CO2 over 100 years (GWP100 = 25) but 86 over 20 years (GWP20 = 86)
Ag_CH4_CO2e_sd<- Ag_CH4_CO2e_sd/a
tmp.SD = data.frame(rasterToPoints(Ag_CH4_CO2e_sd, xy=T, na.rm=T))

tmp<-merge(tmp, tmp.SD, by=c("x", "y"))

#####     
## Plot
tmp$cuts=cut(tmp[,3],breaks=c(cuts= c(seq(0, 5, 0.6), 10 ) )) # dt$cuts=cut(dt[,3],breaks=c(cuts= c(seq(0, 9e3, 11e2), 8e8 ) ))
ggplot(data=tmp) + 
  geom_tile(aes(x=x,y=y,fill=cuts)) + 
  scale_fill_manual(values = c("lightyellow","yellow", "gold",  "orange", "red","red2", "red3", "red4", "purple"), na.value = "white",
                    name = expression("Mg" ~ "ha"^-1~ "yr"^-1))+
  theme_bw() + coord_equal() + theme(panel.grid.major = element_blank()) +
  xlab("Longitude") + ylab("Latitude") + 
  ggtitle(expression("2017-2022 mean CH"[4]*"-CO"[2]*"e" ~ "agricluture soil flux"))


tmp.SD$cuts=cut(tmp.SD$layer,breaks=c(cuts= c(seq(0, 0.2, 0.03), 2) ))
ggplot(data=tmp.SD) + 
  geom_tile(aes(x=x,y=y,fill=cuts)) + 
  scale_fill_manual(values = terrain_hcl(9,  h=c(180,33),l=c(85,95), c=c(65,10), power=c(1/3,1.5) ), #brewer.pal(9, "Purples"), 
                    na.value = "white", name = expression("Mg" ~ "ha"^-1~ "yr"^-1))+ 
  theme_bw() + coord_equal() + theme(panel.grid.major = element_blank()) +
  xlab("Longitude") + ylab("Latitude") + 
  ggtitle(expression("2017-2022 SD CH"[4]*"-CO"[2]*"e" ~ "agricluture soil flux"))


#####     N2O map     ####   
rN2O_mean<-raster("/Users/savannah/Documents/Ch4/Data/Ag_opportunity_cost/_N2O_mean_ag_EDGARv8_2017-2022.tif")
rN2O_mean<-resample(rN2O_mean, rCH4_mean)
Ag_N2O_CO2e<- rN2O_mean*298 # N2O has a global warming potential of 298x CO2 over 100 years (GWP100 = 298) and also over 20 years (GWP20 = 298)
Ag_N2O_CO2e<-Ag_N2O_CO2e/a # Convert from Mg 0.1 degree-1 yr-1 to Mg ha-1 yr-1 by dividing Ag_CO2e raster by the Area of each pixel (in ha) 
t<- data.frame(rasterToPoints(Ag_N2O_CO2e, xy=T, na.rm=T))
tmp<-merge(tmp, t, by=c("x", "y"))

## SD map
rN2O_sd<-raster("/Users/savannah/Documents/Ch4/Data/Ag_opportunity_cost/_N2O_sd_ag_EDGARv8_2017-2022.tif")
rN2O_sd<-resample(rN2O_sd, rCH4_sd)
Ag_N2O_CO2e_sd<- rN2O_sd*298 # GWP 
Ag_N2O_CO2e_sd<- Ag_N2O_CO2e_sd/a
t.SD = data.frame(rasterToPoints(Ag_N2O_CO2e_sd)) #; dt<-subset(dt, layer<2); hist(dt$layer)

tmp<-merge(tmp, t.SD, by=c("x", "y"))

names(tmp)<- c("x", "y", "Ag_CH4_CO2e", "Ag_CH4_CO2e.sd", "Ag_N2O_CO2e", "Ag_N2O_CO2e.sd")

#####     
## Plot
t$cuts=cut(t$layer,breaks=c(cuts=c(seq(0, 2, 0.25), 5 ) )) # c(seq(0, 8, ), 8e8 )
ggplot(data=t) + 
  geom_tile(aes(x=x,y=y,fill=cuts)) + 
  scale_fill_manual(values = c("lightyellow","yellow", "gold",  "orange", "red","red2", "red3", "red4", "purple"),
                    na.value = "white", name = expression("Mg" ~ "ha"^-1~ "yr"^-1))+
  theme_bw() + coord_equal() + theme(panel.grid.major = element_blank()) +
  xlab("Longitude") + ylab("Latitude") + 
  ggtitle(expression("2017-2022 mean N"[2]*"O-CO"[2]*"e" ~ "agricluture soil flux"))


t.SD$cuts=cut(t.SD$layer,breaks=c(cuts= c(seq(0, 0.15, 0.03), 2) )) #set breaks. Max n is 9 for "Reds" palette
ggplot(data=t.SD) + 
  geom_tile(aes(x=x,y=y,fill=cuts)) + 
  scale_fill_manual(values = terrain_hcl(9,  h=c(180,33),l=c(85,95), c=c(65,10), power=c(1/3,1.5) ), #brewer.pal(9, "Purples"),
                    na.value = "white", name = expression("Mg" ~ "ha"^-1~ "yr"^-1))+ 
  theme_bw() + coord_equal() + theme(panel.grid.major = element_blank()) +
  xlab("Longitude") + ylab("Latitude") + 
  ggtitle(expression("2017-2022 SD N"[2]*"O-CO"[2]*"e" ~ "agricluture soil flux"))

#####     Fig 5A & B vOld: Global Ag Net CH4-N2O effect & aggregated biome-specific maps     ####   
# Calculate Net CH4-N2O effect and plot
Ag_net_CH4_N2O <- Ag_CH4_CO2e + Ag_N2O_CO2e
dt = data.frame(rasterToPoints(Ag_net_CH4_N2O))#; dt<-subset(dt, layer<10); hist(dt$layer)

## Add SD 
Ag_net_CH4_N2O_sd <- (Ag_CH4_CO2e_sd^2 + Ag_N2O_CO2e_sd^2)^0.5
dt_sd = data.frame(rasterToPoints(Ag_net_CH4_N2O_sd)) #; dt<-subset(dt, layer<1); hist(dt$layer)

dt<-merge(dt, dt_sd, by=c("x", "y"))

dt$cuts=cut(dt$layer,breaks=c(cuts= c(seq(0, 5, 0.6), 20 ) )) #set breaks. Max n is 9 for "Reds" palette
ggplot(data=dt) + 
  geom_tile(aes(x=x,y=y,fill=cuts)) + 
  scale_fill_manual(values = c("lightyellow","yellow", "gold",  "orange", "red","red2", "red3", "red4", "purple"),
                    na.value = "white", name = expression("CO"[2]*"e" ~ "soil flux (Mg" ~ "ha"^-1~ "yr"^-1*")"))+
  theme_bw() + coord_equal() + theme(panel.grid.major = element_blank()) +
  xlab("Longitude") + ylab("Latitude") + 
  ggtitle(expression("Mean net CH"[4]*"-N"[2]*"O effect from agriculture"))

dt_sd$cuts=cut(dt_sd$layer,breaks=c(cuts= c(seq(0, 0.5, 0.07), 0.8 ) )) #set breaks. Max n is 9 for "Reds" palette
ggplot(data=dt_sd) + 
  geom_tile(aes(x=x,y=y,fill=cuts)) + 
  scale_fill_manual(values = terrain_hcl(9,  h=c(180,33),l=c(85,95), c=c(65,10), power=c(1/3,1.5) ), #brewer.pal(9, "Purples"),
  na.value = "white", name = expression("CO"[2]*"e" ~ "soil flux (Mg" ~ "ha"^-1~ "yr"^-1*")") )+
  theme_bw() + coord_equal() + theme(panel.grid.major = element_blank()) +
  xlab("Longitude") + ylab("Latitude") + 
  ggtitle(expression("SD of net CH"[4]*"-N"[2]*"O effect from agriculture"))


##    Prep: define biome types v2  ######
# Import global biome data set (from 2017 resolve ecoregions) as data table https://resourcewatch.org/data/explore/bio042-Ecoregion-by-Biome?section=Discover&selectedCollection=&zoom=3&lat=0&lng=0&pitch=0&bearing=0&basemap=dark&labels=light&layers=%255B%257B%2522dataset%2522%253A%2522ed1544bb-a092-424e-88c2-8d548f4ef94a%2522%252C%2522opacity%2522%253A1%252C%2522layer%2522%253A%2522eea03513-208b-448a-a543-fa6056d356f3%2522%257D%255D&aoi=&page=1&sort=most-viewed&sortDirection=-1
r<-raster("/Users/savannah/Documents/Ch3/Data/Biomes_ecoregions/Biomes_11km_res.tif") #rgb <- brick(“pathto/rgb.tif”)
plot(r, col=bpy.colors(14)) #plotRGB(rgb) 
# crs_wgs84 <- st_crs(4326) # WGS 84 has EPSG code 4326 # Deprecated Proj.4 representation: +proj=longlat +datum=WGS84 +no_defs 

# Convert to data.table
b<-data.table(as.data.frame(r, xy=T, na.rm=T))
names(b)<-c("Lon", "Lat", "Biome_num")

# Add biome name column based on numeric biome ID
#b$Biome_num<-as.factor(b$Biome_num); levels(b$Biome_num)
b<- b %>%
  mutate(Biome_orig = case_when((Biome_num == 1) ~ 'Tropical and subtropical moist broadleaf forests',
                                (Biome_num == 2) ~ 'Tropical and subtropical dry broadleaf forests',
                                (Biome_num == 3) ~ 'Tropical and subtropical coniferous forests',
                                (Biome_num == 4) ~ 'Temperate broadleaf and mixed forests',
                                (Biome_num == 5) ~ 'Temperate conifer forests',
                                (Biome_num == 6) ~ 'Boreal forests or taiga',
                                (Biome_num == 7) ~ 'Tropical and subtropical grasslands, savannas, and shrublands',
                                (Biome_num == 8) ~ 'Temperate grasslands, savannas, and shrublands',
                                (Biome_num == 9) ~ 'Flooded grasslands and savannas',
                                (Biome_num == 10) ~ 'Montane grasslands and shrublands',
                                (Biome_num == 11) ~ 'Tundra',
                                (Biome_num == 12) ~ 'Mediterranean forests, woodlands, and scrub',
                                (Biome_num == 13) ~ 'Deserts and xeric shrublands',
                                (Biome_num == 14) ~ 'Mangroves',
                                TRUE ~ as.character(NA) )) 

# Re-categorize biomes based on C-P et al (2020).    
b<- b %>%
  mutate(Biome = case_when((Biome_num == 1) ~ 'Subtropical/tropical forest',
                           (Biome_num == 2) ~ 'Subtropical/tropical forest',
                           (Biome_num == 3) ~ 'Subtropical/tropical forest',
                           (Biome_num == 4) ~ 'Temperate broadleaf forest',
                           (Biome_num == 5) ~ 'Temperate conifer forest',
                           (Biome_num == 6) ~ 'Boreal forest',
                           (Biome_num == 7) ~ 'Subtropical/tropical savanna',
                           (Biome_num == 8) ~ NA,
                           (Biome_num == 9) ~ NA, # Left category as-is, Will not include
                           (Biome_num == 10) ~ NA, #'Temperate grassland and savanna',# a bit of a stretch for: 'Montane grasslands and shrublands',
                           (Biome_num == 11) ~ NA,
                           (Biome_num == 12) ~ 'Subtropical/tropical savanna',# a bit of a stretch for: 'Mediterranean forests, woodlands, and scrub',
                           (Biome_num == 13) ~ NA,#'Deserts and xeric shrublands', # Left category as-is, Will not include
                           (Biome_num == 14) ~ NA,#'Mangroves', # Left category as-is, Will not include
                           TRUE ~ as.character(NA) )) 

# Spatial join via KNN with X,Y points: return the k points in b that are nearest to each point in dt
b$rNum<-1:nrow(b)
nn = get.knnx(subset(b, select = c("Lon", "Lat")), 
              subset(dt, select = c("x", "y")), k=1, algorithm="kd_tree") # Lat lon columns only for both data sets

# Join data using index from nn
dt$index<-nn$nn.index # an n x k matrix for the nearest neighbor index
dt$dist<-nn$nn.dist 	# an n x k matrix for the nearest neighbor Euclidean distances.
dt<- merge(dt, b, by.x="index", by.y="rNum")
dt$index<-NULL; dt$dist<-NULL; dt$Lon<-NULL; dt$Lat<-NULL; names(dt); 
colnames(dt)[3] <- "Ag_net_CH4_N2O" # dt$Ag_net_CH4_N2O already in units of Mg ha-1 yr-1
colnames(dt)[4] <- "sd_Ag_net_CH4_N2O" 

#dt$Biome<-as.factor(dt$Biome); levels(dt$Biome); dt %>% plot_frq(Biome)

dt<-merge(dt, tmp, by=c("x", "y"))

fwrite(dt, "/Users/savannah/Documents/Ch3/Data/Ag_opportunity_cost/__dt_Ag_CH4_N2O_xy.csv")

##    Plot: aggregated Ag_net_CH4_N2O by biome   ####### 

dt<-fread("/Users/savannah/Documents/Ch3/Data/Ag_opportunity_cost/__dt_Ag_CH4_N2O_xy.csv")
dt$Ag_net_CH4_N2O[dt$Ag_net_CH4_N2O < 0.1] <- NA # adjust data to only compute mean for values > 0.1 Mg ha-1 yr-1
dt$Ag_net_CH4_N2O[dt$Biome == "Tundra"] <- NA # set Tundra to NA

dt <- dt %>%
  group_by(Biome_orig) %>%
  mutate(mean_Ag_net_CH4_N2O = mean(Ag_net_CH4_N2O, na.rm=T)) 

dt <- dt %>%
  group_by(Biome_orig) %>%
  mutate(SD_Ag_net_CH4_N2O = sd(Ag_net_CH4_N2O, na.rm=T)) # mean(Ag_net_CH4_N2O, na.rm=T) + 

dt %>%
  group_by(Biome) %>%
  summarise(meanAg = mean(Ag_net_CH4_N2O, na.rm=T), n=n()) 

fwrite(dt, "/Users/savannah/Documents/Ch3/Data/Ag_opportunity_cost/__dt_Ag_CH4_N2O_xy.csv")

# Plot aggregated Ag_net_CH4_N2O by biome
#par(mfrow = c(2, 2)) 

df<-subset(dt, select=c("x", "y", "mean_Ag_net_CH4_N2O"))
df$cuts=cut(df$mean_Ag_net_CH4_N2O,breaks=c(cuts= c(seq(0, 5, 0.6), 20 ) )) #set breaks. Max n is 9 for "Reds" palette

ggplot(data=df) + 
  geom_tile(aes(x=x,y=y,fill=cuts)) + 
  scale_fill_manual(values = c("lightyellow","yellow", "gold",  "orange", "red","red2", "red3", "red4", "purple"),
                    na.value = "white", name = expression("CO"[2]*"e" ~ "soil flux (Mg" ~ "ha"^-1~ "yr"^-1*")"))+
  theme_bw() + coord_equal() + theme(panel.grid.major = element_blank()) +
  xlab("Longitude") + ylab("Latitude") + 
  ggtitle(expression("Biome mean net CH"[4]*"-N"[2]*"O effect from ag."))

#mean_Ag_net_CH4_N2O <- rasterFromXYZ(dt[, c("Lon", "Lat", "mean_Ag_net_CH4_N2O")]) 
#plot(mean_Ag_net_CH4_N2O, col=diverge_hcl(8, c=100, l=c(50,90), power=1), main=expression("Mean Net CH"[4]*"-N"[2]*"O effect (CO"[2]*"e" ~ "from agriculture)")) 

df<-subset(dt, select=c("x", "y", "SD_Ag_net_CH4_N2O"))
df$cuts=cut(df$SD_Ag_net_CH4_N2O,breaks=c(cuts= c(seq(0, 120, 20), 250 ) )) #set breaks. Max n is 9 for "Reds" palette

ggplot(data=df) + 
  geom_tile(aes(x=x,y=y,fill=cuts)) + 
  scale_fill_manual(values = terrain_hcl(6,  h=c(180,33),l=c(85,95), c=c(65,10), power=c(1/3,1.5)), #brewer.pal(9, "Purples"), 
                    na.value = "white", name = expression("Mg" ~ "ha"^-1~ "yr"^-1))+ 
  theme_bw() + coord_equal() + theme(panel.grid.major = element_blank()) +
  xlab("Longitude") + ylab("Latitude") + 
  ggtitle(expression("Biome SD net CH"[4]*"-N"[2]*"O effect from ag."))

# difference between net CH4-N2O effect of ecosystem regeneration vs. the predicted net CH4-N2O-CO2 effects for alternate agriculture land use.

################     v2 Alternate ag extent: Higher Ag Net CH4-N2O effect with MODIS IGBP classification     ####   
#' End goal: global map of the Ag_net_CH4_N2O effect with pasture assumed to occupy 26% of "grassland" class from MODIS IGBP classification  
#' 
#' All pixels with non-zero Ag_net_CH4_N2O effect outside of IGBP cropland zones will be re-distributed uniformly 
#' across newly defined pasture lands that cover a smaller area than before (assumed to occupy 26% of IGBP natural grasslands class). 
#' This will result in larger Ag_net_CH4_N2O effects per unit area. 
 
#' LC class 12 Croplands: at least 60% of area is cultivated cropland.
#' LC class 14 Cropland/Natural Veg Mosaics: mosaics of small-scale cultivation 40-60% with natural tree, shrub, or herbaceous vegetation.
#' LC class 10 Grasslands: dominated by herbaceous annuals (<2m). Pasture is assumed to be 26% of grasslands based on FAO estimate also used by Powell group

# From Powell group: "The IGBP land cover system does not differentiate between natural grasslands and cultivated pastures. 
# We split the land area of grasslands between 67°N/S (25) into natural vs. cultivated per grid cell based on the estimate that cultivated pastures covered 26% of global grasslands in 2019 (37)."
# 
# 25. G. Carlsson, & K. Huss-Danell, Nitrogen fixation in perennial forage legumes in the field. Plant Soil 253, 353–372 (2003).
# 37. FAO, FAOSTAT. License: CC BY-NC-SA 3.0 IGO. (2021). https://www.fao.org/faostat/en/#data.
# 49. M. A. Friedl, D. Sulla-Menashe, MCD12Q1 MODIS/Terra+Aqua Land Cover Type Yearly L3 Global 500m SIN Grid V006 [Data set]. NASA EOSDIS Land Processes Distributed Active Archive Center. https://doi.org/10.5067/MODIS/MCD12Q1.006.
# 
# The approach was to treat each of the cells within the latitude bounds as 27% pasture and 73% natural grassland.

#' 1. Spatial join MODIS IGBP with Ag_net_CH4_N2O raster

dt<-fread("/Users/savannah/Documents/Ch3/Data/Ag_opportunity_cost/__dt_Ag_CH4_N2O_xy.csv")
Ag_net_CH4_N2O <- rasterFromXYZ(dt[, c("x", "y", "Ag_net_CH4_N2O")]) 
#writeRaster(Ag_net_CH4_N2O,"/Users/savannah/Documents/Ch3/Data/Ag_opportunity_cost/_Ag_net_CH4_N2O_xy.tif")

IGBP<- raster("/Users/savannah/Documents/Ch3/Data/Ag_opportunity_cost/LandCover_classificaiton/IGBP_MODIS_2021_1km.tif") # 500m MODIS Land Cover Type (MCD12Q1) Version 6.1 data product: ee.ImageCollection("MODIS/061/MCD12Q1")
IGBP<- resample(IGBP, Ag_net_CH4_N2O, method="ngb")
b<-data.table(as.data.frame(IGBP, xy=T, na.rm=T))
names(b)<-c("Lon", "Lat", "LC_num")

b<- b %>%
  mutate(LC_class = case_when((LC_num == 1) ~ 'Evergreen Needleleaf Forests',
                           (LC_num == 2) ~ 'Evergreen Broadleaf Forests',
                           (LC_num == 3) ~ 'Deciduous Needleleaf Forests',
                           (LC_num == 4) ~ 'Deciduous Broadleaf Forests',
                           (LC_num == 5) ~ 'Mixed Forests', 
                           (LC_num == 6) ~ 'Closed Shrublands',
                           (LC_num == 7) ~ 'Open Shrublands',
                           (LC_num == 8) ~ 'Woody Savannas',
                           (LC_num == 9) ~ 'Savannas', 
                           (LC_num == 10) ~ 'Grasslands', 
                           (LC_num == 11) ~ 'Permanent Wetlands',
                           (LC_num == 12) ~ 'Croplands', 
                           (LC_num == 13) ~ 'Urban and Built-up Lands', 
                           (LC_num == 14) ~ 'Cropland/Natural Vegetation Mosaics', 
                           (LC_num == 15) ~ 'Permanent Snow and Ice', 
                           (LC_num == 16) ~ 'Barren', 
                           (LC_num == 17) ~ 'Water Bodies', 
                           TRUE ~ as.character(NA) )) 

# Thought I could directly merge since I resampled raster already (ie. No need to spatial join). But this resulted in empty dt w 0 obs
#dt<- merge(dt, b, by.x=c("x", "y"), by.y=c("Lon", "Lat"))
#b<- merge(b, dt, by=c("Lon", "Lat", "LC_num", "LC_class"), all.x=T, all.y=T)

# Spatial join via KNN with X,Y points: return the k points in dt that are nearest to each point in b
dt$rNum<-1:nrow(dt)
nn = get.knnx(subset(dt, select = c("Lon", "Lat")),
              subset(b, select = c("x", "y")), k=1, algorithm="kd_tree") # Lat lon columns only for both data sets

# Join data using index from nn
b$index<-nn$nn.index # an n x k matrix for the nearest neighbor index
b$dist<-nn$nn.dist 	# an n x k matrix for the nearest neighbor Euclidean distances.
b<- merge(b, b, by.x="index", by.y="rNum")
b$dist<-NULL; b$index<-NULL

b$LC_class<-as.factor(b$LC_class); levels(b$LC_class); #b %>% plot_frq(Biome)
fwrite(b, "/Users/savannah/Documents/Ch3/Data/Ag_opportunity_cost/__dt_Ag_CH4_N2O_xy.csv")

#####     
#' 2. Sum all Ag_net_CH4_N2O pixels in non-cropland areas as classified by MODIS IGBP
dt<-fread("/Users/savannah/Documents/Ch4/Data/Ag_opportunity_cost/__dt_Ag_CH4_N2O_xy.csv")
#dt$Ag_net_CH4_N2O <- dt$Ag_CH4_CO2e + dt$Ag_N2O_CO2e # in [Mg CO2e ha-1 yr-1]

# cropland pixels are MODIS IGBP LC classes 12, 14
notCrop<-subset(dt, LC_num != 12); notCrop<-subset(notCrop, LC_num != 14); 

sum(notCrop$Ag_CH4_CO2e, na.rm = T); sum(dt$Ag_CH4_CO2e, na.rm = T) # CH4: 728340.9, 1311368
sum(notCrop$Ag_N2O_CO2e, na.rm = T); sum(dt$Ag_N2O_CO2e, na.rm = T) # N2O: 101610.5, 185072.5. Sum net effect: 829951.4 in [Mg CO2e ha-1 yr-1]
# sum(notCrop$Ag_net_CH4_N2O, na.rm = T) # Sanity check, same Sum net effect as above: 829951.4 also in [Mg CO2e ha-1 yr-1]

#' 3. Spatially re-distribute the non-cropland ag fluxes across pasture pixels (pasture is assumed to occupy ~50% of each class 10 pixel in latitude bounds +/- 57 degrees)
#'    ~40% pasture in grassland was adjusted from FAO's 26% bc of change in lat bounds
grass57 <- subset(dt, LC_num==10); nrow(grass57) # [1] 674704
grass57 <- subset(grass57, Lat >= -57); grass57 <- subset(grass57, Lat <= 57); nrow(grass57)

# 57 lat: 534367/674704 = 0.7920021. So 79% of all grasslands are in +/- 57 degrees lat

# .26/1 = x/0.7920021 --> x = 0.2346. 
# So new pasture % is: 26 + 23.46 = 49.46 % of all grassland

# Calculate new flux for +/- 57 Lat grassland pixels: 
# divide total non-cropland ag. flux sum by the number of +/- 57 Lat grassland pixels (assumed to be partial-pasture for Fig 5)
grass57_CH4_CO2e<- 100/49.46*sum(notCrop$Ag_CH4_CO2e, na.rm = T) / nrow(grass57) #57 lat: 2.73, 67 lat: 2.483838 [Mg CO2e ha-1 yr-1]  #*100/26 to get accurate per ha flux of pasture-only for Figs 2 and 4
grass57_N2O_CO2e<- 100/49.46*sum(notCrop$Ag_N2O_CO2e, na.rm = T) / nrow(grass57) #57 lat: 0.38, 67 lat: 0.346519 [Mg CO2e ha-1 yr-1] 

# set all nonCropland fluxes to NA
dt<- dt %>% mutate(Ag_CH4_CO2e = case_when(LC_num == 12 | LC_num == 14 ~ Ag_CH4_CO2e,
                                           TRUE ~ NA)) 
dt<- dt %>% mutate(Ag_N2O_CO2e = case_when(LC_num == 12 | LC_num == 14 ~ Ag_N2O_CO2e,
                                           TRUE ~ NA)) 

# Plot cropland-only ghg effects: dt$Ag_net_CH4_N2O <- dt$Ag_CH4_CO2e+dt$Ag_N2O_CO2e

# then replace grassland pixels between +/- 57 degrees Lat with newly calculated pasture flux (from grass57_CH4_CO2e and grass57_N2O_CO2e)
dt<- dt %>% mutate(Ag_CH4_CO2e = case_when(Lat <= 57 & Lat >= -57 & LC_num == 10 ~ grass57_CH4_CO2e,
                                           TRUE ~ Ag_CH4_CO2e)) 
dt<- dt %>% mutate(Ag_N2O_CO2e = case_when(Lat <= 57 & Lat >= -57 & LC_num == 10 ~ grass57_N2O_CO2e,
                                           TRUE ~ Ag_N2O_CO2e)) 

# Plot cropland+pasture ghg effects: dt$Ag_net_CH4_N2O <- dt$Ag_CH4_CO2e+dt$Ag_N2O_CO2e

sum(notCrop$Ag_CH4_CO2e, na.rm = T); sum(dt$Ag_CH4_CO2e, na.rm = T) # CH4: 728340.9, 2055613 (NOT the same as initial sum - sanity check??)
sum(notCrop$Ag_N2O_CO2e, na.rm = T); sum(dt$Ag_N2O_CO2e, na.rm = T) # N2O: 101610.5, 288901.8 (NOT same as original: 185072.5) 

dt$Ag_net_CH4_N2O <- dt$Ag_CH4_CO2e + dt$Ag_N2O_CO2e # in [Mg CO2e ha-1 yr-1]

fwrite(dt, "/Users/savannah/Documents/Ch3/Data/Ag_opportunity_cost/__dt_Ag_CH4_N2O_xy_pastureV2.csv")

#####     v3 post-GCB: exclude rice cultivation areas for CH4     ####   

# Read rice cultivation data
rice<- raster("/Users/savannah/Documents/Ch4/Data/Ag_opportunity_cost/CROPGRIDSv1.08_rice.nc",
              varname = "croparea") # varname can set to one of: harvarea, croparea, qual, set

dt_v2<- fread("/Users/savannah/Documents/Ch4/Data/Ag_opportunity_cost/__dt_Ag_CH4_N2O_xy_pastureV2.csv")

# Load required libraries
library(raster)
library(data.table)

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

dt_v2 %>%
  group_by(Biome) %>%
  summarise(
    #CH4_pred.mean = mean(Ag_CH4_CO2e_noRice, na.rm=TRUE),
    #N2O_pred.mean = mean(Ag_N2O_CO2e, na.rm=TRUE),
    #net_pred.mean = mean(Ag_noRice_CH4_N2O, na.rm=TRUE),
    net_pred.sd = sd(Ag_noRice_CH4_N2O, na.rm=TRUE))

#####     v3 Fig 5A & B: Global Ag Net CH4-N2O effect & aggregated biome-specific maps      #####   
#dt_v2<- fread("/Users/savannah/Documents/Ch4/Data/Ag_opportunity_cost/__dt_Ag_CH4_N2O_xy_pastureV2.csv")
dt<- fread("/Users/savannah/Documents/Ch4/Data/Ag_opportunity_cost/__dt_Ag_CH4_N2O_xy_noRice.csv")
str(dt); summary(dt)

summary(dt$Ag_noRice_CH4_N2O)

dt$cuts=cut(dt$Ag_net_CH4_N2O, breaks= c(cuts= c(0,0.6, 1.2, 1.8, 2.4, 3, 3.6, 4.2, 4.8, 20 ) )) 
              #c(cuts= c(0, seq(2, 6.4, 0.6), 20 ) )) #c(0, 1, 3, 5, 7, 12, 20 ) )) # Hard to get this right . OG:c(0, seq(2, 6.4, 0.6), 20 ) . Mostly red: c(cuts= c(0,seq(1.5, 5, 0.5), 20 ) )) #set breaks. Max n is 9 for "Reds" palette
ggplot(data=dt) + 
  geom_tile(aes(x=x,y=y,fill=cuts)) + 
  scale_fill_manual(values = c("lightyellow","yellow", "gold",  "orange", "red","red2", "red3", "red4", "purple"),
                    na.value = "grey88", name = expression("CO"[2]*"e" ~ "soil flux (Mg" ~ "ha"^-1~ "yr"^-1*")"))+
  theme_bw() + coord_equal() + theme(panel.grid.major = element_blank()) +
  xlab("Longitude") + ylab("Latitude") + 
  ggtitle(expression("Mean net CH"[4]*"-N"[2]*"O effect from agriculture"))

dt %>%
  group_by(Biome) %>%
  summarise(mean_Ag_noRice_CH4_N2O = mean(Ag_noRice_CH4_N2O, na.rm=T), n=n()) 
dt <- dt %>%
  group_by(Biome) %>%
  mutate(mean_Ag_noRice_CH4_N2O = mean(Ag_noRice_CH4_N2O, na.rm=T)) 

# reclassify to NA all pixels in mean_Ag_noRice_CH4_N2O that spatially overlap NA pixels in Ag_noRice_CH4_N2O
#dt$mean_Ag_noRice_CH4_N2O[dt$Ag_noRice_CH4_N2O == NA ] <- NA
#dt$mean_Ag_noRice_CH4_N2O[dt$Biome == ""] <- NA 

df<-subset(dt, select=c("x", "y", "mean_Ag_noRice_CH4_N2O")) 
# ag_CH4_N2O means: c(3.05, 2.53, 2.06, 1.55, 1.87,2.212): 
# 1.4,(trop forest) 1.6,(boreal) 2, (savanna) 2.4, (temp broadleaf), 2.8, (temp conifer), 3.2
df$cuts=cut(df$mean_Ag_noRice_CH4_N2O,breaks=c(cuts= c(1.4, 1.8, 2, 2.2,2.6,3.0) ))#seq(1.4, 3.2, 0.4)  )) #c(seq(0, 3, 0.6), 20) set breaks. Max n is 9 for "Reds" palette

ggplot(data=df) + 
  geom_tile(aes(x=x,y=y,fill=cuts)) + 
  scale_fill_manual(values = c("lightyellow","yellow", "gold",  "orange", "red","red2", "red3", "red4", "purple"),
                    na.value = "gray88", name = expression("CO"[2]*"e" ~ "soil flux (Mg" ~ "ha"^-1~ "yr"^-1*")"))+
  theme_bw() + coord_equal() + theme(panel.grid.major = element_blank()) +
  xlab("Longitude") + ylab("Latitude") + 
  ggtitle(expression("Biome mean net CH"[4]*"-N"[2]*"O effect from ag."))





#####     v4 Fig 5A & B, S3, S7, S8 - vMarch2025: Ag_noRice + Global Ag-adjusted Net CH4-N2O effect    ####   
# Bring in ecosystem regeneration net_CH4_N2O effect results
#b<-fread(paste(root,"Net_ch4_n2o_xy_results_Jan4.csv")) 
b<- fread("/Users/savannah/Documents/Ch4/Data/Data_extraction/global_preds_vMarch2025.csv") #global_preds_vNov24.csv")
#b$Eco_net.CH4.N2O.CO2e<- b$CH4.CO2e.BackTr + b$N2O.CO2e.BackTr
b<-subset(b, select=c(Lon, Lat, CH4.CO2e.BackTr, N2O.CO2e.BackTr, net_effect)) # ,Eco_net.CH4.N2O.CO2e))
str(b)

# Bring in Agriculture data
#dt<-fread("/Users/savannah/Documents/Ch3/Data/Ag_opportunity_cost/__dt_Ag_CH4_N2O_xy_pastureV2.csv")
dt<- fread("/Users/savannah/Documents/Ch4/Data/Ag_opportunity_cost/__dt_Ag_CH4_N2O_xy_noRice.csv")
str(dt)

# Spatial join via KNN with X,Y points: return the k points in b that are nearest to each point in dt
b$rNum<-1:nrow(b)
nn = get.knnx(subset(b, select = c("Lon", "Lat")),
              subset(dt, select = c("Lon", "Lat")), k=1, algorithm="kd_tree") # Lat lon columns only for both data sets

# Join data using index from nn
dt$index<-nn$nn.index # an n x k matrix for the nearest neighbor index
dt$dist<-nn$nn.dist 	# an n x k matrix for the nearest neighbor Euclidean distances.
dt<- merge(dt, b, by.x="index", by.y="rNum")
dt$index<- NULL; #dt$dist<- NULL; dt$X<- NULL; dt$Y<- NULL
hist(dt$dist, breaks = c(seq(0,5,0.5),50 )) 


nn <- subset(dt, select=c("x", "y", "Ag_noRice_CH4_N2O"))
p5a <- ggplot(data=nn) + 
  geom_tile(aes(x=x, y=y, fill=Ag_noRice_CH4_N2O)) +  # Remove the cuts and use continuous values
  scale_fill_gradientn(
    colors = c("white", "khaki1", "gold", "orange", "darkred"),
    values = scales::rescale(c(0, 1, 2, 4)),
    limits = c(0, 4),
    na.value = "gray88",
    name = expression("CO"[2]*"e soil flux (Mg" ~ ha^-1 ~ yr^-1*")")
  ) +
  theme_minimal() + 
  coord_equal() + 
  theme(panel.grid.major = element_blank()) +
  xlab("Longitude") + 
  ylab("Latitude") + 
  ggtitle(expression("Net CH"[4]*"-N"[2]*"O effect from agriculture"))


fig_file <- file.path(figsDir, "figure5a.pdf")
ggsave(fig_file, p5a, width = 12, height = 5)


# dt <- dt %>%
#   group_by(Biome) %>%
#   mutate(Ag_noRice_CH4_N2O_biomeMean = mean(Ag_noRice_CH4_N2O, na.rm = TRUE)) %>%
#   ungroup()

# Revised approach: only rows which =NA for Ag_noRice_CH4_N2O take on biome mean values. All other non-NA Ag ghg flux rows remain the same as original Ag_noRice_CH4_N2O
dt <- dt %>%
  group_by(Biome) %>%
  mutate(Ag_noRice_CH4_N2O_biomeMean =
           case_when(
             is.na(Ag_noRice_CH4_N2O) ~ mean(Ag_noRice_CH4_N2O, na.rm = TRUE), # 0.0000001, #
             TRUE ~ Ag_noRice_CH4_N2O
           )) %>%
  ungroup()

dt <- dt %>%
  group_by(Biome) %>%
  mutate(Ag_noRice_CH4_N2O_biomeMean_mask =
           case_when(
             is.na(net_effect) ~ NA,
             TRUE ~ Ag_noRice_CH4_N2O_biomeMean
           )) %>%
  ungroup()

dt$opp_benefit_CH4_N2O <- dt$net_effect - dt$Ag_noRice_CH4_N2O_biomeMean_mask

par(mfrow = c(1, 1))

df <- subset(dt, select=c("x", "y", "opp_benefit_CH4_N2O"))
p5b<- ggplot(data=df) + 
  geom_tile(aes(x=x, y=y, fill=opp_benefit_CH4_N2O)) +
  scale_fill_gradientn(
    colors = c("darkgreen", "green4", "green3", "palegreen", "white", "khaki1", "gold", "orange"), #, "darkred"
    values = scales::rescale(c(-8, -6, -2, 0, 0.5, 1, 3, 6)),  # Adjusted -1 to 0 and added 0.5
    limits = c(-8, 6),
    na.value = "gray88",
    name = expression("CO"[2]*"e soil flux (Mg" ~ ha^-1 ~ yr^-1*")")
  ) +
  theme_minimal() + 
  coord_equal() + 
  theme(panel.grid.major = element_blank(),
        strip.text.x = element_text(size = 18, colour = "black", face = "bold" )) +
  xlab("Longitude") + 
  ylab("Latitude") + 
  ggtitle(expression("Agriculture-adjusted Net CH"[4]*"-N"[2]*"O effect"),size = 18, colour = "black", face = "bold" )

fig_file <- file.path(figsDir, "figure5b.pdf")
ggsave(fig_file, p5b, width = 12, height = 5)



#####     v4 Fig S7, S8 - vMarch2025: enhanced plotting    ####   

# Figure S7. 
# Mean and SD 2017-2022 CH4 and N2O soil fluxes from agriculture (original EDGAR data)

df<-fread("/Users/savannah/Documents/Ch4/Data/Ag_opportunity_cost/__dt_Ag_CH4_N2O_xy.csv")

p_ch4 <- ggplot(data=df) + 
  geom_tile(aes(x=x, y=y, fill=Ag_CH4_CO2e)) +  
  scale_fill_gradientn(
    colors = c("gray88", "khaki1", "gold", "orange", "darkred", "purple"),
    values = scales::rescale(c(0, 1, 2, 4,6, 8)),
    limits = c(0, 8), na.value = "gray88", # ? not showing up, so changed white to gray88 for 0 values (some of which correspond to NA on land, but existing NAs are ocean and thus I didn't want to lump all 0s into NA category)
    name = expression("CO"[2]*"e soil flux (Mg" ~ ha^-1 ~ yr^-1*")")) +
  theme_minimal() + coord_equal() + theme(panel.grid.major = element_blank()) +
  xlab("Longitude") +   ylab("Latitude") + 
  labs(title = expression("Mean 2017-2022 CH"[4]*"-CO"[2]*"e" ~ "agricluture soil flux"))

p_n2o <- ggplot(data=df) + 
  geom_tile(aes(x=x, y=y, fill=Ag_N2O_CO2e)) +  
  scale_fill_gradientn(
    colors = c("gray88", "khaki1", "gold", "orange", "darkred", "purple"),
    values = scales::rescale(c(0, 1, 2.5)),
    limits = c(0, 2.5), na.value = "gray88",
    name = expression("CO"[2]*"e soil flux (Mg" ~ ha^-1 ~ yr^-1*")")) +
  theme_minimal() + coord_equal() + theme(panel.grid.major = element_blank()) +
  xlab("Longitude") +   ylab("Latitude") + 
  labs(title = expression("Mean 2017-2022 N"[2]*"O-CO"[2]*"e" ~ "agricluture soil flux"))

p_ch4_SD <- ggplot(data=df) + 
  geom_tile(aes(x=x, y=y, fill=Ag_CH4_CO2e.sd)) +  
  scale_fill_gradientn(
    colors = c("gray88", "lightblue2",  "deepskyblue", "deepskyblue2","deepskyblue4"),
    values = scales::rescale(c(0, 0.2, 0.4)),
    limits = c(0, 0.4), na.value = "gray88",
    name = expression("CO"[2]*"e soil flux (Mg" ~ ha^-1 ~ yr^-1*")")) +
  theme_minimal() + coord_equal() + theme(panel.grid.major = element_blank()) +
  xlab("Longitude") +   ylab("Latitude") + 
  labs(title = expression("SD 2017-2022 CH"[4]*"-CO"[2]*"e" ~ "agricluture soil flux"))

p_n2o_SD <- ggplot(data=df) + 
  geom_tile(aes(x=x, y=y, fill=Ag_N2O_CO2e.sd)) +  
  scale_fill_gradientn(
    colors = c("gray88", "lightblue2",  "deepskyblue", "deepskyblue2","deepskyblue4"),
    values = scales::rescale(c(0, 0.05, 0.1)),
    limits = c(0, 0.1), na.value = "gray88",
    name = expression("CO"[2]*"e soil flux (Mg" ~ ha^-1 ~ yr^-1*")")) +
  theme_minimal() + coord_equal() + theme(panel.grid.major = element_blank()) +
  xlab("Longitude") +   ylab("Latitude") + 
  labs(title = expression("SD 2017-2022 N"[2]*"O-CO"[2]*"e" ~ "agricluture soil flux"))

fig_file <- file.path(figsDir, "figS7.pdf")
figS7 <- gridExtra::grid.arrange(p_ch4, p_n2o, p_ch4_SD, p_n2o_SD, 
                                ncol = 2, heights = c(4, 4, 4, 4))

ggsave(fig_file, figS7, width = 18, height = 10)

#####
# Figure S8. Global maps of the net CH4-N2O effect in CO2-equivalent soil flux per unit of agricultural land per year from 
#' (A) the original EDGAR v8 inventory data averaged from 2017-2022 (same as Figure 5A); 
#' (B) a spatial re-distribution of all agricultural CH4 and N2O observations outside of IGBP cropland zones based on our estimate of the actual extent of pasture lands (see methods); 
#' (C) a removal of all CH4 fluxes corresponding to rice cultivation (see methods); and 
#' (D) the final data set with the pasture adjustment, the removal of areas with rice cultivation and NA pixels gap-filled with the biome-specific mean CH4-N2O effect. 
#' We used data shown in (D) in the subtraction to calculate the net CH4-N2O effect of regeneration relative to agriculture. 

# Fig 8 (A)
p <- ggplot(data=df) + 
  geom_tile(aes(x=x, y=y, fill=Ag_net_CH4_N2O)) +  
  scale_fill_gradientn(
    colors = c("white", "khaki1", "gold", "orange", "darkred", "purple"),
    values = scales::rescale(c(0, 1, 2, 4, 20)),
    limits = c(0, 20), na.value = "gray88",
    name = expression("CO"[2]*"e soil flux (Mg" ~ ha^-1 ~ yr^-1*")")) +
  theme_minimal() + coord_equal() + theme(panel.grid.major = element_blank()) +
  xlab("Longitude") +   ylab("Latitude") + 
  labs(title = expression("Net CH"[4]*"-N"[2]*"O effect from agriculture"),
       subtitle="Original EDGAR data")


fig_file <- file.path(figsDir, "figS8a.pdf")
ggsave(fig_file, p, width = 12, height = 5)

# Fig 8 (B)
df<- fread("/Users/savannah/Documents/Ch4/Data/Ag_opportunity_cost/__dt_Ag_CH4_N2O_xy_pastureV2.csv")

p <- ggplot(data=df) + 
  geom_tile(aes(x=x, y=y, fill=Ag_net_CH4_N2O)) +  
  scale_fill_gradientn(
    colors = c("white", "khaki1", "gold", "orange", "darkred", "purple"),
    values = scales::rescale(c(0, 1, 2, 4, 20)),
    limits = c(0, 20), na.value = "gray88",
    name = expression("CO"[2]*"e soil flux (Mg" ~ ha^-1 ~ yr^-1*")")) +
  theme_minimal() + coord_equal() + theme(panel.grid.major = element_blank()) +
  xlab("Longitude") +  ylab("Latitude") + 
  labs(title = expression("Net CH"[4]*"-N"[2]*"O effect from agriculture"),
       subtitle="Pasture-adjusted data")

fig_file <- file.path(figsDir, "figS8b.pdf")
ggsave(fig_file, p, width = 12, height = 5)


#####     Fig S3 - vMarch2025: enhanced plotting    ####   

#' Figure S3. Mean and SD of the original, untransformed greenhouse gas soil flux observations in Mg CO2e ha-1 year-1 for 
#' (A, C) CH4 and (B, D) N2O fluxes. 

b<-fread("/Users/savannah/Documents/Ch4/Data/Data_extraction/__Agg_Nov2024_noPeatland.csv") 

# Re-do Boreal flux obs. vNov2024: Going to assume growing season is half the year not 1/3. 
#vMar2024 already multiplied by 0.33 in 2_prep_stat_mod.R. 
b$CH4.C.converted[b$Biome == "Boreal forest"] <- b$CH4.C.converted[b$Biome == "Boreal forest"]*0.5/0.33
b$N2O.N.converted[b$Biome == "Boreal forest"] <- b$N2O.N.converted[b$Biome == "Boreal forest"]*0.5/0.33

b$CH4.CO2e_Mg_ha<- b$CH4.C.converted*85*1e4*1e-9 # mg CH4-C m-2 y-1 * 85 CO2e/CH4 on 20yr time horizon * 1e4 m2 ha-1 * 1e-9 Mg mg-1
b$N2O.CO2e_Mg_ha<- b$N2O.N.converted*273*1e4*1e-9 # mg N2O-N m-2 y-1 * 273 CO2e/N2O on 20yr time horizon * 1e4 m2 ha-1* 1e-9 Mg mg-1

b <-subset(b, Citation != 174)

nn<- data.frame(b$CH4.CO2e_Mg_ha, b$N2O.CO2e_Mg_ha, b$Biome)
nn$b.Biome<-as.factor(nn$b.Biome); levels(nn$b.Biome)


df<- subset(dt, select=c("x", "y", "Biome"))

library(gridExtra)
library(grid)
library(scales)
library(cowplot)

# Step 1: Calcuye biome-specific statistics from nn data
biome_stats <- nn %>%
  group_by(b.Biome) %>%
  summarize(
    CH4_CO2e_Mg_ha_biomeMean = mean(b.CH4.CO2e_Mg_ha, na.rm = TRUE),
    CH4_CO2e_Mg_ha_biomeSD = sd(b.CH4.CO2e_Mg_ha, na.rm = TRUE),
    N2O_CO2e_Mg_ha_biomeMean = mean(b.N2O.CO2e_Mg_ha, na.rm = TRUE),
    N2O_CO2e_Mg_ha_biomeSD = sd(b.N2O.CO2e_Mg_ha, na.rm = TRUE)
  )

# Step 2: Join these statistics to the spatial df based on Biome
df <- df %>%
  left_join(biome_stats, by = c("Biome" = "b.Biome"))

# Step 3: Create plots for CH4 and N2O means and SDs
# Mean CH4 plot
p_ch4 <- ggplot(data = df) + 
  geom_tile(aes(x = x, y = y, fill = CH4_CO2e_Mg_ha_biomeMean)) +  
  scale_fill_gradientn(
    colors = c("lightgreen","white", "pink", "pink2", "pink3", "deeppink4"),
    values = scales::rescale(c(-0.2, 0,1, 2.5)),
    limits = c(-0.2, 2.5), na.value = "gray88",
    name = expression("CO"[2]*"e soil flux (Mg" ~ ha^-1 ~ yr^-1*")")) +
  theme_minimal() + 
  coord_equal() + 
  theme(panel.grid.major = element_blank()) +
  xlab("Longitude") + 
  ylab("Latitude") + 
  labs(title = expression("Mean CH"[4]*"-CO"[2]*"e" ~ "soil flux"))

# SD CH4 plot
p_ch4_SD <- ggplot(data = df) + 
  geom_tile(aes(x = x, y = y, fill = CH4_CO2e_Mg_ha_biomeSD)) +  
  scale_fill_gradientn(
    colors = c("white", "lightblue2", "deepskyblue", "deepskyblue2", "deepskyblue4"),
    values = scales::rescale(c(0, 2, 5)),
    limits = c(0, 5), na.value = "gray88",
    name = expression("CO"[2]*"e soil flux (Mg" ~ ha^-1 ~ yr^-1*")")) +
  theme_minimal() + 
  coord_equal() + 
  theme(panel.grid.major = element_blank()) +
  xlab("Longitude") + 
  ylab("Latitude") + 
  labs(title = expression("SD CH"[4]*"-CO"[2]*"e" ~ "soil flux"))

# Mean N2O plot
p_n2o <- ggplot(data = df) + 
  geom_tile(aes(x = x, y = y, fill = N2O_CO2e_Mg_ha_biomeMean)) +  
  scale_fill_gradientn(
    colors = c("white", "pink", "pink2", "pink3", "deeppink4"),
    values = scales::rescale(c(0,0.6, 1.2)),
    limits = c(0, 1.2), na.value = "gray88",
    name = expression("CO"[2]*"e soil flux (Mg" ~ ha^-1 ~ yr^-1*")")) +
  theme_minimal() + 
  coord_equal() + 
  theme(panel.grid.major = element_blank()) +
  xlab("Longitude") + 
  ylab("Latitude") + 
  labs(title = expression("Mean N"[2]*"O-CO"[2]*"e" ~ "soil flux"))

# SD N2O plot
p_n2o_SD <- ggplot(data = df) + 
  geom_tile(aes(x = x, y = y, fill = N2O_CO2e_Mg_ha_biomeSD)) +  
  scale_fill_gradientn(
    colors = c("white", "lightblue2", "deepskyblue", "deepskyblue2", "deepskyblue4"),
    values = scales::rescale(c(0, 1, 4)),
    limits = c(0, 4), na.value = "gray88",
    name = expression("CO"[2]*"e soil flux (Mg" ~ ha^-1 ~ yr^-1*")")) +
  theme_minimal() + 
  coord_equal() + 
  theme(panel.grid.major = element_blank()) +
  xlab("Longitude") + 
  ylab("Latitude") + 
  labs(title = expression("SD N"[2]*"O-CO"[2]*"e" ~ "soil flux"))

# Create a common legend for all plots
legend <- get_legend(
  p_ch4 + theme(legend.position = "bottom") +
    guides(fill = guide_colorbar(
      title.position = "top",
      barwidth = 10, barheight = 0.5,
      title.hjust = 0.5
    ))
)

# Combine all plots
final_plot <- grid.arrange( arrangeGrob(p_ch4, p_n2o, p_ch4_SD, p_n2o_SD, ncol = 2),
  legend,heights = c(10, 1),nrow = 2)

# Save the figure with high resolution
# weird white space issues: ggsave("figS3_Global_CH4_N2O_Fluxes.png", final_plot, width = 12, height = 8, dpi = 300)

ggsave("figS3_CH4_mean.png", p_ch4, width = 6, height = 5, dpi = 300)
ggsave("figS3_N2O_mean.png", p_n2o, width = 6, height = 5, dpi = 300)
ggsave("figS3_CH4_SD.png", p_ch4_SD, width = 6, height = 5, dpi = 300)
ggsave("figS3_N2O_SD.png", p_n2o_SD, width = 6, height = 5, dpi = 300)
 

