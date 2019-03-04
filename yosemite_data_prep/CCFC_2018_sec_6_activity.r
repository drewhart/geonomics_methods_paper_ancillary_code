################################################
### STEP 1: Install and load necessary packages 
################################################

#Packages are just special files full of code designed for particular purposes
#These packages will be important forking with the data we need to use to make
#Species Distribution Models


#install packages if needed, then load
required.pkg = c("dismo","raster", "rgeos", "maptools", "rgdal", "jsonlite")
pkgs.not.installed <- required.pkg[!sapply(required.pkg, function(p) require(p, character.only=T))]
if (length(pkgs.not.installed) > 0){
  install.packages(pkgs.not.installed, dependencies=TRUE)
}

library(rgdal)
library(jsonlite)
library(dismo)
library(raster)
library(rgeos)
library(maptools)



#####################################
### STEP 2: Get a CA state shapefile 
#####################################

#A shapefile is just a file containg a bunch of line segments, connected in 'connect-the-dots' fashion,
#that draw the boundaries of some object.
#In this case, they will draw the CA state border
usa = getData("GADM", country = "USA", level = 1)
CA = usa[usa$NAME_1 %in% c("California"),]

#plot the shapefile, for demo
plot(CA)

#write to a shapefile
maptools::writePolyShape(CA, 'CA')

#***FOR OFFLINE USE***
CA = readOGR('CA.shp')



######################################
### STEP 3: Get bioclimate data for CA
######################################

#These are climate variables of relevance to biological questions (e.g. mean annual
#temperature, precipitation of driest or wettest month, etc.) 
#This data comes from: http://www.worldclim.org/bioclim
#Worldclim data was created and is curated by a research group at UC Davis
#The data is in the form of rasters, which are just square grids (think of a chessboard) in which each square 
#(i.e. cell) in the grid is assigned a value for a certain variable (e.g. annual precipitation).
#Digital photographs are saved in precisely the same way, but in the case of geographic rasters, the grids
#come with some information about where they belong on the earth's surface, so that they can be plotted on 
#x-y axes of longitude and latitude respectively, and thus can be used as maps.

#NOTES: 
# 1.) mat = Mean Annual Temperature; ppt = Annual Precipitation
# 2.) the mosaic command just connects our two different regions of climate data along their edges, like puzzle pieces 
# 3.) the crop command just uses the CA border to cookie-cut our climate data, extracting only data within CA

#download worldclim data
worldclim1 = getData('worldclim', var='bio', res=0.5, lon = -100, lat = 35)
worldclim2 = getData('worldclim', var='bio', res=0.5, lon = -130, lat = 35)

#extract the 1st and 12th layers, which according to the Worldclim site are mat and ppt
mat = mosaic(worldclim1[[1]], worldclim2[[1]], fun = 'mean')
ppt = mosaic(worldclim1[[12]], worldclim2[[12]], fun = 'mean')

#write the mat and ppt layers to raster files
raster::writeRaster(mat, filename = 'mat.tif', format = "GTiff")
raster::writeRaster(ppt, filename = 'ppt.tif', format = "GTiff")

#***FOR OFFLINE USE***
#read in the mat and ppt layers from file
mat = raster::raster('mat.tif')
ppt = raster::raster('ppt.tif')

#crop the rasters to the CA shapefile
CA_mat = crop(mat, CA)
CA_mat = CA_mat/10 #NOTE: According to http://www.worldclim.org/formats1, the data ships with units of deg C * 10
names(CA_mat) = "MAT"
CA_ppt = crop(ppt, CA)
names(CA_ppt) = "PPT"

#plot the two rasters, for demo
par(mfrow = c(1,2))
plot(CA_mat, main = 'Mean Annual Temperature (deg Celsius)', xlab = 'lon', ylab = 'lat'); plot(CA, add = T)
plot(CA_ppt, main = 'Mean Annual Precipitation (mm)', xlab = 'lon', ylab = 'lat'); plot(CA, add = T)



##############################################
###STEP 4: Create future climate-scenario data
##############################################

#Ideally, we would download maps of projected future climate variables, as determined by ensembles of 
#climate change projection models such as Dr. Silver discussed in lecture. However, in the interest
#of simplicity here, we will create our own future-scenario data by simply:
#   1.) adding 4 degrees Celsius to our current MAT data
#   2.) multiplying our PPT data by both 0.8 and 1.2 
    #   (i.e. creating two possible future precipitation scenarios, one 20% wetter than the past, one 20% drier)
#Creating these two alternative scenarios is one way of accounting for uncertainty about how climate change
#will affect future climates in CA. Obviously, the future of climate change is uncertain in many ways (e.g. which
#RCP is most representative of how society will respond to climate change? exactly how much warming will a
#given RCP cause? how will that global mean warming affect the diversity of local climates in CA?). But
#the future of CA precipitation is a particularly open question, and this uncertainty is particularly
#important given the relevance of water to ecosystems, especially in ecosystems that experience extreme water
#limitation, such as occur in much of CA.

#create future layers
CA_fut_mat = CA_mat + 4
CA_fut_ppt_wet = CA_ppt * 1.2
CA_fut_ppt_dry = CA_ppt * 0.8



#######################################################
### STEP 5: Download and process species data from GBIF
#######################################################

#These data come from: https://www.gbif.org/
#GBIF is a repository combining data from numerous herbaria and museum collections, as well as from
#iNaturalist and other citizen science platforms.
#The 'geo = T' argument returns only data for specimens that are geolocated, i.e that have a record of the
#latitude and longitude at which they were collected. This way we can use this data to build a Species
#Distribution Model (i.e. a map of all the areas where we predict the species lives).
#We will model Ponderosa pine, the dominant tree species in Sierra Nevada mid-elevation mixed-conifer forest

#NOTES:
  #1. 'pres' is short for presences (presences are just points where the species is recorded as being present)
  #We need to correlate these data with our climate data in order to figure out the range of mean annual temperature
  #and annual precipitation values that allow Ponderosa pine to grow

  #2. Some of these data are citizen science data, collected by lay people on their phones using the
  #iNaturalist app: http://inaturalist.ca/ You can see this from the row of data displayed a few lines down.


#download the data
pres = gbif(genus = 'Pinus', species = 'ponderosa', geo = T, download = T, ext = extent(CA))
#and write to file, for later offline use
write.csv(pres, 'gbif_ponderosa.csv')

#***FOR OFFLINE USE***
#read in data from file
pres = read.csv('gbif_ponderosa.csv')

#print rows containing of iNat data
View(pres[seq(25),c('X', 'basisOfRecord', 'family', 'species', 'lon', 'lat', 'institutionCode', 'occurrenceRemarks')])

#some code to remove problematic data (NAs, no geodetic datum)
pres.all <- subset(pres, !is.na(lon) & !is.na(lat))
unique(pres.all$geodeticDatum)
pres.all <- subset(pres.all, !is.na(geodeticDatum))
coordinates(pres.all) <- c("lon","lat")
proj4string(pres.all) <- proj4string(CA)
pres.all <- pres.all[CA, ]

#subsample to reduce points falling within the same grid cells
pres <- gridSample(pres.all, CA_ppt, n=1)

#save backup, for use in step 12
backup_pres = pres


########################################
###STEP 6: Generate 'pseudoabsence' data
########################################

#'Pseudoabsences' are our best guesses at real absences. It is not crucial that you understand what this
#'means, but in case you're curious: Because we don't have good
#data from numerous locations telling us that Ponderosa pine definitely does NOT grow there, we will instead
#create random points all over our map, in places where we don't have presence data, and use that as a
#guesstimate of absence.
#This brings up yet another type of uncertainty that comes into play in our models: uncertainty about exactly
#where the species IS and where it ISN'T. We do not have perfect, exhaustive data for the former, and we don't
#have any true data for the latter. This does not mean that we cannot create a reasonable model. It simply
#means that there will be some 'noise' in our estimation of the TRUE relationship between climate and species 
#suitability that we are interested in. This is a caveat that we will keep in mind when considering and
#assessing our results. 

#generate pseudoabsence points
pseu = randomPoints(CA_ppt, 1000, pres)
pseu = as.data.frame(pseu)
coordinates(pseu) = c('x', 'y')
proj4string(pseu) = proj4string(CA)
pseu = pseu[CA,]
pseu = as.data.frame(pseu)

#save backup, for use in step 12
backup_pseu = pseu

#plot both presences and pseudoabsences over our other data
plot(CA_ppt[[1]], main = 'Presence and pseudoabsence data, plotted over PPT (mm)')
plot(CA, add = T)
points(pres, col = 'blue', pch = 16, cex = 0.7)
points(pseu, col = 'red', pch = 16, cex = 0.7)




####################################################
###STEP 7: Get the climate values at our data points
####################################################

#Now that we have our datasets prepared, we need to 'extract' our climate data to our presence
#and pseudoabsence points (i.e. for each point, we need to record the temperature and precipitation values on
#our map at that point's location).


#extract ppt and mat data to both presence and pseudoabsence points
ppt.pres = extract(CA_ppt, pres[,c("lon", "lat")])
ppt.pseu = extract(CA_ppt, pseu[,c("x", "y")])
mat.pres = extract(CA_mat, pres[,c("lon", "lat")])
mat.pseu = extract(CA_mat, pseu[,c("x", "y")])




#################################################
###STEP 8: Some final data-processing necessities
#################################################

#Nothing particularly interesting here. Just technical necessities.


#create a column for the dependent variable in our model (where 0 = absence, 1 = presence)
pres <- cbind(pres, PRESENCE=1)
colnames(pres) <- c("x","y","PRESENCE")
pseu <- cbind(pseu, PRESENCE=0)

#combine our presence and pseudoabsence data into a single data.frame
pres.abs.data = rbind(pres, pseu)
mat.data = c(mat.pres, mat.pseu)
ppt.data = c(ppt.pres, ppt.pseu)
reg.data = data.frame(pres.abs.data, MAT = mat.data, PPT = ppt.data)




######################################################
###STEP 9: Create our Species Distribution Model (SDM)
######################################################

#Now we will use a linear regression to relate our species' presence/absence to values of annual precipitation
#and mean temperature.
#A linear regression estimates a 'best fit line', or 'trend line', which is just a linear equation that
#estimates the relationship between our climate variables and a species' presence. Thus, we will estimate b
#and c in:
  #  PRESENCE = b*MAT + c*PPT

#NOTES:
  #We are actually using a special type of linear regression, called a Generalized Linear Model (GLM), just because
  #this allows us to model a dependent variable that is binary (because PRESENCE is 0 if the species is absent, 1 if present)


#use a GLM to create the SDM
form = formula(PRESENCE ~ MAT + PPT)
mod = glm(form, family = binomial(link = 'logit'), data = reg.data)




######################################################################################
###STEP 10: Project our model onto geographic space, under current and future climates
######################################################################################

#Now that we've built a statistical model relating presence/absence of Ponderosa pine to mean annual
#temperature and annual precipitation, we can use that model to predict the species' presence/absence at each
#raster cell, and thus produce estimated maps of both current distribution and possible distributions under
#our future climate scenarios.
#NOTES:
  #1. If you compare our current map, you'll see that it matches up fairly well to other range maps online,
  #such as here: https://en.wikipedia.org/wiki/File:Pinus_ponderosa_subspecies_range_map_1.png. This is
  #pretty neat, given that we used only two basic climatic variables to paint a very broad-stroke picture of
  #Ponderosa's range.


#use the model to project the species' range
  #NOTE: The mask() command isn't necessary, and takes a while, so skip if need be
curr_sdm = predict(raster::stack(CA_mat, CA_ppt), mod, type = 'response')
fut_wet_sdm = predict(raster::stack(CA_fut_mat, CA_fut_ppt_wet), mod, type = 'response')
fut_dry_sdm = predict(raster::stack(CA_fut_mat, CA_fut_ppt_dry), mod, type = 'response')

curr_sdm = mask(curr_sdm, CA)
fut_wet_sdm = mask(fut_wet_sdm, CA)
fut_dry_sdm = mask(fut_dry_sdm, CA)

#set threshold value and extent, to use in mapping the projections
threshold = 0.6
extent <- extent(c(-125, -114, 32, 43))

#create a PDF of the maps
pdf('CCFC_Sec6_Ponderosa_pine_SDMs.pdf')

#map the projections
op <- par(mfrow = c(2,2), oma = c(2,4,2,0) + 0.1, mar = c(2,0,2,1) + 0.1)
plot(curr_sdm>= threshold, main = 'Current', xlab = 'lon', ylab = 'lat', legend = F, ext = extent)
plot(CA, add = T)
plot(fut_wet_sdm >= threshold, main = 'Hotter, wetter', xlab = 'lon', ylab = 'lat', legend = F)
plot(CA, add = T)
plot(fut_dry_sdm >= threshold, main = 'Hotter, drier', xlab = 'lon', ylab = 'lat', legend = F)
plot(CA, add = T)
title('Species Distribution Models of Ponderosa pine,\nfor current climate and future scenarios', outer = T)
par(op)

#turn off device
dev.off()






#############################################################################
###OPTIONAL STEP 11: Plot climatic niche on current and future climate spaces
#############################################################################

#How does the range of temperature and precipitation values defining Ponderosa pine's climatic niche
#compare to the total range of these values in California?
#And how does that relationship change under our future scenarios?
#Compare and discuss the following plots:


#create a PDf to save the climatic niche graphs
pdf('CCFC_Sec6_Ponderosa_pine_climatic_niche_graphs.pdf')

#plot the niche graphs
n_bg_pts = 150000
op <- par(mfrow = c(2,2), oma = c(1,1,3,1) + 0.1, mar = c(4,4,2,1) + 0.1)
plot(sample(CA_mat@data@values, size = n_bg_pts, replace = F), sample(CA_ppt@data@values, size = n_bg_pts, replace = F), col = 'gray', main = 'Current', xlab = 'MAT', ylab = 'PPT', xlim = c(0,30), ylim = c(0,2500))
points(mat.pres, ppt.pres, col = 'blue')
plot(sample(CA_fut_mat@data@values, size = n_bg_pts, replace = F), sample(CA_fut_ppt_wet@data@values, size = n_bg_pts, replace = F), col = 'gray', main = 'Hotter, wetter', xlab = 'MAT', ylab = 'PPT', xlim = c(0,30), ylim = c(0,2500))
points(mat.pres, ppt.pres, col = 'blue')
plot(sample(CA_fut_mat@data@values, size = n_bg_pts, replace = F), sample(CA_fut_ppt_dry@data@values, size = n_bg_pts, replace = F), col = 'gray', main = 'Hotter, drier', xlab = 'MAT', ylab = 'PPT', xlim = c(0,30), ylim = c(0,2500))
points(mat.pres, ppt.pres, col = 'blue')
title('Ponderosa pine climatic niche (blue) on top of CA climate space (gray),\nfor current climate and future scenarios', outer = T)
par(op)

#turn off device
dev.off()




#################################################
#OPTIONAL STEP 12: COMPARE MAT & PPT TO MAT & CWD
#################################################

#read in cwd data
cwd = readRDS('./BCM2014_cwd1981_2010_wy_ave_HST.Rdata')

#reproject the raster
cwd = projectRaster(cwd, crs = crs(mat))

#and crop it
cwd = crop(cwd, CA)

#create future proxies
fut_cwd_dry = cwd + 500
fut_cwd_wet = cwd - 500

#get our backups
p = backup_pres
s = backup_pseu

#plot both presences and pseudoabsences over our other data
plot(cwd [[1]], main = 'Presence and pseudoabsence data, plotted over CWD (mm)')
plot(CA, add = T)
points(p, col = 'blue', pch = 16, cex = 0.7)
points(s, col = 'red', pch = 16, cex = 0.7)

#extract cwd data to both presence and pseudoabsence points
cwd.p = extract(cwd, p[,c("lon", "lat")])
cwd.s = extract(cwd, s[,c("x", "y")])

#create a column for the dependent variable in our model (where 0 = absence, 1 = presence)
p <- cbind(p, PRESENCE=1)
colnames(p) <- c("x","y","PRESENCE")
s <- cbind(s, PRESENCE=0)

#combine our presence and pseudoabsence data into a single data.frame
cwd.pres.abs.data = rbind(p, s)
cwd.data = c(cwd.p, cwd.s)
cwd.reg.data = data.frame(cwd.pres.abs.data, MAT = mat.data, CWD = cwd.data)

#use a GLM to create the SDM
form = formula(PRESENCE ~ MAT + CWD)
mod = glm(form, family = binomial(link = 'logit'), data = cwd.reg.data)

#some code to make the extent and size of the rasters overlap, just for producing the projections
extent(CA_mat) = extent(cwd)
extent(CA_fut_mat) = extent(fut_cwd_wet)
rs_mat = resample(CA_mat, cwd)
rs_fut_mat = resample(CA_fut_mat, cwd)
names(cwd) = c('CWD')
names(fut_cwd_wet) = c('CWD')
names(fut_cwd_dry) = c('CWD')

#use the model to project the species' range
  #NOTE: The mask() command isn't necessary, and takes a while, so skip if need be
cwd_curr_sdm = predict(raster::stack(CA_mat, cwd), mod, type = 'response')
cwd_fut_wet_sdm = predict(raster::stack(CA_fut_mat, fut_cwd_wet), mod, type = 'response')
cwd_fut_dry_sdm = predict(raster::stack(CA_fut_mat, fut_cwd_dry), mod, type = 'response')

cwd_curr_sdm = mask(cwd_curr_sdm, CA)
cwd_fut_wet_sdm = mask(cwd_fut_wet_sdm, CA)
cwd_fut_dry_sdm = mask(cwd_fut_dry_sdm, CA)

#set threshold value and extent, to use in mapping the projections
threshold = 0.6
extent <- extent(c(-125, -114, 32, 43))

#create a PDF of the maps
pdf('CCFC_Sec6_Ponderosa_pine_SDMs_CWD.pdf')

#map the projections
op <- par(mfrow = c(2,2), oma = c(2,4,2,0) + 0.1, mar = c(2,0,2,1) + 0.1)
plot(cwd_curr_sdm>= threshold, main = 'Current', xlab = 'lon', ylab = 'lat', legend = F, ext = extent)
plot(CA, add = T)
plot(cwd_fut_wet_sdm >= threshold, main = 'Hotter, wetter', xlab = 'lon', ylab = 'lat', legend = F)
plot(CA, add = T)
plot(cwd_fut_dry_sdm >= threshold, main = 'Hotter, drier', xlab = 'lon', ylab = 'lat', legend = F)
plot(CA, add = T)
title('Species Distribution Models of Ponderosa pine,\nfor current climate and future scenarios with CWD', outer = T)
par(op)

#turn off device
dev.off()

#create a PDf to save the climatic niche graphs
pdf('CCFC_Sec6_Ponderosa_pine_climatic_niche_graphs_CWD.pdf')

#plot the niche graphs
n_bg_pts = 150000
op <- par(mfrow = c(2,2), oma = c(1,1,3,1) + 0.1, mar = c(4,4,2,1) + 0.1)
plot(sample(CA_mat@data@values, size = n_bg_pts, replace = F), sample(cwd@data@values, size = n_bg_pts, replace = F), col = 'gray', main = 'Current', xlab = 'MAT', ylab = 'CWD', xlim = c(0,30), ylim = c(0,2500))
points(mat.pres, ppt.pres, col = 'blue')
plot(sample(CA_fut_mat@data@values, size = n_bg_pts, replace = F), sample(fut_cwd_wet@data@values, size = n_bg_pts, replace = F), col = 'gray', main = 'Hotter, wetter', xlab = 'MAT', ylab = 'CWD', xlim = c(0,30), ylim = c(0,2500))
points(mat.pres, ppt.pres, col = 'blue')
plot(sample(CA_fut_mat@data@values, size = n_bg_pts, replace = F), sample(fut_cwd_dry@data@values, size = n_bg_pts, replace = F), col = 'gray', main = 'Hotter, drier', xlab = 'MAT', ylab = 'CWD', xlim = c(0,30), ylim = c(0,2500))
points(mat.pres, ppt.pres, col = 'blue')
title('Ponderosa pine climatic niche (blue) on top of CA climate space (gray),\nfor current climate and future scenarios with CWD', outer = T)
par(op)

#turn off device
dev.off()



