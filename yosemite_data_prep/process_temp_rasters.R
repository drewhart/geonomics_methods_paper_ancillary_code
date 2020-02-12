#load packages
library(raster)
library(RColorBrewer)
library(rgdal)
library(jsonlite)
library(dismo)
library(rgeos)
library(maptools)


##################################################################################
# STEP 1: Get series of temp-change rasters, at high spatial res, for study region
##################################################################################

#directory containing the raster files
main_dir = '/home/drew/Desktop/stuff/berk/research/projects/sim/methods_paper/yosemite_data_prep/'
PRISM_data_dir = paste0(main_dir, 'PRISM_tmean_30yr_normal_800M2_annual_bil/')
PRISM_ppt_dir = paste0(main_dir, 'PRISM_ppt_30yr_normal_800M2_annual_bil/')
PRISM_filename = 'PRISM_tmean_30yr_normal_800mM2_annual_bil.bil'
PRISM_ppt_filename = 'PRISM_ppt_30yr_normal_800mM2_annual_bil.bil'

#read in the full raster
prism = raster(paste0(PRISM_data_dir, PRISM_filename))
prism_ppt = raster(paste0(PRISM_ppt_dir, PRISM_ppt_filename))


#set bounding box for simulated region
bbox_x = c(-120.3864, -119.0178)
bbox_y = c(37.2458, 38.6109)
bbox_to_crop = matrix(c(bbox_x, bbox_y), ncol = 2, byrow = T)
bbox_to_crop = extent(bbox_to_crop)
#crop it
prism_crop = crop(prism, bbox_to_crop)
prism_ppt_crop = crop(prism_ppt, bbox_to_crop)
#plot(prism_crop)
#prism_crop
#ncol(prism_crop) == nrow(prism_crop)

#reproject to the projection that all my coarser future data will be in
prism_crop_WGS84 = projectRaster(prism_crop,
                   projectExtent(prism_crop,
                   CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")))
prism_ppt_crop_WGS84 = projectRaster(prism_ppt_crop,
                       projectExtent(prism_ppt_crop,
                       CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")))


#set name of variable in the layer, so that the SDM we're projecting onto
#the raster further down finds the right variable name
prism_crop_WGS84@data@names = 'tmp'
prism_ppt_crop_WGS84@data@names = 'ppt'

#write raster to file
writeRaster(prism_crop_WGS84, 
            'PRISM_tmean_30yr_normal_800M_annual_bil_YOSEMITE_CROP.tif',
            overwrite = T)
writeRaster(prism_ppt_crop_WGS84, 
            'PRISM_ppt_30yr_normal_800M_annual_bil_YOSEMITE_CROP.tif',
            overwrite = T)


#read in annual-mean temperature projections for the region,
#(starting from 2014, because the PRISM normals cover 1980 to 2010)
#overlap with this PRISM raster,
#calculate difference between each cell and this raster's mean within that cell,
#then create new raster for that year
#by adding the fixed value to all PRISM cells within that cell
max_data_dir = paste0(main_dir, 'Cal-Adapt_yosemite_region_tasmax/')
min_data_dir = paste0(main_dir, 'Cal-Adapt_yosemite_region_tasmin/')
max_ppt_data_dir = paste0(main_dir, 'Cal-Adapt_yosemite_region_pptmax/')
min_ppt_data_dir = paste0(main_dir, 'Cal-Adapt_yosemite_region_pptmin/')

max_files = list.files(max_data_dir)
min_files = list.files(min_data_dir)
max_ppt_files = list.files(max_ppt_data_dir)
min_ppt_files = list.files(min_ppt_data_dir)

yrs = seq(2010, 2100, 5)
prism_changes = stack()
prism_ppt_changes = stack()
#make a copy of the prism_crop_WGS84 data, to which I will add yr-on-yr
#changes within each of the Cal-Adapt cells
changed_prism_rast = prism_crop_WGS84
changed_prism_ppt_rast = prism_ppt_crop_WGS84
#create a variable to hold the previous year's CalAdapt rast
prev_caladapt = raster()
prev_caladapt_ppt = raster()
tmp_stack = stack()
ppt_stack = stack()
for (i in seq(length(yrs))){
        yr = yrs[i]
        yr_max = raster(paste0(max_data_dir, max_files[grep(paste0('rcp45-', as.character(yr), '.tif'), max_files)]))
        yr_min = raster(paste0(min_data_dir, min_files[grep(paste0('rcp45-', as.character(yr), '.tif'), min_files)]))
        yr_avg = mean(stack(c(yr_max, yr_min)))
        yr_max_ppt = raster(paste0(max_ppt_data_dir, max_ppt_files[grep(paste0('rcp45-', as.character(yr), '.tif'), max_ppt_files)]))
        yr_min_ppt = raster(paste0(min_ppt_data_dir, min_ppt_files[grep(paste0('rcp45-', as.character(yr), '.tif'), min_ppt_files)]))
        #convert ppt from kg/m^2/s to mm/yr, to reconcile with other datasets
        yr_max_ppt = yr_max_ppt * 60 * 60 * 24 * 365.25
        yr_min_ppt = yr_min_ppt * 60 * 60 * 24 * 365.25
        yr_avg_ppt = mean(stack(c(yr_max_ppt, yr_min_ppt)))

        print(paste('Now processing year', yr))
        #print(yr_avg)
        yr_avg_crop = crop(yr_avg, bbox(prism_crop_WGS84))
        yr_avg_ppt_crop = crop(yr_avg_ppt, bbox(prism_crop_WGS84))
        tmp_stack = stack(tmp_stack, yr_avg_crop)
        ppt_stack = stack(ppt_stack, yr_avg_ppt_crop)

        #just save the raster over the prev_prism_rast variable, if the yr is 2010
        #(i.e. if it's the starting year)
        if (yr == 2010){
               prev_caladapt = yr_avg_crop
               prev_caladapt_ppt = yr_avg_ppt_crop
        } 
        else {
                #get change between previous CalAdapt raster and this one
                change_caladapt = yr_avg_crop - prev_caladapt
                change_caladapt_ppt = yr_avg_ppt_crop - prev_caladapt_ppt
                #print(change_caladapt@data@min)
                #print(change_caladapt@data@max)
                print(1)
                
                #get extent and res for this raster
                ext = bbox(yr_avg_crop)
                ext_ppt = bbox(yr_avg_ppt_crop)
                res_x = (ext[1,2] - ext[1,1]) / ncol(yr_avg_crop)
                res_x_ppt = (ext[1,2] - ext[1,1]) / ncol(yr_avg_ppt_crop)
                res_y = (ext[2,2] - ext[2,1]) / nrow(yr_avg_crop)
                res_y_ppt = (ext[2,2] - ext[2,1]) / nrow(yr_avg_ppt_crop)
                print(2)

                #use extent and res to get seqs of ll and ur corner x and y coords for
                #all cells in this raster
                xlls = seq(ext[1,1], ext[1,2]-res_x, res_x)
                xurs = seq(ext[1,1]+res_x, ext[1,2], res_x)
                ylls = seq(ext[2,1], ext[2,2]-res_y, res_y)
                yurs = seq(ext[2,1]+res_y, ext[2,2], res_y)
                xlls_ppt = seq(ext_ppt[1,1], ext_ppt[1,2]-res_x_ppt, res_x_ppt)
                xurs_ppt = seq(ext_ppt[1,1]+res_x_ppt, ext_ppt[1,2], res_x_ppt)
                ylls_ppt = seq(ext_ppt[2,1], ext_ppt[2,2]-res_y_ppt, res_y_ppt)
                yurs_ppt = seq(ext_ppt[2,1]+res_y_ppt, ext_ppt[2,2], res_y_ppt)
                print(3)
 
                #for each cell, add the diff between the previous CalAdapt raster
                #and this raster for that cell to all cells it overlaps in the PRISM data
                tile_num = 1
                for (col_num in 1:ncol(yr_avg_crop)){
                        for (row_num in 1:nrow(yr_avg_crop)){
                                #print(c(row_num, col_num))
                                this_ext = extent(matrix(c(xlls[col_num], xurs[col_num],
                                                    ylls[row_num], yurs[row_num]),
                                                  ncol = 2, byrow = T))
                                #print(this_ext)
                                cropped_rast = crop(changed_prism_rast, this_ext)
                                cropped_rast = cropped_rast + change_caladapt[row_num, col_num]
                                #if this is the first tile in the mosaic, just save it
                                if (tile_num == 1){
                                        mosaic_rast = cropped_rast
                                }
                                #if it's not the first item in the mosaic, then mosaic it
                                #into the existing mosaic
                                else {
                                        mosaic_rast = mosaic(mosaic_rast, cropped_rast,
                                                             fun = mean)
                                        #print(c(nrow(mosaic_rast), ncol(mosaic_rast)))
                                }
                                
                       }
                tile_num = tile_num + 1 
                }
                changed_prism_rast = mosaic_rast
                print(4)

                #do the same thing for ppt rasters
                tile_num = 1
                for (col_num in 1:ncol(yr_avg_ppt_crop)){
                        for (row_num in 1:nrow(yr_avg_ppt_crop)){
                                print(c(col_num, row_num))
                                #print(c(row_num, col_num))
                                print(4.1)
                                this_ext = extent(matrix(c(xlls_ppt[col_num],
                                                           xurs_ppt[col_num],
                                                           ylls_ppt[row_num],
                                                           yurs_ppt[row_num]),
                                                           ncol = 2, byrow = T))
                                print(4.2)
                                #print(this_ext)
                                cropped_rast = crop(changed_prism_ppt_rast, this_ext)
                                print(4.3)
                                cropped_rast = cropped_rast + change_caladapt_ppt[row_num, col_num]
                                print(4.4)
                                #if this is the first tile in the mosaic, just save it
                                if (tile_num == 1){
                                        mosaic_rast = cropped_rast
                                }
                                #if it's not the first item in the mosaic, then mosaic it
                                #into the existing mosaic
                                else {
                                        mosaic_rast = mosaic(mosaic_rast, cropped_rast,
                                                             fun = mean)
                                        #print(c(nrow(mosaic_rast), ncol(mosaic_rast)))
                                }
                                print(4.5)
                                
                       }
                tile_num = tile_num + 1 
                }
                changed_prism_ppt_rast = mosaic_rast
                print(5)


                #save the changed PRISM raster to the prism_changes stack
                prism_changes = stack(prism_changes, changed_prism_rast)
                prism_ppt_changes = stack(prism_ppt_changes, changed_prism_ppt_rast)

                #replace the old prev_caldapt raster with this one
                prev_caladapt = yr_avg_crop
                prev_caladapt_ppt = yr_avg_ppt_crop
                print(6)
        }
}

#plot start and end rasters, diff between them, and hist of diffs between them
#breaks = seq(-2,25)
#breaks2 = seq(0, 3.5, 0.25)
#tmp_pal = colorRampPalette(c("#2E95EA", "#F73232"))
#tmp_chng_pal = colorRampPalette(c("green", "red"))
#SDM_pal = colorRampPalette(c("red", "green"))
#par(mfrow = c(2,2))
#plot(prism_crop_WGS84, breaks = breaks, col = tmp_pal(length(breaks)))
#plot(prism_changes[[nlayers(prism_changes)]], breaks = breaks, col = tmp_pal(length(breaks)))
#plot(prism_changes[[nlayers(prism_changes)]] - prism_crop_WGS84,
#     breaks = breaks2,
#     col = tmp_chng_pal(length(breaks)))
#hist(prism_changes[[nlayers(prism_changes)]] - prism_crop_WGS84)

#NOTE: There's a strange column on the western edge of the prism_changes
#rasters that contains mostly but not entirely NAs.
#Not sure why, but really not important to figure it out because the
#area I'm choosing to simulate within is arbitrary anyhow.
#So, just subset everything (30-yr normals and their future projections)
#to the same, slightly smaller bbox, then collect them all as a single stack
#called tmp_lyrs
#final_prism = prism_crop_WGS84[15:nrow(prism_crop_WGS84), 15:ncol(prism_crop_WGS84), drop=F]
#final_prism_changes = prism_changes[14:nrow(prism_changes), 15:ncol(prism_changes), drop=F]
#tmp_lyrs = stack(final_prism, final_prism_changes)
tmp_lyrs = stack(crop(prism_crop_WGS84, prism_changes), prism_changes)

#final_prism_ppt = prism_ppt_crop_WGS84[15:nrow(prism_ppt_crop_WGS84), 15:ncol(prism_ppt_crop_WGS84), drop=F]
#ppt_lyrs = prism_ppt_changes[14:nrow(prism_ppt_changes), 15:ncol(prism_ppt_changes), drop=F]
ppt_lyrs = stack(crop(prism_ppt_crop_WGS84, prism_ppt_changes), prism_ppt_changes)
#final_prism_changes_ppt = prism_ppt_changes[14:nrow(prism_ppt_changes), 15:ncol(prism_ppt_changes), drop=F]
#ppt_lyrs = stack(final_prism_ppt, final_prism_changes_ppt)



#####################################
### STEP 2: Get a CA_NV shapefile 
#####################################
#usa = getData("GADM", country = "USA", level = 1)
#CA_NV = usa[usa$NAME_1 %in% c("California", "Nevada"),]
#CA_NV = aggregate(CA_NV, dissolve = T)
#coerce back to SPDF
#CA_NV = SpatialPolygonsDataFrame(CA_NV, data = data.frame(name = c('CA_NV')))

#transform to correct projection
#CA_NV_WGS84 = spTransform(CA_NV, CRS(proj4string(prism_crop_WGS84)))

#plot the shapefile, for demo
#plot(CA_NV)

#write to a shapefile
#writeOGR(CA_NV_WGS84, 'CA_NV_shapefile',
#         layer='CA_NV', driver = 'ESRI Shapefile',
#         overwrite_layer=T)

#***FOR OFFLINE USE***
CA_NV = readOGR(paste0(main_dir, 'CA_NV_shapefile'), 'CA_NV')


######################################
### STEP 3: Get bioclimate data for CA
######################################

#download worldclim data
#worldclim1 = getData('worldclim', var='bio', res=0.5, lon = -100, lat = 35)
#worldclim2 = getData('worldclim', var='bio', res=0.5, lon = -130, lat = 35)

#extract the 1st and 12th layers, which according to the Worldclim site are
#mean annual temp and ppt
#mat = mosaic(worldclim1[[1]], worldclim2[[1]], fun = 'mean')
#ppt = mosaic(worldclim1[[12]], worldclim2[[12]], fun = 'mean')

#reproject
#mat_WGS84 = projectRaster(mat, projectExtent(mat, CRS(proj4string((CA_NV)))))
#ppt_WGS84 = projectRaster(ppt, projectExtent(ppt, CRS(proj4string((CA_NV)))))

#write the mat and ppt layers to raster files
#raster::writeRaster(mat_WGS84, filename = 'bioclim_mat.tif',
#                    format = "GTiff", overwrite=T)
#raster::writeRaster(ppt_WGS84, filename = 'bioclim_ppt.tif',
#                     format = "GTiff", overwrite=T)

#***FOR OFFLINE USE***
#read in the mat and ppt layers from file
mat = raster(paste0(main_dir, 'bioclim_mat.tif'))
ppt = raster(paste0(main_dir, 'bioclim_ppt.tif'))

#crop the rasters to the CA shapefile
CA_NV_mat = crop(mat, CA_NV)
CA_NV_mat = CA_NV_mat/10 #NOTE: According to http://www.worldclim.org/formats1, the data ships with units of deg C * 10
names(CA_NV_mat) = "tmp"

#crop the rasters to the CA shapefile
CA_NV_ppt = crop(ppt, CA_NV)
CA_NV_ppt = CA_NV_ppt/10 #NOTE: According to http://www.worldclim.org/formats1, the data ships with units of deg C * 10
names(CA_NV_ppt) = "ppt"


#plot the raster
#plot(CA_mat, main = 'Mean Annual Temperature (deg Celsius)', xlab = 'lon', ylab = 'lat');
#plot(CA_NV, add = T)



#######################################################
### STEP 5: Download and process species data from GBIF
#######################################################

#download the data
#pres = gbif(genus = 'Sceloporus', species = 'graciosus', geo = T, download = T, ext = extent(CA))
#and write to file, for later offline use
#write.csv(pres, 'GBIF_Sceloporus_graciosus.csv')

#read in data from file
pres_csv = read.csv(paste0(main_dir, 'GBIF_Sceloporus_graciosus.csv'))

#make into a SPDF
pres_copy = pres_csv
coordinates(pres_copy) <- c('lon', 'lat')
proj4string(pres_copy) = CRS("+init=epsg:4326") 
pres_WGS84 = spTransform(pres_copy, CRS(proj4string(prism_crop_WGS84)))

#some code to remove problematic data (NAs, no geodetic datum)
pres_WGS84 <- pres_WGS84[CA_NV, ]

#subsample to reduce points falling within the same grid cells
pres <- gridSample(pres_WGS84, CA_NV_mat, n=1)


########################################
###STEP 6: Generate 'pseudoabsence' data
########################################

#generate pseudoabsence points
pseu = randomPoints(CA_NV_mat, 5*nrow(pres), pres)
pseu = as.data.frame(pseu)
coordinates(pseu) = c('x', 'y')
proj4string(pseu) = proj4string(CA_NV_mat)
pseu = pseu[CA_NV,]
pseu = as.data.frame(pseu)

#plot both presences and pseudoabsences over our other data
#plot(CA_mat[[1]], main = 'Presence and pseudoabsence data, plotted over MAT')
#plot(CA_NV, add = T)
#points(pres, col = 'blue', pch = 16, cex = 0.7)
#points(pseu, col = 'red', pch = 16, cex = 0.7)


####################################################
###STEP 7: Get the climate values at our data points
####################################################

#extract ppt and mat data to both presence and pseudoabsence points
mat.pres = extract(CA_NV_mat, pres[,c("lon", "lat")])
mat.pseu = extract(CA_NV_mat, pseu[,c("x", "y")])
ppt.pres = extract(CA_NV_ppt, pres[,c("lon", "lat")])
ppt.pseu = extract(CA_NV_ppt, pseu[,c("x", "y")])

#multiply ppt by 12, because it is expressed as monthly, not annual
ppt.pres = ppt.pres*12
ppt.pseu = ppt.pseu*12

#################################################
###STEP 8: Some final data-processing necessities
#################################################

#create a column for the dependent variable in our model (where 0 = absence, 1 = presence)
pres <- cbind(pres, presence=1)
colnames(pres) <- c("x","y","presence")
pseu <- cbind(pseu, presence=0)

#combine our presence and pseudoabsence data into a single data.frame
pres.abs.data = rbind(pres, pseu)
mat.data = c(mat.pres, mat.pseu)
ppt.data = c(ppt.pres, ppt.pseu)
reg.data = data.frame(pres.abs.data, tmp = mat.data, ppt=ppt.data)


######################################################
###STEP 9: Create our Species Distribution Model (SDM)
######################################################

#use a GLM to create the SDM
form = formula(presence ~ tmp + ppt)
mod = glm(form, family = binomial(link = 'logit'), data = reg.data)
summary(mod)

########################################################################################
###STEP 10: Project our model onto geographic space, under current and future conditions
########################################################################################

#use the model to project the species' range
#for each raster in our series of temperature-change rasters
SDM_lyrs = stack()
for (lyr_num in 1:nlayers(tmp_lyrs)){
        #get the layer
        tmp_lyr = tmp_lyrs[[lyr_num]] 
        ppt_lyr = ppt_lyrs[[lyr_num]] 
        #set the correct variable name for the SDM prediction further down
        tmp_lyr@data@names = 'tmp'
        ppt_lyr@data@names = 'ppt'
        lyrs = stack(tmp_lyr, ppt_lyr)
        fut_SDM = predict(lyrs, mod, type = 'response')
        SDM_lyrs = stack(SDM_lyrs, fut_SDM)
}
#threshold so that all cells below 0.3 get zeros
#SDM_lyrs[SDM_lyrs < 0.1] = 0

#par(mfrow = c(1,2))
#plot(SDM_lyrs[[1]],
#     breaks = seq(0,1,0.05),
#     col = SDM_pal(length(seq(0,1,0.05))))
#plot(CA_NV, add = T)
#plot(SDM_lyrs[[nlayers(SDM_lyrs)]],
#     breaks = seq(0,1,0.05),
#     col = SDM_pal(length(seq(0,1,0.05))))
#plot(CA_NV, add = T)


###############################################################################
###STEP 10.5: Crop layers to only the section where I don't have missing values
###############################################################################

tmp_lyrs = tmp_lyrs[ , 6:ncol(tmp_lyrs), drop = F]
ppt_lyrs = ppt_lyrs[ , 6:ncol(ppt_lyrs), drop = F]
SDM_lyrs = SDM_lyrs[ , 6:ncol(SDM_lyrs), drop = F]


##################################################################
### STEP 11: Write both series of rasters to their own directories
##################################################################
gnx_rast_dir = paste0(main_dir, 'yosemite_lyrs/')
tmp_dir = paste0(gnx_rast_dir, 'tmp/')
ppt_dir = paste0(gnx_rast_dir, 'ppt/')
sdm_dir = paste0(gnx_rast_dir, 'sdm/')
for (lyr_num in 1:nlayers(SDM_lyrs)){
        tmp_rast = tmp_lyrs[[lyr_num]]
        tmp_rast@extent = round(extent(tmp_rast), 4) 
        ppt_rast = ppt_lyrs[[lyr_num]]
        ppt_rast@extent = round(extent(ppt_rast), 4) 
        sdm_rast = SDM_lyrs[[lyr_num]]
        sdm_rast@extent = round(extent(sdm_rast), 4) 
        yr = as.character(yrs[lyr_num])
        #the first layer will be written to the geonomics_rast_dir,
        #because it will be used as one of the starting layers in the landscape
        #rather than as part of the change series
        if (yr  == '2010'){
                yr = '1980-2010'
                writeRaster(tmp_rast,
                            paste0(gnx_rast_dir, 'tmp_',
                                   yr, '.tif'),
                            overwrite = T)
                writeRaster(ppt_rast,
                            paste0(gnx_rast_dir, 'ppt_',
                                   yr, '.tif'),
                            overwrite = T)
                writeRaster(sdm_rast,
                            paste0(gnx_rast_dir, 'sdm_',
                                   yr, '.tif'),
                            overwrite = T)
        }
        else{
                #set the timestep in the simulation at which each layer
                #should be switched out and used (because this is how
                #geonomics requires that the directory of files be named)
                gnx_tstep = as.character(500 - 1 + 5*lyr_num)
                writeRaster(tmp_rast,
                            paste0(tmp_dir, gnx_tstep, '_tmp_',
                                   yr, '.tif'),
                            overwrite = T)
                writeRaster(ppt_rast,
                            paste0(ppt_dir, gnx_tstep, '_ppt_',
                                   yr, '.tif'),
                            overwrite = T)
                writeRaster(sdm_rast,
                            paste0(sdm_dir, gnx_tstep, '_sdm_',
                                   yr, '.tif'),
                            overwrite = T)
        }
}
print('All done!')

