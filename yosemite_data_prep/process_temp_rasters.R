#load packages
library(raster)

#directory containing the raster files
PRISM_data_dir = './PRISM_tmean_30yr_normal_800M2_annual_bil/'
PRISM_filename = 'PRISM_tmean_30yr_normal_800mM2_annual_bil.bil'

#read in the full raster
prism = raster(paste0(PRISM_data_dir, PRISM_filename))


#set bounding box for simulated region
bbox_x = c(-120.3864, -119.0178)
bbox_y = c(37.2458, 38.6109)
bbox_to_crop = matrix(c(bbox_x, bbox_y), ncol = 2, byrow = T)
bbox_to_crop = extent(bbox_to_crop)
bbox_to_crop
#crop
crop_prism = crop(prism, bbox_to_crop)
#plot(crop_prism)
#crop_prism
#ncol(crop_prism) == nrow(crop_data)

#reproject to the projection that all my coarser future data will be in
crop_prism_WGS84 = projectRaster(crop_data,
                   projectExtent(crop_data,
                   CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")))

#write raster to file
writeRaster(crop_prism_WGS84, 
            'PRISM_tmean_30yr_normal_800M_annual_bil_YOSEMITE_CROP.tif',
            overwrite = T)

#read in annual-mean temperature projections for the region,
#(starting from 2014, because the PRISM normals cover 1980 to 2010)
#overlap with this PRISM raster,
#calculate difference between each cell and this raster's mean within that cell,
#then create new raster for that year
#by adding the fixed value to all PRISM cells within that cell
max_data_dir = 'Cal-Adapt_yosemite_region_tasmax/'
min_data_dir = 'Cal-Adapt_yosemite_region_tasmin/'
max_files = list.files(max_data_dir)
min_files = list.files(min_data_dir)
yrs = seq(2010, 2100, 5)
for (i in seq(length(yrs))){
        yr = yrs[i]
        yr_max = raster(paste0(max_data_dir, max_files[grep(paste0('rcp45-', as.character(yr), '.tif'), max_files)]))
        yr_min = raster(paste0(min_data_dir, min_files[grep(paste0('rcp45-', as.character(yr), '.tif'), min_files)]))
        yr_avg = mean(stack(c(yr_max, yr_min)))
        print(yr)
        print(yr_avg)
        yr_avg_crop = crop(yr_avg, bbox(crop_prism_WGS84))
}


