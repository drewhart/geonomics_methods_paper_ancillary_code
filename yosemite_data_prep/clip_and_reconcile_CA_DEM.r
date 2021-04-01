library(sf)
library(tmap)
library(raster)

dem = raster('../geonomics/geonomics/data/yosemite_demo/CA_DEM.tif')
tmp = raster('../geonomics/geonomics/data/yosemite_demo/yosemite_lyrs/tmp_1980-2010.tif')

#set the bbox values, because this somehow got lost and screwed up when I hackily
#ripped this data out of the annoying-as-hell ArcLayer package that I downloaded it in
#NOTE: values pulled directly from the file's metadata XML file (which I saved in
#this directory)
dem@extent@xmin = -124.482003 
dem@extent@xmax = -114.131170
dem@extent@ymin = 32.528832
dem@extent@ymax = 42.009665

#CRS is already lined up, so just clip, then resample
dem_crop = crop(dem, tmp)
dem_crop_resamp = resample(dem_crop, tmp)

#check the resolutions are the same
res(dem_crop_resamp) == res(tmp)

#then just going to use a linear model to backfill the missing values that were clipped
#out of the original DEM layer because they fall in NV; this is fine because the DEM
#is not being used for analysis at all; rather, it's just being used as the terrain
#onto which the results will be draped
dem_vals = dem_crop_resamp[!is.na(dem_crop_resamp)]
tmp_vals = tmp[!is.na(dem_crop_resamp)]
df = data.frame(alt=dem_vals, tmp=tmp_vals)
mod = lm(alt ~ tmp, data=df)
preds = stats::predict.lm(mod, data.frame(tmp=tmp[is.na(dem_crop_resamp)]))
dem_crop_resamp[is.na(dem_crop_resamp)] = preds

#plot to check
tmap_mode('view')
tm_shape(dem_crop_resamp) + tm_raster(alpha = 0.5) + tm_shape(tmp) + tm_raster(alpha = 0.5)

#now write to new file
writeRaster(dem_crop_resamp, 'yosemite_DEM.tif', 'GTiff', overwrite=T)
