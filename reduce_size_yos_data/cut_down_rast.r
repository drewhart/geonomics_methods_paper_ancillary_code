library(raster)

files = list.files('.')
files = files[grep('.tif$', files)]
print(files)

for (f in files){
    rast = raster(f)
    new = rast[68:nrow(rast), 68:ncol(rast), drop=F]
    new_fname = gsub('.tif', '_90x90.tif', f)
    print(paste('Now writing:', new_fname))
    writeRaster(new, new_fname, 'GTiff')
}
