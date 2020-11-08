##############################################################
# TODO:
# 1. get the env rasters at start and end,
#calculate scaled 3d distance between each pixel before and after change,
# then run spatial regression of phenotype change on env change

# 2. go back and look at draft yosemite fig and work on other parts of it

# 3. generate gif of rotating 3d landscape with phenotype and phenotype change on it
# (figure out hillshading/shadow for that??)

# 4. do same thing for pop density! (still need to KDE it; see py script)

# 5. get drafts of all figs combined and plopped into paper to get Ian's feedback

library(raster)
library(rgdal)
library(rasterVis)
library(rayshader)
library(pals)
library(wesanderson)
library(MASS)
library(nlme)
library(ggplot2)
library(ggthemes)
library(cowplot)

# output controls
save_vid = F
save_figs = F


###################
###################
###### MAPS #######
###################
###################

# read in the phenotype rasters
rasterize_table = function(filename, dem){
    rast = raster(as.matrix(read.table(filename)))
    rast@extent = dem@extent
    rast@crs = dem@crs
    return(rast)
}

# function to normalize a matrix to values between low and high,
# optionally flooring them integers
normalize_low_to_high = function(mat, low, high, floor_it=F){
    minval = min(mat)
    maxval = max(mat)
    out = (mat-minval)/(maxval-minval)
    out = ((high-low) * out)+low
    if (floor_it){
        out = floor(out)
    }
    return(out)
}


# function to get an ixjx3 array of RGB values
# corresponding to the given palette's hex values
# for an ixj array of data
get_rgb_vals = function(arr, pal){
    rgbs = lapply(seq(1, length(pal)), 
                  function(i){return(col2rgb(pal[i]))})
    out_arr = array(0, dim=c(nrow(arr), ncol(arr),3))
    for (i in seq(dim(arr)[1])){
        for (j in seq(dim(arr)[2])){
            out_arr[i, j, ] = rgbs[[arr[i, j]]][,1]
        }
    }
    return(out_arr)
}


# function to make a RGB-contrast ixjx3 array
# from the ixj raster of data, using the given palette
make_rgb_contrast_array = function(rast, pal){
    # convert raster to a matrix
    mat = rayshader::raster_to_matrix(rast)
    # normalize to ints between 1 and palette length
    norm = normalize_low_to_high(mat, 1, length(pal), floor_it=T)
    # convert to a 3-layer array of the RGB values
    rgb_arr = get_rgb_vals(norm, pal)
    # permute to line up correctly for the plot_3D function (I believe?)
    rgb_arr = aperm(rgb_arr, c(2,1,3))
    # create a contrast array
    rgb_contrast = scales::rescale(rgb_arr, to=c(0,1))
    return(rgb_contrast)
}


# function to make a hillshade plot of the given raster and given DEM,
# colored by the given palette, using the given title
make_hillshade_plot = function(rast, dem, pal, title){
    slope = terrain(dem, opt='slope')
    aspect = terrain(dem, opt='aspect')
    hill = hillShade(slope, aspect, angle=65, direction=180)
    # plot the overlay for phenotype change
    plot(hill, main = title, col = grey(1:100/100), legend = FALSE)
    plot(rast, add = TRUE, alpha = .5, col=pal)
}


# function to make a 3d plot of the given raster and given DEM,
# colored by the given palette, using the given title
make_3d_plot = function(rast, dem, pal, title, var_name='no_var', zscale=10,
                        phi=30, theta_start=-45, fov=0,
                        background_color='#000000', title_bar_color='#ffffff',
                        title_color='#000000', snapshot=F, video=F){
    # get the dem as a matrix
    demmat = raster_to_matrix(dem)
    #use some of rayshader's built-in textures on the dem
    # NOTE: FIGURE OUT HOW TO ADD THIS TO GET HILLSHADING? CAN I EVEN DO IT?
    if (F){
    demmat_texture = demmat %>%
          sphere_shade(texture = "desert") %>%
          add_water(detect_water(demmat), color='desert') %>%
          add_shadow(ray_shade(demmat), 0.5) %>%
          add_shadow(ambient_shade(demmat), 0) 
    }
    # get the raster as an RGB-contrast 3-layer array
    rgb_contrast = make_rgb_contrast_array(rast, pal)
    # make the plot 
    plot_3d(rgb_contrast, demmat, windowsize = c(1100,900), zscale = zscale,
            shadowdepth = -50, zoom=0.5, phi=phi, theta=theta_start, fov=fov,
            background = background_color, shadowcolor = "#523E2B")
    if (snapshot){
        render_snapshot(title_text = title, title_bar_color = title_bar_color,
                        title_color = title_color, title_bar_alpha = 1)
    }
    if (video){
        # save 3d rotating video
        angles= seq(0,360,length.out = 1441)[-1]
        for(i in 1:(length(angles)-1)) {
          render_camera(theta=theta_start+angles[i], phi=phi)
          render_snapshot(filename = sprintf("./video_files/yosemite_vid_%s_%i.png", var_name, i),
                          title_text = title, title_bar_color = title_bar_color,
                          title_color = title_color, title_bar_alpha = 1) }
        rgl::rgl.close()
        system(paste0(sprintf("ffmpeg -framerate 60 -i ./video_files/yosemite_vid_%s", var_name),
                      "_%d.png -pix_fmt yuv420p ",
                      sprintf("./video_files/yosemite_%s.mp4", var_name)))
    }
} 

# function to normalize the values in the before and after rasters
# (so that differences calculated as 3d Euclidean distances weight
# each var equally), but to do so using the full range of var values
# across both the b4 and after layers (so that calculated differences
# still reflect the fact that changes in temp are largely increases,
# changes in sdm are largely decreases, ppt is mixed)
get_normed_stacks = function(b4, af){
    for (i in seq(nlayers(b4))){
        vals = c(b4[[i]][,], af[[i]][,])
        minval = min(vals)
        maxval = max(vals)
        b4[[i]] = (b4[[i]] - minval)/(maxval-minval)
        af[[i]] = (af[[i]] - minval)/(maxval-minval)
    }
    return(list(b4, af))
}


# function to make a hillshaded ggplot with the given titile
# and from the given input raster using the layer data
# (which will be labeled in the legend
# using the given display_name, and the colors for which will be 'cols' 
# will be scaled between minval and maxval, with the mean of those two
# vals being used as the contour-line value)
make_gg_hillshade_plot = function(rast, hillrast, display_name,
                                  minval, maxval, cols, title){
  # make both rasters into data.frames
  df = as.data.frame(rast, xy=T)
  hill_df = as.data.frame(hillrast, xy=T)
  colnames(df) = c('x', 'y', 'layer')
  # create the ggplot object
  out = ggplot() +
      # plot the main raster, filling by its value column
      geom_raster(data = df, 
                  aes(x = x, y = y, fill=layer)) +  
      # add the hillshade raster, filling by its value column
      geom_raster(data = hill_df,
                  aes(x = x, y = y, alpha=layer)) +
      # add the contour line at the midway value
      # between the min and max of this variable across time
      geom_contour(data = df, colour='#edede8',
                   breaks=c(minval, mean(c(minval, maxval)), 1.01*maxval),
                   aes(x=x, y=y, z=layer)) +
      # scale the fill using the given colors and display_name
      scale_fill_gradientn(colours=cols, name=display_name,
                           limits=c(minval, maxval*1.01)) +
      # scale the alpha vals of the hillshade layer
      scale_alpha(range = c(0.15, 0.5), guide = "none") +  
      # use the ggthemes map theme
      ggthemes::theme_map() + 
      # tweak theme aspects
      theme(legend.position = 'right',
            legend.background = element_rect(fill = "darkgray"), 
            legend.key = element_rect(fill = "lightblue", color = NA),
            legend.key.size = unit(1.5, "cm"), legend.key.width = unit(0.5,"cm"),
            plot.title = element_text(size = 12),
            #axis.title = element_text(),
            axis.text = element_text(), axis.ticks = element_line()) +
      #xlab("lon") +
      #ylab("lat") + 
      # add the title
      ggtitle(title) +
      coord_quickmap()
  return(out)
}


# function to create mirrored sets of individuals and pad them around the
# real individuals, to take care of edge effects
mirror_pad_inds = function(inds, extent){
    # get original and mirrored x and y coords
    xs = inds$x
    ys = inds$y
    mirr_xs = -xs
    mirr_ys = -ys
    # create a set of inds mirrored to the left
    lpad = inds 
    lpad$x = mirr_xs
    # create a set of inds mirrored to the right
    rpad = inds
    rpad$x = mirr_xs+(2*extent[[1]])
    # create a set of inds mirrored below 
    bpad = inds 
    bpad$y = mirr_ys
    # create a set of inds mirrored above 
    upad = inds
    upad$y = mirr_ys+(2*extent[[2]])
    # concatenate all into new data.frame
    out = rbind(inds, lpad, rpad, bpad, upad) 
    return(out)
}


####################
# Get and prep data:

# TODO: DELETE LO-RES DEM??
# read in the yosemite DEM
dem = raster('yosemite_DEM_90x90.tif')

# read in the before and after kriged phenotypes
b4_z = rasterize_table('b4_cc_krig.txt', dem)
af_z = rasterize_table('af_cc_krig.txt', dem)
# and get the phenotype-change raster
del_z = af_z - b4_z

# downscale everything (for better viz)
factor = 20
# TODO: DELETE LO-RES DEM??
dem = disaggregate(dem, factor, method='bilinear')
b4_z = disaggregate(b4_z, factor, method='bilinear')
af_z = disaggregate(af_z, factor, method='bilinear')
del_z = disaggregate(del_z, factor, method='bilinear')

# read in the hi-res DEM, mosaic, and clip (NOTE: doesn't need reprojection)
hrdem1 = raster('./hi_res_DEM/n37_w120_1arc_v3.tif')
hrdem2 = raster('./hi_res_DEM/n37_w121_1arc_v3.tif')
hrdem3 = raster('./hi_res_DEM/n38_w120_1arc_v3.tif')
hrdem4 = raster('./hi_res_DEM/n38_w121_1arc_v3.tif')
hrdem = mosaic(hrdem1, hrdem2, hrdem3, hrdem4, fun='mean')
hrdem = crop(hrdem, dem)

# get the 3 variables before climate change
b4 = stack(c('./yosemite_env_layers/ppt_1980-2010_90x90.tif', 
             './yosemite_env_layers/tmp_1980-2010_90x90.tif',
             './yosemite_env_layers/sdm_1980-2010_90x90.tif'))

# get the 3 vars after
af = stack(c('./yosemite_env_layers/594_ppt_2100_90x90.tif', 
             './yosemite_env_layers/594_tmp_2100_90x90.tif',
             './yosemite_env_layers/594_sdm_2100_90x90.tif'))

# norm them
normed = get_normed_stacks(b4, af)

# get the Euclidean-distance raster
# (Euclidean distance each pixel traveled in normed env space)
euc = sqrt(sum((normed[[2]] - normed[[1]])^2))


# read in the CSVs of individs' points
b4_inds = read.csv('./b4_cc_individs.csv')
af_inds = read.csv('./af_cc_individs.csv')


# produce 2D KDE plots of population density, as rasters
# NOTE: transpose the KDE matrices because, per the docs,
# x vals are on the rows (not sure why?)
# NOTE: add a margin to the limits, then subset the core of the KDE
# matrix, to attempt to ameliorate the edge effect
# NOTE: normalizing the raster by its own sum,
# then multiplying by total pop size, as I believe this
# gives me cells that can be interpreted as actual pop densities
h=8
real_width = nrow(b4)
margin_width=10
lims=c(-margin_width, real_width+margin_width,
       -margin_width, real_width+margin_width)

mirr_b4_inds = mirror_pad_inds(b4_inds, c(90,90)) 
mirr_af_inds = mirror_pad_inds(af_inds, c(90,90)) 

b4_kde = kde2d(mirr_b4_inds$x, mirr_b4_inds$y,
                        n=real_width + (2*margin_width),
                        h=c(h,h), lims=lims)
b4_kde = b4_kde[[3]][(margin_width+1):(margin_width+real_width),
                     (margin_width+1):(margin_width+real_width)]
b4_kde = raster(t(b4_kde))
b4_kde = nrow(b4_inds) * b4_kde/cellStats(b4_kde, 'sum')
b4_kde@extent = b4@extent
b4_kde@crs = b4@crs
af_kde = kde2d(mirr_af_inds$x, mirr_af_inds$y,
                        n=real_width + (2*margin_width),
                        h=c(h,h), lims=lims)
af_kde = af_kde[[3]][(margin_width+1):(margin_width+real_width),
                     (margin_width+1):(margin_width+real_width)]
af_kde = raster(t(af_kde))
af_kde = nrow(af_inds) * af_kde/cellStats(af_kde, 'sum')
af_kde@extent = af@extent
af_kde@crs = af@crs
  



#################
# Make the plots:

# get hillshade raster
slope = terrain(hrdem, opt='slope')
aspect = terrain(hrdem, opt='aspect')
hill = hillShade(slope, aspect, angle=15, direction=10); 


# plot temperature and phenotype hillshade plots,
# before and after climate change
n_breaks = 100
#cols = coolwarm(n_breaks)
cols = brewer.rdbu(n_breaks*1.5)[125:26]

b4_tmp_hillplot = make_gg_hillshade_plot(b4[[2]], hill, 'tmp (°C)',
                              b4[[2]]@data@min, af[[2]]@data@max,
                              cols, 'temperature before climate change')
af_tmp_hillplot = make_gg_hillshade_plot(af[[2]], hill, 'tmp (°C)',
                              b4[[2]]@data@min, af[[2]]@data@max,
                              cols, 'temperature after climate change')
b4_z_hillplot = make_gg_hillshade_plot(b4_z, hill, '  z',
                              b4_z@data@min, af_z@data@max,
                              cols, 'phenotypes before climate change')
af_z_hillplot = make_gg_hillshade_plot(af_z, hill, '  z',
                              b4_z@data@min, af_z@data@max,
                              cols, 'phenotypes after climate change')

pg_tmp_z = plot_grid(b4_tmp_hillplot, af_tmp_hillplot,
                     b4_z_hillplot, af_z_hillplot, ncol=2)

if (save_figs){
  ggsave(pg_tmp_z, file='b4_af_tmp_z_plot.pdf',
         width=35, height=20, units = "cm", dpi=500)
}


# plot hab suitability and pop density hillshade plots,
# before and after climate change
n_cols = 100
zissou = wes_palette("Zissou1", n_cols, type = "continuous")

min_sdm_val = min(b4[[3]]@data@min, af[[3]]@data@min)
max_sdm_val = max(b4[[3]]@data@max, af[[3]]@data@max)
min_kde_val = min(b4_kde@data@min, af_kde@data@min)
max_kde_val = max(b4_kde@data@max, af_kde@data@max)

b4_sdm_hillplot = make_gg_hillshade_plot(b4[[3]], hill, 'hab',
                              min_sdm_val, max_sdm_val,
                              zissou, 'habitat suitability before climate change')
af_sdm_hillplot = make_gg_hillshade_plot(af[[3]], hill, 'hab',
                              min_sdm_val, max_sdm_val,
                              zissou, 'habitat suitability after climate change')
b4_nt_hillplot = make_gg_hillshade_plot(b4_kde, hill, 'pop dens\n(ind/cell)',
                              min_kde_val, max_kde_val,
                              zissou, 'population density before climate change')
af_nt_hillplot = make_gg_hillshade_plot(af_kde, hill, 'pop dens\n(ind/cell)',
                              min_kde_val, max_kde_val,
                              zissou, 'population density after climate change')
pg_sdm_nt = plot_grid(b4_sdm_hillplot, af_sdm_hillplot,
                      b4_nt_hillplot, af_nt_hillplot, ncol=2)

if (save_figs){
  ggsave(pg_sdm_nt, file='b4_af_sdm_nt_plot.pdf',
         width=35, height=20, units = "cm", dpi=500)
}




# 3d plots
zscale=5
if (save_vid){
  make_3d_plot(b4_z, dem, coolwarm, 'Kriged phenotype before climate change',
               var_name='pheno', zscale=zscale, background_color='#a6a490',
               title_bar_color='#ebebeb', snapshot=F, video=T, phi=22, fov=70)
}

make_3d_plot(del_z, dem, zissou,
             'Difference in kriged phenotypes before and after climate change',
             zscale=zscale, background_color='#000000', title_bar_color='#99cde0')



###############################################################



###############################################################

# get a normalized difference raster, then plot it
diff_inds = (af_kde-b4_kde)/b4_kde
col = colorRampPalette(c("red", "white", "blue"))(255)
plot(diff_inds, col=col)

# investigate relationship between change in pop density and env change
diff_inds_vals = diff_inds[,]
euc_vals = euc[,]
coords = lapply(seq(length(euc_vals)), function(i){return(xyFromCell(diff_inds, i))})
coords = data.frame(matrix(unlist(coords), ncol=2, byrow=T))
colnames(coords) = c('lon', 'lat')
df = data.frame(inds = diff_inds_vals, euc = euc_vals,
                lon = coords$lon, lat = coords$lat, dummy=1)
mod = lme(fixed=inds ~ euc, data=df, random= ~1|dummy,
          correlation=corGaus(1, form = ~ lon + lat))
# plot scatterplot along with model trendline
plot(euc[,], diff_inds[,])



###################
###################
###### ETC. #######
###################
###################

# make the pop-growth plot
nt = read.csv('./Nt_EXTENDED.csv')
# set total burn-in time, then truncate dataset
burn_t = 98
nt$t = nt$t-98
nt = nt[nt$t>=0, ] 
# set ylims
ylims = c(min(nt$Nt)*0.9, max(nt$Nt) * 1.1)

nt_plot = ggplot() +
    geom_line(aes(t, Nt), col='#47b3ab', size=2, data=nt) +
    geom_vline(xintercept= 509, colour = '#fc033d') + 
    geom_vline(xintercept= 594, colour = '#fc033d') + 
    scale_x_continuous(name = '') +
    #scale_x_continuous(name = 't (time steps)') +
    scale_y_continuous(name = '', limits=ylims) + 
    #scale_y_continuous(name = 'Nt (individuals/cell)', limits=ylims) + 
    ggthemes::theme_gdocs() +
    #ggtitle("Population dynamics") +
    theme(axis.title = element_text(size=30), axis.text = element_text(size=25, color='black'))
nt_plot

ggsave(nt_plot, file='Nt_plot.pdf',
       width=40, height=20, units='cm', dpi=500)

# note that mean popsize af/mean popsize b4 = sum(sdm_af)/sum(sdm_b4)
b4sum = sum(b4[[3]][,])
afsum = sum(af[[3]][,])
sdm_ratio = afsum/b4sum
print('sdm ratio')
print(sdm_ratio)
b4meanNt = mean(nt[200:400,'Nt']) 
afmeanNt = mean(nt[600:650,'Nt']) 
Nt_ratio = afmeanNt/b4meanNt
print("Nt ratio")
print(Nt_ratio)


