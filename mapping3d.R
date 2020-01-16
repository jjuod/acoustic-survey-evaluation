library(rgdal)
library(ggplot2)
library(geoviz)
library(rayshader)
library(raster)

# redefine a function from geoviz, because it seems to use wrong zoom level
get_slippy_map_15 <- function(bounding_box, image_source = "stamen", image_type = "watercolor", max_tiles = 30, api_key){
  
  xt_scene <- raster::extent(bounding_box)
  print(xt_scene)
  
  overlay_bbox <-
    sf::st_bbox(c(xmin = xt_scene@xmin,
                  xmax = xt_scene@xmax,
                  ymin = xt_scene@ymin,
                  ymax = xt_scene@ymax),
                crs = sf::st_crs("+proj=longlat +datum=WGS84 +no_defs"))
  
  tile_grid <- slippymath::bbox_to_tile_grid(overlay_bbox, max_tiles = max_tiles)
  
  if(tile_grid$zoom > 15 & image_source == "mapbox" & image_type == "terrain-rgb"){
    message(glue::glue("Zoom level with max_tiles = {max_tiles} is {tile_grid$zoom}. Resetting zoom to 15, which is max for mapbox.terrain-rgb."))
    tile_grid <- slippymath::bbox_to_tile_grid(overlay_bbox, zoom = 15)
  }
  
  
  
  if(image_source=="stamen"){
    if(stringr::str_detect(image_type, "watercolor")){
      query_string <- paste0("http://tile.stamen.com/", image_type, "/{zoom}/{x}/{y}.jpg")
    } else {
      query_string <- paste0("http://tile.stamen.com/", image_type, "/{zoom}/{x}/{y}.png")
    }
    
  } else if (image_source=="mapbox"){
    
    query_string <- paste0("https://api.mapbox.com/v4/mapbox.", image_type, "/{zoom}/{x}/{y}.jpg90",
                           "?access_token=",
                           api_key)
  } else {
    stop(glue::glue("unknown source '{image_source}'"))
  }
  
  #create a temporary dir to hold tiles
  tile_dir <- tempfile(pattern = "map_tiles_")
  dir.create(tile_dir)
  
  images <-
    purrr::pmap(tile_grid$tiles,
                function(x, y, zoom){
                  outfile <- glue::glue("{tile_dir}/{x}_{y}.jpg")
                  curl::curl_download(url = glue::glue(query_string),
                                      destfile = outfile)
                  outfile
                },
                zoom = tile_grid$zoom)
  
  raster_out <- slippymath::compose_tile_grid(tile_grid, images)
  
  unlink(tile_dir, recursive = TRUE)  #kill the temp directory containing tiles
  
  return(raster_out)
}

# with fixed CRS and raise_agl6
add_gps_to_rayshader <- function(raster_input, lat, long, alt, zscale, line_width = 1, colour = "red", alpha = 0.8,
                                 lightsaber = TRUE, clamp_to_ground = FALSE, raise_agl = 0, ground_shadow = FALSE,
                                 as_line = TRUE, point_size = 20, pch=16){
  
  coords <- latlong_to_rayshader_coords(raster_input, lat, long)
  
  distances_x <- coords$x
  
  distances_y <- coords$y
  
  
  if (clamp_to_ground | ground_shadow) {
    
    sp_gps <- sp::SpatialPoints(cbind(long, lat), proj4string = sp::CRS('+init=epsg:4326'))
    
    sp_gps <- sp::spTransform(sp_gps, sp::CRS(as.character(raster::crs(raster_input))))
    
    gps_ground_line <- raster::extract(raster_input, sp_gps)
    
  }
  
  if(clamp_to_ground){
    
    track_altitude <- gps_ground_line
    
  } else {
    
    track_altitude <- alt
  }
  
  if(as_line){
    
    if(!lightsaber){
      rgl::lines3d(
        distances_x,  #lat
        track_altitude / zscale,  #alt
        -distances_y,  #long
        color = colour,
        alpha = alpha,
        lwd = line_width
      )
    } else {
      
      #render track 3 times with transparent & thicker outside
      
      rgl::lines3d(
        distances_x,
        track_altitude / zscale,
        -distances_y,
        color = colour,
        alpha = 0.2,
        lwd = line_width * 6,
        shininess = 25,
        fog = TRUE
      )
      
      rgl::lines3d(
        distances_x,
        track_altitude / zscale,
        -distances_y,
        color = colour,
        alpha = 0.6,
        lwd = line_width * 3,
        shininess = 80,
        fog = TRUE
      )
      
      rgl::lines3d(
        distances_x,
        track_altitude / zscale,
        -distances_y,
        color = lighten(colour),
        alpha = 1,
        lwd = 1,
        shininess = 120
      )
      
    }
    
    if(ground_shadow){
      rgl::lines3d(
        distances_x,
        gps_ground_line / zscale + raise_agl,
        -distances_y,
        color = "black",
        alpha = 0.4,
        lwd = line_width * 2,
        shininess = 25,
        fog = TRUE
      )
    }
  } else {
    # switch shape point vs pch
    if(pch==16){
      rgl::points3d(
            distances_x,  #lat
            track_altitude / zscale + raise_agl,  #alt
            -distances_y,  #long
            color = colour,
            alpha = alpha,
            size = point_size
      )  
    } else {
      rgl::pch3d(pch=pch,
            distances_x,  #lat
            track_altitude / zscale + raise_agl,  #alt
            -distances_y,  #long
            color = colour,
            alpha = alpha,
            size = point_size
      )  
    }
  }
}

# based on add_gps_to_rayshader, but adds text
add_text_to_rayshader <- function(raster_input, lat, long, alt, zscale, text, colour = "red", alpha = 0.8,
                                  clamp_to_ground = FALSE, raise_agl = 0, ground_shadow = FALSE, char_size = 20, cex = rgl::par3d("cex")){
  
  coords <- latlong_to_rayshader_coords(raster_input, lat, long)
  
  distances_x <- coords$x
  
  distances_y <- coords$y
  
  
  if (clamp_to_ground | ground_shadow) {
    
    sp_gps <- sp::SpatialPoints(cbind(long, lat), proj4string = sp::CRS('+init=epsg:4326'))
    
    sp_gps <- sp::spTransform(sp_gps, sp::CRS(as.character(raster::crs(raster_input))))
    
    gps_ground_line <- raster::extract(raster_input, sp_gps)
    
  }
  
  if(clamp_to_ground){
    
    track_altitude <- gps_ground_line
    
  } else {
    
    track_altitude <- alt
  }
  
  rgl::text3d(
    distances_x,  #lat
    track_altitude / zscale + raise_agl,  #alt
    -distances_y,  #long
    text,
    cex = cex,
    color = colour,
    alpha = alpha,
    size = char_size
  )
}

# read in recorder map
gpspos = read.table("ZealandiaPointsOld.txt", sep=" ", nrows=10, h=F)
gpspos$V1 = paste0("Z", substr(gpspos$V1, 2, 2))
ggplot(gpspos, aes(y=V2, x=V3)) + geom_point() + geom_text(aes(label=V1), nudge_y=0.0005) +
  theme(aspect.ratio = 1)
coordinates(gpspos) = c("V3", "V2")
proj4string(gpspos) = CRS("+proj=longlat +datum=WGS84")

# read in satellite image
gpsmap = readGDAL(paste0(outdir, "BQ31_subset.tif"))

# get everything in both projections
gpsposUTM = spTransform(gpspos, CRS(proj4string(gpsmap)))
gpsmapWGS = spTransform(gpsmap, CRS(proj4string(gpspos)))

bbox(gpsmap)
bbox(gpsmapWGS)
bbox(gpspos)
bbox(gpsposUTM)


# download elevation data
mpapikey="pk.eyJ1Ijoiamp1b2QiLCJhIjoiY2p1N3N5YmY1MXh2aDQ1cDhjaXh0ODlzcCJ9.6PGlWLkFEusWqUkOu4a8sQ"

# part copied from slippy_raster, but with manual bounding box:
# download an elevation map for the entire gpsmap location
raster_elev <- get_slippy_map_15(gpsmapWGS, image_source="mapbox", image_type="terrain-rgb", max_tiles=30, api_key=mpapikey)
# delete transformed map, save RAM
rm(gpsmapWGS)
gc()

# project to UTM
raster_elev <- raster::projectRaster(raster_elev, crs = CRS(proj4string(gpsmap)))

# crop to plotting size
mapbbox = extent(c(1745150, 1746000, 5425000, 5426300))
#crs(mapbbox) = CRS(proj4string(gpsmap)) - doesn't work??
raster_elev2 = crop(raster_elev, mapbbox)

# template for crop/resample
overlay_raster <- raster(mapbbox, ncol=850, nrow=1300, crs = CRS(proj4string(gpsmap)))

# resample elevation
raster_elev3 <- raster::resample(raster_elev2, overlay_raster)

# convert/flatten elevation color layers somehow
DEM = -10000 + ((
    raster(raster_elev3, layer = 1) * 256 * 256 +
    raster(raster_elev3, layer = 2) * 256 +
    raster(raster_elev3, layer = 3))
  * 0.1)

# make mountains light
DEMsh = elevation_shade(DEM, elevation_palette=c("#000000", "#ffffff"), png_opacity = 0.9)

# prep for rayshader
elev_mat = matrix(extract(DEM, extent(DEM), method = 'bilinear'),
                  nrow = ncol(DEM), ncol = nrow(DEM))

# crop & resample gpsmap
gpsmapR = stack(paste0(outdir, "BQ31_subset.tif"))
gpsmapR = crop(gpsmapR, mapbbox)
gpsmapR2 <- raster::resample(gpsmapR, overlay_raster)


# weird conversion to PNG
temp_map_image <- tempfile(fileext = ".png")
slippymath::raster_to_png(gpsmapR2, temp_map_image)
gpsmapR2 <- png::readPNG(temp_map_image)
file.remove(temp_map_image)
alpha_layer <- matrix(0.8, nrow = dim(gpsmapR2)[1], ncol = dim(gpsmapR2)[2])
map_image = gpsmapR2
map_image[,,4] = alpha_layer

# produce the rayshader 3d plot
scene <- elev_mat %>%
  sphere_shade(sunangle=270) %>% 
  add_overlay(DEMsh) %>%
  add_overlay(map_image)

plot_3d(scene, elev_mat, solidcolor="grey40")
save_png(scene, "~/Documents/kiwis/scripts/zealandia/elevationmap.png")


# add markers for recorders
scene <- elev_mat %>%
  sphere_shade(sunangle=270) %>% 
  add_overlay(DEMsh) %>%
  add_overlay(map_image)
# remove ZC if using this for 2018 Zealandia data
gpspos_notc = gpspos[gpspos$V1!="ZC",]
plot_3d(scene, elev_mat, solidcolor="grey40", zoom=0.5)
add_gps_to_rayshader(DEM, gpspos_notc$V2, gpspos_notc$V3,
                     zscale=1, clamp_to_ground=T, raise_agl=2,
                     point_size=20, as_line=F, colour="white", alpha=0.9)
#add_text_to_rayshader(DEM, gpspos_notc$V2, gpspos_notc$V3, text=gpspos_notc$V1,
#                      zscale=1, clamp_to_ground=T, raise_agl=30,
#                      char_size=80, colour="white")
add_text_to_rayshader(DEM, gpspos_notc$V2, gpspos_notc$V3, text="x",
                      zscale=1, clamp_to_ground=T, raise_agl=8, alpha=1,
                      char_size=120, cex=1.8, colour="white")
